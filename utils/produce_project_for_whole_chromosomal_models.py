import argparse
import os
import sys
from typing import List, Tuple, Optional

import dotenv
from pydantic import BaseModel

script_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(script_dir, '../src')
sys.path.append(src_dir)

from common.models import (
    Enhancer3dProject, Enhancer3dProjectConfiguration, Enhancer3dProjectDataset,
    Enhancer3dEnhancerAtlasDatasetType, Enhancer3dGencodeAnnotationDatasetType, ChromatinRegion,
    Enhancer3dProjectDatasetMetadata
)

from calculator.models import CalculateDistancesForProjectWorkflowInput
from utils.filesystem_utils import get_bucket_filesystem
from chromatin_model.loaders.packed import load_chromatin_model_ensemble_head_from_filesystem
from chromatin_model.models import ChromatinModelEnsembleHead


class EnsembleArg(BaseModel):
    path: str
    chromosome: str


class Args(BaseModel):
    project_id: str
    project_authors: list[str]
    project_species: list[str]
    project_cell_line: str
    enhancer_atlas_dataset_name: str
    enhancer_atlas_dataset_type: str
    gencode_annotation_dataset_name: str
    gencode_annotation_dataset_type: str
    base_pair_linear_distance_threshold: Optional[int] = None
    enhancer_promoter_pairs_chunk_size: int = 100_000
    output_file: str
    input_segments_file: str
    ensembles: list[EnsembleArg]


def parse_args() -> Args:
    parser = argparse.ArgumentParser(description="Produce project for whole chromosomal models.")
    parser.add_argument(
        "--project_definition_file",
        type=str,
        required=True,
        help="Path to the project definition file (JSON format)."
    )

    args = parser.parse_args()
    project_definition_path = args.project_definition_file
    if not os.path.isfile(project_definition_path):
        raise FileNotFoundError(f"Project definition file not found: {project_definition_path}")

    with open(project_definition_path, 'r') as f:
        args = Args.model_validate_json(f.read())

    return args

# chr1	4932966	7168152
# chr1	8803670	9607446
# -> ChromatinRegion(chromosome=1, start=0, end=4932966),
#    ChromatinRegion(chromosome=1, start=4932966, end=7168152),
#    ChromatinRegion(chromosome=1, start=7168152, end=8803670),
#    ...
def extract_chromosomal_region_from_file(file_path: str, ensemble: EnsembleArg, head: ChromatinModelEnsembleHead) -> List[ChromatinRegion]:
    with open(file_path, 'r') as f:
        lines = [
            line.strip().split('\t')
            for line in f
            if line.strip() and not line.startswith('#')
        ]

    points_in_sequence_for_chromosome = [
        (ensemble.chromosome, int(start), int(end))
        for chr, start, end in lines
        if chr == ensemble.chromosome
    ]

    points_in_sequence_for_chromosome = sorted(points_in_sequence_for_chromosome, key=lambda x: x[1])

    regions = []
    for i, (chr, start, end) in enumerate(points_in_sequence_for_chromosome):
        if i == 0:
            if start > 0:
                regions.append(ChromatinRegion(chromosome=ensemble.chromosome, start=0, end=start))
            regions.append(ChromatinRegion(chromosome=ensemble.chromosome, start=start, end=end))
        else:
            prev_chr, prev_start, prev_end = points_in_sequence_for_chromosome[i - 1]
            if start > prev_end:
                regions.append(ChromatinRegion(chromosome=ensemble.chromosome, start=prev_end, end=start))
            regions.append(ChromatinRegion(chromosome=ensemble.chromosome, start=start, end=end))

        if i == len(points_in_sequence_for_chromosome) - 1:
            if end < head.last_bin:
                regions.append(ChromatinRegion(chromosome=ensemble.chromosome, start=end, end=head.last_bin))

    return regions


def main():
    dotenv.load_dotenv()
    args = parse_args()

    bucket_fs = get_bucket_filesystem()
    model_repository_bucket = os.environ.get("MODEL_REPOSITORY_BUCKET", "model-repository")
    model_heads: List[Tuple[EnsembleArg, ChromatinModelEnsembleHead]] = [
        (ensemble, load_chromatin_model_ensemble_head_from_filesystem(bucket_fs, model_repository_bucket, ensemble.path))
        for ensemble in args.ensembles
    ]

    project = Enhancer3dProject(
        id=args.project_id,
        authors=args.project_authors,
        species=args.project_species,
        cell_lines=[args.project_cell_line]
    )

    configuration = Enhancer3dProjectConfiguration(
        base_pair_linear_distance_threshold=args.base_pair_linear_distance_threshold,
        enhancer_promoter_pairs_chunk_size=args.enhancer_promoter_pairs_chunk_size
    )

    datasets = [
        Enhancer3dProjectDataset(
            ensemble_id=ensemble.path,
            ensemble_region=region,
            enhancer_atlas_dataset_name=args.enhancer_atlas_dataset_name,
            enhancer_atlas_dataset_type=Enhancer3dEnhancerAtlasDatasetType(args.enhancer_atlas_dataset_type),
            gencode_annotation_dataset_name=args.gencode_annotation_dataset_name,
            gencode_annotation_dataset_type=Enhancer3dGencodeAnnotationDatasetType(args.gencode_annotation_dataset_type),
            metadata=Enhancer3dProjectDatasetMetadata(
                cell_line=args.project_cell_line,
                species=args.project_species,
                authors=args.project_authors,
                modeling_pipeline="cudaMMC"
            )
        )
        for ensemble, head in model_heads
        for region in extract_chromosomal_region_from_file(args.input_segments_file, ensemble, head)
    ]

    input = CalculateDistancesForProjectWorkflowInput(
        project=project,
        datasets=datasets,
        configuration=configuration
    )

    with open(args.output_file, 'w') as f:
        f.write(input.model_dump_json(indent=4))


if __name__ == "__main__":
    main()
