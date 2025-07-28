import argparse
from typing import NamedTuple, List
import sys
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(script_dir, '../src')
sys.path.append(src_dir)

from common.models import (
    Enhancer3dProject, Enhancer3dProjectConfiguration, Enhancer3dProjectDataset,
    Enhancer3dEnhancerAtlasDatasetType, Enhancer3dGencodeAnnotationDatasetType, ChromatinRegion,
    Enhancer3dProjectDatasetMetadata
)

from calculator.models import CalculateDistancesForProjectWorkflowInput


class Args(NamedTuple):
    project_id: str
    project_authors: list[str]
    project_species: list[str]
    project_cell_line: str
    ensemble_id_format: str
    ensemble_id_chromosomes: list[int]
    enhancer_atlas_dataset_name: str
    enhancer_atlas_dataset_type: str
    gencode_annotation_dataset_name: str
    gencode_annotation_dataset_type: str
    base_pair_linear_distance_threshold: int | None
    enhancer_promoter_pairs_chunk_size: int
    regions_flank_by: int
    output_file: str
    input_segments_file: str


def parse_args() -> Args:
    parser = argparse.ArgumentParser(description="Produce project for whole chromosomal models.")
    parser.add_argument(
        "--project_id",
        type=str,
        required=True,
        help="Project ID for the Enhancer3D project. This should be a unique identifier.",
    )
    parser.add_argument(
        "--project_authors",
        type=str,
        nargs='*',
        default=[],
        help="List of authors for the project. Provide multiple names separated by spaces.",
    )
    parser.add_argument(
        "--project_species",
        type=str,
        nargs='*',
        default=[],
        help="List of species for the project. Provide multiple species names separated by spaces.",
    )
    parser.add_argument(
        "--project_cell_line",
        type=str,
        help="List of cell lines for the project. Provide multiple cell line names separated by spaces.",
    )

    parser.add_argument(
        "--ensemble_id_format",
        type=str,
        default="{chromosome}",
        help="Format string for the ensemble ID. Use {chromosome} to include the chromosome name.",
    )
    parser.add_argument(
        "--ensemble_id_chromosomes",
        type=int,
        nargs='*',
        default=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22],
        help="List of chromosomes to include in the ensemble IDs. Provide multiple chromosome names separated by spaces.",
    )

    parser.add_argument(
        "--enhancer_atlas_dataset_name",
        type=str,
        default="enhancer_atlas",
        help="Name of the Enhancer Atlas dataset.",
    )
    parser.add_argument(
        "--enhancer_atlas_dataset_type",
        type=str,
        choices=["bed", "tsv_liftovered_ref", "tsv_liftovered_mod"],
        default="bed",
        help="Type of the Enhancer Atlas dataset.",
    )
    parser.add_argument(
        "--gencode_annotation_dataset_name",
        type=str,
        default="gencode_annotation",
        help="Name of the Gencode annotation dataset.",
    )
    parser.add_argument(
        "--gencode_annotation_dataset_type",
        type=str,
        choices=["gtf", "tsv_liftovered_ref", "tsv_liftovered_mod"],
        default="gtf",
        help="Type of the Gencode annotation dataset.",
    )

    parser.add_argument(
        "--base_pair_linear_distance_threshold",
        type=int,
        default=None,
        help="Base pair linear distance threshold for the project configuration. If not specified, no threshold will be applied.",
    )
    parser.add_argument(
        "--enhancer_promoter_pairs_chunk_size",
        type=int,
        default=100000,
        help="Chunk size for processing enhancer-promoter pairs. Default is 100000.",
    )
    parser.add_argument(
        "--regions_flank_by",
        type=int,
        default=0,
        help="Number of base pairs to flank each region by. Default is 0, meaning no flanking.",
    )

    parser.add_argument(
        "--output_file",
        type=str,
        default="project_configuration.json",
        help="Path to the output file where the project configuration will be saved.",
    )

    parser.add_argument(
        "--input_segments_file",
        type=str,
        required=True,
        help="Path to the input segments file containing chromosomal coordinates.",
    )

    args = parser.parse_args()
    return Args(
        project_id=args.project_id,
        project_authors=args.project_authors,
        project_species=args.project_species,
        project_cell_line=args.project_cell_line,
        ensemble_id_format=args.ensemble_id_format,
        ensemble_id_chromosomes=args.ensemble_id_chromosomes,
        enhancer_atlas_dataset_name=args.enhancer_atlas_dataset_name,
        enhancer_atlas_dataset_type=args.enhancer_atlas_dataset_type,
        gencode_annotation_dataset_name=args.gencode_annotation_dataset_name,
        gencode_annotation_dataset_type=args.gencode_annotation_dataset_type,
        base_pair_linear_distance_threshold=args.base_pair_linear_distance_threshold,
        enhancer_promoter_pairs_chunk_size=args.enhancer_promoter_pairs_chunk_size,
        regions_flank_by=args.regions_flank_by,
        output_file=args.output_file,
        input_segments_file=args.input_segments_file
    )


# chr1	4932966	4932966
# chr1	7168152	7168152
# chr1	8803670	8803670
# chr1	9607446	9607446
# -> ChromatinRegion(chromosome=1, start=4932966, end=7168152),
#    ChromatinRegion(chromosome=1, start=7168152, end=8803670),
#    ...
def extract_chromosomal_region_from_file(file_path: str, chromosome: int, flank_by: int = 0) -> List[ChromatinRegion]:
    with open(file_path, 'r') as f:
        lines = [
            line.strip().split('\t')
            for line in f
            if line.strip() and not line.startswith('#')
        ]

    points_in_sequence_for_chromosome = [
        (chromosome, int(start), int(end))
        for chr, start, end in lines
        if chr == f'chr{chromosome}'
    ]

    points_in_sequence_for_chromosome = sorted(points_in_sequence_for_chromosome, key=lambda x: x[1])

    regions = []
    for chr, start, end in points_in_sequence_for_chromosome:
        start = max(0, start - flank_by)
        end = end + flank_by
        regions.append(ChromatinRegion(chromosome=f"chr{chr}", start=start, end=end))

    return regions


def main():
    args = parse_args()
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
            ensemble_id=args.ensemble_id_format.format(chromosome=chromosome),
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
        for chromosome in args.ensemble_id_chromosomes
        for region in extract_chromosomal_region_from_file(args.input_segments_file, chromosome, args.regions_flank_by)
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
