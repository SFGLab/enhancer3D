import argparse
import logging
import os

from chromatin_model import load_chromatin_model_ensemble_from_filesystem
from enhancer3d import ChromatinRegion, Enhancer3dProject, \
    load_enhancer_atlas_dataset_from_filesystem, load_gencode_annotation_dataset_from_filesystem, run_distance_calculation_for_region

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="enh3d")

    parser.add_argument("--project-name", type=str, required=True, default="default-project")
    parser.add_argument("--species", type=str, required=True, default="hg38")
    parser.add_argument("--cell-line", type=str, required=True, default="GM12878")

    parser.add_argument("--mapping-data-path", type=str, required=True, default="./data")
    parser.add_argument("--ensemble-data-path", type=str, required=True)

    parser.add_argument("--enhancer-atlas-dataset-name", type=str, required=True)
    parser.add_argument("--gencode-annotation-dataset-name", type=str, required=True)

    parser.add_argument("--reference-ensemble-name", type=str, required=True)
    parser.add_argument("--reference-ensemble-region", type=str, required=True)

    parser.add_argument("--modification-ensemble-name", type=str, required=True)
    parser.add_argument("--modification-ensemble-region", type=str, required=True)

    parser.add_argument("--output-path", type=str, required=True, default="./output")

    args = parser.parse_args()

    project = Enhancer3dProject(
        name=args.project_name,
        species=args.species,
        cell_line=args.cell_line,
        reference_ensemble_region=ChromatinRegion.from_string(args.reference_ensemble_region),
        modification_ensemble_region=ChromatinRegion.from_string(args.modification_ensemble_region)
    )

    mapping_data_path = args.mapping_data_path
    if not os.path.exists(mapping_data_path):
        raise FileNotFoundError(f"Mapping data path not found: {mapping_data_path}")

    enhancer_atlas_dataset = load_enhancer_atlas_dataset_from_filesystem(
        data_path=mapping_data_path,
        dataset_name=args.enhancer_atlas_dataset_name
    )

    gencode_annotation_dataset = load_gencode_annotation_dataset_from_filesystem(
        data_path=mapping_data_path,
        dataset_name=args.gencode_annotation_dataset_name
    )

    ensemble_data_path = args.ensemble_data_path
    if not os.path.exists(ensemble_data_path):
        raise FileNotFoundError(f"Reference ensemble data path not found: {ensemble_data_path}")

    reference_ensemble = load_chromatin_model_ensemble_from_filesystem(
        data_path=ensemble_data_path,
        model_name=args.reference_ensemble_name
    )

    modification_ensemble = load_chromatin_model_ensemble_from_filesystem(
        data_path=ensemble_data_path,
        model_name=args.modification_ensemble_name
    )

    reference_ensemble_region = ChromatinRegion.from_string(args.reference_ensemble_region)
    modification_ensemble_region = ChromatinRegion.from_string(args.modification_ensemble_region)

    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    logger.info("Running distance calculation for potential enhancer-gene pairs")
    distances_for_potential_enhancer_gene_pairs = run_distance_calculation_for_region(
        project=project,
        enhancer_atlas_dataset=enhancer_atlas_dataset,
        gencode_annotation_dataset=gencode_annotation_dataset,
        reference_ensemble=reference_ensemble,
        modification_ensemble=modification_ensemble
    )

    logger.info(f"Saving the results to {args.output_path}")
    distances_for_potential_enhancer_gene_pairs.to_csv(
        os.path.join(args.output_path, f"{args.project_name}_distances.csv"),
        index=False
    )
