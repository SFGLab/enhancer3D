# enhancer3D

Code repository for the **enhancer3D** resource - a database of 3D chromatin model ensembles and enhancer-promoter distance profiles for archaic (Neanderthal, Denisovan) and modern human genomes (GM12878, H1-ESC, HFFc6). The repository contains the service backend, analysis notebooks, and preprocessed datasets associated with a publication (see [References](#references)).

The live database is accessible at: **https://3dgnome.mini.pw.edu.pl/enhancer3D/**

---

## Research Context

Three-dimensional chromatin organisation shapes gene regulation by modulating spatial proximity between enhancers and promoters. **enhancer3D** provides two complementary analysis modules:

1. **Archaic vs modern comparison** - 3D models for genomic regions flanking archaic-specific structural variants (SVs) from Neanderthal and Denisovan genomes, enabling comparative E-P distance analysis relative to the modern human reference.
2. **Cross-cell-line comparison** - whole-chromosome 3D models for three modern human cell lines, integrated with RNA-seq, ChromHMM chromatin states, and EnhancerAtlas 2.0 annotations to study tissue-specific regulatory architecture.

Models are computed using the **cudaMMC** GPU-accelerated Monte Carlo polymer modelling engine within the **3D-GNOME 3.0** platform, from CTCF ChIA-PET interaction data provided by the 4D Nucleome consortium.

---

## Repository Structure

```
.
├── connector/          # Java/Spark MongoDB connector
├── data/               # Preprocessed datasets (Parquet/JSON)
├── infrastructure/     # Docker, Spark/Livy, Temporal, MongoDB configs
├── playground/         # Exploratory notebooks and generated figures
├── research/           # Core analysis notebooks (reproduce paper results)
├── src/                # Python backend source code
├── utils/              # Project-generation utilities
└── requirements.txt
```

### `connector/`
Java-based Spark connector for MongoDB, used in analytics jobs that query or write enhancer3D model data. Built with Gradle.

- `src/.../MongoConnector.java` - Spark connector implementation
- `src/.../MongoConnectorConfiguration.java` - Connection configuration

### `data/`
Preprocessed datasets used by analysis notebooks and backend services.

| Path | Contents |
|---|---|
| `chromatin_states/` | ChromHMM segmentations for GM12878, H1-ESC, HFFc6 (`.parquet`) |
| `deseq/` | Pairwise DESeq2 differential expression results between cell lines (`.parquet`) |
| `links/` | EnhancerAtlas 2.0 E-P link tables lifted to hg38, per cell line (`.parquet`) |
| `projects/` | JSON modeling project descriptors (8k regional + whole-chromosome ensembles) |

### `infrastructure/`
Docker-based deployment stack for the full enhancer3D backend.

- **Dockerfiles**: `app_api`, `app_calculator`, `app_repacker`, `app_servant`, `spark`, `livy`, `jupyter`
- **`docker-compose.yml`** / `docker-compose.local.yml` / `docker-compose.production.yml` - multi-container orchestration
- **`config/`** - Livy/Spark defaults, MongoDB init scripts, Temporal deployment YAML, Jupyter config

### `playground/`
Exploratory notebooks and intermediate outputs used during analysis development.

Key notebooks:
- `compare_2_cell_line_ep_largest_distances*.ipynb` - E-P distance distribution comparisons between cell lines
- `compare_2_refs_for_links*.ipynb` - Cross-reference and cross-model distance comparisons
- `deseq_exp_2_clean.ipynb`, `deseq_plots_exp1.ipynb` - Secondary DESeq2 visualisations
- `distance_flow_showcase.ipynb`, `model_flow_showcase.ipynb` - Distance-calculation and model workflow demonstrations
- `figs/` - Pre-generated figures (violin, volcano, enrichment, density plots) used in manuscript preparation

### `research/`
**Primary analysis notebooks - start here to reproduce the published results.**

#### `research/enhancer3d/`
HTTP request templates and JSON configs for triggering whole-chromosome modeling jobs via the enhancer3D API (one per cell line: GM12878, H1-ESC, HFFc6).

#### `research/genome_spatial_organization/`
Three-notebook pipeline corresponding to the *Genome Biology* brief report:

| Notebook | Purpose |
|---|---|
| `1_extract_closest_enh_distance_by_gene.ipynb` | Compute per-gene nearest-active-enhancer distances; assign proximity categories (small / mid / large) by chromosome percentile |
| `2_compute_rna_expression_difference_for_cell_lines.ipynb` | Compute pairwise differential expression between cell lines (uses `data/deseq/` or reruns PyDESeq2) |
| `3_compare_ep_distances_to_rna_expression.ipynb` | Integrate proximity changes with log fold changes; reproduce scale-dependent coupling results and pathway enrichment figures |

### `src/`
Python source code for the enhancer3D backend. Entry-point scripts in the root of `src/`:

- `app_api.py` - REST API server for queries and visualisation endpoints
- `app_calculator.py` - Worker executing distance-calculation workflows
- `app_repacker.py` - Worker repacking model outputs into Parquet datasets
- `app_servant.py` - Orchestrator for user-facing operations

Internal modules:

| Module | Role |
|---|---|
| `api/` | Request/response models for the REST API |
| `calculator/` | Temporal activities and workflows for E-P distance computation |
| `chromatin_model/` | Loaders for 3D-GNOME and packed model formats; model packing/unpacking |
| `common/` | Shared Pydantic data models |
| `database/` | Storage abstraction over MinIO (Parquet) and MongoDB |
| `distance_calculation/` | Core Euclidean distance computation and ensemble averaging |
| `repacker/` | Workflows for transforming raw model data into analysis-ready Parquet |
| `servant/` | Orchestration of archaic-modern and cross-cell-line comparison workflows |
| `utils/` | Helpers for filesystem, Mongo, Pandas, Pydantic, Scylla, and Temporal |

### `utils/`
- `produce_project_for_whole_chromosomal_models.py` - Generates project JSON descriptors for whole-chromosome modeling runs (populates `data/projects/`)

---

## Analysis Workflow

To **reproduce the published analyses** without redeploying the backend:

```
1. pip install -r requirements.txt
2. Launch Jupyter (locally or via infrastructure/jupyter.Dockerfile)
3. Run research/genome_spatial_organization/ notebooks in order: 1 -> 2 -> 3
   (precomputed inputs are in data/chromatin_states/, data/links/, data/deseq/)
4. Explore playground/ for supplementary figures and sensitivity analyses
```

To **deploy the full backend** (optional, for database operations and API):
```
cd infrastructure && docker compose -f docker-compose.yml -f docker-compose.local.yml up
```
This starts the application services, Apache Spark + Livy, MongoDB, Temporal workflow engine, and Jupyter. Remember to configure environment variables (e.g. MongoDB credentials) as needed, check .env.example for reference.

---

## Data

### Included in this repository
- `data/chromatin_states/` - ChromHMM annotations for three cell lines
- `data/links/` - E-P link tables (EnhancerAtlas 2.0, hg38)
- `data/deseq/` - Pairwise DESeq2 differential expression results
- `data/projects/` - Modeling project descriptors
- `playground/links/experiment_*/` - Experiment-specific intermediate link tables

### External (download separately)
Full 3D model ensembles, E-P distance tables, and annotation tracks are available at:
**https://3dgnome.mini.pw.edu.pl/download/enhancer3D**

---

## Requirements

```
pip install -r requirements.txt
```

The analysis notebooks additionally require a Jupyter environment. The full backend stack (Spark, Temporal, MongoDB) is defined in `infrastructure/`.

---

## References

1. Wlasnowolski M, Kozlov N, Wojcik M, Jacobs GS, Plewczynski D.
   **enhancer3D: 3D chromatin structures and enhancer-promoter distance profiles for archaic and modern human genomes.**
   *Nucleic Acids Research* 54(D1):D1046-D1052, 2026.
   DOI: [10.1093/nar/gkaf1256](https://doi.org/10.1093/nar/gkaf1256)
