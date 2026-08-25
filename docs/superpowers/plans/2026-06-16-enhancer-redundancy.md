# Enhancer Redundancy Experiment — Implementation Notes

**Goal:** Detect "enhancer redundancy" — a gene's nearest active enhancer changes between cell lines, the old one moves spatially far in c2, a different one is close — and measure whether expression is conserved.

**Spec:** `docs/superpowers/specs/2026-06-16-enhancer-redundancy-design.md` (see the 2026-06-17 revision).

## Design (final — faithful 3D recompute)

| File | Runs where | What it does |
|---|---|---|
| `research/genome_spatial_organization/6_extract_enhancer_redundancy.ipynb` | Spark / cluster | `install_deps(['pandas','numpy','s3fs'])` into the live Livy kernel. **Step 1 (Spark):** nearest active enhancer per gene per cell line (+ gene & enhancer coords). **Step 2 (driver, numpy):** for each (c2, chromosome) load the ensemble from `model-repository` once, map gene-TSS (strand-aware) + c1-nearest-enhancer centre to bins `(pos−first_bin)//res+1`, compute ensemble-mean ‖·‖ → `dist_L1_in_c2`. **Classify (pandas):** per-chromosome tertiles; `nearest_changed`=L1/L2 don't overlap within `PADDING`; redundancy flags (`REQUIRE_SMALL_IN_C1` toggle). Writes `s3://database/enhancer_redundancy/<comparison>/data.csv`. |
| `research/genome_spatial_organization/7_analyze_enhancer_redundancy.ipynb` | Local (`enhancer3D`) | **No storage/model access.** Reads the CSVs + DESeq2, computes frequency per comparison and the **signed log2FC** comparison (redundancy vs non-redundant switchers), both `all` and `measured_only`. CSVs → `export/`, figures → `figs/exp6_*.png`. |

**Spec ↔ implementation deviations (intentional):**
- **Cond 1** (`nearest_enh_c1 != nearest_enh_c2`) → overlap-based (`PADDING`, 5 kb), because enhancer sets differ per cell line so exact `enh_id` `!=` is ~always true.
- **Cond 2** (`c1_nearest region in c2 == large`) → exact, via the 3D recompute. `REQUIRE_SMALL_IN_C1=True` additionally requires it was *small in c1* (the "moved away"/small→large reading); set `False` for the literal spec.
- log2FC comparison contrast = **non-redundant switchers** (old→large in c2, no close replacement).

## How to run

1. **Notebook 6 (cluster) — full rerun, this is the heavy step.** Confirm: `QUERY_IDS` resolve to results parquets with `gene_*`/`enh_*`/coords; `model-repository` holds `<ensemble_id>.coordinates.npy` + `.metadata.json` for all 3 cell lines (`MODEL_REPOSITORY` constant); `s3fs` reaches it (config pulled from the Spark Hadoop conf). Set `REQUIRE_SMALL_IN_C1` as desired. Run all cells → per-comparison `in-window coverage` (should be ≈1.0) + frequencies, then CSVs to S3.
2. **Pull** `s3://database/enhancer_redundancy/*` → `data/whole_chromosomes/enhancer_redundancy/<comparison>/data.csv`.
3. **Notebook 7 (local).** Run → `export/enhancer_redundancy_frequency.csv`, `export/enhancer_redundancy_log2fc_stats.csv`, `figs/exp6_*.png`.

## Environment gotchas (hit & worked around)

- **No numpy/pandas/s3fs on the Livy driver** → `install_deps` (`pip install -t SparkFiles.getRootDirectory()` + add to `sys.path`).
- **`spark.createDataFrame(pandas)` is broken** (pyspark↔JVM mismatch: `PythonSQLUtils` is a JavaPackage) → write with pandas `to_csv(..., storage_options=...)` straight to S3.
- **pyarrow/parquet wedged** on the cluster (duplicate `pandas.period` registration) and **pickle isn't pandas-version-portable** → output is **CSV**; notebook 7 coerces the bool columns back.

## Verified locally (no cluster)

- Both notebooks: JSON valid, all non-magic cells compile.
- Step 2 model-distance: strand-aware TSS→bin + L1-centre→bin + ensemble ‖·‖ checked on a synthetic ensemble (g1=7, g2=8). ✓
- Notebook 7 DESeq2 join + signed-log2FC stats against real `data/deseq/*`. ✓

## Open caveats

- **Driver memory:** Step 2 holds one chromosome ensemble at a time and releases it; if the `whole_all_vs_all_*` models are very large (high model count), watch driver RAM (the per-(c2,chrom) loop already minimises it).
- **`PADDING`** only affects Cond 1 now (the distance is exact); sensitivity-check it.
- Two `export/split_celllines_reactome_terms*.txt` files still carry `git stash pop` conflict markers (notebook-5 outputs, untouched).
