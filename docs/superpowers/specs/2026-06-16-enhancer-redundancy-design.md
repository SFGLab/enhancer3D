# Enhancer Redundancy Experiment — Design

- **Date:** 2026-06-16
- **Status:** Implemented (notebooks `6_*` and `7_*`)
- **Location:** `research/genome_spatial_organization/` (new notebooks `6_*` and `7_*`)
- **Author:** Nik (with Claude)

> **Revision 2026-06-17 (final):** Enhancer sets **differ per cell line** (identical genes). §5 Step 2 is the **faithful 3D recompute** (the spec's heavy step): notebook 6 loads each c2 chromosome ensemble from `model-repository`, maps the gene TSS (strand-aware) and the c1-nearest enhancer's centre to model bins `(pos − first_bin)//resolution + 1`, and computes the ensemble-mean Euclidean distance → `dist_L1_in_c2` (near-full coverage). It runs **on the cluster**: the Livy driver lacks numpy/pandas/s3fs, so notebook 6 installs them into the live kernel (`install_deps` → `pip install -t SparkFiles.getRootDirectory()`), and writes a portable **CSV** per comparison (`spark.createDataFrame` is broken by a pyspark↔JVM version mismatch, and pickle isn't pandas-version-portable). `nearest_changed` (Cond 1) is overlap-based with `PADDING` (default 5 kb), since exact `enh_id` equality is ~always false across the different sets. Redundancy flags use a `REQUIRE_SMALL_IN_C1` toggle: `True` = the small→large reading (old enhancer was small/close in c1 AND large/far in c2); `False` = the literal Cond 2 (just large in c2). **Notebook 7 is local, storage-free** — it reads the CSVs + DESeq2 and reports frequency per comparison and the signed-log2FC effect (redundancy vs non-redundant switchers). *Superseded earlier detours: the local-numpy split, the shared-id lookup join, and the overlap-reuse proxy.* See `docs/superpowers/plans/2026-06-16-enhancer-redundancy.md`.

## 1. Motivation & hypothesis

We study cases where an **evolutionary/regulatory mechanism conserves a gene's transcription level**
despite spatial reorganization of its enhancer landscape between cell lines.

The pattern of interest ("enhancer redundancy / compensation"): the enhancer that was spatially
closest to a gene in cell line `c1` **moves away** (becomes spatially distant) in cell line `c2`,
but a **different** enhancer is now close to the gene in `c2`. If a different enhancer steps in to
keep an active enhancer nearby, expression should stay roughly constant — i.e. the gene's
`|log2FC|` should be **small** relative to genes that lost proximity without a replacement.

This extends the existing `1 → 3` distance-vs-expression pipeline by adding **enhancer identity**
and a **cross-cell-line 3D distance** dimension (the old enhancer's distance measured in the new
cell line's 3D model).

## 2. Inputs

| Input | Source (cluster) | Notes |
|---|---|---|
| Per-(gene, enhancer, ensemble) 3D distances | `s3a://database/results/{query_id}` per cell line (as in notebook 1), or master `s3a://database/distance_calculation` | Has `enh_id`, `enh_chr/start/end`, `gene_chr/start/end/strand`, `avg_dist`, `var_dist`, `region_id`, `ensemble_id`, `has_link`, `cell_line` |
| ChromHMM chromatin states | `s3a://database/chromatin_states` (cluster) / `data/chromatin_states/*.parquet` (local) | `chrom, start, end, name, …` |
| Project → ensemble mapping | `s3a://database/project_configuration` | `project_id → datasets[ensemble_id, ensemble_region, metadata.cell_line]` (see `src/servant/queries.py:53-108`) |
| 3D model ensembles (coordinates) | `model-repository` bucket, by `ensemble_id` | `{ensemble_id}.metadata.json` (+ `first_bin`, `last_bin`, `resolution`) and `{ensemble_id}.coordinates.npy`. Loaded via `load_chromatin_model_ensemble_from_filesystem` (`src/chromatin_model/loaders/packed.py:9`; usage `src/calculator/activities.py:95-99`, buckets at `:61-62`) |
| DESeq2 pairwise results | `data/deseq/{gm12878_vs_h1esc,h1esc_vs_hffc6,gm12878_vs_hffc6}_results.parquet` | `log2FoldChange`, `padj`, `baseMean`; consumed **locally** in notebook 7 |

**Cell lines:** GM12878, H1ESC, HFFC6. **Project filter** (modern human, whole genome):
`whole_all_vs_all_{gm12878,h1esc,hffc6}_fix` (as in notebook 1).

### Input to confirm before implementation
- **Query IDs** for the three cell lines' `results/{query_id}` parquets (notebook 1 uses fixed UUIDs;
  this notebook needs results that retain `dist`/coords columns — confirm the right query IDs or read
  `s3a://database/distance_calculation` directly).
- **`model-repository` bucket** is reachable from the notebook's Spark/driver environment (needed for
  Step 2's recompute). Confirm `MODEL_REPOSITORY_BUCKET` and that the `full_2` per-chromosome ensembles
  for all three cell lines exist there.

## 3. Locked decisions

1. **Compute backend:** Spark on the cluster (mirrors notebook 1). Heavy compute writes a parquet to
   `s3a://database/enhancer_redundancy/…`; local analysis/plots read it back. *(User decision.)*
2. **"Active enhancer" states (broad set, matches notebook 1):**
   `TssA, TssAFlnk, TxFlnk, Tx, TxWk, EnhG, EnhG1, EnhG2, Enh, EnhA1, EnhA2`. An enhancer is *active*
   in a cell line when its locus overlaps a segment in this set. *(User decision.)*
3. **`proximity_category` (small / mid / large):** per-chromosome tertiles of the gene→nearest-active-enhancer
   distance distribution within each cell line (small ≤ 33rd pct, large > 67th pct), reusing notebook 3's
   `add_chromosome_and_quartiles` scheme. The old enhancer's recomputed distance in `c2` is classified against
   **`c2`'s own per-chromosome thresholds**, so `large` there is directly comparable to `proximity_category_c2`.
   *(User decision.)*
4. **"Nearest enhancer changed":** exact `enh_id` string inequality (`nearest_enh_c1 != nearest_enh_c2`).
   *(User decision.)* **Implication (documented, accepted):** because each cell line uses its own
   EnhancerAtlas set and enhancer loci overlap only ~2.5% across cell lines, this condition is satisfied
   for nearly all genes and therefore does **not** materially filter; the two proximity conditions
   (Step 3 below) do the discriminating. The denominator "genes where nearest changed" is thus effectively
   "genes with a defined nearest active enhancer in both `c1` and `c2`".
5. **Expression is the measured OUTCOME, not a filter.** Redundancy cases are defined purely by the three
   enhancer/proximity criteria; we then compare `|log2FC|` of **redundancy genes vs non-redundant switchers**
   (nearest changed + old→`large` in `c2`, but the new nearest enhancer is **not** `small`, i.e. lost proximity
   with no replacement). *(User decision.)*

## 4. Deliverables & structure

Two notebooks, mirroring the existing split (`1`: Spark→S3 heavy compute; `3`: local analysis+plots):

- **`6_extract_enhancer_redundancy.ipynb`** (Spark, cluster) — builds the redundancy table, writes parquet.
- **`7_analyze_enhancer_redundancy.ipynb`** (local) — frequency per comparison + `|log2FC|` comparison + figures.

Comparisons are **directional**: `c1 → c2` asks "did the `c1`-nearest enhancer move away in `c2`?".
The 3 unordered pairs give **6 directed comparisons**. DESeq2 `log2FC` sign is aligned to each contrast's
direction exactly as notebook 3 already flips `h1esc_vs_hffc6`.

## 5. Notebook 6 — pipeline (Spark)

### Step 1 — Nearest active enhancer per gene, per cell line
1. Read per-cell-line distances; filter `project_id ∈ used_projects`, `avg_dist > 0 AND var_dist > 0` (as nb1).
2. Keep `(gene, enh)` pairs where the **enhancer locus** overlaps an active-state segment (decision 2), via the
   EXISTS-join on `chromatin_states` used in notebook 1 (enhancer-side condition; the gene-side active condition
   from nb1 is retained for consistency).
3. Per `(cell_line, gene_id)`, pick `argmin(avg_dist)` with a window function
   (`row_number() over (partition by cell_line, gene_id order by avg_dist asc)` → keep rank 1).
   Output per gene: `nearest_enh_id`, `nearest_avg_dist`, plus the enhancer/gene coordinates and `gene_strand`
   (needed for Step 2 bin mapping).
4. Assign `proximity_category` = per-chromosome tertile of `nearest_avg_dist` within the cell line (decision 3).

Produces `nearest_c1` and `nearest_c2` tables keyed by `gene_id`.

### Step 2 — Old enhancer's distance in c2's model (the heavy step)
For each gene, `L1 = nearest_enh_c1` is a genomic locus (`enh_id = chr:start-end`). Compute its 3D distance to
the gene **in `c2`'s model**, regardless of whether `L1` is active (or present) in `c2`:

1. Resolve the `c2` ensemble for the gene's chromosome from `project_configuration` (the `*_full_2_chr{N}` ensemble);
   load it from `model-repository` (`load_chromatin_model_ensemble_from_filesystem`). Cache one ensemble per
   chromosome (≈21 per cell line).
2. Map positions to model bins using `region_start = ensemble.first_bin` and `ensemble.resolution`
   (formula from `src/distance_calculation/services.py:119-137`):
   - `gene_bin   = (gene_TSS − first_bin) // resolution + 1`, where `gene_TSS = gene_start` if strand `+` else `gene_end`.
   - `enh_bin    = (L1_center − first_bin) // resolution + 1`, where `L1_center = (enh_start + enh_end)//2`.
3. `dist_L1_in_c2 = ensemble.distance_distribution(gene_bin, enh_bin).mean()`
   (`src/chromatin_model/models.py:80-83`; same per-pair machinery as
   `calculate_distances_for_potential_enhancer_gene_pairs`, `services.py:185-217`).
4. Classify `dist_L1_in_c2` against **`c2`'s per-chromosome tertile thresholds** → `proximity_of_L1_in_c2`.
5. Annotate `L1_active_in_c2` (does `L1` overlap an active state in `c2`) for interpretation.

**Edge cases / flags:**
- `L1_center` outside `[first_bin, last_bin]` of the `c2` ensemble → cannot place in the model →
  flag `L1_out_of_c2_window = True` and treat `proximity_of_L1_in_c2 = large` (it is at least as far as the
  window edge; recorded so the analysis can include or exclude these).
- Genes with no defined nearest in either cell line are dropped (no comparison possible).

**Execution shape:** group target `(gene, L1)` rows by `c2` chromosome ensemble so each ensemble loads once;
compute all that chromosome's pairs vectorized over the model stack. This is the costliest stage (raw model
loading + per-pair distance), bounded to the gene set with a defined nearest — far smaller than all-pairs.

> **Optional fast-path / cross-check (not primary):** where some existing `c2` candidate enhancer already
> occupies `enh_bin` for the same gene in `distance_calculation`, its `avg_dist` equals `dist_L1_in_c2` exactly
> (identical bin ⇒ identical coordinates) and can be read by a pure Spark join, skipping the ensemble load for
> those rows. Useful to validate the recompute; coverage is partial so it does not replace it.

### Step 3 — Redundancy classification
Per gene, per directed comparison:
- `nearest_changed = (nearest_enh_c1 != nearest_enh_c2)`  *(exact id; decision 4)*
- `is_redundancy   = nearest_changed AND (proximity_of_L1_in_c2 == 'large') AND (proximity_category_c2 == 'small')`
- `is_nonredundant_switcher = nearest_changed AND (proximity_of_L1_in_c2 == 'large') AND (proximity_category_c2 != 'small')`

### Step 4 — Frequency
Per directed comparison: `frequency = |is_redundancy| / |nearest_changed|`.

### Output table (parquet → `s3a://database/enhancer_redundancy/{comparison}`; pulled to `data/whole_chromosomes/`)
`comparison, c1, c2, gene_id, gene_chr,`
`nearest_enh_c1, avg_dist_c1, proximity_category_c1,`
`nearest_enh_c2, avg_dist_c2, proximity_category_c2,`
`dist_L1_in_c2, proximity_of_L1_in_c2, L1_active_in_c2, L1_out_of_c2_window,`
`nearest_changed, is_redundancy, is_nonredundant_switcher`

## 6. Notebook 7 — outcome analysis (local)

1. Read the redundancy table(s) + `data/deseq/*`; join on `gene_id`; align `log2FC` direction per comparison;
   `abs_log2fc = |log2FoldChange|`.
2. **Frequency per comparison** (recompute/verify Step 4) → table + bar chart → `export/`.
3. **`|log2FC|` comparison:** redundancy genes vs **non-redundant switchers**:
   - Mann-Whitney U (two-sided), median `|log2FC|` per group, n per group.
   - Distribution plots: ECDF + violin/box; joyplot in notebook-3 style.
   - Expected: redundancy genes have **smaller** `|log2FC|`.
4. Optional sanity baseline: also show stable-nearest genes' `|log2FC|` distribution.
5. Exports: figures → `figs/exp6_*.png`, tables → `export/enhancer_redundancy_*.csv`.

## 7. Reused code / precedents
- Active-state filter, project filter, per-gene nearest pattern — `1_extract_closest_enh_distance_by_gene.ipynb`.
- Per-chromosome tertiles, `small→large` change detection, `|log2FC|` joyplot — `3_compare_ep_distances_to_rna_expression.ipynb`.
- Bin mapping + TSS strand logic — `src/distance_calculation/services.py:119-137`.
- 3D distance — `src/chromatin_model/models.py:80-83`; per-pair driver `services.py:185-217`.
- Ensemble loading — `src/chromatin_model/loaders/packed.py:9`; bucket usage `src/calculator/activities.py:61-99`.
- Project/ensemble regions — `s3a://database/project_configuration` via `src/servant/queries.py:53-108`.

## 8. Out of scope (YAGNI)
- No changes to the Temporal/`distance_calculation` production pipeline; the notebook reuses its primitives read-only.
- No new ensembles computed; we only read existing `full_2` models.
- No genome-wide recompute of all enhancer–gene distances — Step 2 is restricted to each gene's single `c1`-nearest enhancer.
- WTC11 (present in `data/projects/`) is excluded; only GM12878/H1ESC/HFFC6 have the DESeq2 + chromatin inputs.

## 9. Risks / open notes
- **Decision 4 implication** (Section 3.4): exact-id "changed" is near-tautological; frequency is driven by the
  proximity criteria. Documented and accepted.
- **Model availability/scale:** per-chromosome ensembles may be large (model_count × bins × 3); confirm driver/executor
  memory. Cache per chromosome and release.
- **Window boundaries:** `L1` just outside the `c2` ensemble window is flagged (`L1_out_of_c2_window`) and counted as
  `large`; report counts so the choice is auditable.
- **Multiple-testing / effect size:** report effect sizes (median diff, rank-biserial) alongside p-values, not p alone.
