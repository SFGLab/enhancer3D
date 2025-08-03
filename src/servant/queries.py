from typing import List

import numpy as np
from pyspark.sql import SparkSession
from pyspark.sql import functions as F, types as T
from pyspark.sql.connect.dataframe import DataFrame
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests

from servant.models import CellLineWithLinksInput, Response, PartialChromatinRegion


@F.udf(T.ArrayType(T.DoubleType()))
def diff(A, B):
    return np.abs(np.array(A) - np.array(B)).tolist()


@F.udf(T.DoubleType())
def var(A):
    return float(np.var(A))


@F.udf(T.DoubleType())
def avg(A):
    return float(np.mean(A))


@F.udf(T.DoubleType())
def mannwhiteneyu(ref, mod):
    result = stats.mannwhitneyu(np.array(ref), np.array(mod), alternative='two-sided')
    return float(result.pvalue)


@F.udf(T.DoubleType())
def bonferroni_correction(pvalues, alpha=0.05):
    reject, pvals_corrected, _, _ = multipletests(pvalues, alpha=alpha, method='bonferroni')
    return float(np.mean(pvals_corrected))


def _projects_by_cell_line(
    spark: SparkSession,
    cell_line: str,
    species: str = None,
    regions: List[PartialChromatinRegion] = None
) -> DataFrame:
    projects_df = (
        spark
        .read
        .json("s3a://database/project_configuration", multiLine=True)
        .select(
            F.col('project_id'),
            F.explode('datasets').alias('datasets')
        )
        .select(
            F.col('project_id'),
            F.col('datasets.ensemble_id').alias('ensemble_id'),
            F.col('datasets.ensemble_region').alias('ensemble_region'),
            F.col('datasets.metadata.cell_line').alias('cell_line'),
            F.col('datasets.metadata.species').alias('species'),
        )
        .where(F.col('cell_line') == cell_line)
    )

    if species:
        projects_df = projects_df.where(F.array_contains(F.col('species'), species))

    if regions and len(regions) > 0:
        regions_schema = T.StructType([
            T.StructField("chromosome", T.StringType(), False),
            T.StructField("start", T.IntegerType(), True),
            T.StructField("end", T.IntegerType(), True)
        ])

        regions_df = spark.createDataFrame(
            [region.model_dump() for region in regions],
            schema=regions_schema
        )

        projects_df = (
            projects_df
            .join(
                F.broadcast(regions_df),
                (projects_df.ensemble_region.chromosome == regions_df.chromosome) &
                (regions_df.start.isNull() | (projects_df.ensemble_region.start <= regions_df.end)) &
                (regions_df.end.isNull() | (projects_df.ensemble_region.end >= regions_df.start)),
                "inner"
            )
            .drop("chromosome", "start", "end")
        )

    return projects_df.alias('projects')


def _distances_for_relevant_projects(
    spark: SparkSession,
    relevant_projects_df: DataFrame,
    regions: List[PartialChromatinRegion],
    gene_ids: List[str]
) -> DataFrame:
    distances_df = (
        spark
        .read
        .parquet("s3a://database/distance_calculation")
        .alias("distances")
        .join(
            relevant_projects_df.select('project_id', 'ensemble_id'),
            ['project_id', 'ensemble_id'],
            how='semi'
        )
        .select(
            'project_id',
            'ensemble_id',
            'region_id',
            'enh_id',
            'enh_chr',
            'enh_start',
            'enh_end',
            'enh_score',
            'avg_dist',
            'var_dist',
            'dist',
            'enh_tSS_distance',
            F.element_at(F.col('project_cell_lines'), 1).alias('cell_line'),
            F.split(F.col('gene_id'), '\.')[0].alias('gene_id'),  # gene_id ENH00001.XXX -> ENH00001
            'gene_chr',
            'gene_start',
            'gene_end',
            'gene_strand',
            'gene_type',
        )
        .where(
            (F.col('gene_type') == 'protein_coding')
            # & (F.col('enh_tSS_distance') < 20_000)
        )
    )

    if gene_ids and len(gene_ids) == 1:
        distances_df = distances_df.where(F.col('gene_id').like(f"%{gene_ids[0]}%"))

    if gene_ids and len(gene_ids) > 1:
        distances_df = distances_df.where(F.col('gene_id').isin(gene_ids))

    if regions and len(regions) > 0:
        regions_schema = T.StructType([
            T.StructField("chromosome", T.StringType(), False),
            T.StructField("start", T.IntegerType(), True),
            T.StructField("end", T.IntegerType(), True)
        ])

        regions_df = spark.createDataFrame(
            [region.model_dump() for region in regions],
            schema=regions_schema
        )

        distances_df = (
            distances_df
            .join(
                F.broadcast(regions_df),
                (distances_df.enh_chr == regions_df.chromosome) &
                (regions_df.start.isNull() | (distances_df.enh_start <= regions_df.end)) &
                (regions_df.end.isNull() | (distances_df.enh_end >= regions_df.start)) &
                (distances_df.gene_chr == regions_df.chromosome) &
                (regions_df.start.isNull() | (distances_df.gene_start <= regions_df.end)) &
                (regions_df.end.isNull() | (distances_df.gene_end >= regions_df.start)),
                "inner"
            )
            .drop("chromosome", "start", "end")
        )

    return distances_df.alias('distances')


def _distances_with_links_for_ensemble_ids(
    spark: SparkSession,
    relevant_projects_df: DataFrame,
    regions: List[PartialChromatinRegion] = None,
    gene_ids: List[str] = None
) -> DataFrame:
    links_df = (
        spark
        .read
        .parquet("s3a://database/links")
        .alias("links")
    )

    distances_df = _distances_for_relevant_projects(spark, relevant_projects_df, regions, gene_ids)
    distances_with_links_df = (
        distances_df
        .repartition(64, 'cell_line', 'gene_id', 'enh_id')
        .join(
            other=links_df.repartition(64, 'cell_line', 'gene_id', 'enh_id'),
            on=["cell_line", "gene_id", "enh_id"],
            how="outer"
        )
        .select(
            F.col('distances.project_id'),
            F.col('distances.ensemble_id'),
            F.col('distances.region_id'),
            F.col('distances.enh_id'),
            F.col('distances.enh_chr'),
            F.col('distances.enh_start'),
            F.col('distances.enh_end'),
            F.col('distances.enh_score'),
            F.col('distances.gene_id'),
            F.col('distances.gene_type'),
            F.col('distances.gene_chr'),
            F.col('distances.gene_start'),
            F.col('distances.gene_end'),
            F.col('distances.gene_strand'),
            F.col('distances.avg_dist'),
            F.col('distances.var_dist'),
            F.col('distances.dist'),
            F.col('distances.enh_tSS_distance'),
            F.col('distances.cell_line'),
            # If has link then True else False
            F.when(F.col('links.gene_id').isNotNull(), True).otherwise(False).alias('has_link'),
        )
    )

    return distances_with_links_df


def query_cell_line_with_links(spark: SparkSession, request_id: str, input: CellLineWithLinksInput) -> Response:
    relevant_projects_df = _projects_by_cell_line(spark, input.cell_line, input.species)
    distances_with_links_df = _distances_with_links_for_ensemble_ids(spark, relevant_projects_df, input.regions, input.gene_ids)

    path = f"s3a://database/results/{request_id}.parquet"
    distances_with_links_df.repartition(1).write.mode("overwrite").parquet(path)

    return Response(
        request_id=request_id,
        path=path
    )


def query_cell_link_cross_comparison(spark: SparkSession, request_id: str, input: CellLineWithLinksInput) -> Response:
    relevant_projects_base_df = _projects_by_cell_line(spark, input.cell_line_base, input.regions, input.species_base)
    relevant_projects_compare_df = _projects_by_cell_line(spark, input.cell_line_compare, input.regions, input.species_compare)

    distances_with_links_base_df = (
        _distances_with_links_for_ensemble_ids(spark, relevant_projects_base_df, input.regions, input.gene_ids)
        .repartition(64, 'gene_id', 'enh_id')
        .alias('distances_base')
    )
    distances_with_links_compare_df = (
        _distances_with_links_for_ensemble_ids(spark, relevant_projects_compare_df, input.regions, input.gene_ids)
        .repartition(64, 'gene_id', 'enh_id')
        .alias('distances_compare')
    )

    cross_comparison_df = (
        distances_with_links_base_df
        .join(
            distances_with_links_compare_df,
            on=['gene_id', 'enh_id'],
            how='inner'
        )
        .select(
            'gene_id',
            'enh_id',
            F.col('distances_base.region_id').alias('region_id_base'),
            F.col('distances_compare.region_id').alias('region_id_compare'),
            F.col('distances_base.gene_chr').alias('gene_chr'),
            F.col('distances_base.gene_start').alias('gene_start'),
            F.col('distances_base.gene_end').alias('gene_end'),
            F.col('distances_base.gene_strand').alias('gene_strand'),
            F.col('distances_base.gene_type').alias('gene_type'),
            F.col('distances_base.enh_chr').alias('enh_chr'),
            F.col('distances_base.enh_start').alias('enh_start'),
            F.col('distances_base.enh_end').alias('enh_end'),
            F.col('distances_base.enh_score').alias('enh_score'),
            diff(F.col('distances_base.dist'), F.col('distances_compare.dist')).alias('diff_dist'),
            var(F.col('distances_base.dist')).alias('var_dist_base'),
            var(F.col('distances_compare.dist')).alias('var_dist_compare'),
            avg(F.col('distances_base.dist')).alias('avg_dist_base'),
            avg(F.col('distances_compare.dist')).alias('avg_dist_compare'),
            mannwhiteneyu(F.col('distances_base.dist'), F.col('distances_compare.dist')).alias('mannwhiteneyu_pvalue')
        )
    )

    path = f"s3a://database/results/{request_id}.parquet"
    cross_comparison_df.repartition(1).write.mode("overwrite").parquet(path)

    return Response(
        request_id=request_id,
        path=path
    )
