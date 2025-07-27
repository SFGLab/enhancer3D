import numpy as np
from pyspark.sql import SparkSession
from pyspark.sql import functions as F, types as T
from pyspark.sql.connect.dataframe import DataFrame
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests

from servant.models import CellLineWithLinksInput, Response


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


def _projects_by_cell_line(spark: SparkSession, cell_line: str) -> DataFrame:
    return (
        spark
        .read
        .json("s3a://database/project_configuration", multiLine=True)
        .select(
            F.col('project_id'),
            F.explode('datasets.ensemble_id').alias('ensemble_id'),
            F.explode('datasets.metadata.cell_line').alias('cell_line')
        )
        .where(F.col('cell_line') == cell_line)
    )


def _distances_with_links_for_ensemble_ids(spark: SparkSession, relevant_projects_df: DataFrame) -> DataFrame:
    links_df = (
        spark
        .read
        .parquet("s3a://database/links")
        .alias("links")
    )

    distances_df = (
        spark
        .read
        .parquet("s3a://database/distance_calculation")
        .alias("distances")
        .join(
            other=relevant_projects_df.alias("relevant_projects"),
            on=F.expr(
                "distances.project_id = relevant_projects.project_id"
                " AND distances.ensemble_id = relevant_projects.ensemble_id"
            ),
            how="inner"
        )
        .select(
            F.col('relevant_projects.project_id').alias('project_id'),
            F.col('relevant_projects.ensemble_id').alias('ensemble_id'),
            'region_id',
            'enh_id',
            'gene_type',
            'avg_dist',
            'var_dist',
            'dist',
            'enh_tSS_distance',
            F.element_at(F.col('project_cell_lines'), 1).alias('cell_line'),
            F.split(F.col('gene_id'), '\.')[0].alias('gene_id'),  # gene_id ENH00001.XXX -> ENH00001
        )
        .where(
            (F.col('gene_type') == 'protein_coding')
            # & (F.col('enh_tSS_distance') < 20_000)
        )
        .alias('distances')
    )

    distances_with_links_df = (
        distances_df
        .join(
            other=links_df,
            on=F.expr(
                "distances.cell_line = links.cell_line"
                " AND distances.gene_id = links.gene_id"
                " AND distances.enh_id = links.enh_id"
            ),
            how="outer"
        )
        .select(
            F.col('distances.project_id'),
            F.col('distances.ensemble_id'),
            F.col('distances.region_id'),
            F.col('distances.enh_id'),
            F.col('distances.gene_type'),
            F.col('distances.avg_dist'),
            F.col('distances.var_dist'),
            F.col('distances.dist'),
            F.col('distances.enh_tSS_distance'),
            # If has link then True else False
            F.when(F.col('links.gene_id').isNotNull(), True).otherwise(False).alias('has_link'),
        )
    )

    return distances_with_links_df


def query_cell_line_with_links(spark: SparkSession, request_id: str, input: CellLineWithLinksInput) -> Response:
    relevant_projects_df = _projects_by_cell_line(spark, input.cell_line)
    distances_with_links_df = _distances_with_links_for_ensemble_ids(spark, relevant_projects_df)

    path = f"s3a://database/results/{request_id}.parquet"
    distances_with_links_df.write.mode("overwrite").parquet(path)

    return Response(
        request_id=request_id,
        path=path
    )
