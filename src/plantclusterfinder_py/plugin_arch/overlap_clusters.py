"""Built-in plugin for computing enzyme statistics.

This plugin demonstrates how to consume PlantClusterFinder results
and produce a derived analysis without modifying core pipeline behavior.
"""

import logging

import polars as pl

from plantclusterfinder_py.plugin_arch.context import ReturnContext, RunContext
from plantclusterfinder_py.plugin_arch.hookspecs import hookimpl

_log = logging.getLogger(__name__)


class OverlapClustersPlugin:
    """A hook specification namespace."""

    expeced_ouptuts = ["cluster_filename"]

    @hookimpl
    def process_results(self, results: pl.DataFrame, context: RunContext) -> ReturnContext:
        """Example plugin logic: process DataFrame to get non-null genes."""
        if "cluster_filename" not in context.output_dict:
            raise ValueError("Setting 'cluster_filename' not in plugin output dictionary")
        _log.debug("Cluster Filename is %s", context.output_dict.get("cluster_filename", ""))
        processed_df = pl.read_csv(
            context.output_dir / context.output_dict["cluster_filename"], separator="\t"
        )
        output_filename = context.output_dir / f"{context.run_name}.non-overlapping_clusters.csv"
        _log.info("Writing output to %s", output_filename)
        processed_df.write_csv(output_filename)
        return ReturnContext(plugin_name="visualization_plugin", output_filename=output_filename)

    def get_no_overlap_clusts(self, clusters: pl.DataFrame) -> pl.DataFrame:
        """Get a list of clusters which do not overlap by combining overlapping clusters.

        Args:
            clusters: The cluster output from PlantClusterFinder.py

        Returns:
            A DataFrame with no overlapping clusters.
                ClusterID - unique ID for cluster (duplicate ID from clusters DataFrame)
                Chromosome - Chromosome ID Cluster is on
                Start bp - Start bp for cluster
                End bp - End bp for cluster
                Gene IDs - list of Gene IDs in cluster
        """
        mclust_expr = (
            pl.col("ClusterID3"),
            # pl.col("Stepsize").first(),
            pl.col("Chromosome").first(),
            pl.col("Start bp").min(),
            pl.col("End bp").max(),
            # union of all genes in mega-cluster
            pl.reduce(function=lambda acc, x: pl.union((acc, x)), exprs=pl.col("Gene IDs")).last(),
        )
        cdf = clusters.with_columns(pl.col("Gene IDs").str.split(","))
        # Get all clusters with overlapping genes to another cluster
        overlapping_clusts = cdf.join(cdf, on="Chromosome", suffix="2").filter(
            pl.col("ClusterID") < pl.col("ClusterID2"),
            (pl.col("Gene IDs").list.set_intersection(pl.col("Gene IDs2")).list.len()) > 0,
        )
        # Get cluster groups
        mega_clusts = (
            overlapping_clusts.join(
                overlapping_clusts, left_on="ClusterID", right_on="ClusterID2", how="left"
            )
            .filter(pl.col("ClusterID_right").is_null())
            .group_by("ClusterID")
            .agg(mclust_expr)
        )
        no_overlap_clusts = cdf.filter(
            pl.col("ClusterID").is_in(overlapping_clusts["ClusterID2"].implode()).not_()
        )
        return no_overlap_clusts.join(mega_clusts, on="ClusterID", how="left").select(
            pl.col("ClusterID"),
            pl.col("Chromosome"),
            pl.coalesce(pl.col("Start bp_right"), pl.col("Start bp")).alias("Start bp"),
            pl.coalesce(pl.col("End bp_right"), pl.col("End bp")).alias("End bp"),
            pl.coalesce(pl.col("Gene IDs_right"), pl.col("Gene IDs").alias("Gene IDs")),
        )
