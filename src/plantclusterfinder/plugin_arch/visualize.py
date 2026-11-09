"""Built-in plugin for computing enzyme statistics.

This plugin demonstrates how to consume PlantClusterFinder results
and produce a derived analysis without modifying core pipeline behavior.
"""

import logging
from collections import namedtuple
from collections.abc import Iterable
from pathlib import Path
from typing import NamedTuple

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pyfaidx
from matplotlib.figure import Figure
from pluggy import PluginValidationError

from plantclusterfinder.plugin_arch.context import ReturnContext, RunContext
from plantclusterfinder.plugin_arch.hookspecs import hookimpl

_log = logging.getLogger(__name__)

# Define font sizes
SIZE_DEFAULT = 14
BAR_HEIGHT = 0.6

scene_colors = ["67f491", "edbf85", "c6d4ff", "2e603f", "3d1f3a", "777FFF"]


class VisualizationPlugin:
    """A hook specification namespace."""

    expected_outputs = None

    @hookimpl
    def process_results(self, results: pl.DataFrame, context: RunContext) -> ReturnContext:
        """Produce visualization of BGCs on chromosomes.

        Args:
            results: DataFrame of gene location data
            context: RunContext object with plugin metadata
        """
        output_path = context.output_dir / f"cluster_location.{context.run_name}.png"
        if context.genome_file is None:
            raise PluginValidationError(
                self, "Cannot visualize clusters without genome_file set in Plugin Context"
            )
        fai_file = context.genome_file.parent / (context.genome_file.name + ".fai")
        if not fai_file.exists():
            # just create the .fai file
            _log.warning("FASTA index (fai) file not found, creating one")
            _ = pyfaidx.Fasta(
                context.genome_file, fai_file, as_raw=False, build_index=True
            )  # not using the pyfaidx result
            if not fai_file.exists():
                raise PluginValidationError(self, "FASTA index for genome file cannot be created.")
        if "cluster_filename" not in context.output_dict:
            raise PluginValidationError(
                self, "Cannot visualize clusters without cluster data (cluster_filename)"
            )
        try:
            clusters = pl.read_csv(context.output_dict["cluster_filename"], separator="\t")
        except Exception as e:
            raise PluginValidationError(
                self, "Error occurred when trying to read cluster data."
            ) from e
        _log.info("Writing output to %s", output_path)
        self.create_graph_for_clusters(fai_file, clusters=clusters, figure_path=output_path)
        return ReturnContext(plugin_name="visualize_bgc_plugin", output_filename=output_path)

    def plot_genomic_regions(
        self,
        clusters: pl.DataFrame,
        chrom_lengths: dict[str, int],
        categories: tuple[str] = ("Predicted BGC",),
    ) -> Figure:
        """Create a plot of the Biosynthetic Gene Clusters in the chromosomes.

        Args:
            clusters: DataFrame of BGCs from PlantClusterFinder
            chrom_lengths: dictionary of chromosome IDs and lengths
            categories: names of categories in input (currently unused)
        """
        plt.rcParams["font.sans-serif"] = ["Roboto", "DejaVu Sans"]
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.size"] = SIZE_DEFAULT
        plt.rcParams["axes.titlesize"] = SIZE_DEFAULT
        plt.rcParams["axes.labelsize"] = SIZE_DEFAULT
        plt.rcParams["xtick.labelsize"] = SIZE_DEFAULT
        plt.rcParams["ytick.labelsize"] = SIZE_DEFAULT
        plt.rcParams["axes.grid.axis"] = "x"
        annotated = []
        # Unique chromosomes and colors for metadata
        chrom_sorted = sorted(chrom_lengths.items(), key=lambda x: -x[1])
        chroms = [k for k, _ in chrom_sorted]  # chromosome names
        max_chrom_length = chrom_sorted[0][1]

        # Define a nice color palette:
        colors = [f"#{scene_colors[x]}" for x in [0, 5, 2]]  # ["#2B2F42", "#8D99AE", "#EF233C"]
        chrom_color = "darkgrey"

        color_map = {cat: colors[i] for i, cat in enumerate(categories)}
        df = clusters

        fig, ax = plt.subplots(figsize=(10, 6))

        # Map chromosome names to y-coordinates
        chrom_y_map = {chrom: i for i, chrom in enumerate(chroms)}

        # 2. Draw background bars for chromosomes
        for row in chrom_sorted:
            chrom, length = row
            y_pos = chrom_y_map[chrom]
            ax.add_patch(
                patches.Rectangle(
                    (0, y_pos - 0.2),
                    int(length),
                    BAR_HEIGHT,
                    color=chrom_color,
                    ec="black",
                    zorder=1,
                )
            )

        df = df.rename({col: col.replace(" ", "_") for col in df.columns}, strict=True)
        # 3. Plot specific regions
        for row_t in df.iter_rows(named=True):
            row_d = dict(row_t)
            true_chromosome: str = str(row_d["Chromosome"])
            start: int = int(row_d["Start_bp"])
            width: int = int(row_d["End_bp"] - start)
            gene_ids: list[str] = row_d["Gene_IDs"]
            y_pos = chrom_y_map[true_chromosome]
            bgc_type = "Predicted BGC"
            bgc_color_type = bgc_type
            known_cluster_bool = False
            high_coexpression_bool = False
            z_order_val = 1
            if known_cluster_bool:
                exp_bgc_genes: dict[str, str] = {}
                cluster_name = [exp_bgc_genes[x] for x in gene_ids if x in exp_bgc_genes][0]
                d = width
                z_order_val = 3
                polygon_shape = np.array(
                    (
                        (start, y_pos + BAR_HEIGHT * 0.8),
                        (start + d / 2, y_pos + 0.05 + BAR_HEIGHT * 0.8),
                        (start + d, y_pos + BAR_HEIGHT * 0.8),
                    )
                )

                if cluster_name not in annotated:
                    triangle = patches.Polygon(
                        polygon_shape,
                        closed=True,
                        color=color_map["Established Cluster"],
                        zorder=4,
                        label=cluster_name,
                    )
                    ax.add_patch(triangle)
                    ax.annotate(
                        cluster_name,
                        (0.5, -1),
                        xycoords=triangle,
                        color=color_map["Established Cluster"],
                        size="x-small",
                        rotation=20.0,
                    )
                    annotated.append(cluster_name)
            if high_coexpression_bool:
                bgc_color_type = "High Coexpression Cluster"
                z_order_val = 4

            ax.add_patch(
                patches.Rectangle(
                    (start, y_pos - (0.2) * 1.15),
                    width,
                    BAR_HEIGHT * 1.1,
                    color=color_map[bgc_color_type],
                    label=bgc_type,
                    zorder=z_order_val,
                )
            )

        # 4. Final Formatting
        # Hide the all but the bottom spines (axis lines)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["top"].set_visible(False)

        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position("left")
        ax.xaxis.set_ticks_position("bottom")

        ax.set_yticks(range(len(chroms)))
        chrom_names = chroms
        ax.set_yticklabels(chrom_names)
        ax.set_xlabel("Position (Mbp)")
        ax.set_title(
            "Biosynthetic Gene Clusters in Arabidopsis as Predicted by PlantClusterFinder.py"
        )
        ax.set_xlim(0, max_chrom_length * 1.05)
        ax.set_ylim(-0.5, len(chroms) - 0.5)
        ax.set_facecolor("dimgrey")

        # Handle Legend (avoid duplicates)
        handles = [patches.Patch(color=color_map[cat], label=cat) for cat in categories]
        ax.legend(handles=handles, title="Cluster Context", loc="lower right")

        plt.tight_layout()
        return fig

    def create_graph_for_clusters(self, fai_file: Path, clusters: pl.DataFrame, figure_path: Path):
        """Create a graph for the BGCs on a chromosome.

        Args:
            fai_file: Location of a .fai file which indicates chromosome sizes
            clusters: Polars DataFrame with Cluster data
            figure_path: Path of figure image output
        """
        # chromosome lengths
        chrom_lengths_df = pl.read_csv(
            fai_file,
            separator="\t",
            has_header=False,
            new_columns=[
                "Chromosome",
                "Length",
                "FASTA_start_index",
                "FASTA_column_width",
                "unknown",
            ],
        )
        chrom_lengths = {
            k: v[0]
            for k, v in chrom_lengths_df.filter(pl.col("Length") > 1e5)
            .select(pl.col("Chromosome"), pl.col("Length"))
            .rows_by_key(key="Chromosome", unique=True)
            .items()
        }
        fig = self.plot_genomic_regions(clusters, chrom_lengths)
        fig.savefig(figure_path)
