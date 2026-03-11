"""Various classes to facilitate the storage of data for cluster finding."""

import logging
import time
from pathlib import Path
from typing import Any

import polars as pl
from pydantic import BaseModel

from plantclusterfinder_py.pgdb_utils import read_in_pgdb_flat

_log = logging.getLogger(__name__)


class ReactionLabels(BaseModel):
    """Information regarding reactions labels from the reaction data input."""

    STClasses: list[str]
    """List of signature and tailoring enzyme classes."""
    SigOrTail: list[str]
    """List of types of enzymes, such as signature, tailoring, etc."""
    Rxn_STClasses: dict
    """Mapping of Reactions to enzyme classes of STClasses."""
    Rxn_SigOrTail: dict
    """Mapping of Reactions to enzyme types of SigOrTail."""


class HardBreaks(BaseModel):
    """Model to store the hard breaks in the sequence as sequence gaps."""

    model_config = {"arbitrary_types_allowed": True}
    """DataFrame containing the
        - chromosome/scaffold
        - identification if it is coding (A) or masked (N) DNA
        - start bp
        - length of masked region
        - intergenic region gap => start and end bp
        corrected for maximal size"""
    SeqGaps: pl.DataFrame | None
    """Number (0-n, default 2) of hypothetical genes that should be inserted.
    Only Used if InsertHypos is set"""
    HypoAmount: int = 2


class AnnotatedMetabolicDomains(BaseModel):
    """Metabolic Domains and their associated Reactions from PGDB."""

    model_config = {"arbitrary_types_allowed": True}
    RxnMetaDom: pl.DataFrame
    """
    DataFrame of all the items stored in the Metabolic Domain File
    Has the columns:
    - Reaction ID (PGDB reactions.dat UNIQUE-ID)
    - Pathway
    - Metabolic Domain
    """
    Rxn2MetaDom: dict
    """Mapping Between Reactions stored in Metabolic Domain file
    and metabolic domains of interest"""


class MclClusterData(BaseModel):
    """Storage for MCL cluster data."""

    model_config = {"arbitrary_types_allowed": True}
    Gene2MCLClusters: dict
    """ Mapping between genes in conversion-ID file and MCL-clustering file. """
    Protein_IDs: list[str] | None
    """ List of protein-IDs stored in protein FASTA file. """
    MCL_Clusters: pl.DataFrame
    """ DataFrame containing the MCL_ID, protein IDs, and Gene Names for each MCL cluster """


class GenesToReactions(BaseModel):
    """Mapping between gene identifiers and reactions. Has all gene and reaction IDs."""

    model_config = {"arbitrary_types_allowed": True}
    GenesAll: pl.DataFrame
    """DataFrame containing all the information
        (Attributes)that were searched for (described in
        genes_attributes, default is {'UNIQUE-ID','PRODUCT', ...
        'ACCESSION-1','ACCESSION-2'}) in genes.dat of the pgdb."""
    ReactionsAll: pl.DataFrame
    """DataFrame containing all the information
        (Attributes) that was searched for (described in
        reactions_attributes, default is {'UNIQUE-ID', ...
        'ENZYMATIC-REACTION', 'EC-NUMBER'}) in reactions.dat
                        of the pgdb."""
    Gene2Rxn: dict
    """Mapping between the genes in GenesAll and reactions in ReactionsAll"""


class Conversion(BaseModel):
    """Conversion maps between gene, transcript, and protein IDs."""

    model_config = {"arbitrary_types_allowed": True}
    """DataFrame containing Gene Name, Transcript Name, and Protein Name"""
    df: pl.DataFrame
    has_pgdb_ids: bool = False
    _cache: dict[str, Any] = {}

    def add_pgdb_id_column(self, pgdb_dir: Path) -> None:
        """Add the PGDB Gene UNIQUE-ID to the Conversion data.

        Args:
            pgdb_dir: Directory of the PGDB data
        """
        if self.has_pgdb_ids:
            return
        self.has_pgdb_ids = True
        proteins_table = (
            read_in_pgdb_flat(
                pgdb_dir / "proteins.dat", attrib=("UNIQUE-ID", "GENE", "COMMON-NAME", "ACCESSION-1")
            )
            .explode("GENE")
            .explode("COMMON-NAME")
        )
        # _log.error(proteins_table.head())
        # time.sleep(1)
        # exit(1)
        # join on the proteins table any way it'll fit, and accept the gene name that comes back
        pg_gene = self.df.join(proteins_table, left_on="Gene Name", right_on="GENE").select(
            pl.col("Gene Name"), pgdb_gene_id1=pl.col("Gene Name")
        )
        pg_prot = self.df.join(proteins_table, left_on="Protein Name", right_on="UNIQUE-ID").select(
            pl.col("Gene Name"), pgdb_gene_id2=pl.col("GENE")
        )
        pg_common_prot = self.df.join(
            proteins_table, left_on="Protein Name", right_on="COMMON-NAME"
        ).select(pl.col("Gene Name"), pgdb_gene_id3=pl.col("GENE"))
        pg_accession_prot = pl.DataFrame(schema={"Gene Name": pl.String, "pgdb_gene_id4": pl.String})
        if isinstance(proteins_table.schema.get("ACCESSION-1", pl.Null), pl.String):
            pg_accession_prot = self.df.join(
                proteins_table, left_on="Transcript Name", right_on="ACCESSION-1"
            ).select(pl.col("Gene Name"), pgdb_gene_id4=pl.col("GENE"))

        pg_ids = pl.concat([pg_gene, pg_prot, pg_common_prot, pg_accession_prot], how="align")
        pg_ids = pg_ids.select(
            pl.col("Gene Name"),
            pl.coalesce(
                (
                    pl.col("pgdb_gene_id1"),
                    pl.col("pgdb_gene_id2"),
                    pl.col("pgdb_gene_id3"),
                    pl.col("pgdb_gene_id4"),
                )
            ).alias("pgdb_gene_id"),
        )

        self.df = self.df.join(pg_ids, on="Gene Name", how="left")

    @property
    def map_gene(self) -> set[str]:
        """Unique list of gene IDs stored in conversion file."""
        if "MapG" in self._cache:
            return self._cache["MapG"]
        self._cache["MapG"] = set(self.df.get_column("Gene Name").unique())
        return self._cache["MapG"]

    @property
    def map_transcript(self) -> set[str]:
        """Unique list of transcript IDs stored in conversion file."""
        if "MapT" in self._cache:
            return self._cache["MapT"]
        self._cache["MapT"] = set(self.df.get_column("Transcript Name").unique())
        return self._cache["MapT"]

    @property
    def map_protein(self) -> set[str]:
        """Unique list of protein IDs stored in conversion file."""
        if "MapP" in self._cache:
            return self._cache["MapP"]
        self._cache["MapP"] = set(self.df.get_column("Protein Name").unique())
        return self._cache["MapP"]

    @property
    def map_gene_to_transcript(self) -> dict:
        """Mapping between MapG and MapT."""
        if "MapG_MapT" in self._cache:
            return self._cache["MapG_MapT"]
        self._cache["MapG_MapT"] = self.df.select(
            pl.col("Gene Name"), pl.col("Transcript Name")
        ).rows_by_key(key="Gene Name")
        return self._cache["MapG_MapT"]

    @property
    def map_gene_to_protein(self) -> dict:
        """Mapping between MapG and MapP."""
        if "MapG_MapP" in self._cache:
            return self._cache["MapG_MapP"]
        self._cache["MapG_MapP"] = self.df.select(
            pl.col("Gene Name"), pl.col("Protein Name")
        ).rows_by_key(key="Gene Name")
        return self._cache["MapG_MapP"]

    @property
    def map_proteins_to_genes(self) -> dict:
        """Mapping between MapG and MapP."""
        if "MapP_MapG" not in self._cache:
            self._cache["MapP_MapG"] = self.df.select(
                pl.col("Gene Name"), pl.col("Protein Name")
            ).rows_by_key(key="Protein Name")
        return self._cache["MapP_MapG"]


class Clusters(BaseModel):
    """Data class for cluster data.

    ClusterMeta, genome background ClusterBackground and step size search information.
    """

    model_config = {"arbitrary_types_allowed": True}
    # cluster information, stepsize, chromosome ID, gene start, gene end, number of enzymes,
    # number of non-enzymes, and number of hypothetical genes
    ClusterMeta: pl.DataFrame
    """DataFrame with cluster information
        - stepsize
        - chromosome ID
        - gene start
        - gene end
        - number of enzymes
        - number of non-enzymes
        - number of hypothetical genes
        - top_cluster (if perform_cutoff is run)"""
    ClusterBackground: pl.DataFrame  # Density information about the genome
    """Density information about the genome
        - window size
        - number of enzymes
        - occurence of window size in genome
        - frequency of seeing a window size of this size with at least as many enzymes in the genome
    """
    FinalStepSize: int
    """step size chosen by criterion"""
    MaxStepSize: int
    """maximal step size computed"""

    def perform_cutoff(self, top_clusters_percent: int) -> list[bool]:
        """---------------------------------------------------------------------------

        Function that (given clusters and a background density information of a
        genome) calculates whether the clusters are within the top n% of enzymes
        dense clusters.

        ---------------------------------------------------------------------------

        Args:
            clusters: Cluster Information (stepsize, chromosome id, gene start,
                    gene end, number of enzymes) and Density information
                    about the genome (window size, number of enzymes.
                    Occurrence of window size with this many enzymes in the genome,
                    occurrence of this window size in the genome. Frequency of seeing
                    a window size of this size with at least as many enzymes in the genome).
            top_clusters_percent: Clusters are labeled as being in the top n% if
                                they are within this top n% of enzyme dense
                                clusters compared to the background.

        Returns:
            List with size Number of Clusters. It is zero if
                the density of the cluster is above or the same as given
                by the input argument top_clusters_percent*100.

        """
        _log.debug("Perform cutoff for clusters in the top %d%%", top_clusters_percent)
        cluster_meta_in_top = (
            self.ClusterMeta.with_row_index("cluster_index")
            .with_columns(window_size=pl.sum_horizontal("enzyme_count", "non_enzyme_gene_count"))
            .join(
                self.ClusterBackground,
                on=["window_size", "enzyme_count"],
                how="left",
            )
            .select(
                pl.col("cluster_index"),
                top_percent_clusters=pl.col("winzyme_frequency").fill_null(0)
                < top_clusters_percent / 100,
            )
        )
        _log.debug("cluster meta in top %s", cluster_meta_in_top.head())
        # prevent re-adding the top_percent_clusters column
        if "top_percent_clusters" in self.ClusterMeta.columns:
            self.ClusterMeta.drop("tep_percent_clusters")
        self.ClusterMeta = self.ClusterMeta.with_row_index("cluster_index").join(
            cluster_meta_in_top, on="cluster_index", how="left"
        )

        top_clusters = cluster_meta_in_top["top_percent_clusters"].to_list()
        _log.debug(
            "Found that %d clusters out of %d are in top",
            len(list(filter(lambda x: x, top_clusters))),
            len(self.ClusterMeta),
        )
        return top_clusters
