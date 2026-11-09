"""Various classes to facilitate the storage of data for cluster finding."""

import logging

import polars as pl
from pydantic import BaseModel

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


class Clusters(BaseModel):
    """Data class for cluster data.

    ClusterMeta, genome background ClusterBackground and step size search information.

    clusters: Cluster Information (stepsize, chromosome id, gene start,
            gene end, number of enzymes) and Density information
            about the genome (window size, number of enzymes.
            Occurrence of window size with this many enzymes in the genome,
            occurrence of this window size in the genome. Frequency of seeing
            a window size of this size with at least as many enzymes in the genome).
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
