"""Gene location data and methods for annotation and analysis."""

import logging
from pathlib import Path
from typing import Any

import polars as pl
from pydantic import BaseModel

from plantclusterfinder_py.cluster_datatypes import (
    Conversion,
    GenesToReactions,
    MclClusterData,
)

_log = logging.getLogger(__name__)


class GeneLoc(BaseModel):
    """--------------------------------------------------------------------------------

    A class for handling genomic data in tabular form using Polars DataFrame. This
    class provides functionality for loading genomic data, sorting genes, and mapping
    gene names to specific gene sets such as Genes, Proteins, or Transcripts.

    The default columns are:
    - gene name
    - gene start (bp)
    - gene end (bp)
    - chromosome name
    - strand
    - hypothetical
    Additionally, the following columns are added as part of clustering
    - gene ID
    - protein ID
    - transcript ID
    - reaction ID(s)

    The class will also hold a list of intergenic regions and hard breaks with the
    columns:
    - chromosome name
    - region start (bp)
    - region end (bp)
    - gene ID preceding
    - gene ID following

    --------------------------------------------------------------------------------

    Raises:
        ValueError: If an invalid mapping type is provided or if the selected gene set is empty.

    """

    model_config = {"arbitrary_types_allowed": True}
    file_path: Path
    """The filepath to the genomic data CSV file."""
    unmasked: bool = True
    """Option to ignore the error if you used unmasked FASTA"""
    is_converted: bool = False
    """Indicates if the Conversion file columns have been added."""
    is_sorted: bool = False
    """Indicates if genes have been sorted by chromosome location."""
    is_annotated: bool = False
    """Indicates if reactions have been annotated."""
    has_hard_breaks: bool = False
    """Indicates if hard breaks have been added"""
    _data: pl.DataFrame
    """The Polars DataFrame that stores the genomic data."""
    _cache: dict = {}

    def __init__(self, **input) -> None:
        """Load data from the specified filepath into the Polars DataFrame.

            Function that reads in TSV files with headers. The headers are used as
            fields in MATLAB structures. These field names are not allowed to have
            certain special characters. Following Characters are currently handled:
            ' ' (space), ')', '(', '#', '/', '-', '.'.
            It is possible that others exist that either should be avoided
            as part of header or the code needs to be adapted.

        Returns:
            **Output**: Polars DataFrame with all the content of the TSV file.

               Expected columns are:
               - gene name
               - gene start (BP)
               - gene end (BP)
               - chromosome name
               - strand
               - hypothetical

        """
        super().__init__(**input)
        _log.debug("Reading TSV file %s to get gene location data.", input["file_path"])
        self._data = pl.read_csv(
            input["file_path"],
            separator="\t",
            schema_overrides={"Strand": pl.datatypes.String},
        )
        if self.data.columns != [
            "Gene Name",
            "Gene Start (bp)",
            "Gene End (bp)",
            "Chromosome Name",
            "Strand",
        ]:
            raise Exception("Gene Location File does not match expected format.")
        if self.data.select(pl.col("Chromosome Name", "Gene Name")).is_duplicated().any():
            _log.warning(
                "There are duplicate gene records, typical of transcript read ingestion."
                + " Please run consolidate_gene_info() to de-duplicate"
            )
        self._data = self.data.with_columns(Hypothetical=pl.lit(False))

        _log.debug("Read %d genes from Gene Location file", self.data.height)

    def __len__(self):
        """Return the length of the DataFrame."""
        return self.data.height

    @property
    def data(self) -> pl.DataFrame:
        """Get access to the underlying DataFrame object."""
        return self._data

    @property
    def columns(self):
        """Get access to columns of DataFrame without needing to pull it directly."""
        return self.data.columns

    def sort_genes(self) -> None:
        """Sort the DataFrame.

        Fist by chromosome index and then conditionally by start or end based on strand.
        """
        _log.debug("Sorting GeneLoc")
        if self.is_sorted:
            _log.debug("GeneLoc is already sorted, returning early")
            return
        self._data = self.data.sort(
            by=[
                pl.col("Chromosome Name"),
                # Column could be "1"/"-1" or "+"/"-", this handles both
                pl.when(pl.col("Strand").str.contains("-", literal=True))
                # when - then Gene End < Gene Start
                .then(pl.col("Gene End (bp)"))
                # when + then Gene Start < Gene end
                .otherwise(pl.col("Gene Start (bp)")),
            ]
        )
        self.is_sorted = True

    def map_to(self, conversion: Conversion, mapping_type: str) -> None:
        """Map gene names to conversion.MapT and add a column 'Map_to_T'.

        This may be deprecated in a future release.

        Args:
            conversion (Conversion): An instance of Conversion
                containing gene sets (MapG, MapP, MapT).
            mapping_type (str): The type of mapping to use ('G', 'P', or 'T').

        Raises:
            ValueError: If the selected gene set is empty or if mapping_type is invalid.

        """
        _log.warning("GeneLoc.map_to() is not stable and may be deprecated in a future release.")

        if not self.is_converted:
            _log.warning("Adding Conversion to Gene Location")
            self.remove_protein_transcript_info(conversion)

        # Create a new column with the count of each gene
        column_name = (
            # Set the column name based on the mapping type
            f"Map_to_{mapping_type}"
        )
        match mapping_type:
            case "G":
                _log.debug("Mapping Gene Location to Gene Names")
                mapping_genes = conversion.map_gene
                mapped_genes = self._data.group_by("Gene Name").len("Map_to_G")
            case "P":
                _log.debug("Mapping Gene Location to Protein Names")
                mapping_genes = conversion.map_protein
                mapped_genes = (
                    self._data.group_by("Gene Name")
                    .n_unique()
                    .select(pl.col("Gene Name"), pl.col("Protein Name").alias(column_name))
                )
            case "T":
                _log.debug("Mapping Gene Location to Transcript Name")
                mapping_genes = conversion.map_transcript
                mapped_genes = (
                    self._data.group_by("Gene Name")
                    .n_unique()
                    .select(
                        pl.col("Gene Name"),
                        pl.col("Transcript Name").alias(column_name),
                    )
                )
            case _:
                _log.error("Mappting type %s is invalid.", mapping_type)
                raise ValueError(
                    "Invalid mapping type! Choose 'G' for Gene, 'P' for Protein,"
                    + " or 'T' for Transcript."
                )

        # Check if the selected gene set is valid
        if not mapping_genes:
            raise ValueError(f"Selected mapping type '{mapping_type}' has no genes to map.")

        self._data = self._data.join(mapped_genes, on="Gene Name")

    def remove_protein_transcript_info(
        self,
        conversion: Conversion,
    ) -> None:
        """---------------------------------------------------------------------------

        Function that adds the Conversion columns to the Gene Location data.

        ---------------------------------------------------------------------------

        Args:
            conversion: Conversion object constructed from conversion file

        """
        if self.is_converted:
            _log.debug("Gene Location already mapped to Conversion, returning early")
            return
        gene_location = self.data.join(conversion.df, on="Gene Name")
        if (
            # This statement returns True if duplicate keys exist in
            # "Gene Name","Protein Name",or "Transcript Name"
            gene_location.select(
                pl.col(["Gene Name", "Protein Name", "Transcript Name"]).is_unique().not_().any()
            )
            .transpose()
            .to_series(0)
            .any()
        ):
            raise Exception(
                "Multiple entries (transcripts) per gene in gtpf file."
                + f"Please clean up entries for {
                    gene_location.group_by('Gene Name').len().filter(pl.col('len') > 1)
                } in your input files.\n",
            )
        self._data = gene_location
        self.is_converted = True

    def __repr__(self) -> str:
        """Return a string representation of the DataFrame."""
        return self.data.__repr__()

    def consolidate_gene_info(self) -> list[str]:
        """---------------------------------------------------------------------------

        Function that collapses multiple entries (e.g. because it originally
        contained transcript or protein information) in the gene location data to
        one entry per gene. It takes the outermost end and start bp of all the
        protein or transcripts of a gene.

        ---------------------------------------------------------------------------

        Returns:
            - **Chromosomes**: A list with unique entries of all chromosome and scaffold
                        names.

        """
        if self.data.select(pl.col("Chromosome Name", "Gene Name")).is_duplicated().any():
            _log.warning("De-duplicating gene data")
            current_columns = self.columns
            current_columns.remove("Gene Start (bp)")
            current_columns.remove("Gene End (bp)")
            self._data = self.data.group_by(pl.col(current_columns)).agg(
                pl.col("Gene Start (bp)").min(), pl.col("Gene End (bp)").max()
            )

        return self.data.select(pl.col("Chromosome Name").unique()).to_series(0).to_list()

    def map_to_conversion_file(self) -> dict:
        """---------------------------------------------------------------------------

        Function that maps genes in the gene location file to the genes in the id
        conversion file.

        ---------------------------------------------------------------------------

        Returns:
            Mapping between the gene location file and the gene-ids
                    in the conversion file.

        """
        if not self.is_converted:
            raise Exception("Cannot map to gene locations before conversion file added")
        _log.debug("Finding the protein names associated with each gene name")
        return self.data.select(pl.col(["Gene Name", "Protein Name"])).rows_by_key(key="Gene Name")

    def remove_non_protein_fasta_gene_locations(
        self,
        conversion: Conversion,
        protein_ids: list[str] | None,
    ) -> dict[Any, Any]:
        """---------------------------------------------------------------------------

        Function that uses the protein fasta file to search for genes that encode
        proteins and then deletes all other entries in the gene location
        variable.

        ---------------------------------------------------------------------------

        Args:
            conversion: Conversion object between genes, transcripts, and proteins
            protein_ids: A cell array of all protein-IDs that are stored in the
                        protein fasta file.

        Returns:
            Mapping between the Gene location file and the gene-IDs in the
                    conversion file.

        """
        if not self.is_converted:
            _log.warning("Implicitly requested adding conversion to Gene Location")
            self.remove_protein_transcript_info(conversion)

        original_gene_count = self.data.height
        if protein_ids is None or len(protein_ids) == 0:
            _log.info("No protein IDs provided to remove, returning early.")
            return self.map_gene_to_protein
        protein_ids_df = pl.DataFrame({"Gene Name": protein_ids})
        self._data = self.data.join(protein_ids_df, on="Gene Name")
        _log.debug(
            f"Gene Location data filtered from {original_gene_count} to {self._data.height} rows"
        )
        return self.map_gene_to_protein

    @property
    def map_gene_to_protein(self):
        """Get mapping dictionary between Gene Name and Protein Name."""
        return self._data.select(pl.col(["Gene Name", "Protein Name"])).rows_by_key("Gene Name")

    def find_intergenic_regions(self) -> pl.DataFrame:
        """---------------------------------------------------------------------------

        Function computes all stretches of DNA that are not encoding at least one
        gene model.

        ---------------------------------------------------------------------------

        Returns:
            **InterChrom**: A DataFrame defining
                - chromosome,
                - start bp,
                - end bp,
                - gene before
                - gene after an intergenic region.

        """
        _log.debug("Finding intergenic regions of GeneLoc")
        if not self.is_sorted:
            _log.warning("Implicitly Requested Gene Sorting")
            self.sort_genes()

        if "InterChrom" in self._cache:
            _log.debug("Returning cached result.")
            return self._cache["InterChrom"]

        ret: list[dict[str, Any]] = []
        previous_gene = None
        """
        Iterate through *sorted* rows and find each region where there is no
            gene located (intergenic region).
        Compare the furthest out gene region in the chromosome to the nearest of the next gene.
         1. If they overlap, the next gene becomes the furthest out gene.
         2. If they don't overlap, an intergenic region is stored.
        """
        # TODO(j-schlesinger): replace with window functions if possible (issue #53) # noqa: FIX002
        # https://github.com/j-schlesinger/PlantClusterFinder.py/issues/53
        for chromosome, gene_start, gene_end, gene_name, strand in self.data.select(
            pl.col(
                [
                    "Chromosome Name",
                    "Gene Start (bp)",
                    "Gene End (bp)",
                    "Gene Name",
                    "Strand",
                ]
            )
        ).iter_rows():
            # in case gene is on - strand
            if "-" in strand:
                global_gene_start = gene_end
                global_gene_end = gene_start
            else:
                global_gene_start = gene_start
                global_gene_end = gene_end
            global_gene_start = min(gene_start, gene_end)
            if not (
                previous_gene is None
                or previous_gene["Gene End"] >= global_gene_start
                or chromosome != previous_gene["Chromosome"]
            ):
                ret.append(
                    {
                        "Chromosome Name": chromosome,
                        "Start": previous_gene["Gene End"] + 1,
                        "End": global_gene_start - 1,
                        "Gene Before": previous_gene["Gene Name"],
                        "Gene After": gene_name,
                    }
                )
            previous_gene = {
                "Chromosome": chromosome,
                "Gene End": global_gene_end,
                "Gene Name": gene_name,
            }

        self._cache["InterChrom"] = pl.DataFrame(ret)
        return self._cache["InterChrom"]

    def annotate_enzyme_info(
        self,
        map_g_gene_pgdb: pl.DataFrame,  # Dict[str, List[str]],
        gene_to_reaction: GenesToReactions,
    ) -> dict[str, bool]:
        """---------------------------------------------------------------------------

        Generates a dictionary for all genes in the gene location
        information and looks up if this gene has a reaction in the pgdb.

        ---------------------------------------------------------------------------

        Args:
            map_g_gene_pgdb: A mapping between Genes in the conversion file and genes
                        stored in the pgdb genes.dat file.
                        Columns are:
                        - pgdb_gene_id
                        - Gene Name
                        - Transcript Name
                        - Protein Name
            gene_to_reaction: Mapping of all genes in the Genes_all to the Reactions in
                Reactions_all

        Returns:
            A dict describing which genes have at least one reaction
                    annotated in the pgdb (True if at least one reaction)

        """
        if not self.is_converted:
            raise Exception("Conversion data not loaded into gene locations")

        if not self.is_annotated:
            genes_with_reactions = map_g_gene_pgdb.with_columns(
                pgdb_reaction_id=pl.col("pgdb_gene_id")
                .replace_strict(gene_to_reaction.Gene2Rxn, default=None)
                .flatten()
            )
            _log.debug(
                "Mapping %d genes with reactions to %d genes",
                len(genes_with_reactions),
                len(self.data),
            )
            # get gene IDs for PGDB data - may come from Protein Name or Gene Name
            self._data = self.data.join(
                map_g_gene_pgdb.select(pl.col("Protein Name"), pl.col("pgdb_gene_id")),
                on="Protein Name",
                how="left",
            )
            self._data = (
                self.data.join(
                    map_g_gene_pgdb.select(pl.col("Gene Name"), pl.col("pgdb_gene_id")),
                    on="Gene Name",
                    how="left",
                )
                .with_columns(
                    pgdb_gene_id=pl.coalesce(pl.col("pgdb_gene_id"), pl.col("pgdb_gene_id_right"))
                )
                .drop("pgdb_gene_id_right")
            )
            self._data = (
                self.data.join(
                    map_g_gene_pgdb.select(pl.col("Transcript Name"), pl.col("pgdb_gene_id")),
                    on="Transcript Name",
                    how="left",
                )
                .with_columns(
                    pgdb_gene_id=pl.coalesce(pl.col("pgdb_gene_id"), pl.col("pgdb_gene_id_right"))
                )
                .drop("pgdb_gene_id_right")
            )
            if self.data.filter(pl.col("pgdb_gene_id").is_not_null()).height == 0:
                _log.error("Cound not map PGDB Gene IDs to Protein Name")
                _log.debug("Map G to Gene PGDB DF: %s", map_g_gene_pgdb)
                _log.debug("Current GeneLocation DF: %s", self.data)
                print(self.data)
                raise Exception("PGDB Protein IDs not mapping to gtpf Protein IDs.")
            # ID genes with any reactions
            self._data = self.data.join(
                genes_with_reactions.select(
                    pl.col("pgdb_gene_id"), pl.col("pgdb_reaction_id").list.first()
                ),
                on="pgdb_gene_id",
                how="left",
            ).with_columns(has_reaction=pl.col("pgdb_reaction_id").is_not_null())
            if self.data.filter(pl.col("has_reaction")).height == 0:
                _log.error("Cound not find PGDB Reaction IDs for any genes")
                _log.debug(
                    "Genes with Reactions DF (filtered for reactions): %s",
                    genes_with_reactions.filter(pl.col("pgdb_reaction_id").is_not_null()),
                )
                _log.debug(
                    "Current GeneLocation DF (filtered for pgdb_gene_id): %s",
                    self.data.filter(pl.col("pgdb_gene_id").is_not_null()),
                )
                raise Exception("Reactions from PGDB could not be mapped to genes.")
            rxn_2_gene = []
            for gene, reactions in gene_to_reaction.Gene2Rxn.items():
                for rxn in reactions:
                    rxn_2_gene.append({"Reaction": rxn[0][0], "Gene": gene})
            rxn_2_gene_df = pl.DataFrame(rxn_2_gene)
            _log.debug("Reaction to Gene DF: %s", rxn_2_gene_df)
            rxn_2_gene_df = rxn_2_gene_df.join(
                gene_to_reaction.ReactionsAll,
                left_on="Reaction",
                right_on="UNIQUE-ID",
                how="left",
            )
            gene_2_ec = (
                rxn_2_gene_df.group_by("Gene")
                .agg(pl.col("EC-NUMBER").flatten())
                .with_columns(
                    pl.when(pl.col("EC-NUMBER").list.first().is_null())
                    .then(None)
                    .otherwise(pl.col("EC-NUMBER"))
                    .alias("EC-NUMBER")
                )
            )
            self._data = self.data.join(
                gene_2_ec,
                left_on="pgdb_gene_id",
                right_on="Gene",
                how="left",
            )
            self.is_annotated = True
        return {
            k: list(v)[0]
            for k, v in self._data.select(pl.col(["Gene Name", "has_reaction"]))
            .rows_by_key("Gene Name", unique=True)
            .items()
        }

    def introduce_hard_breaks(
        self,
        sequence_gaps: pl.DataFrame,  # from annotate_sequencing_gaps()
        max_sequence_gap_size: int | None = None,  # from args
        max_inter_gene_distance: int | None = None,  # from args
        percentile_for_max_inter_gene_distance: float | None = None,  # from args
        max_inter_gene_distance_by_median: int | None = None,  # from args
        hard_break_gaps: list[int] | None = None,  # from args -
    ) -> None:
        """---------------------------------------------------------------------------

        Function that prevents formation of clusters Sequencing gaps that
        exceed max_sequence_gap_size or that are defined not to be crossed by
        hard_break_gaps.

        Creates a set of pseudo-chromosomes in the GeneLoc data which breaks the
        sequence across hard breaks given the gap size criteria.

        ---------------------------------------------------------------------------

        Args:
            sequence_gaps: A DataFrame containing the
                    - chromosome/scaffold,
                    - the identification if it is coding (A) or masked (N) DNA
                    - start BP and
                    - length of a masked region.
                    In addition there is the information stored in with
                    intergenic region the gap is encoded and its start and end
                    bp, corrected for maximal size (= intergenic region).
            max_sequence_gap_size: Maximal masked nucleotide region (in bp) that a cluster
                            is allowed to bridge (-1 or None is inactive).
            max_inter_gene_distance: Maximal intergenic distance that still can be
                            gene cluster.
            percentile_for_max_inter_gene_distance: Prevents clusters to be spanned over
                                            the top n% of intergenic distances.
                                            (-1 or None is inactive)
            max_inter_gene_distance_by_median: Maximal intergenic distance that still can
                                    be spanned by a gene cluster: this number
                                    times the median of the gene sizes.
            hard_break_gaps: A list of numbers. The masked nucleotides sequences of
                            these specific lengths are preventing cluster formation.
                            In other words: specific lengths of gaps that indicate hard breaks.
                            (Default is an empty list []).

        """
        if hard_break_gaps is None:
            hard_break_gaps = []
        if self.has_hard_breaks:
            return
        self.has_hard_breaks = True
        if max_sequence_gap_size is None or max_sequence_gap_size < 1:
            self.has_hard_breaks = True
            _log.warning("Max sequence gap size not set or less than 1.")
        final_max_inter_gene_distance = max_inter_gene_distance

        if (
            percentile_for_max_inter_gene_distance is not None
            and percentile_for_max_inter_gene_distance > 1
        ):
            # set max intergenic distance by percentile
            pctile_set_max_intergenic_distance = (
                self.find_intergenic_regions()
                .select(len=pl.col("End") - pl.col("Start") + 1)
                .quantile(percentile_for_max_inter_gene_distance / 100)
            )["len"][0]
            _log.info(
                "Percentile used for intergenic distance: %ibp",
                percentile_for_max_inter_gene_distance,
            )
            _log.info(
                "Clusters broken when spanning distance larger: %ibp.\n",
                pctile_set_max_intergenic_distance,
            )
            # update max_inter_gene_distance using temporary max
            #
            final_max_inter_gene_distance = (
                pctile_set_max_intergenic_distance
                if final_max_inter_gene_distance is None or final_max_inter_gene_distance < 1
                else min(final_max_inter_gene_distance, pctile_set_max_intergenic_distance)
            )
        if max_inter_gene_distance_by_median is not None and max_inter_gene_distance_by_median > 0:
            median_set_max_intergenic_distance = self.data.select(
                len=pl.col("Gene End (bp)") - pl.col("Gene Start (bp)")
            ).quantile(0.50)["len"][0]
            final_max_inter_gene_distance = (
                median_set_max_intergenic_distance * max_inter_gene_distance_by_median
                if final_max_inter_gene_distance is None or final_max_inter_gene_distance == -1
                else min(
                    final_max_inter_gene_distance,
                    median_set_max_intergenic_distance * max_inter_gene_distance_by_median,
                )
                # TODO(j-schlesinger): Why are these vals are being multiplied # noqa: FIX002, TD003
                # https://github.com/j-schlesinger/PlantClusterFinder.py/issues/54
            )
        # this holds the SeqGaps Hard data
        sequence_gaps = sequence_gaps.with_columns(
            hard_gap=pl.col("seq_gaps_length").is_in(hard_break_gaps) | pl.col("seq_gaps_length")
            > max_sequence_gap_size
        )

        hard_breaks = (
            self.find_intergenic_regions()
            .join_where(  # using join_where as a cheeky filter
                # join intergenic regions to sequence gaps where gap starts in an intergenic region
                sequence_gaps,
                pl.col("Start").is_between(
                    pl.col("start index"),
                    pl.col("end index"),
                ),
            )
            # only get hard gaps
            .filter(pl.col("hard_gap"))
            .select(pl.col(["Chromosome Name", "Start", "End", "Gene Before", "Gene After"]))
            # cumulative count of pseudochromosomes
            .with_columns(pseudochromosome=pl.col("Chromosome Name").cum_count())
        )
        prev_end = 0
        prev_gene = None
        prev_chrom = None
        pseudochromosome = []
        if len(hard_breaks) == 0:
            return
        for h_break in hard_breaks.iter_rows():
            # convert hard breaks into pseudochromosome sections
            (chromosome, start, end, prev_gene, next_gene, pseudochrom_count) = h_break
            if prev_chrom != chromosome:
                prev_chrom = chromosome
                prev_gene = None
            pseudochromosome.append(
                {
                    "pchrom_start": prev_end,
                    "pchrom_end": start,
                    "pchromosome name": f"{chromosome}_{pseudochrom_count - 1}",
                    "Old Chromosome Name": chromosome,
                    "pchrom_start_gene": prev_gene,
                    "pchrom_end_gene": prev_gene,
                }
            )
            prev_end = end
            prev_gene = next_gene
        pseudochrom_df = pl.DataFrame(pseudochromosome)
        # Add pseudochromosomes to genes which are between the break boundaries
        updated_genes = (
            self.data.join_where(
                pseudochrom_df,
                pl.col("Gene Start (bp)") >= pl.col("pchrom_start"),
                pl.col("Gene End (bp)") < pl.col("pchrom_end"),
                pl.col("Chromosome Name") == pl.col("Old Chromosome Name"),
            )
            .with_columns(
                pl.col("pchromosome name").alias("Chromosome Name"),
            )
            .drop(
                [
                    "pchrom_start",
                    "pchrom_end",
                    "pchrom_start_gene",
                    "pchrom_end_gene",
                    "Old Chromosome Name",
                ]
            )
        )
        self._data = self.data.update(updated_genes, on="Gene Name", how="left")

    def insert_genes(self, other_genes: pl.DataFrame) -> None:
        """---

        Insert genes into the GeneLoc data.,
        Modifies in place and returns DataFrame.

        ---

        Args:
            other_genes: a Polars DataFrame containing the other genes columns.
                         **this DataFrame must have the same columns as this object**

        """
        _ = self._data.extend(other_genes)

    def compute_local_duplication(
        self,
        mcl_cluster_data: MclClusterData,
        local_duplication_kb_range: int,
        local_duplication_gene_range: int,
    ) -> pl.DataFrame:
        """---------------------------------------------------------------------------

        Function that computes whether a gene is locally duplicated, and if so,
        which MCL cluster it is associated with.

        ---------------------------------------------------------------------------

        Args:
            mcl_cluster_data : Mapping between the genes in the conversion-ID file and the
                    MCL-clustering file.
            local_duplication_kb_range: Genes that are encoded more than this bp apart, are not
                        called local duplicates.
            local_duplication_gene_range: Number of genes that space two genes from each other.
                        Genes that have more than this many intervening genes
                        between them are not called local duplicates.

        Returns:
            DataFrame with the columns:

            - GeneID: "Gene Name" in GeneLocation
            - is_locally_duplicated: A vector of length of the number of Genes in GeneLoc object,
                being True if the gene was considered to be locally duplicated (a
                gene was within given bp distance and not spaced more than
                with given number of genes from another gene with the same MCL
                clustering given by mcl_cluster_data).
            - ClusterIds: The MCL clusters of the gene (list)

        """
        _log.debug("Computing local duplication")

        if "is_locally_duplicated" in self.data.columns:
            _log.warning("Local duplication already exists, replacing with new inputs")
            self._data.drop("is_locally_duplicated")
            self._data.drop("mcl_clusters")
        # For all clusters, select only those with multiple genes
        possible_ld_clusters = mcl_cluster_data.MCL_Clusters.filter(pl.col("Cluster_Size") > 1)
        # Explode on Gene Name, and select just cluster ID and gene name
        possible_ld_g_w_clust = possible_ld_clusters.explode("Gene_Names").select(
            pl.col("ID").alias("mcl_cluster_id"),
            pl.col("Gene_Names")
            .list.first()
            .alias("Gene Name"),  # should only have one value because we exploded
        )
        # Join on GeneLocation data to get gene start and end point
        possible_ld_genes = possible_ld_g_w_clust.join(
            self.data.with_row_index("gene_index"), on="Gene Name"
        )
        # self-join with cluster criteria
        #  - same mcl_cluster_id
        #  - same Chromosome Name
        #  - different Gene Name
        #  - gene index <= gene duplication range apart
        #  - gene position <= kb duplication range apart
        ld_genes = (
            possible_ld_genes.join_where(
                possible_ld_genes,
                pl.col("mcl_cluster_id") == pl.col("mcl_cluster_id_right"),
                pl.col("Chromosome Name") == pl.col("Chromosome Name_right"),
                pl.col("Gene Name") != pl.col("Gene Name_right"),
                suffix="_right",
            )
            .with_columns(
                intervening_genes=(pl.col("gene_index") - pl.col("gene_index_right")).abs(),
                kb_between_genes=pl.min_horizontal(
                    (pl.col("Gene Start (bp)") - pl.col("Gene End (bp)_right")).abs(),
                    (pl.col("Gene End (bp)") - pl.col("Gene Start (bp)_right")).abs(),
                ),
            )
            .filter(
                pl.col("intervening_genes") <= local_duplication_gene_range,
                pl.col("kb_between_genes") <= local_duplication_kb_range,
            )
            .group_by("Gene Name")
            .agg(
                pl.lit(True).alias("is_locally_duplicated"),
            )
        )

        # Get is_locally_duplicated from join, and use replace_strict to get clusters per gene
        mcl_df = pl.from_dicts(
            [
                {"Gene Name": k[0], "mcl_clusters": v}
                for k, v in mcl_cluster_data.Gene2MCLClusters.items()
            ]
        )
        self._data = (
            self.data.join(mcl_df, on="Gene Name", how="left")
            .join(ld_genes, on="Gene Name", how="left")
            .with_columns(is_locally_duplicated=pl.col("is_locally_duplicated").fill_null(False))
        )

        return self.data

    def write_parquet(self, filename: Path) -> None:
        """Write GeneLoc data as Parquet file."""
        self.data.write_parquet(filename)


# Example usage
if __name__ == "__main__":
    gene_loc = GeneLoc(
        file_path=Path("./csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt"),
        unmasked=False,
    )

    # Get the length of the DataFrame
    print(f"Number of entries: {len(gene_loc)}")
    print("The unsorted dataframe:\n")
    print(gene_loc)

    # Sort the genes
    gene_loc.sort_genes()

    # Display the sorted DataFrame
    print("The sorted dataframe:\n")
    print(gene_loc)
