import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import polars as pl
from pydantic import BaseModel

from plantclusterfinder_py.cluster_datatypes import Conversion, GenesToReactions

log = logging.getLogger(__name__)


class GeneLoc(BaseModel):
    """
    --------------------------------------------------------------------------------

    A class for handling genomic data in tabular form using Polars DataFrame. This
    class provides functionality for loading genomic data, sorting genes, and mapping
    gene names to specific gene sets such as Genes, Proteins, or Transcripts.

    --------------------------------------------------------------------------------

    Attributes:
        file_path (str): The file path to the genomic data CSV file.
        unmasked: Option to ignore the error if you use unmasked DNA files.
        _data (pl.DataFrame): The Polars DataFrame that stores the genomic data.
        _converted: Indicates if the Conversion file columns have been added.
        _sorted: Indicates if genes have been sorted by chromosome location.
        _annotated: Indicates if reactions have been annotated.
    Properties:
        data: Function to access _data attribute.
        columns: data.columns property
        BiomG_G: dictionary that associates gene name and protein name

    Methods:
        sort_genes: Sorts genes by chromosome and position (start or end) depending on strand orientation.
        map_to: Maps gene names to a specified gene set (Gene, Protein, or Transcript) and adds a new column
                with the mapping counts.
        remove_protein_transcript_info: Adds Conversion file columns to GeneLocation data.
        consolidate_gene_info: Returns a summary of gene information.
        map_to_conversion_file: Returns a mapping between conversion file and genelocation data.
        remove_non_protein_fasta_gene_locations: Remove gene locations not found in the protein FASTA file.
        find_intergenic_regions: Finds the region between genes and returns a DataFrame representing that.
        annotate_enzyme_info: Annotates if a gene has a reaction in the PGDB data.
        __repr__: Returns a string representation of the DataFrame for quick inspection.
        __len__: Returns the number of rows (genes) in the DataFrame.

    Raises:
        ValueError: If an invalid mapping type is provided or if the selected gene set is empty.
    """

    model_config = {"arbitrary_types_allowed": True}
    file_path: Path
    unmasked: bool = True
    _converted: bool = False
    _sorted: bool = False
    _annotated: bool = False
    _data: pl.DataFrame
    _cache: dict = dict()

    def __init__(self, **input) -> None:
        """
        Load data from the specified file path into the Polars DataFrame.
            Function that reads in TSV files with headers. The headers are used as
            fields in matlab structures. These field names are not allowed to have
            certain special characters. Follwoing Characters are currently handeled:
            ' ' (space), ')', '(', '#', '/', '-', '.'. It is possible that others
            exist that need to either be avoided as part of header or the code needs
            to be adapted.

        Returns:
            **Output**: Polars dataframe with all the content of the TSV file.

               Expected columns are:
               - gene name
               - gene start (bp)
               - gene end (bp)
               - chromosome name
               - strand
        """
        super().__init__(**input)
        log.debug(f"Reading TSV f{input["file_path"]} to get gene location data.")
        self._data = pl.read_csv(
            input["file_path"],
            separator="\t",
            schema_overrides={"Strand": pl.datatypes.String},
        )
        if self._data.columns != [
            "Gene Name",
            "Gene Start (bp)",
            "Gene End (bp)",
            "Chromosome Name",
            "Strand",
        ]:
            raise Exception("Gene Location File does not match expected format.")
        if (
            self.data.select(pl.col("Chromosome Name", "Gene Name"))
            .is_duplicated()
            .any()
        ):
            log.warning(
                "There are duplicate gene records, typical of transcript read ingestion. Please run consolidate_gene_info() to de-duplicate"
            )

        log.debug(f"Read {self._data.height} genes from Gene Location file")

    def __len__(self):
        """Return the length of the DataFrame."""
        return self._data.height

    @property
    def data(self) -> pl.DataFrame:
        """Get access to the underlying DataFrame object."""
        return self._data

    @property
    def columns(self):
        """Get access to columns of DataFrame without needing to pul it directly"""
        return self._data.columns

    def sort_genes(self) -> None:
        """Sort the DataFrame by chromosome index and then conditionally by start or end based on strand."""
        log.debug("Sorting GeneLoc")
        if self._sorted:
            log.debug("GeneLoc is already sorted, returning early")
            return
        self._data = self.data.sort(
            by=[
                pl.col("Chromosome Name"),
                # Column could be "1"/"-1" or "+"/"-", this handles both
                pl.when(pl.col("Strand").str.contains("-", literal=True))
                .then(pl.col("Gene End (bp)"))
                .otherwise(pl.col("Gene Start (bp)")),
            ]
        )
        self._sorted = True

    def map_to(self, conversion: Conversion, mapping_type: str) -> None:
        """
        Map gene names to conversion.MapT and add a column 'Map_to_T'.

        This may be deprecated in a future release.

        Args:
            conversion (Conversion): An instance of Conversion containing gene sets (MapG, MapP, MapT).
            mapping_type (str): The type of mapping to use ('G', 'P', or 'T').

        Raises:
            ValueError: If the selected gene set is empty or if mapping_type is invalid.
        """
        log.warning(
            "GeneLoc.map_to() is not stable and may be deprecated in a future release."
        )

        if not self._converted:
            log.warning("Adding Conversion to Gene Location")
            self.remove_protein_transcript_info(conversion)

        # Create a new column with the count of each gene
        column_name = (
            # Set the column name based on the mapping type
            f"Map_to_{mapping_type}"
        )
        match mapping_type:
            case "G":
                log.debug("Mapping Gene Location to Gene Names")
                mapping_genes = conversion.MapG
                mapped_genes = self._data.group_by("Gene Name").len("Map_to_G")
            case "P":
                log.debug("Mapping Gene Location to Protein Names")
                mapping_genes = conversion.MapP
                mapped_genes = (
                    self._data.group_by("Gene Name")
                    .n_unique()
                    .select(
                        pl.col("Gene Name"), pl.col("Protein Name").alias(column_name)
                    )
                )
            case "T":
                log.debug("Mapping Gene Location to Transcript Name")
                mapping_genes = conversion.MapT
                mapped_genes = (
                    self._data.group_by("Gene Name")
                    .n_unique()
                    .select(
                        pl.col("Gene Name"),
                        pl.col("Transcript Name").alias(column_name),
                    )
                )
            case _:
                log.error(f"Mappting type {mapping_type} is invalid.")
                raise ValueError(
                    "Invalid mapping type! Choose 'G' for Gene, 'P' for Protein, or 'T' for Transcript."
                )

        # Check if the selected gene set is valid
        if not mapping_genes:
            raise ValueError(
                f"Selected mapping type '{mapping_type}' has no genes to map."
            )

        self._data = self._data.join(mapped_genes, on="Gene Name")

    def remove_protein_transcript_info(
        self,
        conversion: Conversion,
    ) -> None:
        """
        ---------------------------------------------------------------------------

        Function that adds the Conversion columns to the Gene Location data.

        ---------------------------------------------------------------------------

        Args:
            conversion: Conversion object constructed from conversion file
        """
        if self._converted:
            log.debug("Gene Location already mapped to Conversion, returning early")
            return
        vGeneLocation = self.data.join(conversion.df, on="Gene Name")
        if (
            # This statement returns True if duplicate keys exist in "Gene Name","Protein Name",or "Transcript Name"
            vGeneLocation.select(
                pl.col(["Gene Name", "Protein Name", "Transcript Name"])
                .is_unique()
                .not_()
                .any()
            )
            .transpose()
            .to_series(0)
            .any()
        ):
            raise Exception(
                "Gene / Transcript / Protein location file (Gene position file) matches several entries to more than one gene.\n"
                + f"Please clean up entries for {vGeneLocation.group_by(
                    "Gene Name").len()} in your input files.\n",
            )
        self._data = vGeneLocation
        self._converted = True

    def __repr__(self) -> str:
        """Return a string representation of the DataFrame."""
        return self._data.__repr__()

    def consolidate_gene_info(self) -> List[str]:
        """
        ---------------------------------------------------------------------------

        Function that collapses multiple entries (e.g. because it originally
        contained transcript or protein information) in the gene location data to
        one entry per gene. It takes the outermost end and start bp of all the
        protein or transcripts of a gene.

        ---------------------------------------------------------------------------

        Returns:
            - **vUniqueChrom**: A list with unique entries of all chromosome and scaffold
                        names.
        """
        if (
            self.data.select(pl.col("Chromosome Name", "Gene Name"))
            .is_duplicated()
            .any()
        ):
            log.warning("De-duplicating gene data")
            current_columns = self.columns
            current_columns.remove("Gene Start (bp)")
            current_columns.remove("Gene End (bp)")
            self._data = self.data.group_by(pl.col(current_columns)).agg(
                pl.col("Gene Start (bp)").min(), pl.col("Gene End (bp)").max()
            )

        return (
            self.data.select(pl.col("Chromosome Name").unique()).to_series(0).to_list()
        )

    def map_to_conversion_file(self, vMapG: set[str]) -> dict:
        """
        ---------------------------------------------------------------------------

        function that maps genes in the gene location file to the genes in the id
        conversion file.

        ---------------------------------------------------------------------------
        Args:
            vGeneLocation: DataFrame about the gene location information,
                            generated in plantclusterfinder function.
            vMapG: unique list of gene ids stored in the conversion file.

        Returns:
            **vBiomG_G**: mapping between the gene location file and the gene-ids
                    in the conversion file.
        """

        if not self._converted:
            raise Exception("Cannot map to gene locations before conversion file added")
        log.debug("Finding the protein names associated with each gene name")
        return self.data.select(pl.col(["Gene Name", "Protein Name"])).rows_by_key(
            key="Gene Name"
        )

    def remove_non_protein_fasta_gene_locations(
        self,
        conversion: Conversion,
        vBiomG_G: Dict[Any, Any],
        vProtein_IDs: Optional[List[str]],
    ) -> Dict[Any, Any]:
        """
        ---------------------------------------------------------------------------

        Function that uses the protein fasta file to search for genes that encode
        proteins and then deletes all other entries in the gene location
        variable.

        ---------------------------------------------------------------------------

        Args:
            GeneLocation: DataFrame about the gene location information,
                            generated in PlantClusterFinder function.
            conversion: Conversion object between genes, transcripts, and proteins
            BiomG_G: Mapping between the Gene location file and the gene-IDs in the
                    conversion file.
            Protein_IDs: A cell array of all protein-IDs that are stored in the
                        protein fasta file. only active if
                        vRemoveNonProtLocations is 1. Else the array has no rows.

        Returns:
            Tuple containing the following values
            - **GeneLocation2**: DataFrame about the gene location information,
                            generated in PlantClusterFinder function.
            - **BiomG_G**: Mapping between the Gene location file and the gene-IDs in the
                    conversion file.
        """

        if not self._converted:
            log.warning("Implicitly requested adding conversion to Gene Location")
            self.remove_protein_transcript_info(conversion)

        original_gene_count = self._data.height
        if vProtein_IDs is None or len(vProtein_IDs) == 0:
            log.info("No protein IDs provided to remove, returning early.")
            return self.BiomG_G
        protein_ids = pl.DataFrame({"Gene Name": vProtein_IDs})
        self._data = self._data.join(protein_ids, on="Gene Name")
        log.debug(
            f"Gene Location data filtered from {original_gene_count} to {self._data.height} rows"
        )
        return self.BiomG_G

    @property
    def BiomG_G(self):
        """Get mapping dictionary between Gene Name and Protein Name"""
        return self._data.select(pl.col(["Gene Name", "Protein Name"])).rows_by_key(
            "Gene Name"
        )

    def find_intergenic_regions(self) -> pl.DataFrame:
        """
        ---------------------------------------------------------------------------

        Function computes all stretches of DNA that are not encoding at least one
        gene model.

        ---------------------------------------------------------------------------

        Returns:
            **InterChrom**: A DataFrame defining
                - chromosome,
                - start bp,
                - end bp,
                - gene before and
                - gene after a intergenic region.
        """

        if not self._sorted:
            log.warning("Implicitly Requested Gene Sorting")
            self.sort_genes()

        if "InterChrom" in self._cache:
            return self._cache["InterChrom"]

        ret: List[Dict[str, Any]] = []
        last_gene = None
        # Iterate through *sorted* rows and find each region where there is no gene located (intergenic region)
        # Compare the furthest out gene region in the chromosome to the nearest of the next gene
        # 1. if they overlap, the next gene becomes the furthest out gene
        # 2. if they don't overlap, an intergenic region is stored
        # TODO: This is a very slow and expensive way of doing this
        for chromosome, gene_start, gene_end, gene_name, strand in self._data.select(
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
            global_gene_start = min(gene_start, gene_end)
            if not (
                last_gene is None
                or last_gene[1] >= global_gene_start
                or chromosome != last_gene[0]
            ):
                ret.append(
                    {
                        "Chromosome Name": chromosome,
                        "Start": last_gene[1] + 1,
                        "End": global_gene_start - 1,
                        "Gene Before": last_gene[2],
                        "Gene After": gene_name,
                    }
                )
            last_gene = (chromosome, max(gene_start, gene_end), gene_name)

        self._cache["InterChrom"] = pl.DataFrame(ret)
        return self._cache["InterChrom"]

    def annotate_enzyme_info(
        self,
        MapG_Gpgdb: Dict[str, str],
        Gene2Rxn: GenesToReactions,
    ) -> Dict[str, bool]:
        """
        ---------------------------------------------------------------------------

        Generates a dictionary for all genes in the gene location
        information and looks up if this gene has a reaction in the pgdb.

        ---------------------------------------------------------------------------

        Args:
            MapG_Gpgdb: A mapping between Genes in the conversion file and genes
                        stored in the pgdb genes.dat file.
            Gene2Rxn: Mapping of all genes in the Genes_all to the Reactions in
                Reactions_all

        Returns:
            A dict describing which genes have at least one reaction
                    annotated in the pgdb (True if at least one reaction)
        """
        if not self._converted:
            raise Exception("Conversion data not loaded into gene locations")

        genes_with_reactions = [
            MapG_Gpgdb[x[0]][0] for x in Gene2Rxn.Gene2Rxn.items() if len(x[1]) > 0
        ]
        if not self._annotated:
            log.debug(
                f"Mapping {len(genes_with_reactions)} genes with reactions to {
                    len(self._data)} genes"
            )
            self._data = self._data.with_columns(
                has_reaction=pl.col("Protein Name").is_in(genes_with_reactions)
            )
            self._annotated = True
        return {
            k: list(v)[0]
            for k, v in self._data.select(pl.col(["Gene Name", "has_reaction"]))
            .rows_by_key("Gene Name")
            .items()
        }

    def introduce_hard_breaks(
        self,
        # vUniqueChrom: List[str],
        # vInterChrom: pl.DataFrame,
        vSeqGaps: pl.DataFrame,  # from annotate_sequencing_gaps()
        vMaxSeqGapSize: int = -1,  # from args
        vMaxInterGeneDist: int = -1,  # from args
        vPercentileForMaxInterGeneDist: Optional[float] = None,  # from args
        vMaxInterGeneDistByMedian: int = -1,  # from args
        # vBiomG_G: Dict[Any, Any],
        vHardBreakGaps: List[
            int
        ] = [],  # from args - specific lengths of gaps that indicate hard breaks
    ) -> None:
        """
        ---------------------------------------------------------------------------

        Function that prevents formation of clusters Sequencing gaps that
        exceed vMaxSeqGapSize or that are defined not to be crossed by
        vHardBreakGaps.

        ---------------------------------------------------------------------------

        Args:
            vGeneLocation: DataFrame about the gene location information,
                            generated in PlantClusterFinder function.
            vUniqueChrom: A list with unique entries of all chromosome and scaffold
                        names.
            vInterChrom: A DataFrame defining chromosome, start bp, end bp, gene
                        before and gene after a intergenic region.
            vSeqGaps: A string array (cell) DataFrame containing the
                    - chromosome/scaffold,
                    - the identification if it is coding (A) or masked (N) DNA
                    - start bp and
                    - length of a masked region.
                    In addition there is the information stored in with
                    intergenic region the gap is encoded and its start and end
                    bp, corrected for maximal size (= intergenic region).
            vMaxSeqGapSize: Maximal masked nucleotide region (in bp) that a cluster
                            is allowed to bridge (-1 or None is inactive).
            vMaxInterGeneDist: Maximal intergenic distance that still can be
                            gene cluster.
            vPercentileForMaxInterGeneDist: Prevents clusters to be spanned over
                                            the top n% of intergenic distances.
                                            (-1 or None is inactive)
            vMaxInterGeneDistByMedian: Maximal intergenic distance that still can
                                    be spanned by a gene cluster: this number
                                    times the median of the gene sizes.
            vHardBreakGaps: A list of numbers. The masked nucleotides sequences of
                            these specific lengths are preventing cluster formation
                            (Default is an empty list []).
        """
        vGeneLocation = self

        # this holds the SeqGaps_Hard data
        vSeqGaps = vSeqGaps.with_columns(
            hard_gap=pl.col("seq_gaps_length").is_in(vHardBreakGaps)
            | pl.col("seq_gaps_length")
            > vMaxSeqGapSize
        )
        if (
            vPercentileForMaxInterGeneDist is not None
            and vPercentileForMaxInterGeneDist != 1
        ):
            vMaxInterGeneDist_temp = (
                vGeneLocation.find_intergenic_regions()
                .select(len=pl.col("End") - pl.col("Start") + 1)
                .quantile(vPercentileForMaxInterGeneDist / 100)
            )["len"][0]
            log.info(
                "Percentile used for intergenic distance: %ibp",
                vPercentileForMaxInterGeneDist,
            )
            log.info(
                "Clusters broken when spanning distance larger: %ibp.\n",
                vMaxInterGeneDist_temp,
            )
            vMaxInterGeneDist = (
                vMaxInterGeneDist_temp
                if vMaxInterGeneDist is None or vMaxInterGeneDist == -1
                else min(vMaxInterGeneDist, vMaxInterGeneDist_temp)
            )
        if vMaxInterGeneDistByMedian is None or vMaxInterGeneDistByMedian != -1:
            vMaxInterGeneDist_temp = vGeneLocation.data.select(
                len=pl.col("Gene End (bp)") - pl.col("Gene Start (bp)")
            ).quantile(0.50)["len"][0]
            vMaxInterGeneDist = (
                vMaxInterGeneDist_temp * vMaxInterGeneDistByMedian
                if vMaxInterGeneDist is None or vMaxInterGeneDist == -1
                else min(
                    vMaxInterGeneDist,
                    vMaxInterGeneDist_temp * vMaxInterGeneDistByMedian,
                )
            )

        hard_breaks = (
            vGeneLocation.find_intergenic_regions()
            .join_where(  # using join_where as a cheeky filter
                vSeqGaps,
                pl.col("Start").is_between(
                    pl.col("start index"),
                    pl.col("end index"),
                    # | pl.col("End").is_between(
                    # pl.col("start index"), pl.col("end index")
                    # )
                    # | ( pl.lit(vMaxInterGeneDist != -1) & pl.col("Length") > vMaxInterGeneDist)
                ),
            )
            .filter(pl.col("hard_gap"))
            .select(
                pl.col(["Chromosome Name", "Start", "End", "Gene Before", "Gene After"])
            )
            .with_columns(pseudochromosome=pl.col("Chromosome Name").cum_count())
        )
        prev_end = 0
        pseudochromosome = []
        if len(hard_breaks) > 0:
            for h_break in hard_breaks.iter_rows():
                pseudochromosome.append(
                    {
                        "break_start": prev_end,
                        "break_end": h_break[1],
                        "pchromosome name": f"{h_break[0]}_{h_break[5]-1}",
                    }
                )
                prev_end = h_break[2]
            pseudochrom_df = pl.DataFrame(pseudochromosome)
            self._data = self._data.join_where(
                pseudochrom_df,
                pl.col("Start Location (bp)") >= pl.col("break_start"),
                pl.col("End Location (bp)") < pl.col("break_end"),
            ).with_columns(
                pl.col("pchromosome name").alias("Chromosome Name"),
                chromosome_name_old=pl.col("Chromosome Name"),
            )

    def insert_genes(self, other_genes: pl.DataFrame) -> None:
        _ = self._data.extend(other_genes)  # modifies in place and returns DataFrame


# Example usage
if __name__ == "__main__":
    gene_loc = GeneLoc(
        file_path=Path(
            "./csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt"
        ),
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
