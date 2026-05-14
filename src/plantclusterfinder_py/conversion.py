"""Handle the conversion join between PGDB and GFF3 data."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import polars as pl
from pydantic import BaseModel

from plantclusterfinder_py.pgdb_utils import read_in_pgdb_flat

_log = logging.getLogger(__name__)


class Conversion(BaseModel):
    """Conversion maps between gene, transcript, and protein IDs."""

    model_config = {"arbitrary_types_allowed": True}
    """DataFrame containing Gene Name, Transcript Name, and Protein Name"""
    df: pl.DataFrame
    has_pgdb_ids: bool = False
    _cache: dict[str, Any] = {}

    @classmethod
    def read_in_conversion_file(
        cls,
        gene_transcript_protein_mapping_filename: Path,
    ) -> Conversion:
        """---------------------------------------------------------------------------

        Function that reads in the gene, transcript, and protein ID conversion
        file and stores the IDs as unique lists and mapping between genes and
        transcripts aswell as genes and proteins.

        ---------------------------------------------------------------------------

        Args:
            gene_transcript_protein_mapping_filename: TSV file with a header describing
                                                gene, transcript, and protein
                                                identifier and then for all the
                                                identfiers the listings (gene-ID,
                                                transcript-ID, protein-ID).

        Returns:
            Conversion object containing the IDs from the conversion file
        Raises

        """
        _log.debug("Reading in conversion file %s", gene_transcript_protein_mapping_filename)
        conversion_data = pl.read_csv(
            gene_transcript_protein_mapping_filename,
            separator="\t",
            schema_overrides=dict.fromkeys(
                ["Gene Name", "Transcript Name", "Protein Name"], pl.String
            ),
        )
        column_names = conversion_data.columns
        conversion_len = len(conversion_data)
        _log.debug("Read in conversion data, %d rows, columns are %s", conversion_len, column_names)
        conversion = Conversion(df=conversion_data)
        # Validate conversion before output
        if len(conversion.map_gene_to_transcript) < len(conversion.map_gene):
            _log.warning("Some Transcripts match multiple genes in the ID Conversion file")
            transcript_multi = (
                conversion_data.group_by(column_names[1])
                .n_unique()
                .filter(pl.col(column_names[0]) > 1)[column_names[1]]
            )
            for transcript in transcript_multi:
                _log.warning("Please check the transcript %s", transcript)
            _log.error("Please clean the protein ID file")
        if len(conversion.map_gene_to_protein) < len(conversion.map_gene):
            _log.warning("Some Proteins match multiple genes in the ID Conversion file")
            prot_multi = (
                conversion_data.group_by(column_names[2])
                .n_unique()
                .filter(pl.col(column_names[0]) > 1)[column_names[2]]
            )
            for prot in prot_multi:
                _log.warning("Please check the protein %s", prot)
            _log.error("Please clean the protein ID file")
        _log.info(
            "Number of Gene-ID to Transcript-ID mappings: %i", len(conversion.map_gene_to_transcript)
        )
        _log.info("Number of Gene-ID to Protein mappings: %i", len(conversion.map_gene_to_protein))
        return conversion

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
        if len(pg_ids.filter(pl.col("pgdb_gene_id").is_not_null())) > 0:
            self.df = self.df.join(pg_ids, on="Gene Name", how="left").unique()
            return
        _log.warning("Could not map PGDB Gene IDs using proteins.dat, trying genes.dat")
        # I already tested against UNIQUE-ID and PRODUCT of genes.dat
        # as GENE and UNIQUE-ID from proteins.dat.
        # Now I just need to see if my gene or protein name
        # line up with the COMMON-NAME from genes.dat
        genes_table = read_in_pgdb_flat(
            pgdb_dir / "genes.dat", attrib=("UNIQUE-ID", "COMMON-NAME", "ACCESSION-1")
        ).explode("COMMON-NAME")
        _log.debug("Genes table %s", genes_table)
        gg_common_prot = self.df.join(
            genes_table, left_on="Protein Name", right_on="COMMON-NAME"
        ).select(pl.col("Gene Name"), pgdb_gene_id1=pl.col("UNIQUE-ID"))
        gg_common_gene = self.df.join(
            genes_table, left_on="Gene Name", right_on="COMMON-NAME"
        ).select(pl.col("Gene Name"), pgdb_gene_id2=pl.col("UNIQUE-ID"))

        gg_ids = pl.concat([gg_common_prot, gg_common_gene], how="align")
        gg_ids = gg_ids.select(
            pl.col("Gene Name"),
            pl.coalesce(
                (
                    pl.col("pgdb_gene_id1"),
                    pl.col("pgdb_gene_id2"),
                )
            ).alias("pgdb_gene_id"),
        )
        if len(gg_ids.filter(pl.col("pgdb_gene_id").is_not_null())) > 0:
            self.df = self.df.join(gg_ids, on="Gene Name", how="left").unique()
            return
        _log.error("Could not map gene IDs to protein IDs %s", gg_ids)
        raise Exception("Cannot map PGDB Gene IDs to IDs in GLOF file.")

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
