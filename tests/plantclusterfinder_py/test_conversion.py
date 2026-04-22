import polars as pl
import pytest

from tests.plantclusterfinder_py.assets import GENE_COUNT, GTPF_FILE_PATH, PGDB_DIR


def test_read_in_conversion_file():
    """Test read_in_conversion_file.

    The conversion file defines the genes, proteins, and reactions
    that many downstream processes rely on.
    """
    from plantclusterfinder_py.conversion import Conversion

    conversion = Conversion.read_in_conversion_file(GTPF_FILE_PATH)
    assert len(conversion.map_gene) == GENE_COUNT
    assert len(conversion.map_transcript) == GENE_COUNT
    assert len(conversion.map_protein) == GENE_COUNT


@pytest.mark.benchmark
def test_add_pgdb_id_column():
    """Test add_pgdb_id_column from Conversion."""
    from plantclusterfinder_py.conversion import Conversion

    conversion1 = Conversion.read_in_conversion_file(GTPF_FILE_PATH)
    conversion2 = Conversion.read_in_conversion_file(GTPF_FILE_PATH)
    conversion2.add_pgdb_id_column(PGDB_DIR)

    assert len(conversion1.df) == len(conversion2.df)
    assert len(conversion2.df.filter(pl.col("Gene Name").is_null())) == 0
    assert len(conversion2.df.filter(pl.col("Transcript Name").is_null())) == 0
    assert len(conversion2.df.filter(pl.col("Protein Name").is_null())) == 0
    assert len(conversion2.df.filter(pl.col("pgdb_gene_id").is_null())) == 7576
