from pathlib import Path

import pytest


@pytest.fixture
def genelocation_obj():
    """Initialized GeneLocation object for testing."""
    from plantclusterfinder.genelocation import GeneLoc

    return GeneLoc(
        file_path=Path("csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt")
    )


@pytest.fixture
def conversion_obj():
    """Conversion Object."""
    from plantclusterfinder.conversion import Conversion

    pgdb_dir = Path("./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data")

    conversion = Conversion.read_in_conversion_file(
        Path("./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt")
    )
    conversion.add_pgdb_id_column(pgdb_dir)
    return conversion


@pytest.fixture
def t_gene_loc():
    """Initialized GeneLocation object with simple test data for testing."""
    from plantclusterfinder.genelocation import GeneLoc

    return GeneLoc(file_path=Path("tests/data/gene_location_test.tsv"))


@pytest.fixture
def expected_mcl_cluster_obj(conversion_obj):
    """Expected MCL cluster file."""
    from plantclusterfinder.cluster import read_mcl_clustering_file

    return read_mcl_clustering_file(
        Path(
            "./csubellipsoidea/MCL_clustering_expected_results/dump.CsubellipsoideaC_169_227_v2.0.protein.pcf13.mci.I20"
        ),
        Path("./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa"),
        conversion_obj,
        False,
    )


@pytest.fixture(scope="session")
def genelocation_transformed_obj():
    """Initialized and Transformed GeneLocation object for testing.

    This shouldn't need further transformation,
    so it is only created once per session
    """
    from plantclusterfinder.cluster import (
        map_pgdb_gene_ids_to_conversion_file,
        map_pgdb_reactions_to_genes,
        read_in_pgdb_flat,
        read_mcl_clustering_file,
    )
    from plantclusterfinder.conversion import Conversion
    from plantclusterfinder.genelocation import GeneLoc

    pgdb_dir = Path("./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data")
    gene_location = GeneLoc(
        file_path=Path("csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt")
    )
    conversion = Conversion.read_in_conversion_file(
        Path("./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt")
    )
    conversion.add_pgdb_id_column(pgdb_dir)
    gene_location.add_conversion_columns(conversion)
    gene_location.consolidate_gene_info()
    gene_location.sort_genes()
    genes_all = read_in_pgdb_flat(
        pgdb_dir / "genes.dat",
        ("UNIQUE-ID", "ACCESSION-1"),
        unique_keys=("UNIQUE-ID", "ACCESSION-1"),
    )
    gene_info_from_genes_dat = ["ACCESSION-1"]
    ids_to_map = "TP"  # would be from args
    pgdb_conversions_df = map_pgdb_gene_ids_to_conversion_file(
        conversion, genes_all, gene_info_from_genes_dat, ids_to_map
    )
    gene2rxn = map_pgdb_reactions_to_genes(pgdb_dir, conversion)
    gene_location.annotate_enzyme_info(pgdb_conversions_df, gene2rxn)
    mcl_cluster_data = read_mcl_clustering_file(
        Path(
            "./csubellipsoidea/MCL_clustering_expected_results/dump.CsubellipsoideaC_169_227_v2.0.protein.pcf13.mci.I20"
        ),
        Path("./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa"),
        conversion,
        False,  # do not remove proteins
    )
    gene_location.compute_local_duplication(
        mcl_cluster_data,
        100_000,  # kb for local duplication
        10,  # genes for local duplication
    )
    return gene_location
