import random
from pathlib import Path

import polars as pl
import pytest

from tests.plantclusterfinder_py.assets import (
    GAPANALYSIS_EXPECTED_RESULTS,
    GENE_COUNT,
    GLOF_FILE_PATH,
    MCL_EXPECTED_RESULTS_PATH,
    PGDB_DIR,
    PROTEIN_FASTA_PATH,
)


def test_init_genelocation():
    from plantclusterfinder_py.genelocation import GeneLoc

    gene_location = GeneLoc(file_path=GLOF_FILE_PATH)

    assert len(gene_location.data.columns) == 6
    assert len(gene_location) == 9_629


def test_read_in_conversion_to_genepositionfile(genelocation_obj, t_gene_loc, conversion_obj):
    from plantclusterfinder_py.cluster import read_in_conversion_file

    gene_location = genelocation_obj
    conversion = conversion_obj
    expected_columns = [
        "Gene Name",
        "Gene Start (bp)",
        "Gene End (bp)",
        "Chromosome Name",
        "Strand",
        "Hypothetical",
        "Transcript Name",  # added
        "Protein Name",  # added
        "pgdb_gene_id",  # added
    ]
    gene_location.add_conversion_columns(conversion)
    assert len(gene_location.columns) == 9
    assert expected_columns == gene_location.columns
    assert len(gene_location) == 9_629

    test_gene_location = t_gene_loc
    bad_conversion = read_in_conversion_file(Path("./tests/data/bad_conversion.tsv"))
    bad_conversion.df = bad_conversion.df.with_columns(pgdb_gene_id=pl.lit(None))
    bad_conversion.has_pgdb_ids = True
    with pytest.raises(Exception, match="Multiple entries") as _:
        test_gene_location.add_conversion_columns(bad_conversion)


def test_consolidate_gene_info(genelocation_obj, t_gene_loc):
    gene_location = genelocation_obj
    gene_info = gene_location.consolidate_gene_info()
    assert len(gene_info) == 45
    test_gene_location = t_gene_loc
    test_gene_info = test_gene_location.consolidate_gene_info()
    assert test_gene_info == ["A"]


def test_sort_genes(genelocation_obj, t_gene_loc):
    import polars as pl

    gene_location = genelocation_obj
    gene_location.sort_genes()

    assert len(gene_location) == 9629

    test_gene_location = t_gene_loc
    test_gene_location.sort_genes()
    sorted_genes = test_gene_location.data.select(pl.col("Gene Name"))
    expected_gene_order = pl.DataFrame({"Gene Name": ["GENE1", "GENE2", "GENE3", "GENE4"]})
    unexpected_gene_order = pl.DataFrame({"Gene Name": ["GENE2", "GENE1", "GENE4", "GENE3"]})
    assert sorted_genes.equals(expected_gene_order)
    assert not sorted_genes.equals(unexpected_gene_order)


def test_map_genelocation_to_conversion_file(genelocation_obj, conversion_obj):
    # Gene Location Test, possibly redundant (mapping takes place in remove_protein_transcript_info)

    conversion = conversion_obj
    expected_columns = [
        "Gene Name",
        "Gene Start (bp)",
        "Gene End (bp)",
        "Chromosome Name",
        "Strand",
        "Hypothetical",
        "Transcript Name",  # added
        "Protein Name",  # added
        "pgdb_gene_id",  # added
    ]
    gene_location = genelocation_obj
    genelocation_obj.add_conversion_columns(conversion)
    biomg_g = gene_location.map_to_conversion_file()
    assert expected_columns == gene_location.data.columns
    assert len(biomg_g) == 9629


def test_remove_non_protein_fasta_gene_locations(genelocation_obj, conversion_obj):
    conversion = conversion_obj
    gene_location = genelocation_obj
    gene_location.add_conversion_columns(conversion)
    biomg_g = gene_location.map_to_conversion_file()
    ok_prot = list(biomg_g.keys())
    random.shuffle(ok_prot)
    protein_ids = ok_prot[0:100]
    cleaned_gene_location = gene_location.remove_non_protein_fasta_gene_locations(
        conversion, protein_ids
    )
    unclean_gene_location = gene_location.remove_non_protein_fasta_gene_locations(conversion, [])
    assert len(gene_location) == len(unclean_gene_location)
    assert len(cleaned_gene_location) == 100


# Part 9
def test_find_intergenic_regions(genelocation_obj, t_gene_loc):
    import polars as pl

    # Gene Location test
    gene_location = genelocation_obj
    intergenic_regions = gene_location.find_intergenic_regions()
    assert len(intergenic_regions) == GENE_COUNT - 216

    intergenic_regions_2 = t_gene_loc.find_intergenic_regions()
    t_gene_result_expected = pl.DataFrame(
        {
            "Chromosome Name": ["A", "A"],
            "Start": [101, 401],
            "End": [149, 499],
            "Gene Before": ["GENE1", "GENE3"],
            "Gene After": ["GENE2", "GENE4"],
        }
    )
    assert intergenic_regions_2.equals(t_gene_result_expected)


# Part 14
def test_annotate_enzyme_info(genelocation_obj, conversion_obj):
    from plantclusterfinder_py.cluster import (
        map_pgdb_gene_ids_to_conversion_file,
        map_pgdb_reactions_to_genes,
        read_in_pgdb_flat,
    )

    conversion = conversion_obj
    genes_all = read_in_pgdb_flat(
        PGDB_DIR / "genes.dat",
        ("UNIQUE-ID", "ACCESSION-1", "COMMON-NAME"),
        unique_keys=("UNIQUE-ID", "ACCESSION-1"),
    )
    gene_info_from_genes_dat = ["ACCESSION-1", "COMMON-NAME"]
    ids_to_map = "TP"
    pgdb_conversions_df = map_pgdb_gene_ids_to_conversion_file(
        conversion, genes_all, gene_info_from_genes_dat, ids_to_map
    )
    print("PGDB Conversions DataFrame")
    print(pgdb_conversions_df)
    gene_location = genelocation_obj
    gene_location.add_conversion_columns(conversion)
    gene2rxn = map_pgdb_reactions_to_genes(PGDB_DIR, conversion)
    assert not gene_location.is_annotated
    genes_with_reactions = gene_location.annotate_enzyme_info(pgdb_conversions_df, gene2rxn)
    expected_columns = [
        "Gene Name",
        "Gene Start (bp)",
        "Gene End (bp)",
        "Chromosome Name",
        "Strand",
        "Hypothetical",
        "Transcript Name",
        "Protein Name",
        "pgdb_gene_id",  # added
        "pgdb_reaction_id",  # added
        "has_reaction",  # added
        "EC-NUMBER",  # added
    ]
    assert gene_location.is_annotated
    assert expected_columns == gene_location.data.columns
    assert len(genes_with_reactions) == 9629
    assert len(list(filter(lambda x: x, genes_with_reactions.values()))) == 2053


def test_introduce_hard_breaks(genelocation_obj):
    from plantclusterfinder_py.cluster import annotate_sequencing_gaps

    sequencing_gaps_file = GAPANALYSIS_EXPECTED_RESULTS
    sequencing_gaps_fasta = sequencing_gaps_file.with_suffix(
        sequencing_gaps_file.suffix.split("_")[0]
    )
    gene_location = genelocation_obj
    uniq_chr = list(gene_location.data.select("Chromosome Name").unique().to_dict().values())[
        0
    ].to_list()
    sequencing_gaps = annotate_sequencing_gaps(
        sequencing_gaps_fasta, gene_location.find_intergenic_regions(), uniq_chr
    )
    gene_location.introduce_hard_breaks(
        sequencing_gaps,
        max_sequence_gap_size=int(4e4),
        percentile_for_max_inter_gene_distance=80,
    )
    print(gene_location.data.select(pl.col("Chromosome Name").unique()))
    assert gene_location.data.select(pl.col("Chromosome Name").unique()).height >= len(uniq_chr)
    assert gene_location.data.select(pl.col("Chromosome Name").unique()).height >= 45
    assert gene_location.data.select(pl.col("Chromosome Name").unique()).height <= 51


def test_small_gap_introduce_hard_breaks(genelocation_obj):
    from plantclusterfinder_py.cluster import annotate_sequencing_gaps

    sequencing_gaps_file = GAPANALYSIS_EXPECTED_RESULTS
    sequencing_gaps_fasta = sequencing_gaps_file.with_suffix(
        sequencing_gaps_file.suffix.split("_")[0]
    )
    gene_location = genelocation_obj
    uniq_chr = list(gene_location.data.select("Chromosome Name").unique().to_dict().values())[
        0
    ].to_list()
    sequencing_gaps = annotate_sequencing_gaps(
        sequencing_gaps_fasta, gene_location.find_intergenic_regions(), uniq_chr
    )
    gene_location = genelocation_obj
    gene_location.introduce_hard_breaks(
        sequencing_gaps,
        max_sequence_gap_size=2,
        percentile_for_max_inter_gene_distance=80,
    )
    # unpredictable results, probably should be examined
    assert gene_location.data.select(pl.col("Chromosome Name").unique()).height in [46, 144, 157]


def test_insert_genes(genelocation_obj):
    print(genelocation_obj.data.head(10))
    test_genes = [
        (f"TEST{i + 1:02}", 50 + 1000 * i, 50 + 1000 * i, "test_scaffold", "-1", True)
        for i in range(500)
    ]
    gene_location = genelocation_obj.model_copy()
    orig_gene_count = len(gene_location)
    test_gene_df = pl.DataFrame(test_genes, schema=gene_location.data.schema, orient="row")
    print(test_gene_df.head(10))
    gene_location.insert_genes(test_gene_df)
    assert orig_gene_count + 500 == len(gene_location)


def test_compute_local_duplication(genelocation_obj, conversion_obj):
    from plantclusterfinder_py.cluster import read_mcl_clustering_file

    conversion = conversion_obj

    mcl_cluster_data = read_mcl_clustering_file(
        MCL_EXPECTED_RESULTS_PATH,
        PROTEIN_FASTA_PATH,
        conversion,
        False,  # do not remove proteins
    )
    gene_location = genelocation_obj
    gene_location.add_conversion_columns(conversion)
    # find duplicates within 100KB and 10 genes of each gene (defaults)
    local_duplication = gene_location.compute_local_duplication(
        mcl_cluster_data,
        100000,
        10,
    )
    expected_columns = [
        "Gene Name",
        "Gene Start (bp)",
        "Gene End (bp)",
        "Chromosome Name",
        "Strand",
        "Hypothetical",
        "Transcript Name",
        "Protein Name",
        "pgdb_gene_id",
        # "neighbor_genes",
        "mcl_clusters",
        "is_locally_duplicated",
    ]
    print(local_duplication.columns)
    print(local_duplication["mcl_clusters"].head())
    print(local_duplication.filter(pl.col("mcl_clusters").list.len() > 1).height)
    assert expected_columns == local_duplication.columns
    assert local_duplication.filter(pl.col("is_locally_duplicated")).height == 358
    # assert 91 == local_duplication.filter(pl.col("is_locally_duplicated")).height
