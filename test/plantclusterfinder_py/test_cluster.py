import polars as pl

# Part 1


def test_read_in_conversion_file():
    """The conversion file defines the genes, proteins, and reactions
    that many downstream processes rely on."""

    from plantclusterfinder_py.cluster import read_in_conversion_file

    conversion = read_in_conversion_file(
        "./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
    )
    assert len(conversion.MapG) == 9629
    assert len(conversion.MapT) == 9629
    assert len(conversion.MapP) == 9629


# Part 2
def test_read_in_PGDB_Flat():
    """Used by map_pgdb_reactions_to_genes"""

    pgdb_dir = "./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data"
    from plantclusterfinder_py.cluster import read_in_PGDB_flat

    gene_df = read_in_PGDB_flat(
        pgdb_dir + "/genes.dat", ("UNIQUE-ID", "PRODUCT", "ACCESSION-1")
    )
    assert isinstance(gene_df, pl.DataFrame)
    assert len(gene_df) == 2053
    g4s_3028 = gene_df.filter(pl.col("PRODUCT") == ("G4S-3028-MONOMER"))
    assert len(g4s_3028) == 1


def test_match_pgdb_flat():
    """Used by map_pgdb_reactions_to_genes"""
    assert 1 == 0  # TEST INCOMPLETE


def test_map_pgdb_reactions_to_genes():
    assert 1 == 0  # TEST INCOMPLETE


# Part 3
def test_get_metabolic_domains():
    assert 1 == 0  # TEST INCOMPLETE


def test_matchPGDBs_to_MDs():
    assert 1 == 0  # TEST INCOMPLETE


# part 4
def test_label_reactions_with_sig_and_tail():
    pgdb_dir = "./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data"
    reactions_dat = pgdb_dir + "/reactions.dat"
    test_signature_enzymes_tsv = "Inputs/scaffold-tailoring-reactions-05082016.tab"
    from plantclusterfinder_py.cluster import label_reactions_with_sig_and_tail

    labeled_reaction = label_reactions_with_sig_and_tail(
        reactions_dat, test_signature_enzymes_tsv
    )
    assert 1 == 0  # TEST IMCOMPLETE


# part 5
def test_map_PGDBs_to_conversion_file():
    assert 1 == 0  # TEST INCOMPLETE


# part 6
def test_get_MCL_clustering():
    assert 1 == 0  # TEST INCOMPLETE


def test_read_MCL_clustering_File():
    assert 1 == 0  # TEST INCOMPLETE


# part 7
def test_extract_results_with_header():
    assert 1 == 0  # TEST INCOMPLETE


def test_remove_protein_transcript_info_from_genepositionfile():
    assert 1 == 0  # TEST INCOMPLETE


def test_consolidate_gene_info():
    assert 1 == 0  # TEST INCOMPLETE


def test_sort_genes():
    assert 1 == 0  # TEST INCOMPLETE


# Part 8
def test_map_genelocation_to_conversion_file():
    assert 1 == 0  # TEST INCOMPLETE


def test_remove_non_protein_fasta_gene_locations():
    assert 1 == 0  # TEST INCOMPLETE


# Part 9
def test_find_intergenic_regions():
    assert 1 == 0  # TEST INCOMPLETE


# Part 10: No Cluster Functions


# part 11
def test_analyze_plant_cluster_gap_file():
    """Used by get_sequencing_gaps"""
    assert 1 == 0  # TEST INCOMPLETE


def test_run_awk_skripts():
    """Used by get_sequencing_gaps"""
    assert 1 == 0  # TEST INCOMPLETE


def test_get_sequencing_gaps():
    """Used by annotate_sequencing_gaps"""
    assert 1 == 0  # TEST INCOMPLETE


def test_annotate_sequencing_gaps():
    assert 1 == 0  # TEST INCOMPLETE


# part 12


def test_introduce_hard_breaks():
    assert 1 == 0  # TEST INCOMPLETE


# part 13
def test_introduce_hypo_genes():
    assert 1 == 0  # TEST INCOMPLETE


# Part 14
def test_annotate_enyme_info_to_gene_location():
    assert 1 == 0  # TEST INCOMPLETE


# Part 15
def test_compute_local_duplication():
    assert 1 == 0  # TEST INCOMPLETE


# Part 16
def test_compute_background():
    """Used by perform_clustering"""
    assert 1 == 0  # TEST INCOMPLETE


def test_perform_clustering():
    assert 1 == 0  # TEST INCOMPLETE


# Part 17
def test_perform_cutoff():
    assert 1 == 0  # TEST INCOMPLETE


# Part 18
def test_write_gene_output():
    assert 1 == 0  # TEST INCOMPLETE


# Part 19
def test_write_cluster_output():
    assert 1 == 0  # TEST INCOMPLETE
