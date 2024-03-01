# Part 1
def test_read_in_conversion_file():
    """The conversion file defines the genes, proteins, and reactions
       that many downstream processes rely on."""

    from plantclusterfinder_py.cluster import read_in_conversion_file
    conversion = read_in_conversion_file(
        './csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt')
    assert (len(conversion.MapG) == 9629)
    assert (len(conversion.MapT) == 9629)
    assert (len(conversion.MapP) == 9629)


# Part 2
def test_match_pgdb_flat():
    """Used by map_pgdb_reactions_to_genes"""
    pass


def test_read_in_PGDB_Flat():
    """Used by map_pgdb_reactions_to_genes"""
    pass


def test_map_pgdb_reactions_to_genes():
    pass


# Part 3
def test_get_metabolic_domains():
    pass


def test_matchPGDBs_to_MDs():
    pass


# part 4
def test_label_reactions_with_sig_and_tail():
    test_pgdb_file = 'csubellipsoidea/pgdb'
    reactions_dat = test_pgdb_file+'/csubellipsoideacyc/1.0/data/reactions.dat'
    test_signature_enzymes_tsv = 'Inputs/scaffold-tailoring-reactions-05082016.tab'
    pass


# part 5
def test_map_PGDBs_to_conversion_file():
    pass


# part 6
def test_get_MCL_clustering():
    pass


def test_read_MCL_clustering_File():
    pass


# part 7
def test_extract_results_with_header():
    pass


def test_remove_protein_transcript_info_from_genepositionfile():
    pass


def test_consolidate_gene_info():
    pass


def test_sort_genes():
    pass


# Part 8
def test_map_genelocation_to_conversion_file():
    pass


def test_remove_non_protein_fasta_gene_locations():
    pass


# Part 9
def test_find_intergenic_regions():
    pass


# Part 10: No Cluster Functions


# part 11
def test_analyze_plant_cluster_gap_file():
    """Used by get_sequencing_gaps"""
    pass


def test_run_awk_skripts():
    """Used by get_sequencing_gaps"""
    pass


def test_get_sequencing_gaps():
    """Used by annotate_sequencing_gaps"""
    pass


def test_annotate_sequencing_gaps():
    pass

# part 12


def test_introduce_hard_breaks():
    pass


# part 13
def test_introduce_hypo_genes():
    pass


# Part 14
def test_annotate_enyme_info_to_gene_location():
    pass


# Part 15
def test_compute_local_duplication():
    pass


# Part 16
def test_compute_background():
    """Used by perform_clustering"""
    pass


def test_perform_clustering():
    pass


# Part 17
def test_perform_cutoff():
    pass


# Part 18
def test_write_gene_output():
    pass


# Part 19
def test_write_cluster_output():
    pass
