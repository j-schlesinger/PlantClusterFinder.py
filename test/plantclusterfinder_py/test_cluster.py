def test_read_in_conversion_file():
    """The conversion file defines the genes, proteins, and reactions
       that many downstream processes rely on."""

    from plantclusterfinder_py.cluster import read_in_conversion_file
    conversion = read_in_conversion_file(
        './csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt')
    assert (len(conversion.MapG) == 9629)
    assert (len(conversion.MapT) == 9629)
    assert (len(conversion.MapP) == 9629)


def test_label_reactions_with_sig_and_tail():
    test_pgdb_file = 'csubellipsoidea/pgdb'
    reactions_dat = test_pgdb_file+'/csubellipsoideacyc/1.0/data/reactions.dat'
    test_signature_enzymes_tsv = 'Inputs/scaffold-tailoring-reactions-05082016.tab'
    pass


def test_compute_LocalDuplication():
    pass


def test_introduceHypoGenes():
    pass


def test_introduce_HardB():
    pass


def test_annotate_Sequencing_Gaps():
    pass


def test_get_Sequencing_Gaps():
    pass


def test_run_awk_skripts():
    pass


def test_get_metabolic_domains():
    pass


def test_read_MCL_clustering_File():
    pass


def test_get_MCL_clustering():
    pass


def test_map_pgdb_reactions_to_genes():
    pass


def test_match_pgdb_flat():
    pass


def test_matchPGDBs_to_MDs():
    pass


def test_map_PGDBs_to_conversion_file():
    pass


def test_remove_protein_transcript_info_from_genepositionfile():
    pass


def test_consolidate_gene_info():
    pass


def test_sort_genes():
    pass


def test_find_intergenic_regions():
    pass


def test_map_genelocation_to_conversion_file():
    pass


def test_annotate_enyme_info_to_gene_location():
    pass


def test_perform_clustering():
    pass


def test_compute_background():
    pass


def test_perform_cutoff():
    pass


def test_write_Gene_Output():
    pass


def test_write_Cluster_Output():
    pass


def test_extract_results_with_header():
    pass


def test_analyze_PlantClusterGapFile():
    pass


def test_read_in_PGDB_Flat():
    pass


def test_Remove_non_protein_fasta_gene_locations():
    pass
