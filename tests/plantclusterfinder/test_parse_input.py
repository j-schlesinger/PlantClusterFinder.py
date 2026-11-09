from tests.plantclusterfinder.assets import (
    CSUB_GFF3_PATH,
    GENE_COUNT,
    GENOME_FASTA_PATH,
)


def test_extract_genes_from_gff_and_fasta():
    from plantclusterfinder.parse_input import extract_genes_from_gff_and_fasta

    (gtpf_file, glof_file, faa_fasta_file) = extract_genes_from_gff_and_fasta(
        GENOME_FASTA_PATH, CSUB_GFF3_PATH, "Name"
    )

    with gtpf_file.open() as f:
        assert len(f.readlines()) == GENE_COUNT + 1
    with glof_file.open() as f:
        assert len(f.readlines()) == GENE_COUNT + 1
    with faa_fasta_file.open() as f:
        x = sum([1 for x in f if x.startswith(">")])
        assert x == GENE_COUNT


def test_extract_genes_gffutils():
    from plantclusterfinder.parse_input import extract_genes_gffutils

    (gtpf_deque, glof_deque, faa_fasta_deque) = extract_genes_gffutils(
        CSUB_GFF3_PATH, GENOME_FASTA_PATH, "Name", "Name", ".*"
    )
    assert len(gtpf_deque) == GENE_COUNT
    assert len(glof_deque) == GENE_COUNT
    assert len(faa_fasta_deque) == GENE_COUNT


def extract_genes_bcbio():
    from plantclusterfinder.parse_input import extract_genes_bcbio

    (gtpf_deque, glof_deque, faa_fasta_deque) = extract_genes_bcbio(
        CSUB_GFF3_PATH, GENOME_FASTA_PATH, "Name", "Name", ".*"
    )
    assert len(gtpf_deque) == GENE_COUNT
    assert len(glof_deque) == GENE_COUNT
    assert len(faa_fasta_deque) == GENE_COUNT
