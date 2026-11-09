"""Constants used multiple places in testing."""

from pathlib import Path

GENE_COUNT = 9_629
PGDB_GENE_COUNT = 2053
PGDB_DIR = Path("./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data")
CSUB_GFF3_PATH = Path("./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.gene.gff3")
MCL_EXPECTED_RESULTS_PATH = Path(
    "./csubellipsoidea/MCL_clustering_expected_results/dump.CsubellipsoideaC_169_227_v2.0.protein.pcf13.mci.I20"
)

PROTEIN_FASTA_PATH = Path("./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa")
GTPF_FILE_PATH = Path("./csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt")
GENOME_FASTA_PATH = Path("./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.hardmasked.fa")
GLOF_FILE_PATH = Path("csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt")
GAPANALYSIS_EXPECTED_RESULTS = Path(
    "./csubellipsoidea/Gapanalysis_expected_results/CsubellipsoideaC_169_227_v2.0.hardmasked.fa_GAPOutput"
)
