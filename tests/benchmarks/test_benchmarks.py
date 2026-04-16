"""Performance benchmarks for PlantClusterFinder core operations."""

from collections import namedtuple
from pathlib import Path

import polars as pl
import pytest


GLOF_FILE_PATH = Path("csubellipsoidea/glof_CsubellipsoideaC_169_227_v2.0.gene.gff3.txt")
GTPF_FILE_PATH = Path(
    "csubellipsoidea/gtpf_CsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt"
)
PGDB_DIR = Path("./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data")
MCL_EXPECTED_RESULTS_PATH = Path(
    "./csubellipsoidea/MCL_clustering_expected_results/"
    "dump.CsubellipsoideaC_169_227_v2.0.protein.pcf13.mci.I20"
)
PROTEIN_FASTA_PATH = Path(
    "./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa"
)


# --- GeneLoc benchmarks ---


def test_bench_geneloc_init(benchmark):
    """Benchmark loading gene location data from TSV."""
    from plantclusterfinder_py.genelocation import GeneLoc

    benchmark(GeneLoc, file_path=GLOF_FILE_PATH)


def test_bench_geneloc_sort(benchmark):
    """Benchmark sorting genes by chromosome and position."""
    from plantclusterfinder_py.genelocation import GeneLoc

    gene_loc = GeneLoc(file_path=GLOF_FILE_PATH)

    def sort_genes():
        gene_loc.is_sorted = False
        gene_loc.sort_genes()

    benchmark(sort_genes)


def test_bench_geneloc_find_intergenic_regions(benchmark):
    """Benchmark finding intergenic regions across the genome."""
    from plantclusterfinder_py.genelocation import GeneLoc

    gene_loc = GeneLoc(file_path=GLOF_FILE_PATH)
    gene_loc.sort_genes()

    def find_regions():
        gene_loc._cache.pop("InterChrom", None)
        return gene_loc.find_intergenic_regions()

    benchmark(find_regions)


def test_bench_geneloc_consolidate(benchmark):
    """Benchmark consolidating gene info (deduplication)."""
    from plantclusterfinder_py.genelocation import GeneLoc

    gene_loc = GeneLoc(file_path=GLOF_FILE_PATH)
    benchmark(gene_loc.consolidate_gene_info)


# --- Conversion and PGDB benchmarks ---


def test_bench_read_conversion_file(benchmark):
    """Benchmark reading the gene-transcript-protein conversion file."""
    from plantclusterfinder_py.cluster import read_in_conversion_file

    benchmark(read_in_conversion_file, GTPF_FILE_PATH)


def test_bench_read_pgdb_flat(benchmark):
    """Benchmark parsing a PGDB flat file (genes.dat)."""
    from plantclusterfinder_py.pgdb_utils import read_in_pgdb_flat

    benchmark(
        read_in_pgdb_flat,
        PGDB_DIR / "genes.dat",
        ("UNIQUE-ID", "ACCESSION-1"),
        unique_keys=("UNIQUE-ID", "ACCESSION-1"),
    )


def test_bench_read_pgdb_flat_reactions(benchmark):
    """Benchmark parsing the reactions.dat PGDB flat file."""
    from plantclusterfinder_py.pgdb_utils import read_in_pgdb_flat

    benchmark(
        read_in_pgdb_flat,
        PGDB_DIR / "reactions.dat",
        ("UNIQUE-ID", "ENZYMATIC-REACTION", "EC-NUMBER"),
    )


# --- Clustering benchmarks ---


def test_bench_compute_background(benchmark):
    """Benchmark computing background enzyme distribution for clustering."""
    from plantclusterfinder_py.cluster import compute_background
    from plantclusterfinder_py.genelocation import GeneLoc

    gene_location = GeneLoc(
        file_path=Path("./tests/data/background_computation_test.tsv")
    )
    has_reaction = pl.read_csv(
        Path("./tests/data/background_computation_has_reaction.csv")
    )
    gene_location._data = gene_location.data.join(has_reaction, on="Gene Name")

    step_sizes = [1 + x for x in range(4) for _ in range(100)]
    enzyme_counts = [3 for _ in range(400)]
    window_sizes = [3 + x for _ in range(100) for x in range(4)]

    cluster = pl.DataFrame(
        (
            {
                "step_size": x,
                "enzyme_count": y,
                "non_enzyme_gene_count": z - y,
                "test_window_size": z,
            }
            for x, y, z in zip(step_sizes, enzyme_counts, window_sizes, strict=True)
        )
    )

    benchmark(compute_background, cluster, gene_location, pl.DataFrame())


def test_bench_evaluate_criterion(benchmark):
    """Benchmark the criterion evaluation logic for step size selection."""
    from plantclusterfinder_py.cluster import evaluate_criterion

    cluster_tuple = namedtuple(
        "cluster_tuple",
        [
            "chromosome",
            "step_size",
            "start",
            "end",
            "enzyme_count",
            "non_enzyme_gene_count",
            "hypothetical_gene_count",
        ],
    )

    clusters = pl.DataFrame(
        [
            cluster_tuple("TEST", 3, x * 50, (x + 1) * 50 + 12, 6, 12, 0)
            for x in range(10)
        ]
    )

    benchmark(evaluate_criterion, 1, clusters, 3, 24, 3)


def test_bench_get_metabolic_domains(benchmark):
    """Benchmark loading and filtering metabolic domain data."""
    from plantclusterfinder_py.cluster import get_metabolic_domains

    benchmark(
        get_metabolic_domains,
        Path("./Inputs/ReactionMetabolicDomainClassification.txt"),
        ["Specialized Metabolism"],
    )
