import re
from pathlib import Path


def test_regex_search():
    from plantclusterfinder_py.get_positions_of_gap import regex_search

    assert [1, 3, 5] == regex_search("A", "BANANA")
    assert 23 == len(
        regex_search(
            "aoeiu",
            "Now is the time for all good men to come to the aid of their community.",
        )
    )
    assert 17 in regex_search(
        " \n", "Haikus can be fun\nbut sometimes they don't make sense\nrefigerator"
    )
    assert [] == regex_search(
        "259", "You have either chosen, or been chosen to relocate to City 17."
    )


def test_read_fasta():
    from plantclusterfinder_py.get_positions_of_gap import read_fasta

    with open("test/data/test_fasta.fa") as fa_fh:
        fasta_gen = read_fasta(fa_fh)
        seq = {n: s for n, s in fasta_gen}
    assert ">Test Sequence 1" in seq.keys()
    assert ">Test Sequence 2" in seq.keys()
    assert 143 == len(seq[">Test Sequence 1"])
    assert 80 * 16 == len(seq[">Test Sequence 2"])


def test_groupbyseq():
    from plantclusterfinder_py.get_positions_of_gap import (
        groupbyseq,
        read_fasta,
        regex_search,
    )

    output_path = Path("/tmp/test_groupbyseq.csv")  # output
    with open(Path("test/data/test_fasta.fa"), "r") as fa_fh:
        for seq_name, seq_aa in read_fasta(fa_fh):
            fasta_id = re.split(r"[\t\s|]", seq_name)[0].replace(">", "")
            indices = regex_search("N", seq_aa.upper())
            groupbyseq(indices, output_path, "X", fasta_id)
            indices = regex_search(["A", "T", "C", "G", "*"], seq_aa.upper())
            groupbyseq(indices, output_path, "A", fasta_id)
    test_d = dict()
    with open(output_path, "r") as test_output:
        for line in test_output:
            nuc = line[5]
            test_d[nuc] = test_d.get(nuc, 0) + 1
    output_path.unlink()
    assert 18 == test_d["X"]
    assert 230 == test_d["A"]


def test_get_positions_of_gap():
    from plantclusterfinder_py.get_positions_of_gap import get_positions_of_gap

    output_path = Path("/tmp/test_groupbyseq.csv")  # output
    input_path = Path("test/data/test_fasta.fa")
    get_positions_of_gap(input_path, output_path)
    test_d = dict()
    with open(output_path, "r") as test_output:
        for line in test_output:
            nuc = line[5]
            test_d[nuc] = test_d.get(nuc, 0) + 1
    output_path.unlink()
    assert 18 == test_d["N"]
    assert 230 == test_d["A"]
