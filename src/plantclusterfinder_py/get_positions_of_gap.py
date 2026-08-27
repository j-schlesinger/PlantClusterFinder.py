"""Get positions of gaps in FASTA files.

Functions to extract, analyze, and return the gaps seen in FASTA files.
"""

import logging
from argparse import ArgumentParser
from collections import deque
from collections.abc import Generator, Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor
from itertools import groupby
from operator import itemgetter
from pathlib import Path
from typing import cast

import numpy as np

log = logging.getLogger()


def char_run_search(chars: str, input_str: str) -> Iterator[tuple[int, int]]:
    """Search for the positions of characters in a string from an input list of characters.

    Args:
        chars: string containing characters to search for
        input_str: string being searched

    Returns:
        list: indices for the characters input

    """
    str_array = np.frombuffer(b"\x00" + input_str.encode("ascii") + b"\x00", dtype=np.uint8)
    test_list = [ord(x) for x in chars]
    pattern_true = np.isin(str_array, test_list)
    index_list = np.where(np.diff(pattern_true))[0]
    run_starts: list[int] = cast(list[int], ((index_list + 1)[::2]).tolist())
    run_lens: list[int] = cast(list[int], (np.diff(index_list)[::2]).tolist())
    if len(run_lens) == len(run_starts) - 1:
        run_lens.append(1)
    return zip(run_starts, run_lens, strict=True)


# From BioPython
def read_fasta(fp: Iterable[str]) -> Generator[tuple[str, str], None, None]:
    """From BioPython: generator which parses FASTA files into tuples of names and sequences.

    Args:
        fp (file handle): Readable handle for FASTA file.

    Yields:
        (name, sequence): the next name and sequence in the FASTA file.

    """
    name = ""
    seq: deque[str] = deque()
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if len(name) > 0:
                yield name, "".join(seq)
            name = line
            seq.clear()
        else:
            seq.append(line)
    if name:
        yield name, "".join(seq)


def groupbyseq(input_list: list[int], amino_acid: str, chromosome: str) -> list[str]:
    """Write the sequence groups to a file given the inputs.

    NOTE: This function is not used by PlantClusterFinder.py

    Args:
        input_list: List of indices of amino acids from FASTA
        amino_acid: single character representation of amino acid in group
        chromosome: name of chromosome from FASTA file

    Returns:
        Sequence group rows - start index and length of groups

    Raises:
        TypeError: When the data from the indices list undergoes a failed
            type cast. This may happen when trying to coerce the first group
            index to an int, or the length to a string.
        IndexError: When the data from the indices list is corrupted or
            incomplete.

    """
    log.debug("Performing groupbyseq")
    output_str: deque[str] = deque()

    for _, g in groupby(enumerate(input_list), lambda x: x[1] - x[0]):
        cur_group = None
        try:
            cur_group = list(map(itemgetter(1), g))
            start_index = str(cur_group[0] + 1)
            length = str(len(cur_group))
            output_str.extend([chromosome, " ", amino_acid, " ", start_index, " ", length, "\n"])
        except TypeError as e:
            err_group = "ERR GETTING GROUP" if cur_group is None else ",".join(map(str, cur_group))
            log.error(
                "Error %s %s",
                e,
                err_group,
            )
        except IndexError as e:
            err_group = "ERR GETTING GROUP" if cur_group is None else ",".join(map(str, cur_group))
            log.error(
                "Error %s %s",
                e,
                err_group,
            )
    return list(output_str)


def _process_sequence(seq_name_and_seq_aa: tuple[str, str]) -> deque[str]:
    seq_name, seq_aa = seq_name_and_seq_aa
    seq_aa = seq_aa.upper()
    return_rows: deque[str] = deque()
    n_indices = char_run_search("N", seq_aa)
    char = seq_name[1:].split()[0]
    for start, length in n_indices:
        return_rows.extend([char, " N ", str(start), " ", str(length), "\n"])
    a_indices = char_run_search("ATGC*", seq_aa)
    for start, length in a_indices:
        return_rows.extend([char, " A ", str(start), " ", str(length), "\n"])
    return return_rows


def get_positions_of_gap(input_path: Path, output_path: Path, parallelize: int = 1) -> Path:
    """Get the positions of gaps from the input FASTA and save them in an output file.

    Args:
        input_path: Path of input amino acid FASTA file.
        output_path: Path of writable output file. This will be created
            if it does not exist, otherwise it will be appended to.
        parallelize: Parallelize processing of FASTA file (not implemented)

    """
    log.debug("Performing get_positions_of_gap")
    parallelize = max(parallelize, 1)

    def _write_list(fh, input_list: deque[str]) -> None:
        for x in input_list:
            fh.write(x)

    log.debug("Reading AA FASTA from %s", input_path)
    with (
        input_path.open() as fp,
        output_path.open("w") as op,
        ProcessPoolExecutor(max_workers=parallelize) as executor,
    ):
        output_list = list(executor.map(_process_sequence, read_fasta(fp)))
    with output_path.open("w") as op:
        for x in output_list:
            _write_list(op, x)
    return output_path


def main(args=None):
    """Script get_positions_of_gap used to be a standalone AWK application.

    This main function exists for testing and backward compatibility.
    """
    parser = ArgumentParser("PCF_Groupby Seq")
    parser.add_argument(
        "-i",
        "--input_path",
        dest="input_path",
        required=True,
        help="The path of the input FASTA",
    )
    parser.add_argument(
        "-o",
        "--output_path",
        dest="output_path",
        required=True,
        help="The path of the output groups",
    )
    parser.add_argument(
        "-p",
        "--parallel",
        dest="para",
        default=1,
        type=int,
        help="The number of workers to parallelize with",
    )
    args = parser.parse_args(args)
    input_path = Path(args.input_path)
    output_path = Path(args.output_path)

    if output_path.is_file():
        output_path.unlink()
    get_positions_of_gap(input_path, output_path, parallelize=args.para)


if __name__ == "__main__":
    main()
