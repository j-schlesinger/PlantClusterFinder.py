import re
import os
import logging
from argparse import ArgumentParser
from operator import itemgetter
from itertools import groupby
from collections import deque

log = logging.getLogger()


def regex_search(string_pattern_list, string_input):
    """Search for the positions of characters
    in a string from an input list of characters.

    Args:
        string_pattern_list: list of characters to search for
        string_input: string being serched

    Returns:
        list: indices for the characters input
    """
    string_pattern = re.compile("[" + ''.join(string_pattern_list) + "]")
    return [m.start() for m in re.finditer(string_pattern, string_input)]


# From BioPython
def read_fasta(fp):
    """
    From BioPython: generator which parses FASTA files
    into tuples of names and sequences

    Args:
        fp (filehandle): Readable handle for FASTA file

    Yields:
        (name, sequence): the next name and sequence in the FASTA file
    """
    name, seq = None, deque()
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield name, ''.join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield name, ''.join(seq)


def groupbyseq(input_list: list, output_path, amino_acid: str, chromosome: str):
    """
    This function writes the sequence groups to a file given the inputs

    Args:
        input_list: List if indices of of amino acids from FASTA
        output_path: Writable output file
        amino_acid: single character representation of amino acid
        chromosome: name of chromosome from FASTA file

    Returns:
        None

    Raises:
        TypeError: When the data from the indices list undergoes a failed
            type cast. This may happen when trying to coerce the first group
            index to an int, or the length to a string.
        IndexError: When the data from the indices list is corrupted or
            incomplete.
        IOError: When the output_path cannot be written to.
    """
    if os.path.exists(output_path):
        append_write = 'a'  # append if already exists
    else:
        append_write = 'w'  # make a new file if not

    with open(output_path, append_write) as op:
        for k, g in groupby(enumerate(input_list), lambda x: x[1] - x[0]):
            try:
                cur_group = list(map(itemgetter(1), g))
                start_index = str(cur_group[0] + 1)
                length = str(len(cur_group))
                op.write(
                    ' '.join([chromosome, amino_acid, start_index, length]) + "\n")
            except (TypeError or IndexError) as e:
                log.error("Error", e, list(map(itemgetter(1), g)))


def get_positions_of_gap(input_path: str, output_path: str):
    """
    Get the positions of gaps from the input FASTA and save them in an 
        output file.

    Args:
        input_path: Path of input amino acid FASTA file.
        output_path: Path of writable output file. This will be created
            if it does not exist, otherwise it will be appended to.
    """
    with open(input_path, 'r') as fp:
        for seq_name, seq_aa in read_fasta(fp):
            fasta_id = re.split(r'[\t\s|]', seq_name)[0].replace('>', '')
            indices = regex_search("N", seq_aa.upper())
            groupbyseq(indices, output_path, "N", fasta_id)
            indices = regex_search(["A", "T", "C", "G", "*"], seq_aa.upper())
            groupbyseq(indices, output_path, "A", fasta_id)


def main():
    parser = ArgumentParser("PCF_Groupby Seq")
    parser.add_argument('-i', '--input_path', dest="input_path",
                        required=True, help="The path of the input FASTA")
    parser.add_argument('-o', '--output_path',
                        dest="output_path", required=True, help="The path of the output groups")
    args = parser.parse_args()
    input_path = args.input_path
    output_path = args.output_path

    try:
        os.remove(output_path)
    except OSError:
        pass
    get_positions_of_gap(input_path, output_path)


if __name__ == '__main__':
    main()
