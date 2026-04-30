"""Parse Input genome and GFF3 annotation files.

Parse input files and produce gtpf, glof, and proteins FASTA required for PlantClusterFinder.
"""

import logging
import re
from collections import deque
from pathlib import Path
from typing import Any, cast

import gffutils
import polars as pl
from BCBio.GFF import parse as parse_gff
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqRecord import SeqRecord
from pyfaidx import Fasta

_log = logging.getLogger(__name__)

glof_header = (
    "Gene Name",
    "Gene Start (bp)",
    "Gene End (bp)",
    "Chromosome Name",
    "Strand",
)

gtpf_header = ("Gene Name", "Transcript Name", "Protein Name")


def extract_genes_from_gff_and_fasta(
    fasta_file: Path,
    gff3_file: Path,
    gene_name_column: str,
    protein_name_column: str | None = None,
    protein_name_regex: str = ".*",
) -> tuple[Path, Path, Path]:
    r"""Extract genes form GFF and FASTA files.

    Produce the gtpf, glof, and Amino Acid FASTA files from the organism's
     genome FASTA file and GFF3 file.

    Args:
        fasta_file: Path to the genome FASTA file.
        gff3_file: Path the gene GFF3 file.
        gene_name_column: Name column to use in the GFF3 file for genes (e.g. "Name").
        protein_name_column: Name column to use in the GFF3 file for
                proteins/transcripts (mRNA) (e.g. "original_protein_id")
                (Optional: will use gene_name_column if not specified).
        protein_name_regex: a regular expression to isolate the protein name from the column,
                e.g. "(?!\\|)[^|]$" to select everything after the last | character
                (Optional: defaults to ".*").

    Returns:
        Tuple containing the paths of the GTPF, GLOF, and AA FASTA file.

    """
    if not fasta_file.exists():
        raise Exception("Cannot find Fasta File " + str(fasta_file))
    if not gff3_file.exists():
        raise Exception("Cannot find required GFF3 file " + str(gff3_file))
    if not protein_name_column:
        protein_name_column = gene_name_column
    gtpf_filename = Path(fasta_file.parent / f"gtpf_{fasta_file.stem}.tsv")
    aafasta_filename = Path(fasta_file.parent / f"{fasta_file.stem}.protein.faa")
    glof_filename = Path(gff3_file.parent / f"glof_{gff3_file.stem}.tsv")
    retval = (gtpf_filename, glof_filename, aafasta_filename)
    for f in retval:
        try:
            _ = f.open(mode="a")
            f.unlink()
        except Exception as e:
            _log.error("Cannot write to file %s: %s", f, e)
            raise e
    # Attempt to extract GFF3 using BCBio first, but failover to gffutils
    bcbio_ok = False
    gtpf = glof = aafasta = None
    try:
        _log.info("Attempting to extract GFF3 using BCBio.GFF")
        (gtpf, glof, aafasta) = extract_genes_bcbio(
            gff3_file,
            fasta_file,
            gene_name_column,
            protein_name_column,
            protein_name_regex,
        )
        bcbio_ok = True
    except Exception as e:
        _log.warning("BCBio failed to parse GFF3 file; threw the error %s", e)
        _log.info("Attempting to extract GFF3 using gffutils")
    # GFFUtils is more resiliant to structural issues with GFF3 files
    # but it is slower, as it needs to create a SQLite DB before it can run
    if not bcbio_ok:
        try:
            (gtpf, glof, aafasta) = extract_genes_gffutils(
                gff3_file,
                fasta_file,
                gene_name_column,
                protein_name_column,
                protein_name_regex,
            )
        except Exception as e:
            _log.error("GFFUtils also failed to parse GFF3 file; threw the error %s", e)
            raise Exception("Could not parse GFF3 file.") from e
    if gtpf is None or glof is None or aafasta is None:
        raise ValueError("Missing Inputs: GTPF, GLOF, or aaFASTA not created.")
    if len(gtpf) == 0:
        _log.error("No gene-protein pairs found.")
        _log.info("GFF3 File: %s, Attribute of Interest: %s", gff3_file, gene_name_column)
        raise Exception("Could not extract genes from GFF3 file. Manual intevention may be required")
    _log.debug("Writing GTPF %s (%d rows)", gtpf_filename, len(gtpf))
    pl.DataFrame(gtpf, schema=gtpf_header, orient="row").write_csv(gtpf_filename, separator="\t")
    _log.debug("Writing GLOF %s (%d rows)", glof_filename, len(glof))
    pl.DataFrame(glof, schema=glof_header, orient="row").write_csv(glof_filename, separator="\t")
    _log.debug("Writing protein AA FASTA %s (%d records)", aafasta_filename, len(aafasta))
    with aafasta_filename.open(mode="w") as fasta_fh:
        FastaWriter(fasta_fh).write_file(aafasta)
    return retval


def extract_genes_gffutils(
    gff3_file: Path,
    fasta_file: Path,
    gene_name_column: str,
    protein_name_column: str,
    protein_name_regex: str,
) -> tuple[
    deque[tuple[str, str, str]],
    deque[tuple[str, int | None, int | None, str, str]],
    deque[SeqRecord],
]:
    r"""Extract genes from GFF3 file using gffutils.

    Args:
        gff3_file: Path the gene GFF3 file.
        fasta_file: Path to the genome FASTA file.
        gene_name_column: Name column to use in the GFF3 file for genes (e.g. "Name").
        protein_name_column: Name column to use in the GFF3 file for
                proteins/transcripts (mRNA) (e.g. "original_protein_id")
                (Optional: will use gene_name_column if not specified).
        protein_name_regex: a regular expression to isolate the protein name from the column,
                e.g. "(?!\\|)[^|]$" to select everything after the last | character
                (Optional: defaults to ".*").

    Returns:
        Tuple containing data for GTPF, GLOF, and AAFasta.

    """
    glof: deque[tuple[str, int | None, int | None, str, str]] = deque()
    aafasta: deque[SeqRecord] = deque()
    gtpf: deque[tuple[str, str, str]] = deque()
    db_filename = Path(fasta_file.parent / "gff3.db")
    fasta_faidx_obj = Fasta(fasta_file, as_raw=False)
    if not db_filename.exists():
        gff3_data = gffutils.create_db(
            str(gff3_file), str(db_filename), merge_strategy="create_unique"
        )
        _log.debug("Created new gff_db file at %s.", db_filename)
    else:
        gff3_data = gffutils.FeatureDB(db_filename)
        _log.debug("Opened gff_db file at %s.", db_filename)
    no_name = 0
    no_children = 0
    for gene in gff3_data.features_of_type("gene"):
        if not gene.chrom:
            raise ValueError("No chromosome for a gene in GFF3. Invalid data.")
        if gene_name_column not in gene.attributes:
            _log.debug("No %s for gene ID %s, cannot save", gene_name_column, gene.id)
            no_name += 1
            continue
        gene_name: str = cast(str, gene[gene_name_column][0])
        _log.debug("Found gene %s", gene_name)
        longest_child = None
        for child in gff3_data.children(gene.id, featuretype="mRNA"):
            if protein_name_column not in child.attributes:
                _log.debug(
                    f"No {protein_name_column} for gene ID {gene.id}, mRNA {child.id}, cannot save"
                )
                continue
            if longest_child is None or len(child) > len(longest_child):
                longest_child = child
        cds_data: deque[gffutils.Feature] = deque()
        for child in gff3_data.children(gene.id, featuretype="CDS"):
            cds_data.append(child)
        if longest_child is None:
            _log.debug("no child for %s", gene_name)
            no_children += 1
            continue
        longest_child_name = longest_child[protein_name_column][0]
        if protein_name_regex != ".*":
            re_search = re.search(protein_name_regex, longest_child_name)
            if re_search is None:
                _log.error(
                    "Regex %s does not match %s. Cannot extract protein name.",
                    protein_name_regex,
                    longest_child_name,
                )
                raise Exception("Protein Regex not matched in GFF3 file.")
            if len(re_search.groups()) > 0:
                longest_child_name = "".join(re_search.groups())
            else:
                longest_child_name = re_search.group()

        shifted_start = None if not gene.start else gene.start + 1
        gene_end = None if not gene.end else cast(int, gene.end)
        gene_chrom = cast(str, gene.chrom)
        glof.append(
            (
                gene_name,
                shifted_start,
                gene_end,
                gene_chrom,
                gene.strand,
            )
        )
        try:
            gtpf.append(
                (
                    gene_name,
                    longest_child_name,
                    longest_child_name,
                )
            )
        except Exception as e:
            _log.error("failed to append details to gtpf")
            raise e
        aafasta.append(
            SeqRecord(
                Seq("".join(cds.sequence(fasta_faidx_obj) for cds in cds_data)).translate(),
                id=longest_child_name,
            )
        )
    _log.debug(
        "Completed processing %d genes, couldn't process %d genes for not having a name, "
        + "and %d genes for not having CDS/mRNA children",
        len(gtpf),
        no_name,
        no_children,
    )
    # guarantee gff3_data is garbage collected
    del gff3_data
    return (gtpf, glof, aafasta)


def extract_genes_bcbio(
    gff3_file: Path,
    fasta_file: Path,
    gene_name_column: str,
    protein_name_column: str,
    protein_name_regex: str,
) -> tuple[
    deque[tuple[str, str, str]],
    deque[tuple[str, int | None, int | None, str, str]],
    deque[SeqRecord],
]:
    r"""Extract genes from GFF3 file using gffutils.

    Args:
        gff3_file: Path the gene GFF3 file.
        fasta_file: Path to the genome FASTA file.
        gene_name_column: Name column to use in the GFF3 file for genes (e.g. "Name").
        protein_name_column: Name column to use in the GFF3 file for
                proteins/transcripts (mRNA) (e.g. "original_protein_id")
                (Optional: will use gene_name_column if not specified).
        protein_name_regex: a regular expression to isolate the protein name from the column,
                e.g. "(?!\\|)[^|]$" to select everything after the last | character
                (Optional: defaults to ".*").

    Returns:
        Tuple containing data for GTPF, GLOF, and AAFasta.

    """
    glof: deque[tuple[str, int | None, int | None, str, str]] = deque()
    aafasta: deque[SeqRecord] = deque()
    gtpf: deque[tuple[str, str, str]] = deque()
    seq_fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    no_name = 0
    no_children = 0
    with gff3_file.open() as gff_fh:
        # go over every macro feature (chromosome or scaffold)
        for macro_feature in parse_gff(gff_fh, base_dict=seq_fasta_dict):
            # go over every gene or gene-like feature
            for gene in macro_feature.features:
                if gene_name_column not in gene.qualifiers:
                    _log.warning(
                        "Gene %s in Chromosome/Scaffold %s is missing the %s qualifier.",
                        gene.id,
                        macro_feature.id,
                        gene_name_column,
                    )
                    no_name += 1
                    continue
                gene_name = gene.qualifiers.get(gene_name_column, [None]).pop(0)
                if gene_name is None:
                    _log.error("Gene name column %s not found for file %s.", gene_name, gff3_file)
                    raise Exception("Can not find gene name column using bcbio.")
                _log.debug("Found gene %s", gene_name)
                glof.append(
                    (
                        gene_name,
                        int(gene.location.start),
                        int(gene.location.end),
                        macro_feature.id,
                        gene.location.strand,
                    )
                )
                largest_isoform = None
                largest_isoform_len = 0
                for sub_feature in gene.sub_features:
                    if protein_name_column not in sub_feature.qualifiers:
                        continue
                    if len(sub_feature) > largest_isoform_len:
                        largest_isoform = sub_feature
                        largest_isoform_len = len(sub_feature)
                if largest_isoform is None:
                    _log.warning("Could not find transcript for gene %s", gene.id)
                    no_children += 1
                    continue
                longest_child_name = largest_isoform.qualifiers[protein_name_column][0]
                if protein_name_regex != ".*":
                    re_search = re.search(protein_name_regex, longest_child_name)
                    if re_search is None:
                        _log.error(
                            "Regex %s does not match %s. Cannot extract protein name.",
                            protein_name_regex,
                            longest_child_name,
                        )
                        raise Exception("Protein Regex not matched in GFF3 file.")
                    if len(re_search.groups()) > 0:
                        longest_child_name = "".join(re_search.groups())
                    else:
                        longest_child_name = re_search.group()
                aafasta.append(
                    SeqRecord(gene.translate(macro_feature.seq, cds=False), id=longest_child_name)
                )
                gtpf.append(
                    (
                        gene_name,
                        longest_child_name,
                        longest_child_name,
                    )
                )
    # guarantee seq_fasta_dict is garbage collected
    del seq_fasta_dict
    _log.warning("%s genes were ignored for not having attribute %s", no_name, gene_name_column)
    return (gtpf, glof, aafasta)
