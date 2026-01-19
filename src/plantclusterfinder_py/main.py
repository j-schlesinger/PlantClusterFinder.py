import argparse
import logging
import logging.handlers
import math
import os
import queue
import re
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from pathlib import Path
from typing import Any

from scipy.io import loadmat, savemat  # type: ignore

from plantclusterfinder_py.cluster import (
    annotate_sequencing_gaps,
    get_MCL_clustering,
    get_metabolic_domains,
    introduce_hypo_genes,
    label_reactions_with_sig_and_tail,
    map_pgdb_reactions_to_genes,
    map_PGDBs_to_conversion_file,
    match_PGDBs_to_MDs,
    perform_clustering,
    read_in_conversion_file,
    read_MCL_clustering_file,
    write_cluster_output,
    write_gene_output,
)
from plantclusterfinder_py.genelocation import GeneLoc

LOG_FORMAT = "{asctime} | {levelname:<8} | {name} : {message}"


def save(filename: str, *args, caller_locals: dict = {}):
    """
    ---
    Save the state of PlantClusterFinder.py (data only)
    ---
    Args:
        filename: name of file where object is being stored
        caller_locals: result of locals() from calling function
    """
    parser = ArgumentParser(prog="SaveFunction", exit_on_error=False)
    parser.add_argument("-regexp", type=str)
    parser.add_argument("-v7.3", action="store_true", dest="version")
    parsed_args = parser.parse_args(args)
    if isinstance(filename, list):
        filename = "".join(filename)
    if parsed_args.regexp is not None:
        savemat(
            filename,
            {
                x: y
                for x, y in caller_locals.items()
                if re.search(parsed_args.regexp, x)
            },
            long_field_names=True,
        )
    else:
        raise Exception("Does not support non-regexp save functions")


def load(filename: str, *variables) -> dict[str, Any]:
    """
    ---
    Load list of variables from a PlantClusterFinder MATLAB file
    ---
    Args:
        filename: name of file where object is stored
        *variables: variable names to load
    Returns:
        Dict of saved variables loaded from file.
    """
    if isinstance(filename, list):
        filename = "".join(filename)
    if len(variables) == 0:
        return loadmat(filename)
    else:
        return {x: y for x, y in loadmat(filename).items if x in variables}


def _tempsaves(caller_locals: dict) -> dict[str, Any]:
    """
    Perfrom temporary save for PlantClusterFinder between steps
    Args:
        caller_locals: locals() of calling function
    Returns:
        If the existing temp save is further along than the caller,
        the vPartFinished from that, otherwise the caller's vPartFinished
    """
    ret = caller_locals["vPartFinished"]
    temp_filename = caller_locals.get("vGeneLocation_File", "TEMP") + "_tempsave.mat"
    if caller_locals.get("vTempsavesOverwrite", 0) == 1:
        save(
            temp_filename,
            "-regexp",
            "^(?!vRunAfterPart$).",
            "-v7.3",
            caller_locals=caller_locals,
        )
    elif not os.path.isfile(temp_filename):
        save(
            temp_filename,
            "-regexp",
            "^(?!vRunAfterPart$).",
            "-v7.3",
            caller_locals=caller_locals,
        )
    elif os.path.isfile(temp_filename):
        vPartFinished_old_1 = ret
        ret = load(temp_filename, "vPartFinished")["vPartFinished"]
        if vPartFinished_old_1 > ret:
            ret = vPartFinished_old_1
            del vPartFinished_old_1
            save(
                temp_filename,
                "-regexp",
                "^(?!vRunAfterPart$).",
                "-v7.3",
                caller_locals=caller_locals,
            )
    return ret


def main() -> None:
    """
    PlantClusterFinder 1.3
    PlantClusterFinder (PCF) detects metabolic gene clusters in a sequenced
    genome. It uses a gene location file provided by the user (see below) and
    a PGDB created with Pathway Tools as well as further information (see
    below) to identify enzyme-coding genes (metabolic genes) located together
    on a chromosome. Initially only continous stretches of metabolic genes
    lying directly next to each other are allowed. This condition is relaxed
    by iteratively increasing the intervening (non-metabolic) gene size by
    one. Several criteria to select for clusters are provided. In addition to
    this, clusters can be prevented from forming by a section of criteria.
    Details of PCF (version 1.0) can be found in PMID: 28228535.

    The major differences between this version (1.3) and previous versions
    (1.0 and 1.2) are:

      1) Physical breaks of the genome or sequencing gaps of unknown size are
         typically encoded by stretches of Ns in the genome assembly fasta
         file. Previously we inserted N hypothetical genes. This however
         diluted the background of low quality genomes with non-enzymes, and
         hence the likelyhood of a cluster to be classified as top x%% of
         enzyme dense regions was better than in a genome that had good
         quality. In version 1.3 we identify these breaks and prevent
         formation of a cluster over these gaps.
      2) Any sequencing information that is missing is typically hard masked
         with Ns. Previously, any intergenic region affected by at least one
         N was evaluated for its length, and hypothetical genes were inserted
         accordingly (See Schlapfer et al, PMID: 28228535). This led
         sometimes to unrealistic prevention of detecting gene clusters. It
         is unlikely that missing information about a single nucleotide would
         (if it would be known) lead to the finding of multiple gene models.
         Thus here we changed the code to insert 2 hypothetical genes only if
         a strech of unknown sequence is larger than nth percentile of gene
         sequences (set to 5). We also provide the option to NOT insert any
         hypothetical genes all together. Instead we by default we use
         MaxSeqGapSize set to 100000 and MaxInterGeneDistByMedian set to 50
         resulting in similar clusters as in PCF version 1.0.
      3) Large gene poor intergenic regions are present in genomes. In this
         version we provide the option of several parameters to prevent
         clusters from spanning such large gene poor regions.

    Authors: Pascal Schlapfer, December 2017
             Bo Xue, December 2017

    """
    # print(main.__doc__)
    help_epilog = """ Outputs:
         None (Two files are generated: vGeneOutputFile, containing all the
              information about genes and vClusterOutputFile, containing all
              the information about clusters)"""
    # Note, PGDBIdsToMap is only needed here because the pgdb
    # contains protein Ids in place of gene identifiers.

    parser = ArgumentParser(
        prog="PlantClusterFinder",
        description="Find Metabolic Gene Clusters in Plant Genomes",
        epilog=help_epilog,
        formatter_class=RawTextHelpFormatter,
    )
    # Mandatory inputs: all files need full path, order does not matter.
    parser.add_argument(
        "-pgdb",
        required=True,
        type=str,
        help="Give full path of PGDB flat file folder, where the flat files of the species pgdb is stored.",
    )
    parser.add_argument(
        "-rmdf",
        required=True,
        type=str,
        help="""
        TSV file of the metabolic domains of reactions. Check
         if file is outdated (needs to contain all reactions
         that are present in the PGDB that have genes
         annotated to.""",
    )
    parser.add_argument(
        "-md",
        required=True,
        action="extend",
        nargs="+",
        type=str,
        help="""
         A list of Metabolic domains that should be analyzed.
          For example {'Amines and Polyamines Metabolism'; 'Amino Acids '}""",
    )
    parser.add_argument(
        "-psf",
        required=True,
        type=str,
        help="""
        Protein sequence file of the genome.
         Can only contain Protein IDs as header e.g. >protein1-ID""",
    )
    parser.add_argument(
        "-gtpf",
        required=True,
        type=str,
        help="""
        TSV file with a header describing gene, transcript,
         and protein identifier and then for all the
         identfiers the listings (gene-ID, transcript-ID,
         protein-ID). This is used to map information
         regarding transcripts (currently not used) and
         proteins (used for MCL clustering, potentially also
         used in PGDBs) to the genes (e.g. the gene location
         information).""",
    )
    parser.add_argument(
        "-glof",
        required=True,
        type=str,
        help="""
        Gene Location (position) TSV file with gene ID, start bp, end bp,
         chromosome / scaffold, strand (encoded as 1 or -1).
         This file can be generated with a biomart or out of a gff3 file.""",
    )
    parser.add_argument(
        "-dnaf",
        required=True,
        type=str,
        help="""
        Fasta file of the HARDMASKED genome nucleotide sequences.""",
    )
    parser.add_argument(
        "-sitf",
        required=True,
        type=str,
        help="""
        Signature and tailoring reactions TSV file describing the classes
         of enzymes that should be classified as signature or tailoring
         enzymes to identify clusters containing such.""",
    )
    parser.add_argument(
        "-gout",
        required=True,
        type=str,
        help="""
        Path and file name to the gene-Outputfile,
         make sure you have access to the file and folder.
         This is the first results output file that will be generated.""",
    )
    parser.add_argument(
        "-cout",
        required=True,
        type=str,
        help="""
        Path and file name to the cluster-Outputfile,
         make sure you have access to the file and folder.
        This is the second Results output file that will be generated.""",
    )
    # Optional Inputs:
    parser.add_argument(
        "--GenesDat_GeneInfo",
        nargs="+",
        default=["UNIQUE-ID", "ACCESSION-1", "ACCESSION-2"],
        type=str,
        help="""
        Needs to be followed by a list of Attributes in
         genes.dat. Defines where gene ID information
         should be searched for.
         Default is "UNIQUE-ID, ACCESSION-1, ACCESSION-2" """,
    )
    parser.add_argument(
        "--MCLClusterFile",
        type=int,
        choices=[0, 1],
        default=0,
        help="""
        Needs to be followed by either 1 or 0. With a default
         of 0. 1 indicates that the input file is not a
         protein fasta file, but a precomputed MCL
         clustering file (precomputing can be usefull for
         speed). If such a file is used, then nothing but
         gene IDs and tabs and new lines are allowed in this
         file.""",
    )
    parser.add_argument(
        "--HypoGenePercentile",
        type=int,
        default=5,
        help="""
        Needs to be followed by a number. Defines the
         length of sequencing gap that should be
         populated by two hypothetical genes. If this is
         10, then two hypothetical genes are introduced
         as soon as the sequencing gap exceeds the
         (100-10) lower percentile of the background
         size of all genes of the genome (default is 5).""",
    )
    parser.add_argument(
        "--MaxSeqGapSize",
        type=int,
        default=100000,
        help="""
        Needs to be followed by a number. Maximal masked
         nucleotide region (in bp) that a cluster is allowed
         to bridge (default is 100000). """,
    )
    parser.add_argument(
        "--MaxInterGeneDist",
        type=int,
        default=None,
        help="""
        Needs to be followed by a number. Maximal
         intergenic distance ( in bp) that still can be
         crossed by a gene cluster. (default is None, thus
         inactive).""",
    )
    parser.add_argument(
        "--MaxInterGeneDistByMedian",
        type=int,
        default=50,
        help="""
        Needs to be followed by a number. Maximal
         intergenic distance that still can be
         crossed by a gene cluster, here defined
         by this number times the median of the
         gene sizes. (default is 50. Set it to -1
         if it should be inactive).""",
    )
    parser.add_argument(
        "--PercentileForMaxInterGeneDist",
        type=int,
        default=None,
        help="""
        Needs to be followed by a number.
         Calculates the maximal intergenic
         distance that still can be crossed
         by a gene cluster. This number
         determines the percentile to choose
         from all the distances. If this is
         99.9, then a cluster is allowed to
         bridge if the largest intergenic
         distance in the cluster is below
         99.9%% of the background size of all
         intergenic distances of the genome
         (default is -1,thus inactive).""",
    )
    parser.add_argument(
        "--SeqGapSizesChromBreak",
        type=int,
        nargs="+",
        default=[],
        help="""
        Needs to be followed by a list of numbers.
         The masked nucleotides sequences of these
         specific lengths will prevent cluster
         formation (Default is an empty list: []).""",
    )
    parser.add_argument(
        "--PreScreenGaps",
        type=int,
        choices=[0, 1],
        default=0,
        help="""
        Needs to be followed by either 1 or zero. Defines if
         sequencing gaps should be prescreened to delete them
         if they are anyway too small to matter. (default is 0)""",
    )
    parser.add_argument(
        "--OverwriteSeqGapInfo",
        type=int,
        choices=[0, 1],
        default=0,
        help="""
        Needs to be followed by either 1 or 0. Forces
         masked nucleotide analysis for finding
         sequencing gaps (default is 0, no enforcement).""",
    )
    parser.add_argument(
        "--OverwriteMCLClustering",
        type=int,
        choices=[0, 1],
        default=0,
        help="""
        Needs to be followed by either 1 or 0.
         Forces (if set to 1) that the MCL
         clustering is redone even if the result
         files are already present (default is 0).""",
    )
    parser.add_argument(
        "--Verbose",
        type=int,
        default=0,
        help="""
        Needs to be followed by a number (0 is default, 1 gives
         more screen-output, 2 gives exhaustive screen-output).""",
    )
    parser.add_argument(
        "--EnzMinForClusters",
        type=int,
        default=3,
        help="""
        Needs to be followed by a number: Defines the
         minimum number of enyzmes that need to be
         present in a cluster (default is 3).""",
    )
    parser.add_argument(
        "--KbRangeLD",
        type=int,
        default=100000,
        help="""
        Needs to be followed by a number. Genes that are
         separated by more than this size (in bp) are not
         considered as local duplicates. (default is 100000)""",
    )
    parser.add_argument(
        "--GeneRangeLD",
        type=int,
        default=10,
        help="""
        Needs to be followed by a number. Number of
         intervening genes that are allowed to separate two
         local duplicated genes. Two genes that are separated
         by more than this many intervening genes are not
         called local duplicates. (default is 10)""",
    )
    parser.add_argument(
        "--TopPercentClusters",
        type=int,
        default=5,
        help="""
        Needs to be followed by a number. Clusters are
         labeled if they are within this top X%% of
         enzyme dense clusters compared to the
         background. (default is 5)""",
    )
    parser.add_argument(
        "--MinStepSize",
        type=int,
        default=5,
        help="""
        Needs to be followed by a number. Minimal intervening
         gene size that should be computed (default is 5).""",
    )
    parser.add_argument(
        "--MaxStepSize",
        type=int,
        default=20,
        help="""
        Needs to be followed by a number. Maximal intervening
         gene size that should be computed (default is 20).""",
    )
    parser.add_argument(
        "--Criterion",
        type=int,
        default=2,
        choices=[1, 2, 3, 4],
        help="""
        Needs to be followed by a number. Criterion after which
              stepsize (intervening gene size) is chosen (default is 2).
         1: as published in Schlapfer et al 2017 (PMID:28228535)''
                 sum(enzymes in clusters) - sum(non-enzymes in clusters) > 0.
         2: stepsize 3 (because most organisms take this one)
         3: sum(enzymes in clusters) - sum(non-enzymes in clusters, without hypo genes) -2*sum(hypogenes)
         4: sum(clusters in stepsize n) > sum(clusters in stepsize n+1).""",
    )
    parser.add_argument(
        "-para",
        type=int,
        default=0,
        help="""
        How many cpus can be used by the algorithm
         (used for MCL clustering and Sequencing gap search).
         The default is 0 (auto-scale).""",
    )
    parser.add_argument(
        "--UnmaskedDNA",
        type=int,
        default=0,
        help="""
        Needs to be followed by a number. Defines if plant
         Cluster finder should continue despite a non-masked
         genome sequence file was used (default is 0, not
         continue).""",
    )
    parser.add_argument(
        "--PGDBIdsToMap",
        type=str,
        default="G",
        help="""
        Needs to be followed by a string. Can contain the
         Letters G and/or T and/or P. In that case (the pgdb
         is mapped to Gene (G) and/or Transcript (T) and/or
         Protein (P) Ids of the Gene conversion file.
         The default value is 'G'.""",
    )
    parser.add_argument(
        "--RunAfterPart",
        type=int,
        default=0,
        help="""
        Needs to be followed by a number. Sets the parts of
         the algorithm that should be recomputed. (default is 0)""",
    )
    parser.add_argument(
        "--Tempsaves",
        type=int,
        default=0,
        choices=[0, 1],
        help="""Defines if Temporary saves should be made. (default is 0)""",
    )
    parser.add_argument(
        "--TempsavesOverwrite",
        type=int,
        default=0,
        choices=[0, 1],
        help="""Defines if Temporary saves should be overwritten. (Default is 0).""",
    )
    parser.add_argument(
        "--RemoveNonProtLocations",
        type=int,
        default=0,
        choices=[0, 1],
        help="""
        Defines if the protein-fasta file should be
         used to clean out the gene position file
         from sequences that do not show up on the
         protein fasta file, and thus do not have
         had a chance to be predicted to be an
         enzyme. (Default 0)""",
    )
    parser.add_argument(
        "--InsertHypos",
        action="store_true",
        help="""
        Defines if the genome should be populated by
         hypothetical genes in sequencing gaps that are not
         physical breaks. (Default False)""",
    )
    parser.add_argument(
        "--HypoAmount",
        type=int,
        default=0,
        help="""
        Needs to be followed by a number. Defines by how many
         hypothetical genes a sequencing gap that is not a
         physical break should be populated by artificial genes
         to punish such sequencing gaps from becoming clusters.
         (Default 0)""",
    )
    parser.add_argument(
        "--OutputType",
        type=str,
        default="simple",
        choices=["simple", "old", "verbose"],
        help="""
        Defines the format of the Outputfiles: 'old' is using
         the old format to report all clusters of all stepsizes,
         'verbose' defines that top percent clusters of all
         stepsizes are reported. 'simple' just reports the top
         percent clusters of the chosen stepsize (default
         simple).""",
    )
    parser.add_argument(
        "--OutputGeneLocationParquet",
        type=str,
        help="""
        The file location of the parquet file to store the GeneLocation data.
         The file is only written if the argument is set.
        """,
    )
    parser.add_argument(
        "--NoRebuildBlastDB", action="store_false", help=argparse.SUPPRESS
    )
    parser.add_argument("--NoRerunBlastP", action="store_false", help=argparse.SUPPRESS)
    args = parser.parse_args()

    verbose_level = args.Verbose
    target_handler = logging.StreamHandler(sys.stderr)
    target_handler.setFormatter(logging.Formatter(fmt=LOG_FORMAT, style="{"))
    log_queue: queue.Queue[logging.LogRecord] = queue.Queue(-1)
    queue_handler = logging.handlers.QueueHandler(log_queue)
    if verbose_level == 1:
        target_handler.setLevel(logging.INFO)
        logging.getLogger("").setLevel(logging.INFO)
        # queue_handler.setLevel(logging.INFO)
    elif verbose_level >= 2:
        target_handler.setLevel(logging.DEBUG)
        logging.getLogger("").setLevel(logging.DEBUG)
        # queue_handler.setLevel(logging.DEBUG)
    listener = logging.handlers.QueueListener(
        log_queue, target_handler, respect_handler_level=True
    )
    listener.start()
    # Add the queue_handler to the root logger
    logging.getLogger("").addHandler(queue_handler)
    log = logging.getLogger("main")
    if verbose_level == 1:
        log.warning("Log verbosity set to 1 (INFO)")
        log.info("Info Logging Functional")
    elif verbose_level == 2:
        log.warning("Log verbosity set to 2 (DEBUG)")
        log.debug("Debug Logging Functional")
    elif verbose_level > 2:
        log.warning("Maximum verbosity exceeded, defaulting to 2 (DEBUG)")
        log.debug("Debug Logging Functional")

    list_of_MD_to_annotate = args.md
    # Verify is directory
    vPGDB_FlatFileFolder = Path(os.path.dirname(args.pgdb))
    # Verify is file
    vMD_reactions_File = args.rmdf
    vProtein_sequence_FastaFile = args.psf
    vGeneTranscriptProtein_mapping_File = args.gtpf
    vGeneLocation_File = args.glof
    vDNA_FastaFile = Path(args.dnaf)
    vSignatureTailorFile = args.sitf
    # Verify Writable
    gene_output_path = Path(args.gout)
    vClusterOutputFile = Path(args.cout)

    vGenesDat_GeneInfo = args.GenesDat_GeneInfo
    vHypoGenePercentile = args.HypoGenePercentile  # default 5
    vMaxSeqGapSize = args.MaxSeqGapSize
    vMaxInterGeneDist = args.MaxInterGeneDist
    vPercentileForMaxInterGeneDist = args.PercentileForMaxInterGeneDist
    vMaxInterGeneDistByMedian = args.MaxInterGeneDistByMedian
    vSeqGapSizesChromBreak = args.SeqGapSizesChromBreak
    top_percent_clusters = args.TopPercentClusters
    vMinStepSize = args.MinStepSize
    vMaxStepSize = args.MaxStepSize
    vCriterion = args.Criterion
    # Change to Boolean
    vMCLClusterFile = args.MCLClusterFile
    vPreScreenGaps = args.PreScreenGaps
    vOverwriteSeqGapInfo = args.OverwriteSeqGapInfo
    vOverwriteMCLClustering = args.OverwriteMCLClustering
    vUnmaskedDNA = args.UnmaskedDNA
    vTempsaves = args.Tempsaves
    vTempsavesOverwrite = args.TempsavesOverwrite
    vRemoveNonProtLocations = args.RemoveNonProtLocations
    vInsertHypos = args.InsertHypos

    vKbRangeLD = args.KbRangeLD
    vGeneRangeLD = args.GeneRangeLD
    vEnzMinForClusters = args.EnzMinForClusters
    # TODO: validate if computer has `para` CPUs. -1 = auto?
    vPara = args.para  # parallelize
    vPGDBIdsToMap = args.PGDBIdsToMap
    vRunAfterPart = args.RunAfterPart
    vHypoAmount = args.HypoAmount
    vOutputType = args.OutputType
    argsParquetFile = (
        Path(args.OutputGeneLocationParquet)
        if args.OutputGeneLocationParquet is not None
        else None
    )
    # by default, these settings should not be set
    vRebuildBlastDB = not args.NoRebuildBlastDB
    vRerunBlastP = not args.NoRerunBlastP

    error_message: list[str] = []
    if not os.path.isdir(vPGDB_FlatFileFolder):
        error_message.append(
            "Argument -pgdb " + f"{vPGDB_FlatFileFolder} was not a dirctory"
        )
    test_file = {
        vMD_reactions_File: "Metabolic domains reaction file (-rmdf)",
        vProtein_sequence_FastaFile: "Protein Sequence FASTA file (-psf)",
        vGeneTranscriptProtein_mapping_File: "Gene transcript-protein mapping file (-gtpf)",
        vGeneLocation_File: "Gene Location File (-glof)",
        vDNA_FastaFile: "DNA FASTA File (-dnaf)",
        vSignatureTailorFile: "Signature & Tailoring File (-sitf)",
    }
    test_writable_file = {
        gene_output_path: "Gene Output file (-gout)",
        vClusterOutputFile: "Cluster Output file (-cout)",
    }
    log.info("Check files for reading and writing")
    for filename, err in test_file.items():
        if not os.path.isfile(filename):
            error_message.append(f"Could not find {err} at {filename}.")
            continue
        try:
            with open(filename, "r") as f:
                pass
        except FileNotFoundError:
            error_message.append(f"File {filename} is a file that could not be found.")
        except PermissionError:
            error_message.append(
                f"You do not have permission to read the {err}" + " file at {filename}"
            )
        except OSError as e:
            error_message.append(f"Cannot read {err} {filename} because of {e}.")
    test_writable_file = {
        gene_output_path: "Gene Output file (-gout)",
        vClusterOutputFile: "Cluster Output file (-cout)",
    }
    for filename, err in test_writable_file.items():
        if os.path.isfile(filename):
            error_message.append(
                "File already exists at " + f"{filename}, ({err}). Cannot overwrite."
            )
            continue
        try:
            with open(filename, "a") as f:
                pass
        except FileNotFoundError:
            error_message.append(f"File {filename} is a file that could not be found.")
        except PermissionError:
            error_message.append(
                f"You do not have permission to write to the {err}"
                + " file at {filename}"
            )
        except OSError as e:
            error_message.append(
                "Cannot write to " + f"{err} {filename} because of {e}."
            )
    if len(error_message) != 0:
        for message in error_message:
            log.error(message)
        exit(1)

    vPartFinished = 0
    # Save temporary
    # Using scipy.io.savemat to save .mat files
    # Matlab save(filename,)
    if vTempsaves == 1:
        vPartFinished = _tempsaves(locals())["vPartFinished"]
    if os.path.isfile(vGeneLocation_File + "_tempsave.mat") and vRunAfterPart != 0:
        vPartFinished_old = vPartFinished
        vPartFinished = load(
            "".join([vGeneLocation_File, "_tempsave.mat"]), "vPartFinished"
        )["vPartFinished"]
        if vPartFinished >= vRunAfterPart:
            locals().update(load("".join([vGeneLocation_File, "_tempsave.mat"])))
            vPartFinished = vRunAfterPart
        elif vTempsavesOverwrite == 1:
            vPartFinished = vPartFinished_old
            del vPartFinished_old
        else:
            locals().update(load("".join([vGeneLocation_File, "_tempsave.mat"])))
            del vPartFinished_old

    # Part1 Read in Gene-Transcript-Protein File
    if vPartFinished < 1:
        log.info(
            "PART START: Read in file describing gene-ID, transcript-ID, and protein-ID conversion"
        )
        conversion = read_in_conversion_file(vGeneTranscriptProtein_mapping_File)
        vPartFinished = 1
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # Part 2: Map PGDB Gene-ids to reactions
    if vPartFinished < 2:
        log.info("PART START: Read in pgdb flat files")
        # ( vGenes_all, _, vReactions_all, vReactions_all_Att, vG2R,)
        genes_to_reactions = map_pgdb_reactions_to_genes(vPGDB_FlatFileFolder)
        vPartFinished = 2
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # Part 3: Map PGDB reactions to metabolic domains
    if vPartFinished < 3:
        log.info("PART START: Read in and map metabolic domains")
        # (vRxnMetaDom, vRxnMetaDom_Att, vR_MD) =
        annotated_metabolic_domains = get_metabolic_domains(
            vMD_reactions_File, list_of_MD_to_annotate
        )
        # Match PGDBs to MDs
        vRpgdb_Rmd = match_PGDBs_to_MDs(genes_to_reactions, annotated_metabolic_domains)
        #    vG2R.Gene2Rxn, vR_MD, vG2R.ReactionsAll,  # vReactions_all, vRxnMetaDom,)
        vPartFinished = 3
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # Part 4: Label Reactions with signature and tailoring
    if vPartFinished < 4:
        log.info(
            ["PART START: Label reactions with signature and tailoring information"]
        )
        # ( vSTClasses, vSigOrTail, vR_STClasses, vR_SigOrTail,)
        reaction_labels = label_reactions_with_sig_and_tail(
            genes_to_reactions.ReactionsAll, vSignatureTailorFile
        )
        vPartFinished = 4
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # % Map PGDB genes to Maping genes
    if vPartFinished < 5:
        log.info("PART START: Map pgdb genes to conversion file gene-IDs")
        vGenesDat_GeneInfo = [info.strip().upper() for info in vGenesDat_GeneInfo]
        map_genes_to_pgdb_genes = map_PGDBs_to_conversion_file(
            conversion,
            genes_to_reactions.GenesAll,
            vGenesDat_GeneInfo,
            vPGDBIdsToMap,
        )
        vPartFinished = 5
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # % Create MCL clustering mapped to Maping genes
    if vPartFinished < 6:
        if vMCLClusterFile == 0:
            log.info(
                "PART START: Create MCL clustering based on protein file"
                + " and map cluster IDs to conversion gene-IDs"
            )
            # [vG_MCL, vProtein_IDs] =
            vG_MCL = get_MCL_clustering(
                vProtein_sequence_FastaFile,
                conversion,
                vOverwriteMCLClustering,
                vRemoveNonProtLocations,
                vParallelize=vPara,
                rebuildBlastDB=vRebuildBlastDB,
                rerunBlastP=vRerunBlastP,
            )
        else:
            log.info(
                "PART START: Read in MCL clustering and map cluster IDs to conversion gene-IDs"
            )
            vG_MCL = read_MCL_clustering_file(
                Path(f"{vProtein_sequence_FastaFile}.MCLTAB.dump.I20"),
                vProtein_sequence_FastaFile,
                conversion,
                vRemoveNonProtLocations,
            )
        vPartFinished = 6
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # % Read in Gene Position File
    if vPartFinished < 7:
        log.info("PART START: Read in gene position file")
        gene_location = GeneLoc(file_path=vGeneLocation_File, unmasked=vUnmaskedDNA)
        log.debug("Gene Location File Attributes:")
        for i, attribute in enumerate(gene_location.columns):
            log.debug("Attribute %i: %s", i, attribute)
            # Remove protein or transcript information
            log.debug("Remove transcript information")

        gene_location.remove_protein_transcript_info(
            conversion,
        )
        # Consolidate gene information (one entry per gene)
        log.debug("Consolidate gene information")
        vUniqueChrom = gene_location.consolidate_gene_info()
        # Sort genes according to their physical positon on the chromosomes
        log.debug("Sort genes on the chromosomes")
        gene_location.sort_genes()
        vPartFinished = 7
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]
        gene_location = locals().pop("gene_location")

    # % Map Biomartfile to Gene Map
    if vPartFinished < 8:
        log.info("PART START: Map gene location to conversion gene-IDs")
        vBiomG_G = gene_location.map_to_conversion_file(conversion.MapG)

        # If chosen, remove gene position information about genes that are
        # not present in the protein fasta file.
        if vRemoveNonProtLocations == 1:
            log.debug("Remove non-protein coding genes in gene locations.")
            vBiomG_G = gene_location.remove_non_protein_fasta_gene_locations(
                conversion,
                vBiomG_G,
                vG_MCL.Protein_IDs,
            )
        vPartFinished = 8
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

        # % Get places on genome that do not have genes encoded:
    if vPartFinished < 9:
        log.info("PART START: Find intergenic regions")
        vInterChrom = gene_location.find_intergenic_regions()
        vPartFinished = 9
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # % Calculate percentile of gene size that should yield hypo genes
    if vPartFinished < 10:
        log.info(
            "PART START: Calculate size of sequence gap that should be populated by hypothetical genes"
        )
        possibly_none_quantile = (
            gene_location.data["Gene End (bp)"] - gene_location.data["Gene Start (bp)"]
        ).quantile(vHypoGenePercentile / 100.0)
        if possibly_none_quantile is None:
            raise Exception(
                "Sequence gap calculation failed: quantile calculation returned None. Check HypoGenePercentile"
            )
        vGeneLengthHypoMin = math.ceil(possibly_none_quantile)
        vPartFinished = 10
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # % Get sequencing gaps
    if vPartFinished < 11:
        log.info("PART START: Identify sequencing gaps (bases encoded by N)")
        if vGeneLengthHypoMin > 0 and vPreScreenGaps == 1:
            vSeqGaps = annotate_sequencing_gaps(
                vDNA_FastaFile,
                vInterChrom,  # intergenic region
                vUniqueChrom,
                vUnmaskedDNA,
                vPara,
                vOverwriteSeqGapInfo,
            )
        else:
            vSeqGaps = annotate_sequencing_gaps(
                vDNA_FastaFile,
                vInterChrom,
                vUniqueChrom,
                vUnmaskedDNA,
                vPara,
                vOverwriteSeqGapInfo,
            )
        vPartFinished = 11
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

        # % Introduce Hardbreaks (split chromosomes)
    if vPartFinished < 12:
        log.info("PART START: Introduce Hard Breaks")
        if vMaxSeqGapSize is not None:
            log.info(
                "Prevent clusters spanning sequencing gaps longer or equal than %ibp",
                vMaxSeqGapSize,
            )
        if len(vSeqGapSizesChromBreak) != 0:
            log.debug("SeqGapSizesChromBreak %s", vSeqGapSizesChromBreak)
            vNumbers = ", ".join([str(x) for x in vSeqGapSizesChromBreak])
            log.info(
                "Prevent clusters spanning sequencing gaps of length %sbp",
                vNumbers,
            )
        if vMaxInterGeneDist is not None:
            log.info(
                "Prevent clusters spanning intergenic distances longer or equal than %ibp",
                vMaxInterGeneDist,
            )
        if vPercentileForMaxInterGeneDist is not None:
            log.info(
                "Prevent clusters spanning the top %i%% intergenic distances",
                vPercentileForMaxInterGeneDist,
            )
        if vMaxInterGeneDistByMedian is not None:
            log.info(
                "Prevent clusters spanning intergenic distances longer or equal than %i * the median of gene sizes",
                vMaxInterGeneDistByMedian,
            )
            gene_location.introduce_hard_breaks(
                vSeqGaps,
                vMaxSeqGapSize,
                vMaxInterGeneDist,
                vPercentileForMaxInterGeneDist,
                vMaxInterGeneDistByMedian,
                vSeqGapSizesChromBreak,
            )

        vPartFinished = 12
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # % Populate Hypothetical genes
    if vPartFinished < 13:
        log.info(
            "PART START: Introduce Hypothetical genes for gaps longer than %i bp",
            vGeneLengthHypoMin,
        )
        if vInsertHypos:
            if "vHypoAmount" not in locals():
                # This should be impossible as long as args were parsed
                log.error("vHypoAmount does not exist. Needs to be set as an Argument.")
            else:
                hard_breaks2 = introduce_hypo_genes(
                    # hard_breaks1,
                    #  vGeneLengthHypoMin,
                    gene_location,
                    gene_location.find_intergenic_regions(),
                    vInsertHypos,
                    vHypoAmount,
                )

        vPartFinished = 13
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # % Annotate Biomartfile with Enzymes (use Gene-protein mapping)
    if vPartFinished < 14:
        log.info("PART START: Annotate gene location file with enzyme information")
        vBiomG_E = gene_location.annotate_enzyme_info(
            map_genes_to_pgdb_genes,
            genes_to_reactions,
        )
        log.debug(
            "Number of genes in location file identified as enzymes: %i",
            sum(vBiomG_E.values()),
        )
        vPartFinished = 14
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # % Compute LD information
    if vPartFinished < 15:
        log.info("PART START: Compute local gene duplication per gene location")
        # (vG_LD, vG_LD_ClustIDs) =
        gene_location.compute_local_duplication(
            vG_MCL,
            vBiomG_G,
            vKbRangeLD,
            vGeneRangeLD,
        )
        vPartFinished = 15
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # % Perform clustering
    if vPartFinished < 16:
        log.info("PART START: Perform clustering")
        # ( vClusters, vClusterBackground, vFinalStepsize, vStepSize_Max,)
        clusters = perform_clustering(
            gene_location,
            vG_MCL,
            vEnzMinForClusters,
            vKbRangeLD,
            vGeneRangeLD,
            vCriterion,
            vMinStepSize,
            vMaxStepSize,
        )
        vPartFinished = 16
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # % Perform cutoff
    if vPartFinished < 17:
        log.info("PART START: Perform cuttof at %i%%", top_percent_clusters)
        vClusters_Top = clusters.perform_cutoff(top_percent_clusters)
        vPartFinished = 17
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # % Produce Gene Output file
    if vPartFinished < 18:
        log.info("PART START: Produce GeneOutputFile")
        write_gene_output(
            gene_output_path,  # from args
            gene_location,
            clusters,
            top_percent_clusters,  # from args
            list_of_MD_to_annotate,  # from args
            genes_to_reactions,
            map_genes_to_pgdb_genes,
            vBiomG_E,
            vClusters_Top,
            annotated_metabolic_domains,
            vRpgdb_Rmd,
            reaction_labels,
            vOutputType,  # from args
        )
        vPartFinished = 18
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]

    # % Print cluster File
    if vPartFinished < 19:
        log.info("PART START: Produce ClusterOutputFile")
        write_cluster_output(
            vClusterOutputFile,
            gene_location,
            clusters,
            vClusters_Top,
            vOutputType,
        )
        if argsParquetFile is not None:
            gene_location.write_parquet(argsParquetFile)
        vPartFinished = 19
        if vTempsaves == 1:
            vPartFinished = _tempsaves(locals())["vPartFinished"]
    listener.stop()

    def help_content():
        """
        ---------------------------------------------------------------------------
         Function that defines the help output
        -------------------------------------------------------------------------
        Args:
           None

        Returns:
           vCellHelp: A cell array matrix containing the text that should be
                      printed out by the function display_help()
        """

        """
    %Matlab code
    vCellHelp = {'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%';
                 '% ';
                 '% PlantClusterFinder 1.3';
                 '% PlantClusterFinder (PCF) detects metabolic gene clusters in a sequenced';
                 '% genome. It uses a gene location file provided by the user (see below) and';
                 '% a PGDB created with Pathway Tools as well as further information (see';
                 '% below) to identify enzyme-coding genes (metabolic genes) located together';
                 '% on a chromosome. Initially only continous stretches of metabolic genes ';
                 '% lying directly next to each other are allowed. This condition is relaxed';
                 '% by iteratively increasing the intervening (non-metabolic) gene size by';
                 '% one. Several criteria to select for clusters are provided. In addition to';
                 '% this, clusters can be prevented from forming by a section of criteria.';
                 '% Details of PCF (version 1.0) can be found in PMID: 28228535.  ';
                 '% ';
                 '% The major differences between this version (1.3) and previous versions ';
                 '% (1.0 and 1.2) are:';
                 '% ';
                 '%   1) Physical breaks of the genome or sequencing gaps of unknown size are';
                 '%      typically encoded by stretches of Ns in the genome assembly fasta';
                 '%      file. Previously we inserted N hypothetical genes. This however ';
                 '%      diluted the background of low quality genomes with non-enzymes, and ';
                 '%      hence the likelyhood of a cluster to be classified as top x% of ';
                 '%      enzyme dense regions was better than in a genome that had good ';
                 '%      quality. In version 1.3 we identify these breaks and prevent ';
                 '%      formation of a cluster over these gaps. ';
                 '%   2) Any sequencing information that is missing is typically hard masked ';
                 '%      with Ns. Previously, any intergenic region affected by at least one ';
                 '%      N was evaluated for its length, and hypothetical genes were inserted ';
                 '%      accordingly (See Schlapfer et al, PMID: 28228535). This led ';
                 '%      sometimes to unrealistic prevention of detecting gene clusters. It';
                 '%      is unlikely that missing information about a single nucleotide would ';
                 '%      (if it would be known) lead to the finding of multiple gene models. ';
                 '%      Thus here we changed the code to insert 2 hypothetical genes only if ';
                 '%      a strech of unknown sequence is larger than nth percentile of gene ';
                 '%      sequences (set to 5). We also provide the option to NOT insert any ';
                 '%      hypothetical genes all together. Instead we by default we use ';
                 '%      MaxSeqGapSize set to 100000 and MaxInterGeneDistByMedian set to 50 ';
                 '%      resulting in similar clusters as in PCF version 1.0.';
                 '%   3) Large gene poor intergenic regions are present in genomes. In this ';
                 '%      version we provide the option of several parameters to prevent ';
                 '%      clusters from spanning such large gene poor regions.  ';
                 '% ';
                 '% Authors: Pascal Schlapfer, December 2017';
                 '%          Bo Xue, December 2017';
                 '% ';
                 '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%';
                 '% ';
                 '% Mandatory inputs: all files need full path, order does not matter.';
                 '%  ';
                 '%     -pgdb fullpath: Give full path of PGDB flat file folder, where the';
                 '%                     flat files of the species pgdb is stored.';
                 '%     -rmdf fullpath: TSV file of the metabolic domains of reactions. Check';
                 '%                     if file is outdated (needs to contain all reactions';
                 '%                     that are present in the PGDB that have genes';
                 '%                     annotated to).';
                 '%     -md list: A list of Metabolic domains that should be analyzed. For';
                 '%               example {''Amines and Polyamines Metabolism''; ''Amino Acids Metabolism''}';
                 '%     -psf fullpath: Protein sequence file of the genome. Can only contain';
                 '%                    Protein IDs as header e.g. >protein1-ID';
                 '%     -gtpf fullpath: TSV file with a header describing gene, transcript,';
                 '%                     and protein identifier and then for all the';
                 '%                     identfiers the listings (gene-ID, transcript-ID,';
                 '%                     protein-ID). This is used to map information';
                 '%                     regarding transcripts (currently not used) and';
                 '%                     proteins (used for MCL clustering, potentially also';
                 '%                     used in PGDBs) to the genes (e.g. the gene location';
                 '%                     information).';
                 '%     -glof fullpath: TSV file with gene ID, start bp, end bp, chromosome /';
                 '%                     scaffold, strand (encoded as 1 or -1). This file can';
                 '%                     be generated with a biomart or out of a gff3 file.';
                 '%     -dnaf fullpath: Fasta file of the HARDMASKED genome nucleotide';
                 '%                     sequences.';
                 '%     -sitf fullpath: TSV file describing the classes of enzymes that';
                 '%                     should be classified as signature or tailoring';
                 '%                     enzymes to identify clusters containing such.';
                 '%     -gout fullpath: Path and file name to the gene-Outputfile, make sure';
                 '%                     you have access to the file and folder. This is the';
                 '%                     first results output file that will be generated.';
                 '%     -cout fullpath: Path and file name to the cluster-Outputfile, make';
                 '%                     sure you have access to the file and folder. This is';
                 '%                     the second Results output file that will be';
                 '%                     generated.';
                 '%  ';
                 '% Optional Inputs:';
                 '%     ''GenesDat_GeneInfo'': needs to be followed by a list of Attributes in';
                 '%                          genes.dat. Defines where gene ID information';
                 '%                          should be searched for. Default is ''Unique-ID'',';
                 '%                          ''Accession-1'', ''Accession-2''';
                 '%     ''MCLClusterFile'': needs to be followed by either 1 or 0. With a default';
                 '%                       of 0. 1 indicates that the input file is not a';
                 '%                       protein fasta file, but a precomputed MCL';
                 '%                       clustering file (precomputing can be usefull for';
                 '%                       speed). If such a file is used, then nothing but';
                 '%                       gene IDs and tabs and new lines are allowed in this';
                 '%                       file.';
                 '%     ''HypoGenePercentile'': needs to be followed by a number. Defines the';
                 '%                           length of sequencing gap that should be';
                 '%                           populated by two hypothetical genes. If this is';
                 '%                           10, then two hypothetical genes are introduced';
                 '%                           as soon as the sequencing gap exceeds the';
                 '%                           (100-10) lower percentile of the background';
                 '%                           size of all genes of the genome (default is 5).';
                 '%     ''MaxSeqGapSize'': needs to be followed by a number. Maximal masked';
                 '%                      nucleotide region (in bp) that a cluster is allowed';
                 '%                      to bridge (default is 100000). ';
                 '%     ''MaxInterGeneDist'': needs to be followed by a number. Maximal';
                 '%                         intergenic distance ( in bp) that still can be';
                 '%                         crossed by a gene cluster. (default is -1, thus';
                 '%                         inactive).';
                 '%     ''MaxInterGeneDistByMedian'': needs to be followed by a number. Maximal';
                 '%                                 intergenic distance that still can be';
                 '%                                 crossed by a gene cluster, here defined';
                 '%                                 by this number times the median of the';
                 '%                                 gene sizes. (default is 50. Set it to -1';
                 '%                                 if it should be inactive).';
                 '%     ''PercentileForMaxInterGeneDist'': needs to be followed by a number.';
                 '%                                      Calculates the maximal intergenic';
                 '%                                      distance that still can be crossed';
                 '%                                      by a gene cluster. This number';
                 '%                                      determines the percentile to choose';
                 '%                                      from all the distances. If this is';
                 '%                                      99.9, then a cluster is allowed to';
                 '%                                      bridge if the largest intergenic';
                 '%                                      distance in the cluster is below';
                 '%                                      99.9% of the background size of all';
                 '%                                      intergenic distances of the genome';
                 '%                                      (default is -1, thus inactive).';
                 '%     ''SeqGapSizesChromBreak'': needs to be followed by a list of numbers.';
                 '%                              The masked nucleotides sequences of these';
                 '%                              specific lengths will prevent cluster';
                 '%                              formation (Default is an empty list: []).';
                 '%     ''PreScreenGaps'': needs to be followed by either 1 or zero. Defines if';
                 '%                      sequencing gaps should be prescreened to delete them';
                 '%                      if they are anyway too small to matter. (default is';
                 '%                      0)';
                 '%     ''OverwriteSeqGapInfo'': needs to be followed by either 1 or 0. Forces';
                 '%                            masked nucleotide analysis for finding';
                 '%                            sequencing gaps (default is 0, no';
                 '%                            enforcement).';
                 '%     ''OverwriteMCLClustering'': needs to be followed by either 1 or 0.';
                 '%                               Forces (if set to 1) that the MCL';
                 '%                               clustering is redone even if the result';
                 '%                               files are already present (default is 0).';
                 '%     ''Verbose'': needs to be followed by a number (0 is default, 1 gives';
                 '%                more screen-output, 2 gives exhaustive screen-output).';
                 '%     ''EnzMinForClusters'': needs to be followed by a number: Defines the';
                 '%                          minimum number of enyzmes that need to be';
                 '%                          present in a cluster (default is 3).';
                 '%     ''KbRangeLD'': needs to be followed by a number. Genes that are';
                 '%                  separated by more than this size (in bp) are not';
                 '%                  considered as local duplicates. (default is 100000)';
                 '%     ''GeneRangeLD'': needs to be followed by a number. Number of';
                 '%                    intervening genes that are allowed to separate two';
                 '%                    local duplicated genes. Two genes that are separated';
                 '%                    by more than this many intervening genes are not';
                 '%                    called local duplicates. (default is 10)';
                 '%     ''TopPercentClusters'': needs to be followed by a number. Clusters are';
                 '%                           labeled if they are within this top X% of';
                 '%                           enzyme dense clusters compared to the';
                 '%                           background. (default is 5)';
                 '%     ''MinStepSize'': needs to be followed by a number. Minimal intervening';
                 '%                    gene size that should be computed (default is 5).';
                 '%     ''MaxStepSize'': needs to be followed by a number. Maximal intervening';
                 '%                    gene size that should be computed (default is 20).';
                 '%     ''Criterion'': needs to be followed by a number. Criterion after which';
                 '%                  stepsize (intervening gene size) is chosen (default is';
                 '%                  2).';
                 '%         1: as published in Schlapfer et al 2017 (PMID:28228535)';
                 '%             sum(enzymes in clusters) - sum(non-enzymes in clusters) > 0.';
                 '%         2: stepsize 3 (because most organisms take this one)';
                 '%         3: sum(enzymes in clusters) - sum(non-enzymes in clusters, without hypo genes) -2*sum(hypogenes)';
                 '%         4: sum(clusters in stepsize n) > sum(clusters in stepsize n+1).';
                 '%     -help: Display help';
                 '%     -para: needs to be followed by a number: how many cpus can be used by';
                 '%            the algorithm (used for MCL clustering and Sequencing gap';
                 '%            search). The default is 1.';
                 '%     ''UnmaskedDNA'': needs to be followed by a number. Defines if plant';
                 '%                    Cluster finder should continue despite a non-masked';
                 '%                    genome sequence file was used (default is 0, not';
                 '%                    continue).';
                 '%     ''PGDBIdsToMap'': needs to be followed by a string. Can contain the';
                 '%                     Letters G and/or T and/or P. In that case (the pgdb';
                 '%                     is mapped to Gene (G) and/or Transcript (T) and/or';
                 '%                     Protein (P) Ids of the Gene conversion file. The';
                 '%                     default value is ''G''.';
                 '%     ''RunAfterPart'': Needs to be followed by a number. Sets the parts of';
                 '%                     the algorithm that should be recomputed. (default is';
                 '%                     0)';
                 '%     ''Tempsaves'': Defines if Temporary saves should be made. (default is';
                 '%                  0)';
                 '%     ''TempsavesOverwrite'': Defines if Temporary saves should be';
                 '%                           overwritten. (Default is 0).';
                 '%     ''RemoveNonProtLocations'': Defines if the protein-fasta file should be';
                 '%                               used to clean out the gene position file';
                 '%                               from sequences that do not show up on the';
                 '%                               protein fasta file, and thus do not have';
                 '%                               had a chance to be predicted to be an';
                 '%                               enzyme. (Default 0)';
                 '%     ''InsertHypos'': Defines if the genome should be populated by';
                 '%                    hypothetical genes in sequencing gaps that are not';
                 '%                    physical breaks. (Default false)';
                 '%     ''HypoAmount'': needs to be followed by a number. Defines by how many';
                 '%                   hypothetical genes a sequencing gap that is not a';
                 '%                   physical break should be populated by artificial genes';
                 '%                   to punish such sequencing gaps from becoming clusters.';
                 '%                   (Default 0)';
                 '%     ''OutputType'': Defines the format of the Outputfiles: ''old'' is using';
                 '%                   the old format to report all clusters of all stepsizes,';
                 '%                   ''verbose'' defines that top percent clusters of all';
                 '%                   stepsizes are reported. ''simple'' just reports the top';
                 '%                   percent clusters of the chosen stepsize (default';
                 '%                   simple).';
                 '% ';
                 '%  Outputs:';
                 '%     None (Two files are generated: vGeneOutputFile, containing all the';
                 '%          information about genes and vClusterOutputFile, containing all';
                 '%          the information about clusters)';
                 '% ';
                 '%     Shell, standalone (linux, but you can compile a windows version, see';
                 '%                       how to compile new version file):';
                 '%         To use standalone, download matlab runtime from matlab website.';
                 '%         path to your runtime.';
                 '%         ';
                 '%         USE ABSOLUTE PATHS!';
                 '% ';
                 '%         ./run_PlantClusterFinder.sh /share/apps/MATLAB/MATLAB_Runtime/v91 -pgdb "./csubellipsoidea/pgdb/csubellipsoideacyc/1.0/data/" -rmdf "./Inputs/ReactionMetabolicDomainClassification.txt" -md "{''Amines and Polyamines Metabolism''; ''Amino Acids Metabolism''; ''Carbohydrates Metabolism''; ''Cofactors Metabolism''; ''Detoxification Metabolism''; ''Energy Metabolism''; ''Fatty Acids and Lipids Metabolism''; ''Hormones Metabolism''; ''Inorganic Nutrients Metabolism''; ''Nitrogen-Containing Compounds''; ''Nucleotides Metabolism''; ''Phenylpropanoid Derivatives''; ''Polyketides''; ''Primary-Specialized Interface Metabolism''; ''Redox Metabolism''; ''Specialized Metabolism''; ''Sugar Derivatives''; ''Terpenoids''}" -psf "./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.protein.pcf13.fa" -gtpf "./csubellipsoidea/gtpCsubellipsoideaC_169_227_v2.0.annotation_info.txt.txt" -glof "./csubellipsoidea/gloCsubellipsoideaC_169_227_v2.0.gene.gff3.txt" -dnaf "./csubellipsoidea/CsubellipsoideaC_169_227_v2.0.hardmasked.fa" -sitf "./Inputs/scaffold-tailoring-reactions-05082016.tab" -gout "./csubellipsoidea/ARAGene1_3_memex.txt" -cout "./csubellipsoidea/ARAClust1_3_memex.txt" SeqGapSizesChromBreak ''[10000]'' PGDBIdsToMap GTP';
                 '%         % Note, PGDBIdsToMap is only needed here because the pgdb';
                 '%         % contains protein Ids in place of gene identifiers. ';
                 '% '}
    end
    """
        pass

    if __name__ == "__main__":
        main()
