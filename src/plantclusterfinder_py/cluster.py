import logging
import os
import subprocess
from collections import namedtuple
from pathlib import Path
from typing import IO, Any, Dict, Iterator, List, Optional, Tuple

import polars as pl

from plantclusterfinder_py.cluster_datatypes import (
    AnnotatedMetabolicDomains,
    Clusters,
    Conversion,
    GenesToReactions,
    MclClusterData,
    ReactionLabels,
)
from plantclusterfinder_py.genelocation import GeneLoc
from plantclusterfinder_py.get_positions_of_gap import get_positions_of_gap

_log = logging.getLogger(__name__)


def label_reactions_with_sig_and_tail(
    all_reactions_df: pl.DataFrame, signature_tailor_filename: Path
) -> ReactionLabels:
    """
    ---------------------------------------------------------------------------

      Function identifying whether reactions found in the pgdb are of a certain
      type of enzyme and if they are either signature or tailoring enzymes.
      Note, there have to be signature enzymes identified as 'signature' in the
      File signature_tailor_filename.

    ---------------------------------------------------------------------------
    Args:
        all_reactions_df: DataFrame containing all the information
                        (Attributes)that was searched for (described in
                        reactions_attributes argument, default is {'UNIQUE-ID', ...
                        'ENZYMATIC-REACTION', 'EC-NUMBER'};) in reactions.dat
                        of the pgdb.
        signature_tailor_file: TSV file describing the classes of enzymes that
                              should be classified as signature or tailoring
                              enzymes to identify clusters containing such.
    Return:
        ReactionLabels object with the signature and tailoring enzyme classes,
            list of types of enzymes, and mapping between reactions and
            enzyme classes as well as signature or tailoring enzymes
    """
    _log.debug(f"Reading in the signature tailor file {signature_tailor_filename}")
    stdata = pl.read_csv(
        signature_tailor_filename,
        has_header=True,
        separator="\t",
        encoding="ISO-8859-1",
        truncate_ragged_lines=True,
    )
    cols = stdata.columns
    _log.debug(f"read {len(stdata)} rows with columns {cols}")
    if len(stdata) != len(stdata["EC_number or RXN_ID"].unique()):
        _log.warning("signature tailor file has duplicate EC numbers or RXN IDs.")
    stclasses = stdata["Enzyme class"].unique()
    sigortail = stdata["Enzyme role"].unique()
    rxn_class = stdata.select("EC_number or RXN_ID", "Enzyme class").rows_by_key(
        "EC_number or RXN_ID", unique=True
    )
    rxn_sig_or_tail = stdata.select("EC_number or RXN_ID", "Enzyme role").rows_by_key(
        "EC_number or RXN_ID", unique=True
    )
    return ReactionLabels(
        STClasses=list(stclasses),
        SigOrTail=list(sigortail),
        Rxn_STClasses=rxn_class,
        Rxn_SigOrTail=rxn_sig_or_tail,
    )


def introduce_hypo_genes(
    gene_location: GeneLoc,
    sequencing_gaps: pl.DataFrame,
    hypothetical_gene_amount: int = 0,
    hypothetical_gene_gap_length_minimum: int = -1,
) -> GeneLoc:
    """
    ---

      Function that introduces Hypothetical genes if a sequencing gap reached a
      critical size hypothetical_gene_gap_length_minimum.
    ---
    Args:
        gene_location: Object with the gene location, chromosome regions, sequence gaps,
            and minimal length of bp of masked DNA that should be populated by two hypothetical
            genes
        sequence_gaps: DataFrame containing the location and size of sequence gaps.
        hypothetical_gene_gap_length_minimum: A size that defines the minimal length of a
                            sequencing gap that should still be processed. If
                            this is bigger than 0, then all sequencing gaps
                            smaller than this are deleted, unless they are
                            among the list of sequence_gaps, which
                            also needs to be set. Default is -1, not active.
    Optional Args:
        **hypothetical_gene_amount**: integer describing the number of hypothetical genes
            a sequencing gap that is not a physical break should be populated
            by artificial genes.
    Returns:
        **HardBreaks** object containing the gene location information, the
            inter-chromosomal matrix, and the HypoAmount populated.
    """
    _log.debug("Introducing hypothetical genes into sequence")
    # get the 5th and 7th columns of sequence_gaps
    if sequencing_gaps is None:
        raise Exception(
            "Sequencing Gaps must be populated to insert hypothetical genes"
        )
    gap_positions = (
        sequencing_gaps.select(pl.col(["Chromosome Name", "Start", "End"]))
        .with_columns(seq_gaps_length=pl.col("End") - pl.col("Start"))
        .with_columns(
            hypo_to_place=(
                pl.lit(hypothetical_gene_amount)
                if (
                    hypothetical_gene_gap_length_minimum < 0
                    or (
                        pl.col("seq_gaps_length")
                        >= pl.lit(hypothetical_gene_gap_length_minimum)
                    )
                )
                else pl.lit(0)
            )
        )
    )

    # WARNING: This is a slow way to process DataFrame
    ngaps = sequencing_gaps.height
    hypothetical_genes = []
    for i, row in enumerate(gap_positions.iter_rows()):
        _log.debug("Populate hypothetical genes %i of %i", i, ngaps)
        # no hypothetical genes to populate
        hypo_to_place = row[-1]
        if hypo_to_place == 0:
            continue
        place_last_gene_end = row[1]  # gap start position
        space_between_genes = row[3]  # gap length
        for hypo_i in range(hypo_to_place):
            place_to_insert = place_last_gene_end + int(
                space_between_genes / (hypo_to_place + 1) * hypo_i
            )
            new_row = (
                f"Hypothetical_Gene_{hypo_to_place + hypo_i}",
                place_to_insert,
                place_to_insert,
                row[0],  # chromosome/scaffold
                "-1",  # originally strand of real gene at start
                # FIXME: Strand "-" but no clear reason why (issue #27)
                True,
            )
            hypothetical_genes.append(new_row)
    hypo_gene_df = pl.DataFrame(
        hypothetical_genes,
        schema=(
            ("Gene Name", pl.String),
            ("Gene Start (bp)", pl.Int64),
            ("Gene End (bp)", pl.Int64),
            ("Chromosome Name", pl.String),
            ("Strand", pl.String),
            ("Hypothetical", pl.Boolean),
        ),
        orient="row",
    ).with_columns(Hypothetical=pl.lit(True))

    gene_location.insert_genes(hypo_gene_df)
    _log.debug(f"Inserted {len(hypo_gene_df)} hypothetical genes")
    return gene_location


def annotate_sequencing_gaps(
    nucleotide_fasta_filename: Path,
    intergenic_location: pl.DataFrame,
    unique_chromosome_names: List[str],
    is_dna_unmasked: bool = True,
    # parallelize: int = 1,
    overwrite_sequence_gap_info_file: bool = False,
    sequence_gap_sizes_chromosome_breaks: List[int] = [],
    hypothetical_gene_gap_length_minimum: int = -1,
) -> pl.DataFrame:
    """
    ---------------------------------------------------------------------------

      Functions that uses the output of get_sequencing_gaps and reformulates
      it such that it is clear in which intergenic region the sequencing gap
      lies.

    ---------------------------------------------------------------------------

    Args:
        nucleotide_fasta_filename: FASTA file of the masked genome nucleotide sequences.
        intergenic_location: A DataFrame describing
                        - chromosome,
                        - start bp,
                        - end bp,
                        - gene before,
                        - and gene after
                      an intergenic region.
        unique_chromosome_names: A list with unique entries of all chromosome and scaffold
                      names.
    Optional Args:
        - **is_dna_unmasked**: Option to ignore the error if you use unmasked DNA files.
        - **parallelize**: Number of parallel processes to run (deprecated).
        - **hypothetical_gene_gap_length_minimum**: If set to True, the code enforces masked nucleotide
                              analysis for finding sequencing gaps (default is
                              False, no enforcement).
        - **vSeqGapSizesChromBreak**: A list of numbers. The masked nucleotides
                                sequences of these specific lengths are
                                preventing cluster formation (Default is an
                                empty list []). If this is given, then also
                                vGeneLengthHypoMin has to be set. Will lead to
                                prescreening of the sequencing gaps.
        - **vGeneLengthHypoMin**: A size that defines the minimal length of a
                            sequencing gap that should still be processed. If
                            this is bigger than 0, then all sequencing gaps
                            smaller than this are deleted, unless they are
                            among the list of vSeqGapSizesChromBreak, which
                            also needs to be set. Default is -1, not active.
    Returns:
        **DataFrame** containing the
                  - chromosome/scaffold,
                  - the identification if it is coding (A) or masked (N) DNA
                  - start bp and
                  - length of a masked region.
                  In addition there is the information stored in with
                  intergenic region the gap is encoded and its start and end
                  bp, corrected for maximal size (= intergenic region).
    """

    if not (
        hypothetical_gene_gap_length_minimum == -1
        or sequence_gap_sizes_chromosome_breaks == []
    ):
        _log.error(
            "annotate_sequencing_gaps needs both SeqGapSizesChromBreak and GeneLengthHypoMin to be set or to be unset."
        )
        raise Exception("Incompatable parameters in annotate_sequencing_gaps function")
    _log.debug("Annotating sequence gaps")

    gap_output_file = nucleotide_fasta_filename.with_suffix(
        nucleotide_fasta_filename.suffix + "_GAPOutput"
    )
    if overwrite_sequence_gap_info_file or not gap_output_file.exists():
        get_sequencing_gaps(
            nucleotide_fasta_filename,
            1,
            overwrite_sequence_gap_info_file,
            is_dna_unmasked,
        )
    sequencing_gaps = (
        pl.read_csv(
            gap_output_file,
            separator=" ",
            has_header=False,
            new_columns=["scaffold", "amino acid", "start index", "length"],
        )
        .with_columns(
            (pl.col("start index") + pl.col("length") - pl.lit(1)).alias("end index")
        )
        .filter(
            pl.col("amino acid") == pl.lit("N")
        )  # Do not process all non N sequences.
    )
    for bad_chromosome, *_ in (
        sequencing_gaps.filter(pl.col("scaffold").is_in(unique_chromosome_names).not_())
        .select(pl.col("scaffold").unique())
        .iter_rows()
    ):
        _log.warning(
            f"DNA Fasta file '{
                nucleotide_fasta_filename
            }' contains chromosomes or scaffolds not covered with the "
            + f"gene position file: {bad_chromosome}"
        )
    sequencing_gaps = sequencing_gaps.filter(
        pl.col("scaffold").is_in(unique_chromosome_names)
    )  # Do not process sequences that do not lie on a chromosome we want
    sequencing_gaps = (
        # All predicates of join_where() must be true for join to happen
        sequencing_gaps.join_where(
            intergenic_location,
            pl.col("Chromosome Name") == pl.col("scaffold"),
            pl.col("Start")  # intergenic region start
            # intergenic region starts before gap start
            <= pl.col("start index"),
            pl.col("End")  # intergenic region end
            >= pl.col("start index"),  # intergenic region ends after gap start
            pl.col("Start")  # intergenic region start
            <= pl.col("end index"),  # intergenic region starts before gap end
            pl.col("End")  # intergenic region end
            >= pl.col("end index"),  # intergenic region ends after gap end
        )
        .group_by(["scaffold", "start index", "end index"])
        .agg(
            seq_gaps_start=pl.col("Start").min(),
            seq_gaps_end=pl.col("End").max(),
        )
        .with_columns(
            seq_gaps_length=pl.col("seq_gaps_end") - pl.col("seq_gaps_start") + 1
        )
    )
    # find the intergenic region of the gap, if the gap is present within a
    # gene model, get rid of the gap.
    # FIXME: validate logic: (issue #28)
    # only gaps which are in intergenic regions
    # & indicate which intergenic region
    return sequencing_gaps


def get_sequencing_gaps(
    vFile: Path,
    vParallelize: int,
    vOverwriteSeqGapInfo: bool = False,
    vUnmaskedDNA: bool = False,
) -> None:
    """
    ---------------------------------------------------------------------------

      Function that screens a DNA FASTA file for masked stretches of DNA
      (masked by N).

      Creates a "..._GAPOutput" file from the FASTA file, then creates a
      "..._GAPOutput_count.txt" file by running analyze_plant_cluster_gap_file()

    ---------------------------------------------------------------------------
    Args:
        vFile: FASTA file of the masked genome nucleotide sequences.
        vParallelize: Defines if multiple cpus should be used for computation of parts
               of algorithm. Active if larger than 1.
        vOverwriteSeqGapInfo: If set to True, the code enforces masked nucleotide
                              analysis for finding sequencing gaps (default is
                              False, no enforcement).
        vUnmaskedDNA: Option to ignore the error if you use unmasked DNA files.
    Returns:
        None (File is generated with the path of the first input argument,
             modifying the gene name)
    Raises:
        IOError: if file cannot be generated.
    """
    if vOverwriteSeqGapInfo:
        if vFile.exists():
            _log.error("Cannot overwrite existing _GAPOutput with tag set")
            raise Exception("Overwrite Protection Enforcement of Existing Gap File")
    _log.debug(f"Writing sequence gap file {vFile}")
    get_positions_of_gap(
        vFile, vFile.with_suffix(vFile.suffix + "_GAPOutput"), parallelize=vParallelize
    )
    analyze_plant_cluster_gap_file(
        vFile.with_suffix(vFile.suffix + "_GAPOutput"), vUnmaskedDNA
    )


def get_metabolic_domains(
    vGenesMetabolicDomains_File: Path, vMD_to_annotate: str
) -> AnnotatedMetabolicDomains:
    """
    ---------------------------------------------------------------------------

      Function that reads in the metabolic domains of reactions
      and gives a mapping between reactions in these files and the metabolic
      domains that were asked to be annotated (usually a subset of all) by
      vMD_to_annotate.

    ---------------------------------------------------------------------------

    Args:
        vGenesMetabolicDomains_File: TSV of the metabolic domains of reactions.
                                     Check if file is outdated (needs to
                                     contain all reactions that are present in
                                     the PGDB that have genes annotated to).
        vMD_to_annotate: A list of Metabolic domains that should be analyzed.
                         For example {'Amines and Polyamines Metabolism'; ...
                         'Amino Acids Metabolism'

    Returns:
        AnnotatedMetabolicDomains object with all items stored in Metabolic domain file aw well as
            mappings for metabolic domains of interest described in `MD_to_annotate`.
    """
    _log.debug(f"Getting metabolic domains from file {vGenesMetabolicDomains_File}")
    # TODO: Why does this file have \r EOL char? (issue #29)
    rxn_meta_dom = (
        pl.read_csv(vGenesMetabolicDomains_File, separator="\t", eol_char="\r")
        .with_columns(pl.col("Reaction ID").str.strip_chars(" ").alias("Reaction ID"))
        .unpivot(
            [f"Class{n + 1}" for n in range(6)] + ["SM1A"],
            index=["Reaction ID", "Pathway"],
            value_name="Metabolic Domain",
        )
        .filter(pl.col("Metabolic Domain") != "NA")
        .select(pl.col("Reaction ID"), pl.col("Pathway"), pl.col("Metabolic Domain"))
    )

    return AnnotatedMetabolicDomains(
        RxnMetaDom=rxn_meta_dom, Rxn2MetaDom=rxn_meta_dom.rows_by_key(key="Reaction ID")
    )


def read_MCL_clustering_file(
    vFile: Path,
    vProtein_sequence_FastaFile: Path,
    conversion: Conversion,
    vRemoveNonProtLocations: bool,
) -> MclClusterData:
    """
    ---------------------------------------------------------------------------

      Function that reads out the MCL clustering information and gives every
      gene one or several MCL clusters. This is then used for local
      duplication.

    ---------------------------------------------------------------------------

    Args:
        vFile: MCL-clustering output file (I20).
        vProtein_sequence_FastaFile: The FASTA file of the protein sequences.
        conversion: Lists of gene, and protein IDs with mapping from conversion files.
        vRemoveNonProteinLocations: True if protein-IDs should be stored
            in return value.

    Returns:
        MclClusterData object which provides mapping between
        genes in the conversion-ID file and the MCL-clustering file,
        as well as a list of protein-IDs stored in the fasta file.
    """
    _log.debug(f"Reading in MCL Clustering File {vFile}")
    if vRemoveNonProtLocations:
        _protein_init = 300_000
    protein_ids = []
    with open(vProtein_sequence_FastaFile, "r") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                protein_ids.append(line[1:].strip())

    mcl_cluster_data: List[Dict[str, int | str | List[str]]] = []
    clusters_by_p: Dict[str, List[str]] = dict()
    clusters_by_g: Dict[str, List[str]] = dict()

    # Read in MCL clusters, and assign ID based on index
    with open(vFile, "r") as f:
        for i, line in enumerate(f):
            # if i <= init_max:
            mcl_id = f"MCL_{i}"
            prots_in_cluster = line.strip().split()
            # Convert Transcript ID to Gene Name
            genes_in_cluster = [conversion.MapP_MapG[p] for p in prots_in_cluster]
            _log.debug(f"{mcl_id}: {prots_in_cluster}")
            if any([len(x) for x in genes_in_cluster]) == 0:
                _log.error(f"Empty genes cluster in MCL {mcl_id}")
                missing_genes = [
                    p for p in prots_in_cluster if len(conversion.MapP_MapG[p])
                ]
                _log.error(conversion.MapP_MapG)
                raise Exception(f"Genes are missing for proteins {missing_genes}")
            genes_in_cluster = [x[0] for x in genes_in_cluster]
            mcl_cluster_data.append(
                {
                    "ID": mcl_id,
                    "Transcript_IDs": prots_in_cluster,
                    "Gene_Names": genes_in_cluster,
                    "Cluster_Size": len(prots_in_cluster),
                }
            )
            for prot in prots_in_cluster:
                clusters_by_p[prot] = clusters_by_p.get(prot, []) + [mcl_id]
            for gene in genes_in_cluster:
                clusters_by_g[gene] = clusters_by_g.get(gene, []) + [mcl_id]

    proteins_list = list(clusters_by_p.keys())

    return MclClusterData(
        Gene2MCLClusters=clusters_by_g,
        MCL_Clusters=pl.DataFrame(mcl_cluster_data),
        Protein_IDs=proteins_list,
    )


def get_MCL_clustering(
    vProtein_sequence_FastaFile: Path,
    conversion: Conversion,
    vOverwriteMCLClustering: bool,
    vRemoveNonProtLocations: bool,
    vParallelize: int = 0,
    rebuildBlastDB: bool = True,
    rerunBlastP: bool = True,
) -> MclClusterData:
    """
    ---------------------------------------------------------------------------

      Function that uses a protein fasta file of a genome annotation to create
      a blastdb, and then does an *all against all* blast and uses the result to
      do a MCL clustering. This is used for local duplication information.

      MCL algorithm is for **Markov Cluster Algorithm,** which is an unsupervised
      clustering algorithm for graphs. For details, the project is hosted at
      https://micans.org/mcl/

    ---------------------------------------------------------------------------

    Args:
        vProtein_sequence_FastaFile: Protein sequence file of the genome.
                                     Should only contain Protein IDs as
                                     header e.g. >protein1-ID
        Conversion:
            vMapG: Unique list of gene IDs stored in the conversion file.
            vMapP: Unique list of protein IDs stored in the conversion file.
            vMapG_MapP: Mapping between vMapG and vMapP
        vOverwriteMCLClustering: Whether or not to overwrite the MCL clustering file
        vRemoveNonProtLocations,
        vParallelize: Defines if multiple cpus should be used for computation of parts
               of algorithm. Active if larger than 1.
        rebuildBlastDB: If BlastDB exists, should it be rebuilt? (recommended)
        rerunBlastP: if BlastP output exists, should it be rerun?

    Returns:
        MclClusterData object which provides mapping between
        genes in the conversion-ID file and the MCL-clustering file,
        as well as a list of protein-IDs stored in the fasta file.
    """
    _log.debug("Getting MCL clustering by running blastp and mcl")
    int_cpu_count: int = os.cpu_count() or 1
    if vParallelize == 0:
        vParallelize = max(1, int_cpu_count - 1)

    apps = ["makeblastdb", "blastp", "cut", "mcxload", "mcl", "mcxdump"]
    missing_apps: List[str] = []
    app_path: Dict[str, Path] = dict()
    # verify all apps are in $PATH, generate error if so
    for app in apps:
        output: subprocess.CompletedProcess = subprocess.run(
            ["which", app], capture_output=True
        )
        if output.returncode > 0:
            missing_apps.append(app)
        else:
            app_path[app] = Path(str(output.stdout, "utf-8").strip())
    if len(missing_apps) > 0:
        missing_apps_msg = ", ".join(missing_apps)
        raise Exception(
            f"The following applications are not in $PATH: {missing_apps_msg}"
        )

    command_for: Dict[str, List[str]] = dict()
    blastdb_file = f"{vProtein_sequence_FastaFile}.blastDB"
    if not (rebuildBlastDB or os.path.isfile(blastdb_file)):
        command_for["makeblastdb"] = [
            "-in",
            str(vProtein_sequence_FastaFile),
            "-dbtype",
            "prot",
            "-out",
            blastdb_file,
            "-hash_index",
        ]
    else:
        command_for["makeblastdb"] = ["-h"]
    blastp_file = f"{vProtein_sequence_FastaFile}.blastOUT"
    if not (rerunBlastP or os.path.isfile(blastp_file)):
        command_for["blastp"] = [
            "-query",
            str(vProtein_sequence_FastaFile),
            "-db",
            blastdb_file,
            "-evalue",
            "1",
            "-out",
            blastp_file,
            "-outfmt",
            "6 qseqid sseqid evalue bitscore length pident qstart qend sstart send",
        ]
        if vParallelize > 1:
            command_for["blastp"].extend(["-num_threads", str(vParallelize)])
    else:
        command_for["blastp"] = ["-h"]
    command_for["cut"] = [
        "-f",
        "1,2,3",
        blastp_file,
    ]
    cut_output = f"{vProtein_sequence_FastaFile}.MCLIN"
    command_for["mcxload"] = [
        "-abc",
        cut_output,
        "--stream-mirror",
        "--stream-neg-log10",
        "-stream-tf",
        "ceil(200)",
        "-o",
        f"{vProtein_sequence_FastaFile}.MCLMATRIX",
        "-write-tab",
        f"{vProtein_sequence_FastaFile}.MCLTAB",
    ]
    command_for["mcl"] = [
        f"{vProtein_sequence_FastaFile}.MCLMATRIX",
        "-I",
        "2",
        "-o",
        f"{vProtein_sequence_FastaFile}.MCLMATRIX.I20",
    ]
    overwrite_dump_workaround = "" if vOverwriteMCLClustering else f"{os.getpid()}."
    command_for["mcxdump"] = [
        "-icl",
        f"{vProtein_sequence_FastaFile}.MCLMATRIX.I20",
        "-tabr",
        f"{vProtein_sequence_FastaFile}.MCLTAB",
        "-o",
        f"{vProtein_sequence_FastaFile}.{overwrite_dump_workaround}MCLTAB.dump.I20",
    ]

    for app in apps:
        _log.info(f"{app_path[app]} {' '.join(command_for[app])}")
        if app == "cut":
            with open(cut_output, "w") as f:
                subprocess.run(
                    [str(app_path[app])] + command_for[app], check=True, stdout=f
                )
        else:
            subprocess.run(
                [str(app_path[app])] + command_for[app],
                check=True,
            )

    return read_MCL_clustering_file(
        Path(
            f"{vProtein_sequence_FastaFile}.{overwrite_dump_workaround}MCLTAB.dump.I20"
        ),
        vProtein_sequence_FastaFile,
        conversion,
        vRemoveNonProtLocations,
    )


def map_pgdb_reactions_to_genes(
    vPGDB_FlatFileFolder: Path, **kwargs
) -> GenesToReactions:
    """
    ---------------------------------------------------------------------------

      Function reading out information of pgdb flatfiles genes.dat,
      proteins.dat, enzrxns.dat and reactions.dat, mapping genes
      to reactions directly.

    ---------------------------------------------------------------------------

    Args:
        vPGDB_FlatFileFolder: Full path to the folder where the flatfiles of
                              the species pgdb is stored.
    Optional Args:
        - **GenesAttributes** Attributes of genes.dat that should be read out.
                           Default is {'UNIQUE-ID','PRODUCT','ACCESSION-1', ...
                           'ACCESSION-2'}
        - **ProteinsAttributes**: Attributes of proteins.dat that should be read
                               out. Default is {'UNIQUE-ID','CATALYZES', ...
                               'COMPONENTS'}
        - **EnzrxnsAttributes**: Attributes of enzrxns.dat that should be read
                             out. Default is {'UNIQUE-ID','REACTION','ENZYME'}
        - **ReactionsAttributes**: Attributes of reactions.dat that should be read
                               out. Default is {'UNIQUE-ID',
                               'ENZYMATIC-REACTION', 'EC-NUMBER'}
        - **G2P**: Defines which attributes in genes.dat should be matched to which
               attributes in proteins.dat. Default is: ('UNIQUE-ID', 'PRODUCT')
        - **P2P**: Defines which attributes in proteins.dat should be matched to
               which attributes in proteins.dat (Protein complexes). Default
               is: ('UNIQUE-ID', 'COMPONENTS')
        - **P2E**: Defines which attributes in proteins.dat should be matched to
               which attributes in enzrxns.dat. Default is: ('CATALYZES', 'COMPONENTS')
        - **E2R**: Defines which attributes in enzrxns.dat should be matched to
               which attributes in reactions.dat. Default is: ('REACTION', 'UNIQUE-ID')

    Returns:
        **GenesToReactions** object containing
            - GenesAll: the searched-for attributes of all genes in genes.dat of pgdb
            - ReactionsAll: searched-for attributes of all reactions in reactions.dat of pgdb
            - Gene2Rxn: Mapping between records in GenesAll and ReactionsAll
    """
    _log.debug("Mapping PGDB reactions to genes")
    # Attributes of Flatfiles that should be read out
    _genes_attributes = tuple(
        kwargs.get(
            "GenesAttributes", ("UNIQUE-ID", "PRODUCT", "ACCESSION-1", "ACCESSION-2")
        )
    )
    _proteins_attributes = tuple(
        kwargs.get("ProteinsAttributes", ("UNIQUE-ID", "CATALYZES", "COMPONENTS"))
    )
    _enzrxns_attributes = tuple(
        kwargs.get("EnzrxnsAttributes", ("UNIQUE-ID", "REACTION", "ENZYME"))
    )
    _reactions_attributes = kwargs.get(
        "ReactionsAttributes", ("UNIQUE-ID", "ENZYMATIC-REACTION", "EC-NUMBER")
    )
    G2P = kwargs.get("G2P", ("PRODUCT", "UNIQUE-ID"))
    P2P = kwargs.get("P2P", ("UNIQUE-ID", "COMPONENTS"))
    # ("CATALYZES", "COMPONENTS"))
    P2E = kwargs.get("P2E", ("CATALYZES", "UNIQUE-ID"))
    E2R = kwargs.get("E2R", ("REACTION", "UNIQUE-ID"))
    # R2E = kwargs.get("R2E", ("UNIQUE-ID", "REACTION"))
    # this was in the original code, but it's not used
    # R2Pa = ("UNIQUE-ID", "ENZYMATIC-REACTION")

    # Guarantee that UNIQUE-ID is in the attributes
    genes_attributes = (
        _genes_attributes
        if "UNIQUE-ID" in _genes_attributes
        else _genes_attributes + ("UNIQUE-ID",)
    )
    proteins_attributes = (
        _proteins_attributes
        if "UNIQUE-ID" in _proteins_attributes
        else _proteins_attributes + ("UNIQUE-ID",)
    )
    enzrxns_attributes = (
        _enzrxns_attributes
        if "UNIQUE-ID" in _enzrxns_attributes
        else _enzrxns_attributes + ("UNIQUE-ID",)
    )
    reactions_attributes = (
        _reactions_attributes
        if "UNIQUE-ID" in _reactions_attributes
        else _reactions_attributes + ("UNIQUE-ID",)
    )

    # Read in all attributes and store them
    _log.info("Reading in genes.dat")
    genes_all = read_in_pgdb_flat(
        vPGDB_FlatFileFolder / "genes.dat",
        genes_attributes,
        unique_keys=genes_attributes,
    )
    _log.info("Reading in proteins.dat")
    proteins_all = read_in_pgdb_flat(
        vPGDB_FlatFileFolder / "proteins.dat", proteins_attributes
    )
    _log.info("Reading in enzrxns.dat")
    enzrxns_all = read_in_pgdb_flat(
        vPGDB_FlatFileFolder / "enzrxns.dat", enzrxns_attributes
    )
    _log.info("Reading in reactions.dat")
    reactions_all = read_in_pgdb_flat(
        vPGDB_FlatFileFolder / "reactions.dat", reactions_attributes
    )

    # Match flat files to each other
    _log.info("Matching genes to proteins")
    G2_P = match_pgdb_flat(
        genes_all, proteins_all, G2P, primary_id="GENE-ID", secondary_id="PROTIEN-ID"
    )
    _log.info("Matching proteins to protein-complexes")
    P_P = match_pgdb_flat(
        proteins_all,
        proteins_all,
        P2P,
        primary_id="PROTIEN-ID",
        secondary_id="COMPONENTS",
    )

    _log.info("Match proteins and protein complexes to enzymereactions")
    P_E = match_pgdb_flat(
        proteins_all,
        enzrxns_all,
        P2E,
        primary_id="PROTIEN-ID",
        secondary_id="ENZYME-ID",
    )
    _log.info("Matching enzymereactions to reactions")
    E2_R = match_pgdb_flat(
        enzrxns_all,
        reactions_all,
        E2R,
        primary_id="ENZYME-ID",
        secondary_id="REACTION-ID",
    )

    # Match Genes to reactions directly.
    _log.info("Match genes to to reactions")
    # initialize matrix the size of [G2_P, E_R2]
    # gene2rxn = (
    # gene_prot = G2_P.join(P_P, on="PROTIEN-ID", validate="1:1")
    gene_prot = G2_P.join(P_P, on="PROTIEN-ID")
    gene_prot_enz = gene_prot.join(P_E, on="PROTIEN-ID")
    gene_prot_enz_reaction = gene_prot_enz.join(E2_R, on="ENZYME-ID")
    G2_R = gene_prot_enz_reaction.select(pl.col("GENE-ID"), pl.col("REACTION-ID"))
    gene2rxn = G2_R.rows_by_key("GENE-ID")

    return GenesToReactions(
        GenesAll=genes_all, ReactionsAll=reactions_all, Gene2Rxn=gene2rxn
    )


def match_pgdb_flat(
    flat1: pl.DataFrame,
    flat2: pl.DataFrame,
    join_on: tuple[str, str],
    primary_id: str = "UNIQUE-ID",
    secondary_id: str = "UNIQUE-ID_right",
) -> pl.DataFrame:
    """
    ---------------------------------------------------------------------------

      Function that maps two flat files of a pgdb to each other.

    ---------------------------------------------------------------------------

    Args:
        flat1: DataFrame of a read out pgdb Flat file.
        flat2: DataFrame of a read out pgdb Flat file.
        join_on: Definition of which attributes in structure 1 should be
                    matched to which attributes of structure 2

    Returns:
        F1_F2: A mapping between items of the first pgdb file to the second
                one.
    """
    _log.debug(f"Matching PGDB flat files {flat1} and {flat2}")
    if join_on[0] not in flat1.columns:
        _log.error(
            f"Matching column {join_on[0]} not in input data"
            + f" {flat1.columns}, note the other table's columns are {flat2.columns}"
        )
        raise Exception("Column-matching Mismatch")
    if join_on[1] not in flat2.columns:
        _log.error(
            f"Matching column {join_on[1]} not in input data"
            + f" {flat2.columns}, note the other table's columns are {flat1.columns}"
        )
        raise Exception("Column-matching Mismatch")
    if flat1[join_on[0]].dtype == pl.Null:
        _log.error(f"Primary table join ID {[join_on[0]]} is null dtype, cannot join")
        _log.info(f"Columns of flats are are {flat1.columns}, {flat2.columns}")
        raise Exception("Column-matching Mismatch")
    if flat2[join_on[1]].dtype == pl.Null:
        _log.warning(
            f"Cannot map {join_on[0]} to {join_on[1]}"
            + ": dtype is Null, returning Null for second column"
        )
        _log.info(f"Columns of flats are are {flat1.columns}, {flat2.columns}")
        return flat1.select(
            pl.col(join_on[0]).alias(primary_id), pl.lit(None).alias(secondary_id)
        )
    flat2 = flat2.explode(join_on[1]) if flat2[join_on[1]].dtype == pl.List else flat2
    flat1 = flat1.explode(join_on[0]) if flat1[join_on[0]].dtype == pl.List else flat1
    # NOTE: Validation on join can prevent data issues, consider adding
    joined_dbs = flat1.join(
        flat2,
        left_on=join_on[0],
        right_on=join_on[1],
        suffix="_right",
    )
    _log.debug(f"SUCCESSFUL JOIN: {joined_dbs}")
    return joined_dbs.select(
        pl.col("UNIQUE-ID").alias(primary_id), pl.col(join_on[0]).alias(secondary_id)
    )


def read_in_conversion_file(
    vGeneTranscriptProtein_mapping_File: Path,
) -> Conversion:  # , vPara, varargin):
    """
    ---------------------------------------------------------------------------

      Function that reads in the gene, transcript, and protein ID conversion
      file and stores the IDs as unique lists and mapping between genes and
      transcripts aswell as genes and proteins.

    ---------------------------------------------------------------------------
    Args:
        vGeneTranscriptProtein_mapping_File: TSV file with a header describing
                                             gene, transcript, and protein
                                             identifier and then for all the
                                             identfiers the listings (gene-ID,
                                             transcript-ID, protein-ID).

    Returns:
        Conversion object containing the IDs from the conversion file
    Raises
    """
    _log.debug(f"Reading in conversion file {vGeneTranscriptProtein_mapping_File}")
    conversion_data = pl.read_csv(
        vGeneTranscriptProtein_mapping_File,
        separator="\t",
        schema_overrides={
            x: pl.String for x in ["Gene Name", "Transcript Name", "Protein Name"]
        },
    )
    column_names = conversion_data.columns
    conversion_len = len(conversion_data)
    _log.debug(
        "Read in conversion data, %d rows, columns are %s", conversion_len, column_names
    )
    conversion = Conversion(df=conversion_data)
    # Validate conversion before output
    # TODO: are unclean protein ID files an unrecoverable error? (issue #30)
    if len(conversion.MapG_MapT) < len(conversion.MapG):
        _log.warning("Some Transcripts match multiple genes in the ID Conversion file")
        transcript_multi = (
            conversion_data.group_by(column_names[1])
            .n_unique()
            .filter(pl.col(column_names[0]) > 1)[column_names[1]]
        )
        for transcript in transcript_multi:
            _log.warning("Please check the transcript %s", transcript)
        _log.error("Please clean the protein ID file")
    if len(conversion.MapG_MapP) < len(conversion.MapG):
        _log.warning("Some Proteins match multiple genes in the ID Conversion file")
        prot_multi = (
            conversion_data.group_by(column_names[2])
            .n_unique()
            .filter(pl.col(column_names[0]) > 1)[column_names[2]]
        )
        for prot in prot_multi:
            _log.warning("Please check the protein %s", prot)
        _log.error("Please clean the protein ID file")
    _log.info(
        "Number of Gene-ID to Transcript-ID mappings: %i", len(conversion.MapG_MapT)
    )
    _log.info("Number of Gene-ID to Protein mappings: %i", len(conversion.MapG_MapP))
    return conversion


def match_PGDBs_to_MDs(
    Gene2Rxn: GenesToReactions,
    annotated_metabolic_domains: AnnotatedMetabolicDomains,
) -> dict:
    """
    ---------------------------------------------------------------------------

      Function that maps the reactions stored in the pgdbs to the reactions
      stored in the metabolic domain file.

    ---------------------------------------------------------------------------
    Args:
        Gene2Rxn: Mapping of all genes in the Genes_all to the Reactions in
              Reactions_all
        annotated_metabolic_domains: Mapping between all Reactions stored in the Metabolic domain file
               and the metabolic domains that were identified to be of interest
               in vMD_to_annotate.
               And DataFrame of all the items stored in the Metabolic domain file

    Returns:
        vRpgdb_Rmd: A mapping between reactions stored in the pgdb and the
                    reactions given in the metabolic domain file. Used by write_gene_output.
    """
    _log.debug("Matching PGDBs to annotated metabolic domains")

    return (
        Gene2Rxn.ReactionsAll.join(
            annotated_metabolic_domains.RxnMetaDom,
            left_on="UNIQUE-ID",
            right_on="Reaction ID",
        )
        .select(pl.col("UNIQUE-ID").alias("REACTION-ID"), pl.col("Metabolic Domain"))
        .rows_by_key(key="REACTION-ID")
    )


def map_PGDBs_to_conversion_file(
    Conversion: Conversion,
    all_genes_df: pl.DataFrame,
    genes_data_metadata: List[str],
    pgdb_ids_to_map: str,
) -> dict:  # NOTE: this seems to be what creates MapG_Gpgdb
    """
    ---------------------------------------------------------------------------

      Function that maps the conversion file gene IDs to the IDs in the
      genes.dat of the pgdb.

    ---------------------------------------------------------------------------
    Args:
        Conversion: Conversion object from conversion file
        vGenes_all: DataFrame containing all the information
                    (Attributes)that was searched for (described in
                    vGenesAttributes, default is {'UNIQUE-ID','PRODUCT', ...
                    'ACCESSION-1','ACCESSION-2'};) in genes.dat of the pgdb.
        vGenesDat_GeneInfo: List of Attributes of genes.dat. Defines where gene
                            information should be searched for.
        vPGDBIdsToMap: A string that contains G and or T and or P and defines
                       if the entries represent Genes and or transcript ids and
                       or proteins ids of the conversion file.

    Returns:
        A dictionary mapping between Genes in the conversion file and genes
                     stored in the pgdb genes.dat file.
    """
    _log.debug("Mapping PGDBs to conversion file")
    if pgdb_ids_to_map not in "GTP":
        _log.error("Invalid PGDBIdsToMap mapping value, valid values are G, T, and P")
        raise Exception("Invalid Mapping Value")

    gene_map = dict()
    if "G" in pgdb_ids_to_map:
        all_genes_df_genes = all_genes_df.select(
            pl.col("UNIQUE-ID"), pl.col(genes_data_metadata[0])
        ).explode(pl.col(genes_data_metadata[0]))
        gene_map = (
            Conversion.df.join(
                all_genes_df_genes, left_on="Gene Name", right_on=genes_data_metadata[0]
            )
            .select(pl.col("Gene Name"), pl.col("UNIQUE-ID"))
            .rows_by_key(key="UNIQUE-ID")
        )

    transcript_map = dict()
    if "T" in pgdb_ids_to_map:
        all_genes_df_transcript = all_genes_df.select(
            pl.col("UNIQUE-ID"), pl.col(genes_data_metadata[0])
        ).explode(pl.col(genes_data_metadata[0]))
        transcript_map = (
            Conversion.df.join(
                all_genes_df_transcript,
                left_on="Transcript Name",
                right_on=genes_data_metadata[0],
            )
            .select(pl.col("Transcript Name"), pl.col("UNIQUE-ID"))
            .rows_by_key(key="UNIQUE-ID")
        )

    protein_map = dict()
    if "P" in pgdb_ids_to_map:
        all_genes_df_proteins = all_genes_df.select(
            pl.col("UNIQUE-ID"), pl.col(genes_data_metadata[0])
        ).explode(pl.col(genes_data_metadata[0]))
        protein_map = (
            Conversion.df.join(
                all_genes_df_proteins,
                left_on="Protein Name",
                right_on=genes_data_metadata[0],
            )
            .select(pl.col("Protein Name"), pl.col("UNIQUE-ID"))
            .rows_by_key(key="UNIQUE-ID")
        )

    ret_keys = (
        set(gene_map.keys()).union(protein_map.keys()).union(transcript_map.keys())
    )
    # TODO: is there a way to just pass dataframes to result? (issue #31)
    return {
        x: list(gene_map.get(x, []))
        + list(transcript_map.get(x, []))
        + list(protein_map.get(x, []))
        for x in ret_keys
    }


def perform_clustering(
    vGeneLocation: GeneLoc,
    mcl_cluster_data: MclClusterData,  # covered by genelocation
    vEnzMinForClusters: int,
    vKbRangeLD: int,
    vGeneRangeLD: int,
    vCriterion: int,  # TODO: this should be an enum (issue #32)
    vMinStepSize: int = 5,
    vMaxStepSize: int = 20,
) -> Clusters:
    """
    ---------------------------------------------------------------------------

      Function that clusters the genome. It lookes at enzymes that are at max
      spaced by vStepsize genes (checking for 0, 1, ..., vMinStepSize intervening
      genes. It continuous to increase the vStepsize until clusters contain
      more non-enzymes than enzymes. Clusters containing only one group of
      locally duplicated genes are turned down. Clusters with only one reaction
      are turned down as well. Clusters need at least vEnzMinForClusters enzymes.

    ---------------------------------------------------------------------------

    Args:
        vGeneLocation: DataFrame about the gene location information.
        mcl_cluster_data: Mapping between the genes in the conversion-ID file and the
                MCL-clustering file.
        vEnzMinForClusters: Defines the minimum number of enyzmes that need to
                            be present in a cluster.
        vKbRangeLD: Genes that are encoded more than this bp apart, are not
                    called local duplicates.
        vGeneRangeLD: Number of genes that space two genes from each other.
                      Genes that have more than this many intervening genes
                      between them are not called local duplicates.
        vCriterion: Criterion after which step size is chosen
                   - 1: as published in Schlapfer et al 2017
                       sum(enzymes in clusters) - sum(non-enzymes in clusters)
                       > 0.
                   - 2: stepsize 3 (because most organisms take this one)
                   - 3: sum(enzymes in clusters) - sum(non-enzymes in clusters,
                       without hypo genes) -2*sum(hypogenes)
                   - 4: sum(clusters in stepsize n) > sum(clusters in stepsize
                       n+1).
        vMinStepSize: Minimal intervening gene size that should be computed
                      (default is 5).
        vMaxStepSize: Maximal intervening gene size that should be computed
                      (default is 20).

    Returns:
        Cluster object which contains the following information:
        - **Clusters**: Cluster Information (step_size, chromosome, start,
                   end, enzyme_count, non_enzyme_gene_count, hypothetical_gene_count)
        - **ClusterBackground**: Density information about the genome (window size,
                            number of enzymes, Occurrence of windowsize with
                            this many enzymes in the genome, occurence of
                            this windowsize in the genome, Frequency of seeing
                            a windowsize of this size with at least as many
                            enzymes in the genome)
        - **FinalStepsize**: Stepsize that was chosen by the criterion (clusters
                        contain in average more enzymes than non-enzymes).
        - **StepSize_Max**: Maximal stepsize that was computed.
    """
    _log.debug("Performing clustering")
    final_step_size = 0
    continue_search = True  # flag if search should continue
    clusters = []  # array of cluster identifiers
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
        # pl.DataFrame()
    )  # Stepsize, start, end, number of enzymes, number of non-enzyme genes, number of hypos
    clusters_background = (
        pl.DataFrame()
    )  # Size, number of enzymes, count, number of total
    step_size = 1  # step size for search
    # note: hypothetical gene names start with "Hypo", this was originally tracked in a variable

    # used in place of biomG_E
    is_enzyme_criteria = pl.col("has_reaction")

    # get rid of clusters that only contain Local Duplicates

    def ld_is_not_a_problem(cluster_genes: pl.DataFrame, start: int, end: int) -> bool:
        """
        Test if genes between start and end index are all local duplicates.
        Local duplicates must be within a provided range (vKbRangeLD) and
        number of genes (GeneRangeLD).
        """
        if end is None:
            return False
        if cluster_genes.height == 0:
            return False
        # LD is not a problem if genes are over vKBRangeLD apart
        bp_range = cluster_genes.select(
            pl.col("Gene Start (bp)").max().alias("last_gene_start"),
            pl.col("Gene End (bp)").min().alias("first_gene_end"),
        )
        if bp_range["last_gene_start"][0] - bp_range["first_gene_end"][0] > vKbRangeLD:
            return True
        # LD is not a problem if "they are encoded within gene range"
        # line 3387 from matlab - count of non-hypothetical genes
        if (
            cluster_genes.filter(pl.col("Hypothetical").not_()).height
            >= vGeneRangeLD + 2
        ):
            return True
        # LD is not a problem if there are two local duplication groups
        # note: this test is only triggered for clusters of less than vKbRangeLD and vGeneRangeLD size
        # TODO: Validate this column, logic (issue #33)

        mcl_mapper = pl.DataFrame(
            {
                "Gene Name": list(mcl_cluster_data.Gene2MCLClusters.keys()),
                "mcl_clusters": list(mcl_cluster_data.Gene2MCLClusters.values()),
                "mcl_clusters_str": str(
                    list(mcl_cluster_data.Gene2MCLClusters.values())
                ),
            }
        )
        cluster_genes = cluster_genes.join(mcl_mapper, on="Gene Name", how="left")
        unique_mcl_clusters: int = cluster_genes.n_unique("mcl_clusters")
        if unique_mcl_clusters > 1:
            return True
        return False

    gene_count = len(vGeneLocation)
    # is there a way to make this search an expression?
    while step_size <= vMinStepSize or continue_search:  # for each valid step size
        cluster_starts_df = (
            vGeneLocation.data.with_row_index("index")
            .select(
                pl.col("index").alias("start"),
                pl.col("Chromosome Name"),
                pl.col("has_reaction"),
                enzyme_count=pl.col("has_reaction")
                .sum()  # sum counts True as 1, False as 0
                .rolling(index_column="index", period=f"{step_size}i"),
            )
            .filter(
                pl.col("has_reaction")  # start is enzyme
                # & pl.col("enzyme_count") >= 2  # another enzyme in search range
            )
        )
        cluster_starts_count = len(cluster_starts_df)
        _log.debug(
            f"Found {cluster_starts_count} cluster starts for step size {step_size}"
        )
        for i, gene in enumerate(cluster_starts_df.iter_rows()):
            # gene has start index and start chromosome
            start, chrom_start, *_ = gene
            _log.debug(
                "Potential Gene Start %i (index: %i) of %i starts, %i genes (step size %i)",
                i,
                start,
                gene_count,
                cluster_starts_count,
                step_size,
            )
            end: int = start
            if start is None:
                break  # no more enzymes
            # all genes in chrom_start chromosome
            chromosome_criteria = pl.col("Chromosome Name") == chrom_start

            # find the end index
            # - search from end to end+step_size
            # - IF there are enzymes in that range on chromosome
            # - THEN that's the new end and continue the search
            # - ELSE keep established end
            # FIXME: possibly a way to collapse search?
            stop_search_end = False
            while not stop_search_end:
                if end is None:
                    _log.error("type of end is none")
                    raise Exception("Cluster search end was type None")
                end_s = end  # end search start
                # stop search step size
                end_e = min(end + step_size, len(vGeneLocation.data))
                _log.debug(
                    f"Searching for enzyme cluster for index {start} from {end_s} to {
                        end_e
                    }"
                )
                end_df = vGeneLocation.data.with_row_index("index").filter(
                    is_enzyme_criteria
                    & chromosome_criteria
                    & (pl.col("index") > end_s)
                    & (pl.col("index") <= end_e)
                )
                if end_df.height > 0:
                    end = end_df.select(pl.last("index"))["index"][0]
                else:
                    stop_search_end = True

            # .filter( (pl.col("index") > start) & (pl.col("index") <= end))
            cluster_genes = vGeneLocation.data.with_row_index("index")[start:end]
            # get rid of clusters that only have one reaction set
            if len(cluster_genes.filter(is_enzyme_criteria)) < vEnzMinForClusters:
                # abolish cluster and continue search
                continue
            if (
                cluster_genes.filter(pl.col("is_locally_duplicated"))
                .select(pl.col("mcl_clusters"))
                .unique()
                .height
                <= 1
            ):
                # if not ld_is_not_a_problem(cluster_genes, start, end):
                # get rid of clusters that only have one local duplicate
                # criteria for this is outlined in the ld_is_not_a_problem function
                # abolish cluster and continue search
                continue

            # This is the code that writes the cluster into the clusters dataframe
            # what goes in here is described in line 3426 in matlab, but it is flexible
            # based on how you want to use it

            _cluster = cluster_genes
            _chromosome = _cluster["Chromosome Name"][0]

            clusters.append(
                cluster_tuple(
                    _chromosome,
                    step_size,
                    start,
                    end,
                    _cluster.filter(pl.col("has_reaction")).height,  # enzyme_count
                    _cluster.filter(
                        pl.col("has_reaction").not_() & pl.col("Hypothetical").not_()
                    ).height,  # non_enzyme_gene_count
                    _cluster.filter(
                        pl.col("Hypothetical").not_()
                    ).height,  # hypothetical_gene_count
                )
            )

        if len(clusters) > 0:
            clusters_df = pl.DataFrame(clusters)
            _log.debug(f"Head of clusters_df\n{clusters_df.head()}")
            # re-calculate clusters_background to test against criterion
            clusters_background = compute_background(
                clusters_df, vGeneLocation, clusters_background
            )
            match vCriterion:
                case 1:
                    # (Schlapfer et al. 2017 criteria)
                    # sum (enzymes in clusters) - sum(non-enzymes in clusters (incl. hypo genes?)) > 0
                    # TODO: validate criterion (issue #34)
                    enzymes_in_step_size = (
                        clusters_df.filter(pl.col("step_size") == step_size)
                        .select(pl.col("enzyme_count"))
                        .sum()
                    )
                    _log.debug(
                        f"Number of enzymes - number of non-enzymes in clustered genes: {
                            enzymes_in_step_size
                        }",
                    )
                    if step_size == vMaxStepSize:
                        continue_search = False
                        final_step_size = step_size
                        _log.warning(
                            "Maximal step size reached: %i. Consider increasing.",
                            vMaxStepSize,
                        )
                    elif (
                        clusters_df.filter(pl.col("step_size") == step_size).height == 0
                    ):
                        _log.warning("No clusters found for step size %i.", step_size)
                    elif (
                        clusters_df.filter(pl.col("step_size") == step_size).select(
                            (
                                pl.col("enzyme_count").sum()
                                / pl.col("non_enzyme_gene_count").sum()
                            ).alias("enz_over_genes")
                        )["enz_over_genes"][0]
                        < 0.5
                    ):
                        if continue_search:
                            final_step_size = step_size - 1
                        continue_search = False
                case 2:
                    # stepsize 3 (because most organisms take this one)
                    _log.debug(
                        "Number of enzymes - number of non-enzymes (no hypos) - 2* number if hypos in clustered genes: %i - %i - %i",
                        clusters_df["enzyme_count"].count(),
                        -clusters_df["non_enzyme_gene_count"].count(),
                        -2 * clusters_df["hypothetical_gene_count"].count(),
                    )
                    if step_size == vMinStepSize:
                        final_step_size = 3
                        continue_search = False
                case 3:
                    # sum(enzymes in clusters) - sum(non-enzymes in clusters without hypo genes) - 2*sum(hypo genes)
                    if step_size == vMaxStepSize:
                        continue_search = False
                        final_step_size = step_size
                        _log.warning(
                            "Maximal step size reached: %i. Consider increasing.",
                            vMaxStepSize,
                        )
                    elif (
                        clusters_df.filter(pl.col("step_size") == step_size).height == 0
                    ):
                        _log.warning("No clusters found for step size %i.", step_size)
                    elif (
                        clusters_df.filter(pl.col("step_size") == step_size).select(
                            (
                                pl.col("enzyme_count").sum()
                                - (
                                    pl.col("non_enzyme_gene_count").sum()
                                    - pl.col("hypothetical_gene_count").sum()
                                    + 1
                                )
                                - pl.col("enzyme_count").sum()
                                - 2 * pl.col("hypoothetical_gene_count").sum()
                            ).alias("enz_less_genes")
                        )["enz_less_genes"][0]
                        < 0
                    ):
                        if continue_search:
                            final_step_size = step_size - 1
                        continue_search = False
                case 4:
                    # sum(count of clusters in stepsize n) > sum(count of clusters in stepsize n + 1)
                    if step_size != 0:
                        _log.debug(
                            "Number of clusters more than last stepsize %i",
                            clusters_df.filter(
                                pl.col("step_size") > (step_size - 1)
                            ).height,
                        )
                    else:
                        _log.debug(
                            "Number of clusters Stopsize 0: %i",
                            clusters_df.filter(pl.col("step_size") == 0).height,
                        )
                    if step_size == vMaxStepSize:
                        continue_search = False
                        final_step_size = step_size
                        _log.warning(
                            "Maximal step size reached: %i. Consider increasing.",
                            vMaxStepSize,
                        )
                    elif (
                        clusters_df.filter(pl.col("step_size") == step_size).height == 0
                    ):
                        _log.warning("No clusters found for step size %i.", step_size)
                    elif (
                        clusters_df.filter(pl.col("step_size") == step_size - 1).height
                        > clusters_df.filter(pl.col("step_size") == step_size).height
                    ):
                        if continue_search:
                            final_step_size = step_size - 1
                        continue_search = False

                case _:
                    # this should be impossible in normal operation
                    raise Exception(
                        "Invalid search criterion, valid options are 1,2,3 or 4."
                    )
        else:
            final_step_size = step_size - 1
            continue_search = False
        step_size += 1
    if final_step_size == 0:
        final_step_size = step_size - 1
    return Clusters(
        ClusterMeta=pl.DataFrame(clusters),
        ClusterBackground=clusters_background,
        FinalStepSize=final_step_size,
        MaxStepSize=step_size - 1,
    )


def compute_background(
    clusters: pl.DataFrame,
    gene_location: GeneLoc,
    clusters_background: pl.DataFrame,
    # vC_IDs: List[str],
) -> pl.DataFrame:
    """
    ---------------------------------------------------------------------------

      Function that looks at all clusters that were found by the function
      f_perform_clustering(), and computes the background information for the
      window sizes of these clusters (e.g. a cluster has 4 enzymes and two non
      enzymes, then its size is 6. The function goes with a window size 6
      through the genome and checks how many windows exist with 0, 1, 2, ...,
      5, 6 enzymes. It then calculates the probability to see a cluster with
      more or equal than 0 enzymes, more or equal than 1, ... more or equal
      than 6 enzymes. This is then used to define top n% clusters, e.g.
      clusters that are at least as dense as the top 5% of window sizes are
      marked to be of higher interest.

    ---------------------------------------------------------------------------

    Args:
        clusters: Cluster Information (step_size, chromosome id, gene start,
                gene end, number of enzymes, number of non-enzyme genes,
                number of hypothetical genes)
        gene_location: DataFrame about the gene location information,
                generated in PlantClusterFinder function.
        clusters_background: DataFrame of backgrounds for up to vClusters
                max_step_size - 1. (window_size, enzyme_count,
                                    window_occurrence, winzyme_occurrence, winzyme_frequency)
                "Winzyme" is a window with enzyme_count enzymes in it.

    Returns:
        - **vClusterBackground**: Density information about the genome (window size,
                            number of enzymes, Occurrence of windowsize with
                            this many enzymes in the genome, occurrence of
                            this windowsize in the genome, Frequency of seeing
                            a windowsize of this size with at least as many
                            enzymes in the genome)
    """
    # (for a given) window size (number of total genes)
    # (for a given) number of enzymes
    # occurrence of window size in genome
    # frequency of seeing window size with _at least_ this many genes
    _log.debug("Computing background of clusters")

    index_offset = len(gene_location)
    # create chromosome offsets in index to avoid windows across chromosome boundaries
    chromosome_offsets = pl.DataFrame(
        [
            {
                "Chromosome Name": x["Chromosome Name"],
                "index_offset": x["chr_idx"] * index_offset,
            }
            for x in gene_location.data.select(pl.col("Chromosome Name").unique())
            .with_row_index("chr_idx")
            .to_dicts()
        ]
    )
    indexed_gene_location = (
        gene_location.data.with_row_index("idx")
        .join(chromosome_offsets, on="Chromosome Name")
        .with_columns(idx=pl.col("idx") + pl.col("index_offset"))
        .sort(pl.col("idx"))
    )
    # "winzyme" is window x enzyme
    max_window = clusters.select(
        pl.sum_horizontal("enzyme_count", "non_enzyme_gene_count")
        .max()
        .alias("window_size_max")
    )["window_size_max"][0]
    _log.debug(f"Max Window {max_window}")
    # if clusters_background.is_empty():
    # re-calculating background for all windows
    clusters_background = pl.DataFrame(
        [],
        schema=[
            ("window_size", pl.Int64),
            ("enzyme_count", pl.Int64),
            ("winzyme_occurrence", pl.Int32),
            ("window_occurrence", pl.Int32),
            ("winzyme_frequency", pl.Float64),
        ],
    )
    prev_max_window_size = 1
    # else:
    #    prev_max_window_size = clusters_background.select(pl.max("window_size"))[
    #        "window_size"
    #    ][0]
    for window_size in range(prev_max_window_size, max_window + 1):
        background_df_win = (
            indexed_gene_location.select(
                # The first (window_size-1) rows will not be window_size windows
                window_size=pl.len()
                .rolling(index_column="idx", period=f"{window_size}i")
                .cast(pl.Int64),
                enzyme_count=pl.sum("has_reaction")
                .rolling(index_column="idx", period=f"{window_size}i")
                .cast(pl.Int64),
                window_occurrence=1,
            )
            .filter(pl.col("window_size") == window_size)
            .group_by("window_size", "enzyme_count")
            .agg(winzyme_occurrence=pl.sum("window_occurrence"))
            .sort(pl.col("enzyme_count"))
            .with_columns(
                # get the number of winzymes with >= enzyme count
                cumulative_winzyme=pl.col("winzyme_occurrence").cum_sum(reverse=True)
            )
        )
        # windows of size "window_size"
        total_window_occurrence = background_df_win.select(
            pl.sum("winzyme_occurrence").alias("total_window_occurrence")
        )["total_window_occurrence"][0]
        clusters_background.extend(
            background_df_win.with_columns(
                window_occurrence=total_window_occurrence,
                # frequency of windows with enzymes >= enzyme_count
                winzyme_frequency=pl.col("cumulative_winzyme")
                / total_window_occurrence,
            ).select(
                pl.col("window_size"),
                pl.col("enzyme_count"),
                pl.col("winzyme_occurrence"),
                pl.col("window_occurrence"),
                pl.col("winzyme_frequency"),
            )
        )

    return clusters_background


def write_gene_output(
    gene_output_file: Path,
    gene_location: GeneLoc,
    clusters: Clusters,
    top_percent_clusters: int,
    metabolic_domains_to_annotate: List[str],
    Gene2Rxn: GenesToReactions,
    vMapG_Gpgdb: Dict[Any, Any],
    vBiomG_E: Dict[Any, Any],
    clusters_top: List[bool],
    annotated_metabolic_domains: AnnotatedMetabolicDomains,
    vRpgdb_Rmd: Dict[Any, Any],
    reaction_labels: ReactionLabels,
    output_type: str,
) -> None:
    """
    ---------------------------------------------------------------------------

      Function that writes the genefile of PlantClusterFinder: TSV format,
      first it includes the information that was stored in the gene location
      file that was provided to Plantclusterfinder (Gene IDs, Gene Start, Gene
      end, Chromosome information, strand). Then it attaches Enzyme (1 or 0) EC
      information and Reaction information. Then whether the Gene is considered
      locally duplicated and if so, which LD-clusters it's proteins were
      associated with. Then per stepsize that was investigated (intervening
      gene size) the ClusterID where the gene is incorporated and the
      informatino whether the cluster was found to be in the top n% of
      clsuters. This is followed by the information about the metabolic
      domains, and the signature / tailoring information.

    ---------------------------------------------------------------------------

    Args:
        gene_output_file: Path and file name to the gene-Outputfile, make sure
                         you have access to the file and folder. This is the
                         first Results output file that will be generated.
        gene_location: Gene location information, generated in
                        PlantClusterFinder function.
        clusters: Cluster Information (stepsize, chromosome id, gene start, gene end, number of enzymes)
        top_percent_clusters: Clusters are labeled as being in the top n% if
                             they are within this top n% of enzyme dense
                             clusters compared to the background (default is
                             5).
        metabolic_domains_to_annotate: A list of Metabolic domains that should
                         be analyzed. For example {'Amines and Polyamines Metabolism'; ...
                         'Amino Acids Metabolism'
        Gene2Rxn: Mapping of all genes in the Genes_all to the Reactions in
        Reactions_all
        vMapG_Gpgdb: A mapping between Genes in the conversion file and genes
                     stored in the pgdb genes.dat file.
        vBiomG_E: A vector describing which gene has at least one reaction
                  annotated in the pgdb (1, if no reaction 0)
        local_duplication: DataFrame of length of the number of Genes in gene_location,
               being 1 if the gene was considered to be locally duplicated (a
               gene was within vKbRangeLD bp distance and not spaced more than
               with vGeneRangeLD genes from another gene with the same MCL
               clustering given by vG_MCL.
        clusters_Top: List with size Number of Clusters. It is true if
                       the density of the cluster is above or the same as given
                       by the input argument vTopPercentClusters*100.
        annotated_metabolic_domains: Mapping between all Reactions stored in
               the Metabolic domain file and the metabolic domains that were
               identified to be of interest in vMD_to_annotate.
        vRpgdb_Rmd: A mapping between reactions stored in the pgdb and the
                    reactions given in the metabolic domain file
        ReactionLabels: Is a data structure of signature and tailoring enzyme
                    classes.
                    1. Mapping of Reactions in vReactions_all to the enzyme
                      classes of vSTClasses.
                    2. Is a list of types of enzymes, e.g. signature, tailoring
                      and so on.
                    3. Mapping of Reactions in vReactions_all to the enzyme
                      types of vSigOrTail.
        output_type: Defines which type of output is given:
                        "simple" (default setting) gives only top n% for the
                           final stepsize.
                        "verbose" gives all computed stepsize until the final
                           stepsize.
                        "old" gives at least 5 until the final stepsize.

    Returns:
        **None** (File is generated with the path of the first input argument)
    """
    _log.debug("Writing gene output")
    if output_type not in ["old", "simple", "verbose"]:
        raise Exception(
            "Output type not supported, options are 'old' 'simple' or 'verbose'"
        )

    gene_columns = gene_location.data.columns
    _log.debug(f"Columns for gene location data are {gene_columns}")
    _log.debug(f"Head for gene location data is {gene_location.data.head()}")
    # FIXME: This can/should be one select() call (issue #36)
    gene_output = gene_location.data
    # the first several columns are straight from gene_location
    gene_output_cols = [
        pl.col("Gene Name"),
        pl.col("Gene Start (bp)"),
        pl.col("Gene End (bp)"),
        pl.col("Chromosome Name"),
        pl.col("Strand"),
        pl.col("has_reaction").alias("Enzyme"),
        pl.col("EC-NUMBER").list.join(",").alias("EC-NUMBER"),
        pl.col("pgdb_reaction_id").list.join(",").alias("RXN"),
        pl.col("is_locally_duplicated").alias("Locally Duplicated"),
        pl.col("mcl_clusters")
        .list.eval(pl.element().cast(pl.String))  # noqa S307
        # convert mcl_clusters from List[String] to comma-separated list
        .list.join(",")
        .alias("LD_ClusterIds"),
    ]
    # Join Cluster data where cluster_start <= gene index <= cluster_end
    cluster_data = (
        gene_location.data.with_row_index("gene_idx")
        .join_where(
            clusters.ClusterMeta.with_row_index("cluster_idx"),
            (pl.col("start").cast(pl.UInt32) <= pl.col("gene_idx"))
            & (pl.col("gene_idx") <= pl.col("end").cast(pl.UInt32)),
        )
        .group_by(pl.col("Gene Name"), pl.col("step_size"))
    )
    _log.debug(f"Cluster data for gene output {cluster_data.head()}")
    # different output columns depending on output type
    match output_type:
        case "old":
            for i in range(clusters.MaxStepSize):
                # chosen stepsize x
                # top pct clusters
                gene_output = gene_output.join(
                    cluster_data.agg(
                        pl.col("cluster_idx")
                        .implode()
                        .alias(f"cluster_idx_stepsize{i}"),
                        pl.col("top_percent_clusters").sum(),
                    )
                    .filter(pl.col("step_size") == pl.lit(i))
                    .select(
                        pl.col("Gene Name"),
                        pl.col(f"cluster_idx_stepsize{i}"),
                        pl.col("top_percent_clusters").alias(
                            f"top_percent_clusters{i}"
                        ),
                    ),
                    "Gene Name",
                    how="left",
                )
                gene_output_cols.append(
                    pl.col(f"cluster_idx_stepsize{i}")
                    .list.len()
                    .alias(
                        f"{'CHOSEN ' if i == clusters.FinalStepSize else ''}Stepsize {
                            i
                        }"
                    )
                )
                gene_output_cols.append(
                    pl.col(f"top_percent_clusters{i}").alias(
                        f"Top {top_percent_clusters}% {i}"
                    )
                )
        # FIXME: cluster IDs filter columns (vClusters 2 and 3)?
        case "verbose":
            for i in range(clusters.MaxStepSize):
                # chosen stepsize x
                # top pct clusters
                gene_output = gene_output.join(
                    cluster_data.agg(
                        pl.col("cluster_idx")
                        .implode()
                        .alias(f"cluster_idx_stepsize{i}"),
                        pl.col("top_percent_clusters")
                        .sum()
                        .alias(f"top_percent_clusters{i}"),
                    )
                    .filter(pl.col("step_size") == pl.lit(i))
                    .select(
                        pl.col("Gene Name"),
                        pl.col(f"cluster_idx_stepsize{i}"),
                        pl.col(f"top_percent_clusters{i}"),
                    ),
                    "Gene Name",
                    how="left",
                )
                # FIXME: There is no reason I can see to have two loops
                for j in range(clusters.MaxStepSize):
                    gene_output_cols.append(
                        pl.col(f"cluster_idx_stepsize{i}")
                        .list.len()
                        .alias(
                            f"{
                                'CHOSEN ' if i == clusters.FinalStepSize else ''
                            }Stepsize {i},{j}"
                        )
                    )
                    gene_output_cols.append(
                        pl.col(f"top_percent_clusters{i}").alias(
                            f"Top {top_percent_clusters}% {i},{j}"
                        )
                    )
        case "simple":
            # stepsize x
            # top pct clusters
            top_pct_cluster_idx = cluster_data.agg(
                pl.col("cluster_idx")
                .count()
                .alias(f"Stepsize {clusters.FinalStepSize}"),
                pl.col("top_percent_clusters")
                .sum()
                .alias(f"Top {top_percent_clusters}%"),
            ).filter(pl.col("step_size") == pl.lit(clusters.FinalStepSize))
            gene_output = gene_output.join(top_pct_cluster_idx, "Gene Name", how="left")

    # describe metabolic domains to annotate for run
    for md in metabolic_domains_to_annotate:
        gene_output_cols.append(
            pl.col("pgdb_reaction_id")  # for all reaction IDs in gene
            .list.eval(  # noqa s307
                pl.element()
                # for all metabolic domains in gene
                .replace_strict(vRpgdb_Rmd, default=[])
                .list.contains(md)  # is md a metabolic domain?
                .any()  # does md appear once?
            )
            .list.first()
            .fill_null(False)  # convert list[bool] of size 1 to bool
            .alias(md)
        )

    # add the enzyme classes and  if they are signature/tailoring
    gene_output_cols.extend(
        (
            pl.col("pgdb_reaction_id")
            .list.eval(  # noqa s307 -- these "eval" statements don't evaluate code in them
                # get reaction class label for all reactions in gene
                pl.element().replace_strict(reaction_labels.Rxn_STClasses, default=None)
            )
            .list.join(",")  # convert from list[String] to String (for output)
            .alias("Signature Tailoring Enzyme Class"),
            pl.col("pgdb_reaction_id")
            .list.eval(  # noqa s307
                pl.element().replace_strict(reaction_labels.Rxn_SigOrTail, default=None)
            )
            .list.join(",")
            .alias("Signature or Tailoring"),
        )
    )
    # debug: see explain plan for each column
    for i, column in enumerate(gene_output_cols):
        _log.debug(f"Column {i} {column.meta.output_name()}")
        try:
            _log.debug(f" explain: {gene_output.select(column).schema}")
        except Exception as e:
            _log.warning(f" column {i} failed to explain {e}")
    gene_output = gene_output.select(gene_output_cols)
    _log.info(f"Gene Output is {len(gene_output)} rows")
    _log.debug(f"Gene Output head {gene_output.head()}")
    gene_output.write_csv(gene_output_file, separator="\t")


def write_cluster_output(
    cluster_output_file: Path,
    gene_location: GeneLoc,
    clusters: Clusters,
    clusters_in_top_percent: List[bool],
    output_type: str,
) -> None:
    """
    ---------------------------------------------------------------------------

      Function that writes the clusterfile of PlantClusterFinder: TSV format,
      first it includes the stepsize of the cluster (intervening genes), then
      cluster ID, chromosome/scaffold the cluster is located on, start bp, end
      bp, physical size, start gene ID, end gene ID, and all the gene IDs

    ---------------------------------------------------------------------------

    Args:
        cluster_output_file: Path and file name to the cluster-Outputfile, make
                            sure you have access to the file and folder. This
                            is the second Results output file that will be
                            generated.
        gene_location: DataFrame about the gene location information,
                        generated in PlantClusterFinder function.
        clusters: Cluster Information (stepsize, chromosome id, gene start,
                   gene end, number of enzymes)
        clusters_in_top_percent: List with size Number of Clusters x 1. It is zero if
                       the density of the cluster is above or the same as given
                       by the input argument vTopPercentClusters*100.
        output_type: Defines which type of output is given: simple (default
                     setting) gives only top n% for the final stepsize. verbose
                     gives all computed stepsize until the final stepsize. old
                     gives at least 5 until the final stepsize.

    Returns:
        **None** (File is generated with the path of the first input argument)
    """
    _log.debug("Writing cluster output")
    if output_type not in ["old", "simple", "verbose"]:
        raise Exception("Output type not supported, options are 'simple' or 'verbose'")
    indexed_clusters = clusters.ClusterMeta.with_row_index("cluster_id")
    cluster_genes = (
        gene_location.data.with_row_index("gene_idx")
        .join_where(
            indexed_clusters,
            (pl.col("Chromosome Name") == pl.col("chromosome")),
            (pl.col("gene_idx") >= pl.col("start")),
            (pl.col("gene_idx") <= pl.col("end")),
        )
        .group_by("cluster_id")
        .agg(
            pl.col("Gene Name").alias("cluster_genes"),
            pl.col("Gene Start (bp)").min().alias("start_bp"),
            pl.col("Gene End (bp)").min().alias("end_bp"),
        )
    )
    test_data = indexed_clusters.select(
        pl.col("cluster_id"), pl.col("chromosome"), pl.col("start"), pl.col("end")
    )
    _log.debug(f"Clusters Data ID, chromosome, and start and end {test_data}")
    _log.debug(f"Cluster Genes Only {cluster_genes.head()}")
    _log.debug(f"Cluster Genes Rows {len(cluster_genes)}")
    _log.debug(f"Cluster Genes Columns {cluster_genes.columns}")
    cluster_output = indexed_clusters.join(cluster_genes, on="cluster_id")
    _log.debug(f"Cluster output length {len(cluster_output)}")
    _log.debug(f"Cluster Step Sizes {cluster_output['step_size'].describe()}")
    cluster_output = (
        cluster_output.with_columns(
            cluster_length=pl.col("end_bp") - pl.col("start_bp") + pl.lit(1),
            total_gene_count=pl.col("enzyme_count") + pl.col("non_enzyme_gene_count"),
        ).select(
            pl.col("step_size").alias("Stepsize"),
            pl.col("cluster_id").alias("ClusterID"),
            pl.col("top_percent_clusters").alias("Top n% cluster"),
            pl.col("chromosome").alias("Chromosome"),
            pl.col("start_bp").alias("Start bp"),
            pl.col("end_bp").alias("End bp"),
            pl.col("cluster_length").alias("Physical size bp"),
            pl.col("total_gene_count").alias("Number of genes"),
            pl.col("hypothetical_gene_count").alias("Number of Hypothetical genes"),
            pl.col("enzyme_count").alias(
                "number of Enzymes"
            ),  # TODO: capitalize number (changed to match MATLAB)
            pl.col("cluster_genes").list.first().alias("Start Gene"),
            pl.col("cluster_genes").list.last().alias("End Gene"),
            pl.col("cluster_genes").list.join(",").alias("Gene IDs"),
        )
    ).filter(
        (pl.lit(output_type) == pl.lit("old"))
        | (
            (pl.lit(output_type) == pl.lit("verbose"))
            & (pl.col("Stepsize") <= pl.lit(clusters.FinalStepSize))
        )
        | (
            (pl.lit(output_type) == pl.lit("simple"))
            & (pl.col("Stepsize") == pl.lit(clusters.FinalStepSize))
            & pl.col("Top n% cluster")
        )
    )
    _log.debug(f"Cluster output dtypes {cluster_output.dtypes}")
    _log.debug(f"Cluster output {cluster_output}")
    cluster_output.write_csv(cluster_output_file, separator="\t")
    return  # early return to keep the old code


def analyze_plant_cluster_gap_file(
    vFileIn: Path, vUnmaskedDNA: bool, file_output: bool = True
) -> Optional[Tuple[str, Dict[Any, list[Any]]]]:
    """
    ---------------------------------------------------------------------------

      Function that analyzes the masked DNA sequences of a _GAPOutput file, and
      counts the occurence of every size of gap.

    ---------------------------------------------------------------------------

    Args:
        vFileIn: Full path to the _GAPOutput file that whould be analyzed
                 (produced within PlantClusterFinder)
        vUnmaskedDNA: Option to ignore the error if you use unmasked DNA files.
        file_output: Whether output is directed to file or to

    Returns:
        **Optional[str]** (File with list of sequencing gaps numbers is generated with
        ending _count.txt if file_output is set to True)
    """
    if not vFileIn.is_file():
        raise Exception(f"vFileIn {vFileIn} is not a file or does not exist")
    _gap_output = (
        pl.read_csv(
            vFileIn,
            separator=" ",
            has_header=False,
            new_columns=["scaffold", "amino acid", "start index", "length"],
        )
        .filter(pl.col("amino acid") == "N")
        .group_by("length")
        .len()
        .rows_by_key(key="length")
    )
    _f_out = vFileIn.with_suffix(vFileIn.suffix + "_count.txt")
    if not file_output:
        return (str(_f_out), _gap_output)
    with open(_f_out, "w") as f:
        f.write("Gapsize in bp\tNumber of Gaps\n")
        for k, v in sorted(_gap_output.items(), key=lambda x: -list(x)[1][0]):
            f.write(f"{k}\t{list(v)[0]}\n")
    return None


def read_in_pgdb_flat(
    pgdb_filename: Path,
    attrib: tuple[str],
    unique_keys: tuple[str] = ("UNIQUE-ID",),
    pgdb_encoding: str = "latin-1",
) -> pl.DataFrame:
    """
    ---------------------------------------------------------------------------

      Function that reads in a pgdb flat file (e.g. genes.dat, proteins.dat,
      enzrxns.dat, reactions.dat) saving Attributes given by vAttrib.

    ---------------------------------------------------------------------------

    Args:
        pgdb_filename: PGDB Flat file that should be read in.
        attrib: Attributes of the PGDB Flat file that should be saved (all
                 others are discarded.
    Optional Args:
        **unique_keys**: tuple of columns that are unique, i.e. should not be stored as a list

    Returns:
        **DataFrame** structure with the contents of each attribute as a column
    """

    # By encapsulating the parsing logic in a generator,
    #  it opens up the door to potentially deferring reading
    #  or simplifying asynchronous reading
    def _row_gen(
        fh: IO, attrib: Tuple[str], unique_keys: Tuple[str]
    ) -> Iterator[Dict[str, Any]]:
        """
        Generate a dict row from the pgdb flat file
        """
        y: Dict[str, Any] = {att: None for att in attrib}
        prev_key = "COMMENT"
        for line in fh:
            # ignore comment rows, columns
            if line.startswith("#") or line.startswith("COMMENT"):
                if line.startswith("COMMENT - "):
                    prev_key = "COMMENT"
                continue
            # line continues or data row ends
            elif line.startswith("/"):
                # split lines and reset row
                if line.startswith("//"):
                    yield y
                    y = {att: None for att in attrib}
                    prev_key = "COMMENT"
                    continue
                # row ends, extend value from previous row
                elif prev_key == "COMMENT":
                    continue
                else:
                    if not isinstance(y[prev_key], list):
                        y[prev_key] += line[1:]
                        continue
                    last_val = y[prev_key][-1]
                    last_val += line[1:]
                    y[prev_key][-1] = last_val

            # get key, value pair from line
            elif " - " in line:
                key, *value_arr = line.split(" - ")
                prev_key = key
                value: str = (" - ".join(value_arr)).strip()
                if key in y:
                    if y[key] is None:
                        if key == unique_keys[0]:
                            y[key] = value
                        else:
                            y[key] = [value]
                    elif not isinstance(y[key], list):
                        y[key] = [y[key], value]
                    else:
                        y[key].append(value)

    # The PGDB source data doesn't enforce encoding, but Arabiopsis
    with open(pgdb_filename, "r", encoding=pgdb_encoding) as f:
        pgdb_data = pl.from_dicts(
            _row_gen(f, attrib, unique_keys), infer_schema_length=None
        )

    return pgdb_data
