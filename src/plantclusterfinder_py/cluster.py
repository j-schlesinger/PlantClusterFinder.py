"""Functions used in calculating where biosynthetic gene clusters are."""

import logging
import os
import shutil
import subprocess
from collections import namedtuple
from pathlib import Path
from typing import Any

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
from plantclusterfinder_py.pgdb_utils import read_in_pgdb_flat

_log = logging.getLogger(__name__)


def label_reactions_with_sig_and_tail(
    all_reactions_df: pl.DataFrame, signature_tailor_filename: Path
) -> ReactionLabels:
    """---------------------------------------------------------------------------

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
        signature_tailor_filename: TSV file describing the classes of enzymes that
                              should be classified as signature or tailoring
                              enzymes to identify clusters containing such.

    Return:
        ReactionLabels object with the signature and tailoring enzyme classes,
            list of types of enzymes, and mapping between reactions and
            enzyme classes as well as signature or tailoring enzymes

    """
    _log.debug("Reading in the signature tailor file %s", signature_tailor_filename)
    stdata = pl.read_csv(
        signature_tailor_filename,
        has_header=True,
        separator="\t",
        encoding="ISO-8859-1",
        truncate_ragged_lines=True,
    )
    _log.debug("All reactions df is %s", all_reactions_df)
    cols = stdata.columns
    _log.debug("read %d rows with columns %s", len(stdata), cols)
    if len(stdata) != len(stdata["EC_number or RXN_ID"].unique()):
        _log.warning("signature tailor file has duplicate EC numbers or RXN IDs.")
    stclasses = stdata["Enzyme class"].unique()
    sigortail = stdata["Enzyme role"].unique()
    rxn_class_dict = stdata.select("EC_number or RXN_ID", "Enzyme class").rows_by_key(
        "EC_number or RXN_ID", unique=True
    )
    rxn_sig_or_tail_dict = stdata.select("EC_number or RXN_ID", "Enzyme role").rows_by_key(
        "EC_number or RXN_ID", unique=True
    )
    # returns a tuple for each row, this extracts the value the list of Enzyme roles
    rxn_class = {k: v[0] for k, v in rxn_class_dict.items()}
    rxn_sig_or_tail = {k: v[0] for k, v in rxn_sig_or_tail_dict.items()}
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
    """---

      Function that introduces Hypothetical genes if a sequencing gap reached a
      critical size hypothetical_gene_gap_length_minimum.
    ---

    Args:
        gene_location: Object with the gene location, chromosome regions, sequence gaps,
            and minimal length of bp of masked DNA that should be populated by two
            hypothetical genes
        sequencing_gaps: DataFrame containing the location and size of sequence gaps.
        hypothetical_gene_amount: Integer describing the number of hypothetical genes
            a sequencing gap that is not a physical break should be populated
            by artificial genes.
        hypothetical_gene_gap_length_minimum: A size that defines the minimal
            length of a sequencing gap that should still be processed.
            If this is bigger than 0, then all sequencing gaps smaller than
            this are deleted, unless they are among the list of sequence_gaps,
            which also needs to be set. Default is -1, not active.

    Returns:
        **HardBreaks** object containing the gene location information, the
            inter-chromosomal matrix, and the HypoAmount populated.

    """
    _log.debug("Introducing hypothetical genes into sequence")
    # get the 5th and 7th columns of sequence_gaps
    gap_positions = (
        sequencing_gaps.select(pl.col(["Chromosome Name", "Start", "End"]))
        .with_columns(seq_gaps_length=pl.col("End") - pl.col("Start"))
        .with_columns(
            hypo_to_place=pl.when(
                (hypothetical_gene_gap_length_minimum < 0)
                | (pl.col("seq_gaps_length") >= pl.lit(hypothetical_gene_gap_length_minimum))
            )
            .then(pl.lit(hypothetical_gene_amount))
            .otherwise(pl.lit(0))
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
                # a strand is required for data purposes, but not used in any calculation
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
    _log.debug("Inserted %d hypothetical genes", len(hypo_gene_df))
    return gene_location


def annotate_sequencing_gaps(
    nucleotide_fasta_filename: Path,
    intergenic_location: pl.DataFrame,
    unique_chromosome_names: list[str],
    parallelize: int = 1,
    dna_unmasked_ok: bool = True,
    overwrite_sequence_gap_info_file: bool = False,
    sequence_gap_sizes_chromosome_breaks: list[int] | None = None,
    hypothetical_gene_gap_length_minimum: int = -1,
) -> pl.DataFrame:
    """---------------------------------------------------------------------------

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
        unique_chromosome_names: A list with unique entries of all chromosome
            and scaffold dna_unmasked_ok: Option to ignore the error if you use
            unmasked DNA files.
        parallelize: How many threads/cores can this process parallelize across.
        dna_unmasked_ok: Set to true of using unmasked DNA FASTA input.
            (Default is True).
        overwrite_sequence_gap_info_file: Set to true if you want to allow
            overwriting the sequence gap info file.
            (Default is False)
        hypothetical_gene_gap_length_minimum: If set to True, the code enforces
            masked nucleotide analysis for finding sequencing gaps
            (Default is False, no enforcement).
        sequence_gap_sizes_chromosome_breaks: A list of numbers.
            The masked nucleotides sequences of these specific lengths are preventing
            cluster formation If this is given, then
            hypothetical_gene_gap_length_minimum also has to be set. Will lead to
            prescreening of the sequencing gaps.
            (Default is an empty list []).
        hypothetical_gene_gap_length_minimum: A size that defines the minimal
            length of a sequencing gap that should still be processed.
            If this is bigger than 0, then all sequencing gaps smaller than this are
            deleted, unless they are among the list of
            sequence_gap_sizes_chromosome_breaks, which also needs to be set.
            (Default is -1, not active.)

    Returns:
        **DataFrame** containing the
                  - scaffold: chromosome/scaffold,
                  - start index: start bp
                  - end index: end bp and
                  - seq_gaps_start
                  - seq_gaps_end
                  - seq_gaps_length: length of a masked region.
                  In addition there is the information stored in with
                  intergenic region the gap is encoded and its start and end
                  bp, corrected for maximal size (= intergenic region).

    """
    if sequence_gap_sizes_chromosome_breaks is None:
        sequence_gap_sizes_chromosome_breaks = []
    if not (
        hypothetical_gene_gap_length_minimum == -1 or sequence_gap_sizes_chromosome_breaks == []
    ):
        _log.error(
            "annotate_sequencing_gaps needs both SeqGapSizesChromBreak and"
            + "GeneLengthHypoMin to be set or to be unset."
        )
        raise AttributeError("Incompatable parameters in annotate_sequencing_gaps function")
    _log.debug("Annotating sequence gaps")

    gap_output_file = nucleotide_fasta_filename.with_suffix(
        nucleotide_fasta_filename.suffix + "_GAPOutput"
    )
    _log.debug("Using gap output file %s", gap_output_file)
    if overwrite_sequence_gap_info_file or not gap_output_file.exists():
        gap_output_file = get_sequencing_gaps(
            nucleotide_fasta_filename,
            parallelize,
            overwrite_sequence_gap_info_file,
            dna_unmasked_ok,
        )
    sequencing_gaps = (
        pl.read_csv(
            gap_output_file,
            separator=" ",
            has_header=False,
            new_columns=["scaffold", "amino acid", "start index", "length"],
        )
        .with_columns((pl.col("start index") + pl.col("length") - pl.lit(1)).alias("end index"))
        .filter(pl.col("amino acid") == pl.lit("N"))  # Do not process all non N sequences.
    )
    for bad_chromosome, *_ in (
        sequencing_gaps.filter(pl.col("scaffold").is_in(unique_chromosome_names).not_())
        .select(pl.col("scaffold").unique())
        .iter_rows()
    ):
        _log.warning(
            "DNA Fasta file '%s' contains chromosomes or scaffolds not covered with the "
            + f"gene position file: %s",
            nucleotide_fasta_filename,
            bad_chromosome,
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
        .with_columns(seq_gaps_length=pl.col("seq_gaps_end") - pl.col("seq_gaps_start") + 1)
    )
    _log.debug("Found %d sequencing gaps", sequencing_gaps.height)
    _log.debug("Example gaps %s", sequencing_gaps.head())
    # find the intergenic region of the gap, if the gap is present within a
    # gene model, get rid of the gap.
    # logic validated in issue (#28)
    # only gaps which are in intergenic regions
    # & indicate which intergenic region
    return sequencing_gaps


def get_sequencing_gaps(
    sequence_fasta_filename: Path,
    parallelize: int,
    overwrite_sequence_gap_info_file: bool = False,
    dna_unmasked_ok: bool = False,
) -> Path:
    """---------------------------------------------------------------------------

      Function that screens a DNA FASTA file for masked stretches of DNA
      (masked by N).

      Creates a "..._GAPOutput" file from the FASTA file, then creates a
      "..._GAPOutput_count.txt" file by running analyze_plant_cluster_gap_file()

    ---------------------------------------------------------------------------

    Args:
        sequence_fasta_filename: FASTA file of the masked genome nucleotide sequences.
        parallelize: Defines if multiple CPUs should be used for computation of parts
               of algorithm. Active if larger than 1.
        overwrite_sequence_gap_info_file: If set to True, the code enforces masked nucleotide
                              analysis for finding sequencing gaps (default is
                              False, no enforcement).
        dna_unmasked_ok: Option to ignore the error if you use unmasked DNA files.

    Returns:
       Path of gap file
    Raises:
        IOError: if file cannot be generated.
        AttributeError: if attempting to overwrite file

    """
    sequence_gap_filename = sequence_fasta_filename.with_suffix(
        sequence_fasta_filename.suffix + "_GAPOutput"
    )
    if not overwrite_sequence_gap_info_file and sequence_gap_filename.exists():
        _log.error("Cannot overwrite existing _GAPOutput with tag set")
        raise AttributeError("Overwrite Protection Enforcement of Existing Gap File")
    _log.debug("Writing sequence gap file %s", sequence_gap_filename)
    get_positions_of_gap(
        sequence_fasta_filename,
        sequence_gap_filename,
        parallelize=parallelize,
    )
    analyze_plant_cluster_gap_file(
        sequence_gap_filename,
        dna_unmasked_ok,
    )
    return sequence_gap_filename


def get_metabolic_domains(
    genes_metabolic_domains_filename: Path, metabolic_domains_to_annotate: list[str]
) -> AnnotatedMetabolicDomains:
    """---------------------------------------------------------------------------

      Function that reads in the metabolic domains of reactions
      and gives a mapping between reactions in these files and the metabolic
      domains that were asked to be annotated (usually a subset of all) by
      metabolic_domains_to_annotate.

    ---------------------------------------------------------------------------

    Args:
        genes_metabolic_domains_filename: TSV of the metabolic domains of reactions.
                                     Check if file is outdated (needs to
                                     contain all reactions that are present in
                                     the PGDB that have genes annotated to).
        metabolic_domains_to_annotate: A list of Metabolic domains that should be
            analyzed.
            For example {'Amines and Polyamines Metabolism'; 'Amino Acids Metabolism'

    Returns:
        AnnotatedMetabolicDomains object with all items stored in Metabolic domain file
        as well as mappings for metabolic domains of interest described in
        `MD_to_annotate`.

    """
    _log.debug(f"Getting metabolic domains from file {genes_metabolic_domains_filename}")
    # file used to have a \r EOL character
    rxn_meta_dom = (
        pl.read_csv(genes_metabolic_domains_filename, separator="\t")
        .with_columns(pl.col("Reaction ID").str.strip_chars(" ").alias("Reaction ID"))
        .unpivot(
            [f"Class{n + 1}" for n in range(6)] + ["SM1A"],
            index=["Reaction ID", "Pathway"],
            value_name="Metabolic Domain",
        )
        .filter(pl.col("Metabolic Domain") != "NA")
        .filter(pl.col("Metabolic Domain").is_in(metabolic_domains_to_annotate))
        .select(pl.col("Reaction ID"), pl.col("Pathway"), pl.col("Metabolic Domain"))
    )

    return AnnotatedMetabolicDomains(
        RxnMetaDom=rxn_meta_dom, Rxn2MetaDom=rxn_meta_dom.rows_by_key(key="Reaction ID")
    )


def read_mcl_clustering_file(
    mcl_cluster_filename: Path,
    protein_sequence_fasta_filename: Path,
    conversion: Conversion,
    remove_non_protein_locations: bool,
) -> MclClusterData:
    """---------------------------------------------------------------------------

      Function that reads out the MCL clustering information and gives every
      gene one or several MCL clusters. This is then used for local
      duplication.

    ---------------------------------------------------------------------------

    Args:
        mcl_cluster_filename: MCL-clustering output file (I20).
        protein_sequence_fasta_filename: The FASTA file of the protein sequences.
        conversion: Lists of gene, and protein IDs with mapping from conversion files.
        remove_non_protein_locations: True if protein-IDs should be stored
            in return value.

    Returns:
        MclClusterData object which provides mapping between
        genes in the conversion-ID file and the MCL-clustering file,
        as well as a list of protein-IDs stored in the fasta file.

    """
    _log.debug("Reading in MCL Clustering File %s", mcl_cluster_filename)
    protein_ids = []
    with protein_sequence_fasta_filename.open() as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                protein_ids.append(line[1:].strip())

    mcl_cluster_data: list[dict[str, int | str | list[str]]] = []
    clusters_by_p: dict[str, list[str]] = {}
    clusters_by_g: dict[str, list[str]] = {}

    # Read in MCL clusters, and assign ID based on index
    with mcl_cluster_filename.open() as f:
        for i, line in enumerate(f):
            mcl_id = f"MCL_{i}"
            prots_in_cluster = line.strip().split()
            # Convert Transcript ID to Gene Name
            genes_in_cluster = [conversion.map_proteins_to_genes[p] for p in prots_in_cluster]
            if any(len(x) == 0 for x in genes_in_cluster):
                _log.error("Empty genes cluster in MCL %s", mcl_id)
                missing_genes = [
                    p for p in prots_in_cluster if len(conversion.map_proteins_to_genes[p])
                ]
                _log.error(conversion.map_proteins_to_genes)
                raise ValueError(f"Genes are missing for proteins {missing_genes}")
            genes_in_cluster = [x[0] for x in genes_in_cluster]
            if i == 3:
                _log.warning("Genes in Cluster: %s", genes_in_cluster)
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
    _log.debug("Found %d clusters across %d proteins", len(clusters_by_p), len(proteins_list))

    return MclClusterData(
        Gene2MCLClusters=clusters_by_g,
        MCL_Clusters=pl.DataFrame(mcl_cluster_data),
        Protein_IDs=proteins_list,
    )


def get_mcl_clustering(
    protein_sequence_fasta_filename: Path,
    conversion: Conversion,
    overwrite_mcl_clustering_file: bool,
    remove_non_protein_locations: bool,
    parallelize: int = 0,
    rebuild_blast_db: bool = True,
    rerun_blastp: bool = True,
    use_diamond: bool = True,
) -> MclClusterData:
    """---------------------------------------------------------------------------

      Function that uses a protein fasta file of a genome annotation to create
      a blastdb, and then does an *all against all* blast and uses the result to
      do a MCL clustering. This is used for local duplication information.

      MCL algorithm is for **Markov Cluster Algorithm,** which is an unsupervised
      clustering algorithm for graphs. For details, the project is hosted at
      https://micans.org/mcl/

    ---------------------------------------------------------------------------

    Args:
        protein_sequence_fasta_filename: Protein sequence file of the genome.
                                     Should only contain Protein IDs as
                                     header e.g. >protein1-ID
        conversion:
            MapG: Unique list of gene IDs stored in the conversion file.
            MapP: Unique list of protein IDs stored in the conversion file.
            MapG_MapP: Mapping between MapG and MapP.
        overwrite_mcl_clustering_file: Whether or not to overwrite the MCL
            clustering file.
        remove_non_protein_locations: True if protein IDs should be stored in return.
        parallelize: Defines if multiple cpus should be used for computation of parts
               of algorithm. Active if larger than 1.
        rebuild_blast_db: If BlastDB exists, should it be rebuilt? (recommended)
        rerun_blastp: if BlastP output exists, should it be rerun?
        use_diamond: Use Diamond instead of NCBI-BLAST for db and blastp steps

    Returns:
        MclClusterData object which provides mapping between
        genes in the conversion-ID file and the MCL-clustering file,
        as well as a list of protein-IDs stored in the fasta file.

    """
    _log.debug("Getting MCL clustering by running blastp and mcl")
    int_cpu_count: int = os.cpu_count() or 1
    if parallelize == 0:
        parallelize = max(1, int_cpu_count - 1)

    apps = ["makeblastdb", "blastp", "cut", "mcxload", "mcl", "mcxdump"]
    if use_diamond:
        _ = apps.pop(0)
        apps[0] = "diamond"
    missing_apps: list[str] = []
    app_path: dict[str, Path] = {}
    # verify all apps are in $PATH, generate error if so
    for app in apps:
        which_path = shutil.which(app)
        if which_path is None:
            missing_apps.append(app)
        else:
            app_path[app] = Path(which_path)
    if len(missing_apps) > 0:
        missing_apps_msg = ", ".join(missing_apps)
        raise OSError(f"The following applications are not in $PATH: {missing_apps_msg}")

    command_for: dict[str, list[str]] = {}
    blastdb_file = Path(f"{protein_sequence_fasta_filename}.blastDB")
    if use_diamond:
        blastdb_file = Path(f"{protein_sequence_fasta_filename}.diamondDB.dmnd")
    if rebuild_blast_db or not Path(blastdb_file.name + ".dmnd").is_file():
        skip_make_blastdb = False
        if use_diamond:
            command_for["diamond1"] = [
                "makedb",
                "--in",
                str(protein_sequence_fasta_filename),
                "--db",
                str(blastdb_file),
            ]
            if parallelize > 1:
                command_for["diamond1"].extend(["--threads", str(parallelize)])
        else:
            command_for["makeblastdb"] = [
                "-in",
                str(protein_sequence_fasta_filename),
                "-dbtype",
                "prot",
                "-out",
                str(blastdb_file),
                "-hash_index",
            ]
    else:
        _log.warning("Skipping execution of make blast db: output file exists")
        skip_make_blastdb = True
    blastp_file = Path(f"{protein_sequence_fasta_filename}.blastOUT")
    if rerun_blastp or not blastp_file.is_file():
        skip_blastp = False
        if use_diamond:
            command_for["diamond2"] = [
                "blastp",
                "--query",
                str(protein_sequence_fasta_filename),
                "--db",
                str(blastdb_file),
                "--evalue",
                "1",
                "--out",
                str(blastp_file),
                "--outfmt",
                "6",
                "qseqid",
                "sseqid",
                "evalue",
                "bitscore",
                "length",
                "pident",
                "qstart",
                "qend",
                "sstart",
                "send",
            ]
            if parallelize > 1:
                command_for["diamond2"].extend(["--threads", str(parallelize)])
        else:
            command_for["blastp"] = [
                "-query",
                str(protein_sequence_fasta_filename),
                "-db",
                str(blastdb_file),
                "-evalue",
                "1",
                "-out",
                str(blastp_file),
                "-outfmt",
                "6 qseqid sseqid evalue bitscore length pident qstart qend sstart send",
            ]
            if parallelize > 1:
                command_for["blastp"].extend(["-num_threads", str(parallelize)])
    else:
        _log.warning("Skipping execution of blastp: output file exists")
        skip_blastp = True
    command_for["cut"] = [
        "-f",
        "1,2,3",
        str(blastp_file),
    ]
    cut_output = Path(f"{protein_sequence_fasta_filename}.MCLIN")
    command_for["mcxload"] = [
        "-abc",
        str(cut_output),
        "--stream-mirror",
        "--stream-neg-log10",
        "-stream-tf",
        "ceil(200)",
        "-o",
        f"{protein_sequence_fasta_filename}.MCLMATRIX",
        "-write-tab",
        f"{protein_sequence_fasta_filename}.MCLTAB",
    ]
    command_for["mcl"] = [
        f"{protein_sequence_fasta_filename}.MCLMATRIX",
        "-I",
        "2",
        "-o",
        f"{protein_sequence_fasta_filename}.MCLMATRIX.I20",
    ]
    overwrite_dump_workaround = "" if overwrite_mcl_clustering_file else f"{os.getpid()}."
    command_for["mcxdump"] = [
        "-icl",
        f"{protein_sequence_fasta_filename}.MCLMATRIX.I20",
        "-tabr",
        f"{protein_sequence_fasta_filename}.MCLTAB",
        "-o",
        f"{protein_sequence_fasta_filename}.{overwrite_dump_workaround}MCLTAB.dump.I20",
    ]

    if use_diamond:
        if skip_make_blastdb and skip_blastp:
            apps.remove("diamond")
    else:
        if skip_make_blastdb:
            apps.remove("makeblastdb")
        if skip_blastp:
            apps.remove("blastp")
    for app in apps:
        if app == "cut":
            _log.info("executing: %s %s", app_path[app], " ".join(command_for[app]))
            with cut_output.open(mode="w") as f:
                subprocess.run(  # noqa: S603
                    [str(app_path[app])] + command_for[app], check=True, stdout=f
                )
        elif app == "diamond":
            if not skip_make_blastdb:
                _log.info("executing: %s %s", app_path[app], " ".join(command_for["diamond1"]))
                subprocess.run(  # noqa: S603
                    [str(app_path[app])] + command_for[app + "1"],
                    check=True,
                )
            if not skip_blastp:
                _log.info("executing: %s %s", app_path[app], " ".join(command_for["diamond2"]))
                subprocess.run(  # noqa: S603
                    [str(app_path[app])] + command_for[app + "2"],
                    check=True,
                )
        else:
            _log.info("executing: %s %s", app_path[app], " ".join(command_for[app]))
            subprocess.run(  # noqa: S603
                [str(app_path[app])] + command_for[app],
                check=True,
            )

    return read_mcl_clustering_file(
        Path(f"{protein_sequence_fasta_filename}.{overwrite_dump_workaround}MCLTAB.dump.I20"),
        protein_sequence_fasta_filename,
        conversion,
        remove_non_protein_locations,
    )


def map_pgdb_reactions_to_genes(  # noqa: D417
    pgdb_flat_file_folder: Path, conversion: Conversion, **kwargs
) -> GenesToReactions:
    """---------------------------------------------------------------------------

      Function reading out information of pgdb flatfiles genes.dat,
      proteins.dat, enzrxns.dat and reactions.dat, mapping genes
      to reactions directly.

      Get all protein IDs for each gene ID: genes.dat -> proteins.dat (one to many)
      Get all protein IDs in a complex for each protein ID: proteins.dat -> proteins.dat
      Get all reaction IDs for each protein ID: proteins.dat -> enzrxns.dat

    ---------------------------------------------------------------------------

    Args:
        pgdb_flat_file_folder: Full path to the folder where the flatfiles of
                              the species pgdb is stored.
        conversion: Conversion object with gene/transcript IDs to validate against
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
        kwargs.get("GenesAttributes", ("UNIQUE-ID", "PRODUCT", "ACCESSION-1", "ACCESSION-2"))
    )
    _proteins_attributes = tuple(
        kwargs.get("ProteinsAttributes", ("UNIQUE-ID", "CATALYZES", "COMPONENTS", "GENE"))
    )
    _enzrxns_attributes = tuple(kwargs.get("EnzrxnsAttributes", ("UNIQUE-ID", "REACTION", "ENZYME")))
    _reactions_attributes = kwargs.get(
        "ReactionsAttributes", ("UNIQUE-ID", "ENZYMATIC-REACTION", "EC-NUMBER")
    )
    gene_to_protein = kwargs.get("G2P", ("PRODUCT", "UNIQUE-ID"))
    protein_to_protein = kwargs.get("P2P", ("UNIQUE-ID", "COMPONENTS"))
    # ("CATALYZES", "COMPONENTS"))
    protein_to_enzyme = kwargs.get("P2E", ("CATALYZES", "UNIQUE-ID"))
    enzyme_to_reaction = kwargs.get("E2R", ("REACTION", "UNIQUE-ID"))
    # R2E = kwargs.get("R2E", ("UNIQUE-ID", "REACTION"))
    # this was in the original code, but it's not used
    # R2Pa = ("UNIQUE-ID", "ENZYMATIC-REACTION")

    # Guarantee that UNIQUE-ID is in the attributes
    genes_attributes = (
        _genes_attributes if "UNIQUE-ID" in _genes_attributes else _genes_attributes + ("UNIQUE-ID",)
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

    conversion_pgdb_mismatch = False
    # Read in all attributes and store them
    _log.info("Reading in genes.dat")
    genes_all = read_in_pgdb_flat(
        pgdb_flat_file_folder / "genes.dat",
        genes_attributes,
        unique_keys=genes_attributes,
    )
    validate_genes = genes_all.join(conversion.df, left_on="UNIQUE-ID", right_on="Gene Name")
    if len(validate_genes) != len(genes_all):
        _log.error(
            "%d Genes of %d in PGDB file not in conversion file",
            len(genes_all) - len(validate_genes),
            len(genes_all),
        )
        _log.debug("Genes All %s", genes_all.head())
        conversion_pgdb_mismatch = True
    _log.info("Reading in proteins.dat")
    proteins_all = read_in_pgdb_flat(pgdb_flat_file_folder / "proteins.dat", proteins_attributes)
    validate_prots = proteins_all.join(conversion.df, left_on="UNIQUE-ID", right_on="Protein Name")
    if len(validate_prots) != len(proteins_all):
        _log.error(
            "%d Proteins of %d in PGDB file not in conversion file",
            len(proteins_all) - len(validate_prots),
            len(proteins_all),
        )
        _log.debug("Proteins all %s", proteins_all.head())
        conversion_pgdb_mismatch = True
    if conversion_pgdb_mismatch:
        _log.error(
            "GFF3 pre-processing results don't match PGDB data. "
            + "The pre-processing settings may need adjustment"
        )
        _log.debug("Conversion DB head: %s", conversion.df.head())
    _log.info("Reading in enzrxns.dat")
    enzrxns_all = read_in_pgdb_flat(pgdb_flat_file_folder / "enzrxns.dat", enzrxns_attributes)
    _log.info("Reading in reactions.dat")
    reactions_all = read_in_pgdb_flat(pgdb_flat_file_folder / "reactions.dat", reactions_attributes)

    _log.info("Matching proteins to protein-complexes")
    # has GENE (GENE-ID), and enzrxn_id
    components_col_idx = proteins_all.columns.index("COMPONENTS")
    if proteins_all.dtypes[components_col_idx] == pl.Null:
        protein_complexes_df = proteins_all.with_columns(enzrxn_id=pl.col("CATALYZES"))
    else:
        protein_complexes_df = proteins_all.join(
            proteins_all.filter(pl.col("COMPONENTS").is_not_null()).explode(pl.col("COMPONENTS")),
            left_on="UNIQUE-ID",
            right_on="COMPONENTS",
            how="left",
        ).with_columns(
            enzrxn_id=pl.col("CATALYZES").list.concat(
                pl.coalesce(pl.col("CATALYZES_right"), pl.lit([]))
            )
        )
    _log.debug("Protein-Protein Complex DF: %s", protein_complexes_df)

    _log.info("Match proteins and protein complexes to enzyme reactions")
    protein_enz_reaction_df = (
        protein_complexes_df.explode("enzrxn_id")
        .join(enzrxns_all, left_on="enzrxn_id", right_on="UNIQUE-ID")
        .select(pl.col(["GENE", "REACTION"]))
    )
    _log.debug("Protein Enzyme DF: %s", protein_enz_reaction_df)

    gene2rxn = (
        protein_enz_reaction_df.with_columns(pl.col("GENE").list.first().alias("GENE"))
        .group_by(pl.col("GENE"))
        .agg(pl.col("REACTION"))
        .rows_by_key("GENE", unique=True)
    )
    _log.info(
        "map_pgdb_reactions_to_genes Found %d genes, %d with reactions",
        len(genes_all),
        len(
            protein_enz_reaction_df.group_by(pl.col("GENE"))
            .agg(pl.col("REACTION"))
            .filter(pl.col("REACTION").list.len() > 0)
        ),
    )
    _log.debug("Genes ALl %s", genes_all)
    _log.debug("Genes to Reactions DF %s", protein_enz_reaction_df)

    return GenesToReactions(GenesAll=genes_all, ReactionsAll=reactions_all, Gene2Rxn=gene2rxn)


def read_in_conversion_file(
    gene_transcript_protein_mapping_filename: Path,
) -> Conversion:
    """---------------------------------------------------------------------------

      Function that reads in the gene, transcript, and protein ID conversion
      file and stores the IDs as unique lists and mapping between genes and
      transcripts aswell as genes and proteins.

    ---------------------------------------------------------------------------

    Args:
        gene_transcript_protein_mapping_filename: TSV file with a header describing
                                             gene, transcript, and protein
                                             identifier and then for all the
                                             identfiers the listings (gene-ID,
                                             transcript-ID, protein-ID).

    Returns:
        Conversion object containing the IDs from the conversion file
    Raises

    """
    _log.debug("Reading in conversion file %s", gene_transcript_protein_mapping_filename)
    conversion_data = pl.read_csv(
        gene_transcript_protein_mapping_filename,
        separator="\t",
        schema_overrides=dict.fromkeys(["Gene Name", "Transcript Name", "Protein Name"], pl.String),
    )
    column_names = conversion_data.columns
    conversion_len = len(conversion_data)
    _log.debug("Read in conversion data, %d rows, columns are %s", conversion_len, column_names)
    conversion = Conversion(df=conversion_data)
    # Validate conversion before output
    if len(conversion.map_gene_to_transcript) < len(conversion.map_gene):
        _log.warning("Some Transcripts match multiple genes in the ID Conversion file")
        transcript_multi = (
            conversion_data.group_by(column_names[1])
            .n_unique()
            .filter(pl.col(column_names[0]) > 1)[column_names[1]]
        )
        for transcript in transcript_multi:
            _log.warning("Please check the transcript %s", transcript)
        _log.error("Please clean the protein ID file")
    if len(conversion.map_gene_to_protein) < len(conversion.map_gene):
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
        "Number of Gene-ID to Transcript-ID mappings: %i", len(conversion.map_gene_to_transcript)
    )
    _log.info("Number of Gene-ID to Protein mappings: %i", len(conversion.map_gene_to_protein))
    return conversion


def match_pgdbs_to_mds(
    gene_to_reaction: GenesToReactions,
    annotated_metabolic_domains: AnnotatedMetabolicDomains,
) -> dict:
    """---------------------------------------------------------------------------

      Function that maps the reactions stored in the pgdbs to the reactions
      stored in the metabolic domain file.

    ---------------------------------------------------------------------------

    Args:
        gene_to_reaction: Mapping of all genes in the Genes_all to the Reactions in
              Reactions_all.
        annotated_metabolic_domains: Mapping between all Reactions stored in the
            Metabolic domain file and the metabolic domains that were identified
            to be of interest in metabolic_domains_to_annotate.
            And DataFrame of all the items stored in the Metabolic domain file

    Returns:
        A mapping between reactions stored in the pgdb and the
                    reactions given in the metabolic domain file. Used by write_gene_output.

    """
    _log.debug("Matching PGDBs to annotated metabolic domains")

    mds = (
        gene_to_reaction.ReactionsAll.join(
            annotated_metabolic_domains.RxnMetaDom,
            left_on="UNIQUE-ID",
            right_on="Reaction ID",
        )
        .select(pl.col("UNIQUE-ID").alias("REACTION-ID"), pl.col("Metabolic Domain"))
        .group_by(pl.col("REACTION-ID"))
        .agg(pl.col("Metabolic Domain"))
    )
    _log.debug("Metabolic Domains dataframe: %s", mds)
    mds_dict = mds.to_dict(as_series=False)
    return dict(zip(mds_dict["REACTION-ID"], mds_dict["Metabolic Domain"], strict=True))


def map_pgdb_gene_ids_to_conversion_file(
    conversion: Conversion,
    all_genes_df: pl.DataFrame,
    genes_data_metadata: list[str],
    pgdb_ids_to_map: str,
) -> pl.DataFrame:  # NOTE: this is what creates MapG_Gpgdb
    """---------------------------------------------------------------------------

      Function that maps the conversion file gene IDs to the IDs in the
      genes.dat of the pgdb.

    ---------------------------------------------------------------------------

    Args:
        conversion: Conversion object from conversion file
        all_genes_df: DataFrame containing all the information
                    (Attributes)that was searched for (described in
                    GenesAttributes, default is {'UNIQUE-ID','PRODUCT', ...
                    'ACCESSION-1','ACCESSION-2'};) in genes.dat of the pgdb.
        genes_data_metadata: List of Attributes of genes.dat. Defines where gene
                            information should be searched for.
        pgdb_ids_to_map: A string that contains G and or T and or P and defines
                       if the entries represent Genes and or transcript ids and
                       or proteins ids of the conversion file.

    Returns:
        A dictionary mapping between Genes in the conversion file and genes
                     stored in the pgdb genes.dat file.

    """
    _log.debug("Mapping PGDBs to conversion file")
    if pgdb_ids_to_map not in "GTP":
        _log.error("Invalid PGDBIdsToMap mapping value, valid values are G, T, and P")
        raise ValueError("Invalid Mapping Value")

    gene_map = pl.DataFrame(schema={"Gene Name": pl.String, "UNIQUE-ID": pl.String})
    # PGDB GeneID => Gene Name
    # expects just Gene Names
    if "G" in pgdb_ids_to_map:
        _log.debug("PGDB Genes Sample: %s", all_genes_df.head())
        all_genes_df_genes = all_genes_df.select(pl.col("UNIQUE-ID").alias("Gene Name")).unique()
        gene_map = conversion.df.join(all_genes_df_genes, on="Gene Name").select(
            pl.col("Gene Name"), pl.col("Gene Name").alias("UNIQUE-ID")
        )
        _log.debug("Mapping genes, found %d Gene Names to map", gene_map.height)
        _log.debug("Example: %s", gene_map)

    transcript_map = pl.DataFrame(schema={"Transcript Name": pl.String, "UNIQUE-ID": pl.String})
    # Pgdb Protein ID => Gene Name
    # pull from "ACCESSION-1"
    _log.warning("Genes_Data_Metadata %s, ids to map %s", repr(genes_data_metadata), pgdb_ids_to_map)
    if "T" in pgdb_ids_to_map and len(genes_data_metadata) > 0:
        all_genes_df_transcript = (
            all_genes_df.select(pl.col("UNIQUE-ID"), pl.col(genes_data_metadata[0]))
            .explode(pl.col(genes_data_metadata[0]))
            .unique()
        )
        _log.debug("Transcript Map Genes: %s", all_genes_df_transcript.head())
        transcript_map = conversion.df.join(
            all_genes_df_transcript,
            left_on="Transcript Name",
            right_on=genes_data_metadata[0],
        ).select(pl.col("Transcript Name"), pl.col("UNIQUE-ID"))
        _log.debug(
            "Mapping transcript, found %d Transcript Names to map",
            transcript_map.height,
        )

    protein_map = pl.DataFrame(schema={"Protein Name": pl.String, "UNIQUE-ID": pl.String})
    # Pgdb Protein ID => Gene Name
    # pull from PRODUCT
    if "P" in pgdb_ids_to_map and len(genes_data_metadata) > 1:
        all_genes_df_proteins = (
            all_genes_df.select(pl.col("UNIQUE-ID"), pl.col(genes_data_metadata[1]))
            .explode(genes_data_metadata[1])
            .unique()
        )
        protein_map = conversion.df.join(
            all_genes_df_proteins, left_on="Protein Name", right_on=genes_data_metadata[1]
        ).select(pl.col("Protein Name"), pl.col("UNIQUE-ID"))
        _log.debug("Mapping proteins, found %d Protein Names to map", protein_map.height)

    conversion_mapped_df = (
        all_genes_df.join(gene_map, on="UNIQUE-ID", how="left")
        .join(transcript_map, on="UNIQUE-ID", how="left")
        .join(protein_map, on="UNIQUE-ID", how="left")
        .select(
            pl.col("UNIQUE-ID").alias("pgdb_gene_id"),
            pl.col("Gene Name"),
            pl.col("Protein Name"),
            pl.col("Transcript Name"),
        )
    )
    _log.debug("Conversion data mapped: %s", conversion_mapped_df)
    return conversion_mapped_df


def evaluate_criterion(
    criterion: int,
    clusters_df: pl.DataFrame,
    step_size: int,
    maximum_step_size: int,
    minimum_step_size: int = 3,
) -> tuple[int, bool]:
    """Evaluate search continuation criterion for cluster search.

    ---

    Args:
        criterion: A value from 1 to 4 describing search criterion.
                   - 1: as published in Schlapfer et al 2017
                       sum(enzymes in clusters) - sum(non-enzymes in clusters)
                       > 0.
                   - 2: stepsize 3 (because most organisms take this one)
                   - 3: sum(enzymes in clusters) - sum(non-enzymes in clusters,
                       without hypo genes) -2*sum(hypogenes)
                   - 4: sum(clusters in stepsize n) > sum(clusters in stepsize
                       n+1).
        clusters_df: DataFrame of discovered clusters.
        step_size: Current search step size for clusters.
        maximum_step_size: Maximum permitted step size.
        minimum_step_size: Minimum permitted step size.

    Returns:
        Tuple with the final step size, and if the search for clusters should continue

    """
    continue_search = True
    final_step_size = step_size - 1
    match criterion:
        case 1:
            # (Schlapfer et al. 2017 criteria)
            # sum (enzymes in clusters) - sum(non-enzymes in clusters) > 0
            #   (incl. hypo genes?))
            # i.e. more enzymes than non-enzymes (incl hypo) in discovered clusters
            # TODO(j-schlesinger): validate criterion (issue #34)  # noqa: FIX002
            # https://github.com/j-schlesinger/PlantClusterFinder.py/issues/34
            enzymes_in_step_size = (
                clusters_df.filter(pl.col("step_size") == step_size)
                .select(pl.col("enzyme_count"))
                .sum()
            )
            _log.debug(
                "Number of enzymes - number of non-enzymes in clustered genes: %d",
                enzymes_in_step_size,
            )
            if step_size == maximum_step_size:
                continue_search = False
                final_step_size = step_size
                _log.warning(
                    "Maximal step size reached: %i. Consider increasing.",
                    maximum_step_size,
                )
            elif clusters_df.filter(pl.col("step_size") == step_size).height == 0:
                _log.warning("No clusters found for step size %i.", step_size)
            elif (
                clusters_df.filter(pl.col("step_size") == step_size).select(
                    (pl.col("enzyme_count").sum() / pl.col("non_enzyme_gene_count").sum()).alias(
                        "enz_over_genes"
                    )
                )["enz_over_genes"][0]
                < 0.5
            ):
                if continue_search:
                    final_step_size = step_size - 1
                continue_search = False
        case 2:
            # step size 3 (because most organisms take this one)
            _log.debug(
                "Number of enzymes - number of non-enzymes (no hypos) - 2* number"
                + " if hypos in clustered genes: %i - %i - %i",
                clusters_df["enzyme_count"].count(),
                -clusters_df["non_enzyme_gene_count"].count(),
                -2 * clusters_df["hypothetical_gene_count"].count(),
            )
            if step_size >= minimum_step_size:
                final_step_size = 3
                continue_search = False
        case 3:
            # sum(enzymes in clusters) - sum(non-enzymes in clusters without hypo genes)
            #     - 2*sum(hypo genes) < 0
            # i.e. enzymes < non-enzymes + (2x hypo genes)
            if step_size == maximum_step_size:
                continue_search = False
                final_step_size = step_size
                _log.warning(
                    "Maximal step size reached: %i. Consider increasing.",
                    maximum_step_size,
                )
            elif clusters_df.filter(pl.col("step_size") == step_size).height == 0:
                _log.warning("No clusters found for step size %i.", step_size)
            elif (
                clusters_df.filter(pl.col("step_size") == step_size).select(
                    (
                        pl.col("enzyme_count").sum()
                        - (
                            pl.col("non_enzyme_gene_count").sum()
                            - pl.col("hypothetical_gene_count").sum()
                            # + 1
                        )
                        # - pl.col("enzyme_count").sum()
                        - 2 * pl.col("hypothetical_gene_count").sum()
                    ).alias("enz_less_genes")
                )["enz_less_genes"][0]
                < 0
            ):
                if continue_search:
                    final_step_size = step_size - 1
                continue_search = False
        case 4:
            # sum(count of clusters in step size n) >
            #    sum(count of clusters in step size n + 1)
            # i.e. only continue searching while finding more
            #      clusters than the last step size
            if step_size != 0:
                _log.debug(
                    "Number of clusters more than last stepsize %i",
                    clusters_df.filter(pl.col("step_size") > (step_size - 1)).height,
                )
            else:
                _log.debug(
                    "Number of clusters Stopsize 0: %i",
                    clusters_df.filter(pl.col("step_size") == 0).height,
                )
            if step_size == maximum_step_size:
                continue_search = False
                final_step_size = step_size
                _log.warning(
                    "Maximal step size reached: %i. Consider increasing.",
                    maximum_step_size,
                )
            elif clusters_df.filter(pl.col("step_size") == step_size).height == 0:
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
            raise ValueError("Invalid search criterion, valid options are 1,2,3 or 4.")
    return (final_step_size, continue_search)


def perform_clustering(
    gene_location: GeneLoc,
    mcl_cluster_data: MclClusterData,  # covered by genelocation
    minimum_enzyme_for_clusters: int,
    criterion: int,  # TODO(j-schlesinger): this should be an enum (issue #32)  # noqa: FIX002
    # https://github.com/j-schlesinger/PlantClusterFinder.py/issues/32
    minimum_step_size: int = 5,
    maximum_step_size: int = 20,
) -> Clusters:
    """---------------------------------------------------------------------------

      Function that clusters the genome. It looks at enzymes that are at max
      spaced by step_size genes (checking for 0, 1, ..., minimum_step_size intervening
      genes. It continuous to increase the step_size until clusters contain
      more non-enzymes than enzymes. Clusters containing only one group of
      locally duplicated genes are rejected. Clusters with only one reaction
      are rejected as well. Clusters need at least minimum_enzyme_for_clusters enzymes.

    ---------------------------------------------------------------------------

    Args:
        gene_location: DataFrame about the gene location information.
        mcl_cluster_data: Mapping between the genes in the conversion-ID file and the
                MCL-clustering file.
        minimum_enzyme_for_clusters: Defines the minimum number of enyzmes that need to
                            be present in a cluster.
        criterion: Criterion after which step size is chosen
                   - 1: as published in Schlapfer et al 2017
                       sum(enzymes in clusters) - sum(non-enzymes in clusters)
                       > 0.
                   - 2: stepsize 3 (because most organisms take this one)
                   - 3: sum(enzymes in clusters) - sum(non-enzymes in clusters,
                       without hypo genes) -2*sum(hypogenes)
                   - 4: sum(clusters in stepsize n) > sum(clusters in stepsize
                       n+1).
        minimum_step_size: Minimal intervening gene size that should be computed
                      (default is 5).
        maximum_step_size: Maximal intervening gene size that should be computed
                      (default is 20).

    Returns:
        Cluster object which contains the following information:
        - **Clusters**: Cluster Information (step_size, chromosome, start,
                   end, enzyme_count, non_enzyme_gene_count, hypothetical_gene_count).
        - **ClusterBackground**: Density information about the genome (window size,
                            number of enzymes, Occurrence of windowsize with
                            this many enzymes in the genome, occurence of
                            this windowsize in the genome, Frequency of seeing
                            a windowsize of this size with at least as many
                            enzymes in the genome).
        - **FinalStepsize**: Stepsize that was chosen by the criterion (clusters
                        contain in average more enzymes than non-enzymes).
        - **StepSize_Max**: Maximal stepsize that was computed.

    """
    if not (
        gene_location.is_annotated or gene_location.is_converted or not gene_location.has_hard_breaks
    ):
        _log.error("Gene Location does not have required transformation(s)")
        if not gene_location.is_annotated:
            _log.error("- Gene Annotation")
        if not gene_location.is_converted:
            _log.error(" - Conversion file ID mappings")
        if not gene_location.has_hard_breaks:
            _log.error(" - Hard break sections")
        raise AttributeError(
            f"Gene Location {gene_location} is not transformed " + "for cluster identification"
        )
    _log.debug("Performing clustering")
    _log.debug(
        "perform_clustering inputs:\n\t"
        + f"gene_location: {gene_location}\n\t"
        + f"mcl_cluster_data: {mcl_cluster_data.MCL_Clusters}\n\t"
        + f"minimum_enzyme_for_clusters: {minimum_enzyme_for_clusters}\n\t"
        + f"criterion: {criterion}\n\t"
        + f"minimum_step_size: {minimum_step_size}\n\t"
        + f"maximum_step_size: {maximum_step_size}"
    )
    final_step_size = 0
    continue_search = True  # flag if search should continue
    clusters = []  # array of cluster identifiers
    # named tuple to simplify creation of clusterdf
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
    clusters_background = pl.DataFrame()  # Size, number of enzymes, count, number of total
    step_size = 1  # step size for search
    # note: hypothetical gene names start with "Hypo", this fact was originally tracked in a variable

    # used in place of biomG_E
    is_enzyme_criteria = pl.col("has_reaction")

    # get rid of clusters that only contain Local Duplicates

    gene_count = len(gene_location)
    # is there a way to make this search an expression?
    while step_size <= minimum_step_size or continue_search:  # for each valid step size
        cluster_starts_df = (
            gene_location.data.with_row_index("index")
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
        _log.debug(f"Found {cluster_starts_count} cluster starts for step size {step_size}")
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
            stop_search_end = False
            while not stop_search_end:
                end_s = end  # end search start
                # stop search step size
                end_e = min(end + step_size, len(gene_location.data))
                _log.debug(f"Searching for enzyme cluster for index {start} from {end_s} to {end_e}")
                end_df = gene_location.data.with_row_index("index").filter(
                    is_enzyme_criteria
                    & chromosome_criteria
                    & (pl.col("index") > end_s)
                    & (pl.col("index") <= end_e)
                )
                if end_df.height > 0:
                    end = end_df.select(pl.last("index"))["index"][0]
                else:
                    stop_search_end = True

            cluster_genes = gene_location.data.with_row_index("index")[start:end]
            # get rid of clusters that only have one reaction set
            if len(cluster_genes.filter(is_enzyme_criteria)) < minimum_enzyme_for_clusters:
                # abolish cluster and continue search
                continue
            if (
                cluster_genes.filter(pl.col("is_locally_duplicated"))
                .select(pl.col("mcl_clusters"))
                .unique()
                .height
                == 1
            ):
                # get rid of clusters that only have exactly one local duplicate.
                # abolish cluster and continue search
                continue

            # This is the code that writes the cluster into the clusters dataframe
            # what goes in here is described in line 3426 in Matlab, but it is flexible
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
                    _cluster.filter(pl.col("Hypothetical").not_()).height,  # hypothetical_gene_count
                )
            )

        if len(clusters) > 0:
            clusters_df = pl.DataFrame(clusters)
            _log.debug("Head of clusters_df %s", clusters_df.head())
            # recalculate clusters_background to test against criterion
            clusters_background = compute_background(clusters_df, gene_location, clusters_background)
            final_step_size, continue_search = evaluate_criterion(
                criterion, clusters_df, step_size, maximum_step_size, minimum_step_size
            )
        else:
            final_step_size = step_size - 1
            continue_search = False
        step_size += 1
    if final_step_size == 0:
        final_step_size = step_size - 1
    clusters_obj = Clusters(
        ClusterMeta=pl.DataFrame(clusters),
        ClusterBackground=clusters_background,
        FinalStepSize=final_step_size,
        MaxStepSize=step_size - 1,
    )
    _log.info("Found %s clusters.", len(clusters))
    return clusters_obj


def compute_background(
    clusters: pl.DataFrame,
    gene_location: GeneLoc,
    clusters_background: pl.DataFrame,
) -> pl.DataFrame:
    """---------------------------------------------------------------------------

      Function that looks at all clusters that were found by the function
      perform_clustering(), and computes the background information for the
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
        clusters_background: DataFrame of backgrounds for up to clusters
                max_step_size - 1. (window_size, enzyme_count,
                                    window_occurrence, winzyme_occurrence, winzyme_frequency)
                "Winzyme" is a window with enzyme_count enzymes in it.

    Returns:
        - Density information about the genome (window size,
                            number of enzymes. Occurrence of window size with
                            this many enzymes in the genome, occurrence of
                            this window size in the genome, Frequency of seeing
                            a window size of this size with at least as many
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
        pl.sum_horizontal("enzyme_count", "non_enzyme_gene_count").max().alias("window_size_max")
    )["window_size_max"][0]
    _log.debug("Max Window %f", max_window)
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
                winzyme_frequency=pl.col("cumulative_winzyme") / total_window_occurrence,
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
    metabolic_domains_to_annotate: list[str],
    pgdb_to_metabolic_domain_map: dict[Any, Any],
    reaction_labels: ReactionLabels,
    output_type: str,
) -> Path:
    """---------------------------------------------------------------------------

      Function that writes the gene file of PlantClusterFinder: TSV format,
      first it includes the information that was stored in the gene location
      file that was provided to PlantClusterFinder (Gene IDs, Gene Start, Gene
      end, Chromosome information, strand). Then it attaches Enzyme (1 or 0) EC
      information and Reaction information. Then whether the Gene is considered
      locally duplicated and if so, which LD-clusters it's proteins were
      associated with. Then per step size that was investigated (intervening
      gene size) the ClusterID where the gene is incorporated and the
      information whether the cluster was found to be in the top n% of
      clusters. This is followed by the information about the metabolic
      domains, and the signature / tailoring information.

    ---------------------------------------------------------------------------

    Args:
        gene_output_file: Path and file name to the gene-Outputfile, make sure
                         you have access to the file and folder. This is the
                         first Results output file that will be generated.
        gene_location: Gene location information, generated in
                        PlantClusterFinder function.
        clusters: Cluster Information (step size, chromosome id, gene start, gene end,
            number of enzymes).
        top_percent_clusters: Clusters are labeled as being in the top n% if
                             they are within this top n% of enzyme dense
                             clusters compared to the background (default is
                             5).
        metabolic_domains_to_annotate: A list of Metabolic domains that should
                         be analyzed. For example {'Amines and Polyamines Metabolism'; ...
                         'Amino Acids Metabolism'
        annotated_metabolic_domains: Mapping between all Reactions stored in
               the Metabolic domain file and the metabolic domains that were
               identified to be of interest in metabolic_domains_to_annotate.
        pgdb_to_metabolic_domain_map: A mapping between reactions stored in the pgdb and the
                    reactions given in the metabolic domain file
        reaction_labels: Is a data structure of signature and tailoring enzyme
                    classes.
                    1. Mapping of Reactions to the enzyme classes of SigOrTail Classes.
                    2. Is a list of types of enzymes, e.g. signature, tailoring,
                      and so on.
                    3. Mapping of Reactions in Reactions to the enzyme types of SigOrTail.
        output_type: Defines which type of output is given:
                        "simple" (default setting) gives only top n% for the
                           final step size.
                        "verbose" gives all computed step size until the final
                           step size.
                        "old" gives at least 5 until the final step size.

    Returns:
        Path of file generated

    """
    _log.debug("Writing gene output")
    if output_type not in ["old", "simple", "verbose"]:
        raise ValueError("Output type not supported, options are 'old' 'simple' or 'verbose'")

    gene_columns = gene_location.data.columns
    _log.debug("Columns for gene location data are %s", gene_columns)
    _log.debug("Head for gene location data is %s", gene_location.data.head())
    gene_output = gene_location.data
    # the first several columns are straight from gene_location
    _log.debug(
        "PGDB Reaction IDs: %s ",
        gene_output.select(pl.col("pgdb_reaction_id")).filter(
            pl.col("pgdb_reaction_id").is_not_null()
        ),
    )
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
    _log.debug("Cluster data for gene output %s", cluster_data.head())
    # different output columns depending on output type
    match output_type:
        case "old":
            for i in range(clusters.MaxStepSize):
                # chosen step size x
                # top pct clusters
                gene_output = gene_output.join(
                    cluster_data.agg(
                        pl.col("cluster_idx").implode().alias(f"cluster_idx_stepsize{i}"),
                        pl.col("top_percent_clusters").sum(),
                    )
                    .filter(pl.col("step_size") == pl.lit(i))
                    .select(
                        pl.col("Gene Name"),
                        pl.col(f"cluster_idx_stepsize{i}"),
                        pl.col("top_percent_clusters").alias(f"top_percent_clusters{i}"),
                    ),
                    "Gene Name",
                    how="left",
                )
                gene_output_cols.append(
                    pl.col(f"cluster_idx_stepsize{i}")
                    .list.len()
                    .alias(f"{'CHOSEN ' if i == clusters.FinalStepSize else ''}Stepsize {i}")
                )
                gene_output_cols.append(
                    pl.col(f"top_percent_clusters{i}").alias(f"Top {top_percent_clusters}% {i}")
                )
        case "verbose":
            for i in range(clusters.MaxStepSize):
                # chosen step size x
                # top pct clusters
                gene_output = gene_output.join(
                    cluster_data.agg(
                        pl.col("cluster_idx").implode().alias(f"cluster_idx_stepsize{i}"),
                        pl.col("top_percent_clusters").sum().alias(f"top_percent_clusters{i}"),
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
                for j in range(clusters.MaxStepSize):
                    gene_output_cols.append(
                        pl.col(f"cluster_idx_stepsize{i}")
                        .list.len()
                        .alias(f"{'CHOSEN ' if i == clusters.FinalStepSize else ''}Stepsize {i},{j}")
                    )
                    gene_output_cols.append(
                        pl.col(f"top_percent_clusters{i}").alias(
                            f"Top {top_percent_clusters}% {i},{j}"
                        )
                    )
        case "simple":
            # step size x
            # top pct clusters
            top_pct_cluster_idx = cluster_data.agg(
                pl.col("cluster_idx").count().alias(f"Stepsize {clusters.FinalStepSize}"),
                pl.col("top_percent_clusters").sum().alias(f"Top {top_percent_clusters}%"),
            ).filter(pl.col("step_size") == pl.lit(clusters.FinalStepSize))
            gene_output = gene_output.join(top_pct_cluster_idx, "Gene Name", how="left")

    # describe metabolic domains to annotate for run
    for md in metabolic_domains_to_annotate:
        gene_output_cols.append(
            pl.col("pgdb_reaction_id")  # for all reaction IDs in gene
            .list.eval(  # noqa s307
                # for all metabolic domains in gene
                pl.element()
                .replace_strict(pgdb_to_metabolic_domain_map, default=[])
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
        _log.debug("Column %d %s", i, column.meta.output_name())
        try:
            _log.debug(" explain: %s", gene_output.select(column).schema)
        except Exception as e:
            _log.warning(" column %d failed to explain %s", i, e)
    gene_output = gene_output.select(gene_output_cols)
    _log.info("Gene Output is %d rows", len(gene_output))
    _log.debug("Gene Output head %s", gene_output.head())
    gene_output.write_csv(gene_output_file, separator="\t")
    return gene_output_file


def write_cluster_output(
    cluster_output_file: Path,
    gene_location: GeneLoc,
    clusters: Clusters,
    output_type: str,
) -> Path:
    """---------------------------------------------------------------------------

      Function that writes the cluster file of PlantClusterFinder in TSV format.
      First it includes the step size of the cluster (intervening genes), then
      cluster ID, chromosome/scaffold the cluster is located on, start BP, end
      BP, physical size, start gene ID, end gene ID, and all the gene IDs.

    ---------------------------------------------------------------------------

    Args:
        cluster_output_file: Path and file name to the cluster output file.
                            Make sure you have access to the file and folder.
                            This is the second Results output file that will be
                            generated.
        gene_location: DataFrame about the gene location information,
                        generated in PlantClusterFinder function.
        clusters: Cluster Information (step size, chromosome id, gene start,
                   gene end, number of enzymes).
        clusters_in_top_percent: List with size Number of Clusters x 1. It is zero if
                       the density of the cluster is above or the same as given
                       by the input argument top_percent_clusters*100.
        output_type: Defines which type of output is given: simple (default
                     setting) gives only top n% for the final stepsize. verbose
                     gives all computed stepsize until the final stepsize. old
                     gives at least 5 until the final stepsize.

    Returns:
        Path of cluster outpuf file

    """
    _log.debug("Writing cluster output")
    if output_type not in ["old", "simple", "verbose"]:
        raise ValueError("Output type not supported, options are 'simple' or 'verbose'")
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
            pl.col("Gene End (bp)").max().alias("end_bp"),
        )
    )
    test_data = indexed_clusters.select(
        pl.col("cluster_id"), pl.col("chromosome"), pl.col("start"), pl.col("end")
    )
    _log.debug("Clusters Data ID, chromosome, and start and end %s", test_data)
    _log.debug("Cluster Genes Only %s", cluster_genes.head())
    _log.debug("Cluster Genes Rows %d", len(cluster_genes))
    _log.debug("Cluster Genes Columns %s", cluster_genes.columns)
    cluster_output = indexed_clusters.join(cluster_genes, on="cluster_id")
    _log.debug("Cluster output length %d", len(cluster_output))
    _log.debug("Cluster Step Sizes %s", cluster_output["step_size"].describe())
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
            pl.col("enzyme_count").alias("Number of Enzymes"),
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
    _log.debug("Cluster output dtypes %s", cluster_output.dtypes)
    _log.debug("Cluster output %s", cluster_output)
    cluster_output.write_csv(cluster_output_file, separator="\t")
    return cluster_output_file


def analyze_plant_cluster_gap_file(
    plant_cluster_gap_filename: Path,
    dna_unmasked_ok: bool = True,
    file_output: bool = True,
) -> tuple[str, dict[Any, list[Any]]] | Path:
    """---------------------------------------------------------------------------

      Function that analyzes the masked DNA sequences of a _GAPOutput file, and
      counts the occurrence of every size of gap.

    ---------------------------------------------------------------------------

    Args:
        plant_cluster_gap_filename: Full path to the _GAPOutput file that would be analyzed
                 (produced within PlantClusterFinder)
        dna_unmasked_ok: Option to ignore the error if you use unmasked DNA files.
            (Default True)
        file_output: Whether output is directed to file or to

    Returns:
        **Optional[str]** (File with list of sequencing gaps numbers is generated with
        ending _count.txt if file_output is set to True)

    """
    if not plant_cluster_gap_filename.is_file():
        raise OSError(
            f"plant_cluster_gap_filename {plant_cluster_gap_filename}"
            + " is not a file or does not exist"
        )
    _gap_output = (
        pl.read_csv(
            plant_cluster_gap_filename,
            separator=" ",
            has_header=False,
            new_columns=["scaffold", "amino acid", "start index", "length"],
        )
        .filter(pl.col("amino acid") == "N")
        .group_by("length")
        .len()
        .rows_by_key(key="length", unique=True)
    )
    f_out_filename = plant_cluster_gap_filename.with_suffix(
        plant_cluster_gap_filename.suffix + "_count.txt"
    )
    if not file_output:
        return (str(f_out_filename), _gap_output)
    with f_out_filename.open(mode="w") as f:
        f.write("Gapsize in bp\tNumber of Gaps\n")
        for k, v in sorted(_gap_output.items(), key=lambda x: -list(x)[1][0]):
            f.write(f"{k}\t{list(v)[0]}\n")
    return f_out_filename


# read_in_pgdb_flat moved to pgdb_handler.py
