import csv
import logging
import os
import subprocess
from collections import namedtuple
from pathlib import Path
from typing import IO, Any, Dict, Iterable, Iterator, List, Optional, Tuple

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

log = logging.getLogger(__name__)


def label_reactions_with_sig_and_tail(
    vReactions_all: pl.DataFrame, vSignatureTailorFile: Path
) -> ReactionLabels:
    """
    ---------------------------------------------------------------------------

      Function identifying wether reactions found in the pgdb are of a certain
      type of enyzme and if they are either signature or tailoring enzymes.
      Note, there have to be signature enzymes identified as 'signature' in the
      File vSignatureTailorFile.

    ---------------------------------------------------------------------------
    Args:
        vReactions_all: DataFrame containing all the information
                        (Attributes)that was searched for (described in
                        vReactionsAttributes, default is {'UNIQUE-ID', ...
                        'ENZYMATIC-REACTION', 'EC-NUMBER'};) in reactions.dat
                        of the pgdb.
        vSignatureTailorFile: TSV file describing the classes of enzymes that
                              should be classified as signature or tailoring
                              enzymes to identify clusters containing such.
    Return:
        ReactionLabels object with the signature and tailoring enzyme classes,
            list of types of enzymes, and mapping between reactions and
            enzyme classes as well as signature or tailoring enymes
    """
    stdata = pl.read_csv(
        vSignatureTailorFile,
        has_header=True,
        separator="\t",
        encoding="ISO-8859-1",
        truncate_ragged_lines=True,
    )
    cols = stdata.columns
    stclasses = stdata[cols[1]].unique()
    sigortail = stdata[cols[2]].unique()
    rxn_class = stdata.to_pandas().groupby(cols[0])[cols[1]].apply(list).to_dict()
    rxn_sig_or_tail = stdata.to_pandas().groupby(cols[0])[cols[2]].apply(list).to_dict()
    return ReactionLabels(
        STClasses=list(stclasses),
        SigOrTail=list(sigortail),
        Rxn_STClasses=rxn_class,
        Rxn_SigOrTail=rxn_sig_or_tail,
    )


def compute_local_duplication(
    gene_location: GeneLoc,
    mcl_cluster_data: MclClusterData,
    vBiomG_G: Dict[Any, Any],
    vKbRangeLD: int,
    vGeneRangeLD: int,
) -> pl.DataFrame:
    """
    ---------------------------------------------------------------------------

      Function that computes wheather a gene is locally duplicated, and if so,
      which MCL cluster it is associated with.

    ---------------------------------------------------------------------------

    Args:
        gene_location: DataFrame about the gene location information,
                        generated in PlantClusterFinder function.
        mcl_cluster_data : Mapping between the genes in the conversion-ID file and the
                MCL-clustering file.
        vBiomG_G: Mapping between the Gene location file and the gene-IDs in the
                  conversion file.
        vKbRangeLD: Genes that are encoded more than this bp apart, are not
                    called local duplicates.
        vGeneRangeLD: Genes that have more than this many intervening genes
                      between themare not called local duplicates.

    Returns:
        DataFrame with the columns:

        - GeneID: ID of Gene in GeneLocation2
        - is_locally_duplicated: A vector of length of the number of Genes in vGeneLocation2,
               being True if the gene was considered to be locally duplicated (a
               gene was within vKbRangeLD bp distance and not spaced more than
               with vGeneRangeLD genes from another gene with the same MCL
               clustering given by mcl_cluster_data).
        - ClusterIds: The MCL clusters of the gene (list)
    """

    # FIXME: This is glacially slow

    # For each gene G
    # Get genes which appear in the same cluster as key gene
    # mcl_df
    mcl_clusts = set([y for x in mcl_cluster_data.Gene2MCLClusters.values() for y in x])
    mcl_c2g = dict()
    for cluster in mcl_clusts:
        acc = []
        for gene, mcl in mcl_cluster_data.Gene2MCLClusters.items():
            if cluster in mcl:
                acc.append(gene)
        mcl_c2g[cluster] = acc
    mcl_g2g = dict()
    for gene, mcl in mcl_cluster_data.Gene2MCLClusters.items():
        acc = []
        for cluster in mcl:
            acc.extend(mcl_c2g[cluster])
        mcl_g2g[gene] = list(set(acc))
    mcl_df = pl.from_dicts(
        [
            {
                "Gene Name": k,
                "mcl_genes": [x for x in v if x != k],
                "mcl_clusters": mcl_cluster_data.Gene2MCLClusters.get(k, []),
            }
            for k, v in mcl_g2g.items()
        ]
    )

    # For each gene G
    # find genes G_i from -GeneRangeLD to GeneRangeLD (excluding G_0, which would just be G)
    # only consider genes within vKbRangeLd of G
    # nieghbor_genes
    prev_genes: list[pl.Row] = []  # rows from -GeneRangeLD to GeneRangeLD
    next_genes: list[pl.Row] = []
    current_gene = None
    # Collection of genes that surround the key gene
    neighbor_genes: Dict[str, list[str]] = dict()
    for row in gene_location.data.iter_rows():
        next_genes.append(row)
        if len(next_genes) > vGeneRangeLD:
            if current_gene is None:
                current_gene = next_genes.pop(-1)
            else:
                prev_genes.append(current_gene)
                current_gene = next_genes.pop(-1)
        if len(prev_genes) > vGeneRangeLD:
            prev_genes.pop(-1)
        if current_gene is None:
            continue
        prev_genes_filtered: filter[pl.Row] = filter(
            lambda x: (x[2] - current_gene[1]) < vKbRangeLD,
            prev_genes,  # prev_genes[2] is gene end
        )
        next_genes_filtered: filter[pl.Row] = filter(
            lambda x: current_gene[2] - x[1] < vKbRangeLD,
            next_genes,  # next_genes[1] is gene start
        )
        neighbor_genes[current_gene[0]] = list(prev_genes_filtered) + list(
            next_genes_filtered
        )
    neighbor_df = pl.from_dicts(
        [
            {"Gene Name": k, "neighbor_genes": [x[0] for x in v if x != k]}
            for k, v in neighbor_genes.items()
        ]
    )

    # Using the above datasets
    # see if mcl_df genes are in neighbor_genes
    # and mark G as locally duplicated if
    #  any G_i has the same MCL cluster as G
    locally_duplicated = neighbor_df.join(mcl_df, on="Gene Name").with_columns(
        locally_duplicated=pl.col("neighbor_genes")
        .list.set_intersection("mcl_genes")
        .list.len()
        > 0,
    )
    return gene_location.data.join(
        locally_duplicated, on="Gene Name", how="left"
    ).with_columns(pl.col("locally_duplicated").fill_null(False))


def introduce_hypo_genes(
    gene_location: GeneLoc,
    sequencing_gaps: pl.DataFrame,
    # hard_breaks: HardBreaks,
    vHypoAmount: int = 0,
    vGeneLengthHypoMin: int = -1,
) -> GeneLoc:
    """
    ---

      Function that introduces Hypothetical genes if a sequencing gap reached a
      critical size vGeneLengthHypoMin.

    ---
    Args:
        hard_breaks: Object with the gene location, chromosome regions, sequence gaps,
            and minimal length of bp of masked DNA that should be populated by two hypothetical
            genes
        vBiomG_G: Mapping between the Gene location file and the gene-IDs in the
                  conversion file.
        vGeneLengthHypoMin: A size that defines the minimal length of a
                            sequencing gap that should still be processed. If
                            this is bigger than 0, then all sequencing gaps
                            smaller than this are deleted, unless they are
                            among the list of vSeqGapSizesChromBreak, which
                            also needs to be set. Default is -1, not active.
    Optional Args:
        **vHypoAmount**: integer describing the number of hypothetical genes
            a sequencing gap that is not a physical break should be populated
            by artificial genes.
    Returns:
        **HardBreaks** object containing the gene location information, the
            inter-chromosomal matrix, and the HypoAmount populated.
    """
    # get the 5th and 7th columns of SeqGaps
    if sequencing_gaps is None:
        raise Exception(
            "Sequencing Gaps must be populated to insert hypothetical genes"
        )
    gap_positions = (
        sequencing_gaps.select(pl.col(["Chromosome Name", "Start", "End"]))
        .with_columns(seq_gaps_length=pl.col("End") - pl.col("Start"))
        .with_columns(
            hypo_to_place=(
                pl.lit(vHypoAmount)
                if (
                    vGeneLengthHypoMin < 0
                    or (pl.col("seq_gaps_length") >= pl.lit(vGeneLengthHypoMin))
                )
                else pl.lit(0)
            )
        )
    )
    vBreakAfter = []
    for pos in gap_positions.iter_rows():
        vBreakAfter.append(gene_location.find_intergenic_regions())

    # WARNING: This is a slow way to process dataframe
    ngaps = sequencing_gaps.height
    hypothetical_genes = []
    for i, row in enumerate(gap_positions.iter_rows()):
        log.debug("Populate hypothetical genes %i of %i", i, ngaps)
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
                "-1",  # FIXME: Strand - originally strand of start gene, but no clear reason why
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
        ),
        orient="row",
    )

    gene_location.insert_genes(hypo_gene_df)
    return gene_location


def annotate_sequencing_gaps(
    vFile: Path,
    vInterGeneLocation: pl.DataFrame,
    vUniqueChrom: List[str],
    vUnmaskedDNA: bool = True,
    vParallelize: int = 1,
    vOverwriteSeqGapInfo: bool = False,
    vSeqGapSizesChromBreak: List[int] = [],
    vGeneLengthHypoMin: int = -1,
) -> pl.DataFrame:
    """
    ---------------------------------------------------------------------------

      Functions that uses the output of get_sequencing_gaps and reformulates
      it such that it is clear in which intergenic region the sequening gap
      lies.

    ---------------------------------------------------------------------------

    Args:
        vFile: Fasta file of the masked genome nucleotide sequences.
        vInterGeneLocation: A DataFrame describing
                        - chromosome,
                        - start bp,
                        - end bp,
                        - gene before,
                        - and gene after
                      a intergenic region.
        vUniqueChrom: A list with unique entries of all chromosome and scaffold
                      names.
    Optional Args:
        - **vUnmaskedDNA**: Option to ignore the error if you use unmasked DNA files.
        - **vParallelize**: Number of parallel processes to run.
        - **vOverwriteSeqGapInfo**: If set to True, the code enforces masked nucleotide
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

    if not (vGeneLengthHypoMin == -1 or vSeqGapSizesChromBreak == []):
        log.error(
            "annotate_sequencing_gaps needs both SeqGapSizesChromBreak and GeneLengthHypoMin to be set or to be unset."
        )
        raise Exception("Incompatable parameters in annotate_sequencing_gaps function")

    gap_output_file = vFile.with_suffix(vFile.suffix + "_GAPOutput")
    if vOverwriteSeqGapInfo or not gap_output_file.exists():
        get_sequencing_gaps(vFile, 1, vOverwriteSeqGapInfo, vUnmaskedDNA)
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
        sequencing_gaps.filter(pl.col("scaffold").is_in(vUniqueChrom).not_())
        .select(pl.col("scaffold").unique())
        .iter_rows()
    ):
        log.warning(
            f"DNA Fasta file '{
                vFile
            }' contains chromosomes or scaffolds not covered with the "
            + f"gene position file: {bad_chromosome}"
        )
    sequencing_gaps = sequencing_gaps.filter(
        pl.col("scaffold").is_in(vUniqueChrom)
    )  # Do not process sequences that do not lie on a chromosome we want
    sequencing_gaps = (
        # All predicates of join_where() must be true for join to happen
        sequencing_gaps.join_where(
            vInterGeneLocation,
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
    _ = """
        vKillSeqGaps = false(vnvi,1); % Check which sequences should not be processed
        for vi =1:vnvi
            vChromID = find(strcmp(vUniqueChrom,vSeqGaps(vi,1)));
            vInterGeneLocationIDs = ...
                vInterGeneLocation(:,1) == vChromID & ...
                vInterGeneLocation(:,2) <= vSeqGaps_3(vi,1) & ...
                vInterGeneLocation(:,3) >= vSeqGaps_3(vi,1) & ...
                vInterGeneLocation(:,2) <= vSeqGaps_4(vi,1) & ...
                vInterGeneLocation(:,3) >= vSeqGaps_4(vi,1);
            if sum(vInterGeneLocationIDs) == 0
                vKillSeqGaps(vi,1) = 1;
            else
                vInterGeneLocationIDs_find = find(vInterGeneLocationIDs);
                vSeqGaps(vi,5) = cellstr(num2str(vInterGeneLocationIDs_find));
                vSeqGapsStart = max([vInterGeneLocation(
                    vInterGeneLocationIDs_find,2),str2num(char(vSeqGaps(vi,3)))]);
                vSeqGapsEnd = min([vInterGeneLocation(vInterGeneLocationIDs_find,3),str2num(
                    char(vSeqGaps(vi,3))) + str2num(char(vSeqGaps(vi,4))) - 1]);
                vSeqGaps(vi,6) = cellstr(num2str(vSeqGapsStart));
                vSeqGaps(vi,7) = cellstr(num2str(vSeqGapsEnd-vSeqGapsStart+1));
            end
        end
    end"""
    # FIXME: validate logic:
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

      Function that screens a DNA FASTA file for masked streches of DNA
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
            log.error("Cannot overwrite existing _GAPOutput with tag set")
            raise Exception("Overwrite Protection Enforcement of Existing Gap File")
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
    # TODO: Why does this file have \r EOL char?
    rxn_meta_dom = (
        pl.read_csv(vGenesMetabolicDomains_File, separator="\t", eol_char="\r")
        .unpivot(
            [f"Class{n + 1}" for n in range(6)],
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
    if vRemoveNonProtLocations:
        _protein_init = 300_000
    protein_ids = []
    with open(vProtein_sequence_FastaFile, "r") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                protein_ids.append(line[1:].strip())

    mcl_cluster_data = []
    clusters_by_p: Dict[str, List[int]] = dict()
    with open(vFile, "r") as f:
        for i, line in enumerate(f):
            log.debug(f"MCL_Line {i}")
            # if i <= init_max:
            prots_in_cluster = line.strip().split()
            mcl_cluster_data.append(prots_in_cluster)
            for prot in prots_in_cluster:
                clusters_by_p[prot] = clusters_by_p.get(prot, []) + [i]

    clusters_by_g = {
        conversion.MapP_MapG[p][0]: clust for p, clust in clusters_by_p.items()
    }
    proteins_list = list(clusters_by_p.keys())

    return MclClusterData(Gene2MCLClusters=clusters_by_g, Protein_IDs=proteins_list)


def get_MCL_clustering(
    vProtein_sequence_FastaFile: Path,
    conversion: Conversion,
    vOverwriteMCLClustering: bool,
    vRemoveNonProtLocations: bool,
    vParallelize: int = 0,
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

    Returns:
        MclClusterData object which provides mapping between
        genes in the conversion-ID file and the MCL-clustering file,
        as well as a list of protein-IDs stored in the fasta file.
    """
    int_cpu_count: int = os.cpu_count() or 1
    if vParallelize == 0:
        vParallelize = max(1, int_cpu_count - 1)

    apps = ["makeblastdb", "blastp", "cut", "mcxload", "mcl", "mcxdump"]
    missing_apps: List[str] = []
    app_path: Dict[str, Path] = dict()
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
    command_for["makeblastdb"] = [
        "-in",
        str(vProtein_sequence_FastaFile),
        "-dbtype",
        "prot",
        "-out",
        f"{vProtein_sequence_FastaFile}.blastDB",
        "-hash_index",
    ]
    command_for["blastp"] = [
        "-query",
        str(vProtein_sequence_FastaFile),
        "-db",
        f"{vProtein_sequence_FastaFile}.blastDB",
        "-evalue",
        "1",
        "-out",
        f"{vProtein_sequence_FastaFile}.blastOUT",
        "-outfmt",
        "6 qseqid sseqid evalue bitscore length pident qstart qend sstart send",
    ]
    if vParallelize > 1:
        command_for["blastp"].extend(["-num_threads", str(vParallelize)])
    command_for["cut"] = [
        "-f",
        "1,2,3",
        f"{vProtein_sequence_FastaFile}.blastOUT",
    ]
    cut_output = f"{vProtein_sequence_FastaFile}.MCLIN"
    command_for["mcxload"] = [
        "-abc",
        f"{vProtein_sequence_FastaFile}.MCLIN",
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
    overwrite_dump_workaround = (
        ""
        if vOverwriteMCLClustering
        else f"{
        os.getpid()}."
    )
    command_for["mcxdump"] = [
        "-icl",
        f"{vProtein_sequence_FastaFile}.MCLMATRIX.I20",
        "-tabr",
        f"{vProtein_sequence_FastaFile}.MCLTAB",
        "-o",
        f"{vProtein_sequence_FastaFile}.{
            overwrite_dump_workaround}MCLTAB.dump.I20",
    ]

    for app in apps:
        log.info(f"{app_path[app]} {' '.join(command_for[app])}")
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
            f"{vProtein_sequence_FastaFile}.{
                overwrite_dump_workaround}MCLTAB.dump.I20"
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
    log.info("Reading in genes.dat")
    genes_all = read_in_pgdb_flat(
        vPGDB_FlatFileFolder / "genes.dat",
        genes_attributes,
        unique_keys=genes_attributes,
    )
    log.info("Reading in proteins.dat")
    proteins_all = read_in_pgdb_flat(
        vPGDB_FlatFileFolder / "proteins.dat", proteins_attributes
    )
    log.info("Reading in enzrxns.dat")
    enzrxns_all = read_in_pgdb_flat(
        vPGDB_FlatFileFolder / "enzrxns.dat", enzrxns_attributes
    )
    log.info("Reading in reactions.dat")
    reactions_all = read_in_pgdb_flat(
        vPGDB_FlatFileFolder / "reactions.dat", reactions_attributes
    )

    # Match flat files to each other
    log.info("Matching genes to proteins")
    G2_P = match_pgdb_flat(
        genes_all, proteins_all, G2P, primary_id="GENE-ID", secondary_id="PROTIEN-ID"
    )
    log.info("Matching proteins to protein-complexes")
    P_P = match_pgdb_flat(
        proteins_all,
        proteins_all,
        P2P,
        primary_id="PROTIEN-ID",
        secondary_id="COMPONENTS",
    )

    # collect all proteins associated with a primary protein ID
    _ = """
    for i in range(len(P_P)):
        P_P[i, i] = 1
    sum_old = -1
    sum_new = sum(P_P)
    while sum_new != sum_old:
        sum_old = sum_new
        for i in range(len(P_P)):
            ids1 = P_P[i, :] == 1 or P_P[:, i].transpose == 1
            P_P[i, ids1] = 1
            P_P[ids1, i] = 1
        sum_new = sum(P_P)
    """

    log.info("Match proteins and protein complexes to enzymereactions")
    P_E = match_pgdb_flat(
        proteins_all,
        enzrxns_all,
        P2E,
        primary_id="PROTIEN-ID",
        secondary_id="ENZYME-ID",
    )
    log.info("Matching enzymereactions to reactions")
    E2_R = match_pgdb_flat(
        enzrxns_all,
        reactions_all,
        E2R,
        # R2_E = match_pgdb_flat(
        #    reactions_all,
        # enzrxns_all,
        #        R2E,
        primary_id="ENZYME-ID",
        secondary_id="REACTION-ID",
    )

    # Match Genes to reactions directly.
    log.info("Match genes to to reactions")
    # initialize matrix the size of [G2_P, E_R2]
    # gene2rxn = (
    gene_prot = G2_P.join(P_P, on="PROTIEN-ID", validate="1:1")
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
        Flat1: DataFrame of a read out pgdb Flat file.
        Flat2: DataFrame of a read out pgdb Flat file.
        join_on: Definition of which attributes in structure 1 should be
                    matched to which attributes of structure 2

    Returns:
        F1_F2: A mapping between items of the first pgdb file to the second
                one.
    """
    if join_on[0] not in flat1.columns:
        log.error(
            f"Matching column {join_on[0]} not in input data"
            + f" {flat1.columns}, note the other table's columns are {flat2.columns}"
        )
        raise Exception("Column-matching Mismatch")
    if join_on[1] not in flat2.columns:
        log.error(
            f"Matching column {join_on[1]} not in input data"
            + f" {flat2.columns}, note the other table's columns are {flat1.columns}"
        )
        raise Exception("Column-matching Mismatch")
    if flat1[join_on[0]].dtype == pl.Null:
        log.error(f"Primary table join ID {
                  [join_on[0]]} is null dtype, cannot join")
        log.info(f"Columns of flats are are {flat1.columns}, {flat2.columns}")
        raise Exception("Column-matching Mismatch")
    if flat2[join_on[1]].dtype == pl.Null:
        log.warning(
            f"Cannot map {join_on[0]} to {join_on[1]}"
            + ": dtype is Null, returning Null for second column"
        )
        log.info(f"Columns of flats are are {flat1.columns}, {flat2.columns}")
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
    log.debug("SUCCESSFUL JOIN", joined_dbs)
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
    conversion_data = pl.read_csv(
        vGeneTranscriptProtein_mapping_File,
        separator="\t",
        schema_overrides={
            x: pl.String for x in ["Gene Name", "Transcript Name", "Protein Name"]
        },
    )
    column_names = conversion_data.columns
    conversion_len = len(conversion_data)
    log.debug(
        "Read in conversion data, %d rows, columns are %s", conversion_len, column_names
    )
    conversion = Conversion(df=conversion_data)
    # Validate conversion before output
    # TODO: are unclean protein ID files an unrecoverable error?
    if len(conversion.MapG_MapT) < len(conversion.MapG):
        log.warning("Some Transcripts match multiple genes in the ID Conversion file")
        transcript_multi = (
            conversion_data.group_by(column_names[1])
            .n_unique()
            .filter(pl.col(column_names[0]) > 1)[column_names[1]]
        )
        for transcript in transcript_multi:
            log.warning("Please check the transcript %s", transcript)
        log.error("Please clean the protein ID file")
    if len(conversion.MapG_MapP) < len(conversion.MapG):
        log.warning("Some Proteins match multiple genes in the ID Conversion file")
        prot_multi = (
            conversion_data.group_by(column_names[2])
            .n_unique()
            .filter(pl.col(column_names[0]) > 1)[column_names[2]]
        )
        for prot in prot_multi:
            log.warning("Please check the protein %s", prot)
        log.error("Please clean the protein ID file")
    log.info(
        "Number of Gene-ID to Transcript-ID mappings: %i", len(conversion.MapG_MapT)
    )
    log.info("Number of Gene-ID to Protein mappings: %i", len(conversion.MapG_MapP))
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
    vGenes_all: pl.DataFrame,
    vGenesDat_GeneInfo: List[str],
    vPGDBIdsToMap: str,
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
    if vPGDBIdsToMap not in "GTP":
        log.error("Invalid PGDBIdsToMap mapping value, valid values are G, T, and P")
        raise Exception("Invalid Mapping Value")

    gene_map = dict()
    if "G" in vPGDBIdsToMap:
        gene_map = (
            Conversion.df.join(
                vGenes_all, left_on="Gene Name", right_on=vGenesDat_GeneInfo[0]
            )
            .select(pl.col("Gene Name"), pl.col("UNIQUE-ID"))
            .rows_by_key(key="UNIQUE-ID")
        )

    transcript_map = dict()
    if "T" in vPGDBIdsToMap:
        transcript_map = (
            Conversion.df.join(
                vGenes_all, left_on="Transcript Name", right_on=vGenesDat_GeneInfo[0]
            )
            .select(pl.col("Transcript Name"), pl.col("UNIQUE-ID"))
            .rows_by_key(key="UNIQUE-ID")
        )

    protein_map = dict()
    if "P" in vPGDBIdsToMap:
        protein_map = (
            Conversion.df.join(
                vGenes_all, left_on="Protein Name", right_on=vGenesDat_GeneInfo[0]
            )
            .select(pl.col("Protein Name"), pl.col("UNIQUE-ID"))
            .rows_by_key(key="UNIQUE-ID")
        )

    ret_keys = (
        set(gene_map.keys()).union(protein_map.keys()).union(transcript_map.keys())
    )
    # TODO: is there a way to just pass dataframes to result?
    return {
        x: list(gene_map.get(x, []))
        + list(transcript_map.get(x, []))
        + list(protein_map.get(x, []))
        for x in ret_keys
    }


def perform_clustering(
    vGeneLocation: GeneLoc,
    mcl_cluster_data: MclClusterData,  # covered by genelocation
    # vBiomG_E: Dict[Any, Any],# covered by geneLocation
    vEnzMinForClusters: int,
    vKbRangeLD: int,
    vGeneRangeLD: int,
    # vBiomG_G: Dict[Any, Any], #covered by genelocation
    # genes_to_reactions: GenesToReactions, # covered by genelocation
    # vMapG_Gpgdb: Dict[Any, Any], #c covered by genelocation
    vCriterion: int,  # TODO: this should be an enum
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
        vGeneLocation: DataFrame about the gene location information,
                        generated in PlantClusterFinder function.
        vMinStepSize: Minimal intervening gene size that should be computed
                      (default is 5).
        vMaxStepSize: Maximal intervening gene size that should be computed
                      (default is 20).
        vBiomG_E: A dictionary describing which gene has at least one reaction
                  annotated in the pgdb (1, if no reaction 0) (removed?)
       vEnzMinForClusters: Defines the minimum number of enyzmes that need to
                            be present in a cluster.
        vKbRangeLD: Genes that are encoded more than this bp apart, are not
                    called local duplicates.
        vGeneRangeLD: Number of genes that space two genes from each other.
                      Genes that have more than this many intervening genes
                      between themare not called local duplicates.
        mcl_cluster_data: Mapping between the genes in the conversion-ID file and the
                MCL-clustering file.
        vBiomG_G: Mapping between the Gene location file and the gene-IDs in the
                  conversion file.
        Gene2Rxn: Mapping of all genes in the Genes_all to the Reactions in
              Reactions_all
        vMapG_Gpgdb: A mapping between Genes in the conversion file and genes
                     stored in the pgdb genes.dat file.
        vCriterion: Criterion after which stepsize is chosen
                    1: as published in Schlapfer et al 2017
                       sum(enzymes in clusters) - sum(non-enzymes in clusters)
                       > 0.
                    2: stepsize 3 (because most organisms take this one)
                    3: sum(enzymes in clusters) - sum(non-enzymes in clusters,
                       without hypo genes) -2*sum(hypogenes)
                    4: sum(clusters in stepsize n) > sum(clusters in stepsize
                       n+1).

    Returns:
        Cluster object which contains the following information:
        - **Clusters**: Cluster Information (stepsize, chromosome id, gene start,
                   gene end, number of enzymes)
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

    continue_search = True  # flag if search should continue
    clusters = []  # array of cluster identifiers
    cluster_tuple = namedtuple(
        "cluster_tuple",
        [
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
    step_size = 0  # step size for search
    # note: hypothetical gene names start with "Hypo", this was originally tracked in a variable

    # used in place of biomG_E
    is_enzyme_criteria = pl.col("has_reaction")

    # get rid of clusters that only contain Local Duplicates

    def ld_is_not_a_problem(geneLocation: GeneLoc, start: int, end: int) -> bool:
        """
        Test if genes between start and end index are all local duplicates.
        Local duplicates must be within a provided range (vKbRangeLD) and
        number of genes (GeneRangeLD).
        """
        if end is None:
            return False
        cluster_genes = geneLocation.data.with_row_index().filter(
            (pl.col("index") > start) & (pl.col("index") <= end)
        )  # [start:end]
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
            cluster_genes.filter(
                pl.col("Gene Name").str.starts_with("Hpythoteical_gene_").not_()
            ).height
            >= vGeneRangeLD + 2
        ):
            return True
        # LD is not a problem if there are two local duplication groups
        # note: this test is only triggered for clusters of less than vKbRangeLD and vGeneRangeLD size
        # TODO: Validate this column, logic

        mcl_mapper = pl.DataFrame(
            {
                "Gene Name": list(mcl_cluster_data.Gene2MCLClusters.keys()),
                "mcl_clusters": list(mcl_cluster_data.Gene2MCLClusters.values()),
            }
        )
        cluster_genes = cluster_genes.join(mcl_mapper, on="Gene Name", how="left")
        unique_mcl_clusters: int = (
            cluster_genes.group_by("mcl_clusters").n_unique().height
        )
        if unique_mcl_clusters > 1:
            return True
        return False

    while step_size <= vMinStepSize or continue_search:  # for each valid step size
        for i in range(len(vGeneLocation)):  # for each gene
            log.debug("Potential Gene Start %i of %i", i, len(vGeneLocation))
            # index of next enzyme in vGeneLocation.data[i:]
            start: int = (
                vGeneLocation.data.with_row_index()
                .filter((pl.col("index") >= i) & (pl.col("has_reaction")))
                .select(pl.first("index"))["index"][0]
            )
            end: int = start
            abolish_cluster = False  # flag if the cluster should be followed
            if start is None:
                break  # no more enzymes
            # chromosome at start
            chrom_start = vGeneLocation.data["Chromosome Name"][start]
            # all genes in chrom_start chromosome
            chromosome_criteria = pl.col("Chromosome Name") == chrom_start

            # find the end
            stop_search_end = False
            while not stop_search_end:
                if end is None:
                    log.error("type of end is none")
                    raise Exception("Cluster search end was type None")
                end_s = end + 1  # stop search
                # stop search step size
                end_e = min(end + 1 + step_size, len(vGeneLocation.data))
                end_df = vGeneLocation.data.with_row_index().filter(
                    is_enzyme_criteria
                    & chromosome_criteria
                    & (pl.col("index") > end_s)
                    & (pl.col("index") <= end_e)
                )
                if end_df.height > 0:
                    end = end_df.select(pl.last("index"))["index"][0]
                else:
                    stop_search_end = True

            if not abolish_cluster:
                # get rid of clusters that only have one reaction set
                if (
                    len(vGeneLocation.data[start:end].filter(is_enzyme_criteria))
                    < vEnzMinForClusters
                ):
                    abolish_cluster = True
            if not abolish_cluster and not ld_is_not_a_problem(
                vGeneLocation, start, end
            ):
                # get rid of clusters that only have one local duplicate
                # criteria for this is outlined in the ld_is_not_a_problem function
                abolish_cluster = True

            if not abolish_cluster:
                # This is the code that writes the cluster into the clusters dataframe
                # what goes in here is described in line 3426 in matlab, but it is flexible
                # based on how you want to use it

                _cluster = vGeneLocation.data[start:end]

                clusters.append(
                    cluster_tuple(
                        step_size,
                        start,
                        end,
                        _cluster.filter(pl.col("has_reaction")).height,  # enzyme_count
                        _cluster.filter(
                            pl.col("has_reaction").not_()
                            & pl.col("Gene Name")
                            .str.starts_with("Hypothetical_Gene_")
                            .not_()
                        ).height,  # non_enzyme_gene_count
                        _cluster.filter(
                            pl.col("Gene Name")
                            .str.starts_with("Hypothetical_Gene_")
                            .not_()
                        ).height,  # hypothetical_gene_count
                    )
                )

        if len(clusters) > 0:
            clusters_df = pl.DataFrame(clusters)
            # re-calculate clusters_background to test against criterion
            clusters_background = compute_background(
                clusters_df, vGeneLocation, clusters_background
            )
            match vCriterion:
                case 1:
                    # (Schlapfer et al. 2017 criteria)
                    # sum (enzymes in clusters) - sum(non-enzymes in clusters (incl. hypo genes?)) > 0
                    # TODO: validate criterion
                    log.debug(
                        "Number of enzymes - number of non-enzymes in clustered genes: %i",
                        clusters_df.filter(pl.col("step_size") == step_size)
                        .with_columns(pl.col("enzyme_count"))
                        .sum(),
                    )
                    if step_size == vMaxStepSize:
                        continue_search = False
                        final_step_size = step_size
                        log.warning(
                            "Maximal step size reached: %i. Consider increasing.",
                            vMaxStepSize,
                        )
                    elif (
                        clusters_df.filter(pl.col("step_size") == step_size).height == 0
                    ):
                        log.warning("No clusters found for step size %i.", step_size)
                    # FIXME: this criterion is only for clusters in the step_size
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
                    log.debug(
                        "Number of enzymes - number of non-enzymes (no hypos) - 2* number if hypos in clustered genes: %i",
                        clusters_df["enzyme_count"].count()
                        - clusters_df["non_enzyme_gene_count"].count(),
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
                        log.warning(
                            "Maximal step size reached: %i. Consider increasing.",
                            vMaxStepSize,
                        )
                    elif (
                        clusters_df.filter(pl.col("step_size") == step_size).height == 0
                    ):
                        log.warning("No clusters found for step size %i.", step_size)
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
                        # FIXME: Validate this logic
                        if continue_search:
                            final_step_size = step_size - 1
                        continue_search = False
                case 4:
                    # sum(count of clusters in stepsize n) > sum(count of clusters in stepsize n + 1)
                    if step_size != 0:
                        log.debug(
                            "Number of clusters more than last stepsize %i",
                            clusters_df.filter(
                                pl.col("step_size") > (step_size - 1)
                            ).height,
                        )
                    else:
                        log.debug(
                            "Number of clusters Stopsize 0: %i",
                            clusters_df.filter(pl.col("step_size") == 0).height,
                        )
                    if step_size == vMaxStepSize:
                        continue_search = False
                        final_step_size = step_size
                        log.warning(
                            "Maximal step size reached: %i. Consider increasing.",
                            vMaxStepSize,
                        )
                    elif (
                        clusters_df.filter(pl.col("step_size") == step_size).height == 0
                    ):
                        log.warning("No clusters found for step size %i.", step_size)
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
        step_size += 1
    return Clusters(
        ClusterMeta=pl.DataFrame(clusters),
        ClusterBackground=clusters_background,
        # FIXME: there should be no path to final_step_size being unbound
        FinalStepSize=final_step_size or step_size - 1,
        MaxStepSize=step_size - 1,
    )


def compute_background(
    clusters: pl.DataFrame,
    vGeneLocation: GeneLoc,
    clusters_background: pl.DataFrame,
    # vC_IDs: List[str],
) -> pl.DataFrame:
    """
    ---------------------------------------------------------------------------

      Function that lookes at all clusters that were found by the function
      f_perform_clustering(), and computes the background information for the
      window sizes of these clusters (e.g. a cluster has 4 enzymes and two non
      enzymes, then its size is 6. The function goes with a window size 6
      through the genome and checks how many windows exist with 0, 1, 2, ...,
      5, 6 enzymes. It then calculates the probability to see a cluster with
      more or equal than 0 enzymes, more or equal than 1, ... more or equal
      than 6 enzymes. This is then used to define top n% clusters, e.g.
      clusters that are at least as dense as the top 5% of windowsizes are
      marked to be of higher interest.

    ---------------------------------------------------------------------------

    Args:
        vClusters: Cluster Information (stepsize, chromosome id, gene start,
                   gene end, number of enzymes, number of non-enzyme genes,
                   number of hypothetical genes)
        vGeneLocation: DataFrame about the gene location information,
                        generated in PlantClusterFinder function.
        clusters_background: DataFrame of backgrounds for up to vClusters
                max_step_size - 1. (window_size, enzyme_count,
                                    window_occurrence, winzyme_occurrence, winzyme_frequency)
                "Winzyme" is a window with enzyme_count enzymes in it.
        vC_IDs: Defines which Entries in vClusters are Clusters of the current
                Stepsize that was calculated and should be checked if the
                background of that size of cluster (number of genes) was
                already calculated.

    Returns:
        - **vClusterBackground**: Density information about the genome (window size,
                            number of enzymes, Occurrence of windowsize with
                            this many enzymes in the genome, occurrence of
                            this windowsize in the genome, Frequency of seeing
                            a windowsize of this size with at least as many
                            enzymes in the genome)
    """
    # (for a given) window size
    # (for a given) number of enzymes
    # occurrence of window size in genome
    # frequency of seeing window size with at least this many genes

    index_offset = len(vGeneLocation)
    chromosome_offsets = pl.DataFrame(
        [
            {
                "Chromosome Name": x["Chromosome Name"],
                "index_offset": x["chr_idx"] * index_offset,
            }
            for x in vGeneLocation.data.select(pl.col("Chromosome Name").unique())
            .with_row_index("chr_idx")
            .to_dicts()
        ]
    )
    indexed_gene_location = (
        vGeneLocation.data.with_row_index("idx")
        .join(chromosome_offsets, on="Chromosome Name")
        .with_columns(idx=pl.col("idx") + pl.col("index_offset"))
        .sort(pl.col("idx"))
    )
    # "winzyme" is window x enzyme
    max_step_size = clusters.select(pl.col("step_size").max())["step_size"][0]
    if clusters_background.is_empty():
        clusters_background = pl.DataFrame(
            [],
            schema=[
                ("window_size", pl.UInt32),
                ("enzyme_count", pl.UInt32),
                ("winzyme_occurrence", pl.Int32),
                ("window_occurrence", pl.Int32),
                ("winzyme_frequency", pl.Float64),
            ],
        )
        prev_max_step_size = 1
    else:
        prev_max_step_size = clusters_background.select(pl.max("window_size"))[
            "window_size"
        ][0]
    for window_size in range(prev_max_step_size, max_step_size + 1):
        background_df_win = (
            indexed_gene_location.select(
                # The first (window_size-1) rows will not be window_size windows
                window_size=pl.len().rolling(
                    index_column="idx", period=f"{window_size}i"
                ),
                enzyme_count=pl.sum("has_reaction").rolling(
                    index_column="idx", period=f"{window_size}i"
                ),
                window_occurrence=1,
            )
            .filter(pl.col("window_size") == window_size)
            .group_by("window_size", "enzyme_count")
            .agg(winzyme_occurrence=pl.sum("window_occurrence"))
        )
        total_window_occurrence = len(background_df_win)
        clusters_background.extend(
            background_df_win.with_columns(
                window_occurrence=total_window_occurrence,
                winzyme_frequency=pl.col("winzyme_occurrence")
                / total_window_occurrence,
            )
        )

    return clusters_background


def perform_cutoff(clusters: Clusters, vTopPercentClusters: int) -> List[bool]:
    """
    ---------------------------------------------------------------------------

      Function that (given clusters and a background density information of a
      genome) calculates whether the clusters are within the top n% of enzymes
      dense clusters.

    ---------------------------------------------------------------------------

    Args:
        clusters: Cluster Information (stepsize, chromosome id, gene start,
                   gene end, number of enzymes) and
                   Density information about the genome (window size,
                   number of enzymes, Occurrence of windowsize with
                   this many enzymes in the genome, occurence of
                   this windowsize in the genome, Frequency of seeing
                   a windowsize of this size with at least as many
                   enzymes in the genome)
        vTopPercentClusters: Clusters are labeled as being in the top n% if
                             they are within this top n% of enzyme dense
                             clusters compared to the background.

    Returns:
        List with size Number of Clusters. It is zero if
               the density of the cluster is above or the same as given
               by the input argument vTopPercentClusters*100.
    """

    _ = """
        % MATLAB code
        vClusters_Top = false(size(vClusters,1),1);
        for vj = 1:size(vClusters,1)
            vSize = vClusters(vj,3)-vClusters(vj,2) + 1;
            vEnzymeSize = vClusters(vj,4);
            vTempBackground_ID = vClusterBackground(:,1) == vSize & vClusterBackground(:,2) == vEnzymeSize;
            if vClusterBackground(vTempBackground_ID,5) <= vTopPercentClusters/100
                vClusters_Top(vj,1) = 1;
            end
        end
    end"""
    # TODO: This logic should be able to push a modified cluster object
    cluster_meta_in_top = clusters.ClusterMeta.join(
        clusters.ClusterBackground,
        left_on=["step_size", "enzyme_count"],
        right_on=["window_size", "enzyme_count"],
    ).select(in_top=pl.col("winzyme_frequency") >= vTopPercentClusters / 100)
    return cluster_meta_in_top["in_top"].to_list()


def write_gene_output(
    vGeneOutputFile: Path,
    vGeneLocation2: pl.DataFrame,
    clusters: Clusters,
    vTopPercentClusters: List[float],
    vMD_to_annotate: List[str],
    Gene2Rxn: GenesToReactions,
    vMapG_Gpgdb: Dict[Any, Any],
    vBiomG_G: Dict[Any, Any],
    vBiomG_E: Dict[Any, Any],
    local_duplications: pl.DataFrame,
    vClusters_Top: List[bool],
    annotated_metabolic_domains: AnnotatedMetabolicDomains,
    vRpgdb_Rmd: Dict[Any, Any],
    reaction_labels: ReactionLabels,
    vOutputType: str,
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
        vGeneOutputFile: Path and file name to the gene-Outputfile, make sure
                         you have access to the file and folder. This is the
                         first Results output file that will be generated.
        vGeneLocation2: DataFrame about the gene location information,
                        generated in PlantClusterFinder function.
        clusters:
        vTopPercentClusters: Clusters are labeled as being in the top n% if
                             they are within this top n% of enzyme dense
                             clusters compared to the background (default is
                             5).
        vMD_to_annotate: A list of Metabolic domains that should be analyzed.
                         For example {'Amines and Polyamines Metabolism'; ...
                         'Amino Acids Metabolism'
        Gene2Rxn: Mapping of all genes in the Genes_all to the Reactions in
              Reactions_all
        vMapG_Gpgdb: A mapping between Genes in the conversion file and genes
                     stored in the pgdb genes.dat file.
        vBiomG_G: Mapping between the Gene location file and the gene-IDs in the
                  conversion file.
        vReactions_all: DataFrame containing all the information
                        (Attributes)that was searched for (described in
                        vReactionsAttributes, default is {'UNIQUE-ID', ...
                        'ENZYMATIC-REACTION', 'EC-NUMBER'};) in reactions.dat
                        of the pgdb.
        vReactions_all_Att: Attributes of vReactions_all, can be displayed by
                            vReactions_all.(char(vReactions_all_Att(1,n)))
        vBiomG_E: A vector describing which gene has at least one reaction
                  annotated in the pgdb (1, if no reaction 0)
        vG_LD_ClustIDs: The MCL clusters of the gene (comma separated if more
                        than 1 peptide is associated to a gene
        vG_LD: A vector of length of the number of Genes in vGeneLocation2,
               being 1 if the gene was considered to be locally duplicated (a
               gene was within vKbRangeLD bp distance and not spaced more than
               with vGeneRangeLD genes from another gene with the same MCL
               clustering given by vG_MCL.
        vClusters: Cluster Information (stepsize, chromosome id, gene start,
                   gene end, number of enzymes)
        vClusters_Top: Vector with size Number of Clusters x 1. It is zero if
                       the density of the cluster is above or the same as given
                       by the input argument vTopPercentClusters*100.
        vR_MD: Mapping between all Reactions stored in the Metabolic domain file
               and the metabolic domains that were identified to be of interest
               in vMD_to_annotate.
        vRpgdb_Rmd: A mapping between reactions stored in the pgdb and the
                    reactions given in the metabolic domain file
        vSTClasses: Is a list of signature and tailoring enzyme classes.
        vR_STClasses: Mapping of Reactions in vReactions_all to the enzyme
                      classes of vSTClasses.
        vSigOrTail: Is a list of types of enzymes, e.g. signature, tailoring
                    and so on.
        vR_SigOrTail: Mapping of Reactions in vReactions_all to the enzyme
                      types of vSigOrTail.
        vOutputType: Defines which type of output is given: simple (default
                     setting) gives only top n% for the final stepsize. verbose
                     gives all computed stepsize until the final stepsize. old
                     gives at least 5 until the final stepsize.

    Returns:
        **None** (File is generated with the path of the first input argument)
    """
    with vGeneOutputFile.open(mode="w") as fio:
        writer = csv.writer(fio, delimiter="\t")
        # write header
        header = vGeneLocation2.columns[0:5] + [
            "Enzyme",
            "EC",
            "RXN",
            "Locally Duplicated",
            "LD_ClusterIds",
        ]
        if vOutputType not in ["old", "simple", "verbose"]:
            raise Exception(
                "Output type not supported, options are 'old' 'simple' or 'verbose'"
            )
        elif vOutputType == "old":
            for i in range(clusters.MaxStepSize):
                header += [
                    f"{'CHOSEN ' if i == clusters.FinalStepSize else ''}Stepsize {i}",
                    f"{'CHOSEN ' if i == clusters.FinalStepSize else ''}Stepsize {
                        i
                    }, Top {vTopPercentClusters}%",
                ]
        elif vOutputType == "simple":
            header += [f"Stepsize {clusters.FinalStepSize}, Top {vTopPercentClusters}%"]
        elif vOutputType == "verbose":
            for i in range(clusters.MaxStepSize):
                header += [
                    f"Stepsize {i}",
                    f"Stepsize {i}, Top {vTopPercentClusters}%",
                ]
        # metabolic domain info
        header += vMD_to_annotate
        header += ["Signature Tailoring Enzyme Class", "Signature or Tailoring"]

        writer.writerow(header)

        # TODO: write cluster data

    """
        % MATLAB code
        %Gene Location
        fprintf(vFIO, '%s\t%s\t%s\t%s\t%s', vGeneLocation_Att{1,1:5});
        fprintf(vFIO, '\tEnzyme\tEC\tRXN\tLocally Duplicated\tLD_ClusterIds');
        for vi = 0:vFinalStepsize
            fprintf(vFIO, '\tStepsize %i\tStepsize %i, Top %i%%',
                    vi, vi, vTopPercentClusters);
        for vi = 1:size(vMD_to_annotate,1)
            fprintf(vFIO, '\t%s', char(vMD_to_annotate(vi,1)));
        end
        %Signature Tailoring info
        fprintf(vFIO, '\tSignature Tailoring Enzyme Class\tSignature or Tailoring');
        fprintf(vFIO,'\n');

        %% Print Data
        vNvi = size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1);
        for vi = 1:vNvi %Print all genes
            if vVerbose >= 2
                fprintf('Write gene %i of %i\n', vi, vNvi);
            end
            %Print Gene Location
            fprintf(vFIO, '%s', char(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vi,1)));
            fprintf(vFIO, '\t%s', char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(vi,1)));
            fprintf(vFIO, '\t%s', char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vi,1)));
            fprintf(vFIO, '\t%s', char(vGeneLocation2.(char(vGeneLocation_Att(1,6)))(vi,1)));
            fprintf(vFIO, '\t%s', char(vGeneLocation2.(char(vGeneLocation_Att(1,5)))(vi,1)));

            %Compose Gene information (reactions, enzymes, local duplications...)
            vRXN_IDs = sum(Gene2Rxn(sum(vMapG_Gpgdb(vBiomG_G(vi,:),:),1)>0,:),1)>0;
            vRXN_IDs_find = find(vRXN_IDs);
            vECstring = '';
            vRXNstring = '';
            if size(vRXN_IDs_find,2)>0
                vRXNstring = char(
                    vReactions_all.UNIQUE_ID(vRXN_IDs_find(1,1),1));
                if ~isempty(vReactions_all.(char(vReactions_all_Att(1,3))){vRXN_IDs_find(1,1),1})
                    vECstring = char(vReactions_all.(char(vReactions_all_Att(1,3)))(vRXN_IDs_find(1,1),1));
                end
                for vk = 2:size(vReactions_all.(char(vReactions_all_Att(1,3))),2)
                    if ~isempty(vReactions_all.(char(vReactions_all_Att(1,3))){vRXN_IDs_find(1,1),vk})
                        vECstring = [vECstring ', ' char(vReactions_all.(char(vReactions_all_Att(1,3)))(vRXN_IDs_find(1,1),vk))];
                    end
                end
                for vj = 2:size(vRXN_IDs_find,2)
                    vRXNstring = [vRXNstring ', ' char(vReactions_all.UNIQUE_ID(vRXN_IDs_find(1,vj)))];
                    vECstring = [vECstring ', ' char(vReactions_all.(char(vReactions_all_Att(1,3)))(vRXN_IDs_find(1,vj),1))];
                    for vk = 2:size(vReactions_all.(char(vReactions_all_Att(1,3))),2)
                        if ~isempty(vReactions_all.(char(vReactions_all_Att(1,3))){vRXN_IDs_find(1,vj),vk})
                            vECstring = [vECstring ', ' char(vReactions_all.(char(vReactions_all_Att(1,3)))(vRXN_IDs_find(1,vj),vk))];
                        end
                    end
                end
            end
            vOld = size(vRXNstring,2)+1;
            vNew = size(vRXNstring,2);
            while vOld > vNew
                vOld = size(vRXNstring,2);
                vRXNstring = regexprep(vRXNstring, ', , ', ', ');
                vNew = size(vRXNstring,2);
            end
            vOld = size(vECstring,2)+1;
            vNew = size(vECstring,2);
            while vOld > vNew
                vOld = size(vECstring,2);
                vECstring = regexprep(vECstring, ', , ', ', ');
                vNew = size(vECstring,2);
            end
            vRXNstring = regexprep(vRXNstring, ', $', '');
            vRXNstring = regexprep(vRXNstring, '^, ', '');
            vECstring = regexprep(vECstring, ', $', '');
            vECstring = regexprep(vECstring, '^, ', '');
            % Print gene information.
            fprintf(vFIO, '\t%i\t%s\t%s', vBiomG_E(
                vi,1), vECstring, vRXNstring);
            if isempty(vG_LD_ClustIDs{vi,1})
                vG_LD_ClustIDs(vi,1) = cellstr('');
            end
            fprintf(vFIO, '\t%i\t%s', vG_LD(vi,1), char(vG_LD_ClustIDs(vi,1)));
            % Print cluster information
            if strcmp(vOutputType,'old')
                for vj = 0:vStepSize_Max
                    vClustIDs = (vClusters(:,1) == vj) & (vClusters(:,2) <= vi) & (vi <= vClusters(:,3));
                    if sum(vClustIDs)>0
                        fprintf(vFIO, '\t%s\t%i', ['C' num2str(find(vClustIDs))], vClusters_Top(vClustIDs));
                    else
                        fprintf(vFIO, '\t\t');
                    end
                end
            elseif strcmp(vOutputType,'simple')
                vClustIDs = (vClusters(:,1) == vFinalStepsize) & (vClusters(:,2) <= vi) & (vi <= vClusters(:,3));
                if sum(vClustIDs)>0
                    if vClusters_Top(vClustIDs) == 1
                        fprintf(vFIO, '\t%s', ['C' num2str(find(vClustIDs))]);
                    else
                        fprintf(vFIO, '\t');
                    end
                else
                    fprintf(vFIO, '\t');
                end
            elseif strcmp(vOutputType,'verbose')
                for vj = 0:vFinalStepsize
                    vClustIDs = (vClusters(:,1) == vj) & (vClusters(:,2) <= vi) & (vi <= vClusters(:,3));
                    if sum(vClustIDs)>0
                        fprintf(vFIO, '\t%s\t%i', ['C' num2str(find(vClustIDs))], vClusters_Top(vClustIDs));
                    else
                        fprintf(vFIO, '\t\t');
                    end
                end
            end
            vMD_IDs = sum(vR_MD(sum(vRpgdb_Rmd(vRXN_IDs,:),1)>0,:),1)>0;
            for vj = 1:size(vMD_IDs,2)
                fprintf(vFIO, '\t%i',vMD_IDs(1,vj));
            end
            vSTClassesString = '';
            vSTClasses_find = find(sum(vR_STClasses(vRXN_IDs,:),1)>0);
            if size(vSTClasses_find,2)>0
                vSTClassesString = char(vSTClasses(vSTClasses_find(1,1)));
                for vj = 2:size(vSTClasses_find,2)
                    vSTClassesString = [vSTClassesString ', ' char(vSTClasses(vSTClasses_find(1,vj)))];
                end
            end
            vSigOrTailString = '';

            vSIGCol = strcmpi(vSigOrTail,'signature');
            vSIGTailCOls = sum(vR_SigOrTail(vRXN_IDs,:),1)>0;
            if sum(vSIGTailCOls(1,vSIGCol))>0
                vSigOrTailString = 'signature';
            elseif sum(vSIGTailCOls(1,:))>0
                vSigOrTailString = 'tailoring';
            end
            fprintf(vFIO, '\t%s\t%s', vSTClassesString, vSigOrTailString);
            fprintf(vFIO, '\n');
        end
        fclose(vFIO);
    end"""


def write_cluster_output(
    vClusterOutputFile: Path,
    vGeneLocation: pl.DataFrame,
    clusters: Clusters,
    vClusters_Top: List[bool],
    vOutputType: str,
) -> None:
    if vOutputType not in ["old", "simple", "verbose"]:
        raise Exception(
            "Output type not supported, options are 'old' 'simple' or 'verbose'"
        )
    """
    ---------------------------------------------------------------------------

      Function that writes the clusterfile of PlantClusterFinder: TSV format,
      first it includes the stepsize of the cluster (intervening genes), then
      cluster ID, chromosome/scaffold the cluster is located on, start bp, end
      bp, physical size, start gene ID, end gene ID, and all the gene IDs

    ---------------------------------------------------------------------------

    Args:
        vClusterOutputFile: Path and file name to the cluster-Outputfile, make
                            sure you have access to the file and folder. This
                            is the second Results output file that will be
                            generated.
        vClusters: Cluster Information (stepsize, chromosome id, gene start,
                   gene end, number of enzymes)
        vClusters_Top: List with size Number of Clusters x 1. It is zero if
                       the density of the cluster is above or the same as given
                       by the input argument vTopPercentClusters*100.
        vGeneLocation: DataFrame about the gene location information,
                        generated in PlantClusterFinder function.
        vOutputType: Defines which type of output is given: simple (default
                     setting) gives only top n% for the final stepsize. verbose
                     gives all computed stepsize until the final stepsize. old
                     gives at least 5 until the final stepsize.

    Returns:
        **None** (File is generated with the path of the first input argument)
    """

    """
        % MATLAB code
        %Check for output type:
        if ~strcmp(vOutputType,'old') && ~strcmp(vOutputType,'verbose') && ~strcmp(vOutputType,'simple')
            error('OutputType has to either be simple or verbose')
        end
        vFIO = fopen(vClusterOutputFile, 'w');

        %Print Header
        fprintf(vFIO, 'Stepsize\tClusterID\tTop n%% cluster\tChromosome\tStart bp\tEnd bp\tPhysical size bp\tNumber of genes\tNumber of Hypothetical genes\tnumber of Enzymes\tStart Gene\tEnd Gene\tGene IDs\n');

        %Print Data
        vNvi = size(vClusters,1);
        for vi = 1:vNvi
            if strcmp(vOutputType,'old') || ...
                (strcmp(vOutputType,'verbose') && vClusters(vi,1) <= vFinalStepsize) || ...
                (strcmp(vOutputType,'simple') && vClusters(vi,1) == vFinalStepsize && vClusters_Top(vi) == 1)
                if vVerbose >= 2
                    fprintf('Write cluster %i of %i\n', vi, vNvi);
                end
                fprintf(vFIO, '%i\t%s\t%i\t%s', vClusters(vi,1), ['C' num2str(vi)], vClusters_Top(vi,1), char(vGeneLocation2.(char(vGeneLocation_Att(1,6)))(vClusters(vi,2),1)));
                vStartbp = min([str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(vClusters(vi,2):vClusters(vi,3),1)));str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vClusters(vi,2):vClusters(vi,3),1)))]);
                vEndbp = max([str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(vClusters(vi,2):vClusters(vi,3),1)));str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vClusters(vi,2):vClusters(vi,3),1)))]);
                fprintf(vFIO, '\t%i\t%i\t%i',vStartbp,
                        vEndbp, vEndbp-vStartbp+1);
                fprintf(vFIO, '\t%i\t%i\t%i',vClusters(vi,3)-vClusters(vi,2)+1, sum(strncmp(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vClusters(vi,2):vClusters(vi,3)),'Hypothetical_Gene_',18)), vClusters(vi,4));
                vClustersString = char(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vClusters(vi,2),1));
                for vj = vClusters(vi,2)+1:vClusters(vi,3)
                    vClustersString = [vClustersString ', ' char(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vj,1))];
                end
                fprintf(vFIO, '\t%s\t%s\t%s\n',char(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vClusters(vi,2),1)),char(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vClusters(vi,3),1)),vClustersString);
            end
        end
        fclose(vFIO);
    end"""


def analyze_plant_cluster_gap_file(
    vFileIn: Path, vUnmaskedDNA: bool, file_output: bool = True
) -> Optional[Tuple[str, Dict[Any, Iterable[Any]]]]:
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
    pgdb_filename: Path, attrib: tuple[str], unique_keys: tuple[str] = ("UNIQUE-ID",)
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
        Generate a row from the pgdb flat file
        """
        y: Dict[str, Any] = {att: None for att in attrib}
        for line in fh:
            if line.startswith("#"):
                continue
            elif line.startswith("//"):
                yield y
                y = {att: None for att in attrib}
                continue
            elif " - " in line:
                key, *value_arr = line.split(" - ")
                value: str = (" - ".join(value_arr)).strip()
                if key in y:
                    if y[key] is None:
                        if key in unique_keys:
                            y[key] = value
                        else:
                            y[key] = [value]
                    elif not isinstance(y[key], list):
                        y[key] = [y[key], value]
                    else:
                        y[key].append(value)

    with open(pgdb_filename, "r") as f:
        pgdb_data = pl.DataFrame(list(_row_gen(f, attrib, unique_keys)))

    return pgdb_data
