import csv
import logging
import os
import subprocess
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import polars as pl

from plantclusterfinder_py.cluster_datatypes import (
    AnnotatedMetabolicDomains,
    Clusters,
    Conversion,
    GenesToReactions,
    HardBreaks,
    MclClusterData,
    ReactionLabels,
)
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
    vGeneLocation: pl.DataFrame,
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
        vGeneLocation: DataFrame about the gene location information,
                        generated in PlantClusterFinder function.
        mcl_cluster_data : Mapping between the genes in the conversion-ID file and the
                MCL-clustering file.
        vBiomG_G: Mapping between the Gene location file and the gene-IDs in the
                  conversion file.
        vKbRangeLD: Genes that are encoded more than this bp apart, are not
                    called local duplicates.
        vGeneRangeLD: Number of genes that space two genes from each other.
                      Genes that have more than this many intervening genes
                      between themare not called local duplicates.

    Returns:
        DataFrame with the columns:

        - GeneID: ID of Gene in GeneLocation2
        - is_locally_duplicated: A vector of length of the number of Genes in vGeneLocation2,
               being True if the gene was considered to be locally duplicated (a
               gene was within vKbRangeLD bp distance and not spaced more than
               with vGeneRangeLD genes from another gene with the same MCL
               clustering given by vG_MCL.
        - ClusterIds: The MCL clusters of the gene (comma separated if more
                        than 1 peptide is associated to a gene
    """

    """
        % MATLAB code
        vNgenes = size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1);
        vG_LD = false(vNgenes,1);
        vG_LD_ClustIDs = cell(vNgenes,1);
        vGStartbp = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))));
        vGEndbp = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))));
    """
    n_genes = len(vGeneLocation)
    for row in vGeneLocation.iter_rows(named=True):
        log.debug(
            f"Checking local duplication of Gene {row['gene name']} of {len(vGeneLocation)}"
        )
        if list(row.keys())[0].startswith("Hypothetical_Gene_"):
            _start_g = (
                row["gene start"] - vGeneRangeLD - 1
                if row["gene start"] < (vGeneRangeLD - 1)
                else 1
            )
            _end_g = (
                row["gene end"] + vGeneRangeLD + 1
                if (row["gene end"] + vGeneRangeLD + 1) < n_genes
                else n_genes
            )
        """
        for vi = 1:vNgenes
            if vVerbose >= 2
                fprintf(
                    'Checking Local duplication of Gene %i of %i\n', vi, vNgenes);
            end
            if ~strncmp(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vi,1),'Hypothetical_Gene_',18)
                vStartG = vi - vGeneRangeLD - 1;
                vEndG = vi + vGeneRangeLD + 1;
                if vStartG < 1
                    vStartG = 1;
                end
                if vEndG > vNgenes
                    vEndG = vNgenes;
                end
        """
        _ids = np.array((n_genes, 1))

    """
                vIDs = false(vNgenes,1);
                vIDs(vStartG:vEndG,1) = 1;
                vIDs(vi,1) = 0;
                vIDs = vIDs & strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,4))),vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi,1));

                vStartbp = vGEndbp(vi,1)-vKbRangeLD;
                vEndbp = vGStartbp(vi,1)+vKbRangeLD;
                vIDs = vIDs & vGStartbp >= vStartbp & vGEndbp <= vEndbp;
                vIDs_find = find(vIDs);

                vMCL = sum(vG_MCL(sum(vBiomG_G(vi,:),1)>0,:),1)>0;
                vMCL_find = find(vMCL);
                clear vTempMCLString_array
                if size(vMCL_find,2)~=0
                    vTempMCLString_array = cellstr(['MCL_' num2str(vMCL_find(1,1))]);
                end
                for vj = 2:size(vMCL_find,2)
                    vTempMCLString_array = [vTempMCLString_array; cellstr(['MCL_' num2str(vMCL_find(1,1))])];
                end
                if exist('vTempMCLString_array','var')
                    vTempMCLString_array = unique(vTempMCLString_array);
                end
                vTempMCLString = '';
                if exist('vTempMCLString_array','var')
                    if size(vTempMCLString_array,1)~=0
                        vTempMCLString = char(vTempMCLString_array);
                    end
                    for vj = 2:size(vTempMCLString_array,1)
                        vTempMCLString = [vTempMCLString ', ' char(vTempMCLString_array)];
                    end
                end
                vG_LD_ClustIDs(vi,1) = cellstr(vTempMCLString);
                vN_MCL = size(vMCL,2);

                vStop = 0;
                vj = 1;
                while vj <= size(vIDs_find,1) && vStop == 0
                    if sum(vMCL)==0
                        fprintf('Warning: %s has no entry in MCL clustering file.\n',char(vGeneLocation2.(char(vGeneLocation_Att(1,1)))(vi,1)))
                        vG_LD(vi,1) = -1;
                    elseif sum(vMCL == (sum(vG_MCL(sum(vBiomG_G(vIDs_find(vj,1),:),1)>0,:),1)>0)) == vN_MCL
                        vG_LD(vi,1) = 1;
                        vStop = 1;
                    end
                    vj = vj + 1;
                end
            end
        end
    end
    """
    return pl.DataFrame()  # FIXME: needs return value


def introduce_hypo_genes(
    hard_breaks: HardBreaks,
    vBiomG_G: Dict[Any, Any],
    insert_hypos: bool = False,
    vHypoAmount: int = 0,
) -> HardBreaks:
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
        vInsertHypos: If set to True, then Hypothetical genes are introduced.
    Optional Args:
        **vHypoAmount**: integer describing the number of hypothetical genes
            a sequencing gap that is not a physical break should be populated
            by artificial genes.
    Returns:
        **HardBreaks** object containing the gene location information, the
            inter-chromosomal matrix, and the HypoAmount populated.
    """
    if not insert_hypos:
        # why would you even call the function?!?
        return hard_breaks

    """
        % MATLAB code
        vCurrArgin = 8;
        if nargin > vCurrArgin
            vi = 1;
            while vi <= nargin - vCurrArgin
                if strcmp(varargin{1,vi},'HypoAmount')
                    if ischar(varargin{1,vi+1})
                        vHypoAmount = str2num(varargin{1,vi + 1});
                    else
                        vHypoAmount = varargin{1,vi + 1};
                    end
                    vi = vi + 1;
                end
                vi = vi + 1;
            end
        end

        if vInsertHypos > 0
            %Check where and how many Hypothetical genes should be placed.
            if vInsertHypos == 1
                vSeqGaps_5 = str2num(char(vSeqGaps(:,5)));
                vSeqGaps_7 = str2num(char(vSeqGaps(:,7)));
                vHypoToPlace = zeros(size(vSeqGaps_5,1),1);
                vBreakAfter = zeros(size(vSeqGaps_5,1),1);
                if exist('vHypoAmount','var')
                    vNi = size(vSeqGaps,1);
                    for vi = vNi:-1:1
                        if vVerbose >= 2
                            fprintf('Estimate number of gaps: %i of %i\n',
                                    vNi - vi + 1, vNi);
                        end
                        vBreakAfter(vi,1) = vInterChrom(vSeqGaps_5(vi,1),4);

                        if vSeqGaps_7(vi,1) >= vGeneLengthHypoMin
                            vHypoToPlace(vi,1) = vHypoAmount;
                        end
                    end
                else
                    error(
                        'vHypoAmount not defined, needs to be defined for InsertHypos = 1.');
                end
            else
                error('Not other option (>1) programmed.');
            end

            %Sort the Breaks:
            [~, vBreakAfter_sort_ID] = sort(vBreakAfter);

            %Remove several sequencing gaps within one intergenic region.
            vBreakafterkill = false(size(vSeqGaps,1),1);
            vOld_Breakafter_sort_ID = vBreakAfter_sort_ID(1,1);
            for vi = 2:size(vSeqGaps,1)
                if vBreakAfter(vBreakAfter_sort_ID(vi,1),1) == vBreakAfter(vOld_Breakafter_sort_ID,1)
                    vBreakafterkill(vBreakAfter_sort_ID(vi,1),1) = 1;
                    %It is possible that a sequencing gap does not meet the
                    %minimal Hypo length, but the next one does. We populate 2
                    %genes in total once one of them meets the criterion.
                    if vHypoToPlace(vOld_Breakafter_sort_ID,1) < vHypoToPlace(vBreakAfter_sort_ID(vi,1),1)
                        vHypoToPlace(vOld_Breakafter_sort_ID,1) = vHypoToPlace(vBreakAfter_sort_ID(vi,1),1);
                    end
                else
                    vOld_Breakafter_sort_ID = vBreakAfter_sort_ID(vi,1);
                end
            end
            vSeqGaps = vSeqGaps(~vBreakafterkill,:);
    %         vSeqGaps_5 = vSeqGaps_5(~vBreakafterkill,:);
    %         vSeqGaps_7 = vSeqGaps_7(~vBreakafterkill,:);
            vBreakAfter = vBreakAfter(~vBreakafterkill,:);
            vHypoToPlace = vHypoToPlace(~vBreakafterkill,:);

            [~, vBreakAfter_sort_ID] = sort(vBreakAfter);

            %Increase the size of the gene locations.
            for vk = 1:size(vGeneLocation_Att,2)
                vGeneLocation3.(char(vGeneLocation_Att(1,vk))) = cell(size(vGeneLocation2.(char(vGeneLocation_Att(1,vk))),1) + sum(vHypoToPlace), ...
                    size(vGeneLocation2.(char(vGeneLocation_Att(1,vk))),2));
            end

            %Increase the size of all variables that are connected to the gene
            %locations.
            vBiomG_G_3 = logical(spalloc(
                size(vBiomG_G,1) + sum(vHypoToPlace),size(vBiomG_G,2),size(vBiomG_G,2)*2));

            %Populate the new gene location variable with the hypothetical genes.
            vEnd_GeneLocation2 = size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1);
            vEnd_GeneLocation3 = size(vGeneLocation3.(char(vGeneLocation_Att(1,1))),1);
            for vi = size(vSeqGaps,1):-1:1
                if vVerbose >= 2
                    fprintf('Populate hypothetical genes: %i of %i\n',
                            vNi - vi + 1, vNi);
                end
                if vHypoToPlace(vBreakAfter_sort_ID(vi,1),1) ~= 0
                    vStart_GeneLocation2 = vBreakAfter(
                        vBreakAfter_sort_ID(vi,1),1) + 1;
                    vStart_GeneLocation3 = vEnd_GeneLocation3 - \
                        (vEnd_GeneLocation2 - vStart_GeneLocation2);
                    for vl = 1:vHypoToPlace(vBreakAfter_sort_ID(vi,1),1)
                        vPlaceLastGeneEnd = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vStart_GeneLocation2-1,1)));
                        vPlaceBetweenGenes = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(vStart_GeneLocation2,1))) - vPlaceLastGeneEnd - 1;
                        vPlaceToInsert = vPlaceLastGeneEnd + \
                            floor(vPlaceBetweenGenes/(vHypoToPlace(vBreakAfter_sort_ID(vi,1),1)+1)*vl);
                        for vk = 1:size(vGeneLocation_Att,2)
                            vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3:vEnd_GeneLocation3,:) = ...
                            vGeneLocation2.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation2:vEnd_GeneLocation2,:);
                            if vk == 1
                                vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3-vl,1) = cellstr(['Hypothetical_Gene_' num2str(sum(vHypoToPlace(1:vBreakAfter_sort_ID(vi,1)-1,1))+vl)]);
                            elseif vk == 2
                                vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3-vl,1) = cellstr(num2str(vPlaceToInsert));
                            elseif vk == 3
                                vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3-vl,1) = cellstr(num2str(vPlaceToInsert));
                            elseif vk == 4
                                vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3-vl,1) = vGeneLocation2.(char(vGeneLocation_Att(1,vk)))(vBreakAfter(vBreakAfter_sort_ID(vi,1),1),1);
                            elseif vk == 5
                                vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3-vl,1) = vGeneLocation2.(char(vGeneLocation_Att(1,vk)))(vBreakAfter(vBreakAfter_sort_ID(vi,1),1),1);
                            elseif vk == 6
                                vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(vStart_GeneLocation3-vl,1) = vGeneLocation2.(char(vGeneLocation_Att(1,vk)))(vBreakAfter(vBreakAfter_sort_ID(vi,1),1),1);
                            end
                        end
                    end
                    if vVerbose >= 2
                        fprintf('Inserting Hypothetical genes: Copy old %i to %i to new %i to %i\n',
                                vStart_GeneLocation2, vEnd_GeneLocation2, vStart_GeneLocation3, vEnd_GeneLocation3);
                    end
                    vBiomG_G_3(vStart_GeneLocation3:vEnd_GeneLocation3,:) = vBiomG_G(vStart_GeneLocation2:vEnd_GeneLocation2,:);
                    vIDs_to_change = vInterChrom(:,5)>vBreakAfter(vBreakAfter_sort_ID(vi,1),1);
                    vInterChrom(vIDs_to_change,5) = vInterChrom(vIDs_to_change,5) + vHypoToPlace(vBreakAfter_sort_ID(vi,1),1);
                    vIDs_to_change = vInterChrom(:,4)>vBreakAfter(vBreakAfter_sort_ID(vi,1),1);
                    vInterChrom(vIDs_to_change,4) = vInterChrom(vIDs_to_change,4) + vHypoToPlace(vBreakAfter_sort_ID(vi,1),1);
                    vEnd_GeneLocation2 = vStart_GeneLocation2-1;
                    vEnd_GeneLocation3 = vStart_GeneLocation3- \
                        (1+vHypoToPlace(vBreakAfter_sort_ID(vi,1),1));
                end
            end
            vStart_GeneLocation3 = 1;
            vStart_GeneLocation2 = 1;
            vBiomG_G_3(vStart_GeneLocation3:vEnd_GeneLocation3,:) = vBiomG_G(vStart_GeneLocation2:vEnd_GeneLocation2,:);
            for vk = 1:size(vGeneLocation_Att,2)
                vGeneLocation3.(char(vGeneLocation_Att(1,vk)))(1:vEnd_GeneLocation3,:) = ...
                vGeneLocation2.(char(vGeneLocation_Att(1,vk)))(1:vEnd_GeneLocation2,:);
            end
            vGeneLocation2 = vGeneLocation3;
            vBiomG_G = vBiomG_G_3;
        end
    end"""
    return hard_breaks  # FIXME: return value


def introduce_hard_breaks(
    vGeneLocation: pl.DataFrame,
    vUniqueChrom: List[str],
    vInterChrom: pl.DataFrame,
    vSeqGaps: pl.DataFrame,
    vMaxSeqGapSize: int,
    vMaxInterGeneDist: int,
    vPercentileForMaxInterGeneDist: float,
    vMaxInterGeneDistByMedian: int,
    vBiomG_G: Dict[Any, Any],
    vHardBreakGaps: List[int] = [],
) -> HardBreaks:
    """
    ---------------------------------------------------------------------------

      Function that prevents formation of clusters Sequencing gaps that
      exceed vMaxSeqGapSize or that are defined not to be crossed by
      vHardBreakGaps.

    ---------------------------------------------------------------------------

    Args:
        vGeneLocation: DataFrame about the gene location information,
                        generated in PlantClusterFinder function.
        vUniqueChrom: A list with unique entries of all chromosome and scaffold
                      names.
        vInterChrom: A DataFrame defining chromosome, start bp, end bp, gene
                     before and gene after a intergenic region.
        vSeqGaps: A string array (cell) DataFrame containing the
                  - chromosome/scaffold,
                  - the identification if it is coding (A) or masked (N) DNA
                  - start bp and
                  - length of a masked region.
                  In addition there is the information stored in with
                  intergenic region the gap is encoded and its start and end
                  bp, corrected for maximal size (= intergenic region).
        vMaxSeqGapSize: Maximal masked nucleotide region (in bp) that a cluster
                        is allowed to bridge (-1 or None is inactive).
        vMaxInterGeneDist: Maximal intergenic distance that still can be
                           gene cluster.
        vPercentileForMaxInterGeneDist: Prevents clusters to be spanned over
                                        the top n% of intergenic distances.
        vMaxInterGeneDistByMedian: Maximal intergenic distance that still can
                                   be spanned by a gene cluster: this number
                                   times the median of the gene sizes.
        vBiomG_G: Mapping between the Gene location file and the gene-IDs in the
                  conversion file.
        vHardBreakGaps: A list of numbers. The masked nucleotides sequences of
                        these specific lengths are preventing cluster formation
                        (Default is an empty list []).

    Returns:
        HardBreaks object containing the gene location information, the
            inter-chromosomal matrix, and the SeqGaps populated.
    """

    _ = """
        % MATLAB code
        vCurrArgin = 11;
        if nargin > vCurrArgin
            vi = 1;
            while vi <= nargin - vCurrArgin

                vi = vi +1;
            end
        end

        vSeqSize = str2num(char(vSeqGaps(:,end)));
        vSeqHard = false(size(vSeqSize,1),1);
        for vi = 1:size(vHardBreakGaps,2)
            vSeqHard = vSeqHard | vSeqSize == vHardBreakGaps(1,vi);
        end
        if vMaxSeqGapSize ~= -1
            vSeqHard = vSeqHard | vSeqSize >= vMaxSeqGapSize;
        end
        vSeqGaps_Hard = vSeqGaps(vSeqHard,:);
        vSeqGaps = vSeqGaps(~vSeqHard,:);

        %Check maximal intergenic distance via percent
        if vPercentileForMaxInterGeneDist ~= -1
            %Change From Median to Percentile
            %vMultiplierInterGeneDist = median(vInterChrom(:,3)-vInterChrom(:,2)+1);
            %vMaxInterGeneDist_temp = vMultiplierInterGeneDist*vPercentileForMaxInterGeneDist;
            %Use 95% Percentile Cutoff Instead
            vMaxInterGeneDist_temp = prctile(vInterChrom(:,3)-vInterChrom(:,2)+1, vPercentileForMaxInterGeneDist);
            fprintf('Percentile used for intergenic distance: %ibp\nClusters broken when spanning distance larger: %ibp.\n',
                    vPercentileForMaxInterGeneDist, vMaxInterGeneDist_temp);
            if vMaxInterGeneDist == -1
                vMaxInterGeneDist = vMaxInterGeneDist_temp;
            else
                vMaxInterGeneDist = min(
                    [vMaxInterGeneDist,vMaxInterGeneDist_temp]);
            end
        end

        %Check maximal intergenic distance via median gene size
        if vMaxInterGeneDistByMedian ~= -1
            vMaxInterGeneDist_temp = prctile(str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(:,1))) - ...
                                             str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(:,1))), 50);
            if vMaxInterGeneDist == -1
                vMaxInterGeneDist = vMaxInterGeneDist_temp*vMaxInterGeneDistByMedian;
            else
                vMaxInterGeneDist = min(
                    [vMaxInterGeneDist,vMaxInterGeneDist_temp*vMaxInterGeneDistByMedian]);
            end
        end

        %Check intergenic distances
        if vMaxInterGeneDist ~= -1
            vG_start = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(:,1)));
            vG_end = str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(:,1)));
            vAct_Crom = '';
            vNG = size(vGeneLocation2.(char(vGeneLocation_Att(1,4))),1);
            for vi = vNG:-1:2
                if strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi-1,1),vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi,1))
                    if vG_end(vi,1) - vG_start(vi-1,1) >= vMaxInterGeneDist
                        if ~strcmp(vAct_Crom,vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi,1))
                            vChrom_counter = 1;
                            vAct_Crom = vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi,1);
                        else
                            vChrom_counter = vChrom_counter + 1;
                        end
                        vIDs = strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(:,1),vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi,1));
                        vIDs(1:vi-1,1) = 0;
                        vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vIDs,1) = cellstr([char(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vi,1)) '_' num2str(vChrom_counter)]);
                    end
                end
            end
        end

        %Check physical hard breaks
        vGeneLocation2.([char(vGeneLocation_Att(1,4)) '_old']) = vGeneLocation2.(char(vGeneLocation_Att(1,4)));
        vGeneLocation_Att = [vGeneLocation_Att, [char(vGeneLocation_Att(1,4)) '_old']];
        for vi = size(vUniqueChrom,1):-1:1
            vTempHardGaps = vSeqGaps_Hard(strcmp(vSeqGaps_Hard(:,1),vUniqueChrom(vi,1)),:);
            if size(vTempHardGaps,1) > 0
                vTempChrom = cell(size(vTempHardGaps,1)+1,1);
                for vj = 1:size(vTempChrom,1)
                    vTempChrom(vj,1) = cellstr([char(vUniqueChrom(vi,1)) '_' num2str(vj)]);
                end
                for vj = size(vTempHardGaps,1):-1:1
                    vBreakAfter = vInterChrom(
                        str2num(char(vTempHardGaps(vj,5))),4);
                    vTemp_IDs = strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(:,1),vUniqueChrom(vi,1));
                    vTemp_IDs(1:vBreakAfter) = 0;
                    vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vTemp_IDs,1) = vTempChrom(vj + 1,1);
                end
                vTemp_IDs = strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(:,1),vUniqueChrom(vi,1));
                vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vTemp_IDs,1) = vTempChrom(1,1);
            end
        end
    end"""
    return HardBreaks(
        GeneLocation=vGeneLocation, InterChrom=vInterChrom, SeqGaps=vSeqGaps
    )  # FIXME: return value


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

      Functions that uses the output of f_get_Sequencing_Gaps and reformulates
      it such that it is clear in which intergenic region the sequening gap
      lies.

    ---------------------------------------------------------------------------

    Args:
        vFile: Fasta file of the masked genome nucleotide sequences.
        vInterGeneLocation: A DataFrame defining chromosome, start bp, end bp, gene
                     before and gene after a intergenic region.
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
                  - chromosome/scaffold
                  - the identification if it is coding (A) or masked (N) DNA
                  - start bp
                  - length of a masked region
                  In addition there is the information stored in with
                  intergenic region the gap is encoded and its start and end
                  bp, corrected for maximal size (= intergenic region).
    """

    _ = """
        % MATLAB code
        vCurrNargin = 7;
        vGeneLengthHypoMin = -1;
        vSeqGapSizesChromBreak = [];
        vSet = 0;
        if nargin > vCurrNargin
            vi = 1;
            while vi <= nargin - vCurrNargin
                if ischar(varargin{1,vi})
                    if strcmp(varargin{1,vi},'SeqGapSizesChromBreak')
                        vSeqGapSizesChromBreak = varargin{1,vi+1};
                        vSet = vSet + 1;
                        vi = vi + 1;
                    elseif strcmp(varargin{1,vi},'GeneLengthHypoMin')
                        vGeneLengthHypoMin = varargin{1,vi+1};
                        vSet = vSet + 1;
                        vi = vi + 1;
                    end
                end
                vi = vi + 1;
            end
        end
        if vGeneLengthHypoMin ~= -1 && vSet ~= 2
            error('f_annotate_Sequencing_Gaps needs both SeqGapSizesChromBreak and GeneLengthHypoMin to be set or to be unset.');
        end
    """
    if not (vGeneLengthHypoMin == -1 or vSeqGapSizesChromBreak == []):
        log.error(
            "annotate_sequencing_gaps needs both SeqGapSizesChromBreak and GeneLengthHypoMin to be set or to be unset."
        )
        raise Exception("Incompatable parameters in annotate_sequencing_gaps function")
    _ = """
        vMaxInit = 7000000;
        if ispc()
            if ~exist([vFile '_GAPOutput'], 'file') || vOverwriteSeqGapInfo
                f_get_Sequencing_Gaps(
                    vFile, vPara, vVerbose, vOverwriteSeqGapInfo, vUnmaskedDNA)
            end
        else
            [~,vOutput] = system(['ls' ' ' '''' vFile '_GAPOutput''']);
            if strncmp(vOutput,'ls: cannot access', 17) || vOverwriteSeqGapInfo
                f_get_Sequencing_Gaps(
                    vFile, vPara, vVerbose, vOverwriteSeqGapInfo, vUnmaskedDNA)
            end
        end

        vFIO = fopen([vFile '_GAPOutput']);
        vSeqGaps = cell(vMaxInit,4);
        vLine = fgetl(vFIO);
        vL = 1;
        while ischar(vLine)
            %Split the line into four columns using the space-character as
            %delimiter.
            vSplit = regexp(vLine,' ','split');
            if size(vSplit,2) ~= 4
                error('%s is corrupt at line %i. It should contain 4 columns split by spaces, it contains %i.\n\n', [vFile '_GAPOutput'],vL, size(vSplit,2));
            end

            if vGeneLengthHypoMin == -1 || sum(vSeqGapSizesChromBreak == str2num(char(vSplit(1,4)))) > 0 || (vGeneLengthHypoMin ~= -1 && vGeneLengthHypoMin <= str2num(char(vSplit(1,4))))
                if vL > vMaxInit
                    vSeqGaps = [vSeqGaps; cell(1,4)];
                end
                vSeqGaps(vL,:) = vSplit;
                vL = vL + 1;
            end
            vLine = fgetl(vFIO);
        end
        fclose(vFIO);
        vSeqGaps = vSeqGaps(1:(vL-1),:);
    """
    gap_output_file = vFile.with_suffix(vFile.suffix + "_GAPOutput")
    if vOverwriteSeqGapInfo or not gap_output_file.exists():
        get_sequencing_gaps(vFile, 1, vOverwriteSeqGapInfo, vUnmaskedDNA)
    sequencing_gaps = pl.read_csv(
        gap_output_file,
        separator=" ",
        has_header=False,
        new_columns=["scaffold", "amino acid", "start index", "length"],
    )
    # find the intergenic region of the gap, if the gap is present within a
    # gene model, get rid of the gap.
    _ = """
        vSeqGaps = [vSeqGaps, cell(size(vSeqGaps,1),3)];
        vSeqGaps_3 = str2num(char(vSeqGaps(:,3)));
        vSeqGaps_4 = str2num(char(vSeqGaps(:,4)));
        vSeqGaps_4 = vSeqGaps_3 + vSeqGaps_4 - 1;

        vnvi = size(vSeqGaps,1);
        vKillSeqGaps = false(vnvi,1); % Check which sequences should not be processed
        vKillSeqGaps = vKillSeqGaps | ~strcmp(vSeqGaps(:,2),'N'); % Do not process all non N sequences.
        % Do not process sequences that do not lie on a chromosome we want.
        vHasGeneLocChrom = false(vnvi,1);
        for vi = 1:size(vUniqueChrom,1)
            vHasGeneLocChrom = vHasGeneLocChrom | strcmp(vSeqGaps(:,1),vUniqueChrom(vi,1));
        end
        vBadUniqueChrom = unique(vSeqGaps(~vHasGeneLocChrom,1));
        for vi = 1:size(vBadUniqueChrom,1)
            warning('DNA Fasta file ''%s'' contains chromosomes or scaffolds not covered with the gene position file: %s\n',
                    vFile, char(vBadUniqueChrom(vi,1)));
        end
        vKillSeqGaps = vKillSeqGaps | ~vHasGeneLocChrom;
        vSeqGaps = vSeqGaps(~vKillSeqGaps,:);
        vSeqGaps_3 = vSeqGaps_3(~vKillSeqGaps,:);
        vSeqGaps_4 = vSeqGaps_4(~vKillSeqGaps,:);
        vnvi = size(vSeqGaps,1);
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
        vSeqGaps = vSeqGaps(~vKillSeqGaps,:);
    %     vSeqGaps_3 = vSeqGaps_3(~vKillSeqGaps,:);
    %     vSeqGaps_4 = vSeqGaps_4(~vKillSeqGaps,:);
    end"""
    # FIXME: only gaps which are in intergenic regions
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

      Creates a _GAPOutput file from the FASTA file, then creates a
      _GAPOutput_count.txt file by running analyze_plant_cluster_gap_file()
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
            [f"Class{n+1}" for n in range(6)],
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
    vParallelize: int = 1,
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
        "'6 qseqid sseqid evalue bitscore length pident qstart qend sstart send'",
    ]
    if vParallelize > 1:
        command_for["blastp"].extend(["-num_threads", str(vParallelize)])
    command_for["cut"] = [
        "-f",
        "1,2,3",
        f"{vProtein_sequence_FastaFile}.blastOUT",
        ">",
        f"{vProtein_sequence_FastaFile}.MCLIN",
    ]
    command_for["mcxload"] = [
        "-abc",
        f"{vProtein_sequence_FastaFile}.MCLIN",
        "--stream-mirror",
        "--stream-neg-log10",
        "-stream-tf",
        "'ceil(200)'",
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
        log.info(f"{app_path[app]} {' '.join(command_for[app])}")
        subprocess.run(
            str(app_path[app]) + " " + " ".join(command_for[app]),
            shell=True,
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
        log.error(f"Primary table join ID {[join_on[0]]} is null dtype, cannot join")
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
    vGeneLocation: pl.DataFrame,
    vBiomG_E: Dict[Any, Any],
    vEnzMinForClusters: int,
    vKbRangeLD: int,
    vGeneRangeLD: int,
    mcl_cluster_data: MclClusterData,
    vBiomG_G: Dict[Any, Any],
    genes_to_reactions: GenesToReactions,
    vMapG_Gpgdb: Dict[Any, Any],
    vCriterion: int,
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
      aswell. Clusters need at least vEnzMinForClusters enzymes.

    ---------------------------------------------------------------------------

    Args:
        vMinStepSize: Minimal intervening gene size that should be computed
                      (default is 5).
        vMaxStepSize: Maximal intervening gene size that should be computed
                      (default is 20).
        vGeneLocation: DataFrame about the gene location information,
                        generated in PlantClusterFinder function.
        vBiomG_E: A dictionary describing which gene has at least one reaction
                  annotated in the pgdb (1, if no reaction 0)
        vEnzMinForClusters: Defines the minimum number of enyzmes that need to
                            be present in a cluster.
        vKbRangeLD: Genes that are encoded more than this bp apart, are not
                    called local duplicates.
        vGeneRangeLD: Number of genes that space two genes from each other.
                      Genes that have more than this many intervening genes
                      between themare not called local duplicates.
        vG_MCL: Mapping between the genes in the conversion-ID file and the
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
                    2: stepsize 3 (because most organisms take this one
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

    _ = """
        % MATLAB code
        vContinueSearch = 1;
        vClusters = zeros(vMinStepSize*size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1),5); % Stepsize, start, end, number of enzymes, number of hypos;
        vClusterBackground = zeros(vMinStepSize*size(vGeneLocation2.(char(vGeneLocation_Att(1,1))),1),5); % Size, number of enzymes, count, number of total;
        vClusterWrite = 1;
        vClusterBackgroundWrite = 1;
        vStepSize = 0;
        vBiom_Hypo = strncmp(vGeneLocation2.(char(vGeneLocation_Att(1,1))),'Hypo',4);
        while vStepSize <= vMinStepSize || vContinueSearch == 1
            % Get all clusters
            vi = 1;
            while vi <= size(vBiomG_E,1)
                if vVerbose >= 2
                    fprintf('Potential Gene start %i of %i\n',
                            vi, size(vBiomG_E,1));
                end
                vStart = vi - 1 + find(vBiomG_E(vi:end,1),1,'first');
                vEnd = vStart;
                vAbolishCluster = 0;
                if size(vStart,1) == 0
                    vAbolishCluster = 1;
                    vStart = size(vBiomG_E,1) + 1;
                    vEnd = vStart;
                else
                    vChrom_start = vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vStart,1);
                    vBiomG_Chrom = strcmp(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(:,1),vChrom_start);

                    vStopSearchEnd = 0;
                    while vStopSearchEnd == 0
                        vEnds = vEnd+1;
                        vEnde = vEnd+1+vStepSize;
                        if vEnde > size(vBiomG_E,1)
                            vEnde = size(vBiomG_E,1);
                        end
                        if size(find(vBiomG_E(vEnds:vEnde,1) & vBiomG_Chrom(vEnds:vEnde,1),1,'last'),1)>0
                            vEnd = vEnds - 1 + find(vBiomG_E(vEnds:vEnde,1) & vBiomG_Chrom(vEnds:vEnde,1),1,'last');
                        else
                            vStopSearchEnd = 1;
                        end
                    end
                end

                %Cluster needs at least 3 enzymes
                if vAbolishCluster == 0
                    if sum(vBiomG_E(vStart:vEnd)) < vEnzMinForClusters
                        vAbolishCluster = 1;
                    end
                end

                if vAbolishCluster == 0
                    vLDisNOTaProblem = 0;
                    % Get rid of clusters that only contain LD
                    vEIDs = vStart+find(vBiomG_E(vStart:vEnd))-1;
    %                 vEIDs_find = find(vEIDs);
                    % Test if they are encoded within kb range
                    vStart_kb = min([str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(vEIDs,1))); ...
                        str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vEIDs,1)))]);
                    vEnd_kb = max([str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,2)))(vEIDs,1))); ...
                        str2num(char(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vEIDs,1)))]);
                    if vEnd_kb-vStart_kb+1 > vKbRangeLD
                        vLDisNOTaProblem = 1;
                    end
                    % Test if they are encoded witin gene range
                    if vLDisNOTaProblem == 0
                        if sum(~strncmp(vGeneLocation2.(char(vGeneLocation_Att(1,3)))(vStart:vEnd,1),'Hypothetical_Gene_',18)) >= vGeneRangeLD+2
                            vLDisNOTaProblem = 1;
                        end
                    end
                    %Check if there are at least two local duplication groups
                    if vLDisNOTaProblem == 0
                        vMCL = sum(vG_MCL(sum(vBiomG_G(vEIDs(1,1),:),1)>0,:),1)>0;
                        vj = 1;
                        vTwoGroups = 0;
                        while vj <= size(vEIDs,1) && vTwoGroups == 0
                            if sum(vMCL & sum(vG_MCL(sum(vBiomG_G(vEIDs(vj,1),:),1)>0,:),1)==0) > 0
                                vTwoGroups = 1;
                            end
                            vj = vj + 1;
                        end
                        if vTwoGroups == 1
                            vLDisNOTaProblem = 1;
                        end
                    end
                    if vLDisNOTaProblem == 0
                        vAbolishCluster = 1;
                    end
                end
                % Get rid of clusters that only have one reaction set
                if vAbolishCluster == 0
                    vRxns = sum(Gene2Rxn(sum(vMapG_Gpgdb(sum(vBiomG_G(vEIDs(1,1),:),1)>0,:),1)>0,:),1)>0;
                    vj = 1;
                    vTwoGroups = 0;
                    while vj <= size(vEIDs,1) && vTwoGroups == 0
                        if sum(vRxns & sum(Gene2Rxn(sum(vMapG_Gpgdb(sum(vBiomG_G(vEIDs(vj,1),:),1)>0,:),1)>0,:),1)==0) > 0
                            vTwoGroups = 1;
                        end
                        vj = vj + 1;
                    end
                    if vTwoGroups == 0
                        vAbolishCluster = 1;
                    end
                end
                if vAbolishCluster == 0
                    vClusters(vClusterWrite,:) = [vStepSize,vEIDs(1,1),vEIDs(end,1),size(vEIDs,1),sum(vBiom_Hypo(vEIDs(1,1):vEIDs(end,1),1))];
                    vClusterWrite = vClusterWrite + 1;
                end
                vi = vEnd + 1;
            end

            vC_IDs = vClusters(:,1) == vStepSize;
            vC_IDs(vClusterWrite:end,1) = 0;

            % Compute background of clusters
            % Get all sizes of clusters of this round
            [vClusterBackground, vClusterBackgroundWrite] = f_compute_background(
                vClusterBackground, vClusters, vClusterBackgroundWrite, vGeneLocation2, vGeneLocation_Att, vBiomG_E, vC_IDs);

            % Evaluate clusters for criterion (stepsize found?)
            if vVerbose >= 2
                if vCriterion == 1
                    fprintf('Number of enzymes - number of non-enzymes in clustered genes: %i\n', sum(vClusters(
                        vC_IDs,4))-(sum(vClusters(vC_IDs,3)-vClusters(vC_IDs,2)+1)-sum(vClusters(vC_IDs,4))));
                elseif vCriterion == 3
                    fprintf('Number of enzymes - number of non-enzymes (no hypos) - 2* number of hypos in clustered genes: %i\n', sum(vClusters(vC_IDs,4))-(
                        sum(vClusters(vC_IDs,3)-vClusters(vC_IDs,2)+1)-sum(vClusters(vC_IDs,4))-sum(vClusters(vC_IDs,5)))-2*sum(vClusters(vC_IDs,5)))
                elseif vCriterion == 4
                    if vStepSize ~= 0
                        fprintf('Number of clustes more than last stepsize: %i\n', size(vClusters(vC_IDs,1)) - size(vClusters(:,1)==vStepSize-1));
                    else
                        fprintf('Number of clustes Stepsize 0: %i\n',
                                size(vClusters(vC_IDs,1)));
                    end
                end
            end
            if vCriterion == 1
                if vStepSize == vMaxStepSize
                    vContinueSearch = 0;
                    vFinalStepsize = vStepSize;
                    fprintf(
                        'Warning: maximal step size reached: %i. Consider increasing this variable.\n', vMaxStepSize)
                elseif sum(vC_IDs)==0
                    fprintf(
                        'Warning: no clusters found for stepsize %i.\n', vStepSize);
                elseif sum(vClusters(vC_IDs,4))/sum(vClusters(vC_IDs,3)-vClusters(vC_IDs,2)+1)<0.5
                    if vContinueSearch == 1
                        vFinalStepsize = vStepSize-1;
                    end
                    vContinueSearch = 0;
                end
            elseif vCriterion == 2
                if vStepSize == vMinStepSize
                    vFinalStepsize = 3;
                    vContinueSearch = 0;
                end
            elseif vCriterion == 3
                if vStepSize == vMaxStepSize
                    vContinueSearch = 0;
                    vFinalStepsize = vStepSize;
                    fprintf(
                        'Warning: maximal step size reached: %i. Consider increasing this variable.\n', vMaxStepSize)
                elseif sum(vC_IDs)==0
                    fprintf(
                        'Warning: no clusters found for stepsize %i.\n', vStepSize);
                elseif sum(vClusters(vC_IDs,4))-(sum(vClusters(vC_IDs,3)-vClusters(vC_IDs,2)+1)-sum(vClusters(vC_IDs,4))-sum(vClusters(vC_IDs,5)))-2*sum(vClusters(vC_IDs,5)) < 0
                    if vContinueSearch == 1
                        vFinalStepsize = vStepSize-1;
                    end
                    vContinueSearch = 0;
                end
            elseif vCriterion == 4
                if vStepSize == vMaxStepSize
                    vContinueSearch = 0;
                    vFinalStepsize = vStepSize;
                    fprintf(
                        'Warning: maximal step size reached: %i. Consider increasing this variable.\n', vMaxStepSize)
                elseif sum(vC_IDs)==0
                    fprintf(
                        'Warning: no clusters found for stepsize %i.\n', vStepSize);
                elseif size(vClusters(vC_IDs,1)) - size(vClusters(:,1)==vStepSize-1) < 0
                    if vContinueSearch == 1
                        vFinalStepsize = vStepSize-1;
                    end
                    vContinueSearch = 0;
                end
            end
            vStepSize = vStepSize + 1;
        end
        vClusters = vClusters(1:vClusterWrite-1,:);
        vClusterBackground = vClusterBackground(vClusterBackground(:,4) ~= 0,:);
        vStepSize_Max = vStepSize-1;
        fprintf('Stepsize chosen: %i\n', vFinalStepsize);
    end"""
    return Clusters(
        ClusterMeta=pl.DataFrame(),
        ClusterBackground=pl.DataFrame(),
        FinalStepsize=0,
        MaxStepSize=0,
    )  # FIXME: return value


def compute_background(
    vClusters: Clusters,
    vClusterBackgroundWrite,
    vGeneLocation: pl.DataFrame,
    vBiomG_E: List[int],
    vC_IDs: List[str],
) -> tuple[pl.DataFrame, int]:
    """
    ---------------------------------------------------------------------------

      Function that lookes at all clusters that were found by the function
      f_perform_clustering(), and computes the background information for the
      windowsizes of these clusters (e.g. a cluster has 4 enzymes and two non
      enzymes, then its size is 6. The function goes with a window size 6
      through the genome and checks how many windows exist with 0, 1, 2, ...,
      5, 6 enzymes. It then calculates the probability to see a cluster with
      more or equal than 0 enzymes, more or equal than 1, ... more or equal
      than 6 enzymes. This is then used to define top n% clusters, e.g.
      clusters that are at least as dense as the top 5% of windowsizes are
      marked to be of higher interest.

    ---------------------------------------------------------------------------

    Args:
        vClusterBackground: Density information about the genome (window size,
                            number of enzymes, Occurrence of windowsize with
                            this many enzymes in the genome, occurence of
                            this windowsize in the genome, Frequency of seeing
                            a windowsize of this size with at least as many
                            enzymes in the genome)
        vClusters: Cluster Information (stepsize, chromosome id, gene start,
                   gene end, number of enzymes)
        vClusterBackgroundWrite: Defines where the next entry should be written
                                 in vClusterBackground.
        vGeneLocation2: DataFrame about the gene location information,
                        generated in PlantClusterFinder function.
        vBiomG_E: A vector describing which gene has at least one reaction
                  annotated in the pgdb (1, if no reaction 0)
        vC_IDs: Defines which Entries in vClusters are Clusters of the current
                Stepsize that was calculated and should be checked if the
                background of that size of cluster (number of genes) was
                already calculated.

    Returns:
        Tuple with the following information:
        - **vClusterBackground**: Density information about the genome (window size,
                            number of enzymes, Occurrence of windowsize with
                            this many enzymes in the genome, occurence of
                            this windowsize in the genome, Frequency of seeing
                            a windowsize of this size with at least as many
                            enzymes in the genome)
        - **vClusterBackgroundWrite**: Defines where the next entry should be written
                                 in vClusterBackground.
    """

    _ = """
        % MATLAB code
         vBackgroundToCheck = unique(
             vClusters(vC_IDs,3)-vClusters(vC_IDs,2)+1);
         % Check which backgrounds have already been produced and compute new backgrounds
         vi = 1;
         while vi <= size(vBackgroundToCheck,1)
             if sum(vClusterBackground(1:vClusterBackgroundWrite-1,1) == vBackgroundToCheck(vi,1)) == 0
                 vClusterBackground(vClusterBackgroundWrite:vClusterBackgroundWrite+vBackgroundToCheck(vi,1),1) = vBackgroundToCheck(vi,1);
                 vClusterBackground(vClusterBackgroundWrite:vClusterBackgroundWrite+vBackgroundToCheck(vi,1),2) = 0:1:vBackgroundToCheck(vi,1);
                 for vk = 1:size(vBiomG_E,1)-(vBackgroundToCheck(vi,1)-1)
                     if size(unique(vGeneLocation2.(char(vGeneLocation_Att(1,4)))(vk:(vk+vBackgroundToCheck(vi,1)-1),1)),1) == 1
                         vClusterBackground(vClusterBackgroundWrite:vClusterBackgroundWrite+vBackgroundToCheck(vi,1),4) = vClusterBackground(vClusterBackgroundWrite:vClusterBackgroundWrite+vBackgroundToCheck(vi,1),4) + 1;
                         vClusterBackground(vClusterBackgroundWrite+sum(vBiomG_E(vk:vk+vBackgroundToCheck(vi,1)-1,1)),3) = vClusterBackground(vClusterBackgroundWrite+sum(vBiomG_E(vk:vk+vBackgroundToCheck(vi,1)-1,1)),3) + 1;
                     end
                 end
                 for vk = vClusterBackgroundWrite+vBackgroundToCheck(vi,1):-1:vClusterBackgroundWrite
                     vClusterBackground(vk,5) = sum(vClusterBackground(vk:vClusterBackgroundWrite+vBackgroundToCheck(vi,1),3))/vClusterBackground(vk,4);
                 end
                 vClusterBackgroundWrite = vClusterBackgroundWrite + \
                     vBackgroundToCheck(vi,1) + 1;
             end
             vi = vi + 1;
         end
    end"""
    return (pl.DataFrame(), 0)  # FIXME: return values


def perform_cutoff(clusters: Clusters, vTopPercentClusters: int) -> List[float]:
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
    return []  # FIXME: return value


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
    vClusters_Top: List[float],
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
                    f'{"CHOSEN " if i==clusters.FinalStepsize else ""}Stepsize {i}',
                    f'{"CHOSEN " if i==clusters.FinalStepsize else ""}Stepsize {i}, Top {vTopPercentClusters}%',
                ]
        elif vOutputType == "simple":
            header += [f"Stepsize {clusters.FinalStepsize}, Top {vTopPercentClusters}%"]
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
    vClusters_Top: List[float],
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
    def _row_gen(fh, attrib, unique_keys):
        """
        Generate a row from the pgdb flat file
        """
        y = {att: None for att in attrib}
        for line in fh:
            if line.startswith("#"):
                continue
            elif line.startswith("//"):
                yield y
                y = {att: None for att in attrib}
                continue
            elif " - " in line:
                key, *value = line.split(" - ")
                value = (" - ".join(value)).strip()
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
