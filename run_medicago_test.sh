#!/usr/bin/env bash
export MEDICAGODIR=/debian/home/jschlesi/workspace/bioinformatics/medicago/
mkdir -p Outputs
rm ./Outputs/Medicago_Clust_v2.tsv 2>/dev/null
rm ./Outputs/Medicago_Gene_v2.tsv 2>/dev/null
uv run python -m plantclusterfinder_py\
   '--GenesDat_GeneInfo' 'ACCESSION-1' 'PRODUCT'\
   '-md'\
   "Amines and Polyamines Metabolism" "Amino Acids Metabolism"\
   "Carbohydrates Metabolism" "Cofactors Metabolism"\
   "Detoxification Metabolism" "Energy Metabolism"\
   "Fatty Acids and Lipids Metabolism" "Hormones Metabolism"\
   "Inorganic Nutrients Metabolism" "Nitrogen-Containing Compounds"\
   "Nucleotides Metabolism" "Phenylpropanoid Derivatives"\
   "Polyketides" "Primary-Specialized Interface Metabolism"\
   "Redox Metabolism" "Specialized Metabolism" "Sugar Derivatives"\
   "Terpenoids"\
   '-dnaf' "${MEDICAGODIR}/ref_genome/ncbi_dataset/data/GCA_003473485.2/GCA_003473485.2_MtrunA17r5.0-ANR_genomic.fna"\
   '-gff3' "${MEDICAGODIR}/ref_genome/ncbi_dataset/data/GCA_003473485.2/genomic.gff"\
   '-gff3_gene' "locus_tag"\
   '-gff3_prot' "ID"\
   '-gff3_regex' "([^.]+)$"\
   '-pgdb' "${MEDICAGODIR}/mtruncatulacyc/6.0.0.1/data/"\
   '-sitf' './Inputs/scaffold-tailoring-reactions-05082016.tab'\
   '-rmdf' './Inputs/ReactionMetabolicDomainClassification.txt'\
   '-gout' './Outputs/Medicago_Gene_v2.tsv'\
   '-cout' './Outputs/Medicago_Clust_v2.tsv'\
   '--SeqGapSizesChromBreak' 10000 '--PGDBIdsToMap' 'GTP'\
   '--OutputType' 'simple'\
   --Verbose 2\
   "--OutputGeneLocationParquet" "./Outputs/Medicago_gene_loc.parquet"\
   --NoRerunBlastP\
   --NoRebuildBlastDB
