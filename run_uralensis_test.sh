#!/usr/bin/env bash
export LICORICEDIR=/debian/home/jschlesi/workspace/bioinformatics/uralensis
mkdir -p Outputs
rm ./Outputs/Gu_kyotoline_Clust_v2.tsv 2>/dev/null
rm ./Outputs/Gu_kyotoline_Gene_v2.tsv 2>/dev/null
uv run python -m plantclusterfinder_py\
   '-r' 'Gu_kyotoline_Test_Run'\
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
   '-dnaf' "${LICORICEDIR}/ref_genome/Glur_r1.0_hardmask_8chr_wrapped.fna"\
   '-gff3' "${LICORICEDIR}/ref_genome/Glur_r1.01_HC.gff3"\
   '-gff3_gene' "ID"\
   '-gff3_prot' "ID"\
   '-gff3_regex' ".*"\
   '-gtpf' "${LICORICEDIR}/ref_genome/gtpf_Glur_r1.0_hardmask_8chr_wrapped.tsv"\
   '-glof' "${LICORICEDIR}/ref_genome/glof_Glur_r1.01_HC.tsv"\
   '-psf' "${LICORICEDIR}/ref_genome/Glur_r1.01_pep_HC_8chrs.faa"\
   '-pgdb' "${LICORICEDIR}/pgdb/gkyoto_cyccyc/1.0/data/"\
   '-sitf' './Inputs/scaffold-tailoring-reactions-05082016.tab'\
   '-rmdf' './Inputs/ReactionMetabolicDomainClassification.txt'\
   '-gout' './Outputs/Gu_kyotoline_Gene_v2.tsv'\
   '-cout' './Outputs/Gu_kyotoline_Clust_v2.tsv'\
   '--SeqGapSizesChromBreak' 10000 '--PGDBIdsToMap' 'GTP'\
   '--OutputType' 'simple'\
   --Verbose 2\
   "--OutputGeneLocationParquet" "./Outputs/Gu_kyotoline_gene_loc.parquet"\
   --NoRerunBlastP\
   --NoRebuildBlastDB
