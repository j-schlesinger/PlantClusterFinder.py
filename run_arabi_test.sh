#!/usr/bin/env bash
export ARABIDIR=/debian/home/jschlesi/workspace/bioinformatics/Arabidopsis
mkdir -p Outputs
rm ./Outputs/Athaliana_Clust_v2_0.tsv 2>/dev/null
rm ./Outputs/Athaliana_Gene_v2_0.tsv 2>/dev/null
uv run python -m plantclusterfinder_py\
   '-r' 'Arabidopsis_TestRun'\
   '--GenesDat_GeneInfo' 'PRODUCT'\
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
   '-dnaf' "${ARABIDIR}/ref_genome/GCF_000001735.4_TAIR10.1_genomic.fna"\
   '-gff3' "${ARABIDIR}/ref_genome/sorted_and_tidy.gff3"\
   '-gff3_gene' "locus_tag"\
   '-gff3_prot' "orig_protein_id"\
   '-gff3_regex' "[^|]+$"\
   '-gtpf' "${ARABIDIR}/ref_genome/gtpf_GCF_000001735.4_TAIR10.1_genomic.tsv"\
   '-psf'  "${ARABIDIR}/ref_genome/GCF_000001735.4_TAIR10.1_genomic.protein.faa"\
   '-glof' "${ARABIDIR}/ref_genome/glof_sorted_and_tidy.tsv"\
   '-pgdb' "${ARABIDIR}/aracyc/18.0.0/data/"\
   '-sitf' './Inputs/scaffold-tailoring-reactions-05082016.tab'\
   '-rmdf' './Inputs/ReactionMetabolicDomainClassification.txt'\
   '-gout' './Outputs/Athaliana_Gene_v2_0.tsv'\
   '-cout' './Outputs/Athaliana_Clust_v2_0.tsv'\
   '--Criterion' '2'\
   '--SeqGapSizesChromBreak' 10000 '--PGDBIdsToMap' 'GTP'\
   '--OutputType' 'simple'\
   --Verbose 2\
   "--OutputGeneLocationParquet" "./Outputs/Athaliana_gene_loc.parquet"\
   --NoRerunBlastP\
   --NoRebuildBlastDB
