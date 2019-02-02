#!/bin/bash
#set -o errexit

# assigned random user name 
input2=$1
# email address
input3=$2

if [ $# -lt 2 ]; then
  echo 1>&2 "$0: not enough arguments"
  exit 2
elif [ $# -gt 2 ]; then
  echo 1>&2 "$0: too many arguments"
  exit 2
fi

## Set path to phenotype and genotype dataset
pheno_path=/GWD/appbase/projects/RD-TSci-UKB/PheWAS/phenotype500kV3
genotype_path=/GWD/appbase/projects/RD-TSci-UKB/ukb20361/Nov2016/data
output_path=/GWD/appbase/projects/RD-TSci-UKB/PheWAS_output_storage
script_path=/GWD/appbase/projects/RD-TSci-UKB/PheWAS/script500k

HES_Self_cate_QT=Anno_HESelf_QT.csv
HES_Self_cate_BI=Anno_HESelf_BI.csv

#create a user folder for performance 
foldername=$input2
temp_path=$output_path/$foldername

# Set file names for phenotype(s), covariate(s) and the subset of genes
BI_pheno_name=phenotypes_BI.txt
QT_pheno_name=phenotypes_QT.txt
covar_name=covars.txt
input_name=scor_names.txt
set_genes=gene_names.txt
set_snps=snp_names.txt
ind_pheno_name=test.txt

## Set path to PLINK binary, may need >=1.9
#PLINK=/GWD/appbase/projects/GXapp/plink/plink-1.9/plink
Rr=/home/anaconda/anaconda2/bin/R

# Combine all results to make a summary table

$Rr --vanilla --slave --args $temp_path $temp_path/$input_name  $pheno_path < $script_path/combine_SCORE_results.R

# Create a directory for plots= 'plots' 
mkdir -p $temp_path/scor_assoc/plots

# Generate plots for all SNPs across phenotypes
$Rr --vanilla --slave --args $temp_path/scor_assoc < $script_path/generate_SCORE_plots.R
# Generate plots for distribution of risk allele score and effect of risk factor on outcome
#$Rr --vanilla --slave --args $pheno_path $temp_path $pheno_path/covars.txt $temp_path/scor_assoc < $script_path/generate_WscorevsOutcomes_plots.R 

# Clean some temporary files
#mkdir $temp_path/scor_assoc/logs
#mv $temp_path/scor_assoc/*.summary.txt $temp_path/scor_assoc/logs
#rm $temp_path/scor_assoc/${input2}_Bin* $temp_path/scor_assoc/${input2}_QT* $temp_path/scor_assoc/${input2}_SCORsummary* $temp_path/scor_assoc/*.sh
#rm $temp_path/scor_assoc/*summary.* $temp_path/scor_assoc/*.linear $temp_path/scor_assoc/*.logistic 
#rm $temp_path/tmp* $temp_path/input_dta*
cd $temp_path/scor_assoc
zip -r ${input2}.MR.results.zip *phewas.txt plots

echo "Please find the summary results: $temp_path/scor_assoc" | mutt -s "Done SCORE based PheWAS analysis" $input3

