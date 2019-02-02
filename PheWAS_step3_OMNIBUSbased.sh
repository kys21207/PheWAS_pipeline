#!/bin/bash
#set -o errexit

# assigned random user name 
input2=$1
# email address
input3=$2
input1=$3

if [ $# -lt 3 ]; then
  echo 1>&2 "$0: not enough arguments"
  exit 2
elif [ $# -gt 3 ]; then
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
set_genes=gene_names.txt
set_snps=snp_names.txt
ind_pheno_name=test.txt

## Set path to PLINK binary, may need >=1.9
#PLINK=/GWD/appbase/projects/GXapp/plink/plink-1.9/plink
Rr=/home/anaconda/anaconda2/bin/R

# Combine all results to make a summary table

$Rr --vanilla --slave --args $temp_path $pheno_path < $script_path/combine_OMNIBUS_results.R

# Create a directory for plots= 'plots' 
mkdir -p $temp_path/gene_assoc/plots

# Generate plots for all SNPs across phenotypes
$Rr --vanilla --slave --args $temp_path/gene_assoc < $script_path/generate_OMNIBUS_plots.R 
$Rr --vanilla --slave --args $temp_path/gene_assoc < $script_path/generate_univariable_plots.R

# Clean some temporary files
#rm -r $temp_path/tmp_predict $temp_path/predicted_exp/d*
mkdir $temp_path/gene_assoc/logs
mv $temp_path/gene_assoc/*.summary.txt $temp_path/gene_assoc/logs
#rm $temp_path/gene_assoc/${input2}_Bin* $temp_path/gene_assoc/${input2}_QT*  $temp_path/gene_assoc/${input2}.Bin* $temp_path/gene_assoc/${input2}.QT* $temp_path/gene_assoc/*.sh
#rm $temp_path/gene_assoc/*.linear $temp_path/gene_assoc/*.logistic $temp_path/gene_assoc/*_run* $temp_path/gene_assoc/*_prepareExp* $temp_path/gene_assoc/*_Genesummary*
#rm $temp_path/tmp* $temp_path/input_dta* $temp_path/test.log
cd $temp_path/gene_assoc
zip -r ${input2}.predixcan.results.zip *phewas.txt plots

mv $temp_path/gene_assoc $temp_path/gene_assoc${input1}
echo "Please find the summary results: $temp_path/gene_assoc" | mutt -s "Done Gene based PheWAS analysis" ${input3}

