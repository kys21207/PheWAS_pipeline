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
input_gene_snp=gene_snp_names.csv
sub_set_genes=gene_names.txt
sub_set_snps=snp_names.txt
ind_pheno_name=test.txt

## Set path to PLINK binary, may need >=1.9
#PLINK=/GWD/appbase/projects/GXapp/plink/plink-1.9/plink
Rr=/home/anaconda/anaconda2/bin/R

# Combine all results to make a summary table

$Rr --vanilla --slave --args $temp_path $temp_path/$input_gene_snp  $pheno_path < $script_path/combine_SNP_results.R

# Create a directory for plots= 'plots' 
mkdir -p $temp_path/snps_assoc/plots

# Generate plots for all SNPs across phenotypes
$Rr --vanilla --slave --args $temp_path/snps_assoc < $script_path/generate_SNP_plots.R 

# Clean some temporary files
#rm $temp_path/snps_assoc/${input2}_Bin* $temp_path/snps_assoc/${input2}_QT* $temp_path/snps_assoc/${input2}_SNPsummary* $temp_path/snps_assoc/*.sh
mkdir $temp_path/snps_assoc/logs
#mv $temp_path/snps_assoc/*.summary.* $temp_path/snps_assoc/logs
#rm $temp_path/snps_assoc/*.linear $temp_path/snps_assoc/*.logistic $temp_path/snps_assoc/*summary.* 
#rm $temp_path/snps_assoc/${input2}_Bin* $temp_path/snps_assoc/${input2}_QT*  $temp_path/snps_assoc/${input2}.Bin* $temp_path/snps_assoc/${input2}.QT* 
#rm $temp_path/test.* $temp_path/tmp* $temp_path/input_dta* 
cd $temp_path/snps_assoc
zip -r ${input2}.snps.results.zip *phewas.txt plots

echo "Please find the summary results: $temp_path/snp_assoc" | mutt -s "Done SNP based PheWAS analysis" $input3

