#!/bin/bash

#Version 1
#set -o errexit

#select 'snp' or 'gene' or 'score'
input4=$1
#full path for input data
input0=$2
# phenotype selection
input1=$3
# specified user name 
input2=$4
# email address
input3=$5


if [ $# -lt 5 ]; then
  echo 1>&2 "$0: not enough arguments"
  exit 2
elif [ $# -gt 5 ]; then
  echo 1>&2 "$0: too many arguments"
  exit 2
fi

## Set path to phenotype and genotype dataset
pheno_path=/GWD/appbase/projects/RD-TSci-UKB/PheWAS/phenotype500kV3
sampleID_path=/GWD/appbase/projects/RD-TSci-UKB/500k_GWASdownload/DataQC/sample_id
output_path=/GWD/appbase/projects/RD-TSci-UKB/PheWAS_output_storage
script_path=/GWD/appbase/projects/RD-TSci-UKB/PheWAS/script500k
record_path=/GWD/appbase/projects/RD-TSci-UKB/PheWAS/Analysis_recordCenter
predictDB_path=/GWD/appbase/projects/RD-TSci-UKB/PredictDB

sub_sample=impv3
HES_Self_cate_QT=Anno_HESelf_QT.csv
HES_Self_cate_BI=Anno_HESelf_BI.csv

## Set path to PLINK binary, may need >=1.9
#PLINK=/GWD/appbase/projects/GXapp/plink/plink-1.9/plink
Rr=/home/anaconda/anaconda2/bin/R

#create a user folder for performance 
rnum=$(date +%Y%m%d%H%M)
input2=${input2}_${rnum}
temp_path=$output_path/$input2

mkdir -p $temp_path

cp $input0 $temp_path/input_dta0.txt
# remove blank lines 
sed '/^$/d' $temp_path/input_dta0.txt > $temp_path/input_dta.txt
rm $temp_path/input_dta0.txt

$Rr --vanilla --slave --args $input4 $temp_path/input_dta.txt $temp_path < $script_path/checking_inputData.R > $temp_path/test.log

error=$(head -n 1 $temp_path/test.log)
if [ $error == "Error" ] 
then 
  echo "Please check your input file"
  mail -s "Error: please check your input file" $input3
  exit 2
fi   
 
# Set file names for phenotype(s), covariate(s) and the subset of genes
covar_name=covars.txt
input_gene_snp=gene_snp_names.csv
set_genes=gene_names.txt
set_snps=snp_names.txt
ind_pheno_name=test.txt


# Generate a log file for recording 
$Rr --vanilla --slave --args $input4 $input0 $input1 $input3 $temp_path $record_path < $script_path/make_logFile.R

# Generate phenotype files for Qt & Bin
#$Rr --vanilla --slave --args linear $input1 $pheno_path $temp_path < $script_path/Generate_SubPhenotypes.R
#$Rr --vanilla --slave --args logistic $input1 $pheno_path $temp_path < $script_path/Generate_SubPhenotypes.R

# Copy the covariate file into temp directory
#cp $pheno_path/$covar_name $temp_path/$covar_name
#cp $pheno_path/phenotypes_BI.txt $temp_path/
#cp $pheno_path/phenotypes_QT.txt $temp_path/


if [[ $input4 == "snp" ]]
then
  echo "starting to run SNP based PheWAS analysis"

  # Remove all white space or tab from the file in order to check whether the file is empty or not.
  cat $temp_path/$set_snps | tr -d " \t\n\r" >> $temp_path/tmp1.txt
  cat $temp_path/$set_genes | tr -d " \t\n\r" >> $temp_path/tmp2.txt

  # Pull out a subset of variants from bgen files and convert to gen files
  python $script_path/Batch_submission_BGEN2GEN_v7.py ${input2}_snp $temp_path/$set_snps $temp_path/temp $temp_path genotype 

  printf '%s\n' '#!/bin/bash' "$script_path/PheWAS_step2_SNPbased.sh $input2 $input3" > $temp_path/${input2}_run_step2.sh
  chmod 777 $temp_path/${input2}_run_step2.sh
 
  sbatch --dependency=singleton --job-name=${input2}_snp_genotype -e $temp_path/${input2}_run_step2.err -o $temp_path/${input2}_run_step2.out --time=250:00:00 --cpus-per-task=1 --mem=30G $temp_path/${input2}_run_step2.sh


fi


if [[ $input4 == "gene" ]]
then
  echo "starting to run Omnibus Gene based PheWAS analysis"
# Remove all white space or tab from the file in order to check whether the file is empty or not.
cat $temp_path/$set_genes | tr -d " \t\n\r" >> $temp_path/tmp2.txt
mkdir -p $temp_path/predicted_exp
mkdir -p $temp_path/dosages
mkdir -p $temp_path/tmp_predict

# Obtain the list of SNPs from PredictDB 
$Rr --vanilla --slave --args $script_path $predictDB_path $temp_path < $script_path/prepareing_predixcan.R

# Pull out a subset of variants from bgen files and convert to gen files
python $script_path/Batch_submission_BGEN2GEN_v7.py ${input2}_gene $temp_path/$set_snps $temp_path/temp $temp_path genotype 


#    ListJobId=$(eval echo predict_{1..$j} | tr ' ' , )
    printf '%s\n' '#!/bin/bash' "$script_path/PheWAS_step2_OMNIBUSbased.sh $input2 $input3" > $temp_path/${input2}_run_step2.sh
    chmod 777 $temp_path/${input2}_run_step2.sh
    sbatch --dependency=singleton --job-name=${input2}_gene_genotype -e $temp_path/${input2}_run_step2.err -o $temp_path/${input2}_run_step2.out --time=250:00:00 --cpus-per-task=1 --mem=30G $temp_path/${input2}_run_step2.sh


fi


if [[ $input4 == "score" ]]
then
# Pull out a subset of variants from bgen files and convert to gen files
python $script_path/Batch_submission_BGEN2GEN_v7.py ${input2}_score $temp_path/$set_snps $temp_path/temp $temp_path genotype 

  printf '%s\n' '#!/bin/bash' "$script_path/PheWAS_step2_SCORbased.sh $input2 $input3" > $temp_path/${input2}_run_step2.sh
  chmod 777 $temp_path/${input2}_run_step2.sh
  sbatch --dependency=singleton --job-name=${input2}_score_genotype -e $temp_path/${input2}_run_step2.err -o $temp_path/${input2}_run_step2.out --time=250:00:00 --cpus-per-task=1 --mem=30G $temp_path/${input2}_run_step2.sh

fi
