#!/bin/bash

# assigned random user name 
input2=$1
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
sampleID_path=/GWD/appbase/projects/RD-TSci-UKB/500k_GWASdownload/DataQC/sample_id

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
sub_sample=impv3

## Set path to PLINK binary, may need >=1.9
#PLINK=/GWD/appbase/projects/GXapp/plink/plink-1.9/plink
Rr=/home/anaconda/anaconda2/bin/R

cat $temp_path/temp/*.gen > $temp_path/genotype.gen
#1. concatenate the extracted files - for example

# create header for the gen file
echo chr chrpos snp pos otherall effall > $temp_path/gen.header
# cut files and create map file
cut -d " " -f 1-6 $temp_path/genotype.gen > $temp_path/genotype_nohead.map
cat $temp_path/gen.header $temp_path/genotype_nohead.map > $temp_path/genotype.map
rm $temp_path/genotype_nohead.map $temp_path/gen.header

# calculate allele B dosage (allele B, the 2nd allele, is the effect allele) from imputed gene file (prob_AB+2*prob_BB, ie, 2nd column plus 2 times the 3rd column for every 3 columns). replace dosage = "NaN" if sum(probabilities)<0.9
cat $temp_path/genotype.gen | awk '{printf $3 " "; for(i=7;i<=NF;i=i+3) {if ($i+$(i+1)+$(i+2)<0.9) printf "NaN "; else printf $(i+1)+2*$(i+2) " "} printf "\n"}' > $temp_path/tempf0

# then remove end space and transpose
sed 's/ *$//' $temp_path/tempf0 | sed 's/ /,/g' > $temp_path/tempf1
#change the line below to my location for the script .pl
perl $script_path/transpose_comma.pl -i $temp_path/tempf1 -o $temp_path/tempf0
sed 's/,/ /g' $temp_path/tempf0 > $temp_path/tempf1

# merge ID from file impv1.sample
awk '(NR!=2) {print $1}' ${sampleID_path}/${sub_sample}.sample > $temp_path/tempf0
paste $temp_path/tempf0 $temp_path/tempf1 | sed 's/\t/ /g' > $temp_path/genotype.dose
rm $temp_path/tempf0 $temp_path/tempf1 

$Rr --vanilla --slave --args $temp_path < $script_path/finding_missing_snps.R

  if [ -s "$temp_path/tmp1.txt" ]   
  then
    ## run the SNP based association test
    # Create a directory for SNP based association results = 'snps_assoc' 
    mkdir -p $temp_path/snps_assoc

 
   #distribute jobs in terms of # of phenotypes
    numFiles=$(ls $pheno_path/Sub_BI_pheno_* | wc -l )

    for((k=1; k <= numFiles; k++))
      do
        # run a logistic association test
        printf '%s\n' '#!/bin/bash' "$Rr --vanilla --slave --args $pheno_path/Sub_BI_pheno_${k}.txt $pheno_path/$covar_name $temp_path/genotype $temp_path/snps_assoc ${k} < $script_path/dose2binomial_snp.R" > $temp_path/snps_assoc/${input2}_run_BinSNP_${k}.sh
        chmod 777 $temp_path/snps_assoc/${input2}_run_BinSNP_${k}.sh
        #run by GRID
        #qsub -N ${input2}_BinSNP_${k} -q "rhel7" -o $temp_path/snps_assoc -e $temp_path/snps_assoc -l mem_free=20G,hostname="us1salx00635|us1salx00649|us1salx00792|us1salx00791|us1salx00843|us1salx00844|us1salx00845|us1salx00942|us1salx00943" $temp_path/snps_assoc/${input2}_run_BinSNP_${k}.sh  
        sbatch --job-name=${input2}_AssoSNP -e $temp_path/snps_assoc/${input2}_BinSNP_${k}.err -o $temp_path/snps_assoc/${input2}_BinSNP_${k}.out --time=250:00:00 --cpus-per-task=10 --mem=30G $temp_path/snps_assoc/${input2}_run_BinSNP_${k}.sh
        #ListJobId[k]=$(eval echo ${id} : | tr -d 'Submitted batch job') 
      done
    
   # run a linear association test
    printf '%s\n' '#!/bin/bash' "$Rr --vanilla --slave --args $pheno_path/$QT_pheno_name $pheno_path/$covar_name $temp_path/genotype $temp_path/snps_assoc < $script_path/dose2lm_snp.R" > $temp_path/snps_assoc/${input2}_run_QTSNP.sh
    chmod 777 ${temp_path}/snps_assoc/${input2}_run_QTSNP.sh
    #run by GRID
    #qsub -N ${input2}_QTSNP -o $temp_path/snps_assoc -e $temp_path/snps_assoc -q "rhel7" -l mem_free=20G,hostname="us1salx00635|us1salx00649|us1salx00792|us1salx00791|us1salx00843|us1salx00844|us1salx00845|us1salx00942|us1salx00943" $temp_path/snps_assoc/${input2}_run_QTSNP.sh  
    sbatch --job-name=${input2}_AssoSNP -e $temp_path/snps_assoc/${input2}_QTSNP.err -o $temp_path/snps_assoc/${input2}_QTSNP.out --time=250:00:00 --cpus-per-task=10 --mem=30G $temp_path/snps_assoc/${input2}_run_QTSNP.sh

   printf '%s\n' '#!/bin/bash' "$script_path/PheWAS_step3_SNPbased.sh $input2 $input3" > $temp_path/snps_assoc/${input2}_run_SNPsummary.sh
    chmod 777 $temp_path/snps_assoc/${input2}_run_SNPsummary.sh
    
  #qsub -hold_jid $ListJobId -N ${input2}_SNPsummary -o $temp_path/snps_assoc -e $temp_path/snps_assoc -q "rhel7" -l mem_free=20G,hostname="us1salx00635|us1salx00649|us1salx00792|us1salx00791|us1salx00843|us1salx00844|us1salx00845|us1salx00942|us1salx00943" $temp_path/snps_assoc/${input2}_run_SNPsummary.sh
  #sbatch --dependency=afterok:$(echo ${ListJobId[@]} | tr -d ' ') --job-name=${input2}_SNPsummary -e $temp_path/snps_assoc/${input2}_SNPsummary.err -o $temp_path/snps_assoc/${input2}_SNPsummary.out --time=25:00:00 --cpus-per-task=1 --mem=2000 $temp_path/snps_assoc/${input2}_run_SNPsummary.sh
  sbatch --dependency=singleton --job-name=${input2}_AssoSNP -e $temp_path/snps_assoc/${input2}_SNPsummary.err -o $temp_path/snps_assoc/${input2}_SNPsummary.out --time=25:00:00 --cpus-per-task=1 --mem=20G $temp_path/snps_assoc/${input2}_run_SNPsummary.sh

  else
    echo "Please check your input file or you skiped the SNP based assocation."
  fi

