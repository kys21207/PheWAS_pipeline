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
echo chr snp pos otherall effall eaf> $temp_path/gen.header

# calculate allele B dosage (allele B, the 2nd allele, is the effect allele) from imputed gene file (prob_AB+2*prob_BB, ie, 2nd column plus 2 times the 3rd column for every 3 columns). replace dosage = "NaN" if sum(probabilities)<0.9
cat $temp_path/genotype.gen | awk '{printf $1 " " $3 " " $4 " " $5 " " $6 " ";  for(i=7;i<=NF;i=i+3) {if ($i+$(i+1)+$(i+2)<0.9) printf "NaN "; else printf $(i+1)+2*$(i+2) " "} printf "\n"}' > $temp_path/tempf01
# for 
cat $temp_path/tempf01 | awk '{T=0; for(N=6; N<=NF; N++) {T+=$N}; T/=((NF-5)*2); printf $1 " " $2 " " $3 " " $4 " " $5 " " T " ";  for(i=6;i<=NF;i++) {printf $i " "} printf "\n"}' > $temp_path/tempf02
# cut files and create map file
cut -d " " -f 1-6 $temp_path/tempf02 > $temp_path/genotype_nohead.map
cat $temp_path/gen.header $temp_path/genotype_nohead.map > $temp_path/genotype.map
rm $temp_path/genotype_nohead.map $temp_path/gen.header
cat $temp_path/tempf02 | awk '{printf $2 " ";  for(i=7;i<=NF;i++) {printf $i " "} printf "\n"}' > $temp_path/tempf0

# then remove end space and transpose
sed 's/ *$//' $temp_path/tempf0 | sed 's/ /,/g' > $temp_path/tempf1
#change the line below to my location for the script .pl
perl $script_path/transpose_comma.pl -i $temp_path/tempf1 -o $temp_path/tempf0
sed 's/,/ /g' $temp_path/tempf0 > $temp_path/tempf1

# merge ID from file impv1.sample
awk '(NR!=2) {print $1}' ${sampleID_path}/${sub_sample}.sample > $temp_path/tempf0
paste $temp_path/tempf0 $temp_path/tempf1 | sed 's/\t/ /g' > $temp_path/genotype.dose
rm $temp_path/tempf0 $temp_path/tempf1 $temp_path/tempf01 $temp_path/tempf02

$Rr --vanilla --slave --args $temp_path < $script_path/finding_missing_snps.R

#$Rr --vanilla --slave --args $temp_path/score_data.txt $temp_path/genotype $temp_path < $script_path/MR_dosage2Uweightedscore.R
$Rr --vanilla --slave --args $temp_path/score_data.txt $temp_path/genotype $temp_path < $script_path/MR_dosage2Weightedscore.R

    mkdir -p $temp_path/scor_assoc
#    mkdir -p $temp_path/scor_assoc/distr_plots
#    $Rr --vanilla --slave --args $temp_path $temp_path/scor_assoc/distr_plots < $script_path/generate_UWvsWscore_plots.R



    numFiles=$(ls $pheno_path/Sub_BI_pheno_* | wc -l )

    for((k=1; k <= numFiles; k++))
      do
       # run a logistic association test
       printf '%s\n' '#!/bin/bash' "$Rr --vanilla --slave --args $pheno_path/Sub_BI_pheno_${k}.txt $pheno_path/$covar_name $temp_path/genotype $temp_path/scor_assoc $temp_path/scor_names.txt ${k} < $script_path/dose2binomial_score.R" > $temp_path/scor_assoc/${input2}_run_BinSCOR_${k}.sh
       chmod 777 $temp_path/scor_assoc/${input2}_run_BinSCOR_${k}.sh
       #run by GRID
       #qsub -N ${input2}_BinSCOR_${k} -q "rhel7" -o $temp_path/scor_assoc -e $temp_path/scor_assoc -l mem_free=20G,hostname="us1salx00635|us1salx00649|us1salx00792|us1salx00791|us1salx00843|us1salx00844|us1salx00845|us1salx00942|us1salx00943" $temp_path/scor_assoc/${input2}_run_BinSCOR_${k}.sh  
       sbatch --job-name=${input2}_AssoSCOR -e $temp_path/scor_assoc/${input2}_BinSCOR_${k}.err -o $temp_path/scor_assoc/${input2}_BinSCOR_${k}.out --time=250:00:00 --cpus-per-task=10 --mem=30G $temp_path/scor_assoc/${input2}_run_BinSCOR_${k}.sh

    done
    
 
   # run a linear association test
    printf '%s\n' '#!/bin/bash' "$Rr --vanilla --slave --args $pheno_path/$QT_pheno_name $pheno_path/$covar_name $temp_path/genotype $temp_path/scor_assoc $temp_path/scor_names.txt < $script_path/dose2lm_score.R" > $temp_path/scor_assoc/${input2}_run_QTSCOR.sh
    chmod 777 ${temp_path}/scor_assoc/${input2}_run_QTSCOR.sh
    #run by GRID
    #qsub -N ${input2}_QTSCOR -o $temp_path/scor_assoc -e $temp_path/scor_assoc -q "rhel7" -l mem_free=20G,hostname="us1salx00635|us1salx00649|us1salx00792|us1salx00791|us1salx00843|us1salx00844|us1salx00845|us1salx00942|us1salx00943" $temp_path/scor_assoc/${input2}_run_QTSCOR.sh  
    sbatch --job-name=${input2}_AssoSCOR -e $temp_path/scor_assoc/${input2}_QTSCOR.err -o $temp_path/scor_assoc/${input2}_QTSCOR.out --time=250:00:00 --cpus-per-task=10 --mem=30G $temp_path/scor_assoc/${input2}_run_QTSCOR.sh


    printf '%s\n' '#!/bin/bash' "$script_path/PheWAS_step3_SCORbased.sh $input2 $input3" > $temp_path/scor_assoc/${input2}_run_SCORsummary.sh
    chmod 777 $temp_path/scor_assoc/${input2}_run_SCORsummary.sh

  #ListJobId=$(eval echo ${input2}_BinSCOR_{1..$numFiles} ${input2}_QTSCOR | tr ' ' , )
  #qsub -hold_jid $ListJobId -N ${input2}_SCORsummary -o $temp_path/scor_assoc -e $temp_path/scor_assoc -q "rhel7" -l mem_free=20G,hostname="us1salx00635|us1salx00649|us1salx00792|us1salx00791|us1salx00843|us1salx00844|us1salx00845|us1salx00942|us1salx00943" $temp_path/scor_assoc/${input2}_run_SCORsummary.sh
  sbatch --dependency=singleton --job-name=${input2}_AssoSCOR -e $temp_path/scor_assoc/${input2}_SCORsummary.err -o $temp_path/scor_assoc/${input2}_SCORsummary.out --time=25:00:00 --cpus-per-task=1 --mem=20G $temp_path/scor_assoc/${input2}_run_SCORsummary.sh

