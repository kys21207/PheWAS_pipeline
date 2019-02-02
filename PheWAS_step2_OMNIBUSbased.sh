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
output_path=/GWD/appbase/projects/RD-TSci-UKB/PheWAS_output_storage
script_path=/GWD/appbase/projects/RD-TSci-UKB/PheWAS/script500k
predictDB_path=/GWD/appbase/projects/RD-TSci-UKB/PredictDB
sampleID_path=/GWD/appbase/projects/RD-TSci-UKB/500k_GWASdownload/DataQC/sample_id

#create a user folder for performance 
foldername=$input2
temp_path=$output_path/$foldername
expres_path=$temp_path/predicted_exp
 
# Set file names for phenotype(s), covariate(s) and the subset of genes
#BI_pheno_name=phenotypes_BI.txt
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
#Qqsub=/GWD/bioinfo/projects/lsf/SGE/6.2u5/bin/lx24-amd64/qsub

ls $temp_path/temp/chr*.gen > $temp_path/temp/chrs.txt
IFS=$'\n' read -d '' -r -a arr_nameChr < $temp_path/temp/chrs.txt

for (( j=0; j < ${#arr_nameChr[@]}; j++ ))
do
  if [[ ! -s ${arr_nameChr[$j]} ]]; then rm ${arr_nameChr[$j]}; fi
done
ls $temp_path/temp/chr*.gen > $temp_path/temp/chrs.txt
sed -r 's/temp\//\t/g;s/\.gen//g' $temp_path/temp/chrs.txt | awk '{print $2}' > $temp_path/temp/chrs_gen.txt

#readarray arr_nameDB < $temp_path/temp/chrs_gen.txt
IFS=$'\n' read -d '' -r -a arr_nameChr < $temp_path/temp/chrs_gen.txt

for (( j=0; j < ${#arr_nameChr[@]}; j++ ))
 do
   # calculate allele B dosage (allele B, the 2nd allele, is the effect allele) from imputed gene file (prob_AB+2*prob_BB, ie, 2nd column plus 2 times the 3rd column for every 3 columns). replace dosage = "NaN" if sum(probabilities)<0.9
   cat $temp_path/temp/${arr_nameChr[$j]}.gen | awk '{printf $1 " " $3 " " $4 " " $5 " " $6 " ";  for(i=7;i<=NF;i=i+3) {if ($i+$(i+1)+$(i+2)<0.9) printf "NaN "; else printf $(i+1)+2*$(i+2) " "} printf "\n"}' > $temp_path/temp/tempf01
   # for 
   cat $temp_path/temp/tempf01 | awk '{T=0; for(N=6; N<=NF; N++) {T+=$N}; T/=((NF-5)*2); printf $1 " " $2 " " $3 " " $4 " " $5 " " T " ";  for(i=6;i<=NF;i++) {printf $i " "} printf "\n"}' > $temp_path/temp/tempf02
   
   sed 's/0//' $temp_path/temp/tempf02 > $temp_path/dosages/${arr_nameChr[$j]}.dose

   gzip $temp_path/dosages/${arr_nameChr[$j]}.dose

 done  

$Rr --vanilla --slave --args $script_path $predictDB_path ${sampleID_path}/${sub_sample}.fam $temp_path < $script_path/making_predict_express_sh.R
# run the script to generate the predicted expression level for gene(s)
numFiles=$(ls $temp_path/tmp_predict/run_string_*.sh | wc -l )

for((j=1; j <= numFiles; j++))
 do
   chmod 777 $temp_path/tmp_predict/run_string_${j}.sh
#   qsub -N predict_${j} -e $temp_path/tmp_predict -o $temp_path/tmp_predict -q "rhel7" -l mem_free=30G,hostname="us1salx00635|us1salx00649|us1salx00792|us1salx00791|us1salx00843|us1salx00844|us1salx00845" -cwd $temp_path/tmp_predict/run_string_${j}.sh
   sbatch --job-name=${input2}_predict -e $temp_path/tmp_predict/run_string_${j}.err -o $temp_path/tmp_predict/run_string_${j}.out --time=250:00:00 --cpus-per-task=10 --mem=30G $temp_path/tmp_predict/run_string_${j}.sh
 done

if [ -s "$temp_path/tmp2.txt" ] 
then
    ## run the Gene based association test
    # Create a directory for gene based association results = 'gene_assoc' 
    mkdir -p $temp_path/gene_assoc

## select all tissues related to the list of genes 
#$Rr --vanilla --slave --args $pheno_path $temp_path < $script_path/GeneList_filtering.R

    #prepare a expression data and info  
    printf '%s\n' '#!/bin/bash' "$Rr --vanilla --slave --args $predictDB_path  $expres_path $temp_path  < $script_path/prepare_omnibus_gene.R" > $temp_path/gene_assoc/${input2}_run_prepareExp.sh
    chmod 777 $temp_path/gene_assoc/${input2}_run_prepareExp.sh
#    $Qqsub -N ${input2}_prepareExp  -o $temp_path/gene_assoc -e $temp_path/gene_assoc -q "rhel7" -l mem_free=30G,hostname="us1salx00635|us1salx00649|us1salx00792|us1salx00791|us1salx00843|us1salx00844|us1salx00845|us1salx00942|us1salx00943" $temp_path/gene_assoc/${input2}_run_prepareExp.sh  
    id=$(sbatch --dependency=singleton --job-name=${input2}_predict -e $temp_path/gene_assoc/${input2}_prepareExp.err -o $temp_path/gene_assoc/${input2}_prepareExp.out --time=250:00:00 --cpus-per-task=10 --mem=30G $temp_path/gene_assoc/${input2}_run_prepareExp.sh)
    id=$(eval echo ${id}  | tr -d 'Submitted batch job')

   #distribute jobs in terms of # of phenotypes
    numFiles=$(ls $pheno_path/Sub_BI_pheno_* | wc -l )

    for((k=1; k <= numFiles; k++))
     do
        #run a logistic association test
        # Create a sub-directory to store the association results in terms of a specific tissue
        printf '%s\n' '#!/bin/bash' "$Rr --vanilla --slave --args $pheno_path/Sub_BI_pheno_${k}.txt $pheno_path/$covar_name $temp_path $temp_path/gene_assoc $k < $script_path/omnibus4binomial_gene.R" > $temp_path/gene_assoc/${input2}_run_BinGene_${k}.sh
        chmod 777 $temp_path/gene_assoc/${input2}_run_BinGene_${k}.sh
#        $Qqsub -hold_jid ${input2}_prepareExp -N ${input2}_BinGene_${k}  -o $temp_path/gene_assoc -e $temp_path/gene_assoc -q "rhel7" -l mem_free=30G,hostname="us1salx00635|us1salx00649|us1salx00792|us1salx00791|us1salx00843|us1salx00844|us1salx00845|us1salx00942|us1salx00943" $temp_path/gene_assoc/${input2}_run_BinGene_${k}.sh  
        sbatch --dependency=afterany:$(echo ${id[@]} | tr -d ' ') --job-name=${input2}_AssoGene -e $temp_path/gene_assoc/${input2}_BinGene_${k}.err -o $temp_path/gene_assoc/${input2}_BinGene_${k}.out --time=250:00:00 --cpus-per-task=10 --mem=30G $temp_path/gene_assoc/${input2}_run_BinGene_${k}.sh

    done

    #run a linear association test
    # Create a sub-directory to store the association results in terms of a specific tissue
    printf '%s\n' '#!/bin/bash' "$Rr --vanilla --slave --args $pheno_path/$QT_pheno_name $pheno_path/$covar_name $temp_path $temp_path/gene_assoc < $script_path/omnibus4lm_gene.R" > $temp_path/gene_assoc/${input2}_run_QTGene.sh
    chmod 777 $temp_path/gene_assoc/${input2}_run_QTGene.sh
#    $Qqsub -hold_jid ${input2}_prepareExp -N ${input2}_QTGene  -o $temp_path/gene_assoc -e $temp_path/gene_assoc -q "rhel7" -l mem_free=30G,hostname="us1salx00635|us1salx00649|us1salx00792|us1salx00791|us1salx00843|us1salx00844|us1salx00845|us1salx00942|us1salx00943" $temp_path/gene_assoc/${input2}_run_QTGene.sh  
    sbatch --dependency=afterany:$(echo ${id[@]} | tr -d ' ') --job-name=${input2}_AssoGene -e $temp_path/gene_assoc/${input2}_QTGene.err -o $temp_path/gene_assoc/${input2}_QTGene.out --time=250:00:00 --cpus-per-task=10 --mem=30G $temp_path/gene_assoc/${input2}_run_QTGene.sh

    printf '%s\n' '#!/bin/bash' "$script_path/PheWAS_step3_OMNIBUSbased.sh $input2 $input3" > $temp_path/gene_assoc/${input2}_run_OMNIBUSsummary.sh
    chmod 777 $temp_path/gene_assoc/${input2}_run_OMNIBUSsummary.sh
#  ListJobId=$(eval echo ${input2}_BinGene_{1..$numFiles} ${input2}_QTGene | tr ' ' , )
#  $Qqsub -hold_jid $ListJobId -N ${input2}_Genesummary -o $temp_path/gene_assoc -e $temp_path/gene_assoc -q "rhel7" -l mem_free=20G,hostname="us1salx00635|us1salx00649|us1salx00792|us1salx00791|us1salx00843|us1salx00844|us1salx00845|us1salx00942|us1salx00943" $temp_path/gene_assoc/${input2}_run_OMNIBUSsummary.sh
  sbatch --dependency=singleton --job-name=${input2}_AssoGene -e $temp_path/gene_assoc/${input2}_Genesummary.err -o $temp_path/gene_assoc/${input2}_Genesummary.out --time=25:00:00 --cpus-per-task=1 --mem=30G $temp_path/gene_assoc/${input2}_run_OMNIBUSsummary.sh
else
   echo "Please check your input file or you skiped the OMNIBUS_GENE based assocation."
fi
