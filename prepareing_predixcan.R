options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  script.path <-args[1]
  predictDB.path <- args[2]
  work.path <- args[3]
} else {
  stop("Not enough arguments")
}


DB=c("TW_Adipose_Subcutaneous_0.5",
"TW_Adipose_Visceral_Omentum_0.5",
"TW_Adrenal_Gland_0.5",
"TW_Artery_Aorta_0.5",
"TW_Artery_Coronary_0.5",
"TW_Artery_Tibial_0.5",
"TW_Brain_Amygdala_0.5",
"TW_Brain_Anterior_cingulate_cortex_BA24_0.5",
"TW_Brain_Caudate_basal_ganglia_0.5",
"TW_Brain_Cerebellar_Hemisphere_0.5",
"TW_Brain_Cerebellum_0.5",
"TW_Brain_Cortex_0.5",
"TW_Brain_Frontal_Cortex_BA9_0.5",
"TW_Brain_Hippocampus_0.5",
"TW_Brain_Hypothalamus_0.5",
"TW_Brain_Nucleus_accumbens_basal_ganglia_0.5",
"TW_Brain_Putamen_basal_ganglia_0.5",
"TW_Brain_Spinal_cord_cervical_c-1_0.5",
"TW_Brain_Substantia_nigra_0.5",
"TW_Breast_Mammary_Tissue_0.5",
"TW_Cells_EBV-transformed_lymphocytes_0.5",
"TW_Cells_Transformed_fibroblasts_0.5",
"TW_Colon_Sigmoid_0.5",
"TW_Colon_Transverse_0.5",
"TW_Esophagus_Gastroesophageal_Junction_0.5",
"TW_Esophagus_Mucosa_0.5",
"TW_Esophagus_Muscularis_0.5",
"TW_Heart_Atrial_Appendage_0.5",
"TW_Heart_Left_Ventricle_0.5",
"TW_Liver_0.5",
"TW_Lung_0.5",
"TW_Minor_Salivary_Gland_0.5",
"TW_Muscle_Skeletal_0.5",
"TW_Nerve_Tibial_0.5",
"TW_Ovary_0.5",
"TW_Pancreas_0.5",
"TW_Pituitary_0.5",
"TW_Prostate_0.5",
"TW_Skin_Not_Sun_Exposed_Suprapubic_0.5",
"TW_Skin_Sun_Exposed_Lower_leg_0.5",
"TW_Small_Intestine_Terminal_Ileum_0.5",
"TW_Spleen_0.5",
"TW_Stomach_0.5",
"TW_Testis_0.5",
"TW_Thyroid_0.5",
"TW_Uterus_0.5",
"TW_Vagina_0.5",
"TW_Whole_Blood_0.5")

#	 /GWD/appbase/projects/RD-TSci-UKB/500k_GWASdownload/DataQC/sample_id
#	 /GWD/appbase/projects/RD-TSci-PhewasUKBB/ukb20361/Nov2016/predicted_exp_1KG
	 
#	 script.path="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS/script500k"
#	 predictDB.path="/GWD/appbase/projects/RD-TSci-PhewasUKBB/ukb20361/Nov2016/PredictDB_1KG"

     # read the list of genename(s)
     gene.list <- as.vector(as.matrix(read.table(paste(work.path,"/gene_names.txt",sep=""), header=FALSE)))
     # read the gene annotation to convert gene names to ENSEMBEL IDs
	 anno.gene <-read.delim(paste(script.path,"/genenameChr.txt",sep=""),header=T)

  all=NULL
  for(i in DB){
    x=read.csv(paste(predictDB.path,"/",i,".csv",sep=""), header=T)
    x=x[,c("gene","genename","rsid")]
    x$tissue=as.character(i)
    # obtain the list of SNPs according to the gene names from predictDB
    x=subset(x, genename %in% gene.list)
#   if(nrow(x) >0) print(i)
   if(is.null(all)){
      all=x
   }else{
      all=rbind(all,x)
   }
}


tmp=all[!duplicated(all$rsid),]
write.table(tmp$rsid,paste(work.path,"/snp_names.txt",sep=""),row.names=F,quote=F,col.names=F,sep="\t")
tmp=all[!duplicated(all$tissue),]
write.table(tmp$tissue,paste(work.path,"/tissue_names.txt",sep=""),row.names=F,quote=F,col.names=F,sep="\t")
tmp=all[!duplicated(all$gene),c("gene","genename")]
tmp=merge(tmp,anno.gene,by="gene")
tmp=tmp[,c("chrom","gene")]
write.table(tmp,paste(work.path,"/predicted_exp/genelist.txt",sep=""),row.names=F,quote=F,col.names=F,sep="\t")
   


    
 