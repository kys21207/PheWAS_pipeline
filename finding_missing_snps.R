options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  work.path <- args[1]
} else {
  stop("Not enough arguments")
}

#work.path="/GWD/appbase/projects/RD-TSci-UKB/PheWAS_output_storage/Cambridge_snps_batch1_201805301159"
map.list <- read.table(paste(work.path,"genotype.map",sep="/"),header=T)
snp.list <- read.table(paste(work.path,"snp_names.txt",sep="/"),header=F)

missing.snp <- snp.list$V1[!(snp.list$V1 %in% map.list$snp)]
if(length(missing.snp) > 0){
   write.table(missing.snp,paste(work.path,"missing_SNPs.txt",sep="/"),row.names=F,col.names=F,quote=F,sep="\t")
}else{
   sink(paste(work.path,"missing_SNPs.txt",sep="/"))
   print("NO MISSING SNP")
   sink()
}



