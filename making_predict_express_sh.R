options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 4) {
  script.path <-args[1]
  predictDB.path <- args[2]
  sampleID <- args[3]
  work.path <- args[4]
} else {
  stop("Not enough arguments")
}



#predictDB_path="/GWD/appbase/projects/RD-TSci-PhewasUKBB/ukb20361/Nov2016/PredictDB_1KG"
#expres_path="/GWD/appbase/projects/RD-TSci-PhewasUKBB/ukb20361/Nov2016/predicted_exp_1KG"
#fam_path="/GWD/appbase/projects/RD-TSci-PhewasUKBB/ukb20361/Nov2016/dosages_500k"

arry_nameDB <- as.vector(as.matrix(read.table(paste(work.path,"/tissue_names.txt",sep=""),header=F)))

j=1
for (i in arry_nameDB){
   # calculate prediction only
   sink(paste(work.path,"/tmp_predict/run_string_",j,".sh",sep=""))
   cat(paste("#!/bin/bash"),"\n")
   cat(paste("mkdir ",work.path,"/predicted_exp/d",j,sep=""),"\n")
   cat(paste("/GWD/bioinfo/projects/RD-TSci-Software/CB/linuxbrew/bin/python2.7 ",script.path,"/PrediXcan.py --predict --dosages ",work.path,"/dosages --samples ",sampleID," --genelist ",work.path,"/predicted_exp/genelist.txt --weights ",predictDB.path,"/",i,".db --output_dir ",work.path,"/predicted_exp/d",j,sep=""),"\n") 
   cat(paste("gzip ",work.path,"/predicted_exp/d",j,"/predicted_expression.txt",sep=""),"\n") 
   cat(paste("mv ",work.path,"/predicted_exp/d",j,"/predicted_expression.txt.gz ",work.path,"/predicted_exp/",i,".txt.gz",sep=""),"\n")
   sink()
   j=j+1
}