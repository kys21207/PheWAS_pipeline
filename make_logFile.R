 options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 6) {
  analysis.type <- args[1]
  input.dta <- args[2]
  pheno.dta <- args[3]
  email.add <- args[4]
  output.path<-args[5]
  record.path<-args[6]
} else {
  stop("Not enough arguments")
}
#email.add="kys21207@gsk.com"
#analysis.type="snp"
#input.dta="/GWD/appbase/projects/RD-TSci-PhewasUKB/KJ/input_dta.txt"
#pheno.dta="ALL"
#output.path<-"/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/Fanli_201706201359"
#record.path<-"/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS/Analysis_recordCenter"

pipeline.file <- "/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS/script/PheWAS_run.sh"
v <- readLines(pipeline.file)
today<-Sys.Date()
dta<-read.table(input.dta,header=T)
if(analysis.type=="snp"){
  type="SNP based analysis"
}else if(analysis.type=="gene" | analysis.type=="omnibus"){
  type="Gene based analysis"
}else if(analysis.type=="score"){
  dta <- dta[,c("trait","snp")]
  type="MR_Score based analysis"
}
if(pheno.dta=="ALL"){
Binary.phenotypes<- "/GWD/appbase/projects/RD-TSci-UKB/PheWAS/phenotype500kNew/phenotypes_BI.txt"
Quanti.phenotypes<-"/GWD/appbase/projects/RD-TSci-UKB/PheWAS/phenotype500kNew/phenotypes_QT.txt"
}
sink(paste(output.path,"/",sub('\\..*$', '', basename(output.path)),".log",sep=""))
cat("==================== \n")
cat("PheWAS Study Summary \n")
cat("==================== \n")
cat("\n")
cat(paste("USER            : ",Sys.info()[["user"]],"\n",sep=""))
cat(paste("Email           : ",email.add,"\n",sep=""))
cat(paste("DataAsOf        : ",format(today,format="%B %d %Y"),"\n",sep=""))
cat(paste("Pipeline Version: ",gsub("#","",v[[2]]),"\n",sep=""))
cat(paste("Pipeline build  : ",file.info(pipeline.file)$ctime,"\n",sep=""))
cat(paste("Analysis Type   : ",type,"\n",sep=""))
cat("\n")
cat("========================= \n")
cat("List of SNP/Gene interest \n")
cat("========================= \n")
cat("\n")
dta
cat("\n")
cat("====================================== \n")
cat("Phenotypes included in PheWAS analysis \n")
cat("====================================== \n")
cat("\n")
cat(paste(basename(Binary.phenotypes),": ",file.info(Binary.phenotypes)$ctime,"\n",sep=""))
cat(paste(basename(Quanti.phenotypes),": ",file.info(Quanti.phenotypes)$ctime,"\n",sep=""))
cat("\n")
cat("====================================== \n")
cat("Output Link\n")
cat("====================================== \n")
cat("\n")
cat(paste(output.path,"\n",sep=""))
sink()

system(paste("cp ",output.path,"/",sub('\\..*$', '', basename(output.path)),".log ",record.path,sep=""))


