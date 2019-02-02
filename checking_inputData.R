options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  selet.analysis <-args[1]
  input.name <- args[2]
  output.path <- args[3]
} else {
  stop("Not enough arguments")
}

dta <- read.delim(paste(input.name,sep=""),header=T)
names(dta)<-tolower(names(dta))
if(nrow(dta) == 0){
  cat("Error \n")
  cat("no input data")
  stop("no input data")
} else { 
  if(selet.analysis == "snp"){
	if(length(which(names(dta) %in% c("gene","snp"))) < 2 ){
           cat("Error \n")
           cat("check variable names, not enough column(s)")
	   stop("check variable names, not enough column(s)")  
	 }else{
         tmp=dta[,c("gene","snp")]
         if(nrow(tmp) < 1){
           cat("Error \n")
           cat("no input data")
           stop("no input data")
	 }
         write.table(tmp,paste(output.path,"/gene_snp_names.csv",sep=""),row.names=F,col.names=F,quote=F,sep=",")
	 write.table(tmp$snp,paste(output.path,"/snp_names.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
	 write.table(tmp$gene,paste(output.path,"/gene_names.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
        }
   } else if(selet.analysis == "gene" | selet.analysis == "omnibus"){
	if(length(which(names(dta) %in% c("gene"))) < 1){
           cat("Error \n")
           cat("check variable names, not enough column(s), no input data")
	   stop("check variable names, not enough column(s), no input data")  
	 }else{
         tmp=dta[,c("gene")]
         if(length(tmp) < 1){
           cat("Error \n")
           cat("no input data")
           stop("no input data")
         }
	 write.table(tmp,paste(output.path,"/gene_names.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
        }
	
   } else if(selet.analysis == "score"){
	if(length(which(names(dta) %in% c("trait","snp","effect_allele","other_allele","beta","se","effect_af"))) < 7){
           cat("Error \n")
           cat("check variable names, not enough column(s), no input data")
	   stop("check variable names, not enough column(s), no input data")  
	 }else{
         tmp=dta[,c("trait","snp","effect_allele","other_allele","beta","se","effect_af")]
         if(nrow(tmp) < 1){
           cat("Error")
           stop("Error: no input data")
         }
         write.table(tmp$snp,paste(output.path,"/snp_names.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
	 write.table(tmp,paste(output.path,"/score_data.txt",sep=""),row.names=F,quote=F,sep="\t")
        }
  }
  cat("Good")
}
   