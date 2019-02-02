options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  select <- args[1]
  path <- args[2]
} else {
  stop("Not enough arguments")
}

if(select == "snp"){
  dta <- read.table(paste(path,"/snp_names.txt",sep=""),header=F)
}else{
  dta <- read.table(paste(path,"/scor_names.txt",sep=""),header=F)
}
ListNames=as.vector(dta$V1)
if(length(ListNames) < 50) {
  for(i in 1:length(ListNames)){
     sub_names <- ListNames[i]
     if(select == "snp"){
       write.table(sub_names,paste(path,"/snps_assoc/snp_group_",i,".txt",sep=""),row.names=F,col.names=F,quote=F,sep=" ")
     }else{
       write.table(sub_names,paste(path,"/scor_assoc/scor_group_",i,".txt",sep=""),row.names=F,col.names=F,quote=F,sep=" ")
     }       
  }
}else{
  num <- trunc(length(ListNames)/50)
  left <- length(ListNames)-(50*num)
  k=1
  end=0
  for(i in 1:50){
     start <-1+end
	 if(k <= left){
	   end <- num*i + 1*i
	   k <- k+1
	 }else{
	   end <- num + (start-1)
	 }
     sub_names <- ListNames[start:end]
     if(select == "snp"){
       write.table(sub_names,paste(path,"/snps_assoc/snp_group_",i,".txt",sep=""),row.names=F,col.names=F,quote=F,sep=" ")
     }else{
       write.table(sub_names,paste(path,"/scor_assoc/scor_group_",i,".txt",sep=""),row.names=F,col.names=F,quote=F,sep=" ")
     }       
  }
}  
