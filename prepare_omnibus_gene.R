library(data.table)

options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  predictDB.path <- args[1]
  expres.path <- args[2]
  work.path <- args[3]
} else {
  stop("Not enough arguments")
}

#predictDB.path="/GWD/appbase/projects/RD-TSci-UKB/PredictDB"
#expres.path="/GWD/appbase/projects/RD-TSci-UKB/PheWAS_output_storage/KJ_genes_asthma_pathway_201805081431/predicted_exp"
#work.path="/GWD/appbase/projects/RD-TSci-UKB/PheWAS_output_storage/KJ_genes_asthma_pathway_201805081431"


stopwords=c("TW_","_0.5")
tissues <- as.character(as.matrix(read.table(paste(work.path,"/tissue_names.txt",sep=""),header=F)))
genes <- as.character(as.matrix(read.table(paste(work.path,"/gene_names.txt",sep=""),header=F)))
tissues.name=gsub(paste(stopwords,collapse="|"),"",tissues)

merge.exprs <- NULL
info.all <- NULL
j=1
for(i in tissues){
   # read predictDB.csv for annotation
   info <- read.csv(paste(predictDB.path,"/",i,".csv", sep = ""), header=TRUE)
   # convert genenames to gene
   tmp <- subset(info, genename %in% genes,select=c(gene,genename,rsid,r2))
   if(nrow(tmp) > 0){
   aggr.rsNUM <- aggregate(rsid ~ genename , tmp,function(x){y=table(x); toString(unlist(names(y)[y==max(y)]))})
   tmp<- tmp[!duplicated(tmp$gene),]
   tmp$rsid<-NULL
   tmp<-merge(tmp,aggr.rsNUM,by="genename")
   tmp$tissue<-as.character(tissues.name[j])
   columns.genes <- unique(as.vector(as.matrix(subset(info, genename %in% genes,select=c(gene)))))
   exprs <- fread(paste("gunzip -c ",expres.path,"/",i,".txt.gz", sep = ""), header = TRUE)
   # check whether there are genes available in predictDB or not
   columns.match <- na.omit(match(c("FID","IID",columns.genes),names(exprs)))
   exprs <- subset(exprs,,columns.match)
   if(length(columns.match) > 2){
      # subset by pre-defined genes
      names(exprs)<-c("FID","IID",paste0(names(exprs)[-c(1,2)],".",tissues.name[j]))
      if(is.null(merge.exprs)){
	     merge.exprs <- exprs
	  }else{
	     merge.exprs <- merge(merge.exprs,exprs,by=c("FID","IID"))
	  }
       if(is.null(info.all)){
           info.all <- tmp
       }else{
           info.all <- rbind(info.all,tmp)
       }
   }
  j<-j+1
  }
}

write.table(merge.exprs,paste(work.path,"/merge.exprs.txt",sep=""),row.names=F,quote=F,sep="\t")
write.table(info.all,paste(work.path,"/all.info.txt",sep=""),row.names=F,quote=F,sep="\t")

