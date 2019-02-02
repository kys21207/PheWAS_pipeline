
options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  output.path <- args[1]
  anno.path <-args[2]
} else {
  stop("Not enough arguments")
}


#binding two data frames that don't have the same set of columns
sbind = function(x, y, fill=NA) {
    sbind.fill = function(d, cols){ 
        for(c in cols)
            d[[c]] = fill
        d
    }

    x = sbind.fill(x, setdiff(names(y),names(x)))
    y = sbind.fill(y, setdiff(names(x),names(y)))

    rbind(x, y)
}


anno.QT <- read.csv(paste(anno.path,"/Anno_HESelf_QT.csv",sep=""))
anno.BI <- read.csv(paste(anno.path,"/Anno_HESelf_BI.csv",sep=""))

# please check the column names
anno.QT <-anno.QT[, c("pheno","PHENOTYPE","Category")]
anno.BI <-anno.BI[, c("pheno","PHENOTYPE","Category")]

anno.QT <- subset(anno.QT, !is.na(pheno))
anno.BI <- subset(anno.BI, !is.na(pheno))

  file1 <- list.files(paste(output.path,"/gene_assoc/",sep=""), pattern="omnibus.*\\.logistic$")
  arg1 <- identical(file1, character(0))
   if(!arg1){
    BI.result<-NULL
    for(i in file1){
      vari <-read.table(paste(output.path,"/gene_assoc/",i,sep=""),header=T)
      vari <- vari[,c("GeneName","N","CASES","CONTS","DF","LRT","PLRT","DFavg","LRTavg","PLRTavg","PHENOTYPE")] 
      names(vari)=c("GENENAME","N","CASES","CONTS","DF","LRT","P","DFavg","LRTavg","Pavg","pheno")
      vari <- vari[vari$CASES > 200,]
      vari$TYPE1 <- "Binary"
      vari$TYPE2 <- "GENE"
      if(is.null(BI.result)){
        BI.result <- vari
      }else{
        BI.result <- rbind(BI.result,vari)
      }
    }
   }

   vari <-read.table(paste(output.path,"/gene_assoc/omnibus.linear",sep=""),header=T)
   vari <- vari[,c("GeneName","N","DF","SS","PVALUE","DFavg","SSavg","PVALUEavg","PHENOTYPE")] 
   names(vari)=c("GENENAME","N","DF","SS","P","DFavg","SSavg","Pavg","pheno")
   vari$TYPE1 <- "QTL"
   vari$TYPE2 <- "GENE"
   QT.result <- vari

  file1 <- list.files(paste(output.path,"/gene_assoc/",sep=""), pattern="*ind.*\\.logistic$")
  arg1 <- identical(file1, character(0))
  if(!arg1){
   BI.ind.result<-NULL
   for(i in file1){
     vari <-read.delim(paste(output.path,"/gene_assoc/",i,sep=""),header=T)
     vari <- vari[,c("GENENAME","TISSUE","r2","PHENOTYPE","CASES","CONTS","OR","BETA","SE","L95","U95","PVALUE","rsid")] 
     names(vari)=c("GENENAME","Tissue","r2","pheno","CASES","CONTS","OR","BETA","SE","L95","U95","P","rsIDs")
     vari <- vari[vari$CASES > 20,]
     vari$TYPE1 <- "Binary"
     vari$TYPE2 <- "GENE"
     if(is.null(BI.ind.result)){
       BI.ind.result <- vari
     }else{
       BI.ind.result <- rbind(BI.ind.result,vari)
     }
   }
     BI.ind.result <- merge(BI.ind.result,anno.BI, by="pheno")
  }


  file1 <- list.files(paste(output.path,"/gene_assoc/",sep=""), pattern="*ind.linear$")
  arg2 <- identical(file1, character(0))
  if(!arg2){

   vari <-read.delim(paste(output.path,"/gene_assoc/univariable.ind.linear",sep=""),header=T)
   vari <- vari[,c("GENENAME","TISSUE","r2","PHENOTYPE","N","BETA","SE","L95","U95","PVALUE","rsid")] 
   names(vari)=c("GENENAME","Tissue","r2","pheno","N","BETA","SE","L95","U95","P","rsIDs")
   vari$TYPE1 <- "QTL"
   vari$TYPE2 <- "GENE"
   QT.ind.result <- vari
   QT.ind.result <- merge(QT.ind.result,anno.QT, by="pheno")

 }

   
  BI.result <- merge(BI.result,anno.BI, by="pheno")
  QT.result <- merge(QT.result,anno.QT, by="pheno")
  combin.result <- sbind(BI.result,QT.result)

  if( !arg1 & !arg2){
     combin.ind.result <- sbind(BI.ind.result,QT.ind.result)
   }else if( !arg1 & arg2){
     combin.ind.result <- BI.ind.result
   }else if(arg1 & !arg2){
     combin.ind.result <- QT.ind.result
   }

  for(i in unique(combin.result$GENENAME)){
     sub.result <- subset(combin.result, GENENAME %in% as.character(i))
     sub.result$FDR <- p.adjust(sub.result$P,method="BH")
     sub.result$FDRavg <- p.adjust(sub.result$Pavg,method="BH")
     write.table(sub.result,paste(output.path,"/gene_assoc/",as.character(i),".omnibus.summary.txt",sep=""),row.names=F,quote=F,sep="\t")
     sub.result <- sub.result[,!(names(sub.result) %in% c("TYPE1","TYPE2"))]
     write.table(sub.result,paste(output.path,"/gene_assoc/",as.character(i),".omnibus.phewas.txt",sep=""),row.names=F,quote=F,sep="\t")

     sub.ind.result <- subset(combin.ind.result, GENENAME %in% as.character(i))
     sub.ind.result$FDR <- p.adjust(sub.ind.result$P,method="BH")
     write.table(sub.ind.result,paste(output.path,"/gene_assoc/",as.character(i),".ind.univariable.summary.txt",sep=""),row.names=F,quote=F,sep="\t")
     sub.ind.result <- sub.ind.result[,!(names(sub.ind.result) %in% c("TYPE1","TYPE2"))]
     write.table(sub.ind.result,paste(output.path,"/gene_assoc/",as.character(i),".ind.univariable.phewas.txt",sep=""),row.names=F,quote=F,sep="\t")

  }

 




