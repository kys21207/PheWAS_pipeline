
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
      #print(i)
      vari <-read.table(paste(output.path,"/gene_assoc/",i,sep=""),header=T)
      vari <- vari[,c("GeneName","N","CASES","CONTS","DF","LRT","PLRT","DFavg","LRTavg","PLRTavg","BETA","SE","OR","PWALDavg","PHENOTYPE")] 
      names(vari)=c("GENENAME","N","CASES","CONTS","DF","LRT","P","DFavg","LRTavg","Pavg","BETA","SE","OR","PWavg","pheno")
      vari <- vari[vari$CASES > 200,]
      if(nrow(vari) > 0){
      vari$TYPE1 <- "Binary"
      vari$TYPE2 <- "GENE"
      if(is.null(BI.result)){
        BI.result <- vari
      }else{
        BI.result <- rbind(BI.result,vari)
      }
     }
    }
   }

   vari <-read.table(paste(output.path,"/gene_assoc/omnibus.linear",sep=""),header=T)
   vari <- vari[,c("GeneName","N","DF","SS","PVALUE","DFavg","SSavg","PVALUEavg","BETA","SE","PWALDavg","PHENOTYPE")] 
   names(vari)=c("GENENAME","N","DF","SS","P","DFavg","SSavg","Pavg","BETA","SE","PWavg","pheno")
   vari$TYPE1 <- "QTL"
   vari$TYPE2 <- "GENE"
   QT.result <- vari

  file1 <- list.files(paste(output.path,"/gene_assoc/",sep=""), pattern="*sig.*\\.logistic$")
  arg1 <- identical(file1, character(0))
  if(!arg1){
   BI.sig.result<-NULL
   for(i in file1){
     print(i)
     vari <-read.delim(paste(output.path,"/gene_assoc/",i,sep=""),header=T)
     vari <- vari[,c("GENENAME","TISSUE","r2","PHENOTYPE","CASES","CONTS","OR","BETA","SE","L95","U95","PVALUE","rsid")] 
     names(vari)=c("GENENAME","Tissue","r2","pheno","CASES","CONTS","OR","BETA","SE","L95","U95","P","rsIDs")
     vari <- vari[vari$CASES > 200,]
     if(nrow(vari) > 0){
     vari$TYPE1 <- "Binary"
     vari$TYPE2 <- "GENE"
     if(is.null(BI.sig.result)){
       BI.sig.result <- vari
     }else{
       BI.sig.result <- rbind(BI.sig.result,vari)
     }
    } 
  }
     BI.sig.result <- merge(BI.sig.result,anno.BI, by="pheno")
   }


  file1 <- list.files(paste(output.path,"/gene_assoc/",sep=""), pattern="*sig.linear$")
  arg2 <- identical(file1, character(0))
  if(!arg2){

   vari <-read.delim(paste(output.path,"/gene_assoc/univariable.sig.linear",sep=""),header=T)
   vari <- vari[,c("GENENAME","TISSUE","r2","PHENOTYPE","N","BETA","SE","L95","U95","PVALUE","rsid")] 
   names(vari)=c("GENENAME","Tissue","r2","pheno","N","BETA","SE","L95","U95","P","rsIDs")
   vari$TYPE1 <- "QTL"
   vari$TYPE2 <- "GENE"
   QT.sig.result <- vari
   QT.sig.result <- merge(QT.sig.result,anno.QT, by="pheno")

 }

   
  BI.result <- merge(BI.result,anno.BI, by="pheno")
  QT.result <- merge(QT.result,anno.QT, by="pheno")
  combin.result <- sbind(BI.result,QT.result)

  if( !arg1 & !arg2){
     combin.sig.result <- sbind(BI.sig.result,QT.sig.result)
   }else if( !arg1 & arg2){
     combin.sig.result <- BI.sig.result
   }else if(arg1 & !arg2){
     combin.sig.result <- QT.sig.result
   }

  for(i in unique(combin.result$GENENAME)){
     sub.result <- subset(combin.result, GENENAME %in% as.character(i))
     sub.result$FDR <- p.adjust(sub.result$P,method="BH")
     sub.result$FDRavg <- p.adjust(sub.result$Pavg,method="BH")
     write.table(sub.result,paste(output.path,"/gene_assoc/",as.character(i),".omnibus.summary.txt",sep=""),row.names=F,quote=F,sep="\t")
     sub.result <- sub.result[,!(names(sub.result) %in% c("TYPE1","TYPE2"))]
     write.table(sub.result,paste(output.path,"/gene_assoc/",as.character(i),".omnibus.phewas.txt",sep=""),row.names=F,quote=F,sep="\t")

 #    sub.sig.result <- subset(combin.sig.result, GENENAME %in% as.character(i))
 #    sub.sig.result$FDR <- p.adjust(sub.sig.result$P,method="BH")
 #    write.table(sub.sig.result,paste(output.path,"/gene_assoc/",as.character(i),".sig.univariable.summary.txt",sep=""),row.names=F,quote=F,sep="\t")
 #    sub.sig.result <- sub.sig.result[,!(names(sub.sig.result) %in% c("TYPE1","TYPE2"))]
 #    write.table(sub.sig.result,paste(output.path,"/gene_assoc/",as.character(i),".sig.univariable.phewas.txt",sep=""),row.names=F,quote=F,sep="\t")

  }

 




