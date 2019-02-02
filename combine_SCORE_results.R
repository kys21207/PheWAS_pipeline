
options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  output.path <- args[1]
  score.name <- args[2]
  anno.path <-args[3]
} else {
  stop("Not enough arguments")
}


#output.path="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/LY_MRtest_201708091235"
#score.name="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/LY_MRtest_201708091235/scor_names.txt"
#anno.path="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS/phenotype"

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

 score.names=read.table(score.name, header=FALSE)

 names(score.names)<-c("riskfactor")

# pheno.data <- read.table(pheno.name,header=T,sep="\t")

 names.for.gene <- as.vector(unique(na.omit(score.names$riskfactor)))
 
  # combine all single variant results across the binary traits
  binary.gene.files <- list.files(paste(output.path,"/scor_assoc/",sep=""), pattern="*[.]logistic$")
  binary.gene.arg <- identical(binary.gene.files, character(0))
  if(!binary.gene.arg){
     n.binary.gene <- length(binary.gene.files)

 for(i in 1:n.binary.gene){
#print(i)
# chr gene pos otherall effall analysed.Freq1 analysed.Rsq N BETA SE PVALUE L95 U95 PHENOTYPE

   vari <-read.table(paste(output.path,"/scor_assoc/",binary.gene.files[i],sep=""),header=T,colClasses=c("character",rep("numeric",11),"character"))
   vari <- vari[,c("trait_score","N","CASES","CONTS","OR","BETA","SE","L95","U95","PWALD","PHENOTYPE")] 
   names(vari)=c("riskfactor","N","CASES","CONTS","OR","BETA","SE","L95","U95","P","pheno")
   vari <- vari[vari$CASES > 200,]
   vari$TYPE1 <- "Binary"
   vari$TYPE2 <- "SCORE"
   all.binary.gene <- vari
   for(k in names.for.gene){
     if(i==1){
      write.table(all.binary.gene[all.binary.gene$riskfactor == k & !is.na(all.binary.gene$riskfactor),],paste(output.path,"/scor_assoc/",gsub("/","_",k),".BIsummary.txt",sep=""),row.names=F,quote=F,sep="\t")
     } else {
      write.table(all.binary.gene[all.binary.gene$riskfactor == k & !is.na(all.binary.gene$riskfactor),],paste(output.path,"/scor_assoc/",gsub("/","_",k),".BIsummary.txt",sep=""),row.names=F,col.names=F,quote=F,append=TRUE,sep="\t")
     }
   }

 }
}

#print(dim(all.binary.snp))


# combine all results for single variants across the QTL traits

qtl.gene.files <- list.files(paste(output.path,"/scor_assoc/",sep=""), pattern="*[.]linear$")
qtl.gene.arg <- identical(qtl.gene.files, character(0))
if(!qtl.gene.arg){

  n.qtl.gene <- length(qtl.gene.files)

# combine all single variant results 
 for(i in 1:n.qtl.gene){
#print(i)
 
#chr snp pos otherall effall analysed.Freq1 analysed.Rsq N BETA SE PVALUE L95 U95 PHENOTYPE

   vari <-read.table(paste(output.path,"/scor_assoc/",qtl.gene.files[i],sep=""),header=T,colClasses=c("character",rep("numeric",6),"character"))
   vari <- vari[,c("trait_score","N","BETA","SE","L95","U95","PVALUE","PHENOTYPE")] 
   names(vari)=c("riskfactor","N","BETA","SE","L95","U95","P","pheno")
   vari$TYPE1 <- "QTL"
   vari$TYPE2 <- "SCORE"
   all.qtl.gene <- vari
#   all.qtl.snp <- merge(all.qtl.gene,gene.snp.names,by="riskfactor",all.x=T,all.y=F) 
   for(k in names.for.gene){
     if(i==1){
      write.table(all.qtl.gene[all.qtl.gene$riskfactor == k & !is.na(all.qtl.gene$riskfactor),],paste(output.path,"/scor_assoc/",gsub("/","_",k),".QTsummary.txt",sep=""),row.names=F,quote=F,sep="\t")
     } else {
      write.table(all.qtl.gene[all.qtl.gene$riskfactor == k & !is.na(all.qtl.gene$riskfactor),],paste(output.path,"/scor_assoc/",gsub("/","_",k),".QTsummary.txt",sep=""),row.names=F,col.names=F,quote=F,append=TRUE,sep="\t")
     }
   }
}
}


# combine all snp results
gene.BIfiles <- list.files(paste(output.path,"/scor_assoc/",sep=""), pattern="*[.]BIsummary[.]txt$")
gene.QTfiles <- list.files(paste(output.path,"/scor_assoc/",sep=""), pattern="*[.]QTsummary[.]txt$")

#check whether a file is empty or not
n.cell=1
exclude.cell=NULL
for(k in gene.BIfiles){
  if( file.info(paste(output.path,"/scor_assoc/",k,sep=""))$size <= 200){
    exclude.cell=c(exclude.cell,n.cell)
  }
   n.cell=n.cell+1 
}
  if(!is.null(exclude.cell)){
    gene.BIfiles=gene.BIfiles[-exclude.cell]
  }
   gene.BInames <- as.character(unlist(lapply( strsplit(gene.BIfiles,"[.]"),function(x)x[[1]])))

n.cell=1
exclude.cell=NULL
for(k in gene.QTfiles){
  if( file.info(paste(output.path,"/scor_assoc/",k,sep=""))$size <= 200){
    exclude.cell=c(exclude.cell,n.cell)
  }
   n.cell=n.cell+1 
}
  if(!is.null(exclude.cell)){
   gene.QTfiles=gene.QTfiles[-exclude.cell]
  }
   gene.QTnames <- as.character(unlist(lapply( strsplit(gene.QTfiles,"[.]"),function(x)x[[1]])))

#find the available gene names in both binary and QT traits
intersect.names <- intersect(gene.BInames,gene.QTnames)
diff.gene.BInames <- setdiff(gene.BInames,gene.QTnames)
diff.gene.QTnames <- setdiff(gene.QTnames,gene.BInames)

anno.QT <- read.csv(paste(anno.path,"/Anno_HESelf_QT.csv",sep=""))
anno.BI <- read.csv(paste(anno.path,"/Anno_HESelf_BI.csv",sep=""))

# please check the column names
anno.QT <-anno.QT[, c("pheno","PHENOTYPE","Category")]
anno.BI <-anno.BI[, c("pheno","PHENOTYPE","Category")]

anno.QT <- subset(anno.QT, !is.na(pheno))
anno.BI <- subset(anno.BI, !is.na(pheno))

for(i in intersect.names){
 # print(i)
  BI.file <-read.table(paste(output.path,"/scor_assoc/",i,".BIsummary.txt",sep=""),header=T,sep="\t")
    
  BI.result <- merge(BI.file,anno.BI, by="pheno")
  QT.file <-read.table(paste(output.path,"/scor_assoc/",i,".QTsummary.txt",sep=""),header=T,sep="\t")
    
  QT.result <- merge(QT.file,anno.QT, by="pheno")
  combin.result <- sbind(BI.result,QT.result)
  combin.result$FDR <- p.adjust(combin.result$P,method="BH")
  write.table(combin.result,paste(output.path,"/scor_assoc/",gsub("/","_",i),".summary.txt",sep=""),row.names=F,quote=F,sep="\t")
  combin.result <- combin.result[,!(names(combin.result) %in% c("TYPE1","TYPE2"))]
  write.table(combin.result,paste(output.path,"/scor_assoc/",gsub("/","_",i),".phewas.txt",sep=""),row.names=F,quote=F,sep="\t")
}

BI.gene.arg <- identical(diff.gene.BInames, character(0))
if(!BI.gene.arg){
for(i in diff.gene.BInames){
 # print(i)
  BI.file <-read.table(paste(output.path,"/scor_assoc/",i,".BIsummary.txt",sep=""),header=T,sep="\t")
  BI.result <- merge(BI.file,anno.BI, by="pheno")

  BI.result$FDR <- p.adjust(BI.result$P,method="BH")
  write.table(BI.result,paste(output.path,"/scor_assoc/",gsub("/","_",i),".summary.txt",sep=""),row.names=F,quote=F,sep="\t")
  BI.result <- BI.result[,!(names(BI.result) %in% c("TYPE1","TYPE2"))]
  write.table(BI.result,paste(output.path,"/scor_assoc/",gsub("/","_",i),".phewas.txt",sep=""),row.names=F,quote=F,sep="\t")
}
}

QT.gene.arg <- identical(diff.gene.QTnames, character(0))
if(!QT.gene.arg){
for(i in diff.gene.QTnames){
 # print(i)
  QT.file <-read.table(paste(output.path,"/scor_assoc/",i,".QTsummary.txt",sep=""),header=T,sep="\t")
  QT.result <- merge(QT.file,anno.QT, by="pheno")

  QT.result$FDR <- p.adjust(QT.result$P,method="BH")
  write.table(QT.result,paste(output.path,"/scor_assoc/",gsub("/","_",i),".summary.txt",sep=""),row.names=F,quote=F,sep="\t")
  QT.result <- QT.result[,!(names(QT.result) %in% c("TYPE1","TYPE2"))]
  write.table(QT.result,paste(output.path,"/scor_assoc/",gsub("/","_",i),".phewas.txt",sep=""),row.names=F,quote=F,sep="\t")
}
}

 




