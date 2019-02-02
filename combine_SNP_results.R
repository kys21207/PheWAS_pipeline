
options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  output.path <- args[1]
  gene.snp.name <- args[2]
  anno.path <-args[3]
} else {
  stop("Not enough arguments")
}


#output.path="/GWD/appbase/projects/RD-TSci-UKB/PheWAS_output_storage/Bisect_asthma_less_201804191316"
#gene.snp.name="/GWD/appbase/projects/RD-TSci-UKB/PheWAS_output_storage/Bisect_asthma_less_201804191316/gene_snp_names.csv"
#anno.path="/GWD/appbase/projects/RD-TSci-UKB/PheWAS/KJ/phenotype1"

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

 gene.snp.names=read.csv(gene.snp.name, header=FALSE)
 gene.snp.names=data.frame(lapply(gene.snp.names, trimws), stringsAsFactors = FALSE)

 if(all(grepl("rs",na.omit(gene.snp.names$V1)))){
   names(gene.snp.names)<-c("SNP","GENENAME")
 }else{
   names(gene.snp.names)<-c("GENENAME","SNP")
 }
# pheno.data <- read.table(pheno.name,header=T,sep="\t")

 names.for.gene <- as.vector(unique(na.omit(gene.snp.names$GENENAME)))
 gene.snp.names <- na.omit(gene.snp.names)
 names.for.snp <-  as.vector(unique(na.omit(gene.snp.names$GENENAME)))
 
  # combine all single variant results across the binary traits
  binary.snp.files <- list.files(paste(output.path,"/snps_assoc/",sep=""), pattern="*[.]logistic$")
  binary.snp.arg <- identical(binary.snp.files, character(0))
  if(!binary.snp.arg){
     n.binary.snp <- length(binary.snp.files)

 for(i in 1:n.binary.snp){
print(i)
# chr snp pos otherall effall analysed.Freq1 analysed.Rsq N BETA SE PVALUE L95 U95 PHENOTYPE

   vari <-read.table(paste(output.path,"/snps_assoc/",binary.snp.files[i],sep=""),header=T)
#   vari <-vari[,c(-5,-11)]
   vari <- vari[,c("chr","snp","pos","otherall","effall","analysed.Freq1","analysed.Rsq","N","CASES","CONTS","OR","BETA","SE","L95","U95","PWALD","PHENOTYPE")] 
   names(vari)=c("CHR","SNP","BP","otherall","effall","analysed.Freq1","Rsq","N","CASES","CONTS","OR","BETA","SE","L95","U95","P","pheno")
   vari <- vari[vari$CASES > 200,]
   if(nrow(vari) > 0){
   vari$TYPE1 <- "Binary"
   vari$TYPE2 <- "SNP"
   all.binary.snp <- vari
   all.binary.snp <- merge(all.binary.snp,gene.snp.names,by="SNP",all.x=T,all.y=F) 
   for(k in names.for.snp){
     if(i==1){
      write.table(all.binary.snp[all.binary.snp$GENENAME == k & !is.na(all.binary.snp$GENENAME),],paste(output.path,"/snps_assoc/",k,".BIsummary.txt",sep=""),row.names=F,quote=F,sep="\t")
     } else {
      write.table(all.binary.snp[all.binary.snp$GENENAME == k & !is.na(all.binary.snp$GENENAME),],paste(output.path,"/snps_assoc/",k,".BIsummary.txt",sep=""),row.names=F,col.names=F,quote=F,append=TRUE,sep="\t")
     }
   }
   # SNP analysis only without gene names
   if(i==1){
     write.table(all.binary.snp[which(is.na(all.binary.snp$GENENAME)),],paste(output.path,"/snps_assoc/nogene.BIsummary.txt",sep=""),row.names=F,quote=F,sep="\t")
   } else {
     write.table(all.binary.snp[which(is.na(all.binary.snp$GENENAME)),],paste(output.path,"/snps_assoc/nogene.BIsummary.txt",sep=""),row.names=F,col.names=F,quote=F,append=TRUE,sep="\t")
   }
  }
 }
}

#print(dim(all.binary.snp))


# combine all results for single variants across the QTL traits

qtl.snp.files <- list.files(paste(output.path,"/snps_assoc/",sep=""), pattern="*[.]linear$")
qtl.snp.arg <- identical(qtl.snp.files, character(0))
if(!qtl.snp.arg){

  n.qtl.snp <- length(qtl.snp.files)

# combine all single variant results 
 for(i in 1:n.qtl.snp){
#print(i)
 
#chr snp pos otherall effall analysed.Freq1 analysed.Rsq N BETA SE PVALUE L95 U95 PHENOTYPE

   vari <-read.table(paste(output.path,"/snps_assoc/",qtl.snp.files[i],sep=""),header=T)
#   pheno.name <- as.character(unlist(lapply( strsplit(qtl.snp.files[i],"_"),function(x)x[[1]])))
   vari <- vari[,c("chr","snp","pos","otherall","effall","analysed.Freq1","analysed.Rsq","N","BETA","SE","L95","U95","PVALUE","PHENOTYPE")] 
   names(vari)=c("CHR","SNP","BP","otherall","effall","analysed.Freq1","Rsq","N","BETA","SE","L95","U95","P","pheno")
   vari$TYPE1 <- "QTL"
   vari$TYPE2 <- "SNP"
   all.qtl.snp <- vari
   all.qtl.snp <- merge(all.qtl.snp,gene.snp.names,by="SNP",all.x=T,all.y=F) 
   for(k in names.for.snp){
     if(i==1){
      write.table(all.qtl.snp[all.qtl.snp$GENENAME == k & !is.na(all.qtl.snp$GENENAME),],paste(output.path,"/snps_assoc/",k,".QTsummary.txt",sep=""),row.names=F,quote=F,sep="\t")
     } else {
      write.table(all.qtl.snp[all.qtl.snp$GENENAME == k & !is.na(all.qtl.snp$GENENAME),],paste(output.path,"/snps_assoc/",k,".QTsummary.txt",sep=""),row.names=F,col.names=F,quote=F,append=TRUE,sep="\t")
     }
   }
   # SNP analysis only without gene names
   if(i==1){
     write.table(all.qtl.snp[which(is.na(all.qtl.snp$GENENAME)),],paste(output.path,"/snps_assoc/nogene.QTsummary.txt",sep=""),row.names=F,quote=F,sep="\t")
   } else {
     write.table(all.qtl.snp[which(is.na(all.qtl.snp$GENENAME)),],paste(output.path,"/snps_assoc/nogene.QTsummary.txt",sep=""),row.names=F,col.names=F,quote=F,append=TRUE,sep="\t")
   }
}
}


# combine all snp results
snps.BIfiles <- list.files(paste(output.path,"/snps_assoc/",sep=""), pattern="*[.]BIsummary[.]txt$")
snps.QTfiles <- list.files(paste(output.path,"/snps_assoc/",sep=""), pattern="*[.]QTsummary[.]txt$")

#check whether a file is empty or not
n.cell=1
exclude.cell=NULL
for(k in snps.BIfiles){
  if( file.info(paste(output.path,"/snps_assoc/",k,sep=""))$size <= 200){
    exclude.cell=c(exclude.cell,n.cell)
  }
   n.cell=n.cell+1 
}
  if(!is.null(exclude.cell)){
    snps.BIfiles=snps.BIfiles[-exclude.cell]
  }
   snps.BInames <- as.character(unlist(lapply( strsplit(snps.BIfiles,"[.]"),function(x)x[[1]])))

n.cell=1
exclude.cell=NULL
for(k in snps.QTfiles){
  if( file.info(paste(output.path,"/snps_assoc/",k,sep=""))$size <= 200){
    exclude.cell=c(exclude.cell,n.cell)
  }
   n.cell=n.cell+1 
}
  if(!is.null(exclude.cell)){
    snps.QTfiles=snps.QTfiles[-exclude.cell]
  }
   snps.QTnames <- as.character(unlist(lapply( strsplit(snps.QTfiles,"[.]"),function(x)x[[1]])))

#find the available gene names in both binary and QT traits
intersect.names <- intersect(snps.BInames,snps.QTnames)
diff.snps.BInames <- setdiff(snps.BInames,snps.QTnames)
diff.snps.QTnames <- setdiff(snps.QTnames,snps.BInames)

anno.QT <- read.csv(paste(anno.path,"/Anno_HESelf_QT.csv",sep=""))
anno.BI <- read.csv(paste(anno.path,"/Anno_HESelf_BI.csv",sep=""))

# please check the column names
anno.QT <-anno.QT[, c("pheno","PHENOTYPE","Category")]
anno.BI <-anno.BI[, c("pheno","PHENOTYPE","Category")]

anno.QT <- subset(anno.QT, !is.na(pheno))
anno.BI <- subset(anno.BI, !is.na(pheno))

for(i in intersect.names){
 # print(i)
  QT.file <-read.table(paste(output.path,"/snps_assoc/",i,".QTsummary.txt",sep=""),header=T,sep="\t")
  QT.file$otherall<-ifelse(QT.file$otherall=="TRUE","T",as.character(QT.file$otherall))
  QT.file$effall<-ifelse(QT.file$effall=="TRUE","T",as.character(QT.file$effall))
  BI.file <-read.table(paste(output.path,"/snps_assoc/",i,".BIsummary.txt",sep=""),header=T,sep="\t")
  BI.file$otherall<-ifelse(BI.file$otherall=="TRUE","T",as.character(BI.file$otherall))
  BI.file$effall<-ifelse(BI.file$effall=="TRUE","T",as.character(BI.file$effall))
    
  QT.result <- merge(QT.file,anno.QT, by="pheno")
  BI.result <- merge(BI.file,anno.BI, by="pheno")


  combin.result <- sbind(BI.result,QT.result)
#  combin.result$FDR <- p.adjust(combin.result$P,method="BH")
  write.table(combin.result,paste(output.path,"/snps_assoc/",i,".summary.txt",sep=""),row.names=F,quote=F,sep="\t")
#create results by individual SNP
 for(j in unique(combin.result$SNP)){
   snp.result<-combin.result[combin.result$SNP==j,]
   snp.result$FDR <- p.adjust(snp.result$P,method="BH")
   snp.result <- snp.result[,!(names(snp.result) %in% c("TYPE1","TYPE2"))]
   write.table(snp.result,paste(output.path,"/snps_assoc/",j,".",i,".phewas.txt",sep=""),row.names=F,quote=F,sep="\t")
 }
}

BI.snps.arg <- identical(diff.snps.BInames, character(0))
if(!BI.snps.arg){
for(i in diff.snps.BInames){
 # print(i)
  BI.file <-read.table(paste(output.path,"/snps_assoc/",i,".BIsummary.txt",sep=""),header=T,sep="\t")
  BI.file$otherall<-ifelse(BI.file$otherall=="TRUE","T",as.character(BI.file$otherall))
  BI.file$effall<-ifelse(BI.file$effall=="TRUE","T",as.character(BI.file$effall))
  BI.result <- merge(BI.file,anno.BI, by="pheno")

#  BI.result$FDR <- p.adjust(BI.result$P,method="BH")
  write.table(BI.result,paste(output.path,"/snps_assoc/",i,".summary.txt",sep=""),row.names=F,quote=F,sep="\t")
 for(j in unique(BI.result$SNP)){
   snp.result<-BI.result[BI.result$SNP==j,]
   snp.result$FDR <- p.adjust(snp.result$P,method="BH")
   snp.result <- snp.result[,!(names(snp.result) %in% c("TYPE1","TYPE2"))]
   write.table(snp.result,paste(output.path,"/snps_assoc/",j,".",i,".phewas.txt",sep=""),row.names=F,quote=F,sep="\t")
 }
}
}

QT.snps.arg <- identical(diff.snps.QTnames, character(0))
if(!QT.snps.arg){
for(i in diff.snps.QTnames){
 # print(i)
  QT.file <-read.table(paste(output.path,"/snps_assoc/",i,".QTsummary.txt",sep=""),header=T,sep="\t")
  QT.file$otherall<-ifelse(QT.file$otherall=="TRUE","T",as.character(QT.file$otherall))
  QT.file$effall<-ifelse(QT.file$effall=="TRUE","T",as.character(QT.file$effall))
  QT.result <- merge(QT.file,anno.QT, by="pheno")

#  QT.result$FDR <- p.adjust(QT.result$P,method="BH")
  write.table(QT.result,paste(output.path,"/snps_assoc/",i,".summary.txt",sep=""),row.names=F,quote=F,sep="\t")
 for(j in unique(QT.result$SNP)){
   snp.result<-QT.result[QT.result$SNP==j,]
   snp.result$FDR <- p.adjust(snp.result$P,method="BH")
   snp.result <- snp.result[,!(names(snp.result) %in% c("TYPE1","TYPE2"))]
   write.table(snp.result,paste(output.path,"/snps_assoc/",j,".",i,".phewas.txt",sep=""),row.names=F,quote=F,sep="\t")
 }
}
}

 

 




