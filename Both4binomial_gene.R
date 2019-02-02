library(data.table)
library(GenABEL)

options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 5) {
  phenos.name <- args[1]
  covars.name <- args[2]
  work.path <- args[3]
  output.path <- args[4]
  numFile <- args[5]
} else {
  stop("Not enough arguments")
}
out.signif <- 6

#work.path <- "/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/test500_201709132026"
#phenos.name<-"/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/test500_201709132026/Sub_BI_pheno_1.txt"
#covars.name<-"/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/test500_201709132026/covars.txt"
#output.path<-"/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/test500_201709132026/gene_assoc"

# read all phenotypes
phenos <- fread(phenos.name, header = TRUE)
# read covariates 
covars <- fread(covars.name, header= TRUE)
# read predictDB.csv for annotation

merge.exprs <- fread(paste(work.path,"/merge.exprs.txt",sep=""),header=T)
info.all <- read.delim(paste(work.path,"/all.info.txt",sep=""),header=T)

info.genes <- info.all[!duplicated(info.all$gene),c("gene","genename")]

#write.table(merge.exprs,paste(output.path,"/exprs.txt",sep=""),row.names=F,quote=F,sep="\t")
stopifnot(length(unique(phenos$FID)) == nrow(phenos))

phenos.include <-  phenos$FID %in% merge.exprs$FID
covars.include <-  covars$FID %in% merge.exprs$FID

cat("Analysing ", sum(phenos.include), "/", nrow(phenos),
    " subjects with non-missing phenotypes, covariates, and genotypes\n", sep = "")
phenos <- phenos[phenos.include, , drop = FALSE]
covars <- covars[covars.include, , drop = FALSE]
end.covar <- ncol(covars)+1 

phenos <- merge(covars, phenos, by=c("FID","IID"))

stopifnot(nrow(phenos) != 0)

cat("Using ", length(phenos$FID), "/", nrow(merge.exprs),
    " genotyped subjects\n", sep = "")

check<-0
check2<-0
end.covar=16
for(i in end.covar:ncol(phenos)){

pheno <- subset(phenos,,c(1,2,i,3:(end.covar-1)))
# exclude missing subjects
pheno <- na.omit(pheno)
pheno.exclude <- pheno[[3]] != -9 & !is.na(pheno[[3]])
pheno <- pheno[pheno.exclude,,drop = FALSE]
merge.expr <- merge.exprs[merge.exprs$FID %in% pheno$FID,,drop=FALSE]
pheno <- pheno[pheno$FID %in% merge.expr$FID,,drop=FALSE]
merge.expr <- subset(merge.expr[match(pheno$FID, merge.expr$FID),],, 3:ncol(merge.expr))

cat("Analysis model: glm( ", names(pheno)[3], " ~ ",
    paste(names(pheno)[-3:-1], collapse = " + "), " )\n", sep = "")

pheno1 <- as.numeric(data.matrix(subset(pheno, , 3)))

covars <- data.matrix(subset(pheno, , -3:-1))
if (ncol(covars) > 0) {
  m0 <- glm(pheno1 ~ covars, family = binomial)
} else {
  m0 <- glm(pheno1 ~ 1, family = binomial)
}

#cat("Analysing", nrow(info), "genes\n")
if (ncol(covars) > 0) { # keep if statement outside vapply
  assoc <- as.data.frame(t(vapply(1:nrow(info.genes), function(idx) {
    return(tryCatch({
      geno1<-data.matrix(subset(merge.expr, , names(merge.expr)[grepl(info.genes$gene[idx],names(merge.expr))]))
      # Clean data by means of winsorization, i.e., by shrinking outlying observations to the border of the main part of the data.
      if(is.null(ncol(geno1))){
          geno2<-geno1
      } else {
          geno2<-rowMeans(geno1)
      }
      m1 <- suppressWarnings(glm(pheno1 ~ geno1 + covars, family = binomial))
      lrt1 <- max(m0$deviance - m1$deviance, 0)
      df1 <- ifelse(is.null(ncol(geno1)),1,ncol(geno1))
      m2 <- suppressWarnings(glm(pheno1 ~ geno2 + covars, family = binomial))
      lrt2 <- max(m0$deviance - m2$deviance, 0)
      n <- length(na.omit(residuals(m1))) # seemingly no direct way to extract
      c(n, df1, lrt1, lrt2)
    },
                    error = function(e) return(c(NA, NA, NA, NA))))},
                                  c(N = 0, DF = 0, LRT = 0., LRTavg = 0.))))
}
# calculate P-values for LRT
assoc$CASES <- sum(pheno1==1)
assoc$CONTS <- sum(pheno1==0)
assoc$PLRT <- with(assoc, signif(pchisq(LRT, DF, lower.tail = FALSE), out.signif))
assoc$LRT <- signif(assoc$LRT, out.signif)
assoc$PLRTavg <- with(assoc, signif(pchisq(LRTavg, df = 1, lower.tail = FALSE), out.signif))
assoc$DFavg <- 1
assoc$LRTavg <- signif(assoc$LRTavg, out.signif)
assoc$GeneName <- as.character(info.genes$genename)
assoc$PHENOTYPE<-names(pheno)[3]

sel <- !is.na(assoc$PLRT)
cat("Reporting valid association results for ",
    sum(sel), "/", nrow(assoc), " genes\n", sep = "")

if(check == 0){
   write.table(assoc[sel, , drop = FALSE],
            file = paste(output.path,"/omnibus.",numFile,".logistic",sep=""),
            row.names = FALSE, quote = FALSE)
   check<-1
}else{
   write.table(assoc[sel, , drop = FALSE],
            file = paste(output.path,"/omnibus.",numFile,".logistic",sep=""),
            row.names = FALSE,col.names = FALSE, quote = FALSE, append=TRUE)
}

sig.output<-NULL
sig.pvalue<-0.05/(ncol(phenos)-end.covar)
sig.assoc<-assoc
if(nrow(sig.assoc) > 0 & sig.assoc$CASES[1] > 20){
   for(j in 1:nrow(sig.assoc)){
      assoc <- as.data.frame(t(vapply(names(merge.expr)[grepl(as.character(info.genes[info.genes$genename == sig.assoc$GeneName[j],"gene"]),names(merge.expr))], function(idx) {
      return(tryCatch({
      geno1 <- data.matrix(subset(merge.expr, , idx))
      # Clean data by means of winsorization, i.e., by shrinking outlying observations to the border of the main part of the data.
#      geno1 <- winsorize(geno1)
      m1 <- suppressWarnings(glm(pheno1 ~ geno1 + covars, family = binomial))
      beta <- coef(m1)["geno1"]
      se <- sqrt(vcov(m1)["geno1", "geno1"])
      n <- length(na.omit(residuals(m1))) # seemingly no direct way to extract
      c(n, beta, se)
    },
                    error = function(e) return(c(NA, NA, NA))))},
                                  c(N = 0, BETA = 0., SE = 0.))))
   
assoc$CASES <- sum(pheno1==1)
assoc$CONTS <- sum(pheno1==0)
assoc$BETA <- signif(assoc$BETA, out.signif)
assoc$OR <- signif(exp(assoc$BETA),3)
assoc$SE <- signif(assoc$SE, out.signif)
assoc$L95 <- signif(exp(assoc$BETA-1.96*assoc$SE),3)
assoc$U95 <- signif(exp(assoc$BETA+1.96*assoc$SE),3)
assoc$PVALUE <- with(assoc, signif(pnorm(-abs(BETA)/SE)*2, out.signif))
assoc$TISSUE<-unlist(lapply( strsplit(names(merge.expr)[grepl(as.character(info.genes[info.genes$genename == sig.assoc$GeneName[j],"gene"]),names(merge.expr))],"[.]"),function(x)x[[3]]))
assoc$GENENAME<-sig.assoc$GeneName[j]
assoc$PHENOTYPE<-names(pheno)[3]

sig.output<-merge(assoc,info.all,by.x=c("GENENAME","TISSUE"),by.y=c("genename","tissue"))
if(check2 == 0){
   write.table(sig.output,
            file = paste(output.path,"/univariable.sig.",numFile,".logistic",sep=""),
            row.names = FALSE, quote = FALSE,sep="\t")
   check2<-1
}else{
   write.table(sig.output,
            file = paste(output.path,"/univariable.sig.",numFile,".logistic",sep=""),
            row.names = FALSE,col.names = FALSE, quote = FALSE, append=TRUE,sep="\t")
}
  } #for
} #if


}
