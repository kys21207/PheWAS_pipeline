library(data.table)
library(GenABEL)

options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 4) {
  phenos.name <- args[1]
  covars.name <- args[2]
  work.path <- args[3]
  output.path <- args[4]
} else {
  stop("Not enough arguments")
}
out.signif <- 6


#phenos.name="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/test4/phenotypes_QT.txt"
#covars.name="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/test4/covars.txt"
#work.path="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/test4"
#output.path="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/test4/gene_assoc"

# read all phenotypes
phenos <- fread(phenos.name, header = TRUE)
# read covariates 
covars <- fread(covars.name, header= TRUE)
# read predictDB.csv for annotation

merge.exprs <- fread(paste(work.path,"/merge.exprs.txt",sep=""),header=T)
info.all <- read.delim(paste(work.path,"/all.info.txt",sep=""),header=T)

info.genes <- info.all[!duplicated(info.all$gene),c("gene","genename")]

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
for(i in end.covar:ncol(phenos)){

pheno <- subset(phenos,,c(1,2,i,3:(end.covar-1)))
# exclude missing subjects
pheno <- na.omit(pheno)
pheno.exclude <- pheno[[3]] != -9 & !is.na(pheno[[3]])
pheno <- pheno[pheno.exclude,,drop = FALSE]
merge.expr <- merge.exprs[merge.exprs$FID %in% pheno$FID,,drop=FALSE]
pheno <- pheno[pheno$FID %in% merge.expr$FID,,drop=FALSE]
merge.expr <- subset(merge.expr[match(pheno$FID, merge.expr$FID),],, 3:ncol(merge.expr))

cat("Analysis model: lm( ", names(pheno)[3], " ~ ",
    paste(names(pheno)[-3:-1], collapse = " + "), " )\n", sep = "")

pheno1 <- as.numeric(data.matrix(subset(pheno, , 3)))
#Inverse normal transformation
pheno1 <- rntransform(pheno1)
covars <- data.matrix(subset(pheno, , -3:-1))
if (ncol(covars) > 0) {
  m0 <- lm(pheno1 ~ covars)
} else {
  m0 <- lm(pheno1 ~ 1)
}
cat("Analysing", nrow(info.genes), "genes\n")
if (ncol(covars) > 0) { # keep if statement outside vapply
  assoc <- as.data.frame(t(vapply(1:nrow(info.genes), function(idx) {
    return(tryCatch({
      geno1<-data.matrix(subset(merge.expr, , names(merge.expr)[grepl(info.genes$gene[idx],names(merge.expr))]))
      if(is.null(ncol(geno1))){
          geno2<-geno1
      } else {
          geno2<-rowMeans(geno1)
      }
      m1 <- suppressWarnings(lm(pheno1 ~ geno1+covars))
      m2 <- suppressWarnings(lm(pheno1 ~ geno2+covars))
      m01 <- anova(m0,m1)
      m02 <- anova(m0,m2)
      df01 <- m01[2,"Df"]
      ss01 <- m01[2,"Sum of Sq"]
      p01 <- m01[2,"Pr(>F)"]
      df02 <- m02[2,"Df"]
      ss02 <- m02[2,"Sum of Sq"]
      p02 <- m02[2,"Pr(>F)"]
      c(length(na.omit(m1$residuals)), df01, ss01, p01, df02, ss02, p02) # hacky way to get sample size
    },
                    error = function(e) return(c(NA, NA, NA, NA, NA, NA, NA))))},
                                  c(N = 0, DF = 0,  SS = 0., PVALUE = 0., DFavg = 0,  SSavg = 0., PVALUEavg = 0.))))
}
assoc$GeneName <- as.character(info.genes$genename)
assoc$PHENOTYPE<-names(pheno)[3]

sel <- !is.na(assoc$PVALUE)
cat("Reporting valid association results for ",
    sum(sel), "/", nrow(assoc), " genes\n", sep = "")

if(check == 0){
   write.table(assoc[sel, , drop = FALSE],
            file = paste(output.path,"/omnibus.linear",sep=""),
            row.names = FALSE, quote = FALSE)
   check<-1
}else{
   write.table(assoc[sel, , drop = FALSE],
            file = paste(output.path,"/omnibus.linear",sep=""),
            row.names = FALSE,col.names = FALSE, quote = FALSE, append=TRUE)
}

sig.output<-NULL
sig.pvalue<-0.05/(ncol(phenos)-end.covar)
sig.assoc<-assoc
if(nrow(sig.assoc) > 0){
   for(j in 1:nrow(sig.assoc)){
      assoc <- as.data.frame(t(vapply(names(merge.expr)[grepl(as.character(info.genes[info.genes$genename == sig.assoc$GeneName[j],"gene"]),names(merge.expr))], function(idx) {
      return(tryCatch({
      geno1<-data.matrix(subset(merge.expr, , idx))
      m1 <- suppressWarnings(lm(pheno1 ~ geno1 + covars))
      beta <- coef(m1)["geno1"]
      se <- sqrt(vcov(m1)["geno1", "geno1"])
      c(length(na.omit(m1$residuals)), beta, se) # hacky way to get sample size
    },
                    error = function(e) return(c(NA, NA, NA))))},
                                  c(N = 0, BETA = 0., SE = 0.))))
   assoc$PVALUE <- with(assoc, signif(pnorm(-abs(BETA)/SE)*2, out.signif))
   assoc$BETA <- signif(assoc$BETA, out.signif)
   assoc$SE <- signif(assoc$SE, out.signif)
   assoc$L95 <- signif((assoc$BETA-1.96*assoc$SE),out.signif)
   assoc$U95 <- signif((assoc$BETA+1.96*assoc$SE),out.signif)
   assoc$TISSUE<-unlist(lapply( strsplit(names(merge.expr)[grepl(as.character(info.genes[info.genes$genename == sig.assoc$GeneName[j],"gene"]),names(merge.expr))],"[.]"),function(x)x[[3]]))
   assoc$GENENAME<-sig.assoc$GeneName[j]
   assoc$PHENOTYPE<-names(pheno)[3]

   sig.output<-merge(assoc,info.all,by.x=c("GENENAME","TISSUE"),by.y=c("genename","tissue"))

   if(check2 == 0){
      write.table(sig.output,
            file = paste(output.path,"/univariable.sig.linear",sep=""),
            row.names = FALSE, quote = FALSE,sep="\t")
      check2<-1
   }else{
      write.table(sig.output,
            file = paste(output.path,"/univariable.sig.linear",sep=""),
            row.names = FALSE,col.names = FALSE, quote = FALSE, append=TRUE,sep="\t")
   }
  }
 }

}
