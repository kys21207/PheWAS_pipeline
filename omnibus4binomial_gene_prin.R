library(data.table)
library(GenABEL)

options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 4) {
  phenos.name <- args[1]
  covars.name <- args[2]
  work.path <- args[3]
  output.path <- args[4]
  numFile <- args[5]
} else {
  stop("Not enough arguments")
}
out.signif <- 6

#work.path <- "/GWD/appbase/projects/RD-TSci-UKB/PheWAS_output_storage/KJ_genes_adult_pathway_201805081431"
#phenos.name<-"/GWD/appbase/projects/RD-TSci-UKB/PheWAS/phenotype500kV3/Sub_BI_pheno_1.txt"
#covars.name<-"/GWD/appbase/projects/RD-TSci-UKB/PheWAS/phenotype500kV3/covars.txt"
#output.path<-"/GWD/appbase/projects/RD-TSci-UKB/PheWAS_output_storage/KJ_genes_adult_pathway_201805081431/gene_assoc"

pheno <-read.table(paste(work.path,"/pheno.txt",sep=""),header=T)
# read all phenotypes
phenos <- fread(phenos.name, header = TRUE)
# read covariates 
covars <- fread(covars.name, header= TRUE)
# read predictDB.csv for annotation

merge.exprs <- fread(paste(work.path,"/gene_exp_prin2.txt",sep=""),header=T)
merge.exprs$aver.exp <- rowSums(as.matrix(merge.exprs[,-c(1,2)]),na.rm = T)
merge.exprs <- subset(merge.exprs,,c("FID","IID","aver.exp"))
#info.all <- read.delim(paste(work.path,"/all.info.txt",sep=""),header=T)

info.genes <- read.delim(paste(work.path,"/gene_names.txt",sep=""),header=T)

#write.table(merge.exprs,paste(output.path,"/exprs.txt",sep=""),row.names=F,quote=F,sep="\t")
stopifnot(length(unique(phenos$FID)) == nrow(phenos))

phenos.include <-  phenos$FID %in% merge.exprs$FID
covars.include <-  covars$FID %in% merge.exprs$FID

cat("Analysing ", sum(phenos.include), "/", nrow(phenos),
    " subjects with non-missing phenotypes, covariates, and genotypes\n", sep = "")
phenos <- phenos[phenos.include, , drop = FALSE]
covars <- covars[covars.include, , drop = FALSE]

phenos <- subset(phenos,,!(names(phenos) %in% names(pheno)[3]))
pheno.exclude <- pheno[[3]] != -9 & !is.na(pheno[[3]])
pheno <- pheno[pheno.exclude,,drop = FALSE]
merge.exprs <- data.matrix(subset(merge.exprs[match(pheno$FID, merge.exprs$FID),],, 3))
phenos.names <-names(phenos)
phenos <-data.matrix(subset(phenos[match(pheno$FID, phenos$FID),],,-2:-1))
covars <-data.matrix(subset(covars[match(pheno$FID, covars$FID),],,-2:-1))

cat("Analysis model: glm( ", names(pheno)[3], " ~ ",
    paste(names(pheno)[-3:-1], collapse = " + "), " )\n", sep = "")

pheno1 <- as.numeric(data.matrix(subset(pheno, , 3)))

#covars <- data.matrix(subset(covars, , -2:-1))
#if (ncol(covars) > 0) {
#  m0 <- glm(pheno1 ~ covars+merge.exprs, family = binomial)
#} else {
#  m0 <- glm(pheno1 ~ 1, family = binomial)
#}

#cat("Analysing", nrow(info), "genes\n")
if (ncol(covars) > 0) { # keep if statement outside vapply
  assoc <- as.data.frame(t(vapply(1:ncol(phenos), function(idx) {
 #  assoc <- as.data.frame(t(vapply(1:10, function(idx) {
    return(tryCatch({
      phen1<-phenos[,idx]
      m1 <- suppressWarnings(glm(pheno1 ~ phen1 + merge.exprs + phen1:merge.exprs + covars, family = binomial,na.action=na.omit))
      beta1 <- coef(m1)["phen1"]
      se1 <- sqrt(vcov(m1)["phen1", "phen1"])
      beta2 <- coef(m1)["merge.exprs"]
      se2 <- sqrt(vcov(m1)["merge.exprs", "merge.exprs"])
      beta3 <- coef(m1)["phen1:merge.exprs"]
      se3 <- sqrt(vcov(m1)["phen1:merge.exprs", "phen1:merge.exprs"])
      n <- length(na.omit(residuals(m1))) # seemingly no direct way to extract
      c(n, beta1, se1, beta2, se2, beta3, se3)
    },
                    error = function(e) return(c(NA, NA, NA, NA, NA, NA, NA))))},
                                  c(N = 0, BETA_pheno = 0., SE_pheno = 0., BETA_exp = 0., SE_exp = 0.,BETA_int = 0., SE_int = 0.))))
}
# calculate P-values for LRT
assoc$CASES <- sum(pheno1==1)
assoc$CONTS <- sum(pheno1==0)
assoc$BETA_pheno <- signif(assoc$BETA_pheno, out.signif)
assoc$OR_pheno <- signif(exp(assoc$BETA_pheno),3)
assoc$SE_pheno <- signif(assoc$SE_pheno, out.signif)
assoc$P_pheno <- with(assoc, signif(pnorm(-abs(BETA_pheno)/SE_pheno)*2, out.signif))

assoc$BETA_exp <- signif(assoc$BETA_exp, out.signif)
assoc$OR_exp <- signif(exp(assoc$BETA_exp),3)
assoc$SE_exp <- signif(assoc$SE_exp, out.signif)
assoc$P_exp <- with(assoc, signif(pnorm(-abs(BETA_exp)/SE_exp)*2, out.signif))

assoc$BETA_int <- signif(assoc$BETA_int, out.signif)
assoc$OR_int <- signif(exp(assoc$BETA_int),3)
assoc$SE_int <- signif(assoc$SE_int, out.signif)
assoc$P_int <- with(assoc, signif(pnorm(-abs(BETA_int)/SE_int)*2, out.signif))

assoc$PhenoName <- as.character(phenos.names[-2:-1])
assoc$PHENOTYPE<-names(pheno)[3]

sel <- !is.na(assoc$P_pheno)
cat("Reporting valid association results for ",
    sum(sel), "/", nrow(assoc), " genes\n", sep = "")

write.table(assoc[sel, , drop = FALSE],
            file = paste(output.path,"/omnibus.",numFile,".logistic",sep=""),
            row.names = FALSE, quote = FALSE)
