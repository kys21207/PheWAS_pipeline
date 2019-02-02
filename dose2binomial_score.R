library(data.table)
options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 6) {
  phenos.name <- args[1]
  covar.name <- args[2]
  genos.name <- args[3]
  output.path <- args[4]
  list.scors <- args[5]
  group.num <- args[6]
} else {
  stop("Not enough arguments")
}
out.signif <- 6

#phenos.name="phenotypes_BI.txt"
#genos.name="genotype"
#covar.name="covars.txt"
#output.path="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/test2/snps_assoc"

phenos <- fread(phenos.name,header = TRUE)
covar <- fread(covar.name,header = TRUE)
dosage <- fread(paste(genos.name, "scor", sep = "."), header = TRUE)

#need the following scripts for GRID 
listSCORs <- read.table(list.scors,header=F)
info <- listSCORs
listSCORs <- as.vector(listSCORs$V1)
names(info)<-"trait_score"
dosage <- subset(dosage, ,c("ID_1",listSCORs))

phenos <- merge(covar ,phenos,by=c("FID","IID"))
end.covar.col <- ncol(covar)+1

stopifnot(length(unique(phenos$FID)) == nrow(phenos))

#phenos.include <- apply(!is.na(phenos), 1, all) & phenos$FID %in% dosage$ID_1
phenos.include <- phenos$FID %in% dosage$ID_1
cat("Analysing ", sum(phenos.include), "/", nrow(phenos),
    " subjects with non-missing phenotypes, covariates, and genotypes\n", sep = "")
phenos <- phenos[phenos.include, , drop = FALSE]

cat("Using ", length(phenos.include), "/", nrow(dosage),
    " genotyped subjects\n", sep = "")
#stopifnot(nrow(dose) == nrow(phenos))
#stopifnot(ncol(dosage) == nrow(info))
stopifnot(!any(is.na(dosage)))


check <- 0
#ptm <- proc.time()

for(i in end.covar.col:ncol(phenos)){
pheno <- subset(phenos,,c(1,2,i,3:(end.covar.col-1)))
# exclude missing subjects
pheno <- na.omit(pheno)
pheno.exclude <- pheno[[3]] != -9 & !is.na(pheno[[3]])
pheno <- pheno[pheno.exclude,,drop = FALSE]
dose <- dosage[dosage$ID_1 %in% pheno$FID,,drop=FALSE]
dose <- dose[match(pheno$FID, dose$ID_1), 2:ncol(dose)]
stopifnot(nrow(dose) == nrow(pheno))
stopifnot(ncol(dose) == nrow(info))

cat("Analysis model: glm( ", names(pheno)[3], " ~ ",
    paste(names(pheno)[-3:-1], collapse = " + "), " )\n", sep = "")
	
pheno1 <- data.matrix(subset(pheno, , 3))
if(sum(pheno1==1) >= 20){
covars <- data.matrix(subset(pheno, , -3:-1))
if (ncol(covars) > 0) {
  m0 <- glm(pheno1 ~ covars, family = binomial)
} else {
  m0 <- glm(pheno1 ~ 1, family = binomial)
}

cat("Analysing", nrow(info), "variants\n")
if (ncol(covars) > 0) { # keep if statement outside vapply
  assoc <- as.data.frame(t(vapply(1:nrow(info), function(idx) {
    return(tryCatch({
      geno1 <- data.matrix(subset(dose, , idx))
      m1 <- suppressWarnings(glm(pheno1 ~ geno1 + covars, family = binomial))
      lrt <- max(m0$deviance - m1$deviance, 0)
      beta <- coef(m1)["geno1"]
      se <- sqrt(vcov(m1)["geno1", "geno1"])
      n <- length(na.omit(residuals(m1))) # seemingly no direct way to extract
      c(n, beta, se, lrt)
    },
                    error = function(e) return(c(NA, NA, NA, NA))))},
                                  c(N = 0, BETA = 0., SE = 0., LRT = 0.))))
} else {
  assoc <- as.data.frame(t(vapply(1:nrow(info), function(idx) {
    return(tryCatch({
      geno1 <- data.matrix(subset(dose, , idx))
      m1 <- suppressWarnings(glm(pheno1 ~ geno1, family = binomial))
      lrt <- max(m0$deviance - m1$deviance, 0)
      beta <- coef(m1)["geno1"]
      se <- sqrt(vcov(m1)["geno1", "geno1"])
      n <- length(na.omit(residuals(m1))) # seemingly no direct way to extract
      c(n, beta, se, lrt)
    },
                    error = function(e) return(c(NA, NA, NA, NA))))},
                                  c(N = 0, BETA = 0., SE = 0., LRT = 0.))))
}

# calculate P-values before rounding BETA, SE, LRT
assoc$CASES <- sum(pheno1 == 1)
assoc$CONTS <- sum(pheno1 ==0)
assoc$PWALD <- with(assoc, signif(pnorm(-abs(BETA)/SE)*2, out.signif))
assoc$PLRT <- with(assoc, signif(pchisq(LRT, df = 1, lower.tail = FALSE), out.signif))
assoc$BETA <- signif(assoc$BETA, out.signif)
assoc$OR <- signif(exp(assoc$BETA),3)
assoc$SE <- signif(assoc$SE, out.signif)
assoc$LRT <- signif(assoc$LRT, out.signif)
assoc$L95 <- signif(exp(assoc$BETA-1.96*assoc$SE),3)
assoc$U95 <- signif(exp(assoc$BETA+1.96*assoc$SE),3)
assoc$PHENOTYPE<-names(pheno)[3]

sel <- !is.na(assoc$LRT)
cat("Reporting valid association results for ",
    sum(sel), "/", nrow(assoc), " variants\n", sep = "")
if(check == 0){
   write.table(cbind(info[sel,, drop = FALSE], assoc[sel, , drop = FALSE]),
            file = paste(output.path,"/scor_",group.num,".logistic",sep=""),
            row.names = FALSE, quote = FALSE)
   check <- 1
}else{
   write.table(cbind(info[sel,, drop = FALSE], assoc[sel, , drop = FALSE]),
            file = paste(output.path,"/scor_",group.num,".logistic",sep=""),col.names=FALSE,
            row.names = FALSE, quote = FALSE,append=TRUE)
}
#print(proc.time()-ptm)

}
}
