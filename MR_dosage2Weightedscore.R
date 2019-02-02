options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  input.name <- args[1]
  genos.name <- args[2]
  output.path <- args[3]
} else {
  stop("Not enough arguments")
}
#input.name="/GWD/appbase/projects/RD-TSci-UKB/PheWAS_output_storage/Cambridge_MR_batch2_201805100818/score_data.txt"
#genos.name="/GWD/appbase/projects/RD-TSci-UKB/PheWAS_output_storage/Cambridge_MR_batch2_201805100818/genotype"
#output.path="/GWD/appbase/projects/RD-TSci-UKB/PheWAS_output_storage/Cambridge_MR_batch2_201805100818"

input.dta <- read.delim(input.name,header=T)
info <- read.table(paste(genos.name, "map", sep = "."), header = TRUE, colClasses=c("numeric","character","numeric","character","character","numeric"))
dosage <- read.table(paste(genos.name, "dose", sep = "."), header = TRUE)
#info: chr snp pos otherall effall maf
#input.dta: trait,snp,chr,pos,effect_allele,other_allele,beta,se,pvalue,eaf,gene
#exclude input SNPs which are not available in UK Biobank
#colnames(input.dta)<-c(
input.dta$alleles<-toupper(paste(input.dta$effect_allele,input.dta$other_allele,sep=""))
input.dta$snp<-as.character(input.dta$snp)
info$alleles.info<-toupper(paste(info$effall,info$otherall,sep=""))
info<-subset(info,,c("snp","alleles.info","otherall","effall","eaf"))
names(info)[names(info) == "eaf"] <- "eaf.info"
new.dosage<-NULL

updated.input<-NULL
for(j in unique(input.dta$trait)){
print(j)
temp.dta<-subset(input.dta, trait %in% j)
temp.dta<-merge(temp.dta,info,by="snp")
if(nrow(temp.dta) > 0){
for(i in 1:nrow(temp.dta)){
  switch(temp.dta$alleles[i],
     AC={
       switch(temp.dta$alleles.info[i],
	 CA = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         GT = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         AC = { temp.dta$beta[i] <- temp.dta$beta[i] },
         TG = { temp.dta$beta[i] <- temp.dta$beta[i] }		 
	   )},
     AG={
       switch(temp.dta$alleles.info[i],
	 GA = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         CT = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         AG = { temp.dta$beta[i] <- temp.dta$beta[i] },
         TC = { temp.dta$beta[i] <- temp.dta$beta[i] }		 
	   )},
	 TC={
       switch(temp.dta$alleles.info[i],
	 CT = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         GA = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         TC = { temp.dta$beta[i] <- temp.dta$beta[i] },
         AG = { temp.dta$beta[i] <- temp.dta$beta[i] }		 
	   )},
     TG={
       switch(temp.dta$alleles.info[i],
	 GT = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         CA = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         TG = { temp.dta$beta[i] <- temp.dta$beta[i] },
         AC = { temp.dta$beta[i] <- temp.dta$beta[i] }		 
	   )},
     GA={
       switch(temp.dta$alleles.info[i],
	 AG = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         TC = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         GA = { temp.dta$beta[i] <- temp.dta$beta[i] },
         CT = { temp.dta$beta[i] <- temp.dta$beta[i] }		 
	   )},
     CT={
       switch(temp.dta$alleles.info[i],
	 TC = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         AG = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         CT = { temp.dta$beta[i] <- temp.dta$beta[i] },
         GA = { temp.dta$beta[i] <- temp.dta$beta[i] }		 
	   )},
     GT={
       switch(temp.dta$alleles.info[i],
	 TG = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         AC = { temp.dta$beta[i] <- temp.dta$beta[i]*(-1) },
         GT = { temp.dta$beta[i] <- temp.dta$beta[i] },
         CA = { temp.dta$beta[i] <- temp.dta$beta[i] }	 
	   )},
     DI={
	   temp.dta$beta[i] <- ifelse(nchar(temp.dta$effall[i]) > 1 & nchar(temp.dta$otherall[i]) == 1, temp.dta$beta[i]*(-1), temp.dta$beta[i]) 
	   },
     ID={
	   temp.dta$beta[i] <- ifelse(nchar(temp.dta$effall[i]) == 1 & nchar(temp.dta$otherall[i]) > 1, temp.dta$beta[i]*(-1), temp.dta$beta[i]) 
	   },
     AT={
       switch(temp.dta$alleles.info[i],
	 TA = { 
		     temp.dta$beta[i]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), temp.dta$beta[i]*(-1),temp.dta$beta[i]) 
		 },
         AT = { 
		     temp.dta$beta[i]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), temp.dta$beta[i]*(-1),temp.dta$beta[i]) 
		 }
	   )},
     GC={
       switch(temp.dta$alleles.info[i],
	 CG = { 
		     temp.dta$beta[i]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), temp.dta$beta[i]*(-1),temp.dta$beta[i]) 
		 },
         GC = { 
		     temp.dta$beta[i]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), temp.dta$beta[i]*(-1),temp.dta$beta[i]) 
		 }
	   )},
     TA={
       switch(temp.dta$alleles.info[i],
	 AT = { 
		     temp.dta$beta[i]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), temp.dta$beta[i]*(-1),temp.dta$beta[i]) 
		 },
         TA = { 
		     temp.dta$beta[i]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), temp.dta$beta[i]*(-1),temp.dta$beta[i]) 
		 }
	   )},
     CG={
       switch(temp.dta$alleles.info[i],
	 GC = { 
		     temp.dta$beta[i]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), temp.dta$beta[i]*(-1),temp.dta$beta[i]) 
		 },
         CG = { 
		     temp.dta$beta[i]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), temp.dta$beta[i]*(-1),temp.dta$beta[i]) 
		 }
	   )}
    )
}
if(is.null(updated.input)){
   updated.input<-temp.dta
}else{
   updated.input<-rbind(updated.input,temp.dta)
}

dosage.sub<-subset(dosage,,c("ID_1",temp.dta$snp))
if(i > 1){	   
dosage.rsName<-colnames(dosage.sub[,-1])
temp.dta<-temp.dta[match(temp.dta$snp,dosage.rsName),]
#dosage.sub$score <- as.matrix(dosage.sub[,-1]) %*% as.vector((temp.dta$beta/sum(temp.dta$beta))*nrow(temp.dta))
dosage.sub$score <- as.matrix(dosage.sub[,-1]) %*% as.vector(temp.dta$beta)
}else{
dosage.sub$score <- dosage.sub[,2]*temp.dta$beta
}
dosage.sub<-subset(dosage.sub,,c("ID_1","score"))
names(dosage.sub)[names(dosage.sub) == "score"] <- paste("score_",j,sep="")

if(is.null(new.dosage)){
   new.dosage<-dosage.sub
}else{
   new.dosage<-merge(new.dosage,dosage.sub,by="ID_1")
}
}
}
write.table(updated.input,paste(output.path,"/updated.data.txt",sep=""),row.names=F,quote=F,sep="\t")
	 
write.table(new.dosage,paste(output.path,"/genotype.scor",sep=""),row.names=F,quote=F,sep="\t")
write.table(colnames(new.dosage)[-1],paste(output.path,"/scor_names.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")

