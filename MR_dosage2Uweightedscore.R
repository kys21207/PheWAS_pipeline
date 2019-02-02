library(data.table)
options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  input.name <- args[1]
  genos.name <- args[2]
  output.path <- args[3]
} else {
  stop("Not enough arguments")
}
#input.name="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/MR/MR_test_SNP_lists_for_phewas.txt"
#genos.name="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/MR/genotype"
#output.path="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/MR"

input.dta <- read.delim(input.name,header=T)
info <- read.table(paste(genos.name, "map", sep = "."), header = TRUE, colClasses=c("numeric","character","numeric","character","character","numeric"))
dosage <- read.table(paste(genos.name, "dose", sep = "."), header = TRUE)
#dosage1<-dosage
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
for(j in unique(input.dta$trait)){
temp.dta<-subset(input.dta, trait %in% j)
temp.dta<-merge(temp.dta,info,by="snp")
for(i in 1:nrow(temp.dta)){
  switch(temp.dta$alleles[i],
     AC={
       switch(temp.dta$alleles.info[i],
	 CA = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] },
         GT = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] }
	   )},
     AG={
       switch(temp.dta$alleles.info[i],
	 GA = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] },
         CT = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] }
	   )},
     TC={
       switch(temp.dta$alleles.info[i],
	 CT = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] },
         GA = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] }
	   )},
     TG={
       switch(temp.dta$alleles.info[i],
	 GT = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] },
         CA = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] }
	   )},
     GA={
       switch(temp.dta$alleles.info[i],
	 AG = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] },
         TC = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] }
	   )},
     CT={
       switch(temp.dta$alleles.info[i],
	 TC = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] },
         AG = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] }
	   )},
     GT={
       switch(temp.dta$alleles.info[i],
	 TG = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] },
         AC = { dosage[,temp.dta$snp[i]]<-2-dosage[,temp.dta$snp[i]] }
      	   )},
     DI={
	   dosage[,temp.dta$snp[i]] <- ifelse(nchar(temp.dta$effall[i]) > 1 & nchar(temp.dta$otherall[i]) == 1, 2-dosage[,temp.dta$snp[i]], dosage[,temp.dta$snp[i]]) 
	   },
     ID={
	   dosage[,temp.dta$snp[i]] <- ifelse(nchar(temp.dta$effall[i]) == 1 & nchar(temp.dta$otherall[i]) > 1, 2-dosage[,temp.dta$snp[i]], dosage[,temp.dta$snp[i]]) 
	   },
     AT={
       switch(temp.dta$alleles.info[i],
	 TA = { 
	     dosage[,temp.dta$snp[i]]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), 2-dosage[,temp.dta$snp[i]],dosage[,temp.dta$snp[i]])
         },
         AT = { 
	     dosage[,temp.dta$snp[i]]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), 2-dosage[,temp.dta$snp[i]],dosage[,temp.dta$snp[i]])
         }
	   )},
     GC={
       switch(temp.dta$alleles.info[i],
	 CG = { 
	     dosage[,temp.dta$snp[i]]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), 2-dosage[,temp.dta$snp[i]],dosage[,temp.dta$snp[i]]) 
		 },
         GC = { 
	     dosage[,temp.dta$snp[i]]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), 2-dosage[,temp.dta$snp[i]],dosage[,temp.dta$snp[i]]) 
		 }
	   )},
     TA={
       switch(temp.dta$alleles.info[i],
	 AT = { 
	     dosage[,temp.dta$snp[i]]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), 2-dosage[,temp.dta$snp[i]],dosage[,temp.dta$snp[i]]) 
		 },
         TA = { 
	     dosage[,temp.dta$snp[i]]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), 2-dosage[,temp.dta$snp[i]],dosage[,temp.dta$snp[i]]) 
		 }
	   )},
     CG={
       switch(temp.dta$alleles.info[i],
	 GC = { 
	     dosage[,temp.dta$snp[i]]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), 2-dosage[,temp.dta$snp[i]],dosage[,temp.dta$snp[i]]) 
		 },
         CG = { 
	     dosage[,temp.dta$snp[i]]<-ifelse((!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] <= 0.45 & temp.dta$eaf.info[i] >= 0.55) | (!is.na(temp.dta$eaf[i]) & temp.dta$eaf[i] >= 0.55 & temp.dta$eaf.info[i] <= 0.45), 2-dosage[,temp.dta$snp[i]],dosage[,temp.dta$snp[i]])
		 }
	   )}
    )
}
#if(is.null(updated.input)){
#   updated.input<-temp.dta
#}else{
#   updated.input<-rbind(updated.input,temp.dta)
#}

dosage.sub<-subset(dosage,,c("ID_1",temp.dta$snp))
#dosage1.sub <-subset(dosage1,,c("ID_1",temp.dta$snp))
   
dosage.sub$Uscore <- rowSums(dosage.sub[,-1], na.rm = TRUE)
dosage.sub<-subset(dosage.sub,,c("ID_1","Uscore"))
names(dosage.sub)[names(dosage.sub) == "Uscore"] <- paste("Uscore_",j,sep="")

if(is.null(new.dosage)){
   new.dosage<-dosage.sub
}else{
   new.dosage<-merge(new.dosage,dosage.sub,by="ID_1")
}

}
#write.table(updated.input,"test.txt",row.names=F,quote=F,sep="\t")
	 
write.table(new.dosage,paste(output.path,"/unweighted_genotype.scor",sep=""),row.names=F,quote=F,sep="\t")
write.table(colnames(new.dosage)[-1],paste(output.path,"/scor_names.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")

