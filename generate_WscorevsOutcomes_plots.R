options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 4) {
  phenos.path <- args[1]
  work.path <- args[2]
  covar.name <-args[3]
  output.path <- args[4]
} else {
  stop("Not enough arguments")
}

library(ggplot2)
library(data.table)
library(GenABEL)

#phenos.path="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/test1_201708011446"
#output.path="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/test1_201708011446/scor_assoc"
#work.path="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/test1_201708011446"
#covar.name="/GWD/appbase/projects/RD-TSci-PhewasUKB/PheWAS_output_storage/test1_201708011446/covars.txt"

covar <- fread(covar.name,header = TRUE)

covars <- names(covar)[-c(1:2)]

types<-c("QTL","Binary")
output.files <- list.files(paste(output.path,"/",sep=""), pattern=paste("*[.]phewas[.]txt$",sep=""))
output.arg <- identical(output.files, character(0))
if(!output.arg){
  for(i in output.files){
     dta <- read.delim(paste(output.path,"/",i,sep=""),header=T)
     sig.p.snp <- 0.05/nrow(dta)
#    dta$riskfactor<-dta$SCORENAME
     score.name<-unique(dta$riskfactor[!is.na(dta$riskfactor)])
	 list.pheno <- dta[dta$P <= sig.p.snp,]
	 if(nrow(list.pheno) > 0){
            scores=fread(paste(work.path,"/genotype.scor",sep=""),header=T)
	    scores=subset(scores, ,c("ID_1",as.character(score.name)))
#       d2=fread(paste(work.path,"/unweighted_genotype.scor",sep=""),header=T)
            binary.pheno <- subset(list.pheno, TYPE1 %in% "Binary")
	    if(nrow(binary.pheno) > 0){
	      pheno.dta <- fread(paste(phenos.path,"/phenotypes_BI.txt",sep=""),header = TRUE)
              pheno.dta <- subset(pheno.dta, ,c("FID",as.character(binary.pheno$pheno)))
              all.dta <- merge(scores,pheno.dta,by.x="ID_1",by.y="FID")
              all.dta <- merge(all.dta,covar,by.x="ID_1",by.y="FID")		  
          for(k in binary.pheno$pheno){
             sub.dta<-na.omit(subset(all.dta,,c("ID_1",as.character(score.name),as.character(k),as.character(covars))))
             sub.dta<-sub.dta[get(as.character(k)) != -9]
             breaks = c(seq(min(sub.dta[,score.name,with=F]),max(sub.dta[,score.name,with=F]),by =(max(sub.dta[,score.name,with=F])-min(sub.dta[,score.name,with=F]))/10))
             sub.dta[, bin :=cut(get(as.character(score.name)),breaks,labels = 1:10)]
             sub.dta[, interv :=cut(get(as.character(score.name)),breaks)]
             df0=as.data.frame(sub.dta[,.(.N),by=bin])
			 names(df0)=c("bin","N.x")
             df1=as.data.frame(sub.dta[,.(.N),by=interv])
			 df1=cbind(df0,df1)
             df1$prop <-100*df1$N/sum(df1$N)

             df2=as.data.frame(sub.dta[,list(coef=summary(glm(get(as.character(k)) ~ get(as.character(score.name))+age+gender+PC1+PC2+PC3+PC4+PC5,family = binomial))$coef[2],se=summary(glm(get(as.character(k)) ~ get(as.character(score.name))+age+gender+PC1+PC2+PC3+PC4+PC5,family = binomial))$coef[4]),by=interv])
             df2=na.omit(merge(df1,df2,by="interv"))
             df2$lower=exp(df2$coef-1.96*df2$se)
             df2$upper=exp(df2$coef+1.96*df2$se)
             df2$OR=exp(df2$coef)
             df2<-df2[order(df2$bin),]
             names(df2)<-c("Risk.Allele.Score","bin","X1","x2","Proportion(%)","coef","se","Lower","Upper","OR")
 #            tabular((Risk.Allele.Score +1)~(Proportion=1)+OR+Lower+Upper,data=df2)
             df2<-df2[,c("Risk.Allele.Score","Proportion(%)","OR","Lower","Upper")]            		 
             write.table(df2,paste(output.path,"/plots/",k,".scoreDist.txt",sep=""),row.names=F,quote=F,sep="\t")
             
			 
		  }
      }
       qt.pheno <- subset(list.pheno, TYPE1 %in% "QTL")
	   if(nrow(qt.pheno) > 0){
	      pheno.dta <- fread(paste(phenos.path,"/phenotypes_QT.txt",sep=""),header = TRUE)
              pheno.dta <- subset(pheno.dta, ,c("FID",as.character(qt.pheno$pheno)))
              all.dta <- merge(scores,pheno.dta,by.x="ID_1",by.y="FID")		  
              all.dta <- merge(all.dta,covar,by.x="ID_1",by.y="FID")		  
          for(k in qt.pheno$pheno){
             sub.dta<-na.omit(subset(all.dta,,c("ID_1",as.character(score.name),as.character(k),as.character(covars))))
             breaks = c(seq(min(sub.dta[,score.name,with=F]),max(sub.dta[,score.name,with=F]),by =(max(sub.dta[,score.name,with=F])-min(sub.dta[,score.name,with=F]))/10))
             sub.dta[, bin :=cut(get(as.character(score.name)),breaks,labels = 1:10)]
             sub.dta[, interv :=cut(get(as.character(score.name)),breaks)]
             sub.dta[, residual:=summary(lm(get(as.character(k)) ~ get(as.character(score.name))+age+gender+PC1+PC2+PC3+PC4+PC5))$residuals]
             sub.dta[, trait := rntransform(residual)]
             df0=as.data.frame(sub.dta[,.(.N),by=bin])
			 names(df0)=c("bin","N.x")
             df1=as.data.frame(sub.dta[,.(.N),by=interv])
			 df1=cbind(df0,df1)
             df1$prop <-100*df1$N/sum(df1$N)
                          

             df2=as.data.frame(sub.dta[,list(coef=summary(lm(trait ~ get(as.character(score.name))))$coefficients[2],se=summary(lm(trait ~ get(as.character(score.name))))$coefficients[4]),by=interv])
             df2=merge(df1,df2,by="interv")

             file.out <-paste(output.path,"/plots/",k,".scoreDist.pdf",sep="")
             pdf(file=file.out, height=8.5, width= 11)
 
             df1<-na.omit(df1[order(df1$bin),])
             df2<-na.omit(df2[order(df2$bin),])

             par(mar = c(5,5,2,5))
  
             mp<-with(df1, barplot(prop, col="red3", 
                 ylab="Proportion (%)",xlab="Risk allele scores",
                 ylim=c(0,max(df1$prop)),axes=F))
             text(mp, par('usr')[3], labels = df1$interv, srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
             axis(2)
#             na.omit(all)
             par(new = T)
             with(df2, plot(coef, axes=F, xlab=NA, ylab=NA))
             with(df2, points(bin,coef))
             segments(as.numeric(df2$bin),as.numeric(df2$coef-df2$se),as.numeric(df2$bin),as.numeric(df2$coef+df2$se))
             epsilon=0.02
             segments(as.numeric(df2$bin)-epsilon,as.numeric(df2$coef-df2$se),as.numeric(df2$bin)+epsilon,as.numeric(df2$coef-df2$se))
             segments(as.numeric(df2$bin)-epsilon,as.numeric(df2$coef+df2$se),as.numeric(df2$bin)+epsilon,as.numeric(df2$coef+df2$se))
             axis(side = 4)
             mtext(side = 4, line = 3, paste("Beta of ",as.character(k)))
             dev.off()
		  }
	 }
	  
   }
  }
 }
