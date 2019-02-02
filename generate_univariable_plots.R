library(ggplot2)

options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  output.path <- args[1]
} else {
  stop("Not enough arguments")
}

prepare_dta.forestPlot_QT <- function(dta,gene.name,p_cut,num_cut){

x=subset(dta, P != 0 & P < p_cut)
if(nrow(x) > 0){
x=x[!duplicated(x),]
x=x[order(x$P),]
if(nrow(x) < num_cut) {
  x=x[1:nrow(x),]
}else {
  x=x[1:num_cut,]
}
cluster=as.data.frame(aggregate(P ~ Tissue, x,function(x){min(x)}))
cluster$rank=rank(-cluster$P)
cluster$P=NULL
x=merge(x,cluster,by="Tissue")
x=x[order(rank(x$Tissue),-x$P),]
x=x[order(x$rank),]

label <- factor(sprintf("%s [ %5.0f]",paste(x$PHENOTYPE,(max(x$rank)-x$rank+1),sep="."),x$N))
mean <- x$BETA
lower<- x$L95
upper<- x$U95
pvalue <- -log10(x$P)
pvalue <-ifelse(pvalue > 20, 20,pvalue)
title <-paste("Plot of ",x$GENENAME[1],sep="")

df <- data.frame(label, mean, lower, upper, pvalue, x$Tissue, title)
df$style="plain"
#df$just=1
#insert categories
df$ind <- seq_len(nrow(df))
row.max.category=as.data.frame(aggregate(ind ~ x.Tissue, df,function(x){max(x)}))
for(i in 1:nrow(row.max.category)){
    df <- rbind(df,data.frame(label = row.max.category$x.Tissue[i], mean = NA,lower = NA, upper=NA, pvalue=NA, x.Tissue=NA, title=NA, style="bold",ind=row.max.category$ind[i]+0.5))
}
df <- df[order(df$ind),]
# reverses the factor level ordering for labels after coord_flip()
df$label <- factor(df$label, levels=df$label)
df
}else{
df<-NULL
}

}

forestPlot_QT<-function(df){
 title <- df$title[1]
plot<-ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
        theme_bw()+
        ggtitle(title) + 
        geom_pointrange() + 
		geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) +
		geom_point(aes(size = pvalue, colour=pvalue)) + 
        scale_colour_gradientn(name="-log10(pvalue)",colors = 
		c("darkblue","lightblue","green","yellow","red"),
                         breaks=c(20,10,8,6,4,2), labels = c("+20-10","10-8","8-6","6-4","4-2","<2"), limits=c(0,20)) +
		scale_size(name="-log10(pvalue)",breaks=c(20,10,8,6,4,2), labels = c("+20-10","10-8","8-6","6-4","4-2","<2"), limits=c(0,20)) +
        guides(color=guide_legend(), size=guide_legend()) +
        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
 		theme( axis.text.y = element_text(face = df$style, color = 'gray30', size = 10, hjust=0, angle=0)
		#axis.text.y = element_text(face = rev(Style), color = 'gray30', size = 14, hjust=0, angle=0),
#                    axis.title.x = element_text(size = 20, color = 'gray30', vjust = 0),
                   ) +
        xlab("") + ylab("Beta (95% CI)")  
plot
}


prepare_dta.forestPlot_BI <- function(dta,gene.name,p_cut,num_cut){
#p_cut: pvalue threshold to show
#num_cut: # of traits to show
x=subset(dta, P < p_cut & P != 0)
if(nrow(x) >0){
x=x[!duplicated(x),]
x=x[order(x$P),]
if(nrow(x) < num_cut) {
  x=x[1:nrow(x),]
}else {
  x=x[1:num_cut,]
}
cluster=as.data.frame(aggregate(P ~ Tissue, x,function(x){min(x)}))
cluster$rank=rank(-cluster$P)
cluster$P=NULL
x=merge(x,cluster,by="Tissue")
x=x[order(rank(x$Tissue),-x$P),]
x=x[order(x$rank),]

label <- factor(sprintf("%s [%5.0f, %5.0f]",paste(x$PHENOTYPE,(max(x$rank)-x$rank+1),sep="."),x$CASES,x$CONT))
OR <- x$OR
lower<- x$L95
upper<- x$U95
pvalue <- -log10(x$P)
pvalue <-ifelse(pvalue > 20, 20,pvalue)
#title <-paste("Plot of ",x$snp[1]," in ",x$gene[1],sep="")
title <-paste("Plot of ",x$GENENAME[1],sep="")

df <- data.frame(label, OR, lower, upper, pvalue, x$Tissue, title)
df$style="plain"
#insert categories
df$ind <- seq_len(nrow(df))
row.max.category=as.data.frame(aggregate(ind ~ x.Tissue, df,function(x){max(x)}))
for(i in 1:nrow(row.max.category)){
    df <- rbind(df,data.frame(label = row.max.category$x.Tissue[i], OR = NA,lower = NA, upper=NA, pvalue=NA, x.Tissue=NA, title=NA, style="bold",ind=row.max.category$ind[i]+0.5))
}
df <- df[order(df$ind),]
df$label <- factor(df$label, levels=df$label)
df
}else{
df<-NULL
}
}


forestPlot_BI <- function(df){
 title <- df$title[1]
 ticks<-c(seq(.1, 1, by =.1), seq(0, 10, by =1), seq(10, 100, by =10))
plot<-ggplot(data=df, aes(x=label, y=OR, ymin=lower, ymax=upper)) +
        ggtitle(title) + 
        theme_bw() + # use a white background
        geom_pointrange() + 
		geom_errorbar(aes(ymin=lower, ymax=upper), width=.1) +
		geom_point(aes(size = pvalue, colour=pvalue)) + 
        scale_colour_gradientn(name="-log10(pvalue)",colors = 
		c("darkblue","lightblue","green","yellow","red"),
                         breaks=c(20,10,8,6,4,2), labels = c("+20-10","10-8","8-6","6-4","4-2","<2"), limits=c(0,20)) +
		scale_size(name="-log10(pvalue)",breaks=c(20,10,8,6,4,2), labels = c("+20-10","10-8","8-6","6-4","4-2","<2"), limits=c(0,20)) +
#		theme( axis.text.x = element_text(face = rev(Style), color = 'gray30', size = 14, hjust=0, angle=0),
		#axis.text.y = element_text(face = rev(Style), color = 'gray30', size = 14, hjust=0, angle=0),
#                    axis.title.x = element_text(size = 20, color = 'gray30', vjust = 0),
#                   ) +
#        geom_text(aes(label = gsub("\\s2", "", label), y = 0), hjust = 0) +
        guides(color=guide_legend(), size=guide_legend()) +
        scale_y_log10(breaks=ticks, labels = ticks) +
        geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
  		theme( axis.text.y = element_text(face = df$style, color = 'gray30', size = 10, hjust=0, angle=0)
		#axis.text.y = element_text(face = rev(Style), color = 'gray30', size = 14, hjust=0, angle=0),
#                    axis.title.x = element_text(size = 20, color = 'gray30', vjust = 0),
                   ) +
       xlab("") + ylab("OR (95% CI)") 
plot
}


types<-c("QTL","Binary")

output.files <- list.files(paste(output.path,"/",sep=""), pattern=paste("*univariable.summary.txt$",sep=""))
output.arg <- identical(output.files, character(0))
if(!output.arg){
  for(i in output.files){
     dta <- read.delim(paste(output.path,"/",i,sep=""),header=T)
           gene.dta <- subset(dta, TYPE2 %in% "GENE" & !is.na(P))
           if(nrow(gene.dta) != 0){
           sig.p.gene <- 0.05/nrow(gene.dta)
           sig.fdr<- 0.05
		   p_cut<-0.05
		   num_cut <- 50
           gene.name <- as.vector(unique(gene.dta$GENENAME))
           # generate data for plots for gene x SNP combinations
	    for(k in types){
           gene.dta.type <- subset(gene.dta, TYPE1 %in% k)
	     if(k == "QTL"){
   	       df<-prepare_dta.forestPlot_QT(gene.dta.type,gene.name,p_cut,num_cut)
               #png(filename=paste(output.path,"/plots/",as.character(j),"_",gene.name,".FPQ.png",sep=""),width=600,height=450,res=100)
                 if(!is.null(df)){ 
                 pdf(file=paste(output.path,"/plots/univariable.",gene.name,".FPQ.pdf",sep=""), height=8.5, width= 11)
                 plot<-forestPlot_QT(df)
                 print(plot)
                 dev.off()
                 }
             }else if(k == "Binary"){
               df<-prepare_dta.forestPlot_BI(gene.dta.type,gene.name,p_cut,num_cut)
                if(!is.null(df)){ 
                 pdf(file=paste(output.path,"/plots/univariable.",gene.name,".FPB.pdf",sep=""), height=8.5, width= 11)
                 plot<-forestPlot_BI(df)
                 print(plot)
                 dev.off()
             }
             }		
            }	   
			
		  }
         }
 
}
