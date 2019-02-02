library(ggplot2)

options(echo = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  output.path <- args[1]
} else {
  stop("Not enough arguments")
}

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

draw_plot <- function(data.plot,sig.p,sig.fdr,cutoff.OR,cutoff.qtl,g,v){
  data.plot.logp.outside <- subset(data.plot, logp > 50)
  data.plot.logp <- subset(data.plot, logp <= 50)
  if(nrow(data.plot.logp.outside) > 0){
     par(xpd = TRUE, mar =c(5.5,4.5,9,5), mfrow = c(4,1))
  }else{
     par(mar =c(5.5,4.5,4.5,5), mfrow = c(4,1))
  }

  data.plot.logp.outside <- data.plot.logp.outside[order(data.plot.logp.outside$logp),]
  ylim <-  c(0, max(c(-1*log10(sig.p)+2, data.plot.logp$logp), na.rm = T))
  plot( data.plot.logp$x, data.plot.logp$logp, xlab = "", xlim=c(0,max(data.plot.logp$x)+1),col = data.plot.logp$color, cex = data.plot.logp$cex, pch = data.plot.logp$pch,xaxt = "n",  ylab =expression(SNP: -log[10](italic(p))), ylim = ylim)
  abline(h = -1*log10(sig.p), col = "grey",xpd=F)
  pos <- aggregate(data.plot.logp$x, by = list(data.plot.logp$TAG), median)
  text(pos[[2]]-0.1, par("usr")[3] - 0.5, labels = pos[[1]], srt = -45, pos = 4, xpd = TRUE,  cex = 0.8)  
  #axis(1, at = pos[[2]], labels = pos[[1]], las =1, cex = 1.2 )
  data.plot.logp.sig <- data.plot.logp[!is.na(data.plot.logp$logp) &data.plot.logp$logp >= -1*log10(sig.p), ]
  max.p <- max(data.plot.logp$logp,na.rm=T)
 # print(data.plot.logp.sig)
  if(!is.null(data.plot.logp.sig) & nrow(data.plot.logp.sig)>0) {
	    text(data.plot.logp.sig$x, data.plot.logp.sig$logp, labels = data.plot.logp.sig$pheno,cex = 0.8, pos = data.plot.logp.sig$pos, offset = 0.25)  
  }
  if(nrow(data.plot.logp.outside) > 0){
    title(paste(v,sprintf("[ %s, %1.1f%% ]",data.plot.logp$effa[1],as.numeric(data.plot.logp$freq[1])*100)," in ",g), line=8)
    for(i in 1:nrow(data.plot.logp.outside)){
       points(data.plot.logp.outside$x[i], max(c(-1*log10(sig.p)+2, data.plot.logp$logp))+(i*1.3+2),col = data.plot.logp.outside$color[i], cex = data.plot.logp.outside$cex[i], pch = data.plot.logp.outside$pch[i])
       text(data.plot.logp.outside$x[i]+1, max(c(-1*log10(sig.p)+2, data.plot.logp$logp))+(i*1.3+2), labels = paste(data.plot.logp.outside$pheno[i],"(",round(data.plot.logp.outside$logp[i]),")",sep=""),cex = 0.8, pos = data.plot.logp.outside$pos[i], offset = 0.25)
    }
  }else{
    title(paste(v,sprintf("[ %s, %1.1f%% ]",data.plot.logp$effa[1],as.numeric(data.plot.logp$freq[1])*100)," in ",g), line=3)
  }
  ylim <-  c(0, max(c(-1*log10(sig.p)+2, data.plot.logp$logp), na.rm = T))
  plot( data.plot.logp$x, data.plot.logp$logp, xlab = "", xlim=c(0,max(data.plot.logp$x)+1),col = data.plot.logp$color, cex = data.plot.logp$cex, pch = data.plot.logp$pch,xaxt = "n",  ylab =expression(SNP: -log[10](italic(p))), ylim = ylim)
  abline(h = -1*log10(sig.p), col = "grey",xpd=F)
  pos <- aggregate(data.plot.logp$x, by = list(data.plot.logp$TAG), median)
  text(pos[[2]]-0.1, par("usr")[3] - 0.5, labels = pos[[1]], srt = -45, pos = 4, xpd = TRUE,  cex = 0.8)  
  if(nrow(data.plot.logp.outside) > 0){
    title(paste(v,sprintf("[ %s, %1.1f%% ]",data.plot.logp$effa[1],as.numeric(data.plot.logp$freq[1])*100)," in ",g), line=8)
    for(i in 1:nrow(data.plot.logp.outside)){
       points(data.plot.logp.outside$x[i], max(c(-1*log10(sig.p)+2, data.plot.logp$logp))+(i*1.3+2),col = data.plot.logp.outside$color[i], cex = data.plot.logp.outside$cex[i], pch = data.plot.logp.outside$pch[i])
#       text(data.plot.logp.outside$x[i]+1, max(c(-1*log10(sig.p)+2, data.plot.logp$logp))+(i*1.3+2), labels = paste(data.plot.logp.outside$pheno[i],"(",round(data.plot.logp.outside$logp[i]),")",sep=""),cex = 0.8, pos = data.plot.logp.outside$pos[i], offset = 0.25)
    }
  }else{
    title(paste(v,sprintf("[ %s, %1.1f%% ]",data.plot.logp$effa[1],as.numeric(data.plot.logp$freq[1])*100)," in ",g), line=3)
  }

 
  data.plot.logFDR.outside <- subset(data.plot, logp > 50)
  data.plot.logFDR <- subset(data.plot, logp <= 50)

  data.plot.logFDR.outside <- data.plot.logFDR.outside[order(data.plot.logFDR.outside$logFDR),]
  ylim <- c(0, max(c(-1*log10(sig.fdr)+2, data.plot.logFDR$logFDR), na.rm = T))
  plot( data.plot.logFDR$x, data.plot.logFDR$logFDR, xlab = "", xlim=c(0,max(data.plot.logFDR$x)+1), 
        col = data.plot.logFDR$color, cex = data.plot.logFDR$cex, pch = data.plot.logFDR$pch,
        xaxt = "n",  ylab =expression(-log[10](italic(FDR))),  ylim = ylim)
  abline(h = -1*log10(sig.fdr), col = "grey",xpd=F)
  pos <- aggregate(data.plot.logFDR$x, by = list(data.plot.logFDR$TAG), median)
  text(pos[[2]]-0.1, par("usr")[3] - 0.2 , labels = pos[[1]], srt = -45, pos = 4, xpd = TRUE,  cex = 0.8)  
  #axis(1, at = pos[[2]], labels = pos[[1]], las =1, cex = 1.2 )
  data.plot.logFDR.sig <- data.plot.logFDR[!is.na(data.plot.logFDR$logFDR) &data.plot.logFDR$logFDR >= -1*log10(sig.fdr), ]
  max.p <- max(data.plot.logFDR$logFDR,na.rm=T)
 # print(data.plot.logFDR.sig)
  if(!is.null(data.plot.logFDR) & nrow(data.plot.logFDR.sig)>0) {
   data.plot.logFDR.sig1 <- subset(data.plot.logFDR.sig, type2 == "SNP" )
   text(data.plot.logFDR.sig$x, data.plot.logFDR.sig$logFDR, labels = data.plot.logFDR.sig$pheno,cex = 0.8, pos = data.plot.logFDR.sig$pos, offset = 0.25)  
  }
  if(nrow(data.plot.logFDR.outside) > 0){
    title(paste(v,sprintf("[ %s, %1.1f%% ]",data.plot.logFDR$effa[1],as.numeric(data.plot.logFDR$freq[1])*100)," in ",g), line=8)
    for(i in 1:nrow(data.plot.logFDR.outside)){
       points(data.plot.logFDR.outside$x[i], max(c(-1*log10(sig.fdr)+2, data.plot.logFDR$logFDR))+(i*1.3+2),col = data.plot.logFDR.outside$color[i], cex = data.plot.logFDR.outside$cex[i], pch = data.plot.logFDR.outside$pch[i])
       text(data.plot.logFDR.outside$x[i]+1, max(c(-1*log10(sig.fdr)+2, data.plot.logFDR$logFDR))+(i*1.3+2), labels = paste(data.plot.logFDR.outside$pheno[i],"(",round(data.plot.logFDR.outside$logFDR[i]),")",sep=""),cex = 0.8, pos = data.plot.logFDR.outside$pos[i], offset = 0.25)
    }
  }else{
    title(paste(v,sprintf("[ %s, %1.1f%% ]",data.plot.logFDR$effa[1],as.numeric(data.plot.logFDR$freq[1])*100)," in ",g), line=3)
  }

  ylim <- c(0, max(c(-1*log10(sig.fdr)+2, data.plot.logFDR$logFDR), na.rm = T))
  plot( data.plot.logFDR$x, data.plot.logFDR$logFDR, xlab = "", xlim=c(0,max(data.plot.logFDR$x)+1), 
        col = data.plot.logFDR$color, cex = data.plot.logFDR$cex, pch = data.plot.logFDR$pch,
        xaxt = "n",  ylab =expression(-log[10](italic(FDR))), ylim = ylim)
  abline(h = -1*log10(sig.fdr), col = "grey",xpd=F)
  pos <- aggregate(data.plot.logFDR$x, by = list(data.plot.logFDR$TAG), median)
  text(pos[[2]]-0.1, par("usr")[3] - 0.2 , labels = pos[[1]], srt = -45, pos = 4, xpd = TRUE,  cex = 0.8)  
  if(nrow(data.plot.logFDR.outside) > 0){
    title(paste(v,sprintf("[ %s, %1.1f%% ]",data.plot.logFDR$effa[1],as.numeric(data.plot.logFDR$freq[1])*100)," in ",g), line=8)
    for(i in 1:nrow(data.plot.logFDR.outside)){
       points(data.plot.logFDR.outside$x[i], max(c(-1*log10(sig.fdr)+2, data.plot.logFDR$logFDR))+(i*1.3+2),col = data.plot.logFDR.outside$color[i], cex = data.plot.logFDR.outside$cex[i], pch = data.plot.logFDR.outside$pch[i])
#       text(data.plot.logFDR.outside$x[i]+1, max(c(-1*log10(sig.fdr)+2, data.plot.logFDR$logFDR))+(i*1.3+2), labels = paste(data.plot.logFDR.outside$pheno[i],"(",round(data.plot.logFDR.outside$logFDR[i]),")",sep=""),cex = 0.8, pos = data.plot.logFDR.outside$pos[i], offset = 0.25)
    }
  }else{
    title(paste(v,sprintf("[ %s, %1.1f%% ]",data.plot.logFDR$effa[1],as.numeric(data.plot.logFDR$freq[1])*100)," in ",g), line=3)
  }

  if(unique(data.plot.logp$type2[!is.na(data.plot.logp$type2)]) == "SNP"){
    add_legend("topright",c(paste(cutoff.OR, "(OR) or ", cutoff.qtl, "(MD)",sep = ""),"SNP: OR<1 or MD<0","SNP: OR>1 or MD>0"),cex = 0.8,pch = c(rep(16, length(cutoff.OR)),1, 16),pt.cex = c(1:length(cutoff.OR),1,1)+0.6, bty = "n", xjust = 0, y.intersp = 1,title = "Effect Magnitude",title.adj=0, ncol = 2)
  }else if(unique(data.plot.logp$type2[!is.na(data.plot.logp$type2)]) == "SCORE"){
    add_legend("topright",c(paste(cutoff.OR, "(OR) or ", cutoff.qtl, "(MD)",sep = ""),"SCORE: OR<1 or MD<0","SCORE: OR>1 or MD>0"),cex = 0.8,pch = c(rep(16, length(cutoff.OR)),1, 16),pt.cex = c(1:length(cutoff.OR),1,1)+0.6, bty = "n", xjust = 0, y.intersp = 1,title = "Effect Magnitude",title.adj=0, ncol = 2)
  } 
}



prepare_dta.forestPlot_QT <- function(dta,snp.name,gene.name,p_cut,num_cut){

x=subset(dta, SNP %in% snp.name  & P < p_cut & P != 0)
if(nrow(x) > 0){
x=x[!duplicated(x),]
x=x[order(x$P),]
if(nrow(x) < num_cut) {
  x=x[1:nrow(x),]
}else {
  x=x[1:num_cut,]
}
cluster=as.data.frame(aggregate(P ~ Category, x,function(x){min(x)}))
cluster$rank=rank(-cluster$P)
cluster$P=NULL
x=merge(x,cluster,by="Category")
x=x[order(rank(x$Category),-x$P),]
x=x[order(x$rank),]
#x$BETA <- ifelse(x$analysed.Freq1 > 0.5, -1*x$BETA,x$BETA)
#x$L95 <- ifelse(x$analysed.Freq1 > 0.5, -1*x$L95,x$L95)
#x$U95 <- ifelse(x$analysed.Freq1 > 0.5, -1*x$U95,x$U95)
x$freq <- x$analysed.Freq1 
x$effa<-as.character(x$effall)

label <- factor(sprintf("%s [ %5.0f]",x$PHENOTYPE,x$N))
mean <- x$BETA
lower<- x$L95
upper<- x$U95
pvalue <- -log10(x$P)
pvalue <-ifelse(pvalue > 20, 20,pvalue)
title <-paste("Plot of ",x$SNP[1]," ", sprintf("[ %s, %1.1f%% ]",x$effa[1],x$freq[1]*100)," in ",x$GENENAME[1],sep="")

df <- data.frame(label, mean, lower, upper, pvalue, x$Category, title)
df$style="plain"
#df$just=1
#insert categories
df$ind <- seq_len(nrow(df))
row.max.category=as.data.frame(aggregate(ind ~ x.Category, df,function(x){max(x)}))
for(i in 1:nrow(row.max.category)){
    df <- rbind(df,data.frame(label = row.max.category$x.Category[i], mean = NA,lower = NA, upper=NA, pvalue=NA, x.Category=NA, title=NA, style="bold",ind=row.max.category$ind[i]+0.5))
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


prepare_dta.forestPlot_BI <- function(dta,snp.name,gene.name,p_cut,num_cut){
#p_cut: pvalue threshold to show
#num_cut: # of traits to show
x=subset(dta, SNP %in% snp.name & P < p_cut & P != 0)
if(nrow(x) >0){
x=x[!duplicated(x),]
x=x[order(x$P),]
if(nrow(x) < num_cut) {
  x=x[1:nrow(x),]
}else {
  x=x[1:num_cut,]
}
cluster=as.data.frame(aggregate(P ~ Category, x,function(x){min(x)}))
cluster$rank=rank(-cluster$P)
cluster$P=NULL
x=merge(x,cluster,by="Category")
x=x[order(rank(x$Category),-x$P),]
x=x[order(x$rank),]
#x$OR <- ifelse(x$analysed.Freq1 > 0.5, 1/x$OR,x$OR)
#x$OR <- exp(x$BETA)
#x$L95 <- ifelse(x$analysed.Freq1 > 0.5, 1/x$L95,x$L95)
#x$U95 <- ifelse(x$analysed.Freq1 > 0.5, 1/x$U95,x$U95)
x$freq <- x$analysed.Freq1
x$effa<-as.character(x$effall)


label <- factor(sprintf("%s [%5.0f, %5.0f]",x$PHENOTYPE,x$CASES,x$CONT))
OR <- x$OR
lower<- x$L95
upper<- x$U95
pvalue <- -log10(x$P)
pvalue <-ifelse(pvalue > 20, 20,pvalue)
#title <-paste("Plot of ",x$snp[1]," in ",x$gene[1],sep="")
title <-paste("Plot of ",x$SNP[1]," ", sprintf("[ %s, %1.1f%% ]",x$effa[1],x$freq[1]*100)," in ",x$GENENAME[1],sep="")

df <- data.frame(label, OR, lower, upper, pvalue, x$Category, title)
df$style="plain"
#insert categories
df$ind <- seq_len(nrow(df))
row.max.category=as.data.frame(aggregate(ind ~ x.Category, df,function(x){max(x)}))
for(i in 1:nrow(row.max.category)){
    df <- rbind(df,data.frame(label = row.max.category$x.Category[i], OR = NA,lower = NA, upper=NA, pvalue=NA, x.Category=NA, title=NA, style="bold",ind=row.max.category$ind[i]+0.5))
}
df <- df[order(df$ind),]
df$label <- factor(df$label, levels=df$label)
df
}else{
df<-NULL
}
}


forestPlot_BI <- function(df){
title<-df$title[1]
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


prepare_dta.plot <- function(sub.dta,cols,cutoff.OR,cutoff.qtl){
           #TYPE1: QTL or Binary
           #TYPE2: SNP or GENE
           sub.dta <- sub.dta[,c("Category","PHENOTYPE","TYPE1","TYPE2","P","FDR","BETA","analysed.Freq1","effall")]
           sub.dta <- subset(sub.dta, Category != "")
#           sub.dta$TISSUE <- strtrim(sub.dta$TISSUE,6)
           #threshold.gene
           
           sub.dta$TAG = strtrim(gsub(" ","_",sub.dta$Category),20)
           sub.dta$pheno = strtrim(gsub(" ","_",sub.dta$PHENOTYPE),20)
#           sub.dta$pheno = ifelse(!is.na(sub.dta$TISSUE),paste(sub.dta$PHENOTYPE,sub.dta$TISSUE,sep=":"),as.character(sub.dta$PHENOTYPE))
           sub.dta$logp <- -1*log10(sub.dta$P)
           sub.dta$logFDR <- -1*log10(sub.dta$FDR)
           sub.dta$freq <-sub.dta$analysed.Freq1 
           sub.dta$effa <-as.character(sub.dta$effall)

           data.plot<- sub.dta[c("TAG", "pheno","TYPE1","TYPE2","logp","logFDR","BETA")]
           names(data.plot) <- c("TAG", "pheno","type1","type2","logp","logFDR","effect")

           count <- table(data.plot$TAG)
           fillmore <- names(count)[count<25]
           for( k in fillmore) {
              for( j in (count[k]+1):20)
                 data.plot[nrow(data.plot)+1,]<- c(k,NA,NA,NA, NA, NA, NA, NA, NA)  
           }
           data.plot$effect <- as.numeric(data.plot$effect)
           data.plot$logp <- as.numeric(data.plot$logp)  
           data.plot$logFDR <- as.numeric(data.plot$logFDR)
           data.plot<- data.plot[order(data.plot$TAG, -1* data.plot$logp), ]
           data.plot$pos <- 4              
           data.plot$x <- 1:nrow(data.plot)
           #assign size by effect size
           data.plot$cex <- NA
           index<- (!is.na(data.plot$type1)) & data.plot$type1 != "QTL"
           data.plot$cex[index]<- c(1:length(cutoff.OR))[as.numeric(cut(abs(data.plot$effect[index]),breaks = c(log(cutoff.OR), Inf), include.lowest = T, right = F))]+0.6
           index<- (!is.na(data.plot$type1)) &data.plot$type1 == "QTL"
           data.plot$cex[index]<- c(1:length(cutoff.qtl))[as.numeric(cut(abs(data.plot$effect[index]),breaks = c(cutoff.qtl, Inf), include.lowest = T, right = F))]+0.6

           ##assign point by effect direction
           data.plot$pch <- c(1,16)[ifelse(as.numeric(data.plot$effect)>=0, 2, 1)]
  
           #assign color by phenotype TAG
           data.plot$color <- cols[as.numeric(as.factor(data.plot$TAG))]
		   data.plot
}

types<-c("QTL","Binary")
rainbowCols<-c("#FF0000FF", "#FFAA00FF", "#AAFF00FF", "#00FF00FF", "#00FFAAFF", "#00AAFFFF","#0000FFFF", "#AA00FFFF", "#FF00AAFF")

output.files <- list.files(paste(output.path,"/",sep=""), pattern=paste("*[.]summary[.]txt$",sep=""))
output.arg <- identical(output.files, character(0))
if(!output.arg){
  for(i in output.files){
     dta <- read.delim(paste(output.path,"/",i,sep=""),header=T)
     dta <- subset(dta,!(Category %in% c("Brain imaging","Procedures")))
           dta$otherall<-ifelse(as.character(dta$otherall) %in% "TRUE","T",as.character(dta$otherall))
           dta$effall<-ifelse(as.character(dta$effall) %in% "TRUE","T",as.character(dta$effall))

           snp.dta <- subset(dta, TYPE2 %in% "SNP" & !is.na(P))
           if(nrow(snp.dta) != 0){
          
           sig.fdr<- 0.05
           cols <- rep(rainbowCols, 30)
           cutoff.OR<- 1:3
           cutoff.qtl<- c(0, 0.5, 1)
		   p_cut<-0.05
		   num_cut <- 20
           SNP.names <- unique(snp.dta$SNP[!is.na(snp.dta$SNP)]) 
		   gene.name <- unlist(lapply( strsplit(i,"[.]"),function(x)x[[1]]))
           # generate data for plots for gene x SNP combinations
           for(j in SNP.names){
	     sub.snp.dta <- subset(snp.dta, SNP %in% j)
             sig.p.snp <- 0.05/nrow(sub.snp.dta)
             sub.snp.dta$FDR <- p.adjust(sub.snp.dta$P,method="BH")
 	    for(k in types){
             sub.snp.dta.type <- subset(sub.snp.dta, TYPE1 %in% k)			
	     if(k == "QTL"){
   	       df<-prepare_dta.forestPlot_QT(sub.snp.dta.type,as.character(j),gene.name,p_cut,num_cut)
               #png(filename=paste(output.path,"/plots/",as.character(j),"_",gene.name,".FPQ.png",sep=""),width=600,height=450,res=100)
                 if(!is.null(df)){ 
                 pdf(file=paste(output.path,"/plots/",as.character(j),".",gene.name,".FPQ.pdf",sep=""), height=8.5, width= 11)
                 plot<-forestPlot_QT(df)
                 print(plot)
                 dev.off()
                 }
             }else if(k == "Binary"){
               df<-prepare_dta.forestPlot_BI(sub.snp.dta.type,as.character(j),gene.name,p_cut,num_cut)
               #png(filename=paste(output_path,"/plots/",x$snp[1],"_",x$GENENAME[1],".FPB.png",sep=""),width=650,height=450,res=100)
                if(!is.null(df)){ 
                 pdf(file=paste(output.path,"/plots/",as.character(j),".",gene.name,".FPB.pdf",sep=""), height=8.5, width= 11)
                 plot<-forestPlot_BI(df)
                 print(plot)
                 dev.off()
                }
             }
            }			   
	     data.plot <- prepare_dta.plot(sub.snp.dta,cols,cutoff.OR,cutoff.qtl)
             data.plot[mapply(is.infinite, data.plot)] <- 300

             v <- as.character(j)
 #            gene.name <- unlist(lapply( strsplit(i,"[.]"),function(x)x[[1]]))
             g <- as.character(gene.name)

             file.out <-paste(output.path,"/plots/",v,".",g,".MP.pdf",sep="")
             pdf(file=file.out, height=16, width= 11)
		       if(sum(!is.na(sub.snp.dta$GENENAME) != 0)){
			     g <- unique(sub.snp.dta$GENENAME)
			   }else{
		         g <- "N/A"	 
			   }
	     draw_plot(data.plot,sig.p.snp,sig.fdr,cutoff.OR,cutoff.qtl,g,v)	
             dev.off()
	
		  }
         }
   }
}

