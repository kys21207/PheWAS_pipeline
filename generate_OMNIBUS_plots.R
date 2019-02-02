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


draw_plot <- function(data.plot,sig.p,sig.fdr,cutoff.OR,cutoff.qtl,v){
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
    title(paste(v), line=8)
    for(i in 1:nrow(data.plot.logp.outside)){
       points(data.plot.logp.outside$x[i], max(c(-1*log10(sig.p)+2, data.plot.logp$logp))+(i*1.3+2),col = data.plot.logp.outside$color[i], cex = data.plot.logp.outside$cex[i], pch = data.plot.logp.outside$pch[i])
       text(data.plot.logp.outside$x[i]+1, max(c(-1*log10(sig.p)+2, data.plot.logp$logp))+(i*1.3+2), labels = paste(data.plot.logp.outside$pheno[i],"(",round(data.plot.logp.outside$logp[i]),")",sep=""),cex = 0.8, pos = data.plot.logp.outside$pos[i], offset = 0.25)
    }
  }else{
    title(paste(v), line=3)
  }
  ylim <-  c(0, max(c(-1*log10(sig.p)+2, data.plot.logp$logp), na.rm = T))
  plot( data.plot.logp$x, data.plot.logp$logp, xlab = "", xlim=c(0,max(data.plot.logp$x)+1),col = data.plot.logp$color, cex = data.plot.logp$cex, pch = data.plot.logp$pch,xaxt = "n",  ylab =expression(SNP: -log[10](italic(p))), ylim = ylim)
  abline(h = -1*log10(sig.p), col = "grey",xpd=F)
  pos <- aggregate(data.plot.logp$x, by = list(data.plot.logp$TAG), median)
  text(pos[[2]]-0.1, par("usr")[3] - 0.5, labels = pos[[1]], srt = -45, pos = 4, xpd = TRUE,  cex = 0.8)  
  if(nrow(data.plot.logp.outside) > 0){
    title(paste(v), line=8)
    for(i in 1:nrow(data.plot.logp.outside)){
       points(data.plot.logp.outside$x[i], max(c(-1*log10(sig.p)+2, data.plot.logp$logp))+(i*1.3+2),col = data.plot.logp.outside$color[i], cex = data.plot.logp.outside$cex[i], pch = data.plot.logp.outside$pch[i])
#       text(data.plot.logp.outside$x[i]+1, max(c(-1*log10(sig.p)+2, data.plot.logp$logp))+(i*1.3+2), labels = paste(data.plot.logp.outside$pheno[i],"(",round(data.plot.logp.outside$logp[i]),")",sep=""),cex = 0.8, pos = data.plot.logp.outside$pos[i], offset = 0.25)
    }
  }else{
    title(paste(v), line=3)
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
    title(paste(v), line=8)
    for(i in 1:nrow(data.plot.logFDR.outside)){
       points(data.plot.logFDR.outside$x[i], max(c(-1*log10(sig.fdr)+2, data.plot.logFDR$logFDR))+(i*1.3+2),col = data.plot.logFDR.outside$color[i], cex = data.plot.logFDR.outside$cex[i], pch = data.plot.logFDR.outside$pch[i])
       text(data.plot.logFDR.outside$x[i]+1, max(c(-1*log10(sig.fdr)+2, data.plot.logFDR$logFDR))+(i*1.3+2), labels = paste(data.plot.logFDR.outside$pheno[i],"(",round(data.plot.logFDR.outside$logFDR[i]),")",sep=""),cex = 0.8, pos = data.plot.logFDR.outside$pos[i], offset = 0.25)
    }
  }else{
    title(paste(v), line=3)
  }

  ylim <- c(0, max(c(-1*log10(sig.fdr)+2, data.plot.logFDR$logFDR), na.rm = T))
  plot( data.plot.logFDR$x, data.plot.logFDR$logFDR, xlab = "", xlim=c(0,max(data.plot.logFDR$x)+1), 
        col = data.plot.logFDR$color, cex = data.plot.logFDR$cex, pch = data.plot.logFDR$pch,
        xaxt = "n",  ylab =expression(-log[10](italic(FDR))), ylim = ylim)
  abline(h = -1*log10(sig.fdr), col = "grey",xpd=F)
  pos <- aggregate(data.plot.logFDR$x, by = list(data.plot.logFDR$TAG), median)
  text(pos[[2]]-0.1, par("usr")[3] - 0.2 , labels = pos[[1]], srt = -45, pos = 4, xpd = TRUE,  cex = 0.8)  
  if(nrow(data.plot.logFDR.outside) > 0){
    title(paste(v), line=8)
    for(i in 1:nrow(data.plot.logFDR.outside)){
       points(data.plot.logFDR.outside$x[i], max(c(-1*log10(sig.fdr)+2, data.plot.logFDR$logFDR))+(i*1.3+2),col = data.plot.logFDR.outside$color[i], cex = data.plot.logFDR.outside$cex[i], pch = data.plot.logFDR.outside$pch[i])
#       text(data.plot.logFDR.outside$x[i]+1, max(c(-1*log10(sig.fdr)+2, data.plot.logFDR$logFDR))+(i*1.3+2), labels = paste(data.plot.logFDR.outside$pheno[i],"(",round(data.plot.logFDR.outside$logFDR[i]),")",sep=""),cex = 0.8, pos = data.plot.logFDR.outside$pos[i], offset = 0.25)
    }
  }else{
    title(paste(v), line=3)
  }

#  if(unique(data.plot.logp$type2[!is.na(data.plot.logp$type2)]) == "SNP"){
#    add_legend("topright",c(paste(cutoff.OR, "(OR) or ", cutoff.qtl, "(MD)",sep = ""),"SNP: OR<1 or MD<0","SNP: OR>1 or MD>0"),cex = 0.8,pch = c(rep(16, length(cutoff.OR)),1, 16),pt.cex = c(1:length(cutoff.OR),1,1)+0.6, bty = "n", xjust = 0, y.intersp = 1,title = "Effect Magnitude",title.adj=0, ncol = 2)
#  }else if(unique(data.plot.logp$type2[!is.na(data.plot.logp$type2)]) == "SCORE"){
#    add_legend("topright",c(paste(cutoff.OR, "(OR) or ", cutoff.qtl, "(MD)",sep = ""),"SCORE: OR<1 or MD<0","SCORE: OR>1 or MD>0"),cex = 0.8,pch = c(rep(16, length(cutoff.OR)),1, 16),pt.cex = c(1:length(cutoff.OR),1,1)+0.6, bty = "n", xjust = 0, y.intersp = 1,title = "Effect Magnitude",title.adj=0, ncol = 2)
#  } 
}


prepare_dta.plot <- function(sub.dta,cols,cutoff.OR,cutoff.qtl){
           #TYPE1: QTL or Binary
           #TYPE2: SNP or GENE or SCORE
           sub.dta <- sub.dta[,c("Category","PHENOTYPE","TYPE1","TYPE2","P","FDR","BETA")]
           sub.dta <- subset(sub.dta, Category != "")
#           sub.dta$TISSUE <- strtrim(sub.dta$TISSUE,6)
           #threshold.gene
           
           sub.dta$TAG = strtrim(gsub(" ","_",sub.dta$Category),20)
           sub.dta$pheno = strtrim(gsub(" ","_",sub.dta$PHENOTYPE),20)
           sub.dta$logp <- -1*log10(sub.dta$P)
           sub.dta$logFDR <- -1*log10(sub.dta$FDR)
           data.plot<- sub.dta[c("TAG", "pheno","TYPE1","TYPE2","logp","logFDR","BETA")]
           names(data.plot) <- c("TAG", "pheno","type1","type2","logp","logFDR","effect")

           count <- table(data.plot$TAG)
           fillmore <- names(count)[count<25]
           for( k in fillmore) {
              for( j in (count[k]+1):20)
                 data.plot[nrow(data.plot)+1,]<- c(k,NA,NA,NA, NA, NA, NA)  
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
           data.plot$cex[index]<- c(1:length(cutoff.OR))[as.numeric(cut(abs(data.plot$effect[index]),breaks = c(log(cutoff.OR), Inf), include.lowest = T, right = F))]+0.5
           index<- (!is.na(data.plot$type1)) &data.plot$type1 == "QTL"
           data.plot$cex[index]<- c(1:length(cutoff.qtl))[as.numeric(cut(abs(data.plot$effect[index]),breaks = c(cutoff.qtl, Inf), include.lowest = T, right = F))]+0.5

           ##assign point by effect direction
           data.plot$pch <- c(1,16)[ifelse(as.numeric(data.plot$effect)>=0, 2, 1)]
  
           #assign color by phenotype TAG
           data.plot$color <- cols[as.numeric(as.factor(data.plot$TAG))]
		   data.plot
}

types<-c("QTL","Binary")
rainbowCols<-c("#FF0000FF", "#FFAA00FF", "#AAFF00FF", "#00FF00FF", "#00FFAAFF", "#00AAFFFF","#0000FFFF", "#AA00FFFF", "#FF00AAFF")

output.files <- list.files(paste(output.path,"/",sep=""), pattern=paste("*omnibus[.]summary[.]txt$",sep=""))
output.arg <- identical(output.files, character(0))
if(!output.arg){
  for(i in output.files){
     dta <- read.delim(paste(output.path,"/",i,sep=""),header=T)
     dta <- subset(dta,!(Category %in% c("Brain imaging","Procedures")))
     gene.dta <- subset(dta, TYPE2 %in% "GENE" & !is.na(P))
     if(nrow(gene.dta) != 0){
        gene.dta$FDR <- p.adjust(gene.dta$P,method="BH")
        sig.p.snp <- 0.05/nrow(gene.dta)
        sig.fdr<- 0.05
        cols <- rep(rainbowCols, 30)
        cutoff.OR<- 1:3
        cutoff.qtl<- c(0, 0.5, 1)
	    p_cut<-0.05
		num_cut <- 20
        sub.gene.dta <- gene.dta[,c("Category","PHENOTYPE","TYPE1","TYPE2","BETA","P","FDR")]
#        sub.gene.dta$BETA <- 1		
	data.plot <- prepare_dta.plot(sub.gene.dta,cols,cutoff.OR,cutoff.qtl)
        data.plot[mapply(is.infinite, data.plot)] <- NA
        v <- as.character(unique(gene.dta$GENENAME))
        file.out <-paste(output.path,"/plots/",v,".omnibus.MP.pdf",sep="")
        pdf(file=file.out, height=16, width= 11)
        draw_plot(data.plot,sig.p.snp,sig.fdr,cutoff.OR,cutoff.qtl,v)	
        dev.off()
        sub.gene.dta <- gene.dta[,c("Category","PHENOTYPE","TYPE1","TYPE2","BETA","Pavg","FDRavg")]
#        sub.gene.dta$BETA <- 1	
        names(sub.gene.dta)<-c("Category","PHENOTYPE","TYPE1","TYPE2","P","FDR","BETA")		
        data.plot <- prepare_dta.plot(sub.gene.dta,cols,cutoff.OR,cutoff.qtl)
        data.plot[mapply(is.infinite, data.plot)] <- 300
        file.out <-paste(output.path,"/plots/",v,".omnibus.avg.MP.pdf",sep="")
        pdf(file=file.out, height=16, width= 11)
        draw_plot(data.plot,sig.p.snp,sig.fdr,cutoff.OR,cutoff.qtl,v)	
        dev.off()
     }
  }
}
