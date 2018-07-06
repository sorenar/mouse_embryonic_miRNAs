
library(ggplot2)
require(gridExtra)

#reading the gene and miRNA clusters with parsing to match the targetfinding method
miRNAclusters = read.delim("miRNA_clusters_16.txt",header = FALSE, stringsAsFactors = FALSE)
colnames(miRNAclusters) = c("Symbol","cluster.number")

mRNAclusters = read.delim("mRNA_clusters_30.txt",stringsAsFactors = FALSE,header = FALSE)
colnames(mRNAclusters) = c("Symbol","cluster.number")

mRNAtpms <- read.table("gene.tpm.modified.txt", header=T, row.names = "Symbol")

miRNAtpms <- read.csv("CPMs.TMMnormalized.modified.csv", header=T, row.names = "Symbol")
miRNAtpms <- miRNAtpms + 0.1 # pseudo counts
miRNAtpms[,colnames(miRNAtpms) %in% paste0("t",seq(1,48))] = NA  # dealing with the missing data labeled as t1 to t48

# averaging the experessions for each cluster at each of the time-stage points
miRNA.tpms = matrix(0, ncol = 16,nrow = 96)  
for(i in 1:length(unique(miRNAclusters$cluster.number))){
  for(j in 1:96){
    miRNA.tpms[j,i] = median(apply(miRNAtpms[miRNAclusters[which(miRNAclusters$cluster.number == i),1],(2*j-1):(2*j)],1,function(x) mean(as.numeric(x))),na.rm = TRUE)
  }
}

mRNA.tpms = matrix(0, ncol = 30,nrow = 96)
for(i in 1:30){
  for(j in 1:96){
    mRNA.tpms[j,i] = median(apply(mRNAtpms[mRNAclusters[which(mRNAclusters$cluster.number == i),1],(2*j-1):(2*j)],1,function(x) mean(as.numeric(x))),na.rm = TRUE)
  }
}


# loading the correlation matrix (corMat), tissue specificity (tissue.specificity) and the cluster means of mRNA and miRNA expressions (mRNA.tpms and miRNA.tpms)
load("~/Documents/MortazaviLab/Integration/FinalCodes/correlations_16miRC.RData")

#loading the target prediction matrix (PVals and targets)
load("~/Documents/MortazaviLab/Integration/FinalCodes/predictionTargets.RData")

# selecting the interaction with significant enrichment
sigInteractions = arrayInd(which(PVals < 0.05/30/length(unique(miRNAclusters$cluster.number))),c(length(unique(miRNAclusters$cluster.number)),30))

samples <- read.csv("miRNA_samples.all.modified.csv", header=TRUE, stringsAsFactors = F)
Tissues <- factor(samples$tissue[seq(1,192,2)],levels=c("forebrain","midbrain","hindbrain","neural.tube",
                                                        "cranioface","liver",
                                                        "heart","kidney","lung","intestine"
                                                        ,"stomach","limb"))
Stages <- rep(c(10.5,11.5,12.5,13.5,14.5,15.5,16.5,18.5),12)
colors = cbind(c("forebrain","#800000"),c("midbrain","#B22222"),c("hindbrain","#DC143C"),c("neural.tube","#FF8C00"),c("cranioface","gold"),c("liver", "#9370DB"),c("heart","#31a354"),c("kidney","#0000FF"),c("lung","green2"),c("intestine","#1E90FF"),c("stomach","#40E0D0"),c("limb","#EE82EE"))
colnames(colors) = colors[1,]
colors = colors[,unique(Tissues)]

# iterating through the significant interactions with negative correlation to plot the miRNA-mRNA expression side-by-side                               
for(i in 1:length(unique(miRNAclusters$cluster.number))){
  for(j in 1:length(unique(mRNAclusters$cluster.number))){
    if(PVals[i,j]<0.05/30/length(unique(miRNAclusters$cluster.number)) && corMat[i,j]<0){
      df <- data.frame(tissue = Tissues, time= Stages,miRNA = miRNA.tpms[,i],mRNA = mRNA.tpms[,j])
      df <- na.omit(df)
      p <- ggplot(df, aes(x=time, y=miRNA, group=tissue, color=tissue)) + 
        geom_line() +
        geom_line(data = df[df$tissue %in% Tissues[tissue.specificity[,i]==1],],size = 1.5) +
        scale_y_log10()+
        labs(y="Median profile (logCPM)",x="Developmental stage",title = paste0("miRNA Cluster",i," (",length(unique(miRNAclusters[miRNAclusters$cluster.number == i,1])) ," miRNAs)"))
      
      p <- p + theme(text = element_text(size = 20))
      p <- p + theme_bw()
      p <- p + scale_color_manual(values = colors[2,], guide = F)
      p <- p + theme(axis.text=element_text(size=9),
                     axis.title=element_text(size=12),plot.title = element_text(size=12) )
      p <- p + scale_x_discrete(limit = c(10.5,11.5, 12.5, 13.5, 14.5, 15.5, 16.5,18.5),
                                labels = c("e10.5","e11.5","e12.5", "e13.5", "e14.5", "e15.5", "e16.5","P0"), expand = c(0.05, 0.05))
      print(p)
      assign("mifig1",p)
      
      
      
      p <- ggplot(df, aes(x=time, y=mRNA, group=tissue, color=tissue)) + 
        geom_line() +
        geom_line(data = df[df$tissue %in% Tissues[tissue.specificity[,i]==1],],size = 1.5)+
        scale_y_log10()+
        labs(y="Median profile (logCPM)",x="Developmental stage",title = paste0("mRNA Cluster",j," (",length(unique(mRNAclusters[mRNAclusters$cluster.number == j,1])) ," genes)"))
      p <- p + theme(text = element_text(size = 20))
      p <- p + theme_bw()
      p <- p + scale_color_manual(values = colors[2,], guide = F)
      p <- p + theme(axis.text=element_text(size=9),
                     axis.title=element_text(size=12),plot.title = element_text(size=12) )
      p <- p + scale_x_discrete(limit = c(10.5,11.5, 12.5, 13.5, 14.5, 15.5, 16.5,18.5),
                                labels = c("e10.5","e11.5","e12.5", "e13.5", "e14.5", "e15.5", "e16.5","P0"), expand = c(0.05, 0.05))
      print(p)
      assign("mfig1",p)
      grid.arrange(mifig1,mfig1, ncol=2)
      dev.print(pdf, paste0('miR',i,'-','gene',j,'_log.pdf'),width = 7.5, height = 3)
    }
  }
}

######################### GO analysis
geneMap = read.delim("StringTie_geneNames_CommonNames.txt",stringsAsFactors = FALSE,header = FALSE)
colnames(geneMap) = c("Symbol","Ensambl","Name")

# writting the input files for the metascape analysis at "http://metascape.org/gp/index.html#/main/step1"
for(n in 1:dim(sigInteractions)[1]){
  targetlist = targets[[sigInteractions[n,1]]][[sigInteractions[n,2]]][-1,]
  mRNAs = unique(targetlist$Symbol)
  write.table(unique(geneMap[which(geneMap$Symbol %in% mRNAs),3]),paste0("target",sigInteractions[n,1],"-",sigInteractions[n,2],".csv"),sep = ",",quote = FALSE,row.names = FALSE,col.names = FALSE)
}

# Once the GO analysis was run on the Metascape server the result table was saved with the following format and read for enrichment analysis
# "miC#-gC#_metascape_result.xlsx"                              
library(readxl)

# using the output of the Metascape analysis to plot the enriched terms 
for(n in 1:dim(sigInteractions)[1]){
  miR = sigInteractions[n,1]
  gene = sigInteractions[n,2]
  test = read_excel(paste0(miR,"-",gene,"_metascape_result.xlsx"), sheet = 2)
  test = test[!duplicated(test$Description),]
  test = test[order(test$LogP,decreasing = TRUE),]
  test$Description = factor(test$Description,levels = test$Description)
  test = test[order(test$LogP,decreasing = FALSE),]
  test$colors[test$LogP < -4] = "#dadaeb"
  test$colors[test$LogP < -6] = "#9e9ac8"
  test$colors[test$LogP < -12] = "#54278f"
  png(paste0(miR,"-",gene,"_metascape_GO1.png"), width=9, height=6, units="in", res=300)
  ggplot(test[1:15,], aes(x=Description, y=-LogP)) +
    geom_bar(stat='identity',fill = rev(test$colors[1:15]),color="black") +
    scale_y_continuous(limits=c(0,17))+
    coord_flip()+
    theme(axis.text = element_text(size = 10),axis.title=element_text(size=12))
  dev.off()
}



