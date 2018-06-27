
# reading the expression data and the cluster info
miRNAclusters = read.delim("miRNA_clusters_16.txt",header = FALSE, stringsAsFactors = FALSE)
colnames(miRNAclusters) = c("Symbol","cluster.number")

mRNAclusters = read.delim("mRNA_clusters_30.txt",stringsAsFactors = FALSE,header = FALSE) # these clusters are with stringTie ids: MSTRG#
colnames(mRNAclusters) = c("Symbol","cluster.number")

mRNAtpms <- read.table("gene.tpm.modified.txt", header=T, row.names = "Symbol")

miRNAtpms <- read.csv("CPMs.TMMnormalized.modified.csv", header=T, row.names = "Symbol")
miRNAtpms <- miRNAtpms + 0.1
miRNAtpms[,colnames(miRNAtpms) %in% paste0("t",seq(1,48))] = NA

# Reading the desin matrix and sample file
edesignfile <- read.csv("miRNA.maSigPro.design.csv", header=TRUE, stringsAsFactors = F)
samples <- read.csv("miRNA_samples.all.modified.csv", header=TRUE, stringsAsFactors = F)

colors = cbind(c("forebrain","#800000"),c("midbrain","#B22222"),c("hindbrain","#DC143C"),c("neural.tube","#FF8C00"),c("cranioface","gold"),c("liver", "#9370DB"),c("heart","#31a354"),c("kidney","#0000FF"),c("lung","green2"),c("intestine","#1E90FF"),c("stomach","#40E0D0"),c("limb","#EE82EE"))
colnames(colors) = colors[1,]

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

tissues = unique(edesignfile$tissue)
stdevMat = as.data.frame(matrix(0, max(miRNAclusters$cluster.number),length(tissues)),colnames = tissues)
for(i in 1:length(unique(miRNAclusters$cluster.number))){
  for(j in 1:length(tissues)){
    s = (which(samples$tissue == tissues[j]))
    s = round((s[seq(1,length(s),2)]+1)/2)
    stdevMat[i,j] = sd(miRNA.tpms[s,i],na.rm = TRUE)
  }
}


library("gplots")
png("tissueSpecificity.png", width=5, height=5, units="in", res=300)
heatmap.2(apply(stdevMat,2,function(x) as.numeric(x)),scale = "row",Rowv = FALSE, Colv = FALSE,colsep = 1:12,rowsep = 1:max(miRNAclusters$cluster.number),ColSideColors = c(colors[2,as.character(unique(tissues))]), sepwidth = c(0.02,0.05), sepcolor = "purple",na.color = "grey",margin=c(8, 4),xlab = "Tissues", ylab = "Clusters",main = "Tissue specificity matrix",col = bluered,tracecol = NA,key = FALSE,labCol = as.character(unique(tissues)))
dev.off()



# defining tissue specificity matrix
scaled_stdev = scale(t(stdevMat),center = TRUE)
scaled_stdev[scaled_stdev>=0]<- 1
scaled_stdev[scaled_stdev<0] <- 0
tissue.specificity = c()
for(i in 1:12){
  tissue.specificity = c(tissue.specificity, rep(scaled_stdev[i,],8))
}
tissue.specificity = matrix(tissue.specificity,nrow =8*12,ncol = length(unique(miRNAclusters$cluster.number)),byrow = TRUE)


tissue.specificity[tissue.specificity == 1] = TRUE
tissue.specificity[is.na(tissue.specificity)] = FALSE
corMat = matrix(0, ncol = 30,nrow = length(unique(miRNAclusters$cluster.number)))
for(i in 1:length(unique(miRNAclusters$cluster.number))){
  X = miRNA.tpms[as.logical(tissue.specificity[,i]),i]
  corMat[i,] = cor(X,mRNA.tpms[as.logical(tissue.specificity[,i]),],use = "pairwise.complete.obs")
}


png("SelectiveCorrelation.png",  width=9, height=6, units="in", res=300)
heatmap.2(apply(corMat,2,function(x) as.numeric(x)),margins = c(5,5),Rowv = FALSE, Colv = FALSE,trace="none",dendrogram = 'none',density.info = "none",colsep = 1:30,rowsep = 1:length(unique(miRNAclusters$cluster.number)), sepwidth = c(0.02,0.05), sepcolor = "black",na.color = "grey",xlab = "mRNA Clusters", ylab = "miRNA Clusters",col = colorpanel(5,"#f03b20","#bcbddc","#dadaeb"),tracecol = NA,main = "Partial correlation of miRNA-mRNA clusters",labCol = paste0("gC",1:30),labRow = paste0("miC",1:length(unique(miRNAclusters$cluster.number))),key = TRUE,keysize = 1,key.par=list(mar=c(3.5,0,3,0)),key.xtickfun=function() {
  breaks <- parent.frame()$breaks
  return(list(
    at=parent.frame()$scale01(c(breaks[1],
                                breaks[length(breaks)])),
    labels=c(as.character(breaks[1]),
             as.character(breaks[length(breaks)]))
  ))
},key.xlab = "Correlation Coefficient",key.title = "Pearson Correlation")
dev.off()

save(tissue.specificity,corMat,miRNA.tpms,mRNA.tpms,file = "correlations_16miRC.RData")
