# miRNAtap compiles the predicted targets of the miRNAs for 
library(miRNAtap)
library(gplots)

#reading the mRNA and miRNA clusters
miRNAclusters = read.delim("miRNA_clusters_16.txt",header = FALSE, stringsAsFactors = FALSE)
colnames(miRNAclusters) = c("Symbol","cluster.number")

mRNAclusters = read.delim("mRNA_clusters_30.txt",stringsAsFactors = FALSE,header = FALSE)
colnames(mRNAclusters) = c("Symbol","cluster.number")
# Reading the Ensembl to ID mapping file
entrezMap = read.delim("gencode.vM10.metadata.EntrezGene",stringsAsFactors = FALSE,header = FALSE)
colnames(entrezMap) = c("Symbol","Entrez.ID")

# geneMap has the Ensambl ids for the mRNAs
geneMap = read.delim("StringTie_geneNames_CommonNames.txt",stringsAsFactors = FALSE,header = FALSE)
colnames(geneMap) = c("Symbol","Ensambl","Name")


z = merge(mRNAclusters,geneMap, by = "Symbol")
mRNAclusters = merge(z,entrezMap, by.x = 3,by.y = 1)

# targets is a list with the size of miC# x gC which holds the target list for each miCxgC interaction
targets = list()
for (i in 1:length(unique(miRNAclusters$cluster.number))){
  targets[[i]] = list()
  for (l in 1:length(unique(mRNAclusters$cluster.number))){
    targets[[i]][[l]] = data.frame(matrix(ncol = 8))
    colnames(targets[[i]][[l]]) = c("Entrez.ID", "miR.Symbol", "rank_product","rank_final","Ensambl","Symbol","cluster.number","Name")
  }
  
  # going through the miRNAs in each cluster to get the predicted targets and compile them
  for (mirna in miRNAclusters[miRNAclusters$cluster.number == i,1]) {
    # miRNAtap uses 5 different resources for target prediction
    p <- getPredictedTargets(mirna, species = 'mmu',method = 'geom',min_src = 3,synonyms = TRUE)
    if(is.null(p)){ 
      print(paste("jumping for no prediction for:",mirna))
      next
    }
    p = as.data.frame(p)
    p$rank_final =  p$rank_final/max(p$rank_final)
    r = cbind(rownames(p),mirna,p[,c("rank_product","rank_final")])
    colnames(r) = c("Entrez.ID", "miR.Symbol","rank_product","rank_final")
    
    # distributing the targets of each miRNA between different mRNA clusters
    for (j in 1:length(unique(mRNAclusters$cluster.number))){
      mRNAcluster = mRNAclusters[mRNAclusters$cluster.number == j,]
      c = merge(r,mRNAcluster,by = "Entrez.ID")
      if(!all(isEmpty(c))){
        targets[[i]][[j]] = rbind(targets[[i]][[j]], c)
      }
      
    }
    
  }
  
}
# counting the unique number of target mRNAs in each interactions
mRNAt = data.frame(matrix(ncol = length(unique(mRNAclusters$cluster.number)),nrow = length(unique(miRNAclusters$cluster.number))))
mRNAn = c()
for (i in 1:length(unique(miRNAclusters$cluster.number))){
  for (j in 1:length(unique(mRNAclusters$cluster.number))){
    mRNAn[j] = length(mRNAclusters[mRNAclusters$cluster.number == j,1]) 
    mRNAt[i,j] = length(unique(targets[[i]][[j]]$Entrez.ID))
  }
}

# calculating the contingency table and the PValues for the Chi-square test

PVals = data.frame(matrix(ncol = length(unique(mRNAclusters$cluster.number)),nrow = length(unique(miRNAclusters$cluster.number))))
for (i in 1:length(unique(miRNAclusters$cluster.number))){
  for (j in 1:length(unique(mRNAclusters$cluster.number))){
    
    contigencyTable = rbind(c(mRNAt[i,j],sum(mRNAt[i,])-mRNAt[i,j]),c(mRNAn[j]-mRNAt[i,j],sum(mRNAn[-j])-sum(mRNAt[i,])+mRNAt[i,j]))
    test = chisq.test(contigencyTable)
    PVals[i,j] = test$p.value
    
  }
}


# plotting the -log10 of the P-value 

x = apply(-log10(PVals),2,function(x) as.numeric(as.character(x)))
breaks = c(0,4.2,6,8,10,12)
png("chisq-Pvals.png", width=9, height=6, units="in", res=300)
par(mar=c(4,4,1,1))
rownames(x) = paste0("miC",1:length(unique(miRNAclusters$cluster.number)))
colnames(x) = paste0("gC",1:length(unique(mRNAclusters$cluster.number)))
heatmap.2(x,Rowv = NA,Colv = NA,trace="none",dendrogram = 'none',density.info = "none",colsep = 0:length(unique(mRNAclusters$cluster.number)),rowsep = 0:length(unique(miRNAclusters$cluster.number)),sepcolor = "black", col = colorpanel(5,"#dadaeb","#9e9ac8","#54278f"),key = TRUE,xlab = "mRNA Clusters", ylab = "miRNA Clusters",key.ylab = "P-Values(-log10)",keysize = 1,key.par=list(mar=c(3.5,0,3,0)),key.xtickfun=function() {
  breaks <- parent.frame()$breaks
  return(list(
    at=parent.frame()$scale01(c(breaks[1],
                                breaks[length(breaks)])),
    labels=c(as.character(breaks[1]),
             as.character(breaks[length(breaks)]))
  ))
},key.xlab = "P-Values",main = "Target Enrichment (Chi-sq)")
dev.off() 

save(targets,mRNAt,mRNAt,PVals,file = "predictionTargets.RData")
