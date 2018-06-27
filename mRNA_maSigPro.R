# Running maSigPro on miRNAs in hpc
library(maSigPro)
library(ggplot2)


edesignfile <- read.csv("mRNA.maSigPro.design.csv", header=TRUE, stringsAsFactors = F)


## make experimental design
time <- edesignfile$time
replicates <- edesignfile$reps
forebrain <- edesignfile$forebrain
midbrain <- edesignfile$midbrain
hindbrain <- edesignfile$hindbrain
cranioface <- edesignfile$cranioface
neural.tube <- edesignfile$neural.tube
liver <- edesignfile$liver
heart <- edesignfile$heart
limb <- edesignfile$limb
stomach <- edesignfile$stomach
intestine <- edesignfile$intestine
kidney <- edesignfile$kidney
lung <- edesignfile$lung

samples <- edesignfile$Sample
ss.edesign <- cbind(time,replicates,forebrain, midbrain,
                    hindbrain, cranioface, neural.tube,
                    liver, heart
                    , limb, stomach,
                    intestine, kidney, lung)
rownames(ss.edesign) <- paste(edesignfile$Sample)
design <- make.design.matrix(ss.edesign, degree=3)



## read the miRNA file and filter the miRNAs with zero variation across the experiments
data.counts = read.delim("gene.counts.txt", row.names="Symbol")
# remove the spike ins
data.counts = data.counts[,colnames(data.counts) %in% edesignfile$Sample]
data.filtered = data.counts[apply(data.counts,1,var,na.rm = TRUE) != 0,]
data.t <- data.filtered + 0.1


datap <- p.vector(data.t, design = design, counts = FALSE) # find the significant genes
datat <- T.fit(datap, alfa = 0.01) # find the significant differences

# saving the maSigpro regressed data
save(datat,file = "mRNA_maSigPro_TMM.RData")

get<-get.siggenes(datat,rsq=0.7, vars="all") # get the list of significant genes
x <- see.genes(get$sig.genes,k=30)$cut
y <- data.frame(x)
write.table(y,"mRNA_clusters_30.txt", sep="\t", quote=FALSE, col.names=FALSE)

colors = cbind(c("forebrain","#800000"),c("midbrain","#B22222"),c("hindbrain","#DC143C"),c("neural.tube","#FF8C00"),c("cranioface","gold"),c("liver", "#9370DB"),c("heart","#31a354"),c("kidney","#0000FF"),c("lung","green2"),c("intestine","#1E90FF"),c("stomach","#40E0D0"),c("limb","#EE82EE"))
colors = t(colors)
rownames(colors) = colors[,1]
samples <- read.csv("mRNA_samples.all.modified.csv", header=TRUE, stringsAsFactors = F)
tpms <- read.delim("gene.tpm.modified.txt", header=T, row.names = "Symbol")
tpms <- tpms + 0.1

tpms[,colnames(tpms) %in% paste0("t",seq(1,48))] = NA

# plotting the miRNA cluster profiles
for (clust in 1:30){
  cluster <- rownames(y)[y$x == clust]
  cluster.data <- tpms[cluster,]
  med <- apply(cluster.data, 2, median)
  
  tissues <- c()
  medians <- c()
  std.err <- c()
  
  for (tissue in unique(edesignfile$tissue)){
    tissue.medians <- med[samples$tissue == tissue]
    list.split <- split(tissue.medians, ceiling(seq_along(tissue.medians)/2))
    f.d <- c()
    se <- c()
    for (i in list.split){
      f.d <- c(f.d, mean(i))
      se <- c(se, sd(i)/sqrt(2))
    }
    tissues <- c(tissues , rep(tissue,8))
    medians <- c(medians, f.d)
    std.err <- c(std.err, se)
  }
  Tissues <- Tissue <- factor(tissues, 
                              levels=c("forebrain","midbrain","hindbrain","neural.tube",
                                       "cranioface","liver",
                                       "heart","kidney","lung","intestine"
                                       ,"stomach","limb"))
  df <- data.frame(tissue = Tissues, time=rep(c(10.5,11.5,12.5,13.5,14.5,15.5,16.5,18.5),12),
                   medians = medians, se = std.err)
  df <- na.omit(df)
  
  clust.name <- plot_num <- paste("cluster",clust,": ",length(cluster), " genes",sep = "")
  
  p <- ggplot(df, aes(x=time, y=medians, group=tissue, color = tissue)) + 
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=medians-se, ymax=medians+se), width=.8,
                  position=position_dodge(0.0))+ 
    #labs(y="",x="",title = clust.name) 
    labs(y="Median profile (TPM)",x="Developmental stage",title = clust.name)
  
  p <- p + theme(text = element_text(size = 20))
  #p <- p + geom_line(aes(size = 0))
  #p <- p + theme_bw()
  
  p <- p + theme_bw()
  p <- p + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.line = element_line(colour = "black"))
  
  p <- p + scale_color_manual(values = as.character(colors[unique(df$tissue),2]), guide = F)
  
  p <- p + theme(axis.text=element_text(size=20),
                 axis.title=element_text(size=20, face = "bold"),plot.title = element_text(size=25) )
  p <- p + scale_x_discrete(limit = c(10.5,11.5, 12.5, 13.5, 14.5, 15.5, 16.5,18.5),
                            labels = c("e10.5","e11.5","e12.5", "e13.5", "e14.5", "e15.5", "e16.5","P0"), expand = c(0.05, 0.05))
  print(p)
  plot_num <- paste("cluster",clust,".pdf",sep="")
  assign(paste("fig",clust,sep=""),p)
}

require(gridExtra)
grid.arrange(fig1,fig2,fig3,fig4,fig5,fig6,fig7,fig8,fig9,fig10,fig11,fig12,
             fig13,fig14,fig15,fig16,fig17,fig18,fig19,fig20,fig21,fig22,fig23,
             fig24,fig25,fig26,fig27,fig28,fig29,fig30,
             ncol=4)
dev.print(pdf, paste0('miRNA_clusters_',16,'.pdf'),width = 50, height = 60)
