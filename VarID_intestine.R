## VarID is part of the RaceID 0.1.4 package.
## This script contains the code to reproduce the anlysis of intestinal epithelial single-cell RNA-seq data
## presented in Figures 4 of Grun XXX
## Data from Haber, A. L. et al. A single-cell survey of the small intestinal epithelium. Nature 551, 333–339 (2017)


require(RaceID)

## load input data
inputData <- readRDS("inputData_intestine.rds")


## Run VarID on mouse intestinal epitheial cells (see vignette for details: vignette("RaceID")).
## The results depend on a random seed. Although a seed argument is integrated, the same seed might lead to slightly
## different results depending on the architecture.


## discard replicates from female mice to avoid strong batch effect
types <- sub("_.+","",colnames(inputData$wt))
f <- ! types %in% c("B1","B2") 
sN <- SCseq(inputData$wt[,f])
sN <- filterdata(sN,mintotal=1000,CGenes=c("Mki67",rownames(inputData$wt)[grep("Rp(l|s)|^Gm\\d",rownames(inputData$wt))]))
expData  <- getExpData(sN)
res   <- pruneKnn(expData,large=TRUE,regNB=FALSE,knn=10,alpha=10,no_cores=5)
cl    <- graphCluster(res,pvalue=0.01)
probs <- transitionProbs(res,cl)

## compute noise from corrected variance
x <- as.matrix(sN@expdata)[sN@genes,colnames(sN@ndata)]
noise <- compNoise(x,res,regNB=FALSE,pvalue=0.01,no_cores=5)

sN  <- updateSC(sN,res=res,cl=cl,noise=noise,flo=.1)
sN <- comptsne(sN)
sN <- compumap(sN)

## expression UMAP
plotmap(sN,um=TRUE)

## map of transition probabilities
plotTrProbs(sN,probs,tp=.5,prthr=0,cthr=0,um=TRUE)

## marker expression dot plot
genes <- c("Lgr5","Hmgb2","Cox6a1","Krt19","Reg3b","Fabp6","Apoa1","Alpi","Dclk1","Sst","Gcg","Chgb","Sox4","Dll1","Defa24","Lyz1","Muc2","Agr2","Zg16","Mki67")
##initialize cluster order (depends on random seed)
clusterV <- sort(unique(sN@cpart))
fractDotPlot(sN, genes, cluster=clusterV, zsc=TRUE,cap=3)

## background model for inference of pruned network
plotBackVar(res)

## baseline noise model
plotNoiseModel(noise)
plotNoiseModel(noise,corrected=TRUE)


## compute noise from the residuals of a negative binomial regression for comparison
x <- as.matrix(sN@expdata)[sN@genes,colnames(sN@ndata)]
noiseR <- compNoise(x,res,regNB=TRUE,pvalue=0.01,no_cores=5)
plotRegNB(expData,noiseR,"(Intercept)")
plotRegNB(expData,noiseR,"beta")
plotRegNB(expData,noiseR,"theta")
## no residual dependence between variance and mean of Pearson's residuals has remained (could be used alternatively)
plotPearsonRes(noiseR,log=TRUE)

## Identify differentially variable genes with increased variability in the Lgr5+ intestinal stem cell cluster.
## The cluster number depends on the random seed. Select the cluster with highest expression of Lgr5.
plotexpmap(sN,"Lgr5",um=TRUE)
Lgr5HiCluster <- 11
ngenes <- diffNoisyGenes(noise,cl,Lgr5HiCluster,no_cores=5)
## differentially variable genes derived based on variability computed from the negative binomial regression
ngenesR <- diffNoisyGenes(noiseR,cl,Lgr5HiCluster,no_cores=5)

## Compute differentially expressed genes in the Lgr5-high cluster.
dgenes <- clustdiffgenes(sN,Lgr5HiCluster)
## Apply p-value and fold change cutoff
gV <- rownames(ngenes)[ngenes$pvalue < .001 & ngenes$log2FC > log2(1.25) ]
gR <- rownames(ngenesR)[ngenesR$pvalue < .001 & ngenesR$log2FC > log2(1.25) ]
gD <- rownames(dgenes)[dgenes$padj < .001 & dgenes$fc > 1.25 ]

## overlap of genes with enhanced variability computed with the corrected variance- and the negative binomial regression-based method
require(VennDiagram)
## shut down open graphic devices
dev.off()
draw.pairwise.venn(length(gV),length(gR),length(intersect(gV,gR)))

## overlap of up-regulated genes and genes with enhanced variability
## shut down open graphic devices
dev.off()
draw.pairwise.venn(length(gV),length(gD),length(intersect(gV,gD)))


## plot heatmap of gene expression and gene expression variability for the top 50 genes with increased variability in the Kit-high cluster
##initialize cluster order (depends on random seed)
clusterV <- sort(unique(sN@cpart))
## plot expression and store the order of rows
ph <- plotmarkergenes(sN,genes=head(gV,50),cl=clusterV,noise=FALSE,zsc=FALSE)
## plot variability keeping the same ordering of genes
plotmarkergenes(sN,genes=ph$tree_row$labels[ ph$tree_row$order ],cl=clusterV,noise=TRUE, cluster_rows=FALSE)

## plot expression of Foxa3, Hopx, Sox4, Tox3
plotexpmap(sN,"Foxa3",um=TRUE,log=TRUE)
plotexpmap(sN,"Hopx",um=TRUE,log=TRUE)
plotexpmap(sN,"Sox4",um=TRUE,log=TRUE)
plotexpmap(sN,"Tox3",um=TRUE,log=TRUE)

## plot variability of Gata1 and Mpo
plotexpmap(sN,"Foxa3",noise=TRUE,um=TRUE,log=TRUE)
plotexpmap(sN,"Hopx",noise=TRUE,um=TRUE,log=TRUE)
plotexpmap(sN,"Sox4",noise=TRUE,um=TRUE,log=TRUE)
plotexpmap(sN,"Tox3",noise=TRUE,um=TRUE,log=TRUE)


## transcription factor network analysis, TF annotation from AnimalTFDB 3.0

tfg <- intersect(gV,inputData$tf$Symbol)

require(parallel)
require(GENIE3)
require(doParallel)
no_cores <- detectCores()
g3 <- GENIE3(getfdata(sN)[tfg,],nCores=no_cores - 2,verbose=TRUE)

gF <- g3
for ( i in 1:ncol(gF) ){ gF[is.na(gF[,i]) | is.nan(gF[,i]),i] <- 0 }
thr <- .08
gF <- gF * ( gF > thr )
require(igraph)
gG <- graph_from_adjacency_matrix(t(gF), mode = "directed", diag = FALSE, weighted = TRUE)
plot.igraph(gG, edge.arrow.size=0.05, vertex.size = 15, vertex.label.cex = .5)

