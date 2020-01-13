## VarID is part of the RaceID (>= v0.1.4) package.
## This script contains the code to reproduce the anlysis of hematopoietic bone marrow single-cell RNA-seq data
## presented in Figures 1-3 of Grün D (2019) Revealing Dynamics of Gene Expression Variability in Cell State Space. Nature Methods 17(1):45-49.  doi: 10.1038/s41592-019-0632-3
## Data from Tusi, B. K. et al. Population snapshots predict early haematopoietic and erythroid hierarchies. Nature 555, 54–60 (2018)

## The analysis in this script was done with RaceID v0.1.6!! For letter versions cluster numbers and parameter settings potentially change.

require(RaceID)

## load input data
inputData <- readRDS("inputData_hematopoiesis.rds")


## Run VarID on normal bone marrow progenitors (see vignette for details: vignette("RaceID")).
## The results depend on a random seed. Although a seed argument is integrated, the same seed might lead to slightly
## different results depending on the architecture.

sN <- SCseq(inputData$wt)
sN <- filterdata(sN,mintotal=1000,CGenes=rownames(inputData$wt)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(inputData$wt))])
expData  <- getExpData(sN)

res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=10,alpha=10,no_cores=5,seed=12345)
cl    <- graphCluster(res,pvalue=0.01)
probs <- transitionProbs(res,cl)

## compute noise from corrected variance
x     <- as.matrix(sN@expdata)[sN@genes,colnames(sN@ndata)]
noise <- compNoise(x,res,regNB=FALSE,pvalue=0.01,no_cores=5)

sN <- updateSC(sN,res=res,cl=cl,noise=noise,flo=.1)
sN <- comptsne(sN)
sN <- compumap(sN)

## expression UMAP
plotmap(sN,um=TRUE)

## map of transition probabilities
plotTrProbs(sN,probs,tp=.5,prthr=0.01,cthr=0,um=TRUE)

## marker expression dot plot
genes <- c("Kit","Flt3","Dntt","Ebf1","Cd19","Lmo4","Ms4a2","Ear10","Cd74","Irf8","Mpo","Elane","Ngp","Mpl","Pf4","Car1","Gata1","Hbb.bs","Ptgfrn","Mki67")
##initialize cluster order (depends on random seed)
clusterV <- sort(unique(sN@cpart))
fractDotPlot(sN, genes, cluster=clusterV,zsc=TRUE,cap=3)

## background model for inference of pruned network
plotBackVar(res)

## baseline noise model
plotNoiseModel(noise)
plotNoiseModel(noise,corrected=TRUE)


## Identify differentially variable genes with increased variability in the Kit+ multipotent progenitor cluster.
## The cluster number depends on the random seed. Select the cluster with highest expression of Kit.
plotexpmap(sN,"Kit",um=TRUE)
kitHiCluster <- 19
ngenes <- diffNoisyGenes(noise,cl,kitHiCluster,no_cores=5)
## Compute differentially expressed genes in the Kit-high cluster.
dgenes <- clustdiffgenes(sN,kitHiCluster)
## Apply p-value and fold change cutoff
gV <- rownames(ngenes)[ngenes$pvalue < .001 & ngenes$log2FC > log2(1.25) ]
gD <- rownames(dgenes)[dgenes$padj < .001 & dgenes$fc > 1.25 ]

## overlap of up-regulated genes and genes with enhanced variability
require(VennDiagram)
## shut down open graphic devices
dev.off()
draw.pairwise.venn(length(gD),length(gV),length(intersect(gD,gV)))


## plot heatmap of gene expression and gene expression variability for the top 50 genes with increased variability in the Kit-high cluster
##initialize cluster order (depends on random seed)
clusterV <- sort(unique(sN@cpart))
## plot expression and store the order of rows
ph <- plotmarkergenes(sN,genes=head(gV,50),cl=clusterV,noise=FALSE,zsc=FALSE)
## plot variability keeping the same ordering of genes
plotmarkergenes(sN,genes=ph$tree_row$labels[ ph$tree_row$order ],cl=clusterV,noise=TRUE, cluster_rows=FALSE)

## plot expression of Gata1 and Mpo
plotexpmap(sN,"Gata1",um=TRUE,log=TRUE)
plotexpmap(sN,"Mpo",um=TRUE,log=TRUE)

## plot variability of Gata1 and Mpo
plotexpmap(sN,"Gata1",noise=TRUE,um=TRUE,log=TRUE)
plotexpmap(sN,"Mpo",noise=TRUE,um=TRUE,log=TRUE)

## compute noise from the residuals of a negative binomial regression for comparison
x      <- as.matrix(sN@expdata)[sN@genes,colnames(sN@ndata)]
noiseR <- compNoise(x,res,regNB=TRUE,pvalue=0.01,no_cores=5)
plotRegNB(expData,noiseR,"(Intercept)")
plotRegNB(expData,noiseR,"beta")
plotRegNB(expData,noiseR,"theta")
## residual dependence between variance and mean of Pearson's residuals has remained (therefore not a good option for this dataset)
plotPearsonRes(noiseR,log=TRUE)


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



## analysis of noise dynamics
## initialize with clusters on neutrophil branch (depends on random seed of initial analysis)
neutroClusters <- c(19,3,7,4)

n <- names(sN@cpart)[sN@cpart %in% neutroClusters]

## run StemID analysis
sl <- SCseq(sN@expdata[,n])
sl <- filterdata(sl,mintotal=1000,CGenes=rownames(data)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(data))])
sl <- compdist(sl)
sl <- clustexp(sl,samp=1000,cln=6)
sl <- findoutliers(sl,outminc=2)
sl <- compumap(sl)
sl <- comptsne(sl,perplexity=30)

ltr <- Ltree(sl)
ltr <- compentropy(ltr)
ltr <- projcells(ltr,cthr=5,nmode=TRUE,knn=3)
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr,pthr=0.05)

plotgraph(ltr)

## select trajectory from lineage tree; clusters depend on initial input (which was determined by a random seed)
clusterTrajectory <- c(6,5,10,4,2,7,1)
n <- cellsfromtree(ltr,clusterTrajectory)$f

k <- getfdata(sN,n=n)
u <- sN@noise[rownames(k),n]
y <- sl@cpart

## analyze pseudo-temporal expression profiles with FateID. See vignette for details: vignette("FateID").
require(FateID)
## plot expression profile of neutrophil marker
plotexpression(k,y=y,g="Elane",n=n, col=rep("lightgrey",length(sN@fcol)),cex=1,alpha=1)

## discard non-expressed genes
fs   <- filterset(k,n=n,minexpr=3,minnumber=1)

## create self-organizing map (SOM) of pseudo-temporal expression profiles
s1d  <- getsom(fs,nb=200,alpha=.5)
ps   <- procsom(s1d,corthr=.85,minsom=10)

## plot SOM of pseudo-temporal expression profiles highlighting the original clusters
y    <- sN@cpart[n] 
plotheatmap(ps$all.z, xpart=y[n], xcol=sN@fcol, ypart=ps$nodes, xgrid=FALSE, ygrid=TRUE, xlab=TRUE)


## create self-organizing map of pseudo-temporal variability profiles for the same genes using the same order of cells
fn <- sN@noise[rownames(fs),colnames(fs)]
s1dn <- getsom(fn,nb=200,alpha=.5)
pn   <- procsom(s1dn,corthr=.85,minsom=10)

## plot SOM of pseudo-temporal varibility profiles highlighting the original clusters
plotheatmap(pn$all.z, xpart=y[n], xcol=sN@fcol, ypart=pn$nodes, xgrid=FALSE, ygrid=TRUE, xlab=TRUE)

## plot noise and expression profiles for a given module i, e.g. i=2

u <- sN@noise[rownames(k),n]
w <- k[,n]
w <- t(t(w)/apply(w,1,sum))
u <- t(t(u)/apply(u,1,sum))

i <- 2
g <- names(pn$nodes)[pn$nodes==i]
## noise
plotexpression(u,y=sN@cpart[n] ,g=g,n=n,name=paste("Module:",i),col=rep("lightgrey",length(sN@fcol)),cex=1,alpha=1)
## expression
plotexpression(w,y=sN@cpart[n] ,g=g,n=n,name=paste("Module:",i),col=rep("lightgrey",length(sN@fcol)),cex=1,alpha=1)


## Reactome Pathway analysis of module i
require(ReactomePA)
require(org.Mm.eg.db)
require(mGSZ)
react <- function(set,universe=NULL,pv=0.05){
    conv <- AnnotationDbi::toTable(org.Mm.egSYMBOL)
    ac <- aggregate(as.numeric(conv$gene_id),by=list(conv$symbol),min)
    k <- ac$x
    names(k) <- ac$Group.1
    gene_id <- as.character(k[set[set %in% names(k)]])
    if ( !is.null(universe)){
        universe_id <- as.character(k[universe[universe %in% names(k)]])
        enrichPathway(gene=gene_id,organism = "mouse",pvalueCutoff=pv, readable=T, universe=universe_id)
    }else{
        enrichPathway(gene=gene_id,organism = "mouse",pvalueCutoff=pv, readable=T)
    }
}
i <- 3
universe <- rownames(getfdata(sN))
raM <-react( names(pn$nodes)[pn$nodes == i], universe=universe, pv=0.05)

## modify Description for convenient plotting
raMplot <- raM
raMplot@result$Description <- sapply(raM@result$Description,FUN=function(x){substr(x,1,50)})      
barplot(raMplot, showCategory=30)




## Co-analysis of normal and EPO-stimulated bone marrow cells

dataA <- cbind(inputData$wt,inputData$EPO)

sA <- SCseq(dataA)
sA <- filterdata(sA,mintotal=1000,CGenes=rownames(dataA)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(dataA))])
expDataA <- getExpData(sA)
resA     <- pruneKnn(expDataA,knn=10,alpha=10,no_cores=5)
clA      <- graphCluster(resA)
probsA   <- transitionProbs(resA,clA)

x <- as.matrix(sA@expdata)[sA@genes,colnames(sA@ndata)]
noiseA <- compNoise(x,resA,no_cores=5)
sA    <- updateSC(sA,res=resA,cl=clA,noise=noiseA,flo=.1)
sA    <- compumap(sA)

plotmap(sA,um=TRUE)

## plot samples in UMAP
types <- sub("B.+","",colnames(sA@ndata))
plotsymbolsmap(sA,types,um=TRUE)

## compute genes with enhanced variability in EPO-stimulated versus normal erythrocyte progenitors
## initialize clusters (depends on random seed)
normCluster <- 15
epoCluster <- 17
ngenesA <- diffNoisyGenes(noiseA,clA,epoCluster,normCluster,no_cores=5)
gA <- rownames(ngenesA)[ngenesA$pvalue < .001 & ngenesA$log2FC > log2(1.25) ]
nhA <- head(gA,50)

## plot heatmap of variability and expression for the top 50 genes ordered by increased variability between EPO-stimulated and normal erythrocyte progenitors
## initialize cluster order (depends on random seed)
clusterA <- c(4,14,16,9,5,15,17,3,11,13,18)
## plot expression and store the order of rows
ph <- plotmarkergenes(sA,genes=nhA,cl=clusterA,samples=types,noise=FALSE,zsc=FALSE)
## plot variability keeping the same ordering of genes
plotmarkergenes(sA,genes=ph$tree_row$labels[ ph$tree_row$order ],cl=clusterA,samples=types,noise=TRUE,cluster_rows=FALSE)

## compute differentially expressed genes between EPO-stimulated and normal erythrocyte progenitors
x <- diffexpnb(getfdata(sA),A=names(sA@cpart)[sA@cpart == normCluster],B=names(sA@cpart)[sA@cpart == epoCluster])
f <- x$res$padj < .001 & x$res$log2FoldChange > log2(1.25)
plotdiffgenesnb(x)
gDE <- rownames(x$res)[f]


require(VennDiagram)
dev.off()
draw.pairwise.venn(length(gDE),length(gA),length(intersect(gDE,gA)))

## perform Reactome Pathway analysis of differentially variable genes
require(ReactomePA)
require(org.Mm.eg.db)
require(mGSZ)
react <- function(set,universe=NULL,pv=0.05){
    conv <- AnnotationDbi::toTable(org.Mm.egSYMBOL)
    ac <- aggregate(as.numeric(conv$gene_id),by=list(conv$symbol),min)
    k <- ac$x
    names(k) <- ac$Group.1
    gene_id <- as.character(k[set[set %in% names(k)]])
    if ( !is.null(universe)){
        universe_id <- as.character(k[universe[universe %in% names(k)]])
        enrichPathway(gene=gene_id,organism = "mouse",pvalueCutoff=pv, readable=T, universe=universe_id)
    }else{
        enrichPathway(gene=gene_id,organism = "mouse",pvalueCutoff=pv, readable=T)
    }
}

## differentially variable genes
z <-  react( gA, universe=universe,pv=0.05)
universe <- rownames(getfdata(sA))
barplot(z, showCategory=30)
## differentially expressed genes
zd <-  react( gDE, universe=universe,pv=0.05)
barplot(zd, showCategory=30)

