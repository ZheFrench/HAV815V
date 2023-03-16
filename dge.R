#################################################################
#
# date: April 02, 2021
# platform: Ubuntu 10.04
# R.version : 4.2.1
# author: Villemin Jean-Philippe
# team: Bioinformatique et biologie des systèmes du cancer : J. Colinge 
# Institute : IRCM
#
# dge.R
# Usage : 
# 
# dge.R
# 
# Description : 
#
#
#################################################################

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(glue))

expression <- read.csv("./data/UnknowDataset.csv",sep="\t",stringsAsFactors=F)
rownames(expression) <- expression$X
expression$X <- NULL

colnames(expression)

conditions <- gsub(names(expression),pattern="_rep..?",replacement="",perl=T)


good <- apply(expression,1,function(x) sum(x>5))>=3
sum(good)
expression <- expression[good,]

dge <- DGEList(expression,genes=rownames(expression),group=factor(conditions))
rownames(dge$counts) <- rownames(expression)

dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge, verbose=TRUE)

# Another way to get normalized count is to use sweep function
# and to divide by the calcNormFactors
# Should be able to retrieve this from JC code somewhere.
cpm(dge)
#The edgeR documentation advises to use cpm or rpkm to export normalized values out of edgeR.
# cpm uses TMM normalization factors automatically.
# cpm stands for counts per million

write.table(cpm(dge),file="./data/norm-counts.txt",col.names =NA,quote=FALSE,row.names=TRUE,sep="\t")


dge$group <- factor(conditions)
#groups <- c(rep(levels(dge$group)[1],3),rep(levels(dge$group)[2],3))
phenotypes <- data.frame(SAMPLES=names(expression),GROUPS=conditions)
write.table(phenotypes,file="./data/phenotypes.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

design <- model.matrix(~0+dge$group) # no intercept #x0 = 1,  force model throught the origin

comp <- c("A-B")
cm <- makeContrasts(contrasts=comp,levels=dge$group)

y     <- estimateDisp(dge,design,robust=T)

fit.y <- glmFit(y,design)
lrt   <- glmLRT(fit.y,contrast = cm)

# No filter
sel.r <- topTags(lrt,p.value=1,adjust.method="BH",n = nrow(y$counts))

write.table(sel.r$table,file="./data/all-gene-ratios.txt",quote=F,row.names=F,sep="\t")
write.table(sel.r$table$genes,file="./data/all-genes.txt",quote=F,row.names=FALSE,col.names=FALSE,sep="\t")

# Filtering 

max.FDR <- 0.01
min.FC  <- 1.5
max.readCount <- 10

sel.r <- sel.r[sel.r$table$FDR<=max.FDR & abs(sel.r$table[[2]])>=min.FC & apply(dge$counts[rownames(sel.r),],1,sum)>=max.readCount,]

sel.rows<- rownames(sel.r)
sel     <- cbind(sel.r$table,dge$counts[sel.rows,])
colnames(sel)
sel.twt <- sel[order(sel$logFC),]

write.table(sel.r$table,file="./data/filtered-gene-ratios.txt",quote=F,row.names=FALSE,sep="\t")
write.table(sel.r$table$genes,file="./data/filtered-genes.txt",quote=F,row.names=FALSE,col.names=FALSE,sep="\t")

