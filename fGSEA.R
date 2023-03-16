#################################################################
#
# date: March 16, 2023
# platform: Ubuntu 18.04
# R.version : 4.2.2
# author: Villemin Jean-Philippe
# team: Bioinformatique et biologie des systèmes du cancer : J. Colinge 
# Institute : IRCM
#
# fGSEA.R
# Usage : 
# 
# fGSEA.R
# 
# Description : 
#
# GSEA can be done with clusterProfiler or fgsea package.
# Both offer complementary functions.
#
##############################################################################

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
#suppressPackageStartupMessages(library(reshape2))

suppressPackageStartupMessages(library(ggplot2))

suppressPackageStartupMessages(library(glue))

suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(fgsea))

library(viridis)

# USEFUL LINKS
#https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/
#https://github.com/ctlab/fgsea
#https://stephenturner.github.io/deseq-to-fgsea/

# ==========   Settings =============================

base.dir <- "./data//"

set.seed(423)
#------------------

#------------------------------------------------------------
#   Read directory with files previously created with logFC and apply fGSEA
#-------------------------------------------------------------------------
file.gmt  <- "./data/h.all.v7.2.symbols.gmt" # From MsiGdb
#file.gmt  <- "./data/c5.go.mf.v2023.1.Hs.symbols.gmt" # From MsiGdb also
# you can use different references.

asbolutepath2file <-  "./data/all-gene-ratios.txt" 

# Two ways to read gmt (two functions from fgsea or clusterProfiler)
h.All     <- gmtPathways(file.gmt) 
h.All.bis <- read.gmt(file.gmt)

# OUTPUT
final.file.padj <- glue("A-fgsea.padj.txt")
final.file.nes  <- glue("A-fgsea.nes.txt")

# fread (from data.table)
dataframe.expression <- fread(asbolutepath2file,data.table=F)

dataframe.expression <- subset(dataframe.expression,select=c(genes,logFC))

dataframe.expression <- dataframe.expression[order(dataframe.expression$logFC),]

# ClusterProfiler need a decreasing order...
dataframe.expression.decreasing <- dataframe.expression[order(dataframe.expression$logFC, decreasing = TRUE),]
ranks_decreasing <- deframe(dataframe.expression.decreasing)

ranks <- deframe(dataframe.expression)

# But genes used are the same so I used it to plot with heatplot function of clusterProfiler the leading edge genes contained in the signature I am interested in
egmt2 <- GSEA(ranks_decreasing, TERM2GENE = h.All.bis, by = "fgsea",verbose=TRUE ,nPermSimple = 10000 ,minGSSize  = 10, maxGSSize  = 325 , eps = 0,  pvalueCutoff = 1)
# You dont need the whole object to be written
write.table(egmt2, file=glue("./A-full-gsea-clusterprofiler.tsv"),quote=F,row.names=F,sep="\t")

# fgseaMultilevel is from fgsea package
# not using fgseaMultilevel...a really tiny difference in NES score due to the fact it just use fgsea

# Yeah I do it again (I know...)
# Because I use then plotEnrichment function
fgseaRes     <- fgseaMultilevel(pathways=h.All, stats=ranks,eps=0, nPermSimple = 10000 ,minSize  = 10, maxSize  = 325)


# Clearly there is something with EMT. ( but I would have expect a positive NES)
# Check the doc of the package and read about to sort the gene etc to correct...
for (pathway in names(h.All)){
  if (pathway %in% c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")){
    print(pathway)
    #print(h.All[[pathway]])
    p<- plotEnrichment(h.All[[pathway]], ranks) + labs(title=pathway)
    png(file=glue("{pathway}_plotEnrichment.png"))
    print(p)
    dev.off()
    
  }  
}
  
  fgseaResTidy <- fgseaRes %>% as_tibble() # it's like a dataframe

  fwrite( filter(fgseaResTidy ,padj <= 1), file=glue("A-full-fgsea.tsv"), sep="\t", sep2=c("", "/", ""))

  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]

  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]

  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
png(file=glue("A-top-global.png"),width=900)
plotGseaTable(h.All[topPathways], ranks, fgseaRes, gseaParam=0.5) 
dev.off()


#  You can also do a gene ontology enrichment analysis
# -----------------------------------------------------

file.gmt  <- "./data/c5.go.mf.v2023.1.Hs.symbols.gmt" # From MsiGdb
GO_ALL <-  read.gmt(file.gmt)  
GO_ALLterm2gene <- GO_ALL[, c("term", "gene")]

# Use your DE Gene List
genes   <- fread(glue("./data/filtered-genes.txt"),data.table=F)

# list to character
genes <- unlist(genes, use.names=FALSE)

dataframe <- as.data.frame(genes)
dataframe$geneslist  <- "genes.list1"

genes.enriched.GO_ALL <- enricher(genes, 
  TERM2GENE=GO_ALLterm2gene, TERM2NAME=NA,
  pvalueCutoff = 0.05,pAdjustMethod = "BH",
  minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.05)

# I let you adapt the size/scale of the plots
png(file=glue("./data/barplot.png"),width=900,height=1+dim(genes.enriched.GO_ALL)[[1]]*50)
barplot(genes.enriched.GO_ALL  ,  drop = TRUE,  showCategory = 50,  title = "GO", font.size = 12) +scale_fill_viridis(option="viridis",direction = -1)
dev.off()

png(file=glue("./data/heatplot.png"),width=900,height=1+dim(genes.enriched.GO_ALL)[[1]]*50)
heatplot(genes.enriched.GO_ALL)+ theme(  axis.text.y = element_text( size = 12))
dev.off()