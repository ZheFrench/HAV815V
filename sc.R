library(Seurat)

library(glue)
library(ggplot2)

dir_path <- "/home/jp/Desktop/TD.JC/PBMC-10x/"

seurat.object <- CreateSeuratObject(
     counts = Read10X( data.dir = dir_path ,  gene.column = 2),
     min.cells = 3, min.features = 200
)

seurat.object

seurat.object          <-  PercentageFeatureSet(object = seurat.object, pattern = "^MT-", col.name = "mitoRatio")
seurat.object          <-  PercentageFeatureSet(object = seurat.object, pattern = "^RP[SL]",col.name = "riboRatio")


#print(seurat.object@assays$RNA@counts)
#print(seurat.object@meta.data)

write.table(seurat.object@meta.data ,file = "meta.tsv", col.names=TRUE,row.names =  FALSE, quote = FALSE,sep = "\t")
write.table(seurat.object@assays$RNA@counts ,file = "counts.tsv", col.names=TRUE,row.names =  TRUE, quote = FALSE,sep = "\t")

colnames(seurat.object@meta.data)
names(seurat.object@assays)


# Normalisation
seurat.object <- SCTransform(seurat.object,  verbose = FALSE)

seurat.object <- RunPCA(seurat.object, assay = "SCT", verbose = FALSE)
seurat.object <- FindNeighbors(seurat.object, reduction = "pca", dims = 1:30) #, k.param=40
seurat.object <- FindClusters(seurat.object, verbose = FALSE)
seurat.object <- RunUMAP(seurat.object, reduction = "pca", dims = 1:30)
seurat.object <- RunTSNE(seurat.object, dims=1:10)



png("nCount.png",width=600,height=300)
FeaturePlot(seurat.object, features = "nCount_RNA") + theme(legend.position = "right")
dev.off()

png("nFeature.png",width=600,height=300)
FeaturePlot(seurat.object, features = "nFeature_RNA") + theme(legend.position = "right") 
dev.off()


bench.dir <- "."

x <- c(0,1)
for (cluster in x) {
    genes.cl <- FindMarkers(seurat.object,ident.1 = cluster, min.pct = 0.25)
    write.table(genes.cl ,file = glue("{bench.dir}/{cluster}.tsv"), col.names=NA,row.names =  TRUE, quote = FALSE,sep = "\t")
}

p1 <- DimPlot(seurat.object, reduction = "umap", label = TRUE)
png("dimplot_umap.png",width=600,height=600)
p1
dev.off()

p1 <- DimPlot(seurat.object, reduction = "tsne", label = TRUE)
png("dimplot_tsnee.png",width=600,height=600)
p1
dev.off()


features <- c("CD19","CD27","CD40","CXCR5")
   
p3 <- VlnPlot(seurat.object, features = features )
png("VlnPlot.genes.png",width=900,height=900)
p3
dev.off()

p4 <-FeaturePlot(seurat.object, features = features)
png("FeaturePlot.genes.png",width=900,height=900)
p4
dev.off()

save(seurat.object,file="seurat.object.rda")
load("seurat.object.rda")

# head(seurat.object@assays$SCT@counts)
dataframe <-   as.data.frame(as.matrix(seurat.object@assays$SCT@counts))

