#################################################
#These are stuffs from the first session with me
#################################################

###############
# Exercice 1 :
###############

# keep gene names
rownames(selection) <- selection[,1]
# remove this column because now it is in rownames if you want to keep gene names
selection[,1] <- NULL

mat <- data.matrix(selection[,13:44])

# Check the type of object you are playing with
#class(mat)
#typeof(mat)

# Write to file
write.table(mat,file="mat.txt", col.names=NA,row.names = TRUE, sep="\t", quote=FALSE)

# Webtools from the broad to explain hierachical clustering and kmeans.
# the different parameters (linkage, correlation or distance)
#https://software.broadinstitute.org/morpheus/

###############
# Exercice 2 : 
###############

# Complete correction with the use of cuttree...
mat <- data.matrix(selection[, 14:45])
nmat <- t(scale(t(mat)))

d.samples <- dist(t(nmat))
h.samples <- hclust(d.samples, method="ward.D")

d.genes <- dist(nmat) 
h.genes <- hclust(d.genes, method="ward.D")

rbcols <- redblue(100)

clusters.gene <- cutree(h.genes,4)
# should have explain more these two lines
cols.gene <- rainbow(4)
sample.rows <- cols.gene[clusters.gene] 

clusters <- cutree(h.samples,7)
cols <- rainbow(7)
sample.cols <- cols[clusters]

#  ColSideColors RowSideColors to add color label for columns and rows.
heatmap(nmat,scale="none",col=rbcols,Colv=d.samples,ColSideColors=sample.cols,RowSideColors=sample.rows)

# Library ComplexHeatmap based on ggplot, you get nicer heatmap with a lot of functionalities
# https://jokergoo.github.io/ComplexHeatmap-reference/book/index.htmls

###############
# Exercice 3  :
###############
# To get the list of object you can retrieve from the list returned by applying prcomp or tsne 
# ex: scale values with "scale" for pca or "x" for component
names(pca) 
names(rt) 
write.table(pca$x,file="mat.txt", col.names=NA,row.names = TRUE, sep="\t", quote=FALSE)

# how to install pacakge
install.packages("gplots")
