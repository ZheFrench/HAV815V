#' ---
#' title: "Go Parse David Output"
#' author: "Jean-Philippe Villemin"
#' date: "Fev 2, 2021"

#' ---

#' ## Description
#' Parse David ouput and plot barplot with -log10 PValcorrected.
#' 
#' ## Usage
#' **Rscript GoParseDavid.R **
#'
#' ## Notes
#' None.
#'
library(knitr)
library(reshape)
library(plyr)
library(ggplot2)
library(tidyr)
library(stringr)
require(gridExtra)

######################################################
#########         MAIN                     ###########
######################################################

file  <- "./data/DAVID.txt"
filename <- file
type  <- "GOTERM_BP_DIRECT"
pval  <- 0.0000001
color <- "red"
y     <- 1

cutoff<-pval

mydata <- read.csv(file=file, header=TRUE, sep="\t",stringsAsFactors=FALSE)

mydata1 <- data.frame(do.call('rbind', strsplit(as.character(mydata[,2]),'~',fixed=TRUE)))
clean <- cbind(mydata1, mydata  )

clean <- clean[clean$Benjamini <= cutoff,]
clean <- clean[clean$Category ==  type,]

colnames(clean)[6]
colnames(clean)[2] <-"GO"
colnames(clean)[6] <-"percent"
head(clean,5)

filename = gsub(".txt", "", filename)
final=paste0(c(filename, type,cutoff,"barplot.svg"),collapse=".")


clean$log10Benjamini <- -log10(clean$Benjamini)
clean_order <- clean[order(clean$Benjamini),]
clean_order$GO <- factor(clean_order$GO, levels = clean_order$GO[order(clean_order$log10Benjamini)])

head(clean_order)
nrow <- nrow(clean_order)

plot1 <- ggplot(clean_order, aes(x=GO, y=log10Benjamini)) +
 geom_bar(width=0.8,stat='identity', color=color,fill="#D7D7D7",linewidth=1)+
 geom_text(data=clean_order, aes(x = GO, y =  y , label =GO,hjust="left"),size=nrow*0.4) +
 geom_hline(yintercept = -log10(cutoff), color='black',linewidth=1) +
 coord_flip() +
 theme_bw(base_size=2)+
 theme(text = element_text(size = 9),
       axis.text = element_text(size = 12),
       axis.text.y=element_blank(),
       )+
 xlab(type)  


print(final)
# Modify width, and size of geom_text to create something clean depending of the number of terms.
svg(file=final,width=12,height=nrow*0.5)
plot1
dev.off()

