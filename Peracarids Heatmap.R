rm(list=ls(all=TRUE))
setwd("~/GitHub/Peracarids")
Abu <- read.csv("~/GitHub/Peracarids/DataHM.csv")

# Load Libraries
library(gplots)  # for heatmap.2
library (vegan)  # for hierachical clustering
library(RColorBrewer)

# Obtain matrix
row.names(Abu) <- Abu$Site
Abu <- Abu[, -1]
Abu.prop <- Abu/rowSums(Abu)

# Determine the maximum relative abundance for each column
maxab <- apply(Abu.prop, 2, max)
head(maxab)

# Remove the genera with less than 5% as their maximum relative abundance
n1 <- names(which(maxab < 0.05))
Abu.prop.1 <- Abu.prop[, -which(names(Abu.prop) %in% n1)]

# Calculate the Bray-Curtis dissimilarity matrix on the full dataset:
data.dist <- vegdist(Abu, method = "bray")

# Do average linkage hierarchical clustering. 
row.clus <- hclust(data.dist, "aver") 
plot(row.clus)


# make the heatmap 
Abu.prop.1[Abu.prop.1==0] = NA # for cero cells to appear in white

heatmap.2(as.matrix(Abu.prop.1), Rowv = as.dendrogram(row.clus), Colv=FALSE,
          col = colorRampPalette(c("yellow", "red")), margins = c(11, 6),
          trace = "none", density.info = "none", xlab = "Taxa", ylab = "Site+Depth", 
          lhei = c(2, 8)) 