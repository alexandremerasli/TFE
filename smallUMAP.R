## 

library("smallvis")
library("dplyr")


##

setwd("/Users/alexandre/Documents/TFE/Workspace/Synchrony_data/")
condition <- read.csv('condition.csv',sep=";",header = FALSE)
conditionStr <- read.csv('condition.csv',sep=";",header = FALSE)
condition <- unlist(condition, use.names=FALSE)
conditionStr <- unlist(conditionStr, use.names=FALSE)
subjects <- 1:length(condition)

# Using a custom epoch_callback
for (i in subjects) {
  print(i)
  if (condition[i]==0) {
    conditionStr[i] <- "NA"  
  } else {
    conditionStr[i] <- "SSA" 
  }
}
uniq_spec <- unique(condition)
colors <- rainbow(length(uniq_spec))
names(colors) <- uniq_spec
umap_plot <- function(x) {
  plot(x, col = rgb(condition==1,0,condition==0),pch=16)
  text(x[c(1:26),1],x[c(1:26),2],subjects-1,pos = 4,cex = 0.7)
}

# ## Euclidean metric on initial distance matrix
#
# distance_matrix <- read.csv('dist_EEG.csv',sep=";",header = FALSE)
# distance_matrix <- matrix(unlist(distance_matrix, use.names=FALSE),nrow=26,ncol=26)
# set.seed(42)
# points <- smallvis(distance_matrix, method = "umap", perplexity = 3, eta = 0.01, epoch_callback = umap_plot,max_iter = 5000)
#
#
#
# ## Euclidean metric on normed distance matrix
#
# distance_matrix <- read.csv('dist_EEG_umap.csv',sep=";",header = FALSE)
# distance_matrix <- matrix(unlist(distance_matrix, use.names=FALSE),nrow=26,ncol=26)
# set.seed(42)
# points <- smallvis(distance_matrix, method = "umap", perplexity = 3, eta = 0.01, epoch_callback = umap_plot)


## Correlation metric on study matrix

library(ape)

study_matrix <- read.csv('ISC_EEG.csv',sep=",",header = FALSE)
study_matrix <- matrix(unlist(study_matrix, use.names=FALSE),nrow=26,ncol=26)
# #set.seed(42) 
correlation <- cor(study_matrix, use = "pairwise.complete.obs", method = "pearson")
correlation <- 1 - correlation
pcoa_points <- pcoa(correlation, correction="none", rn=NULL)$vectors
points <- smallvis(pcoa_points, method = "umap", perplexity = 3, eta = 0.01, epoch_callback = umap_plot)

# points <- smallvis(1-study_matrix/3, method = "umap", perplexity = 3, eta = 0.01, epoch_callback = umap_plot)

#distance_matrix <- read.csv('dist_EEG_test.csv',sep=";",header = FALSE)
#
#
## UMAP on initial distance matrix

# pcoa_points <- read.csv('dist_EEG_coordinates.csv',sep=";",header = FALSE)
# pcoa_points <- matrix(unlist(pcoa_points, use.names=FALSE),nrow=26,ncol=25)
# set.seed(42)
# distance <- dist(pcoa_points,method = "euclidean")
# print(distance)
# points <- smallvis(pcoa_points, method = "umap", perplexity = 13, eta = 0.01, epoch_callback = umap_plot)


# ## Euclidean distances on correlation metric on study matrix (small dataset)
# 
# study_matrix <- read.csv('ISC_EEG.csv',sep=",",header = FALSE)
# study_matrix <- matrix(unlist(study_matrix, use.names=FALSE),nrow=26,ncol=26)
# set.seed(42)
# 
# correlation <- cor(study_matrix, use = "pairwise.complete.obs", method = "pearson")
# correlation <- 1 - correlation
# points <- smallvis(correlation, method = "umap", perplexity = 3, eta = 0.01, epoch_callback = umap_plot)
# 

# ## Comparing with approximate UMAP (comparison not good for the moment)
# 
# library(umap)
# distance_matrix <- read.csv('dist_EEG_basic.csv',sep=";",header = FALSE)
# distance_matrix <- matrix(unlist(distance_matrix, use.names=FALSE),nrow=26,ncol=26)
# study_matrix <- read.csv('ISC_EEG.csv',sep=",",header = FALSE)
# study_matrix <- matrix(unlist(study_matrix, use.names=FALSE),nrow=26,ncol=26)
# set.seed(42)
# 
# correlation <- cor(study_matrix, use = "pairwise.complete.obs", method = "pearson")
# correlation <- 1 - correlation
# custom.config = umap.defaults
# custom.config$n_neighbors = 3
# custom.config$n_epochs = 5000
# custom.config$metric = "euclidean"
# custom.config$min_dist = 0.00001
# 
# pointsUMAP <- umap(distance_matrix,custom.config)
# 
# umap_plot(pointsUMAP$layout)
