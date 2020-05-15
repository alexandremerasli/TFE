## 

library("dplyr")


##

setwd("/Users/alexandre/Documents/TFE/Workspace/Synchrony_data/")
condition <- read.csv('condition.csv',sep=";",header = FALSE)
conditionStr <- read.csv('condition.csv',sep=";",header = FALSE)
condition <- unlist(condition, use.names=FALSE)
conditionStr <- unlist(conditionStr, use.names=FALSE)
N <- length(condition)
subjects <- 1:N

# Using a custom epoch_callback
for (i in subjects) {
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
  par(pty="s")
  plot(x, col = rgb(condition==1,0,condition==0),pch=16,asp = 1)
  text(x[c(1:26),1],x[c(1:26),2],subjects-1,pos = 4,cex = 0.7)
}

library(multiview)

## MVMDS metric euclidean
distance_matrix_EEG <- read.csv('dist_EEG_basic.csv',sep=";",header = FALSE)
distance_matrix_EEG <- matrix(unlist(distance_matrix_EEG,use.names=FALSE),nrow=26,ncol=26)

distance_matrix_EDA<- read.csv('dist_EDA_basic.csv',sep=";",header = FALSE)
distance_matrix_EDA <- matrix(unlist(distance_matrix_EDA,use.names=FALSE),nrow=26,ncol=26)

distance_matrix_IBI <- read.csv('dist_IBI_basic.csv',sep=";",header = FALSE)
distance_matrix_IBI <- matrix(unlist(distance_matrix_IBI,use.names=FALSE),nrow=26,ncol=26)

#set.seed(42) 
points <- mvmds(list(dist(distance_matrix_EEG),dist(distance_matrix_EDA),dist(distance_matrix_IBI)), k = 2)
umap_plot(points)

## MVMDS directly on distance matrix
distance_matrix_EEG <- read.csv('dist_EEG_coordinates.csv',sep=";",header = FALSE)
distance_matrix_EEG <- matrix(unlist(distance_matrix_EEG,use.names=FALSE),nrow=26,ncol=26)

distance_matrix_EDA<- read.csv('dist_EDA_coordinates.csv',sep=";",header = FALSE)
distance_matrix_EDA <- matrix(unlist(distance_matrix_EDA,use.names=FALSE),nrow=26,ncol=26)

distance_matrix_IBI <- read.csv('dist_IBI_coordinates.csv',sep=";",header = FALSE)
distance_matrix_IBI <- matrix(unlist(distance_matrix_IBI,use.names=FALSE),nrow=26,ncol=26)

#set.seed(42) 
points <- mvmds(list(dist(distance_matrix_EEG),dist(distance_matrix_EDA),dist(distance_matrix_IBI)), k = 2)
umap_plot(points)

library(multiview)

# ## MVTSNE metric euclidean
# distance_matrix_EEG <- read.csv('dist_EEG_basic.csv',sep=";",header = FALSE)
# distance_matrix_EEG <- matrix(unlist(distance_matrix_EEG,use.names=FALSE),nrow=26,ncol=26)
# 
# distance_matrix_EDA<- read.csv('dist_EDA_basic.csv',sep=";",header = FALSE)
# distance_matrix_EDA <- matrix(unlist(distance_matrix_EDA,use.names=FALSE),nrow=26,ncol=26)
# 
# distance_matrix_IBI <- read.csv('dist_IBI_basic.csv',sep=";",header = FALSE)
# distance_matrix_IBI <- matrix(unlist(distance_matrix_IBI,use.names=FALSE),nrow=26,ncol=26)
# 
# #set.seed(42) 
# points <- mvtsne(list(dist(distance_matrix_EEG),dist(distance_matrix_EDA),dist(distance_matrix_IBI)), k = 2,perplexity = 10)
# umap_plot(points$embedding)
# 
# ## MVTSNE directly on distance matrix
# distance_matrix_EEG <- read.csv('dist_EEG_coordinates.csv',sep=";",header = FALSE)
# distance_matrix_EEG <- matrix(unlist(distance_matrix_EEG,use.names=FALSE),nrow=26,ncol=26)
# 
# distance_matrix_EDA<- read.csv('dist_EDA_coordinates.csv',sep=";",header = FALSE)
# distance_matrix_EDA <- matrix(unlist(distance_matrix_EDA,use.names=FALSE),nrow=26,ncol=26)
# 
# distance_matrix_IBI <- read.csv('dist_IBI_coordinates.csv',sep=";",header = FALSE)
# distance_matrix_IBI <- matrix(unlist(distance_matrix_IBI,use.names=FALSE),nrow=26,ncol=26)
# 
# #set.seed(42)
# points <- mvtsne(list(dist(distance_matrix_EEG),dist(distance_matrix_EDA),dist(distance_matrix_IBI)), k = 2,perplexity = 5)
# umap_plot(points$embedding)
# 
# 
# ## MVSC metric euclidean
# distance_matrix_EEG <- read.csv('dist_EEG_basic.csv',sep=";",header = FALSE)
# distance_matrix_EEG <- matrix(unlist(distance_matrix_EEG,use.names=FALSE),nrow=26,ncol=26)
# 
# distance_matrix_EDA<- read.csv('dist_EDA_basic.csv',sep=";",header = FALSE)
# distance_matrix_EDA <- matrix(unlist(distance_matrix_EDA,use.names=FALSE),nrow=26,ncol=26)
# 
# distance_matrix_IBI <- read.csv('dist_IBI_basic.csv',sep=";",header = FALSE)
# distance_matrix_IBI <- matrix(unlist(distance_matrix_IBI,use.names=FALSE),nrow=26,ncol=26)
# 
# #set.seed(42) 
# points <- mvsc(list(dist(distance_matrix_EEG),dist(distance_matrix_EDA),dist(distance_matrix_IBI)), k = 2)
# umap_plot(points$evectors)
# 
# ## MVTSNE directly on distance matrix
# distance_matrix_EEG <- read.csv('dist_EEG_coordinates.csv',sep=";",header = FALSE)
# distance_matrix_EEG <- matrix(unlist(distance_matrix_EEG,use.names=FALSE),nrow=26,ncol=26)
# 
# distance_matrix_EDA<- read.csv('dist_EDA_coordinates.csv',sep=";",header = FALSE)
# distance_matrix_EDA <- matrix(unlist(distance_matrix_EDA,use.names=FALSE),nrow=26,ncol=26)
# 
# distance_matrix_IBI <- read.csv('dist_IBI_coordinates.csv',sep=";",header = FALSE)
# distance_matrix_IBI <- matrix(unlist(distance_matrix_IBI,use.names=FALSE),nrow=26,ncol=26)
# 
# #set.seed(42)
# points <- mvsc(list(dist(distance_matrix_EEG),dist(distance_matrix_EDA),dist(distance_matrix_IBI)), k = 2)
# umap_plot(points$evectors)
