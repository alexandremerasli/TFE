## 

library("smallvis")
library("dplyr")


##

setwd("/Users/alexandre/Documents/TFE/Workspace/Synchrony_data/")
distance_matrix <- read.csv('dist_EEG.csv',sep=";",header = FALSE)
distance_matrix <- read.csv('dist_EEG_umap.csv',sep=";",header = FALSE)
distance_matrix <- read.csv('dist_EEG_test.csv',sep=";",header = FALSE)
condition <- read.csv('condition.csv',sep=";",header = FALSE)
conditionStr <- read.csv('condition.csv',sep=";",header = FALSE)
distance_matrix <- matrix(unlist(distance_matrix, use.names=FALSE),nrow=26,ncol=26)
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

points <- smallvis(distance_matrix, method = "umap", perplexity = 3, eta = 0.01, epoch_callback = umap_plot)
