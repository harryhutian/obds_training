library(ggplot2)
library(tidyverse)
library(cowplot)
library(umap)

# Perform PCA. How many principal components do you think you should keep for follow up analysis?

log_counts <- read.csv("data/logcounts.csv", row.names = 1)

cell_metadata <- read.csv("data/cell_metadata.csv", row.names = 1)

# convert data.frame into matrix
log_counts_matrix <- as.matrix(log_counts)
# check dimension of matrix before PCA, know what to expect
dim(log_counts_matrix)

# perform PCA
log_counts_pca <- prcomp(t(log_counts_matrix),
                        center = TRUE, scale. = TRUE)

plot(log_counts_pca)
View(log_counts_pca$x)

# plot pca

ggplot(as.data.frame(log_counts_pca$x), aes(x= PC1, y=PC2))+
    geom_point()
screeplot(log_counts_pca, npcs = 20, type = 'lines')


# from screeplot you might take botton of elbow
# plot % variance
log_counts_pca$sdev

# sdev is sqr root of varicance so we square it

scree_pca_sdv <- data.frame(
    var = log_counts_pca$sdev**2,
    pc = seq_along(log_counts_pca$sdev)) %>%
    mutate(percent_var = var/sum(var)*100,
           cumulative_var = cumsum(percent_var))

# visualize the cummulative variance
ggplot(scree_pca_sdv[1:50,], aes(x= pc, y=cumulative_var))+
    geom_point()

# visualize the percent variance
ggplot(scree_pca_sdv[1:50,], aes(x= pc, y=percent_var))+
    geom_point()

# combine the cell_metadata with log_counts, and plot it
pca_df <- as.data.frame(log_counts_pca$x)
pca_df <- rownames_to_column(pca_df)
cell_metadata <- rownames_to_column(cell_metadata)

pca_df <- pca_df %>%
    left_join(cell_metadata)

# Use the experimental metadata to visualise which cell types visually cluster together in principal component space.
ggplot(pca_df, aes(x= PC1, y=PC2, colour = Status))+
    geom_point()

ggplot(pca_df, aes(x= PC1, y=PC2, colour = Time))+
    geom_point()

ggplot(pca_df, aes(x= PC1, y=PC2, colour = Infection))+
    geom_point()

ggplot(pca_df) +
    geom_density(aes(PC1, fill = Status), color = "black",
                 alpha = 0.5)+
    facet_grid(Time~Infection)+
    theme_cowplot()

str(log_counts_pca)
View(log_counts_pca$rotation)

topgenes_PC1 <- log_counts_pca$rotation[,1]
View(topgenes_PC1)
top10_genes_PC1 <- sort(topgenes_PC1, decreasing = TRUE)
View(top10_genes_PC1)

gene_metadata <- read.csv("data/gene_metadata.csv",
                          row.names = 1)

# clustering the data
kmeans <- kmeans(t(log_counts_matrix), centers = 4)

# subsetting, and adding the new column
pca_df$cluster <- as.factor(kmeans$cluster[pca_df$rowname])


ggplot(pca_df, aes(x= PC1, y=PC2,
                   colour = cluster, shape = Time))+
    geom_point()
#For good clustering we want small sum(withinss) and large betweenss , so this ratio we want to be as large as possible.

kmeans$withinss
kmeans$betweenss

candidate_k = 2:10

km <- sapply(candidate_k, function(i){
    km <- kmeans(t(log_counts_matrix),
           centers = i)
    sum(km$withinss)}
    )

kmeans_sums <- data.frame(sum_withinss = km, k = candidate_k)
ggplot(kmeans_sums, aes(x= k, y= sum_withinss))+
    geom_point()+


gg_cluster <- ggplot(
                pca_df, aes(x= PC1, y=PC2,
                colour = cluster))+
    geom_point()

gg_time <- ggplot(
    pca_df, aes(x= PC1, y=PC2,
                colour = Time))+
    geom_point()

gg_infection <- ggplot(
    pca_df, aes(x= PC1, y=PC2,
                colour = Infection))+
    geom_point()


plot_grid(gg_cluster,gg_time, gg_infection)

# umap
umap_dc <- umap(log_counts_pca$x)

umap_coord <- as.data.frame(umap_dc$layout)

umap_coord <- cbind(umap_coord,cell_metadata)

ggplot(umap_coord,aes(x = V1, y = V2, color = Time))+
    geom_point()+
    theme_cowplot()
