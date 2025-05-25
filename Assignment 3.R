# Assignment 3

library(tidyverse)
library(ggplot2)
library(ggExtra)
library(dplyr)
library(GGally)
library(parallel)
library(cluster)
library(factoextra)

###################################################
########### Exploratory Data Analysis #############

#------------(a) Data Description_________________#

data <- read.table("STA4026_Assignment_Clustering.txt", header = FALSE)

colnames(data) <- c("V1", "V2") # we can decide of the variables later

og_size <- dim(data) # 2 variables with 5000 observations

data_types <-  sapply(data, class) # data type - both integers

missing_vals <-  colSums(is.na(data)) # 0 missing values in both columns

duplicates <- any(duplicated(data)) # false - no duplicates

infinite_vals <- any(sapply(data, function(x) any(is.infinite(x)))) # false - no infinite values

data_summary <- summary(data)

#----------(c) Exploratory Analysis---------------#

pairPlot <- ggplot(data, aes(x=V1, y=V2)) +
  geom_point(alpha = 0.5) +
  labs(x = expression(V[1]), y = expression(V[2])) + 
  theme_minimal()
ggMarginal(pairPlot, type = "densigram")


#-----(d) Exploratory Analysis of Distance--------#

dist <- dist(data, method = "euclidean")
dist_vals <- as.vector(dist)

ggplot(data.frame(Distance = dist_vals), aes(x = Distance)) +
  geom_histogram(aes(y = after_stat(count)),
    bins = 30, fill = "lightskyblue", color = "black") +
  geom_density(aes(y = after_stat(count) * (max(dist_vals) - min(dist_vals))/30),
    alpha = 0.2, color = "red") +
  labs(x = "Distance",
       y = "Frequency") +
  theme_minimal()

#----------(e) Outlier Identification-------------#

data_clean <- data |>
  mutate(dist_from_median = sqrt((V1 - median(V1))^2 + (V2 - median(V2))^2)) |>  mutate(observation = row_number())

outliers <-  data_clean |>
  arrange(desc(dist_from_median)) |>
  head(10)

data_cleaned <- anti_join(data, outliers)

ggplot(data, aes(V1, V2)) +
  geom_point() +
  geom_point(data = outliers, color = "red", size = 2) 

#----------(f) Correlation Analysis---------------#
cor_mat <- cor(data_cleaned$V1, data_cleaned$V2)


###################################################
############ Hyper-parameter Tuning ###############

#--------------(a) Selecting K--------------------#

# Range of K values
K <- 2:20

# compute distance matrix once
dists <- dist(data_cleaned)

### ---- K-MEANS ---- ###
sil_widths_kmeans <- numeric(length(K))

for (i in seq_along(K)) {
  k <- K[i]
  km <- kmeans(data_cleaned, centers = k, nstart =100)
  sil <- silhouette(km$cluster, dists)
  sil_widths_kmeans[i] <- mean(sil[, 3])
}

# Find the best K (max silhouette)
best_k_kmeans <- K[which.max(sil_widths_kmeans)]

# Plot for k-means
plot(K, sil_widths_kmeans, type = 'b', pch = 19, col = 'blue',
     xlab = "Number of Clusters K", ylab = "Average Silhouette Score",
     main = "")
abline(v = best_k_kmeans, lty = 2, col = "darkblue")

### ---- K-MEDOIDS (CLARA) ---- ###
sil_widths_medoid <- numeric(length(K))

for (i in seq_along(K)) {
  k <- K[i]
  km <- clara(data_cleaned, k = k, metric = "euclidean", pamLike = TRUE, samples = 50)
  sil <- silhouette(km$clustering, dists)
  sil_widths_medoid[i] <- mean(sil[, 3])
}

# Best K for medoids
best_k_medoids <- K[which.max(sil_widths_medoid)]

# Plot for k-medoids
plot(K, sil_widths_medoid, type = 'b', pch = 19, col = 'red',
     xlab = "Number of Clusters K", ylab = "Average Silhouette Score",
     main = "")
abline(v = best_k_medoids, lty = 2, col = "darkred")

#----------(b) Initialisation Sensitivity---------#

set.seed(2025)  # reproducibility

iter <- 50 # number of iterations

# Define k-means optimal K function
opt_kmeans <- function(data, k_range) {
  sil_width <- sapply(k_range, function(k) {
    km <- kmeans(data, centers = k, nstart = 50, iter.max = 50)
    sil <- silhouette(km$cluster, dist(data))
    mean(sil[, "sil_width"])
  })
  k_range[which.max(sil_width)]
}

# Define k-medoids optimal K function
opt_kmed <- function(data, k_range) {
  sil_width <- sapply(k_range, function(k) {
    cl <- clara(data, k, metric = "euclidean", samples = 50, pamLike = TRUE)
    cl$silinfo$avg.width
  })
  k_range[which.max(sil_width)]
}

# Create cluster for parallel processing
cl <- makeCluster(detectCores() - 1)

# Export variables and functions to the cluster
clusterExport(cl, c("opt_kmeans", "opt_kmed", "data_cleaned", "k_vals", "silhouette", "dist", "kmeans", "clara"))

# Parallel runs for k-means
opt_kmeans_pll <- parSapply(cl, 1:iter, function(i) opt_kmeans(data_cleaned, K))

# Parallel runs for k-medoids
opt_kmed_pll <- parSapply(cl, 1:iter, function(i) opt_kmed(data_cleaned, K))

stopCluster(cl)  # Stop cluster

# Frequency tables
kmeans_freq <- table(opt_kmeans_pll)
kmedoids_freq <- table(opt_kmed_pll)

# Plot k-means results #Optimal Number of Clusters (K*) - k-means (Parallel 100 runs)
HPT_kmeans <- ggplot(data.frame(K = names(kmeans_freq), Frequency = as.integer(kmeans_freq)), aes(x = K, y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "",
       x = "Optimal number of clusters (K*)",
       y = "Frequency") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Plot k-medoids results Optimal Number of Clusters (K*) - k-medoids (Parallel 100 runs)
HPT_kmeds <- ggplot(data.frame(K = names(kmedoids_freq), Frequency = as.integer(kmedoids_freq)), aes(x = K, y = Frequency)) +
  geom_bar(stat = "identity", fill = "#CC0000") +
  labs(title = "",
       x = "Optimal number of clusters (K*)",
       y = "Frequency") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Save plots (optional)
ggsave("HPT_kmeans.png", HPT_kmeans)
ggsave("HPT_kmeds.png", HPT_kmeds)

#---------(c) Increasing Initialisations----------#

cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("opt_kmeans", "data_cleaned", "k_vals", "silhouette", "dist", "kmeans"))
opt_kmeans_pll <- parSapply(cl, 1:200, function(i) opt_kmeans(data_cleaned, K))
stopCluster(cl)

parallel_k_freq <- table(opt_kmeans_pll)

#Distribution of Optimal K from 200 Parallel Initialisations
Kmean_200<-ggplot(data.frame(K = names(parallel_k_freq), Frequency = as.integer(parallel_k_freq)), aes(x = K, y = Frequency)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  labs(title = "",
       x = "Optimal number of clusters (K*)",
       y = "Frequency") +
  theme_minimal()

# Create cluster with number of cores minus one
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("opt_kmed", "data_cleaned", "k_vals", "clara"))
opt_kmed_pll <- parSapply(cl, 1:200, function(i) {
  opt_kmed(data_cleaned, K)
})

# Stop the cluster after finishing
stopCluster(cl)

# Distribution of Optimal K from 300 Parallel Initialisations (k-medoids)
opt_kmed_pll_freq <- table(opt_kmed_pll)

Kmed_200<-ggplot(data.frame(K = names(opt_kmed_pll_freq), Frequency = as.integer(opt_kmed_pll_freq)), aes(x = K, y = Frequency)) +
  geom_bar(stat = "identity", fill = "#CC0000") +
  labs(title = "",
       x = "Optimal number of clusters (K*)",
       y = "Frequency") +
  theme_minimal()

ggsave("HPT_kmeans_300.png", Kmean_300)
ggsave("HPT_kmeds_300.png", Kmed_300)

#----------(d) Selecting K using Gap--------------#
# Computing the gap statistic
set.seed(2025)
km_gap_stat <- clusGap(data_cleaned, FUN = kmeans, K.max = 20, B = 50)

kmed_gap_stat <- clusGap(data_cleaned, FUN = clara, K.max = 20, B = 50)
# Find the index of max gap value
best_k_kmeans <- which.max(km_gap_stat$Tab[,"gap"])
best_k_medoids <- which.max(kmed_gap_stat$Tab[,"gap"])


km_gap<-plot(km_gap_stat, main="", xlab="K", col="steelblue", pch=16, ylim = c(0.29,0.6))
kmed_gap<-plot(kmed_gap_stat, main="", xlab="K", col="#CC0000", pch=16, ylim = c(0.29,0.6))


plot(km_gap_stat, main="", xlab="K", col="steelblue", pch=16, ylim=c(0.29,0.6))
abline(v = best_k_kmeans, lty=2, col="steelblue")

plot(kmed_gap_stat, main="", xlab="K", col="#CC0000", pch=16, ylim=c(0.29,0.6))
abline(v = best_k_medoids, lty=2, col="#CC0000")

ggsave("km_gap.jpeg", km_gap)
ggsave("kmed_gap.jpeg", kmed_gap)

best_k_kmeans <- km_gap_stat$Tab[which.max(km_gap_stat$Tab[,"gap"]), "k"]
best_k_medoids <- kmed_gap_stat$Tab[which.max(kmed_gap_stat$Tab[,"gap"]), "k"]

###################################################
################Cluster Analysis###################
#---------(a) Silhouette Score Analysis-----------#

#................ K MEANS 15......................#

set.seed(2025)
colours15 <- c("#3300CC","#008080","#FF0000", "#FF9900","#0000FF", "#9900FF",
               "#FFFF33", "#339900", "#FF99FF", "#0099CC", "#99FF00",
               "#003399","#FFCC00", "#FF0066","#00CCCC")

km15 <- kmeans(data_cleaned, centers = 15, nstart = 100, iter.max = 50) # k-means with 15 cluster

data_cleaned$Cluster <- as.factor(km15$cluster)

# Compute silhouette widths
sil15 <- silhouette(km15$cluster, dist(data_cleaned))
data_cleaned$SilWidth = sil15[, 'sil_width']

plot(sil15, col = colours15, border = NA, main = "")

data_cleaned$Cluster <- as.factor(km15$cluster)

# Scatter plot of the clustered data
cluster_plot <- ggplot(data_cleaned, aes(x = V1, y = V2, color = Cluster)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "",
       x = "V1",
       y = "V2") +
  scale_color_manual(values = colours15) # Applying custom colors

# Print the plot
print(cluster_plot)

#................ K MEANS 16.....................#

set.seed(2025)
colours16 <- c("#008080","#FF0000", "#FF9900","#0000FF", "#9900FF",
                 "#FFFF33", "#339900", "#FF99FF", "#0099CC", "#99FF00",
                 "#003399","#FFCC00", "#FF0066","#00CCCC","#3300CC","#999999")

km16 <- kmeans(data_cleaned, centers = 16, nstart = 100, iter.max = 50) # k-means with 16 cluster

data_cleaned$Cluster16 <- as.factor(km16$cluster)

# Compute silhouette widths
sil16 <- silhouette(km16$cluster, dist(data_cleaned))

plot(sil16, col = colours16, border = NA, main = "")

data_cleaned$Cluster16 <- as.factor(km16$cluster)

# Scatter plot of the clustered data
cluster_plot <- ggplot(data_cleaned, aes(x = V1, y = V2, color = Cluster16)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "",
       x = "V1",
       y = "V2") +
  scale_color_manual(values = colours16) # Applying custom colors

# Print the plot
print(cluster_plot)

#................K MEDOIDS 14.....................#

set.seed(2025)
colours14 <- c("#3300CC","#008080","#FF0000", "#FF9900","#0000FF", "#9900FF",
               "#FFFF33", "#339900", "#FF99FF", "#0099CC", "#99FF00",
               "#003399","#FFCC00", "#FF0066")

cl14 <- clara(data_cleaned, 14, metric= "euclidean", samples = 50, pamLike=T)

sil14c <- silhouette(cl14$clustering, dist(data_cleaned))

plot(sil14c, col = colours14, border = NA, main = "")

data_cleaned$Cluster14c <- as.factor(cl14$clustering)

# Scatter plot of the clustered data_cleaned
cluster_plot <- ggplot(data_cleaned, aes(x = V1, y = V2, color = Cluster14c)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "",
       x = "V1",
       y = "V2",
       color="Cluster") +
  scale_color_manual(values = colours14) # Applying custom colors

# Print the plot
print(cluster_plot)

#................K MEDOIDS 20.....................#

set.seed(2025)
colours20 <- c(
  "#008080", "#FF0000", "#FF9900", "#0000FF", "#9900FF",
  "#FFFF33", "#339900", "#FF99FF", "#0099CC", "#99FF00",
  "#003399", "#FFCC00", "#FF0066", "#00CCCC", "#CC0066", "#3300CC",
  "#666699", "#FF6600", "#66FF66", "#993300"
)


cl20 <- clara(data_cleaned, 20, metric= "euclidean", samples = 50, pamLike=T)

sil20c <- silhouette(cl20$clustering, dist(data_cleaned))

plot(sil20c, col = colours20, border = NA, main = "")

data_cleaned$Cluster20c <- as.factor(cl20$clustering)

# Scatter plot of the clustered data
cluster_plot <- ggplot(data_cleaned, aes(x = V1, y = V2, color = Cluster20c)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "",
       x = "V1",
       y = "V2",
       color="Cluster") +
  scale_color_manual(values = colours20) # Applying custom colors

# Print the plot
print(cluster_plot)

#--------------(b) Outlier analysis---------------#
# locate outliers on cluster and silhouette plots
set.seed(2025)

km_with_outliers <- kmeans(data, centers = 15, nstart = 100, iter.max = 50) # k-means with 15 cluster on data with outliers

outlier_indices <- c(461, 340, 2540, 4722, 442, 486, 501, 2317, 535, 4723)

data$Cluster <- as.factor(km_with_outliers$cluster)

# Compute silhouette widths
sil_with_outliers <- silhouette(km_with_outliers$cluster, dist(data))
data$SilWidth = sil_with_outliers[, 'sil_width']

plot(sil_with_outliers, col = colours15, border = NA, main = "")

# Scatter plot of the clustered data
cluster_plot <- ggplot(data, aes(x = V1, y = V2, color = Cluster)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(125000, 875000), 
             color = "red", 
             linetype = "dotted",
             linewidth = 0.8) +
  theme_minimal() +
  labs(title = "",
       x = "V1",
       y = "V2") +
  scale_color_manual(values = colours15) +
  geom_point(data[outlier_indices, ],
             mapping = aes(x = V1, y=V2, col = "black"),
             col ="black", fill = "black", size = 2)

# Print the plot
print(cluster_plot)

##### compare original cluster with outliers to clustering without outliers
common_cols <- intersect(colnames(data), colnames(data_cleaned))
data_subset <- data[, common_cols]
data_cleaned_subset <- data_cleaned[, common_cols]

# original cluster 
set.seed(2025)
km_with_outliers <- kmeans(data_subset, centers = 15, nstart = 100, iter.max = 50)

# cluster without outliers
km_wo_outliers <- kmeans(data_cleaned_subset, centers = 15, nstart = 100, iter.max = 50)

# compare cluster centers
all_centers <- rbind(
  data.frame(km_with_outliers$centers, Type = "With Outliers"),
  data.frame(km_wo_outliers$centers, Type = "Without Outliers")
)

ggplot() +
  geom_point(data = data, aes(x = V1, y = V2), color = "gray20", alpha = 0.3) +
  geom_point(data = all_centers, aes(x = V1, y = V2, color = Type), size = 4) +
  scale_color_manual(values = c("With Outliers" = "green", "Without Outliers" = "purple")) +
  labs(x = "V1", y = "V2") + 
  geom_point(data[outlier_indices, ],
             mapping = aes(x = V1, y=V2, col = "black"),
             col ="black", fill = "black", size = 2) +
  theme_minimal()

#--------------(c) Post-Processing----------------#

km15 <- kmeans(data_cleaned, 15, nstart = 100, iter.max = 50)

# silhouette scores
sil15 <- silhouette(km15$cluster, dist(data_cleaned))
negative_indices <- which(sil15[, 3] < 0) # negative silhouette points

# reassign neg score obs, to closest different cluster (based on centroid diff)
centroids <- km15$centers
updated_clusters <- km15$cluster

for (i in negative_indices) {
  curr_cluster <- km15$cluster[i]
  dist <- sqrt(rowSums((t(t(centroids) - unlist(data_cleaned[i, ]))^2))) # calc dist to all centroids
  alt_clusters <- order(dist) # get dist to other clusters (asc - closest to furthest)
  alt_clusters <- alt_clusters[alt_clusters != curr_cluster] # ignore distance to its current cluster
  
  # only reassing if a better cluster exists
  if (length(alt_clusters) > 0) {
    updated_clusters[i] <- alt_clusters[1] # assign to closest diff cluster
  }
}

# creating a new data frame with updated clusters
new_df <- cbind(data_cleaned, cluster = as.factor(updated_clusters))

# plot the clusters with the re-assigned obs in black
ggplot(new_df, aes(x = V1, y = V2, color = cluster)) +
  geom_point(alpha = 0.7) +
  geom_point(data = new_df[negative_indices, ], 
             aes(x = V1, y = V2), 
             color = "black", size = 3) +
  theme_minimal() +
  scale_color_manual(values = colours15) +
  labs(x = "V1", y = "V2", color = "Clusters (Updated)")

# silhouette comparisons
# Original silhouette 
original_sil <- sil15
plot(original_sil, col = colours15, border = NA, main = "Original Silhouette")

# Updated silhouette after reassignment:
updated_sil <- silhouette(updated_clusters, dist(data_cleaned))
plot(updated_sil, col = colours15, border = NA, main = "")

# Compare average silhouette widths
cat("Average silhouette width before:", mean(original_sil[, 3]), "\n")
cat("Average silhouette width after: ", mean(updated_sil[, 3]), "\n")

# summary table of negative-score points
results_table <- data.frame(
  Index = negative_indices,
  Original_Silhouette_Score = round(sil15[negative_indices, 3], 4),
  Original_Cluster = sil15[negative_indices, 1],
  Reassigned_Cluster = sil15[negative_indices, 2],
  New_Silhouette_Score = round(updated_sil[negative_indices, 3], 4)
)






