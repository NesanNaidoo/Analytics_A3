# Assignment 3

library(tidyverse)
library(ggplot2)
library(ggExtra)
library(dplyr)

###################################################
########### Exploratory Data Analysis #############

#-------------Data Description___________________#

data <- read.table("STA4026_Assignment_Clustering.txt", header = FALSE)

colnames(data) <- c("V1", "V2") # we can decide of the variables later

og_size <- dim(data) # 2 variables with 5000 observations

data_types <-  sapply(data, class) # data type - both integers

missing_vals <-  colSums(is.na(data)) # 0 missing values in both columns

duplicates <- any(duplicated(data)) # false - no duplicates

infinite_vals <- any(sapply(data, function(x) any(is.infinite(x)))) # false - no infinite values

data_summary <- summary(data)

#------------Exploratory Analysis----------------#

pairPlot <- ggplot(data, aes(x=V1, y=V2)) +
  geom_point(alpha = 0.5) +
  labs(x = expression(V[1]), y = expression(V[2])) + 
  theme_minimal()
ggMarginal(pairPlot, type = "densigram")





#-------Exploratory Analysis of Distance--------#

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

#------------Outlier Identification------------#

data_clean <- data |>
  mutate(dist_from_median = sqrt((V1 - median(V1))^2 + (V2 - median(V2))^2)) |>  mutate(observation = row_number())

outliers <-  data_clean |>
  arrange(desc(dist_from_median)) |>
  head(10)

data_cleaned <- anti_join(data, outliers)

ggplot(data, aes(V1, V2)) +
  geom_point() +
  geom_point(data = outliers, color = "red", size = 2) 

#------------Correlation Analysis---------------#
cor_mat <- cor(data_cleaned$V1, data_cleaned$V2)




###Alternative EDA
library(GGally)
ggpairs(data_cleaned, title=NULL)+ #pair plot
  theme_bw() + # Start with a black-and-white theme
  theme(
    plot.background = element_blank(), # Remove background
    panel.background = element_blank(), # Remove panel background
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    strip.background = element_blank(), # Remove background from facet labels
    plot.margin = unit(c(1, 1, 1, 1), "lines") # Adjust plot margin
  )


#pairwise distances
distances = dist(data_cleaned, method = "euclidean", diag = FALSE, upper = FALSE)

hist(distances, breaks = 30, freq = FALSE, col = "pink" , main = "", ylim = c(0, 2.5e-06))
lines(density(distances), col = "blue", lwd = 2)


# Hyperparamter Tuning
##### (a) Selecting K#####

library(cluster)
library(factoextra)
library(ggplot2)

# Range of K values
K <- 2:20

# compute distance matrix once
dists <- dist(data_cleaned)

### ---- K-MEANS ---- ###
sil_widths_kmeans <- numeric(length(K))

for (i in seq_along(K)) {
  k <- K[i]
  km <- kmeans(data_cleaned, centers = k, nstart =100)
  ss <- silhouette(km$cluster, dists)
  sil_widths_kmeans[i] <- mean(ss[, 3])
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
  ss <- silhouette(km$clustering, dists)
  sil_widths_medoid[i] <- mean(ss[, 3])
}

# Best K for medoids
best_k_medoids <- K[which.max(sil_widths_medoid)]

# Plot for k-medoids
plot(K, sil_widths_medoid, type = 'b', pch = 19, col = 'red',
     xlab = "Number of Clusters K", ylab = "Average Silhouette Score",
     main = "")
abline(v = best_k_medoids, lty = 2, col = "darkred")





#### (b) Initialisation Sensitivity #####

library(parallel)
library(ggplot2)

set.seed(2025)  # reproducibility

k_vals <- 2:20
num_iter <- 50

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
opt_kmeans_pll <- parSapply(cl, 1:num_iter, function(i) opt_kmeans(data_cleaned, k_vals))

# Parallel runs for k-medoids
opt_kmed_pll <- parSapply(cl, 1:num_iter, function(i) opt_kmed(data_cleaned, k_vals))

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



#### (c) Increasing Initialisations ####

library(parallel)

cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("opt_kmeans", "data_cleaned", "k_vals", "silhouette", "dist", "kmeans"))
opt_kmeans_pll <- parSapply(cl, 1:200, function(i) opt_kmeans(data_cleaned, k_vals))
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
  opt_kmed(data_cleaned, k_vals)
})

# Stop the cluster after finishing
stopCluster(cl)

# Distribution of Optimal K from 300 Parallel Initialisations (k-medoids)
opt_kmed_pll_freq <- table(opt_kmed_pll)

library(ggplot2)
Kmed_200<-ggplot(data.frame(K = names(opt_kmed_pll_freq), Frequency = as.integer(opt_kmed_pll_freq)), aes(x = K, y = Frequency)) +
  geom_bar(stat = "identity", fill = "#CC0000") +
  labs(title = "",
       x = "Optimal number of clusters (K*)",
       y = "Frequency") +
  theme_minimal()

ggsave("HPT_kmeans_300.png", Kmean_300)
ggsave("HPT_kmeds_300.png", Kmed_300)



#### Selecting K using Gap ####
# Computing the gap statistic
set.seed(2025)
gap_stat_km <- clusGap(data_cleaned, FUN = kmeans, K.max = 20, B = 50)

gap_stat_kmed <- clusGap(data_cleaned, FUN = clara, K.max = 20, B = 50)
# Find the index of max gap value
best_k_kmeans <- which.max(gap_stat_km$Tab[,"gap"])
best_k_medoids <- which.max(gap_stat_kmed$Tab[,"gap"])


km_gap<-plot(gap_stat_km, main="", xlab="K", col="steelblue", pch=16, ylim = c(0.29,0.6))
kmed_gap<-plot(gap_stat_kmed, main="", xlab="K", col="#CC0000", pch=16, ylim = c(0.29,0.6))


plot(gap_stat_km, main="", xlab="K", col="steelblue", pch=16, ylim=c(0.29,0.6))
abline(v = best_k_kmeans, lty=2, col="steelblue")

plot(gap_stat_kmed, main="", xlab="K", col="#CC0000", pch=16, ylim=c(0.29,0.6))
abline(v = best_k_medoids, lty=2, col="#CC0000")





ggsave("km_gap.jpeg", km_gap)
ggsave("kmed_gap.jpeg", kmed_gap)

best_k_kmeans <- gap_stat_km$Tab[which.max(gap_stat_km$Tab[,"gap"]), "k"]
best_k_medoids <- gap_stat_kmed$Tab[which.max(gap_stat_kmed$Tab[,"gap"]), "k"]







#Recommendations for Q3
# Use k =15,K=16  for k means
# Use K= 20, k=14 for k-mediods

##############################################################################
#### QUESTION 3 ####
library(dplyr)
set.seed(10)

#### K MEANS 15
my_colors <- c("#3300CC","#008080","#FF0000", "#FF9900","#0000FF", "#9900FF",
               "#FFFF33", "#339900", "#FF99FF", "#0099CC", "#99FF00",
               "#003399","#FFCC00", "#FF0066","#00CCCC")

km <- kmeans(data, centers = 15, nstart = 50) # k-means with 15 cluster

data$Cluster <- as.factor(km$cluster)

# Compute silhouette widths
sil <- silhouette(km$cluster, dist(data))
data$SilWidth = sil[, 'sil_width']

plot(sil, col = my_colors, border = NA, main = "")

data$Cluster <- as.factor(km$cluster)

# Scatter plot of the clustered data
cluster_plot <- ggplot(data, aes(x = X1, y = X2, color = Cluster)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "",
       x = "X1",
       y = "X2") +
  scale_color_manual(values = my_colors) # Applying custom colors

# Print the plot
print(cluster_plot)

####### K MEANS 14#########

set.seed(10)
my_colors14 <- c("#008080","#FF0000", "#FF9900","#0000FF", "#9900FF",
                 "#FFFF33", "#339900", "#FF99FF", "#0099CC", "#99FF00",
                 "#003399","#FFCC00", "#FF0066","#00CCCC")

km14 <- kmeans(data, centers = 14, nstart = 50) # k-means with 14 cluster

data$Cluster14 <- as.factor(km14$cluster)

# Compute silhouette widths
sil14 <- silhouette(km14$cluster, dist(data))

plot(sil14, col = my_colors14, border = NA, main = "")

data$Cluster14 <- as.factor(km14$cluster)

# Scatter plot of the clustered data
cluster_plot <- ggplot(data, aes(x = X1, y = X2, color = Cluster14)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "",
       x = "X1",
       y = "X2") +
  scale_color_manual(values = my_colors14) # Applying custom colors

# Print the plot
print(cluster_plot)

###### K MEDOIDS 14 ######

set.seed(10)

cl14 <- clara(data, 14, metric= "euclidean", samples = 50, pamLike=T)

sil14c <- silhouette(cl14$clustering, dist(data))

plot(sil14c, col = my_colors14, border = NA, main = "")

data$Cluster14c <- as.factor(cl14$clustering)

# Scatter plot of the clustered data
cluster_plot <- ggplot(data, aes(x = X1, y = X2, color = Cluster14c)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "",
       x = "X1",
       y = "X2",
       color="Cluster") +
  scale_color_manual(values = my_colors14) # Applying custom colors

# Print the plot
print(cluster_plot)
###### K MEDOIDS 16#########

my_colors16 <- c("#008080","#FF0000", "#FF9900","#0000FF", "#9900FF",
                 "#FFFF33", "#339900", "#FF99FF", "#0099CC", "#99FF00",
                 "#003399","#FFCC00", "#FF0066","#00CCCC", "#CC0066", "#3300CC")

set.seed(10)

cl16 <- clara(data, 16, metric= "euclidean", samples = 50, pamLike=T)

sil16c <- silhouette(cl16$clustering, dist(data))

plot(sil16c, col = my_colors16, border = NA, main = "")

data$Cluster16c <- as.factor(cl16$clustering)

# Scatter plot of the clustered data
cluster_plot <- ggplot(data, aes(x = X1, y = X2, color = Cluster16c)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "",
       x = "X1",
       y = "X2",
       color="Cluster") +
  scale_color_manual(values = my_colors16) # Applying custom colors

# Print the plot
print(cluster_plot)
###### Outlier analysis##############
 
set.seed(10)
kmm1 <- kmeans(clean_data, 15, nstart = 50)

my_colors16 <- c("#008080","#FF0000", "#FF9900","#0000FF", "#9900FF",
                 "#FFFF33", "#339900", "#FF99FF", "#0099CC", "#99FF00",
                 "#003399","#FFCC00", "#FF0066","#00CCCC", "#CC0066")

outlier_indices <- c(1201, 2282, 2317, 2322, 2382, 2479, 2540, 4722, 4723, 4809)

df_ggplot <- cbind(clean_data, clusters = as.factor(kmm1$cluster))

ggplot(df_ggplot, aes(x = X1, y=X2, col = clusters)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  scale_color_manual(values = my_colors16) +
  labs(x = "X1", y = "X2", color = "clusters") +
  geom_point(df_ggplot[outlier_indices, ],
             mapping = aes(x = 'X1', y='X2', col = "black"),
             col ="black", shape = 22, fill = "black", size = 2)
######Post Processing:########


ll <- silhouette(kmm1$cluster, dist(clean_data))

negative_indices <- which(ll[, 3] < 0)
negative_scores <- silhouette_scores[which(silhouette_scores$sil_width < 0 ), ]
ll[, 1][negative_indices] <- ll[, 2][negative_indices]

plot(ll, col = my_colors16, border=NA, main="", cex = 1)

new_df <- cbind(data, cluster = as.factor(ll[, 1]))

ggplot(new_df, aes(x = X1, y=X2, col = cluster)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  scale_color_manual(values = my_colors16) +
  labs(x = "X1", y = "X2", color = "clusters")