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

duplicates <- any(duplicated(data)) # flase - no duplicates

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
    bins = 50, fill = "lightskyblue", color = "black") +
  geom_density(aes(y = after_stat(count) * (max(dist_vals) - min(dist_vals))/50),
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
  geom_point(data = outliers, color = "red", size = 3) 

#------------Correlation Analysis---------------#
cor_mat <- cor(data_cleaned$V1, data_cleaned$V2)
