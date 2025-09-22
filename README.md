# Analytics_A3

## Overview

This project provides a detailed cluster analysis of a two-dimensional numerical dataset. It employs and compares two prominent clustering algorithms: K-Means and K-Medoids (using the CLARA implementation). The analysis includes a thorough exploratory data analysis (EDA), robust hyperparameter tuning for selecting the optimal number of clusters, an investigation into initialization sensitivity, and post-processing techniques to refine the resulting clusters.

## Features

*   **Exploratory Data Analysis (EDA):** Comprehensive initial analysis including data profiling, visualization with marginal density plots, distance distribution analysis, and outlier identification.
*   **Hyperparameter Tuning:** Systematic selection of the optimal number of clusters (K) using both the average Silhouette Score and the Gap Statistic.
*   **Initialization Sensitivity:** A parallelized approach to run clustering with multiple initializations (up to 200 runs), ensuring the stability and robustness of the chosen K.
*   **Comparative Cluster Analysis:** Implementation of both K-Means and K-Medoids (CLARA) algorithms to compare their performance and results on the dataset.
*   **In-depth Evaluation:** Detailed evaluation of cluster quality using Silhouette plots and visualizations of the clustered data points.
*   **Outlier Impact Assessment:** Analysis of how outliers affect clustering results by comparing models trained on the raw vs. cleaned dataset.
*   **Post-Processing:** A method to improve cluster cohesion by reassigning points with negative silhouette scores to their nearest neighboring cluster and evaluating the improvement.

## Dataset

The analysis is performed on the `STA4026_Assignment_Clustering.txt` dataset.

*   **Dimensions:** 5000 observations and 2 numerical features (V1, V2).
*   **Format:** A text file with two space-separated columns.
*   **Quality:** The dataset contains no missing or duplicate values. A small number of outliers are identified and handled during the analysis.

## Analysis Workflow

The core logic is contained in the `NDXNES005_CHTDOM001_Analytics_2025_A3.R` script. The workflow is as follows:

1.  **Data Loading and EDA:** The dataset is loaded, and an initial exploratory analysis is performed to understand its characteristics. Outliers are identified, and a cleaned version of the dataset is created for the main analysis.
2.  **Hyperparameter Tuning:** The script determines the optimal number of clusters (K) by iterating through a range of values (2-20). It evaluates this range using the Silhouette and Gap Statistic methods for both K-Means and K-Medoids.
3.  **Sensitivity Analysis:** To ensure the chosen K is robust, the tuning process is repeated in parallel, and the frequency of each resulting optimal K is plotted.
4.  **Clustering and Evaluation:** Using the identified optimal K values, both K-Means and K-Medoids algorithms are applied to the cleaned data. The results are visualized through cluster scatter plots and evaluated using Silhouette plots.
5.  **Outlier & Post-Processing Analysis:** The impact of outliers is assessed by running the final K-Means model on the original (uncleaned) data. A post-processing step is also implemented to reassign poorly clustered points and measure the improvement in the average silhouette score.

## Requirements

The analysis requires R and the following packages. You can install them by running this command in your R console:

```r
install.packages(c("tidyverse", "ggplot2", "ggExtra", "dplyr", "GGally", "parallel", "cluster", "factoextra"))
```

## How to Run

1.  Clone or download this repository to your local machine.
2.  Set your R working directory to the repository's root folder. Ensure the data file `STA4026_Assignment_Clustering.txt` is present.
3.  Open an R or RStudio session and execute the analysis script:

    ```r
    source("NDXNES005_CHTDOM001_Analytics_2025_A3.R")
    ```
4.  The script will perform the full analysis, printing results to the console and generating plots. Visualizations may also be saved as PNG or JPEG files in the working directory.
