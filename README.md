# Clustering Compositional Data using Dirichlet Mixture Model

This repository provides the R implementation for the paper **"[Clustering compositional data using Dirichlet mixture model](https://doi.org/10.1371/journal.pone.0268438)"** by Samyajoy Pal and Christian Heumann. The method introduced in the paper proposes a clustering algorithm using the Dirichlet mixture model (DMM) without transforming compositional data. The repository includes code for fitting the model, predicting clusters, and retrieving estimated parameters.

## Abstract

A model-based clustering method for compositional data is explored in this article. Most methods for compositional data analysis require some kind of transformation. The proposed method builds a mixture model using Dirichlet distribution  which works with the unit sum constraint. The mixture model uses a hard EM algorithm with some modification to overcome the problem of fast convergence with empty clusters. This work includes a rigorous simulation study to evaluate the performance of the proposed method over varied dimensions, number of clusters, and overlap. The performance of the model is also compared with other popular clustering algorithms often used for compositional data analysis (e.g. KMeans, Gaussian mixture model (GMM) Gaussian Mixture Model with Hard EM (Hard GMM), partition around medoids (PAM), Clustering Large Applications based on Randomized Search (CLARANS), Density-Based Spatial Clustering of Applications with Noise (DBSCAN) etc.) for simulated data as well as two real data problems coming from the business and marketing domain and physical science domain, respectively. The study has shown promising results exploiting different distributional patterns of compositional data.

## Table of Contents

- [Requirements](#requirements)
- [Usage](#usage)
- [Methods](#methods)
- [License](#license)

## Requirements

- R version 3.6 or higher
- R libraries:
  - `DirichletReg`
  - `ggplot2`

You can install these libraries by running the following command in R:

```R
install.packages(c("DirichletReg", "ggplot2"))
```

## usage

1. Clone the repository:
    ```bash
    git clone git@github.com:samyajoypal/HardDMM.git
    cd harddmm
    ```
2. Load the R script into your R environment:
    ```R
    source("DirichletMixtureModel.R")
    ```

3. Example Usage:
    ```R
    # Load required libraries
    library(DirichletReg)

    # Generate sample data from a Dirichlet distribution
    sample_data <- rbind(
      rdirichlet(500, c(30, 20, 10)),
      rdirichlet(100, c(10, 20, 30)),
      rdirichlet(300, c(15, 15, 15))
    )

    # Initialize and fit the model
    dmm <- DirichletMixtureModel$new(data = sample_data, k = 3)
    dmm$fit()

    # Get estimated parameters
    params <- dmm$get_params()
    print(params$alphas)
    print(params$weights)

    # Predict clusters for new data
    predicted_clusters <- dmm$predict(sample_data)
    print(predicted_clusters)

    ```

## Methods

The model employs the following core methods:

1. `fit()`: Fits the Dirichlet Mixture Model to the data using a modified Hard EM algorithm.
2. `predict(new_data)`: Predicts the clusters for new compositional data based on the fitted model.
3. `get_params()`: Returns the estimated model parameters (`alphas` and `pi`).

## License

This repository is licensed under the MIT License. See the `LICENSE` file for more details.
