library(DirichletReg)
set.seed(123)  # Set seed to 123

# Define the Dirichlet Mixture Model Class
DirichletMixtureModel <- setRefClass(
  "DirichletMixtureModel",
  fields = list(
    data = "matrix",
    k = "numeric",
    alphas_hats = "matrix",
    pi_hats = "numeric",
    clusters = "matrix"
  ),
  
  methods = list(
    
    # Initialization Method
    initialize = function(data, k) {
      "Initialize with data and number of clusters."
      .self$data <- as.matrix(data)
      .self$k <- k
      params <- .self$initial_params(.self$data, .self$k)
      .self$alphas_hats <- params$alpha
      .self$pi_hats <- params$pi
      .self$clusters <- matrix(nrow=nrow(.self$data), ncol=1) # Empty cluster assignment
    },
    
    # Method for initializing parameters
    initial_params = function(data, k, c=60) {
      "Initializes the parameters (alpha and pi) using K-Means clustering."
      kmeans_result <- kmeans(data, centers = k)
      return(list('alpha' = kmeans_result$centers * c, 'pi' = as.numeric(table(kmeans_result$cluster)/nrow(data))))
    },
    
    # Compute log-likelihood
    log_likelihood = function(x, pi, alpha) {
      "Computes the log-likelihood for the Dirichlet mixture model."
      pi_fj <- matrix(nrow=nrow(alpha), ncol=nrow(x))
      for (j in 1:nrow(alpha)) {
        pi_fj[j, ] <- pi[j] * ddirichlet(x, alpha[j, ])
      }
      return(sum(log(colSums(pi_fj))))
    },
    
    # Compute gamma_ij
    gamma_ij = function(x, pi, alpha) {
      "Computes the responsibility matrix gamma_ij for each data point and cluster."
      pij_fj <- matrix(nrow=nrow(x), ncol=nrow(alpha))
      for (j in 1:nrow(alpha)) {
        pij_fj[, j] <- pi[j] * ddirichlet(x, alpha[j, ])
      }
      return(pij_fj / rowSums(pij_fj))
    },
    
    # Invert the digamma function
    inv_digamma = function(x, tol = 1e-16, max_iter = 100) {
      "Uses Newton-Raphson method to invert the digamma function."
      M <- as.numeric(x >= -2.22)
      y <- M * (exp(x) + 0.5) + (1 - M) * (-1 / (x - digamma(1)))
      for (iter in 1:max_iter) {
        y_new <- y - (digamma(y) - x) / trigamma(y)
        if (abs(y_new - y) < tol) {
          break
        }
        y <- y_new
      }
      return(y)
    },
    
    # Fit the Dirichlet Mixture Model
    fit = function(epsilon = 1e-4, max_iter = 1000) {
      "Fits the Dirichlet Mixture Model using the EM algorithm."
      loglik_old <- .self$log_likelihood(.self$data, .self$pi_hats, .self$alphas_hats)
      log_diff <- epsilon + 1
      iters <- 0
      
      while (log_diff > epsilon && iters < max_iter) {
        # E-step: Compute responsibilities (gamma_ij)
        gammas <- .self$gamma_ij(.self$data, .self$pi_hats, .self$alphas_hats)
        
        # M-step: Update the mixture weights (pi)
        .self$pi_hats <- colSums(gammas) / nrow(.self$data)
        
        # Assign clusters based on responsibilities
        .self$clusters <- matrix(apply(gammas, 1, which.max), ncol = 1)
        cluster_group <- table(factor(.self$clusters, levels = 1:.self$k))
        
        # Update alpha for each cluster
        for (i in 1:.self$k) {
          if (cluster_group[i] == 0) {
            .self$alphas_hats[i, ] <- .self$alphas_hats[i, ] # No update for empty clusters
          } else {
            cluster_indices <- which(.self$clusters == i)
            cluster_data_points <- as.matrix(.self$data[cluster_indices, ])
            psi_a <- digamma(sum(.self$alphas_hats[i, ])) + colMeans(log(cluster_data_points))
            .self$alphas_hats[i, ] <- sapply(psi_a, .self$inv_digamma)
          }
        }
        
        # Compute log-likelihood and check for convergence
        loglik_new <- .self$log_likelihood(.self$data, .self$pi_hats, .self$alphas_hats)
        log_diff <- abs(loglik_old - loglik_new)
        loglik_old <- loglik_new
        iters <- iters + 1
      }
      
      return(list('alphas' = .self$alphas_hats, 'weights' = .self$pi_hats, 'clusters' = .self$clusters))
    },
    
    # Predict clusters for new data
    predict = function(new_data) {
      "Predicts the cluster for new data using the fitted model."
      new_data <- as.matrix(new_data)
      gammas <- .self$gamma_ij(new_data, .self$pi_hats, .self$alphas_hats)
      return(apply(gammas, 1, which.max))
    },
    
    # Get model parameters (alphas and pi)
    get_params = function() {
      "Returns the fitted parameters (alphas and pi)."
      return(list('alphas' = .self$alphas_hats, 'weights' = .self$pi_hats))
    }
  )
)

# Example Usage
sample_data <- rbind(
  rdirichlet(500, c(30, 20, 10)),
  rdirichlet(100, c(0.5, 20, 30)),
  rdirichlet(300, c(5, 5, 5))
)

dmm <- DirichletMixtureModel$new(data = sample_data, k = 3)
dmm$fit()
params <- dmm$get_params()

print(params$alphas)
print(params$weights)
predicted_clusters <- dmm$predict(sample_data)
print(predicted_clusters)

