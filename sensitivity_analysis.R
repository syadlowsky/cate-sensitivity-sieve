#  Copyright 2019, 2020 Steve Yadlowsky
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

library(glmnet)
library(caret)
library(grf)

source("fast_polym.R")

estimate_theta <- function(y, train_x, test_x, gamma = 0, regularization=0.05, degree=1) {
   n <- length(y)

  if (is.null(dim(train_x)) | nrow(train_x) < 4) {
    stop("Need more data")
  }

  if (is.null(dim(test_x))) {
    set_flag <- T
    test_x <- matrix(c(test_x, test_x), nrow=2)
  } else {
    set_flag <- F
  }

  tryCatch({
  poly_model <- fast_polym(train_x, degree=degree)
  train_x <- data.frame(predict(poly_model, train_x))
  test_x <- data.frame(predict(poly_model, test_x))
  train_x <- model.matrix( ~ ., train_x)
  test_x <- model.matrix( ~ ., test_x)
  }, error=function(...) {
    print(train_x)
    print(test_x)
    print(...)
    stop(...)
  })

  if (set_flag) {
    test_x <- matrix(test_x[1,], nrow=1)
  }

  Gamma <- exp(gamma)
  Gamma_weight <- sqrt(Gamma) - 1
  theta <- rep(mean(y), n)
  weights <- rep(1, n)
  weights_last <- weights
  d <- ncol(train_x)

  beta <- rep(0, d)
  for (i in 1:15) {
    x_weight <- weights * train_x
    xtx <- crossprod(x_weight)
    diag(xtx) <- diag(xtx) + regularization
    cond_num <- rcond(xtx)
    if (is.na(cond_num) | cond_num < 1e-8) {
      break
    }
    beta <- solve(xtx, crossprod(x_weight, weights*y))
    theta <- train_x %*% beta
    weights_last <- weights
    weights <- (Gamma_weight) * as.numeric(theta > y) + 1
    if (all(weights==weights_last)) {
      break
    }
  }
  thetas <- as.vector(test_x %*% beta)
  return(thetas)
}

estimate_nu <- function(y, train_x, test_x, regularization=0.05, degree=1) {
  n <- length(y)

  if (is.null(dim(train_x)) | nrow(train_x) < 4) {
    stop("Need more data")
  }

  if (is.null(dim(test_x))) {
    set_flag <- T
    test_x <- matrix(c(test_x, test_x), nrow=2)
  } else {
    set_flag <- F
  }

  tryCatch({
  poly_model <- fast_polym(train_x, degree=degree)
  train_x <- data.frame(predict(poly_model, train_x))
  test_x <- data.frame(predict(poly_model, test_x))
  train_x <- model.matrix( ~ ., train_x)
  test_x <- model.matrix( ~ ., test_x)
  }, error=function(...) {
    print("nu")
    print(train_x)
    print(test_x)
    print(...)
    stop(...)
  })

  if (set_flag) {
    test_x <- matrix(test_x[1,], nrow=1)
  }

  weights <- rep(1, n)
  weights_last <- weights
  alpha <- rep(1.0/n, n)
  d <- ncol(train_x)
  beta <- rep(0, d)
  for (i in 1:10) {
    x_weight <- weights * train_x
    xtx <- crossprod(x_weight)
    diag(xtx) <- diag(xtx) + regularization
    cond_num <- rcond(xtx)
    if (is.na(cond_num) | cond_num < 1e-8) {
      break
    }
    beta <- solve(xtx, crossprod(x_weight, weights*y))
    score <- train_x %*% beta
    weights_last <- weights
    weights <- as.vector(exp(-score/2)/(1+exp(-score)))
    if (any(is.na(weights)) | any(weights < 1e-20)) {
      break
    }
    if (all(abs(weights-weights_last)/weights_last < 1e-3)) {
      break
    }
  }
  nu <- as.vector(test_x %*% beta)
  nu <- 1/(1+exp(-nu))
  return(nu)
}

out_of_sample_residuals <- function(w, y, x, test, k=10, np=F) {
  train_mu_hat <- rep(0, length(y))
  if (np) {
    test_mu_hat <- rep(0, nrow(test))
  } else {
    model_1 <- cv.glmnet(x[w==1,], y[w==1], keep=T, alpha=0)
    train_mu_hat[w==1] <- model_1$fit.preval[, model_1$lambda == model_1$lambda.min]

    test_mu_hat <- predict(model_1, test, s="lambda.min")
  }
  return(list(train=train_mu_hat, test=test_mu_hat))
}

score <- function(y, w, theta_hat, mu_hat, e_hat, weights, cap=0.1) {
  prop_weights <- (1 - e_hat) / e_hat
  prop_weights <- pmin(prop_weights, cap*(sum(prop_weights) - prop_weights))
  psi <- ifelse(w==1, y + weights * (y - theta_hat - mu_hat) * prop_weights,
                theta_hat + mu_hat)
  psi
}

compute_score_internal <- function(y, w, x, gamma, k, fit_hyperparams, folds, e_hat, np=F, cap=0.05) {
  num_s <- length(fit_hyperparams)
  n <- nrow(x)

  data <- data.frame(y, w, x)
  
  thetas <- matrix(0, nrow=n, ncol=num_s)
  nus <- matrix(0, nrow=n, ncol=num_s)
  mus <- matrix(0, nrow=n, ncol=k)
  errors <- matrix(0, nrow=n, ncol=num_s)
  scores <- rep(0, n)
  adj_gamma <- (exp(gamma) - 1)

  for (fold_iter in 1:k) {
    fold <- folds[[fold_iter]]
    train <- data[-fold,]
    test <- data[fold,]
    train_x <- as.matrix(x[-fold,])
    test_x <- as.matrix(x[fold,])

    mu_hat <- out_of_sample_residuals(w[-fold], y[-fold], train_x, test_x, np=np)
    training_residuals = y[-fold] - mu_hat$train
    test_mu <- mu_hat$test
    mus[-fold, fold_iter] <- training_residuals
    mus[fold, fold_iter] <- test_mu

    training_residuals <- training_residuals[train$w==1]
    train_x <- train_x[train$w==1,]
    train <- train[train$w == 1,]

    errors[1,] <- 0
    cv_folds <- createFolds(train$y <= median(train$y), k=10)
    for (cv_fold in cv_folds) {
    for (s_iter in 0:(num_s-1)) {
      test_theta <- estimate_theta(training_residuals[-cv_fold],
                                   train_x[-cv_fold,],
                                   train_x[cv_fold,],
                                   regularization=fit_hyperparams[[1 + s_iter]]$regularization,
                                   gamma=gamma,
                                   degree=fit_hyperparams[[1 + s_iter]]$degree)

      weights <- 1 + adj_gamma * (training_residuals[cv_fold] < test_theta)

      errors[1, s_iter+1] <- errors[1, s_iter+1] + sum(train$w[cv_fold] * weights * (training_residuals[cv_fold] - test_theta)^2)
    }
    }
    best_error_theta <- which.min(errors[1,]) - 1
    for (s_iter in 0:(num_s-1)) {
      test_theta <- estimate_theta(training_residuals,
                                   train_x,
                                   test_x,
                                   regularization=fit_hyperparams[[1 + s_iter]]$regularization,
                                   gamma=gamma,
                                   degree=fit_hyperparams[[1 + s_iter]]$degree)


      weights <- 1 + adj_gamma * (test$y < test_theta + test_mu)

      if (s_iter == best_error_theta) {
        thetas[fold, 1] <- test_theta
      }
    }
  }
  best_error_theta <- 0
  
  for (fold_iter in 1:k) {
    fold <- folds[[fold_iter]]
    train <- data[-fold,]
    test <- data[fold,]
    train_x <- as.matrix(x[-fold,])
    test_x <- as.matrix(x[fold,])

    training_residuals <- mus[-fold, fold_iter]
    test_residuals <- test$y - mus[fold, fold_iter]

    theta_hat <- thetas[,best_error_theta+1]
    training_theta <- theta_hat[-fold][train$w==1]
    test_theta <- theta_hat[fold]

    training_residuals <- training_residuals[train$w==1]
    training_residuals_below_theta <- training_residuals < training_theta
    train_x <- train_x[train$w==1,]
    train <- train[train$w==1,]

    errors[1,] <- 0
    cv_folds <- createFolds(training_residuals_below_theta, k=10)
    for (cv_fold in cv_folds) {
    for (s_iter in 0:(num_s-1)) {
      test_residuals_below_theta <- test_residuals < test_theta
      test_nu <- estimate_nu(training_residuals_below_theta[-cv_fold],
                             train_x[-cv_fold,],
                             train_x[cv_fold,],
                             regularization=fit_hyperparams[[1 + s_iter]]$regularization,
                             degree=fit_hyperparams[[1 + s_iter]]$degree)

      errors[1, s_iter+1] <- errors[1, s_iter+1] + sum(training_residuals_below_theta[cv_fold] - test_nu)^2
    }
    }
    best_error_nu <- which.min(errors[1,]) - 1

    for (s_iter in 0:(num_s-1)) {
      test_residuals_below_theta <- test_residuals < test_theta
      test_nu <- estimate_nu(training_residuals_below_theta,
                             train_x,
                             test_x,
                             regularization=fit_hyperparams[[1 + s_iter]]$regularization,
                             degree=fit_hyperparams[[1 + s_iter]]$degree)

      if (s_iter == best_error_nu) {
        nus[fold, 1] <- test_nu
      }
    }
  }
  best_error_nu <- 0

  for (fold_iter in 1:k) {
    fold <- folds[[fold_iter]]
    train <- data[-fold,]
    test <- data[fold,]
    train_x <- as.matrix(x[-fold,])
    test_x <- as.matrix(x[fold,])

    training_residuals <- mus[-fold, fold_iter]
    test_mu <- mus[fold, fold_iter]
    
    theta_hat <- thetas[,best_error_theta+1]
    training_theta <- theta_hat[-fold]
    test_theta <- theta_hat[fold]

    nu_hat <- nus[,best_error_nu+1]
    training_nu <- nu_hat[-fold]
    test_nu <- nu_hat[fold]

    weights <- (1 + adj_gamma * (test$y < test_theta + test_mu)) / (1 + adj_gamma * test_nu)

    scores[fold] <- score(test$y, test$w, test_theta, test_mu, e_hat[fold], weights = weights, cap=cap)
  }

  return(scores)
}

estimate_propensity <- function(w, x=NULL, folds=NULL, np=F, hyperparams=list(min.node.size=4,num.trees=800,honesty=F,mtry=1.0,alpha=0.1,imbalance.penalty=2.0)) {
  if (np) {
   propensity_model <- grf::regression_forest(x, w,
                                              min.node.size = hyperparams$min.node.size,
				       num.trees=hyperparams$num.trees,
				       honesty=F,
							 mtry = hyperparams$mtry,
							 alpha = hyperparams$alpha,
							 imbalance.penalty = hyperparams$imbalance.penalty)
    fit <- predict(propensity_model)$predictions
    return(fit)
  } else {
    k <- length(folds)
    fold_assignment <- rep(0, length(w))
    for (fold_num in 1:k) {
      fold_assignment[folds[[fold_num]]] = fold_num
    }
    propensity_model <- cv.glmnet(x, w, family="binomial", foldid=fold_assignment, keep=T)
    return(propensity_model$fit.preval[, propensity_model$lambda == propensity_model$lambda.min])
  }
}


sensitivity_analysis <- function(y, w, x, gamma = 0, k=5, fit_hyperparams=list(), propensity_hyperparams=NULL, np=F, np.prop=F) {
  folds <- createFolds(w, k = k)

  if (!is.null(propensity_hyperparams)){
    e_hat <- estimate_propensity(w, x, folds, np.prop, hyperparams=propensity_hyperparams)
  } else {
    e_hat <- estimate_propensity(w, x, folds, np.prop)
  }

  score_1 <- compute_score_internal(y, w, x, gamma, k,fit_hyperparams,  folds, e_hat, np=np)
  score_0 <- -compute_score_internal(-y, 1-w, x, gamma, k,fit_hyperparams, folds, 1-e_hat, np=np)
  score_1_u <- -compute_score_internal(-y, w, x, gamma, k,fit_hyperparams, folds, e_hat, np=np)
  score_0_u <- compute_score_internal(y, 1-w, x, gamma, k,fit_hyperparams, folds, 1-e_hat, np=np)

  tau_hat_lower <- mean(score_1 - score_0)
  tau_hat_upper <- mean(score_1_u - score_0_u)
  tau_hat_lower_se <- sqrt(var(score_1 - score_0) / length(score_0))
  tau_hat_upper_se <- sqrt(var(score_1_u - score_0_u) / length(score_1_u))
  return(list(lower=tau_hat_lower, upper=tau_hat_upper,
              lower.se = tau_hat_lower_se, upper.se = tau_hat_upper_se))
}
