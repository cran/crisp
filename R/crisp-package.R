#' crisp: A package for fitting a model that partitions the covariate space into blocks in a data-adaptive way.
#'
#' This package is called crisp for "Convex Regression with Interpretable Sharp Partitions",
#' which considers the problem of predicting an outcome variable on the basis of two covariates,
#' using an interpretable yet non-additive model. CRISP partitions the covariate space into
#' blocks in a data-adaptive way, and fits a mean model within each block.
#' Unlike other partitioning methods, CRISP is fit using a non-greedy approach by solving a
#' convex optimization problem, resulting in low-variance fits. More details are provided
#' in Petersen, A., Simon, N., and Witten, D. (2016). Convex Regression with Interpretable
#' Sharp Partitions. Journal of Machine Learning Research, 17(94): 1-31 <http://jmlr.org/papers/volume17/15-344/15-344.pdf>.
#'
#' The main functions are: (1)\code{\link{crisp}} and (2)\code{\link{crispCV}}. The first function
#' \code{\link{crisp}} fits CRISP for a sequence of tuning parameters and provides the fits
#' for this entire sequence of tuning parameters. The second function \code{\link{crispCV}} considers
#' a sequence of tuning parameters and provides the fits, but also returns the optimal tuning parameter,
#' as chosen using K-fold cross-validation.
#'
#' @examples
#' \dontrun{
#' #general example illustrating all functions
#' #see specific function help pages for details of using each function
#'
#' #generate data (using a very small 'n' for illustration purposes)
#' set.seed(1)
#' data <- sim.data(n = 15, scenario = 2)
#' #plot the mean model for the scenario from which we generated data
#' plot(data)
#'
#' #fit model for a range of tuning parameters, i.e., lambda values
#' #lambda sequence is chosen automatically if not specified
#' crisp.out <- crisp(X = data$X, y = data$y)
#' #or fit model and select lambda using 2-fold cross-validation
#' #note: use larger 'n.fold' (e.g., 10) in practice
#' crispCV.out <- crispCV(X = data$X, y = data$y, n.fold = 2)
#'
#' #summarize all of the fits
#' summary(crisp.out)
#' #or just summarize a single fit
#' #we examine the fit with an index of 25. that is, lambda of
#' crisp.out$lambda.seq[25]
#' summary(crisp.out, lambda.index = 25)
#' #lastly, we can summarize the fit chosen using cross-validation
#' summary(crispCV.out)
#' #and also plot the cross-validation error
#' plot(summary(crispCV.out))
#' #the lambda chosen by cross-validation is also available using
#' crispCV.out$lambda.cv
#'
#' #plot the estimated relationships between two predictors and outcome
#' #do this for a specific fit
#' plot(crisp.out, lambda.index = 25)
#' #or for the fit chosen using cross-validation
#' plot(crispCV.out)
#'
#' #we can make predictions for a covariate matrix with new observations
#' #new.X with 20 observations
#' new.data <- sim.data(n = 20, scenario = 2)
#' new.X <- new.data$X
#' #these will give the same predictions:
#' yhat1 <- predict(crisp.out, new.X = new.X, lambda.index = crispCV.out$index.cv)
#' yhat2 <- predict(crispCV.out, new.X = new.X)
#' }
#' @docType package
#' @name crisp-package
#' @aliases crisp-package
NULL
