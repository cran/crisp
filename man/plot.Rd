% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CRISP_functions.R
\name{plot}
\alias{plot}
\alias{plot.crisp}
\alias{plot.crispCV}
\title{Plots Fit from \code{\link{crisp}} or \code{\link{crispCV}}.}
\usage{
\method{plot}{crisp}(x, lambda.index, title = NULL, x1lab = NULL,
  x2lab = NULL, min = NULL, max = NULL, cex.axis = 1, cex.lab = 1,
  color1 = "seagreen1", color2 = "steelblue1", color3 = "darkorchid4",
  ...)

\method{plot}{crispCV}(x, title = NULL, x1lab = NULL, x2lab = NULL,
  min = NULL, max = NULL, cex.axis = 1, cex.lab = 1,
  color1 = "seagreen1", color2 = "steelblue1", color3 = "darkorchid4",
  ...)
}
\arguments{
\item{x}{An object of class \code{crisp} or \code{crispCV}, which result from running the \code{\link{crisp}} or \code{\link{crispCV}} functions, respectively.}

\item{lambda.index}{The index for the desired value of lambda, i.e., \code{x$lambda.seq[lambda.index]}.}

\item{title}{The title of the plot. By default, the value of lambda is noted.}

\item{x1lab}{The axis label for the first covariate. By default, it is "X1".}

\item{x2lab}{The axis label for the second covariate. By default, it is "X2".}

\item{min, max}{The minimum and maximum y-values, respectively, to use when plotting the fit. By default, they are chosen to be the minimum and maximum of all of the fits, i.e., the minimum and maximum of \code{unlist(x$M.hat.list)}.}

\item{cex.axis}{The magnification to be used for axis annotation relative to the current setting of \code{cex}.}

\item{cex.lab}{The magnification to be used for x and y labels relative to the current setting of \code{cex}.}

\item{color1, color2, color3}{The colors to use to create the color gradient for plotting the response values. At least the first two must be specified, or the defaults of \code{"seagreen1"}, \code{"steelblue1"}, and \code{"darkorchid4"} will be used.}

\item{...}{Additional arguments to be passed, which are ignored in this function.}
}
\value{
None.
}
\description{
This function plots fit of the class \code{crispCV}, or class \code{crisp} with a user-specified tuning parameter.
}
\examples{
\dontrun{
#See ?'crisp-package' for a full example of how to use this package

#generate data (using a very small 'n' for illustration purposes)
set.seed(1)
data <- sim.data(n = 15, scenario = 2)

#fit model for a range of tuning parameters, i.e., lambda values
#lambda sequence is chosen automatically if not specified
crisp.out <- crisp(X = data$X, y = data$y)
#or fit model and select lambda using 2-fold cross-validation
#note: use larger 'n.fold' (e.g., 10) in practice
crispCV.out <- crispCV(X = data$X, y = data$y, n.fold = 2)

#plot the estimated relationships between two predictors and outcome
#do this for a specific fit
plot(crisp.out, lambda.index = 25)
#or for the fit chosen using cross-validation
plot(crispCV.out)
}
}

