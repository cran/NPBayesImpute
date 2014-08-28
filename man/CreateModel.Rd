\name{CreateModel}
\alias{CreateModel}
\title{
Create and initialize the Rcpp_Lcm model object
}
\description{
CreateModel creates and initializes an Rcpp_Lcm \link{Rcpp_Lcm-class} object for non-parametric multiple imputation of discrete multivariate categorical data with structural zeros.
}
\usage{
CreateModel(X, MCZ, K, Nmax, aalpha, balpha)
}
\arguments{
  \item{X}{a data frame with the dataset with missing values. All variables must be unordered factors.}
  \item{MCZ}{a dataframe with the definition of the structural zeros. Placeholder components are represented with NAs. Variables in MCZ must be factors with the same levels as X. Rows do not need to define disjoint regions of the contingency table. See Manrique-Vallier and Reiter (2014) for details of the definition of structural zeros.}
  \item{K}{the maximum number of mixture components.}
  \item{Nmax}{An upper truncation limit for the augmented sample size.}
  \item{aalpha}{the hyper parameter 'a' for alpha in stick-breaking prior distribution.}
  \item{balpha}{the hyper parameter 'b' for alpha in stick-breaking prior distribution.}
}
\details{
This should be the first function one should call to use the library. The returned model object will be referenced in all subsequent calls. }

\value{

CreateModel returns an Rcpp_lcm object. The returned model object will be referenced in all subsequent calls.
}

\references{
  Manrique-Vallier, D. and Reiter, J.P. (2013), "Bayesian Estimation of Discrete Multivariate Latent Structure Models with Structural Zeros", JCGS. 
  \cr
  Manrique-Vallier, D. and Reiter, J.P. (2014), "Bayesian Multiple Imputation for Large-Scale Categorical Data with Structural Zeros", Survey Methodology.
}
\examples{
\donttest{
require(NPBayesImpute)
data('NYexample')

#create the model
model <- CreateModel(X,MCZ,50,200000,0.25,0.25)

#run 1 burnins, 20 mcmc iterations and thin every 2 iterations
model$Run(1,20,2)

#retrieve parameters from the final iteration
result <- model$snapshot

#convert ImputedX matrix to dataframe, using proper factors/names etc.
ImputedX <- GetDataFrame(result$ImputedX,X)
#View(ImputedX)
}
}
