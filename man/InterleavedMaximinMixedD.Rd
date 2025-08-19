\name{InterleavedMaximinDMixVars}
\alias{InterleavedMaximinDMixVars}
\alias{ILMmDMixVarsAlg6}
\alias{ILMmDMixVarsAlg8}
\alias{ILMmDMixVarsAlg9}
\title{Interleaved lattice-based maximin distance designs for mixed continuous, ordinal, and binary variables}
\description{
Generates an interleaved lattice-based maximin distance design for mixed continuous, ordinal, and binary variables. 
}
\usage{
InterleavedMaximinDMixVars(p,n,discrete_dims,ordinal_levels,weight=rep(1,p));
ILMmDMixVarsAlg6(p,n,discrete_dims,ordinal_levels,weight=rep(1,p));
ILMmDMixVarsAlg8(p,n,discrete_dims,ordinal_levels,weight=rep(1,p),N=1,ReturnAll=0);
ILMmDMixVarsAlg9(p,n,discrete_dims,ordinal_levels,weight=rep(1,p),pfrom=8);
}
\arguments{
\item{p}{Number of dimensions, must be an integer greater than one.}
\item{n}{Targeted number of points, must be an integer greater than one.}
\item{discrete_dims}{Index of discrete variables in 1:p.}
\item{ordinal_levels}{Allowable levels of ordinal variables by a matrix.}
\item{weight}{Optional, the weights used in the distance measure, higher for more important variable.}
\item{N}{Number of results saved.}
\item{ReturnAll}{Whether or not outputing all designs with the highest separation distance, regardless of their sample size.}
\item{pfrom}{Number of dimensions to supplement from.}
}
\details{
These functions generate an interleaved lattice-based maximin distance design for mixed continuous, ordinal, and binary variables in p dimensions and at most n points, following the algorithms provided in the paper "Maximin distance designs for mixed continuous, ordinal, and binary variables". Function InterleavedMaximinDMixVars uses the recommended algorithm provided in the paper. Functions ILMmDMixVarsAlg6, ILMmDMixVarsAlg8, and ILMmDMixVarsAlg9 use Algorithms 6, 8, and 9, respectively. For ILMmDMixVarsAlg6, p must be no greater than 5. For ILMmDMixVarsAlg8, p must be no greater than 8. 
}
\value{
The value returned from the function is a list containing the following components:
\item{Design}{The generated design.}
\item{SeparationDistance}{The separation distance of the generated design.}
\item{m}{The actual number of points of the generated design.}
\item{DesignTransformed}{The generated design that is transformed to the rectangular design space given the weights.}
\item{weight}{The weight used in the distance measure, higher for more important variable.}
\item{s_vector}{The numbers of distinct levels of the generated design.}
\item{L01}{The base design.}
}
\references{
Lan, Hui and He, Xu (2025). "maximin distance designs for mixed continuous, ordinal, and binary variables", arXiv:2507.23405v1. 
}
\examples{
ordinal_levels_inpaper <- list(
  Q1 = c(0, 1),
  Q2 = c(0, 0.5, 1),
  Q3 = c(0, 0.1, 0.3, 0.6, 1),
  Q4 = c(0, 0.25, 0.5, 0.75, 1),
  Q5 = c(0, 0.2, 0.3, 0.5, 0.7, 1),
  Q6 = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  Q7 = c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1),
  Q8 = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
)

discrete_dims <- c(2,4)
ordinal_levels <- list(
  v1 = ordinal_levels_inpaper$Q6,
  v2 = ordinal_levels_inpaper$Q7
)
res <- InterleavedMaximinDMixVars(p=4,n=50,discrete_dims,ordinal_levels,weight=(3/4)^(0:3));
res$Design

discrete_dims <- c()
ordinal_levels <- c()
res <- InterleavedMaximinDMixVars(p=3,n=50,discrete_dims,ordinal_levels,weight=rep(1,3));
res$Design

discrete_dims <- c(1,2,3)
ordinal_levels <- list(
  v1 = ordinal_levels_inpaper$Q8,
  v2 = ordinal_levels_inpaper$Q3,
  v3 = ordinal_levels_inpaper$Q8
)
res <- ILMmDMixVarsAlg6(p=3,n=30,discrete_dims,ordinal_levels);
res$Design

discrete_dims <- c(3,4,5)
ordinal_levels <- list(
  v1 = ordinal_levels_inpaper$Q1,
  v2 = ordinal_levels_inpaper$Q5,
  v3 = ordinal_levels_inpaper$Q8
)
res <- ILMmDMixVarsAlg8(p=5,n=40,discrete_dims,ordinal_levels);
res$Design

discrete_dims <- c(4:10)
ordinal_levels <- list(
  v1 = ordinal_levels_inpaper$Q1,
  v2 = ordinal_levels_inpaper$Q1,
  v3 = ordinal_levels_inpaper$Q1,
  v4 = ordinal_levels_inpaper$Q2,
  v5 = ordinal_levels_inpaper$Q3,
  v6 = ordinal_levels_inpaper$Q4,
  v7 = ordinal_levels_inpaper$Q5
)
res <- ILMmDMixVarsAlg9(p=10,n=60,discrete_dims,ordinal_levels,pfrom=5);
res$Design
}
