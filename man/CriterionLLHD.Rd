\name{CriterionLLHD}
\alias{WSL}
\alias{WPL}
\alias{WDL}
\alias{bWSL}
\alias{bWFL}
\alias{bWDL}
\alias{WSL2}
\alias{WFL2}
\alias{WDL2}
\title{Criterion value of a lattice-based Latin hypercube design}
\description{
Computes the criterion value of a lattice-based Latin hypercube design. 
}
\usage{
WSL(n,v);
WPL(n,v);
WDL(n,v);
bWSL(n,v);
bWFL(n,v);
bWDL(n,v);
WSL2(n,v);
WFL2(n,v);
WDL2(n,v);
}
\arguments{
\item{n}{Number of points.}
\item{v}{Generator vector of the lattice-based Latin hypercube design.}
}
\details{
These functions compute the criterion value of a lattice-based Latin hypercube design. 
WSL, WPL, and WDL give the wrap-around reciprocal separation distance, wrap-around projective separation distance, and wrap-around discrepancy, respectively, for lattice-based Latin hypercube designs. 
bWSL, bWFL, and bWDL give the bivariate wrap-around reciprocal separation distance, bivariate wrap-around fill distance, and bivariate wrap-around discrepancy, respectively, for lattice-based Latin hypercube designs. 
WSL2, WFL2, and WDL2 give the wrap-around reciprocal separation distance, wrap-around fill distance, and wrap-around discrepancy, respectively, for two-dimensional lattice-based Latin hypercube designs. 
}
\value{
The value returned from the function gives the criterion value of the lattice-based Latin hypercube design.
}
\references{
Xu He, Junpeng Gong, and Zhaohui Li (2025) "optimized and regularly repeated lattice-based Latin hypercube designs for large-scale computer experiment", <arXiv:2506.04582v1>
}
\examples{
WSL(100,c(3,7,13));
WPL(100,c(3,7,13));
WDL(100,c(3,7,13));
WSL2(100,c(3,7));
WSL2(100,c(3,13));
WSL2(100,c(13,7));
WFL2(100,c(3,7));
WDL2(101,c(3,7));
bWSL(100,c(3,7,13));
bWFL(100,c(3,7,13));
bWDL(101,c(3,7,13));
}
