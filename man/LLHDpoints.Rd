\name{LLHDpoints}
\alias{LLHDpoints}
\title{Points of lattice-based Latin hypercube designs}
\description{
Generates the points of a lattice-based Latin hypercube design. 
}
\usage{
LLHDpoints(n,v,delta); 
}
\arguments{
\item{n}{Number of points, must be an integer greater than one.}
\item{v}{Generator vector.}
\item{delta}{Translation parameter.}
}
\details{
This function gives the points of a lattice-based Latin hypercube design, following the formulas provided in the paper "optimized and regularly repeated lattice-based Latin hypercube designs for large-scale computer experiment".  
}
\value{
The value returned from the function gives the design matrix for the lattice-based Latin hypercube design. 
}
\references{
Xu He, Junpeng Gong, and Zhaohui Li (2025) "optimized and regularly repeated lattice-based Latin hypercube designs for large-scale computer experiment", <arXiv:2506.04582v1>}
\examples{
set.seed(85)
v <- LLHD(n=20,d=3,criterion="WS")
design <- LLHDpoints(n=20,v=v,delta=sample(x=1:20,size=3,replace=TRUE))
design
}
