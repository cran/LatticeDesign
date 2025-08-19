\name{LLHD}
\alias{LLHD}
\title{Optimal generator vector for Lattice-based Latin hypercube designs}
\description{
Generates the optimal generator vector for lattice-based Latin hypercube designs. 
}
\usage{
LLHD(n,d,criterion="WS",T=1000,nstart=max(c(floor(T/(length(CoprimeNumbers(n))*d*5)),1)));  
}
\arguments{
\item{n}{Number of points, must be an integer greater than one.}
\item{d}{Number of dimensions, must be an integer greater than one.}
\item{criterion}{Optional, the criterion used in optimization, which can be "WS", "WP", "WD", "bWS", "bWF", or "bWD".}
\item{T}{Optional, number of iterations.}
\item{nstart}{Optional, number of random starts.}
}
\details{
This function generates the optimal generator vector for lattice-based Latin hypercube designs, following the algorithms provided in the paper "optimized and regularly repeated lattice-based Latin hypercube designs for large-scale computer experiment".  
}
\value{
The value returned from the function gives a generator vector for lattice-based Latin hypercube designs. 
}
\references{
Xu He, Junpeng Gong, and Zhaohui Li (2025) "optimized and regularly repeated lattice-based Latin hypercube designs for large-scale computer experiment", <arXiv:2506.04582v1>}
\examples{
set.seed(85)
v <- LLHD(n=20,d=3,criterion="WS")
design <- LLHDpoints(n=20,v=v,delta=sample(x=1:20,size=3,replace=TRUE))
design
}
