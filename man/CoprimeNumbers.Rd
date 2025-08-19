\name{CoprimeNumbers}
\alias{CoprimeNumbers}
\title{Coprime mumbers}
\description{
Find all potive numbers that is coprime with n and lower than n/2. 
}
\usage{
CoprimeNumbers(n); 
}
\arguments{
\item{n}{A positive integer.}
}
\details{
This function finds all potive numbers that is coprime with n and lower than n/2.
}
\value{
The value returned from the function gives the vector containing all potive numbers that are coprime with n and lower than n/2. 
}
\references{
Xu He, Junpeng Gong, and Zhaohui Li (2025) "optimized and regularly repeated lattice-based Latin hypercube designs for large-scale computer experiment", <arXiv:2506.04582v1>}
\examples{
CoprimeNumbers(20);
}
