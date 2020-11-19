\name{InterleavedMinimaxD}
\alias{InterleavedMinimaxD}
\title{Interleaved lattice-based minimax distance designs}
\description{
Generates an interleaved lattice-based minimax distance design. 
}
\usage{
InterleavedMinimaxD(p,n,maxdissimilarity=2*p);
}
\arguments{
\item{p}{Number of dimensions, must be an integer between 2 and 8.}
\item{n}{Targeted number of points, must be an integer greater than one.}
\item{maxdissimilarity}{Optional, the maximum dissimilarity allowed for the number of levels.}
}
\details{
These functions generate an interleaved lattice-based minimax distance design in p dimensions and at most n points, following the algorithm provided in the paper "Interleaved lattice-based minimax distance designs".   
}
\value{
The value returned from the function is a list containing the following components:
\item{Design}{The generated design.}
\item{TargetFillDistance}{The target fill distance, an estimate of the fill distance.}
\item{ActualSize}{The actual number of points of the generated design.}
\item{s_vector}{The numbers of distinct levels of the generated design.}
\item{L01}{The base design.}
}
\references{
He, Xu (2017). "Interleaved lattice-based minimax distance designs", \emph{Biometrika}, 104(3): 713-725.
}
\examples{
InterleavedMinimaxD(p=2,n=20);
}
