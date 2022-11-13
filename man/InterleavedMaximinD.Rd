\name{InterleavedMaximinD}
\alias{InterleavedMaximinD}
\alias{InterleavedMaximinDAlg1}
\alias{InterleavedMaximinDAlg2}
\alias{InterleavedMaximinDAlg3}
\title{Interleaved lattice-based maximin distance designs}
\description{
Generates an interleaved lattice-based maximin distance design. 
}
\usage{
InterleavedMaximinD(p,n,weight=rep(1,p));
InterleavedMaximinDAlg1(p,n,weight=rep(1,p));
InterleavedMaximinDAlg2(p,n,weight=rep(1,p));
InterleavedMaximinDAlg3(p,n,weight=rep(1,p));
}
\arguments{
\item{p}{Number of dimensions, must be an integer greater than one.}
\item{n}{Targeted number of points, must be an integer greater than one.}
\item{weight}{Optional, the weights used in the distance measure, higher for more important variable.}
}
\details{
This function generates an interleaved lattice-based maximin distance design in p dimensions and at least n points, following the algorithms provided in the paper "Interleaved lattice-based maximin distance designs". Function InterleavedMaximinD uses the recommended algorithm provided in the paper. Functions InterleavedMaximinDAlg1, InterleavedMaximinDAlg2, and InterleavedMaximinDAlg3 use Algorithm 1, 2, and 3, respectively. For InterleavedMaximinDAlg1, p must be no greater than 5. For InterleavedMaximinDAlg3, p must be greater than 8. 
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
He, Xu (2019). "Interleaved lattice-based maximin distance designs", \emph{Biometrika}, 106(2): 453-464. 
}
\examples{
InterleavedMaximinD(p=3,n=10,weight=rep(1,3));
InterleavedMaximinDAlg1(p=3,n=10);
InterleavedMaximinDAlg2(p=6,n=10);
InterleavedMaximinDAlg3(p=9,n=257);
}
