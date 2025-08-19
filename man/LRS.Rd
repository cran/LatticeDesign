\name{LRS}
\alias{LRS}
\title{Vertexes of a polytope giving halfspace definition}
\description{
Computes the radius, widths, and vertexes of a polytope giving halfspace definition. 
The program is a R shell of LRS (v.5.1a with lrsmp.h), a reverse search vertex enumeration program/CH package in C which is developed by David Avis <http://cgm.cs.mcgill.ca/~avis/C/lrs.html>. 
Consider the problem of Ax<=b, where A is an n*p matrix, x is a p-vector, and b is an n-vector. Please make sure that the solution of x is nonempty and bounded. Then the nonequalities give the halfspace definition of a polytope. Also make sure that A and b are rational numbers. 
}
\usage{
LRS(numerator,denominator);
}
\arguments{
\item{numerator}{The numerators of cbind(b,A), an n*(p+1) matrix of integer numbers.}
\item{denominator}{The denominators of cbind(b,A), an n*(p+1) matrix of integer numbers.}
}
\details{
This function computes the radius, widths, and vertexes of a polytope giving halfspace definition. 
It is used in constructing interleaved lattice-based minimax distance designs. 
Currently only tested when the maximum values of numerators and denominators are below 2^20. 
If the nonequalities are not defined by rational numbers, round-up to small rational numbers is needed before calling the function. 
The computation is slow for large p but very fast for slow p. Avoid redundant nonequalities may accelerate the calculation. 
}
\value{
The value returned from the function is a list containing the following components:
\item{Radius}{The maximum L2 distance of vertexes to the origin.}
\item{MaxValue}{The maximum k-dimensional value of the vertexes, for k from 1 to p.}
\item{Vertexes}{The vertexes of the polytope.}
}
\references{
Avis, David. LRS, \url{http://cgm.cs.mcgill.ca/~avis/C/lrs.html}.

He, Xu (2017). "Interleaved lattice-based minimax distance designs", \emph{Biometrika}, 104(3): 713-725.
}
\examples{
num = matrix(0,5,3)
den = matrix(1,5,3)
num[1,2] = -1;  den[1,2] = 2; 
num[1,1] = 1;   den[1,1] = 8;
num[2,3] = -1;
num[2,1] = 1;   den[2,1] = 2;
num[3,2] = -1;  den[3,2] = 4;
num[3,3] = -1;  den[3,3] = 2;
num[3,1] = 5;   den[3,1] = 32;
num[4,2] = 1; 
num[4,1] = 0;
num[5,3] = 1;
num[5,1] = 0; 
LRS(num,den)
}
\seealso{
\code{\link{InterleavedMinimaxD}}.
}