\name{RSPD}
\alias{RSPD}
\title{Rotated sphere packing designs}
\description{
Generates a rotated sphere packing design. 
}
\usage{
RSPD(p=2,n,rotation="magic",w=100)
}
\arguments{
\item{p}{Number of dimensions, must be an integer greater than one.}
\item{n}{Number of points, must be a positive integer.}
\item{rotation}{Optional, whether to use the magic rotation matrix (for p=2, recommended) or random rotation matrices.}
\item{w}{Number of rotation matrices to try, fixed to 1 when p=2 and rotation="magic".}
}
\details{
This function generates a rotated sphere packing design. 
}
\value{
The value returned from the function is a list containing the following components:
\item{Design}{The generated design.}
\item{generator}{The generator matrix.}
\item{rotation}{The rotation matrix.}
\item{delta}{The value of parameter delta.}
\item{Theta}{The value of parameter Theta.}
\item{l}{The value of parameter l.}
\item{FillDistance}{The fill distance of the design for the nonboundary region.}
}
\references{
He, Xu (2017). "Rotated sphere packing designs", \emph{Journal of the American Statistical Association}, 112(520): 1612-1622. 
}
\examples{
RSPD(p=2,n=50,rotation="magic",w=100)
}
