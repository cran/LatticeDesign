\name{SlicedRSPD}
\alias{SlicedRSPD}
\title{Sliced rotated sphere packing designs by partitioning a design}
\description{
Generates a sliced rotated sphere packing design by partitioning one rotated sphere packing design. 
}
\usage{
SlicedRSPD(p=2,n,rotation="magic",w=100)
}
\arguments{
\item{p}{Number of dimensions, must be an integer greater than one.}
\item{n}{Number of points, must be a positive integer.}
\item{rotation}{Optional, whether to use magic rotation matrices (for p=2, recommended) or random rotation matrices.}
\item{w}{Number of rotation matrices to try.}
}
\details{
This function generates a rotated sphere packing design and the slice indexes of points. 
}
\value{
The value returned from the function is a list containing the following components:
\item{Design}{The generated design.}
\item{slices}{The slice indexes of design points.}
\item{generator}{The generator matrix.}
\item{rotation}{The rotation matrix.}
\item{delta}{The value of parameter delta.}
\item{Theta}{The value of parameter Theta.}
\item{l}{The value of parameter l.}
\item{FillDistance}{The fill distance of the design for the nonboundary region.}
}
\references{
He, Xu (2019). "Sliced rotated sphere packing designs", \emph{Technometrics}, 61(1): 66-76. 
}
\examples{
SlicedRSPD(p=2,n=50,rotation="magic",w=100)
}
