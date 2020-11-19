\name{AdaptiveRSPD}
\alias{AdaptiveRSPD}
\title{Sliced rotated sphere packing designs by enlarging a design}
\description{
Generates a sliced rotated sphere packing design by enlarging one rotated sphere packing design. 
}
\usage{
AdaptiveRSPD(p=2,n,w=100)
}
\arguments{
\item{p}{Number of dimensions, must be an integer greater than one.}
\item{n}{Number of points of the small design, must be a positive integer.}
\item{w}{Number of rotation matrices to try.}
}
\details{
This function generates a small rotated sphere packing design and the candidate points for enlarging it. 
}
\value{
The value returned from the function is a list containing the following components:
\item{Design}{The generated design.}
\item{candidates}{The candidate points to add.}
\item{generator}{The generator matrix.}
\item{rotation}{The rotation matrix.}
\item{delta}{The value of parameter delta.}
\item{Theta}{The value of parameter Theta.}
\item{l}{The value of parameter l.}
\item{FillDistance}{The fill distance of the design for the nonboundary region.}
}
\references{
He, Xu (2018). "Sliced rotated sphere packing designs", \emph{Technometrics}, 61(1): 66-76.
}
\examples{
AdaptiveRSPD(p=2,n=50,w=100)
}
