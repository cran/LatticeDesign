\name{DPMPD}
\alias{DPMPD}
\title{Densest packing-based maximum projection designs}
\description{
Generates a densest packing-based maximum projection design. 
}
\usage{
DPMPD(p,n,rotation="magic",w=100)
}
\arguments{
\item{p}{Number of dimensions, must be an integer greater than one and no higher than eight.}
\item{n}{Number of points, must be an integer greater than one.}
\item{rotation}{Optional, whether to use magic rotation matrices (for p=2,3,4,6,8, recommended) or random rotation matrices.}
\item{w}{Number of rotation matrices to try.}
}
\details{
This function generates a densest packing-based maximum projection design in two to eight dimensions. 
For p=2,4,8 with rotation="magic", the designs are generated following the Biometrika paper "Lattice-based designs possessing quasi-optimal separation distance on all projections". 
For p=3,6 with rotation="magic", the designs are generated following the arXiv paper "Lattice-based designs with quasi-uniform projections". 
For other p or rotation!="magic", the designs are generated from random rotations. 
}
\value{
The value returned from the function is a list containing the following components:
\item{Design}{The generated design.}
\item{ProjectiveSeparationDistance}{The projective separation distance of the generated design, from one-dimensional projections to the unprojected design.}
}
\references{
He, Xu (2021). "Lattice-based designs possessing quasi-optimal separation distance on all projections", \emph{Biometrika}, 108(2): 443-454.

He, Xu (2018). "Lattice-based designs with quasi-uniform projections", arXiv:1709.02062v2. 
}
\examples{
DPMPD(p=4,n=200,w=100)
}
