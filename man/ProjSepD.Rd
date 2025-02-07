\name{ProjSepD}
\alias{ProjSepD}
\title{Projective separation distance of a design}
\description{
Computes the projective separation distance of a design. 
}
\usage{
ProjSepD(design);
}
\arguments{
\item{design}{The experimental design, must be a matrix whose rows indicate experimental runs.}
}
\details{
This function computes the squared projective separation distance of a design.  
}
\value{
The value returned from the function gives the squared one-dimensional, two-dimensional, ..., (p-1)-dimensional projective separation distances, and the unprojected separation distance, where p is the number of dimensions of the design.
}
\references{
He, Xu (2021). "Lattice-based designs possessing quasi-optimal separation distance on all projections", \emph{Biometrika}, 108(2): 443-454.
}
\examples{
design = rbind(1:3,c(41,1.2,1.3),c(5.4,5.48,5.7),c(4.3,2.3,2));
ProjSepD(design);
}
