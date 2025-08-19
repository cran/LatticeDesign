\name{InterleavedFillSepD}
\alias{InterleavedFillSepD}
\title{Interleaved lattice-based fill and separation distance designs}
\description{
Generates an interleaved lattice-based designs with low fill and high separation distance properties. 
}
\usage{
InterleavedFillSepD(p,n,w=rep(1,p),pfrom=p,a=1/2,nmin=floor(n*.8),
  nmax=ceiling(n*1.2),coefF=-4,coefS=1,msC=0,NL=10,NP=100,NJ=10,NS=100);
}
\arguments{
\item{p}{Number of dimensions.}
\item{n}{Targeted number of points, must be an integer greater than one.}
\item{w}{Optional, weights of the dimensions.}
\item{pfrom}{Optional, number of dimensions designs are generated and supplemented from, no more than p and no more than 8.}
\item{a}{Optional, translation parameter, with a=0 for uniform design and a=1 for pushing the points to the boundary.}
\item{nmin}{Optional, minimal acceptable number of points, no less than n.}
\item{nmax}{Optional, maximal acceptable number of points, no greater than n.}
\item{coefF}{Optional, coefficient of r_F in the criterion.}
\item{coefS}{Optional, coefficient of r_S in the criterion.}
\item{msC}{Optional, maximal allowed sum of projuct weights for pairs of binary aliased dimensions.}
\item{NL}{Optional, maximum number of lattices to try.}
\item{NP}{Optional, maximum number of dimension permutations to try.}
\item{NJ}{Optional, maximum number of discretization choices for each continuous lattice and s combination.}
\item{NS}{Optional, maximum number of choices to supplement from.}
}
\details{
This function generate an interleaved lattice-based designs with low fill and high separation distance properties in p dimensions and around n points, following the algorithm provided in the paper "Efficient Kriging using interleaved lattice-based designs with low fill and high separation distance properties".   
}
\value{
The value returned from the function is a matrix containing the generated design. Remark that no qualified design might be found if (a) both nmax and nmin are to close to n, (b) both n and msC are small, or (c) both pfrom and msC are small while p-pfrom is big. 
}
\references{
He, Xu (2024). "Efficient Kriging using interleaved lattice-based designs with low fill and high separation distance properties", \emph{SIAM/ASA Journal on Uncertainty Quantification}, 12(4): 1113-1134. 
}
\examples{
InterleavedFillSepD(p=2,n=20);
}
