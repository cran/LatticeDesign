\name{LatticeDesign-package}
\alias{LatticeDesign-package}
\docType{package}
\title{LatticeDesign package}
\description{
Generate lattice-based space-filling designs with fill or separation distance properties. 
These include interleaved lattice-based minimax distance designs, 
interleaved lattice-based maximin distance designs, 
interleaved lattice-based designs with both fill and separation distance properties, 
(sliced) rotated sphere packing designs, 
and densest packing-based maximum projections designs. 
}
\details{
\tabular{ll}{
Package: \tab LatticeDesign\cr
Type: \tab Package\cr
Version: \tab 3.0-1\cr
Date: \tab 2025-2-6\cr
License: \tab LGPL-2.1\cr
}
Important functions in this package are: 
\code{\link{InterleavedMinimaxD}} generates an interleaved lattice-based minimax distance design. 
\code{\link{InterleavedMaximinD}} generates an interleaved lattice-based maximin distance design. 
\code{\link{InterleavedFillSepD}} generates an interleaved lattice-based design with both fill and separation distance properties. 
\code{\link{DPMPD}} generates a densest packing-based maximum projection design. 
\code{\link{RSPD}} generates a rotated sphere packing design. 
\code{\link{SlicedRSPD}} generates a sliced rotated sphere packing design by partitioning one rotated sphere packing design. 
\code{\link{AdaptiveRSPD}} generates a sliced rotated sphere packing design by enlarging one rotated sphere packing design. 

All those functions generate space-filling designs with fill or separation distance properties. 
Such designs are useful for accurate emulation of computer experiments, fitting nonparametric models and resource allocation. 
They are constructed from lattices, i.e., sets of points with group structures. 

\code{\link{RSPD}} and \code{\link{DPMPD}} generate designs in two to eight dimensions 
with both unprojected and projective distance properties. 
Such designs are desirable when possibly the output value is insensitive to some variables. 
\code{\link{DPMPD}} can be seen as an upgrade of \code{\link{RSPD}} using new magic rotation matrices. 
Another distinction is that \code{\link{RSPD}} generates designs with better unprojected fill distance for nonboundary regions 
while \code{\link{DPMPD}} generates designs with better unprojected separation distance. 
\code{\link{RSPD}} and \code{\link{DPMPD}} construct designs by rescaling, rotating, translating and extracting 
the points of the lattice with asymptotically optimal fill and separation distance, respectively. 

\code{\link{SlicedRSPD}} and \code{\link{AdaptiveRSPD}} generate sliced rotated sphere packing designs, 
i.e., a rotated sphere packing design that can be partitioned into several smaller rotated sphere packing designs. 
\code{\link{SlicedRSPD}} partitions one rotated sphere packing design. 
The generated designs are useful for computer experiments with a categorical variable, 
computer experiments from multiply resources and model validation. 
Alternatively, \code{\link{AdaptiveRSPD}} enlarges a smaller rotated sphere packing design, 
which is useful for adaptive design of computer experiments. 

\code{\link{InterleavedMinimaxD}} generates designs in two to eight dimensions with low fill distance.  
\code{\link{InterleavedMaximinD}} generates designs with high separation distance. 
\code{\link{InterleavedFillSepD}} generates designs with both properties. 
\code{\link{InterleavedMaximinD}} and \code{\link{InterleavedFillSepD}} allow users to specify the relative importance of variables 
and is applicable to problems with any number of variables. 
These designs are useful for accurate emulation of computer experiments when 
the variables are almost equally important in predicting the output value 
or relatively accurate a priori guess on the variable importance is available.
On the other hand, these designs are poor in projective distance properties 
and are thus not recommended when the output value is insensitive to many unknown variables. 
}
\author{
Maintainer: Xu He <hexu@amss.ac.cn>
}
\references{
He, Xu (2017). "Rotated sphere packing designs", \emph{Journal of the American Statistical Association}, 112(520): 1612-1622. 

He, Xu (2017). "Interleaved lattice-based minimax distance designs", \emph{Biometrika}, 104(3): 713-725.

He, Xu (2018). "Lattice-based designs with quasi-uniform projections", arXiv:1709.02062v2. 

He, Xu (2019). "Interleaved lattice-based maximin distance designs", \emph{Biometrika}, 106(2): 453-464. 

He, Xu (2019). "Sliced rotated sphere packing designs", \emph{Technometrics}, 61(1): 66-76. 

He, Xu (2021). "Lattice-based designs possessing quasi-optimal separation distance on all projections", \emph{Biometrika}, 108(2): 443-454.

He, Xu (2024). "Efficient Kriging using interleaved lattice-based designs with low fill and high separation distance properties", \emph{SIAM/ASA Journal on Uncertainty Quantification}, 12(4): 1113-1134. 
}

\keyword{ package }

