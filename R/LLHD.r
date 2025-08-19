#### Lattice-based Latin hypercube designs


r <- function(x) { return( x-floor(x) ); }  # apply to vectors as well

imods <- function(x,z) {  # Find reminder of x divided by z, x and z can be vectors.
	return(x-z*floor(x/z))
}

primes <- function(n) {  # find all primes no larger than n that divide n. 
	primess = NULL
	thetry = 2
	on = n
	while(thetry <= sqrt(n+.1)) {
		if(imods(n,thetry)==0)  {  primess = c(primess,thetry);  while(round(n/thetry)==n/thetry)  n = n/thetry;  }
		thetry = thetry+1 
	}
	if(n>thetry) primess = c(primess,n)
	return(primess[primess<on])
}

CoprimeNumbers <- function(n) {  # Find the choices of v, i.e., numbers coprime with n and lower than n/2. 
	NCP <- rep(1,floor(n/2))
	for(j in primes(n))  NCP[seq(from=j,to=floor(n/2),by=j)] = 0
	return( (1:floor(n/2))[NCP==1] )
}

WSL2 <- function(n,v) {  
	.C(WSL2_c, PACKAGE="LatticeDesign", as.integer(n), as.integer(v), as.double(0) )[[3]]
}

WFL2 <- function(n,v) {  
	.C(WFL2_c, PACKAGE="LatticeDesign", as.integer(n), as.integer(v), as.double(0) )[[3]]
}

WDL2 <- function(n,v) {  
	.C(WDL2_c, PACKAGE="LatticeDesign", as.integer(n), as.integer(v), as.double(0) )[[3]]
}

WSL <- function(n,v) {  
	.C(WSL_c, PACKAGE="LatticeDesign", as.integer(n), as.integer(length(v)), as.integer(v), as.double(0) )[[4]]
}

WPL <- function(n,v) {  
	.C(WPL_c, PACKAGE="LatticeDesign", as.integer(n), as.integer(length(v)), as.integer(v), as.double(0) )[[4]]
}

WDL <- function(n,v) {  
	.C(WDL_c, PACKAGE="LatticeDesign", as.integer(n), as.integer(length(v)), as.integer(v), as.double(0) )[[4]]
}

bWSL <- function(n,v) {  
	.C(bWSL_c, PACKAGE="LatticeDesign", as.integer(n), as.integer(length(v)), as.integer(v), as.double(0) )[[4]]
}

bWFL <- function(n,v) {  
	.C(bWFL_c, PACKAGE="LatticeDesign", as.integer(n), as.integer(length(v)), as.integer(v), as.double(0) )[[4]]
}

bWDL <- function(n,v) {  
	.C(bWDL_c, PACKAGE="LatticeDesign", as.integer(n), as.integer(length(v)), as.integer(v), as.double(0) )[[4]]
}


## Construct an LLHD. 
LLHD <- function(n,d,criterion="WS",T=1000,nstart= max(c( floor(T/(length(CoprimeNumbers(n))*d*5)), 1)) ){  
	vchoices <- CoprimeNumbers(n)
	if(d>length(vchoices)) return(sample(c(vchoices,LLHD(n=n,d=d-length(vchoices),criterion=criterion,T=T,nstart=nstart))))
	if(d==length(vchoices)) return(sample(vchoices))
	if(d==1) return(sample(vchoices,1))
	cXbest <- Inf
	ttnew <- floor(T/nstart*(0:(nstart-1)))+1
	csns <- rep(0,d)
	for(tt in 1:T){  
		if(sum(tt==ttnew)>0) {
			v <- sample(vchoices,d,replace=FALSE)
			if(criterion=="WS") cX <- WSL(n,v)
			if(criterion=="WP") cX <- WPL(n,v)
			if(criterion=="WD") cX <- WDL(n,v)
			if(criterion=="bWS") {
				cs <- matrix(0,d,d); cX <- 0; 
				for(j in 1:(d-1)) for(k in (j+1):d) { cs[j,k] <- WSL2(n,v[c(j,k)]); cX <- cX + cs[j,k]; }
			}
			if(criterion=="bWF") {
				cs <- matrix(0,d,d); cX <- 0; 
				for(j in 1:(d-1)) for(k in (j+1):d) { cs[j,k] <- WFL2(n,v[c(j,k)]); cX <- cX + cs[j,k]; }
			}
			if(criterion=="bWD") {
				cs <- matrix(0,d,d); cX <- 0; 
				for(j in 1:(d-1)) for(k in (j+1):d) { cs[j,k] <- WDL2(n,v[c(j,k)]); cX <- cX + cs[j,k]; }
			}
			if(cX<cXbest)  { vbest <- v; cXbest <- cX; }
		} 
		else { 
			k <- sample(1:d,1)
			vtry <- v
			vtry[k] <- sample(setdiff(vchoices,v),1)
			if(criterion=="WS") cXtry <- WSL(n,vtry)
			if(criterion=="WP") cXtry <- WPL(n,vtry)
			if(criterion=="WD") cXtry <- WDL(n,vtry)
			if(criterion=="bWS") {
				cXtry <- cX 
				if(k>1) for(j in 1:(k-1)) { csn <- WSL2(n,vtry[c(j,k)]); cXtry <- cXtry-cs[j,k]+csn; csns[j] <- csn; }
				if(k<d) for(j in (k+1):d) { csn <- WSL2(n,vtry[c(k,j)]); cXtry <- cXtry-cs[k,j]+csn; csns[j] <- csn; }
			}
			if(criterion=="bWF") {
				cXtry <- cX 
				if(k>1) for(j in 1:(k-1)) { csn <- WFL2(n,vtry[c(j,k)]); cXtry <- cXtry-cs[j,k]+csn; csns[j] <- csn; }
				if(k<d) for(j in (k+1):d) { csn <- WFL2(n,vtry[c(k,j)]); cXtry <- cXtry-cs[k,j]+csn; csns[j] <- csn; }
			}
			if(criterion=="bWD") {
				cXtry <- cX 
				if(k>1) for(j in 1:(k-1)) { csn <- WDL2(n,vtry[c(j,k)]); cXtry <- cXtry-cs[j,k]+csn; csns[j] <- csn; }
				if(k<d) for(j in (k+1):d) { csn <- WDL2(n,vtry[c(k,j)]); cXtry <- cXtry-cs[k,j]+csn; csns[j] <- csn; }
			}
			if(cXtry<cX)  { 
				v <- vtry; cX <- cXtry; 
				if(criterion=="bWS") if(k>1) for(j in 1:(k-1)) cs[j,k] = csns[j]
				if(criterion=="bWS") if(k<d) for(j in (k+1):d) cs[k,j] = csns[j]
				if(criterion=="bWF") if(k>1) for(j in 1:(k-1)) cs[j,k] = csns[j]
				if(criterion=="bWF") if(k<d) for(j in (k+1):d) cs[k,j] = csns[j]
				if(criterion=="bWD") if(k>1) for(j in 1:(k-1)) cs[j,k] = csns[j]
				if(criterion=="bWD") if(k<d) for(j in (k+1):d) cs[k,j] = csns[j]
			}
			if(cXtry<cXbest)  { vbest <- vtry; cXbest <- cXtry; }
		}
	}
	return(vbest)
}


LLHDpoints <- function(n,v,delta) {  # Generate the design points of an LLHD
	X <- matrix(0,n,length(v))
	for(k in 1:length(v)) X[,k] = r( (0:(n-1))/n*v[k] +delta[k]/n + 1/2/n )
	return(X)
}
