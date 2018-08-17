## Rotated sphere packing designs and sliced rotated sphere packing designs
## Author: Xu He, hexu@amss.ac.cn

mod = function(a,b) { return(a-b*floor(a/b)) }

fgr2 <- function(D) {
	dis <- (dist(D[, 1]))^2
	for (j in 2:dim(D)[2]) dis <- dis * (dist(D[, j]))^2
	DIST <- as.matrix(dis)
	fn <- sum(1/dis)
	lfn <- log(fn)
	return(lfn)
}
  
addfrom <- function(p,tryF,ballr2) {
	sumF <- sum(tryF) 
	if( sum(tryF^2)*(p+1)/p - sumF^2/p > ballr2 ) return(0)
	dd2m <- rep(0,2*p+2)
	dd2m[1:p] <- -2*(p+1)*tryF +2*sumF +p
	dd2m[(p+1):(2*p)] <- 2*(p+1)*tryF -2*sumF +p
	dd2m[2*p+1] <- -2*sumF +p
	dd2m[2*p+2] <- 2*sumF +p
	minv <- 0
	mind <- 0
	for(s in 1:(2*p+2)) if(dd2m[s]<minv) { minv <- dd2m[s]; mind <- s; }
	return(mind)
}

addAm <- function(ballr,mm) {
	if(mm==1) FE <- (-floor(ballr)):floor(ballr)   
	if(mm>1) { 
		FE = cbind(0, addAm(sqrt(ballr^2),mm-1));
		if(floor(ballr)>0) for(s in 1:floor(ballr)) { FEE <- addAm(sqrt(ballr^2-s^2),mm-1); FE <- rbind( FE, cbind(rep(s,length(FEE)/(mm-1)),FEE) ); }
		if(floor(ballr)>0) for(s in (-1):(-floor(ballr))) { FEE <- addAm(sqrt(ballr^2-s^2),mm-1); FE <- rbind( FE, cbind(rep(s,length(FEE)/(mm-1)),FEE) ); } 
	}
	return(FE)
}


RSPD <- function(p=2, n, rotation="magic", w=100){

	if(!p>=2|!p==round(p)) stop("p must be an integer greater than one.") 
	if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.") 
	if(!w>=1|!w==round(w)) stop("w must be a positive integer.") 

	G <- matrix(-1,p,p)
	diag(G) <- p-sqrt(p+1)
	G <- G/(sqrt(p+1)-1)/sqrt(p)
	etamin <- sqrt(p+1)/sqrt(2*p)
	detG <- (p+1)^((p-1)/2)/p^{p/2}
	rhoc <- sqrt((p+2)/12)
	Omega <- pi^(p/2)/gamma(p/2+1)
	Theta <- Omega*rhoc^p/detG
	l <- (n*detG)^(1/p)
	
	nsd <- 0
	for(r in (p+1)/p^seq(1.7,.99,-0.05)) {
		ballr <- l/2 * sqrt(p) + rhoc*r 
		ballr2 <- ballr^2
		N <- Omega * ballr^p / detG *1.1 + 1000*p
		
		FE <- matrix(0,N,p)
		nrowFE <- 1
		srowFE <- 1
		while(srowFE<=nrowFE ) {  
			tF <- FE[srowFE,1:p]
			maxF <- max(tF)
			minF <- min(tF)
			for(s in 1:p) if(tF[s]>=max(tF)-1) {
				tryF <- tF
				tryF[s] <- tryF[s]+1
				if(addfrom(p,tryF,ballr2)==s) { 
					nrowFE <- nrowFE +1
					FE[nrowFE,1:p] <- tryF
				}
			}
			for(s in 1:p) if(tF[s]<=min(tF)+1) {
				tryF <- tF
				tryF[s] <- tryF[s]-1
				if(addfrom(p,tryF,ballr2)==p+s) { 
					nrowFE <- nrowFE +1
					FE[nrowFE,1:p] <- tryF
				}
			}
			tryF <- tF+1
			if(addfrom(p,tryF,ballr2)==2*p+1) { 
				nrowFE <- nrowFE +1
				FE[nrowFE,1:p] <- tryF
			}
			tryF <- tF-1
			if(addfrom(p,tryF,ballr2)==2*p+2) { 
				nrowFE <- nrowFE +1
				FE[nrowFE,1:p] <- tryF
			}
			srowFE <- srowFE +1
			if(nrowFE>=dim(FE)[1]-2*p-2) FE=rbind(FE,matrix(0,N,p))
		}
		FE <- FE[1:nrowFE,]
		
		minfgr2 <- 0  
		for(ll in 1:(w*2)) {
			
			CPair <- matrix(0,p*(p-1)/2,2)
			row <- 1
			for(i in 1:(p-1)) for(j in (i+1):p) { CPair[row,1] <- i; CPair[row,2] <- j; row <- row+1; }
			R <- diag(p)
			if(p>2|rotation!='magic') for(a in 1:((p*(p-1)/2))) {
				alpha <- runif(1,0,2*pi)
				thepair <- CPair[a,]
				W <- diag(p)
				W[thepair[1],thepair[1]] <- cos(alpha); W[thepair[2],thepair[2]] <- cos(alpha); 
				W[thepair[1],thepair[2]] <- sin(alpha); W[thepair[2],thepair[1]] <- -sin(alpha); 
				R <- R%*%W
			}
			GR <- G %*% R
			E <- FE %*% GR 

			isE <- FALSE
			isB <- FALSE
			isL <- FALSE
			epsilonV <- matrix(0,3,p)
			for(i in 1:20) {
				epsilon <- runif(p,-rhoc*r,rhoc*r)
				D <- E
				for(j in 1:p) { if(length(D)<=p) break; D <- D[ abs(D[,j]-epsilon[j])<l/2 ,]; }
				if(length(D)<=p) next; 
				if(dim(D)[1]==n) { isE <- TRUE; epsilonV[1,] <- epsilon; break; }
				if(isL==FALSE & dim(D)[1]<n) { isL <- TRUE; epsilonV[2,] <- epsilon; }
				if(isB==FALSE & dim(D)[1]>n) { isB <- TRUE; epsilonV[3,] <- epsilon; }
				if(isL==TRUE & isB==TRUE) break;
			}
			try1 <- i
			if(try1==20) next; 
			try2 <- 0
			try3 <- 0
			
			if(isE==FALSE) for(j in 1:(p-1)) {
				try2 <- j
				epsilon <- c(epsilonV[2,1:j],epsilonV[3,(j+1):p])
				D <- E
				for(j in 1:p) { if(length(D)<=p) break;  D <- D[ abs(D[,j]-epsilon[j])<l/2 ,]; }
				if(length(D)<=p) next; 
				if(dim(D)[1]==n) { isE <- TRUE; epsilonV[1,] <- epsilon; break; }
				if(dim(D)[1]<n) { epsilonV[2,] <- epsilon; break; }
				if(dim(D)[1]>n) { epsilonV[3,] <- epsilon; }
			}

			if(isE==FALSE) for(k in 1:20) { 
				try3 <- k
				epsilon <- (epsilonV[2,]+epsilonV[3,])/2
				D <- E
				for(j in 1:p) { if(length(D)<=p) break; D <- D[ abs(D[,j]-epsilon[j])<l/2 ,]; }
				if(length(D)<=p) next; 
				if(dim(D)[1]==n) { isE <- TRUE; epsilonV[1,] <- epsilon; break; }
				if(dim(D)[1]<n) { epsilonV[2,] <- epsilon; }
				if(dim(D)[1]>n) { epsilonV[3,] <- epsilon; }
			}
			if(try3==20) next;
			
			nsd <- nsd+1
			for(j in 1:p) D[,j] <- (D[,j]-epsilon[j]) / l + .5

			fgr2D <- fgr2(D)
			if(nsd==1|fgr2D<minfgr2) {
				Dopt <- D
				minfgr2 <- fgr2D
				Ropt <- R
				epsilonopt <- epsilon
			}
			
			if(nsd>=w) break; 
		}
		
		if(nsd>=w) break; 
	}
		
	D <- Dopt

	list(Design=Dopt,generator=G,rotation=Ropt,delta=-epsilonopt,Theta=Theta,l=l,FillDistance=rhoc/l)
}


AdaptiveRSPD <- function(p=2, n, w=100){

	if(!p>=2|!p==round(p)) stop("p must be an integer greater than one.") 
	if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.") 
	if(!w>=1|!w==round(w)) stop("w must be a positive integer.") 

	G <- matrix(-(sqrt(p+1)+1)/sqrt(2)/p,p,p)
	diag(G) <- diag(G) + sqrt(2)/2
	detG <- abs(det(G))
	intm <- floor((p+1)/2)
	rhoc <- 1/2 * sqrt(2*intm*(p+1-intm)/(p+1))
	Omega <- pi^(p/2)/gamma(p/2+1)
	Theta <- Omega*rhoc^p/detG
	l <- (n*detG)^(1/p)
	
	nsd=0
	for(r in (p+1)/p^seq(1.7,.99,-0.05))
	{
		ballr <- l/2 * sqrt(p) + rhoc*r 
		FE <- addAm(ballr*sqrt(2),p)
		FEdis <- rep(0,dim(FE)[1])
		for(i in 1:dim(FE)[1]) { FEdis[i] <- sum(FE[i,])^2/2 + sum(FE[i,]^2)/2; } 
		FE <- FE[ FEdis <= ballr^2, ]
		
		minfgr2=10^10
		for(ll in 1:(w*2)) {
			
			CPair <- matrix(0,p*(p-1)/2,2)
			row <- 1
			for(i in 1:(p-1)) for(j in (i+1):p) { CPair[row,1] <- i; CPair[row,2] <- j; row <- row+1; }
			R <- diag(p)
			for(a in 1:((p*(p-1)/2))) {
				alpha <- runif(1,0,2*pi)
				thepair <- CPair[a,]
				W <- diag(p)
				W[thepair[1],thepair[1]] <- cos(alpha); W[thepair[2],thepair[2]] <- cos(alpha); 
				W[thepair[1],thepair[2]] <- sin(alpha); W[thepair[2],thepair[1]] <- -sin(alpha); 
				R <- R%*%W
			}
			GR <- G %*% R
			E <- FE %*% GR 

			isE <- FALSE
			isB <- FALSE
			isL <- FALSE
			epsilonV <- matrix(0,3,p)
			for(i in 1:20) {
				epsilon <- runif(p,-rhoc*r,rhoc*r)
				D <- E
				for(j in 1:p) { if(length(D)<=p) break; D <- D[ abs(D[,j]-epsilon[j])<l/2 ,]; }
				if(length(D)<=p) next; 
				if(dim(D)[1]==n) { isE <- TRUE; epsilonV[1,] <- epsilon; break; }
				if(isL==FALSE & dim(D)[1]<n) { isL <- TRUE; epsilonV[2,] <- epsilon; }
				if(isB==FALSE & dim(D)[1]>n) { isB <- TRUE; epsilonV[3,] <- epsilon; }
				if(isL==TRUE & isB==TRUE) break;
			}
			try1 <- i
			if(try1==20) next; 
			try2 <- 0
			try3 <- 0
			
			if(isE==FALSE) for(j in 1:(p-1)) {
				try2 <- j
				epsilon <- c(epsilonV[2,1:j],epsilonV[3,(j+1):p])
				D <- E
				for(j in 1:p) { if(length(D)<=p) break;  D <- D[ abs(D[,j]-epsilon[j])<l/2 ,]; }
				if(length(D)<=p) next; 
				if(dim(D)[1]==n) { isE <- TRUE; epsilonV[1,] <- epsilon; break; }
				if(dim(D)[1]<n) { epsilonV[2,] <- epsilon; break; }
				if(dim(D)[1]>n) { epsilonV[3,] <- epsilon; }
			}

			if(isE==FALSE) for(k in 1:20) { 
				try3 <- k
				epsilon <- (epsilonV[2,]+epsilonV[3,])/2
				D <- E
				for(j in 1:p) { if(length(D)<=p) break; D <- D[ abs(D[,j]-epsilon[j])<l/2 ,]; }
				if(length(D)<=p) next; 
				if(dim(D)[1]==n) { isE <- TRUE; epsilonV[1,] <- epsilon; break; }
				if(dim(D)[1]<n) { epsilonV[2,] <- epsilon; }
				if(dim(D)[1]>n) { epsilonV[3,] <- epsilon; }
			}
			if(try3==20) next;
			
			nsd=nsd+1
			for(j in 1:p) D[,j] <- (D[,j]-epsilon[j]) / l + .5

			fgr2D=fgr2(D)
			if(nsd==1|fgr2D<minfgr2) {
				Dopt <- D
				minfgr2 <- fgr2D
				Ropt <- R
				epsilonopt <- epsilon
			}
			
			if(nsd>=w) break; 
		}
		
		if(nsd>=w) break; 
	}
		
	D0=Dopt
	R=Ropt
	epsilon=epsilonopt*sqrt(2+2/p)

	n1=n*(p+1)
	
	G <- matrix(-1,p,p)
	diag(G) <- p-sqrt(p+1)
	G <- G/(sqrt(p+1)-1)/sqrt(p)
	etamin <- sqrt(p+1)/sqrt(2*p)
	detG <- (p+1)^((p-1)/2)/p^{p/2}
	rhoc <- sqrt((p+2)/12)
	Omega <- pi^(p/2)/gamma(p/2+1)
	Theta <- Omega*rhoc^p/detG
	l <- (n1*detG)^(1/p)

	ballr <- l/2 * sqrt(p) + sqrt(sum(epsilon^2)) 
	ballr2 <- ballr^2
	N <- Omega * ballr^p / detG *1.1 + 1000*p
	
	FE <- matrix(0,N,p)
	nrowFE <- 1
	srowFE <- 1
	while(srowFE<=nrowFE ) {  
		tF <- FE[srowFE,1:p]
		maxF <- max(tF)
		minF <- min(tF)
		for(s in 1:p) if(tF[s]>=max(tF)-1) {
			tryF <- tF
			tryF[s] <- tryF[s]+1
			if(addfrom(p,tryF,ballr2)==s) { 
				nrowFE <- nrowFE +1
				FE[nrowFE,1:p] <- tryF
			}
		}
		for(s in 1:p) if(tF[s]<=min(tF)+1) {
			tryF <- tF
			tryF[s] <- tryF[s]-1
			if(addfrom(p,tryF,ballr2)==p+s) { 
				nrowFE <- nrowFE +1
				FE[nrowFE,1:p] <- tryF
			}
		}
		tryF <- tF+1
		if(addfrom(p,tryF,ballr2)==2*p+1) { 
			nrowFE <- nrowFE +1
			FE[nrowFE,1:p] <- tryF
		}
		tryF <- tF-1
		if(addfrom(p,tryF,ballr2)==2*p+2) { 
			nrowFE <- nrowFE +1
			FE[nrowFE,1:p] <- tryF
		}
		srowFE <- srowFE +1
		if(nrowFE>=dim(FE)[1]-2*p-2) FE=rbind(FE,matrix(0,N,p))
	}
	FE <- FE[1:nrowFE,]
		
	GR <- G %*% Ropt
	E <- FE %*% GR 
	candidates <- cbind(E,FE,matrix(1,dim(FE)[1],1))
	for(j in 1:p) { candidates <- candidates[ abs(candidates[,j]-epsilon[j])<l/2 ,]; }
	for(j in 1:p) candidates[,j] <- (candidates[,j]-epsilon[j]) / l + .5

	for(i in 1:dim(candidates)[1]) if(mod(sum(candidates[i,(p+1):(2*p)]),p+1)==0) candidates[i,2*p+1]=0; 
	candidates <- candidates[ candidates[,2*p+1]==1, 1:p]
	
	list(Design=Dopt,candidates=candidates,generator=G,rotation=Ropt,delta=-epsilonopt,Theta=Theta,l=l,FillDistance=rhoc/l)
}


SlicedRSPD <- function(p=2, n, rotation="magic", w=100){

	if(!p>=2|!p==round(p)) stop("p must be an integer greater than one.") 
	if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.") 
	if(!w>=1|!w==round(w)) stop("w must be a positive integer.") 

	G <- matrix(-1,p,p)
	diag(G) <- p-sqrt(p+1)
	G <- G/(sqrt(p+1)-1)/sqrt(p)
	etamin <- sqrt(p+1)/sqrt(2*p)
	detG <- (p+1)^((p-1)/2)/p^{p/2}
	rhoc <- sqrt((p+2)/12)
	Omega <- pi^(p/2)/gamma(p/2+1)
	Theta <- Omega*rhoc^p/detG
	l <- (n*detG)^(1/p)
	
	nsd <- 0
	for(r in (p+1)/p^seq(1.7,.99,-0.05)) {
		ballr <- l/2 * sqrt(p) + rhoc*r 
		ballr2 <- ballr^2
		N <- Omega * ballr^p / detG *1.1 + 1000*p
		
		FE <- matrix(0,N,p)
		nrowFE <- 1
		srowFE <- 1
		while(srowFE<=nrowFE ) {  
			tF <- FE[srowFE,1:p]
			maxF <- max(tF)
			minF <- min(tF)
			for(s in 1:p) if(tF[s]>=max(tF)-1) {
				tryF <- tF
				tryF[s] <- tryF[s]+1
				if(addfrom(p,tryF,ballr2)==s) { 
					nrowFE <- nrowFE +1
					FE[nrowFE,1:p] <- tryF
				}
			}
			for(s in 1:p) if(tF[s]<=min(tF)+1) {
				tryF <- tF
				tryF[s] <- tryF[s]-1
				if(addfrom(p,tryF,ballr2)==p+s) { 
					nrowFE <- nrowFE +1
					FE[nrowFE,1:p] <- tryF
				}
			}
			tryF <- tF+1
			if(addfrom(p,tryF,ballr2)==2*p+1) { 
				nrowFE <- nrowFE +1
				FE[nrowFE,1:p] <- tryF
			}
			tryF <- tF-1
			if(addfrom(p,tryF,ballr2)==2*p+2) { 
				nrowFE <- nrowFE +1
				FE[nrowFE,1:p] <- tryF
			}
			srowFE <- srowFE +1
			if(nrowFE>=dim(FE)[1]-2*p-2) FE=rbind(FE,matrix(0,N,p))
		}
		FE <- FE[1:nrowFE,]
		
		minfgr2=10^10
		for(ll in 1:(w*2)) {
			
			CPair <- matrix(0,p*(p-1)/2,2)
			row <- 1
			for(i in 1:(p-1)) for(j in (i+1):p) { CPair[row,1] <- i; CPair[row,2] <- j; row <- row+1; }
			R <- diag(p)
			if(p>2|rotation!='magic') for(a in 1:((p*(p-1)/2))) {
				alpha <- runif(1,0,2*pi)
				thepair <- CPair[a,]
				W <- diag(p)
				W[thepair[1],thepair[1]] <- cos(alpha); W[thepair[2],thepair[2]] <- cos(alpha); 
				W[thepair[1],thepair[2]] <- sin(alpha); W[thepair[2],thepair[1]] <- -sin(alpha); 
				R <- R%*%W
			}
			GR <- G %*% R
			E <- FE %*% GR 

			isE <- FALSE
			isB <- FALSE
			isL <- FALSE
			epsilonV <- matrix(0,3,p)
			for(i in 1:20) {
				epsilon <- runif(p,-rhoc*r,rhoc*r)
				D <- E
				for(j in 1:p) { if(length(D)<=p) break; D <- D[ abs(D[,j]-epsilon[j])<l/2 ,]; }
				if(length(D)<=p) next; 
				if(dim(D)[1]==n) { isE <- TRUE; epsilonV[1,] <- epsilon; break; }
				if(isL==FALSE & dim(D)[1]<n) { isL <- TRUE; epsilonV[2,] <- epsilon; }
				if(isB==FALSE & dim(D)[1]>n) { isB <- TRUE; epsilonV[3,] <- epsilon; }
				if(isL==TRUE & isB==TRUE) break;
			}
			try1 <- i
			if(try1==20) next; 
			try2 <- 0
			try3 <- 0
			
			if(isE==FALSE) for(j in 1:(p-1)) {
				try2 <- j
				epsilon <- c(epsilonV[2,1:j],epsilonV[3,(j+1):p])
				D <- E
				for(j in 1:p) { if(length(D)<=p) break;  D <- D[ abs(D[,j]-epsilon[j])<l/2 ,]; }
				if(length(D)<=p) next; 
				if(dim(D)[1]==n) { isE <- TRUE; epsilonV[1,] <- epsilon; break; }
				if(dim(D)[1]<n) { epsilonV[2,] <- epsilon; break; }
				if(dim(D)[1]>n) { epsilonV[3,] <- epsilon; }
			}

			if(isE==FALSE) for(k in 1:20) { 
				try3 <- k
				epsilon <- (epsilonV[2,]+epsilonV[3,])/2
				D <- E
				for(j in 1:p) { if(length(D)<=p) break; D <- D[ abs(D[,j]-epsilon[j])<l/2 ,]; }
				if(length(D)<=p) next; 
				if(dim(D)[1]==n) { isE <- TRUE; epsilonV[1,] <- epsilon; break; }
				if(dim(D)[1]<n) { epsilonV[2,] <- epsilon; }
				if(dim(D)[1]>n) { epsilonV[3,] <- epsilon; }
			}
			if(try3==20) next;
			
			nsd <- nsd+1
			for(j in 1:p) D[,j] <- (D[,j]-epsilon[j]) / l + .5

			fgr2D <- fgr2(D)
			if(nsd==1|fgr2D<minfgr2) {
				Dopt <- D
				minfgr2 <- fgr2D
				Ropt <- R
				epsilonopt <- epsilon
			}
			
			if(nsd>=w) break; 
		}
		
		if(nsd>=w) break; 
	}
		
	D <- Dopt
	DD = D 
	for(j in 1:p) DD[,j] <- (D[,j]- .5)*l +epsilonopt[j]
	DDD = round( DD %*% solve(Ropt) %*% solve(G) )
	slices = DDD[,1] 
	for(i in 1:dim(DDD)[1]) slices[i] = mod(sum(DDD[i,]),p+1)

	list(Design=Dopt,slices=slices,generator=G,rotation=Ropt,delta=-epsilonopt,Theta=Theta,l=l,FillDistance=rhoc/l)
}

