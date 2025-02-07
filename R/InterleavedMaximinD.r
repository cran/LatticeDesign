#### Interleaved lattice-based maximin distance designs
#### Author: Xu He, hexu@amss.ac.cn

BAdd <- function(a,b) {  # Addition in Z_2, a and b are vectors
	ss <- a+b 
	for(i in 1:length(a)) if(ss[i]>=2) ss[i] <- ss[i]-2
	ss
}

BAdds <- function(a,b) {  # Addition in Z_2, a is a matrix and b is a vector
	ssm <- matrix(0,dim(a)[1],dim(a)[2])
	for(i in 1:dim(a)[1]) {
		ssm[i,] <- a[i,]+b
		for(j in 1:dim(a)[2]) if(ssm[i,j]>=2) ssm[i,j] <- ssm[i,j]-2
	}
	ssm
}

BAddss <- function(a,b) {  # Addition in Z_2, a and b are matrices 
	ssm <- matrix(0,0,dim(a)[2])
	for(i in 1:dim(a)[1]) for(j in 1:dim(b)[1]) {
		ss <- a[i,]+b[j,]
		for(j in 1:dim(a)[2]) if(ss[j]>=2) ss[j] <- ss[j]-2
		ssm <- rbind(ssm,ss)
	}
	ssm
}

AddToGroup <- function(a,b) {  # Add elements of b into the group of a, a is a matrix and b is a vector
	NewGenerator <- 1
	aa <- a
	if(dim(a)[1]==0) { aa <- matrix(0,1,length(b)); for(j in 1:length(b)) aa[1,j] <- b[j]; NewGenerator <- 0; }
	if(dim(a)[1]>0) for(i in 1:dim(a)[1]) if(sum(a[i,]==b)==dim(a)[2]) NewGenerator <- 0
	if(NewGenerator==1) aa <- rbind(a,b,BAdds(a,b))
	aa
}

SiteBinary <- function(a) {  # Converge a binary vector into an integer representing the site of it in the StatusVector
	sum(a*2^(0:(length(a)-1)))
}

r_v <- function(G) {   # compute r_vector(G), the sample sizes for the lattice for s_k being 2 or 1.
	p <- dim(G)[1] 
	is.rep <- rep(0,p)
	for(i in 1:p) if(sum(G[i,])==1) is.rep <- is.rep+G[i,]
	p23 <- p-sum(is.rep)
	r_vector <- rep(0,2^p23)
	L01alt <- NULL
	if(p23>0) {
		Galt <- G[,is.rep==0]
		L01alt <- matrix(0,1,p23)
		for(i in 1:p) if(sum(Galt[i,])>1) if(max(Galt[i,])<2) { 
			Ladd <- L01alt; for(j in 1:p23) Ladd[,j] <- Ladd[,j]+Galt[i,j]; L01alt <- rbind(L01alt,Ladd); }
		L01alt <- L01alt-floor(L01alt/2)*2
		for(i in 1:dim(L01alt)[1]) {
			rtoadd <- 0
			for(j in 1:p23) if(L01alt[i,j]==0) rtoadd <- c(rtoadd,rtoadd+2^(j-1))
			r_vector[rtoadd+1] <- r_vector[rtoadd+1] + 1
		}
	}
	list(r_vector,is.rep,L01alt)
}

Sep_Known <- function(L01alt,is.rep,s_vector,weight) {  # Compute separation distance based on known G (and thus known L01alt).
	p <- length(is.rep)
	SepSquare <- sum(weight^2)
	if(sum(is.rep)>0) SepSquare <- min(weight[is.rep==1]/(s_vector[is.rep==1]-1))^2
	if(sum(is.rep)<p) {
		for(j in 1:p) if(is.rep[j]==0) if(s_vector[j]>2) 
			SepSquare <- min(SepSquare,(2*weight[j]/(s_vector[j]-1))^2)
		for(j in 1:dim(L01alt)[2]) L01alt[,j] <- L01alt[,j] * weight[is.rep==0][j]
		for(i in 2:dim(L01alt)[1]) SepSquare <- min(SepSquare,sum((L01alt[i,]/(s_vector[is.rep==0]-1))^2))
	}
	sqrt(SepSquare)
}

m_exact_unknown <- function(G,s_vector) {   # Compute m(G,s_vector), the sample size for the lattice based on G in \prod_k\{0,s_k-1\}. 
	p <- length(s_vector) 
	is.rep <- rep(0,p)
	for(i in 1:p) if(sum(G[i,])==1) is.rep <- is.rep+G[i,]
	p23 <- p-sum(is.rep)
	if(p23==0) m <- prod(s_vector)
	if(p23>0) {
		Galt <- G[,is.rep==0]
		L01alt <- matrix(0,1,p23)
		for(i in 1:p23) if(sum(Galt[i,])>1) if(max(Galt[i,])<2) { 
			Ladd <- L01alt; for(j in 1:p23) Ladd[,j] <- Ladd[,j]+Galt[i,j]; L01alt <- rbind(L01alt,Ladd); }
		L01alt <- L01alt-floor(L01alt/2)*2
		
		r_vector <- rep(0,2^p23)
		for(i in 1:dim(L01alt)[1]) {
			rtoadd <- 0
			for(j in 1:p23) if(L01alt[i,j]==0) rtoadd <- c(rtoadd,rtoadd+2^(j-1))
			r_vector[rtoadd+1] <- r_vector[rtoadd+1] + 1
		}

		salt <- s_vector[is.rep==0]
		aalt <- floor(salt/2)
		balt <- salt-aalt*2
		m=0
		rtomul <- c(aalt[1],balt[1])
		if(p23>1) for(j in 2:p23) rtomul <- c(rtomul*aalt[j],rtomul*balt[j]) 
		m <- sum(r_vector*rtomul)
		for(j in 1:p) if(is.rep[j]==1) m <- m*s_vector[j]
	}
	m
}

m_exact_known <- function(G,r_vector,s_vector) {   # Compute m(G,s_vector), the sample size for the lattice based on G in \prod_k\{0,s_k-1\}. 
	p <- length(s_vector) 
	is.rep <- rep(0,p)
	for(i in 1:p) if(sum(G[i,])==1) is.rep <- is.rep+G[i,]
	p23 <- p-sum(is.rep)
	if(p23==0) m=prod(s_vector)
	if(p23>0) {
		salt <- s_vector[is.rep==0]
		aalt <- floor(salt/2)
		balt <- salt-aalt*2
		m=0
		rtomul <- c(aalt[1],balt[1])
		if(p23>1) for(j in 2:p23) rtomul <- c(rtomul*aalt[j],rtomul*balt[j]) 
		m <- sum(r_vector*rtomul)
		for(j in 1:p) if(is.rep[j]==1) m <- m*s_vector[j]
	}
	m
}

m_exact_fromL01 <- function(L01alt,is.rep,s_vector) {   # Compute m(L01alt,s_vector), the sample size for the lattice based design in \prod_k\{0,s_k-1\}. 
	p <- length(is.rep)
	p23 <- p-sum(is.rep)
	r_vector <- rep(0,2^p23)
	if(p23>0) {
		for(i in 1:dim(L01alt)[1]) {
			rtoadd <- 0
			for(j in 1:p23) if(L01alt[i,j]==0) rtoadd <- c(rtoadd,rtoadd+2^(j-1))
			r_vector[rtoadd+1] <- r_vector[rtoadd+1] + 1
		}
	}
	if(p23==0) m=prod(s_vector)
	if(p23>0) {
		salt <- s_vector[is.rep==0]
		aalt <- floor(salt/2)
		balt <- salt-aalt*2
		m=0
		rtomul <- c(aalt[1],balt[1])
		if(p23>1) for(j in 2:p23) rtomul <- c(rtomul*aalt[j],rtomul*balt[j]) 
		m <- sum(r_vector*rtomul)
		for(j in 1:p) if(is.rep[j]==1) m <- m*s_vector[j]
	}
	m
}

CheckArrange <- function(p,p3,VectorListCannot,VectorListYes)  {
  # Return 0 if cannot arrange; return number of (nonzero) cosets if can arrange. 
  # VectorListYes must have group structure; VectorListCannot needs not. 
  # VectorListYes can have zero rows; VectorListCannot needs at least one row. 
	StatusList <- rep(-1,2^p-1)
	if(dim(VectorListCannot)[1]>0) for(i in 1:dim(VectorListCannot)[1]) StatusList[SiteBinary(VectorListCannot[i,])] <- 2^p3
	VectorListArrange <- rep(list(matrix(0,0,p)),2^p3-1)
	if(dim(VectorListCannot)[1]>0) for(i in 1:dim(VectorListCannot)[1]) { 
		VectorToArrange <- VectorListCannot[i,]
		if(StatusList[SiteBinary(VectorToArrange)]==0) return(0)
		if(StatusList[SiteBinary(VectorToArrange)]!=2^p3) next 
		j=1; while(j<=2^p3-1) {
			AddHere <- 1 
			if(dim(VectorListArrange[[j]])[1]==0) break
			VectorListToYes <- BAdds(VectorListArrange[[j]],VectorToArrange)
			StatusListNew <- StatusList
			VectorListNewYes <- VectorListYes
			for(ii in 1:dim(VectorListToYes)[1]) {
				temp <- AddToGroup(VectorListNewYes,VectorListToYes[ii,])
				if(dim(temp)[1]>dim(VectorListNewYes)[1]) {
					VectorListNewYes <- temp
				}
			}
			for(iii in 1:dim(VectorListNewYes)[1]) {
				if(StatusListNew[SiteBinary(VectorListNewYes[iii,])]>0) { AddHere <- 0; break; }
				StatusListNew[SiteBinary(VectorListNewYes[iii,])] = 0
			}
			for(jj in 1:(2^p3-1)) if(dim(VectorListArrange[[jj]])[1]>0) { 
				if(AddHere==0) break 
				VectorListToCoset <- BAddss(VectorListArrange[[jj]],VectorListNewYes) 
				for(iii in 1:dim(VectorListToCoset)[1]) {
					OldStatus <- StatusListNew[SiteBinary(VectorListToCoset[iii,])]
					if(OldStatus!=jj & OldStatus!=2^p3 & OldStatus!=-1) { AddHere <- 0; break; }
					StatusListNew[SiteBinary(VectorListToCoset[iii,])] <- jj
				}
			}
			if(AddHere==1) break
			j=j+1
		}
		if(j==2^p3) return(0)
		VectorListArrange[[j]] <- rbind(VectorListArrange[[j]],VectorToArrange)
		StatusList[SiteBinary(VectorToArrange)] <- j
		if(dim(VectorListArrange[[j]])[1]>1) {
			VectorListYes <- VectorListNewYes
			for(jj in 1:(2^p3-1)) if(dim(VectorListArrange[[jj]])[1]>0) { 
				VectorListToCoset <- BAddss(VectorListArrange[[jj]],VectorListYes) 
				for(iii in 1:dim(VectorListToCoset)[1]) {
					if(StatusList[SiteBinary(VectorListToCoset[iii,])]!=jj) {
						VectorListArrange[[jj]] <- rbind(VectorListArrange[[jj]],VectorListToCoset[iii,])
						StatusList[SiteBinary(VectorListToCoset[iii,])] <- jj
					}
				}
			}
			StatusList <- StatusListNew
		}
	}
	for(j in 1:(2^p3-1)) if(dim(VectorListArrange[[j]])[1]==0) return(j-1)
	return(2^p3-1) 
}


InterleavedMaximinDAlg1 <- function(p,n,weight=rep(1,p)) {
	
	if(!p>=2|!p<=5|!p==round(p)) stop("p must be an integer greater than one and no greater than five.") 
	if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.") 
	if(!length(weight)==p) stop("weight must be a vector with length p.") 
	if(p==2) GsName <- data(GeneratorMatrices2, envir=environment())
	if(p==3) GsName <- data(GeneratorMatrices3, envir=environment())
	if(p==4) GsName <- data(GeneratorMatrices4, envir=environment())
	if(p==5) GsName <- data(GeneratorMatrices5, envir=environment())
	Gs <- get(GsName)
	Gsn <- dim(Gs)[1]/p
	SepBest <- 0
	GindexBest <- 0
	s_vectorBest <- rep(0,p)
	mBest <- 0
	
	for(Gindex in 1:Gsn) {
		G <- Gs[(Gindex*p-p+1):(Gindex*p),]
		temp <- r_v(G)
		r_vector <- temp[[1]]
		is.rep <- temp[[2]]
		L01alt <- temp[[3]]
		p3 <- sum(G==2)
		
		s_vector <- rep(2,p)
		s_vector[p] = max(c( ceiling(n*2^p3/prod(ceiling(s_vector[-p]/2)*2) /2)*2-1, 2 ))  
		tobreak <- 0
		while(tobreak==0) {
			m <- prod(ceiling(s_vector[is.rep==0]/2)*2)*prod(s_vector[is.rep==1])/2^p3  # A quick upper bound of m
			if(m>=n) m <- m_exact_known(G,r_vector,s_vector)
			if(m<n) s_vector[p] <- s_vector[p]+1
			if(m>=n) {
				Sep <- Sep_Known(L01alt,is.rep,s_vector,weight)
				if(Sep>SepBest) {
					SepBest <- Sep
					GindexBest <- Gindex
					s_vectorBest <- s_vector
					mBest <- m
					p3Best <- p3
					p1Best <- sum(is.rep)
					if(is.null(L01alt)) L01alt <- matrix(0,1,0)
					L01altBest <- L01alt
					is.repBest <- is.rep
				}

				for(kk in p:1) if(s_vector[kk]>2) break; 
				if(kk<=1) { tobreak<-1; break; }
				kk <- kk-1
				while(kk>0) {   # Check if the kk dimension has no chance
					if(is.rep[kk]==1) if(weight[kk]/s_vector[kk]>=SepBest) break
					if(is.rep[kk]==0) if(weight[kk]/s_vector[kk]*2>=SepBest) break
					kk <- kk-1; 
				}
				if(kk==0) { tobreak<-1; break; }
				s_vector[kk] <- s_vector[kk]+1 
				s_vector[(kk+1):p] <- 2
				s_vector[p] = max(c( ceiling(n*2^p3/prod(ceiling(s_vector[-p]/2)*2) /2)*2-1, 2 ))  
			}
		}
	}

	L01Best <- matrix(0,dim(L01altBest)[1],p)
	L01Best[,is.repBest==0] <- L01altBest
	Design <- L01Best
	L01BestS <- L01Best
	for(j in p:1) {
		k <- j
		if(is.repBest[k]==0) if(s_vectorBest[j]>3) for(i in 2:floor(s_vectorBest[j]/2)) {
			DesignAdd <- L01Best 
			DesignAdd[,j] <- DesignAdd[,j]+(i-1)*2
			Design <- rbind(Design,DesignAdd)
		}
		if(is.repBest[k]==0) if(floor(s_vectorBest[j]/2)*2 != s_vectorBest[j]) {
			DesignAdd <- L01Best[L01Best[,j]==0,,drop=FALSE] 
			DesignAdd[,j] <- DesignAdd[,j]+s_vectorBest[j]-1
			Design <- rbind(Design,DesignAdd)
		}	
		if(is.repBest[k]==1) for(i in 2:s_vectorBest[j]) {
			DesignAdd <- L01Best 
			DesignAdd[,j] <- DesignAdd[,j]+(i-1)
			Design <- rbind(Design,DesignAdd)
		}
		L01Best <- Design
	}
	for(j in 1:p) Design[,j] <- Design[,j]/(s_vectorBest[j]-1)
	DesignTransformed <- Design
	for(j in 1:p) DesignTransformed[,j] <- DesignTransformed[,j]*weight[j]
	
	return(list(Design=Design,SeparationDistance=SepBest,m=mBest,DesignTransformed=DesignTransformed,weight=weight,s_vector=s_vectorBest,L01=L01BestS))
}


InterleavedMaximinDAlg2 <- function(p,n,weight=rep(1,p)) {

	if(!p>=2|!p==round(p)) stop("p must be an integer greater than one.") 
	if(!p<=8) warning("The program can be extremely slow for p>8.")
	if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.") 
	if(!length(weight)==p) stop("weight must be a vector with length p.") 

	SepBest <- 0
	
	## Exact solution for p3=0. 
	
	p3 <- 0 
	s_vector <- rep(2,p)
	m <- prod(s_vector)
	while(m<n) {
		SepEachNext <- weight/s_vector
		s_vector[order(-SepEachNext)[1]] <- s_vector[order(-SepEachNext)[1]]+1
		m <- prod(s_vector)
	}
	SepEach <- weight/(s_vector-1)
	SepUpperBound <- min(SepEach)
	if(SepUpperBound>SepBest) {
		SepBest <- SepUpperBound
		mBest <- m
		p3Best <- p3
		p1Best <- p
		s_vectorBest <- s_vector
		NCannotBest <- 0
		DimOrderBest <- 1:p
		is.repBest <- rep(1,p)
		L01altBest <- matrix(0,1,0)
		VectorGeneratorAltBest <- matrix(0,0,0)
	}
	
	## Fast solution for p3=1, all alt., not guaranteed best s_vector. 
	
	p3 <- 1
	s_vector <- rep(2,p)
	m <- ceiling(prod(s_vector)/2)
	while(m<n) {
		SepEachNext <- weight/s_vector
		s_vector[order(-SepEachNext)[1]] <- s_vector[order(-SepEachNext)[1]]+1
		m <- ceiling(prod(s_vector)/2)
	}
	SepEach <- weight/(s_vector-1)
	SepUpperBound <- min(c( sqrt(min(SepEach)^2+sort(SepEach)[2]^2), min(SepEach) ))
	if(SepUpperBound>SepBest) {
		SepBest <- SepUpperBound
		mBest <- m
		p3Best <- p3
		p1Best <- 0
		s_vectorBest <- s_vector
		NCannotBest <- sum(SepEach<SepUpperBound)
		DimOrderBest <- 1:p
		is.repBest <- rep(0,p)
		L01altBest <- matrix(0,1,p)
		for(j in p:1) { L01altBestAdd <- L01altBest; L01altBestAdd[,j]<-L01altBestAdd[,j]+1; L01altBest <- rbind(L01altBest,L01altBestAdd); }
		temp <- rep(0,2^p)
		for(j in 1:p) temp <- temp + L01altBest[,j]
		L01altBest <- L01altBest[floor(temp/2)*2==temp,]
		VectorGeneratorAltBest <- cbind(rep(1,p-1),diag(p-1))
	}
	
	  ## Search from p3=1. Only consider best Sep for given p3 and s_vector. 
	for(p3 in 1:(p-1)){
		s_vector <- rep(2,p)
		s_vector[p] = max(c( ceiling(n*2^p3/prod(ceiling(s_vector[-p]/2)*2) /2)*2-1, 2 ))  # A quick upper bound of m: prod(ceiling(s_vector[-p]/2)*2)/2^p3
		tobreak <- 0
		while(tobreak==0) {  # m<n: increase last s; IsChance==0: increase other s and set last as 2. 
			SepEach <- weight/(s_vector-1)
			DimMustAlt <- (SepEach<SepBest)
			if(sum(DimMustAlt)<=p3) {  # Force some dimensions to be Alt. 
				DimSeq <- s_vector + (floor(s_vector/2)*2==s_vector)*max(s_vector) + (SepEach<SepBest)*max(s_vector)*2
				DimMustAlt[order(DimSeq)[1:(p3+1-sum(DimMustAlt))]] = TRUE
			}  
			m <- ceiling(prod(s_vector[DimMustAlt])/2) 
			if(p3>=2) m <- m - prod(floor(s_vector[DimMustAlt]/2)*2)/2 + prod(floor(s_vector[DimMustAlt]/2)*2)/2^p3
			if(sum(DimMustAlt)<p) if(sum(DimMustAlt)>0) m <- m * prod(s_vector[DimMustAlt==0])  # Another upper bound of m, good for low p3. p3=1 exact; p3=2 very good. 
			if(m<n) { s_vector[p] <- s_vector[p]+1; next; }

			IsChance <- 1  # Is the choice of p3 and s_vector has a chance to be SepBest? 0 stands for no chance.
			SepUpperBound <- sqrt(p) 
			if(max(s_vector)>2) SepUpperBound <- min(SepEach[s_vector>2]*2)
			SepUpperBound <- min(c(SepUpperBound,sqrt(sum(SepEach^2))))
			if(SepUpperBound<=SepBest) IsChance <- 0

			if(IsChance==1) {
				DimOrder <- order(SepEach)
				wdsOrdered <- weight[DimOrder] / (s_vector[DimOrder]-1) # weight divided by (s-1) ordered.
				VectorInConsideration <- cbind(diag(p),wdsOrdered)
				colnames(VectorInConsideration) <- NULL
				VectorListCannot <- matrix(0,0,p)

				  # Add vectors that must be added. 
				VectorToCannotIndex <- order(VectorInConsideration[,p+1])[1]
				VectorToCannot <- VectorInConsideration[VectorToCannotIndex,1:p]
				while(sqrt(sum((VectorToCannot*wdsOrdered)^2))<=SepBest+10^-14) {
					for(kkk in 1:p) if(VectorToCannot[kkk]>0) break
					if(kkk-1>0) for(i in 1:(kkk-1)) {
						TempVector <- VectorToCannot
						TempVector[i] <- TempVector[i]+1
						TempVector <- c(TempVector,sqrt(sum((TempVector*wdsOrdered)^2)))
						VectorInConsideration <- rbind(VectorInConsideration,TempVector)
					}
					VectorInConsideration <- VectorInConsideration[-VectorToCannotIndex,,drop=FALSE]
					VectorListCannot <- rbind(VectorListCannot,VectorToCannot) 
					VectorToCannotIndex <- order(VectorInConsideration[,p+1])[1]
					VectorToCannot <- VectorInConsideration[VectorToCannotIndex,1:p]
				}
			}
					
			  # Increase nonrep dimensions if necessary. 
			if(IsChance==1) {
				is.rep <- rep(1,p) 
				for(j in 1:p) is.rep[j] <- (sum(VectorListCannot[,j])==0)
				while(sum(is.rep)>p-p3-1) {
					VectorToCannotIndex <- 1
					VectorToCannot <- VectorInConsideration[VectorToCannotIndex,1:p]
					for(kkk in 1:p) if(VectorToCannot[kkk]>0) break
					if(kkk-1>0) for(i in 1:(kkk-1)) {
						TempVector <- VectorToCannot
						TempVector[i] <- TempVector[i]+1
						TempVector <- c(TempVector,sqrt(sum((TempVector*wdsOrdered)^2)))
						VectorInConsideration <- rbind(VectorInConsideration,TempVector)
					}
					VectorInConsideration <- VectorInConsideration[-VectorToCannotIndex,,drop=FALSE]
					VectorListCannot <- rbind(VectorListCannot,VectorToCannot) 
					for(j in 1:p) is.rep[j] <- (sum(VectorListCannot[,j])==0)
				}
			}
			  
			  # Check if the vectors can be arranged (check if u is large enough). 
			if(IsChance==1) { 
				VectorListMust <- matrix(0,0,p)
				ncosets <- CheckArrange(p,p3,VectorListCannot,VectorListMust)
				IsChance <- (ncosets>0) 
			}
			
			  # Now check subproblems for each p1
			if(IsChance==1) { 
				Op1 <- sum(is.rep)
				for(p1aim in Op1:0) {
					if(p1aim!=Op1) {	# Increase nonrep dimensions by one. 
						VectorToCannotIndex <- 1
						VectorToCannot <- VectorInConsideration[VectorToCannotIndex,1:p]
						for(kkk in 1:p) if(VectorToCannot[kkk]>0) break
						if(kkk-1>0) for(i in 1:(kkk-1)) {
							TempVector <- VectorToCannot
							TempVector[i] <- TempVector[i]+1
							TempVector <- c(TempVector,sqrt(sum((TempVector*wdsOrdered)^2)))
							VectorInConsideration <- rbind(VectorInConsideration,TempVector)
						}
						VectorInConsideration <- VectorInConsideration[-VectorToCannotIndex,,drop=FALSE]
						VectorListCannot <- rbind(VectorListCannot,VectorToCannot) 
						for(j in 1:p) is.rep[j] <- (sum(VectorListCannot[,j])==0)
					}
					p1 <- sum(is.rep) 

					Sep <- SepUpperBound
					if(p1>0) Sep <- min(c( Sep, sort(SepEach)[p+1-p1] ))
					YesAlready <- 0
					VectorListNo <- VectorListCannot[,!is.rep,drop=FALSE]
					for(i in dim(VectorListNo)[1]:1) if(sum(VectorListNo[i,])==0) VectorListNo <- VectorListNo[-i,,drop=FALSE]
					VectorListYes <- matrix(0,0,p-p1)
					VectorGeneratorAlt <- VectorListYes
					VectorListInfo <- VectorInConsideration
					for(j in p:1) if(is.rep[j]==1) VectorListInfo <- VectorListInfo[,-j,drop=FALSE]
					for(i in dim(VectorListInfo)[1]:1) if(sum(VectorListInfo[i,1:(p-p1)])==0) VectorListInfo <- VectorListInfo[-i,,drop=FALSE]
					while(dim(VectorListYes)[1]<2^(p-p1-p3)-1) {
						VectorToNoIndex <- order(VectorListInfo[,p-p1+1])[1]
						SepN <- min(VectorListInfo[,p-p1+1])
						VectorToNo <- VectorListInfo[VectorToNoIndex,1:(p-p1)]
						for(kkk in 1:(p-p1)) if(VectorToNo[kkk]>0) break
						if(kkk-1>0) for(i in 1:(kkk-1)) {
							TempVector <- VectorToNo
							TempVector[i] <- TempVector[i]+1
							TempVector <- c(TempVector,sqrt(sum((TempVector*wdsOrdered[!is.rep])^2)))
							VectorListInfo <- rbind(VectorListInfo,TempVector)
						}
						VectorListInfo <- VectorListInfo[-VectorToNoIndex,,drop=FALSE]
						NeedTo <- 1
						if(dim(VectorListYes)[1]>0) for(i in 1:dim(VectorListYes)[1]) if(sum(VectorToNo==VectorListYes[i,])==p-p1) NeedTo <- 0
						if(dim(VectorListYes)[1]==0) for(i in 1:dim(VectorListNo)[1]) if(sum(VectorToNo==VectorListNo[i,])==p-p1) NeedTo <- 0
						if(dim(VectorListYes)[1]>0) for(i in 1:dim(VectorListNo)[1]) for(j in 1:dim(VectorListYes)[1]) if(sum(VectorToNo==BAdd(VectorListNo[i,],VectorListYes[j,]))==p-p1) NeedTo <- 0
						if(NeedTo==1) {
							CanNo <- ( CheckArrange(p-p1,p3,rbind(VectorListNo,VectorToNo),VectorListYes) >0 )
							if(CanNo) VectorListNo <- rbind(VectorListNo,VectorToNo)
							if(!CanNo) { 
								VectorListYes <- AddToGroup(VectorListYes,VectorToNo)
								VectorGeneratorAlt <- rbind(VectorGeneratorAlt,VectorListYes)
								if(YesAlready==0) {
									Sep <- min(c( Sep, SepN ))
									YesAlready <- 1
								}
							}
						}
					}
					L01alt <- rbind(rep(0,p-p1),VectorListYes)  ## VectorGenerator and L01 not updated. 
					m <- m_exact_fromL01(L01alt,is.rep,s_vector[DimOrder])
					if(m>=n) {  
						SepBest <- Sep
						mBest <- m
						p3Best <- p3
						p1Best <- p1
						s_vectorBest <- s_vector
						NCannotBest <- dim(VectorListCannot)[1]
						DimOrderBest <- DimOrder
						is.repBest <- is.rep
						L01altBest <- L01alt
						VectorGeneratorAltBest <- VectorGeneratorAlt
						if(p1==0) IsChance <- 0
					}
					if(m<n) { s_vector[p]<-s_vector[p]+1; break; }
				}
			}

			if(IsChance==0) {
				for(kk in p:1) if(s_vector[kk]>2) break; 
				if(kk==1) { tobreak<-1; break; }
				s_vector[kk-1] <- s_vector[kk-1]+1 
				while(kk>1) {   # Check if the kk-1 dimension has no chance
					if(weight[kk-1]/(s_vector[kk-1]-1)*2>=SepBest) break
					kk <- kk-1; 
					if(kk>1) s_vector[kk-1] <- s_vector[kk-1]+1 
				}
				if(kk==1) { tobreak<-1; break; }
				s_vector[kk:p] <- 2
				s_vector[p] = max(c( ceiling(n*2^p3/prod(ceiling(s_vector[-p]/2)*2) /2)*2-1, 2 ))  
			}
		}
	}
	
## Generate the best design	

	VectorGeneratorBest <- matrix(0,dim(VectorGeneratorAltBest)[1],p)
	VectorGeneratorBest[,DimOrderBest[is.repBest==0]] <- VectorGeneratorAltBest
	for(j in p:1) if(is.repBest[j]==1) { temp<-rep(0,p); temp[DimOrderBest[j]]<-1; VectorGeneratorBest <- rbind(temp,VectorGeneratorBest); }
	L01Best <- matrix(0,dim(L01altBest)[1],p)
	L01Best[,DimOrderBest[is.repBest==0]] <- L01altBest
	Design <- L01Best
	L01BestS <- L01Best
	for(j in p:1) {
		for(jj in 1:p) if(DimOrderBest[jj]==j) break
		if(is.repBest[jj]==0) if(s_vectorBest[j]>3) for(i in 2:floor(s_vectorBest[j]/2)) {
			DesignAdd <- L01Best 
			DesignAdd[,j] <- DesignAdd[,j]+(i-1)*2
			Design <- rbind(Design,DesignAdd)
		}
		if(is.repBest[jj]==0) if(floor(s_vectorBest[j]/2)*2 != s_vectorBest[j]) {
			DesignAdd <- L01Best[L01Best[,j]==0,,drop=FALSE] 
			DesignAdd[,j] <- DesignAdd[,j]+s_vectorBest[j]-1
			Design <- rbind(Design,DesignAdd)
		}	
		if(is.repBest[jj]==1) for(i in 2:s_vectorBest[j]) {
			DesignAdd <- L01Best 
			DesignAdd[,j] <- DesignAdd[,j]+(i-1)
			Design <- rbind(Design,DesignAdd)
		}
		L01Best <- Design
	}
	for(j in 1:p) Design[,j] <- Design[,j]/(s_vectorBest[j]-1)
	DesignTransformed <- Design
	for(j in 1:p) DesignTransformed[,j] <- DesignTransformed[,j]*weight[j]

	return(list(Design=Design,SeparationDistance=SepBest,m=mBest,DesignTransformed=DesignTransformed,weight=weight,s_vector=s_vectorBest,L01=L01BestS))
}


InterleavedMaximinDAlg3 <- function(p,n,weight=rep(1,p)) {

	if(!p>=9|!p==round(p)) stop("p must be an integer greater than eight.") 
	if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.") 
	if(!length(weight)==p) stop("weight must be a vector with length p.") 

	SepBest <- 0
	WeightOrder <- order(-weight)
	weightOrdered <- weight[WeightOrder]
	Dpart <- InterleavedMaximinDAlg2(8,n,weightOrdered[1:8])
	s_vectorOrdered <- c(Dpart$s_vector,rep(2,p-8))
	L01 <- Dpart$L01
	for(j in 8:1) if(max(L01[,j])==0) { L01Add <- L01; L01Add[,j] <- 1; L01 <- rbind(L01,L01Add); }

	for(pp in 8:(p-1)){
		L01List <- cbind(L01,rep(-1,dim(L01)[1]),rep(0,dim(L01)[1]))
		for(j in 1:pp) L01List[,pp+2] <- L01List[,pp+2] + (L01List[,j]*weightOrdered[j]/(s_vectorOrdered[j]-1))^2
		L01List[,pp+2] <- sqrt(L01List[,pp+2])
		VectorListYes <- matrix(0,0,pp)
		VectorOrder <- order(L01List[,pp+2])
		L01List[VectorOrder[1],pp+1] <- 0
		L01List[VectorOrder[2],pp+1] <- 1
		if(dim(L01List)[1]>2) for(i in 3:dim(L01List)[1]) if(L01List[VectorOrder[i],pp+1]==-1) {
			L01List[VectorOrder[i],pp+1] <- 1
			VectorListYes <- AddToGroup(VectorListYes,BAdd(L01List[VectorOrder[i],1:pp],L01List[VectorOrder[2],1:pp]))
			for(ii in 1:dim(L01List)[1]) for(iii in 1:dim(VectorListYes)[1]) 
				if(sum(L01List[ii,1:pp]==VectorListYes[iii,])==pp) 
					L01List[ii,pp+1] <- 0
		}
		L01 <- L01List[,1:(pp+1)]
	}
	
	L01BestS <- L01
	p1 <- 0
	for(j in 1:p) {
		for(i in 2:dim(L01BestS)[1]) if(sum(L01BestS[i,-j])==0) {
			p1 <- p1+1
			L01BestS <- L01BestS[L01BestS[,j]==0,,drop=FALSE]
			break
		}
	}
	
	L01List[,p+1] <- 0
	for(j in 1:p) L01List[,p+1] <- L01List[,p+1] + (L01List[,j]*weightOrdered[j]/(s_vectorOrdered[j]-1))^2
	L01List[,p+1] <- sqrt(L01List[,p+1])
	mBest <- Dpart$m
	L01Best <- matrix(0,dim(L01)[1],p) 
	for(j in 1:p) L01Best[,WeightOrder[j]] <- L01[,j] 
	s_vectorBest <- rep(0,p)
	for(j in 1:p) s_vectorBest[WeightOrder[j]] <- s_vectorOrdered[j]
	Design <- L01Best
	for(j in p:1) {
		if(s_vectorBest[j]>3) for(i in 2:floor(s_vectorBest[j]/2)) {
			DesignAdd <- L01Best 
			DesignAdd[,j] <- DesignAdd[,j]+(i-1)*2
			Design <- rbind(Design,DesignAdd)
		}
		if(floor(s_vectorBest[j]/2)*2 != s_vectorBest[j]) {
			DesignAdd <- L01Best[L01Best[,j]==0,,drop=FALSE] 
			DesignAdd[,j] <- DesignAdd[,j]+s_vectorBest[j]-1
			Design <- rbind(Design,DesignAdd)
		}	
		L01Best <- Design
	}
	for(j in 1:p) Design[,j] <- Design[,j]/(s_vectorBest[j]-1)
	DesignTransformed <- Design
	for(j in 1:p) DesignTransformed[,j] <- DesignTransformed[,j]*weight[j]

	SepBest <- min(L01List[-1,p+1])
	if(max(s_vectorBest)>2) SepBest <- min(c(SepBest,2*weight[s_vectorBest>2]/(s_vectorBest[s_vectorBest>2]-1)))

	return(list(Design=Design,SeparationDistance=SepBest,m=mBest,DesignTransformed=DesignTransformed,weight=weight,s_vector=s_vectorBest,L01=L01BestS))
}


InterleavedMaximinD <- function(p,n,weight=rep(1,p)) {

	if(!p>=2|!p==round(p)) stop("p must be an integer greater than one and no greater than five.") 
	if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.") 
	if(!length(weight)==p) stop("weight must be a vector with length p.") 

	if(p<=5) return(InterleavedMaximinDAlg1(p,n,weight))
	if(p<=8) return(InterleavedMaximinDAlg2(p,n,weight))
	return(InterleavedMaximinDAlg3(p,n,weight))
}
