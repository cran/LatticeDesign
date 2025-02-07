## Densest packing-based maximum projection designs
## Author: Xu He, hexu@amss.ac.cn


ProjSepD <- function(design) {
	.C(ProjSep, PACKAGE="LatticeDesign", 
	as.double(design), as.integer(dim(design)[2]), as.integer(dim(design)[1]), as.double(rep(0,dim(design)[2])) )[[4]] 
}


DPMPD <- function(p=2, n, rotation="magic", w=100){

	if(!p>=2|!p<=8|!p==round(p)) stop("p must be an integer greater than one and no greater than eight.") 
	if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.") 
	if(!w>=1|!w==round(w)|!w<=10000) stop("w must be a positive integer greater no greater than 10000.") 

	#### Generate the unrotated candidate sets. 

	detG <- 2
	if(p==2) detG <- 2*sqrt(3)
	if(p==6) detG <- sqrt(3)
	if(p==7) detG <- sqrt(2)
	if(p==8) detG <- 1
	r <- (n*detG)^(1/p)/2*sqrt(p)+1

	ss <- floor( r )
	FE <- matrix(0,ss+1,2)
	FE[,2] = 0:ss
	FE[,1] = FE[,2]^2
	for(j in 2:p) {  # Generate all nonnegative integer points with radius <= r
		nFE <- matrix(0,0,j+1)
		for(k in 0:ss) { oFE = cbind(FE,rep(k,dim(FE)[1])); oFE[,1] = oFE[,1]+k^2; nFE <- rbind(nFE,oFE); }
		FE <- nFE[nFE[,1]<= floor( r^2 ), ]
	} 
	FE[,1]=0; for(j in 2:(p+1)) FE[,1]=FE[,1]+FE[,j]; FE = FE[ floor(FE[,1]/2)*2==FE[,1] ,]; ## Delete half
	
	if(p>=6) {  ## Add half while reflecting to negative first dimension
		FE1 = FE[ FE[,2]>0, ]
		FE1[,2] = -FE1[,2]
		FE2 = FE[ FE[,2]==1, ]
		FE2[,2] = -FE2[,2]
		FE3 = rbind(FE,FE2)
		FE4 = FE3 +0.5
		FE5 = FE +0.5
		FE5[,2] = -FE5[,2]-1
		FE <- rbind(FE,FE1,FE4,FE5)
	}
	if(p<6) {  ## Reflect to negative first dimension
		nFE = FE[FE[,2]>0,]
		nFE[,2] <- -nFE[,2]
		FE <- rbind(FE,nFE)
	}
	
	if(p==2) { FE[,3] <- FE[,3]*sqrt(3); }  ## Rescale the last dimension
	if(p==6) { FE[,7] <- FE[,7]*sqrt(3); }
	if(p==7) { FE[,8] <- FE[,8]*sqrt(2); }
	
	FE[,1] <- 0; 
	for(j in 2:(p+1)) FE[,1] <- FE[,1] +FE[,j]^2
	FE <- FE[FE[,1]<=r^2,] 
	FE <- FE[,-1]

	for(j in 2:p) {  ## Reflect to negative second-to-pth dimensions 
		nFE = FE[FE[,j]>0,]
		nFE[,j] <- -nFE[,j]
		nFE[,1] <- -nFE[,1]
		FE <- rbind(FE,nFE)
	}


	## Prepare rotation for p=4 and 8

	mv = 25
	if(w>100) mv = 55
	vlist <- matrix(0,0,3)
	for(v1 in 1:mv) for(v2 in v1:mv) {
		vlist = rbind(vlist,c(v1,v2,v1^2+v2^2))
	}

	vvlist <- matrix(0,0,6)
	for(i in 1:dim(vlist)[1]) for(j in 1:dim(vlist)[1]) {
		ratio = vlist[j,3]/vlist[i,3]
		if(ratio>1) if(floor(sqrt(ratio))!=sqrt(ratio)) if(floor(ratio)*vlist[i,3]==vlist[j,3]) vvlist = rbind(vvlist,c(vlist[i,1:2],vlist[j,1:2],ratio,0)) 
	}
	vvlist[,6] = abs( (vvlist[,3]-vvlist[,1]*sqrt(vvlist[,5])) / sqrt( (vvlist[,3]-vvlist[,1]*sqrt(vvlist[,5]))^2 + (vvlist[,4]-vvlist[,2]*sqrt(vvlist[,5]))^2 ) )
	j=2
	while(j<=dim(vvlist)[1]) {
		for(k in 1:j) {
			if(k<j) if(abs(vvlist[j,6]-vvlist[k,6])<10^-10) { vvlist = vvlist[-j,]; break; }
		}
		if(k==j) j=j+1
	}
	for(i in 1:dim(vvlist)[1]) vvlist[i,6] = max(vvlist[i,1:4])
	
	vvlist1 = vvlist[vvlist[,5]==2,]
	vvlist2 = vvlist[vvlist[,5]==5,]
	vvlist3 = vvlist[vvlist[,5]==13,]

	if(w<=100) {
		vvlist1 = vvlist1[ order(vvlist1[,6])[1:10], ]
		vvlist2 = vvlist2[ order(vvlist2[,6])[1:10], ]
		vvlist3 = vvlist3[ order(vvlist3[,6])[1:10], ]
	}
	


	#### Try w rotation matrices
	
	ress = matrix(0,w,p+1)  # Minimum of separation distance of 1,...,p dimensional projections, overall score. 
	l <- (n*detG)^(1/p)
	maxscore <- -10^10
	
	for(ll in 1:w) {
		
		#### Generate rotation matrix 
		
		if(rotation!="magic"|p==5|p==7) {  ## Random rotation
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
		}
		
		if(rotation=="magic"&p==2) {  ## Magic rotation matrices
			q1 <- 3
			bb1 <- 20; bb2 <- 4; 
			if(w>100) {  bb1 <- floor(w^(2/3)); bb2 <- floor(w^(1/3)); }
			u1 = ceiling(ll/bb1)
			u2 = ceiling( ( ll - (u1-1)*bb1 ) / bb2 )
			u3 = ll - (u1-1)*bb1 - (u2-1)*bb2 -floor(bb2/2)
			if(u3^2==3*u1^2+u2^2) u3 = bb2 -floor(bb2/2)+1
			R = rbind( c(u1*sqrt(q1)+u3,-u2), c(u2,u1*sqrt(q1)+u3) )
			for(j in 1:2) R[,j] = R[,j]/sqrt(sum(R[,j]^2))
			R1 = rbind( c(sqrt(3)+1,sqrt(3)-1), c(-sqrt(3)+1,sqrt(3)+1) ) /2/sqrt(2) 
			R = R1 %*% R
		}
		
		if(rotation=="magic"&p==3) {  ## Magic rotation matrices
			qlist = 2:(w+10)
			for(k in 2:20) qlist = qlist[ floor(qlist/k^3)*k^3!=qlist ]
			bb1 <- 10; bb2 <- 55; 
			if(w>100) bb2 <- floor((w+10)/2)
			q1 = qlist[ll]; 
			if(ll>bb1) q1 = qlist[ll-bb1]/8
			if(ll>bb2) q1 = qlist[ll-bb2]/27
			R3=matrix(0,3,3)
			R3[1,1]=R3[3,2]=R3[2,3] = 1-q1
			R3[2,1]=R3[1,2]=R3[3,3] = q1-q1^(1/3)
			R3[3,1]=R3[2,2]=R3[1,3] = (1-q1)*(q1^(1/3)+q1^(2/3))
			R <- R3 / sqrt(sum(R3[,3]*R3[,3]))
		}
		
		if(rotation=="magic"&p==4) {  ## Magic rotation matrices
			bb1 <- 10; 
			if(w>100) bb1 <- floor(sqrt(w))
			if(w>bb1*76) bb1 <- ceiling(w/76)
			rot1 = vvlist1[ceiling(ll/bb1),]
			GRR1 = rbind( rot1[c(3,1)], rot1[c(4,2)] )
			q1 = rot1[5]
			R1 <- GRR1 %*% rbind( c(1,1), c(-sqrt(q1),sqrt(q1)) )
			for(j in 1:2) R1[,j] = R1[,j]/sqrt(sum(R1[,j]^2))
			
			rot2 = vvlist2[ll-ceiling(ll/bb1)*bb1+bb1,]
			GRR2 = rbind( rot2[c(3,1)], rot2[c(4,2)] )
			q2 = rot2[5]
			R2 <- GRR2 %*% rbind( c(1,1), c(-sqrt(q2),sqrt(q2)) )
			for(j in 1:2) R2[,j] = R2[,j]/sqrt(sum(R2[,j]^2))
			
			R <- rbind( cbind(R2[1,1]*R1,R2[1,2]*R1), cbind(R2[2,1]*R1,R2[2,2]*R1) )
		}
		
		if(rotation=="magic"&p==6) {  ## Magic rotation matrices
			bb1 <- 50; bb2 <- 25; 
			if(w>100)  { bb2 <- floor((w/5)^(1/3)*5); bb1 <- floor((w/5)^(2/3)*5); }
			q1 <- 3
			u1 = ceiling(ll/bb1)
			u2 = ceiling( ( ll - (ceiling(ll/bb1)-1)*bb1 ) / bb2 )
			u3 = ceiling( ( ll - (ceiling(ll/bb2)-1)*bb2 ) / 5 ) - ceiling(bb2/10)
			if(u3==0) u3 = ceiling(bb2/5) - ceiling(bb2/10) + 1
			if(w<=100) if(u1==2&u2==2&abs(u3)==2) u1 = 3
			R2 = rbind( c(u1*sqrt(q1)+u3,-u2), c(u2,u1*sqrt(q1)+u3) )
			for(j in 1:2) R2[,j] = R2[,j]/sqrt(sum(R2[,j]^2))

			qlist = 2:100
			for(k in 2:4) qlist = qlist[ floor(qlist/k^3)*k^3!=qlist ]
			q1 = qlist[ ll-(ceiling(ll/5)-1)*5 ]; 
			R3=matrix(0,3,3)
			R3[1,1]=R3[3,2]=R3[2,3] = 1-q1
			R3[2,1]=R3[1,2]=R3[3,3] = q1-q1^(1/3)
			R3[3,1]=R3[2,2]=R3[1,3] = (1-q1)*(q1^(1/3)+q1^(2/3))
			R3 <- R3 / sqrt(sum(R3[,3]*R3[,3]))
			
			R <- rbind( cbind(R3[1,1]*R2,R3[1,2]*R2,R3[1,3]*R2), cbind(R3[2,1]*R2,R3[2,2]*R2,R3[2,3]*R2), cbind(R3[3,1]*R2,R3[3,2]*R2,R3[3,3]*R2) )
		}
		
		if(rotation=="magic"&p==8) {  ## Magic rotation matrices
			bb2 <- floor(w^(1/3)); bb1 <- floor(w^(2/3)); 

			rot1 = vvlist1[ ceiling(ll/bb1) ,]
			GRR1 = rbind( rot1[c(3,1)], rot1[c(4,2)] )
			q1 = rot1[5]
			R1 <- GRR1 %*% rbind( c(1,1), c(-sqrt(q1),sqrt(q1)) )
			for(j in 1:2) R1[,j] = R1[,j]/sqrt(sum(R1[,j]^2))
			
			rot2 = vvlist2[ ceiling( (ll-ceiling(ll/bb1)*bb1+bb1) / bb2 ) ,]
			GRR2 = rbind( rot2[c(3,1)], rot2[c(4,2)] )
			q2 = rot2[5]
			R2 <- GRR2 %*% rbind( c(1,1), c(-sqrt(q2),sqrt(q2)) )
			for(j in 1:2) R2[,j] = R2[,j]/sqrt(sum(R2[,j]^2))

			rot3 = vvlist3[ ll-ceiling(ll/bb2)*bb2+bb2 ,]
			GRR3 = rbind( rot3[c(3,1)], rot3[c(4,2)] )
			q3 = rot3[5]
			R3 <- GRR3 %*% rbind( c(1,1), c(-sqrt(q3),sqrt(q3)) )
			for(j in 1:2) R3[,j] = R3[,j]/sqrt(sum(R3[,j]^2))
			
			R <- rbind( cbind(R2[1,1]*R1,R2[1,2]*R1), cbind(R2[2,1]*R1,R2[2,2]*R1) )
			R <- rbind( cbind(R3[1,1]*R,R3[1,2]*R), cbind(R3[2,1]*R,R3[2,2]*R) )
		}
		
		#### Find delta and thus epsilon and design

		E <- FE %*% R 

		isE <- FALSE
		isB <- FALSE
		isL <- FALSE
		epsilonV <- matrix(0,3,p)
		for(i in 1:100) {
			while(1) {
				delta <- runif(p,-1,1) 
				if( sum( sort(abs(delta))[(p-1):p] ) > 1) next
				if(p==2) delta[p] <- delta[p]*sqrt(3)
				if(p==6) delta[p] <- delta[p]*sqrt(3)
				if(p==7) delta[p] <- delta[p]*sqrt(2)
				if( sum( (abs(delta)+l/2)^2 ) <= r^2) break
			}
			epsilon <- delta %*% R 
			design <- E
			for(j in 1:p) { if(length(design)<=p) break; design <- design[ abs(design[,j]-epsilon[j])<l/2 ,]; }
			if(length(design)<=p) next; 
			if(dim(design)[1]==n) { isE <- TRUE; epsilonV[1,] <- epsilon; break; }
			if(isL==FALSE & dim(design)[1]<n) { isL <- TRUE; epsilonV[2,] <- epsilon; }
			if(isB==FALSE & dim(design)[1]>n) { isB <- TRUE; epsilonV[3,] <- epsilon; }
			if(isL==TRUE & isB==TRUE) break;
		}
		
		if(isE==FALSE) for(j in 1:(p-1)) {
			epsilon <- c(epsilonV[2,1:j],epsilonV[3,(j+1):p])
			design <- E
			for(j in 1:p) { if(length(design)<=p) break;  design <- design[ abs(design[,j]-epsilon[j])<l/2 ,]; }
			if(length(design)<=p) next; 
			if(dim(design)[1]==n) { isE <- TRUE; epsilonV[1,] <- epsilon; break; }
			if(dim(design)[1]<n) { epsilonV[2,] <- epsilon; break; }
			if(dim(design)[1]>n) { epsilonV[3,] <- epsilon; }
		}

		if(isE==FALSE) for(k in 1:100) { 
			epsilon <- (epsilonV[2,]+epsilonV[3,])/2
			design <- E
			for(j in 1:p) { if(length(design)<=p) break; design <- design[ abs(design[,j]-epsilon[j])<l/2 ,]; }
			if(length(design)<=p) next; 
			if(dim(design)[1]==n) { isE <- TRUE; epsilonV[1,] <- epsilon; break; }
			if(dim(design)[1]<n) { epsilonV[2,] <- epsilon; }
			if(dim(design)[1]>n) { epsilonV[3,] <- epsilon; }
		}
		
		for(j in 1:p) design[,j] <- (design[,j]-epsilon[j]) / l + .5
		
		
		#### Evaluate design 
		
		ress[ll,1:p] = ProjSepD(design)
		ress[ll,p+1] = sum( log(ress[ll,1:p]*n^(1/(1:p))) * c(1/2,2:p) )
		if(ll==1|ress[ll,p+1]>maxscore) {
			maxscore <- ress[ll,p+1]
			designBest <- design
		}

	}
	
	for(j in 1:p) designBest[,j] = designBest[,j]-(min(designBest[,j])+max(designBest[,j])-1)/2
	return(list(Design=designBest,ProjectiveSeparationDistance=sqrt(ProjSepD(designBest))))
}


