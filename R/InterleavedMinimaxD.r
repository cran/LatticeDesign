
LRS <- function(numerator,denominator) {  # A program to compute the vertexes of a polytope 
	n = dim(numerator)[1]  # number of constraints
	p = dim(numerator)[2]  # number of dimensions plus one
	OutMv = rep(-1,p) 
	MaxOut = 1024 * 4 * (p-1)  # Maximum length for output. A better formula shall be 2^(p-1) * 32
	
	while(1) {
		Out = rep(-1.5,MaxOut)
		res = .C(lrs_XuHe, PACKAGE="LatticeDesign", 
		as.integer(n), as.integer(p), as.integer(numerator), as.integer(denominator), 
		  as.double(Out), as.double(OutMv), as.integer(MaxOut)) 
		if(res[[7]] > 0) break 
		MaxOut = MaxOut * 4
	}
	Vertexes = matrix(res[[5]][1:res[[7]]], ncol = p-1, byrow = TRUE) 
	list(Radius=res[[6]][1],MaxValue=res[[6]][2:p],Vertexes=Vertexes)
}


selects <- function(p,n,ps,expair,coefC,maxdissimilarity,minndp) {  # A function to find eligible s_vector's
	maxnrow = floor( 4^(log(n)/log(10)) * 100 )  # Maximum length for output. 1024 should be enough for cases with n<=100. 
	while(1){
		lssst = matrix(0,p+2,maxnrow)
		res = .C(selectlss, PACKAGE="LatticeDesign", 
		as.integer(p), as.integer(n), as.integer(ps), as.integer(dim(expair)[1]), as.integer(expair), as.integer(coefC), 
		  as.double(maxdissimilarity), as.double(minndp), as.integer(lssst), as.integer(maxnrow)) 
		if(res[[10]]<maxnrow-1) break
		maxnrow = maxnrow * 4
	}
	if(res[[10]]>0) return( matrix( res[[9]][1:(res[[10]]*(p+2))], ncol = p+2, byrow = TRUE ) )
	if(res[[10]]==0) return( matrix( 0, 0, p+2 ) )
}


gcd = function(a,b) {  # greatest common divider of two integers. 
	a = ceiling(abs(a))
	b = ceiling(abs(b))
	while(b!=0) {
		r = a %% b
		a = b
		b = r
	}
	return(a) 
}


setGmML = setClass("GmML", slots = c(G="matrix", ps="numeric", coefD="matrix", coefC="matrix", Evor="matrix", expair="matrix", coin="matrix", is.rep="numeric"))
  ## A class for information on interleaved lattices. 
  ## coefD: the 2^(ps[2]+ps[3]) * (pp+1) matrix. The first pp columns are the product points. The last column gives the coset index. 
  #### Defined if ps[2]+ps[3]>0. 
  ## coefC: the matrix to record the number of points for each coset.  
  #### Defined if ps[2]+ps[3]>0. 
  #### The 2^ps[3] columns for cosets. 
  #### The 2^(ps[2]+ps[3]) rows for scenarios. 
  #### The first row for all dimensions being 0 or 1. 
  #### The 2,4,6,... rows for the last dimension being 0. 
  #### The 3,4,7,8,... rows for the second last dimension being 0. 
  #### Minus one an after decomposition, the elements with one shall be 0 for corresponding dimensions. 
  ## expair: the pairs of dimensions that are exchangeable, indexing on the alt dimensions. 
  #### Defined if ps[2]+ps[3]>1. 
  ## is.rep: 1 for repetitive dimensions and 0 for alternative dimensions. 
  ## coin: 1 for coincident dimension pairs

initG1 = function(pp) {  # Diagonal generator matrix with ones. 
	GmML2 = setGmML(G=diag(pp))
	GmML2@ps = c(pp,0,0)

	GmML2
}

initG2 = function(pp) {  # Diagonal generator matrix with twos. 
	GmML2 = setGmML(G=diag(pp)*2)
	GmML2@ps = c(0,0,pp)

	GmML2@coefD = matrix(0,2^pp,pp+1) 
	is = rep(0,pp)
	for(i in 1:(2^pp)) {
		GmML2@coefD[i,] = c(is,i-1) 
		is[pp]=is[pp]+1
		if(pp>1) for(j in pp:2) if(is[j]==2) { is[j]=0; is[j-1]=is[j-1]+1; }
	}

	GmML2
}

inducG1 = function(GmML1) {  # Generator matrix induced with one more replicative dimension.
	GmML2 = GmML1
	pp = dim(GmML1@G)[1]+1
	GmML2@G = matrix(0,pp,pp)
	GmML2@G[1,1]=1
	GmML2@G[2:pp,2:pp]=GmML1@G
	GmML2@ps[1]=GmML1@ps[1]+1
	
	GmML2
}

inducG = function(GmML1,addv) {  # Generator matrix induced with one more alternative dimension.
	GmML2 = GmML1
	pp = dim(GmML1@G)[1]+1
	GmML2@G = matrix(0,pp,pp)
	GmML2@G[1,]=addv
	GmML2@G[2:pp,2:pp]=GmML1@G
	GmML2@ps[2]=GmML1@ps[2]+1
	
	coefD1=cbind(rep(0,2^(GmML1@ps[2]+GmML1@ps[3])),GmML1@coefD)
	coefD2=coefD1
	coefD2[,1]=1
	for(j in 2:pp) if(addv[j]==1) coefD2[,j]=1-coefD2[,j]
	GmML2@coefD=rbind(coefD1,coefD2)

	GmML2
}

newG = function(H) {  # Arbitrary generator matrix. Used for ps[2]>=2. 
	pp = dim(H)[1]
	ps = rep(0,3)
	for(i in 1:pp) if(prod(H[i,]==diag(pp)[i,])) ps[1]=ps[1]+1
	for(i in pp:1) if(prod(H[i,]==diag(pp)[i,]*2)) ps[3]=ps[3]+1
	ps[2]=pp-ps[1]-ps[3]
	if(ps[1]==pp) GmML2 = initG1(pp) 
	if(ps[1]>0&ps[1]<pp) {
		GmML1 = newG(H[2:pp,2:pp,drop=FALSE])
		GmML2 = inducG1(GmML1)
	}
	if(ps[1]==0&ps[2]>0) {
		GmML1 = newG(H[2:pp,2:pp,drop=FALSE])
		GmML2 = inducG(GmML1,H[1,])
	}
	if(ps[1]==0&ps[2]==0) {
		GmML2 = initG2(pp)
	}
	
	GmML2
}

supp = function(GmML1) {  # Supplement generating information.  
	pp=dim(GmML1@G)[1]

	if(GmML1@ps[3]>0) {
		GmML1@coefC=matrix(0,2^(GmML1@ps[2]+GmML1@ps[3]),2^GmML1@ps[3])
		for(i in 1:(2^(GmML1@ps[2]+GmML1@ps[3]))) {
			corres=0
			for(j in (GmML1@ps[2]+GmML1@ps[3]):1) { if(GmML1@coefD[i,j]==0) corres=c(corres,corres+2^(GmML1@ps[2]+GmML1@ps[3]-j)); }
			GmML1@coefC[ corres+1, GmML1@coefD[i,GmML1@ps[2]+GmML1@ps[3]+1]+1 ] = GmML1@coefC[ corres+1, GmML1@coefD[i,GmML1@ps[2]+GmML1@ps[3]+1]+1 ]+1
		}
	}
	
	expair = matrix(0,0,2)
	if(GmML1@ps[2]==1) {
		H = GmML1@G
		for(j1 in 1:(GmML1@ps[2]+GmML1@ps[3]-1)) for(j2 in (j1+1):(GmML1@ps[2]+GmML1@ps[3])) if(H[GmML1@ps[1]+1,GmML1@ps[1]+j1]==H[GmML1@ps[1]+1,GmML1@ps[1]+j2]) expair=rbind(expair,c(j1,j2))
	}
	if(GmML1@ps[2]>1) {
		H = GmML1@G[(GmML1@ps[1]+1):(GmML1@ps[1]+GmML1@ps[2]),(GmML1@ps[1]+1):pp]
		for(j1 in 1:(dim(H)[2]-1)) for(j2 in (j1+1):dim(H)[2]) {
			if(j2<=dim(H)[1]) if(prod(H[j1,(dim(H)[1]+1):dim(H)[2]]==H[j2,(dim(H)[1]+1):dim(H)[2]])==1) expair=rbind(expair,c(j1,j2))
			if(j1 >dim(H)[1]) if(prod(H[,j1]==H[,j2])==1) expair=rbind(expair,c(j1,j2))
			if(j1<=dim(H)[1]&j2>dim(H)[1]) if(H[j1,j2]==1) if(sum(H[j1,])==2|sum(H[,j2])==1) expair=rbind(expair,c(j1,j2))
		}
	}
	GmML1@expair <- expair

	Evor = diag(pp)
	if(GmML1@ps[2]+GmML1@ps[3]>0) for(j in (GmML1@ps[1]+1):pp) Evor[j,j]=2
	if(GmML1@ps[3]>0) {
		D0 = GmML1@coefD[ GmML1@coefD[,GmML1@ps[2]+GmML1@ps[3]+1]==0 , 1:(GmML1@ps[2]+GmML1@ps[3]) ] 
		D0 = D0[-1,,drop=FALSE]
		j1=1 
		while(1) {
			if(j1>dim(D0)[1]) break
			todelete=FALSE
			for(j2 in 1:dim(D0)[1]) if(j2!=j1) if(prod((D0[j1,]-D0[j2,])>=0)) todelete=TRUE
			if(todelete==TRUE) D0=D0[-j1,] else j1=j1+1
		}
		if(GmML1@ps[1]>0) Evor = rbind(Evor, cbind(matrix(0,dim(D0)[1],GmML1@ps[1]),D0)) else Evor = rbind(Evor,D0)
	}
	GmML1@Evor <- Evor
	
	coin = matrix(1,pp,pp)  # 1: the columns are coincident
	GM1 = GmML1@G[1:(pp-GmML1@ps[3]),,drop=FALSE]
	for(i1 in 1:(pp-1)) for(i2 in (i1+1):pp) { 
		for(j in 1:dim(GM1)[1])  if(GM1[j,i1]!=GM1[j,i2])  {  coin[i1,i2]=coin[i2,i1]=0;  break;  }
	}
	for(i1 in 1:pp)  coin[i1,i1]=0 
	GmML1@coin = coin

	is.rep = rep(1,pp) 
	for(j in 1:pp) if(sum(GmML1@G[j,])>1) is.rep[j] = 0
	GmML1@is.rep = is.rep

	GmML1
}


InterleavedMinimaxD = function(p,n,maxdissimilarity=2*p)
{
	if(!p>=2|!p<=8|!p==round(p)) stop("p must be an integer greater than one and no greater than eight.") 
	if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.") 

	GMs <- NULL
	if(p==2) GsName <- data(GMs2, envir=environment())
	if(p==3) GsName <- data(GMs3, envir=environment())
	if(p==4) GsName <- data(GMs4, envir=environment())
	if(p==5) GsName <- data(GMs5, envir=environment())
	if(p==6) GsName <- data(GMs6, envir=environment())
	if(p==7) GsName <- data(GMs7, envir=environment())
	if(p==8) GsName <- data(GMs8, envir=environment())
#	GMs <- get(GsName)
	CRofP=p  # The best covering radius of the transformed design so far for the given interleaved lattice. 
	
	for(indexG in 1:(dim(GMs)[1]/p)){  # Consider one generator matrix
		GM = matrix(0,p,p)
		for(j in 1:p) GM[,j]=GMs[((indexG-1)*p+1):(indexG*p),j]
		theGmML = newG(GM)
		ps = theGmML@ps
		theGmML = supp(theGmML)

		alllss = selects(p,n,ps,theGmML@expair,theGmML@coefC,maxdissimilarity,0)  # Find all eligible s_vectors
		if(dim(alllss)[1]>0) for(choicelss in 1:dim(alllss)[1]) {
			lss=alllss[choicelss,1:p]
			nd=alllss[choicelss,p+1]
			thecoef=alllss[choicelss,p+2]+1
			if(max(lss)>=2^15)  next;  # have to add this patch to avoid error. Big integer causes problems. 
			
			TCR = lss
			if(ps[2]+ps[3]>0) TCR[(ps[1]+1):p] = lss[(ps[1]+1):p]-1
			Anum = theGmML@Evor
			Aden = matrix(1,dim(Anum)[1],dim(Anum)[2]) 
			for(k in 1:p) Aden[,k] = TCR[k] 
			bnum = Anum[,1]^2
			bden = Aden[,1]^2*2
			gcds = rep(0,length(bnum))
			tobreak = FALSE
			for(k in 2:p) {
				for(i in 1:length(bnum))  gcds[i] = gcd(bden[i],Aden[i,k]^2*2)
				snum = bnum*(Aden[,k]^2*2/gcds) + Anum[,k]^2*(bden/gcds)
				sden = (bden/gcds)*Aden[,k]^2*2
				if(max(c(snum,sden))>2^30) { tobreak = TRUE; break; }  # make sure later computed numbers are below 2^31-1. 
				for(i in 1:length(bnum))  gcds[i] = gcd(snum[i],sden[i])
				bnum = snum/gcds
				bden = sden/gcds
				if(max(c(bnum,bden))>2^20) { tobreak = TRUE; break; }  # make sure later computed numbers are below 2^31-1. 
			}
			if(tobreak==TRUE)  next;  # have to add this patch to avoid error. Big integer causes problems. 
			num = cbind(bnum,-Anum)
			den = cbind(bden,Aden)
			num = rbind(num, cbind(matrix(0,p,1),diag(p)) )
			den = rbind(den, matrix(1,p,p+1) )
			AA = LRS(num,den)
			CRofD = AA$Radius
			VorMr = AA$MaxValue
			VorV = AA$Vertexes
			TCRnew = lss+1 - 2*VorMr*TCR
			Sradiuses=rep(0,dim(VorV)[1])
			for(i in 1:length(Sradiuses)) Sradiuses[i]=sum((VorV[i,]*TCR/TCRnew)^2)
			CRofD = sqrt(max(Sradiuses))  # The best covering radius of the transformed design for the M and lss. 
			if(CRofD<CRofP) { 
				CRofP = CRofD; 
				theGmMLbest = theGmML; lssbest = lss; ndbest=nd; coefbest = thecoef; 
				TCRnewbest = TCRnew; TCRbest = TCR; VorMrbest=VorMr; 
			}
		}
	}

	theGmML = theGmMLbest
	ps = theGmML@ps
	if(ps[3]>0) {
		Ebase = theGmML@coefD[ theGmML@coefD[,ps[2]+ps[3]+1]==coefbest-1 , 1:(ps[2]+ps[3]) ]
		if(ps[1]>0) Ebase = cbind( matrix(0,dim(Ebase)[1],ps[1]) , Ebase )
	}
	if(ps[3]==0) Ebase = matrix(0,1,p) 

	Ethis = Ebase
	if(ps[3]>0) for(k in p:(ps[1]+1)) {
		Eadded = Ethis
		if(lssbest[k]>=4) for(ii in 1:(floor(lssbest[k]/2)-1)) {
			Eadd = Ethis
			Eadd[,k] = Ethis[,k]+2*ii
			Eadded = rbind(Eadded,Eadd)
		}
		if(lssbest[k]>=3) if(floor(lssbest[k]/2)*2<lssbest[k]) {
			Eadd = Ethis[Ethis[,k]==0,,drop=FALSE]
			Eadd[,k] = Eadd[,k]+2*floor(lssbest[k]/2)
			Eadded = rbind(Eadded,Eadd)
		}
		Ethis = Eadded
	}	
	if(ps[1]>0) for(k in ps[1]:1) {
		Eadded = Ethis
		if(lssbest[k]>1) for(ii in 1:(lssbest[k]-1)) {
			Eadd = Ethis
			Eadd[,k] = Ethis[,k]+ii
			Eadded = rbind(Eadded,Eadd)
		}
		Ethis = Eadded
	}

	Dbest = Ethis
	for(k in 1:p) Dbest[,k]=(Dbest[,k]+1)/TCRnewbest[k]-VorMrbest[k]*TCRbest[k]/TCRnewbest[k]; 

	return(list(Design=Dbest,TargetFillDistance=CRofP,ActualSize=ndbest,s_vector=lssbest,L01=Ebase))
}


