
tn <- function(lsss,coefC) {   # compute the number of points with regard to the alt dimensions. lsss gives s for the alt dimensions. 
	if(dim(coefC)[2]==0) return(prod(lsss))
	oneortwo = (lsss>2)
	counts = rep(0,dim(coefC)[2])
	if(sum(oneortwo)==0)  counts = coefC[ sum((2-lsss)*2^((length(lsss)-1):0))+1 ,]
	if(sum(oneortwo)!=0)  { 
		for(j in 1:length(lsss)) if(lsss[j]>2) break
		thedim = j 
		thetwos = floor((lsss[j]-1)/2)
		lsss1 = lsss
		lsss1[thedim]=2
		lsss2 = lsss
		lsss2[thedim]=lsss[thedim]-thetwos*2
		counts = thetwos*tn(lsss1,coefC) + tn(lsss2,coefC)
	}
	
	counts
}

RhoCcal <- function(sthis,coin) {  # Calculate RhoC, root averaged square of correlations
	sthiseven = ( floor(sthis/2)*2 == sthis )
	ss = sqrt(3/(sthis^2-1)) * sthiseven
	RhoC = ss %*% coin %*% ss /length(sthis) /(length(sthis)-1)
	return(RhoC)
}

RhoFcalappro <- function(A) {  # Calculate approximate fill distance with numeric inputs
	bestse = Inf
	mden = 100
	if(dim(A)[2]>=8) mden = floor(400/dim(A)[2])
	for(den in 1:mden) {
		se = sum(abs(round(A*den)-A*den))/den
		if(se<bestse) {
			Anum = round(A*den)
			bestden = den
			bestse = se
		}
	}
	for(k in dim(Anum)[2]:1) if(prod(Anum[,k]==0)==1) Anum = Anum[,-k,drop=FALSE]
	for(i in dim(Anum)[1]:1) if(prod(Anum[i,]==0)==1) Anum = Anum[-i,,drop=FALSE]
	bnum = Anum[,1]^2
	Aden = matrix(bestden,dim(Anum)[1],dim(Anum)[2])
	bden = rep(bestden^2*2,dim(Anum)[1])
	if(dim(Anum)[2]>1) for(k in 2:dim(Anum)[2]) bnum = bnum + Anum[,k]^2
	num = cbind(bnum,-Anum)
	den = cbind(bden,Aden)
	num = rbind(num, cbind(matrix(0,dim(Anum)[2],1),diag(dim(Anum)[2])) )
	den = rbind(den, matrix(1,dim(Anum)[2],dim(Anum)[2]+1) )
	AAp = LRS(num,den)
	RhoFp = AAp$Radius
	RhoFp
}

RhoFcal <- function(p,Anum,TCR,wusenum) {  # Calculate fill distance with integer inputs
	for(k in 1:p) Anum[,k] = Anum[,k] *wusenum[k]
	for(k in dim(Anum)[2]:1) if(prod(Anum[,k]==0)==1) { Anum = Anum[,-k,drop=FALSE]; TCR = TCR[-k]; p = p-1; }
	for(i in dim(Anum)[1]:1) if(prod(Anum[i,]==0)==1) Anum = Anum[-i,,drop=FALSE]
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
		if(max(c(bnum,bden))>2^15) { tobreak = TRUE; break; }  # make sure later computed numbers are below 2^31-1. 
	}
	if(tobreak==TRUE)  return(RhoFcalappro(Anum%*%diag(1/TCR)));  # have to add this patch to avoid error. Big integer causes problems. 
	num = cbind(bnum,-Anum)
	den = cbind(bden,Aden)
	num = rbind(num, cbind(matrix(0,p,1),diag(p)) )
	den = rbind(den, matrix(1,p,p+1) )
	AA = LRS(num,den)
	RhoF = AA$Radius
	RhoF
}

Gcharateristic <- function(Evor,sthis,wnum) {  # compute a feature statistic of the (L,s,w). Those with the same value are same in r_F and m. 
	p = dim(Evor)[2]
	row1 = rep(1,dim(Evor)[1]) 
	for(k in 1:p) row1 = row1 * (Evor[,k]<=1)
	v = Evor[row1==1,,drop=FALSE]
	ndims = matrix(0,p,p)
	for(i in 1:dim(v)[1]) ndims[sum(v[i,,drop=FALSE]),v[i,,drop=FALSE]==1] = ndims[sum(v[i,,drop=FALSE]),v[i,,drop=FALSE]==1]+1
	dimscore = rbind(wnum,sthis,ndims)

	dimperm = do.call(order, as.data.frame(-t(dimscore)))
	v = v[,dimperm,drop=FALSE]
	s = sthis[dimperm]
	ww = wnum[dimperm]
	dimscore = dimscore[,dimperm,drop=FALSE]
	dimranks = 1:p
	for(k in 2:p) if(identical(dimscore[,k],dimscore[,k-1])) dimranks[k] = dimranks[k-1]
	
	G = matrix(0,0,p)
	L01 = matrix(0,0,p)
	vscore = cbind(matrix(0,dim(v)[1],2+p))
	for(k in 1:p) vscore[,2] = vscore[,2]+v[,k,drop=FALSE]
	for(k in 1:p) vscore[,3+p-dimranks[k]] = vscore[,3+p-dimranks[k]]+v[,k,drop=FALSE]
	
	while(min(vscore[,1])==0) {
		rowopt = do.call(order, as.data.frame(vscore))[1]
		thevor = v[rowopt,]
		G = rbind(G,thevor)
		dimscore = rbind(dimscore,thevor) 
		L01add = L01
		for(k in 1:p) L01add[,k] = ((L01[,k]+thevor[k])==1)
		L01 = rbind(L01,thevor,L01add,drop=FALSE)
		vscore[rowopt,1] = 1
		if(dim(L01add)[1]>0) for(i in 1:dim(L01add)[1]) for(j in 1:dim(v)[1]) if(prod(L01add[i,]==v[j,,drop=FALSE])==1) 
			vscore[j,1] = 1

		dimperm = do.call(order, as.data.frame(-t(dimscore)))
		v = v[,dimperm,drop=FALSE]
		dimscore = dimscore[,dimperm,drop=FALSE]
		dimranks = 1:p
		for(k in 2:p) if(identical(dimscore[,k],dimscore[,k-1])) dimranks[k] = dimranks[k-1]

		G = G[,dimperm,drop=FALSE]
		L01 = L01[,dimperm,drop=FALSE]
		vscore[,3:(p+2)] = 0
		for(k in 1:p) vscore[,3+p-dimranks[k]] = vscore[,3+p-dimranks[k]]+v[,k,drop=FALSE]
	}
	return(c(ww,s,G%*%2^((p-1):0)))
}


FS <- function(p=2,contchoices,n=100,w=rep(1,p),wsupp=NULL,nmin=floor(n*.8),nmax=ceiling(n*1.2),coefF=-4,coefS=1,msC=0,NL=10,NP=100,NJ=10) {
	# the function to generate interleaved lattice-based designs with both fill and separation distance properties. 

	pp = p+length(wsupp)
	coefn = (coefF+coefS)/p
	coefnpp = (coefF+coefS)/pp
	bestse = Inf
	for(wd in 1:100) {
		se = sum(abs(round(w*wd)-w*wd))/wd
		if(se<bestse) {
			wnum = round(w*wd)
			wden = wd
			bestse = se
		}
	}

	GMs <- NULL
	if(p==2) GsName <- data(GMs2, envir=environment())
	if(p==3) GsName <- data(GMs3, envir=environment())
	if(p==4) GsName <- data(GMs4, envir=environment())
	if(p==5) GsName <- data(GMs5, envir=environment())
	if(p==6) GsName <- data(GMs6, envir=environment())
	if(p==7) GsName <- data(GMs7, envir=environment())
	if(p==8) GsName <- data(GMs8, envir=environment())
#	GMs <- get(GsName)

	for(indexG in 1:(dim(GMs)[1]/p)) { 
		GM = matrix(0,p,p)
		for(j in 1:p) GM[,j]=GMs[((indexG-1)*p+1):(indexG*p),j]
		theGmML = newG(GM)
		theGmML = supp(theGmML)
		if(indexG==1)  {  
			GMlist = list(theGmML);  
		}  
		if(indexG>1)  {  
			GMlist[[indexG]] = theGmML;  
		}
	}
	
	sdis = matrix(0,1,p)
	for(k in p:1) {
		sdisadd = sdis
		sdisadd[,k] = 1
		sdis = rbind(sdis,sdisadd)
	}
	sdisbase = sdis
	sdep = matrix(0,dim(sdis)[1],dim(sdis)[1])
	for(i in 1:(dim(sdis)[1]-1)) for(j in (i+1):dim(sdis)[1]) {
		if(prod(sdis[i,]<=sdis[j,])==1) { sdep[i,j] = -1; sdep[j,i] =  1; }
		if(prod(sdis[i,]>=sdis[j,])==1) { sdep[i,j] =  1; sdep[j,i] = -1; }
	}
	sdepbase = sdep
	wori = w; 
	wnumori = wnum; 

	contchoices[,3] = contchoices[,6]^coefF*contchoices[,7]^coefS*100^coefn
	if(length(wsupp)>0)  contchoices[,3] = (n/prod(w)/prod(wsupp))^coefnpp * 
	  sqrt( ( contchoices[,6]*(n/100)^(-1/p)*prod(w)^(1/p) )^2 + sum((wsupp/4)^2) )^coefF *
	  sqrt( ( contchoices[,7]*(n/100)^(-1/p)*prod(w)^(1/p) )^2 + sum((wsupp/4)^2) )^coefS  
	criterionopt = 0
	juniors = matrix(0,0,9+p*3)
	criterionhatbest = 0
	if(length(wsupp)>0) criterionhatbestbound = 0
				
	for(indexc in order(-contchoices[,3])[1:min(c(NL,dim(contchoices)[1]))]) {
	  thechoice = contchoices[ indexc, ]
	  if(length(wsupp)==0) if(thechoice[3]<criterionhatbest) break
	  if(length(wsupp) >0) if(thechoice[3]<criterionhatbestbound) break
	  indexG = thechoice[1]
	  theGmML = GMlist[[indexG]]
	  Plist = matrix(0,0,p*3-theGmML@ps[3])
	  Flist = matrix(0,0,p*3-theGmML@ps[3]+1)
	  
	  for(indexP in 1:NP) {
		if(min(wnum)==max(wnum)) if(indexP>1) break
		set.seed(indexP); 
		oriP = sample(1:p); 
		if(indexP==1) oriP = 1:p
		w = wori[oriP]
		wnum = wnumori[oriP]
		if(min(wnum)!=max(wnum)) { 
			Pch = Gcharateristic(Evor=theGmML@Evor,wnum=wnum,sthis=rep(0,p))
			Pcheck = Plist
			for(kk in 1:(p*3-theGmML@ps[3]))  Pcheck = Pcheck[ Pcheck[,kk]==Pch[kk], ,drop=FALSE] 
			if(dim(Pcheck)[1]>0) next   
			Plist = rbind(Plist, Pch)
		}
		countJ = 0
		scont = thechoice[ (8):(p+7) ] * w
		scont = scont * ( n/(prod(scont)/2^theGmML@ps[3]) )^(1/p)
		seq2 = scont<=2
		while(min(scont)<2) {
			seq2 = scont<=2+1e-15
			scont[seq2] = 2
			if(sum(!seq2)==0) break
			scont[!seq2] = scont[!seq2] * ( n/(prod(scont)/2^theGmML@ps[3]) )^(1/sum(!seq2))
		}
		
		sdis = sdisbase
		for(k in 1:p) sdis[,k] = sdis[,k] + floor(scont[k])
		if(dim(theGmML@expair)[1]>0) for(ii in 1:dim(theGmML@expair)[1])
			if( wnum[theGmML@expair[ii,1]] == wnum[theGmML@expair[ii,2]] )  for(i in 1:dim(sdis)[1]) 
				if( sdis[i,theGmML@expair[ii,1]] > sdis[i,theGmML@expair[ii,2]] ) 
					sdis[i,theGmML@expair[ii,1:2]] = sdis[i,theGmML@expair[ii,2:1]]
		dup = duplicated(sdis)
		sdis = sdis[ dup==0, ,drop=FALSE]
		sdep = sdepbase[ dup==0, dup==0, drop=FALSE]
				
		mhats = sdis[,1]
		for(k in 2:p) mhats = mhats * sdis[,k]
		mhats = mhats / 2^theGmML@ps[3]
		msign = sign(mhats-n)
		sCs = rep(0,dim(sdis)[1])
		for(i in 1:(p-1)) for(j in (i+1):p) if(theGmML@coin[i,j]==1) 
			sCs = sCs + (sdis[,i]==2)*(sdis[,j]==2)*w[i]*w[j]
		
		isOpen = (sCs<=msC)*3  # 0, not sC; 1, exceeded m; 2, tried and proper m; 3, waiting; 4, ready. 
		sdep[,!isOpen] = 0
		for(k in 1:dim(sdep)[1]) if(isOpen[k]==0) sdep[sdep[,k]==-1,k] = 0
		for(k in 1:dim(sdep)[1]) if(msign[k]>=0) sdep[sdep[,k]==-1,k] = 0
		for(k in 1:dim(sdep)[1]) if(msign[k]<=0) sdep[sdep[,k]== 1,k] = 0
		for(k in 1:dim(sdep)[1]) if(isOpen[k]>0) if(prod(sdep[k,]==0)==1) isOpen[k] = 4
		serr = abs(sdis[,1]-scont[1])
		for(k in 2:p) serr = serr + abs(sdis[,k]-scont[k])
	
		while(sum(isOpen==4)>0 & countJ<NJ) {
			thetry = order(serr+(isOpen!=4)*2*p)[1]
			sthis = sdis[thetry,]
			sC = sCs[thetry]/wden^2

			Fcheck = Flist
			Fch = Gcharateristic(Evor=theGmML@Evor,wnum=wnum,sthis=sthis)
			for(kk in 1:(p*3-theGmML@ps[3]))  Fcheck = Fcheck[ Fcheck[,kk]==Fch[kk], ,drop=FALSE] 
			if(dim(Fcheck)[1]>0)  RhoF = Fcheck[1,p*3-theGmML@ps[3]+1]
			if(dim(Fcheck)[1]==0) { 
				Evoruse = theGmML@Evor[,wnum>0,drop=FALSE]
				for(ii in dim(Evoruse)[1]:1) if(sum(Evoruse[ii,])==0) Evoruse = Evoruse[-ii,]
				RhoF = RhoFcal(p=sum(wnum>0), Anum=Evoruse, TCR=sthis[wnum>0]*wden, wnum[wnum>0])
				Flist = rbind(Flist, c(Fch,RhoF) )

				ms = tn(sthis,theGmML@coefC)
				thecoset = order(abs(ms/n-1)*2+(ms>n))[1]
				m = ms[thecoset]
				thetran = theGmML@coefD[ theGmML@coefD[,p+1]+1==thecoset, 1:p, drop=FALSE][1,]

				mhat = prod(sthis)/2^theGmML@ps[3]
				RhoS = 0 
				for(kk in 1:p) RhoS = RhoS + (theGmML@Evor[(p+1):dim(theGmML@Evor)[1],kk]/sthis[kk]*w[kk])^2
				RhoS = sqrt(min(RhoS)) /2
				if(max(sthis)>2) RhoS = min(c(RhoS,w[sthis>2]/sthis[sthis>2]))
				RhoFUB = sqrt( RhoF^2 + sum((wsupp/2)^2) )
				RhoSLB = RhoS
				criterionhat = RhoF^coefF*RhoS^coefS*(mhat/prod(w))^coefn
				criterion = RhoF^coefF*RhoS^coefS*(m/prod(w))^coefn
				if(length(wsupp)>0) criterionhatbound = (mhat/prod(w)/prod(wsupp))^coefnpp * RhoFUB^coefF * RhoSLB^coefS  
				Gv = thetran 
				for(ii in 1:(p-theGmML@ps[3])) Gv = Gv + theGmML@G[ii,] * 2^ii

				juniors = rbind( juniors, c(indexG, w,sthis,criterion,criterionhat,RhoF,RhoS,mhat,m,thecoset,sC,Gv) )
				if(m<=nmax) if(m>=nmin) if(criterionhat>criterionhatbest) criterionhatbest = criterionhat
				if(m<=nmax) if(m>=nmin) if(length(wsupp)>0) if(criterionhatbound>criterionhatbestbound) criterionhatbestbound = criterionhatbound
			}
			
			if(mhat<n) if(max(ms)>=nmin) { 
				isOpen[thetry] = 2; 
				sdep[,thetry] =  0; 
				for(ii in 1:dim(sdis)[1]) if(isOpen[ii]==3) if(prod(sdep[ii,]==0)==1) isOpen[ii]=4; 
			}
			if(mhat<n) if(max(ms)<nmin) { 
				isOpen[thetry] = 1; 
				isOpen[ (sdep[,thetry]==-1) * (isOpen==3) == 1 ] = 1; 
			}
			if(mhat>n) if(min(ms)<=nmax) { 
				isOpen[thetry] = 2; 
				sdep[,thetry] =  0; 
				for(ii in 1:dim(sdis)[1]) if(isOpen[ii]==3) if(prod(sdep[ii,]==0)==1) isOpen[ii]=4; 
			}
			if(mhat>n) if(min(ms)>nmax) { 
				isOpen[thetry] = 1; 
				isOpen[ (sdep[,thetry]== 1) * (isOpen==3) == 1 ] = 1; 
			}
			if(mhat==n) isOpen[thetry] = 2 
			countJ = countJ+1 
		}
	  }
	}

	return(juniors)
}


FSsupp <- function(p=2,n=100,wsupp,juniors,nmin=floor(n*.8),nmax=ceiling(n*1.2),coefF=-4,coefS=1,msC=0,NS=100) {
	# the function to supplement columns to interleaved lattice-based designs with both fill and separation distance properties. 

	pp = p+length(wsupp)
	coefn = (coefF+coefS)/pp
	seniors = juniors
	seniors = seniors[ seniors[,2*p+9]<=max(c(msC, min(seniors[,2*p+9]))), ,drop=FALSE]
	seniors = seniors[ seniors[,2*p+7]<=max(c(nmax,min(seniors[,2*p+7]))), ,drop=FALSE]
	seniors = seniors[ seniors[,2*p+7]>=min(c(nmin,max(seniors[,2*p+7]))), ,drop=FALSE]
	seniors = unique(seniors)
	seniors = seniors[ order(abs(seniors[,2*p+7]-n)), ,drop=FALSE]
	seniors = seniors[ order(-seniors[,2*p+3]), ,drop=FALSE]
	prodw = prod(juniors[2:(p+1)])
	crbound = (seniors[,2*p+6]/prodw/prod(wsupp))^coefn * 
	  sqrt( seniors[,2*p+4]^2 + sum((wsupp/4)^2) )^coefF *
	  sqrt( seniors[,2*p+5]^2 + sum((wsupp/4)^2) )^coefS  
	finals = matrix(0,0,9+pp*3)
	
	GMs <- NULL
	if(p==2) GsName <- data(GMs2, envir=environment())
	if(p==3) GsName <- data(GMs3, envir=environment())
	if(p==4) GsName <- data(GMs4, envir=environment())
	if(p==5) GsName <- data(GMs5, envir=environment())
	if(p==6) GsName <- data(GMs6, envir=environment())
	if(p==7) GsName <- data(GMs7, envir=environment())
	if(p==8) GsName <- data(GMs8, envir=environment())
#	GMs <- get(GsName)
	
	criterionhatbest = 0
	for(indexs in 1:min(c(NS,dim(seniors)[1]))) {
		if(crbound[indexs] < criterionhatbest) next;
		thechoice = seniors[ indexs, ]
		indexG = thechoice[1]
		s = thechoice[(p+2):(2*p+1)]
		w = thechoice[2:(p+1)]
		mhat = thechoice[p*2+6]
		m = thechoice[p*2+7]
		thecoset = thechoice[p*2+8]
		sC = thechoice[p*2+9]
		
		Gv = thechoice[(2*p+10):(3*p+9)]
		thetran = Gv - floor(Gv/2)*2
		Gv = ( Gv - thetran )/2
		G01 = matrix(0,0,p) 
		while(max(Gv)>0) {
			G01 = rbind(G01, Gv - floor(Gv/2)*2) 
			Gv = ( Gv - G01[dim(G01)[1],] )/2
		}
		Ebase = matrix(0,1,dim(G01)[2])
		for(i in 1:dim(G01)[1]) {
			Ebasenew = Ebase 
			for(j in 1:dim(G01)[2]) if(G01[i,j]==1) Ebasenew[,j] = 1-Ebasenew[,j]
			Ebase = rbind(Ebase,Ebasenew)
		}
		L01 = Ebase
		
		while(length(w)<pp) {
			isOpen = rep(1,dim(L01)[1])
			isOpen[1] = 0
			Ebase = matrix(0,1,dim(L01)[2])
			G01 = matrix(0,0,dim(L01)[2])
			dist2 = rep(0,dim(L01)[1])
			for(j in 1:dim(G01)[2]) dist2 = dist2 + (L01[,j]*w[j]/s[j])^2
			
			while(sum(isOpen)>0) {
				therow = order(dist2 - max(dist2)*isOpen*2)[1]
				Ebasenew = Ebase 
				for(j in 1:dim(L01)[2]) if(L01[therow,j]==1) Ebasenew[,j] = 1-Ebasenew[,j]
				for(j in 1:dim(L01)[1]) if(isOpen[j]) for(i in 1:dim(Ebasenew)[1]) if(prod(Ebasenew[i,]==L01[j,])==1)
					isOpen[j] = 0
				Ebase = rbind(Ebase,Ebasenew)
				G01 = rbind(G01,L01[therow,])
			}
			
			sCadd = rep(0,2^dim(G01)[1]) 
			sCadd[1] = Inf  # only one level for the new column, the worst scenario. 
			for(j in 1:dim(G01)[2]) if(s[j]==2) sCadd[sum(G01[,j]*2^((dim(G01)[1]-1):0))+1] = sCadd[sum(G01[,j]*2^((dim(G01)[1]-1):0))+1] + w[j]*wsupp[length(w)-p+1]
			indexadd = order(-sCadd)[2^dim(G01)[1]]
			sC = sC + sCadd[indexadd]
			theadd = rep(0,dim(G01)[1]) 
			for(j in 1:dim(G01)[1]) theadd[j] = floor(((indexadd-1)%%2^(dim(G01)[1]-j+1))/2^(dim(G01)[1]-j))
			
			G01 = cbind(G01,theadd)
			s = c(s,2)
			w = c(w,wsupp[length(w)-p+1])
			thetran = c(thetran,0)
			Ebase = matrix(0,1,dim(G01)[2])
			for(i in 1:dim(G01)[1]) {
				Ebasenew = Ebase 
				for(j in 1:dim(G01)[2]) if(G01[i,j]==1) Ebasenew[,j] = 1-Ebasenew[,j]
				Ebase = rbind(Ebase,Ebasenew)
			}
			L01 = Ebase
		}
		
		Evor = diag(pp)*2
		D0 = L01[-1,,drop=FALSE]
		j1=1 
		while(1) {
			if(j1>dim(D0)[1]) break
			todelete=FALSE
			for(j2 in 1:dim(D0)[1]) if(j2!=j1) if(prod((D0[j1,]-D0[j2,])>=0)) todelete=TRUE
			if(todelete==TRUE) D0=D0[-j1,] else j1=j1+1
		}
		Evor = rbind(Evor,D0)

		bestse = Inf
		for(wd in 1:100) {
			se = sum(abs(round(w*wd)-w*wd))/wd
			if(se<bestse) {
				wnum = round(w*wd)
				wden = wd
				bestse = se
			}
		}
		
		Evoruse = Evor
		for(ii in pp:1) if(wnum[ii]==0) Evoruse = Evoruse[-ii,]
		RhoF = RhoFcal(p=sum(wnum>0), Anum=Evoruse[,wnum>0,drop=FALSE], TCR=s[wnum>0]*wden, wnum[wnum>0])
		RhoS = 0 
		for(kk in 1:pp) RhoS = RhoS + (Evor[(pp+1):dim(Evor)[1],kk]/s[kk]*w[kk])^2
		RhoS = sqrt(min(RhoS)) /2
		if(max(s)>2) RhoS = min(c(RhoS,w[s>2]/s[s>2]))

		criterionhat = RhoF^coefF*RhoS^coefS*(mhat/prod(w))^coefn
		criterion = RhoF^coefF*RhoS^coefS*(m/prod(w))^coefn
		
		Gv = thetran 
		for(ii in 1:dim(G01)[1]) Gv = Gv + G01[ii,] * 2^ii
		
		finals = rbind( finals, c(indexG, w,s,criterion,criterionhat,RhoF,RhoS,mhat,m,thecoset,sC,Gv) )
		if(criterionhat>criterionhatbest) criterionhatbest = criterionhat
	}
	
	return(finals)
}


FSE <- function(p,thechoice) {  # Generate integer-valued design (levels 0 to s-1) based on information
	dimord = order(-thechoice[2:(p+1)])
	sbest = thechoice[(p+2):(p*2+1)]

	Gv = thechoice[(2*p+10):(3*p+9)]
	thetran = Gv - floor(Gv/2)*2
	Gv = ( Gv - thetran )/2
	G01 = matrix(0,0,p) 
	while(max(Gv)>0) {
		G01 = rbind(G01, Gv - floor(Gv/2)*2) 
		Gv = ( Gv - G01[dim(G01)[1],] )/2
	}
	Ebase = matrix(0,1,dim(G01)[2])
	for(i in 1:dim(G01)[1]) {
		Ebasenew = Ebase 
		for(j in 1:dim(G01)[2]) if(G01[i,j]==1) Ebasenew[,j] = 1-Ebasenew[,j]
		Ebase = rbind(Ebase,Ebasenew)
	}
	for(j in 1:p) if(thetran[j]==1) Ebase[,j] = 1-Ebase[,j]

	Ethis = Ebase
	for(k in p:1) {
		Eadded = Ethis
		if(sbest[k]>=4) for(ii in 1:(floor(sbest[k]/2)-1)) {
			Eadd = Ethis
			Eadd[,k] = Ethis[,k]+2*ii
			Eadded = rbind(Eadded,Eadd)
		}
		if(sbest[k]>=3) if(floor(sbest[k]/2)*2<sbest[k]) {
			Eadd = Ethis[Ethis[,k]==0,,drop=FALSE]
			Eadd[,k] = Eadd[,k]+2*floor(sbest[k]/2)
			Eadded = rbind(Eadded,Eadd)
		}
		Ethis = Eadded
	}	

	Ethis = Ethis[,dimord,drop=FALSE]
	return(Ethis)
}


ILdesign <- function(Ethis,a=1/2) {  # Generate scaled design based on integer-valued design (levels 0 to s-1)
	Dbest = Ethis
	for(k in 1:dim(Ethis)[2]) {
		ss = max(Ethis[,k])+1
		Dbest[,k] = Dbest[,k]/(ss-a)+(1-a)/(ss-a)/2
	}
	return(Dbest)
}


InterleavedFillSepD <- function(p,n,w=rep(1,p),pfrom=p,a=1/2,nmin=floor(n*.8),nmax=ceiling(n*1.2),coefF=-4,coefS=1,msC=0,NL=10,NP=100,NJ=10,NS=100) {
	# the main function to generate interleaved lattice-based designs, the discrete optimization part. 

	if(!n>=2|!n==round(n)) stop("n must be an integer greater than one.") 
	if(!p>=2|!p==round(p)) stop("p must be an integer greater than one.") 
	if(!pfrom>=2|!pfrom==round(pfrom)|!pfrom<=p) stop("pfrom must be an integer greater than one and no higher than p.") 
	if(pfrom>8) { pfrom=8; warning("pfrom adjusted to 8."); }
	if(!length(w)==p|!prod(w>0)) stop("w must be a positive p-vector.") 
	if(!a>=0|!a<=1) stop("a must be no less than zero and no greater than one.") 
	if(!nmin<=n|!nmax>=n) stop("nmin must not be greater than n and nmax must not be less than n.") 
	if(!coefF<=0) stop("coefF must be non-positive.") 
	if(!coefS>=0) stop("coefS must be non-negative.") 
	if(!msC>=0) stop("msC must be non-negative.") 
	if(!NL>=1|!NL==round(NL)) stop("NL must be an positive integer.") 
	if(!NP>=1|!NP==round(NP)) stop("NP must be an positive integer.") 
	if(!NJ>=1|!NJ==round(NJ)) stop("NJ must be an positive integer.") 
	if(!NS>=1|!NS==round(NS)) stop("NS must be an positive integer.") 
	
	wori = w
	w = sort(w,decreasing=TRUE)

	contchoices <- NULL
	if(pfrom==2) GsName <- data(CCs2, envir=environment())
	if(pfrom==3) GsName <- data(CCs3, envir=environment())
	if(pfrom==4) GsName <- data(CCs4, envir=environment())
	if(pfrom==5) GsName <- data(CCs5, envir=environment())
	if(pfrom==6) GsName <- data(CCs6, envir=environment())
	if(pfrom==7) GsName <- data(CCs7, envir=environment())
	if(pfrom==8) GsName <- data(CCs8, envir=environment())
#	contchoices <- get(GsName)	# indexG, trial, rhoB, cthick, cdensity, fill, sep, s
	
	if(pfrom < p) {
		juniors = FS(p=pfrom,contchoices,n=n,w=w[1:pfrom],wsupp=w[(pfrom+1):p],nmin=nmin,nmax=nmax,coefF=coefF,coefS=coefS,msC=msC,NL=NL,NP=NP,NJ=NJ)
		thefinals = FSsupp(p=pfrom,n=n,wsupp=w[(pfrom+1):p],juniors=juniors,nmin=nmin,nmax=nmax,coefF=coefF,coefS=coefS,msC=msC,NS=NS)
		thefinals = cbind(thefinals,1:dim(thefinals)[1])
	}
	if(pfrom == p) {
		thefinals = FS(p=p,contchoices,n=n,w=w,wsupp=NULL,nmin=nmin,nmax=nmax,coefF=coefF,coefS=coefS,msC=msC,NL=NL,NP=NP,NJ=NJ)
		thefinals = cbind(thefinals,1:dim(thefinals)[1])
	}
	thefinals = thefinals[ thefinals[,2*p+9]<=max(c(msC, min(thefinals[,2*p+9]))), ,drop=FALSE ]
	thefinals = thefinals[ thefinals[,2*p+7]<=max(c(nmax,min(thefinals[,2*p+7]))), ,drop=FALSE ]
	thefinals = thefinals[ thefinals[,2*p+7]>=min(c(nmin,max(thefinals[,2*p+7]))), ,drop=FALSE ]
	thefinals = thefinals[ order(abs(thefinals[,2*p+7]-n)), ,drop=FALSE ]
	thefinals = thefinals[ order(-thefinals[,2*p+3]),    ,drop=FALSE ]
	
	E = FSE(p=p,thechoice=thefinals[1,])
	D = ILdesign(Ethis=E,a=a)
	return(D[,rank(-wori,ties.method="first"),drop=FALSE])
}


