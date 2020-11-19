#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int tnp2(int* counts, int* lsss, int* coefC, int nnrow, int nncol, int ll) {
	int i,j; 
	int sumoneortwo=0; for(j=0;j<ll;j++) if(lsss[j]>2) sumoneortwo = 1; 
	if(sumoneortwo==0) { 
		int therow=0; 
		int tomul; for(j=0;j<ll;j++) { tomul=1; for(i=0;i<ll-j;i++) tomul*=2; therow += (2-lsss[j])*tomul; }
		for(i=0;i<nncol;i++) counts[i] = coefC[therow+i*nnrow]; 
	}
	if(sumoneortwo) { 
		for(j=0;j<ll;j++) if(lsss[j]>2) break; 
		int thedim = j;
		int thetwos = (lsss[j]-1)/2; 
		int lsss1[ll]; for(j=0;j<ll;j++) lsss1[j]=lsss[j]; lsss1[thedim]=2; 
		int lsss2[ll]; for(j=0;j<ll;j++) lsss2[j]=lsss[j]; lsss2[thedim]=lsss[thedim]-thetwos*2; 
		int counts1[nncol]; tnp2(counts1, lsss1, coefC, nnrow, nncol, ll); 
		int counts2[nncol]; tnp2(counts2, lsss2, coefC, nnrow, nncol, ll); 
		for(i=0;i<nncol;i++) counts[i]=thetwos*counts1[i] + counts2[i]; 
	}
	return 0; 
}

int tnp1(int* counts, int* lsss, int* coefC, int nnrow, int nncol, int ll) {
	int i,j; 
	int sumoneortwo=0; for(j=0;j<ll;j++) if(lsss[j]>2) sumoneortwo = 1; 
	if(sumoneortwo==0) { 
		int therow=1; 
		int tomul; for(j=0;j<ll;j++) { tomul=1; for(i=0;i<ll-j;i++) tomul*=2; therow += (2-lsss[j])*tomul; }
		for(i=0;i<nncol;i++) counts[i] = coefC[therow+i*nnrow]; 
	}
	if(sumoneortwo) { 
		for(j=0;j<ll;j++) if(lsss[j]>2) break; 
		int thedim = j;
		int thetwos = (lsss[j]-1)/2; 
		int lsss1[ll]; for(j=0;j<ll;j++) lsss1[j]=lsss[j]; lsss1[thedim]=2; 
		int lsss2[ll]; for(j=0;j<ll;j++) lsss2[j]=lsss[j]; lsss2[thedim]=lsss[thedim]-thetwos*2; 
		int counts1[nncol]; tnp1(counts1, lsss1, coefC, nnrow, nncol, ll); 
		int counts2[nncol]; tnp1(counts2, lsss2, coefC, nnrow, nncol, ll); 
		for(i=0;i<nncol;i++) counts[i]=thetwos*counts1[i] + counts2[i]; 
	}
	return 0; 
}


void selectlss(int* p_, int* n_, int* ps_, int* expair_nrow_, int* expair, int* coefC, double* maxdissimilarity_, double* minndp_, int* lssm, int* maxnrow_) {
	int p = *p_;
	int n = *n_;
	int ps1 = ps_[0];
	int ps2 = ps_[1];
	int ps3 = ps_[2];
	int expair_nrow = *expair_nrow_; 
	double maxdissimilarity = *maxdissimilarity_;
	double minndp = *minndp_; 
	int minnd = n*(minndp+.0001/n); 
	int maxnrow = *maxnrow_; // the maximum number of rows for lssm. 
	int therow=0;  // the index for the row of lssm to write. 

	int i,j; 
	double lmean = pow((double)n,1/(double)p);
	double lgoal[p]; for(j=0;j<p;j++) lgoal[j]=lmean;
	if(ps2>0) for(j=ps1;j<p;j++) lgoal[j]=lmean*pow(2.0,(double)ps3/(ps2+ps3)); 
	int llow[p]; for(j=0;j<p;j++) llow[j]=lgoal[j]/2;
	for(j=0;j<p;j++) if(llow[j]<1) llow[j]=1; 
	if(ps3>0) for(j=ps1;j<p;j++) if(llow[j]<2) llow[j]=2;  
	int lhigh[p]; for(j=0;j<p;j++) lhigh[j]=lgoal[j]*2+1;
	int lss[p+1]; for(j=0;j<p;j++) lss[j]=llow[j];  // s, followed by 

	double lreal[p]; 
	double td; 
	int tobreak; 
	int highestlastl, lowestn, thecoef, nd;
	int nncol=1;  if(ps3>0) for(j=0;j<ps3;j++) nncol *= 2; 
	int nnrow=1;  if(ps2+ps3>0) for(j=0;j<ps2+ps3;j++) nnrow *= 2; 
	int counts2[nncol], counts1[nncol], twos[nncol], ones[nncol], lastl[nncol], ns[nncol]; 
			
	lss[p-2]=lss[p-2]-1; 
	while(1) {
		lss[p-2] = lss[p-2] +1; 
		if(p>2) for(j=p-2;j>=1;j--) if(lss[j]>lhigh[j]) { lss[j]=llow[j]; lss[j-1]=lss[j-1]+1; } 
		if(lss[0]>lhigh[0]) break; 
		lss[p-1]=lhigh[p-1]; 
	// check exchange without the last dimension
		tobreak=0;  
		if(ps1>1) for(j=0;j<ps1-1;j++) if(lss[j+1]<lss[j]) tobreak=1; 
		if(minnd==n) if(ps1>1) for(j=0;j<ps1-1;j++) if((n/lss[j])*lss[j]!=n) tobreak=1; 
		if(tobreak==1) continue; 
		if(ps2+ps3>1) if(expair_nrow>0) for(i=0;i<expair_nrow;i++) if(expair[i+expair_nrow]+ps1<p) 
			if(lss[expair[i+expair_nrow]+ps1-1]<lss[expair[i]+ps1-1]) { tobreak=1; break; }
		if(tobreak==1) continue; 
	// check dissimilarity without the last dimension
		for(j=0;j<p-1;j++) lreal[j] = lss[j]; 
		if(ps3>0) for(j=ps1;j<p-1;j++) lreal[j] = lreal[j] / pow(2.0,(double)ps3/(ps2+ps3));   
		lmean=1; for(j=0;j<p-1;j++) lmean=lmean*lreal[j];
		lmean=pow(lmean,1/((double)p-1));
		td=0; for(j=0;j<p-1;j++) td=td+fabs( log(lreal[j]/lmean) );  
		if(td > log(maxdissimilarity)) continue;
	// select the coefficient and compute the lss[p-1] 
		if(ps3>0) {  
			tnp2(counts2, lss+ps1, coefC, nnrow, nncol, ps2+ps3-1);
			if(ps1>0) for(i=0;i<nncol;i++) for(j=0;j<ps1;j++) counts2[i] *= lss[j]; 
			tobreak=1; for(i=0;i<nncol;i++) if(counts2[i]<=n) tobreak=0; 
			if(tobreak==1) continue; 
			tnp1(counts1, lss+ps1, coefC, nnrow, nncol, ps2+ps3-1);
			if(ps1>0) for(i=0;i<nncol;i++) for(j=0;j<ps1;j++) counts1[i] *= lss[j]; 
			for(i=0;i<nncol;i++) twos[i] = n/counts2[i];
			for(i=0;i<nncol;i++) ones[i] = (n-twos[i]*counts2[i] >= counts1[i]);
			for(i=0;i<nncol;i++) lastl[i] = twos[i]*2+ones[i]; 
			for(i=0;i<nncol;i++) ns[i] = twos[i]*counts2[i]+ones[i]*counts1[i]; 
			highestlastl = 1; 
			lowestn = n+1; 
			thecoef = -1; 
			for(i=0;i<nncol;i++) if(lastl[i]>=2) {
				if(lastl[i]>highestlastl) { highestlastl=lastl[i]; lowestn=ns[i]; thecoef=i; }
				if(lastl[i]==highestlastl&&ns[i]<lowestn) { highestlastl=lastl[i]; lowestn=ns[i]; thecoef=i; }
			}
			lss[p-1] = lastl[thecoef]; 
			nd = twos[thecoef]*counts2[thecoef] + ones[thecoef]*counts1[thecoef];  
			if(nd<minnd) continue; 
			lss[p]=nd; 
			lss[p+1]=thecoef; 
		}
		if(ps3==0) {
			thecoef = 1;
			lss[p-1] = n; 
			for(j=0;j<p-1;j++) lss[p-1] /= lss[j]; 
			nd = 1; 
			for(j=0;j<p;j++) nd *= lss[j]; 
			if(lss[p-1]<lss[p-2]) continue; 
			if(nd<minnd) continue; 
			lss[p]=nd; 
			lss[p+1]=thecoef; 
		}
	// check exchange with the last dimension
		tobreak=0; 
		if(ps2+ps3>1) if(expair_nrow>0) for(i=0;i<expair_nrow;i++) 
			if(lss[expair[i+expair_nrow]+ps1-1]<lss[expair[i]+ps1-1]) { tobreak=1; break; }
		if(tobreak==1) continue; 
	// check dissimilarity with the last dimension
		for(j=0;j<p;j++) lreal[j] = lss[j]; 
		if(ps3>0) for(j=ps1;j<p;j++) lreal[j] = lreal[j] / pow(2.0,(double)ps3/(ps2+ps3));   
		lmean=1; for(j=0;j<p;j++) lmean=lmean*lreal[j];
		lmean=pow(lmean,1/((double)p)); 
		td=0; for(j=0;j<p;j++) td=td+fabs( log(lreal[j]/lmean) );  
		if(td > log(maxdissimilarity)) continue;

	// save results 
		for(j=0;j<p+2;j++) lssm[therow*(p+2)+j]=lss[j];
		therow += 1; 
		if(therow==maxnrow) therow -= 1;
	}

	*maxnrow_ = therow; 
	return; 
}
