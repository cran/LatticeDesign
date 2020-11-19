#include <stdio.h>
#include <stdlib.h>

int compare_double(const void * a, const void * b)
{
    return ( *(double*)a - *(double*)b >0 );
}

void ProjSep(double* design, int* p_, int* n_, double* minds) {
	int i,j,k;
	int p = *p_;
	int n = *n_;

	double* dists = (double*)malloc(sizeof(double) *p*n*(n-1)/2);

	for(k=0;k<p;k++) {
		for(i=0;i<n-1;i++) for(j=i+1;j<n;j++) {
			dists[ (i*n-i*(i+1)/2 +j-i-1)*p +k ] = ( design[k*n+i] - design[k*n+j] )*( design[k*n+i] - design[k*n+j] );
		}
	}
	for(i=0;i<n*(n-1)/2;i++) {
		qsort(dists+i*p,p,sizeof(double),compare_double);
		for(k=1;k<p;k++) dists[i*p+k] = dists[i*p+k-1] + dists[i*p+k]; 
	}
	
	for(j=0;j<p;j++) minds[j]=j*100+100; 
	for(k=0;k<p;k++) for(i=0;i<n*(n-1)/2;i++) if(dists[i*p+k]<minds[k]) minds[k]=dists[i*p+k];
		
	free(dists);
	
	return; 
}

