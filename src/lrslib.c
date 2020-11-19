/* lrslib.c     library code for lrs                     */

/* modified by Gary Roumanis for multithread plrs compatability */
/* truncate needs mod to supress last pivot */
/* need to add a test for non-degenerate pivot step in reverse I guess */
/* Copyright: David Avis 2005,2011 avis@cs.mcgill.ca         */

/* This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA  02110-1335, USA.
 */



#include <stdio.h>
#include <string.h>
#include "lrslib.h"
#include <math.h>

/* Globals; these need to be here, rather than lrslib.h, so they are
   not multiply defined. */

int* lrs_inv; int* lrs_idv; 
double* lrs_o; 
double* mv; 
double sumradius2;
int lrs_therow, lrs_thecol, lrs_theout, lrs_n, lrs_m, lrs_MO; 
int lrs_indexcol, isOut, MaxOutputLength; 

//FILE *lrs_ifp;			/* input file pointer       */
//FILE *lrs_ofp;			/* output file pointer      */

static unsigned long dict_count, dict_limit, cache_tries, cache_misses;

/* Variables and functions global to this file only */
static long lrs_checkpoint_seconds = 0;

static long lrs_global_count = 0;	/* Track how many lrs_dat records are 
					   allocated */

static lrs_dat_p *lrs_global_list[MAX_LRS_GLOBALS + 1];

static lrs_dic *new_lrs_dic (long m, long d, long m_A);


static void cache_dict (lrs_dic ** D_p, lrs_dat * global, long i, long j);
static long check_cache (lrs_dic ** D_p, lrs_dat * global, long *i_p, long *j_p);
static void save_basis (lrs_dic * D, lrs_dat * Q);

static void pushQ (lrs_dat * global, long m, long d, long m_A);

#ifdef TIMES
static void ptimes ();
static double get_time();
#endif


/*******************************/
/* signals handling            */
/*******************************/
#ifdef SIGNALS
static void checkpoint ();
static void die_gracefully ();
static void setup_signals ();
static void timecheck ();
#endif

/*******************************/
/* functions  for external use */
/*******************************/

void prat_XuHe (lrs_mp Nin, lrs_mp Din)	/*reduce and print Nin/Din  */
{
	rattodouble(Nin,Din,lrs_o+lrs_MO);
	if(mv[lrs_indexcol+1]<lrs_o[lrs_MO]) mv[lrs_indexcol+1] = lrs_o[lrs_MO];
	sumradius2 += pow(lrs_o[lrs_MO],2.0);
	lrs_indexcol += 1;
	lrs_MO += 1;
	if(lrs_MO==MaxOutputLength) { lrs_MO=0; isOut=1; } 
	if(lrs_indexcol==lrs_n-1) { 
		lrs_indexcol=0;
		if(sumradius2>mv[0]) mv[0] = sumradius2;
		sumradius2 = 0; 
	}
}

/*
void pmp_XuHe (lrs_mp a)	// print the long precision integer a 
{
	lrs_onv[lrs_MO] = (int) a[1];
	lrs_odv[lrs_MO] = 1;
	if((int)a[0]<0) lrs_onv[lrs_MO] = -lrs_onv[lrs_MO];
	lrs_MO += 1;
}
*/

void lrs_XuHe(int* nhrep_, int* ndim_, int* input_numerator, int* input_denominator, 
  double* output_vertex, double* mv_, int* MaxOutput_) {
	// results: the maximum value of each dimension, followed by the covering radius. 
	
// Step 0: Initialize 

	lrs_inv = input_numerator;  // by columns
	lrs_idv = input_denominator;  // by columns
	lrs_o   = output_vertex;  // by rows
	lrs_therow = 0; 
	lrs_thecol = 0;
	lrs_m = *nhrep_;  // number of input rows
	lrs_n = *ndim_;  // number of input columns
	lrs_indexcol = 0;  // index for the site for writing outputs
	lrs_MO = 0;  // index for the site for writing outputs
	MaxOutputLength = *MaxOutput_;  // the number of rows of outputs, maximum allowed as inputs and actual as outputs
	isOut = 0;  // whether has go beyond the maximum output. 
	sumradius2 = 0; // sum of squared radius for the current vertex. 
	mv = mv_;  // the maximum radius followed by the maximum value of the dimension 
	
	lrs_dic *P;			/* structure for holding current dictionary and indices */
	lrs_dat *Q;			/* structure for holding static problem data            */
	lrs_mp_vector output;		/* holds one line of output; ray,vertex,facet,linearity */
	lrs_mp_matrix Lin;		/* holds input linearities if any are found             */
	long col;			/* output column index for dictionary                   */
	long startcol = 0;
	long prune = FALSE;		/* if TRUE, getnextbasis will prune tree and backtrack  */

//	lrs_mp_init (ZERO, stdin, stdout)) 
//	lrs_ifp = stdin;
//	lrs_ofp = stdout;
	lrs_record_digits = 0;
	long dec_digits = DEFAULT_DIGITS;
	lrs_digits = DEC2DIG (dec_digits);	/* max permitted no. of digits   */
	if (lrs_digits > MAX_DIGITS)  lrs_digits = MAX_DIGITS;

	lrs_global_count = 0;
	lrs_checkpoint_seconds = 0;

	
/* Step 1: Allocate lrs_dat, lrs_dic and set up the problem                      */

	Q = lrs_alloc_dat ("LRS globals");	/* allocate and init structure for static problem data */

	Q->hull = FALSE;
	stringcpy (Q->fname, "VorVerCal");
	Q->m = (long) *nhrep_;  // number of input rows
	Q->n = (long) *ndim_;  // number of input columns (p+1)

	P = lrs_alloc_dic (Q);	/* allocate and initialize lrs_dic                     */

	lrs_read_dic (P, Q);  /* read remainder of input to setup P and Q            */

//	lrs_theout = 0; 
	output = lrs_alloc_mp_vector (Q->n);	/* output holds one line of output from dictionary     */


/* Step 2: Find a starting cobasis from default of specified order               */
/*         P is created to hold  active dictionary data and may be cached        */
/*         Lin is created if necessary to hold linearity space                   */
/*         Print linearity space if any, and retrieve output from first dict.    */

	lrs_getfirstbasis (&P, Q, &Lin, FALSE);


	if(!Q->restart)
		for (col = startcol; col < Q->nredundcol; col++)	/* print linearity space               */
			lrs_printoutput (Q, Lin[col]);	/* Array Lin[][] holds the coeffs.     */


  prune=lrs_checkbound(P,Q);
  do
  {
//2015.6.5   after maxcobases reached, generate subtrees that have not been enumerated
    if ((Q->maxcobases > 0) &&  (Q->count[2] >=Q->maxcobases))
    {
       prune=TRUE;
    }     // if Q-> maxcobases...

    for (col = 0; col <= P->d; col++)          /* print output vertex/ray if any */
	    if (lrs_getsolution (P, Q, output, col))
	        lrs_printoutput (Q, output);
  }   while (!Q->lponly && lrs_getnextbasis (&P, Q, prune));  // do ...
 
   lrs_clear_mp_vector (output, Q->n);
   lrs_free_dic (P,Q);           /* deallocate lrs_dic */
   lrs_free_dat (Q);             /* deallocate lrs_dat */
 
	mv[0] = sqrt(mv[0]); 
	*MaxOutput_ = lrs_MO; 
	if(isOut==1) *MaxOutput_ = -1;  // Output length is not enough.
	return;
}




/*******************/
/* lrs_printoutput */
/*******************/
void 
lrs_printoutput (lrs_dat * Q, lrs_mp_vector output)
{

if (Q->countonly)
   return;

  long i;

  if (Q->hull || zero (output[0]))	/*non vertex */
    {
      for (i = 0; i < Q->n; i++)
	prat_XuHe (output[i], output[0]);

    }
  else
    {				/* vertex   */
      for (i = 1; i < Q->n; i++)
	prat_XuHe (output[i], output[0]);  // /*reduce and print output[i]/output[0]  */
    }
}
/**************************/
/* end of lrs_printoutput */
/**************************/

long 
lrs_getsolution (lrs_dic * P, lrs_dat * Q, lrs_mp_vector output, long col)
   /* check if column indexed by col in this dictionary */
   /* contains output                                   */
   /* col=0 for vertex 1....d for ray/facet             */
{

	
  long j;			/* cobasic index     */

  lrs_mp_matrix A = P->A;
  long *Row = P->Row;

  if (col == ZERO)		/* check for lexmin vertex */
    return lrs_getvertex (P, Q, output);

/*  check for rays: negative in row 0 , positive if lponly */

  if (Q->lponly)
    {
      if (!positive (A[0][col]))
	return FALSE;
    }

  else if (!negative (A[0][col]))
    return FALSE;

/*  and non-negative for all basic non decision variables */

  j = Q->lastdv + 1;
  while (j <= P->m && !negative (A[Row[j]][col]))
    j++;

  if (j <= P->m)
    return FALSE;

  if (Q->geometric || Q->allbases || lexmin (P, Q, col) || Q->lponly)

    return lrs_getray (P, Q, col, Q->n, output);

  return FALSE;			/* no more output in this dictionary */

}				/* end of lrs_getsolution */


/***********************************/
/* allocate and initialize lrs_dat */
/***********************************/
lrs_dat *
lrs_alloc_dat (const char *name)
{
  lrs_dat *Q;
  long i;

//  if (lrs_global_count >= MAX_LRS_GLOBALS)  return Q;
  
  Q = (lrs_dat *) malloc (sizeof (lrs_dat));
  if (Q == NULL)
    return Q;			/* failure to allocate */

  lrs_global_list[lrs_global_count] = Q;
  Q->id = lrs_global_count;
  lrs_global_count++;
  Q->name=(char *) CALLOC ((unsigned) strlen(name)+1, sizeof (char));
  strcpy(Q->name,name); 

/* initialize variables */
  Q->m = 0L;
  Q->n = 0L;
  Q->inputd = 0L;
  Q->deepest = 0L;
  Q->nlinearity = 0L;
  Q->nredundcol = 0L;
  Q->runs = 0L;
  Q->subtreesize=MAXD;
  Q->seed = 1234L;
  Q->totalnodes = 0L;
  for (i = 0; i < 10; i++)
    {
      Q->count[i] = 0L;
      Q->cest[i] = 0.0;
      if(i < 5)
	Q->startcount[i] = 0L;
    }
  Q->count[2] = 1L;		/* basis counter */
  Q->startcount[2] = 0L;	/* starting basis counter */
/* initialize flags */
  Q->allbases = FALSE;
  Q->bound = FALSE;            /* upper/lower bound on objective function given */
  Q->countonly = FALSE;        /* produce the usual output */
  Q->debug = FALSE;
  Q->frequency = 0L;
  Q->dualdeg = FALSE;          /* TRUE if dual degenerate starting dictionary */
  Q->geometric = FALSE;
  Q->getvolume = FALSE;
  Q->homogeneous = TRUE;
  Q->polytope = FALSE;
  Q->hull = FALSE;
  Q->incidence = FALSE;
  Q->lponly = FALSE;
  Q->maxdepth = MAXD;
  Q->mindepth = -MAXD;
  Q->maxoutput = 0L;
  Q->maxcobases = 0L;     /* after maxcobases have been found unexplored subtrees reported */
  Q->nash = FALSE;
  Q->nonnegative = FALSE;
  Q->printcobasis = FALSE;
  Q->printslack = FALSE;
  Q->truncate = FALSE;          /* truncate tree when moving from opt vertex        */
  Q->verbose=FALSE;
  Q->voronoi = FALSE;
  Q->maximize = FALSE;		/*flag for LP maximization                          */
  Q->minimize = FALSE;		/*flag for LP minimization                          */
  Q->restart = FALSE;		/* TRUE if restarting from some cobasis             */
  Q->givenstart = FALSE;	/* TRUE if a starting cobasis is given              */
  Q->strace = -1L;		/* turn on  debug at basis # strace */
  Q->etrace = -1L;		/* turn off debug at basis # etrace */

  Q->saved_flag = 0;		/* no cobasis saved initially, db */
  lrs_alloc_mp (Q->Nvolume);
  lrs_alloc_mp (Q->Dvolume);
  lrs_alloc_mp (Q->sumdet);
  lrs_alloc_mp (Q->saved_det);
  lrs_alloc_mp (Q->boundn);
  lrs_alloc_mp (Q->boundd);
  itomp (ZERO, Q->Nvolume);
  itomp (ONE, Q->Dvolume);
  itomp (ZERO, Q->sumdet);
/* 2012.6.1 */
  Q->unbounded = FALSE;

  return Q;
}				/* end of allocate and initialize lrs_dat */

/****************************/
/* set up lrs_dic structure */
/****************************/
long 
lrs_read_dic (lrs_dic * P, lrs_dat * Q)
/* read constraint matrix and set up problem and dictionary  */

{
  lrs_mp Temp, mpone, mpten;
  lrs_mp_vector oD;		/* Denom for objective function */

  long i, j;

/* assign local variables to structures */

  lrs_mp_matrix A;
  lrs_mp_vector Gcd, Lcm;
  long hull = Q->hull;
  long m, d;

  lrs_alloc_mp(Temp); lrs_alloc_mp(mpone);
  lrs_alloc_mp(Tempn); lrs_alloc_mp(Tempd); lrs_alloc_mp(mpten);
  A = P->A;
  m = Q->m;
  d = Q->inputd;
  Gcd = Q->Gcd;
  Lcm = Q->Lcm;

  oD = lrs_alloc_mp_vector (d);

  itomp (ONE, mpone);
  itomp(10L,mpten);
  itomp (ONE, A[0][0]);
  itomp (ONE, Lcm[0]);
  itomp (ONE, Gcd[0]);

  for (i = 1; i <= m; i++)	/* read in input matrix row by row                 */ 
  {
    itomp (ONE, Lcm[i]);	/* Lcm of denominators */
    itomp (ZERO, Gcd[i]);	/* Gcd of numerators */
    for (j = hull; j <= d; j++)	/* hull data copied to cols 1..d */
	{
	  if (readrat (A[i][j], A[0][j]))
	    lcm (Lcm[i], A[0][j]);	/* update lcm of denominators */
	  copy (Temp, A[i][j]);
	  gcd (Gcd[i], Temp);	/* update gcd of numerators   */
	}
 
    if (hull)
	{
	  itomp (ZERO, A[i][0]);	/*for hull, we have to append an extra column of zeroes */
	  if (!one (A[i][1]) || !one (A[0][1]))		/* all rows must have a one in column one */
	    Q->polytope = FALSE;
	}
    if (!zero (A[i][hull]))	/* for H-rep, are zero in column 0     */
	  Q->homogeneous = FALSE;	/* for V-rep, all zero in column 1     */

    storesign (Gcd[i], POS);
    storesign (Lcm[i], POS);

    if (mp_greater (Gcd[i], mpone) || mp_greater (Lcm[i], mpone))
	  for (j = 0; j <= d; j++)
	  {
	    exactdivint (A[i][j], Gcd[i], Temp);	/*reduce numerators by Gcd  */
	    mulint (Lcm[i], Temp, Temp);	/*remove denominators */
	    exactdivint (Temp, A[0][j], A[i][j]);	/*reduce by former denominator */
	  }
  }				/* end of for i=       */

  lrs_clear_mp(Temp); lrs_clear_mp(mpone);
  lrs_clear_mp(Tempn); lrs_clear_mp(Tempd); lrs_clear_mp(mpten);
  lrs_clear_mp_vector (oD,d);
  return TRUE;

}

 		/* end of if(voronoi)     */

/* In lrs_getfirstbasis and lrs_getnextbasis we use D instead of P */
/* since the dictionary P may change, ie. &P in calling routine    */

#define D (*D_p)

long 
lrs_getfirstbasis (lrs_dic ** D_p, lrs_dat * Q, lrs_mp_matrix * Lin, long no_output)
/* gets first basis, FALSE if none              */
/* P may get changed if lin. space Lin found    */
/* no_output is TRUE supresses output headers   */
{
  long i, j, k;

/* assign local variables to structures */

  lrs_mp_matrix A;
  long *B, *C, *Col;
  long *inequality;
  long *linearity;
  long hull = Q->hull;
  long m, d, lastdv, nlinearity, nredundcol;

  lrs_alloc_mp(Temp); lrs_alloc_mp(scale);

  if (Q->lponly)
    no_output = TRUE;
  m = D->m;
  d = D->d;

  lastdv = Q->lastdv;
  nredundcol = 0L;		/* will be set after getabasis        */
  nlinearity = Q->nlinearity;	/* may be reset if new linearity read or in getabasis*/
  linearity = Q->linearity;

  A = D->A;
  B = D->B;
  C = D->C;
  Col = D->Col;
  inequality = Q->inequality;

  if (Q->runs > 0)		/* arrays for estimator */
    {
      Q->isave = (long *) CALLOC ((unsigned) (m * d), sizeof (long));
      Q->jsave = (long *) CALLOC ((unsigned) (m * d), sizeof (long));
    }
/* default is to look for starting cobasis using linearies first, then     */
/* filling in from last rows of input as necessary                         */
/* linearity array is assumed sorted here                                  */
/* note if restart/given start inequality indices already in place         */
/* from nlinearity..d-1                                                    */
  for (i = 0; i < nlinearity; i++)	/* put linearities first in the order */
    inequality[i] = linearity[i];

  k = 0;			/* index for linearity array   */

  if (Q->givenstart)
    k = d;
  else
    k = nlinearity;
  for (i = m; i >= 1; i--)
    {
      j = 0;
      while (j < k && inequality[j] != i)
	j++;			/* see if i is in inequality  */
      if (j == k)
	inequality[k++] = i;
    }
/* for voronoi convert to h-description using the transform                  */
/* a_0 .. a_d-1 -> (a_0^2 + ... a_d-1 ^2)-2a_0x_0-...-2a_d-1x_d-1 + x_d >= 0 */
/* note constant term is stored in column d, and column d-1 is all ones      */
/* the other coefficients are multiplied by -2 and shifted one to the right  */
  if (!Q->maximize && !Q->minimize)
    for (j = 0; j <= d; j++)
      itomp (ZERO, A[0][j]);

/* Now we pivot to standard form, and then find a primal feasible basis       */
/* Note these steps MUST be done, even if restarting, in order to get         */
/* the same index/inequality correspondance we had for the original prob.     */
/* The inequality array is used to give the insertion order                   */
/* and is defaulted to the last d rows when givenstart=FALSE                  */

  if(Q->nonnegative) 
   {
/* no need for initial pivots here, labelling already done */
     Q->lastdv = d;
     Q->nredundcol = 0;
   }
  else
  {
     if (!getabasis (D, Q, inequality))
          return FALSE;
/* bug fix 2009.12.2 */
     nlinearity=Q->nlinearity;   /*may have been reset if some lins are redundant*/
  }

  nredundcol = Q->nredundcol;
  lastdv = Q->lastdv;
  d = D->d;


/* Reset up the inequality array to remember which index is which input inequality */
/* inequality[B[i]-lastdv] is row number of the inequality with index B[i]              */
/* inequality[C[i]-lastdv] is row number of the inequality with index C[i]              */

  for (i = 1; i <= m; i++)
    inequality[i] = i;
  if (nlinearity > 0)		/* some cobasic indices will be removed */
    {
      for (i = 0; i < nlinearity; i++)	/* remove input linearity indices */
	inequality[linearity[i]] = 0;
      k = 1;			/* counter for linearities         */
      for (i = 1; i <= m - nlinearity; i++)
	{
	  while (k <= m && inequality[k] == 0)
	    k++;		/* skip zeroes in corr. to linearity */
	  inequality[i] = inequality[k++];
	}
    }
				/* end if linearity */

  if (nredundcol > 0)
    {
      const unsigned int Qn = Q->n;
      *Lin = lrs_alloc_mp_matrix (nredundcol, Qn);

      for (i = 0; i < nredundcol; i++)
	{
		

	  if (!(Q->homogeneous && Q->hull && i == 0))	/* skip redund col 1 for homog. hull */
	    {
		
	      lrs_getray (D, Q, Col[0], D->C[0] + i - hull, (*Lin)[i]);		/* adjust index for deletions */
	    }

	  if (!removecobasicindex (D, Q, 0L))
	    {
	      lrs_clear_mp_matrix (*Lin, nredundcol, Qn);
	      return FALSE;
	    }
	}
	

    }				/* end if nredundcol > 0 */

/* Do dual pivots to get primal feasibility */
  if (!primalfeasible (D, Q))
    {
      return FALSE;
    }

/* Now solve LP if objective function was given */
  if (Q->maximize || Q->minimize)
    {
      Q->unbounded = !lrs_solvelp (D, Q, Q->maximize);
      if (Q->lponly)		
        {
          lrs_clear_mp(Temp); lrs_clear_mp(scale);
          return TRUE;
        }
      else                         /* check to see if objective is dual degenerate */
       {
	  j = 1;
	  while (j <= d && !zero (A[0][j]))
	    j++;
	  if (j <= d)
	    Q->dualdeg = TRUE;
	}
    }
  else
/* re-initialize cost row to -det */
    {
      for (j = 1; j <= d; j++)
	{
	  copy (A[0][j], D->det);
	  storesign (A[0][j], NEG);
	}

      itomp (ZERO, A[0][0]);	/* zero optimum objective value */
    }

/* reindex basis to 0..m if necessary */
/* we use the fact that cobases are sorted by index value */

  while (C[0] <= m)
    {
      i = C[0];
      j = inequality[B[i] - lastdv];
      inequality[B[i] - lastdv] = inequality[C[0] - lastdv];
      inequality[C[0] - lastdv] = j;
      C[0] = B[i];
      B[i] = i;
      reorder1 (C, Col, ZERO, d);
    }

  if (Q->restart)
    {
      if (!restartpivots (D, Q))
	return FALSE;
      D->lexflag = lexmin (D, Q, ZERO);		/* see if lexmin basis */
   }

/* Check to see if necessary to resize */
  if (Q->inputd > d)
    *D_p = resize (D, Q);

  lrs_clear_mp(Temp); lrs_clear_mp(scale);
  return TRUE;
}
/********* end of lrs_getfirstbasis  ***************/


/*****************************************/
/* getnextbasis in reverse search order  */
/*****************************************/


long 
lrs_getnextbasis (lrs_dic ** D_p, lrs_dat * Q, long backtrack)
	 /* gets next reverse search tree basis, FALSE if none  */
	 /* switches to estimator if maxdepth set               */
	 /* backtrack TRUE means backtrack from here            */

{
  /* assign local variables to structures */
  long i = 0L, j = 0L;
  long m = D->m;
  long d = D->d;
  long saveflag;
  long cob_est=0;     /* estimated number of cobases in subtree from current node */


  if (backtrack && D->depth == 0)
    return FALSE;                       /* cannot backtrack from root      */

  if (Q->maxoutput > 0 && Q->count[0]+Q->count[1]-Q->hull >= Q->maxoutput)
     return FALSE;                      /* output limit reached            */

  while ((j < d) || (D->B[m] != m))	/*main while loop for getnextbasis */
    {
      if (D->depth >= Q->maxdepth)
        {
          if (Q->runs > 0 && !backtrack )     /*get an estimate of remaining tree */
            {

//2015.2.9 do iterative estimation backtracking when estimate is small
                 saveflag=Q->printcobasis;
                 Q->printcobasis=FALSE;
                 cob_est=lrs_estimate (D, Q);
                 Q->printcobasis=saveflag;
                 if(cob_est <= Q->subtreesize) /* stop iterative estimation */
                  {
                    if(cob_est > 0)   /* when zero we are at a leaf */
                       {  ; //lrs_printcobasis(D,Q,ZERO);
                       }
                    backtrack=TRUE;
                  }

            }
            else    // either not estimating or we are backtracking

              if (!backtrack && !Q->printcobasis) 
                 if(!lrs_leaf(D,Q))    /* 2015.6.5 cobasis returned if not a leaf */
                      ; //lrs_printcobasis(D,Q,ZERO);

            backtrack = TRUE;

 
	     if (Q->maxdepth == 0 && cob_est <= Q->subtreesize)	/* root estimate only */
		 { return FALSE;	} /* no nextbasis  */
       }     // if (D->depth >= Q->maxdepth)


/*      if ( Q->truncate && negative(D->A[0][0]))*/   /* truncate when moving from opt. vertex */
/*          backtrack = TRUE;    2011.7.14 */

      if (backtrack)		/* go back to prev. dictionary, restore i,j */
	{
	  backtrack = FALSE;

	  if (check_cache (D_p, Q, &i, &j))
	    {
			;
		}
	  else
	    {
	      D->depth--;
	      selectpivot (D, Q, &i, &j);
	      pivot (D, Q, i, j);
	      update (D, Q, &i, &j);	/*Update B,C,i,j */
	    }
	  j++;			/* go to next column */
	}			/* end of if backtrack  */
      if (D->depth < Q->mindepth) 		{ break;	}	/*return new dictionary */

      /* try to go down tree */
/* 2011.7.14 patch */
      while ((j < d) && (!reverse (D, Q, &i, j) || (Q->truncate && Q->minratio[D->m]==1)))
	j++;
      if (j == d )
	backtrack = TRUE;
      else
	/*reverse pivot found */
	{
	  cache_dict (D_p, Q, i, j);
	  /* Note that the next two lines must come _after_ the 
	     call to cache_dict */

	  D->depth++;
	  if (D->depth > Q->deepest)
	    Q->deepest++;

	  pivot (D, Q, i, j);
	  update (D, Q, &i, &j);	/*Update B,C,i,j */

	  D->lexflag = lexmin (D, Q, ZERO);	/* see if lexmin basis */
	  Q->count[2]++;
	  Q->totalnodes++;

	  save_basis (*D_p, Q);
	  if (Q->strace == Q->count[2])
	    Q->debug = TRUE;
	  if (Q->etrace == Q->count[2])
	    Q->debug = FALSE;
		{ return TRUE;	}	/*return new dictionary */
	}


    }				/* end of main while loop for getnextbasis */
  return FALSE;			/* done, no more bases */
}				/*end of lrs_getnextbasis */

/*************************************/
/* print out one line of output file */
/*************************************/
long 
lrs_getvertex (lrs_dic * P, lrs_dat * Q, lrs_mp_vector output)
/*Print out current vertex if it is lexmin and return it in output */
/* return FALSE if no output generated  */
{
  long i;
  long ind;			/* output index                                  */
  long ired;			/* counts number of redundant columns            */
/* assign local variables to structures */
  long *redundcol = Q->redundcol;
  long *count = Q->count;
  long hull;
  long lexflag;

  hull = Q->hull;
  lexflag = P->lexflag;
  if (lexflag || Q->allbases)
    ++(Q->count[1]);
#ifdef PLRS
        // do not print vertex again in PLRS at root
	if(P->depth == Q->mindepth ){
		return FALSE;
	}

#else
	//If we are at minimum depth and not at root do not print vertex
	if(P->depth == Q->mindepth && Q->mindepth != 0){
		return FALSE;
	}
#endif

  linint (Q->sumdet, 1, P->det, 1);
  if (Q->getvolume)
   {
    updatevolume (P, Q);
    if(Q->verbose)   /* this will print out a triangulation */
	; //lrs_printcobasis(P,Q,ZERO);
   }


  /*print cobasis if printcobasis=TRUE and count[2] a multiple of frequency */
  /* or for lexmin basis, except origin for hull computation - ugly!        */

  if (Q->printcobasis)
    if ((lexflag && !hull)  || ((Q->frequency > 0) && (count[2] == (count[2] / Q->frequency) * Q->frequency)))
	if(P->depth != Q->mindepth || Q->mindepth == 0) //Don't print cobasis if this is a restart cobasis
		; //lrs_printcobasis(P,Q,ZERO);

  if (hull)
    return FALSE;		/* skip printing the origin */

  if (!lexflag && !Q->allbases && !Q->lponly)	/* not lexmin, and not printing forced */
    return FALSE;


  /* copy column 0 to output */

  i = 1;
  ired = 0;
  copy (output[0], P->det);

  for (ind = 1; ind < Q->n; ind++)	/* extract solution */

    if ((ired < Q->nredundcol) && (redundcol[ired] == ind))
      /* column was deleted as redundant */
      {
	itomp (ZERO, output[ind]);
	ired++;
      }
    else
      /* column not deleted as redundant */
      {
	getnextoutput (P, Q, i, ZERO, output[ind]);
	i++;
      }

  reducearray (output, Q->n);
  if (lexflag && one(output[0]))
      ++Q->count[4];               /* integer vertex */

  return TRUE;
}				/* end of lrs_getvertex */

long 
lrs_getray (lrs_dic * P, lrs_dat * Q, long col, long redcol, lrs_mp_vector output)
/*Print out solution in col and return it in output   */
/*redcol =n for ray/facet 0..n-1 for linearity column */
/*hull=1 implies facets will be recovered             */
/* return FALSE if no output generated in column col  */
{
  long i;
  long ind;			/* output index                                  */
  long ired;			/* counts number of redundant columns            */
/* assign local variables to structures */
  long *redundcol = Q->redundcol;
  long *count = Q->count;
  long hull = Q->hull;
  long n = Q->n;

#ifdef PLRS
        // do not print vertex again in PLRS at root
	if(P->depth == Q->mindepth ){
		return FALSE;
	}

#else
	//If we are at minimum depth and not at origin do not print ray
	if(P->depth == Q->mindepth && Q->mindepth != 0){
		return FALSE;
	}
#endif

  if (redcol == n)
    {
      ++count[0];
      if (Q->printcobasis)
	if(P->depth != Q->mindepth || Q->mindepth == 0) //Don't print cobasis if this is a restart cobasis
		; //lrs_printcobasis(P,Q,col);
    }

  i = 1;
  ired = 0;

  for (ind = 0; ind < n; ind++)	/* print solution */
    {
      if (ind == 0 && !hull)	/* must have a ray, set first column to zero */
	itomp (ZERO, output[0]);

      else if ((ired < Q->nredundcol) && (redundcol[ired] == ind))
	/* column was deleted as redundant */
	{
	  if (redcol == ind)	/* true for linearity on this cobasic index */
	    /* we print reduced determinant instead of zero */
	    copy (output[ind], P->det);
	  else
	    itomp (ZERO, output[ind]);
	  ired++;
	}
      else
	/* column not deleted as redundant */
	{
	  getnextoutput (P, Q, i, col, output[ind]);
	  i++;
	}

    }
  reducearray (output, n);
/* printslack for rays: 2006.10.10 */
/* printslack inequality indices  */

  return TRUE;
}				/* end of lrs_getray */


void 
getnextoutput (lrs_dic * P, lrs_dat * Q, long i, long col, lrs_mp out)
      /* get A[B[i][col] and copy to out */
{
  long row;
  long m = P->m;
  long d = P->d;
  long lastdv = Q->lastdv;
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *Row = P->Row;
  long j;

  if (i == d && Q->voronoi)
    return;			/* skip last column if voronoi set */

  row = Row[i];

  if (Q->nonnegative) 	/* if m+i basic get correct value from dictionary          */
                        /* the slack for the inequality m-d+i contains decision    */
                        /* variable x_i. We first see if this is in the basis      */
                        /* otherwise the value of x_i is zero, except for a ray    */
                        /* when it is one (det/det) for the actual column it is in */
    {
	for (j = lastdv+ 1; j <= m; j++)
           {
             if ( Q->inequality[B[j]-lastdv] == m-d+i )
	          {
                    copy (out, A[Row[j]][col]);
                    return;
                  }
           }
/* did not find inequality m-d+i in basis */
      if ( i == col )
            copy(out,P->det);
      else
            itomp(ZERO,out);
	
    }
  else
    copy (out, A[row][col]);

}				/* end of getnextoutput */


/************************/
/*  Estimation function */
/************************/
long  
lrs_estimate (lrs_dic * P, lrs_dat * Q)
		   /*returns estimate of subtree size (no. cobases) from current node    */
		   /*current node is not counted.                   */
		   /*cest[0]rays [1]vertices [2]bases [3]volume     */
                   /*    [4] integer vertices                       */
{

  lrs_mp_vector output;		/* holds one line of output; ray,vertex,facet,linearity */
  lrs_mp Nvol, Dvol;		/* hold volume of current basis */
  long estdepth = 0;		/* depth of basis/vertex in subtree for estimate */
  long i = 0, j = 0, k, nchild, runcount, col;
  double prod = 0.0;
  double cave[] =
  {0.0, 0.0, 0.0, 0.0, 0.0};
  double nvertices, nbases, nrays, nvol, nivertices;
  long rays = 0;
  double newvol = 0.0;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *isave = Q->isave;
  long *jsave = Q->jsave;
  double *cest = Q->cest;
  long d = P->d;
  lrs_alloc_mp(Nvol); lrs_alloc_mp(Dvol);
/* Main Loop of Estimator */

  output = lrs_alloc_mp_vector (Q->n);	/* output holds one line of output from dictionary     */

  for (runcount = 1; runcount <= Q->runs; runcount++)
    {				/* runcount counts number of random probes */
      j = 0;
      nchild = 1;
      prod = 1;
      nvertices = 0.0;
      nbases = 0.0;
      nrays = 0.0;
      nvol = 0.0;
      nivertices =0.0;

      while (nchild != 0)	/* while not finished yet */
	{

	  nchild = 0;
	  while (j < d)
	    {
	      if (reverse (P, Q, &i, j))
		{
		  isave[nchild] = i;
		  jsave[nchild] = j;
		  nchild++;
		}
	      j++;
	    }

	  if (estdepth == 0 && nchild == 0)
	    {
	      cest[0] = cest[0] + rays;		/* may be some rays here */
              lrs_clear_mp(Nvol); lrs_clear_mp(Dvol);
              lrs_clear_mp_vector(output, Q->n);
	      return(0L);		/*subtree is a leaf */
	    }

	  prod = prod * nchild;
	  nbases = nbases + prod;

	  if (nchild > 0)	/*reverse pivot found choose random child */
	    {
	      k = myrandom (Q->seed, nchild);
	      Q->seed = myrandom (Q->seed, 977L);
	      i = isave[k];
	      j = jsave[k];
	      estdepth++;
	      Q->totalnodes++;	/* calculate total number of nodes evaluated */
	      pivot (P, Q, i, j);
	      update (P, Q, &i, &j);	/*Update B,C,i,j */
	      if (lexmin (P, Q, ZERO))	/* see if lexmin basis for vertex */
                {
		  nvertices = nvertices + prod;
                                                    /* integer vertex estimate */
                  if( lrs_getvertex(P,Q,output))
                      { --Q->count[1];
                       if  (one(output[0] ))
                         { --Q->count[4];
                           nivertices = nivertices + prod;
                         }
                      }
                }

	      rays = 0;
	      for (col = 1; col <= d; col++)
		if (negative (A[0][col]) && (lrs_ratio (P, Q, col) == 0) && lexmin (P, Q, col))
		  rays++;
	      nrays = nrays + prod * rays;	/* update ray info */

	      if (Q->getvolume)
		{
		  rescaledet (P, Q, Nvol, Dvol);	/* scales determinant in case input rational */
		  rattodouble (Nvol, Dvol, &newvol);
		  nvol = nvol + newvol * prod;	/* adjusts volume for degree              */
		}
	      j = 0;
	    }
	}
      cave[0] = cave[0] + nrays;
      cave[1] = cave[1] + nvertices;
      cave[2] = cave[2] + nbases;
      cave[3] = cave[3] + nvol;
      cave[4] = cave[4] + nivertices;

/*  backtrack to root and do it again */

      while (estdepth > 0)
	{
	  estdepth = estdepth - 1;
	  selectpivot (P, Q, &i, &j);
	  pivot (P, Q, i, j);
	  update (P, Q, &i, &j);	/*Update B,C,i,j */
	  j++;
	}

    }				/* end of for loop on runcount */

  j=(long) cave[2]/Q->runs;

//2015.2.9   Do not update totals if we do iterative estimating and subtree is too big
  if(Q->subtreesize == 0  || j <= Q->subtreesize )
       for (i = 0; i < 5; i++)
           cest[i] = cave[i] / Q->runs + cest[i];

  
  lrs_clear_mp(Nvol); lrs_clear_mp(Dvol);
  lrs_clear_mp_vector(output, Q->n);
  return(j);
}				/* end of lrs_estimate  */


/*********************************/
/* Internal functions            */
/*********************************/
/* Basic Dictionary functions    */
/******************************* */

long 
reverse (lrs_dic * P, lrs_dat * Q, long *r, long s)
/*  find reverse indices  */
/* TRUE if B[*r] C[s] is a reverse lexicographic pivot */
{
  long i, j, row, col;

/* assign local variables to structures */
  long *B = P->B;
  long *C = P->C;
  lrs_mp_matrix A = P->A;
  long *Row = P->Row;
  long *Col = P->Col;
  long d = P->d;

  col = Col[s];
  if (!negative (A[0][col]))
    {
      Q->minratio[P->m]=0;  /* 2011.7.14 */
      return (FALSE);
    }

  *r = lrs_ratio (P, Q, col);
  if (*r == 0)			/* we have a ray */
    {
      Q->minratio[P->m]=0;  /* 2011.7.14 */
      return (FALSE);
    }

  row = Row[*r];

/* check cost row after "pivot" for smaller leaving index    */
/* ie. j s.t.  A[0][j]*A[row][col] < A[0][col]*A[row][j]     */
/* note both A[row][col] and A[0][col] are negative          */

  for (i = 0; i < d && C[i] < B[*r]; i++)
    if (i != s)
      {
	j = Col[i];
	if (positive (A[0][j]) || negative (A[row][j]))		/*or else sign test fails trivially */
	  if ((!negative (A[0][j]) && !positive (A[row][j])) ||
	      comprod (A[0][j], A[row][col], A[0][col], A[row][j]) == -1)
	    {			/*+ve cost found */
              Q->minratio[P->m]=0;  /* 2011.7.14 */

	      return (FALSE);
	    }
      }
  return (TRUE);
}				/* end of reverse */

long 
selectpivot (lrs_dic * P, lrs_dat * Q, long *r, long *s)
/* select pivot indices using lexicographic rule   */
/* returns TRUE if pivot found else FALSE          */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s] */
{
  long j, col;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Col = P->Col;
  long d = P->d;

  *r = 0;
  *s = d;
  j = 0;

/*find positive cost coef */
  while ((j < d) && (!positive (A[0][Col[j]])))
    j++;

  if (j < d)			/* pivot column found! */
    {
      *s = j;
      col = Col[j];

      /*find min index ratio */
      *r = lrs_ratio (P, Q, col);
      if (*r != 0)
	return (TRUE);		/* unbounded */
    }
  return (FALSE);
}				/* end of selectpivot        */
/******************************************************* */

void 
pivot (lrs_dic * P, lrs_dat * Q, long bas, long cob)
		     /* Qpivot routine for array A              */
		     /* indices bas, cob are for Basis B and CoBasis C    */
		     /* corresponding to row Row[bas] and column       */
		     /* Col[cob]   respectively                       */
{
  long r, s;
  long i, j;
  lrs_mp Ns, Nt, Ars;
/* assign local variables to structures */

  lrs_mp_matrix A = P->A;
  long *Row = P->Row;
  long *Col = P->Col;
  long d, m_A;

  lrs_alloc_mp(Ns); lrs_alloc_mp(Nt); lrs_alloc_mp(Ars);

  d = P->d;
  m_A = P->m_A;
  Q->count[3]++;    /* count the pivot */

  r = Row[bas];
  s = Col[cob];

  copy (Ars, A[r][s]);
  storesign (P->det, sign (Ars));	/*adjust determinant to new sign */

  for (i = 0; i <= m_A; i++)
    if (i != r)
      for (j = 0; j <= d; j++)
	if (j != s)
	  {
/*          A[i][j]=(A[i][j]*Ars-A[i][s]*A[r][j])/P->det; */

	    mulint (A[i][j], Ars, Nt);
	    mulint (A[i][s], A[r][j], Ns);
	    decint (Nt, Ns);
	    exactdivint (Nt, P->det, A[i][j]);
	  }			/* end if j ....  */

  if (sign (Ars) == POS)
    {
      for (j = 0; j <= d; j++)	/* no need to change sign if Ars neg */
	/*   A[r][j]=-A[r][j];              */
	if (!zero (A[r][j]))
	  changesign (A[r][j]);
    }				/* watch out for above "if" when removing this "}" ! */
  else
    for (i = 0; i <= m_A; i++)
      if (!zero (A[i][s]))
	changesign (A[i][s]);

  /*  A[r][s]=P->det;                  */

  copy (A[r][s], P->det);		/* restore old determinant */
  copy (P->det, Ars);
  storesign (P->det, POS);		/* always keep positive determinant */

/* set the new rescaled objective function value */

  mulint (P->det, Q->Lcm[0], P->objden);
  mulint (Q->Gcd[0], A[0][0], P->objnum);

  if (!Q->maximize)
        changesign (P->objnum);
  if (zero (P->objnum))
        storesign (P->objnum, POS);
  else
	reduce (P->objnum,P->objden);

  lrs_clear_mp(Ns); lrs_clear_mp(Nt); lrs_clear_mp(Ars);
}				/* end of pivot */

long 
primalfeasible (lrs_dic * P, lrs_dat * Q)
/* Do dual pivots to get primal feasibility */
/* Note that cost row is all zero, so no ratio test needed for Dual Bland's rule */
{
  long primalinfeasible = TRUE;
  long i, j;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Row = P->Row;
  long *Col = P->Col;
  long m, d, lastdv;
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

/*temporary: try to get new start after linearity */

  while (primalinfeasible)
    {
      i=lastdv+1;
      while (i <= m && !negative (A[Row[i]][0]) )
        i++;
      if (i <= m )
	{
	      j = 0;		/*find a positive entry for in row */
	      while (j < d && !positive (A[Row[i]][Col[j]]))
		j++;
	      if (j >= d)
		return (FALSE);	/* no positive entry */
	      pivot (P, Q, i, j);
	      update (P, Q, &i, &j);
        }
      else
         primalinfeasible = FALSE;
    }				/* end of while primalinfeasibile */
  return (TRUE);
}				/* end of primalfeasible */


long 
lrs_solvelp (lrs_dic * P, lrs_dat * Q, long maximize)
/* Solve primal feasible lp by Dantzig`s rule and lexicographic ratio test */
/* return TRUE if bounded, FALSE if unbounded                              */
{
  long i, j;
/* assign local variables to structures */
  long d = P->d;

  while (dan_selectpivot (P, Q, &i, &j))
    {
      Q->count[3]++;
      pivot (P, Q, i, j);
      update (P, Q, &i, &j);	/*Update B,C,i,j */
    }

  if (j < d && i == 0)		/* selectpivot gives information on unbounded solution */
    {
      return FALSE;
    }
  return TRUE;
}				/* end of lrs_solvelp  */

long 
getabasis (lrs_dic * P, lrs_dat * Q, long order[])
/* Pivot Ax<=b to standard form */
/*Try to find a starting basis by pivoting in the variables x[1]..x[d]        */
/*If there are any input linearities, these appear first in order[]           */
/* Steps: (a) Try to pivot out basic variables using order                    */
/*            Stop if some linearity cannot be made to leave basis            */
/*        (b) Permanently remove the cobasic indices of linearities           */
/*        (c) If some decision variable cobasic, it is a linearity,           */
/*            and will be removed.                                            */
{
  long i, j, k;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long *linearity = Q->linearity;
  long *redundcol = Q->redundcol;
  long m, d, nlinearity;
  long nredundcol = 0L;		/* will be calculated here */
  m = P->m;
  d = P->d;
  nlinearity = Q->nlinearity;

	for (j = 0l; j < m; j++)
	{
		i = 0l;
		while (i <= m && B[i] != d + order[j])
			i++;			/* find leaving basis index i */
		if (j < nlinearity && i > m)	/* cannot pivot linearity to cobasis */
		{
	  		return FALSE;
		}
		if (i <= m)
		{			/* try to do a pivot */
	  		k = 0l;
	  		while (C[k] <= d && zero (A[Row[i]][Col[k]])){
	    			k++;
			}
	  		if (C[k] <= d)
	    		{
				
	      			pivot (P, Q, i, k);
	      			update (P, Q, &i, &k);
    			}
	  		else if (j < nlinearity)
	    		{			/* cannot pivot linearity to cobasis */
				if (zero (A[Row[i]][0]))
				{
			  		linearity[j]=0l;
				}
		      		else
				{
			  		return FALSE;
				}
	    		}
			

		}
	}

/* update linearity array to get rid of redundancies */
  i = 0;
  k = 0;			/* counters for linearities         */
  while (k < nlinearity)
    {
      while (k < nlinearity && linearity[k] == 0)
	k++;
      if (k < nlinearity)
	linearity[i++] = linearity[k++];
    }

  nlinearity = i;
/* bug fix, 2009.6.27 */    Q->nlinearity = i;

/* column dependencies now can be recorded  */
/* redundcol contains input column number 0..n-1 where redundancy is */
  k = 0;
  while (k < d && C[k] <= d)
    {
      if (C[k] <= d){		/* decision variable still in cobasis */
	redundcol[nredundcol++] = C[k] - Q->hull;	/* adjust for hull indices */
		
	}
      k++;
    }

/* now we know how many decision variables remain in problem */
  Q->nredundcol = nredundcol;
  Q->lastdv = d - nredundcol;

/* Remove linearities from cobasis for rest of computation */
/* This is done in order so indexing is not screwed up */

  for (i = 0; i < nlinearity; i++)
    {				/* find cobasic index */
      k = 0;
      while (k < d && C[k] != linearity[i] + d)
	k++;
      if (k >= d)
	{
	  return FALSE;
	}
      if (!removecobasicindex (P, Q, k))
	return FALSE;
      d = P->d;
    }
/* set index value for first slack variable */

/* Check feasability */
  if (Q->givenstart)
    {
      i = Q->lastdv + 1;
      while (i <= m && !negative (A[Row[i]][0]))
	i++;
      if (i <= m); 
	}
  return TRUE;
}				/*  end of getabasis */

long 
removecobasicindex (lrs_dic * P, lrs_dat * Q, long k)
/* remove the variable C[k] from the problem */
/* used after detecting column dependency    */
{
  long i, j, cindex, deloc;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Col = P->Col;
  long m, d;
  m = P->m;
  d = P->d;

  cindex = C[k];		/* cobasic index to remove              */
  deloc = Col[k];		/* matrix column location to remove     */

  for (i = 1; i <= m; i++)	/* reduce basic indices by 1 after index */
    if (B[i] > cindex)
      B[i]--;
  
  for (j = k; j < d; j++)	/* move down other cobasic variables    */
    {
      C[j] = C[j + 1] - 1;	/* cobasic index reduced by 1           */
      Col[j] = Col[j + 1];
    }

  if (deloc != d)               
    {
  /* copy col d to deloc */
      for (i = 0; i <= m; i++)
        copy (A[i][deloc], A[i][d]);

  /* reassign location for moved column */
      j = 0;
      while (Col[j] != d)
        j++;
      
    Col[j] = deloc;
  }

  P->d--;
  return TRUE;
}				/* end of removecobasicindex */

lrs_dic *
resize (lrs_dic * P, lrs_dat * Q)
	/* resize the dictionary after some columns are deleted, ie. inputd>d */
	/* a new lrs_dic record is created with reduced size, and items copied over */
{
  lrs_dic *P1;			/* to hold new dictionary in case of resizing */

  long i, j;
  long m, d, m_A;


  m = P->m;
  d = P->d;
  m_A = P->m_A;

/* get new dictionary record */

  P1 = new_lrs_dic (m, d, m_A);

/* copy data from P to P1    */
  P1->i = P->i;
  P1->j = P->j;
  P1->depth = P->depth;
  P1->m = P->m;
  P1->d = P1->d_orig = d;
  P1->lexflag = P->lexflag;
  P1->m_A = P->m_A;
  copy (P1->det, P->det);
  copy (P1->objnum, P->objnum);
  copy (P1->objden, P->objden);

  for (i = 0; i <= m; i++)
    {
      P1->B[i] = P->B[i];
      P1->Row[i] = P->Row[i];
    }
  for (i = 0; i <= m_A; i++)
    {
      for (j = 0; j <= d; j++)
	copy (P1->A[i][j], P->A[i][j]);
    }


  for (j = 0; j <= d; j++)
    {
      P1->Col[j] = P->Col[j];
      P1->C[j] = P->C[j];
    }

  lrs_free_dic (P,Q);

/* Reassign cache pointers */

  Q->Qhead = P1;
  Q->Qtail = P1;
  P1->next = P1;
  P1->prev = P1;

  return P1;

}
/********* resize                    ***************/


long 
restartpivots (lrs_dic * P, lrs_dat * Q)
/* facet contains a list of the inequalities in the cobasis for the restart */
/* inequality contains the relabelled inequalities after initialization     */
{
  long i, j, k;
  long *Cobasic;		/* when restarting, Cobasic[j]=1 if j is in cobasis */
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long *inequality = Q->inequality;
  long *facet = Q->facet;
  long nlinearity = Q->nlinearity;
  long m, d, lastdv;
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  Cobasic = (long *) CALLOC ((unsigned) m + d + 2, sizeof (long));

  /* set Cobasic flags */
  for (i = 0; i < m + d + 1; i++)
    Cobasic[i] = 0;
  for (i = 0; i < d; i++)	/* find index corresponding to facet[i] */
    {
      j = 1;
      while (facet[i + nlinearity] != inequality[j])
	j++;
      Cobasic[j + lastdv] = 1;
    }

  /* Note that the order of doing the pivots is important, as */
  /* the B and C vectors are reordered after each pivot       */

/* Suggested new code from db starts */
  i=m;
  while (i>d){
    while(Cobasic[B[i]]){
      k = d - 1;
      while ((k >= 0) && (zero (A[Row[i]][Col[k]]) || Cobasic[C[k]])){
       k--;
      }
      if (k >= 0)  {
           /*db asks: should i really be modified here? (see old code) */
           /*da replies: modifying i only makes is larger, and so      */
           /*the second while loop will put it back where it was       */
           /*faster (and safer) as done below                          */
       long  ii=i;
       pivot (P, Q, ii, k);
       update (P, Q, &ii, &k);
      } else {
       free(Cobasic);
       return FALSE;
      }
    }
    i--;
  }
/* Suggested new code from db ends */

  if (lexmin (P, Q, ZERO))
    --Q->count[1];		/* decrement vertex count if lexmin */
/* check restarting from a primal feasible dictionary               */
  for (i = lastdv + 1; i <= m; i++)
    if (negative (A[Row[i]][0]))
      {
        free(Cobasic);
	return FALSE;
      }
  free(Cobasic);
  return TRUE;

}				/* end of restartpivots */


long 
lrs_ratio (lrs_dic * P, lrs_dat * Q, long col)	/*find lex min. ratio */
		  /* find min index ratio -aig/ais, ais<0 */
		  /* if multiple, checks successive basis columns */
		  /* recoded Dec 1997                     */
{
  long i, j, comp, ratiocol, basicindex, start, nstart, cindex, bindex;
  long firstime;		/*For ratio test, true on first pass,else false */
  lrs_mp Nmin, Dmin;
  long degencount, ndegencount;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *Row = P->Row;
  long *Col = P->Col;
  long *minratio = Q->minratio;
  long m, d, lastdv;

  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;


  nstart=0;
  ndegencount=0;
  degencount = 0;
  minratio[P->m]=1;   /*2011.7.14 non-degenerate pivot flag */

  for (j = lastdv + 1; j <= m; j++)
    {
      /* search rows with negative coefficient in dictionary */
      /*  minratio contains indices of min ratio cols        */
      if (negative (A[Row[j]][col]))
       {
	minratio[degencount++] = j;
        if(zero (A[Row[j]][0]))
          minratio[P->m]=0;   /*2011.7.14 degenerate pivot flag */
       }
    }				/* end of for loop */
  if (degencount == 0)
    return (degencount);	/* non-negative pivot column */

  lrs_alloc_mp(Nmin); lrs_alloc_mp(Dmin);
  ratiocol = 0;			/* column being checked, initially rhs */
  start = 0;			/* starting location in minratio array */
  bindex = d + 1;		/* index of next basic variable to consider */
  cindex = 0;			/* index of next cobasic variable to consider */
  basicindex = d;		/* index of basis inverse for current ratio test, except d=rhs test */
  while (degencount > 1)	/*keep going until unique min ratio found */
    {
      if (B[bindex] == basicindex)	/* identity col in basis inverse */
	{
	  if (minratio[start] == bindex)
	    /* remove this index, all others stay */
	    {
	      start++;
	      degencount--;
	    }
	  bindex++;
	}
      else
	/* perform ratio test on rhs or column of basis inverse */
	{
	  firstime = TRUE;
	  /*get next ratio column and increment cindex */
	  if (basicindex != d)
	    ratiocol = Col[cindex++];
	  for (j = start; j < start + degencount; j++)
	    {
	      i = Row[minratio[j]];	/* i is the row location of the next basic variable */
	      comp = 1;		/* 1:  lhs>rhs;  0:lhs=rhs; -1: lhs<rhs */
	      if (firstime)
		firstime = FALSE;	/*force new min ratio on first time */
	      else
		{
		  if (positive (Nmin) || negative (A[i][ratiocol]))
		    {
		      if (negative (Nmin) || positive (A[i][ratiocol]))
			comp = comprod (Nmin, A[i][col], A[i][ratiocol], Dmin);
		      else
			comp = -1;
		    }

		  else if (zero (Nmin) && zero (A[i][ratiocol]))
		    comp = 0;

		  if (ratiocol == ZERO)
		    comp = -comp;	/* all signs reversed for rhs */
		}
	      if (comp == 1)
		{		/*new minimum ratio */
		  nstart = j;
		  copy (Nmin, A[i][ratiocol]);
		  copy (Dmin, A[i][col]);
		  ndegencount = 1;
		}
	      else if (comp == 0)	/* repeated minimum */
		minratio[nstart + ndegencount++] = minratio[j];

	    }			/* end of  for (j=start.... */
	  degencount = ndegencount;
	  start = nstart;
	}			/* end of else perform ratio test statement */
      basicindex++;		/* increment column of basis inverse to check next */
    }				/*end of while loop */
  lrs_clear_mp(Nmin); lrs_clear_mp(Dmin);
  return (minratio[start]);
}				/* end of ratio */



long 
lexmin (lrs_dic * P, lrs_dat * Q, long col)
  /*test if basis is lex-min for vertex or ray, if so TRUE */
  /* FALSE if a_r,g=0, a_rs !=0, r > s          */
{
/*do lexmin test for vertex if col=0, otherwise for ray */
  long r, s, i, j;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long m = P->m;
  long d = P->d;

  for (i = Q->lastdv + 1; i <= m; i++)
    {
      r = Row[i];
      if (zero (A[r][col]))	/* necessary for lexmin to fail */
	for (j = 0; j < d; j++)
	  {
	    s = Col[j];
	    if (B[i] > C[j])	/* possible pivot to reduce basis */
	      {
		if (zero (A[r][0]))	/* no need for ratio test, any pivot feasible */
		  {
		    if (!zero (A[r][s]))
		      return (FALSE);
		  }
		else if (negative (A[r][s]) && ismin (P, Q, r, s))
		  {
		    return (FALSE);
		  }
	      }			/* end of if B[i] ... */
	  }
    }
  return (TRUE);
}				/* end of lexmin */

long 
ismin (lrs_dic * P, lrs_dat * Q, long r, long s)
/*test if A[r][s] is a min ratio for col s */
{
  long i;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long m_A = P->m_A;

  for (i = 1; i <= m_A; i++)
    if ((i != r) &&
	negative (A[i][s]) && comprod (A[i][0], A[r][s], A[i][s], A[r][0]))
      {
	return (FALSE);
      }

  return (TRUE);
}

void 
update (lrs_dic * P, lrs_dat * Q, long *i, long *j)
 /*update the B,C arrays after a pivot */
 /*   involving B[bas] and C[cob]           */
{

  long leave, enter;
/* assign local variables to structures */
  long *B = P->B;
  long *C = P->C;
  long *Row = P->Row;
  long *Col = P->Col;
  long m = P->m;
  long d = P->d;

  leave = B[*i];
  enter = C[*j];
  B[*i] = enter;
  reorder1 (B, Row, *i, m + 1);
  C[*j] = leave;
  reorder1 (C, Col, *j, d);
/* restore i and j to new positions in basis */
  for (*i = 1; B[*i] != enter; (*i)++);		/*Find basis index */
  for (*j = 0; C[*j] != leave; (*j)++);		/*Find co-basis index */
}				/* end of update */

long 
lrs_degenerate (lrs_dic * P, lrs_dat * Q)
/* TRUE if the current dictionary is primal degenerate */
/* not thoroughly tested   2000/02/15                  */
{
  long i;
  long *Row;

  lrs_mp_matrix A = P->A;
  long d = P->d;
  long m = P->m;

  Row = P->Row;

  for (i = d + 1; i <= m; i++)
    if (zero (A[Row[i]][0]))
      return TRUE;

  return FALSE;
}


/*********************************************************/
/*                 Miscellaneous                         */
/******************************************************* */

void 
reorder (long a[], long range)
/*reorder array in increasing order with one misplaced element */
{
  long i, temp;
  for (i = 0; i < range - 1; i++)
    if (a[i] > a[i + 1])
      {
	temp = a[i];
	a[i] = a[i + 1];
	a[i + 1] = temp;
      }
  for (i = range - 2; i >= 0; i--)
    if (a[i] > a[i + 1])
      {
	temp = a[i];
	a[i] = a[i + 1];
	a[i + 1] = temp;
      }

}				/* end of reorder */

void 
reorder1 (long a[], long b[], long newone, long range)
/*reorder array a in increasing order with one misplaced element at index newone */
/*elements of array b are updated to stay aligned with a */
{
  long temp;
  while (newone > 0 && a[newone] < a[newone - 1])
    {
      temp = a[newone];
      a[newone] = a[newone - 1];
      a[newone - 1] = temp;
      temp = b[newone];
      b[newone] = b[newone - 1];
      b[--newone] = temp;
    }
  while (newone < range - 1 && a[newone] > a[newone + 1])
    {
      temp = a[newone];
      a[newone] = a[newone + 1];
      a[newone + 1] = temp;
      temp = b[newone];
      b[newone] = b[newone + 1];
      b[++newone] = temp;
    }
}				/* end of reorder1 */


void 
rescaledet (lrs_dic * P, lrs_dat * Q, lrs_mp Vnum, lrs_mp Vden)
  /* rescale determinant to get its volume */
  /* Vnum/Vden is volume of current basis  */
{
  lrs_mp gcdprod;		/* to hold scale factors */
  long i;
/* assign local variables to structures */
  long *B = P->B;
  long *C = P->C;
  long m, d, lastdv;

  lrs_alloc_mp(gcdprod);
  m = P->m;
  d = P->d;
  lastdv = Q->lastdv;

  itomp (ONE, gcdprod);
  itomp (ONE, Vden);
  for (i = 0; i < d; i++)
    if (B[i] <= m)
      {
	mulint (Q->Gcd[Q->inequality[C[i] - lastdv]], gcdprod, gcdprod);
	mulint (Q->Lcm[Q->inequality[C[i] - lastdv]], Vden, Vden);
      }
  mulint (P->det, gcdprod, Vnum);
  reduce (Vnum, Vden);
  lrs_clear_mp(gcdprod);
}				/* end rescaledet */

void 
rescalevolume (lrs_dic * P, lrs_dat * Q, lrs_mp Vnum, lrs_mp Vden)
/* adjust volume for dimension */
{
  lrs_mp temp, dfactorial;
/* assign local variables to structures */
  long lastdv = Q->lastdv;

  lrs_alloc_mp(temp); lrs_alloc_mp(dfactorial);

  /*reduce Vnum by d factorial  */
  getfactorial (dfactorial, lastdv);
  mulint (dfactorial, Vden, Vden);
  if (Q->hull && !Q->homogeneous)
    {				/* For hull option multiply by d to correct for lifting */
      itomp (lastdv, temp);
      mulint (temp, Vnum, Vnum);
    }

  reduce (Vnum, Vden);
  lrs_clear_mp(temp); lrs_clear_mp(dfactorial);
}


void 
updatevolume (lrs_dic * P, lrs_dat * Q)		/* rescale determinant and update the volume */
{
  lrs_mp tN, tD, Vnum, Vden;
  lrs_alloc_mp(tN); lrs_alloc_mp(tD); lrs_alloc_mp(Vnum); lrs_alloc_mp(Vden);
  rescaledet (P, Q, Vnum, Vden);
  copy (tN, Q->Nvolume);
  copy (tD, Q->Dvolume);
  linrat (tN, tD, ONE, Vnum, Vden, ONE, Q->Nvolume, Q->Dvolume);
  lrs_clear_mp(tN); lrs_clear_mp(tD); lrs_clear_mp(Vnum); lrs_clear_mp(Vden);

}				/* end of updatevolume */



/***************************************************/
/* Routines for redundancy checking                */
/***************************************************/

long 
checkredund (lrs_dic * P, lrs_dat * Q)
/* Solve primal feasible lp by least subscript and lex min basis method */
/* to check redundancy of a row in objective function                   */
/* returns TRUE if redundant, else FALSE                                */
{
  lrs_mp Ns, Nt;
  long i, j;
  long r, s;

/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Row, *Col;
  long d = P->d;

  lrs_alloc_mp(Ns); lrs_alloc_mp(Nt);
  Row = P->Row;
  Col = P->Col;

  while (selectpivot (P, Q, &i, &j))
    {
      Q->count[2]++;

/* sign of new value of A[0][0]            */
/* is      A[0][s]*A[r][0]-A[0][0]*A[r][s] */

      r = Row[i];
      s = Col[j];

      mulint (A[0][s], A[r][0], Ns);
      mulint (A[0][0], A[r][s], Nt);

      if (mp_greater (Ns, Nt))
        {
          lrs_clear_mp(Ns); lrs_clear_mp(Nt);
	  return FALSE;		/* non-redundant */
        }

      pivot (P, Q, i, j);
      update (P, Q, &i, &j);	/*Update B,C,i,j */

    }

  lrs_clear_mp(Ns); lrs_clear_mp(Nt);

  return !(j < d && i == 0);	/* unbounded is also non-redundant */


}				/* end of checkredund  */

long 
checkcobasic (lrs_dic * P, lrs_dat * Q, long index)
/* TRUE if index is cobasic and nonredundant                         */
/* FALSE if basic, or degen. cobasic, where it will get pivoted out  */

{

/* assign local variables to structures */

  lrs_mp_matrix A = P->A;
  long *C, *Row, *Col;
  long d = P->d;
  long m = P->m;
  long i = 0;
  long j = 0;
  long s;

  C = P->C;
  Row = P->Row;
  Col = P->Col;

  while ((j < d) && C[j] != index)
    j++;

  if (j == d)
    return FALSE;		/* not cobasic index */

/* index is cobasic */

/* not debugged for new LOC 
   s=LOC[index];
 */
  s = Col[j];
  i = Q->lastdv + 1;

  while ((i <= m) &&
	 (zero (A[Row[i]][s]) || !zero (A[Row[i]][0])))
    i++;

  if (i > m)
    {
      return TRUE;
    }

  pivot (P, Q, i, j);
  update (P, Q, &i, &j);	/*Update B,C,i,j */

  return FALSE;			/*index is no longer cobasic */
}				/* end of checkcobasic */

long 
checkindex (lrs_dic * P, lrs_dat * Q, long index)
/* 0 if index is non-redundant inequality */
/* 1 if index is redundant     inequality */
/* 2 if index is input linearity          */
/*NOTE: row is returned all zero if redundant!! */
{
  long i, j;

  lrs_mp_matrix A = P->A;
  long *Row = P->Row;
  long *B = P->B;
  long d = P->d;
  long m = P->m;

/* each slack index must be checked for redundancy */
/* if in cobasis, it is pivoted out if degenerate */
/* else it is non-redundant                       */

  if (checkcobasic (P, Q, index))
    return ZERO;

/* index is basic   */
/* not debugged for new LOC 
   i=LOC[index];
 */
  j = 1;
  while ((j <= m) && (B[j] != index))
    j++;

  i = Row[j];

  /* copy row i to cost row, and set it to zero */

  for (j = 0; j <= d; j++)
    {
      copy (A[0][j], A[i][j]);
      changesign (A[0][j]);
      itomp (ZERO, A[i][j]);
    }


  if (checkredund (P, Q))
    return ONE;

/* non-redundant, copy back and change sign */

  for (j = 0; j <= d; j++)
    {
      copy (A[i][j], A[0][j]);
      changesign (A[i][j]);
    }

  return ZERO;

}				/* end of checkindex */


/***************************************************************/
/*                                                             */
/*     Routines for caching, allocating etc.                   */
/*                                                             */
/***************************************************************/

/* From here mostly Bremner's handiwork */

static void
cache_dict (lrs_dic ** D_p, lrs_dat * global, long i, long j)
{

  if (dict_limit > 1)
    {
      /* save row, column indicies */
      (*D_p)->i = i;
      (*D_p)->j = j;

/* Make a new, blank spot at the end of the queue to copy into */ 

      pushQ (global, (*D_p)->m, (*D_p)->d, (*D_p)->m_A);

      copy_dict (global, global->Qtail, *D_p);	/* Copy current dictionary */
    }
  *D_p = global->Qtail;

}

void 
copy_dict (lrs_dat * global, lrs_dic * dest, lrs_dic * src)
{
  long m = src->m;
  long m_A = src->m_A;        /* number of rows in A */
  long d = src->d;
  long r,s;

#ifdef GMP 
  for ( r=0;r<=m_A;r++)
    for( s=0;s<=d;s++)
       copy(dest->A[r][s],src->A[r][s]);

#else            /* fast copy for MP and LRSLONG arithmetic */
  /* Note that the "A" pointer trees need not be copied, since they
     always point to the same places within the corresponding space
*/
/* I wish I understood the above remark. For the time being, do it the easy way for Nash */
  if(global->nash)
  {
  for ( r=0;r<=m_A;r++)
    for( s=0;s<=d;s++)
       copy(dest->A[r][s],src->A[r][s]);
  }
  else
  memcpy (dest->A[0][0], (global->Qtail->prev)->A[0][0],
          (d + 1) * (lrs_digits + 1) * (m_A + 1) * sizeof (long));

#endif

  dest->i = src->i;
  dest->j = src->j;
  dest->m = m;
  dest->d = d;
  dest->m_A  = src->m_A;

  dest->depth = src->depth;
  dest->lexflag = src->lexflag;

  copy (dest->det, src->det);
  copy (dest->objnum, src->objnum);
  copy (dest->objden, src->objden);

  memcpy (dest->B, src->B, (m + 1) * sizeof (long));
  memcpy (dest->C, src->C, (d + 1) * sizeof (long));
  memcpy (dest->Row, src->Row, (m + 1) * sizeof (long));
  memcpy (dest->Col, src->Col, (d + 1) * sizeof (long));
}

/* 
 * pushQ(lrs_dat *globals,m,d):
 * this routine ensures that Qtail points to a record that 
 * may be copied into.
 *
 * It may create a new record, or it may just move the head pointer
 * forward so that know that the old record has been overwritten.
 */

static void
pushQ (lrs_dat * global, long m, long d ,long m_A)
{

  if ((global->Qtail->next) == global->Qhead)
    {
      /* the Queue is full */
      if (dict_count < dict_limit)
	{
	  /* but we are allowed to create more */
	  lrs_dic *p;

	  p = new_lrs_dic (m, d, m_A);

	  if (p)
	    {
	      /* we successfully created another record */

	      p->next = global->Qtail->next;
	      (global->Qtail->next)->prev = p;
	      (global->Qtail->next) = p;
	      p->prev = global->Qtail;

	      dict_count++;
	      global->Qtail = p;
	    }
	  else
	    {
	      /* virtual memory exhausted. bummer */
	      global->Qhead = global->Qhead->next;
	      global->Qtail = global->Qtail->next;
	    }
	}
      else
	{
	  /*
	   * user defined limit reached. start overwriting the
	   * beginning of Q 
	   */
	  global->Qhead = global->Qhead->next;
	  global->Qtail = global->Qtail->next;
	}
    }

  else
    {
      global->Qtail = global->Qtail->next;
    }
}

lrs_dic *
lrs_getdic(lrs_dat *Q)
/* create another dictionary for Q without copying any values */
/* derived from lrs_alloc_dic,  used by nash.c                */
{
lrs_dic *p;

  long m;

  m = Q->m;

/* nonnegative flag set means that problem is d rows "bigger"     */
/* since nonnegative constraints are not kept explicitly          */


  if(Q->nonnegative)
    m = m+Q->inputd;

  p = new_lrs_dic (m, Q->inputd, Q->m);
  if (!p)
    return NULL;

  p->next = p;
  p->prev = p;
  Q->Qhead = p;
  Q->Qtail = p;

  return p;
}


#define NULLRETURN(e) if (!(e)) return NULL;

static lrs_dic *
new_lrs_dic (long m, long d, long m_A)
{
  lrs_dic *p;

  NULLRETURN (p = (lrs_dic *) malloc (sizeof (lrs_dic)));


  NULLRETURN (p->B = (long int*) calloc ((m + 1), sizeof (long)));
  NULLRETURN (p->Row = (long int*) calloc ((m + 1), sizeof (long)));

  NULLRETURN (p->C =  (long int*) calloc ((d + 1), sizeof (long)));
  NULLRETURN (p->Col = (long int*) calloc ((d + 1), sizeof (long)));

#ifdef GMP
  lrs_alloc_mp(p->det);
  lrs_alloc_mp(p->objnum);
  lrs_alloc_mp(p->objden);
#endif

  p->d_orig=d;
  p->A=lrs_alloc_mp_matrix(m_A,d);


  return p;
}

void 
lrs_free_dic (lrs_dic * P, lrs_dat *Q)
{
/* do the same steps as for allocation, but backwards */
/* gmp variables cannot be cleared using free: use lrs_clear_mp* */
  lrs_dic *P1;

/* repeat until cache is empty */

  do
  {
    /* I moved these here because I'm not certain the cached dictionaries
       need to be the same size. Well, it doesn't cost anything to be safe. db */

  long d = P->d_orig;
  long m_A = P->m_A;

  lrs_clear_mp_matrix (P->A,m_A,d);

/* "it is a ghastly error to free something not assigned my malloc" KR167 */
/* so don't try: free (P->det);                                           */

  lrs_clear_mp (P->det);      
  lrs_clear_mp (P->objnum);      
  lrs_clear_mp (P->objden);      

  free (P->Row);
  free (P->Col);
  free (P->C);
  free (P->B);

/* go to next record in cache if any */
  P1 =P->next;
  free (P);
  P=P1;

  }  while (Q->Qhead != P );
}

void
lrs_free_dat ( lrs_dat *Q )
{
   long m=Q->m;

/* most of these items were allocated in lrs_alloc_dic */

  lrs_clear_mp_vector (Q->Gcd,m);
  lrs_clear_mp_vector (Q->Lcm,m);

  lrs_clear_mp (Q->sumdet);
  lrs_clear_mp (Q->Nvolume);
  lrs_clear_mp (Q->Dvolume);
  lrs_clear_mp (Q->saved_det);
  lrs_clear_mp (Q->boundd);
  lrs_clear_mp (Q->boundn);

  free (Q->inequality);
  free (Q->facet);
  free (Q->redundcol);
  free (Q->linearity);
  free (Q->minratio);
  free (Q->temparray);

  free (Q->name);  
  free (Q->saved_C);

  lrs_global_count--;

  free(Q);
}


static long
check_cache (lrs_dic ** D_p, lrs_dat * global, long *i_p, long *j_p)
{
/* assign local variables to structures */



  cache_tries++;

  if (global->Qtail == global->Qhead)
    {
      /* Q has only one element */
      cache_misses++;
      return 0;

    }
  else
    {
      global->Qtail = global->Qtail->prev;

      *D_p = global->Qtail;

      *i_p = global->Qtail->i;
      *j_p = global->Qtail->j;
      return 1;
    }
}


lrs_dic *
lrs_alloc_dic (lrs_dat * Q)
/* allocate and initialize lrs_dic */
{

  lrs_dic *p;
  long i, j;
  long m, d, m_A;

  if (Q->hull)                       /* d=col dimension of A */
    Q->inputd = Q->n;                /* extra column for hull */
  else
    Q->inputd = Q->n - 1;

  m = Q->m;
  d = Q->inputd;
  m_A = m;   /* number of rows in A */

/* nonnegative flag set means that problem is d rows "bigger"     */
/* since nonnegative constraints are not kept explicitly          */

  if(Q->nonnegative)
    m = m+d;

  p = new_lrs_dic (m, d, m_A);
  if (!p)
    return NULL;

  p->next = p;
  p->prev = p;
  Q->Qhead = p;
  Q->Qtail = p;


  dict_count = 1;
  dict_limit = 50; 
  cache_tries = 0;
  cache_misses = 0;

/* Initializations */

  p->d = p->d_orig = d;
  p->m = m;
  p->m_A  = m_A;
  p->depth = 0L;
  p->lexflag = TRUE;
  itomp (ONE, p->det);
  itomp (ZERO, p->objnum);
  itomp (ONE, p->objden);

/*m+d+1 is the number of variables, labelled 0,1,2,...,m+d  */
/*  initialize array to zero   */
  for (i = 0; i <= m_A; i++)
    for (j = 0; j <= d; j++)
      itomp (ZERO, p->A[i][j]);

  Q->inequality = (long int*) CALLOC ((m + 1), sizeof (long));
  if (Q->nlinearity == ZERO)   /* linearity may already be allocated */
      Q->linearity  = (long int*) CALLOC ((m + 1), sizeof (long));

  Q->facet =  (long int*) CALLOC ((unsigned) d + 1, sizeof (long));
  Q->redundcol = (long int*) CALLOC ((d + 1), sizeof (long));
  Q->minratio = (long int*) CALLOC ((m + 1), sizeof (long));
                         /*  2011.7.14  minratio[m]=0 for degen =1 for nondegen pivot*/
  Q->temparray = (long int*) CALLOC ((unsigned) d + 1, sizeof (long));

  Q->inequality[0] = 2L;
  Q->Gcd = lrs_alloc_mp_vector(m);
  Q->Lcm = lrs_alloc_mp_vector(m);
  Q->saved_C = (long int*) CALLOC (d + 1, sizeof (long));

  Q->lastdv = d;      /* last decision variable may be decreased */
                      /* if there are redundant columns          */

/*initialize basis and co-basis indices, and row col locations */
/*if nonnegative, we label differently to avoid initial pivots */
/* set basic indices and rows */
 if(Q->nonnegative)
  for (i = 0; i <= m; i++)
    {
      p->B[i] = i;
      if (i <= d )
          p->Row[i]=0; /* no row for decision variables */
      else 
          p->Row[i]=i-d;
    }
 else
   for (i = 0; i <= m; i++)
    {
      if (i == 0 )
          p->B[0]=0;
      else
          p->B[i] = d + i;
      p->Row[i] = i;
    }

  for (j = 0; j < d; j++)
    {
      if(Q->nonnegative)
          p->C[j] = m+j+1;
      else
          p->C[j] = j + 1;
      p->Col[j] = j + 1;
    }
  p->C[d] = m + d + 1;
  p->Col[d] = 0;
  return p;
}				/* end of lrs_alloc_dic */


/* 
   this routine makes a copy of the information needed to restart, 
   so that we can guarantee that if a signal is received, we 
   can guarantee that nobody is messing with it.
   This as opposed to adding all kinds of critical regions in 
   the main line code.

   It is also used to make sure that in case of overflow, we
   have a valid cobasis to restart from.
 */
static void
save_basis (lrs_dic * P, lrs_dat * Q)
{
  int i;
/* assign local variables to structures */
  long *C = P->C;
  long d;
  d = P->d;

  Q->saved_flag = 1;

  for (i = 0; i < 3; i++)
    Q->saved_count[i] = Q->count[i];

  for (i = 0; i < d + 1; i++)
    Q->saved_C[i] = C[i];

  copy (Q->saved_det, P->det);

  Q->saved_d = P->d;
  Q->saved_depth = P->depth;
}

/* digits overflow is a call from lrs_mp package */


/* Routines based on lp_solve */

void
lrs_set_row(lrs_dic *P, lrs_dat *Q, long row, long num[], long den[], long ineq)
/* convert to lrs_mp then call lrs_set_row */
{
 lrs_mp_vector Num, Den;
 long d;
 long j;
 
  d = P->d;

  Num=lrs_alloc_mp_vector(d+1);
  Den=lrs_alloc_mp_vector(d+1);

  for (j=0;j<=d;j++)
  {
  itomp(num[j],Num[j]);
  itomp(den[j],Den[j]);
  }
 
  lrs_set_row_mp(P,Q,row,Num,Den,ineq);

  lrs_clear_mp_vector(Num,d+1);
  lrs_clear_mp_vector(Den,d+1);

}

void
lrs_set_row_mp(lrs_dic *P, lrs_dat *Q, long row, lrs_mp_vector num, lrs_mp_vector den, long ineq)
/* set row of dictionary using num and den arrays for rational input */
/* ineq = 1 (GE)   - ordinary row  */
/*      = 0 (EQ)   - linearity     */
{
  lrs_mp Temp, mpone;
  lrs_mp_vector oD;             /* denominator for row  */

  long i, j;

/* assign local variables to structures */

  lrs_mp_matrix A;
  lrs_mp_vector Gcd, Lcm;
  long hull;
  long m, d;
  lrs_alloc_mp(Temp); lrs_alloc_mp(mpone);
  hull = Q->hull;
  A = P->A;
  m = P->m;
  d = P->d;
  Gcd = Q->Gcd;
  Lcm = Q->Lcm;

  oD = lrs_alloc_mp_vector (d);
  itomp (ONE, mpone);
  itomp (ONE, oD[0]);

  i=row;
  itomp (ONE, Lcm[i]);      /* Lcm of denominators */
  itomp (ZERO, Gcd[i]);     /* Gcd of numerators */
  for (j = hull; j <= d; j++)       /* hull data copied to cols 1..d */
        {
          copy( A[i][j],num[j-hull]);
          copy(oD[j],den[j-hull]);
          if (!one(oD[j]))
            lcm (Lcm[i], oD[j]);      /* update lcm of denominators */
          copy (Temp, A[i][j]);
          gcd (Gcd[i], Temp);   /* update gcd of numerators   */
        }

  if (hull)
        {
          itomp (ZERO, A[i][0]);        /*for hull, we have to append an extra column of zeroes */
          if (!one (A[i][1]) || !one (oD[1]))         /* all rows must have a one in column one */
            Q->polytope = FALSE;
        }
  if (!zero (A[i][hull]))   /* for H-rep, are zero in column 0     */
        Q->homogeneous = FALSE; /* for V-rep, all zero in column 1     */

  storesign (Gcd[i], POS);
  storesign (Lcm[i], POS);
  if (mp_greater (Gcd[i], mpone) || mp_greater (Lcm[i], mpone))
        for (j = 0; j <= d; j++)
          {
            exactdivint (A[i][j], Gcd[i], Temp);        /*reduce numerators by Gcd  */
            mulint (Lcm[i], Temp, Temp);        /*remove denominators */
            exactdivint (Temp, oD[j], A[i][j]);       /*reduce by former denominator */
          }

  if ( ineq == EQ )        /* input is linearity */
     {
      Q->linearity[Q->nlinearity]=row;
      Q->nlinearity++;
     }

/* 2010.4.26   Set Gcd and Lcm for the non-existant rows when nonnegative set */


  if(Q->nonnegative && row==m)
      for(j=1;j<=d;j++)
         { itomp (ONE, Lcm[m+j]);
           itomp (ONE, Gcd[m+j]);
         }


  lrs_clear_mp_vector (oD,d);
  lrs_clear_mp(Temp); lrs_clear_mp(mpone);
}          /* end of lrs_set_row_mp */

void 
lrs_set_obj(lrs_dic *P, lrs_dat *Q, long num[], long den[], long max)
{
  long i;

  if (max == MAXIMIZE)
       Q->maximize=TRUE;
  else
       {
       Q->minimize=TRUE;
       for(i=0;i<=P->d;i++)
         num[i]=-num[i];
       }

  lrs_set_row(P,Q,0L,num,den,GE);
}

void
lrs_set_obj_mp(lrs_dic *P, lrs_dat *Q, lrs_mp_vector num, lrs_mp_vector den, long max)
{
  long i;

  if (max == MAXIMIZE)
       Q->maximize=TRUE;
  else
       {
       Q->minimize=TRUE;
       for(i=0;i<=P->d;i++)
         changesign(num[i]);
       }

  lrs_set_row_mp(P,Q,0L,num,den,GE);
}


long
lrs_solve_lp(lrs_dic *P, lrs_dat *Q)
/* user callable function to solve lp only */
{
  lrs_mp_matrix Lin;		/* holds input linearities if any are found             */
  long col;

  Q->lponly = TRUE;

  if (!lrs_getfirstbasis (&P, Q, &Lin, FALSE))
    return FALSE;

/* There may have been column redundancy                */
/* If so the linearity space is obtained and redundant  */
/* columns are removed. User can access linearity space */
/* from lrs_mp_matrix Lin dimensions nredundcol x d+1   */

  for (col = 0; col < Q->nredundcol; col++)	/* print linearity space               */
    lrs_printoutput (Q, Lin[col]);   	        /* Array Lin[][] holds the coeffs.     */

  return TRUE;
} /* end of lrs_solve_lp */


long 
dan_selectpivot (lrs_dic * P, lrs_dat * Q, long *r, long *s)
/* select pivot indices using dantzig simplex method             */
/* largest coefficient with lexicographic rule to avoid cycling  */
/* Bohdan Kaluzny's handiwork                                    */
/* returns TRUE if pivot found else FALSE                        */
/* pivot variables are B[*r] C[*s] in locations Row[*r] Col[*s]  */
{
  long j,k,col;
  lrs_mp coeff;
/* assign local variables to structures */
  lrs_mp_matrix A = P->A;
  long *Col = P->Col;
  long d = P->d;

  lrs_alloc_mp (coeff);
  *r = 0;
  *s = d;
  j = 0;
  k = 0;
 
  itomp(0,coeff);
/*find positive cost coef */
  while (k < d) 
     {
       if(mp_greater(A[0][Col[k]],coeff))
        {
          j = k;
          copy(coeff,A[0][Col[j]]);
	}
      k++;
     }

  if (positive(coeff))			/* pivot column found! */
    {
      *s = j;
      col = Col[j];

      /*find min index ratio */
      *r = lrs_ratio (P, Q, col);
      if (*r != 0)
        {
        lrs_clear_mp(coeff);
	return (TRUE);		/* unbounded */
        }
    }
  lrs_clear_mp(coeff);
  return (FALSE);
}				/* end of dan_selectpivot        */


long
lrs_checkbound(lrs_dic *P, lrs_dat *Q)
{
/* check bound on objective and return TRUE if exceeded */

  if(!Q->bound)
     return FALSE;

  if( Q->maximize && comprod(Q->boundn,P->objden,P->objnum,Q->boundd) == 1 )
       {
         return TRUE;
       }
  if( Q->minimize && comprod(Q->boundn,P->objden,P->objnum,Q->boundd) == -1 )
       {
         return TRUE;
       }
  return FALSE;
}


long
lrs_leaf(lrs_dic *P, lrs_dat *Q)
{
/* check if current dictionary is a leaf of reverse search tree */
  long    col=0;
  long    tmp=0;

  while (col < P->d && !reverse(P,Q,&tmp,col))
                 col++;
  if(col < P->d) 
     return 0;            /* dictionary is not a leaf */
  else
     return 1;
}




long 
readrat (lrs_mp Na, lrs_mp Da)
 /* read a rational or integer and convert to lrs_mp with base BASE */
 /* returns true if denominator is not one                      */
 /* returns 999 if premature end of file                        */
{
  if( lrs_therow == lrs_m ) return(999L);
  itomp( (long)lrs_inv[lrs_thecol*lrs_m+lrs_therow], Na ); 
  int den = lrs_idv[lrs_thecol*lrs_m+lrs_therow]; 
  itomp( (long)lrs_idv[lrs_thecol*lrs_m+lrs_therow], Da ); 
  lrs_thecol++;
  if( lrs_thecol == lrs_n ) { lrs_thecol=0; lrs_therow++; }
  if( den == 1) return(FALSE); 
  else return(TRUE);
}
