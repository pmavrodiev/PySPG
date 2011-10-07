/*
COMPILE WITH:

gcc -std=gnu99 -I/usr/share/R/include -fpic  -O3 -c wisdom_from_R.c -o wisdom_from_R.o -fopenmp -lgomp
gcc -std=gnu99 -shared -o wisdom_from_R.so wisdom_from_R.o -L/usr/lib/R/lib -lR -lgomp -O3
*/
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <omp.h>

/*The C representation of R matrices is a one dimensional vector.
  Example:

  The following R matrix:
   
    A    D
    B    E
    C    F

    will be stored as the one dimensional vector: |A|B|C|D|E|F
*/

#define ESTIMATES(I_,J_) REAL(_estimates)[I_ + J_*N]
#define Z(I_,J_) REAL(_z)[I_ + J_*alpha_length]
#define ZZ(I_,J_) REAL(_zz)[I_ + J_*alpha_length]
#define ZZZ(I_,J_) REAL(_zzz)[I_ + J_*alpha_length]


/*'z' stores the collective error for each alpha, beta combination
  'estimates' is self explanatory
  'zz' stores the group diversity for each alpha, beta combination
  'zzz' stores the WOC indicator for each alpha, beta combination
*/

SEXP compute_dt(SEXP _dt, SEXP _estimates,SEXP _alpha_start,SEXP _alpha_step, SEXP _alpha_length, SEXP _beta_start, SEXP _beta_step, SEXP _beta_length , SEXP _z, SEXP _zz, SEXP _zzz,SEXP _Dn, SEXP _Truth) {

  double alpha = 0;
  double beta  = 0;
  double dt, Dn, Truth;
  int i,j,t,k,p,alpha_length,beta_length;
  double random_number, differential, current_mean,current_diversity;
  int T,N,A,B;
  SEXP Rdim;
  int checked;
  int start;
  double *sortedEstimates = NULL;

  /* INPUT CHECK */
  if (!isNumeric(_dt) || !isNumeric(_alpha_start) || !isNumeric(_alpha_step) || !isNumeric(_alpha_length) || \
      !isNumeric(_beta_start) || !isNumeric(_beta_step) || !isNumeric(_beta_length) || !isNumeric(_Dn) || \
      !isNumeric(_Truth)) 
    error("one of the required inputs is not a double");
  

  if (!isMatrix(_estimates))
    error("a matrix is required for estimates in compute_dt");
  
  if (!isMatrix(_z))
    error("a matrix is required for z in compute_dt");  
  
  if (!isMatrix(_zz))
    error("a matrix is required for zz in compute_dt");

  if (!isMatrix(_zzz))
    error("a matrix is required for zzz in compute_dt");
  
  _estimates = coerceVector(_estimates,REALSXP);
  _z = coerceVector(_z,REALSXP);
  _zz = coerceVector(_zz,REALSXP);
  _zzz = coerceVector(_zzz,REALSXP);
  _alpha_length = coerceVector(_alpha_length,INTSXP);
  _beta_length = coerceVector(_beta_length,INTSXP);

  dt = REAL(_dt)[0];
  Dn = REAL(_Dn)[0];
  Truth = REAL(_Truth)[0];

  alpha_length = INTEGER(_alpha_length)[0];
  beta_length = INTEGER(_beta_length)[0];

  
  Rdim = getAttrib(_z, R_DimSymbol);
  A = INTEGER(Rdim)[0]; 
  B = INTEGER(Rdim)[1];
  if (A != alpha_length)
    error("matrix z has wrong number of rows %d. Expected %d",A,alpha_length);
  if (B != beta_length)
    error("matrix z has wrong number of cols %d. Expected %d",B,beta_length);
  
  Rdim = getAttrib(_zz, R_DimSymbol);
  A = INTEGER(Rdim)[0]; 
  B = INTEGER(Rdim)[1];
  if (A != alpha_length)
    error("matrix zz has wrong number of rows %d. Expected %d",A,alpha_length);
  if (B != beta_length)
    error("matrix zz has wrong number of cols %d. Expected %d",B,beta_length);

  Rdim = getAttrib(_zzz, R_DimSymbol);
  A = INTEGER(Rdim)[0]; 
  B = INTEGER(Rdim)[1];
  if (A != alpha_length)
    error("matrix zzz has wrong number of rows %d. Expected %d",A,alpha_length);
  if (B != beta_length)
    error("matrix zzz has wrong number of cols %d. Expected %d",B,beta_length);
  
  /*******************/ 

  /*get the dimensions of _estimates */
  Rdim = getAttrib(_estimates, R_DimSymbol);
  N = INTEGER(Rdim)[0]; 
  T = INTEGER(Rdim)[1]; //T is already T/deltat

  sortedEstimates = (double *) malloc(N*sizeof(double));
  memset(sortedEstimates,0.0,N*sizeof(double));
  if (sortedEstimates == NULL) Rprintf("Memory cannot be allocated\n");
  /**********************************/
  Rprintf("%d,%d, %f\n",N,T,Truth);

  GetRNGstate();
  for (i=0; i<alpha_length; i++) {
    Rprintf("i=%d(%d)\n",i,alpha_length);
    alpha = REAL(_alpha_start)[0] + i*REAL(_alpha_step)[0];
    for (j=0; j<beta_length; j++) {
      Rprintf("j=%d(%d)\n",j,beta_length);
      beta = REAL(_beta_start)[0] + j*REAL(_beta_step)[0];
      for (t=1; t<T; t++) {
	  current_mean = 0;
	  checked = 0;
#pragma omp parallel for default(none) shared(N,_estimates,alpha,dt,current_mean,checked,Dn,beta,t,T) \
  private(random_number,k,p,differential)
	  	  
	  for (k=0; k<N; k++) {
	    random_number = rnorm(0,0.001);
	    /*calculate current mean*/
	    if (!checked) {
            #pragma omp critical
	      {
		  if (!checked) {
		    for (p=0; p<N; p++) 
		      current_mean += ESTIMATES(p,(t-1));
		    
		    current_mean /= N;		
		    checked = 1;
		  }
	      }
	    }	   
	  
	    differential = alpha*dt*(current_mean-ESTIMATES(k,(t-1))) + beta*dt*(ESTIMATES(k,0)-ESTIMATES(k,(t-1))) + sqrt(2*Dn*dt)*random_number;
	  
	    ESTIMATES(k,t) = ESTIMATES(k,(t-1)) + differential;	 
	    if (ESTIMATES(k,t) < 0) Rprintf("%f, %f\n",random_number, ESTIMATES(k,t));
			    
	  }
      }
      //update z here
      /*calculate the mean of the log of the estimates at the last time step*/
      current_mean = 0;
      //combining the WOC calculation with this for cycle for efficiency
      for (p=0; p<N; p++) {
	current_mean += log(ESTIMATES(p,(T-1)));
	sortedEstimates[p] = ESTIMATES(p,(T-1));
      }
      current_mean /= N;
      Z(i,j) = pow(log(Truth)-current_mean,2);
      /*calculate the group diversity*/
      current_diversity=0;
      for (p=0; p<N; p++)
	current_diversity += pow((log(ESTIMATES(p,(T-1))) - current_mean),2);
      current_diversity /= N;
      ZZ(i,j) = current_diversity;
      /*calculate the WOC*/
      R_qsort(sortedEstimates,1,N);
      start = ceil(N/2);
      while (start > 1) {
	if ((sortedEstimates[start-1] <= Truth) && (sortedEstimates[N-start] >= Truth)) {
	  ZZZ(i,j)= (double) start;
	  break;
	}	 	
	start--;
      }
    }
  }
  PutRNGstate();
  free(sortedEstimates);
  return R_NilValue;
}
