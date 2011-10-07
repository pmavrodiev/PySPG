#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>


#define YOUT(T_, I_)  REAL(_yout)[T_ + I_*(T+1)]
double *uppertriangular;


/*compute the weight matrix for each k_i*/
/*
  weightmatrix - pointer to an upper triangular matrix
  currentestimates - SEXP pointer to the estimates at the current time step
  term - the term of the RK4 method, i.e. which k do we want to compute - k0 (i.e. Euler) or k1,k2,k3,k4.
  k1,k2,k3,k4 - pointers to matrices that would store the corresponding k terms
  N - total number of agents
  a - the value of alpha
 */
void weight(double *weightmatrix, SEXP currentestimates, short term,double *k1, double *k2, double *k3, double *k4, int *N,double *a) {
  int i,j;
  double w=0,xi,xj;
  int index=0;
  
  
  //  currentestimates = coerceVector(currentestimates,REALSXP);

  for (i=0; i < *N; i++) {
    xi=REAL(currentestimates)[i];
    for (j=i+1; j < *N; j++) {
      xj=REAL(currentestimates)[j];
      /*calculate the weight that i gives to j*/
      switch (term) {
      case 0:
	w= 1.0 / (1+exp(fabs(xj-xi) / (*a)));
	break;
      case 1:
	w= 1.0 / (1+exp(fabs(xj-xi + 0.5*(k1[j]-k1[i])) / (*a)));	
	break;
      case 2:
	w= 1.0 / (1+exp(fabs(xj-xi+0.5*(k2[j]-k2[i])) / (*a)));	
	break;
      case 3:
	w= 1.0 / (1+exp(fabs(xj-xi+k3[j]-k3[i]) / (*a)));	
	break;
      default:
	error("weight(): wrong value for \"term\" parameter - %d",term);
      }
      /*cache the weight into the weight matrix*/
      index=((2*(*N)-i-1)*i) / 2 - i + j-1; 
      weightmatrix[index] = w;
    }
  }
}

/*
  _y0 - vector of initial estimates
  _N - number of agents
  _h - step size for rk4
  _alpha - social influence
  _beta - individual conviction
  _C - constant
  _Dn - diffusion coefficient for noise
  _yout - matrix to store the solution for all agents. size is (T+1)x(N+1)
  _T - number of time steps
  _ywoc - WOC indicator for the given values of alpha and beta
  _Truth - true value
*/

SEXP myrk4 (SEXP _y0 ,SEXP _N,SEXP _h, SEXP _alpha, SEXP _beta, SEXP _C, SEXP _Dn, SEXP _yout, SEXP _T, SEXP _ywoc ,SEXP _Truth) {
  
  int i,j,t,T,N,index,start;
  double a,b,C,Dn,h,*k1,*k2,*k3,*k4,*sorted_estimates;
  double normconst=0,K1,K2,K3,K4,noise,xi,xj;
  SEXP _y;
  SEXP Rdim;
  /*input check*/
  if (!isNumeric(_C) || !isNumeric(_Dn)) 
    error("one of the required inputs is not a double");

  if (!isMatrix(_yout))
    error("a matrix is required for yout in myrk4");
  
  
  _y0 = coerceVector(_y0,REALSXP);
  _ywoc = coerceVector(_ywoc,REALSXP);
  _yout = coerceVector(_yout,REALSXP);
  _Truth = coerceVector(_Truth,REALSXP);
  _N = coerceVector(_N,INTSXP);
  _T = coerceVector(_T,INTSXP);
  N = INTEGER(_N)[0];
  T = INTEGER(_T)[0];
  a = REAL(_alpha)[0];
  b = REAL(_beta)[0];
  C = REAL(_C)[0];
  Dn = REAL(_Dn)[0];
  h = REAL(_h)[0];

  Rdim=getAttrib(_yout, R_DimSymbol);
  i=INTEGER(Rdim)[0];
  j=INTEGER(Rdim)[1];
  if (i != (T+1))
    error("_yout has wrong number of rows %d. Expected %d\n",i,T+1);
  if (j != (N+1))
    error("_yout has wrong number of cols %d. Expected %d\n",j,N+1);

  _y = allocVector(REALSXP,N);
  copyVector(_y,_y0); /*destination,source*/
  /*allocate memory*/
  uppertriangular = (double *) calloc(N*(N-1)/2,sizeof(double));
  if (uppertriangular == NULL) error("myrk4(): failed to allocate memory of %f bytes to weightmatrix",(N*(N-1)/2)*sizeof(double));
  k1 = (double *) calloc(N,sizeof(double));
  if (k1 == NULL) error("myrk4(): failed to allocate memory of %f bytes to k1",N*sizeof(double));
  k2 = (double *) calloc(N,sizeof(double));
  if (k2 == NULL) error("myrk4(): failed to allocate memory of %f bytes to k2",N*sizeof(double));
  k3 = (double *) calloc(N,sizeof(double));
  if (k3 == NULL) error("myrk4(): failed to allocate memory of %f bytes to k3",N*sizeof(double));
  k4 = (double *) calloc(N,sizeof(double));
  if (k4 == NULL) error("myrk4(): failed to allocate memory of %f bytes to k4",N*sizeof(double));
  sorted_estimates = (double *) calloc(N,sizeof(double));
  if (sorted_estimates == NULL)
    error("myrk4(): failed to allocate memory of %f bytes to sorted_estimates",N*sizeof(double));
  GetRNGstate();
  
  for (t=1; t<=T; t++) {
    /*CALCULATE THE WEIGHT MATRIX AT LEVEL K1*/
    weight(uppertriangular,_y,0,k1,k2,k3,k4,&N,&a);
    /*calculate the k1 vector for all agents*/
    for (i=0; i < N; i++) {
      xi=REAL(_y)[i];
      /*calculate the normalization constant for agent i*/
      normconst = 0.0;
      for (j=0; j < N; j++) 
	normconst += 1.0 / (1+exp(fabs(REAL(_y)[j]-xi) / a));
      normconst = pow(normconst,-1);
      /*fill the k1 vector with social influence only*/
      K1=0;
      for (j=0; j < N; j++) {
	if (j != i) {
	  xj = REAL(_y)[j];
	  index = ((2*N-i-1)*i) / 2 - i + j-1; 
	  if (i>j) index = ((2*N-j-1)*j) / 2 - j + i-1;
	  K1 += (uppertriangular[index]*(xj-xi));//sum up the weights
	}
      }
      K1 *= normconst;
      /*add to the k1 vector individual conviction*/
      K1 += b*(REAL(_y0)[i]-xi);
      K1 *= h;
      /*add to the k1 vector noise*/
      noise = sqrt(2*Dn*h)*rnorm(0,0.001);
      K1 += noise;
      k1[i] = K1;
    }
    
    /********   DONE WITH K1   ***************************/
    /*CALCULATE THE WEIGHT MATRIX AT LEVEL K2*/
    weight(uppertriangular,_y,1,k1,k2,k3,k4,&N,&a);
    /*calculate the k1 vector for all agents*/
    for (i=0; i < N; i++) {
      /*calculate the normalization constant for agent i*/
      normconst = 0.0;
      for (j=0; j < N; j++) 
	normconst += 1.0 / (1+exp(fabs(REAL(_y)[j]-REAL(_y)[i]+0.5*(k1[j]-k1[i])) / a));
      normconst = pow(normconst,-1);
      /*fill the k2 vector with social influence only*/
      for (j=0; j < N; j++) {
	if (j != i) {
	  index=((2*N-i-1)*i) / 2 - i + j-1; 
	  if (i>j) index = ((2*N-j-1)*j) / 2 - j + i-1;
	  
	  k2[i] += uppertriangular[index]*(REAL(_y)[j]-REAL(_y)[i]+0.5*(k1[j]-k1[i])); //sum up the weights
	}
      }
      k2[i] *= normconst;
      /*add to the k2 vector individual conviction*/
      k2[i] += b*(REAL(_y0)[i]-REAL(_y)[i]-0.5*k1[i]);
      k2[i] *= h;
      /*add to the k2 vector noise*/
      k2[i] += sqrt(2*Dn*h)*rnorm(0,0.001);
    }
  /********   DONE WITH K2   ***************************/
  /*CALCULATE THE WEIGHT MATRIX AT LEVEL K3*/
    weight(uppertriangular,_y,2,k1,k2,k3,k4,&N,&a);
    /*calculate the k3 vector for all agents*/
    for (i=0; i < N; i++) {
      /*calculate the normalization constant for agent i*/
      normconst = 0.0;
      for (j=0; j < N; j++) 
  	normconst += 1.0 / (1+exp(fabs(REAL(_y)[j]-REAL(_y)[i]+0.5*(k2[j]-k2[i])) / a));
      normconst = pow(normconst,-1);
      /*fill the k3 vector with social influence only*/
      for (j=0; j < N; j++) {
  	if (j != i) {
  	  index=((2*N-i-1)*i) / 2 - i + j-1; 
  	  if (i>j) index = ((2*N-j-1)*j) / 2 - j + i-1;
  	  
  	  k3[i] += uppertriangular[index]*(REAL(_y)[j]-REAL(_y)[i]+0.5*(k2[j]-k2[i])); //sum up the weights
  	}
      }
      k3[i] *= normconst;
      /*add to the k3 vector individual conviction*/
      k3[i] += b*(REAL(_y0)[i]-REAL(_y)[i]-0.5*k2[i]);
      k3[i] *= h;
      /*add to the k3 vector noise*/
      k3[i] += sqrt(2*Dn*h)*rnorm(0,0.001);  
    }
    /********   DONE WITH K3  ***************************/
    /*CALCULATE THE WEIGHT MATRIX AT LEVEL K4*/
    weight(uppertriangular,_y,3,k1,k2,k3,k4,&N,&a);
    /*calculate the k4 vector for all agents*/
    for (i=0; i < N; i++) {
      /*calculate the normalization constant for agent i*/
      normconst = 0.0;
      for (j=0; j < N; j++) 
  	normconst += 1.0 / (1+exp(fabs(REAL(_y)[j]-REAL(_y)[i]+k3[j]-k3[i]) / a));
      normconst = pow(normconst,-1);
      /*fill the k4 vector with social influence only*/
      for (j=0; j < N; j++) {
  	if (j != i) {
  	  index=((2*N-i-1)*i) / 2 - i + j-1; 
  	  if (i>j) index = ((2*N-j-1)*j) / 2 - j + i-1;
  	  
  	  k4[i] += uppertriangular[index]*(REAL(_y)[j]-REAL(_y)[i]+k3[j]-k3[i]); //sum up the weights
  	}
      }
      k4[i] *= normconst;
      /*add to the k4 vector individual conviction*/
      k4[i] += b*(REAL(_y0)[i]-REAL(_y)[i]-k3[i]);
      k4[i] *= h;
      /*add to the k4 vector noise*/
      k4[i] += sqrt(2*Dn*h)*rnorm(0,0.001);  
    }
    /********   DONE WITH K4  ***************************/
    //update _y
    for (i=0; i<N; i++) {
      REAL(_y)[i] += (1.0 / 6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
      YOUT(t,(i+1)) = REAL(_y)[i];
      if (t==T)
	sorted_estimates[i]=REAL(_y)[i];
    }
    
  } //end for (t=1;t<=T;T++) 

  /*calculate WOC*/
  R_qsort(sorted_estimates,1,N);
  start = ceil(N/2);
  while (start > 1) {
    if ((sorted_estimates[start-1] <= REAL(_Truth)[0]) && (sorted_estimates[N-start] >= REAL(_Truth)[0])) {
      REAL(_ywoc)[0] = (double) start;
      break;
    }	 	
    start--;
  }
  /***************/
  
  free(uppertriangular);
  free(k1);
  free(k2);
  free(k3);
  free(k4);
  free(sorted_estimates);
  
  PutRNGstate();
  return R_NilValue;
}
 






