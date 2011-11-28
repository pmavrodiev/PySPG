/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Fri Oct 7 2011
    copyright            : (C) 2011 by Pavlin Mavrodiev
    email                : pmavrodiev@ethz.ch
 ***************************************************************************/


#include "base.h"

#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <algorithm>

#include <dranxor.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>


using namespace std;
using namespace CTGlobal;

double *estimates; //store estimates at time t
double *initial_estimates; //store initial estimates

/*read "filename" from base.ct containing the initial estimates*/
void read_initial_estimates() {
  using namespace CTGlobal;
  FILE *inFile=NULL;
  double estimate=0;
  int counter = 0;
  estimates = (double *) calloc(N,sizeof(double));
  if (estimates==NULL)
    perror("read_initial_estimates(): failed to allocate memory to estimates\n");
  initial_estimates = (double *) calloc(N,sizeof(double));
  if (initial_estimates==NULL)
    perror("read_initial_estimates(): failed to allocate memory to initial_estimates\n");
   inFile = fopen(filename.c_str(),"r");
  if (inFile == NULL) perror("read_initial_estimates():");
  
  while ((fscanf(inFile,"%lf",&estimate) == 1) && (counter<N)) {
    estimates[counter] = estimate;
    initial_estimates[counter++] = estimate;
  }
  
  fclose(inFile);
}


void weigth(double **weigthmatrix, double **currentestimates, short term, double **k1, double **k2, double **k3) {
  
  int i,j;
  double w=0,xi,xj;
  int index=0;
  
  for (i=0; i < N; i++) {
    xi = (*currentestimates)[i];
    for (j=i+1; j < N; j++) {
      xj = (*currentestimates)[j];
      /*calculate the weight that i gives to j*/
      switch (term) {
      case 0:
	w= 1.0 / (1+exp(fabs(xj-xi) / (alpha)));
	break;
      case 1:
	w= 1.0 / (1+exp(fabs(xj-xi + 0.5*((*k1)[j]-(*k1)[i])) / (alpha)));	
	break;
      case 2:
	w= 1.0 / (1+exp(fabs(xj-xi+0.5*((*k2)[j]-(*k2)[i])) / (alpha)));	
	break;
      case 3:
	w= 1.0 / (1+exp(fabs(xj-xi+(*k3)[j]-(*k3)[i]) / (alpha)));	
	break;
      default:
	perror("weight(): wrong value for \"term\" parameter, must be one of 0,1,2,3");
      }
      /*cache the weight into the weight matrix*/
      index=((2*(N)-i-1)*i) / 2 - i + j-1; 
      (*weigthmatrix)[index] = w;
    }
  }
}


/*** AGGREGATE REGIMR ***/
void run_aggr() {
  /* Euler: x_i(t+1)=x_i(t)+deltat*delta_i_t */
  double delta_i_t=0,mean=0; 
  /*loop over the whole time period*/
  for (long i=0; i<t; i++){
    /*simulate the dynamics of all agents*/
    for (unsigned j=0; j<N; j++) {
      mean=gsl_stats_mean(estimates,1,N);
      double randnum = rt_rand_gaussian()*0.001;
      delta_i_t = deltat*(alpha*(mean-estimates[j])+beta*(initial_estimates[j]-estimates[j]))+sqrt(2*Dn*deltat)*randnum;
      estimates[j] += delta_i_t;
    }
  }
}

/***FULL INFO REGIME ***/
void run_full() {
  double *weigth_matrix = (double *) calloc(N*(N-1)/2,sizeof(double));
  double *k1 = (double *) calloc(N,sizeof(double));
  double *k2 = (double *) calloc(N,sizeof(double));
  double *k3 = (double *) calloc(N,sizeof(double));
  double *k4 = (double *) calloc(N,sizeof(double));
  double xi,xj,normconst,K1,K2,K3,K4,noise;
  int index;

  if (weigth_matrix == NULL) 
    perror("run_full(): failed to allocate memory to weigth matrix\n");
  if (k1 == NULL) 
    perror("run_full(): failed to allocate memory to weigth matrix\n");  
  if (k2 == NULL) 
    perror("run_full(): failed to allocate memory to weigth matrix\n");  
  if (k3 == NULL) 
    perror("run_full(): failed to allocate memory to weigth matrix\n");
  if (k4 == NULL) 
    perror("run_full(): failed to allocate memory to weigth matrix\n");
  
  for (long ii=0; ii<t; ii++) {
    /*calculate the weigth matrix at level k1*/
    weigth(&weigth_matrix,&estimates,0,&k1,&k2,&k3);
    /*calculate the k1 vector for all agents*/
    for (unsigned i=0; i<N; i++) {
      xi=estimates[i];
      /*calculate the normalization constant for agent i*/
      normconst=0.0;
      for (unsigned j=0; j<N; j++) 
	normconst += 1.0 / (1+exp(fabs(estimates[j]-xi)/ alpha));
      normconst = pow(normconst,-1);
      /*fill the k1 vector with social influence only*/
      K1=0;
      for (unsigned j=0; j<N; j++) {
	if (j != i) {
	  xj=estimates[j];
	  index = ((2*N-i-1)*i) / 2 - i + j-1;
	  if (i>j) index = ((2*N-j-1)*j) / 2 - j + i-1;
	  K1 += weigth_matrix[index]*(xj-xi);
	}
      }
      K1 *= normconst;
      /*add to the k1 vector individual conviction*/
      K1 += beta*(initial_estimates[i]-xi);
      K1 *= deltat;
      /*add to the k1 vector noise*/
      noise = sqrt(2*Dn*deltat)*rt_rand_gaussian()*0.001;
      K1 += noise;
      k1[i] = K1;
    }
    /*done with k1*/
   /*calculate the weigth matrix at level k2*/
    weigth(&weigth_matrix,&estimates,1,&k1,&k2,&k3);
    /*calculate the k2 vector for all agents*/
    for (unsigned i=0; i<N; i++) {
      xi=estimates[i];
      /*calculate the normalization constant for agent i*/
      normconst=0.0;
      for (unsigned j=0; j<N; j++) 
	normconst += 1.0 / (1+exp(fabs(estimates[j]-xi+0.5*(k1[j]-k1[i]))/ alpha));
      normconst = pow(normconst,-1);
      /*fill the k2 vector with social influence only*/
      for (unsigned j=0; j<N; j++) {
	if (j != i) {
	  xj=estimates[j];
	  index = ((2*N-i-1)*i) / 2 - i + j-1;
	  if (i>j) index = ((2*N-j-1)*j) / 2 - j + i-1;
	  k2[i] += weigth_matrix[index]*(estimates[j]-estimates[i]+0.5*(k1[j]-k1[i]));
	}
      }
      k2[i] *= normconst;
      /*add to the k2 vector individual conviction*/
      k2[i] += beta*(initial_estimates[i]-estimates[i]-0.5*k1[i]);
      k2[i] *= deltat;
      /*add to the k2 vector noise*/
      k2[i] += sqrt(2*Dn*deltat)*rt_rand_gaussian()*0.001;
    }
    /*done with k2*/
    
   /*calculate the weigth matrix at level k3*/
    weigth(&weigth_matrix,&estimates,2,&k1,&k2,&k3);
    /*calculate the k3 vector for all agents*/
    for (unsigned i=0; i<N; i++) {
      xi=estimates[i];
      /*calculate the normalization constant for agent i*/
      normconst=0.0;
      for (unsigned j=0; j<N; j++) 
	normconst += 1.0 / (1+exp(fabs(estimates[j]-xi+0.5*(k2[j]-k2[i])) / alpha ));
      normconst = pow(normconst,-1);
      /*fill the k3 vector with social influence only*/
      for (unsigned j=0; j<N; j++) {
	if (j != i) {
	  xj=estimates[j];
	  index = ((2*N-i-1)*i) / 2 - i + j-1;
	  if (i>j) index = ((2*N-j-1)*j) / 2 - j + i-1;
	  k3[i] += weigth_matrix[index]*(estimates[j]-estimates[i]+0.5*(k2[j]-k2[i]));
	}
      }
      k3[i] *= normconst;
      /*add to the k3 vector individual conviction*/
      k3[i] += beta*(initial_estimates[i]-estimates[i]-0.5*k2[i]);
      k3[i] *= deltat;
      /*add to the k3 vector noise*/
      k3[i] += sqrt(2*Dn*deltat)*rt_rand_gaussian()*0.001;
    }
    /*done with k3*/
   /*calculate the weigth matrix at level k4*/
    weigth(&weigth_matrix,&estimates,3,&k1,&k2,&k3);
    /*calculate the k4 vector for all agents*/
    for (unsigned i=0; i<N; i++) {
      xi=estimates[i];
      /*calculate the normalization constant for agent i*/
      normconst=0.0;
      for (unsigned j=0; j<N; j++) 
	normconst += 1.0 / (1+exp(fabs(estimates[j]-xi+k3[j]-k3[i])/alpha));
      normconst = pow(normconst,-1);
      /*fill the k4 vector with social influence only*/
      for (unsigned j=0; j<N; j++) {
	if (j != i) {
	  xj=estimates[j];
	  index = ((2*N-i-1)*i) / 2 - i + j-1;
	  if (i>j) index = ((2*N-j-1)*j) / 2 - j + i-1;
	  k4[i] += weigth_matrix[index]*(estimates[j]-estimates[i]+k3[j]-k3[i]);
	}
      }
      k4[i] *= normconst;
      /*add to the k2 vector individual conviction*/
      k4[i] += beta*(initial_estimates[i]-estimates[i]-k3[i]);
      k4[i] *= deltat;
      /*add to the k2 vector noise*/
      k4[i] += sqrt(2*Dn*deltat)*rt_rand_gaussian()*0.001;
    }
    /*done with k4*/

    /* update estimates*/
    for (unsigned i=0; i<N; i++) {
      estimates[i] += (1.0 / 6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    }

  } //end ii
  free(weigth_matrix);
  free(k1);
  free(k2);
  free(k3);
  free(k4);

}

int main(int argc, char *argv[])
{
 
  initialize_program( argc,  argv );
  read_initial_estimates();
  double collective_error = 0;
  double group_diversity = 0;
  double woc=0;
  
  if (regime == "AGGR")
    run_aggr();
  else if (regime == "FULL")
    run_full();
  else 
    {
      printf("main(): Unrecognized regime %s\n",regime.c_str());
      perror("");
    }
  /*calculate the final collective error and group diversity*/
  /*log of the estimates. how much I miss R's vectorization */
  double *dummy = (double *) malloc(N*sizeof(double));
  for (unsigned j=0; j<N; j++) 
    dummy[j] = log(estimates[j]);

  double mean=gsl_stats_mean(dummy,1,N);
  collective_error = pow(lnTruth-mean,2);
  group_diversity =  gsl_stats_variance_m(dummy,1,N,mean);
  free(dummy);

  /*calculate the final WOC indicator*/
  sort(&estimates[0], &estimates[N]);
  int start = ceil(N/2);
  while (start > 1) {
    if ((estimates[start-1] <= exp(lnTruth)) && (estimates[N-start] >= exp(lnTruth))) {
      woc = (double) start;
      break;
    }	 	
    start--;
  }
  /***************/

  std::cout << collective_error << "\t";
  std::cout << group_diversity << "\t";
  std::cout << woc << std::endl;

  free(estimates);
  free(initial_estimates);  
  return EXIT_SUCCESS;
  
}

================================================================================================
  
  /* loop over the whole time period*/
  for (long i=0; i<t; i++) {
    /*simulate the dynamics of all agents*/
    for (unsigned j=0; j<N; j++) {
      delta_i_t = 
    }
  }
  
  /*calculate the final collective error and group diversity*/
  /*log of the estimates. how much I miss R's vectorization */
  double *dummy = (double *) malloc(N*sizeof(double));
  for (unsigned j=0; j<N; j++) 
    dummy[j] = log(estimates[j]);

  double mean=gsl_stats_mean(dummy,1,N);
  collective_error = pow(lnTruth-mean,2);
  group_diversity =  gsl_stats_variance_m(dummy,1,N,mean);
  free(dummy);

  /*calculate the final WOC indicator*/
  sort(&estimates[0], &estimates[N]);
  int start = ceil(N/2);
  while (start > 1) {
    if ((estimates[start-1] <= exp(lnTruth)) && (estimates[N-start] >= exp(lnTruth))) {
      woc = (double) start;
      break;
    }	 	
    start--;
  }
  /***************/
