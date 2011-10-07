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

#include <dranxor.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>


using namespace std;

double *estimates; //store estimates at time t
double *initial_estimates; //store initial estimates

/*read "filename" from base.ct containing the initial estimates*/
void read_initial_estimates() {
  using namespace CTGlobal;
  FILE *inFile=NULL;
  double estimate=0;
  int counter = 0;
  estimates = (double *) malloc(N*sizeof(double));
  initial_estimates = (double *) malloc(N*sizeof(double));
  inFile = fopen(filename.c_str(),"r");
  if (inFile == NULL) perror("read_initial_estimates():");
  
  while ((fscanf(inFile,"%lf",&estimate) == 1) && (counter<N)) {
    estimates[counter] = estimate;
    initial_estimates[counter++] = estimate;
  }
  
  fclose(inFile);
}



int main(int argc, char *argv[])
{
  using namespace CTGlobal;

  initialize_program( argc,  argv );
  read_initial_estimates();
  double collective_error = 0;
  double group_diversity = 0;
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

  /*calculate the final collective error and group diversity*/
  /*log of the estimates. how much I miss R's vectorization */
  double *dummy = (double *) malloc(N*sizeof(double));
  for (unsigned j=0; j<N; j++) 
    dummy[j] = log(estimates[j]);

  mean=gsl_stats_mean(dummy,1,N);
  collective_error = pow(lnTruth-mean,2);
  group_diversity =  gsl_stats_variance_m(dummy,1,N,mean);
  free(dummy);
  
  std::cout << collective_error << "\t";
  std::cout << group_diversity << std::endl;

  free(estimates);
  free(initial_estimates);  
  return EXIT_SUCCESS;
  
}


