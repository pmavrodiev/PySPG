/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Fri Nov 25 2011
    copyright            : (C) 2011 by Pavlin Mavrodiev
    email                : pmavrodiev@ethz.ch
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ***************************************************************************/


#include "base.h"

#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <unistd.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/timeb.h>
#include <iomanip>

using namespace std;
using namespace CTGlobal;

gsl_rng *r; //random number generator
gsl_vector *estimates = NULL; //store estimates at time t
gsl_vector *ranks = NULL; //store initial estimates
FILE *inFile=NULL; //read in the initial estimates

/*gracefully exit by closing all open file handles, freeing memory, etc.*/
void grace_exit() {
  if (estimates != NULL) gsl_vector_free(estimates);
  if (ranks != NULL) gsl_vector_free(ranks);
  if (inFile != NULL) fclose(inFile);
  if (r != NULL)  gsl_rng_free(r);
  exit(1);
}

/*Compute the ranks of the population in the current time period.
  Modified from http://www.gnu.org/s/gsl/manual/html_node/Computing-the-rank.html*/
void compute_ranks () {
  size_t i;
  size_t n = estimates->size;
  gsl_permutation * perm = gsl_permutation_alloc(n);
  if (perm == NULL) grace_exit();
  gsl_permutation * rank = gsl_permutation_alloc(n);
  if (rank == NULL) {
    gsl_permutation_free(perm);
    grace_exit();
  }
  /*store the estimates relative to the truth , i.e. yi=xi-Truth*/
  gsl_vector *relative_estimates = gsl_vector_calloc(N);
  if (relative_estimates == NULL) {
    gsl_permutation_free (perm);
    gsl_permutation_free (rank);
    grace_exit();
  }
  gsl_vector_memcpy(relative_estimates,estimates);
  gsl_vector_add_constant(relative_estimates,-lnTruth);
  /*probably there is a more elegant way to store abs values in-place*/
  for (i=0; i<N; i++) {
    double el = gsl_vector_get(relative_estimates,i); 
    gsl_vector_set(relative_estimates,i,fabs(el));
  }
  /**************************************************************/

  gsl_sort_vector_index (perm, relative_estimates);
  gsl_permutation_inverse (rank, perm);
  
  for (i = 0; i < n; i++) {
    double vi = gsl_vector_get(relative_estimates, i);
    gsl_vector_set(ranks,i,(rank->data[i]) / (double) N);
  }
  
  gsl_permutation_free (perm);
  gsl_permutation_free (rank);
  gsl_vector_free(relative_estimates);
 }


/*read "filename" from base.ct containing the initial estimates*/
int read_initial_estimates() {
  using namespace CTGlobal;  
  double estimate=0;
  int counter = 0;

  estimates = gsl_vector_calloc(N); //init and set all elements to 0
  if (estimates == NULL) {
    GSL_ERROR("read_initial_estimates(): failed to allocate memory to estimates",GSL_ETOL);
    *f_dynamics_out<<"read_initial_estimates(): failed to allocate memory to estimates"<<endl;
    grace_exit();
  }

  ranks = gsl_vector_calloc(N); //init and set all elements to 0
  if (ranks == NULL) {
    GSL_ERROR("read_initial_estimates(): failed to allocate memory to ranks",GSL_ETOL);
    *f_dynamics_out<<"read_initial_estimates(): failed to allocate memory to ranks"<<endl;
    grace_exit();
  }
  
  inFile = fopen(filename.c_str(),"r");
  if (inFile == NULL) {
    perror(getcwd(NULL, 100));
    grace_exit();
  }
  
  if (gsl_vector_fscanf (inFile,estimates)) {//TODO: handle mismatch between N and the length of the input file
    GSL_ERROR("read_initial_estimates(): failed to read input file", GSL_ETOL);
    *f_dynamics_out<<"read_initial_estimates(): failed to read input file"<<endl;
    grace_exit();
 }
  // for (int i=0; i<N; i++)
  //   gsl_vector_set(estimates,i,opinions[i]);
  
  fclose(inFile);

  return 0;
}


int main(int argc, char *argv[]) {
 
  initialize_program( argc,  argv );
  read_initial_estimates();
  double collective_error = 0.0;
  double group_diversity = 0.0;
  double woc=0;
  /* Euler: x_i(t+1)=x_i(t)+deltat*delta_i_t */
  double delta_i_t=0.0,randnum=0.0,D=0.0,mean=0.0;
  struct timeb tp;

  r = gsl_rng_alloc(gsl_rng_taus);
  if (r == NULL) {
    cerr<<"main(): Not enough memory to initialize random generator"<<endl;
    *f_dynamics_out<<"main(): Not enough memory to initialize random generator"<<endl;
    return EXIT_FAILURE;
  }
  ftime(&tp);
  //seed the rng with different value for each realization
  gsl_rng_set(r,tp.time+tp.millitm); 
  
  double *dummy = (double *) malloc(N*sizeof(double));
  if (dummy == NULL) {
    *f_dynamics_out<<"main(): Not enough memory to initialize dummy"<<endl;
    grace_exit();
  } 
  
  /*estimates of agent with rank 1 in the beginning*/
  int rank_one=gsl_vector_min_index(ranks);

  /*compute initial collective error and group diversity*/
  /* LINEAR ESTIMATES CONSIDERED AT THIS POINT */
  for (unsigned j=0; j<N; j++) 
    //dummy[j] = log(gsl_vector_get(estimates,j));   
    dummy[j] = gsl_vector_get(estimates,j);   
  mean=gsl_stats_mean(dummy,1,N);
  collective_error = pow(lnTruth-mean,2); //if linear estimates then lnTruth should be the "linear" truth
  group_diversity = gsl_stats_variance_m(dummy,1,N,mean);
  pid_t pid = getpid(); //output the process id to distinguish multiple workers executing the same parameter combination
  /*process id \t starting time \t collective error \t group diversity \t estimates of agent with lowest rank \t the rank itself*/
  *f_dynamics_out<<pid<<"\t"<<0<<"\t"<<collective_error<<"\t"<<group_diversity<<"\t"<<gsl_vector_get(estimates,rank_one)<<"\t"<<gsl_vector_get(ranks,rank_one)<<"\t"<<gsl_vector_get(estimates,gsl_vector_min_index(ranks))<<"\t"<<gsl_vector_get(ranks,gsl_vector_min_index(ranks))<<endl;
  /*output all estimates*/
  for (int j=0; j<N; j++)
    *f_dynamics_out<<gsl_vector_get(estimates,j)<<endl;

  /*loop over the whole time period*/
  for (int i = 0; i < t; i++) {
    /*simulate the dynamics of all agents*/
    compute_ranks();
    for (unsigned j=0; j<N; j++) {
      try {
	randnum = gsl_ran_gaussian(r,1.0);
      }
      catch (bad_alloc &ba){
	cerr<<"main():bad_alloc caught: "<<ba.what()<<endl;
	*f_dynamics_out<<"main():bad_alloc caught: "<<ba.what()<<endl;		
	return EXIT_FAILURE;
      }
      D = (Dmax*Dmin*exp(eta*gsl_vector_get(ranks,j))) / (Dmax+Dmin*(exp(eta*gsl_vector_get(ranks,j))-1));
	//D = W / (1.0 + (W-1.0)*exp(-1.0*eta*gsl_vector_get(ranks,j)));
      delta_i_t = D*sqrt(deltat)*randnum;
      gsl_vector_set(estimates,j,gsl_vector_get(estimates,j)+delta_i_t);
    }
    /*save the collective error and group diversity at the end of the current time step*/
    /* LINEAR ESTIMATES CONSIDERED AT THIS POINT */
    for (unsigned j=0; j<N; j++) 
      //dummy[j] = log(gsl_vector_get(estimates,j));   
      dummy[j] = gsl_vector_get(estimates,j);   
    
    mean=gsl_stats_mean(dummy,1,N);
    collective_error = pow(lnTruth-mean,2); //if linear estimates then lnTruth should be the "linear" truth
    group_diversity = gsl_stats_variance_m(dummy,1,N,mean);
    pid_t pid = getpid(); //output the process id to distinguish multiple workers executing the same parameter combination
    *f_dynamics_out<<pid<<"\t"<<i+1<<"\t"<<collective_error<<"\t"<<group_diversity<<"\t"<<gsl_vector_get(estimates,rank_one)<<"\t"<<gsl_vector_get(ranks,rank_one)<<"\t"<<gsl_vector_get(estimates,gsl_vector_min_index(ranks))<<"\t"<<gsl_vector_get(ranks,gsl_vector_min_index(ranks))<<endl;
    /*output all estimates*/
    for (int jj=0; jj<N; jj++)
      *f_dynamics_out<<gsl_vector_get(estimates,jj)<<endl;

    /*********************************************/
  }
  
  
  gsl_vector_free(estimates);
  gsl_vector_free(ranks);  
  gsl_rng_free(r);
  free(dummy);
  //  std::cout << collective_error << "\t";
  //  std::cout << group_diversity << "\t";
  //  std::cout << woc << std::endl;
  return EXIT_SUCCESS;
  
}


