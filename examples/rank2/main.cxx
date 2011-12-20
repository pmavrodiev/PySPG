/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Tue Dec 20 2011
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
#include <string.h>
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
    grace_exit();
  }

  ranks = gsl_vector_calloc(N); //init and set all elements to 0
  if (ranks == NULL) {
    GSL_ERROR("read_initial_estimates(): failed to allocate memory to ranks",GSL_ETOL);
    grace_exit();
  }
  
  inFile = fopen(filename.c_str(),"r");
  if (inFile == NULL) {
    perror(getcwd(NULL, 100));
    grace_exit();
  }
  
  if (gsl_vector_fscanf (inFile,estimates)) {//TODO: handle mismatch between N and the length of the input file
    GSL_ERROR("read_initial_estimates(): failed to read input file", GSL_ETOL);
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
  double median_fpt=-1.0,mean_fpt=-1.0,mode_fpt=-1.0;
  struct timeb tp;

  r = gsl_rng_alloc(gsl_rng_taus);
  if (r == NULL) {
    cerr<<"main(): Not enough memory to initialize random generator"<<endl;
    return EXIT_FAILURE;
  }
  ftime(&tp);
  //seed the rng with different value for each realization
  gsl_rng_set(r,tp.time+tp.millitm); 
  
  double *dummy = (double *) malloc(N*sizeof(double));
  if (dummy == NULL) {
    grace_exit();
  } 
  

  bool mean_done=false,median_done=false,mode_done=false;
  double time_step=0;
  /*loop over the whole time period*/
  while (!mean_done || !median_done || !mode_done) {
    time_step++;
    /*simulate the dynamics of all agents for one time step*/
    compute_ranks();
    for (unsigned j=0; j<N; j++) {
      try {
	randnum = gsl_ran_gaussian(r,1.0);
      }
      catch (bad_alloc &ba){
	cerr<<"main():bad_alloc caught: "<<ba.what()<<endl;
	return EXIT_FAILURE;
      }
      D = (Dmax*Dmin*exp(eta*gsl_vector_get(ranks,j))) / (Dmax+Dmin*(exp(eta*gsl_vector_get(ranks,j))-1));
	//D = W / (1.0 + (W-1.0)*exp(-1.0*eta*gsl_vector_get(ranks,j)));
      delta_i_t = D*sqrt(deltat)*randnum;
      gsl_vector_set(estimates,j,gsl_vector_get(estimates,j)+delta_i_t);
    }
    /*calculate the FPT of the mean*/
    if (!mean_done) {
      memcpy(dummy,gsl_vector_ptr(estimates,0),N*sizeof(double));
      //dummy = gsl_vector_ptr(estimates,0);
      mean_fpt=fabs(gsl_stats_mean(dummy,1,N)-lnTruth);
      if (mean_fpt <= fabs(Truth_bracket*lnTruth)) {
	mean_done=true;
	mean_fpt=time_step;
      }
    }
    /*calculate the FPT of the median*/
    if (!median_done) {
      gsl_sort_vector(estimates);
      //dummy = gsl_vector_ptr(estimates,0);
      memcpy(dummy,gsl_vector_ptr(estimates,0),N*sizeof(double));
      median_fpt = fabs(gsl_stats_median_from_sorted_data(dummy,1,N) - lnTruth);
      if (median_fpt <= fabs(Truth_bracket*lnTruth)) {
	median_done=true;
	median_fpt=time_step;
      }
    }
    /*calculate the FPT of the mode*/
    if (!mode_done) {
      mode_done=true;
      /*TODO*/
    }
    /*********************************************/
  } //end while (!mean_done || !median_done || !mode_done) {
  
  
  gsl_vector_free(estimates);
  gsl_vector_free(ranks);  
  gsl_rng_free(r);
  free(dummy);
  std::cout << mean_fpt << "\t";
  std::cout << median_fpt << "\t";
  std::cout << mode_fpt << std::endl;
  return EXIT_SUCCESS;
  
}


