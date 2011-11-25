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
#include <algorithm>

#include <dranxor.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>

using namespace std;
using namespace CTGlobal;

gsl_vector *estimates=NULL; //store estimates at time t
gsl_vector *initial_estimates=NULL; //store initial estimates

unsigned int *ranks=NULL; //vector storing the ranks of all agents
FILE *inFile=NULL; //read in the initial estimates

/*gracefully exit by closing all open file handles, freeing memory, etc.*/
void grace_exit() {
  if (estimates != NULL) gsl_vector_free(estimates);
  if (initial_estimates != NULL) gsl_vector_free(initial_estimates);
  if (ranks != NULL) free(ranks);
  if (inFile != NULL) fclose(inFile);
  exit(1);
}

/*Compute the ranks of the population in the current time period.
  Taken from http://www.gnu.org/s/gsl/manual/html_node/Computing-the-rank.html*/
void compute_ranks (gsl_vector *v) {
  size_t i;
  size_t n = v->size;
  gsl_permutation * perm = gsl_permutation_alloc(n);
  gsl_permutation * rank = gsl_permutation_alloc(n);
  
  gsl_sort_vector_index (perm, v);
  gsl_permutation_inverse (rank, perm);
  
  for (i = 0; i < n; i++) {
      double vi = gsl_vector_get(v, i);
      printf ("element = %d, value = %g, rank = %d\n",
	      i, vi, rank->data[i]);
    }
  
  gsl_permutation_free (perm);
  gsl_permutation_free (rank);
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

  initial_estimates = gsl_vector_calloc(N); //init and set all elements to 0
  if (initial_estimates == NULL) {
    GSL_ERROR("read_initial_estimates(): failed to allocate memory to initial_estimates",GSL_ETOL);
    grace_exit();
  }

  inFile = fopen(filename.c_str(),"r");
  if (inFile == NULL) {
    perror("read_initial_estimates():");
    grace_exit();
  }

  if (gsl_vector_fscanf (inFile,estimates)) {//TODO: handle mismatch between N and the length of the input file
    GSL_ERROR("read_initial_estimates(): failed to read input file", GSL_ETOL);
    grace_exit();
  }
  /*
  while ((fscanf(inFile,"%lf",&estimate) == 1) && (counter<N)) {
    gsl_vector_set(estimates,counter) = estimates;
    gsl_vector_set(initial_estimates,counter++) = estimates;
  }  

  if (counter<N) { //fewers estimates in the input file than we expected from base.ct
    printf("read_initial_estimates(): Expected %d estimates from base.ct but got only %d from %s\n",N,counter,filename.c_str());
    grace_exit();
  }
  */
  fclose(inFile);

  return 0;
}


int main(int argc, char *argv[]) {
 
  initialize_program( argc,  argv );
  read_initial_estimates();
  double collective_error = 0;
  double group_diversity = 0;
  double woc=0;
  /* Euler: x_i(t+1)=x_i(t)+deltat*delta_i_t */
  double delta_i_t=0,mean=0;


  for (int i = 0; i < N; i++)
         {
	   std::cout <<  gsl_vector_get (estimates, i)<<std::endl;
         }
  //std::cout << collective_error << "\t";
  //std::cout << group_diversity << "\t";
  //std::cout << woc << std::endl;

  free(estimates);
  free(initial_estimates);  
  return EXIT_SUCCESS;
  
}


