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

//double opinions[100] = {
// 0.0304541776502582,
// 0.0136051394294611,
// 0.0431026788469159, 
// 0.0605170431191805, 
// 0.011067750504134 ,
// 0.0123585058069386,
// 0.0305099463621518, 
// 0.179329529796854 ,
// 0.0624516059507075,
// 0.0662676616284229,
// 0.0445737518843613, 
// 0.235273227961589 ,
// 0.0710027174747392,
// 0.0389383684356751,
// 0.30726454570021  ,
// 0.0405609127605168,
// 0.0485079932751849,
// 0.0231962544002043,
// 0.0287594489617644, 
// 0.0578672118132589,
// 0.0529006379372462,
// 0.0195822215113274,
// 0.0293770396670783, 
// 0.023973263512744 ,
// 0.196233825684164 ,
// 0.0425385136072795,
// 0.021894506916699 ,
// 0.0956668339153118,
// 0.0562371521953241,
// 0.178281837257854 ,
// 0.722914354338746 ,
// 0.0883513854648337,
// 0.0433088352447022,
// 0.101146502869907 ,
// 0.0290858883672253, 
// 0.063047393861378 ,
// 0.0273716476006422,
// 0.126017691629181 ,
// 0.070285464978148 ,
// 0.0337161632913864,
// 0.0645509997419784,
// 0.0839399295342982,
// 0.17739755094298  ,
// 0.0184636065766418,
// 0.0671706471526621,
// 0.189824475734004 ,
// 0.0453995879831098, 
// 0.102495605076897 ,
// 0.0245657534029934, 
// 0.127452014970281 ,
// 0.0220622423118579, 
// 0.291387561598921 ,
// 0.0313344527111527,
// 0.0149628851051261,
// 0.0288653540859111, 
// 0.125381098283336 ,
// 0.0373500505674803, 
// 0.0178544061050245,
// 0.0400423507876844, 
// 0.101597615754887 ,
// 0.0376217833863892, 
// 0.100050688294814 ,
// 0.10080589070722  ,
// 0.0773597761027292,
// 0.187270932584715 ,
// 0.0193999835567671,
// 0.237642467610782 ,
// 0.046352992699674 ,
// 0.0832581837984352,
// 0.0300979035345081,
// 0.0799021102129099, 
// 0.0714401599496476, 
// 0.119930945307331 ,
// 0.0212521500799866,
// 0.104738191502093 ,
// 0.0263640462519334,
// 0.0364401719754576,
// 0.123691373513805 ,
// 0.0428754301860053, 
// 0.087501122959982 ,
// 0.121910296112429 ,
// 0.0407545780593138, 
// 0.0353168397324526, 
// 0.0301868349490974,
// 0.0616621016369282, 
// 0.151612794798738 ,
// 0.236236472476187 ,
// 0.03073581317959  ,
// 0.0723721958602008, 
// 0.0167675631313543, 
// 0.0178052424514191, 
// 0.0251066017538328,
// 0.0156064193125263, 
// 0.0166002961228833, 
// 0.341721199821824 ,
// 0.0721522746676428, 
// 0.0280623148701852,
// 0.0363378751421454,
// 0.077544698307166 ,
// 0.0488149539027285
// };

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
/*TODO-FINISH THIS*/
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
  gsl_vector_add_constant(relative_estimates,-(exp(lnTruth)+50.0));
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
  
  /*shift the distribution to the right to avoid agents going to negative opinions*/
  gsl_vector_add_constant(estimates,50.0);
  
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
      D = W / (1.0 + (W-1.0)*exp(-1.0*eta*gsl_vector_get(ranks,j)));
      delta_i_t = D*sqrt(deltat)*randnum;
      gsl_vector_set(estimates,j,gsl_vector_get(estimates,j)+delta_i_t);
    }
    /*save the collective error and group diversity at the end of the current time step*/
    for (unsigned j=0; j<N; j++) 
      dummy[j] = log(gsl_vector_get(estimates,j));   
    
    mean=gsl_stats_mean(dummy,1,N);
  
    collective_error = pow(lnTruth-mean,2);
    group_diversity = gsl_stats_variance_m(dummy,1,N,mean);
    *f_dynamics_out<<i<<"\t"<<collective_error<<"\t"<<group_diversity<<endl;
    /*********************************************/
  }
  
  
  gsl_vector_free(estimates);
  gsl_vector_free(ranks);  
  gsl_rng_free(r);
  free(dummy);
  std::cout << collective_error << "\t";
  std::cout << group_diversity << "\t";
  std::cout << woc << std::endl;
  return EXIT_SUCCESS;
  
}


