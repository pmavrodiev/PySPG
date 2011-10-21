/* main.cxx
 * 
 * Copyright (C) 2011 Pavlin Mavrodiev
 * Email: pmavrodiev@ethz.ch
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */


#include "base.h"

#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <algorithm>

//#include <dranxor.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/timeb.h>
#include <iomanip>

using namespace std;
using namespace CTGlobal;

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
//void run_aggr(int *size) {
  /* Euler: x_i(t+1)=x_i(t)+deltat*delta_i_t */
//  double delta_i_t=0,mean=0; 
  /*loop over the whole time period*/
//  for (long i=0; i<T; i++){
    /*simulate the dynamics of all agents*/
//   for (unsigned j=0; j<*size; j++) {
//     mean=gsl_stats_mean(estimates,1,*size);
//     double randnum = rt_rand_gaussian()*0.001;
//     delta_i_t = deltat*(alpha*(mean-estimates[j])+beta*(initial_estimates[j]-estimates[j]))+sqrt(2*Dn*deltat)*randnum;
//     estimates[j] += delta_i_t;
//   }
// }
//}

/***FULL INFO REGIME ***/
////void run_full(int *size) {
////  double *weigth_matrix = (double *) calloc((*size)*((*size)-1)/2,sizeof(double));
////  double *k1 = (double *) calloc((*size),sizeof(double));
////  double *k2 = (double *) calloc((*size),sizeof(double));
////  double *k3 = (double *) calloc((*size),sizeof(double));
////  double *k4 = (double *) calloc((*size),sizeof(double));
////  double xi,xj,normconst,K1,K2,K3,K4,noise;
////  int index;
////
////  if (weigth_matrix == NULL) 
////    perror("run_full(): failed to allocate memory to weigth matrix\n");
////  if (k1 == NULL) 
////    perror("run_full(): failed to allocate memory to weigth matrix\n");  
////  if (k2 == NULL) 
////    perror("run_full(): failed to allocate memory to weigth matrix\n");  
////  if (k3 == NULL) 
////    perror("run_full(): failed to allocate memory to weigth matrix\n");
////  if (k4 == NULL) 
////    perror("run_full(): failed to allocate memory to weigth matrix\n");
////  
////  for (long ii=0; ii<t; ii++) {
////    /*calculate the weigth matrix at level k1*/
////    weigth(&weigth_matrix,&estimates,0,&k1,&k2,&k3);
////    /*calculate the k1 vector for all agents*/
////    for (unsigned i=0; i<(*size); i++) {
////      xi=estimates[i];
////      /*calculate the normalization constant for agent i*/
////      normconst=0.0;
////      for (unsigned j=0; j<(*size); j++) 
////	normconst += 1.0 / (1+exp(fabs(estimates[j]-xi)/ alpha));
////      normconst = pow(normconst,-1);
////      /*fill the k1 vector with social influence only*/
////      K1=0;
////      for (unsigned j=0; j<(*size); j++) {
////	if (j != i) {
////	  xj=estimates[j];
////	  index = ((2*(*size)-i-1)*i) / 2 - i + j-1;
////	  if (i>j) index = ((2*(*size)-j-1)*j) / 2 - j + i-1;
////	  K1 += weigth_matrix[index]*(xj-xi);
////	}
////      }
////      K1 *= normconst;
////      /*add to the k1 vector individual conviction*/
////      K1 += beta*(initial_estimates[i]-xi);
////      K1 *= deltat;
////      /*add to the k1 vector noise*/
////      noise = sqrt(2*Dn*deltat)*rt_rand_gaussian()*0.001;
////      K1 += noise;
////      k1[i] = K1;
////    }
////    /*done with k1*/
////   /*calculate the weigth matrix at level k2*/
////    weigth(&weigth_matrix,&estimates,1,&k1,&k2,&k3);
////    /*calculate the k2 vector for all agents*/
////    for (unsigned i=0; i<(*size); i++) {
////      xi=estimates[i];
////      /*calculate the normalization constant for agent i*/
////      normconst=0.0;
////      for (unsigned j=0; j<(*size); j++) 
////	normconst += 1.0 / (1+exp(fabs(estimates[j]-xi+0.5*(k1[j]-k1[i]))/ alpha));
////      normconst = pow(normconst,-1);
////      /*fill the k2 vector with social influence only*/
////      for (unsigned j=0; j<(*size); j++) {
////	if (j != i) {
////	  xj=estimates[j];
////	  index = ((2*(*size)-i-1)*i) / 2 - i + j-1;
////	  if (i>j) index = ((2*(*size)-j-1)*j) / 2 - j + i-1;
////	  k2[i] += weigth_matrix[index]*(estimates[j]-estimates[i]+0.5*(k1[j]-k1[i]));
////	}
////      }
////      k2[i] *= normconst;
////      /*add to the k2 vector individual conviction*/
////      k2[i] += beta*(initial_estimates[i]-estimates[i]-0.5*k1[i]);
////      k2[i] *= deltat;
////      /*add to the k2 vector noise*/
////      k2[i] += sqrt(2*Dn*deltat)*rt_rand_gaussian()*0.001;
////    }
////    /*done with k2*/
////    
////   /*calculate the weigth matrix at level k3*/
////    weigth(&weigth_matrix,&estimates,2,&k1,&k2,&k3);
////    /*calculate the k3 vector for all agents*/
////    for (unsigned i=0; i<(*size); i++) {
////      xi=estimates[i];
////      /*calculate the normalization constant for agent i*/
////      normconst=0.0;
////      for (unsigned j=0; j<(*size); j++) 
////	normconst += 1.0 / (1+exp(fabs(estimates[j]-xi+0.5*(k2[j]-k2[i])) / alpha ));
////      normconst = pow(normconst,-1);
////      /*fill the k3 vector with social influence only*/
////      for (unsigned j=0; j<(*size); j++) {
////	if (j != i) {
////	  xj=estimates[j];
////	  index = ((2*(*size)-i-1)*i) / 2 - i + j-1;
////	  if (i>j) index = ((2*(*size)-j-1)*j) / 2 - j + i-1;
////	  k3[i] += weigth_matrix[index]*(estimates[j]-estimates[i]+0.5*(k2[j]-k2[i]));
////	}
////      }
////      k3[i] *= normconst;
////      /*add to the k3 vector individual conviction*/
////      k3[i] += beta*(initial_estimates[i]-estimates[i]-0.5*k2[i]);
////      k3[i] *= deltat;
////      /*add to the k3 vector noise*/
////      k3[i] += sqrt(2*Dn*deltat)*rt_rand_gaussian()*0.001;
////    }
////    /*done with k3*/
////   /*calculate the weigth matrix at level k4*/
////    weigth(&weigth_matrix,&estimates,3,&k1,&k2,&k3);
////    /*calculate the k4 vector for all agents*/
////    for (unsigned i=0; i<(*size); i++) {
////      xi=estimates[i];
////      /*calculate the normalization constant for agent i*/
////      normconst=0.0;
////      for (unsigned j=0; j<(*size); j++) 
////	normconst += 1.0 / (1+exp(fabs(estimates[j]-xi+k3[j]-k3[i])/alpha));
////      normconst = pow(normconst,-1);
////      /*fill the k4 vector with social influence only*/
////      for (unsigned j=0; j<(*size); j++) {
////	if (j != i) {
////	  xj=estimates[j];
////	  index = ((2*(*size)-i-1)*i) / 2 - i + j-1;
////	  if (i>j) index = ((2*(*size)-j-1)*j) / 2 - j + i-1;
////	  k4[i] += weigth_matrix[index]*(estimates[j]-estimates[i]+k3[j]-k3[i]);
////	}
////      }
////      k4[i] *= normconst;
////      /*add to the k2 vector individual conviction*/
////      k4[i] += beta*(initial_estimates[i]-estimates[i]-k3[i]);
////      k4[i] *= deltat;
////      /*add to the k2 vector noise*/
////      k4[i] += sqrt(2*Dn*deltat)*rt_rand_gaussian()*0.001;
////    }
////    /*done with k4*/
////
////    /* update estimates*/
////    for (unsigned i=0; i<(*size); i++) {
////      estimates[i] += (1.0 / 6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
////    }
////
////  } //end ii
////  free(weigth_matrix);
////  free(k1);
////  free(k2);
////  free(k3);
////  free(k4);
////
////}
////
int main(int argc, char *argv[])
{
 
  initialize_program( argc,  argv );
  double *estimates, *estimates_group1, *estimates_group2;
  gsl_rng *r;
  time_t seconds;
  struct timeb tp;
  //empty element at the end for the 'copy' algo below + 0.5 for correct floor rounding to unsigned int
  //for some reason n2 will not equal its expected value from parameters_gsplit.dat. e.g. n2=0.01 results in some 
  //trailing digits at the very high bits, i.e. n2=0.0100000000007893, which casted to int without care gives
  //n2=0.02 almost always.
  unsigned int size_group1 = static_cast<unsigned int> ((1.0-n2)*N+1.0+0.5); 
  unsigned int size_group2 =  static_cast<unsigned int> (n2*N+1.0+0.5);
 
  /*******************************/
  /*****INIT THE POPULATION*******/
  r = gsl_rng_alloc(gsl_rng_taus);
  if (r == NULL) {
    cerr<<"main(): Not enough memory to initialize random generator"<<endl;
    return EXIT_FAILURE;
  }
  
  ftime(&tp);
  //seconds=time(NULL);
  *f_dynamics_out <<tp.time+tp.millitm<<endl;
  //seed the rng with different value for each realization
  gsl_rng_set(r,tp.time+tp.millitm);
 
 
  //generate N estimate
  try {
    estimates_group1 = new double[size_group1]; 
    estimates_group2 = new double[size_group2];
  }
  catch (bad_alloc &ba){
    cerr<<"main():bad_alloc caught: "<<ba.what()<<endl;
    return EXIT_FAILURE;
  }
  

  for (int i=0; i < size_group1-1; i++) 
    estimates_group1[i] = gsl_ran_gaussian(r,sigma1);
  for (int i=0; i < size_group2-1; i++)
    estimates_group2[i] = gsl_ran_gaussian(r,sigma2);
  /******************************/
  /******************************/
  

  /*PM: This if-then hell can be probably made more gracious*/
  if (control_case == "c") {
    try {
      estimates = new double[N+1];
    }
    catch (bad_alloc &ba){
      cerr<<"main():bad_alloc caught for estimates: "<<ba.what()<<endl;
      return EXIT_FAILURE;
    }
    //merge the two groups
    copy(estimates_group1,estimates_group1 + size_group1, estimates);
    copy(estimates_group2,estimates_group2 + size_group2, estimates + size_group1 -1);
    //don't need them anymore
    delete [] estimates_group1;
    delete [] estimates_group2;
    //save the initial estimates
    for (int i=0; i < N; i++) 
      *f_dynamics_out <<estimates[i]<<endl;
    

  } //end if



  else if (control_case == "t") {
  } //end else if
  else {
    cerr<<"main(): Invalid value for 'Case' "<<control_case<<". Possible values are 'c' or 't'"<<endl;
    return EXIT_FAILURE;
  }


  /*  if (regime == "AGGR")
    run_aggr(&N);
  else if (regime == "FULL")
    run_full(&N);
  else 
    {
      printf("main(): Unrecognized regime %s\n",regime.c_str());
      perror("");
    }
  */

  gsl_rng_free(r);
  return EXIT_SUCCESS;
  
}


