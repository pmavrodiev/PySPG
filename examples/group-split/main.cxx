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

gsl_rng *r; //random number generator

/*
 The main interaction routine. Parameters:
 whole_group - a pointer to a pointer to all estimates
 regime - AGGR or FULL. Currently FULL is not implemented
 group1 - a pointer to a pointer to the estimates of group 1
 group2 - a pointer to a pointer to the estimates of group 2
 size_group1 - the size of group1 + 1. The +1 is needed for the copy algorithm
 size_group2 - the size of group2 + 1. The +1 is needed for the copy algorithm 
 
 */
int interact(double **whole_group, string *regime, double **group1, double **group2, unsigned int *size_group1, unsigned int *size_group2) {

  double *initial_estimates, *initial_group1, *initial_group2;
  double *estimates, *estimates_group1,*estimates_group2;//simpify notation
  bool gsplit_flag = 0; //0 - no splitting, 1 - splitting
  /* save the initial estimates for calculating individual conviction */
  if (whole_group == NULL) { //group splitting   
    estimates_group1=*group1;
    estimates_group2=*group2;
    initial_group1 = new double[*size_group1];
    initial_group2 = new double[*size_group2];
    copy(estimates_group1, estimates_group1 + *size_group1,initial_group1);
    copy(estimates_group2, estimates_group2 + *size_group2,initial_group2);
    gsplit_flag = 1;
  }
  else if ((group1 == NULL) && (group2 == NULL))  { //no group splitting   
    estimates = *whole_group;
    initial_estimates = new double[N+1];
    copy(estimates, estimates+N,initial_estimates);
    gsplit_flag = 0;    
  }
  else {
    cerr<<"interact(): Something is wrong with the default parameters. Check if interact() is invoked correctly"<<endl;
    *f_dynamics_out<<"interact(): Something is wrong with the default parameters. Check if interact() is invoked correctly"<<endl;
    return EXIT_FAILURE;
  } //end else
  /********************************************************************/
  if (*regime == "AGGR") {   
    double delta_i_t = 0, mean = 0;
    if (!gsplit_flag) { //no splitting    
      /* Euler: x_i(t+1)=x_i(t)+deltat*delta_i_t */
      /* loop over the whole time period */
      for (int i=0; i<T; i++) {
	/* simulate the dynamics of all agents */
	mean=gsl_stats_mean(estimates,1,N); //average of all estimates at current time step
	for (int j = 0; j < N; j++) {
	  double randnum = gsl_ran_gaussian(r,0.001);
	  delta_i_t = deltat*(alpha*(mean-estimates[j])+beta*(initial_estimates[j]-estimates[j]))+sqrt(2*Dn*deltat)*randnum;
	  estimates[j] += delta_i_t;
	}
      }
      /*save the final estimates*/
      for (int j=0; j<N; j++)
	*f_dynamics_out<<estimates[j]<<endl;
      delete [] initial_estimates;
    } //end if
    else if (gsplit_flag) { //splitting
      /*loop over the whole time period*/
      for (int i=0; i<T; i++) {
	/*simulate the dynamics of group 1. if the group size is 1, no need to do anything. remember group size is +1 at the end for the copy algorithm*/
	mean=gsl_stats_mean(estimates_group1,1,(*size_group1-1));
	for (int j=0; (j < (*size_group1-1) && (*size_group1) != 2) ; j++) {
	  double randnum = gsl_ran_gaussian(r,0.001);
	  delta_i_t = deltat*(alpha*(mean-estimates_group1[j])+beta*(initial_group1[j]-estimates_group1[j]))+sqrt(2*Dn*deltat)*randnum;
	  estimates_group1[j] += delta_i_t;
	}
	/*simulate the dynamics of group 2. if the group size is 1, no need to do anything. remember group size is +1 at the end for the copy algorithm*/
	mean=gsl_stats_mean(estimates_group2,1,(*size_group2-1));
	for (int j=0; (j < (*size_group2-1) && (*size_group2) != 2) ; j++) {
	  double randnum = gsl_ran_gaussian(r,0.001);
	  delta_i_t = deltat*(alpha*(mean-estimates_group2[j])+beta*(initial_group2[j]-estimates_group2[j]))+sqrt(2*Dn*deltat)*randnum;
	  estimates_group2[j] += delta_i_t;
	}
      } //end first interaction step
      /*save the final state of the two groups*/
      for (int i=0; i < (*size_group1-1); i++)
	*f_dynamics_out<<";\t"<<estimates_group1[i]<<endl;
      for (int i=0; i < (*size_group2-1); i++)
	*f_dynamics_out<<"-\t"<<estimates_group2[i]<<endl;
      /*let the two representative agents from both groups interact
        TODO: probably they would reach consensus much faster than T, so
        try with T/2 or T/4 */
      double rep_agent_group1 = gsl_stats_mean(estimates_group1,1,(*size_group1-1));
      double rep_agent_group1_initial = rep_agent_group1;
      double rep_agent_group2 = gsl_stats_mean(estimates_group2,1,(*size_group2-1));
      double rep_agent_group2_initial = rep_agent_group2;
      for (int i=0; i<T; i++) {
	mean = (rep_agent_group1 + rep_agent_group2) / 2;
	//first representative agent
	double randnum = gsl_ran_gaussian(r,0.001);
	delta_i_t = deltat*(alpha*(mean-rep_agent_group1)+beta*(rep_agent_group1_initial-rep_agent_group1))+sqrt(2*Dn*deltat)*randnum;
	rep_agent_group1 += delta_i_t;
	//second representative agent
	randnum = gsl_ran_gaussian(r,0.001);
	delta_i_t = deltat*(alpha*(mean-rep_agent_group2)+beta*(rep_agent_group2_initial-rep_agent_group2))+sqrt(2*Dn*deltat)*randnum;
	rep_agent_group2 += delta_i_t;
      } //end of second state interaction
      //save the final state of the two representative agents
      *f_dynamics_out<<";\t"<<rep_agent_group1<<endl;
      *f_dynamics_out<<"-\t"<<rep_agent_group2<<endl;
      delete [] initial_group1;
      delete [] initial_group2;
    } //end else if
    else {
      if (initial_group1 != NULL)
	delete [] initial_group1;
      if (initial_group2 != NULL)
	delete [] initial_group2;      
      if (initial_estimates != NULL)
	delete [] initial_estimates;
      cerr<< "Error in interact(): unexpected value for gsplit_flag"<<endl;
      *f_dynamics_out<< "Error in interact(): unexpected value for gsplit_flag"<<endl;
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  } //end AGGR regime
  else if (*regime == "FULL") {
    return EXIT_SUCCESS;
  } //end FULL regime
  else {
    cerr<< "Regime "<<*regime<<" not implemented"<<endl;
    *f_dynamics_out<< "Regime "<<*regime<<" not implemented"<<endl;
    return EXIT_FAILURE;
  }
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
//void run_aggr(int *size) {
  /* Euler: x_i(t+1)=x_i(t)+deltat*delta_i_t */
//  double delta_i_t=0,mean=0; 
  /*loop over the whole time period*/
//  for (long i=0; i<T; i++){
    /*simulate the dynamics of all agents*/
//     mean=gsl_stats_mean(estimates,1,*size);
//   for (unsigned j=0; j<*size; j++) {
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
    *f_dynamics_out<<"main(): Not enough memory to initialize random generator"<<endl;
    return EXIT_FAILURE;
  }
  
  ftime(&tp);
  //seed the rng with different value for each realization
  gsl_rng_set(r,tp.time+tp.millitm); 
  //generate N estimate
  try {
    estimates_group1 = new double[size_group1]; 
    estimates_group2 = new double[size_group2];
    for (int i=0; i < size_group1-1; i++) 
      estimates_group1[i] = gsl_ran_gaussian(r,sigma1);
    for (int i=0; i < size_group2-1; i++)
      estimates_group2[i] = gsl_ran_gaussian(r,sigma2);

  }
  catch (bad_alloc &ba){
    cerr<<"main():bad_alloc caught: "<<ba.what()<<endl;
    *f_dynamics_out<<"main():bad_alloc caught: "<<ba.what()<<endl;		
    return EXIT_FAILURE;
  }
  /******************************/
  /******************************/  

  /*PM: This if-then hell can be probably made more gracious*/
  if (control_case == "c") {
    try {
      estimates = new double[N+1];
      //merge the two groups
	copy(estimates_group1,estimates_group1 + size_group1, estimates);
	copy(estimates_group2,estimates_group2 + size_group2, estimates + size_group1 -1);
      //don't need them anymore
      delete [] estimates_group1;
      delete [] estimates_group2;
    }
    catch (bad_alloc &ba){
      cerr<<"main():bad_alloc caught for estimates: "<<ba.what()<<endl;
      *f_dynamics_out<<"main():bad_alloc caught for estimates: "<<ba.what()<<endl;
      return EXIT_FAILURE;
    }
    //save the initial estimates
    for (int i=0; i < N; i++) 
      *f_dynamics_out <<estimates[i]<<endl;
    
    /*start the interaction according to 'regime'*/
    int rv = interact(&estimates,&regime,NULL,NULL,NULL,NULL); 

    delete [] estimates;
  } //end if control_case == "c"
  else if (control_case == "t") {
    //save the initial estimates
    for (int i=0; i < size_group1-1; i++)
      *f_dynamics_out<<";\t"<<estimates_group1[i]<<endl;
    for (int i=0; i < size_group2-1; i++)
      *f_dynamics_out<<"-\t"<<estimates_group2[i]<<endl;

    int rv = interact(NULL,&regime,&estimates_group1,&estimates_group2,&size_group1,&size_group2);
    //don't need them anymore
    delete [] estimates_group1;
    delete [] estimates_group2;
    
  } //end else if control_case == "t"
  else {
    cerr<<"main(): Invalid value for 'Case' "<<control_case<<". Possible values are 'c' or 't'"<<endl;
    *f_dynamics_out<<"main(): Invalid value for 'Case' "<<control_case<<". Possible values are 'c' or 't'"<<endl;
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


