#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <math.h>

int main(int argc, char **argv) {
  FILE *inFile=NULL;
  double estimate=0;
  inFile = fopen("wisdom.in","r");
  while (fscanf(inFile,"%lf",&estimate) == 1) {
    printf("%.18lg\n",estimate);
  }
  
  return 0;
}
