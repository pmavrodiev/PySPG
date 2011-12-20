//:::~ 
//:::~ File automaticamente generado por ctt.py  NO EDITAR
//:::~ 

#ifndef __CTTBASE_H
#define __CTTBASE_H
#include <cstdlib>
#include <string>
#include <iostream>
#include <ctime>
#include <sys/time.h>
#include <fstream>
#include <dranxor.h>

//:::~ includes found in backends

#define VERSION_NUMBER ""
#define RELEASE_DATE ""


//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::~ Definiciones de las variables seteables desde el infile
//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++

namespace CTGlobal
{

  extern std::string filein /* = "input.rank2" */ ;

  extern std::string prog_name /*  */ ;

  extern int verbosityLevel /* = 0*/;
  extern bool quietRun /* = false*/;

//:::~  time step for the simulation
extern  double deltat /* =  0.01 */ ;
//:::~  number of agents
extern  int N /* =  100 */ ;
//:::~  the linear (or logarithm of) the Truth
extern  double lnTruth /* =  6 */ ;
//:::~  the maximum diffusion for the noise term
extern  double Dmax /* =  3.0 */ ;
//:::~  the minimum diffusion for the noise term
extern  double Dmin /* =  0.01 */ ;
//:::~  sensitivity of agents to their ranks
extern  double eta /* =  5.0 */ ;
//:::~  number of realizations per W,eta pair
extern  int B /* =  100 */ ;
//:::~  initial opinions
extern std::string filename /* =  "/home/pmavrodiev/run/rank.in" */ ;
//:::~  neighbourhood around the truth
extern  double Truth_bracket /* =  0.15 */ ;
//:::~  semilla de los numeros aleatorios
extern  long randomseed /* =  0 */ ;

// from backends

//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::~ Parse de las variables
//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++
void input_variables(std::istream &fin);

//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::~ Asking for output
//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++

void query_output(std::string type, std::string msg , std::string opening_closing = "<>",int indent=0);

//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::~ Small Help
//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++

void help_available();
//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::~ Option reading
//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++
void initialize_program(int &argc, char** & argv);
};

#endif

