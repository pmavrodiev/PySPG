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
#include <sstream>

#define VERSION_NUMBER ""
#define RELEASE_DATE ""


//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::~ Definiciones de las variables seteables desde el infile
//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++
//:::~ +++++++++++++++++++++++++++++++++++++++++++++++++++++++

namespace CTGlobal
{

  extern std::string filein /* = "input.rank" */ ;

  extern std::string prog_name /*  */ ;

  extern int verbosityLevel /* = 0*/;
  extern bool quietRun /* = false*/;

//:::~  sets whether to store the dynamics
extern bool  store_dynamics /* =  false */ ;
//:::~  file name to store the dynamics
extern std::string store_dynamics_filename /* =  "dynamics.out" */ ;
//:::~  time step for the simulation
extern  double deltat /* =  0.01 */ ;
//:::~  number of time steps in the simulation
extern  long t /* =  300 */ ;
//:::~  number of agents
extern  int N /* =  100 */ ;
//:::~  the logarithm of the Truth
extern  double lnTruth /* =  6.62 */ ;
//:::~  maximum diffusion for the noise term of the agent ranked last
extern  double W /* =  5.0 */ ;
//:::~  sensitivity of agents to their ranks
extern  double eta /* =  5.0 */ ;
//:::~  number of realizations per W,eta pair
extern  int R /* =  100 */ ;
//:::~  initial opinions
extern std::string filename /* =  "/home/pmavrodiev/run/rank.in" */ ;
//:::~  semilla de los numeros aleatorios
extern  long randomseed /* =  0 */ ;

// from backends
extern  std::fstream *f_dynamics_out  ;

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

