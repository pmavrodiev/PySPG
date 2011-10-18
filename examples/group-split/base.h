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

  extern std::string filein /* = "input.wisdom" */ ;

  extern std::string prog_name /*  */ ;

  extern int verbosityLevel /* = 0*/;
  extern bool quietRun /* = false*/;

//:::~  information regime
extern std::string regime /* =  "AGGR" */ ;
//:::~  strength of social influence
extern  double alpha /* =  1.0 */ ;
//:::~  strength of individual conviction
extern  double beta /* =  1.0 */ ;
//:::~  time step for the simulation
extern  double deltat /* =  0.01 */ ;
//:::~  number of time steps in the simulation
extern  long t /* =  300 */ ;
//:::~  initial estimates (fig 4 top-left) are stored here
extern std::string filename /* =  "wisdom.in" */ ;
//:::~  number of agents
extern  int N /* =  100 */ ;
//:::~  noise diffusion coefficient
extern  double Dn /* =  0.5 */ ;
//:::~  the logarithm of the Truth
extern  double lnTruth /* =  -2.0 */ ;
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

