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

  extern std::string filein /* = "input.gsplit" */ ;

  extern std::string prog_name /*  */ ;

  extern int verbosityLevel /* = 0*/;
  extern bool quietRun /* = false*/;

//:::~  sets whether to store the dynamics
extern bool  store_dynamics /* =  false */ ;
//:::~  file name to store the dynamics
extern std::string store_dynamics_filename /* =  "dynamics.out" */ ;
//:::~  information regime
extern std::string regime /* =  "AGGR" */ ;
//:::~  (c)ontrol case (no splitting) or (t)est case (splitting)
extern std::string control_case /* =  "c" */ ;
//:::~  std of group one
extern  double sigma1 /* =  1.0 */ ;
//:::~  std of group two
extern  double sigma2 /* =  1.0 */ ;
//:::~  strength of social influence
extern  double alpha /* =  1.0 */ ;
//:::~  strength of individual conviction
extern  double beta /* =  1.0 */ ;
//:::~  time step for the simulation
extern  double deltat /* =  0.01 */ ;
//:::~  number of time steps for interaction
extern  int T /* =  300 */ ;
//:::~  number of agents
extern  int N /* =  100 */ ;
//:::~  number of realizations
extern  int M /* =  100 */ ;
//:::~  relative size of the second group
extern  double n2 /* =  0.01 */ ;
//:::~  noise diffusion coefficient
extern  double Dn /* =  0.5 */ ;
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

