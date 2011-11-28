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

  extern std::string filein /* = "input.random_walk" */ ;

  extern std::string prog_name /*  */ ;

  extern int verbosityLevel /* = 0*/;
  extern bool quietRun /* = false*/;

//:::~  model type
extern std::string model /* =  "UNBIASED" */ ;
//:::~  diffusion coefficient
extern  double D /* =  1.0 */ ;
//:::~  drift coefficient
extern  double drift /* =  0.0 */ ;
//:::~  filtering time equalts this times ts
extern  long filter_timesteps /* =  0 */ ;
//:::~  fnumber of timesteps in the simulation
extern  long simulation_timesteps /* =  100 */ ;
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

