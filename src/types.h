#pragma once

//! NOTE: 13 JUL 23
//! Updated MPI compilers (e.g., OpenMPI) don'e use C++ MPI bindings
//! anymore, which has broken all of the "MPI:: " commands.
//! We need to redefine the MPI namespace and a few typedefs here.
//! Much of the updated code has been taken from Niklas Donocik (toxciq-math)
//! which can be found at https://github.com/siloutil/aikef

//! LL updated all MPI commands on 13 JUL 23.
//! Relevant changed files are:
//  CBlk_FieldCom.cpp
//  CBlk_MPI_PartCom.cpp
//  CBlk_Part-Com.cpp
//  CHybrid_BlockCom.cpp
//  CHybrid_FileCom.cpp
//  CHybrid_Init_IonNeutral.cpp
//  CHybrid_IonProfiles.cpp
//  CHybrid_MPI_Init.cpp
//  CHybrid_MPI_Massload.cpp
//  defines.h
//  utils_GN_MPI.cpp

#include <mpi.h>

namespace MPI {
    typedef MPI_Comm Intracomm;
    typedef MPI_Request Request;
    typedef MPI_Datatype Datatype;
}


//! As "int" is 4 BYTE on all todays architechturs (which
//! results in a Range of [-2e9;+2e.9]), use int and avoid
//! using long.

//! All variable types are replaced with own types for 
//! the following reasons:
//! 1) If an error accurs, by setting all types to double precision
//!    an memory overflow can be easyly excluded
//! 2) Different systems use different precison for the same type 
//! 	(eg. int 2byte or 4 byte).
//! 3) Plotting File Stream is set to single Precesion by default.
//!    This can be easily replaced if necessary.

//! integer datatypes
#define MPI_INT64  MPI_LONG_LONG
typedef long long  INT64;

#define MPI_INT32  MPI_INT
typedef int  INT32;

//! floating point variables:
#define MPI_D_REAL  MPI_DOUBLE
typedef double  D_REAL;

#define MPI_F_REAL  MPI_FLOAT
typedef float  F_REAL;

#define MPI_PARTICLE_REAL  MPI_DOUBLE
typedef double  PARTICLE_REAL;

#define MPI_WEIGHT_REAL  MPI_DOUBLE
typedef double  WEIGHT_REAL;

//! always use float for FILE_REAL, else
//! silo functions will return errors
#define MPI_FILE_REAL  MPI_FLOAT
typedef float  FILE_REAL;