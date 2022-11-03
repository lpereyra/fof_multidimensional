/* file variables.h
 * brief declares global variables.
 *
 * This file declares all global variables. 
 * Further variables should be added here, and declared as 'extern'. 
 * The actual existence of these variables is provided by the file 'variables.c'. 
 * To produce 'variables.c' from 'variables.h', do the following:
 *
 *    - Erase all #define's, typedef's, and enum's
 *    - add #include "variables.h", delete the #ifndef VARIABLES_H conditional
 *    - delete all keywords 'extern'
 *    - delete all struct definitions enclosed in {...}, e.g.
 *      "extern struct cosmoparam {....} cp;"
 *      becomes "struct cosmoparam cp;"
 */

#ifndef VARIABLES_H
#define VARIABLES_H

#ifndef NPARTMIN
  #define NPARTMIN 2
#endif

#ifndef FOF
  #define FOF 0.2
#endif

#define NDIM 7

/* If defined, the type variable 
 * is set to "double", otherwise to FLOAT 
 */
#ifdef PRECDOUBLE
typedef double type_real;
#else
typedef float type_real;
#endif

/* Precision del codigo (enteros) */
#ifdef LONGIDS
typedef unsigned long long type_int;
#else
typedef unsigned int type_int;
#endif

extern struct cosmoparam
{
  double   lbox  ;  /* Boxsize [Kpc / h]                        */
  type_int npart ;  /* Particle number                          */
} cp;

/* Input simulation files */
extern struct SnapST
{
  int nfiles;
  char root[200], name[200];
} snap;

extern struct gridst
{
  unsigned long ngrid;
  unsigned long nobj;
  type_int *icell;
} grid;

/* This structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern struct particle_data 
{
  type_real pos[NDIM];
  type_int  id;
} *P;

extern char message[200];

#endif
