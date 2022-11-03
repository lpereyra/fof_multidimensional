/*  file allocate.c
 *  Routines for allocating particle
*/

#include <stdlib.h>
#include <stdio.h>
#include "variables.h"
#include "allocate.h"

/*
 *  Allocates structure with the position, 
 *  velocity, ids of the particles and
 *  an auxiliary array for FoF groups.
 */ 
extern int allocate_particles(struct particle_data **Q, const type_int size)
{

  *Q = NULL;

  *Q = (struct particle_data *) malloc(size*sizeof(struct particle_data));

  if(!*Q) 
  {
    fprintf(stderr, "cannot allocate particles\n" );
    return(0);
  }    

  return ( 1 );
}

/*
 *  Deallocates structure with the position, 
 *  velocity, ids of the particles and
 *  an auxiliary array for FoF groups.
 */ 
extern void free_particles(struct particle_data **Q)
{
  if(*Q) free(*Q);
}

/*
 *  Reallocates structure with the position, 
 *  velocity, ids of the particles and
 *  an auxiliary array for FoF groups.
 */ 
extern int reallocate_particles(struct particle_data **Q, const type_int size)
{

  *Q = (struct particle_data *) realloc(*Q,size*sizeof(struct particle_data));

  if(!*Q) 
  {
    fprintf(stderr, "cannot reallocate particles\n" );
    return(0);
  }    

  return ( 1 );
}


