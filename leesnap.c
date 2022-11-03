/* file leesnap.c
 * Routine for read 
 * GADGET snapshot 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "variables.h"
#include "allocate.h"
#include "leesnap.h"
#include "colores.h"


void read(void)
{
  type_int i, idim;

  cp.lbox = 1.0;
  cp.npart = 1000;
  
  if(!allocate_particles(&P, cp.npart))  exit(1);

  for(i=0; i<cp.npart; i++)
  {
    P[i].id = i;
    for(idim=0;idim<NDIM;idim++)
  	  P[i].pos[idim]   = (type_real)drand48()*cp.lbox;
  }
 
  return;
}
