#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "variables.h"
#include "grid.h"

extern void grid_init(void)
{
  unsigned long long nalloc = pow(grid.ngrid, NDIM) + 1;
  fprintf(stdout,"Allocating %.5f gb\n",
          (double)(nalloc*sizeof(type_int))/1024.0/1024.0/1024.0);
  fflush(stdout);

  grid.icell	= (type_int *) malloc(nalloc*sizeof(type_int));
  assert(grid.icell != NULL);
  for(unsigned long i = 0; i<=nalloc; i++)
    grid.icell[i] = cp.npart;

  return;
}

/*  This function brings the particles 
 *  into the same order as
 *  the sorted auxiliary. 
 */
static void reorder_particles(type_int *tmp_Id)
{ 
  type_int i;

  for(i=0;i<cp.npart;i++)
  {
    if(tmp_Id[i] != i)
    {

      struct particle_data P_source = P[i];

      type_int  idsource = tmp_Id[i];
      type_int  dest     = tmp_Id[i];

      while(1)
      {
        struct particle_data P_save = P[dest];

	      type_int idsave = tmp_Id[dest];

        P[dest] = P_source;

        tmp_Id[dest] = idsource;

        if(dest == i)  break;

        P_source = P_save;

        idsource = idsave;
  	    dest = idsource;
      } // close while
    }  // close if

  } // close for

  return;
}

extern void grid_build(void)
{
  unsigned long long nalloc = pow(grid.ngrid, NDIM);
  type_int i, j, idim;
  type_int *tmp_Id; // Auxiliary array
  long long ix, ibox;
  double fac = (double)grid.ngrid/(double)cp.lbox ;
  fprintf(stdout,"Building Grid..... Ngrid = %lu\n",grid.ngrid);
  tmp_Id = (type_int *) malloc(grid.nobj*sizeof(type_int));
  assert(tmp_Id != NULL);

  for(i = 0; i < grid.nobj; i++)
  {
  
    ibox = 0;
    for(idim=0; idim<NDIM; idim++)
    {
      ibox *= grid.ngrid;
      ix = (long long)((double)P[i].pos[idim]*fac);
      ix = (ix >= (long long)grid.ngrid) ? (long long)grid.ngrid-1 : ( (ix<0) ? 0 : ix );
      ibox += ix;
    }

    assert(ibox>=0 && ibox<nalloc);

    tmp_Id[i] = grid.icell[ibox];
    grid.icell[ibox] = i;
  }

  // Sorted tmp_Id

  j = 0;
  for(ibox = 0; ibox < nalloc; ibox++)
  {      
    i = grid.icell[ibox];
    grid.icell[ibox] = j;

    while(i != cp.npart)
    {
      unsigned long tmp = tmp_Id[i];
      tmp_Id[i] = j;
      j++;	
      i = tmp;
    }
  }

  //Ghost cell
  grid.icell[nalloc] = cp.npart;

  reorder_particles(tmp_Id);

  free(tmp_Id);
}


extern void grid_free(void)
{
  if(grid.icell!=NULL) free(grid.icell);

  return;
}
