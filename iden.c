#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include "variables.h"
#include "allocate.h"
#include "leesnap.h"
#include "grid.h"
#include "iden.h"
#include "io.h"

#define DIV_CEIL(x,y) (x+y-1)/y

static struct iden_st iden;
static struct temporary Temp;
static type_int *gr;
#ifdef LOCK
  static omp_lock_t *lock;
#endif
 
static inline type_int Raiz(type_int i, type_int * restrict ar)
{
  if(i != ar[i])
    ar[i] = Raiz(ar[i],ar);
 
  return ar[i];
}

static inline void Unir(type_int u, type_int v, type_int * restrict ar)
{
 
  type_int z;

  while(ar[u] != ar[v])
  { 
      if(ar[u] < ar[v])
      {
#ifdef LOCK
          if(u == ar[u])
          {
            omp_set_lock(&(lock[u]));
            z = 0;

            if(u == ar[u])
            {
              ar[u] = ar[v];  
              z = 1;
            } 

            omp_unset_lock(&(lock[u]));
            if(z==1) break;             

          }
#else
          if(u == ar[u])
          {
            ar[u] = ar[v];  
            break;             
          }
#endif
          
          z = ar[u];   
          ar[u] = ar[v];
          u = z;

      }else{
#ifdef LOCK
          if(v == ar[v])
          {
            omp_set_lock(&(lock[v]));
            z = 0;

            if(v == ar[v])
            {
              ar[v] = ar[u];   
              z = 1;
            }

            omp_unset_lock(&(lock[v]));
            if(z == 1) break;            
          }
#else
          if(v == ar[v])
          {
            ar[v] = ar[u];   
            break;             
          }
#endif

          z = ar[v];   
          ar[v] = ar[u];   
          v = z;

      }
  }

}

static void test(const type_int ic, long long *counters)
{
  type_int i, idim;
  type_real xx, dd;
  long long ix, ibox;

  ibox = 0;
  for(ix=0; ix<NDIM; ix++)
  {
    ibox *= grid.ngrid;
    ibox += counters[ix];
  }

  for(i=grid.icell[ibox];i<grid.icell[ibox+1];i++)
  {
    if(ic<i)
    {
      dd = 0.0;
      for(idim=0;idim<NDIM;idim++)
      {
        xx = P[i].pos[idim] - P[ic].pos[idim];
        dd += xx*xx;
      }

      if(dd < iden.r0)
      {
        Unir(ic,i,gr);
      }
    } // cierra el if
  }

  return;
}

static void loopn_recurse(const type_int ic, long long *min, long long *max, long long *counters, type_int n, void (*func)(const type_int, long long*)) 
{
  long long i;

  for(i = min[n]; i <= max[n]; ++i) 
  {
    counters[n] = i;

    if (NDIM - n > 1) {
      loopn_recurse(ic, min, max, counters, n + 1, func);
    } else {
      func(ic, counters);
    }
  }
}

static void loopn(const type_int ic, long long *min, long long *max, void (*func)(const type_int, long long*)) 
{
  loopn_recurse(ic, min, max, (long long[NDIM]){}, 0, func);
}

static void busv(const type_int ic)
{
  type_int idim;
  double fac = (double)grid.ngrid/(double)cp.lbox;
  long long ix, ixs[NDIM], ixe[NDIM];

  for(idim=0;idim<NDIM;idim++)
  {
    ix        = (long)(P[ic].pos[idim]*fac);
    ixs[idim] = ix - 1;
    ixe[idim] = ix + 1;
    ixs[idim] = (ixs[idim] >= (long long)grid.ngrid) ? (long long)grid.ngrid-1 : ( (ixs[idim]<0) ? 0 : ixs[idim]);
    ixe[idim] = (ixe[idim] >= (long long)grid.ngrid) ? (long long)grid.ngrid-1 : ( (ixe[idim]<0) ? 0 : ixe[idim]);
  }
  
  loopn(ic, ixs, ixe, test);

}

static void linkedlist(type_int * restrict ar)
{
  type_int i, g;
  
  Temp.ll = (type_int *) calloc(iden.nobj,sizeof(type_int));

  iden.ngrupos = 0;
  for(i=0;i<iden.nobj;i++)
  {
    ar[i] = Raiz(i,ar);
    assert(ar[i]>=i);
    if(Temp.ll[ar[i]] < NPARTMIN)
    {
      Temp.ll[ar[i]]++;
      if(Temp.ll[ar[i]]==NPARTMIN)
      { 
        iden.ngrupos++;
        Temp.ll[ar[i]] = NPARTMIN + iden.ngrupos;
      }
    }
  }

  iden.ngrupos++;  // SUMA UNO;

  Temp.head   = (type_int *) malloc(iden.ngrupos*sizeof(type_int));
  Temp.npgrup = (type_int *) malloc(iden.ngrupos*sizeof(type_int));

  for(i=0;i<iden.ngrupos;i++)
  {
    Temp.head[i]   = iden.nobj;
    Temp.npgrup[i] =  0;
  }

  for(i=0;i<iden.nobj;i++)
  {
    if(Temp.ll[ar[i]]>NPARTMIN)
    {
      ar[i] = Temp.ll[ar[i]] - NPARTMIN;
    }else{
      ar[i] = 0;
    }

    g = ar[i];

    #ifdef DEBUG
    assert((g >= 0) && (g < iden.ngrupos));
    #endif
    Temp.ll[i] = Temp.head[g];
    Temp.head[g] = i;
    Temp.npgrup[g]++;
  }

  return;
}

static void Write_Groups(const type_int niv)
{
  type_int i,j,k,npar,gn;
  char filename[200];
  FILE *pfout;

  i = iden.ngrupos-1; // LE RESTO UNO POR EL GRUPO 0 PARA ESCRIBIR EN EL ARCHIVO

  ///////////////////////////////////////////////////////
  sprintf(filename,"fof.bin");
  pfout=fopen(filename,"w");
  fwrite(&i,sizeof(type_int),1,pfout);
  //////////////////////////////////////////////////////

  npar = gn = 0;

  for(i=1;i<iden.ngrupos;i++)
  {
    j = 0;
    k = Temp.head[i];

    fwrite(&i,sizeof(type_int),1,pfout);
    fwrite(&Temp.npgrup[i],sizeof(type_int),1,pfout);    

    while(k != iden.nobj)
    {
      fwrite(&P[k].id,sizeof(type_int),1,pfout);
      k = Temp.ll[k];
      j++;
    }
    
    assert(j == Temp.npgrup[i]);

    npar+=j;
    gn++;
  }

  assert(gn == (iden.ngrupos-1));
  fclose(pfout);

  fprintf(stdout,"num de grupos %u num de particulas en grupos %u\n",gn,npar);
  fflush(stdout);

  return;
}

extern void identification(void)
{
  type_int i, tid;

  iden.nobj = cp.npart;
  iden.r0  = FOF * cp.lbox /  pow(cp.npart,1.0/(double)NDIM);

  fprintf(stdout,"Linking length = %f\n",iden.r0);
  
  gr = (type_int *) malloc(iden.nobj*sizeof(type_int));
  #ifdef LOCK
    lock = (omp_lock_t *) malloc(iden.nobj*sizeof(omp_lock_t));
  #endif
  for(i=0;i<iden.nobj;i++)
  {
    #ifdef LOCK
      omp_init_lock(&(lock[i]));
    #endif
    gr[i] = i;
  }

  grid.ngrid = (long)(cp.lbox/iden.r0);

  if(grid.ngrid > NGRIDMAX)
  {
    fprintf(stdout,"Using NGRIDMAX = %d\n",NGRIDMAX);
    grid.ngrid = NGRIDMAX;
  }
  
  grid.nobj = iden.nobj;
  grid_init();
  grid_build();

  iden.r0 *= iden.r0;

  #ifdef NTHREADS
  omp_set_dynamic(0);
  omp_set_num_threads(NTHREADS);
  #endif
    
  fprintf(stdout,"Comienza identificacion.....\n");
  fprintf(stdout,"Running on %d threads\n",NTHREADS);
  fflush(stdout);

  #ifdef LOCK
    #pragma omp parallel default(none) private(tid,i) \
    shared(P,iden,cp,lock,stdout)   
  #else
    #pragma omp parallel default(none) private(tid,i) \
    shared(P,iden,cp,stdout)   
  #endif
  {
    tid = omp_get_thread_num(); 

    for(i = tid*DIV_CEIL(iden.nobj,NTHREADS);
    i<(tid==NTHREADS-1 ? iden.nobj : (tid+1)*DIV_CEIL(iden.nobj,NTHREADS));
    i++)
    {
     
      if(i%100==0) fprintf(stdout,"%u %u %u %.4f\n",tid,i,iden.nobj,(float)i/(float)iden.nobj);

      #pragma omp task
      {
        busv(i);
      }
    }

  }  /****** TERMINA SECCION PARALELA ****************************/

  fprintf(stdout,"Sale del paralelo\n"); fflush(stdout);

  #ifdef LOCK
  for(i=0;i<iden.nobj;i++) 
    omp_destroy_lock(&(lock[i]));
  free(lock);
  #endif

  linkedlist(gr);
  Write_Groups(0);

  free(gr);
  free(Temp.ll);
  free(Temp.head);
  free(Temp.npgrup);
  grid_free();

  return;
}

