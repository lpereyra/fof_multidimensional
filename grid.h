#ifndef GRID_H
#define GRID_H

#ifndef NGRIDMAX
  #define NGRIDMAX 512
#endif

extern void grid_init(void);
extern void grid_build(void);
extern void grid_free(void);

#endif
