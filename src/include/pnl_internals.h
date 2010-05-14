#ifndef _PNL_INTERNALS_H
#define _PNL_INTERNALS_H

#define FREE(x) if (x != NULL) { free (x); x = NULL; }
#define MALLOC_DOUBLE(n) malloc ((n) * sizeof (double))
#define MALLOC_INT(n) malloc ((n) * sizeof (int))

extern double pnl_d1mach (int i);
extern double pnl_dlamch (char *);


#endif /* _PNL_INTERNALS_H */
