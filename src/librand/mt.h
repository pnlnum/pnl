#ifndef _MT_H
#define _MT_H 

#define N 624

typedef struct 
{
  unsigned long mt[N];
  int mti;
} mt_state;

extern void pnl_mt_set_seed(mt_state *state, unsigned long int s);
extern unsigned long pnl_mt_genrand (mt_state *state);
extern double pnl_mt_genrand_double (mt_state *state);

#endif /* _MT_H */
