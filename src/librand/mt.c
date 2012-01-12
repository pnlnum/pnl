/*
 * Written and (C) by Jérôme Lelong <jerome.lelong@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as 
 * published by  the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License  along with this program.  If not, see
 * <http://www.gnu.org/licenses/>. 
 *
 * UPDATE:
 *  Original implementation modified to ensure unsigned long
 *  is 32-bit even on 64-bit machines.
 *
 * This implementation is based on the following
 * 
 * A C-program for MT19937: Integer version (1999/10/28)         
 *  genrand() generates one pseudorandom unsigned integer (32bit)
 * which is uniformly distributed among 0 to 2^32-1  for each    
 * call. sgenrand(seed) sets initial values to the working area  
 * of 624 words. Before genrand(), sgenrand(seed) must be        
 * called once. (seed is any 32-bit integer.)                    
 *   Coded by Takuji Nishimura, considering the suggestions by   
 * Topher Cooper and Marc Rieffel in July-Aug. 1997.             
 * 
 * This library is free software; you can redistribute it and/or  
 * modify it under the terms of the GNU Library General Public    
 * License as published by the Free Software Foundation; either   
 * version 2 of the License, or (at your option) any later        
 * version.                                                       
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.           
 * See the GNU Library General Public License for more details.   
 * You should have received a copy of the GNU Library General     
 * Public License along with this library; if not, write to the   
 * Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  
 * 02111-1307  USA                                                
 * 
 * Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura.
 * Copyright (C) 2009 Makoto Matsumoto, Takuji Nishimura and      
 * Mutsuo Saito.                                                  
 * When you use this, send an email to: matumoto@math.keio.ac.jp  
 * with an appropriate reference to your work.                    
 * 
 * REFERENCE                                                      
 * M. Matsumoto and T. Nishimura,                                 
 * "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform 
 * Pseudo-Random Number Generator",                               
 * ACM Transactions on Modeling and Computer Simulation,          
 * Vol. 8, No. 1, January 1998, pp 3--30.                         
 */

#include "pnl/pnl_random.h"


/* #define PNL_MT_N 624, in mt.h */
#define PNL_MT_M 397
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */
#define BIT32_MASK 0xffffffffUL


/* Initializing the array with a seed */
void pnl_mt_sseed(mt_state *state, unsigned long int s)
{
  int i;

  if ( s == 0 ) s = 4357; /* this is the default seed */

  state->mt[0] = s & BIT32_MASK;

  for ( i=1 ; i<PNL_MT_N ; i++ ) 
    {
      state->mt[i] = (1812433253UL * (state->mt[i-1]  ^ (state->mt[i-1] >> 30))) + i;
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array mt[].                        */
      state->mt[i] &= BIT32_MASK;
    }
  state->mti = i;
}

#define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)
/* Draws the next number in the sequence. 
 * It returns an unsigned long to which we apply a 32-bit
 * mask to make sure that the returned value is 32-bit long
 * even on 64-bit machine
 */
unsigned long pnl_mt_genrand (mt_state *state)
{
  unsigned long y;

  if (state->mti >= PNL_MT_N) 
    { /* generate PNL_MT_N words at one time */
      int kk;

      for ( kk=0 ; kk<PNL_MT_N-PNL_MT_M ; kk++ ) 
        {
          y = (state->mt[kk] & UPPER_MASK) | (state->mt[kk+1] & LOWER_MASK);
          state->mt[kk] = state->mt[kk+PNL_MT_M] ^ (y >> 1) ^ MAGIC(y);
        }
      for (; kk<PNL_MT_N-1 ; kk++) 
        {
          y = (state->mt[kk] & UPPER_MASK) | (state->mt[kk+1] & LOWER_MASK);
          state->mt[kk] = state->mt[kk+(PNL_MT_M-PNL_MT_N)] ^ (y >> 1) ^ MAGIC(y);
        }
      y = (state->mt[PNL_MT_N-1] & UPPER_MASK) | (state->mt[0] & LOWER_MASK);
      state->mt[PNL_MT_N-1] = state->mt[PNL_MT_M-1] ^ (y >> 1) ^ MAGIC(y);

      state->mti = 0;
    }

  y = state->mt[state->mti];
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  state->mti++;
  return y;
}

double pnl_mt_genrand_double (mt_state *state)
{
  /* normalisation by 2^32  */
  return pnl_mt_genrand (state) / 4294967296.0;
}

