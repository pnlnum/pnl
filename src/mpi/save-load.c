
/************************************************************************/
/* Copyright Jérôme Lelong <jerome.lelong@gmail.com>                    */
/*                                                                      */
/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU Lesser General Public License as       */
/* published by the Free Software Foundation, either version 3 of the   */
/* License, or (at your option) any later version.                      */
/*                                                                      */
/* This program is distributed in the hope that it will be useful, but  */
/* WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    */
/* Lesser General Public License for more details.                      */
/*                                                                      */
/* You should have received a copy of the GNU Lesser General Public     */
/* License along with this program.  If not, see                        */
/* <http://www.gnu.org/licenses/>.                                      */
/************************************************************************/

#include <pnl/pnl_random.h>
#include <pnl/pnl_mpi.h>




/**
 * Creates an array of PnlRng from a file
 * 
 * @param str The name of the file containing the generators
 * @param n the number of generators to be read
 * 
 * @return an array of n rng
 */
PnlRng ** pnl_rng_create_from_file (char *str, int n)
{
  PnlRng **rngtab;
  int i, count, size, pos;
  FILE *xdr;
  char *buf;


  rngtab = malloc (n * sizeof(PnlRng *));
  
  size = 0;
  xdr = fopen (str, "rb");
  while (fgetc(xdr) != EOF) { size++; }
  buf = malloc (size);
  rewind (xdr);
  fread (buf, sizeof(char), size,  xdr);
  fclose (xdr);
  pos = 0;

  /*
   * Read  the total number of generators in the file
   */
  MPI_Unpack (buf, size, &pos, &count, 1, MPI_INT, MPI_COMM_WORLD);
  if ( n > count )
    {
      printf ("File %s only contains %d generators\n", str, count);
      free (buf); return NULL;
    }
  for ( i=0 ; i<n ; i++ )
    {
      rngtab[i] = pnl_rng_new();      
      pnl_object_mpi_unpack (PNL_OBJECT(rngtab[i]), buf, size, &pos, MPI_COMM_WORLD);
    }
  free (buf);
  return rngtab;
}

/**
 * Save an array of rng to a file
 * 
 * @param rngtab an array of rng
 * @param n the number of rng to save, must be less or equal than the number
 * of elements of rngtab
 * @param str the name of the file to which saving rngtab
 * 
 * @return an MPI error code
 */
int pnl_rng_save_to_file (PnlRng **rngtab, int n, char *str)
{
  int i, size, pos, info;
  char *buf;
  FILE *xdr;

  /*
   * Writes the total number of generators
   */
  xdr = fopen (str, "wb");
  size = 0; pos = 0;
  if ((info = MPI_Pack_size (1, MPI_INT, MPI_COMM_WORLD, &size))) return info;
  if ((buf = malloc (size)) == NULL) return MPI_ERR_BUFFER;
  if ((info = MPI_Pack (&n, 1, MPI_INT, buf, size, &pos, MPI_COMM_WORLD))) return info;
  fwrite (buf, sizeof(char), pos, xdr);
  free (buf);

  /*
   * Writes each generator
   */
  for ( i=0 ; i<n ; i++ )
    {
      pos = 0;
      if ((info = pnl_object_mpi_pack_size (PNL_OBJECT(rngtab[i]), MPI_COMM_WORLD, &size))) return info;
      if ((buf = malloc (size)) == NULL) return MPI_ERR_BUFFER;
      if ((info = pnl_object_mpi_pack (PNL_OBJECT(rngtab[i]), buf, size, &pos, MPI_COMM_WORLD))) return info;
      fwrite (buf, sizeof(char), pos, xdr);
      free (buf);
    }
  fclose (xdr);
  return MPI_SUCCESS;
}
