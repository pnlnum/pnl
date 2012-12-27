
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
#include <pnl/pnl_list.h>
#include <pnl/pnl_mpi.h>


/**
 * Create an array of PnlRng from a file
 *
 * @param str The name of the file containing the generators
 * @param n the number of generators to be read
 *
 * @return an array of n rng
 */
PnlRng** pnl_rng_create_from_file (char *str, int n)
{
  PnlRng **rngtab;
  int i, size, pos;
  FILE *xdr;
  char *buf;

  size = 0;
  
  if ( (xdr = fopen (str, "rb")) == NULL ) { return NULL; }

  
  while (fgetc(xdr) != EOF) { size++; }
  buf = malloc (size);
  rewind (xdr);
  fread (buf, sizeof(char), size,  xdr);
  fclose (xdr);
  pos = 0;

  rngtab = malloc (n * sizeof(PnlRng *));
  for ( i=0 ; i<n ; i++ )
    {
      if ( pos >= size )
        {
          printf ("File %s only contains %d generators\n", str, i);
          break;
        }
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
  int i;
  FILE *stream;

  if ( (stream = fopen (str, "wb")) == NULL ) return MPI_ERR_BUFFER;

  for ( i=0 ; i<n ; i++ )
    {
      int info;
      if ((info = pnl_object_save (PNL_OBJECT(rngtab[i]), stream))) return info;
    }
  fclose (stream);
  return MPI_SUCCESS;
}

/**
 * Save a PnlObject into a stream
 *
 * @param obj a PnlObject
 * @param stream the stream obtained when opening a file. The file should
 * opened as a binary file
 *
 * @return an error code equal to MPI_SUCCESS when the function succeeds
 */
int pnl_object_save (PnlObject *obj, FILE *stream)
{
  int size, pos, info;
  char *buf;

  pos = 0; size = 0;
  if ((info = pnl_object_mpi_pack_size (PNL_OBJECT(obj), MPI_COMM_WORLD, &size))) return info;
  if ((buf = malloc (size)) == NULL) return MPI_ERR_BUFFER;
  if ((info = pnl_object_mpi_pack (PNL_OBJECT(obj), buf, size, &pos, MPI_COMM_WORLD)))
    {
      free (buf); return info;
    }
  info = fwrite (buf, sizeof(char), pos, stream);
  free (buf);
  if ( info < pos ) return MPI_ERR_BUFFER;
  return MPI_SUCCESS;
}

/**
 * Load an object from a string buffer
 *
 * @param buf a string buffer
 * @param bufsize the size of the buffer in bytes
 * @param pos (input/output) current position in buf
 * @return a PnlObject
 */
static PnlObject* load_from_buf (char *buf, int bufsize, int *pos)
{
  PnlObject *O;
  int parent_id, id, savepos;

  savepos = *pos;
  if (MPI_Unpack(buf, bufsize, pos, &parent_id, 1, MPI_INT, MPI_COMM_WORLD)
      || MPI_Unpack(buf, bufsize, pos, &id, 1, MPI_INT, MPI_COMM_WORLD))
    {
      printf ("Cannot read any correct PnlTypes in pnl_object_load\n");
      return NULL;
    }
  O = pnl_object_create(id);
  *pos = savepos;
  pnl_object_mpi_unpack(O, buf, bufsize, pos, MPI_COMM_WORLD);
  return O;
}

/**
 * Load a object from a stream
 *
 * @param stream the stream obtained when opening a file. The file should
 * opened as a binary file
 *
 * @return a PnlObject or NULL if the stream was empty or it did not
 * contain any PnlObject
 */
PnlObject* pnl_object_load (FILE *stream)
{
  PnlObject *O;
  long pos_in_stream;
  int  bufsize, pos;
  char *buf;

  pos_in_stream = ftell(stream);
  bufsize = 0;
  while (fgetc(stream) != EOF) { bufsize++; }
  if ( bufsize == 0 ) return NULL;
  buf = malloc (bufsize);
  fseek(stream, pos_in_stream, SEEK_SET);
  fread (buf, 1, bufsize,  stream);

  pos = 0;
  O = load_from_buf (buf, bufsize, &pos);
  free (buf);
  fseek(stream, pos_in_stream + pos, SEEK_SET);
  return O;
}

/**
 * Load objects from a stream and stores them into a PnlList
 *
 * @param stream the stream obtained when opening a file. The file should
 * opened as a binary file
 *
 * @return a PnlList or NULL if the stream was empty or it did not
 * contain any PnlObjects
 */
PnlList* pnl_object_load_into_list (FILE *stream)
{

  PnlList *L;
  long     pos_in_stream;
  int      bufsize, pos;
  char    *buf;


  pos_in_stream = ftell(stream);
  bufsize = 0;
  while (fgetc(stream) != EOF) { bufsize++; }
  if ( bufsize == 0 ) return NULL;
  buf = malloc (bufsize);
  fseek(stream, pos_in_stream, SEEK_SET);
  fread (buf, 1, bufsize,  stream);

  pos = 0;
  L = pnl_list_new ();
  while ( pos < bufsize )
    {
      PnlObject *O;
      if ((O = load_from_buf (buf, bufsize, &pos)) == NULL)
        {
          free (buf); break;
        }
      pnl_list_insert_last (L, O);
    }
  free (buf);
  fseek(stream, pos_in_stream + pos, SEEK_SET);
  if (L->len == 0)
    {
      pnl_list_free(&L); L=NULL;
    }
  return L;
}

