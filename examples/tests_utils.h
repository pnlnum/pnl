#ifndef _TESTS_H
#define _TESTS_H

#include <string.h>
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_sp_matrix.h"

/**
 * \defgroup PnlTest a unit test framework
 */
 /* @{  */

struct tst_list_t
{
  char *label;
  void (*func)();
};

#define MAKE_ENUM(f) {#f, f}

#ifdef DATA_DIR
  #define DATA_FILE(filename) DATA_DIR "/" filename
#endif


typedef struct tst_list_t tst_list;

extern void run_all_test (tst_list *l);
extern void menu_test (tst_list *l);

/* New unit test framework */
extern int verbose;
extern int pnl_test_eq (double x, double y, double relerr, const char *label, const char *fmt, ...);
extern int pnl_test_eq_rel (double x, double y, double relerr, const char *label, const char *fmt, ...);
extern int pnl_test_eq_abs (double x, double y, double abserr, const char *label, const char *fmt, ...);
extern int pnl_test_vect_eq_abs(const PnlVect *X, const PnlVect *Y, double relerr, const char *str, const char *fmt, ...);
extern int pnl_test_vect_complex_eq_abs (const PnlVectComplex *X, const PnlVectComplex *Y, double abserr, const char *str, const char *fmt, ...);
extern int pnl_test_mat_eq(const PnlMat *X, const PnlMat *Y, double relerr, const char *str, const char *fmt, ...);
extern int pnl_test_mat_eq_abs(const PnlMat *X, const PnlMat *Y, double relerr, const char *str, const char *fmt, ...);
extern int pnl_test_mat_int_eq(const PnlMatInt *X, const PnlMatInt *Y, const char *str, const char *fmt, ...);
extern int pnl_test_sp_mat_int_eq(const PnlSpMatInt *X, const PnlSpMatInt *Y, const char *str, const char *fmt, ...);
extern int pnl_test_mat_complex_eq_abs (const PnlMatComplex *X, const PnlMatComplex *Y, double abserr, const char *str, const char *fmt, ...);
extern int pnl_test_hmat_eq(const PnlHmat *X, const PnlHmat *Y, double relerr, const char *str, const char *fmt, ...);
extern int pnl_test_hmat_eq_rel(const PnlHmat *X, const PnlHmat *Y, double relerr, const char *str, const char *fmt, ...);
extern int pnl_test_hmat_eq_abs(const PnlHmat *X, const PnlHmat *Y, double relerr, const char *str, const char *fmt, ...);

extern void pnl_test_init(int argc, char **argv);
extern int pnl_test_finalize(const char *str);
extern void pnl_test_set_ok (const char *str);
extern void pnl_test_set_fail (const char *str, double res, double expected);
extern void pnl_test_set_fail0(const char *str);
extern int pnl_test_is_verbose ();

/* @} */

#endif /* TESTS_H */

