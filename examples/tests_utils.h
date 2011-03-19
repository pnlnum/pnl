#ifndef _TESTS_H
#define _TESTS_H

struct tst_list_t {
    char *label;
    void (*func)();
};

#define MAKE_ENUM(f)  {#f, f}

typedef struct tst_list_t tst_list;


extern void run_all_test (tst_list *l);
extern void menu_test (tst_list *l);

/* New unit test framework */
extern int verbose;
extern int pnl_test_eq (double x, double y, double relerr, const char *label, const char *fmt, ...);
extern int pnl_test_eq_rel (double x, double y, double relerr, const char *label, const char *fmt, ...);
extern int pnl_test_eq_abs (double x, double y, double abserr, const char *label, const char *fmt, ...);
extern int pnl_test_mat_eq(const PnlMat *X, const PnlMat *Y, double relerr, const char *str, const char *fmt, ...);
extern int pnl_test_mat_eq_rel(const PnlMat *X, const PnlMat *Y, double relerr, const char *str, const char *fmt, ...);
extern int pnl_test_mat_eq_abs(const PnlMat *X, const PnlMat *Y, double relerr, const char *str, const char *fmt, ...);

extern void pnl_test_init();
extern int pnl_test_finalize();


#endif /* TESTS_H */

