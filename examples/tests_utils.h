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
extern int pnl_eq_rel (double x, double y, double relerr, const char *label, const char *fmt, ...);
extern int pnl_eq_abs (double x, double y, double abserr, const char *label, const char *fmt, ...);
extern void pnl_init_tests ();
extern int pnl_finalize_tests ();


#endif /* TESTS_H */

