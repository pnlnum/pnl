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

#endif /* TESTS_H */

