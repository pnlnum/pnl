#ifndef _TESTS_H
#define _TESTS_H

struct tst_list_t {
    int id;
    char *label;
    void (*func)();
};

#define MAKE_ENUM(id, f)  {id, #f, f}
#define NULL_INT -1

typedef struct tst_list_t tst_list;


extern void run_all_test (tst_list *l);
extern void run_all_test (tst_list *l);

#endif /* TESTS_H */

