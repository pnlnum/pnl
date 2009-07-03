#ifndef _TESTS_H
#define _TESTS_H

struct list_t {
    int id;
    char *label;
    void (*func)();
};

#define MAKE_ENUM(id, f)  {id, #f, f}
#define NULL_INT -1

typedef struct list_t list;


#endif /* TESTS_H */

