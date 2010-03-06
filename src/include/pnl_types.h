#ifndef _PNL_TYPES_H
#define _PNL_TYPES_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */



#define NULLINT -1
/* when label=NULL, key must be set to NULLINT */
typedef struct
{
    const char * label;
    int          key;
} enum_member;

/* the array of enum_member must end with {NULL, NULLINT} */
typedef struct 
{
    unsigned      size;     /*!< size in bytes of an enum member */
    enum_member * members;  /*!< a pointer to the first member of the enum */
    unsigned      count;    /*!< number of members in the array */
} enum_members;

typedef struct 
{
    int            value;
    enum_members * members;
} enumeration;

#define DEFINE_ENUM(Name, Members)  enum_members Name = { sizeof(Members[0]), &Members[0], sizeof(Members)/sizeof(Members[0])-1 };

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_TYPES_H */
