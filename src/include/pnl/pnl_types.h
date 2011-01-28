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
} PnlEnum;


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* _PNL_TYPES_H */
