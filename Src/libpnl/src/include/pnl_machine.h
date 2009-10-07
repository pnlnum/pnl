#ifndef MACHINE_H
#define MACHINE_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "config.h"

/* Define  C2F entry point conversion */
#if defined(WTU)
#if defined(USE_SHARP_SIGN)
#define C2F(name) name##_
#else
#define C2F(name) name/**/_
#endif
#else
#define C2F(name) name
#endif


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif
