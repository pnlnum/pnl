
/**************************************************/
/* Windows specific macros (not common with Unix) */
/**************************************************/
#include "configwin_specific.h"


/***********************************************/
/* Macros which are common to Windows and Unix */
/***********************************************/

/* Home directory of Premia */
#define PREMIADIR ""

/* Unix Home directory of Premia */
#define PREMIA_HOME_DIR ""

/* directory where to find the doc of Premia */
#define PREMIA_MAN_DIR ""

/* directory where to find the pdf and html doc of Premia */
#define PREMIA_MAN_PDF_DIR ""

/* directory where to find tex files of Premia */
#define PREMIA_MAN_TEX_DIR ""

/* Home directory of Premia */
#define PREMIA_SRC_DIR ""

/* Runnning on Cygwin */
#undef _CYGWIN

/* Running on a Win32 system */
#define _WIN32 1

/* Specifiy the program to use to view pdf files */
#define PDF_VIEWER_PROG "acroread"

/* #undef LAUNCH_USE_OPEN */

/* Define if have tgamma (true gamma) */
#undef HAVE_TGAMMA

/* Define if have lgamma (log gamma) */
#undef HAVE_LGAMMA

/* Define if erfl (C99) */
#undef HAVE_ERFL

/* Define if header complex.h exists (C99) */
#undef HAVE_COMPLEX_H

#if !defined(HAVE_EXP10)
#define exp10(x) pow((double) 10.0,x)
#else 
/* missing in some math.h : */
extern double exp10(double);
#endif

/* do not check indices before access */
#define PREMIA_RANGE_CHECK_OFF 1

#undef HAVE_INLINE 

