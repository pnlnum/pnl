/* to be put first, because optype.h redefines OUT IN already defined in windows.h */
#ifdef _MSC_VER
#include <windows.h> // Usefull for LPSTR
#undef IN
#undef OUT
#else
#include <unistd.h> // Usefull for getcwd()
#include <dirent.h>
#endif

#include <sys/types.h>

#include "optype.h"
#include "var.h"
#include "premia_obj.h"
#include "method.h"
#include "test.h"
#include "timeinfo.h"
#include "tools.h"
#include "ftools.h"
#include "pnl_random.h"
#include "config.h"

/* Models Equities*/
extern Model BS1D_model;
extern Model BSDISDIV1D_model;
extern Model BS2D_model;
extern Model BSND_model;
extern Model DUP1D_model;
extern Model MER1D_model;
extern Model HW1D_model;
extern Model HES1D_model;
extern Model TIMEHES1D_model;
extern Model DOUBLEHES1D_model;
extern Model SCOTT1D_model;
extern Model STEIN1D_model;
extern Model BNS_model;
extern Model DPS_model;
extern Model FPS1D_model;
extern Model FPS2D_model;
extern Model WISHART2D_model;
extern Model MERHES1D_model;
extern Model PUREJUMP1D_model;
extern Model VARIANCEGAMMA1D_model;
extern Model NIG1D_model;
extern Model KOU1D_model;
extern Model RSKOU1D_model;
extern Model TEMPEREDSTABLE1D_model;
extern Model RSTEMPEREDSTABLE1D_model;
extern Model CGMY1D_model;
extern Model JUMP1D_model;
extern Model NONPAR1D_model;
extern Model VARSWAP3D_model;

#ifdef XNSP 
extern Model NSPTEST_model;
#endif 

Model  *models_e[]=
  {
    &BS1D_model,
    &BSDISDIV1D_model,
    &BS2D_model,
    &BSND_model,
    &DUP1D_model,
    &MER1D_model,
    &PUREJUMP1D_model,
    &VARIANCEGAMMA1D_model,
    &NIG1D_model,
    &KOU1D_model,
    &RSKOU1D_model,
    &TEMPEREDSTABLE1D_model,
    &RSTEMPEREDSTABLE1D_model,
    &CGMY1D_model,
    &HW1D_model,
    &HES1D_model,
    &TIMEHES1D_model,
    &DOUBLEHES1D_model,
    &SCOTT1D_model,
    &STEIN1D_model,
    &FPS1D_model,
    &FPS2D_model,
    &WISHART2D_model,
    &MERHES1D_model,
    &BNS_model,
    &DPS_model,
    &NONPAR1D_model,
    &VARSWAP3D_model,
#ifdef XNSP 
    &NSPTEST_model,
#endif 
    NULL
  };

/* Models Energy*/
extern Model JUMP1D_model;

Model  *models_g[]=
  {
    &JUMP1D_model,
    NULL
  };

/* Models Inflation*/
extern Model  JarrowYildirim1D_model;
extern Model INFLATION_LMM_HESTON1D_model;
Model  *models_f[]=
  {
    & JarrowYildirim1D_model,
    & INFLATION_LMM_HESTON1D_model,
    NULL
  };

/* Models Insurance*/
extern Model BSCIR2D_model;

Model  *models_a[]=
  {
    &BSCIR2D_model,
    NULL
  };


/* Models Interest Rates*/
extern Model Vasicek1D_model;
extern Model Cir1D_model;
extern Model HullWhite1D_model;
extern Model CirPP1D_model;
extern Model BlackKarasinski1D_model;
extern Model SG1D_model;
extern Model LRSHJM1D_model;
extern Model BharChiarella1D_model;
extern Model HK1D_model;
extern Model LMM1D_model;
extern Model LMM_HESTON1D_model;
extern Model LMM_JUMP1D_model;
extern Model HullWhite2D_model;
extern Model Cir2D_model;
extern Model QTSM2D_model;
extern Model Affine3D_model;

Model *models_i[]=
  {
    &Vasicek1D_model,
    &Cir1D_model,
    &HullWhite1D_model,
    &CirPP1D_model,
    &BlackKarasinski1D_model,
    &SG1D_model,
    &LRSHJM1D_model,
    &BharChiarella1D_model,
    &HK1D_model,
    &LMM1D_model,
    &LMM_HESTON1D_model,
    &LMM_JUMP1D_model,
    &HullWhite2D_model,
    &Cir2D_model,
    &QTSM2D_model,
    &Affine3D_model,
    NULL

  };

/* Models Credit*/
extern Model CirPP2D_model;
extern Model COPULA_model;
extern Model DYNAMIC_CDO_model;

Model *models_c[]=
  {
    &CirPP2D_model,
    &COPULA_model,
    &DYNAMIC_CDO_model,
    NULL
  };


extern Family STD_family;
extern Family LIM_family;
extern Family LIMDISC_family;
extern Family DOUBLIM_family;
extern Family PAD_family;
extern Family VOL_family;
extern Family STD2D_family;
extern Family STDND_family;
extern Family STDi_family;
extern Family STDf_family;
extern Family STDg_family;
extern Family STDc_family;
extern Family STDNDc_family;
extern Family STDa_family;

Family *families_e[]=
  {
    &STD_family,
    &LIM_family,
    &LIMDISC_family,
    &DOUBLIM_family,
    &PAD_family,
    &VOL_family,
    &STD2D_family,
    &STDND_family,
    NULL
  };

Family *families_i[]=
  {
    &STDi_family,
    NULL
  };

Family *families_f[]=
  {
    &STDf_family,
    NULL
  };

Family *families_g[]=
  {
    &STDg_family,
    NULL
  };

Family *families_c[]=
  {
    &STDc_family,
    &STDNDc_family,
    NULL
  };

Family *families_a[]=
  {
    &STDa_family,
    NULL
  };

extern  Pricing BS1D_STD_pricing;
extern  Pricing BS1D_LIM_pricing;
extern  Pricing BS1D_LIMDISC_pricing;
extern  Pricing BS1D_PAD_pricing;
extern  Pricing BS1D_DOUBLIM_pricing;
extern  Pricing BSDISDIV1D_STD_pricing;
extern  Pricing BS2D_STD2D_pricing;
extern  Pricing BSND_STDND_pricing;
extern  Pricing DUP1D_STD_pricing;
extern  Pricing MER1D_STD_pricing;
extern  Pricing MER1D_LIM_pricing;
extern  Pricing MER1D_PAD_pricing;
extern  Pricing HW1D_STD_pricing;
extern  Pricing HES1D_STD_pricing;
extern  Pricing HES1D_LIM_pricing;
extern  Pricing HES1D_PAD_pricing;
extern  Pricing HES1D_VOL_pricing;
extern  Pricing DOUBLEHES1D_STD_pricing;
extern  Pricing DOUBLEHES1D_VOL_pricing;
extern  Pricing TIMEHES1D_STD_pricing;
extern  Pricing SCOTT1D_STD_pricing;
extern  Pricing STEIN1D_STD_pricing;
extern  Pricing BNS_STD_pricing;
extern  Pricing DPS_STD_pricing;
extern  Pricing FPS1D_STD_pricing;
extern  Pricing FPS2D_STD_pricing;
extern  Pricing WISHART2D_STD2D_pricing;
extern  Pricing MERHES1D_STD_pricing;
extern  Pricing MERHES1D_LIM_pricing;
extern  Pricing MERHES1D_PAD_pricing;
extern  Pricing PUREJUMP1D_PAD_pricing;
extern  Pricing VARIANCEGAMMA1D_STD_pricing;
extern  Pricing VARIANCEGAMMA1D_LIM_pricing;
extern  Pricing VARIANCEGAMMA1D_PAD_pricing;
extern  Pricing NIG1D_STD_pricing;
extern  Pricing NIG1D_LIM_pricing;
extern  Pricing NIG1D_PAD_pricing;
extern  Pricing KOU1D_STD_pricing;
extern  Pricing KOU1D_LIM_pricing;
extern  Pricing KOU1D_PAD_pricing;
extern  Pricing RSKOU1D_STD_pricing;
extern  Pricing RSKOU1D_LIM_pricing;
extern  Pricing TEMPEREDSTABLE1D_STD_pricing;
extern  Pricing TEMPEREDSTABLE1D_LIM_pricing;
extern  Pricing TEMPEREDSTABLE1D_VOL_pricing;
extern  Pricing RSTEMPEREDSTABLE1D_STD_pricing;
extern  Pricing RSTEMPEREDSTABLE1D_LIM_pricing;
extern  Pricing CGMY1D_STD_pricing;
extern  Pricing CGMY1D_PAD_pricing;
extern  Pricing NONPAR1D_VOL_pricing;
extern  Pricing VARSWAP3D_STD_pricing;

Pricing *pricings_e[]=
  {
    &BS1D_STD_pricing,
    &BS1D_LIM_pricing,
    &BS1D_LIMDISC_pricing,
    &BS1D_PAD_pricing,
    &BS1D_DOUBLIM_pricing,
    &BSDISDIV1D_STD_pricing,
    &BS2D_STD2D_pricing,
    &BSND_STDND_pricing,
    &PUREJUMP1D_PAD_pricing,
    &DUP1D_STD_pricing,
    &MER1D_STD_pricing,
    &MER1D_LIM_pricing,
    &MER1D_PAD_pricing,
    &HW1D_STD_pricing,
    &HES1D_STD_pricing,
    &HES1D_LIM_pricing,
    &HES1D_PAD_pricing,
    &HES1D_VOL_pricing,
    &DOUBLEHES1D_STD_pricing,
    &DOUBLEHES1D_VOL_pricing,
    &TIMEHES1D_STD_pricing,
    &SCOTT1D_STD_pricing,
    &STEIN1D_STD_pricing,
    &BNS_STD_pricing,
    &DPS_STD_pricing,
    &FPS1D_STD_pricing,
    &FPS2D_STD_pricing,
    &WISHART2D_STD2D_pricing,
    &MERHES1D_STD_pricing,
    &MERHES1D_LIM_pricing,
    &MERHES1D_PAD_pricing,
    &VARIANCEGAMMA1D_STD_pricing,
    &VARIANCEGAMMA1D_LIM_pricing,
    &VARIANCEGAMMA1D_PAD_pricing,
    &NIG1D_STD_pricing,
    &NIG1D_LIM_pricing,
    &NIG1D_PAD_pricing,
    &KOU1D_STD_pricing,
    &KOU1D_LIM_pricing,
    &KOU1D_PAD_pricing,
    &RSKOU1D_STD_pricing,
    &RSKOU1D_LIM_pricing,
    &TEMPEREDSTABLE1D_STD_pricing,
    &TEMPEREDSTABLE1D_LIM_pricing,
    &TEMPEREDSTABLE1D_VOL_pricing,
    &RSTEMPEREDSTABLE1D_STD_pricing,
    &RSTEMPEREDSTABLE1D_LIM_pricing,
    &CGMY1D_STD_pricing,
    &CGMY1D_PAD_pricing,
    &NONPAR1D_VOL_pricing,
    &VARSWAP3D_STD_pricing,
    NULL
  };

extern  Pricing Vasicek1D_STDi_pricing;
extern  Pricing Cir1D_STDi_pricing;
extern  Pricing HullWhite1D_STDi_pricing;
extern  Pricing CirPP1D_STDi_pricing;
extern  Pricing BlackKarasinski1D_STDi_pricing;
extern  Pricing SG1D_STDi_pricing;
extern  Pricing LRSHJM1D_STDi_pricing;
extern  Pricing BharChiarella1D_STDi_pricing;
extern  Pricing HK1D_STDi_pricing;
extern  Pricing LMM1D_STDi_pricing;
extern  Pricing LMM_HESTON1D_STDi_pricing;
extern  Pricing LMM_JUMP1D_STDi_pricing;
extern  Pricing HullWhite2D_STDi_pricing;
extern  Pricing Cir2D_STDi_pricing;
extern  Pricing QTSM2D_STDi_pricing;
extern  Pricing Affine3D_STDi_pricing;

Pricing *pricings_i[]=
  {
    &Vasicek1D_STDi_pricing,
    &Cir1D_STDi_pricing,
    &HullWhite1D_STDi_pricing,
    &CirPP1D_STDi_pricing,
    &BlackKarasinski1D_STDi_pricing,
    &SG1D_STDi_pricing,
    &LRSHJM1D_STDi_pricing,
    &BharChiarella1D_STDi_pricing,
    &HK1D_STDi_pricing,
    &LMM1D_STDi_pricing,
    &LMM_HESTON1D_STDi_pricing,
    &LMM_JUMP1D_STDi_pricing,
    &HullWhite2D_STDi_pricing,
    &Cir2D_STDi_pricing,
    &QTSM2D_STDi_pricing,
    &Affine3D_STDi_pricing,
    NULL
  };

extern  Pricing CirPP2D_STDc_pricing;
extern  Pricing COPULA_STDNDc_pricing;
extern  Pricing DYNAMIC_CDO_STDNDc_pricing;



Pricing *pricings_c[]=
  {
    &CirPP2D_STDc_pricing,
    &COPULA_STDNDc_pricing,
    &DYNAMIC_CDO_STDNDc_pricing,
    NULL
  };

extern  Pricing JUMP1D_STDg_pricing;
Pricing *pricings_g[]=
  {
    &JUMP1D_STDg_pricing, 
    NULL
  };

extern  Pricing BSCIR2D_STDa_pricing;
Pricing *pricings_a[]=
  {
    &BSCIR2D_STDa_pricing, 
    NULL
  };

extern  Pricing JarrowYildirim1D_STDf_pricing;
extern  Pricing INFLATION_LMM_HESTON1D_STDf_pricing;
Pricing *pricings_f[]=
  {
    &JarrowYildirim1D_STDf_pricing,
    &INFLATION_LMM_HESTON1D_STDf_pricing,
    NULL
  };

PremiaAsset premia_assets[] =
  {
    {"equity", models_e, families_e, pricings_e, 'e'},
    {"interest", models_i, families_i, pricings_i, 'i'},
    {"credit", models_c, families_c, pricings_c, 'c'},
#if !(defined(PremiaCurrentVersion) && PremiaCurrentVersion < (2007+2))
    {"inflation", models_f, families_f, pricings_f, 'f'},
    {"energy", models_g, families_g, pricings_g, 'g'},
#endif
#if !(defined(PremiaCurrentVersion) && PremiaCurrentVersion < (2010+2))
    {"insurance", models_a, families_a, pricings_a, 'a'},
#endif
    {NULL, NULL, NULL, NULL},
  };


int _interactive_call=0;

extern TimeInfo computation_time_info;

FILE *out_stream=NULL;

char premiasrcdir[MAX_PATH_LEN];
char premia_data_dir[MAX_PATH_LEN];
char premiapersodir[MAX_PATH_LEN];
char premiamandir[MAX_PATH_LEN];

char PREMIA_OUT[MAX_PATH_LEN]="";
char GNUPLOT_DAT[MAX_PATH_LEN]="";
char TITLES_TEX[MAX_PATH_LEN]="";
char GNUPLOT_SCREEN_PLT[MAX_PATH_LEN]="";
char GNUPLOT_FILE_PLT[MAX_PATH_LEN]="";
char GNU_TEX[MAX_PATH_LEN]="";
char PREMIA_LOG[MAX_PATH_LEN]="";
char SESSION_LOG[MAX_PATH_LEN]="";

#ifdef _WIN32
const char *path_sep = "\\";
#else
const char *path_sep = "/";
#endif

/**********************************************************************
In a string (str_from), this function replaces all occurences
of a substring (str_find) with a smaller string (str_replace)
**********************************************************************/
void replace_smaller_substring_all(char *str_from, const char *str_find, const char *str_replace)
{
  char *ptr_begin;

  size_t len_find, len_replace;

  len_find    = strlen(str_find);
  len_replace = strlen(str_replace);

  if (len_replace <= len_find)
    {
      while ((ptr_begin=strstr(str_from,str_find)) != NULL)
        {
          memmove(ptr_begin,               str_replace,            len_replace);
          memmove(ptr_begin + len_replace, ptr_begin + len_find, strlen(ptr_begin) + len_find + 1); // +1 is to include the '\0' at the end
        }
    }
}


/**********************************************************************
Copies the full path of the current executable into exec_dir
argv[0] must be passed to the function
It is advised to call this function at the beginning of main()
  so that the directory from where the user launched the executable
  still corresponds to the current work directory returned by getcwd.
**********************************************************************/
void get_exec_directory(char *exec_dir, const char *argv0)
{
  char *tail_char;

#ifdef _WIN32
  /* On windows argv[0] always represents the full path to the executable */
  strcpy(exec_dir, argv0);
#else
  char current_work_dir[MAX_PATH_LEN]="";
  char exec_path[MAX_PATH_LEN]="";
  char str_to_replace[4]="";

  /* If argv[0] starts with "./", we remove it */
  if ((*argv0 == '.') && (*(argv0+1) == *path_sep)) strcpy(exec_path, (argv0+2));
  else strcpy(exec_path, argv0);

  /* We replace the string "/./" by "/" anywhere in argv[0]
     This is usefull if we want later on to move up and down in the path */
  sprintf(str_to_replace, "%c.%c", *path_sep, *path_sep);
  replace_smaller_substring_all(exec_path, str_to_replace, path_sep);

  if (getcwd(current_work_dir, MAX_PATH_LEN) == NULL) 
    {
      perror("Error from getcwd in function get_exec_directory");
      abort();
    }

  /* Usually getcwd returns a path with no "/" at the end
     but in case it does we remove it (by replacing it with '\0').
     It is at least the case when the current work directory is the root (just "/") */
  if (strlen(current_work_dir) >= 1)
    {
      tail_char = current_work_dir + strlen(current_work_dir) - 1;
      if (*tail_char == *path_sep) *tail_char = '\0';
    }

  /* If argv[0] is a relative path (not starting with "/")
     we add in front the current work directory */
  if (*exec_path != *path_sep) sprintf(exec_dir, "%s%c%s", current_work_dir, *path_sep, exec_path);
  /* If argv[0] is already the full path. Note that on some Unix shells,
     argv[0] will contain the full path if you have typed ".." in the command */
  else sprintf(exec_dir, "%s", exec_path);
#endif

  /* We move the end of the string to the last "/" so to remove the executable's name */
  if ((tail_char = strrchr (exec_dir, *path_sep)) != NULL) *tail_char = '\0';
}


/**********************************************************************
Sets global variables corresponding to files locations
**********************************************************************/
static void premia_set_global_file_vars()
{
  sprintf(PREMIA_OUT,"premia.out");
  sprintf(GNUPLOT_DAT,"%s%sgnupremia.dat",premiapersodir,path_sep);
  sprintf(TITLES_TEX,"titles.tex");
  sprintf(GNUPLOT_SCREEN_PLT,"%s%sgnutoscreen.plt",premiapersodir,path_sep);
  sprintf(GNUPLOT_FILE_PLT,"%s%sgnutofile.plt",premiapersodir,path_sep);
  sprintf(GNU_TEX,"%s%sgnupremia.tex",premiapersodir,path_sep);
  sprintf(PREMIA_LOG,"%s%spremia.log",premiapersodir,path_sep);
  sprintf(SESSION_LOG,"%s%ssession.log",premiapersodir,path_sep);      
}


/**********************************************************************
Sets global variables corresponding to directories and files locations
This is a new version of premia_set_global_vars(), the difference
  being that premia_dir is self determined (at execution time),
  instead of being determined from ./configure
**********************************************************************/
void premia_self_set_global_vars(char *premia_dir)
{
  char *tail_char;

  /* we go up one directory (the executable is in the "Unix" directory) to get the Premia root dir */
  if ((tail_char = strrchr (premia_dir, *path_sep)) != NULL) *tail_char = '\0';

  sprintf(premiamandir,    "%s%c%s%c%s", premia_dir, *path_sep, "man", *path_sep, "pdf_html");
  sprintf(premiasrcdir,    "%s%c%s",     premia_dir, *path_sep, "Src");
  sprintf(premiapersodir,  "%s%c%s",     premia_dir, *path_sep, "Unix");
  sprintf(premia_data_dir, "%s%c%s",     premia_dir, *path_sep, "data");

  premia_set_global_file_vars();
}


/**********************************************************************
Sets global variables corresponding to directories and files locations
This is still usefull for the NSP version of Premia
**********************************************************************/
int premia_set_global_vars()
{
#if defined(__MINGW32__) || defined(__APPLE__) || defined(__CYGWIN__) || defined(netbsd) || defined(freebsd) || defined (linux)
  char *nsp_dir;
  DIR *dir;
  nsp_dir = NULL;
#endif

#ifdef _MSC_VER
  LPSTR tail;
      
  GetModuleFileName(GetModuleHandle(NULL),premiapersodir, 256);
  printf("---> %s\n", premiapersodir); 

  /* remove executable name */
  if ((tail = strrchr (premiapersodir, '\\')) != (LPSTR) NULL)
    *tail = '\0';
  /* remove last directory component */
  if ((tail = strrchr (premiapersodir, '\\')) != (LPSTR) NULL)
    *tail = '\0';

  printf("---> %s\n", premiapersodir); 
  strcpy(premia_data_dir, premiapersodir);
  
  if ((tail = strrchr (premia_data_dir, '\\')) != (LPSTR) NULL)
    *tail = '\0';
  
  sprintf(premiamandir, "%s\\man\\pdf_html", premia_data_dir );
  sprintf(premia_data_dir, "%s\\data", premiapersodir);
#else /* _MSC_VER */

  sprintf(premiamandir, "%s", PREMIA_MAN_PDF_DIR);
  sprintf(premiasrcdir, "%s", PREMIA_SRC_DIR);
  sprintf(premiapersodir, "%s", PREMIA_HOME_DIR);
  sprintf(premia_data_dir, "%s%sdata", PREMIADIR, path_sep);
#endif /* _MSC_VER */

#if defined(__MINGW32__) || defined(__APPLE__) || defined(__CYGWIN__) || defined(netbsd) || defined(freebsd) || defined (linux)
  if ((dir = opendir (premiamandir)) == NULL)
    {
      if ((nsp_dir = getenv ("NSP")) == NULL) return FAIL;
      sprintf(premiamandir, "%s%spremia-bin%sman%spdf_html", nsp_dir, path_sep, path_sep, path_sep);
    }
  else closedir (dir);


  if ((dir = opendir (premia_data_dir)) == NULL)
    {
      if (nsp_dir == NULL && (nsp_dir = getenv ("NSP")) == NULL) return FAIL;
      sprintf(premia_data_dir, "%s%spremia-bin%sdata", nsp_dir, path_sep, path_sep);
    }
  else  closedir (dir);

#endif
  premia_set_global_file_vars();
  return OK;
}
