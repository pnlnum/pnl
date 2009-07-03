#ifndef _PREMIA_OBJ
#define _PREMIA_OBJ

#include "optype.h"

extern Model **premia_models[];
extern Family **premia_families[];
extern Pricing **premia_pricings[];


extern Model *models_e[];
extern Model *models_i[];
extern Model *models_g[];
extern Model *models_f[];
extern Model *models_c[];

extern Family *families_e[];
extern Family *families_i[];
extern Family *families_g[];
extern Family *families_f[];
extern Family *families_c[];

extern Pricing *pricings_e[];
extern Pricing *pricings_i[];
extern Pricing *pricings_g[];
extern Pricing *pricings_f[];
extern Pricing *pricings_c[];


typedef struct PremiaAsset {
  const char *name;
  Model **models;
  Family **families;
  Pricing **pricings;
  char label;
} PremiaAsset;

extern PremiaAsset premia_assets[];

/* set to 1 when premia is called from the command line */
extern int _interactive_call;

extern char premiasrcdir[MAX_PATH_LEN];
extern char premia_data_dir[MAX_PATH_LEN];
extern char premiapersodir[MAX_PATH_LEN];
extern const char *path_sep;

extern char premiamandir[MAX_PATH_LEN];


extern void get_exec_directory(char *exec_dir, const char *argv0);
extern void premia_self_set_global_vars(char *premia_dir);

extern int premia_set_global_vars();
#endif /* _PREMIA_OBJ */
