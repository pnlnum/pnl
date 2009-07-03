
/* to be put first, because optype.h redefines OUT IN
 already defined in windows.h */
#if defined(_WIN32) && !defined(_CYGWIN)
#include <windows.h>
#undef IN
#undef OUT
#endif


#include "optype.h"
#include "var.h"
#include "premia_obj.h"
#include "method.h"
#include "test.h"
#include "timeinfo.h"

#include "tools.h"
#include "ftools.h"
#include "error_msg.h"
#include "pnl_random.h"

#include "config.h"

extern TimeInfo computation_time_info;

extern FILE *out_stream;
extern char premiasrcdir[MAX_PATH_LEN];
extern char premiadir[MAX_PATH_LEN];

extern const char *path_sep;
extern char premiamandir[MAX_PATH_LEN];

extern char PREMIA_OUT[MAX_PATH_LEN];
extern char GNUPLOT_DAT[MAX_PATH_LEN];
extern char TITLES_TEX[MAX_PATH_LEN];
extern char GNUPLOT_SCREEN_PLT[MAX_PATH_LEN];
extern char GNUPLOT_FILE_PLT[MAX_PATH_LEN];
extern char GNU_TEX[MAX_PATH_LEN];
extern char PREMIA_LOG[MAX_PATH_LEN];
extern char SESSION_LOG[MAX_PATH_LEN];



static int premia_interactive_menu(  Planning          *pt_plan,
                                     Model             **models,
                                     Family            **families,
                                     Pricing           **pricings,
                                     int               user)
{

  Model*            pt_model;
  Option*           pt_option;
  Pricing*          pt_pricing;
  PricingMethod*    pt_method;
  DynamicTest*      pt_test;
  PricingMethod*    pt_methods_available[MAX_METHODS];

  if (OutputFile(&out_stream)!=OK) return WRONG;
  if (SelectModel(user,pt_plan,models,families,pricings,&pt_model)!=OK) return WRONG;
  if (SelectOption(user,pt_plan,families,pt_model,pricings,&pt_option)!=OK) return WRONG;
  if (SelectPricing(user,pt_model,pt_option,pricings,&pt_pricing)!=OK) return WRONG;

  while(1){
    if (SelectMethod(user,pt_plan,pt_pricing,pt_option,pt_model,&pt_method)!=OK)
      return FAIL;
    if (SelectTest(user,pt_plan,pt_pricing,pt_option,pt_model,pt_method,&pt_test)!=OK)
      return FAIL;
    if (GetTimeInfo(user,pt_plan,&computation_time_info)!=OK)
      return FAIL;
          
    if ((pt_plan->Action=='p')||
        (( pt_plan->Action=='t')&&
         (GetTest(user,pt_plan,pt_pricing,pt_option,pt_test)==OK))){
      (void)ShowPlanning(NAMEONLYTOFILE,pt_plan);
      (void)Action(pt_model,pt_option,pt_pricing,
                   pt_method,pt_test,NAMEONLYTOFILE,pt_plan,&computation_time_info);
      Fprintf(TOSCREEN,"\nComputing...\n");
      Iterate(pt_plan,&(pt_plan->Par[0]),0,pt_plan->Action,pt_model,pt_option,
              pt_pricing,pt_method,pt_test,TOFILE,&computation_time_info);
      pt_methods_available[pt_plan->NumberOfMethods]=pt_method;
      if (pt_plan->Action=='t' || MoreAction(&(pt_plan->NumberOfMethods))==FAIL)
        break;
      else
        free_premia_method(pt_method);
    }
  }

  fclose(out_stream);

  if ((pt_plan->Action=='p') && (pt_plan->VarNumber>0))
    (void)BuildGnuStuff(pt_plan,pt_model,pt_option,pt_pricing,pt_methods_available);

  if (pt_plan->Action=='t')
    {
      (void)FreeTest(pt_test); 
      (void)BuildGnuStuffTest(pt_model,pt_option,pt_pricing,pt_method,pt_test); 
    }
  free_premia_model(pt_model);
  free_premia_option(pt_option);
  free_premia_method(pt_method);
  return OK;
}


static int premia_treat_input_file ( Planning          *pt_plan,
                                     Model             **models,
                                     Family            **families,
                                     Pricing           **pricings,
                                     int               user,
                                     char              FileRead[MAX_LINE][MAX_CHAR_LINE])
{

  Model*            pt_model;
  Option*           pt_option;
  Pricing*          pt_pricing;
  PricingMethod*    pt_method;
  DynamicTest*      pt_test;
  PricingMethod*    pt_methods_available[MAX_METHODS];
  
  pt_plan->Action=FChooseAction(FileRead);
  if (OutputFile(&out_stream)!=OK) return FAIL;
  if (FSelectModel(FileRead,user,pt_plan,models,&pt_model)!=OK) return FAIL;
  if (FSelectOption(FileRead,user,pt_plan,families,pt_model,pricings,&pt_option)!=OK) return FAIL;
  if (FSelectPricing(FileRead,user,pt_model,pt_option,pricings,&pt_pricing)!=OK) return FAIL;
  
  while(1){
    if (FSelectMethod(FileRead,user,pt_plan,pt_pricing,pt_option,pt_model,&pt_method)!=OK) return FAIL;
    if (FSelectTest(FileRead,user,pt_plan,pt_pricing,pt_option,pt_model,pt_method,&pt_test)!=OK)
      return FAIL;
    if (FGetTimeInfo(FileRead,user,pt_plan,&computation_time_info)!=OK) return FAIL;
            
    if ((pt_plan->Action=='p')||(( pt_plan->Action=='t'))){
      (void)ShowPlanning(NAMEONLYTOFILE,pt_plan);
      (void)Action(pt_model,pt_option,pt_pricing,
                   pt_method,pt_test,NAMEONLYTOFILE,pt_plan,&computation_time_info);
      
      Fprintf(TOSCREEN,"\nComputing...\n");
      Iterate(pt_plan,&(pt_plan->Par[0]),0,pt_plan->Action,pt_model,
              pt_option,pt_pricing,pt_method,pt_test,TOFILE,&computation_time_info);
      pt_methods_available[pt_plan->NumberOfMethods]=pt_method;
      if (pt_plan->Action=='t' || FMoreAction(FileRead,&(pt_plan->NumberOfMethods))==FAIL)
        break;
      else
        free_premia_method(pt_method);
    }
  } 
  fclose(out_stream);
  if ((pt_plan->Action=='p') && (pt_plan->VarNumber>0))
    (void)BuildGnuStuff(pt_plan,pt_model,pt_option,pt_pricing,pt_methods_available); 

  if (pt_plan->Action=='t')
    {
      (void)FreeTest(pt_test); 
      (void)BuildGnuStuffTest(pt_model,pt_option,pt_pricing,pt_method,pt_test); 
    }
  
  free_premia_model(pt_model);
  free_premia_option(pt_option);
  free_premia_method(pt_method);
  return OK;
}

int main(int argc, char *argv[])
{
  Planning         *pt_plan=NULL;
  int              user, i;
  char             FileRead[MAX_LINE][MAX_CHAR_LINE];
  char             product;
  char             exec_dir[MAX_PATH_LEN]="";
  
  _interactive_call=1;
  
  if ((pt_plan = malloc (sizeof(Planning)))==NULL)
    {
      printf("Memory Allocation error\n"); exit(1);
    }

  get_exec_directory(exec_dir,argv[0]);
  premia_self_set_global_vars(exec_dir);

  if( argc==2 ){
    ReadInputFile(argv[1],FileRead);
    InputMode(&user);
    WellcomeMsg(user);
    if ((InitErrorMsg()==OK)&&(InitVar()==OK))
      {
        ResetPlanning(pt_plan);
        product=FChooseProduct(FileRead);
        i = 0;
        while (premia_assets[i].name != NULL)
          {
            if (product == premia_assets[i].label)
              {
                premia_treat_input_file(pt_plan, premia_assets[i].models,  premia_assets[i].families, 
                                        premia_assets[i].pricings,user, FileRead);
                break;
              }
            i++;
          }
        if (premia_assets[i].name == NULL) return FAIL;
      }
                        
  }else{
    InputMode(&user);
    WellcomeMsg(user);
    if ((InitErrorMsg()==OK)&&(InitVar()==OK)){
      do{
        ResetPlanning(pt_plan);
        product = ChooseProduct();
        i = 0;
        while (premia_assets[i].name != NULL)
          {
            if (product == premia_assets[i].label)
              {
                pt_plan->Action = ChooseAction (product);
                premia_interactive_menu(pt_plan, premia_assets[i].models,  premia_assets[i].families, 
                                        premia_assets[i].pricings,user);
                break;
              }
            i++;
          }
      }while (NextSession(pt_plan,pt_plan->Action,user)==OK);
    }
  }
  
  (void)ExitVar();
  free(pt_plan);
  return OK;
}
