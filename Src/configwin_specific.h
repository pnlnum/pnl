/**************************************************/
/* Windows specific macros (not common with Unix) */
/**************************************************/

//The dllexport and dllimport storage-class attributes are Microsoft-specific 
//extensions to the C and C++ languages.
//They enable you to export and import functions, data, and objects to and from a DLL.
#ifdef PREMIA_EXPORTS
#define PREMIA_API __declspec(dllexport)
#else
#define PREMIA_API __declspec(dllimport)
#endif

#pragma warning ( disable : 4217 4049 )

