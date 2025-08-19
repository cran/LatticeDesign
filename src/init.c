#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void ProjSep(void *, void *, void *, void *);
extern void lrs_XuHe (void*, void*, void*, void*, void*, void*, void*);
extern void selectlss(void*, void*, void*, void*, void*, void*, void*, void*, void*, void*);
extern void WSL2_c(void*, void*, void*);
extern void WFL2_c(void*, void*, void*);
extern void WDL2_c(void*, void*, void*);
extern void WSL_c(void*, void*, void*, void*);
extern void WPL_c(void*, void*, void*, void*);
extern void WDL_c(void*, void*, void*, void*);
extern void bWSL_c(void*, void*, void*, void*);
extern void bWFL_c(void*, void*, void*, void*);
extern void bWDL_c(void*, void*, void*, void*);

static const R_CMethodDef CEntries[] = {
    {"ProjSep", (DL_FUNC) &ProjSep, 4},
    {"lrs_XuHe", (DL_FUNC) &lrs_XuHe, 7},
    {"selectlss", (DL_FUNC) &selectlss, 10},
    {"WSL2_c", (DL_FUNC) &WSL2_c, 3},
    {"WFL2_c", (DL_FUNC) &WFL2_c, 3},
    {"WDL2_c", (DL_FUNC) &WDL2_c, 3},
    {"WSL_c", (DL_FUNC) &WSL_c, 4},
    {"WPL_c", (DL_FUNC) &WPL_c, 4},
    {"WDL_c", (DL_FUNC) &WDL_c, 4},
    {"bWSL_c", (DL_FUNC) &bWSL_c, 4},
    {"bWFL_c", (DL_FUNC) &bWFL_c, 4},
    {"bWDL_c", (DL_FUNC) &bWDL_c, 4},
   {NULL, NULL, 0}
};

void R_init_LatticeDesign(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
