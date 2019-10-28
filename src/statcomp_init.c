#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void HVG_C(void *, void *, void *, void *);
extern void HVG_penetrable_C(void *, void *, void *, void *, void *);
extern void ordinal_pattern_loop(void *, void *, void *, void *, void *, void *);
extern void weighted_ordinal_pattern_loop(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"HVG_C",                         (DL_FUNC) &HVG_C,                         4},
    {"HVG_penetrable_C",              (DL_FUNC) &HVG_penetrable_C,              5},
    {"ordinal_pattern_loop",          (DL_FUNC) &ordinal_pattern_loop,          6},
    {"weighted_ordinal_pattern_loop", (DL_FUNC) &weighted_ordinal_pattern_loop, 5},
    {NULL, NULL, 0}
};

void R_init_statcomp(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
