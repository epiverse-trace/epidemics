#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  The following symbols/expressions for .NAME have been omitted

    _epidemics_model_default_internal
    _epidemics_model_diphtheria_internal
    _epidemics_model_vacamole_internal

  Most likely possible values need to be added below.
*/

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void seirv_model_initmod_desolve(void *);
extern void seirv_model_rhs_dde(void *);
extern void seirv_model_rhs_desolve(void *);

/* .Call calls */
extern SEXP seirv_model_contents(SEXP);
extern SEXP seirv_model_create(SEXP);
extern SEXP seirv_model_initial_conditions(SEXP, SEXP);
extern SEXP seirv_model_metadata(SEXP);
extern SEXP seirv_model_rhs_r(SEXP, SEXP, SEXP);
extern SEXP seirv_model_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP seirv_model_set_user(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"seirv_model_initmod_desolve", (DL_FUNC) &seirv_model_initmod_desolve, 1},
    {"seirv_model_rhs_dde",         (DL_FUNC) &seirv_model_rhs_dde,         1},
    {"seirv_model_rhs_desolve",     (DL_FUNC) &seirv_model_rhs_desolve,     1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"seirv_model_contents",           (DL_FUNC) &seirv_model_contents,           1},
    {"seirv_model_create",             (DL_FUNC) &seirv_model_create,             1},
    {"seirv_model_initial_conditions", (DL_FUNC) &seirv_model_initial_conditions, 2},
    {"seirv_model_metadata",           (DL_FUNC) &seirv_model_metadata,           1},
    {"seirv_model_rhs_r",              (DL_FUNC) &seirv_model_rhs_r,              3},
    {"seirv_model_set_initial",        (DL_FUNC) &seirv_model_set_initial,        4},
    {"seirv_model_set_user",           (DL_FUNC) &seirv_model_set_user,           2},
    {NULL, NULL, 0}
};

void R_init_epidemics(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
