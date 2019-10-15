#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _Ares_Ares_mix_Combined_Independent(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Ares_mix_Combined_LogLikelihood(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Ares_mix_Combined_Mstep(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Ares_mix_Combined_Mstep_Mu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Ares_mix_Combined_Mstep2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Ares_mix_Joint(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Ares_mix_SNP_Gradient(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Ares_mix_SNP_LogLikelihood(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Ares_mix_SNP_Mstep(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Min_Alpha_Ares_mix_SNP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Min_Mu_Ares_mix_Combined(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Min_Mu_Ares_mix_SNP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Min_Mu_Binomial(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Min_Mu_Binomial_repara(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Min_Omega_Ares_mix_SNP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Min_Omega_Binomial(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Min_Omega_Binomial_repara(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Mixed_Binomial_LogLikelihood(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Mixed_Binomial_LogLikelihood_newton(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Mixed_Binomial_LogLikelihood_repara(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Solve_Ares_mix_Combined(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Solve_Ares_mix_Combined_Independent(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Solve_Ares_mix_Combined_Mu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Solve_Ares_mix_SNP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Solve_mixed_Binomial_Goldensection(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Solve_mixed_Binomial_newton_Goldensection(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _Ares_Solve_mixed_Binomial_repara(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_Ares_Ares_mix_Combined_Independent",             (DL_FUNC) &_Ares_Ares_mix_Combined_Independent,             13},
    {"_Ares_Ares_mix_Combined_LogLikelihood",           (DL_FUNC) &_Ares_Ares_mix_Combined_LogLikelihood,           16},
    {"_Ares_Ares_mix_Combined_Mstep",                   (DL_FUNC) &_Ares_Ares_mix_Combined_Mstep,                   13},
    {"_Ares_Ares_mix_Combined_Mstep_Mu",                (DL_FUNC) &_Ares_Ares_mix_Combined_Mstep_Mu,                13},
    {"_Ares_Ares_mix_Combined_Mstep2",                  (DL_FUNC) &_Ares_Ares_mix_Combined_Mstep2,                  13},
    {"_Ares_Ares_mix_Joint",                            (DL_FUNC) &_Ares_Ares_mix_Joint,                            13},
    {"_Ares_Ares_mix_SNP_Gradient",                     (DL_FUNC) &_Ares_Ares_mix_SNP_Gradient,                      7},
    {"_Ares_Ares_mix_SNP_LogLikelihood",                (DL_FUNC) &_Ares_Ares_mix_SNP_LogLikelihood,                 9},
    {"_Ares_Ares_mix_SNP_Mstep",                        (DL_FUNC) &_Ares_Ares_mix_SNP_Mstep,                         7},
    {"_Ares_Min_Alpha_Ares_mix_SNP",                    (DL_FUNC) &_Ares_Min_Alpha_Ares_mix_SNP,                    11},
    {"_Ares_Min_Mu_Ares_mix_Combined",                  (DL_FUNC) &_Ares_Min_Mu_Ares_mix_Combined,                  18},
    {"_Ares_Min_Mu_Ares_mix_SNP",                       (DL_FUNC) &_Ares_Min_Mu_Ares_mix_SNP,                       11},
    {"_Ares_Min_Mu_Binomial",                           (DL_FUNC) &_Ares_Min_Mu_Binomial,                           10},
    {"_Ares_Min_Mu_Binomial_repara",                    (DL_FUNC) &_Ares_Min_Mu_Binomial_repara,                    10},
    {"_Ares_Min_Omega_Ares_mix_SNP",                    (DL_FUNC) &_Ares_Min_Omega_Ares_mix_SNP,                    11},
    {"_Ares_Min_Omega_Binomial",                        (DL_FUNC) &_Ares_Min_Omega_Binomial,                        10},
    {"_Ares_Min_Omega_Binomial_repara",                 (DL_FUNC) &_Ares_Min_Omega_Binomial_repara,                 10},
    {"_Ares_Mixed_Binomial_LogLikelihood",              (DL_FUNC) &_Ares_Mixed_Binomial_LogLikelihood,               8},
    {"_Ares_Mixed_Binomial_LogLikelihood_newton",       (DL_FUNC) &_Ares_Mixed_Binomial_LogLikelihood_newton,        8},
    {"_Ares_Mixed_Binomial_LogLikelihood_repara",       (DL_FUNC) &_Ares_Mixed_Binomial_LogLikelihood_repara,        8},
    {"_Ares_Solve_Ares_mix_Combined",                   (DL_FUNC) &_Ares_Solve_Ares_mix_Combined,                   12},
    {"_Ares_Solve_Ares_mix_Combined_Independent",       (DL_FUNC) &_Ares_Solve_Ares_mix_Combined_Independent,       12},
    {"_Ares_Solve_Ares_mix_Combined_Mu",                (DL_FUNC) &_Ares_Solve_Ares_mix_Combined_Mu,                12},
    {"_Ares_Solve_Ares_mix_SNP",                        (DL_FUNC) &_Ares_Solve_Ares_mix_SNP,                         6},
    {"_Ares_Solve_mixed_Binomial_Goldensection",        (DL_FUNC) &_Ares_Solve_mixed_Binomial_Goldensection,         6},
    {"_Ares_Solve_mixed_Binomial_newton_Goldensection", (DL_FUNC) &_Ares_Solve_mixed_Binomial_newton_Goldensection,  6},
    {"_Ares_Solve_mixed_Binomial_repara",               (DL_FUNC) &_Ares_Solve_mixed_Binomial_repara,                6},
    {NULL, NULL, 0}
};

void R_init_Ares(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
