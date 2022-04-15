// Automatically generated, editing not advised.
#ifndef R_GLMNET_H
#define R_GLMNET_H
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("glmnet", String)
#else
#define _(String) (String)
#endif

extern "C" {
SEXP _glmnet_elnet_exp(SEXP kaSEXP, SEXP parmSEXP, SEXP xSEXP, SEXP ySEXP, SEXP wSEXP, SEXP jdSEXP, SEXP vpSEXP, SEXP clSEXP, SEXP neSEXP, SEXP nxSEXP, SEXP nlamSEXP, SEXP flminSEXP, SEXP ulamSEXP, SEXP thrSEXP, SEXP isdSEXP, SEXP intrSEXP, SEXP maxitSEXP, SEXP pbSEXP, SEXP lmuSEXP, SEXP a0SEXP, SEXP caSEXP, SEXP iaSEXP, SEXP ninSEXP, SEXP rsqSEXP, SEXP almSEXP, SEXP nlpSEXP, SEXP jerrSEXP);
SEXP _glmnet_spelnet_exp(SEXP kaSEXP, SEXP parmSEXP, SEXP xSEXP, SEXP ySEXP, SEXP wSEXP, SEXP jdSEXP, SEXP vpSEXP, SEXP clSEXP, SEXP neSEXP, SEXP nxSEXP, SEXP nlamSEXP, SEXP flminSEXP, SEXP ulamSEXP, SEXP thrSEXP, SEXP isdSEXP, SEXP intrSEXP, SEXP maxitSEXP, SEXP pbSEXP, SEXP lmuSEXP, SEXP a0SEXP, SEXP caSEXP, SEXP iaSEXP, SEXP ninSEXP, SEXP rsqSEXP, SEXP almSEXP, SEXP nlpSEXP, SEXP jerrSEXP);
SEXP _glmnet_lognet_exp(SEXP parmSEXP, SEXP xSEXP, SEXP ySEXP, SEXP gSEXP, SEXP jdSEXP, SEXP vpSEXP, SEXP clSEXP, SEXP neSEXP, SEXP nxSEXP, SEXP nlamSEXP, SEXP flminSEXP, SEXP ulamSEXP, SEXP thrSEXP, SEXP isdSEXP, SEXP intrSEXP, SEXP maxitSEXP, SEXP koptSEXP, SEXP pbSEXP, SEXP lmuSEXP, SEXP a0SEXP, SEXP caSEXP, SEXP iaSEXP, SEXP ninSEXP, SEXP nulldevSEXP, SEXP devSEXP, SEXP almSEXP, SEXP nlpSEXP, SEXP jerrSEXP);
SEXP _glmnet_splognet_exp(SEXP parmSEXP, SEXP xSEXP, SEXP ySEXP, SEXP gSEXP, SEXP jdSEXP, SEXP vpSEXP, SEXP clSEXP, SEXP neSEXP, SEXP nxSEXP, SEXP nlamSEXP, SEXP flminSEXP, SEXP ulamSEXP, SEXP thrSEXP, SEXP isdSEXP, SEXP intrSEXP, SEXP maxitSEXP, SEXP koptSEXP, SEXP pbSEXP, SEXP lmuSEXP, SEXP a0SEXP, SEXP caSEXP, SEXP iaSEXP, SEXP ninSEXP, SEXP nulldevSEXP, SEXP devSEXP, SEXP almSEXP, SEXP nlpSEXP, SEXP jerrSEXP);
SEXP _glmnet_fishnet_exp(SEXP parmSEXP, SEXP xSEXP, SEXP ySEXP, SEXP gSEXP, SEXP wSEXP, SEXP jdSEXP, SEXP vpSEXP, SEXP clSEXP, SEXP neSEXP, SEXP nxSEXP, SEXP nlamSEXP, SEXP flminSEXP, SEXP ulamSEXP, SEXP thrSEXP, SEXP isdSEXP, SEXP intrSEXP, SEXP maxitSEXP, SEXP pbSEXP, SEXP lmuSEXP, SEXP a0SEXP, SEXP caSEXP, SEXP iaSEXP, SEXP ninSEXP, SEXP nulldevSEXP, SEXP devSEXP, SEXP almSEXP, SEXP nlpSEXP, SEXP jerrSEXP);
SEXP _glmnet_spfishnet_exp(SEXP parmSEXP, SEXP xSEXP, SEXP ySEXP, SEXP gSEXP, SEXP wSEXP, SEXP jdSEXP, SEXP vpSEXP, SEXP clSEXP, SEXP neSEXP, SEXP nxSEXP, SEXP nlamSEXP, SEXP flminSEXP, SEXP ulamSEXP, SEXP thrSEXP, SEXP isdSEXP, SEXP intrSEXP, SEXP maxitSEXP, SEXP pbSEXP, SEXP lmuSEXP, SEXP a0SEXP, SEXP caSEXP, SEXP iaSEXP, SEXP ninSEXP, SEXP nulldevSEXP, SEXP devSEXP, SEXP almSEXP, SEXP nlpSEXP, SEXP jerrSEXP);
SEXP _glmnet_multelnet_exp(SEXP parmSEXP, SEXP xSEXP, SEXP ySEXP, SEXP wSEXP, SEXP jdSEXP, SEXP vpSEXP, SEXP clSEXP, SEXP neSEXP, SEXP nxSEXP, SEXP nlamSEXP, SEXP flminSEXP, SEXP ulamSEXP, SEXP thrSEXP, SEXP isdSEXP, SEXP jsdSEXP, SEXP intrSEXP, SEXP maxitSEXP, SEXP pbSEXP, SEXP lmuSEXP, SEXP a0SEXP, SEXP caSEXP, SEXP iaSEXP, SEXP ninSEXP, SEXP rsqSEXP, SEXP almSEXP, SEXP nlpSEXP, SEXP jerrSEXP);
SEXP _glmnet_multspelnet_exp(SEXP parmSEXP, SEXP xSEXP, SEXP ySEXP, SEXP wSEXP, SEXP jdSEXP, SEXP vpSEXP, SEXP clSEXP, SEXP neSEXP, SEXP nxSEXP, SEXP nlamSEXP, SEXP flminSEXP, SEXP ulamSEXP, SEXP thrSEXP, SEXP isdSEXP, SEXP jsdSEXP, SEXP intrSEXP, SEXP maxitSEXP, SEXP pbSEXP, SEXP lmuSEXP, SEXP a0SEXP, SEXP caSEXP, SEXP iaSEXP, SEXP ninSEXP, SEXP rsqSEXP, SEXP almSEXP, SEXP nlpSEXP, SEXP jerrSEXP);
SEXP _glmnet_get_int_parms(SEXP smlSEXP, SEXP epsSEXP, SEXP bigSEXP, SEXP mnlamSEXP, SEXP rsqmaxSEXP, SEXP pminSEXP, SEXP exmxSEXP, SEXP itraceSEXP);
SEXP _glmnet_get_int_parms2(SEXP epsnrSEXP, SEXP mxitnrSEXP);
SEXP _glmnet_chg_fract_dev(SEXP argSEXP);
SEXP _glmnet_chg_dev_max(SEXP argSEXP);
SEXP _glmnet_chg_min_flmin(SEXP argSEXP);
SEXP _glmnet_chg_big(SEXP argSEXP);
SEXP _glmnet_chg_min_lambdas(SEXP irgSEXP);
SEXP _glmnet_chg_min_null_prob(SEXP argSEXP);
SEXP _glmnet_chg_max_exp(SEXP argSEXP);
SEXP _glmnet_chg_itrace(SEXP irgSEXP);
SEXP _glmnet_chg_bnorm(SEXP argSEXP, SEXP irgSEXP);
SEXP _glmnet_get_bnorm(SEXP argSEXP, SEXP irgSEXP);
SEXP _glmnet_chg_epsnr(SEXP argSEXP);
SEXP _glmnet_chg_mxitnr(SEXP irgSEXP);
SEXP _glmnet_wls_exp(SEXP alm0SEXP, SEXP almcSEXP, SEXP alphaSEXP, SEXP mSEXP, SEXP noSEXP, SEXP niSEXP, SEXP xSEXP, SEXP rSEXP, SEXP xvSEXP, SEXP vSEXP, SEXP intrSEXP, SEXP juSEXP, SEXP vpSEXP, SEXP clSEXP, SEXP nxSEXP, SEXP thrSEXP, SEXP maxitSEXP, SEXP aSEXP, SEXP aintSEXP, SEXP gSEXP, SEXP iaSEXP, SEXP iySEXP, SEXP izSEXP, SEXP mmSEXP, SEXP ninoSEXP, SEXP rsqcSEXP, SEXP nlpSEXP, SEXP jerrSEXP);
SEXP _glmnet_spwls_exp(SEXP alm0SEXP, SEXP almcSEXP, SEXP alphaSEXP, SEXP mSEXP, SEXP noSEXP, SEXP niSEXP, SEXP xSEXP, SEXP xmSEXP, SEXP xsSEXP, SEXP rSEXP, SEXP xvSEXP, SEXP vSEXP, SEXP intrSEXP, SEXP juSEXP, SEXP vpSEXP, SEXP clSEXP, SEXP nxSEXP, SEXP thrSEXP, SEXP maxitSEXP, SEXP aSEXP, SEXP aintSEXP, SEXP gSEXP, SEXP iaSEXP, SEXP iySEXP, SEXP izSEXP, SEXP mmSEXP, SEXP ninoSEXP, SEXP rsqcSEXP, SEXP nlpSEXP, SEXP jerrSEXP);
SEXP storePB(SEXP);

void F77_SUB(coxnet)(
		     double *parm,
		     int *no,
		     int *ni,
		     double *x,
		     double *y,
		     double *d,
		     double *o,
		     double *w,
		     int *jd,
		     double *vp,
		     double *cl,
		     int *ne,
		     int *nx,
		     int *nlam,
		     double *flmin,
		     double *ulam,
		     double *thr,
		     int *maxit,
		     int *isd,
		     int *lmu,
		     double *ca,
		     int *ia,
		     int *nin,
		     double *dev0,
		     double *fdev,
		     double *alm,
		     int *nlp,
		     int *jerrc
		     );

void F77_SUB(loglike)(
		      int *no,
		      int *ni,
		      double *x,
		      double *y,
		      double *d,
		      double *g,
		      double *w,
		      int *nlam,
		      double *a,
		      double *flog,
		      int *jerr
		      );
} // end extern "C"

static const R_CallMethodDef CallEntries[] = {
    {"_glmnet_elnet_exp", (DL_FUNC) &_glmnet_elnet_exp, 27},
    {"_glmnet_spelnet_exp", (DL_FUNC) &_glmnet_spelnet_exp, 27},
    {"_glmnet_lognet_exp", (DL_FUNC) &_glmnet_lognet_exp, 28},
    {"_glmnet_splognet_exp", (DL_FUNC) &_glmnet_splognet_exp, 28},
    {"_glmnet_fishnet_exp", (DL_FUNC) &_glmnet_fishnet_exp, 28},
    {"_glmnet_spfishnet_exp", (DL_FUNC) &_glmnet_spfishnet_exp, 28},
    {"_glmnet_multelnet_exp", (DL_FUNC) &_glmnet_multelnet_exp, 27},
    {"_glmnet_multspelnet_exp", (DL_FUNC) &_glmnet_multspelnet_exp, 27},
    {"_glmnet_get_int_parms", (DL_FUNC) &_glmnet_get_int_parms, 8},
    {"_glmnet_get_int_parms2", (DL_FUNC) &_glmnet_get_int_parms2, 2},
    {"_glmnet_chg_fract_dev", (DL_FUNC) &_glmnet_chg_fract_dev, 1},
    {"_glmnet_chg_dev_max", (DL_FUNC) &_glmnet_chg_dev_max, 1},
    {"_glmnet_chg_min_flmin", (DL_FUNC) &_glmnet_chg_min_flmin, 1},
    {"_glmnet_chg_big", (DL_FUNC) &_glmnet_chg_big, 1},
    {"_glmnet_chg_min_lambdas", (DL_FUNC) &_glmnet_chg_min_lambdas, 1},
    {"_glmnet_chg_min_null_prob", (DL_FUNC) &_glmnet_chg_min_null_prob, 1},
    {"_glmnet_chg_max_exp", (DL_FUNC) &_glmnet_chg_max_exp, 1},
    {"_glmnet_chg_itrace", (DL_FUNC) &_glmnet_chg_itrace, 1},
    {"_glmnet_chg_bnorm", (DL_FUNC) &_glmnet_chg_bnorm, 2},
    {"_glmnet_get_bnorm", (DL_FUNC) &_glmnet_get_bnorm, 2},
    {"_glmnet_chg_epsnr", (DL_FUNC) &_glmnet_chg_epsnr, 1},
    {"_glmnet_chg_mxitnr", (DL_FUNC) &_glmnet_chg_mxitnr, 1},
    {"_glmnet_wls_exp", (DL_FUNC) &_glmnet_wls_exp, 28},
    {"_glmnet_spwls_exp", (DL_FUNC) &_glmnet_spwls_exp, 30},
    {"storePB", (DL_FUNC) &storePB, 1},
    {NULL, NULL, 0}
};

#define FDEF(name)  {#name, (DL_FUNC) &F77_SUB(name), sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}

static R_NativePrimitiveArgType coxnet_t[] = {
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(spelnet)(
		      int *ka,
		      double *parm,
		      int *no,
		      int *ni,
		      double *x,
		      int *ix,
		      int *jx,
		      double *y,
		      double *w,
		      int *jd,
		      double *vp,
		      double *cl,
		      int *ne,
		      int *nx,
		      int *nlam,
		      double *flmin,
		      double *ulam,
		      double *thr,
		      int *isd,
		      int *intr,
		      int *maxit,
		      int *lmu,
		      double *a0,
		      double *ca,
		      int *ia,
		      int *nin,
		      double *rsq,
		      double *alm,
		      int *nlp,
		      int *jerr
		      );

static R_NativePrimitiveArgType loglike_t[] = {
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP
};

static R_FortranMethodDef fMethods[] = {
  FDEF(coxnet) ,
  FDEF(loglike) ,
  {NULL, NULL, 0}
};

void R_init_glmnet(DllInfo *dll){
  R_registerRoutines(dll, NULL, CallEntries, fMethods, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

#endif
