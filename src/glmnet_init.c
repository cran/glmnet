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

#define FDEF(name)  {#name, (DL_FUNC) &F77_SUB(name), sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}
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

static R_NativePrimitiveArgType spelnet_t[] = {
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
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
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(elnet)(
		    int *ka,
		    double *parm,
		    int *no,
		    int *ni,
		    double *x,
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

static R_NativePrimitiveArgType elnet_t[] = {
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
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
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(fishnet)(
		      double *parm,
		      int *no,
		      int *ni,
		      double *x,
		      double *y,
		      double *g,
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
		      double *dev0,
		      double *dev,
		      double *alm,
		      int *nlp,
		      int *jerr
		      );

static R_NativePrimitiveArgType fishnet_t[] = {
  REALSXP,
  INTSXP,
  INTSXP,
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
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(spfishnet)(
			double *parm,
			int *no,
			int *ni,
			double *x,
			int *ix,
			int *jx,
			double *y,
			double *g,
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
			double *dev0,
			double *dev,
			double *alm,
			int *nlp,
			int *jerr
			);

static R_NativePrimitiveArgType spfishnet_t[] = {
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
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
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(splognet)(
		       double *parm,
		       int *no,
		       int *ni,
		       int *nc,
		       double *x,
		       int *ix,
		       int *jx,
		       double *y,
		       double *g,
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
		       int *kopt,
		       int *lmu,
		       double *a0,
		       double *ca,
		       int *ia,
		       int *nin,
		       double *dev0,
		       double *dev,
		       double *alm,
		       int *nlp,
		       int *jerr
		       );

static R_NativePrimitiveArgType splognet_t[] = {
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
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
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(lognet)(
		     double *parm,
		     int *no,
		     int *ni,
		     int *nc,
		     double *x,
		     double *y,
		     double *g,
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
		     int *kopt,
		     int *lmu,
		     double *a0,
		     double *ca,
		     int *ia,
		     int *nin,
		     double *dev0,
		     double *dev,
		     double *alm,
		     int *nlp,
		     int *jerr
		     );

static R_NativePrimitiveArgType lognet_t[] = {
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
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
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(multspelnet)(
			  double *parm,
			  int *no,
			  int *ni,
			  int *nr,
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
			  int *jsd,
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

static R_NativePrimitiveArgType multspelnet_t[] = {
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  INTSXP,
  INTSXP,
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
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP
};
void F77_SUB(multelnet)(
			double *parm,
			int *no,
			int *ni,
			int *nr,
			double *x,
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
			int *jsd,
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

static R_NativePrimitiveArgType multelnet_t[] = {
  REALSXP,
  INTSXP,
  INTSXP,
  INTSXP,
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
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  INTSXP
};
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
void F77_SUB(get_int_parms)(
			    double *sml,
			    double *eps,
			    double *big,
			    int *mnlam,
			    double *rsqmax,
			    double *pmin,
			    double *exmx
			    );

static R_NativePrimitiveArgType get_int_parms_t[] = {
  REALSXP,
  REALSXP,
  REALSXP,
  INTSXP,
  REALSXP,
  REALSXP,
  REALSXP
};
void F77_SUB(chg_fract_dev)(
			    double *fdev
			    );

static R_NativePrimitiveArgType chg_fract_dev_t[] = {
  REALSXP
};
void F77_SUB(chg_dev_max)(
			  double *devmax
			  );

static R_NativePrimitiveArgType chg_dev_max_t[] = {
  REALSXP
};
void F77_SUB(chg_min_flmin)(
			    double *eps
			    );

static R_NativePrimitiveArgType chg_min_flmin_t[] = {
  REALSXP
};
void F77_SUB(chg_big)(
		      double *big
		      );

static R_NativePrimitiveArgType chg_big_t[] = {
  REALSXP
};
void F77_SUB(chg_min_lambdas)(
			      int *mnlam
			      );

static R_NativePrimitiveArgType chg_min_lambdas_t[] = {
  INTSXP
};
void F77_SUB(chg_min_null_prob)(
				double *pmin
				);

static R_NativePrimitiveArgType chg_min_null_prob_t[] = {
  REALSXP
};
void F77_SUB(chg_max_exp)(
			  double *exmx
			  );

static R_NativePrimitiveArgType chg_max_exp_t[] = {
  REALSXP
};
void F77_SUB(chg_bnorm)(
			double *prec,
			int *mxit
			);

static R_NativePrimitiveArgType chg_bnorm_t[] = {
  REALSXP,
  INTSXP
};
void F77_SUB(get_bnorm)(
			double *arg,
			int *irg
			);

static R_NativePrimitiveArgType get_bnorm_t[] = {
  REALSXP,
  INTSXP
};

static R_FortranMethodDef fMethods[] = {
  FDEF(coxnet) ,
  FDEF(spelnet) ,
  FDEF(elnet) ,
  FDEF(fishnet) ,
  FDEF(spfishnet) ,
  FDEF(splognet) ,
  FDEF(lognet) ,
  FDEF(multspelnet) ,
  FDEF(multelnet) ,
  FDEF(loglike) ,
  FDEF(get_int_parms) ,
  FDEF(chg_fract_dev) ,
  FDEF(chg_dev_max) ,
  FDEF(chg_min_flmin) ,
  FDEF(chg_big) ,
  FDEF(chg_min_lambdas) ,
  FDEF(chg_min_null_prob) ,
  FDEF(chg_max_exp) ,
  FDEF(chg_bnorm) ,
  FDEF(get_bnorm) ,
  {NULL, NULL, 0}
};

void R_init_glmnet(DllInfo *dll){
  R_registerRoutines(dll, NULL, NULL, fMethods, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
#endif
