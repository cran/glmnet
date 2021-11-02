#pragma once

#ifndef GLMNETPP_MOCK_LEGACY

extern "C" {

void wls_(
double *alm0,
double *almc,
double *alpha,
int *m,
int *no,
int *ni,
double *x,
double *r,
double *v,
int *intr,
int *ju,
double *vp,
double *cl,
int *nx,
double *thr,
int *maxit,
double *a,
double *aint,
double *g,
int *ia,
int *iy,
int *iz,
int *mm,
int *nino,
double *rsqc,
int *nlp,
int *jerr
);

void elnet_(
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

void spelnet_(
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

void elnet1_(
double *beta,
int *ni,
int *ju,
double *vp,
double *cl,
double *g,
int *no,
int *ne,
int *nx,
double *x,
int *nlam,
double *flmin,
double *ulam,
double *thr,
int *maxit,
double *xv,  
int *lmu,
double *ao,
int *ia,
int *kin,
double *rsqo,
double *almo,
int *nlp,
int *jerr
);

void elnet2_(
double *beta,
int *ni,
int *ju,
double *vp,
double *cl,
double *y,
int *no,
int *ne,
int *nx,
double *x,
int *nlam,
double *flmin,
double *ulam,
double *thr,
int *maxit,
double *xv,  
int *lmu,
double *ao,
int *ia,
int *kin,
double *rsqo,
double *almo,
int *nlp,
int *jerr
);

void spelnet1_(
double *beta,
int *ni,
double *g,
int *no,
double *w,
int *ne,
int *nx,
double *x,
int *ix,
int *jx,
int *ju,
double *vp,
double *cl,
int *nlam,
double *flmin,
double *ulam,
double *thr,
int *maxit,
double *xm,
double *xs,
double *xv,  
int *lmu,
double *ao,
int *ia,
int *kin,
double *rsqo,
double *almo,
int *nlp,
int *jerr
);

void spelnet2_(
double *beta,
int *ni,
double *gy,
double *w,
int *no,
int *ne,
int *nx,
double *x,
int *ix,
int *jx,
int *ju,
double *vp,
double *cl,
int *nlam,
double *flmin,
double *ulam,
double *thr,
int *maxit,
double *xm,
double *xs,
double *xv,  
int *lmu,
double *ao,
int *ia,
int *kin,
double *rsqo,
double *almo,
int *nlp,
int *jerr
);

void standard_(
int *no,
int *ni,
double *x,
double *y,
double *w,
int *isd,
int *intr,
int *ju,
double *g,
double *xm,
double *xs,
double *ym,
double *ys,
double *xv,
int *jerr
);

void standard1_(
int *no,
int *ni,
double *x,
double *y,
double *w,
int *isd,
int *intr,
int *ju,
double *xm,
double *xs,
double *ym,
double *ys,
double *xv,
int *jerr
);

void spstandard_(
int *no,
int *ni,
double *x,
int *ix,
int *jx,
double *y,
double *w,
int *ju,
int *isd,
int *intr,
double *g,
double *xm,
double *xs,
double *ym,
double *ys,
double *xv,
int *jerr
);

void spstandard1_(
int *no,
int *ni,
double *x,
int *ix,
int *jx,
double *y,
double *w,
int *ju,
int *isd,
int *intr,
double *xm,
double *xs,
double *ym,
double *ys,
double *xv,
int *jerr
);

void chkvars_(
int *no,
int *ni,
double *x,
int *ju
);

void spchkvars_(
int *no,
int *ni,
double *x,
int *ix,
int *ju
);

void get_int_parms_(
double *sml,
double *eps,
double *big,
int *mnlam,
double *rsqmax,
double *pmin,
double *exmx,
int *itrace
);

void chg_fract_dev_(
double *fdev
);

void chg_dev_max_(
double *devmax
);

void chg_min_flmin_(
double *eps
);

void chg_big_(
double *big
);

void chg_min_lambdas_(
int *mnlam
);

void chg_min_null_prob_(
double *pmin
);

void chg_max_exp_(
double *exmx
);

void chg_itrace_(
int *itrace
);

void chg_bnorm_(
double *prec,
int *mxit
);

void chg_epsnr_(
double *epsnr
);
 
void chg_mxitnr_(
int *mxitnr
);
 
void setpb_(
int *val
);

void get_int_parms_(
double* sml,
double* eps,
double* big,
int* mnlam,
double* rsqmax,
double* pmin,
double* exmx,
int* itrace
);

} // end extern "C"

#else

inline void setpb_(int *m) {}
inline void elnet1_(
    double* beta,
    int* ni,
    int* ju,
    double* vp,
    double* cl,
    double* g,
    int* no,
    int* ne,
    int* nx,
    double* x,
    int* nlam,
    double* flmin,
    double* ulam,
    double* thr,
    int* maxit,
    double* xv,
    int* lmu,
    double* ao,
    int* ia,
    int* kin,
    double* rsqo,
    double* almo,
    int* nlp,
    int* jerr
) {}

inline void elnet2_(
    double* beta,
    int* ni,
    int* ju,
    double* vp,
    double* cl,
    double* y,
    int* no,
    int* ne,
    int* nx,
    double* x,
    int* nlam,
    double* flmin,
    double* ulam,
    double* thr,
    int* maxit,
    double* xv,
    int* lmu,
    double* ao,
    int* ia,
    int* kin,
    double* rsqo,
    double* almo,
    int* nlp,
    int* jerr
) {}

#endif
