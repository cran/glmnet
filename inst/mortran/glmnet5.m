"
c
c                          newGLMnet (5/12/14)
c                            
c                        
c                 Elastic net with squared-error loss
c
c dense predictor matrix:
c
c call elnet(ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
c            intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call spelnet(ka,parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
c             isd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse row format
c
c
c other inputs:
c
c   ka = algorithm flag
c      ka=1 => covariance updating algorithm
c      ka=2 => naive algorithm
c   parm = penalty member index (0 <= parm <= 1)
c        = 0.0 => ridge
c        = 1.0 => lasso
c   no = number of observations
c   ni = number of predictor variables
c   y(no) = response vector (overwritten)
c   w(no)= observation weights (overwritten)
c   jd(jd(1)+1) = predictor variable deletion flag
c      jd(1) = 0  => use all variables
c      jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
c   vp(ni) = relative penalties for each predictor variable
c      vp(j) = 0 => jth variable unpenalized
c   cl(2,ni) = interval constraints on coefficient values (overwritten)
c      cl(1,j) = lower bound for jth coefficient value (<= 0.0)
c      cl(2,j) = upper bound for jth coefficient value (>= 0.0)
c   ne = maximum number of variables allowed to enter largest model
c        (stopping criterion)
c   nx = maximum number of variables allowed to enter all models
c        along path (memory allocation, nx > ne).
c   nlam = (maximum) number of lamda values
c   flmin = user control of lamda values (>=0)
c      flmin < 1.0 => minimum lamda = flmin*(largest lamda value)
c      flmin >= 1.0 => use supplied lamda values (see below)
c   ulam(nlam) = user supplied lamda values (ignored if flmin < 1.0)
c   thr = convergence threshold for each lamda solution.
c      iterations stop when the maximum reduction in the criterion value
c      as a result of each parameter update over a single pass
c      is less than thr times the null criterion value.
c      (suggested value, thr=1.0e-5)
c   isd = predictor variable standarization flag:
c      isd = 0 => regression on original predictor variables
c      isd = 1 => regression on standardized predictor variables
c      Note: output solutions always reference original
c            variables locations and scales.
c   intr = intercept flag
c      intr = 0/1 => don't/do include intercept in model
c   maxit = maximum allowed number of passes over the data for all lambda
c      values (suggested values, maxit = 100000)
c
c output:
c
c   lmu = actual number of lamda values (solutions)
c   a0(lmu) = intercept values for each solution
c   ca(nx,lmu) = compressed coefficient values for each solution
c   ia(nx) = pointers to compressed coefficients
c   nin(lmu) = number of compressed coefficients for each solution
c   rsq(lmu) = R**2 values for each solution
c   alm(lmu) = lamda values corresponding to each solution
c   nlp = actual number of passes over the data for all lamda values
c   jerr = error flag:
c      jerr  = 0 => no error
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 10000 => maxval(vp) <= 0.0
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
c
c
c
c least-squares utility routines:
c
c
c uncompress coefficient vectors for all solutions:
c
c call solns(ni,nx,lmu,ca,ia,nin,b)
c
c input:
c
c    ni,nx = input to elnet
c    lmu,ca,ia,nin = output from elnet
c
c output:
c
c    b(ni,lmu) = all elnet returned solutions in uncompressed format
c
c
c uncompress coefficient vector for particular solution:
c
c call uncomp(ni,ca,ia,nin,a)
c
c input:
c
c    ni = total number of predictor variables
c    ca(nx) = compressed coefficient values for the solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for the solution
c
c output:
c
c    a(ni) =  uncompressed coefficient vector
c             referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor matrix:
c
c call modval(a0,ca,ia,nin,n,x,f);
c
c input:
c
c    a0 = intercept
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    n = number of predictor vectors (observations)
c    x(n,ni) = full (uncompressed) predictor matrix
c
c output:
c
c    f(n) = model predictions
c
c
c evaluate linear model from compressed coefficients and
c compressed predictor matrix:
c
c call cmodval(a0,ca,ia,nin,x,ix,jx,n,f);
c
c input:
c
c    a0 = intercept
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    x, ix, jx = predictor matrix in compressed sparse row format
c    n = number of predictor vectors (observations)
c
c output:
c
c    f(n) = model predictions
c
c
c
c
c                           Multiple response
c                  elastic net with squared-error loss
c
c dense predictor matrix:
c
c call multelnet(parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
c                jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call multspelnet(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
c             isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse row format
c
c other inputs:
c
c   nr = number of response variables
c   y(no,nr) = response data matrix (overwritten)
c   jsd = response variable standardization flag
c      jsd = 0 => regression using original response variables
c      jsd = 1 => regression using standardized response variables
c      Note: output solutions always reference original
c            variables locations and scales.
c   all other inputs same as elnet/spelnet above
c
c output:
c
c   a0(nr,lmu) = intercept values for each solution
c   ca(nx,nr,lmu) = compressed coefficient values for each solution
c   all other outputs same as elnet/spelnet above
c   (jerr = 90000 => bounds adjustment non convergence)
c
c
c
c multiple response least-squares utility routines:
c
c
c uncompress coefficient matrix for all solutions:
c
c call multsolns(ni,nx,nr,lmu,ca,ia,nin,b)
c
c input:
c
c    ni,nx,nr = input to multelnet
c    lmu,ca,ia,nin = output from multelnet
c
c output:
c
c    b(ni,nr,lmu) = all multelnet returned solutions in uncompressed format
c
c
c uncompress coefficient matrix for particular solution:
c
c call multuncomp(ni,nr,nx,ca,ia,nin,a)
c
c input:
c
c    ni,nr,nx = input to multelnet
c    ca(nx,nr) = compressed coefficient values for the solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for the solution
c
c output:
c
c    a(ni,nr) =  uncompressed coefficient matrix
c             referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor matrix:
c
c call multmodval(nx,nr,a0,ca,ia,nin,n,x,f);
c
c input:
c
c    nx,nr = input to multelnet
c    a0(nr) = intercepts
c    ca(nx,nr) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    n = number of predictor vectors (observations)
c    x(n,ni) = full (uncompressed) predictor matrix
c
c output:
c
c    f(nr,n) = model predictions
c
c
c evaluate linear model from compressed coefficients and
c compressed predictor matrix:
c
c call multcmodval(nx,nr,a0,ca,ia,nin,x,ix,jx,n,f);
c
c input:
c
c    nx,nr = input to multelnet
c    a0(nr) = intercepts
c    ca(nx,nr) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    x, ix, jx = predictor matrix in compressed sparse row format
c    n = number of predictor vectors (observations)
c
c output:
c
c    f(nr,n) = model predictions
c
c   
c
c
c          Symmetric binomial/multinomial logistic elastic net
c
c
c dense predictor matrix:
c
c call lognet (parm,no,ni,nc,x,y,o,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
c              intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call splognet (parm,no,ni,nc,x,ix,jx,y,o,jd,vp,cl,ne,nx,nlam,flmin,
c      ulam,thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse row format
c
c
c other inputs:
c
c   parm,no,ni,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,maxit
c    = same as elnet above.
c   
c   nc = number of classes (distinct outcome values)
c        nc=1 => binomial two-class logistic regression
c            (all output references class 1)
c   y(no,max(2,nc)) = number of each class at each design point
c      entries may have fractional values or all be zero (overwritten)
c   o(no,nc) = observation off-sets for each class
c   kopt = optimization flag
c      kopt = 0 => Newton-Raphson (recommended)
c      kpot = 1 => modified Newton-Raphson (sometimes faster)
c      kpot = 2 => nonzero coefficients same for each class (nc > 1)
c
c
c output:
c
c   lmu,ia,nin,alm,nlp = same as elent above
c
c   a0(nc,lmu) = intercept values for each class at each solution
c   ca(nx,nc,lmu) = compressed coefficient values for each class at
c                each solution
c   dev0 = null deviance (intercept only model)
c   fdev(lmu) = fraction of devience explained by each solution
c   jerr = error flag
c      jerr = 0  => no error
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 8000 + k => null probability < 1.0e-5 for class k
c         jerr = 9000 + k => null probability for class k
c                            > 1.0 - 1.0e-5
c         jerr = 10000 => maxval(vp) <= 0.0
c         jerr = 90000 => bounds adjustment non convergence
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
c         jerr = -20000-k => max(p*(1-p)) < 1.0e-6 at kth lamda value.
c    o(no,nc) = training data values for last (lmu_th) solution linear
c               combination.
c
c
c
c logistic/multinomial utilitity routines:
c
c
c uncompress coefficient vectors for all solutions:
c
c call lsolns(ni,nx,nc,lmu,ca,ia,nin,b)
c
c input:
c
c    ni,nx,nc = input to lognet
c    lmu,ca,ia,nin = output from lognet
c
c output:
c
c    b(ni,nc,lmu) = all lognet returned solutions in uncompressed format
c
c
c uncompress coefficient vector for particular solution:
c
c call luncomp(ni,nx,nc,ca,ia,nin,a)
c
c input:
c
c    ni, nx, nc = same as above
c    ca(nx,nc) = compressed coefficient values (for each class)
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients
c
c output:
c
c    a(ni,nc) =  uncompressed coefficient vectors
c                 referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor vectors:
c
c call lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans);
c
c input:
c
c    nt = number of observations
c    x(nt,ni) = full (uncompressed) predictor vectors
c    nc, nx = same as above
c    a0(nc) = intercepts
c    ca(nx,nc) = compressed coefficient values (for each class)
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients
c
c output:
c
c ans(nc,nt) = model predictions
c
c
c evaluate linear model from compressed coefficients and
c compressed predictor matrix:
c
c call lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f);
c
c input:
c
c    nc, nx = same as above
c    a0(nc) = intercept
c    ca(nx,nc) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    x, ix, jx = predictor matrix in compressed sparse row format
c    n = number of predictor vectors (observations)
c
c output:
c
c    f(nc,n) = model predictions
c
c
c
c
c                        Poisson elastic net
c
c
c dense predictor matrix:                      
c
c call fishnet (parm,no,ni,x,y,o,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
c               isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c sparse predictor matrix:
c
c call spfishnet (parm,no,ni,x,ix,jx,y,o,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
c               isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c    x, ix, jx = predictor data matrix in compressed sparse row format
c
c other inputs:
c
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   parm,no,ni,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,maxit
c    = same as elnet above
c
c output:
c
c   lmu,a0,ca,ia,nin,alm = same as elnet above
c   dev0,fdev = same as lognet above
c   nlp = total number of passes over predictor variables
c   jerr = error flag
c      jerr = 0  => no error
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 8888 => negative response count y values
c         jerr = 9999 => no positive observations weights
c         jerr = 10000 => maxval(vp) <= 0.0
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
c    o(no) = training data values for last (lmu_th) solution linear
c            combination.
c
c 
c Poisson utility routines:
c
c
c same as elnet above:
c
c    call solns(ni,nx,lmu,ca,ia,nin,b)
c    call uncomp(ni,ca,ia,nin,a)
c    call modval(a0,ca,ia,nin,n,x,f);
c    call cmodval(a0,ca,ia,nin,x,ix,jx,n,f);
c
c compute deviance for given uncompressed data and set of uncompressed
c solutions
c
c call deviance(no,ni,x,y,o,w,nsol,a0,a,flog,jerr)
c
c input:
c
c   no = number of observations
c   ni = number of predictor variables
c   x(no,ni) = predictor data matrix flat file
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   w(no)= observation weights
c   nsol = number of solutions
c   a0(nsol) = intercept for each solution
c   a(ni,nsol) = solution coefficient vectors (uncompressed)
c
c output:
c
c   flog(nsol) = respective deviance values minus null deviance
c   jerr = error flag - see above
c
c
c compute deviance for given compressed data and set of uncompressed solutions
c
c call spdeviance(no,ni,x,ix,jx,y,o,w,nsol,a0,a,flog,jerr)
c
c input:
c
c   no = number of observations
c   ni = number of predictor variables
c   x, ix, jx = predictor data matrix in compressed sparse row format
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   w(no)= observation weights
c   nsol = number of solutions
c   a0(nsol) = intercept for each solution
c   a(ni,nsol) = solution coefficient vectors (uncompressed)
c
c output
c
c   flog(nsol) = respective deviance values minus null deviance
c   jerr = error flag - see above
c
c
c compute deviance for given compressed data and compressed solutions
c
c call cspdeviance(no,x,ix,jx,y,o,w,nx,lmu,a0,ca,ia,nin,flog,jerr)
c
c input:
c
c   no = number of observations
c   x, ix, jx = predictor data matrix in compressed sparse row format
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   w(no)= observation weights
c   nx = input to spfishnet
c   lmu,a0(lmu),ca(nx,lmu),ia(nx),nin(lmu) = output from spfishnet
c
c output
c
c   flog(lmu) = respective deviance values minus null deviance
c   jerr = error flag - see above
c
c
c
c          Elastic net with Cox proportional hazards model
c
c
c dense predictor matrix:
c
c call coxnet (parm,no,ni,x,y,d,o,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
c              maxit,isd,lmu,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c input:
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c   y(no) = observation times
c   d(no) = died/censored indicator
c       d(i)=0.0 => y(i) = censoring time
c       d(i)=1.0 => y(i) = death time
c   o(no) = observation off-sets
c   parm,no,ni,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,maxit
c                = same as fishnet above
c
c output:
c
c   lmu,ca,ia,nin,dev0,fdev,alm,nlp = same as fishnet above
c   jerr = error flag
c      jerr = 0  => no error - output returned
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 8888 => all observations censored (d(i)=0.0)
c         jerr = 9999 => no positive observations weights
c         jerr = 10000 => maxval(vp) <= 0.0
c         jerr = 20000, 30000 => initialization numerical error
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
c         jerr = -30000-k => numerical error at kth lambda value
c    o(no) = training data values for last (lmu_th) solution linear
c            combination.
c
c
c
c coxnet utility routines:
c
c
c same as elnet above:
c
c    call solns(ni,nx,lmu,ca,ia,nin,b)
c    call uncomp(ni,ca,ia,nin,a)
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor matrix:
c
c call cxmodval(ca,ia,nin,n,x,f);
c
c input:
c
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    n = number of predictor vectors (observations)
c    x(n,ni) = full (uncompressed) predictor matrix
c
c output:
c
c    f(n) = model predictions
c
c
c compute log-likelihood for given data set and vectors of coefficients
c
c call loglike(no,ni,x,y,d,o,w,nvec,a,flog,jerr)
c
c input:
c
c   no = number of observations
c   ni = number of predictor variables
c   x(no,ni) = predictor data matrix flat file
c   y(no) = observation times
c   d(no) = died/censored indicator
c       d(i)=0.0 => y(i) = censoring time
c       d(i)=1.0 => y(i) = death time
c   o(no) = observation off-sets
c   w(no)= observation weights
c   nvec = number of coefficient vectors
c   a(ni,nvec) = coefficient vectors (uncompressed)
c
c output
c
c   flog(nvec) = respective log-likelihood values
c   jerr = error flag - see coxnet above
c
c
c
c
c                Changing internal parameter values
c
c
c call chg_fract_dev(fdev)
c   fdev = minimum fractional change in deviance for stopping path
c      default = 1.0e-5
c
c call chg_dev_max(devmax)
c   devmax = maximum fraction of explained deviance for stopping path
c      default = 0.999
c
c call chg_min_flmin(eps)
c   eps = minimum value of flmin (see above). default= 1.0e-6
c
c call chg_big(big)
c   big = large floating point number. default = 9.9e35
c
c call chg_min_lambdas(mnlam)
c   mnlam = minimum number of path points (lambda values) allowed
c      default = 5
c
c call chg_min_null_prob(pmin)
c   pmin = minimum null probability for any class. default = 1.0e-9
c
c call chg _max_exp(exmx)
c   exmx = maximum allowed exponent. default = 250.0
c
c call chg_bnorm(prec,mxit)
c   prec = convergence threshold for multi response bounds adjustment
c          solution. default = 1.0e-10.
c   mxit = maximum iterations for multiresponse bounds adjustment solution
c          default = 100.
c
c
c             Obtain current internal parameter values
c
c call get_int_parms(fdev,eps,big,mnlam,devmax,pmin,exmx)
c call get_bnorm(prec,mxit);
c
c
c
"
subroutine get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx);
data sml0,eps0,big0,mnlam0,rsqmax0,pmin0,exmx0
   /1.0e-5,1.0e-6,9.9e35,5,0.999,1.0e-9,250.0/;
sml=sml0; eps=eps0; big=big0; mnlam=mnlam0; rsqmax=rsqmax0;
pmin=pmin0; exmx=exmx0;
return;
entry chg_fract_dev(arg); sml0=arg; return;
entry chg_dev_max(arg); rsqmax0=arg; return;
entry chg_min_flmin(arg); eps0=arg; return;
entry chg_big(arg); big0=arg; return;
entry chg_min_lambdas(irg); mnlam0=irg; return;
entry chg_min_null_prob(arg); pmin0=arg; return;
entry chg_max_exp(arg); exmx0=arg; return;
end;
subroutine elnet
 (ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,maxit,
   lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
real x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam),cl(2,ni);
real ulam(nlam),a0(nlam),rsq(nlam),alm(nlam);
integer jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: vq;
%mortran
if maxval(vp).le.0.0 < jerr=10000; return;>
allocate(vq(1:ni),stat=jerr); if(jerr.ne.0) return;
vq=max(0.0,vp); vq=vq*ni/sum(vq);
if ka.eq.1 <
   call elnetu
    (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,maxit,
      lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
>
else <
   call elnetn
    (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,maxit,
      lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
>
deallocate(vq);
return;
end;
subroutine elnetu
 (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,maxit,
   lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
real x(no,ni),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni);
real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam);
integer jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam
      integer, dimension (:), allocatable :: ju
%mortran
allocate(g(1:ni),stat=jerr);
allocate(xm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xs(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ju(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xv(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(vlam(1:nlam),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
call chkvars(no,ni,x,ju);
if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0;
if maxval(ju).le.0 < jerr=7777; return;>
call standard(no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr);
if(jerr.ne.0) return;
cl=cl/ys; if isd.gt.0 < <j=1,ni; cl(:,j)=cl(:,j)*xs(j);>>
if(flmin.ge.1.0) vlam=ulam/ys;
call elnet1(parm,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxit,xv,
   lmu,ca,ia,nin,rsq,alm,nlp,jerr);
if(jerr.gt.0) return;
<k=1,lmu; alm(k)=ys*alm(k); nk=nin(k);
   <l=1,nk; ca(l,k)=ys*ca(l,k)/xs(ia(l));> a0(k)=0.0;
   if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)));
>
deallocate(xm,xs,g,ju,xv,vlam);
return;
end;
subroutine standard (no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr);
real x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni); integer ju(ni);
%fortran
      real, dimension (:), allocatable :: v
%mortran
allocate(v(1:no),stat=jerr); if(jerr.ne.0) return;
w=w/sum(w); v=sqrt(w);
if intr.eq.0 < ym=0.0; y=v*y;
   ys=sqrt(dot_product(y,y)-dot_product(v,y)**2); y=y/ys;
   <j=1,ni; if(ju(j).eq.0) next; xm(j)=0.0; x(:,j)=v*x(:,j);
      xv(j)=dot_product(x(:,j),x(:,j));
      if isd.ne.0 < xbq=dot_product(v,x(:,j))**2; vc=xv(j)-xbq;
         xs(j)=sqrt(vc); x(:,j)=x(:,j)/xs(j); xv(j)=1.0+xbq/vc;
      >
      else < xs(j)=1.0;>
   >
>
else <
   <j=1,ni; if(ju(j).eq.0) next;
      xm(j)=dot_product(w,x(:,j)); x(:,j)=v*(x(:,j)-xm(j));
      xv(j)=dot_product(x(:,j),x(:,j)); if(isd.gt.0) xs(j)=sqrt(xv(j));
   >
   if isd.eq.0 < xs=1.0;>
   else <
      <j=1,ni; if(ju(j).eq.0) next; x(:,j)=x(:,j)/xs(j);>
      xv=1.0;
   >
   ym=dot_product(w,y); y=v*(y-ym); ys=sqrt(dot_product(y,y)); y=y/ys;
>
g=0.0; <j=1,ni; if(ju(j).ne.0) g(j)=dot_product(y,x(:,j));>
deallocate(v);
return;
end;
subroutine elnet1 (beta,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,ulam,thr,maxit,xv,
   lmu,ao,ia,kin,rsqo,almo,nlp,jerr);
real vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(nlam),xv(ni);
real cl(2,ni);
integer ju(ni),ia(nx),kin(nlam);
%fortran
      real, dimension (:), allocatable :: a,da
      integer, dimension (:), allocatable :: mm
      real, dimension (:,:), allocatable :: c
      allocate(c(1:ni,1:nx),stat=jerr)
%mortran
call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx);
allocate(a(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(da(1:ni),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
bta=beta; omb=1.0-bta;
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
rsq=0.0; a=0.0; mm=0; /nlp,nin/=0; iz=0; mnl=min(mnlam,nlam);
<m=1,nlam;
   if flmin.ge.1.0 < alm=ulam(m);>
   elseif m.gt.2 < alm=alm*alf;>
   elseif m.eq.1 < alm=big;>
   else < alm=0.0;
      <j=1,ni; if(ju(j).eq.0) next; if(vp(j).le.0.0) next;
         alm=max(alm,abs(g(j))/vp(j));
      >
      alm=alf*alm/max(bta,1.0e-3);
   >
   dem=alm*omb; ab=alm*bta; rsq0=rsq; jz=1;
   loop < if(iz*jz.ne.0) go to :b:; nlp=nlp+1; dlx=0.0;
      <k=1,ni; if(ju(k).eq.0) next;
         ak=a(k); u=g(k)+ak*xv(k); v=abs(u)-vp(k)*ab; a(k)=0.0;
         if(v.gt.0.0)
            a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)));
         if(a(k).eq.ak) next;
         if mm(k).eq.0 < nin=nin+1; if(nin.gt.nx) exit;
            <j=1,ni; if(ju(j).eq.0) next;
               if mm(j).ne.0 < c(j,nin)=c(k,mm(j)); next;>
               if j.eq.k < c(j,nin)=xv(j); next;>
               c(j,nin)=dot_product(x(:,j),x(:,k));
            >
            mm(k)=nin; ia(nin)=k;
         >
         del=a(k)-ak; rsq=rsq+del*(2.0*g(k)-del*xv(k));
         dlx=max(xv(k)*del**2,dlx);
         <j=1,ni; if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del;>
      >
      if(dlx.lt.thr) exit; if(nin.gt.nx) exit;
      if nlp.gt.maxit < jerr=-m; return;>
      :b: iz=1; da(1:nin)=a(ia(1:nin));
      loop < nlp=nlp+1; dlx=0.0;
         <l=1,nin; k=ia(l); ak=a(k); u=g(k)+ak*xv(k); v=abs(u)-vp(k)*ab;
            a(k)=0.0;
            if(v.gt.0.0)
               a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)));
            if(a(k).eq.ak) next;
            del=a(k)-ak; rsq=rsq+del*(2.0*g(k)-del*xv(k));
            dlx=max(xv(k)*del**2,dlx);
            <j=1,nin; g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del;>
         >
         if(dlx.lt.thr) exit; if nlp.gt.maxit < jerr=-m; return;>
      >
      da(1:nin)=a(ia(1:nin))-da(1:nin);
      <j=1,ni; if(mm(j).ne.0) next;
         if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin));
      >
      jz=0;      
   >
   if nin.gt.nx < jerr=-10000-m;  exit;>
   if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin)); kin(m)=nin;
   rsqo(m)=rsq; almo(m)=alm; lmu=m;
   if(m.lt.mnl) next; if(flmin.ge.1.0) next;
   me=0; <j=1,nin; if(ao(j,m).ne.0.0) me=me+1;> if(me.gt.ne) exit;
   if(rsq-rsq0.lt.sml*rsq) exit; if(rsq.gt.rsqmax) exit;
>
deallocate(a,mm,c,da);
return;
end;
subroutine elnetn (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
   intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
real vp(ni),x(no,ni),y(no),w(no),ulam(nlam),cl(2,ni);
real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam);
integer jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: xm,xs,xv,vlam
      integer, dimension (:), allocatable :: ju
%mortran
allocate(xm(1:ni),stat=jerr);
allocate(xs(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ju(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xv(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(vlam(1:nlam),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
call chkvars(no,ni,x,ju);
if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0;
if maxval(ju).le.0 < jerr=7777; return;>
call standard1(no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr);
if(jerr.ne.0) return;
cl=cl/ys; if isd.gt.0 < <j=1,ni; cl(:,j)=cl(:,j)*xs(j);>>
if(flmin.ge.1.0) vlam=ulam/ys;
call elnet2(parm,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,vlam,thr,maxit,xv,
   lmu,ca,ia,nin,rsq,alm,nlp,jerr);
if(jerr.gt.0) return;
<k=1,lmu; alm(k)=ys*alm(k); nk=nin(k);
   <l=1,nk; ca(l,k)=ys*ca(l,k)/xs(ia(l));> a0(k)=0.0;
   if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)));
>
deallocate(xm,xs,ju,xv,vlam);
return;
end;
subroutine standard1 (no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr);
real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni); integer ju(ni);
%fortran
      real, dimension (:), allocatable :: v
%mortran
allocate(v(1:no),stat=jerr); if(jerr.ne.0) return;
w=w/sum(w); v=sqrt(w);
if intr.eq.0 < ym=0.0; y=v*y;
   ys=sqrt(dot_product(y,y)-dot_product(v,y)**2); y=y/ys;
   <j=1,ni; if(ju(j).eq.0) next; xm(j)=0.0; x(:,j)=v*x(:,j);
      xv(j)=dot_product(x(:,j),x(:,j));
      if isd.ne.0 < xbq=dot_product(v,x(:,j))**2; vc=xv(j)-xbq;
         xs(j)=sqrt(vc); x(:,j)=x(:,j)/xs(j); xv(j)=1.0+xbq/vc;
      >
      else < xs(j)=1.0;>
   >
   go to :out:;
>
<j=1,ni; if(ju(j).eq.0) next;
   xm(j)=dot_product(w,x(:,j)); x(:,j)=v*(x(:,j)-xm(j));
   xv(j)=dot_product(x(:,j),x(:,j)); if(isd.gt.0) xs(j)=sqrt(xv(j));
>
if isd.eq.0 < xs=1.0;>
else < <j=1,ni; if(ju(j).eq.0) next; x(:,j)=x(:,j)/xs(j);>
   xv=1.0;
>
ym=dot_product(w,y); y=v*(y-ym); ys=sqrt(dot_product(y,y)); y=y/ys;
:out:deallocate(v);
return;
end;
subroutine elnet2(beta,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,ulam,thr,maxit,xv,
   lmu,ao,ia,kin,rsqo,almo,nlp,jerr);
real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(nlam),xv(ni);
real cl(2,ni);
integer ju(ni),ia(nx),kin(nlam);
%fortran
      real, dimension (:), allocatable :: a,g
      integer, dimension (:), allocatable :: mm,ix
%mortran
call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx);
allocate(a(1:ni),stat=jerr);
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(g(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ix(1:ni),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
bta=beta; omb=1.0-bta; ix=0;
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
rsq=0.0; a=0.0; mm=0; /nlp,nin/=0; iz=0; mnl=min(mnlam,nlam); alm=0.0;
<j=1,ni; if(ju(j).eq.0) next; g(j)=abs(dot_product(y,x(:,j)));>
<m=1,nlam; alm0=alm;
   if flmin.ge.1.0 < alm=ulam(m);>
   elseif m.gt.2 < alm=alm*alf;>
   elseif m.eq.1 < alm=big;>
   else < alm0=0.0;
      <j=1,ni; if(ju(j).eq.0) next; if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j));>
      alm0=alm0/max(bta,1.0e-3); alm=alf*alm0;
   >
   dem=alm*omb; ab=alm*bta; rsq0=rsq; jz=1;
   tlam=bta*(2.0*alm-alm0);
   <k=1,ni; if(ix(k).eq.1) next; if(ju(k).eq.0) next;
      if(g(k).gt.tlam*vp(k)) ix(k)=1;
   >
   loop < if(iz*jz.ne.0) go to :b:;
      :again:nlp=nlp+1; dlx=0.0;
      <k=1,ni; if(ix(k).eq.0) next; gk=dot_product(y,x(:,k));
         ak=a(k); u=gk+ak*xv(k); v=abs(u)-vp(k)*ab; a(k)=0.0;
         if(v.gt.0.0)
            a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)));
         if(a(k).eq.ak) next;
         if mm(k).eq.0 < nin=nin+1; if(nin.gt.nx) exit;
            mm(k)=nin; ia(nin)=k;
         >
         del=a(k)-ak; rsq=rsq+del*(2.0*gk-del*xv(k));
         y=y-del*x(:,k); dlx=max(xv(k)*del**2,dlx);         
      >
      if(nin.gt.nx) exit;
      if dlx.lt.thr < ixx=0;
         <k=1,ni; if(ix(k).eq.1) next; if(ju(k).eq.0) next;
            g(k)=abs(dot_product(y,x(:,k)));
            if g(k).gt.ab*vp(k) < ix(k)=1; ixx=1;>
         >
         if(ixx.eq.1) go to :again:;
         exit;
      >
      if nlp.gt.maxit < jerr=-m; return;> 
      :b: iz=1;
      loop < nlp=nlp+1; dlx=0.0;
         <l=1,nin; k=ia(l); gk=dot_product(y,x(:,k));
            ak=a(k); u=gk+ak*xv(k); v=abs(u)-vp(k)*ab; a(k)=0.0;
            if(v.gt.0.0)
               a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)));
            if(a(k).eq.ak) next;
            del=a(k)-ak; rsq=rsq+del*(2.0*gk-del*xv(k));
            y=y-del*x(:,k); dlx=max(xv(k)*del**2,dlx);
         >
         if(dlx.lt.thr) exit; if nlp.gt.maxit < jerr=-m; return;>
      >
      jz=0;
   >
   if nin.gt.nx < jerr=-10000-m;  exit;>
   if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin)); kin(m)=nin;
   rsqo(m)=rsq; almo(m)=alm; lmu=m;
   if(m.lt.mnl) next; if(flmin.ge.1.0) next;
   me=0; <j=1,nin; if(ao(j,m).ne.0.0) me=me+1;> if(me.gt.ne) exit;
   if(rsq-rsq0.lt.sml*rsq) exit; if(rsq.gt.rsqmax) exit;
>
deallocate(a,mm,g,ix);
return;
end;
subroutine chkvars(no,ni,x,ju);
real x(no,ni); integer ju(ni);
<j=1,ni; ju(j)=0; t=x(1,j);
   <i=2,no; if(x(i,j).eq.t) next; ju(j)=1; exit;>
>
return;
end;
subroutine uncomp(ni,ca,ia,nin,a);
real ca(*),a(ni); integer ia(*);
a=0.0; if(nin.gt.0) a(ia(1:nin))=ca(1:nin);
return;
end;
subroutine modval(a0,ca,ia,nin,n,x,f);
real ca(nin),x(n,*),f(n); integer ia(nin);
f=a0; if(nin.le.0) return;
<i=1,n; f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)));>
return;
end;
subroutine spelnet
 (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,
   maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
real x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni);
real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam);
integer ix(*),jx(*),jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: vq;
%mortran
if maxval(vp).le.0.0 < jerr=10000; return;>
allocate(vq(1:ni),stat=jerr); if(jerr.ne.0) return;
vq=max(0.0,vp); vq=vq*ni/sum(vq);
if ka.eq.1 <
   call spelnetu
    (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,isd,
      intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
>
else <
   call spelnetn
    (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,
      maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
>
deallocate(vq);
return;
end;
subroutine spelnetu
 (parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,intr,
   maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
real x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni);
real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam);
integer ix(*),jx(*),jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam
      integer, dimension (:), allocatable :: ju
%mortran
allocate(g(1:ni),stat=jerr);
allocate(xm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xs(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ju(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xv(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(vlam(1:nlam),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
call spchkvars(no,ni,x,ix,ju);
if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0;
if maxval(ju).le.0 < jerr=7777; return;>
call spstandard(no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys,xv,jerr);
if(jerr.ne.0) return;
cl=cl/ys; if isd.gt.0 < <j=1,ni; cl(:,j)=cl(:,j)*xs(j);>>
if(flmin.ge.1.0) vlam=ulam/ys;
call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vlam,thr,maxit,
   xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr);
if(jerr.gt.0) return;
<k=1,lmu; alm(k)=ys*alm(k); nk=nin(k);
   <l=1,nk; ca(l,k)=ys*ca(l,k)/xs(ia(l));> a0(k)=0.0;
   if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)));
>
deallocate(xm,xs,g,ju,xv,vlam);
return;
end;
subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys,xv,jerr);
real x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni); integer ix(*),jx(*),ju(ni);
w=w/sum(w);
if intr.eq.0 < ym=0.0;
   ys=sqrt(dot_product(w,y**2)-dot_product(w,y)**2); y=y/ys;
   <j=1,ni; if(ju(j).eq.0) next; xm(j)=0.0; jb=ix(j); je=ix(j+1)-1; 
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2);
      if isd.ne.0 < xbq=dot_product(w(jx(jb:je)),x(jb:je))**2; vc=xv(j)-xbq;
         xs(j)=sqrt(vc); xv(j)=1.0+xbq/vc;
      >
      else < xs(j)=1.0;>
   >
>
else <
   <j=1,ni; if(ju(j).eq.0) next;
      jb=ix(j); je=ix(j+1)-1; xm(j)=dot_product(w(jx(jb:je)),x(jb:je));
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2;
      if(isd.gt.0) xs(j)=sqrt(xv(j));
   >
   if isd.eq.0 < xs=1.0;> else < xv=1.0;>
   ym=dot_product(w,y); y=y-ym; ys=sqrt(dot_product(w,y**2)); y=y/ys;
>
g=0.0;
<j=1,ni; if(ju(j).eq.0) next; jb=ix(j); je=ix(j+1)-1;
   g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j);
>
return;
end;
subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,ulam,
   thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr);
real g(ni),vp(ni),x(*),ulam(nlam),w(no);
real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni),cl(2,ni);
integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam);
%fortran
      real, dimension (:), allocatable :: a,da
      integer, dimension (:), allocatable :: mm
      real, dimension (:,:), allocatable :: c
      allocate(c(1:ni,1:nx),stat=jerr)
%mortran
call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx);
allocate(a(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(da(1:ni),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
bta=beta; omb=1.0-bta;
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
rsq=0.0; a=0.0; mm=0; /nlp,nin/=0; iz=0; mnl=min(mnlam,nlam);
<m=1,nlam;
   if flmin.ge.1.0 < alm=ulam(m);>
   elseif m.gt.2 < alm=alm*alf;>
   elseif m.eq.1 < alm=big;>
   else < alm=0.0;
      <j=1,ni; if(ju(j).eq.0) next; if(vp(j).le.0.0) next;
         alm=max(alm,abs(g(j))/vp(j));
      >
      alm=alf*alm/max(bta,1.0e-3);
   >
   dem=alm*omb; ab=alm*bta; rsq0=rsq; jz=1;
   loop < if(iz*jz.ne.0) go to :b:; nlp=nlp+1; dlx=0.0;
      <k=1,ni; if(ju(k).eq.0) next;
         ak=a(k); u=g(k)+ak*xv(k); v=abs(u)-vp(k)*ab; a(k)=0.0;
         if(v.gt.0.0)
            a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)));
         if(a(k).eq.ak) next;
         if mm(k).eq.0 < nin=nin+1; if(nin.gt.nx) exit;
            <j=1,ni; if(ju(j).eq.0) next;
               if mm(j).ne.0 < c(j,nin)=c(k,mm(j)); next;>
               if j.eq.k < c(j,nin)=xv(j); next;>
               c(j,nin)=
                  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k));
            >
            mm(k)=nin; ia(nin)=k;
         >
         del=a(k)-ak; rsq=rsq+del*(2.0*g(k)-del*xv(k));
         dlx=max(xv(k)*del**2,dlx);
         <j=1,ni; if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del;>
      >
      if(dlx.lt.thr) exit; if(nin.gt.nx) exit;
      if nlp.gt.maxit < jerr=-m; return;>
      :b: iz=1; da(1:nin)=a(ia(1:nin));
      loop < nlp=nlp+1; dlx=0.0;
         <l=1,nin; k=ia(l);
            ak=a(k); u=g(k)+ak*xv(k); v=abs(u)-vp(k)*ab; a(k)=0.0;
            if(v.gt.0.0)
               a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)));
            if(a(k).eq.ak) next;
            del=a(k)-ak; rsq=rsq+del*(2.0*g(k)-del*xv(k));
            dlx=max(xv(k)*del**2,dlx);
            <j=1,nin; g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del;>
         >
         if(dlx.lt.thr) exit; if nlp.gt.maxit < jerr=-m; return;>
      >
      da(1:nin)=a(ia(1:nin))-da(1:nin);
      <j=1,ni; if(mm(j).ne.0) next;
         if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin));
      >
      jz=0;
   >
   if nin.gt.nx < jerr=-10000-m;  exit;>
   if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin)); kin(m)=nin;
   rsqo(m)=rsq; almo(m)=alm; lmu=m;
   if(m.lt.mnl) next; if(flmin.ge.1.0) next;
   me=0; <j=1,nin; if(ao(j,m).ne.0.0) me=me+1;> if(me.gt.ne) exit;
   if(rsq-rsq0.lt.sml*rsq) exit; if(rsq.gt.rsqmax) exit;
>
deallocate(a,mm,c,da);
return;
end;
subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,
   thr,isd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
real x(*),vp(ni),y(no),w(no),ulam(nlam),cl(2,ni);
real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam);
integer ix(*),jx(*),jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: xm,xs,xv,vlam
      integer, dimension (:), allocatable :: ju
%mortran
allocate(xm(1:ni),stat=jerr);
allocate(xs(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ju(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xv(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(vlam(1:nlam),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
call spchkvars(no,ni,x,ix,ju);
if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0;
if maxval(ju).le.0 < jerr=7777; return;>
call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,xv,jerr);
if(jerr.ne.0) return;
cl=cl/ys; if isd.gt.0 < <j=1,ni; cl(:,j)=cl(:,j)*xs(j);>>
if(flmin.ge.1.0) vlam=ulam/ys;
call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vlam,thr,maxit,
   xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr);
if(jerr.gt.0) return;
<k=1,lmu; alm(k)=ys*alm(k); nk=nin(k);
   <l=1,nk; ca(l,k)=ys*ca(l,k)/xs(ia(l));> a0(k)=0.0;
   if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)));
>
deallocate(xm,xs,ju,xv,vlam);
return;
end;
subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,xv,jerr);
real x(*),y(no),w(no),xm(ni),xs(ni),xv(ni); integer ix(*),jx(*),ju(ni);
w=w/sum(w);
if intr.eq.0 < ym=0.0;
   ys=sqrt(dot_product(w,y**2)-dot_product(w,y)**2); y=y/ys;
   <j=1,ni; if(ju(j).eq.0) next; xm(j)=0.0; jb=ix(j); je=ix(j+1)-1;
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2);
      if isd.ne.0 < xbq=dot_product(w(jx(jb:je)),x(jb:je))**2; vc=xv(j)-xbq;
         xs(j)=sqrt(vc); xv(j)=1.0+xbq/vc;
      >
      else < xs(j)=1.0;>
   >
   return;
>
<j=1,ni; if(ju(j).eq.0) next;
   jb=ix(j); je=ix(j+1)-1; xm(j)=dot_product(w(jx(jb:je)),x(jb:je));
   xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2;
   if(isd.gt.0) xs(j)=sqrt(xv(j));
>
if isd.eq.0 < xs=1.0;> else < xv=1.0;>
ym=dot_product(w,y); y=y-ym; ys=sqrt(dot_product(w,y**2)); y=y/ys;
return;
end;
subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,ulam,
   thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr);
real y(no),w(no),x(*),vp(ni),ulam(nlam),cl(2,ni);
real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni);
integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam);
%fortran
      real, dimension (:), allocatable :: a,g
      integer, dimension (:), allocatable :: mm,iy
%mortran
call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx);
allocate(a(1:ni),stat=jerr);
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(g(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(iy(1:ni),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
bta=beta; omb=1.0-bta; alm=0.0; iy=0;
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
rsq=0.0; a=0.0; mm=0; o=0.0; /nlp,nin/=0; iz=0; mnl=min(mnlam,nlam);
<j=1,ni; if(ju(j).eq.0) next;
   jb=ix(j); je=ix(j+1)-1;
   g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j));
>
<m=1,nlam; alm0=alm;
   if flmin.ge.1.0 < alm=ulam(m);>
   elseif m.gt.2 < alm=alm*alf;>
   elseif m.eq.1 < alm=big;>
   else < alm0=0.0;
      <j=1,ni; if(ju(j).eq.0) next; if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j));>
      alm0=alm0/max(bta,1.0e-3); alm=alf*alm0;
   >   
   dem=alm*omb; ab=alm*bta; rsq0=rsq; jz=1;
   tlam=bta*(2.0*alm-alm0);
   <k=1,ni; if(iy(k).eq.1) next; if(ju(k).eq.0) next;
      if(g(k).gt.tlam*vp(k)) iy(k)=1;
   >   
   loop < if(iz*jz.ne.0) go to :b:;
      :again:nlp=nlp+1; dlx=0.0;
      <k=1,ni; if(iy(k).eq.0) next; jb=ix(k); je=ix(k+1)-1;
         gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k);
         ak=a(k); u=gk+ak*xv(k); v=abs(u)-vp(k)*ab; a(k)=0.0;
         if(v.gt.0.0)
            a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)));
         if(a(k).eq.ak) next;
         if mm(k).eq.0 < nin=nin+1; if(nin.gt.nx) exit;
            mm(k)=nin; ia(nin)=k;
         >
         del=a(k)-ak; rsq=rsq+del*(2.0*gk-del*xv(k));
         y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k);
         o=o+del*xm(k)/xs(k); dlx=max(xv(k)*del**2,dlx);         
      >
      if(nin.gt.nx) exit;
      if dlx.lt.thr < ixx=0;
         <j=1,ni; if(iy(j).eq.1) next; if(ju(j).eq.0) next;
            jb=ix(j); je=ix(j+1)-1;
            g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j));
            if g(j).gt.ab*vp(j) < iy(j)=1; ixx=1;>
         >
         if(ixx.eq.1) go to :again:;
         exit;
      >
      if nlp.gt.maxit < jerr=-m; return;>
      :b: iz=1;
      loop < nlp=nlp+1; dlx=0.0;
         <l=1,nin; k=ia(l); jb=ix(k); je=ix(k+1)-1;
            gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k);
            ak=a(k); u=gk+ak*xv(k); v=abs(u)-vp(k)*ab; a(k)=0.0;
            if(v.gt.0.0)
               a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*dem)));
            if(a(k).eq.ak) next;
            del=a(k)-ak; rsq=rsq+del*(2.0*gk-del*xv(k));
            y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k);
            o=o+del*xm(k)/xs(k); dlx=max(xv(k)*del**2,dlx);
         >
         if(dlx.lt.thr) exit; if nlp.gt.maxit < jerr=-m; return;>            
      >
      jz=0;
   >
   if nin.gt.nx < jerr=-10000-m;  exit;>
   if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin)); kin(m)=nin;
   rsqo(m)=rsq; almo(m)=alm; lmu=m;
   if(m.lt.mnl) next; if(flmin.ge.1.0) next;
   me=0; <j=1,nin; if(ao(j,m).ne.0.0) me=me+1;> if(me.gt.ne) exit;
   if(rsq-rsq0.lt.sml*rsq) exit; if(rsq.gt.rsqmax) exit;
>
deallocate(a,mm,g,iy);
return;
end;
subroutine spchkvars(no,ni,x,ix,ju);
real x(*); integer ix(*),ju(ni);
<j=1,ni; ju(j)=0; jb=ix(j); nj=ix(j+1)-jb; if(nj.eq.0) next;
   je=ix(j+1)-1;
   if nj.lt.no < <i=jb,je; if(x(i).eq.0.0) next; ju(j)=1; exit;>>
   else < t=x(jb); <i=jb+1,je; if(x(i).eq.t) next; ju(j)=1; exit;>>
>
return;
end;
subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f);
real ca(*),x(*),f(n); integer ia(*),ix(*),jx(*);
f=a0;
<j=1,nin; k=ia(j); kb=ix(k); ke=ix(k+1)-1;
   f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke);
>
return;
end;
function row_prod(i,j,ia,ja,ra,w);
integer ia(*),ja(*); real ra(*),w(*);
row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),
        ia(i+1)-ia(i),ia(j+1)-ia(j),w);
return;
end;
function dot(x,y,mx,my,nx,ny,w);
real x(*),y(*),w(*); integer mx(*),my(*);
/i,j/=1; s=0.0;
loop <
   until mx(i).ge.my(j) < i=i+1; if(i.gt.nx) go to :done:;>
   if(mx(i).eq.my(j)) go to :equal:;
   until my(j).ge.mx(i) < j=j+1; if(j.gt.ny) go to :done:;>
   if(mx(i).eq.my(j)) go to :equal:; next;
   :equal:s=s+w(mx(i))*x(i)*y(j);
   i=i+1; if(i.gt.nx) exit; j=j+1; if(j.gt.ny) exit;
>
:done: dot=s;
return;
end;
subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
   isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr);
real x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam);
real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni);
integer jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: xm,xs,ww,vq,xv
      integer, dimension (:), allocatable :: ju
%mortran
if maxval(vp).le.0.0 < jerr=10000; return;>
allocate(ww(1:no),stat=jerr);
allocate(ju(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(vq(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xm(1:ni),stat=ierr); jerr=jerr+ierr;
if kopt.eq.2 < allocate(xv(1:ni),stat=ierr); jerr=jerr+ierr;>
if isd.gt.0 < allocate(xs(1:ni),stat=ierr); jerr=jerr+ierr;>
if(jerr.ne.0) return;
call chkvars(no,ni,x,ju);
if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0;
if maxval(ju).le.0 < jerr=7777; return;>
vq=max(0.0,vp); vq=vq*ni/sum(vq);
<i=1,no; ww(i)=sum(y(i,:)); if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i);>
sw=sum(ww); ww=ww/sw;
if nc.eq.1 < call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs);
   if isd.gt.0 < <j=1,ni; cl(:,j)=cl(:,j)*xs(j);>>
   call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,
         thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr);
>
elseif kopt.eq.2 < call multlstandard1(no,ni,x,ww,ju,isd,intr,xm,xs,xv);
   if isd.gt.0 < <j=1,ni; cl(:,j)=cl(:,j)*xs(j);>>
   call multlognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,thr,
         intr,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr);
>
else < call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs);
   if isd.gt.0 < <j=1,ni; cl(:,j)=cl(:,j)*xs(j);>>
   call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,thr,
         isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr);
>
if(jerr.gt.0) return; dev0=2.0*sw*dev0;
<k=1,lmu; nk=nin(k);
   <ic=1,nc; if isd.gt.0 < <l=1,nk; ca(l,ic,k)=ca(l,ic,k)/xs(ia(l));>>
      if intr.eq.0 < a0(ic,k)=0.0;>
      else < a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)));>
   >
>
deallocate(ww,ju,vq,xm); if(isd.gt.0) deallocate(xs);
if(kopt.eq.2) deallocate(xv);
return;
end;
subroutine lstandard1 (no,ni,x,w,ju,isd,intr,xm,xs);
real x(no,ni),w(no),xm(ni),xs(ni); integer ju(ni);
if intr.eq.0 <
   <j=1,ni; if(ju(j).eq.0) next; xm(j)=0.0;
      if isd.ne.0 < vc=dot_product(w,x(:,j)**2)-dot_product(w,x(:,j))**2;
         xs(j)=sqrt(vc); x(:,j)=x(:,j)/xs(j);
      >
   >
   return;
>
<j=1,ni; if(ju(j).eq.0) next;
   xm(j)=dot_product(w,x(:,j)); x(:,j)=x(:,j)-xm(j);
   if isd.gt.0 < xs(j)=sqrt(dot_product(w,x(:,j)**2)); x(:,j)=x(:,j)/xs(j);>
>
return;
end;
subroutine multlstandard1 (no,ni,x,w,ju,isd,intr,xm,xs,xv);
real x(no,ni),w(no),xm(ni),xs(ni),xv(ni); integer ju(ni);
if intr.eq.0 <
   <j=1,ni; if(ju(j).eq.0) next; xm(j)=0.0;
      xv(j)=dot_product(w,x(:,j)**2);
      if isd.ne.0 < xbq=dot_product(w,x(:,j))**2; vc=xv(j)-xbq;
         xs(j)=sqrt(vc); x(:,j)=x(:,j)/xs(j); xv(j)=1.0+xbq/vc;
      >
   >
   return;
>
<j=1,ni; if(ju(j).eq.0) next;
   xm(j)=dot_product(w,x(:,j)); x(:,j)=x(:,j)-xm(j);
   xv(j)=dot_product(w,x(:,j)**2);
   if isd.gt.0 < xs(j)=sqrt(xv(j)); x(:,j)=x(:,j)/xs(j); xv(j)=1.0;>
>
return;
end;
subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin,ulam,shri,
    isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr);
real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni);
real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam);
integer ju(ni),m(nx),kin(nlam);
%fortran
      real, dimension (:), allocatable :: b,bs,v,r,xv,q,ga
      integer, dimension (:), allocatable :: mm,ixx
%mortran
call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx);
allocate(b(0:ni),stat=jerr);
allocate(xv(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ga(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(bs(0:ni),stat=ierr); jerr=jerr+ierr;
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ixx(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(r(1:no),stat=ierr); jerr=jerr+ierr;
allocate(v(1:no),stat=ierr); jerr=jerr+ierr;
allocate(q(1:no),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
fmax=log(1.0/pmin-1.0); fmin=-fmax; vmin=(1.0+pmin)*pmin*(1.0-pmin);
bta=parm; omb=1.0-bta;
q0=dot_product(w,y); if q0.le.pmin < jerr=8001; return;>
if q0.ge.1.0-pmin < jerr=9001; return;>
if(intr.eq.0.0) q0=0.5;
ixx=0; al=0.0; bz=0.0; if(intr.ne.0) bz=log(q0/(1.0-q0));
if nonzero(no,g).eq.0 < vi=q0*(1.0-q0); b(0)=bz; v=vi*w;
   r=w*(y-q0); q=q0; xmz=vi; dev1=-(bz*q0+log(1.0-q0));
>
else < b(0)=0.0;
   if intr.ne.0 < b(0)=azero(no,y,g,w,jerr); if(jerr.ne.0) return;>
   q=1.0/(1.0+exp(-b(0)-g)); v=w*q*(1.0-q); r=w*(y-q); xmz=sum(v);
   dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)));
>
if kopt.gt.0 <
   if isd.gt.0.and.intr.ne.0 < xv=0.25;>
   else < <j=1,ni; if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2);>>
>
dev0=dev1;
<i=1,no; if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i));
   if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i));
>
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
m=0; mm=0; /nlp,nin/=0; mnl=min(mnlam,nlam); bs=0.0; b(1:ni)=0.0;
shr=shri*dev0;
<j=1,ni; if(ju(j).eq.0) next; ga(j)=abs(dot_product(r,x(:,j)));>
<ilm=1,nlam; al0=al;
   if flmin.ge.1.0 < al=ulam(ilm);>
   elseif ilm.gt.2 < al=al*alf;>
   elseif ilm.eq.1 < al=big;>
   else < al0=0.0;
      <j=1,ni; if(ju(j).eq.0) next; if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j));>
      al0=al0/max(bta,1.0e-3); al=alf*al0;
   >
   al2=al*omb; al1=al*bta; tlam=bta*(2.0*al-al0);
   <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1;
   >
   :again:continue;
   loop < bs(0)=b(0); if(nin.gt.0) bs(m(1:nin))=b(m(1:nin));
      if kopt.eq.0 <
         <j=1,ni; if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2);>
      >
      loop < nlp=nlp+1; dlx=0.0;
         <k=1,ni; if(ixx(k).eq.0) next;
            bk=b(k); gk=dot_product(r,x(:,k));
            u=gk+xv(k)*b(k); au=abs(u)-vp(k)*al1;
            if au.le.0.0 < b(k)=0.0;>
            else <
               b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)));
            >
            d=b(k)-bk; if(abs(d).le.0.0) next; dlx=max(dlx,xv(k)*d**2);
            r=r-d*v*x(:,k);
            if mm(k).eq.0 < nin=nin+1; if(nin.gt.nx) exit;
               mm(k)=nin; m(nin)=k;
            >
         >
         if(nin.gt.nx) exit;
         d=0.0; if(intr.ne.0) d=sum(r)/xmz;
         if d.ne.0.0 < b(0)=b(0)+d; dlx=max(dlx,xmz*d**2); r=r-d*v;>
         if(dlx.lt.shr) exit; if nlp.gt.maxit < jerr=-ilm; return;>
         loop < nlp=nlp+1; dlx=0.0;
            <l=1,nin; k=m(l); bk=b(k);
               gk=dot_product(r,x(:,k));
               u=gk+xv(k)*b(k); au=abs(u)-vp(k)*al1;
               if au.le.0.0 < b(k)=0.0;>
               else <
                  b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)));
               >
               d=b(k)-bk; if(abs(d).le.0.0) next; dlx=max(dlx,xv(k)*d**2);
               r=r-d*v*x(:,k);
            >
            d=0.0; if(intr.ne.0) d=sum(r)/xmz;
            if d.ne.0.0 < b(0)=b(0)+d; dlx=max(dlx,xmz*d**2); r=r-d*v;>
            if(dlx.lt.shr) exit; if nlp.gt.maxit < jerr=-ilm; return;>
         >
      >
      if(nin.gt.nx) exit;
      <i=1,no; fi=b(0)+g(i);
         if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)));
         if fi.lt.fmin < q(i)=0.0;> elseif fi.gt.fmax < q(i)=1.0;>
         else < q(i)=1.0/(1.0+exp(-fi));> 
      >
      v=w*q*(1.0-q); xmz=sum(v); if(xmz.le.vmin) exit; r=w*(y-q);
      if xmz*(b(0)-bs(0))**2.lt.shr < ix=0;
         <j=1,nin; k=m(j);
            if(xv(k)*(b(k)-bs(k))**2.lt.shr) next; ix=1; exit;
         >
         if ix.eq.0 <
            <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
               ga(k)=abs(dot_product(r,x(:,k)));
               if ga(k).gt.al1*vp(k) < ixx(k)=1; ix=1;>
            >
            if(ix.eq.1) go to :again:;
            exit;
         >
      >      
   >
   if nin.gt.nx < jerr=-10000-ilm;  exit;>
   if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin)); kin(ilm)=nin;
   a0(ilm)=b(0); alm(ilm)=al; lmu=ilm;
   devi=dev2(no,w,y,q,pmin);
   dev(ilm)=(dev1-devi)/dev0; if(xmz.le.vmin) exit;
   if(ilm.lt.mnl) next; if(flmin.ge.1.0) next;
   me=0; <j=1,nin; if(a(j,ilm).ne.0.0) me=me+1;> if(me.gt.ne) exit;
   if(dev(ilm).gt.devmax) exit; if(dev(ilm)-dev(ilm-1).lt.sml) exit;
>
g=log(q/(1.0-q));
deallocate(b,bs,v,r,xv,q,mm,ga,ixx);
return;
end;
function dev2(n,w,y,p,pmin);
real w(n),y(n),p(n);
pmax=1.0-pmin; s=0.0;
<i=1,n; pi=min(max(pmin,p(i)),pmax);
   s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi));
>
dev2=s;
return;
end;
function azero(n,y,g,q,jerr);
parameter(eps=1.0e-7);
real y(n),g(n),q(n);
%fortran
      real, dimension (:), allocatable :: e,p,w
%mortran
allocate(e(1:n),stat=jerr);
allocate(p(1:n),stat=ierr); jerr=jerr+ierr;
allocate(w(1:n),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
az=0.0; e=exp(-g); qy=dot_product(q,y); p=1.0/(1.0+e);
loop < w=q*p*(1.0-p);
   d=(qy-dot_product(q,p))/sum(w); az=az+d; if(abs(d).lt.eps) exit;
   ea0=exp(-az); p=1.0/(1.0+ea0*e);
>
azero=az;
deallocate(e,p,w);
return;
end;
subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin,ulam,shri,
    isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr);
real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam);
real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni);
integer ju(ni),m(nx),kin(nlam);
%fortran
      real, dimension (:,:), allocatable :: q
      real, dimension (:), allocatable :: sxp,sxpl
      real, dimension (:), allocatable :: di,v,r,ga
      real, dimension (:,:), allocatable :: b,bs,xv
      integer, dimension (:), allocatable :: mm,is,ixx
      allocate(b(0:ni,1:nc),stat=jerr)
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr     
%mortran
call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx); exmn=-exmx;
allocate(r(1:no),stat=ierr); jerr=jerr+ierr;
allocate(v(1:no),stat=ierr); jerr=jerr+ierr;
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(is(1:max(nc,ni)),stat=ierr); jerr=jerr+ierr;
allocate(sxp(1:no),stat=ierr); jerr=jerr+ierr;
allocate(sxpl(1:no),stat=ierr); jerr=jerr+ierr;
allocate(di(1:no),stat=ierr); jerr=jerr+ierr;
allocate(ga(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ixx(1:ni),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
pmax=1.0-pmin; emin=pmin/pmax; emax=1.0/emin;
pfm=(1.0+pmin)*pmin; pfx=(1.0-pmin)*pmax; vmin=pfm*pmax;
bta=parm; omb=1.0-bta; dev1=0.0; dev0=0.0;
<ic=1,nc; q0=dot_product(w,y(:,ic));
   if q0.le.pmin < jerr =8000+ic; return;>
   if q0.ge.1.0-pmin < jerr =9000+ic; return;>
   if intr.eq.0 < q0=1.0/nc; b(0,ic)=0.0;>
   else < b(0,ic)=log(q0); dev1=dev1-q0*b(0,ic);>
   b(1:ni,ic)=0.0;
>
if(intr.eq.0) dev1=log(float(nc)); ixx=0; al=0.0;
if nonzero(no*nc,g).eq.0 <
   b(0,:)=b(0,:)-sum(b(0,:))/nc; sxp=0.0;
   <ic=1,nc; q(:,ic)=exp(b(0,ic)); sxp=sxp+q(:,ic);>
>
else < <i=1,no; g(i,:)=g(i,:)-sum(g(i,:))/nc;> sxp=0.0;
   if intr.eq.0 < b(0,:)=0.0;>
   else < call kazero(nc,no,y,g,w,b(0,:),jerr); if(jerr.ne.0) return;>
   dev1=0.0;
   <ic=1,nc; q(:,ic)=b(0,ic)+g(:,ic);
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic)); 
      q(:,ic)=exp(q(:,ic)); sxp=sxp+q(:,ic);
   >
   sxpl=w*log(sxp); <ic=1,nc; dev1=dev1+dot_product(y(:,ic),sxpl);>
>
<ic=1,nc; <i=1,no; if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic));>>
dev0=dev0+dev1;
if kopt.gt.0 <
   if isd.gt.0.and.intr.ne.0 < xv=0.25;>
   else < <j=1,ni; if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2);>>
>
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
m=0; mm=0; nin=0; nlp=0; mnl=min(mnlam,nlam); bs=0.0; shr=shri*dev0;
ga=0.0;
<ic=1,nc; r=w*(y(:,ic)-q(:,ic)/sxp);
   <j=1,ni; if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))));>
>
<ilm=1,nlam; al0=al;
   if flmin.ge.1.0 < al=ulam(ilm);>
   elseif ilm.gt.2 < al=al*alf;>
   elseif ilm.eq.1 < al=big;>
   else < al0=0.0;
      <j=1,ni; if(ju(j).eq.0) next; if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j));>
      al0=al0/max(bta,1.0e-3); al=alf*al0;
   >
   al2=al*omb; al1=al*bta; tlam=bta*(2.0*al-al0);
   <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1;
   >
   :again:continue;
   loop < /ix,jx/=0; ig=0;
      <ic=1,nc; bs(0,ic)=b(0,ic);
         if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic);
         xmz=0.0;
         <i=1,no; pic=q(i,ic)/sxp(i);
            if pic.lt.pfm < pic=0.0; v(i)=0.0;>
            elseif pic.gt.pfx < pic=1.0; v(i)=0.0;>
            else < v(i)=w(i)*pic*(1.0-pic); xmz=xmz+v(i);>
            r(i)=w(i)*(y(i,ic)-pic);
         >
         if(xmz.le.vmin) next; ig=1;
         if kopt.eq.0 <
           <j=1,ni; if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2);>
         >     
         loop < nlp=nlp+1; dlx=0.0;
            <k=1,ni; if(ixx(k).eq.0) next;
               bk=b(k,ic); gk=dot_product(r,x(:,k));
               u=gk+xv(k,ic)*b(k,ic); au=abs(u)-vp(k)*al1;
               if au.le.0.0 < b(k,ic)=0.0;>
               else <
                  b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/
                     (xv(k,ic)+vp(k)*al2)));
               >
               d=b(k,ic)-bk; if(abs(d).le.0.0) next;
               dlx=max(dlx,xv(k,ic)*d**2); r=r-d*v*x(:,k);
               if mm(k).eq.0 < nin=nin+1;
                  if nin.gt.nx < jx=1; exit;>
                  mm(k)=nin; m(nin)=k;
               >
            >
            if(jx.gt.0) exit;
            d=0.0; if(intr.ne.0) d=sum(r)/xmz;
            if d.ne.0.0 < b(0,ic)=b(0,ic)+d; dlx=max(dlx,xmz*d**2); r=r-d*v;>
            if(dlx.lt.shr) exit;
            if nlp.gt.maxit < jerr=-ilm; return;>
            loop < nlp=nlp+1; dlx=0.0;
               <l=1,nin; k=m(l); bk=b(k,ic);
                  gk=dot_product(r,x(:,k));
                  u=gk+xv(k,ic)*b(k,ic); au=abs(u)-vp(k)*al1;
                  if au.le.0.0 < b(k,ic)=0.0;>
                  else <
                     b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/
                        (xv(k,ic)+vp(k)*al2)));
                  >
                  d=b(k,ic)-bk; if(abs(d).le.0.0) next;
                  dlx=max(dlx,xv(k,ic)*d**2); r=r-d*v*x(:,k);
               >
               d=0.0; if(intr.ne.0) d=sum(r)/xmz;
               if d.ne.0.0 < b(0,ic)=b(0,ic)+d;
                  dlx=max(dlx,xmz*d**2); r=r-d*v;
               >
               if(dlx.lt.shr) exit; if nlp.gt.maxit < jerr=-ilm; return;>                
            >
         >
         if(jx.gt.0) exit;
         if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1;
         if ix.eq.0 <
            <j=1,nin; k=m(j); 
               if xv(k,ic)*(b(k,ic)-bs(k,ic))**2.gt.shr < ix=1; exit;>               
            >
         >       
         <i=1,no; fi=b(0,ic)+g(i,ic);
            if(nin.gt.0)
              fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)));
            fi=min(max(exmn,fi),exmx); sxp(i)=sxp(i)-q(i,ic); 
            q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i));
            sxp(i)=sxp(i)+q(i,ic);
         >     
      >
      s=-sum(b(0,:))/nc; b(0,:)=b(0,:)+s; di=s;
      <j=1,nin; l=m(j);
         if vp(l).le.0.0 < s=sum(b(l,:))/nc;>
         else < s=elc(parm,nc,cl(:,l),b(l,:),is);>
         b(l,:)=b(l,:)-s; di=di-s*x(:,l);
      > 
      di=exp(di); sxp=sxp*di; <ic=1,nc; q(:,ic)=q(:,ic)*di;>
      if(jx.gt.0) exit; if(ig.eq.0) exit;
      if ix.eq.0 <
         <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next; ga(k)=0.0;>
         <ic=1,nc; r=w*(y(:,ic)-q(:,ic)/sxp);
            <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
               ga(k)=max(ga(k),abs(dot_product(r,x(:,k))));
            >
         >
         <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
            if ga(k).gt.al1*vp(k) < ixx(k)=1; ix=1;>
         >
         if(ix.eq.1) go to :again:;
         exit;
      >
   > 
   if jx.gt.0 < jerr=-10000-ilm;  exit;> devi=0.0;
   <ic=1,nc;
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic); a0(ic,ilm)=b(0,ic);
      <i=1,no; if (y(i,ic).le.0.0) next;
         devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i));
      >
   >
   kin(ilm)=nin; alm(ilm)=al; lmu=ilm;
   dev(ilm)=(dev1-devi)/dev0; if(ig.eq.0) exit;
   if(ilm.lt.mnl) next; if(flmin.ge.1.0) next;
   if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne) exit;
   if(dev(ilm).gt.devmax) exit; if(dev(ilm)-dev(ilm-1).lt.sml) exit;
>
g=log(q); <i=1,no; g(i,:)=g(i,:)-sum(g(i,:))/nc;>
deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx);
return;
end;
subroutine kazero(kk,n,y,g,q,az,jerr);
parameter(eps=1.0e-7);
real y(n,kk),g(n,kk),q(n),az(kk);
%fortran
      real, dimension (:), allocatable :: s
      real, dimension (:,:), allocatable :: e
      allocate(e(1:n,1:kk),stat=jerr)
%mortran
allocate(s(1:n),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
az=0.0; e=exp(g); <i=1,n; s(i)=sum(e(i,:));>
loop < dm=0.0;
   <k=1,kk; /t,u/=0.0;
      <i=1,n; pik=e(i,k)/s(i);
         t=t+q(i)*(y(i,k)-pik); u=u+q(i)*pik*(1.0-pik);
      >
      d=t/u; az(k)=az(k)+d; ed=exp(d); dm=max(dm,abs(d));
      <i=1,n; z=e(i,k); e(i,k)=z*ed; s(i)=s(i)-z+e(i,k);>
   >
> until dm.lt.eps;
az=az-sum(az)/kk;
deallocate(e,s);
return;
end;
function elc(parm,n,cl,a,m);
real a(n),cl(2); integer m(n);
fn=n; am=sum(a)/fn;
if parm.eq.0.0.or.n.eq.2 < elc=am; go to :chk:;>
<i=1,n; m(i)=i;> call psort7(a,m,1,n);
if a(m(1)).eq.a(m(n)) < elc=a(1); go to :chk:;>
if mod(n,2).eq.1 < ad=a(m(n/2+1));>
else < ad=0.5*(a(m(n/2+1))+a(m(n/2)));>
if parm.eq.1.0 < elc=ad; go to :chk:;>
b1=min(am,ad); b2=max(am,ad); k2=1;
until a(m(k2)).gt.b1 < k2=k2+1;> k1=k2-1;
until a(m(k2)).ge.b2 < k2=k2+1;>
r=parm/((1.0-parm)*fn); is=0; sm=n-2*(k1-1);
<k=k1,k2-1; sm=sm-2.0; s=r*sm+am;
   if s.gt.a(m(k)).and.s.le.a(m(k+1)) < is=k; exit;>
>
if is.ne.0 < elc=s; go to :chk:;> r2=2.0*r; s1=a(m(k1)); am2=2.0*am;
cri=r2*sum(abs(a-s1))+s1*(s1-am2); elc=s1;
<k=k1+1,k2; s=a(m(k)); if(s.eq.s1) next;
   c=r2*sum(abs(a-s))+s*(s-am2);
   if c.lt.cri < cri=c; elc=s;> s1=s;
>
:chk:elc=max(maxval(a-cl(2)),min(minval(a-cl(1)),elc));
return;
end;
function nintot(ni,nx,nc,a,m,nin,is);
real a(nx,nc); integer m(nx),is(ni);
is=0; nintot=0;
<ic=1,nc; <j=1,nin; k=m(j); if(is(k).ne.0) next;
   if(a(j,ic).eq.0.0) next; is(k)=k; nintot=nintot+1;
>>   
return;
end;
subroutine luncomp(ni,nx,nc,ca,ia,nin,a);
real ca(nx,nc),a(ni,nc); integer ia(nx);
a=0.0; 
<ic=1,nc; if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic);>
return;
end;
subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans);
real a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt); integer ia(nx);
<i=1,nt; <ic=1,nc; ans(ic,i)=a0(ic);
   if(nin.gt.0)
      ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1:nin)));
>>
return;
end;
subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,cl,ne,nx,nlam,flmin,
   ulam,thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr);
real x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam);
real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni);
integer ix(*),jx(*),jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: xm,xs,ww,vq,xv
      integer, dimension (:), allocatable :: ju
%mortran
if maxval(vp).le.0.0 < jerr=10000; return;>
allocate(ww(1:no),stat=jerr);
allocate(ju(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(vq(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xs(1:ni),stat=ierr); jerr=jerr+ierr;
if kopt.eq.2 < allocate(xv(1:ni),stat=ierr); jerr=jerr+ierr;>
if(jerr.ne.0) return;
call spchkvars(no,ni,x,ix,ju);
if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0;
if maxval(ju).le.0 < jerr=7777; return;>
vq=max(0.0,vp); vq=vq*ni/sum(vq);
<i=1,no; ww(i)=sum(y(i,:)); if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i);>
sw=sum(ww); ww=ww/sw;
if nc.eq.1 < call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs);
   if isd.gt.0 < <j=1,ni; cl(:,j)=cl(:,j)*xs(j);>>
   call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,cl,ne,nx,nlam,
      flmin,ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,
      alm,nlp,jerr);
>
elseif kopt.eq.2 <
   call multsplstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs,xv);
   if isd.gt.0 < <j=1,ni; cl(:,j)=cl(:,j)*xs(j);>>
   call multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,
      ulam,thr,intr,maxit,xv,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr);
>
else < call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs);
   if isd.gt.0 < <j=1,ni; cl(:,j)=cl(:,j)*xs(j);>>
   call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,
      ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,
      ia,nin,dev0,dev,alm,nlp,jerr);
>
if(jerr.gt.0) return; dev0=2.0*sw*dev0;
<k=1,lmu; nk=nin(k);
   <ic=1,nc; if isd.gt.0 < <l=1,nk; ca(l,ic,k)=ca(l,ic,k)/xs(ia(l));>>
      if intr.eq.0 < a0(ic,k)=0.0;>
      else < a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)));>
   >
>
deallocate(ww,ju,vq,xm,xs); if(kopt.eq.2) deallocate(xv);
return;
end;
subroutine multsplstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs,xv);
real x(*),w(no),xm(ni),xs(ni),xv(ni); integer ix(*),jx(*),ju(ni);
if intr.eq.0 < 
   <j=1,ni; if(ju(j).eq.0) next; xm(j)=0.0; jb=ix(j); je=ix(j+1)-1;
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2);
      if isd.ne.0 < xbq=dot_product(w(jx(jb:je)),x(jb:je))**2; vc=xv(j)-xbq;
         xs(j)=sqrt(vc); xv(j)=1.0+xbq/vc;
      >
      else < xs(j)=1.0;>
   >
   return;
>
<j=1,ni; if(ju(j).eq.0) next; jb=ix(j); je=ix(j+1)-1;
   xm(j)=dot_product(w(jx(jb:je)),x(jb:je));
   xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2;
   if isd.gt.0 < xs(j)=sqrt(xv(j)); xv(j)=1.0;>
>
if(isd.eq.0) xs=1.0;
return;
end;
subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs);
real x(*),w(no),xm(ni),xs(ni); integer ix(*),jx(*),ju(ni);
if intr.eq.0 <
   <j=1,ni; if(ju(j).eq.0) next; xm(j)=0.0; jb=ix(j); je=ix(j+1)-1;
      if isd.ne.0 <
         vc=dot_product(w(jx(jb:je)),x(jb:je)**2)
            -dot_product(w(jx(jb:je)),x(jb:je))**2;
         xs(j)=sqrt(vc);
      >
      else < xs(j)=1.0;>
   >
   return;
>
<j=1,ni; if(ju(j).eq.0) next; jb=ix(j); je=ix(j+1)-1;
   xm(j)=dot_product(w(jx(jb:je)),x(jb:je));
   if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2);
>
if(isd.eq.0) xs=1.0;
return;
end;
subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,nlam,
   flmin,ulam,shri,isd,intr,maxit,kopt,xb,xs,
   lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr);
real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni);
real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam);
real xb(ni),xs(ni); integer ix(*),jx(*),ju(ni),m(nx),kin(nlam);
%fortran
      real, dimension (:), allocatable :: xm,b,bs,v,r,sc,xv,q,ga
      integer, dimension (:), allocatable :: mm,ixx
%mortran
call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx);
allocate(b(0:ni),stat=jerr);
allocate(xm(0:ni),stat=ierr); jerr=jerr+ierr;
allocate(xv(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(bs(0:ni),stat=ierr); jerr=jerr+ierr;
allocate(ga(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ixx(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(q(1:no),stat=ierr); jerr=jerr+ierr;
allocate(r(1:no),stat=ierr); jerr=jerr+ierr;
allocate(v(1:no),stat=ierr); jerr=jerr+ierr;
allocate(sc(1:no),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
fmax=log(1.0/pmin-1.0); fmin=-fmax; vmin=(1.0+pmin)*pmin*(1.0-pmin);
bta=parm; omb=1.0-bta;
q0=dot_product(w,y); if q0.le.pmin < jerr=8001; return;>
if q0.ge.1.0-pmin < jerr=9001; return;>
if(intr.eq.0) q0=0.5; bz=0.0; if(intr.ne.0) bz=log(q0/(1.0-q0));
if nonzero(no,g).eq.0 < vi=q0*(1.0-q0); b(0)=bz; v=vi*w;
   r=w*(y-q0); q=q0; xm(0)=vi; dev1=-(bz*q0+log(1.0-q0));
>
else < b(0)=0.0;
   if intr.ne.0 < b(0)=azero(no,y,g,w,jerr); if(jerr.ne.0) return;>
   q=1.0/(1.0+exp(-b(0)-g)); v=w*q*(1.0-q); r=w*(y-q); xm(0)=sum(v);
   dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)));
>
if kopt.gt.0 <
   if isd.gt.0.and.intr.ne.0 < xv=0.25;>
   else <
      <j=1,ni; if(ju(j).eq.0) next; jb=ix(j); je=ix(j+1)-1;
         xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2);
      >
   >
>
b(1:ni)=0.0; dev0=dev1;
<i=1,no; if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i));
   if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i));
>
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
m=0; mm=0; nin=0; /o,svr/=0.0; mnl=min(mnlam,nlam); bs=0.0; /nlp,nin/=0;
shr=shri*dev0; al=0.0; ixx=0;
<j=1,ni; if(ju(j).eq.0) next;
   jb=ix(j); je=ix(j+1)-1; jn=ix(j+1)-ix(j);
   sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o;
   gj=dot_product(sc(1:jn),x(jb:je));
   ga(j)=abs((gj-svr*xb(j))/xs(j));
>
<ilm=1,nlam; al0=al;
   if flmin.ge.1.0 < al=ulam(ilm);>
   elseif ilm.gt.2 < al=al*alf;>
   elseif ilm.eq.1 < al=big;>
   else < al0=0.0;
      <j=1,ni; if(ju(j).eq.0) next; if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j));>
      al0=al0/max(bta,1.0e-3); al=alf*al0;
   >
   al2=al*omb; al1=al*bta; tlam=bta*(2.0*al-al0);
   <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1;
   >
   :again:continue;   
   loop <  bs(0)=b(0); if(nin.gt.0) bs(m(1:nin))=b(m(1:nin));
      <j=1,ni; if(ixx(j).eq.0) next;
         jb=ix(j); je=ix(j+1)-1; jn=ix(j+1)-ix(j);
         sc(1:jn)=v(jx(jb:je));
         xm(j)=dot_product(sc(1:jn),x(jb:je));
         if kopt.eq.0 <
            xv(j)=dot_product(sc(1:jn),x(jb:je)**2);
            xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2;
         >
      >
      loop < nlp=nlp+1; dlx=0.0;
         <k=1,ni; if(ixx(k).eq.0) next;
            jb=ix(k); je=ix(k+1)-1; jn=ix(k+1)-ix(k); bk=b(k);
            sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o;
            gk=dot_product(sc(1:jn),x(jb:je));
            gk=(gk-svr*xb(k))/xs(k);
            u=gk+xv(k)*b(k); au=abs(u)-vp(k)*al1;
            if au.le.0.0 < b(k)=0.0;>
            else <
               b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)));
            >
            d=b(k)-bk; if(abs(d).le.0.0) next; dlx=max(dlx,xv(k)*d**2);
            if mm(k).eq.0 < nin=nin+1; if(nin.gt.nx) exit;
               mm(k)=nin; m(nin)=k; sc(1:jn)=v(jx(jb:je));
               xm(k)=dot_product(sc(1:jn),x(jb:je));
            >
            r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k);
            o=o+d*(xb(k)/xs(k));
            svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k);
         >
         if(nin.gt.nx) exit;
         d=0.0; if(intr.ne.0) d=svr/xm(0);
         if d.ne.0.0 < b(0)=b(0)+d; dlx=max(dlx,xm(0)*d**2); r=r-d*v;
            svr=svr-d*xm(0);
         >
         if(dlx.lt.shr) exit;
         if nlp.gt.maxit < jerr=-ilm; return;>
         loop < nlp=nlp+1; dlx=0.0;
            <l=1,nin; k=m(l); jb=ix(k); je=ix(k+1)-1;
               jn=ix(k+1)-ix(k); bk=b(k);
               sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o;
               gk=dot_product(sc(1:jn),x(jb:je));
               gk=(gk-svr*xb(k))/xs(k);
               u=gk+xv(k)*b(k); au=abs(u)-vp(k)*al1;
               if au.le.0.0 < b(k)=0.0;>
               else <
                  b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)));
               >
               d=b(k)-bk; if(abs(d).le.0.0) next; dlx=max(dlx,xv(k)*d**2);
               r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k);
               o=o+d*(xb(k)/xs(k));
               svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k);
            >
            d=0.0; if(intr.ne.0) d=svr/xm(0);
            if d.ne.0.0 < b(0)=b(0)+d; dlx=max(dlx,xm(0)*d**2); r=r-d*v;
               svr=svr-d*xm(0);
            >
            if(dlx.lt.shr) exit;
            if nlp.gt.maxit < jerr=-ilm; return;>
         >
      >
      if(nin.gt.nx) exit;
      sc=b(0); b0=0.0;
      <j=1,nin; l=m(j); jb=ix(l); je=ix(l+1)-1;
         sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l);
         b0=b0-b(l)*xb(l)/xs(l);
      >
      sc=sc+b0;
      <i=1,no; fi=sc(i)+g(i);
         if fi.lt.fmin < q(i)=0.0;> elseif fi.gt.fmax < q(i)=1.0;>
         else < q(i)=1.0/(1.0+exp(-fi));>
      >
      v=w*q*(1.0-q); xm(0)=sum(v); if(xm(0).lt.vmin) exit;
      r=w*(y-q); svr=sum(r); o=0.0;
      if xm(0)*(b(0)-bs(0))**2.lt.shr < kx=0;
         <j=1,nin; k=m(j);
            if(xv(k)*(b(k)-bs(k))**2.lt.shr) next; kx=1; exit;
         >
         if kx.eq.0 <
            <j=1,ni; if(ixx(j).eq.1) next; if(ju(j).eq.0) next;
               jb=ix(j); je=ix(j+1)-1; jn=ix(j+1)-ix(j);
               sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o;
               gj=dot_product(sc(1:jn),x(jb:je));
               ga(j)=abs((gj-svr*xb(j))/xs(j));
               if ga(j).gt.al1*vp(j) < ixx(j)=1; kx=1;>
            >
            if(kx.eq.1) go to :again:;     
            exit;
         >
      >
   >
   if nin.gt.nx < jerr=-10000-ilm;  exit;>
   if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin)); kin(ilm)=nin;
   a0(ilm)=b(0); alm(ilm)=al; lmu=ilm;
   devi=dev2(no,w,y,q,pmin);
   dev(ilm)=(dev1-devi)/dev0;
   if(ilm.lt.mnl) next; if(flmin.ge.1.0) next;
   me=0; <j=1,nin; if(a(j,ilm).ne.0.0) me=me+1;> if(me.gt.ne) exit;
   if(dev(ilm).gt.devmax) exit; if(dev(ilm)-dev(ilm-1).lt.sml) exit;
   if(xm(0).lt.vmin) exit;
>
g=log(q/(1.0-q));
deallocate(xm,b,bs,v,r,sc,xv,q,mm,ga,ixx);
return;
end;
subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,nlam,flmin,
   ulam,shri,isd,intr,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr);
real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni);
real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni);
integer ix(*),jx(*),ju(ni),m(nx),kin(nlam);
%fortran
      real, dimension (:,:), allocatable :: q
      real, dimension (:), allocatable :: sxp,sxpl
      real, dimension (:), allocatable :: sc,xm,v,r,ga
      real, dimension (:,:), allocatable :: b,bs,xv
      integer, dimension (:), allocatable :: mm,is,iy
      allocate(b(0:ni,1:nc),stat=jerr)
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr      
%mortran
call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx); exmn=-exmx;
allocate(xm(0:ni),stat=ierr); jerr=jerr+ierr;
allocate(r(1:no),stat=ierr); jerr=jerr+ierr;
allocate(v(1:no),stat=ierr); jerr=jerr+ierr;
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ga(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(iy(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(is(1:max(nc,ni)),stat=ierr); jerr=jerr+ierr;
allocate(sxp(1:no),stat=ierr); jerr=jerr+ierr;
allocate(sxpl(1:no),stat=ierr); jerr=jerr+ierr;
allocate(sc(1:no),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
pmax=1.0-pmin; emin=pmin/pmax; emax=1.0/emin;
pfm=(1.0+pmin)*pmin; pfx=(1.0-pmin)*pmax; vmin=pfm*pmax;
bta=parm; omb=1.0-bta; dev1=0.0; dev0=0.0;
<ic=1,nc; q0=dot_product(w,y(:,ic));
   if q0.le.pmin < jerr =8000+ic; return;>
   if q0.ge.1.0-pmin < jerr =9000+ic; return;>
   if(intr.eq.0) q0=1.0/nc;
   b(1:ni,ic)=0.0; b(0,ic)=0.0;
   if intr.ne.0 < b(0,ic)=log(q0); dev1=dev1-q0*b(0,ic);>
>
if(intr.eq.0) dev1=log(float(nc)); iy=0; al=0.0;
if nonzero(no*nc,g).eq.0 <
   b(0,:)=b(0,:)-sum(b(0,:))/nc; sxp=0.0;
   <ic=1,nc; q(:,ic)=exp(b(0,ic)); sxp=sxp+q(:,ic);>
>
else < <i=1,no; g(i,:)=g(i,:)-sum(g(i,:))/nc;> sxp=0.0;
   if intr.eq.0 < b(0,:)=0.0;>
   else < call kazero(nc,no,y,g,w,b(0,:),jerr); if(jerr.ne.0) return;>
   dev1=0.0;
   <ic=1,nc; q(:,ic)=b(0,ic)+g(:,ic);
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic));
      q(:,ic)=exp(q(:,ic)); sxp=sxp+q(:,ic);
   >
   sxpl=w*log(sxp); <ic=1,nc; dev1=dev1+dot_product(y(:,ic),sxpl);>
>
<ic=1,nc; <i=1,no; if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic));>>
dev0=dev0+dev1;
if kopt.gt.0 <
   if isd.gt.0.and.intr.ne.0 < xv=0.25;>
   else <
      <j=1,ni; if(ju(j).eq.0) next; jb=ix(j); je=ix(j+1)-1;
         xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2);
      >
   >
>
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
m=0; mm=0; nin=0; nlp=0; mnl=min(mnlam,nlam); bs=0.0; svr=0.0; o=0.0;
shr=shri*dev0; ga=0.0;
<ic=1,nc; v=q(:,ic)/sxp; r=w*(y(:,ic)-v); v=w*v*(1.0-v);
   <j=1,ni; if(ju(j).eq.0) next;
      jb=ix(j); je=ix(j+1)-1; jn=ix(j+1)-ix(j);
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je));
      gj=dot_product(sc(1:jn),x(jb:je));
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j));
   >
>
<ilm=1,nlam; al0=al;
   if flmin.ge.1.0 < al=ulam(ilm);>
   elseif ilm.gt.2 < al=al*alf;>
   elseif ilm.eq.1 < al=big;>
   else < al0=0.0;     
      <j=1,ni; if(ju(j).eq.0) next; if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j));>
      al0=al0/max(bta,1.0e-3); al=alf*al0;
   >
   al2=al*omb; al1=al*bta; tlam=bta*(2.0*al-al0);
   <k=1,ni; if(iy(k).eq.1) next; if(ju(k).eq.0) next;
      if(ga(k).gt.tlam*vp(k)) iy(k)=1;
   >
   :again:continue;
   loop < /ixx,jxx/=0; ig=0;
      <ic=1,nc; bs(0,ic)=b(0,ic);
         if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic);
         xm(0)=0.0; svr=0.0; o=0.0;
         <i=1,no; pic=q(i,ic)/sxp(i);
            if pic.lt.pfm <pic=0.0; v(i)=0.0;>
            elseif pic.gt.pfx < pic=1.0; v(i)=0.0;>
            else < v(i)=w(i)*pic*(1.0-pic); xm(0)=xm(0)+v(i);>
            r(i)=w(i)*(y(i,ic)-pic); svr=svr+r(i);
         >
         if(xm(0).le.vmin) next; ig=1;
         <j=1,ni; if(iy(j).eq.0) next;
            jb=ix(j); je=ix(j+1)-1;
            xm(j)=dot_product(v(jx(jb:je)),x(jb:je));
            if kopt.eq.0 <
               xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2);
               xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2;
            >
         >
         loop < nlp=nlp+1; dlx=0.0;
            <k=1,ni; if(iy(k).eq.0) next;
               jb=ix(k); je=ix(k+1)-1; jn=ix(k+1)-ix(k); bk=b(k,ic);
               sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je));
               gk=dot_product(sc(1:jn),x(jb:je));
               gk=(gk-svr*xb(k))/xs(k);
               u=gk+xv(k,ic)*b(k,ic); au=abs(u)-vp(k)*al1;
               if au.le.0.0 < b(k,ic)=0.0;>
               else <
                  b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/
                  (xv(k,ic)+vp(k)*al2)));
               >
               d=b(k,ic)-bk; if(abs(d).le.0.0) next;
               dlx=max(dlx,xv(k,ic)*d**2);
               if mm(k).eq.0 < nin=nin+1;
                  if nin.gt.nx < jxx=1; exit;>
                  mm(k)=nin; m(nin)=k;
                  xm(k)=dot_product(v(jx(jb:je)),x(jb:je));
               >
               r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k);
               o=o+d*(xb(k)/xs(k));
               svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k);
            >
            if(jxx.gt.0) exit;
            d=0.0; if(intr.ne.0) d=svr/xm(0);
            if d.ne.0.0 < b(0,ic)=b(0,ic)+d; dlx=max(dlx,xm(0)*d**2);
               r=r-d*v; svr=svr-d*xm(0);
            >
            if(dlx.lt.shr) exit; if nlp.gt.maxit < jerr=-ilm; return;>
            loop < nlp=nlp+1; dlx=0.0;
               <l=1,nin; k=m(l); jb=ix(k); je=ix(k+1)-1;
                  jn=ix(k+1)-ix(k); bk=b(k,ic);
                  sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je));
                  gk=dot_product(sc(1:jn),x(jb:je));
                  gk=(gk-svr*xb(k))/xs(k);
                  u=gk+xv(k,ic)*b(k,ic); au=abs(u)-vp(k)*al1;
                  if au.le.0.0 < b(k,ic)=0.0;>
                  else <
                     b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/
                     (xv(k,ic)+vp(k)*al2)));
                  >
                  d=b(k,ic)-bk; if(abs(d).le.0.0) next;
                  dlx=max(dlx,xv(k,ic)*d**2);
                  r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k);
                  o=o+d*(xb(k)/xs(k));
                  svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k);
               >
               d=0.0; if(intr.ne.0) d=svr/xm(0);
               if d.ne.0.0 < b(0,ic)=b(0,ic)+d; dlx=max(dlx,xm(0)*d**2);
                  r=r-d*v; svr=svr-d*xm(0);
               >
               if(dlx.lt.shr) exit; if nlp.gt.maxit < jerr=-ilm; return;>                
            >
         >
         if(jxx.gt.0) exit;
         if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1;
         if ixx.eq.0 <
            <j=1,nin; k=m(j);
               if xv(k,ic)*(b(k,ic)-bs(k,ic))**2.gt.shr < ixx=1; exit;>               
            >
         >
         sc=b(0,ic)+g(:,ic); b0=0.0;
         <j=1,nin; l=m(j); jb=ix(l); je=ix(l+1)-1;
            sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l);
            b0=b0-b(l,ic)*xb(l)/xs(l);
         >
         sc=min(max(exmn,sc+b0),exmx);     
         sxp=sxp-q(:,ic);
         q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp);
         sxp=sxp+q(:,ic);
      >
      s=-sum(b(0,:))/nc; b(0,:)=b(0,:)+s; sc=s; b0=0.0;
      <j=1,nin; l=m(j);
         if vp(l).le.0.0 < s=sum(b(l,:))/nc;>
         else < s=elc(parm,nc,cl(:,l),b(l,:),is);>
         b(l,:)=b(l,:)-s;
         jb=ix(l); je=ix(l+1)-1;
         sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l);
         b0=b0+s*xb(l)/xs(l);
      >
      sc=sc+b0; sc=exp(sc); sxp=sxp*sc;  <ic=1,nc; q(:,ic)=q(:,ic)*sc;>
      if(jxx.gt.0) exit; if(ig.eq.0) exit;
      if ixx.eq.0 <
         <j=1,ni; if(iy(j).eq.1) next; if(ju(j).eq.0) next; ga(j)=0.0;>
         <ic=1,nc; v=q(:,ic)/sxp; r=w*(y(:,ic)-v); v=w*v*(1.0-v);
            <j=1,ni; if(iy(j).eq.1) next; if(ju(j).eq.0) next;
               jb=ix(j); je=ix(j+1)-1; jn=ix(j+1)-ix(j);
               sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je));
               gj=dot_product(sc(1:jn),x(jb:je));
               ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j));
            >
         >
         <k=1,ni; if(iy(k).eq.1) next; if(ju(k).eq.0) next;
            if ga(k).gt.al1*vp(k) < iy(k)=1; ixx=1;>
         >
         if(ixx.eq.1) go to :again:;
         exit;
      >
   > 
   if jxx.gt.0 < jerr=-10000-ilm; exit;> devi=0.0;
   <ic=1,nc;
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic); a0(ic,ilm)=b(0,ic);
      <i=1,no; if (y(i,ic).le.0.0) next;
         devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i));
      >
   >
   kin(ilm)=nin; alm(ilm)=al; lmu=ilm;
   dev(ilm)=(dev1-devi)/dev0; if(ig.eq.0) exit;
   if(ilm.lt.mnl) next; if(flmin.ge.1.0) next;
   if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne) exit;
   if(dev(ilm).gt.devmax) exit; if(dev(ilm)-dev(ilm-1).lt.sml) exit;
>
g=log(q); <i=1,no; g(i,:)=g(i,:)-sum(g(i,:))/nc;>
deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy);
return;
end;
subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f);
real a0(nc),ca(nx,nc),x(*),f(nc,n); integer ia(*),ix(*),jx(*);
<ic=1,nc; f(ic,:)=a0(ic);>
<j=1,nin; k=ia(j); kb=ix(k); ke=ix(k+1)-1;
   <ic=1,nc; f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke);>
>
return;
end;
subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
   maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr);
real x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam);
real ca(nx,nlam),dev(nlam),alm(nlam),cl(2,ni);
integer jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: xs,ww,vq
      integer, dimension (:), allocatable :: ju
%mortran
if maxval(vp).le.0.0 < jerr=10000; return;>
allocate(ww(1:no),stat=jerr);
allocate(ju(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(vq(1:ni),stat=ierr); jerr=jerr+ierr;
if isd.gt.0 < allocate(xs(1:ni),stat=ierr); jerr=jerr+ierr;>
if(jerr.ne.0) return;
call chkvars(no,ni,x,ju);
if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0;
if maxval(ju).le.0 < jerr=7777; return;>
vq=max(0.0,vp); vq=vq*ni/sum(vq);
ww=max(0.0,w); sw=sum(ww);
if sw.le.0.0 < jerr=9999; return;> ww=ww/sw;
call cstandard(no,ni,x,ww,ju,isd,xs);
if isd.gt.0 < <j=1,ni; cl(:,j)=cl(:,j)*xs(j);>>
call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,thr,
         isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr);
if(jerr.gt.0) return; dev0=2.0*sw*dev0;
if isd.gt.0 < <k=1,lmu; nk=nin(k); ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk));>>
deallocate(ww,ju,vq); if(isd.gt.0) deallocate(xs);
return;
end;
subroutine cstandard (no,ni,x,w,ju,isd,xs);
real x(no,ni),w(no),xs(ni); integer ju(ni);
<j=1,ni; if(ju(j).eq.0) next;
   xm=dot_product(w,x(:,j)); x(:,j)=x(:,j)-xm;
   if isd.gt.0 < xs(j)=sqrt(dot_product(w,x(:,j)**2)); x(:,j)=x(:,j)/xs(j);>
>
return;
end;
subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,cl,ne,nx,nlam,flmin,ulam,cthri,
    isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr);
real x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam);
real ao(nx,nlam),dev(nlam),alm(nlam),cl(2,ni);
integer ju(ni),m(nx),kin(nlam);
%fortran
      real, dimension (:), allocatable :: w,dk,v,xs,wr,a,as,f,dq
      real, dimension (:), allocatable :: e,uu,ga
      integer, dimension (:), allocatable :: jp,kp,mm,ixx
%mortran
call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx);
sml=sml*100.0; devmax=devmax*0.99/0.999;
allocate(e(1:no),stat=jerr);
allocate(uu(1:no),stat=ierr); jerr=jerr+ierr;
allocate(f(1:no),stat=ierr); jerr=jerr+ierr;
allocate(w(1:no),stat=ierr); jerr=jerr+ierr;
allocate(v(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(a(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(as(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xs(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ga(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ixx(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(jp(1:no),stat=ierr); jerr=jerr+ierr;
allocate(kp(1:no),stat=ierr); jerr=jerr+ierr;
allocate(dk(1:no),stat=ierr); jerr=jerr+ierr;
allocate(wr(1:no),stat=ierr); jerr=jerr+ierr;
allocate(dq(1:no),stat=ierr); jerr=jerr+ierr;
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0)go to :done:;
call groups(no,y,d,q,nk,kp,jp,t0,jerr);
if(jerr.ne.0) go to :done:; alpha=parm;
oma=1.0-alpha; nlm=0; ixx=0; al=0.0;
dq=d*q; call died(no,nk,dq,kp,jp,dk);
a=0.0; f(1)=0.0; fmax=log(huge(f(1))*0.1);
if nonzero(no,g).ne.0 < f=g-dot_product(q,g);
   e=q*exp(sign(min(abs(f),fmax),f));
> 
else < f=0.0; e=q;>
r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu);
rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0); dev0=rr;
<i=1,no; if y(i).lt.t0.or.q(i).le.0.0 < /w(i),wr(i)/=0.0;>>  
call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu);
if(jerr.ne.0) go to :done:;
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
m=0; mm=0; /nlp,nin/=0; mnl=min(mnlam,nlam); as=0.0; cthr=cthri*dev0;
<j=1,ni; if(ju(j).eq.0) next; ga(j)=abs(dot_product(wr,x(:,j)));>
<ilm=1,nlam; al0=al;
   if flmin.ge.1.0 < al=ulam(ilm);>
   elseif ilm.gt.2 < al=al*alf;>
   elseif ilm.eq.1 < al=big;>
   else < al0=0.0;
      <j=1,ni; if(ju(j).eq.0) next; if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j));>
      al0=al0/max(parm,1.0e-3); al=alf*al0;
   >
   sa=alpha*al; omal=oma*al; tlam=alpha*(2.0*al-al0);
   <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1;
   >
   :again:continue;
   loop < if(nin.gt.0) as(m(1:nin))=a(m(1:nin));
      call vars(no,ni,x,w,ixx,v);
      loop < nlp=nlp+1; dli=0.0;
         <j=1,ni; if(ixx(j).eq.0) next;
            u=a(j)*v(j)+dot_product(wr,x(:,j));
            if abs(u).le.vp(j)*sa < at=0.0;>
            else < at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/
               (v(j)+vp(j)*omal)));
            >
            if at.ne.a(j) < del=at-a(j); a(j)=at; dli=max(dli,v(j)*del**2);
               wr=wr-del*w*x(:,j); f=f+del*x(:,j);
               if mm(j).eq.0 < nin=nin+1; if(nin.gt.nx) exit;
                  mm(j)=nin; m(nin)=j;
               >
            >
         >
         if(nin.gt.nx) exit; if(dli.lt.cthr) exit;
         if nlp.gt.maxit < jerr=-ilm; return;>
         loop < nlp=nlp+1;  dli=0.0;
            <l=1,nin; j=m(l);
               u=a(j)*v(j)+dot_product(wr,x(:,j));
               if abs(u).le.vp(j)*sa < at=0.0;>
               else < at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/
                  (v(j)+vp(j)*omal)));
               >
               if at.ne.a(j) < del=at-a(j); a(j)=at; dli=max(dli,v(j)*del**2);
                  wr=wr-del*w*x(:,j); f=f+del*x(:,j);
               >
            >
            if(dli.lt.cthr) exit; if nlp.gt.maxit < jerr=-ilm; return;>
         >
      >
      if(nin.gt.nx) exit;
      e=q*exp(sign(min(abs(f),fmax),f));
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu);
      if jerr.ne.0 < jerr=jerr-ilm;  go to :done:;>
      ix=0;
      <j=1,nin; k=m(j);
         if(v(k)*(a(k)-as(k))**2.lt.cthr) next; ix=1; exit;>
      if ix.eq.0 <
         <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
            ga(k)=abs(dot_product(wr,x(:,k)));
            if ga(k).gt.sa*vp(k) < ixx(k)=1; ix=1;>
         >
         if(ix.eq.1) go to :again:;
         exit;
      >
   >
   if nin.gt.nx < jerr=-10000-ilm;  exit;>
   if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin)); kin(ilm)=nin;
   alm(ilm)=al; lmu=ilm;
   dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr;
   if(ilm.lt.mnl) next; if(flmin.ge.1.0) next;
   me=0; <j=1,nin; if(ao(j,ilm).ne.0.0) me=me+1;> if(me.gt.ne) exit;
   if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml) exit;
   if(dev(ilm).gt.devmax) exit;
>
g=f;
:done: deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx);
return;
end;
subroutine cxmodval(ca,ia,nin,n,x,f);
real ca(nin),x(n,*),f(n); integer ia(nin);
f=0.0; if(nin.le.0) return;
<i=1,n; f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)));>
return;
end;
subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr);
real y(no),d(no),q(no); integer jp(no),kp(*"nk");
<j=1,no; jp(j)=j;> call psort7(y,jp,1,no);
nj=0; <j=1,no; if(q(jp(j)).le.0.0) next; nj=nj+1; jp(nj)=jp(j);>
if nj.eq.0 < jerr=20000; return;>
j=1; until d(jp(j)).gt.0.0 < j=j+1;> until j.gt.nj;
if j.ge.nj-1 < jerr=30000; return;>
t0=y(jp(j)); j0=j-1;
if j0.gt.0 <
   until y(jp(j0)).lt.t0 < j0=j0-1;> until j0.eq.0;
   if j0.gt.0 < nj=nj-j0; <j=1,nj; jp(j)=jp(j+j0);>>
>
jerr=0; nk=0; yk=t0; j=2;
loop <
   until d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk < j=j+1;> until j.gt.nj;
   nk=nk+1; kp(nk)=j-1; if(j.gt.nj) exit;
   if j.eq.nj < nk=nk+1; kp(nk)=nj; exit;>
   yk=y(jp(j)); j=j+1;
>
return;
end;
subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u);
real d(no),dk(nk),wr(no),w(no);
real e(no),u(no),b,c; integer kp(nk),jp(no);
call usk(no,nk,kp,jp,e,u);
b=dk(1)/u(1); c=dk(1)/u(1)**2; jerr=0;
<j=1,kp(1); i=jp(j);
   w(i)=e(i)*(b-e(i)*c); if w(i).le.0.0 < jerr=-30000; return;>
   wr(i)=d(i)-e(i)*b;
>
<k=2,nk; j1=kp(k-1)+1; j2=kp(k);
   b=b+dk(k)/u(k); c=c+dk(k)/u(k)**2;
   <j=j1,j2; i=jp(j);
      w(i)=e(i)*(b-e(i)*c); if w(i).le.0.0 < jerr=-30000; return;>
      wr(i)=d(i)-e(i)*b;
   >
>
return;
end;
subroutine vars(no,ni,x,w,ixx,v);
real x(no,ni),w(no),v(ni); integer ixx(ni);
<j=1,ni; if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2);>
return;
end;
subroutine died(no,nk,d,kp,jp,dk);
real d(no),dk(nk); integer kp(nk),jp(no);
dk(1)=sum(d(jp(1:kp(1))));
<k=2,nk; dk(k)=sum(d(jp((kp(k-1)+1):kp(k))));>
return;
end;
subroutine usk(no,nk,kp,jp,e,u);
real e(no),u(nk),h; integer kp(nk),jp(no);
h=0.0;
<k=nk,1,-1; j2=kp(k);
   j1=1; if(k.gt.1) j1=kp(k-1)+1;
   <j=j2,j1,-1; h=h+e(jp(j));>
   u(k)=h;
>
return;
end;
function risk(no,ni,nk,d,dk,f,e,kp,jp,u);
real d(no),dk(nk),f(no);
integer kp(nk),jp(no); real e(no),u(nk),s;
call usk(no,nk,kp,jp,e,u); u=log(u);
risk=dot_product(d,f)-dot_product(dk,u);
return;
end;
subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr);
real x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(nlam);
%fortran
      real, dimension (:), allocatable :: dk,f,xm,dq,q
      real, dimension (:), allocatable :: e,uu
      integer, dimension (:), allocatable :: jp,kp
%mortran
allocate(e(1:no),stat=jerr);
allocate(q(1:no),stat=ierr); jerr=jerr+ierr;
allocate(uu(1:no),stat=ierr); jerr=jerr+ierr;
allocate(f(1:no),stat=ierr); jerr=jerr+ierr;
allocate(dk(1:no),stat=ierr); jerr=jerr+ierr;
allocate(jp(1:no),stat=ierr); jerr=jerr+ierr;
allocate(kp(1:no),stat=ierr); jerr=jerr+ierr;
allocate(dq(1:no),stat=ierr); jerr=jerr+ierr;
allocate(xm(1:ni),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) go to :done:;
q=max(0.0,w); sw=sum(q);
if sw.le.0.0 < jerr=9999; go to :done:;>
call groups(no,y,d,q,nk,kp,jp,t0,jerr);
if(jerr.ne.0) go to :done:; fmax=log(huge(e(1))*0.1);
dq=d*q; call died(no,nk,dq,kp,jp,dk); gm=dot_product(q,g)/sw;
<j=1,ni; xm(j)=dot_product(q,x(:,j))/sw;>
<lam=1,nlam;
   <i=1,no; f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm));
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)));
   >
   flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu);
>
:done: deallocate(e,uu,dk,f,jp,kp,dq);
return;
end;
subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
   isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr);
real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam);
real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni);
integer jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: xm,xs,ww,vq
      integer, dimension (:), allocatable :: ju
%mortran
if maxval(vp).le.0.0 < jerr=10000; return;>
if minval(y).lt.0.0 < jerr=8888; return;>
allocate(ww(1:no),stat=jerr);
allocate(ju(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(vq(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xm(1:ni),stat=ierr); jerr=jerr+ierr;
if isd.gt.0 < allocate(xs(1:ni),stat=ierr); jerr=jerr+ierr;>
if(jerr.ne.0) return;
call chkvars(no,ni,x,ju);
if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0;
if maxval(ju).le.0 < jerr=7777; go to :done:;>
vq=max(0.0,vp); vq=vq*ni/sum(vq);
ww=max(0.0,w); sw=sum(ww); if sw.le.0.0 < jerr=9999; go to :done:;>
ww=ww/sw;
call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs);
if isd.gt.0 < <j=1,ni; cl(:,j)=cl(:,j)*xs(j);>>
call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,thr,
      isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr);
if(jerr.gt.0) go to :done:; dev0=2.0*sw*dev0;
<k=1,lmu; nk=nin(k);
   if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk));
   if intr.eq.0 < a0(k)=0.0;>
   else < a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)));>
>
:done:deallocate(ww,ju,vq,xm); if(isd.gt.0) deallocate(xs);
return;
end;
subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,cl,ne,nx,nlam,flmin,ulam,shri,
    isd,intr,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr);
real x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam);
real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni);
integer ju(ni),m(nx),kin(nlam);
%fortran
      real, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga
      integer, dimension (:), allocatable :: mm,ixx
%mortran
call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx); sml=sml*10.0;
allocate(a(1:ni),stat=jerr);
allocate(as(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(t(1:no),stat=ierr); jerr=jerr+ierr;
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ga(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ixx(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(wr(1:no),stat=ierr); jerr=jerr+ierr;
allocate(v(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(w(1:no),stat=ierr); jerr=jerr+ierr;
allocate(f(1:no),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
bta=parm; omb=1.0-bta;
t=q*y; yb=sum(t); fmax=log(huge(bta)*0.1);
if nonzero(no,g).eq.0 <
   if intr.ne.0 < w=q*yb;  az=log(yb); f=az; dv0=yb*(az-1.0);>
   else < w=q; az=0.0; f=az; dv0=-1.0;>
>
else < w=q*exp(sign(min(abs(g),fmax),g)); v0=sum(w);
   if intr.ne.0 < eaz=yb/v0; w=eaz*w; az=log(eaz); f=az+g;
      dv0=dot_product(t,g)-yb*(1.0-az);
   >
   else < az=0.0; f=g; dv0=dot_product(t,g)-v0;>
>
a=0.0; as=0.0; wr=t-w; v0=1.0; if(intr.ne.0) v0=yb; dvr=-yb;
<i=1,no; if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i));> dvr=dvr-dv0; dev0=dvr;
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
m=0; mm=0; /nlp,nin/=0; mnl=min(mnlam,nlam); shr=shri*dev0; ixx=0; al=0.0;
<j=1,ni; if(ju(j).eq.0) next; ga(j)=abs(dot_product(wr,x(:,j)));>
<ilm=1,nlam; al0=al;
   if flmin.ge.1.0 < al=ulam(ilm);>
   elseif ilm.gt.2 < al=al*alf;>
   elseif ilm.eq.1 < al=big;>
   else < al0=0.0;
      <j=1,ni; if(ju(j).eq.0) next; if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j));>
      al0=al0/max(bta,1.0e-3); al=alf*al0;
   >
   al2=al*omb; al1=al*bta; tlam=bta*(2.0*al-al0);
   <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1;
   >
   :again:continue;
   loop <
      az0=az; if(nin.gt.0) as(m(1:nin))=a(m(1:nin));
      <j=1,ni; if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2);>
      loop < nlp=nlp+1; dlx=0.0;
         <k=1,ni; if(ixx(k).eq.0) next; ak=a(k);           
            u=dot_product(wr,x(:,k))+v(k)*ak; au=abs(u)-vp(k)*al1;
            if au.le.0.0 < a(k)=0.0;>
            else <
               a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)));
            >
            if(a(k).eq.ak) next; d=a(k)-ak; dlx=max(dlx,v(k)*d**2);
            wr=wr-d*w*x(:,k); f=f+d*x(:,k);
            if mm(k).eq.0 < nin=nin+1; if(nin.gt.nx) exit;
               mm(k)=nin; m(nin)=k;
            >
         >
         if(nin.gt.nx) exit;
         if intr.ne.0 < d=sum(wr)/v0;
            az=az+d; dlx=max(dlx,v0*d**2); wr=wr-d*w; f=f+d;
         >
         if(dlx.lt.shr) exit; if nlp.gt.maxit < jerr=-ilm; return;>
         loop < nlp=nlp+1; dlx=0.0;
            <l=1,nin; k=m(l); ak=a(k);
               u=dot_product(wr,x(:,k))+v(k)*ak; au=abs(u)-vp(k)*al1;
               if au.le.0.0 < a(k)=0.0;>
               else <
                  a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)));
               >
               if(a(k).eq.ak) next; d=a(k)-ak; dlx=max(dlx,v(k)*d**2);
               wr=wr-d*w*x(:,k); f=f+d*x(:,k)
            >
            if intr.ne.0 < d=sum(wr)/v0; az=az+d;
               dlx=max(dlx,v0*d**2); wr=wr-d*w; f=f+d;
            >
            if(dlx.lt.shr) exit; if nlp.gt.maxit < jerr=-ilm; return;>
         >
      >
      if(nin.gt.nx) exit;
      w=q*exp(sign(min(abs(f),fmax),f)); v0=sum(w); wr=t-w;
      if v0*(az-az0)**2.lt.shr < ix=0;
         <j=1,nin; k=m(j);
            if(v(k)*(a(k)-as(k))**2.lt.shr) next; ix=1; exit;
         >
         if ix.eq.0 <
            <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
               ga(k)=abs(dot_product(wr,x(:,k)));
               if ga(k).gt.al1*vp(k) < ixx(k)=1; ix=1;>
            >
            if(ix.eq.1) go to :again:;
            exit;
         >
      >    
   >
   if nin.gt.nx < jerr=-10000-ilm;  exit;>
   if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin)); kin(ilm)=nin;
   a0(ilm)=az; alm(ilm)=al; lmu=ilm;
   dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr;
   if(ilm.lt.mnl) next; if(flmin.ge.1.0) next;
   me=0; <j=1,nin; if(ca(j,ilm).ne.0.0) me=me+1;> if(me.gt.ne) exit;
   if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml) exit;
   if(dev(ilm).gt.devmax) exit;
>
g=f;
:done:deallocate(t,w,wr,v,a,f,as,mm,ga,ixx);
return;
end;
function nonzero(n,v);
real v(n);
nonzero=0; <i=1,n; if v(i).ne.0.0 < nonzero=1; return;>>
return;
end;
subroutine solns(ni,nx,lmu,a,ia,nin,b);
real a(nx,lmu),b(ni,lmu); integer ia(nx),nin(lmu);
<lam=1,lmu; call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam));>
return;
end;
subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b);
real a(nx,nc,lmu),b(ni,nc,lmu); integer ia(nx),nin(lmu);
<lam=1,lmu; call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam));>
return;
end;
subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr);
real x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam);
%fortran
      real, dimension (:), allocatable :: w
%mortran
if minval(y).lt.0.0 < jerr=8888; return;>
allocate(w(1:no),stat=jerr); if(jerr.ne.0) return;
w=max(0.0,q); sw=sum(w); if sw.le.0.0 < jerr=9999; go to :done:;>
yb=dot_product(w,y)/sw; fmax=log(huge(y(1))*0.1);
<lam=1,nlam; s=0.0;
   <i=1,no; if(w(i).le.0.0) next;
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:));
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)));
   >
   flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s);
>
:done: deallocate(w);
return;
end;
subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,cl,ne,nx,nlam,flmin,
   ulam,thr,isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr);
real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni);
real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam);
integer ix(*),jx(*),jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: xm,xs,ww,vq
      integer, dimension (:), allocatable :: ju
%mortran
if maxval(vp).le.0.0 < jerr=10000; return;>
if minval(y).lt.0.0 < jerr=8888; return;>
allocate(ww(1:no),stat=jerr);
allocate(ju(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(vq(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xs(1:ni),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
call spchkvars(no,ni,x,ix,ju);
if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0;
if maxval(ju).le.0 < jerr=7777; go to :done:;>
vq=max(0.0,vp); vq=vq*ni/sum(vq);
ww=max(0.0,w); sw=sum(ww); if sw.le.0.0 < jerr=9999; go to :done:;>
ww=ww/sw;
call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs);
if isd.gt.0 < <j=1,ni; cl(:,j)=cl(:,j)*xs(j);>>
call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,thr,
      isd,intr,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr);
if(jerr.gt.0) go to :done:; dev0=2.0*sw*dev0;
<k=1,lmu; nk=nin(k);
   if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk));
   if intr.eq.0 < a0(k)=0.0;>
   else < a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)));>
>
:done:deallocate(ww,ju,vq,xm,xs);
return;
end;
subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,cl,ne,nx,nlam,flmin,ulam,
   shri,isd,intr,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr);
real x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),xs(ni);
real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni);
integer ix(*),jx(*),ju(ni),m(nx),kin(nlam);
%fortran
      real, dimension (:), allocatable :: qy,t,w,wr,v,a,as,xm,ga
      integer, dimension (:), allocatable :: mm,ixx
%mortran
call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx); sml=sml*10.0;
allocate(a(1:ni),stat=jerr);
allocate(as(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(t(1:no),stat=ierr); jerr=jerr+ierr;
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ga(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ixx(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(wr(1:no),stat=ierr); jerr=jerr+ierr;
allocate(v(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(w(1:no),stat=ierr); jerr=jerr+ierr;
allocate(qy(1:no),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
bta=parm; omb=1.0-bta; fmax=log(huge(bta)*0.1);
qy=q*y; yb=sum(qy);
if nonzero(no,g).eq.0 < t=0.0;
   if intr.ne.0 < w=q*yb; az=log(yb); uu=az;
      xm=yb*xb; dv0=yb*(az-1.0);
   >
   else < w=q; xm=0.0; /uu,az/=0.0; dv0=-1.0;>
>
else < w=q*exp(sign(min(abs(g),fmax),g)); ww=sum(w); t=g;
   if intr.ne.0 < eaz=yb/ww;
      w=eaz*w; /az,uu/=log(eaz); dv0=dot_product(qy,g)-yb*(1.0-az);
   >
   else < /uu,az/=0.0; dv0=dot_product(qy,g)-ww;>
   <j=1,ni; if(ju(j).eq.0) next; jb=ix(j); je=ix(j+1)-1;
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je));
   >
>   
tt=yb*uu; ww=1.0; if(intr.ne.0) ww=yb; wr=qy-q*(yb*(1.0-uu)); a=0.0; as=0.0;
dvr=-yb;
<i=1,no; if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i));> dvr=dvr-dv0; dev0=dvr;
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
m=0; mm=0; /nlp,nin/=0; mnl=min(mnlam,nlam); shr=shri*dev0; al=0.0; ixx=0;
<j=1,ni; if(ju(j).eq.0) next;
   jb=ix(j); je=ix(j+1)-1;
   ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))
         -uu*(xm(j)-ww*xb(j))-xb(j)*tt)/xs(j);
>
<ilm=1,nlam; al0=al;
   if flmin.ge.1.0 < al=ulam(ilm);>
   elseif ilm.gt.2 < al=al*alf;>
   elseif ilm.eq.1 < al=big;>
   else < al0=0.0;
      <j=1,ni; if(ju(j).eq.0) next; if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j));>
      al0=al0/max(bta,1.0e-3); al=alf*al0;
   >
   al2=al*omb; al1=al*bta; tlam=bta*(2.0*al-al0);
   <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1;
   >
   :again:continue;
   loop <
      az0=az; if(nin.gt.0) as(m(1:nin))=a(m(1:nin));
      <j=1,ni; if(ixx(j).eq.0) next; jb=ix(j); je=ix(j+1)-1;
         xm(j)=dot_product(w(jx(jb:je)),x(jb:je));
         v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)
                -2.0*xb(j)*xm(j)+ww*xb(j)**2)/xs(j)**2;
      >   
      loop <
         nlp=nlp+1; dlx=0.0;
         <k=1,ni; if(ixx(k).eq.0) next; jb=ix(k); je=ix(k+1)-1; ak=a(k);           
            u=(dot_product(wr(jx(jb:je)),x(jb:je))
               -uu*(xm(k)-ww*xb(k))-xb(k)*tt)/xs(k)+v(k)*ak;
            au=abs(u)-vp(k)*al1;
            if au.le.0.0 < a(k)=0.0;>
            else <
               a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)));
            >
            if(a(k).eq.ak) next;
            if mm(k).eq.0 < nin=nin+1; if(nin.gt.nx) exit;
               mm(k)=nin; m(nin)=k;
            >            
            d=a(k)-ak; dlx=max(dlx,v(k)*d**2); dv=d/xs(k);
            wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je);
            t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je);
            uu=uu-dv*xb(k); tt=tt-dv*xm(k);
         >
         if(nin.gt.nx) exit;
         if intr.ne.0 < d=tt/ww-uu;
            az=az+d; dlx=max(dlx,ww*d**2); uu=uu+d;
         >
         if(dlx.lt.shr) exit; if nlp.gt.maxit < jerr=-ilm; return;>
         loop < nlp=nlp+1; dlx=0.0;
            <l=1,nin; k=m(l);
               jb=ix(k); je=ix(k+1)-1; ak=a(k);
               u=(dot_product(wr(jx(jb:je)),x(jb:je))
                  -uu*(xm(k)-ww*xb(k))-xb(k)*tt)/xs(k)+v(k)*ak;
               au=abs(u)-vp(k)*al1;
               if au.le.0.0 < a(k)=0.0;>
               else <
                  a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)));
               >
               if(a(k).eq.ak) next; d=a(k)-ak; dlx=max(dlx,v(k)*d**2);
               dv=d/xs(k); wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je);
               t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je);
               uu=uu-dv*xb(k); tt=tt-dv*xm(k);
            >
            if intr.ne.0 < d=tt/ww-uu; az=az+d;
               dlx=max(dlx,ww*d**2); uu=uu+d;
            >
            if(dlx.lt.shr) exit; if nlp.gt.maxit < jerr=-ilm; return;>
         >
      >
      if(nin.gt.nx) exit;
      euu=exp(sign(min(abs(uu),fmax),uu));
      w=euu*q*exp(sign(min(abs(t),fmax),t)); ww=sum(w);
      wr=qy-w*(1.0-uu); tt=sum(wr);
      if ww*(az-az0)**2.lt.shr < kx=0;
         <j=1,nin; k=m(j);
            if(v(k)*(a(k)-as(k))**2.lt.shr) next; kx=1; exit;
         >
         if kx.eq.0 <
            <j=1,ni; if(ixx(j).eq.1) next; if(ju(j).eq.0) next;
               jb=ix(j); je=ix(j+1)-1;
               xm(j)=dot_product(w(jx(jb:je)),x(jb:je));
               ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))
                  -uu*(xm(j)-ww*xb(j))-xb(j)*tt)/xs(j);
               if ga(j).gt.al1*vp(j) < ixx(j)=1; kx=1;>
            >
            if(kx.eq.1) go to :again:;
            exit;
         >
      >
   >
   if nin.gt.nx < jerr=-10000-ilm;  exit;>
   if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin)); kin(ilm)=nin;
   a0(ilm)=az; alm(ilm)=al; lmu=ilm;
   dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr;
   if(ilm.lt.mnl) next; if(flmin.ge.1.0) next;
   me=0; <j=1,nin; if(ca(j,ilm).ne.0.0) me=me+1;> if(me.gt.ne) exit;
   if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml) exit;
   if(dev(ilm).gt.devmax) exit;
>
g=t+uu;
:done:deallocate(t,w,wr,v,a,qy,xm,as,mm,ga,ixx);
return;
end;
subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr);
real x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam);
integer ix(*),jx(*);
%fortran
      real, dimension (:), allocatable :: w,f
%mortran
if minval(y).lt.0.0 < jerr=8888; return;>
allocate(w(1:no),stat=jerr);
allocate(f(1:no),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
w=max(0.0,q); sw=sum(w); if sw.le.0.0 < jerr=9999; go to :done:;>
yb=dot_product(w,y)/sw; fmax=log(huge(y(1))*0.1);
<lam=1,nlam; f=a0(lam);
   <j=1,ni; if(a(j,lam).eq.0.0) next; jb=ix(j); je=ix(j+1)-1;
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je);
   >
   f=f+g;
   s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)));
   flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s);
>
:done: deallocate(w,f);
return;
end;
subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,jerr);
real x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(nlam);
integer ix(*),jx(*),nin(nlam),ia(nx);
%fortran
      real, dimension (:), allocatable :: w,f
%mortran
if minval(y).lt.0.0 < jerr=8888; return;>
allocate(w(1:no),stat=jerr);
allocate(f(1:no),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
w=max(0.0,q); sw=sum(w); if sw.le.0.0 < jerr=9999; go to :done:;>
yb=dot_product(w,y)/sw; fmax=log(huge(y(1))*0.1);
<lam=1,nlam; f=a0(lam);
   <k=1,nin(lam); j=ia(k); jb=ix(j); je=ix(j+1)-1;
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je);
   >
   f=f+g;
   s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)));
   flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s);
>
:done: deallocate(w,f);
return;
end;
subroutine multelnet
 (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,jsd,intr,maxit,
   lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
real x(no,ni),y(no,nr),w(no),vp(ni),ca(nx,nr,nlam);
real ulam(nlam),a0(nr,nlam),rsq(nlam),alm(nlam),cl(2,ni);
integer jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: vq;
%mortran
if maxval(vp).le.0.0 < jerr=10000; return;>
allocate(vq(1:ni),stat=jerr); if(jerr.ne.0) return;
vq=max(0.0,vp); vq=vq*ni/sum(vq);
call multelnetn(parm,no,ni,nr,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,isd,
   jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
deallocate(vq);
return;
end;
subroutine multelnetn (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
   isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
real vp(ni),x(no,ni),y(no,nr),w(no),ulam(nlam),cl(2,ni);
real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam);
integer jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: xm,xs,xv,ym,ys
      integer, dimension (:), allocatable :: ju
      real, dimension (:,:,:), allocatable :: clt
      allocate(clt(1:2,1:nr,1:ni),stat=jerr);
%mortran
allocate(xm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xs(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ym(1:nr),stat=ierr); jerr=jerr+ierr;
allocate(ys(1:nr),stat=ierr); jerr=jerr+ierr;
allocate(ju(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xv(1:ni),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
call chkvars(no,ni,x,ju);
if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0;
if maxval(ju).le.0 < jerr=7777; return;>
call multstandard1(no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym,ys,xv,ys0,jerr);
if(jerr.ne.0) return;
<j=1,ni; <k=1,nr; <i=1,2; clt(i,k,j)=cl(i,j);>>>
if isd.gt.0 < <j=1,ni; <k=1,nr; <i=1,2; clt(i,k,j)=clt(i,k,j)*xs(j);>>>>
if jsd.gt.0 < <j=1,ni; <k=1,nr; <i=1,2; clt(i,k,j)=clt(i,k,j)/ys(k);>>>>
call multelnet2(parm,ni,nr,ju,vp,clt,y,no,ne,nx,x,nlam,flmin,ulam,thr,maxit,xv,
   ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr);
if(jerr.gt.0) return;
<k=1,lmu;  nk=nin(k);
   <j=1,nr;
      <l=1,nk; ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l));>
      if intr.eq.0 < a0(j,k)=0.0;>
      else < a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)));>
   >
>
deallocate(xm,xs,ym,ys,ju,xv,clt);
return;
end;
subroutine multstandard1
   (no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym,ys,xv,ys0,jerr);
real x(no,ni),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),ys(nr);
integer ju(ni);
%fortran
      real, dimension (:), allocatable :: v
%mortran
allocate(v(1:no),stat=jerr); if(jerr.ne.0) return;
w=w/sum(w); v=sqrt(w);
if intr.eq.0 <
   <j=1,ni; if(ju(j).eq.0) next; xm(j)=0.0; x(:,j)=v*x(:,j);
      z=dot_product(x(:,j),x(:,j));
      if isd.gt.0 < xbq=dot_product(v,x(:,j))**2; vc=z-xbq;
         xs(j)=sqrt(vc); x(:,j)=x(:,j)/xs(j); xv(j)=1.0+xbq/vc;
      >
      else < xs(j)=1.0; xv(j)=z;>
   >
   ys0=0.0;
   <j=1,nr; ym(j)=0.0; y(:,j)=v*y(:,j);
      z=dot_product(y(:,j),y(:,j));
      if jsd.gt.0 < u=z-dot_product(v,y(:,j))**2; ys0=ys0+z/u;
         ys(j)=sqrt(u); y(:,j)=y(:,j)/ys(j);
      >
      else < ys(j)=1.0; ys0=ys0+z;>
   >
   go to :out:;
>   
<j=1,ni; if(ju(j).eq.0) next;
   xm(j)=dot_product(w,x(:,j)); x(:,j)=v*(x(:,j)-xm(j));
   xv(j)=dot_product(x(:,j),x(:,j)); if(isd.gt.0) xs(j)=sqrt(xv(j));
>
if isd.eq.0 < xs=1.0;>
else < <j=1,ni; if(ju(j).eq.0) next; x(:,j)=x(:,j)/xs(j);>
   xv=1.0;
>
ys0=0.0;
<j=1,nr;
   ym(j)=dot_product(w,y(:,j)); y(:,j)=v*(y(:,j)-ym(j));
   z=dot_product(y(:,j),y(:,j));
   if jsd.gt.0 < ys(j)=sqrt(z); y(:,j)=y(:,j)/ys(j);>
   else < ys0=ys0+z;>
>
if jsd.eq.0 < ys=1.0;> else < ys0=nr;>
:out:deallocate(v);
return;
end;
subroutine multelnet2(beta,ni,nr,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,ulam,thri,
   maxit,xv,ys0,lmu,ao,ia,kin,rsqo,almo,nlp,jerr);
real vp(ni),y(no,nr),x(no,ni),ulam(nlam),ao(nx,nr,nlam);
real rsqo(nlam),almo(nlam),xv(ni),cl(2,nr,ni);
integer ju(ni),ia(nx),kin(nlam);
%fortran
      real, dimension (:), allocatable :: g,gk,del,gj
      integer, dimension (:), allocatable :: mm,ix,isc
      real, dimension (:,:), allocatable :: a
      allocate(a(1:nr,1:ni),stat=jerr)
%mortran
call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx);
allocate(gj(1:nr),stat=ierr); jerr=jerr+ierr;
allocate(gk(1:nr),stat=ierr); jerr=jerr+ierr;
allocate(del(1:nr),stat=ierr); jerr=jerr+ierr;
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(g(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ix(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(isc(1:nr),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
bta=beta; omb=1.0-bta; ix=0; thr=thri*ys0/nr;
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
rsq=ys0; a=0.0; mm=0; /nlp,nin/=0; iz=0; mnl=min(mnlam,nlam); alm=0.0;
<j=1,ni; if(ju(j).eq.0) next; g(j)=0.0;
   <k=1,nr; g(j)=g(j)+dot_product(y(:,k),x(:,j))**2;>
   g(j)=sqrt(g(j));
>
<m=1,nlam; alm0=alm;
   if flmin.ge.1.0 < alm=ulam(m);>
   elseif m.gt.2 < alm=alm*alf;>
   elseif m.eq.1 < alm=big;>
   else < alm0=0.0;
      <j=1,ni; if(ju(j).eq.0) next;
         if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j));
      >
      alm0=alm0/max(bta,1.0e-3); alm=alf*alm0;
   >
   dem=alm*omb; ab=alm*bta; rsq0=rsq; jz=1;
   tlam=bta*(2.0*alm-alm0);
   <k=1,ni; if(ix(k).eq.1) next; if(ju(k).eq.0) next;
      if(g(k).gt.tlam*vp(k)) ix(k)=1;
   >
   loop < if(iz*jz.ne.0) go to :b:;
      :again:nlp=nlp+1; dlx=0.0;
      <k=1,ni; if(ix(k).eq.0) next; gkn=0.0;
         <j=1,nr; gj(j)=dot_product(y(:,j),x(:,k));
            gk(j)=gj(j)+a(j,k)*xv(k); gkn=gkn+gk(j)**2
         >
         gkn=sqrt(gkn); u=1.0-ab*vp(k)/gkn; del=a(:,k);
         if u.le.0.0 < a(:,k)=0.0;>
         else < a(:,k)=gk*(u/(xv(k)+dem*vp(k)));
            call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),
               dem*vp(k),ab*vp(k),a(:,k),isc,jerr);
            if(jerr.ne.0) return;
         >
         del=a(:,k)-del; if(maxval(abs(del)).le.0.0) next;
         <j=1,nr; rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k));
            y(:,j)=y(:,j)-del(j)*x(:,k); dlx=max(dlx,xv(k)*del(j)**2);
         >
         if mm(k).eq.0 < nin=nin+1; if(nin.gt.nx) exit;
            mm(k)=nin; ia(nin)=k;
         >
      >
      if(nin.gt.nx) exit;
      if dlx.lt.thr < ixx=0;
         <k=1,ni; if(ix(k).eq.1) next; if(ju(k).eq.0) next; g(k)=0.0;
            <j=1,nr; g(k)=g(k)+dot_product(y(:,j),x(:,k))**2;>
            g(k)=sqrt(g(k));
            if g(k).gt.ab*vp(k) < ix(k)=1; ixx=1;>
         >
         if(ixx.eq.1) go to :again:;
         exit;
      >
      if nlp.gt.maxit < jerr=-m; return;>
      :b: iz=1;
      loop < nlp=nlp+1; dlx=0.0;
         <l=1,nin; k=ia(l); gkn=0.0;
            <j=1,nr; gj(j)=dot_product(y(:,j),x(:,k));
               gk(j)=gj(j)+a(j,k)*xv(k); gkn=gkn+gk(j)**2
            >
            gkn=sqrt(gkn); u=1.0-ab*vp(k)/gkn; del=a(:,k);
            if u.le.0.0 < a(:,k)=0.0;>
            else < a(:,k)=gk*(u/(xv(k)+dem*vp(k)));
               call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),
                  dem*vp(k),ab*vp(k),a(:,k),isc,jerr);
               if(jerr.ne.0) return;
            >
            del=a(:,k)-del; if(maxval(abs(del)).le.0.0) next;
            <j=1,nr; rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k));
               y(:,j)=y(:,j)-del(j)*x(:,k); dlx=max(dlx,xv(k)*del(j)**2);
            >
         >
         if(dlx.lt.thr) exit; if nlp.gt.maxit < jerr=-m; return;>
      >
      jz=0;
   >
   if nin.gt.nx < jerr=-10000-m;  exit;>
   if nin.gt.0 < <j=1,nr; ao(1:nin,j,m)=a(j,ia(1:nin));>>
   kin(m)=nin;
   rsqo(m)=1.0-rsq/ys0; almo(m)=alm; lmu=m;
   if(m.lt.mnl) next; if(flmin.ge.1.0) next;
   me=0; <j=1,nin; if(ao(j,1,m).ne.0.0) me=me+1;> if(me.gt.ne) exit;
   if(rsq0-rsq.lt.sml*rsq) exit; if(rsqo(m).gt.rsqmax) exit;
>
deallocate(a,mm,g,ix,del,gj,gk);
return;
end;
subroutine chkbnds(nr,gk,gkn,xv,cl,al1,al2,a,isc,jerr);
real gk(nr),cl(2,nr),a(nr); integer isc(nr);
kerr=0; al1p=1.0+al1/xv; al2p=al2/xv; isc=0;
gsq=gkn**2; asq=dot_product(a,a); usq=0.0;
loop < vmx=0.0;
   <k=1,nr; v=max(a(k)-cl(2,k),cl(1,k)-a(k));
      if v.gt.vmx < vmx=v; kn=k;>
   >
   if(vmx.le.0.0) exit; if(isc(kn).ne.0) exit;
   gsq=gsq-gk(kn)**2; g=sqrt(gsq)/xv;
   if(a(kn).lt.cl(1,kn)) u=cl(1,kn); if(a(kn).gt.cl(2,kn)) u=cl(2,kn);
   usq=usq+u**2;
   if usq.eq.0.0 < b=max(0.0,(g-al2p)/al1p);>
   else < b0=sqrt(asq-a(kn)**2);
      b=bnorm(b0,al1p,al2p,g,usq,kerr); if(kerr.ne.0) exit;
   >
   asq=usq+b**2; if asq.le.0.0 < a=0.0; exit;>
   a(kn)=u; isc(kn)=1; f=1.0/(xv*(al1p+al2p/sqrt(asq)));
   <j=1,nr; if(isc(j).eq.0) a(j)=f*gk(j);>
>
if(kerr.ne.0) jerr=kerr;
return;
end;
subroutine chkbnds1(nr,gk,gkn,xv,cl1,cl2,al1,al2,a,isc,jerr);
real gk(nr),a(nr); integer isc(nr);
kerr=0; al1p=1.0+al1/xv; al2p=al2/xv; isc=0;
gsq=gkn**2; asq=dot_product(a,a); usq=0.0;
loop < vmx=0.0;
   <k=1,nr; v=max(a(k)-cl2,cl1-a(k));
      if v.gt.vmx < vmx=v; kn=k;>
   >
   if(vmx.le.0.0) exit; if(isc(kn).ne.0) exit;
   gsq=gsq-gk(kn)**2; g=sqrt(gsq)/xv;
   if(a(kn).lt.cl1) u=cl1; if(a(kn).gt.cl2) u=cl2;
   usq=usq+u**2;
   if usq.eq.0.0 < b=max(0.0,(g-al2p)/al1p);>
   else < b0=sqrt(asq-a(kn)**2);
      b=bnorm(b0,al1p,al2p,g,usq,kerr); if(kerr.ne.0) exit;
   >
   asq=usq+b**2; if asq.le.0.0 < a=0.0; exit;>
   a(kn)=u; isc(kn)=1; f=1.0/(xv*(al1p+al2p/sqrt(asq)));
   <j=1,nr; if(isc(j).eq.0) a(j)=f*gk(j);>
>
if(kerr.ne.0) jerr=kerr;
return;
end;
function bnorm(b0,al1p,al2p,g,usq,jerr);
data thr,mxit /1.0e-10,100/;
b=b0; zsq=b**2+usq; if zsq.le.0.0 < bnorm=0.0; return;>
z=sqrt(zsq); f=b*(al1p+al2p/z)-g; jerr=0;
<it=1,mxit;  b=b-f/(al1p+al2p*usq/(z*zsq));
   zsq=b**2+usq; if zsq.le.0.0 < bnorm=0.0; return;>
   z=sqrt(zsq); f=b*(al1p+al2p/z)-g;
   if(abs(f).le.thr) exit; if b.le.0.0 < b=0.0; exit;>
>
bnorm=b; if(it.ge.mxit) jerr=90000;
return;
entry chg_bnorm(arg,irg); thr=arg; mxit=irg; return;
entry get_bnorm(arg,irg); arg=thr; irg=mxit; return;
end;
subroutine multsolns(ni,nx,nr,lmu,a,ia,nin,b);
real a(nx,nr,lmu),b(ni,nr,lmu); integer ia(nx),nin(lmu);
<lam=1,lmu; call multuncomp(ni,nr,nx,a(1,1,lam),ia,nin(lam),b(1,1,lam));>
return;
end;
subroutine multuncomp(ni,nr,nx,ca,ia,nin,a);
real ca(nx,nr),a(ni,nr); integer ia(nx);
a=0.0;
if nin.gt.0 < <j=1,nr; a(ia(1:nin),j)=ca(1:nin,j);>>
return;
end;
subroutine multmodval(nx,nr,a0,ca,ia,nin,n,x,f);
real a0(nr),ca(nx,nr),x(n,*),f(nr,n); integer ia(nx);
<i=1,n; f(:,i)=a0;> if(nin.le.0) return;
<i=1,n; <j=1,nr; f(j,i)=f(j,i)+dot_product(ca(1:nin,j),x(i,ia(1:nin)));>>
return;
end;
subroutine multspelnet
 (parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
   jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
real x(*),y(no,nr),w(no),vp(ni),ulam(nlam),cl(2,ni);
real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam);
integer ix(*),jx(*),jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: vq;
%mortran
if maxval(vp).le.0.0 < jerr=10000; return;>
allocate(vq(1:ni),stat=jerr); if(jerr.ne.0) return;
vq=max(0.0,vp); vq=vq*ni/sum(vq);
call multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,
   ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
deallocate(vq);
return;
end;
subroutine multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flmin,
   ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr);
real x(*),vp(ni),y(no,nr),w(no),ulam(nlam),cl(2,ni);
real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam);
integer ix(*),jx(*),jd(*),ia(nx),nin(nlam);
%fortran
      real, dimension (:), allocatable :: xm,xs,xv,ym,ys
      integer, dimension (:), allocatable :: ju
      real, dimension (:,:,:), allocatable :: clt
      allocate(clt(1:2,1:nr,1:ni),stat=jerr)
%mortran
allocate(xm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xs(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ym(1:nr),stat=ierr); jerr=jerr+ierr;
allocate(ys(1:nr),stat=ierr); jerr=jerr+ierr;
allocate(ju(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(xv(1:ni),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
call spchkvars(no,ni,x,ix,ju);
if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0;
if maxval(ju).le.0 < jerr=7777; return;>
call multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,
   xm,xs,ym,ys,xv,ys0,jerr);
if(jerr.ne.0) return;
<j=1,ni; <k=1,nr; <i=1,2; clt(i,k,j)=cl(i,j);>>>
if isd.gt.0 < <j=1,ni; <k=1,nr; <i=1,2; clt(i,k,j)=clt(i,k,j)*xs(j);>>>>
if jsd.gt.0 < <j=1,ni; <k=1,nr; <i=1,2; clt(i,k,j)=clt(i,k,j)/ys(k);>>>>
call multspelnet2(parm,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,clt,nlam,flmin,
   ulam,thr,maxit,xm,xs,xv,ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr);
if(jerr.gt.0) return;
<k=1,lmu; nk=nin(k);
   <j=1,nr;
      <l=1,nk; ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l));>
      if intr.eq.0 < a0(j,k)=0.0;>
      else < a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)));>
   >
>
deallocate(xm,xs,ym,ys,ju,xv,clt);
return;
end;
subroutine multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,
   xm,xs,ym,ys,xv,ys0,jerr);
real x(*),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),ys(nr);
integer ix(*),jx(*),ju(ni);
w=w/sum(w);
if intr.eq.0 <
   <j=1,ni; if(ju(j).eq.0) next; xm(j)=0.0; jb=ix(j); je=ix(j+1)-1;
      z=dot_product(w(jx(jb:je)),x(jb:je)**2);
      if isd.gt.0 < xbq=dot_product(w(jx(jb:je)),x(jb:je))**2; vc=z-xbq;
         xs(j)=sqrt(vc); xv(j)=1.0+xbq/vc;
      >
      else < xs(j)=1.0; xv(j)=z;>
   >
   ys0=0.0;   
   <j=1,nr; ym(j)=0.0; z=dot_product(w,y(:,j)**2);
      if jsd.gt.0 < u=z-dot_product(w,y(:,j))**2; ys0=ys0+z/u;
         ys(j)=sqrt(u); y(:,j)=y(:,j)/ys(j);
      >
      else < ys(j)=1.0; ys0=ys0+z;>
   >
   return;
>
<j=1,ni; if(ju(j).eq.0) next;
   jb=ix(j); je=ix(j+1)-1; xm(j)=dot_product(w(jx(jb:je)),x(jb:je));
   xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2;
   if(isd.gt.0) xs(j)=sqrt(xv(j));
>
if isd.eq.0 < xs=1.0;> else < xv=1.0;>
ys0=0.0;
<j=1,nr;
   ym(j)=dot_product(w,y(:,j)); y(:,j)=y(:,j)-ym(j);
   z=dot_product(w,y(:,j)**2);
   if jsd.gt.0 < ys(j)=sqrt(z); y(:,j)=y(:,j)/ys(j);>
   else < ys0=ys0+z;>
>
if jsd.eq.0 < ys=1.0;> else < ys0=nr;>
return;
end;
subroutine multspelnet2(beta,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,
   ulam,thri,maxit,xm,xs,xv,ys0,lmu,ao,ia,kin,rsqo,almo,nlp,jerr);
real y(no,nr),w(no),x(*),vp(ni),ulam(nlam),cl(2,nr,ni);
real ao(nx,nr,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni);
integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam);
%fortran
      real, dimension (:), allocatable :: g,gj,gk,del,o
      integer, dimension (:), allocatable :: mm,iy,isc
      real, dimension (:,:), allocatable :: a
      allocate(a(1:nr,1:ni),stat=jerr)
%mortran
call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx);
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(g(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(gj(1:nr),stat=ierr); jerr=jerr+ierr;
allocate(gk(1:nr),stat=ierr); jerr=jerr+ierr;
allocate(del(1:nr),stat=ierr); jerr=jerr+ierr;
allocate(o(1:nr),stat=ierr); jerr=jerr+ierr;
allocate(iy(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(isc(1:nr),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
bta=beta; omb=1.0-bta; alm=0.0; iy=0; thr=thri*ys0/nr;
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
rsq=ys0; a=0.0; mm=0; o=0.0; /nlp,nin/=0; iz=0; mnl=min(mnlam,nlam);
<j=1,ni; if(ju(j).eq.0) next; jb=ix(j); je=ix(j+1)-1; g(j)=0.0;
   <k=1,nr;
      g(j)=g(j)+(dot_product(y(jx(jb:je),k),w(jx(jb:je))*x(jb:je))/xs(j))**2;
   >
   g(j)=sqrt(g(j));
>
<m=1,nlam; alm0=alm;
   if flmin.ge.1.0 < alm=ulam(m);>
   elseif m.gt.2 < alm=alm*alf;>
   elseif m.eq.1 < alm=big;>
   else < alm0=0.0;
      <j=1,ni; if(ju(j).eq.0) next;
         if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j));
      >
      alm0=alm0/max(bta,1.0e-3); alm=alf*alm0;
   >
   dem=alm*omb; ab=alm*bta; rsq0=rsq; jz=1;
   tlam=bta*(2.0*alm-alm0);
   <k=1,ni; if(iy(k).eq.1) next; if(ju(k).eq.0) next;
      if(g(k).gt.tlam*vp(k)) iy(k)=1;
   >
   loop < if(iz*jz.ne.0) go to :b:;
      :again:nlp=nlp+1; dlx=0.0;
      <k=1,ni; if(iy(k).eq.0) next; jb=ix(k); je=ix(k+1)-1; gkn=0.0;
         <j=1,nr;
            gj(j)=dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(k);
            gk(j)=gj(j)+a(j,k)*xv(k); gkn=gkn+gk(j)**2;
         >
         gkn=sqrt(gkn); u=1.0-ab*vp(k)/gkn; del=a(:,k);
         if u.le.0.0 < a(:,k)=0.0;>
         else < a(:,k)=gk*(u/(xv(k)+dem*vp(k)));
            call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),
               dem*vp(k),ab*vp(k),a(:,k),isc,jerr);
            if(jerr.ne.0) return;
         >
         del=a(:,k)-del; if(maxval(abs(del)).le.0.0) next;
         if mm(k).eq.0 < nin=nin+1; if(nin.gt.nx) exit;
            mm(k)=nin; ia(nin)=k;
         >
         <j=1,nr; rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k));
            y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k);
            o(j)=o(j)+del(j)*xm(k)/xs(k); dlx=max(xv(k)*del(j)**2,dlx);
         >
      >
      if(nin.gt.nx) exit;
      if dlx.lt.thr < ixx=0;
         <j=1,ni; if(iy(j).eq.1) next; if(ju(j).eq.0) next;
            jb=ix(j); je=ix(j+1)-1; g(j)=0.0;
            <k=1,nr; g(j)=g(j)+
             (dot_product(y(jx(jb:je),k)+o(k),w(jx(jb:je))*x(jb:je))/xs(j))**2;
            >
            g(j)=sqrt(g(j));
            if g(j).gt.ab*vp(j) < iy(j)=1; ixx=1;>
         >
         if(ixx.eq.1) go to :again:;
         exit;
      >
      if nlp.gt.maxit < jerr=-m; return;>
      :b: iz=1;
      loop < nlp=nlp+1; dlx=0.0;
         <l=1,nin; k=ia(l); jb=ix(k); je=ix(k+1)-1; gkn=0.0;
            <j=1,nr; gj(j)=
                 dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(k);
               gk(j)=gj(j)+a(j,k)*xv(k); gkn=gkn+gk(j)**2;
            >
            gkn=sqrt(gkn); u=1.0-ab*vp(k)/gkn; del=a(:,k);
            if u.le.0.0 < a(:,k)=0.0;>
            else < a(:,k)=gk*(u/(xv(k)+dem*vp(k)));
               call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),
                  dem*vp(k),ab*vp(k),a(:,k),isc,jerr);
               if(jerr.ne.0) return;
            >
            del=a(:,k)-del; if(maxval(abs(del)).le.0.0) next;
            <j=1,nr; rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k));
               y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k);
               o(j)=o(j)+del(j)*xm(k)/xs(k); dlx=max(xv(k)*del(j)**2,dlx);
            >
         >
         if(dlx.lt.thr) exit; if nlp.gt.maxit < jerr=-m; return;>
      >
      jz=0;
   >
   if nin.gt.nx < jerr=-10000-m;  exit;>
   if nin.gt.0 < <j=1,nr; ao(1:nin,j,m)=a(j,ia(1:nin));>>
   kin(m)=nin;
   rsqo(m)=1.0-rsq/ys0; almo(m)=alm; lmu=m;
   if(m.lt.mnl) next; if(flmin.ge.1.0) next;
   me=0; <j=1,nin; if(ao(j,1,m).ne.0.0) me=me+1;> if(me.gt.ne) exit;
   if(rsq0-rsq.lt.sml*rsq) exit; if(rsqo(m).gt.rsqmax) exit;
>
deallocate(a,mm,g,iy,gj,gk,del,o);
return;
end;
subroutine multlognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin,ulam,
   shri,intr,maxit,xv,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr);
real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),cl(2,ni);
real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),xv(ni);
integer ju(ni),m(nx),kin(nlam);
%fortran
      real, dimension (:,:), allocatable :: q,r,b,bs
      real, dimension (:), allocatable :: sxp,sxpl,ga,gk,del
      integer, dimension (:), allocatable :: mm,is,ixx,isc
      allocate(b(0:ni,1:nc),stat=jerr)
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(r(1:no,1:nc),stat=ierr); jerr=jerr+ierr;
%mortran
call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx); exmn=-exmx;
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(is(1:max(nc,ni)),stat=ierr); jerr=jerr+ierr;
allocate(sxp(1:no),stat=ierr); jerr=jerr+ierr;
allocate(sxpl(1:no),stat=ierr); jerr=jerr+ierr;
allocate(ga(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ixx(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(gk(1:nc),stat=ierr); jerr=jerr+ierr;
allocate(del(1:nc),stat=ierr); jerr=jerr+ierr;
allocate(isc(1:nc),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
pmax=1.0-pmin; emin=pmin/pmax; emax=1.0/emin;
bta=parm; omb=1.0-bta; dev1=0.0; dev0=0.0;
<ic=1,nc; q0=dot_product(w,y(:,ic));
   if q0.le.pmin < jerr =8000+ic; return;>
   if q0.ge.pmax < jerr =9000+ic; return;>
   if intr.eq.0 < q0=1.0/nc; b(0,ic)=0.0;>
   else < b(0,ic)=log(q0); dev1=dev1-q0*b(0,ic);>
   b(1:ni,ic)=0.0;
>
if(intr.eq.0) dev1=log(float(nc)); ixx=0; al=0.0;
if nonzero(no*nc,g).eq.0 <
   b(0,:)=b(0,:)-sum(b(0,:))/nc; sxp=0.0;
   <ic=1,nc; q(:,ic)=exp(b(0,ic)); sxp=sxp+q(:,ic);>
>
else < <i=1,no; g(i,:)=g(i,:)-sum(g(i,:))/nc;> sxp=0.0;
   if intr.eq.0 < b(0,:)=0.0;>
   else < call kazero(nc,no,y,g,w,b(0,:),jerr); if(jerr.ne.0) return;>
   dev1=0.0;
   <ic=1,nc; q(:,ic)=b(0,ic)+g(:,ic);
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic));
      q(:,ic)=exp(q(:,ic)); sxp=sxp+q(:,ic);
   >
   sxpl=w*log(sxp); <ic=1,nc; dev1=dev1+dot_product(y(:,ic),sxpl);>
>
<ic=1,nc; <i=1,no; if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic));>>
dev0=dev0+dev1;
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
m=0; mm=0; nin=0; nlp=0; mnl=min(mnlam,nlam); bs=0.0; shr=shri*dev0;
ga=0.0;
<ic=1,nc; r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp);
   <j=1,ni; if(ju(j).ne.0) ga(j)=ga(j)+dot_product(r(:,ic),x(:,j))**2;>
>
ga=sqrt(ga);
<ilm=1,nlam; al0=al;
   if flmin.ge.1.0 < al=ulam(ilm);>
   elseif ilm.gt.2 < al=al*alf;>
   elseif ilm.eq.1 < al=big;>
   else < al0=0.0;
      <j=1,ni; if(ju(j).eq.0) next; if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j));>
      al0=al0/max(bta,1.0e-3); al=alf*al0;
   >
   al2=al*omb; al1=al*bta; tlam=bta*(2.0*al-al0);
   <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1;
   >
   :again:continue;
   loop < /ix,jx,kx/=0; t=0.0;
      <ic=1,nc; t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp));>
      if t.lt.eps < kx=1; exit;> t=2.0*t; alt=al1/t; al2t=al2/t;
      <ic=1,nc;
         bs(0,ic)=b(0,ic); if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic);
         r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t;
         d=0.0; if(intr.ne.0) d=sum(r(:,ic));
         if d.ne.0.0 <
            b(0,ic)=b(0,ic)+d; r(:,ic)=r(:,ic)-d*w; dlx=max(dlx,d**2);
         >
      >
      loop < nlp=nlp+nc; dlx=0.0;
         <k=1,ni; if(ixx(k).eq.0) next; gkn=0.0;
            <ic=1,nc; gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k);
               gkn=gkn+gk(ic)**2;
            >
            gkn=sqrt(gkn); u=1.0-alt*vp(k)/gkn; del=b(k,:);
            if u.le.0.0 < b(k,:)=0.0;>
            else < b(k,:)=gk*(u/(xv(k)+vp(k)*al2t));
               call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),
                  cl(2,k),vp(k)*al2t,alt*vp(k),b(k,:),isc,jerr);
               if(jerr.ne.0) return;
             >
            del=b(k,:)-del; if(maxval(abs(del)).le.0.0) next;
            <ic=1,nc; dlx=max(dlx,xv(k)*del(ic)**2);
               r(:,ic)=r(:,ic)-del(ic)*w*x(:,k);
            >
            if mm(k).eq.0 < nin=nin+1;
               if nin.gt.nx < jx=1; exit;>
               mm(k)=nin; m(nin)=k;
            >
         >
         if(jx.gt.0) exit; if(dlx.lt.shr) exit;
         if nlp.gt.maxit < jerr=-ilm; return;>
         loop < nlp=nlp+nc; dlx=0.0;
            <l=1,nin; k=m(l); gkn=0.0;
               <ic=1,nc; gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k);
                  gkn=gkn+gk(ic)**2;
               >
               gkn=sqrt(gkn); u=1.0-alt*vp(k)/gkn; del=b(k,:);
               if u.le.0.0 < b(k,:)=0.0;>
               else < b(k,:)=gk*(u/(xv(k)+vp(k)*al2t));
                  call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),
                     cl(2,k),vp(k)*al2t,alt*vp(k),b(k,:),isc,jerr);
                  if(jerr.ne.0) return;
               >
               del=b(k,:)-del; if(maxval(abs(del)).le.0.0) next;
               <ic=1,nc; dlx=max(dlx,xv(k)*del(ic)**2);
                  r(:,ic)=r(:,ic)-del(ic)*w*x(:,k);
               >
            >
            if(dlx.lt.shr) exit; if nlp.gt.maxit < jerr=-ilm; return;>

         >
      >
      if(jx.gt.0) exit;
      <ic=1,nc;
         if((b(0,ic)-bs(0,ic))**2.gt.shr) ix=1;
         if ix.eq.0 <
            <j=1,nin; k=m(j);
               if xv(k)*(b(k,ic)-bs(k,ic))**2.gt.shr < ix=1; exit;>

            >
         >
         <i=1,no; fi=b(0,ic)+g(i,ic);
            if(nin.gt.0)
              fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)));
            fi=min(max(exmn,fi),exmx); sxp(i)=sxp(i)-q(i,ic);
            q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i));
            sxp(i)=sxp(i)+q(i,ic);
         >
      >
      s=-sum(b(0,:))/nc; b(0,:)=b(0,:)+s;
      if(jx.gt.0) exit;
      if ix.eq.0 <
         <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next; ga(k)=0.0;>
         <ic=1,nc; r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp);
            <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
               ga(k)=ga(k)+dot_product(r(:,ic),x(:,k))**2;
            >
         >
         ga=sqrt(ga);
         <k=1,ni; if(ixx(k).eq.1) next; if(ju(k).eq.0) next;
            if ga(k).gt.al1*vp(k) < ixx(k)=1; ix=1;>
         >
         if(ix.eq.1) go to :again:;
         exit;
      >
   >
   if kx.gt.0 < jerr=-20000-ilm;  exit;>
   if jx.gt.0 < jerr=-10000-ilm;  exit;> devi=0.0;
   <ic=1,nc;
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic); a0(ic,ilm)=b(0,ic);
      <i=1,no; if (y(i,ic).le.0.0) next;
         devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i));
      >
   >
   kin(ilm)=nin; alm(ilm)=al; lmu=ilm;
   dev(ilm)=(dev1-devi)/dev0;
   if(ilm.lt.mnl) next; if(flmin.ge.1.0) next;
   me=0; <j=1,nin; if(a(j,1,ilm).ne.0.0) me=me+1;> if(me.gt.ne) exit;
   if(dev(ilm).gt.devmax) exit; if(dev(ilm)-dev(ilm-1).lt.sml) exit;
>
g=log(q); <i=1,no; g(i,:)=g(i,:)-sum(g(i,:))/nc;>
deallocate(sxp,b,bs,r,q,mm,is,ga,ixx,gk,del,sxpl);
return;
end;
subroutine multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,nlam,
   flmin,ulam,shri,intr,maxit,xv,xb,xs,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr);
real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni),xv(ni);
real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni);
integer ix(*),jx(*),ju(ni),m(nx),kin(nlam);
%fortran
      real, dimension (:,:), allocatable :: q,r,b,bs
      real, dimension (:), allocatable :: sxp,sxpl,ga,gk,del,sc,svr
      integer, dimension (:), allocatable :: mm,is,iy,isc
      allocate(b(0:ni,1:nc),stat=jerr)
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr
      allocate(r(1:no,1:nc),stat=ierr); jerr=jerr+ierr
%mortran
call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx); exmn=-exmx;
allocate(mm(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(ga(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(gk(1:nc),stat=ierr); jerr=jerr+ierr;
allocate(del(1:nc),stat=ierr); jerr=jerr+ierr;
allocate(iy(1:ni),stat=ierr); jerr=jerr+ierr;
allocate(is(1:max(nc,ni)),stat=ierr); jerr=jerr+ierr;
allocate(sxp(1:no),stat=ierr); jerr=jerr+ierr;
allocate(sxpl(1:no),stat=ierr); jerr=jerr+ierr;
allocate(svr(1:nc),stat=ierr); jerr=jerr+ierr;
allocate(sc(1:no),stat=ierr); jerr=jerr+ierr;
allocate(isc(1:nc),stat=ierr); jerr=jerr+ierr;
if(jerr.ne.0) return;
pmax=1.0-pmin; emin=pmin/pmax; emax=1.0/emin;
bta=parm; omb=1.0-bta; dev1=0.0; dev0=0.0;
<ic=1,nc; q0=dot_product(w,y(:,ic));
   if q0.le.pmin < jerr =8000+ic; return;>
   if q0.ge.pmax < jerr =9000+ic; return;>
   b(1:ni,ic)=0.0;
   if intr.eq.0 < q0=1.0/nc; b(0,ic)=0.0;>
   else < b(0,ic)=log(q0); dev1=dev1-q0*b(0,ic);>
>
if(intr.eq.0) dev1=log(float(nc)); iy=0; al=0.0;
if nonzero(no*nc,g).eq.0 <
   b(0,:)=b(0,:)-sum(b(0,:))/nc; sxp=0.0;
   <ic=1,nc; q(:,ic)=exp(b(0,ic)); sxp=sxp+q(:,ic);>
>
else < <i=1,no; g(i,:)=g(i,:)-sum(g(i,:))/nc;> sxp=0.0;
   if intr.eq.0 < b(0,:)=0.0;>
   else < call kazero(nc,no,y,g,w,b(0,:),jerr); if(jerr.ne.0) return;>
   dev1=0.0;
   <ic=1,nc; q(:,ic)=b(0,ic)+g(:,ic);
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic));
      q(:,ic)=exp(q(:,ic)); sxp=sxp+q(:,ic);
   >
   sxpl=w*log(sxp); <ic=1,nc; dev1=dev1+dot_product(y(:,ic),sxpl);>
>
<ic=1,nc; <i=1,no; if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic));>>
dev0=dev0+dev1;
if flmin.lt.1.0 < eqs=max(eps,flmin); alf=eqs**(1.0/(nlam-1));>
m=0; mm=0; nin=0; nlp=0; mnl=min(mnlam,nlam); bs=0.0;
shr=shri*dev0; ga=0.0;
<ic=1,nc; r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp); svr(ic)=sum(r(:,ic));
   <j=1,ni; if(ju(j).eq.0) next;
      jb=ix(j); je=ix(j+1)-1;
      gj=dot_product(r(jx(jb:je),ic),x(jb:je));
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2;
   >
>
ga=sqrt(ga);
<ilm=1,nlam; al0=al;
   if flmin.ge.1.0 < al=ulam(ilm);>
   elseif ilm.gt.2 < al=al*alf;>
   elseif ilm.eq.1 < al=big;>
   else < al0=0.0;
      <j=1,ni; if(ju(j).eq.0) next; if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j));>
      al0=al0/max(bta,1.0e-3); al=alf*al0;
   >
   al2=al*omb; al1=al*bta; tlam=bta*(2.0*al-al0);
   <k=1,ni; if(iy(k).eq.1) next; if(ju(k).eq.0) next;
      if(ga(k).gt.tlam*vp(k)) iy(k)=1;
   >
   :again:continue;
   loop < /ixx,jxx,kxx/=0; t=0.0;
      <ic=1,nc; t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp));>
      if t.lt.eps < kxx=1; exit;> t=2.0*t; alt=al1/t; al2t=al2/t;
      <ic=1,nc; bs(0,ic)=b(0,ic); if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic);
         r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t; svr(ic)=sum(r(:,ic));
         if intr.ne.0 < b(0,ic)=b(0,ic)+svr(ic); r(:,ic)=r(:,ic)-svr(ic)*w;
            dlx=max(dlx,svr(ic)**2);
         >
      >
      loop < nlp=nlp+nc; dlx=0.0;
         <k=1,ni; if(iy(k).eq.0) next;
            jb=ix(k); je=ix(k+1)-1; del=b(k,:); gkn=0.0;
            <ic=1,nc;
              u=(dot_product(r(jx(jb:je),ic),x(jb:je))-svr(ic)*xb(k))/xs(k);
               gk(ic)=u+del(ic)*xv(k); gkn=gkn+gk(ic)**2;
            >
            gkn=sqrt(gkn); u=1.0-alt*vp(k)/gkn;
            if u.le.0.0 < b(k,:)=0.0;>
            else <
               b(k,:)=gk*(u/(xv(k)+vp(k)*al2t));
               call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),
                  vp(k)*al2t,alt*vp(k),b(k,:),isc,jerr);
               if(jerr.ne.0) return;
            >
            del=b(k,:)-del; if(maxval(abs(del)).le.0.0) next;
            <ic=1,nc; dlx=max(dlx,xv(k)*del(ic)**2);
               r(jx(jb:je),ic)=r(jx(jb:je),ic)
                  -del(ic)*w(jx(jb:je))*(x(jb:je)-xb(k))/xs(k);
            >
            if mm(k).eq.0 < nin=nin+1;
               if nin.gt.nx < jxx=1; exit;>
               mm(k)=nin; m(nin)=k;
            >
         >
         if(jxx.gt.0) exit;
         if(dlx.lt.shr) exit; if nlp.gt.maxit < jerr=-ilm; return;>
         loop < nlp=nlp+nc; dlx=0.0;
            <l=1,nin; k=m(l); jb=ix(k); je=ix(k+1)-1; del=b(k,:); gkn=0.0;
               <ic=1,nc;
                  u=(dot_product(r(jx(jb:je),ic),x(jb:je))
                     -svr(ic)*xb(k))/xs(k);
                  gk(ic)=u+del(ic)*xv(k); gkn=gkn+gk(ic)**2;
               >
               gkn=sqrt(gkn); u=1.0-alt*vp(k)/gkn;
               if u.le.0.0 < b(k,:)=0.0;>
               else <
                  b(k,:)=gk*(u/(xv(k)+vp(k)*al2t));
                  call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),
                     vp(k)*al2t,alt*vp(k),b(k,:),isc,jerr);
                  if(jerr.ne.0) return;
               >
               del=b(k,:)-del; if(maxval(abs(del)).le.0.0) next;
               <ic=1,nc; dlx=max(dlx,xv(k)*del(ic)**2);
                  r(jx(jb:je),ic)=r(jx(jb:je),ic)
                     -del(ic)*w(jx(jb:je))*(x(jb:je)-xb(k))/xs(k);
               >
            >
            if(dlx.lt.shr) exit; if nlp.gt.maxit < jerr=-ilm; return;>

         >
      >
      if(jxx.gt.0) exit;
      <ic=1,nc;
         if((b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1;
         if ixx.eq.0 <
            <j=1,nin; k=m(j);
               if xv(k)*(b(k,ic)-bs(k,ic))**2.gt.shr < ixx=1; exit;>

            >
         >
         sc=b(0,ic)+g(:,ic); b0=0.0;
         <j=1,nin; l=m(j); jb=ix(l); je=ix(l+1)-1;
            sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l);
            b0=b0-b(l,ic)*xb(l)/xs(l);
         >
         sc=min(max(exmn,sc+b0),exmx);
         sxp=sxp-q(:,ic);
         q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp);
         sxp=sxp+q(:,ic);
      >
      s=sum(b(0,:))/nc; b(0,:)=b(0,:)-s;
      if(jxx.gt.0) exit;
      if ixx.eq.0 <
         <j=1,ni; if(iy(j).eq.1) next; if(ju(j).eq.0) next; ga(j)=0.0;>
         <ic=1,nc; r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp);
            <j=1,ni; if(iy(j).eq.1) next; if(ju(j).eq.0) next;
               jb=ix(j); je=ix(j+1)-1;
               gj=dot_product(r(jx(jb:je),ic),x(jb:je));
               ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2;
            >
         >
         ga=sqrt(ga);
         <k=1,ni; if(iy(k).eq.1) next; if(ju(k).eq.0) next;
            if ga(k).gt.al1*vp(k) < iy(k)=1; ixx=1;>
         >
         if(ixx.eq.1) go to :again:;
         exit;
      >
   >
   if kxx.gt.0 < jerr=-20000-ilm;  exit;>
   if jxx.gt.0 < jerr=-10000-ilm; exit;> devi=0.0;
   <ic=1,nc;
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic); a0(ic,ilm)=b(0,ic);
      <i=1,no; if (y(i,ic).le.0.0) next;
         devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i));
      >
   >
   kin(ilm)=nin; alm(ilm)=al; lmu=ilm;
   dev(ilm)=(dev1-devi)/dev0;
   if(ilm.lt.mnl) next; if(flmin.ge.1.0) next;
   me=0; <j=1,nin; if(a(j,1,ilm).ne.0.0) me=me+1;> if(me.gt.ne) exit;
   if(dev(ilm).gt.devmax) exit; if(dev(ilm)-dev(ilm-1).lt.sml) exit;
>
g=log(q); <i=1,no; g(i,:)=g(i,:)-sum(g(i,:))/nc;>
deallocate(sxp,b,bs,r,q,mm,is,sc,ga,iy,gk,del,sxpl);
return;
end;
%fortran
      subroutine psort7 (v,a,ii,jj)
c
c     puts into a the permutation vector which sorts v into
c     increasing order. the array v is not modified.
c     only elements from ii to jj are considered.
c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements
c
c     this is a modification of cacm algorithm #347 by r. c. singleton,
c     which is a modified hoare quicksort.
c
      dimension a(jj),v(jj),iu(20),il(20)
      integer t,tt
      integer a
      real v
      m=1
      i=ii
      j=jj
 10   if (i.ge.j) go to 80
 20   k=i
      ij=(j+i)/2
      t=a(ij)
      vt=v(t)
      if (v(a(i)).le.vt) go to 30
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      vt=v(t)
 30   l=j
      if (v(a(j)).ge.vt) go to 50
      a(ij)=a(j)
      a(j)=t
      t=a(ij)
      vt=v(t)
      if (v(a(i)).le.vt) go to 50
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      vt=v(t)
      go to 50
 40   a(l)=a(k)
      a(k)=tt
 50   l=l-1
      if (v(a(l)).gt.vt) go to 50
      tt=a(l)
      vtt=v(tt)
 60   k=k+1
      if (v(a(k)).lt.vt) go to 60
      if (k.le.l) go to 40
      if (l-i.le.j-k) go to 70
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 90
 70   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 90
 80   m=m-1
      if (m.eq.0) return
      i=il(m)
      j=iu(m)
 90   if (j-i.gt.10) go to 20
      if (i.eq.ii) go to 10
      i=i-1
 100  i=i+1
      if (i.eq.j) go to 80
      t=a(i+1)
      vt=v(t)
      if (v(a(i)).le.vt) go to 100
      k=i
 110  a(k+1)=a(k)
      k=k-1
      if (vt.lt.v(a(k))) go to 110
      a(k+1)=t
      go to 100
      end
%%
