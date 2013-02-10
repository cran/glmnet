c
c                          newGLMnet (2/07/13)
c
c
c                 Elastic net with squared-error loss
c
c dense predictor matrix:
c
c call elnet(ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,maxit,
c            lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call spelnet(ka,parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
c             isd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
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
c                jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call multspelnet(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
c             isd,jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
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
c              maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call splognet (parm,no,ni,nc,x,ix,jx,y,o,jd,vp,cl,ne,nx,nlam,flmin,
c             ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse row format
c
c
c other inputs:
c
c   parm,no,ni,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,maxit
c    = same as elnet above.
c
c   nc = number of classes (distinct outcome values)
c        nc=1 => binomial two-class logistic regression
c            (all output references class 1)
c   y(no,max(2,nc)) = number of each class at each design point(overwritten)
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
c               isd,maxit,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c sparse predictor matrix:
c
c call spfishnet (parm,no,ni,x,ix,jx,y,o,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,
c               isd,maxit,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c    x, ix, jx = predictor data matrix in compressed sparse row format
c
c other inputs:
c
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   parm,no,ni,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,maxit
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
c   pmin = minimum null probability for any class. default = 1.0e-5
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
      subroutine get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)          769
      data sml0,eps0,big0,mnlam0,rsqmax0,pmin0,exmx0  /1.0e-5,1.0e-6,9.9    771 
     *e35,5,0.999,1.0e-5,250.0/
      sml=sml0                                                              771
      eps=eps0                                                              771
      big=big0                                                              771
      mnlam=mnlam0                                                          771
      rsqmax=rsqmax0                                                        772
      pmin=pmin0                                                            772
      exmx=exmx0                                                            773
      return                                                                774
      entry chg_fract_dev(arg)                                              774
      sml0=arg                                                              774
      return                                                                775
      entry chg_dev_max(arg)                                                775
      rsqmax0=arg                                                           775
      return                                                                776
      entry chg_min_flmin(arg)                                              776
      eps0=arg                                                              776
      return                                                                777
      entry chg_big(arg)                                                    777
      big0=arg                                                              777
      return                                                                778
      entry chg_min_lambdas(irg)                                            778
      mnlam0=irg                                                            778
      return                                                                779
      entry chg_min_null_prob(arg)                                          779
      pmin0=arg                                                             779
      return                                                                780
      entry chg_max_exp(arg)                                                780
      exmx0=arg                                                             780
      return                                                                781
      end                                                                   782
      subroutine elnet  (ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,u    785 
     *lam,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam),cl(2,ni)                 786
      real ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)                          787
      integer jd(*),ia(nx),nin(nlam)                                        788
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10021                                     791
      jerr=10000                                                            791
      return                                                                791
10021 continue                                                              792
      allocate(vq(1:ni),stat=jerr)                                          792
      if(jerr.ne.0) return                                                  793
      vq=max(0.0,vp)                                                        793
      vq=vq*ni/sum(vq)                                                      794
      if(ka .ne. 1)goto 10041                                               795
      call elnetu  (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,    798 
     *isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                            799
10041 continue                                                              800
      call elnetn (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,i    803 
     *sd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                              804
10031 continue                                                              804
      deallocate(vq)                                                        805
      return                                                                806
      end                                                                   807
      subroutine elnetu  (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ula    810 
     *m,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)                  811
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         812
      integer jd(*),ia(nx),nin(nlam)                                        813
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           818
      allocate(xm(1:ni),stat=ierr)                                          818
      jerr=jerr+ierr                                                        819
      allocate(xs(1:ni),stat=ierr)                                          819
      jerr=jerr+ierr                                                        820
      allocate(ju(1:ni),stat=ierr)                                          820
      jerr=jerr+ierr                                                        821
      allocate(xv(1:ni),stat=ierr)                                          821
      jerr=jerr+ierr                                                        822
      allocate(vlam(1:nlam),stat=ierr)                                      822
      jerr=jerr+ierr                                                        823
      if(jerr.ne.0) return                                                  824
      call chkvars(no,ni,x,ju)                                              825
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  826
      if(maxval(ju) .gt. 0)goto 10071                                       826
      jerr=7777                                                             826
      return                                                                826
10071 continue                                                              827
      call standard(no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)               828
      if(jerr.ne.0) return                                                  829
      cl=cl/ys                                                              829
      if(isd .le. 0)goto 10091                                              829
10100 do 10101 j=1,ni                                                       829
      cl(:,j)=cl(:,j)*xs(j)                                                 829
10101 continue                                                              829
10102 continue                                                              829
10091 continue                                                              830
      if(flmin.ge.1.0) vlam=ulam/ys                                         831
      call elnet1(parm,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxi    833 
     *t,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  834
10110 do 10111 k=1,lmu                                                      834
      alm(k)=ys*alm(k)                                                      834
      nk=nin(k)                                                             835
10120 do 10121 l=1,nk                                                       835
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          835
10121 continue                                                              836
10122 continue                                                              836
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         837
10111 continue                                                              838
10112 continue                                                              838
      deallocate(xm,xs,g,ju,xv,vlam)                                        839
      return                                                                840
      end                                                                   841
      subroutine standard (no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)        842
      real x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                  842
      integer ju(ni)                                                        843
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           846
      if(jerr.ne.0) return                                                  847
      w=w/sum(w)                                                            847
      v=sqrt(w)                                                             848
10130 do 10131 j=1,ni                                                       848
      if(ju(j).eq.0)goto 10131                                              849
      xm(j)=dot_product(w,x(:,j))                                           849
      x(:,j)=v*(x(:,j)-xm(j))                                               850
      xv(j)=dot_product(x(:,j),x(:,j))                                      850
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        851
10131 continue                                                              852
10132 continue                                                              852
      if(isd .ne. 0)goto 10151                                              852
      xs=1.0                                                                852
      goto 10161                                                            853
10151 continue                                                              854
10170 do 10171 j=1,ni                                                       854
      if(ju(j).eq.0)goto 10171                                              854
      x(:,j)=x(:,j)/xs(j)                                                   854
10171 continue                                                              855
10172 continue                                                              855
      xv=1.0                                                                856
10161 continue                                                              857
10141 continue                                                              857
      ym=dot_product(w,y)                                                   857
      y=v*(y-ym)                                                            857
      ys=sqrt(dot_product(y,y))                                             857
      y=y/ys                                                                857
      g=0.0                                                                 858
10180 do 10181 j=1,ni                                                       858
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                             858
10181 continue                                                              859
10182 continue                                                              859
      deallocate(v)                                                         860
      return                                                                861
      end                                                                   862
      subroutine elnet1 (beta,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,ulam,t    864 
     *hr,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      real vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    865 
     *nlam),xv(ni)
      real cl(2,ni)                                                         866
      integer ju(ni),ia(nx),kin(nlam)                                       867
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)                874
      allocate(a(1:ni),stat=ierr)                                           874
      jerr=jerr+ierr                                                        875
      allocate(mm(1:ni),stat=ierr)                                          875
      jerr=jerr+ierr                                                        876
      allocate(da(1:ni),stat=ierr)                                          876
      jerr=jerr+ierr                                                        877
      if(jerr.ne.0) return                                                  878
      bta=beta                                                              878
      omb=1.0-bta                                                           879
      if(flmin .ge. 1.0)goto 10201                                          879
      eqs=max(eps,flmin)                                                    879
      alf=eqs**(1.0/(nlam-1))                                               879
10201 continue                                                              880
      rsq=0.0                                                               880
      a=0.0                                                                 880
      mm=0                                                                  880
      nlp=0                                                                 880
      nin=nlp                                                               880
      iz=0                                                                  880
      mnl=min(mnlam,nlam)                                                   881
10210 do 10211 m=1,nlam                                                     882
      if(flmin .lt. 1.0)goto 10231                                          882
      alm=ulam(m)                                                           882
      goto 10221                                                            883
10231 if(m .le. 2)goto 10241                                                883
      alm=alm*alf                                                           883
      goto 10221                                                            884
10241 if(m .ne. 1)goto 10251                                                884
      alm=big                                                               884
      goto 10261                                                            885
10251 continue                                                              885
      alm=0.0                                                               886
10270 do 10271 j=1,ni                                                       886
      if(ju(j).eq.0)goto 10271                                              886
      if(vp(j).le.0.0)goto 10271                                            887
      alm=max(alm,abs(g(j))/vp(j))                                          888
10271 continue                                                              889
10272 continue                                                              889
      alm=alf*alm/max(bta,1.0e-3)                                           890
10261 continue                                                              891
10221 continue                                                              891
      dem=alm*omb                                                           891
      ab=alm*bta                                                            891
      rsq0=rsq                                                              891
      jz=1                                                                  892
10280 continue                                                              892
10281 continue                                                              892
      if(iz*jz.ne.0) go to 10290                                            892
      nlp=nlp+1                                                             892
      dlx=0.0                                                               893
10300 do 10301 k=1,ni                                                       893
      if(ju(k).eq.0)goto 10301                                              894
      ak=a(k)                                                               894
      u=g(k)+ak*xv(k)                                                       894
      v=abs(u)-vp(k)*ab                                                     894
      a(k)=0.0                                                              896
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d    897 
     *em)))
      if(a(k).eq.ak)goto 10301                                              898
      if(mm(k) .ne. 0)goto 10321                                            898
      nin=nin+1                                                             898
      if(nin.gt.nx)goto 10302                                               899
10330 do 10331 j=1,ni                                                       899
      if(ju(j).eq.0)goto 10331                                              900
      if(mm(j) .eq. 0)goto 10351                                            900
      c(j,nin)=c(k,mm(j))                                                   900
      goto 10331                                                            900
10351 continue                                                              901
      if(j .ne. k)goto 10371                                                901
      c(j,nin)=xv(j)                                                        901
      goto 10331                                                            901
10371 continue                                                              902
      c(j,nin)=dot_product(x(:,j),x(:,k))                                   903
10331 continue                                                              904
10332 continue                                                              904
      mm(k)=nin                                                             904
      ia(nin)=k                                                             905
10321 continue                                                              906
      del=a(k)-ak                                                           906
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      907
      dlx=max(xv(k)*del**2,dlx)                                             908
10380 do 10381 j=1,ni                                                       908
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                               908
10381 continue                                                              909
10382 continue                                                              909
10301 continue                                                              910
10302 continue                                                              910
      if(dlx.lt.thr)goto 10282                                              910
      if(nin.gt.nx)goto 10282                                               911
      if(nlp .le. maxit)goto 10401                                          911
      jerr=-m                                                               911
      return                                                                911
10401 continue                                                              912
10290 continue                                                              912
      iz=1                                                                  912
      da(1:nin)=a(ia(1:nin))                                                913
10410 continue                                                              913
10411 continue                                                              913
      nlp=nlp+1                                                             913
      dlx=0.0                                                               914
10420 do 10421 l=1,nin                                                      914
      k=ia(l)                                                               914
      ak=a(k)                                                               914
      u=g(k)+ak*xv(k)                                                       914
      v=abs(u)-vp(k)*ab                                                     915
      a(k)=0.0                                                              917
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d    918 
     *em)))
      if(a(k).eq.ak)goto 10421                                              919
      del=a(k)-ak                                                           919
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      920
      dlx=max(xv(k)*del**2,dlx)                                             921
10430 do 10431 j=1,nin                                                      921
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                  921
10431 continue                                                              922
10432 continue                                                              922
10421 continue                                                              923
10422 continue                                                              923
      if(dlx.lt.thr)goto 10412                                              923
      if(nlp .le. maxit)goto 10451                                          923
      jerr=-m                                                               923
      return                                                                923
10451 continue                                                              924
      goto 10411                                                            925
10412 continue                                                              925
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                      926
10460 do 10461 j=1,ni                                                       926
      if(mm(j).ne.0)goto 10461                                              927
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))            928
10461 continue                                                              929
10462 continue                                                              929
      jz=0                                                                  930
      goto 10281                                                            931
10282 continue                                                              931
      if(nin .le. nx)goto 10481                                             931
      jerr=-10000-m                                                         931
      goto 10212                                                            931
10481 continue                                                              932
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 932
      kin(m)=nin                                                            933
      rsqo(m)=rsq                                                           933
      almo(m)=alm                                                           933
      lmu=m                                                                 934
      if(m.lt.mnl)goto 10211                                                934
      if(flmin.ge.1.0)goto 10211                                            935
      me=0                                                                  935
10490 do 10491 j=1,nin                                                      935
      if(ao(j,m).ne.0.0) me=me+1                                            935
10491 continue                                                              935
10492 continue                                                              935
      if(me.gt.ne)goto 10212                                                936
      if(rsq-rsq0.lt.sml*rsq)goto 10212                                     936
      if(rsq.gt.rsqmax)goto 10212                                           937
10211 continue                                                              938
10212 continue                                                              938
      deallocate(a,mm,c,da)                                                 939
      return                                                                940
      end                                                                   941
      subroutine elnetn (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam    943 
     *,thr,isd,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no),w(no),ulam(nlam),cl(2,ni)                  944
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         945
      integer jd(*),ia(nx),nin(nlam)                                        946
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          951
      allocate(xs(1:ni),stat=ierr)                                          951
      jerr=jerr+ierr                                                        952
      allocate(ju(1:ni),stat=ierr)                                          952
      jerr=jerr+ierr                                                        953
      allocate(xv(1:ni),stat=ierr)                                          953
      jerr=jerr+ierr                                                        954
      allocate(vlam(1:nlam),stat=ierr)                                      954
      jerr=jerr+ierr                                                        955
      if(jerr.ne.0) return                                                  956
      call chkvars(no,ni,x,ju)                                              957
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  958
      if(maxval(ju) .gt. 0)goto 10511                                       958
      jerr=7777                                                             958
      return                                                                958
10511 continue                                                              959
      call standard1(no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)                960
      if(jerr.ne.0) return                                                  961
      cl=cl/ys                                                              961
      if(isd .le. 0)goto 10531                                              961
10540 do 10541 j=1,ni                                                       961
      cl(:,j)=cl(:,j)*xs(j)                                                 961
10541 continue                                                              961
10542 continue                                                              961
10531 continue                                                              962
      if(flmin.ge.1.0) vlam=ulam/ys                                         963
      call elnet2(parm,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,vlam,thr,maxi    965 
     *t,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  966
10550 do 10551 k=1,lmu                                                      966
      alm(k)=ys*alm(k)                                                      966
      nk=nin(k)                                                             967
10560 do 10561 l=1,nk                                                       967
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          967
10561 continue                                                              968
10562 continue                                                              968
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         969
10551 continue                                                              970
10552 continue                                                              970
      deallocate(xm,xs,ju,xv,vlam)                                          971
      return                                                                972
      end                                                                   973
      subroutine standard1 (no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)         974
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)                        974
      integer ju(ni)                                                        975
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           978
      if(jerr.ne.0) return                                                  979
      w=w/sum(w)                                                            979
      v=sqrt(w)                                                             980
10570 do 10571 j=1,ni                                                       980
      if(ju(j).eq.0)goto 10571                                              981
      xm(j)=dot_product(w,x(:,j))                                           981
      x(:,j)=v*(x(:,j)-xm(j))                                               982
      xv(j)=dot_product(x(:,j),x(:,j))                                      982
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        983
10571 continue                                                              984
10572 continue                                                              984
      if(isd .ne. 0)goto 10591                                              984
      xs=1.0                                                                984
      goto 10601                                                            985
10591 continue                                                              985
10610 do 10611 j=1,ni                                                       985
      if(ju(j).eq.0)goto 10611                                              985
      x(:,j)=x(:,j)/xs(j)                                                   985
10611 continue                                                              986
10612 continue                                                              986
      xv=1.0                                                                987
10601 continue                                                              988
10581 continue                                                              988
      ym=dot_product(w,y)                                                   988
      y=v*(y-ym)                                                            988
      ys=sqrt(dot_product(y,y))                                             988
      y=y/ys                                                                989
      deallocate(v)                                                         990
      return                                                                991
      end                                                                   992
      subroutine elnet2(beta,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,ulam,th    994 
     *r,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    995 
     *nlam),xv(ni)
      real cl(2,ni)                                                         996
      integer ju(ni),ia(nx),kin(nlam)                                       997
      real, dimension (:), allocatable :: a,g                                   
      integer, dimension (:), allocatable :: mm,ix                              
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1002
      allocate(a(1:ni),stat=jerr)                                          1003
      allocate(mm(1:ni),stat=ierr)                                         1003
      jerr=jerr+ierr                                                       1004
      allocate(g(1:ni),stat=ierr)                                          1004
      jerr=jerr+ierr                                                       1005
      allocate(ix(1:ni),stat=ierr)                                         1005
      jerr=jerr+ierr                                                       1006
      if(jerr.ne.0) return                                                 1007
      bta=beta                                                             1007
      omb=1.0-bta                                                          1007
      ix=0                                                                 1008
      if(flmin .ge. 1.0)goto 10631                                         1008
      eqs=max(eps,flmin)                                                   1008
      alf=eqs**(1.0/(nlam-1))                                              1008
10631 continue                                                             1009
      rsq=0.0                                                              1009
      a=0.0                                                                1009
      mm=0                                                                 1009
      nlp=0                                                                1009
      nin=nlp                                                              1009
      iz=0                                                                 1009
      mnl=min(mnlam,nlam)                                                  1009
      alm=0.0                                                              1010
10640 do 10641 j=1,ni                                                      1010
      if(ju(j).eq.0)goto 10641                                             1010
      g(j)=abs(dot_product(y,x(:,j)))                                      1010
10641 continue                                                             1011
10642 continue                                                             1011
10650 do 10651 m=1,nlam                                                    1011
      alm0=alm                                                             1012
      if(flmin .lt. 1.0)goto 10671                                         1012
      alm=ulam(m)                                                          1012
      goto 10661                                                           1013
10671 if(m .le. 2)goto 10681                                               1013
      alm=alm*alf                                                          1013
      goto 10661                                                           1014
10681 if(m .ne. 1)goto 10691                                               1014
      alm=big                                                              1014
      goto 10701                                                           1015
10691 continue                                                             1015
      alm0=0.0                                                             1016
10710 do 10711 j=1,ni                                                      1016
      if(ju(j).eq.0)goto 10711                                             1016
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           1016
10711 continue                                                             1017
10712 continue                                                             1017
      alm0=alm0/max(bta,1.0e-3)                                            1017
      alm=alf*alm0                                                         1018
10701 continue                                                             1019
10661 continue                                                             1019
      dem=alm*omb                                                          1019
      ab=alm*bta                                                           1019
      rsq0=rsq                                                             1019
      jz=1                                                                 1020
      tlam=bta*(2.0*alm-alm0)                                              1021
10720 do 10721 k=1,ni                                                      1021
      if(ix(k).eq.1)goto 10721                                             1021
      if(ju(k).eq.0)goto 10721                                             1022
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       1023
10721 continue                                                             1024
10722 continue                                                             1024
10730 continue                                                             1024
10731 continue                                                             1024
      if(iz*jz.ne.0) go to 10290                                           1025
10740 continue                                                             1025
      nlp=nlp+1                                                            1025
      dlx=0.0                                                              1026
10750 do 10751 k=1,ni                                                      1026
      if(ix(k).eq.0)goto 10751                                             1026
      gk=dot_product(y,x(:,k))                                             1027
      ak=a(k)                                                              1027
      u=gk+ak*xv(k)                                                        1027
      v=abs(u)-vp(k)*ab                                                    1027
      a(k)=0.0                                                             1029
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1030 
     *em)))
      if(a(k).eq.ak)goto 10751                                             1031
      if(mm(k) .ne. 0)goto 10771                                           1031
      nin=nin+1                                                            1031
      if(nin.gt.nx)goto 10752                                              1032
      mm(k)=nin                                                            1032
      ia(nin)=k                                                            1033
10771 continue                                                             1034
      del=a(k)-ak                                                          1034
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1035
      y=y-del*x(:,k)                                                       1035
      dlx=max(xv(k)*del**2,dlx)                                            1036
10751 continue                                                             1037
10752 continue                                                             1037
      if(nin.gt.nx)goto 10732                                              1038
      if(dlx .ge. thr)goto 10791                                           1038
      ixx=0                                                                1039
10800 do 10801 k=1,ni                                                      1039
      if(ix(k).eq.1)goto 10801                                             1039
      if(ju(k).eq.0)goto 10801                                             1040
      g(k)=abs(dot_product(y,x(:,k)))                                      1041
      if(g(k) .le. ab*vp(k))goto 10821                                     1041
      ix(k)=1                                                              1041
      ixx=1                                                                1041
10821 continue                                                             1042
10801 continue                                                             1043
10802 continue                                                             1043
      if(ixx.eq.1) go to 10740                                             1044
      goto 10732                                                           1045
10791 continue                                                             1046
      if(nlp .le. maxit)goto 10841                                         1046
      jerr=-m                                                              1046
      return                                                               1046
10841 continue                                                             1047
10290 continue                                                             1047
      iz=1                                                                 1048
10850 continue                                                             1048
10851 continue                                                             1048
      nlp=nlp+1                                                            1048
      dlx=0.0                                                              1049
10860 do 10861 l=1,nin                                                     1049
      k=ia(l)                                                              1049
      gk=dot_product(y,x(:,k))                                             1050
      ak=a(k)                                                              1050
      u=gk+ak*xv(k)                                                        1050
      v=abs(u)-vp(k)*ab                                                    1050
      a(k)=0.0                                                             1052
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1053 
     *em)))
      if(a(k).eq.ak)goto 10861                                             1054
      del=a(k)-ak                                                          1054
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1055
      y=y-del*x(:,k)                                                       1055
      dlx=max(xv(k)*del**2,dlx)                                            1056
10861 continue                                                             1057
10862 continue                                                             1057
      if(dlx.lt.thr)goto 10852                                             1057
      if(nlp .le. maxit)goto 10881                                         1057
      jerr=-m                                                              1057
      return                                                               1057
10881 continue                                                             1058
      goto 10851                                                           1059
10852 continue                                                             1059
      jz=0                                                                 1060
      goto 10731                                                           1061
10732 continue                                                             1061
      if(nin .le. nx)goto 10901                                            1061
      jerr=-10000-m                                                        1061
      goto 10652                                                           1061
10901 continue                                                             1062
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1062
      kin(m)=nin                                                           1063
      rsqo(m)=rsq                                                          1063
      almo(m)=alm                                                          1063
      lmu=m                                                                1064
      if(m.lt.mnl)goto 10651                                               1064
      if(flmin.ge.1.0)goto 10651                                           1065
      me=0                                                                 1065
10910 do 10911 j=1,nin                                                     1065
      if(ao(j,m).ne.0.0) me=me+1                                           1065
10911 continue                                                             1065
10912 continue                                                             1065
      if(me.gt.ne)goto 10652                                               1066
      if(rsq-rsq0.lt.sml*rsq)goto 10652                                    1066
      if(rsq.gt.rsqmax)goto 10652                                          1067
10651 continue                                                             1068
10652 continue                                                             1068
      deallocate(a,mm,g,ix)                                                1069
      return                                                               1070
      end                                                                  1071
      subroutine chkvars(no,ni,x,ju)                                       1072
      real x(no,ni)                                                        1072
      integer ju(ni)                                                       1073
10920 do 10921 j=1,ni                                                      1073
      ju(j)=0                                                              1073
      t=x(1,j)                                                             1074
10930 do 10931 i=2,no                                                      1074
      if(x(i,j).eq.t)goto 10931                                            1074
      ju(j)=1                                                              1074
      goto 10932                                                           1074
10931 continue                                                             1075
10932 continue                                                             1075
10921 continue                                                             1076
10922 continue                                                             1076
      return                                                               1077
      end                                                                  1078
      subroutine uncomp(ni,ca,ia,nin,a)                                    1079
      real ca(*),a(ni)                                                     1079
      integer ia(*)                                                        1080
      a=0.0                                                                1080
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                  1081
      return                                                               1082
      end                                                                  1083
      subroutine modval(a0,ca,ia,nin,n,x,f)                                1084
      real ca(nin),x(n,*),f(n)                                             1084
      integer ia(nin)                                                      1085
      f=a0                                                                 1085
      if(nin.le.0) return                                                  1086
10940 do 10941 i=1,n                                                       1086
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      1086
10941 continue                                                             1087
10942 continue                                                             1087
      return                                                               1088
      end                                                                  1089
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam   1092 
     *,flmin,ulam,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)                     1093
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1094
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1095
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10961                                    1098
      jerr=10000                                                           1098
      return                                                               1098
10961 continue                                                             1099
      allocate(vq(1:ni),stat=jerr)                                         1099
      if(jerr.ne.0) return                                                 1100
      vq=max(0.0,vp)                                                       1100
      vq=vq*ni/sum(vq)                                                     1101
      if(ka .ne. 1)goto 10981                                              1102
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,u   1105 
     *lam,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10991                                                           1106
10981 continue                                                             1107
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,ul   1110 
     *am,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10991 continue                                                             1111
10971 continue                                                             1111
      deallocate(vq)                                                       1112
      return                                                               1113
      end                                                                  1114
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,f   1117 
     *lmin,ulam,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)                     1118
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1119
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1120
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                          1125
      allocate(xm(1:ni),stat=ierr)                                         1125
      jerr=jerr+ierr                                                       1126
      allocate(xs(1:ni),stat=ierr)                                         1126
      jerr=jerr+ierr                                                       1127
      allocate(ju(1:ni),stat=ierr)                                         1127
      jerr=jerr+ierr                                                       1128
      allocate(xv(1:ni),stat=ierr)                                         1128
      jerr=jerr+ierr                                                       1129
      allocate(vlam(1:nlam),stat=ierr)                                     1129
      jerr=jerr+ierr                                                       1130
      if(jerr.ne.0) return                                                 1131
      call spchkvars(no,ni,x,ix,ju)                                        1132
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1133
      if(maxval(ju) .gt. 0)goto 11011                                      1133
      jerr=7777                                                            1133
      return                                                               1133
11011 continue                                                             1134
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,jerr)      1135
      if(jerr.ne.0) return                                                 1136
      cl=cl/ys                                                             1136
      if(isd .le. 0)goto 11031                                             1136
11040 do 11041 j=1,ni                                                      1136
      cl(:,j)=cl(:,j)*xs(j)                                                1136
11041 continue                                                             1136
11042 continue                                                             1136
11031 continue                                                             1137
      if(flmin.ge.1.0) vlam=ulam/ys                                        1138
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla   1140 
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1141
11050 do 11051 k=1,lmu                                                     1141
      alm(k)=ys*alm(k)                                                     1141
      nk=nin(k)                                                            1142
11060 do 11061 l=1,nk                                                      1142
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1142
11061 continue                                                             1143
11062 continue                                                             1143
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                        1144
11051 continue                                                             1145
11052 continue                                                             1145
      deallocate(xm,xs,g,ju,xv,vlam)                                       1146
      return                                                               1147
      end                                                                  1148
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,j   1149 
     *err)
      real x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                     1149
      integer ix(*),jx(*),ju(ni)                                           1150
      w=w/sum(w)                                                           1151
11070 do 11071 j=1,ni                                                      1151
      if(ju(j).eq.0)goto 11071                                             1152
      jb=ix(j)                                                             1152
      je=ix(j+1)-1                                                         1152
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1153
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1154
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1155
11071 continue                                                             1156
11072 continue                                                             1156
      if(isd .ne. 0)goto 11091                                             1156
      xs=1.0                                                               1156
      goto 11101                                                           1156
11091 continue                                                             1156
      xv=1.0                                                               1156
11101 continue                                                             1157
11081 continue                                                             1157
      ym=dot_product(w,y)                                                  1157
      y=y-ym                                                               1157
      ys=sqrt(dot_product(w,y**2))                                         1157
      y=y/ys                                                               1157
      g=0.0                                                                1158
11110 do 11111 j=1,ni                                                      1158
      if(ju(j).eq.0)goto 11111                                             1158
      jb=ix(j)                                                             1158
      je=ix(j+1)-1                                                         1159
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)           1160
11111 continue                                                             1161
11112 continue                                                             1161
      return                                                               1162
      end                                                                  1163
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm   1165 
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      real g(ni),vp(ni),x(*),ulam(nlam),w(no)                              1166
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni),cl(2,n   1167 
     *i)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1168
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1175
      allocate(a(1:ni),stat=ierr)                                          1175
      jerr=jerr+ierr                                                       1176
      allocate(mm(1:ni),stat=ierr)                                         1176
      jerr=jerr+ierr                                                       1177
      allocate(da(1:ni),stat=ierr)                                         1177
      jerr=jerr+ierr                                                       1178
      if(jerr.ne.0) return                                                 1179
      bta=beta                                                             1179
      omb=1.0-bta                                                          1180
      if(flmin .ge. 1.0)goto 11131                                         1180
      eqs=max(eps,flmin)                                                   1180
      alf=eqs**(1.0/(nlam-1))                                              1180
11131 continue                                                             1181
      rsq=0.0                                                              1181
      a=0.0                                                                1181
      mm=0                                                                 1181
      nlp=0                                                                1181
      nin=nlp                                                              1181
      iz=0                                                                 1181
      mnl=min(mnlam,nlam)                                                  1182
11140 do 11141 m=1,nlam                                                    1183
      if(flmin .lt. 1.0)goto 11161                                         1183
      alm=ulam(m)                                                          1183
      goto 11151                                                           1184
11161 if(m .le. 2)goto 11171                                               1184
      alm=alm*alf                                                          1184
      goto 11151                                                           1185
11171 if(m .ne. 1)goto 11181                                               1185
      alm=big                                                              1185
      goto 11191                                                           1186
11181 continue                                                             1186
      alm=0.0                                                              1187
11200 do 11201 j=1,ni                                                      1187
      if(ju(j).eq.0)goto 11201                                             1187
      if(vp(j).le.0.0)goto 11201                                           1188
      alm=max(alm,abs(g(j))/vp(j))                                         1189
11201 continue                                                             1190
11202 continue                                                             1190
      alm=alf*alm/max(bta,1.0e-3)                                          1191
11191 continue                                                             1192
11151 continue                                                             1192
      dem=alm*omb                                                          1192
      ab=alm*bta                                                           1192
      rsq0=rsq                                                             1192
      jz=1                                                                 1193
11210 continue                                                             1193
11211 continue                                                             1193
      if(iz*jz.ne.0) go to 10290                                           1193
      nlp=nlp+1                                                            1193
      dlx=0.0                                                              1194
11220 do 11221 k=1,ni                                                      1194
      if(ju(k).eq.0)goto 11221                                             1195
      ak=a(k)                                                              1195
      u=g(k)+ak*xv(k)                                                      1195
      v=abs(u)-vp(k)*ab                                                    1195
      a(k)=0.0                                                             1197
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1198 
     *em)))
      if(a(k).eq.ak)goto 11221                                             1199
      if(mm(k) .ne. 0)goto 11241                                           1199
      nin=nin+1                                                            1199
      if(nin.gt.nx)goto 11222                                              1200
11250 do 11251 j=1,ni                                                      1200
      if(ju(j).eq.0)goto 11251                                             1201
      if(mm(j) .eq. 0)goto 11271                                           1201
      c(j,nin)=c(k,mm(j))                                                  1201
      goto 11251                                                           1201
11271 continue                                                             1202
      if(j .ne. k)goto 11291                                               1202
      c(j,nin)=xv(j)                                                       1202
      goto 11251                                                           1202
11291 continue                                                             1203
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))       1205
11251 continue                                                             1206
11252 continue                                                             1206
      mm(k)=nin                                                            1206
      ia(nin)=k                                                            1207
11241 continue                                                             1208
      del=a(k)-ak                                                          1208
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1209
      dlx=max(xv(k)*del**2,dlx)                                            1210
11300 do 11301 j=1,ni                                                      1210
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                              1210
11301 continue                                                             1211
11302 continue                                                             1211
11221 continue                                                             1212
11222 continue                                                             1212
      if(dlx.lt.thr)goto 11212                                             1212
      if(nin.gt.nx)goto 11212                                              1213
      if(nlp .le. maxit)goto 11321                                         1213
      jerr=-m                                                              1213
      return                                                               1213
11321 continue                                                             1214
10290 continue                                                             1214
      iz=1                                                                 1214
      da(1:nin)=a(ia(1:nin))                                               1215
11330 continue                                                             1215
11331 continue                                                             1215
      nlp=nlp+1                                                            1215
      dlx=0.0                                                              1216
11340 do 11341 l=1,nin                                                     1216
      k=ia(l)                                                              1217
      ak=a(k)                                                              1217
      u=g(k)+ak*xv(k)                                                      1217
      v=abs(u)-vp(k)*ab                                                    1217
      a(k)=0.0                                                             1219
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1220 
     *em)))
      if(a(k).eq.ak)goto 11341                                             1221
      del=a(k)-ak                                                          1221
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1222
      dlx=max(xv(k)*del**2,dlx)                                            1223
11350 do 11351 j=1,nin                                                     1223
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                 1223
11351 continue                                                             1224
11352 continue                                                             1224
11341 continue                                                             1225
11342 continue                                                             1225
      if(dlx.lt.thr)goto 11332                                             1225
      if(nlp .le. maxit)goto 11371                                         1225
      jerr=-m                                                              1225
      return                                                               1225
11371 continue                                                             1226
      goto 11331                                                           1227
11332 continue                                                             1227
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                     1228
11380 do 11381 j=1,ni                                                      1228
      if(mm(j).ne.0)goto 11381                                             1229
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))           1230
11381 continue                                                             1231
11382 continue                                                             1231
      jz=0                                                                 1232
      goto 11211                                                           1233
11212 continue                                                             1233
      if(nin .le. nx)goto 11401                                            1233
      jerr=-10000-m                                                        1233
      goto 11142                                                           1233
11401 continue                                                             1234
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1234
      kin(m)=nin                                                           1235
      rsqo(m)=rsq                                                          1235
      almo(m)=alm                                                          1235
      lmu=m                                                                1236
      if(m.lt.mnl)goto 11141                                               1236
      if(flmin.ge.1.0)goto 11141                                           1237
      me=0                                                                 1237
11410 do 11411 j=1,nin                                                     1237
      if(ao(j,m).ne.0.0) me=me+1                                           1237
11411 continue                                                             1237
11412 continue                                                             1237
      if(me.gt.ne)goto 11142                                               1238
      if(rsq-rsq0.lt.sml*rsq)goto 11142                                    1238
      if(rsq.gt.rsqmax)goto 11142                                          1239
11141 continue                                                             1240
11142 continue                                                             1240
      deallocate(a,mm,c,da)                                                1241
      return                                                               1242
      end                                                                  1243
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flm   1245 
     *in,ulam,  thr,isd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),vp(ni),y(no),w(no),ulam(nlam),cl(2,ni)                     1246
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1247
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1248
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         1253
      allocate(xs(1:ni),stat=ierr)                                         1253
      jerr=jerr+ierr                                                       1254
      allocate(ju(1:ni),stat=ierr)                                         1254
      jerr=jerr+ierr                                                       1255
      allocate(xv(1:ni),stat=ierr)                                         1255
      jerr=jerr+ierr                                                       1256
      allocate(vlam(1:nlam),stat=ierr)                                     1256
      jerr=jerr+ierr                                                       1257
      if(jerr.ne.0) return                                                 1258
      call spchkvars(no,ni,x,ix,ju)                                        1259
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1260
      if(maxval(ju) .gt. 0)goto 11431                                      1260
      jerr=7777                                                            1260
      return                                                               1260
11431 continue                                                             1261
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,jerr)       1262
      if(jerr.ne.0) return                                                 1263
      cl=cl/ys                                                             1263
      if(isd .le. 0)goto 11451                                             1263
11460 do 11461 j=1,ni                                                      1263
      cl(:,j)=cl(:,j)*xs(j)                                                1263
11461 continue                                                             1263
11462 continue                                                             1263
11451 continue                                                             1264
      if(flmin.ge.1.0) vlam=ulam/ys                                        1265
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla   1267 
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1268
11470 do 11471 k=1,lmu                                                     1268
      alm(k)=ys*alm(k)                                                     1268
      nk=nin(k)                                                            1269
11480 do 11481 l=1,nk                                                      1269
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1269
11481 continue                                                             1270
11482 continue                                                             1270
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                        1271
11471 continue                                                             1272
11472 continue                                                             1272
      deallocate(xm,xs,ju,xv,vlam)                                         1273
      return                                                               1274
      end                                                                  1275
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,je   1276 
     *rr)
      real x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)                           1276
      integer ix(*),jx(*),ju(ni)                                           1277
      w=w/sum(w)                                                           1278
11490 do 11491 j=1,ni                                                      1278
      if(ju(j).eq.0)goto 11491                                             1279
      jb=ix(j)                                                             1279
      je=ix(j+1)-1                                                         1279
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1280
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1281
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1282
11491 continue                                                             1283
11492 continue                                                             1283
      if(isd .ne. 0)goto 11511                                             1283
      xs=1.0                                                               1283
      goto 11521                                                           1283
11511 continue                                                             1283
      xv=1.0                                                               1283
11521 continue                                                             1284
11501 continue                                                             1284
      ym=dot_product(w,y)                                                  1284
      y=y-ym                                                               1284
      ys=sqrt(dot_product(w,y**2))                                         1284
      y=y/ys                                                               1285
      return                                                               1286
      end                                                                  1287
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm   1289 
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      real y(no),w(no),x(*),vp(ni),ulam(nlam),cl(2,ni)                     1290
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)          1291
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1292
      real, dimension (:), allocatable :: a,g                                   
      integer, dimension (:), allocatable :: mm,iy                              
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1297
      allocate(a(1:ni),stat=jerr)                                          1298
      allocate(mm(1:ni),stat=ierr)                                         1298
      jerr=jerr+ierr                                                       1299
      allocate(g(1:ni),stat=ierr)                                          1299
      jerr=jerr+ierr                                                       1300
      allocate(iy(1:ni),stat=ierr)                                         1300
      jerr=jerr+ierr                                                       1301
      if(jerr.ne.0) return                                                 1302
      bta=beta                                                             1302
      omb=1.0-bta                                                          1302
      alm=0.0                                                              1302
      iy=0                                                                 1303
      if(flmin .ge. 1.0)goto 11541                                         1303
      eqs=max(eps,flmin)                                                   1303
      alf=eqs**(1.0/(nlam-1))                                              1303
11541 continue                                                             1304
      rsq=0.0                                                              1304
      a=0.0                                                                1304
      mm=0                                                                 1304
      o=0.0                                                                1304
      nlp=0                                                                1304
      nin=nlp                                                              1304
      iz=0                                                                 1304
      mnl=min(mnlam,nlam)                                                  1305
11550 do 11551 j=1,ni                                                      1305
      if(ju(j).eq.0)goto 11551                                             1306
      jb=ix(j)                                                             1306
      je=ix(j+1)-1                                                         1307
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1308
11551 continue                                                             1309
11552 continue                                                             1309
11560 do 11561 m=1,nlam                                                    1309
      alm0=alm                                                             1310
      if(flmin .lt. 1.0)goto 11581                                         1310
      alm=ulam(m)                                                          1310
      goto 11571                                                           1311
11581 if(m .le. 2)goto 11591                                               1311
      alm=alm*alf                                                          1311
      goto 11571                                                           1312
11591 if(m .ne. 1)goto 11601                                               1312
      alm=big                                                              1312
      goto 11611                                                           1313
11601 continue                                                             1313
      alm0=0.0                                                             1314
11620 do 11621 j=1,ni                                                      1314
      if(ju(j).eq.0)goto 11621                                             1314
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           1314
11621 continue                                                             1315
11622 continue                                                             1315
      alm0=alm0/max(bta,1.0e-3)                                            1315
      alm=alf*alm0                                                         1316
11611 continue                                                             1317
11571 continue                                                             1317
      dem=alm*omb                                                          1317
      ab=alm*bta                                                           1317
      rsq0=rsq                                                             1317
      jz=1                                                                 1318
      tlam=bta*(2.0*alm-alm0)                                              1319
11630 do 11631 k=1,ni                                                      1319
      if(iy(k).eq.1)goto 11631                                             1319
      if(ju(k).eq.0)goto 11631                                             1320
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       1321
11631 continue                                                             1322
11632 continue                                                             1322
11640 continue                                                             1322
11641 continue                                                             1322
      if(iz*jz.ne.0) go to 10290                                           1323
10740 continue                                                             1323
      nlp=nlp+1                                                            1323
      dlx=0.0                                                              1324
11650 do 11651 k=1,ni                                                      1324
      if(iy(k).eq.0)goto 11651                                             1324
      jb=ix(k)                                                             1324
      je=ix(k+1)-1                                                         1325
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1326
      ak=a(k)                                                              1326
      u=gk+ak*xv(k)                                                        1326
      v=abs(u)-vp(k)*ab                                                    1326
      a(k)=0.0                                                             1328
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1329 
     *em)))
      if(a(k).eq.ak)goto 11651                                             1330
      if(mm(k) .ne. 0)goto 11671                                           1330
      nin=nin+1                                                            1330
      if(nin.gt.nx)goto 11652                                              1331
      mm(k)=nin                                                            1331
      ia(nin)=k                                                            1332
11671 continue                                                             1333
      del=a(k)-ak                                                          1333
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1334
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1335
      o=o+del*xm(k)/xs(k)                                                  1335
      dlx=max(xv(k)*del**2,dlx)                                            1336
11651 continue                                                             1337
11652 continue                                                             1337
      if(nin.gt.nx)goto 11642                                              1338
      if(dlx .ge. thr)goto 11691                                           1338
      ixx=0                                                                1339
11700 do 11701 j=1,ni                                                      1339
      if(iy(j).eq.1)goto 11701                                             1339
      if(ju(j).eq.0)goto 11701                                             1340
      jb=ix(j)                                                             1340
      je=ix(j+1)-1                                                         1341
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1342
      if(g(j) .le. ab*vp(j))goto 11721                                     1342
      iy(j)=1                                                              1342
      ixx=1                                                                1342
11721 continue                                                             1343
11701 continue                                                             1344
11702 continue                                                             1344
      if(ixx.eq.1) go to 10740                                             1345
      goto 11642                                                           1346
11691 continue                                                             1347
      if(nlp .le. maxit)goto 11741                                         1347
      jerr=-m                                                              1347
      return                                                               1347
11741 continue                                                             1348
10290 continue                                                             1348
      iz=1                                                                 1349
11750 continue                                                             1349
11751 continue                                                             1349
      nlp=nlp+1                                                            1349
      dlx=0.0                                                              1350
11760 do 11761 l=1,nin                                                     1350
      k=ia(l)                                                              1350
      jb=ix(k)                                                             1350
      je=ix(k+1)-1                                                         1351
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1352
      ak=a(k)                                                              1352
      u=gk+ak*xv(k)                                                        1352
      v=abs(u)-vp(k)*ab                                                    1352
      a(k)=0.0                                                             1354
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1355 
     *em)))
      if(a(k).eq.ak)goto 11761                                             1356
      del=a(k)-ak                                                          1356
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1357
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1358
      o=o+del*xm(k)/xs(k)                                                  1358
      dlx=max(xv(k)*del**2,dlx)                                            1359
11761 continue                                                             1360
11762 continue                                                             1360
      if(dlx.lt.thr)goto 11752                                             1360
      if(nlp .le. maxit)goto 11781                                         1360
      jerr=-m                                                              1360
      return                                                               1360
11781 continue                                                             1361
      goto 11751                                                           1362
11752 continue                                                             1362
      jz=0                                                                 1363
      goto 11641                                                           1364
11642 continue                                                             1364
      if(nin .le. nx)goto 11801                                            1364
      jerr=-10000-m                                                        1364
      goto 11562                                                           1364
11801 continue                                                             1365
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1365
      kin(m)=nin                                                           1366
      rsqo(m)=rsq                                                          1366
      almo(m)=alm                                                          1366
      lmu=m                                                                1367
      if(m.lt.mnl)goto 11561                                               1367
      if(flmin.ge.1.0)goto 11561                                           1368
      me=0                                                                 1368
11810 do 11811 j=1,nin                                                     1368
      if(ao(j,m).ne.0.0) me=me+1                                           1368
11811 continue                                                             1368
11812 continue                                                             1368
      if(me.gt.ne)goto 11562                                               1369
      if(rsq-rsq0.lt.sml*rsq)goto 11562                                    1369
      if(rsq.gt.rsqmax)goto 11562                                          1370
11561 continue                                                             1371
11562 continue                                                             1371
      deallocate(a,mm,g,iy)                                                1372
      return                                                               1373
      end                                                                  1374
      subroutine spchkvars(no,ni,x,ix,ju)                                  1375
      real x(*)                                                            1375
      integer ix(*),ju(ni)                                                 1376
11820 do 11821 j=1,ni                                                      1376
      ju(j)=0                                                              1376
      jb=ix(j)                                                             1376
      nj=ix(j+1)-jb                                                        1376
      if(nj.eq.0)goto 11821                                                1377
      je=ix(j+1)-1                                                         1378
      if(nj .ge. no)goto 11841                                             1378
11850 do 11851 i=jb,je                                                     1378
      if(x(i).eq.0.0)goto 11851                                            1378
      ju(j)=1                                                              1378
      goto 11852                                                           1378
11851 continue                                                             1378
11852 continue                                                             1378
      goto 11861                                                           1379
11841 continue                                                             1379
      t=x(jb)                                                              1379
11870 do 11871 i=jb+1,je                                                   1379
      if(x(i).eq.t)goto 11871                                              1379
      ju(j)=1                                                              1379
      goto 11872                                                           1379
11871 continue                                                             1379
11872 continue                                                             1379
11861 continue                                                             1380
11831 continue                                                             1380
11821 continue                                                             1381
11822 continue                                                             1381
      return                                                               1382
      end                                                                  1383
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                         1384
      real ca(*),x(*),f(n)                                                 1384
      integer ia(*),ix(*),jx(*)                                            1385
      f=a0                                                                 1386
11880 do 11881 j=1,nin                                                     1386
      k=ia(j)                                                              1386
      kb=ix(k)                                                             1386
      ke=ix(k+1)-1                                                         1387
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                             1388
11881 continue                                                             1389
11882 continue                                                             1389
      return                                                               1390
      end                                                                  1391
      function row_prod(i,j,ia,ja,ra,w)                                    1392
      integer ia(*),ja(*)                                                  1392
      real ra(*),w(*)                                                      1393
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(   1395 
     *i),ia(j+1)-ia(j),w)
      return                                                               1396
      end                                                                  1397
      function dot(x,y,mx,my,nx,ny,w)                                      1398
      real x(*),y(*),w(*)                                                  1398
      integer mx(*),my(*)                                                  1399
      i=1                                                                  1399
      j=i                                                                  1399
      s=0.0                                                                1400
11890 continue                                                             1400
11891 continue                                                             1400
11900 continue                                                             1401
11901 if(mx(i).ge.my(j))goto 11902                                         1401
      i=i+1                                                                1401
      if(i.gt.nx) go to 11910                                              1401
      goto 11901                                                           1402
11902 continue                                                             1402
      if(mx(i).eq.my(j)) go to 11920                                       1403
11930 continue                                                             1403
11931 if(my(j).ge.mx(i))goto 11932                                         1403
      j=j+1                                                                1403
      if(j.gt.ny) go to 11910                                              1403
      goto 11931                                                           1404
11932 continue                                                             1404
      if(mx(i).eq.my(j)) go to 11920                                       1404
      goto 11891                                                           1405
11920 continue                                                             1405
      s=s+w(mx(i))*x(i)*y(j)                                               1406
      i=i+1                                                                1406
      if(i.gt.nx)goto 11892                                                1406
      j=j+1                                                                1406
      if(j.gt.ny)goto 11892                                                1407
      goto 11891                                                           1408
11892 continue                                                             1408
11910 continue                                                             1408
      dot=s                                                                1409
      return                                                               1410
      end                                                                  1411
      subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,cl,ne,nx,nlam,flmin,u   1413 
     *lam,thr,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)             1414
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni)         1415
      integer jd(*),ia(nx),nin(nlam)                                       1416
      real, dimension (:), allocatable :: xm,xs,ww,vq,xv                        
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 11951                                    1420
      jerr=10000                                                           1420
      return                                                               1420
11951 continue                                                             1421
      allocate(ww(1:no),stat=jerr)                                         1422
      allocate(ju(1:ni),stat=ierr)                                         1422
      jerr=jerr+ierr                                                       1423
      allocate(vq(1:ni),stat=ierr)                                         1423
      jerr=jerr+ierr                                                       1424
      allocate(xm(1:ni),stat=ierr)                                         1424
      jerr=jerr+ierr                                                       1425
      if(kopt .ne. 2)goto 11971                                            1425
      allocate(xv(1:ni),stat=ierr)                                         1425
      jerr=jerr+ierr                                                       1425
11971 continue                                                             1426
      if(isd .le. 0)goto 11991                                             1426
      allocate(xs(1:ni),stat=ierr)                                         1426
      jerr=jerr+ierr                                                       1426
11991 continue                                                             1427
      if(jerr.ne.0) return                                                 1428
      call chkvars(no,ni,x,ju)                                             1429
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1430
      if(maxval(ju) .gt. 0)goto 12011                                      1430
      jerr=7777                                                            1430
      return                                                               1430
12011 continue                                                             1431
      vq=max(0.0,vp)                                                       1431
      vq=vq*ni/sum(vq)                                                     1432
12020 do 12021 i=1,no                                                      1432
      ww(i)=sum(y(i,:))                                                    1432
      y(i,:)=y(i,:)/ww(i)                                                  1432
12021 continue                                                             1432
12022 continue                                                             1432
      sw=sum(ww)                                                           1432
      ww=ww/sw                                                             1433
      if(nc .ne. 1)goto 12041                                              1433
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             1434
      if(isd .le. 0)goto 12061                                             1434
12070 do 12071 j=1,ni                                                      1434
      cl(:,j)=cl(:,j)*xs(j)                                                1434
12071 continue                                                             1434
12072 continue                                                             1434
12061 continue                                                             1435
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,cl,ne,nx,nlam,fl   1437 
     *min,ulam,  thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      goto 12031                                                           1438
12041 if(kopt .ne. 2)goto 12081                                            1438
      call multlstandard1(no,ni,x,ww,ju,isd,xm,xs,xv)                      1439
      if(isd .le. 0)goto 12101                                             1439
12110 do 12111 j=1,ni                                                      1439
      cl(:,j)=cl(:,j)*xs(j)                                                1439
12111 continue                                                             1439
12112 continue                                                             1439
12101 continue                                                             1440
      call multlognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,   1442 
     *ulam,thr,  isd,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 12121                                                           1443
12081 continue                                                             1443
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             1444
      if(isd .le. 0)goto 12141                                             1444
12150 do 12151 j=1,ni                                                      1444
      cl(:,j)=cl(:,j)*xs(j)                                                1444
12151 continue                                                             1444
12152 continue                                                             1444
12141 continue                                                             1445
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam   1447 
     *,thr,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
12121 continue                                                             1448
12031 continue                                                             1448
      if(jerr.gt.0) return                                                 1448
      dev0=2.0*sw*dev0                                                     1449
12160 do 12161 k=1,lmu                                                     1449
      nk=nin(k)                                                            1450
12170 do 12171 ic=1,nc                                                     1450
      if(isd .le. 0)goto 12191                                             1450
12200 do 12201 l=1,nk                                                      1450
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1450
12201 continue                                                             1450
12202 continue                                                             1450
12191 continue                                                             1451
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1452
12171 continue                                                             1453
12172 continue                                                             1453
12161 continue                                                             1454
12162 continue                                                             1454
      deallocate(ww,ju,vq,xm)                                              1454
      if(isd.gt.0) deallocate(xs)                                          1455
      if(kopt.eq.2) deallocate(xv)                                         1456
      return                                                               1457
      end                                                                  1458
      subroutine lstandard1 (no,ni,x,w,ju,isd,xm,xs)                       1459
      real x(no,ni),w(no),xm(ni),xs(ni)                                    1459
      integer ju(ni)                                                       1460
12210 do 12211 j=1,ni                                                      1460
      if(ju(j).eq.0)goto 12211                                             1461
      xm(j)=dot_product(w,x(:,j))                                          1461
      x(:,j)=x(:,j)-xm(j)                                                  1462
      if(isd .le. 0)goto 12231                                             1462
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 1462
      x(:,j)=x(:,j)/xs(j)                                                  1462
12231 continue                                                             1463
12211 continue                                                             1464
12212 continue                                                             1464
      return                                                               1465
      end                                                                  1466
      subroutine multlstandard1 (no,ni,x,w,ju,isd,xm,xs,xv)                1467
      real x(no,ni),w(no),xm(ni),xs(ni),xv(ni)                             1467
      integer ju(ni)                                                       1468
12240 do 12241 j=1,ni                                                      1468
      if(ju(j).eq.0)goto 12241                                             1469
      xm(j)=dot_product(w,x(:,j))                                          1469
      x(:,j)=x(:,j)-xm(j)                                                  1470
      xv(j)=dot_product(w,x(:,j)**2)                                       1471
      if(isd .le. 0)goto 12261                                             1471
      xs(j)=sqrt(xv(j))                                                    1471
      x(:,j)=x(:,j)/xs(j)                                                  1471
      xv(j)=1.0                                                            1471
12261 continue                                                             1472
12241 continue                                                             1473
12242 continue                                                             1473
      return                                                               1474
      end                                                                  1475
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin,u   1477 
     *lam,shri,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)           1478
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1479
      integer ju(ni),m(nx),kin(nlam)                                       1480
      real, dimension (:), allocatable :: b,bs,v,r,xv,q,ga                      
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               1485
      allocate(b(0:ni),stat=jerr)                                          1486
      allocate(xv(1:ni),stat=ierr)                                         1486
      jerr=jerr+ierr                                                       1487
      allocate(ga(1:ni),stat=ierr)                                         1487
      jerr=jerr+ierr                                                       1488
      allocate(bs(0:ni),stat=ierr)                                         1488
      jerr=jerr+ierr                                                       1489
      allocate(mm(1:ni),stat=ierr)                                         1489
      jerr=jerr+ierr                                                       1490
      allocate(ixx(1:ni),stat=ierr)                                        1490
      jerr=jerr+ierr                                                       1491
      allocate(r(1:no),stat=ierr)                                          1491
      jerr=jerr+ierr                                                       1492
      allocate(v(1:no),stat=ierr)                                          1492
      jerr=jerr+ierr                                                       1493
      allocate(q(1:no),stat=ierr)                                          1493
      jerr=jerr+ierr                                                       1494
      if(jerr.ne.0) return                                                 1495
      fmax=log(1.0/pmin-1.0)                                               1495
      fmin=-fmax                                                           1495
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1496
      bta=parm                                                             1496
      omb=1.0-bta                                                          1497
      q0=dot_product(w,y)                                                  1497
      if(q0 .gt. pmin)goto 12281                                           1497
      jerr=8001                                                            1497
      return                                                               1497
12281 continue                                                             1498
      if(q0 .lt. 1.0-pmin)goto 12301                                       1498
      jerr=9001                                                            1498
      return                                                               1498
12301 continue                                                             1499
      ixx=0                                                                1499
      al=0.0                                                               1499
      bz=log(q0/(1.0-q0))                                                  1500
      if(nonzero(no,g) .ne. 0)goto 12321                                   1500
      vi=q0*(1.0-q0)                                                       1500
      b(0)=bz                                                              1500
      v=vi*w                                                               1501
      r=w*(y-q0)                                                           1501
      q=q0                                                                 1501
      xmz=vi                                                               1501
      dev1=-(bz*q0+log(1.0-q0))                                            1502
      goto 12331                                                           1503
12321 continue                                                             1503
      b(0)=azero(no,y,g,w,jerr)                                            1503
      if(jerr.ne.0) return                                                 1504
      q=1.0/(1.0+exp(-b(0)-g))                                             1504
      v=w*q*(1.0-q)                                                        1504
      r=w*(y-q)                                                            1504
      xmz=sum(v)                                                           1505
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1506
12331 continue                                                             1507
12311 continue                                                             1507
      if(kopt .le. 0)goto 12351                                            1508
      if(isd .le. 0)goto 12371                                             1508
      xv=0.25                                                              1508
      goto 12381                                                           1509
12371 continue                                                             1509
12390 do 12391 j=1,ni                                                      1509
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1509
12391 continue                                                             1509
12392 continue                                                             1509
12381 continue                                                             1510
12361 continue                                                             1510
12351 continue                                                             1511
      dev0=dev1                                                            1512
12400 do 12401 i=1,no                                                      1512
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1513
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1514
12401 continue                                                             1515
12402 continue                                                             1515
      if(flmin .ge. 1.0)goto 12421                                         1515
      eqs=max(eps,flmin)                                                   1515
      alf=eqs**(1.0/(nlam-1))                                              1515
12421 continue                                                             1516
      m=0                                                                  1516
      mm=0                                                                 1516
      nlp=0                                                                1516
      nin=nlp                                                              1516
      mnl=min(mnlam,nlam)                                                  1516
      bs=0.0                                                               1516
      b(1:ni)=0.0                                                          1517
      shr=shri*dev0                                                        1518
12430 do 12431 j=1,ni                                                      1518
      if(ju(j).eq.0)goto 12431                                             1518
      ga(j)=abs(dot_product(r,x(:,j)))                                     1518
12431 continue                                                             1519
12432 continue                                                             1519
12440 do 12441 ilm=1,nlam                                                  1519
      al0=al                                                               1520
      if(flmin .lt. 1.0)goto 12461                                         1520
      al=ulam(ilm)                                                         1520
      goto 12451                                                           1521
12461 if(ilm .le. 2)goto 12471                                             1521
      al=al*alf                                                            1521
      goto 12451                                                           1522
12471 if(ilm .ne. 1)goto 12481                                             1522
      al=big                                                               1522
      goto 12491                                                           1523
12481 continue                                                             1523
      al0=0.0                                                              1524
12500 do 12501 j=1,ni                                                      1524
      if(ju(j).eq.0)goto 12501                                             1524
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1524
12501 continue                                                             1525
12502 continue                                                             1525
      al0=al0/max(bta,1.0e-3)                                              1525
      al=alf*al0                                                           1526
12491 continue                                                             1527
12451 continue                                                             1527
      al2=al*omb                                                           1527
      al1=al*bta                                                           1527
      tlam=bta*(2.0*al-al0)                                                1528
12510 do 12511 k=1,ni                                                      1528
      if(ixx(k).eq.1)goto 12511                                            1528
      if(ju(k).eq.0)goto 12511                                             1529
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1530
12511 continue                                                             1531
12512 continue                                                             1531
10740 continue                                                             1532
12520 continue                                                             1532
12521 continue                                                             1532
      bs(0)=b(0)                                                           1532
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1533
      if(kopt .ne. 0)goto 12541                                            1534
12550 do 12551 j=1,ni                                                      1534
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                       1534
12551 continue                                                             1535
12552 continue                                                             1535
12541 continue                                                             1536
12560 continue                                                             1536
12561 continue                                                             1536
      nlp=nlp+1                                                            1536
      dlx=0.0                                                              1537
12570 do 12571 k=1,ni                                                      1537
      if(ixx(k).eq.0)goto 12571                                            1538
      bk=b(k)                                                              1538
      gk=dot_product(r,x(:,k))                                             1539
      u=gk+xv(k)*b(k)                                                      1539
      au=abs(u)-vp(k)*al1                                                  1540
      if(au .gt. 0.0)goto 12591                                            1540
      b(k)=0.0                                                             1540
      goto 12601                                                           1541
12591 continue                                                             1542
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          1543
12601 continue                                                             1544
12581 continue                                                             1544
      d=b(k)-bk                                                            1544
      if(abs(d).le.0.0)goto 12571                                          1544
      dlx=max(dlx,xv(k)*d**2)                                              1545
      r=r-d*v*x(:,k)                                                       1546
      if(mm(k) .ne. 0)goto 12621                                           1546
      nin=nin+1                                                            1546
      if(nin.gt.nx)goto 12572                                              1547
      mm(k)=nin                                                            1547
      m(nin)=k                                                             1548
12621 continue                                                             1549
12571 continue                                                             1550
12572 continue                                                             1550
      if(nin.gt.nx)goto 12562                                              1551
      d=sum(r)/xmz                                                         1552
      if(d .eq. 0.0)goto 12641                                             1552
      b(0)=b(0)+d                                                          1552
      dlx=max(dlx,xmz*d**2)                                                1552
      r=r-d*v                                                              1552
12641 continue                                                             1553
      if(dlx.lt.shr)goto 12562                                             1553
      if(nlp .le. maxit)goto 12661                                         1553
      jerr=-ilm                                                            1553
      return                                                               1553
12661 continue                                                             1554
12670 continue                                                             1554
12671 continue                                                             1554
      nlp=nlp+1                                                            1554
      dlx=0.0                                                              1555
12680 do 12681 l=1,nin                                                     1555
      k=m(l)                                                               1555
      bk=b(k)                                                              1556
      gk=dot_product(r,x(:,k))                                             1557
      u=gk+xv(k)*b(k)                                                      1557
      au=abs(u)-vp(k)*al1                                                  1558
      if(au .gt. 0.0)goto 12701                                            1558
      b(k)=0.0                                                             1558
      goto 12711                                                           1559
12701 continue                                                             1560
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          1561
12711 continue                                                             1562
12691 continue                                                             1562
      d=b(k)-bk                                                            1562
      if(abs(d).le.0.0)goto 12681                                          1562
      dlx=max(dlx,xv(k)*d**2)                                              1563
      r=r-d*v*x(:,k)                                                       1564
12681 continue                                                             1565
12682 continue                                                             1565
      d=sum(r)/xmz                                                         1566
      if(d .eq. 0.0)goto 12731                                             1566
      b(0)=b(0)+d                                                          1566
      dlx=max(dlx,xmz*d**2)                                                1566
      r=r-d*v                                                              1566
12731 continue                                                             1567
      if(dlx.lt.shr)goto 12672                                             1567
      if(nlp .le. maxit)goto 12751                                         1567
      jerr=-ilm                                                            1567
      return                                                               1567
12751 continue                                                             1568
      goto 12671                                                           1569
12672 continue                                                             1569
      goto 12561                                                           1570
12562 continue                                                             1570
      if(nin.gt.nx)goto 12522                                              1571
12760 do 12761 i=1,no                                                      1571
      fi=b(0)+g(i)                                                         1572
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1573
      if(fi .ge. fmin)goto 12781                                           1573
      q(i)=0.0                                                             1573
      goto 12771                                                           1573
12781 if(fi .le. fmax)goto 12791                                           1573
      q(i)=1.0                                                             1573
      goto 12801                                                           1574
12791 continue                                                             1574
      q(i)=1.0/(1.0+exp(-fi))                                              1574
12801 continue                                                             1575
12771 continue                                                             1575
12761 continue                                                             1576
12762 continue                                                             1576
      v=w*q*(1.0-q)                                                        1576
      xmz=sum(v)                                                           1576
      if(xmz.le.vmin)goto 12522                                            1576
      r=w*(y-q)                                                            1577
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 12821                           1577
      ix=0                                                                 1578
12830 do 12831 j=1,nin                                                     1578
      k=m(j)                                                               1579
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 12831                           1579
      ix=1                                                                 1579
      goto 12832                                                           1580
12831 continue                                                             1581
12832 continue                                                             1581
      if(ix .ne. 0)goto 12851                                              1582
12860 do 12861 k=1,ni                                                      1582
      if(ixx(k).eq.1)goto 12861                                            1582
      if(ju(k).eq.0)goto 12861                                             1583
      ga(k)=abs(dot_product(r,x(:,k)))                                     1584
      if(ga(k) .le. al1*vp(k))goto 12881                                   1584
      ixx(k)=1                                                             1584
      ix=1                                                                 1584
12881 continue                                                             1585
12861 continue                                                             1586
12862 continue                                                             1586
      if(ix.eq.1) go to 10740                                              1587
      goto 12522                                                           1588
12851 continue                                                             1589
12821 continue                                                             1590
      goto 12521                                                           1591
12522 continue                                                             1591
      if(nin .le. nx)goto 12901                                            1591
      jerr=-10000-ilm                                                      1591
      goto 12442                                                           1591
12901 continue                                                             1592
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1592
      kin(ilm)=nin                                                         1593
      a0(ilm)=b(0)                                                         1593
      alm(ilm)=al                                                          1593
      lmu=ilm                                                              1594
      devi=dev2(no,w,y,q,pmin)                                             1595
      dev(ilm)=(dev1-devi)/dev0                                            1595
      if(xmz.le.vmin)goto 12442                                            1596
      if(ilm.lt.mnl)goto 12441                                             1596
      if(flmin.ge.1.0)goto 12441                                           1597
      me=0                                                                 1597
12910 do 12911 j=1,nin                                                     1597
      if(a(j,ilm).ne.0.0) me=me+1                                          1597
12911 continue                                                             1597
12912 continue                                                             1597
      if(me.gt.ne)goto 12442                                               1598
      if(dev(ilm).gt.devmax)goto 12442                                     1598
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12442                             1599
12441 continue                                                             1600
12442 continue                                                             1600
      g=log(q/(1.0-q))                                                     1601
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                                  1602
      return                                                               1603
      end                                                                  1604
      function dev2(n,w,y,p,pmin)                                          1605
      real w(n),y(n),p(n)                                                  1606
      pmax=1.0-pmin                                                        1606
      s=0.0                                                                1607
12920 do 12921 i=1,n                                                       1607
      pi=min(max(pmin,p(i)),pmax)                                          1608
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1609
12921 continue                                                             1610
12922 continue                                                             1610
      dev2=s                                                               1611
      return                                                               1612
      end                                                                  1613
      function azero(n,y,g,q,jerr)                                         1614
      parameter(eps=1.0e-7)                                                1615
      real y(n),g(n),q(n)                                                  1616
      real, dimension (:), allocatable :: e,p,w                                 
      allocate(e(1:n),stat=jerr)                                           1620
      allocate(p(1:n),stat=ierr)                                           1620
      jerr=jerr+ierr                                                       1621
      allocate(w(1:n),stat=ierr)                                           1621
      jerr=jerr+ierr                                                       1622
      if(jerr.ne.0) return                                                 1623
      az=0.0                                                               1623
      e=exp(-g)                                                            1623
      qy=dot_product(q,y)                                                  1623
      p=1.0/(1.0+e)                                                        1624
12930 continue                                                             1624
12931 continue                                                             1624
      w=q*p*(1.0-p)                                                        1625
      d=(qy-dot_product(q,p))/sum(w)                                       1625
      az=az+d                                                              1625
      if(abs(d).lt.eps)goto 12932                                          1626
      ea0=exp(-az)                                                         1626
      p=1.0/(1.0+ea0*e)                                                    1627
      goto 12931                                                           1628
12932 continue                                                             1628
      azero=az                                                             1629
      deallocate(e,p,w)                                                    1630
      return                                                               1631
      end                                                                  1632
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin   1634 
     *,ulam,shri,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam)              1635
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni)          1636
      integer ju(ni),m(nx),kin(nlam)                                       1637
      real, dimension (:,:), allocatable :: q                                   
      real, dimension (:), allocatable :: sxp,sxpl                              
      real, dimension (:), allocatable :: di,v,r,ga                             
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is,ixx                          
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               1648
      exmn=-exmx                                                           1649
      allocate(r(1:no),stat=ierr)                                          1649
      jerr=jerr+ierr                                                       1650
      allocate(v(1:no),stat=ierr)                                          1650
      jerr=jerr+ierr                                                       1651
      allocate(mm(1:ni),stat=ierr)                                         1651
      jerr=jerr+ierr                                                       1652
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1652
      jerr=jerr+ierr                                                       1653
      allocate(sxp(1:no),stat=ierr)                                        1653
      jerr=jerr+ierr                                                       1654
      allocate(sxpl(1:no),stat=ierr)                                       1654
      jerr=jerr+ierr                                                       1655
      allocate(di(1:no),stat=ierr)                                         1655
      jerr=jerr+ierr                                                       1656
      allocate(ga(1:ni),stat=ierr)                                         1656
      jerr=jerr+ierr                                                       1657
      allocate(ixx(1:ni),stat=ierr)                                        1657
      jerr=jerr+ierr                                                       1658
      if(jerr.ne.0) return                                                 1659
      pmax=1.0-pmin                                                        1659
      emin=pmin/pmax                                                       1659
      emax=1.0/emin                                                        1660
      pfm=(1.0+pmin)*pmin                                                  1660
      pfx=(1.0-pmin)*pmax                                                  1660
      vmin=pfm*pmax                                                        1661
      bta=parm                                                             1661
      omb=1.0-bta                                                          1661
      dev1=0.0                                                             1661
      dev0=0.0                                                             1662
12940 do 12941 ic=1,nc                                                     1662
      q0=dot_product(w,y(:,ic))                                            1663
      if(q0 .gt. pmin)goto 12961                                           1663
      jerr =8000+ic                                                        1663
      return                                                               1663
12961 continue                                                             1664
      if(q0 .lt. 1.0-pmin)goto 12981                                       1664
      jerr =9000+ic                                                        1664
      return                                                               1664
12981 continue                                                             1665
      b(0,ic)=log(q0)                                                      1665
      dev1=dev1-q0*b(0,ic)                                                 1665
      b(1:ni,ic)=0.0                                                       1666
12941 continue                                                             1667
12942 continue                                                             1667
      ixx=0                                                                1667
      al=0.0                                                               1668
      if(nonzero(no*nc,g) .ne. 0)goto 13001                                1669
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1669
      sxp=0.0                                                              1670
13010 do 13011 ic=1,nc                                                     1670
      q(:,ic)=exp(b(0,ic))                                                 1670
      sxp=sxp+q(:,ic)                                                      1670
13011 continue                                                             1671
13012 continue                                                             1671
      goto 13021                                                           1672
13001 continue                                                             1672
13030 do 13031 i=1,no                                                      1672
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1672
13031 continue                                                             1672
13032 continue                                                             1672
      sxp=0.0                                                              1673
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1673
      if(jerr.ne.0) return                                                 1674
      dev1=0.0                                                             1675
13040 do 13041 ic=1,nc                                                     1675
      q(:,ic)=b(0,ic)+g(:,ic)                                              1676
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             1677
      q(:,ic)=exp(q(:,ic))                                                 1677
      sxp=sxp+q(:,ic)                                                      1678
13041 continue                                                             1679
13042 continue                                                             1679
      sxpl=w*log(sxp)                                                      1679
13050 do 13051 ic=1,nc                                                     1679
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  1679
13051 continue                                                             1680
13052 continue                                                             1680
13021 continue                                                             1681
12991 continue                                                             1681
13060 do 13061 ic=1,nc                                                     1681
13070 do 13071 i=1,no                                                      1681
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1681
13071 continue                                                             1681
13072 continue                                                             1681
13061 continue                                                             1682
13062 continue                                                             1682
      dev0=dev0+dev1                                                       1683
      if(kopt .le. 0)goto 13091                                            1684
      if(isd .le. 0)goto 13111                                             1684
      xv=0.25                                                              1684
      goto 13121                                                           1685
13111 continue                                                             1685
13130 do 13131 j=1,ni                                                      1685
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1685
13131 continue                                                             1685
13132 continue                                                             1685
13121 continue                                                             1686
13101 continue                                                             1686
13091 continue                                                             1687
      if(flmin .ge. 1.0)goto 13151                                         1687
      eqs=max(eps,flmin)                                                   1687
      alf=eqs**(1.0/(nlam-1))                                              1687
13151 continue                                                             1688
      m=0                                                                  1688
      mm=0                                                                 1688
      nin=0                                                                1688
      nlp=0                                                                1688
      mnl=min(mnlam,nlam)                                                  1688
      bs=0.0                                                               1688
      shr=shri*dev0                                                        1689
      ga=0.0                                                               1690
13160 do 13161 ic=1,nc                                                     1690
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1691
13170 do 13171 j=1,ni                                                      1691
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))           1691
13171 continue                                                             1692
13172 continue                                                             1692
13161 continue                                                             1693
13162 continue                                                             1693
13180 do 13181 ilm=1,nlam                                                  1693
      al0=al                                                               1694
      if(flmin .lt. 1.0)goto 13201                                         1694
      al=ulam(ilm)                                                         1694
      goto 13191                                                           1695
13201 if(ilm .le. 2)goto 13211                                             1695
      al=al*alf                                                            1695
      goto 13191                                                           1696
13211 if(ilm .ne. 1)goto 13221                                             1696
      al=big                                                               1696
      goto 13231                                                           1697
13221 continue                                                             1697
      al0=0.0                                                              1698
13240 do 13241 j=1,ni                                                      1698
      if(ju(j).eq.0)goto 13241                                             1698
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1698
13241 continue                                                             1699
13242 continue                                                             1699
      al0=al0/max(bta,1.0e-3)                                              1699
      al=alf*al0                                                           1700
13231 continue                                                             1701
13191 continue                                                             1701
      al2=al*omb                                                           1701
      al1=al*bta                                                           1701
      tlam=bta*(2.0*al-al0)                                                1702
13250 do 13251 k=1,ni                                                      1702
      if(ixx(k).eq.1)goto 13251                                            1702
      if(ju(k).eq.0)goto 13251                                             1703
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1704
13251 continue                                                             1705
13252 continue                                                             1705
10740 continue                                                             1706
13260 continue                                                             1706
13261 continue                                                             1706
      ix=0                                                                 1706
      jx=ix                                                                1706
      ig=0                                                                 1707
13270 do 13271 ic=1,nc                                                     1707
      bs(0,ic)=b(0,ic)                                                     1708
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1709
      xmz=0.0                                                              1710
13280 do 13281 i=1,no                                                      1710
      pic=q(i,ic)/sxp(i)                                                   1711
      if(pic .ge. pfm)goto 13301                                           1711
      pic=0.0                                                              1711
      v(i)=0.0                                                             1711
      goto 13291                                                           1712
13301 if(pic .le. pfx)goto 13311                                           1712
      pic=1.0                                                              1712
      v(i)=0.0                                                             1712
      goto 13321                                                           1713
13311 continue                                                             1713
      v(i)=w(i)*pic*(1.0-pic)                                              1713
      xmz=xmz+v(i)                                                         1713
13321 continue                                                             1714
13291 continue                                                             1714
      r(i)=w(i)*(y(i,ic)-pic)                                              1715
13281 continue                                                             1716
13282 continue                                                             1716
      if(xmz.le.vmin)goto 13271                                            1716
      ig=1                                                                 1717
      if(kopt .ne. 0)goto 13341                                            1718
13350 do 13351 j=1,ni                                                      1718
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                    1718
13351 continue                                                             1719
13352 continue                                                             1719
13341 continue                                                             1720
13360 continue                                                             1720
13361 continue                                                             1720
      nlp=nlp+1                                                            1720
      dlx=0.0                                                              1721
13370 do 13371 k=1,ni                                                      1721
      if(ixx(k).eq.0)goto 13371                                            1722
      bk=b(k,ic)                                                           1722
      gk=dot_product(r,x(:,k))                                             1723
      u=gk+xv(k,ic)*b(k,ic)                                                1723
      au=abs(u)-vp(k)*al1                                                  1724
      if(au .gt. 0.0)goto 13391                                            1724
      b(k,ic)=0.0                                                          1724
      goto 13401                                                           1725
13391 continue                                                             1726
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   1728 
     *)
13401 continue                                                             1729
13381 continue                                                             1729
      d=b(k,ic)-bk                                                         1729
      if(abs(d).le.0.0)goto 13371                                          1730
      dlx=max(dlx,xv(k,ic)*d**2)                                           1730
      r=r-d*v*x(:,k)                                                       1731
      if(mm(k) .ne. 0)goto 13421                                           1731
      nin=nin+1                                                            1732
      if(nin .le. nx)goto 13441                                            1732
      jx=1                                                                 1732
      goto 13372                                                           1732
13441 continue                                                             1733
      mm(k)=nin                                                            1733
      m(nin)=k                                                             1734
13421 continue                                                             1735
13371 continue                                                             1736
13372 continue                                                             1736
      if(jx.gt.0)goto 13362                                                1737
      d=sum(r)/xmz                                                         1738
      if(d .eq. 0.0)goto 13461                                             1738
      b(0,ic)=b(0,ic)+d                                                    1738
      dlx=max(dlx,xmz*d**2)                                                1738
      r=r-d*v                                                              1738
13461 continue                                                             1739
      if(dlx.lt.shr)goto 13362                                             1740
      if(nlp .le. maxit)goto 13481                                         1740
      jerr=-ilm                                                            1740
      return                                                               1740
13481 continue                                                             1741
13490 continue                                                             1741
13491 continue                                                             1741
      nlp=nlp+1                                                            1741
      dlx=0.0                                                              1742
13500 do 13501 l=1,nin                                                     1742
      k=m(l)                                                               1742
      bk=b(k,ic)                                                           1743
      gk=dot_product(r,x(:,k))                                             1744
      u=gk+xv(k,ic)*b(k,ic)                                                1744
      au=abs(u)-vp(k)*al1                                                  1745
      if(au .gt. 0.0)goto 13521                                            1745
      b(k,ic)=0.0                                                          1745
      goto 13531                                                           1746
13521 continue                                                             1747
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   1749 
     *)
13531 continue                                                             1750
13511 continue                                                             1750
      d=b(k,ic)-bk                                                         1750
      if(abs(d).le.0.0)goto 13501                                          1751
      dlx=max(dlx,xv(k,ic)*d**2)                                           1751
      r=r-d*v*x(:,k)                                                       1752
13501 continue                                                             1753
13502 continue                                                             1753
      d=sum(r)/xmz                                                         1754
      if(d .eq. 0.0)goto 13551                                             1754
      b(0,ic)=b(0,ic)+d                                                    1755
      dlx=max(dlx,xmz*d**2)                                                1755
      r=r-d*v                                                              1756
13551 continue                                                             1757
      if(dlx.lt.shr)goto 13492                                             1757
      if(nlp .le. maxit)goto 13571                                         1757
      jerr=-ilm                                                            1757
      return                                                               1757
13571 continue                                                             1758
      goto 13491                                                           1759
13492 continue                                                             1759
      goto 13361                                                           1760
13362 continue                                                             1760
      if(jx.gt.0)goto 13272                                                1761
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                            1762
      if(ix .ne. 0)goto 13591                                              1763
13600 do 13601 j=1,nin                                                     1763
      k=m(j)                                                               1764
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 13621                1764
      ix=1                                                                 1764
      goto 13602                                                           1764
13621 continue                                                             1765
13601 continue                                                             1766
13602 continue                                                             1766
13591 continue                                                             1767
13630 do 13631 i=1,no                                                      1767
      fi=b(0,ic)+g(i,ic)                                                   1769
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1770
      fi=min(max(exmn,fi),exmx)                                            1770
      sxp(i)=sxp(i)-q(i,ic)                                                1771
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    1772
      sxp(i)=sxp(i)+q(i,ic)                                                1773
13631 continue                                                             1774
13632 continue                                                             1774
13271 continue                                                             1775
13272 continue                                                             1775
      s=-sum(b(0,:))/nc                                                    1775
      b(0,:)=b(0,:)+s                                                      1775
      di=s                                                                 1776
13640 do 13641 j=1,nin                                                     1776
      l=m(j)                                                               1777
      if(vp(l) .gt. 0.0)goto 13661                                         1777
      s=sum(b(l,:))/nc                                                     1777
      goto 13671                                                           1778
13661 continue                                                             1778
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                     1778
13671 continue                                                             1779
13651 continue                                                             1779
      b(l,:)=b(l,:)-s                                                      1779
      di=di-s*x(:,l)                                                       1780
13641 continue                                                             1781
13642 continue                                                             1781
      di=exp(di)                                                           1781
      sxp=sxp*di                                                           1781
13680 do 13681 ic=1,nc                                                     1781
      q(:,ic)=q(:,ic)*di                                                   1781
13681 continue                                                             1782
13682 continue                                                             1782
      if(jx.gt.0)goto 13262                                                1782
      if(ig.eq.0)goto 13262                                                1783
      if(ix .ne. 0)goto 13701                                              1784
13710 do 13711 k=1,ni                                                      1784
      if(ixx(k).eq.1)goto 13711                                            1784
      if(ju(k).eq.0)goto 13711                                             1784
      ga(k)=0.0                                                            1784
13711 continue                                                             1785
13712 continue                                                             1785
13720 do 13721 ic=1,nc                                                     1785
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1786
13730 do 13731 k=1,ni                                                      1786
      if(ixx(k).eq.1)goto 13731                                            1786
      if(ju(k).eq.0)goto 13731                                             1787
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                          1788
13731 continue                                                             1789
13732 continue                                                             1789
13721 continue                                                             1790
13722 continue                                                             1790
13740 do 13741 k=1,ni                                                      1790
      if(ixx(k).eq.1)goto 13741                                            1790
      if(ju(k).eq.0)goto 13741                                             1791
      if(ga(k) .le. al1*vp(k))goto 13761                                   1791
      ixx(k)=1                                                             1791
      ix=1                                                                 1791
13761 continue                                                             1792
13741 continue                                                             1793
13742 continue                                                             1793
      if(ix.eq.1) go to 10740                                              1794
      goto 13262                                                           1795
13701 continue                                                             1796
      goto 13261                                                           1797
13262 continue                                                             1797
      if(jx .le. 0)goto 13781                                              1797
      jerr=-10000-ilm                                                      1797
      goto 13182                                                           1797
13781 continue                                                             1797
      devi=0.0                                                             1798
13790 do 13791 ic=1,nc                                                     1799
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1799
      a0(ic,ilm)=b(0,ic)                                                   1800
13800 do 13801 i=1,no                                                      1800
      if(y(i,ic).le.0.0)goto 13801                                         1801
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1802
13801 continue                                                             1803
13802 continue                                                             1803
13791 continue                                                             1804
13792 continue                                                             1804
      kin(ilm)=nin                                                         1804
      alm(ilm)=al                                                          1804
      lmu=ilm                                                              1805
      dev(ilm)=(dev1-devi)/dev0                                            1805
      if(ig.eq.0)goto 13182                                                1806
      if(ilm.lt.mnl)goto 13181                                             1806
      if(flmin.ge.1.0)goto 13181                                           1807
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 13182             1808
      if(dev(ilm).gt.devmax)goto 13182                                     1808
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13182                             1809
13181 continue                                                             1810
13182 continue                                                             1810
      g=log(q)                                                             1810
13810 do 13811 i=1,no                                                      1810
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1810
13811 continue                                                             1811
13812 continue                                                             1811
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                           1812
      return                                                               1813
      end                                                                  1814
      subroutine kazero(kk,n,y,g,q,az,jerr)                                1815
      parameter(eps=1.0e-7)                                                1816
      real y(n,kk),g(n,kk),q(n),az(kk)                                     1817
      real, dimension (:), allocatable :: s                                     
      real, dimension (:,:), allocatable :: e                                   
      allocate(e(1:n,1:kk),stat=jerr)                                           
      allocate(s(1:n),stat=ierr)                                           1822
      jerr=jerr+ierr                                                       1823
      if(jerr.ne.0) return                                                 1824
      az=0.0                                                               1824
      e=exp(g)                                                             1824
13820 do 13821 i=1,n                                                       1824
      s(i)=sum(e(i,:))                                                     1824
13821 continue                                                             1825
13822 continue                                                             1825
13830 continue                                                             1825
13831 continue                                                             1825
      dm=0.0                                                               1826
13840 do 13841 k=1,kk                                                      1826
      t=0.0                                                                1826
      u=t                                                                  1827
13850 do 13851 i=1,n                                                       1827
      pik=e(i,k)/s(i)                                                      1828
      t=t+q(i)*(y(i,k)-pik)                                                1828
      u=u+q(i)*pik*(1.0-pik)                                               1829
13851 continue                                                             1830
13852 continue                                                             1830
      d=t/u                                                                1830
      az(k)=az(k)+d                                                        1830
      ed=exp(d)                                                            1830
      dm=max(dm,abs(d))                                                    1831
13860 do 13861 i=1,n                                                       1831
      z=e(i,k)                                                             1831
      e(i,k)=z*ed                                                          1831
      s(i)=s(i)-z+e(i,k)                                                   1831
13861 continue                                                             1832
13862 continue                                                             1832
13841 continue                                                             1833
13842 continue                                                             1833
      if(dm.lt.eps)goto 13832                                              1833
      goto 13831                                                           1834
13832 continue                                                             1834
      az=az-sum(az)/kk                                                     1835
      deallocate(e,s)                                                      1836
      return                                                               1837
      end                                                                  1838
      function elc(parm,n,cl,a,m)                                          1839
      real a(n),cl(2)                                                      1839
      integer m(n)                                                         1840
      fn=n                                                                 1840
      am=sum(a)/fn                                                         1841
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 13881                       1841
      elc=am                                                               1841
      go to 13890                                                          1841
13881 continue                                                             1842
13900 do 13901 i=1,n                                                       1842
      m(i)=i                                                               1842
13901 continue                                                             1842
13902 continue                                                             1842
      call psort7(a,m,1,n)                                                 1843
      if(a(m(1)) .ne. a(m(n)))goto 13921                                   1843
      elc=a(1)                                                             1843
      go to 13890                                                          1843
13921 continue                                                             1844
      if(mod(n,2) .ne. 1)goto 13941                                        1844
      ad=a(m(n/2+1))                                                       1844
      goto 13951                                                           1845
13941 continue                                                             1845
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       1845
13951 continue                                                             1846
13931 continue                                                             1846
      if(parm .ne. 1.0)goto 13971                                          1846
      elc=ad                                                               1846
      go to 13890                                                          1846
13971 continue                                                             1847
      b1=min(am,ad)                                                        1847
      b2=max(am,ad)                                                        1847
      k2=1                                                                 1848
13980 continue                                                             1848
13981 if(a(m(k2)).gt.b1)goto 13982                                         1848
      k2=k2+1                                                              1848
      goto 13981                                                           1848
13982 continue                                                             1848
      k1=k2-1                                                              1849
13990 continue                                                             1849
13991 if(a(m(k2)).ge.b2)goto 13992                                         1849
      k2=k2+1                                                              1849
      goto 13991                                                           1850
13992 continue                                                             1850
      r=parm/((1.0-parm)*fn)                                               1850
      is=0                                                                 1850
      sm=n-2*(k1-1)                                                        1851
14000 do 14001 k=k1,k2-1                                                   1851
      sm=sm-2.0                                                            1851
      s=r*sm+am                                                            1852
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 14021                   1852
      is=k                                                                 1852
      goto 14002                                                           1852
14021 continue                                                             1853
14001 continue                                                             1854
14002 continue                                                             1854
      if(is .eq. 0)goto 14041                                              1854
      elc=s                                                                1854
      go to 13890                                                          1854
14041 continue                                                             1854
      r2=2.0*r                                                             1854
      s1=a(m(k1))                                                          1854
      am2=2.0*am                                                           1855
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    1855
      elc=s1                                                               1856
14050 do 14051 k=k1+1,k2                                                   1856
      s=a(m(k))                                                            1856
      if(s.eq.s1)goto 14051                                                1857
      c=r2*sum(abs(a-s))+s*(s-am2)                                         1858
      if(c .ge. cri)goto 14071                                             1858
      cri=c                                                                1858
      elc=s                                                                1858
14071 continue                                                             1858
      s1=s                                                                 1859
14051 continue                                                             1860
14052 continue                                                             1860
13890 continue                                                             1860
      elc=max(maxval(a-cl(2)),min(minval(a-cl(1)),elc))                    1861
      return                                                               1862
      end                                                                  1863
      function nintot(ni,nx,nc,a,m,nin,is)                                 1864
      real a(nx,nc)                                                        1864
      integer m(nx),is(ni)                                                 1865
      is=0                                                                 1865
      nintot=0                                                             1866
14080 do 14081 ic=1,nc                                                     1866
14090 do 14091 j=1,nin                                                     1866
      k=m(j)                                                               1866
      if(is(k).ne.0)goto 14091                                             1867
      if(a(j,ic).eq.0.0)goto 14091                                         1867
      is(k)=k                                                              1867
      nintot=nintot+1                                                      1868
14091 continue                                                             1868
14092 continue                                                             1868
14081 continue                                                             1869
14082 continue                                                             1869
      return                                                               1870
      end                                                                  1871
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             1872
      real ca(nx,nc),a(ni,nc)                                              1872
      integer ia(nx)                                                       1873
      a=0.0                                                                1874
14100 do 14101 ic=1,nc                                                     1874
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            1874
14101 continue                                                             1875
14102 continue                                                             1875
      return                                                               1876
      end                                                                  1877
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      1878
      real a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                             1878
      integer ia(nx)                                                       1879
14110 do 14111 i=1,nt                                                      1879
14120 do 14121 ic=1,nc                                                     1879
      ans(ic,i)=a0(ic)                                                     1881
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   1882 
     *:nin)))
14121 continue                                                             1882
14122 continue                                                             1882
14111 continue                                                             1883
14112 continue                                                             1883
      return                                                               1884
      end                                                                  1885
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,cl,ne,nx,nlam   1887 
     *,flmin,  ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp
     *,jerr)
      real x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)                 1888
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni)         1889
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1890
      real, dimension (:), allocatable :: xm,xs,ww,vq,xv                        
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 14141                                    1894
      jerr=10000                                                           1894
      return                                                               1894
14141 continue                                                             1895
      allocate(ww(1:no),stat=jerr)                                         1896
      allocate(ju(1:ni),stat=ierr)                                         1896
      jerr=jerr+ierr                                                       1897
      allocate(vq(1:ni),stat=ierr)                                         1897
      jerr=jerr+ierr                                                       1898
      allocate(xm(1:ni),stat=ierr)                                         1898
      jerr=jerr+ierr                                                       1899
      allocate(xs(1:ni),stat=ierr)                                         1899
      jerr=jerr+ierr                                                       1900
      if(kopt .ne. 2)goto 14161                                            1900
      allocate(xv(1:ni),stat=ierr)                                         1900
      jerr=jerr+ierr                                                       1900
14161 continue                                                             1901
      if(jerr.ne.0) return                                                 1902
      call spchkvars(no,ni,x,ix,ju)                                        1903
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1904
      if(maxval(ju) .gt. 0)goto 14181                                      1904
      jerr=7777                                                            1904
      return                                                               1904
14181 continue                                                             1905
      vq=max(0.0,vp)                                                       1905
      vq=vq*ni/sum(vq)                                                     1906
14190 do 14191 i=1,no                                                      1906
      ww(i)=sum(y(i,:))                                                    1906
      y(i,:)=y(i,:)/ww(i)                                                  1906
14191 continue                                                             1906
14192 continue                                                             1906
      sw=sum(ww)                                                           1906
      ww=ww/sw                                                             1907
      if(nc .ne. 1)goto 14211                                              1907
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     1908
      if(isd .le. 0)goto 14231                                             1908
14240 do 14241 j=1,ni                                                      1908
      cl(:,j)=cl(:,j)*xs(j)                                                1908
14241 continue                                                             1908
14242 continue                                                             1908
14231 continue                                                             1909
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,cl,ne,n   1912 
     *x,nlam,  flmin,ulam,thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0
     *,dev,  alm,nlp,jerr)
      goto 14201                                                           1913
14211 if(kopt .ne. 2)goto 14251                                            1913
      call multsplstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs,xv)              1914
      if(isd .le. 0)goto 14271                                             1914
14280 do 14281 j=1,ni                                                      1914
      cl(:,j)=cl(:,j)*xs(j)                                                1914
14281 continue                                                             1914
14282 continue                                                             1914
14271 continue                                                             1915
      call multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nl   1917 
     *am,flmin,  ulam,thr,isd,maxit,xv,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,a
     *lm,nlp,jerr)
      goto 14291                                                           1918
14251 continue                                                             1918
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     1919
      if(isd .le. 0)goto 14311                                             1919
14320 do 14321 j=1,ni                                                      1919
      cl(:,j)=cl(:,j)*xs(j)                                                1919
14321 continue                                                             1919
14322 continue                                                             1919
14311 continue                                                             1920
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,f   1922 
     *lmin,  ulam,thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm
     *,nlp,jerr)
14291 continue                                                             1923
14201 continue                                                             1923
      if(jerr.gt.0) return                                                 1923
      dev0=2.0*sw*dev0                                                     1924
14330 do 14331 k=1,lmu                                                     1924
      nk=nin(k)                                                            1925
14340 do 14341 ic=1,nc                                                     1925
      if(isd .le. 0)goto 14361                                             1925
14370 do 14371 l=1,nk                                                      1925
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1925
14371 continue                                                             1925
14372 continue                                                             1925
14361 continue                                                             1926
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1927
14341 continue                                                             1928
14342 continue                                                             1928
14331 continue                                                             1929
14332 continue                                                             1929
      deallocate(ww,ju,vq,xm,xs)                                           1929
      if(kopt.eq.2) deallocate(xv)                                         1930
      return                                                               1931
      end                                                                  1932
      subroutine multsplstandard2(no,ni,x,ix,jx,w,ju,isd,xm,xs,xv)         1933
      real x(*),w(no),xm(ni),xs(ni),xv(ni)                                 1933
      integer ix(*),jx(*),ju(ni)                                           1934
14380 do 14381 j=1,ni                                                      1934
      if(ju(j).eq.0)goto 14381                                             1934
      jb=ix(j)                                                             1934
      je=ix(j+1)-1                                                         1935
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1936
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1937
      if(isd .le. 0)goto 14401                                             1937
      xs(j)=sqrt(xv(j))                                                    1937
      xv(j)=1.0                                                            1937
14401 continue                                                             1938
14381 continue                                                             1939
14382 continue                                                             1939
      if(isd.eq.0) xs=1.0                                                  1940
      return                                                               1941
      end                                                                  1942
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,xm,xs)                1943
      real x(*),w(no),xm(ni),xs(ni)                                        1943
      integer ix(*),jx(*),ju(ni)                                           1944
14410 do 14411 j=1,ni                                                      1944
      if(ju(j).eq.0)goto 14411                                             1944
      jb=ix(j)                                                             1944
      je=ix(j+1)-1                                                         1945
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1946
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   1947 
     *)**2)
14411 continue                                                             1948
14412 continue                                                             1948
      if(isd.eq.0) xs=1.0                                                  1949
      return                                                               1950
      end                                                                  1951
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,nl   1953 
     *am,  flmin,ulam,shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,
     *alm,nlp,jerr)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)               1954
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1955
      real xb(ni),xs(ni)                                                   1955
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1956
      real, dimension (:), allocatable :: xm,b,bs,v,r,sc,xv,q,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               1961
      allocate(b(0:ni),stat=jerr)                                          1962
      allocate(xm(0:ni),stat=ierr)                                         1962
      jerr=jerr+ierr                                                       1963
      allocate(xv(1:ni),stat=ierr)                                         1963
      jerr=jerr+ierr                                                       1964
      allocate(bs(0:ni),stat=ierr)                                         1964
      jerr=jerr+ierr                                                       1965
      allocate(ga(1:ni),stat=ierr)                                         1965
      jerr=jerr+ierr                                                       1966
      allocate(mm(1:ni),stat=ierr)                                         1966
      jerr=jerr+ierr                                                       1967
      allocate(ixx(1:ni),stat=ierr)                                        1967
      jerr=jerr+ierr                                                       1968
      allocate(q(1:no),stat=ierr)                                          1968
      jerr=jerr+ierr                                                       1969
      allocate(r(1:no),stat=ierr)                                          1969
      jerr=jerr+ierr                                                       1970
      allocate(v(1:no),stat=ierr)                                          1970
      jerr=jerr+ierr                                                       1971
      allocate(sc(1:no),stat=ierr)                                         1971
      jerr=jerr+ierr                                                       1972
      if(jerr.ne.0) return                                                 1973
      fmax=log(1.0/pmin-1.0)                                               1973
      fmin=-fmax                                                           1973
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1974
      bta=parm                                                             1974
      omb=1.0-bta                                                          1975
      q0=dot_product(w,y)                                                  1975
      if(q0 .gt. pmin)goto 14431                                           1975
      jerr=8001                                                            1975
      return                                                               1975
14431 continue                                                             1976
      if(q0 .lt. 1.0-pmin)goto 14451                                       1976
      jerr=9001                                                            1976
      return                                                               1976
14451 continue                                                             1976
      bz=log(q0/(1.0-q0))                                                  1977
      if(nonzero(no,g) .ne. 0)goto 14471                                   1977
      vi=q0*(1.0-q0)                                                       1977
      b(0)=bz                                                              1977
      v=vi*w                                                               1978
      r=w*(y-q0)                                                           1978
      q=q0                                                                 1978
      xm(0)=vi                                                             1978
      dev1=-(bz*q0+log(1.0-q0))                                            1979
      goto 14481                                                           1980
14471 continue                                                             1980
      b(0)=azero(no,y,g,w,jerr)                                            1980
      if(jerr.ne.0) return                                                 1981
      q=1.0/(1.0+exp(-b(0)-g))                                             1981
      v=w*q*(1.0-q)                                                        1981
      r=w*(y-q)                                                            1981
      xm(0)=sum(v)                                                         1982
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1983
14481 continue                                                             1984
14461 continue                                                             1984
      if(kopt .le. 0)goto 14501                                            1985
      if(isd .le. 0)goto 14521                                             1985
      xv=0.25                                                              1985
      goto 14531                                                           1986
14521 continue                                                             1987
14540 do 14541 j=1,ni                                                      1987
      if(ju(j).eq.0)goto 14541                                             1987
      jb=ix(j)                                                             1987
      je=ix(j+1)-1                                                         1988
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          1989
14541 continue                                                             1990
14542 continue                                                             1990
14531 continue                                                             1991
14511 continue                                                             1991
14501 continue                                                             1992
      b(1:ni)=0.0                                                          1992
      dev0=dev1                                                            1993
14550 do 14551 i=1,no                                                      1993
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1994
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1995
14551 continue                                                             1996
14552 continue                                                             1996
      if(flmin .ge. 1.0)goto 14571                                         1996
      eqs=max(eps,flmin)                                                   1996
      alf=eqs**(1.0/(nlam-1))                                              1996
14571 continue                                                             1997
      m=0                                                                  1997
      mm=0                                                                 1997
      nin=0                                                                1997
      o=0.0                                                                1997
      svr=o                                                                1997
      mnl=min(mnlam,nlam)                                                  1997
      bs=0.0                                                               1997
      nlp=0                                                                1997
      nin=nlp                                                              1998
      shr=shri*dev0                                                        1998
      al=0.0                                                               1998
      ixx=0                                                                1999
14580 do 14581 j=1,ni                                                      1999
      if(ju(j).eq.0)goto 14581                                             2000
      jb=ix(j)                                                             2000
      je=ix(j+1)-1                                                         2000
      jn=ix(j+1)-ix(j)                                                     2001
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2002
      gj=dot_product(sc(1:jn),x(jb:je))                                    2003
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      2004
14581 continue                                                             2005
14582 continue                                                             2005
14590 do 14591 ilm=1,nlam                                                  2005
      al0=al                                                               2006
      if(flmin .lt. 1.0)goto 14611                                         2006
      al=ulam(ilm)                                                         2006
      goto 14601                                                           2007
14611 if(ilm .le. 2)goto 14621                                             2007
      al=al*alf                                                            2007
      goto 14601                                                           2008
14621 if(ilm .ne. 1)goto 14631                                             2008
      al=big                                                               2008
      goto 14641                                                           2009
14631 continue                                                             2009
      al0=0.0                                                              2010
14650 do 14651 j=1,ni                                                      2010
      if(ju(j).eq.0)goto 14651                                             2010
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2010
14651 continue                                                             2011
14652 continue                                                             2011
      al0=al0/max(bta,1.0e-3)                                              2011
      al=alf*al0                                                           2012
14641 continue                                                             2013
14601 continue                                                             2013
      al2=al*omb                                                           2013
      al1=al*bta                                                           2013
      tlam=bta*(2.0*al-al0)                                                2014
14660 do 14661 k=1,ni                                                      2014
      if(ixx(k).eq.1)goto 14661                                            2014
      if(ju(k).eq.0)goto 14661                                             2015
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2016
14661 continue                                                             2017
14662 continue                                                             2017
10740 continue                                                             2018
14670 continue                                                             2018
14671 continue                                                             2018
      bs(0)=b(0)                                                           2018
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                2019
14680 do 14681 j=1,ni                                                      2019
      if(ixx(j).eq.0)goto 14681                                            2020
      jb=ix(j)                                                             2020
      je=ix(j+1)-1                                                         2020
      jn=ix(j+1)-ix(j)                                                     2021
      sc(1:jn)=v(jx(jb:je))                                                2022
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 2023
      if(kopt .ne. 0)goto 14701                                            2024
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              2025
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                2026
14701 continue                                                             2027
14681 continue                                                             2028
14682 continue                                                             2028
14710 continue                                                             2028
14711 continue                                                             2028
      nlp=nlp+1                                                            2028
      dlx=0.0                                                              2029
14720 do 14721 k=1,ni                                                      2029
      if(ixx(k).eq.0)goto 14721                                            2030
      jb=ix(k)                                                             2030
      je=ix(k+1)-1                                                         2030
      jn=ix(k+1)-ix(k)                                                     2030
      bk=b(k)                                                              2031
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2032
      gk=dot_product(sc(1:jn),x(jb:je))                                    2033
      gk=(gk-svr*xb(k))/xs(k)                                              2034
      u=gk+xv(k)*b(k)                                                      2034
      au=abs(u)-vp(k)*al1                                                  2035
      if(au .gt. 0.0)goto 14741                                            2035
      b(k)=0.0                                                             2035
      goto 14751                                                           2036
14741 continue                                                             2037
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          2038
14751 continue                                                             2039
14731 continue                                                             2039
      d=b(k)-bk                                                            2039
      if(abs(d).le.0.0)goto 14721                                          2039
      dlx=max(dlx,xv(k)*d**2)                                              2040
      if(mm(k) .ne. 0)goto 14771                                           2040
      nin=nin+1                                                            2040
      if(nin.gt.nx)goto 14722                                              2041
      mm(k)=nin                                                            2041
      m(nin)=k                                                             2041
      sc(1:jn)=v(jx(jb:je))                                                2042
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 2043
14771 continue                                                             2044
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2045
      o=o+d*(xb(k)/xs(k))                                                  2046
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2047
14721 continue                                                             2048
14722 continue                                                             2048
      if(nin.gt.nx)goto 14712                                              2049
      d=svr/xm(0)                                                          2050
      if(d .eq. 0.0)goto 14791                                             2050
      b(0)=b(0)+d                                                          2050
      dlx=max(dlx,xm(0)*d**2)                                              2050
      r=r-d*v                                                              2050
14791 continue                                                             2051
      svr=svr-d*xm(0)                                                      2051
      if(dlx.lt.shr)goto 14712                                             2052
      if(nlp .le. maxit)goto 14811                                         2052
      jerr=-ilm                                                            2052
      return                                                               2052
14811 continue                                                             2053
14820 continue                                                             2053
14821 continue                                                             2053
      nlp=nlp+1                                                            2053
      dlx=0.0                                                              2054
14830 do 14831 l=1,nin                                                     2054
      k=m(l)                                                               2054
      jb=ix(k)                                                             2054
      je=ix(k+1)-1                                                         2055
      jn=ix(k+1)-ix(k)                                                     2055
      bk=b(k)                                                              2056
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2057
      gk=dot_product(sc(1:jn),x(jb:je))                                    2058
      gk=(gk-svr*xb(k))/xs(k)                                              2059
      u=gk+xv(k)*b(k)                                                      2059
      au=abs(u)-vp(k)*al1                                                  2060
      if(au .gt. 0.0)goto 14851                                            2060
      b(k)=0.0                                                             2060
      goto 14861                                                           2061
14851 continue                                                             2062
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          2063
14861 continue                                                             2064
14841 continue                                                             2064
      d=b(k)-bk                                                            2064
      if(abs(d).le.0.0)goto 14831                                          2064
      dlx=max(dlx,xv(k)*d**2)                                              2065
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2066
      o=o+d*(xb(k)/xs(k))                                                  2067
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2068
14831 continue                                                             2069
14832 continue                                                             2069
      d=svr/xm(0)                                                          2070
      if(d .eq. 0.0)goto 14881                                             2070
      b(0)=b(0)+d                                                          2070
      dlx=max(dlx,xm(0)*d**2)                                              2070
      r=r-d*v                                                              2070
14881 continue                                                             2071
      svr=svr-d*xm(0)                                                      2071
      if(dlx.lt.shr)goto 14822                                             2072
      if(nlp .le. maxit)goto 14901                                         2072
      jerr=-ilm                                                            2072
      return                                                               2072
14901 continue                                                             2073
      goto 14821                                                           2074
14822 continue                                                             2074
      goto 14711                                                           2075
14712 continue                                                             2075
      if(nin.gt.nx)goto 14672                                              2076
      sc=b(0)                                                              2076
      b0=0.0                                                               2077
14910 do 14911 j=1,nin                                                     2077
      l=m(j)                                                               2077
      jb=ix(l)                                                             2077
      je=ix(l+1)-1                                                         2078
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      2079
      b0=b0-b(l)*xb(l)/xs(l)                                               2080
14911 continue                                                             2081
14912 continue                                                             2081
      sc=sc+b0                                                             2082
14920 do 14921 i=1,no                                                      2082
      fi=sc(i)+g(i)                                                        2083
      if(fi .ge. fmin)goto 14941                                           2083
      q(i)=0.0                                                             2083
      goto 14931                                                           2083
14941 if(fi .le. fmax)goto 14951                                           2083
      q(i)=1.0                                                             2083
      goto 14961                                                           2084
14951 continue                                                             2084
      q(i)=1.0/(1.0+exp(-fi))                                              2084
14961 continue                                                             2085
14931 continue                                                             2085
14921 continue                                                             2086
14922 continue                                                             2086
      v=w*q*(1.0-q)                                                        2086
      xm(0)=sum(v)                                                         2086
      if(xm(0).lt.vmin)goto 14672                                          2087
      r=w*(y-q)                                                            2087
      svr=sum(r)                                                           2087
      o=0.0                                                                2088
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 14981                         2088
      kx=0                                                                 2089
14990 do 14991 j=1,nin                                                     2089
      k=m(j)                                                               2090
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 14991                           2090
      kx=1                                                                 2090
      goto 14992                                                           2091
14991 continue                                                             2092
14992 continue                                                             2092
      if(kx .ne. 0)goto 15011                                              2093
15020 do 15021 j=1,ni                                                      2093
      if(ixx(j).eq.1)goto 15021                                            2093
      if(ju(j).eq.0)goto 15021                                             2094
      jb=ix(j)                                                             2094
      je=ix(j+1)-1                                                         2094
      jn=ix(j+1)-ix(j)                                                     2095
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2096
      gj=dot_product(sc(1:jn),x(jb:je))                                    2097
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      2098
      if(ga(j) .le. al1*vp(j))goto 15041                                   2098
      ixx(j)=1                                                             2098
      kx=1                                                                 2098
15041 continue                                                             2099
15021 continue                                                             2100
15022 continue                                                             2100
      if(kx.eq.1) go to 10740                                              2101
      goto 14672                                                           2102
15011 continue                                                             2103
14981 continue                                                             2104
      goto 14671                                                           2105
14672 continue                                                             2105
      if(nin .le. nx)goto 15061                                            2105
      jerr=-10000-ilm                                                      2105
      goto 14592                                                           2105
15061 continue                                                             2106
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                2106
      kin(ilm)=nin                                                         2107
      a0(ilm)=b(0)                                                         2107
      alm(ilm)=al                                                          2107
      lmu=ilm                                                              2108
      devi=dev2(no,w,y,q,pmin)                                             2109
      dev(ilm)=(dev1-devi)/dev0                                            2110
      if(ilm.lt.mnl)goto 14591                                             2110
      if(flmin.ge.1.0)goto 14591                                           2111
      me=0                                                                 2111
15070 do 15071 j=1,nin                                                     2111
      if(a(j,ilm).ne.0.0) me=me+1                                          2111
15071 continue                                                             2111
15072 continue                                                             2111
      if(me.gt.ne)goto 14592                                               2112
      if(dev(ilm).gt.devmax)goto 14592                                     2112
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 14592                             2113
      if(xm(0).lt.vmin)goto 14592                                          2114
14591 continue                                                             2115
14592 continue                                                             2115
      g=log(q/(1.0-q))                                                     2116
      deallocate(xm,b,bs,v,r,sc,xv,q,mm,ga,ixx)                            2117
      return                                                               2118
      end                                                                  2119
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,n   2121 
     *lam,flmin,  ulam,shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev
     *,alm,nlp,jerr)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni)    2122
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni)          2123
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2124
      real, dimension (:,:), allocatable :: q                                   
      real, dimension (:), allocatable :: sxp,sxpl                              
      real, dimension (:), allocatable :: sc,xm,v,r,ga                          
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is,iy                           
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2135
      exmn=-exmx                                                           2136
      allocate(xm(0:ni),stat=ierr)                                         2136
      jerr=jerr+ierr                                                       2137
      allocate(r(1:no),stat=ierr)                                          2137
      jerr=jerr+ierr                                                       2138
      allocate(v(1:no),stat=ierr)                                          2138
      jerr=jerr+ierr                                                       2139
      allocate(mm(1:ni),stat=ierr)                                         2139
      jerr=jerr+ierr                                                       2140
      allocate(ga(1:ni),stat=ierr)                                         2140
      jerr=jerr+ierr                                                       2141
      allocate(iy(1:ni),stat=ierr)                                         2141
      jerr=jerr+ierr                                                       2142
      allocate(is(1:max(nc,ni)),stat=ierr)                                 2142
      jerr=jerr+ierr                                                       2143
      allocate(sxp(1:no),stat=ierr)                                        2143
      jerr=jerr+ierr                                                       2144
      allocate(sxpl(1:no),stat=ierr)                                       2144
      jerr=jerr+ierr                                                       2145
      allocate(sc(1:no),stat=ierr)                                         2145
      jerr=jerr+ierr                                                       2146
      if(jerr.ne.0) return                                                 2147
      pmax=1.0-pmin                                                        2147
      emin=pmin/pmax                                                       2147
      emax=1.0/emin                                                        2148
      pfm=(1.0+pmin)*pmin                                                  2148
      pfx=(1.0-pmin)*pmax                                                  2148
      vmin=pfm*pmax                                                        2149
      bta=parm                                                             2149
      omb=1.0-bta                                                          2149
      dev1=0.0                                                             2149
      dev0=0.0                                                             2150
15080 do 15081 ic=1,nc                                                     2150
      q0=dot_product(w,y(:,ic))                                            2151
      if(q0 .gt. pmin)goto 15101                                           2151
      jerr =8000+ic                                                        2151
      return                                                               2151
15101 continue                                                             2152
      if(q0 .lt. 1.0-pmin)goto 15121                                       2152
      jerr =9000+ic                                                        2152
      return                                                               2152
15121 continue                                                             2153
      b(1:ni,ic)=0.0                                                       2153
      b(0,ic)=log(q0)                                                      2153
      dev1=dev1-q0*b(0,ic)                                                 2154
15081 continue                                                             2155
15082 continue                                                             2155
      iy=0                                                                 2155
      al=0.0                                                               2156
      if(nonzero(no*nc,g) .ne. 0)goto 15141                                2157
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         2157
      sxp=0.0                                                              2158
15150 do 15151 ic=1,nc                                                     2158
      q(:,ic)=exp(b(0,ic))                                                 2158
      sxp=sxp+q(:,ic)                                                      2158
15151 continue                                                             2159
15152 continue                                                             2159
      goto 15161                                                           2160
15141 continue                                                             2160
15170 do 15171 i=1,no                                                      2160
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2160
15171 continue                                                             2160
15172 continue                                                             2160
      sxp=0.0                                                              2161
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 2161
      if(jerr.ne.0) return                                                 2162
      dev1=0.0                                                             2163
15180 do 15181 ic=1,nc                                                     2163
      q(:,ic)=b(0,ic)+g(:,ic)                                              2164
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             2165
      q(:,ic)=exp(q(:,ic))                                                 2165
      sxp=sxp+q(:,ic)                                                      2166
15181 continue                                                             2167
15182 continue                                                             2167
      sxpl=w*log(sxp)                                                      2167
15190 do 15191 ic=1,nc                                                     2167
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  2167
15191 continue                                                             2168
15192 continue                                                             2168
15161 continue                                                             2169
15131 continue                                                             2169
15200 do 15201 ic=1,nc                                                     2169
15210 do 15211 i=1,no                                                      2169
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               2169
15211 continue                                                             2169
15212 continue                                                             2169
15201 continue                                                             2170
15202 continue                                                             2170
      dev0=dev0+dev1                                                       2171
      if(kopt .le. 0)goto 15231                                            2172
      if(isd .le. 0)goto 15251                                             2172
      xv=0.25                                                              2172
      goto 15261                                                           2173
15251 continue                                                             2174
15270 do 15271 j=1,ni                                                      2174
      if(ju(j).eq.0)goto 15271                                             2174
      jb=ix(j)                                                             2174
      je=ix(j+1)-1                                                         2175
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        2176
15271 continue                                                             2177
15272 continue                                                             2177
15261 continue                                                             2178
15241 continue                                                             2178
15231 continue                                                             2179
      if(flmin .ge. 1.0)goto 15291                                         2179
      eqs=max(eps,flmin)                                                   2179
      alf=eqs**(1.0/(nlam-1))                                              2179
15291 continue                                                             2180
      m=0                                                                  2180
      mm=0                                                                 2180
      nin=0                                                                2180
      nlp=0                                                                2180
      mnl=min(mnlam,nlam)                                                  2180
      bs=0.0                                                               2180
      svr=0.0                                                              2180
      o=0.0                                                                2181
      shr=shri*dev0                                                        2181
      ga=0.0                                                               2182
15300 do 15301 ic=1,nc                                                     2182
      v=q(:,ic)/sxp                                                        2182
      r=w*(y(:,ic)-v)                                                      2182
      v=w*v*(1.0-v)                                                        2183
15310 do 15311 j=1,ni                                                      2183
      if(ju(j).eq.0)goto 15311                                             2184
      jb=ix(j)                                                             2184
      je=ix(j+1)-1                                                         2184
      jn=ix(j+1)-ix(j)                                                     2185
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2186
      gj=dot_product(sc(1:jn),x(jb:je))                                    2187
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2188
15311 continue                                                             2189
15312 continue                                                             2189
15301 continue                                                             2190
15302 continue                                                             2190
15320 do 15321 ilm=1,nlam                                                  2190
      al0=al                                                               2191
      if(flmin .lt. 1.0)goto 15341                                         2191
      al=ulam(ilm)                                                         2191
      goto 15331                                                           2192
15341 if(ilm .le. 2)goto 15351                                             2192
      al=al*alf                                                            2192
      goto 15331                                                           2193
15351 if(ilm .ne. 1)goto 15361                                             2193
      al=big                                                               2193
      goto 15371                                                           2194
15361 continue                                                             2194
      al0=0.0                                                              2195
15380 do 15381 j=1,ni                                                      2195
      if(ju(j).eq.0)goto 15381                                             2195
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2195
15381 continue                                                             2196
15382 continue                                                             2196
      al0=al0/max(bta,1.0e-3)                                              2196
      al=alf*al0                                                           2197
15371 continue                                                             2198
15331 continue                                                             2198
      al2=al*omb                                                           2198
      al1=al*bta                                                           2198
      tlam=bta*(2.0*al-al0)                                                2199
15390 do 15391 k=1,ni                                                      2199
      if(iy(k).eq.1)goto 15391                                             2199
      if(ju(k).eq.0)goto 15391                                             2200
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      2201
15391 continue                                                             2202
15392 continue                                                             2202
10740 continue                                                             2203
15400 continue                                                             2203
15401 continue                                                             2203
      ixx=0                                                                2203
      jxx=ixx                                                              2203
      ig=0                                                                 2204
15410 do 15411 ic=1,nc                                                     2204
      bs(0,ic)=b(0,ic)                                                     2205
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          2206
      xm(0)=0.0                                                            2206
      svr=0.0                                                              2206
      o=0.0                                                                2207
15420 do 15421 i=1,no                                                      2207
      pic=q(i,ic)/sxp(i)                                                   2208
      if(pic .ge. pfm)goto 15441                                           2208
      pic=0.0                                                              2208
      v(i)=0.0                                                             2208
      goto 15431                                                           2209
15441 if(pic .le. pfx)goto 15451                                           2209
      pic=1.0                                                              2209
      v(i)=0.0                                                             2209
      goto 15461                                                           2210
15451 continue                                                             2210
      v(i)=w(i)*pic*(1.0-pic)                                              2210
      xm(0)=xm(0)+v(i)                                                     2210
15461 continue                                                             2211
15431 continue                                                             2211
      r(i)=w(i)*(y(i,ic)-pic)                                              2211
      svr=svr+r(i)                                                         2212
15421 continue                                                             2213
15422 continue                                                             2213
      if(xm(0).le.vmin)goto 15411                                          2213
      ig=1                                                                 2214
15470 do 15471 j=1,ni                                                      2214
      if(iy(j).eq.0)goto 15471                                             2215
      jb=ix(j)                                                             2215
      je=ix(j+1)-1                                                         2216
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             2217
      if(kopt .ne. 0)goto 15491                                            2218
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       2219
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          2220
15491 continue                                                             2221
15471 continue                                                             2222
15472 continue                                                             2222
15500 continue                                                             2222
15501 continue                                                             2222
      nlp=nlp+1                                                            2222
      dlx=0.0                                                              2223
15510 do 15511 k=1,ni                                                      2223
      if(iy(k).eq.0)goto 15511                                             2224
      jb=ix(k)                                                             2224
      je=ix(k+1)-1                                                         2224
      jn=ix(k+1)-ix(k)                                                     2224
      bk=b(k,ic)                                                           2225
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2226
      gk=dot_product(sc(1:jn),x(jb:je))                                    2227
      gk=(gk-svr*xb(k))/xs(k)                                              2228
      u=gk+xv(k,ic)*b(k,ic)                                                2228
      au=abs(u)-vp(k)*al1                                                  2229
      if(au .gt. 0.0)goto 15531                                            2229
      b(k,ic)=0.0                                                          2229
      goto 15541                                                           2230
15531 continue                                                             2231
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   2233 
     *)
15541 continue                                                             2234
15521 continue                                                             2234
      d=b(k,ic)-bk                                                         2234
      if(abs(d).le.0.0)goto 15511                                          2235
      dlx=max(dlx,xv(k,ic)*d**2)                                           2236
      if(mm(k) .ne. 0)goto 15561                                           2236
      nin=nin+1                                                            2237
      if(nin .le. nx)goto 15581                                            2237
      jxx=1                                                                2237
      goto 15512                                                           2237
15581 continue                                                             2238
      mm(k)=nin                                                            2238
      m(nin)=k                                                             2239
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             2240
15561 continue                                                             2241
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2242
      o=o+d*(xb(k)/xs(k))                                                  2243
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2244
15511 continue                                                             2245
15512 continue                                                             2245
      if(jxx.gt.0)goto 15502                                               2246
      d=svr/xm(0)                                                          2247
      if(d .eq. 0.0)goto 15601                                             2247
      b(0,ic)=b(0,ic)+d                                                    2247
      dlx=max(dlx,xm(0)*d**2)                                              2248
      r=r-d*v                                                              2248
      svr=svr-d*xm(0)                                                      2249
15601 continue                                                             2250
      if(dlx.lt.shr)goto 15502                                             2250
      if(nlp .le. maxit)goto 15621                                         2250
      jerr=-ilm                                                            2250
      return                                                               2250
15621 continue                                                             2251
15630 continue                                                             2251
15631 continue                                                             2251
      nlp=nlp+1                                                            2251
      dlx=0.0                                                              2252
15640 do 15641 l=1,nin                                                     2252
      k=m(l)                                                               2252
      jb=ix(k)                                                             2252
      je=ix(k+1)-1                                                         2253
      jn=ix(k+1)-ix(k)                                                     2253
      bk=b(k,ic)                                                           2254
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2255
      gk=dot_product(sc(1:jn),x(jb:je))                                    2256
      gk=(gk-svr*xb(k))/xs(k)                                              2257
      u=gk+xv(k,ic)*b(k,ic)                                                2257
      au=abs(u)-vp(k)*al1                                                  2258
      if(au .gt. 0.0)goto 15661                                            2258
      b(k,ic)=0.0                                                          2258
      goto 15671                                                           2259
15661 continue                                                             2260
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   2262 
     *)
15671 continue                                                             2263
15651 continue                                                             2263
      d=b(k,ic)-bk                                                         2263
      if(abs(d).le.0.0)goto 15641                                          2264
      dlx=max(dlx,xv(k,ic)*d**2)                                           2265
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2266
      o=o+d*(xb(k)/xs(k))                                                  2267
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2268
15641 continue                                                             2269
15642 continue                                                             2269
      d=svr/xm(0)                                                          2270
      if(d .eq. 0.0)goto 15691                                             2270
      b(0,ic)=b(0,ic)+d                                                    2270
      dlx=max(dlx,xm(0)*d**2)                                              2271
      r=r-d*v                                                              2271
      svr=svr-d*xm(0)                                                      2272
15691 continue                                                             2273
      if(dlx.lt.shr)goto 15632                                             2273
      if(nlp .le. maxit)goto 15711                                         2273
      jerr=-ilm                                                            2273
      return                                                               2273
15711 continue                                                             2274
      goto 15631                                                           2275
15632 continue                                                             2275
      goto 15501                                                           2276
15502 continue                                                             2276
      if(jxx.gt.0)goto 15412                                               2277
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                         2278
      if(ixx .ne. 0)goto 15731                                             2279
15740 do 15741 j=1,nin                                                     2279
      k=m(j)                                                               2280
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 15761                2280
      ixx=1                                                                2280
      goto 15742                                                           2280
15761 continue                                                             2281
15741 continue                                                             2282
15742 continue                                                             2282
15731 continue                                                             2283
      sc=b(0,ic)+g(:,ic)                                                   2283
      b0=0.0                                                               2284
15770 do 15771 j=1,nin                                                     2284
      l=m(j)                                                               2284
      jb=ix(l)                                                             2284
      je=ix(l+1)-1                                                         2285
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   2286
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            2287
15771 continue                                                             2288
15772 continue                                                             2288
      sc=min(max(exmn,sc+b0),exmx)                                         2289
      sxp=sxp-q(:,ic)                                                      2290
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          2291
      sxp=sxp+q(:,ic)                                                      2292
15411 continue                                                             2293
15412 continue                                                             2293
      s=-sum(b(0,:))/nc                                                    2293
      b(0,:)=b(0,:)+s                                                      2293
      sc=s                                                                 2293
      b0=0.0                                                               2294
15780 do 15781 j=1,nin                                                     2294
      l=m(j)                                                               2295
      if(vp(l) .gt. 0.0)goto 15801                                         2295
      s=sum(b(l,:))/nc                                                     2295
      goto 15811                                                           2296
15801 continue                                                             2296
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                     2296
15811 continue                                                             2297
15791 continue                                                             2297
      b(l,:)=b(l,:)-s                                                      2298
      jb=ix(l)                                                             2298
      je=ix(l+1)-1                                                         2299
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         2300
      b0=b0+s*xb(l)/xs(l)                                                  2301
15781 continue                                                             2302
15782 continue                                                             2302
      sc=sc+b0                                                             2302
      sc=exp(sc)                                                           2302
      sxp=sxp*sc                                                           2302
15820 do 15821 ic=1,nc                                                     2302
      q(:,ic)=q(:,ic)*sc                                                   2302
15821 continue                                                             2303
15822 continue                                                             2303
      if(jxx.gt.0)goto 15402                                               2303
      if(ig.eq.0)goto 15402                                                2304
      if(ixx .ne. 0)goto 15841                                             2305
15850 do 15851 j=1,ni                                                      2305
      if(iy(j).eq.1)goto 15851                                             2305
      if(ju(j).eq.0)goto 15851                                             2305
      ga(j)=0.0                                                            2305
15851 continue                                                             2306
15852 continue                                                             2306
15860 do 15861 ic=1,nc                                                     2306
      v=q(:,ic)/sxp                                                        2306
      r=w*(y(:,ic)-v)                                                      2306
      v=w*v*(1.0-v)                                                        2307
15870 do 15871 j=1,ni                                                      2307
      if(iy(j).eq.1)goto 15871                                             2307
      if(ju(j).eq.0)goto 15871                                             2308
      jb=ix(j)                                                             2308
      je=ix(j+1)-1                                                         2308
      jn=ix(j+1)-ix(j)                                                     2309
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2310
      gj=dot_product(sc(1:jn),x(jb:je))                                    2311
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2312
15871 continue                                                             2313
15872 continue                                                             2313
15861 continue                                                             2314
15862 continue                                                             2314
15880 do 15881 k=1,ni                                                      2314
      if(iy(k).eq.1)goto 15881                                             2314
      if(ju(k).eq.0)goto 15881                                             2315
      if(ga(k) .le. al1*vp(k))goto 15901                                   2315
      iy(k)=1                                                              2315
      ixx=1                                                                2315
15901 continue                                                             2316
15881 continue                                                             2317
15882 continue                                                             2317
      if(ixx.eq.1) go to 10740                                             2318
      goto 15402                                                           2319
15841 continue                                                             2320
      goto 15401                                                           2321
15402 continue                                                             2321
      if(jxx .le. 0)goto 15921                                             2321
      jerr=-10000-ilm                                                      2321
      goto 15322                                                           2321
15921 continue                                                             2321
      devi=0.0                                                             2322
15930 do 15931 ic=1,nc                                                     2323
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          2323
      a0(ic,ilm)=b(0,ic)                                                   2324
15940 do 15941 i=1,no                                                      2324
      if(y(i,ic).le.0.0)goto 15941                                         2325
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           2326
15941 continue                                                             2327
15942 continue                                                             2327
15931 continue                                                             2328
15932 continue                                                             2328
      kin(ilm)=nin                                                         2328
      alm(ilm)=al                                                          2328
      lmu=ilm                                                              2329
      dev(ilm)=(dev1-devi)/dev0                                            2329
      if(ig.eq.0)goto 15322                                                2330
      if(ilm.lt.mnl)goto 15321                                             2330
      if(flmin.ge.1.0)goto 15321                                           2331
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 15322             2332
      if(dev(ilm).gt.devmax)goto 15322                                     2332
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 15322                             2333
15321 continue                                                             2334
15322 continue                                                             2334
      g=log(q)                                                             2334
15950 do 15951 i=1,no                                                      2334
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2334
15951 continue                                                             2335
15952 continue                                                             2335
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy)                      2336
      return                                                               2337
      end                                                                  2338
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  2339
      real a0(nc),ca(nx,nc),x(*),f(nc,n)                                   2339
      integer ia(*),ix(*),jx(*)                                            2340
15960 do 15961 ic=1,nc                                                     2340
      f(ic,:)=a0(ic)                                                       2340
15961 continue                                                             2341
15962 continue                                                             2341
15970 do 15971 j=1,nin                                                     2341
      k=ia(j)                                                              2341
      kb=ix(k)                                                             2341
      ke=ix(k+1)-1                                                         2342
15980 do 15981 ic=1,nc                                                     2342
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    2342
15981 continue                                                             2343
15982 continue                                                             2343
15971 continue                                                             2344
15972 continue                                                             2344
      return                                                               2345
      end                                                                  2346
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,cl,ne,nx,nlam,flmin,   2348 
     *ulam,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam)              2349
      real ca(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)                        2350
      integer jd(*),ia(nx),nin(nlam)                                       2351
      real, dimension (:), allocatable :: xs,ww,vq                              
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16001                                    2355
      jerr=10000                                                           2355
      return                                                               2355
16001 continue                                                             2356
      allocate(ww(1:no),stat=jerr)                                         2357
      allocate(ju(1:ni),stat=ierr)                                         2357
      jerr=jerr+ierr                                                       2358
      allocate(vq(1:ni),stat=ierr)                                         2358
      jerr=jerr+ierr                                                       2359
      if(isd .le. 0)goto 16021                                             2359
      allocate(xs(1:ni),stat=ierr)                                         2359
      jerr=jerr+ierr                                                       2359
16021 continue                                                             2360
      if(jerr.ne.0) return                                                 2361
      call chkvars(no,ni,x,ju)                                             2362
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2363
      if(maxval(ju) .gt. 0)goto 16041                                      2363
      jerr=7777                                                            2363
      return                                                               2363
16041 continue                                                             2364
      vq=max(0.0,vp)                                                       2364
      vq=vq*ni/sum(vq)                                                     2365
      ww=max(0.0,w)                                                        2365
      sw=sum(ww)                                                           2366
      if(sw .gt. 0.0)goto 16061                                            2366
      jerr=9999                                                            2366
      return                                                               2366
16061 continue                                                             2366
      ww=ww/sw                                                             2367
      call cstandard(no,ni,x,ww,ju,isd,xs)                                 2368
      if(isd .le. 0)goto 16081                                             2368
16090 do 16091 j=1,ni                                                      2368
      cl(:,j)=cl(:,j)*xs(j)                                                2368
16091 continue                                                             2368
16092 continue                                                             2368
16081 continue                                                             2369
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,   2371 
     *thr,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 2371
      dev0=2.0*sw*dev0                                                     2372
      if(isd .le. 0)goto 16111                                             2372
16120 do 16121 k=1,lmu                                                     2372
      nk=nin(k)                                                            2372
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   2372
16121 continue                                                             2372
16122 continue                                                             2372
16111 continue                                                             2373
      deallocate(ww,ju,vq)                                                 2373
      if(isd.gt.0) deallocate(xs)                                          2374
      return                                                               2375
      end                                                                  2376
      subroutine cstandard (no,ni,x,w,ju,isd,xs)                           2377
      real x(no,ni),w(no),xs(ni)                                           2377
      integer ju(ni)                                                       2378
16130 do 16131 j=1,ni                                                      2378
      if(ju(j).eq.0)goto 16131                                             2379
      xm=dot_product(w,x(:,j))                                             2379
      x(:,j)=x(:,j)-xm                                                     2380
      if(isd .le. 0)goto 16151                                             2380
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 2380
      x(:,j)=x(:,j)/xs(j)                                                  2380
16151 continue                                                             2381
16131 continue                                                             2382
16132 continue                                                             2382
      return                                                               2383
      end                                                                  2384
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,cl,ne,nx,nlam,flmin,   2386 
     *ulam,cthri,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam)              2387
      real ao(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)                        2388
      integer ju(ni),m(nx),kin(nlam)                                       2389
      real, dimension (:), allocatable :: w,dk,v,xs,wr,a,as,f,dq                
      real, dimension (:), allocatable :: e,uu,ga                               
      integer, dimension (:), allocatable :: jp,kp,mm,ixx                       
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2395
      sml=sml*100.0                                                        2395
      devmax=devmax*0.99/0.999                                             2396
      allocate(e(1:no),stat=jerr)                                          2397
      allocate(uu(1:no),stat=ierr)                                         2397
      jerr=jerr+ierr                                                       2398
      allocate(f(1:no),stat=ierr)                                          2398
      jerr=jerr+ierr                                                       2399
      allocate(w(1:no),stat=ierr)                                          2399
      jerr=jerr+ierr                                                       2400
      allocate(v(1:ni),stat=ierr)                                          2400
      jerr=jerr+ierr                                                       2401
      allocate(a(1:ni),stat=ierr)                                          2401
      jerr=jerr+ierr                                                       2402
      allocate(as(1:ni),stat=ierr)                                         2402
      jerr=jerr+ierr                                                       2403
      allocate(xs(1:ni),stat=ierr)                                         2403
      jerr=jerr+ierr                                                       2404
      allocate(ga(1:ni),stat=ierr)                                         2404
      jerr=jerr+ierr                                                       2405
      allocate(ixx(1:ni),stat=ierr)                                        2405
      jerr=jerr+ierr                                                       2406
      allocate(jp(1:no),stat=ierr)                                         2406
      jerr=jerr+ierr                                                       2407
      allocate(kp(1:no),stat=ierr)                                         2407
      jerr=jerr+ierr                                                       2408
      allocate(dk(1:no),stat=ierr)                                         2408
      jerr=jerr+ierr                                                       2409
      allocate(wr(1:no),stat=ierr)                                         2409
      jerr=jerr+ierr                                                       2410
      allocate(dq(1:no),stat=ierr)                                         2410
      jerr=jerr+ierr                                                       2411
      allocate(mm(1:ni),stat=ierr)                                         2411
      jerr=jerr+ierr                                                       2412
      if(jerr.ne.0)go to 11910                                             2413
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2414
      if(jerr.ne.0) go to 11910                                            2414
      alpha=parm                                                           2415
      oma=1.0-alpha                                                        2415
      nlm=0                                                                2415
      ixx=0                                                                2415
      al=0.0                                                               2416
      dq=d*q                                                               2416
      call died(no,nk,dq,kp,jp,dk)                                         2417
      a=0.0                                                                2417
      f(1)=0.0                                                             2417
      fmax=log(huge(f(1))*0.1)                                             2418
      if(nonzero(no,g) .eq. 0)goto 16171                                   2418
      f=g-dot_product(q,g)                                                 2419
      e=q*exp(sign(min(abs(f),fmax),f))                                    2420
      goto 16181                                                           2421
16171 continue                                                             2421
      f=0.0                                                                2421
      e=q                                                                  2421
16181 continue                                                             2422
16161 continue                                                             2422
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 2423
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         2423
      dev0=rr                                                              2424
16190 do 16191 i=1,no                                                      2424
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 16211                   2424
      w(i)=0.0                                                             2424
      wr(i)=w(i)                                                           2424
16211 continue                                                             2424
16191 continue                                                             2425
16192 continue                                                             2425
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2426
      if(jerr.ne.0) go to 11910                                            2427
      if(flmin .ge. 1.0)goto 16231                                         2427
      eqs=max(eps,flmin)                                                   2427
      alf=eqs**(1.0/(nlam-1))                                              2427
16231 continue                                                             2428
      m=0                                                                  2428
      mm=0                                                                 2428
      nlp=0                                                                2428
      nin=nlp                                                              2428
      mnl=min(mnlam,nlam)                                                  2428
      as=0.0                                                               2428
      cthr=cthri*dev0                                                      2429
16240 do 16241 j=1,ni                                                      2429
      if(ju(j).eq.0)goto 16241                                             2429
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2429
16241 continue                                                             2430
16242 continue                                                             2430
16250 do 16251 ilm=1,nlam                                                  2430
      al0=al                                                               2431
      if(flmin .lt. 1.0)goto 16271                                         2431
      al=ulam(ilm)                                                         2431
      goto 16261                                                           2432
16271 if(ilm .le. 2)goto 16281                                             2432
      al=al*alf                                                            2432
      goto 16261                                                           2433
16281 if(ilm .ne. 1)goto 16291                                             2433
      al=big                                                               2433
      goto 16301                                                           2434
16291 continue                                                             2434
      al0=0.0                                                              2435
16310 do 16311 j=1,ni                                                      2435
      if(ju(j).eq.0)goto 16311                                             2435
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2435
16311 continue                                                             2436
16312 continue                                                             2436
      al0=al0/max(parm,1.0e-3)                                             2436
      al=alf*al0                                                           2437
16301 continue                                                             2438
16261 continue                                                             2438
      sa=alpha*al                                                          2438
      omal=oma*al                                                          2438
      tlam=alpha*(2.0*al-al0)                                              2439
16320 do 16321 k=1,ni                                                      2439
      if(ixx(k).eq.1)goto 16321                                            2439
      if(ju(k).eq.0)goto 16321                                             2440
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2441
16321 continue                                                             2442
16322 continue                                                             2442
10740 continue                                                             2443
16330 continue                                                             2443
16331 continue                                                             2443
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2444
      call vars(no,ni,x,w,ixx,v)                                           2445
16340 continue                                                             2445
16341 continue                                                             2445
      nlp=nlp+1                                                            2445
      dli=0.0                                                              2446
16350 do 16351 j=1,ni                                                      2446
      if(ixx(j).eq.0)goto 16351                                            2447
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2448
      if(abs(u) .gt. vp(j)*sa)goto 16371                                   2448
      at=0.0                                                               2448
      goto 16381                                                           2449
16371 continue                                                             2449
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o   2451 
     *mal)))
16381 continue                                                             2452
16361 continue                                                             2452
      if(at .eq. a(j))goto 16401                                           2452
      del=at-a(j)                                                          2452
      a(j)=at                                                              2452
      dli=max(dli,v(j)*del**2)                                             2453
      wr=wr-del*w*x(:,j)                                                   2453
      f=f+del*x(:,j)                                                       2454
      if(mm(j) .ne. 0)goto 16421                                           2454
      nin=nin+1                                                            2454
      if(nin.gt.nx)goto 16352                                              2455
      mm(j)=nin                                                            2455
      m(nin)=j                                                             2456
16421 continue                                                             2457
16401 continue                                                             2458
16351 continue                                                             2459
16352 continue                                                             2459
      if(nin.gt.nx)goto 16342                                              2459
      if(dli.lt.cthr)goto 16342                                            2460
      if(nlp .le. maxit)goto 16441                                         2460
      jerr=-ilm                                                            2460
      return                                                               2460
16441 continue                                                             2461
16450 continue                                                             2461
16451 continue                                                             2461
      nlp=nlp+1                                                            2461
      dli=0.0                                                              2462
16460 do 16461 l=1,nin                                                     2462
      j=m(l)                                                               2463
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2464
      if(abs(u) .gt. vp(j)*sa)goto 16481                                   2464
      at=0.0                                                               2464
      goto 16491                                                           2465
16481 continue                                                             2465
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o   2467 
     *mal)))
16491 continue                                                             2468
16471 continue                                                             2468
      if(at .eq. a(j))goto 16511                                           2468
      del=at-a(j)                                                          2468
      a(j)=at                                                              2468
      dli=max(dli,v(j)*del**2)                                             2469
      wr=wr-del*w*x(:,j)                                                   2469
      f=f+del*x(:,j)                                                       2470
16511 continue                                                             2471
16461 continue                                                             2472
16462 continue                                                             2472
      if(dli.lt.cthr)goto 16452                                            2472
      if(nlp .le. maxit)goto 16531                                         2472
      jerr=-ilm                                                            2472
      return                                                               2472
16531 continue                                                             2473
      goto 16451                                                           2474
16452 continue                                                             2474
      goto 16341                                                           2475
16342 continue                                                             2475
      if(nin.gt.nx)goto 16332                                              2476
      e=q*exp(sign(min(abs(f),fmax),f))                                    2477
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2478
      if(jerr .eq. 0)goto 16551                                            2478
      jerr=jerr-ilm                                                        2478
      go to 11910                                                          2478
16551 continue                                                             2479
      ix=0                                                                 2480
16560 do 16561 j=1,nin                                                     2480
      k=m(j)                                                               2481
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 16561                           2481
      ix=1                                                                 2481
      goto 16562                                                           2481
16561 continue                                                             2482
16562 continue                                                             2482
      if(ix .ne. 0)goto 16581                                              2483
16590 do 16591 k=1,ni                                                      2483
      if(ixx(k).eq.1)goto 16591                                            2483
      if(ju(k).eq.0)goto 16591                                             2484
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2485
      if(ga(k) .le. sa*vp(k))goto 16611                                    2485
      ixx(k)=1                                                             2485
      ix=1                                                                 2485
16611 continue                                                             2486
16591 continue                                                             2487
16592 continue                                                             2487
      if(ix.eq.1) go to 10740                                              2488
      goto 16332                                                           2489
16581 continue                                                             2490
      goto 16331                                                           2491
16332 continue                                                             2491
      if(nin .le. nx)goto 16631                                            2491
      jerr=-10000-ilm                                                      2491
      goto 16252                                                           2491
16631 continue                                                             2492
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               2492
      kin(ilm)=nin                                                         2493
      alm(ilm)=al                                                          2493
      lmu=ilm                                                              2494
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   2495
      if(ilm.lt.mnl)goto 16251                                             2495
      if(flmin.ge.1.0)goto 16251                                           2496
      me=0                                                                 2496
16640 do 16641 j=1,nin                                                     2496
      if(ao(j,ilm).ne.0.0) me=me+1                                         2496
16641 continue                                                             2496
16642 continue                                                             2496
      if(me.gt.ne)goto 16252                                               2497
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 16252              2498
      if(dev(ilm).gt.devmax)goto 16252                                     2499
16251 continue                                                             2500
16252 continue                                                             2500
      g=f                                                                  2501
11910 continue                                                             2501
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)              2502
      return                                                               2503
      end                                                                  2504
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 2505
      real ca(nin),x(n,*),f(n)                                             2505
      integer ia(nin)                                                      2506
      f=0.0                                                                2506
      if(nin.le.0) return                                                  2507
16650 do 16651 i=1,n                                                       2507
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      2507
16651 continue                                                             2508
16652 continue                                                             2508
      return                                                               2509
      end                                                                  2510
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         2511
      real y(no),d(no),q(no)                                               2511
      integer jp(no),kp(*)                                                 2512
16660 do 16661 j=1,no                                                      2512
      jp(j)=j                                                              2512
16661 continue                                                             2512
16662 continue                                                             2512
      call psort7(y,jp,1,no)                                               2513
      nj=0                                                                 2513
16670 do 16671 j=1,no                                                      2513
      if(q(jp(j)).le.0.0)goto 16671                                        2513
      nj=nj+1                                                              2513
      jp(nj)=jp(j)                                                         2513
16671 continue                                                             2514
16672 continue                                                             2514
      if(nj .ne. 0)goto 16691                                              2514
      jerr=20000                                                           2514
      return                                                               2514
16691 continue                                                             2515
      j=1                                                                  2515
16700 continue                                                             2515
16701 if(d(jp(j)).gt.0.0)goto 16702                                        2515
      j=j+1                                                                2515
      if(j.gt.nj)goto 16702                                                2515
      goto 16701                                                           2516
16702 continue                                                             2516
      if(j .lt. nj-1)goto 16721                                            2516
      jerr=30000                                                           2516
      return                                                               2516
16721 continue                                                             2517
      t0=y(jp(j))                                                          2517
      j0=j-1                                                               2518
      if(j0 .le. 0)goto 16741                                              2519
16750 continue                                                             2519
16751 if(y(jp(j0)).lt.t0)goto 16752                                        2519
      j0=j0-1                                                              2519
      if(j0.eq.0)goto 16752                                                2519
      goto 16751                                                           2520
16752 continue                                                             2520
      if(j0 .le. 0)goto 16771                                              2520
      nj=nj-j0                                                             2520
16780 do 16781 j=1,nj                                                      2520
      jp(j)=jp(j+j0)                                                       2520
16781 continue                                                             2520
16782 continue                                                             2520
16771 continue                                                             2521
16741 continue                                                             2522
      jerr=0                                                               2522
      nk=0                                                                 2522
      yk=t0                                                                2522
      j=2                                                                  2523
16790 continue                                                             2523
16791 continue                                                             2523
16800 continue                                                             2524
16801 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 16802                     2524
      j=j+1                                                                2524
      if(j.gt.nj)goto 16802                                                2524
      goto 16801                                                           2525
16802 continue                                                             2525
      nk=nk+1                                                              2525
      kp(nk)=j-1                                                           2525
      if(j.gt.nj)goto 16792                                                2526
      if(j .ne. nj)goto 16821                                              2526
      nk=nk+1                                                              2526
      kp(nk)=nj                                                            2526
      goto 16792                                                           2526
16821 continue                                                             2527
      yk=y(jp(j))                                                          2527
      j=j+1                                                                2528
      goto 16791                                                           2529
16792 continue                                                             2529
      return                                                               2530
      end                                                                  2531
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     2532
      real d(no),dk(nk),wr(no),w(no)                                       2533
      real e(no),u(no),b,c                                                 2533
      integer kp(nk),jp(no)                                                2534
      call usk(no,nk,kp,jp,e,u)                                            2535
      b=dk(1)/u(1)                                                         2535
      c=dk(1)/u(1)**2                                                      2535
      jerr=0                                                               2536
16830 do 16831 j=1,kp(1)                                                   2536
      i=jp(j)                                                              2537
      w(i)=e(i)*(b-e(i)*c)                                                 2537
      if(w(i) .gt. 0.0)goto 16851                                          2537
      jerr=-30000                                                          2537
      return                                                               2537
16851 continue                                                             2538
      wr(i)=d(i)-e(i)*b                                                    2539
16831 continue                                                             2540
16832 continue                                                             2540
16860 do 16861 k=2,nk                                                      2540
      j1=kp(k-1)+1                                                         2540
      j2=kp(k)                                                             2541
      b=b+dk(k)/u(k)                                                       2541
      c=c+dk(k)/u(k)**2                                                    2542
16870 do 16871 j=j1,j2                                                     2542
      i=jp(j)                                                              2543
      w(i)=e(i)*(b-e(i)*c)                                                 2543
      if(w(i) .gt. 0.0)goto 16891                                          2543
      jerr=-30000                                                          2543
      return                                                               2543
16891 continue                                                             2544
      wr(i)=d(i)-e(i)*b                                                    2545
16871 continue                                                             2546
16872 continue                                                             2546
16861 continue                                                             2547
16862 continue                                                             2547
      return                                                               2548
      end                                                                  2549
      subroutine vars(no,ni,x,w,ixx,v)                                     2550
      real x(no,ni),w(no),v(ni)                                            2550
      integer ixx(ni)                                                      2551
16900 do 16901 j=1,ni                                                      2551
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                        2551
16901 continue                                                             2552
16902 continue                                                             2552
      return                                                               2553
      end                                                                  2554
      subroutine died(no,nk,d,kp,jp,dk)                                    2555
      real d(no),dk(nk)                                                    2555
      integer kp(nk),jp(no)                                                2556
      dk(1)=sum(d(jp(1:kp(1))))                                            2557
16910 do 16911 k=2,nk                                                      2557
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  2557
16911 continue                                                             2558
16912 continue                                                             2558
      return                                                               2559
      end                                                                  2560
      subroutine usk(no,nk,kp,jp,e,u)                                      2561
      real e(no),u(nk),h                                                   2561
      integer kp(nk),jp(no)                                                2562
      h=0.0                                                                2563
16920 do 16921 k=nk,1,-1                                                   2563
      j2=kp(k)                                                             2564
      j1=1                                                                 2564
      if(k.gt.1) j1=kp(k-1)+1                                              2565
16930 do 16931 j=j2,j1,-1                                                  2565
      h=h+e(jp(j))                                                         2565
16931 continue                                                             2566
16932 continue                                                             2566
      u(k)=h                                                               2567
16921 continue                                                             2568
16922 continue                                                             2568
      return                                                               2569
      end                                                                  2570
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             2571
      real d(no),dk(nk),f(no)                                              2572
      integer kp(nk),jp(no)                                                2572
      real e(no),u(nk),s                                                   2573
      call usk(no,nk,kp,jp,e,u)                                            2573
      u=log(u)                                                             2574
      risk=dot_product(d,f)-dot_product(dk,u)                              2575
      return                                                               2576
      end                                                                  2577
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 2578
      real x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(nlam)          2579
      real, dimension (:), allocatable :: dk,f,xm,dq,q                          
      real, dimension (:), allocatable :: e,uu                                  
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          2585
      allocate(q(1:no),stat=ierr)                                          2585
      jerr=jerr+ierr                                                       2586
      allocate(uu(1:no),stat=ierr)                                         2586
      jerr=jerr+ierr                                                       2587
      allocate(f(1:no),stat=ierr)                                          2587
      jerr=jerr+ierr                                                       2588
      allocate(dk(1:no),stat=ierr)                                         2588
      jerr=jerr+ierr                                                       2589
      allocate(jp(1:no),stat=ierr)                                         2589
      jerr=jerr+ierr                                                       2590
      allocate(kp(1:no),stat=ierr)                                         2590
      jerr=jerr+ierr                                                       2591
      allocate(dq(1:no),stat=ierr)                                         2591
      jerr=jerr+ierr                                                       2592
      allocate(xm(1:ni),stat=ierr)                                         2592
      jerr=jerr+ierr                                                       2593
      if(jerr.ne.0) go to 11910                                            2594
      q=max(0.0,w)                                                         2594
      sw=sum(q)                                                            2595
      if(sw .gt. 0.0)goto 16951                                            2595
      jerr=9999                                                            2595
      go to 11910                                                          2595
16951 continue                                                             2596
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2597
      if(jerr.ne.0) go to 11910                                            2597
      fmax=log(huge(e(1))*0.1)                                             2598
      dq=d*q                                                               2598
      call died(no,nk,dq,kp,jp,dk)                                         2598
      gm=dot_product(q,g)/sw                                               2599
16960 do 16961 j=1,ni                                                      2599
      xm(j)=dot_product(q,x(:,j))/sw                                       2599
16961 continue                                                             2600
16962 continue                                                             2600
16970 do 16971 lam=1,nlam                                                  2601
16980 do 16981 i=1,no                                                      2601
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       2602
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        2603
16981 continue                                                             2604
16982 continue                                                             2604
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          2605
16971 continue                                                             2606
16972 continue                                                             2606
11910 continue                                                             2606
      deallocate(e,uu,dk,f,jp,kp,dq)                                       2607
      return                                                               2608
      end                                                                  2609
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,cl,ne,nx,nlam,flmin,u   2611 
     *lam,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    2612
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)               2613
      integer jd(*),ia(nx),nin(nlam)                                       2614
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 17001                                    2618
      jerr=10000                                                           2618
      return                                                               2618
17001 continue                                                             2619
      if(minval(y) .ge. 0.0)goto 17021                                     2619
      jerr=8888                                                            2619
      return                                                               2619
17021 continue                                                             2620
      allocate(ww(1:no),stat=jerr)                                         2621
      allocate(ju(1:ni),stat=ierr)                                         2621
      jerr=jerr+ierr                                                       2622
      allocate(vq(1:ni),stat=ierr)                                         2622
      jerr=jerr+ierr                                                       2623
      allocate(xm(1:ni),stat=ierr)                                         2623
      jerr=jerr+ierr                                                       2624
      if(isd .le. 0)goto 17041                                             2624
      allocate(xs(1:ni),stat=ierr)                                         2624
      jerr=jerr+ierr                                                       2624
17041 continue                                                             2625
      if(jerr.ne.0) return                                                 2626
      call chkvars(no,ni,x,ju)                                             2627
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2628
      if(maxval(ju) .gt. 0)goto 17061                                      2628
      jerr=7777                                                            2628
      go to 11910                                                          2628
17061 continue                                                             2629
      vq=max(0.0,vp)                                                       2629
      vq=vq*ni/sum(vq)                                                     2630
      ww=max(0.0,w)                                                        2630
      sw=sum(ww)                                                           2630
      if(sw .gt. 0.0)goto 17081                                            2630
      jerr=9999                                                            2630
      go to 11910                                                          2630
17081 continue                                                             2631
      ww=ww/sw                                                             2632
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             2633
      if(isd .le. 0)goto 17101                                             2633
17110 do 17111 j=1,ni                                                      2633
      cl(:,j)=cl(:,j)*xs(j)                                                2633
17111 continue                                                             2633
17112 continue                                                             2633
17101 continue                                                             2634
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,t   2636 
     *hr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11910                                            2636
      dev0=2.0*sw*dev0                                                     2637
17120 do 17121 k=1,lmu                                                     2637
      nk=nin(k)                                                            2638
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2639
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2640
17121 continue                                                             2641
17122 continue                                                             2641
11910 continue                                                             2641
      deallocate(ww,ju,vq,xm)                                              2641
      if(isd.gt.0) deallocate(xs)                                          2642
      return                                                               2643
      end                                                                  2644
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,cl,ne,nx,nlam,flmin,u   2646 
     *lam,shri,  isd,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)                    2647
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)               2648
      integer ju(ni),m(nx),kin(nlam)                                       2649
      real, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga                    
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2653
      sml=sml*10.0                                                         2654
      allocate(a(1:ni),stat=jerr)                                          2655
      allocate(as(1:ni),stat=ierr)                                         2655
      jerr=jerr+ierr                                                       2656
      allocate(t(1:no),stat=ierr)                                          2656
      jerr=jerr+ierr                                                       2657
      allocate(mm(1:ni),stat=ierr)                                         2657
      jerr=jerr+ierr                                                       2658
      allocate(ga(1:ni),stat=ierr)                                         2658
      jerr=jerr+ierr                                                       2659
      allocate(ixx(1:ni),stat=ierr)                                        2659
      jerr=jerr+ierr                                                       2660
      allocate(wr(1:no),stat=ierr)                                         2660
      jerr=jerr+ierr                                                       2661
      allocate(v(1:ni),stat=ierr)                                          2661
      jerr=jerr+ierr                                                       2662
      allocate(w(1:no),stat=ierr)                                          2662
      jerr=jerr+ierr                                                       2663
      allocate(f(1:no),stat=ierr)                                          2663
      jerr=jerr+ierr                                                       2664
      if(jerr.ne.0) return                                                 2665
      bta=parm                                                             2665
      omb=1.0-bta                                                          2666
      t=q*y                                                                2666
      yb=sum(t)                                                            2666
      fmax=log(huge(bta)*0.1)                                              2667
      if(nonzero(no,g) .ne. 0)goto 17141                                   2667
      w=q*yb                                                               2667
      az=log(yb)                                                           2667
      f=az                                                                 2667
      dv0=yb*(log(yb)-1.0)                                                 2667
      goto 17151                                                           2668
17141 continue                                                             2668
      w=q*exp(sign(min(abs(g),fmax),g))                                    2668
      v0=sum(w)                                                            2668
      eaz=yb/v0                                                            2669
      w=eaz*w                                                              2669
      az=log(eaz)                                                          2669
      f=az+g                                                               2670
      dv0=dot_product(t,g)-yb*(1.0-az)                                     2671
17151 continue                                                             2672
17131 continue                                                             2672
      a=0.0                                                                2672
      as=0.0                                                               2672
      wr=t-w                                                               2672
      v0=yb                                                                2672
      dvr=-yb                                                              2673
17160 do 17161 i=1,no                                                      2673
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               2673
17161 continue                                                             2673
17162 continue                                                             2673
      dvr=dvr-dv0                                                          2673
      dev0=dvr                                                             2674
      if(flmin .ge. 1.0)goto 17181                                         2674
      eqs=max(eps,flmin)                                                   2674
      alf=eqs**(1.0/(nlam-1))                                              2674
17181 continue                                                             2675
      m=0                                                                  2675
      mm=0                                                                 2675
      nlp=0                                                                2675
      nin=nlp                                                              2675
      mnl=min(mnlam,nlam)                                                  2675
      shr=shri*dev0                                                        2675
      ixx=0                                                                2675
      al=0.0                                                               2676
17190 do 17191 j=1,ni                                                      2676
      if(ju(j).eq.0)goto 17191                                             2676
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2676
17191 continue                                                             2677
17192 continue                                                             2677
17200 do 17201 ilm=1,nlam                                                  2677
      al0=al                                                               2678
      if(flmin .lt. 1.0)goto 17221                                         2678
      al=ulam(ilm)                                                         2678
      goto 17211                                                           2679
17221 if(ilm .le. 2)goto 17231                                             2679
      al=al*alf                                                            2679
      goto 17211                                                           2680
17231 if(ilm .ne. 1)goto 17241                                             2680
      al=big                                                               2680
      goto 17251                                                           2681
17241 continue                                                             2681
      al0=0.0                                                              2682
17260 do 17261 j=1,ni                                                      2682
      if(ju(j).eq.0)goto 17261                                             2682
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2682
17261 continue                                                             2683
17262 continue                                                             2683
      al0=al0/max(bta,1.0e-3)                                              2683
      al=alf*al0                                                           2684
17251 continue                                                             2685
17211 continue                                                             2685
      al2=al*omb                                                           2685
      al1=al*bta                                                           2685
      tlam=bta*(2.0*al-al0)                                                2686
17270 do 17271 k=1,ni                                                      2686
      if(ixx(k).eq.1)goto 17271                                            2686
      if(ju(k).eq.0)goto 17271                                             2687
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2688
17271 continue                                                             2689
17272 continue                                                             2689
10740 continue                                                             2690
17280 continue                                                             2690
17281 continue                                                             2690
      az0=az                                                               2691
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2692
17290 do 17291 j=1,ni                                                      2692
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                        2692
17291 continue                                                             2693
17292 continue                                                             2693
17300 continue                                                             2693
17301 continue                                                             2693
      nlp=nlp+1                                                            2693
      dlx=0.0                                                              2694
17310 do 17311 k=1,ni                                                      2694
      if(ixx(k).eq.0)goto 17311                                            2694
      ak=a(k)                                                              2695
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2695
      au=abs(u)-vp(k)*al1                                                  2696
      if(au .gt. 0.0)goto 17331                                            2696
      a(k)=0.0                                                             2696
      goto 17341                                                           2697
17331 continue                                                             2698
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           2699
17341 continue                                                             2700
17321 continue                                                             2700
      if(a(k).eq.ak)goto 17311                                             2700
      d=a(k)-ak                                                            2700
      dlx=max(dlx,v(k)*d**2)                                               2701
      wr=wr-d*w*x(:,k)                                                     2701
      f=f+d*x(:,k)                                                         2702
      if(mm(k) .ne. 0)goto 17361                                           2702
      nin=nin+1                                                            2702
      if(nin.gt.nx)goto 17312                                              2703
      mm(k)=nin                                                            2703
      m(nin)=k                                                             2704
17361 continue                                                             2705
17311 continue                                                             2706
17312 continue                                                             2706
      if(nin.gt.nx)goto 17302                                              2706
      d=sum(wr)/v0                                                         2707
      az=az+d                                                              2707
      dlx=max(dlx,v0*d**2)                                                 2707
      wr=wr-d*w                                                            2707
      f=f+d                                                                2708
      if(dlx.lt.shr)goto 17302                                             2708
      if(nlp .le. maxit)goto 17381                                         2708
      jerr=-ilm                                                            2708
      return                                                               2708
17381 continue                                                             2709
17390 continue                                                             2709
17391 continue                                                             2709
      nlp=nlp+1                                                            2709
      dlx=0.0                                                              2710
17400 do 17401 l=1,nin                                                     2710
      k=m(l)                                                               2710
      ak=a(k)                                                              2711
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2711
      au=abs(u)-vp(k)*al1                                                  2712
      if(au .gt. 0.0)goto 17421                                            2712
      a(k)=0.0                                                             2712
      goto 17431                                                           2713
17421 continue                                                             2714
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           2715
17431 continue                                                             2716
17411 continue                                                             2716
      if(a(k).eq.ak)goto 17401                                             2716
      d=a(k)-ak                                                            2716
      dlx=max(dlx,v(k)*d**2)                                               2717
      wr=wr-d*w*x(:,k)                                                     2717
      f=f+d*x(:,k)                                                         2719
17401 continue                                                             2719
17402 continue                                                             2719
      d=sum(wr)/v0                                                         2719
      az=az+d                                                              2719
      dlx=max(dlx,v0*d**2)                                                 2719
      wr=wr-d*w                                                            2719
      f=f+d                                                                2720
      if(dlx.lt.shr)goto 17392                                             2720
      if(nlp .le. maxit)goto 17451                                         2720
      jerr=-ilm                                                            2720
      return                                                               2720
17451 continue                                                             2721
      goto 17391                                                           2722
17392 continue                                                             2722
      goto 17301                                                           2723
17302 continue                                                             2723
      if(nin.gt.nx)goto 17282                                              2724
      w=q*exp(sign(min(abs(f),fmax),f))                                    2724
      v0=sum(w)                                                            2724
      wr=t-w                                                               2725
      if(v0*(az-az0)**2 .ge. shr)goto 17471                                2725
      ix=0                                                                 2726
17480 do 17481 j=1,nin                                                     2726
      k=m(j)                                                               2727
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 17481                            2727
      ix=1                                                                 2727
      goto 17482                                                           2728
17481 continue                                                             2729
17482 continue                                                             2729
      if(ix .ne. 0)goto 17501                                              2730
17510 do 17511 k=1,ni                                                      2730
      if(ixx(k).eq.1)goto 17511                                            2730
      if(ju(k).eq.0)goto 17511                                             2731
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2732
      if(ga(k) .le. al1*vp(k))goto 17531                                   2732
      ixx(k)=1                                                             2732
      ix=1                                                                 2732
17531 continue                                                             2733
17511 continue                                                             2734
17512 continue                                                             2734
      if(ix.eq.1) go to 10740                                              2735
      goto 17282                                                           2736
17501 continue                                                             2737
17471 continue                                                             2738
      goto 17281                                                           2739
17282 continue                                                             2739
      if(nin .le. nx)goto 17551                                            2739
      jerr=-10000-ilm                                                      2739
      goto 17202                                                           2739
17551 continue                                                             2740
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2740
      kin(ilm)=nin                                                         2741
      a0(ilm)=az                                                           2741
      alm(ilm)=al                                                          2741
      lmu=ilm                                                              2742
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               2743
      if(ilm.lt.mnl)goto 17201                                             2743
      if(flmin.ge.1.0)goto 17201                                           2744
      me=0                                                                 2744
17560 do 17561 j=1,nin                                                     2744
      if(ca(j,ilm).ne.0.0) me=me+1                                         2744
17561 continue                                                             2744
17562 continue                                                             2744
      if(me.gt.ne)goto 17202                                               2745
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17202              2746
      if(dev(ilm).gt.devmax)goto 17202                                     2747
17201 continue                                                             2748
17202 continue                                                             2748
      g=f                                                                  2749
11910 continue                                                             2749
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                                2750
      return                                                               2751
      end                                                                  2752
      function nonzero(n,v)                                                2753
      real v(n)                                                            2754
      nonzero=0                                                            2754
17570 do 17571 i=1,n                                                       2754
      if(v(i) .eq. 0.0)goto 17591                                          2754
      nonzero=1                                                            2754
      return                                                               2754
17591 continue                                                             2754
17571 continue                                                             2755
17572 continue                                                             2755
      return                                                               2756
      end                                                                  2757
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               2758
      real a(nx,lmu),b(ni,lmu)                                             2758
      integer ia(nx),nin(lmu)                                              2759
17600 do 17601 lam=1,lmu                                                   2759
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        2759
17601 continue                                                             2760
17602 continue                                                             2760
      return                                                               2761
      end                                                                  2762
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           2763
      real a(nx,nc,lmu),b(ni,nc,lmu)                                       2763
      integer ia(nx),nin(lmu)                                              2764
17610 do 17611 lam=1,lmu                                                   2764
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             2764
17611 continue                                                             2765
17612 continue                                                             2765
      return                                                               2766
      end                                                                  2767
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               2768
      real x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)       2769
      real, dimension (:), allocatable :: w                                     
      if(minval(y) .ge. 0.0)goto 17631                                     2772
      jerr=8888                                                            2772
      return                                                               2772
17631 continue                                                             2773
      allocate(w(1:no),stat=jerr)                                          2773
      if(jerr.ne.0) return                                                 2774
      w=max(0.0,q)                                                         2774
      sw=sum(w)                                                            2774
      if(sw .gt. 0.0)goto 17651                                            2774
      jerr=9999                                                            2774
      go to 11910                                                          2774
17651 continue                                                             2775
      yb=dot_product(w,y)/sw                                               2775
      fmax=log(huge(y(1))*0.1)                                             2776
17660 do 17661 lam=1,nlam                                                  2776
      s=0.0                                                                2777
17670 do 17671 i=1,no                                                      2777
      if(w(i).le.0.0)goto 17671                                            2778
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          2779
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      2780
17671 continue                                                             2781
17672 continue                                                             2781
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2782
17661 continue                                                             2783
17662 continue                                                             2783
11910 continue                                                             2783
      deallocate(w)                                                        2784
      return                                                               2785
      end                                                                  2786
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,cl,ne,nx,nlam   2788 
     *,flmin,  ulam,thr,isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr
     *)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)               2789
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2790
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           2791
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 17691                                    2795
      jerr=10000                                                           2795
      return                                                               2795
17691 continue                                                             2796
      if(minval(y) .ge. 0.0)goto 17711                                     2796
      jerr=8888                                                            2796
      return                                                               2796
17711 continue                                                             2797
      allocate(ww(1:no),stat=jerr)                                         2798
      allocate(ju(1:ni),stat=ierr)                                         2798
      jerr=jerr+ierr                                                       2799
      allocate(vq(1:ni),stat=ierr)                                         2799
      jerr=jerr+ierr                                                       2800
      allocate(xm(1:ni),stat=ierr)                                         2800
      jerr=jerr+ierr                                                       2801
      allocate(xs(1:ni),stat=ierr)                                         2801
      jerr=jerr+ierr                                                       2802
      if(jerr.ne.0) return                                                 2803
      call spchkvars(no,ni,x,ix,ju)                                        2804
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2805
      if(maxval(ju) .gt. 0)goto 17731                                      2805
      jerr=7777                                                            2805
      go to 11910                                                          2805
17731 continue                                                             2806
      vq=max(0.0,vp)                                                       2806
      vq=vq*ni/sum(vq)                                                     2807
      ww=max(0.0,w)                                                        2807
      sw=sum(ww)                                                           2807
      if(sw .gt. 0.0)goto 17751                                            2807
      jerr=9999                                                            2807
      go to 11910                                                          2807
17751 continue                                                             2808
      ww=ww/sw                                                             2809
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     2810
      if(isd .le. 0)goto 17771                                             2810
17780 do 17781 j=1,ni                                                      2810
      cl(:,j)=cl(:,j)*xs(j)                                                2810
17781 continue                                                             2810
17782 continue                                                             2810
17771 continue                                                             2811
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,flmi   2813 
     *n,ulam,thr,  isd,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jer
     *r)
      if(jerr.gt.0) go to 11910                                            2813
      dev0=2.0*sw*dev0                                                     2814
17790 do 17791 k=1,lmu                                                     2814
      nk=nin(k)                                                            2815
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2816
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2817
17791 continue                                                             2818
17792 continue                                                             2818
11910 continue                                                             2818
      deallocate(ww,ju,vq,xm,xs)                                           2819
      return                                                               2820
      end                                                                  2821
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,cl,ne,nx,nlam   2823 
     *,flmin,ulam,  shri,isd,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,alm,nl
     *p,jerr)
      real x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),xs(ni)          2824
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)               2825
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2826
      real, dimension (:), allocatable :: qy,t,w,wr,v,a,as,xm,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2830
      sml=sml*10.0                                                         2831
      allocate(a(1:ni),stat=jerr)                                          2832
      allocate(as(1:ni),stat=ierr)                                         2832
      jerr=jerr+ierr                                                       2833
      allocate(t(1:no),stat=ierr)                                          2833
      jerr=jerr+ierr                                                       2834
      allocate(mm(1:ni),stat=ierr)                                         2834
      jerr=jerr+ierr                                                       2835
      allocate(ga(1:ni),stat=ierr)                                         2835
      jerr=jerr+ierr                                                       2836
      allocate(ixx(1:ni),stat=ierr)                                        2836
      jerr=jerr+ierr                                                       2837
      allocate(wr(1:no),stat=ierr)                                         2837
      jerr=jerr+ierr                                                       2838
      allocate(v(1:ni),stat=ierr)                                          2838
      jerr=jerr+ierr                                                       2839
      allocate(xm(1:ni),stat=ierr)                                         2839
      jerr=jerr+ierr                                                       2840
      allocate(w(1:no),stat=ierr)                                          2840
      jerr=jerr+ierr                                                       2841
      allocate(qy(1:no),stat=ierr)                                         2841
      jerr=jerr+ierr                                                       2842
      if(jerr.ne.0) return                                                 2843
      bta=parm                                                             2843
      omb=1.0-bta                                                          2843
      fmax=log(huge(bta)*0.1)                                              2844
      qy=q*y                                                               2844
      yb=sum(qy)                                                           2845
      if(nonzero(no,g) .ne. 0)goto 17811                                   2845
      w=q*yb                                                               2845
      az=log(yb)                                                           2845
      uu=az                                                                2846
      xm=yb*xb                                                             2846
      t=0.0                                                                2846
      dv0=yb*(log(yb)-1.0)                                                 2847
      goto 17821                                                           2848
17811 continue                                                             2848
      w=q*exp(sign(min(abs(g),fmax),g))                                    2848
      ww=sum(w)                                                            2848
      eaz=yb/ww                                                            2849
      w=eaz*w                                                              2849
      az=log(eaz)                                                          2849
      uu=az                                                                2849
      t=g                                                                  2849
      dv0=dot_product(qy,g)-yb*(1.0-az)                                    2850
17830 do 17831 j=1,ni                                                      2850
      if(ju(j).eq.0)goto 17831                                             2850
      jb=ix(j)                                                             2850
      je=ix(j+1)-1                                                         2851
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2852
17831 continue                                                             2853
17832 continue                                                             2853
17821 continue                                                             2854
17801 continue                                                             2854
      tt=yb*uu                                                             2854
      ww=yb                                                                2854
      wr=qy-q*(yb*(1.0-uu))                                                2854
      a=0.0                                                                2854
      as=0.0                                                               2855
      dvr=-yb                                                              2856
17840 do 17841 i=1,no                                                      2856
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             2856
17841 continue                                                             2856
17842 continue                                                             2856
      dvr=dvr-dv0                                                          2856
      dev0=dvr                                                             2857
      if(flmin .ge. 1.0)goto 17861                                         2857
      eqs=max(eps,flmin)                                                   2857
      alf=eqs**(1.0/(nlam-1))                                              2857
17861 continue                                                             2858
      m=0                                                                  2858
      mm=0                                                                 2858
      nlp=0                                                                2858
      nin=nlp                                                              2858
      mnl=min(mnlam,nlam)                                                  2858
      shr=shri*dev0                                                        2858
      al=0.0                                                               2858
      ixx=0                                                                2859
17870 do 17871 j=1,ni                                                      2859
      if(ju(j).eq.0)goto 17871                                             2860
      jb=ix(j)                                                             2860
      je=ix(j+1)-1                                                         2861
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   2863 
     *)-xb(j)*tt)/xs(j)
17871 continue                                                             2864
17872 continue                                                             2864
17880 do 17881 ilm=1,nlam                                                  2864
      al0=al                                                               2865
      if(flmin .lt. 1.0)goto 17901                                         2865
      al=ulam(ilm)                                                         2865
      goto 17891                                                           2866
17901 if(ilm .le. 2)goto 17911                                             2866
      al=al*alf                                                            2866
      goto 17891                                                           2867
17911 if(ilm .ne. 1)goto 17921                                             2867
      al=big                                                               2867
      goto 17931                                                           2868
17921 continue                                                             2868
      al0=0.0                                                              2869
17940 do 17941 j=1,ni                                                      2869
      if(ju(j).eq.0)goto 17941                                             2869
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2869
17941 continue                                                             2870
17942 continue                                                             2870
      al0=al0/max(bta,1.0e-3)                                              2870
      al=alf*al0                                                           2871
17931 continue                                                             2872
17891 continue                                                             2872
      al2=al*omb                                                           2872
      al1=al*bta                                                           2872
      tlam=bta*(2.0*al-al0)                                                2873
17950 do 17951 k=1,ni                                                      2873
      if(ixx(k).eq.1)goto 17951                                            2873
      if(ju(k).eq.0)goto 17951                                             2874
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2875
17951 continue                                                             2876
17952 continue                                                             2876
10740 continue                                                             2877
17960 continue                                                             2877
17961 continue                                                             2877
      az0=az                                                               2878
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2879
17970 do 17971 j=1,ni                                                      2879
      if(ixx(j).eq.0)goto 17971                                            2879
      jb=ix(j)                                                             2879
      je=ix(j+1)-1                                                         2880
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2881
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   2883 
     *b(j)**2)/xs(j)**2
17971 continue                                                             2884
17972 continue                                                             2884
17980 continue                                                             2884
17981 continue                                                             2884
      nlp=nlp+1                                                            2885
      dlx=0.0                                                              2886
17990 do 17991 k=1,ni                                                      2886
      if(ixx(k).eq.0)goto 17991                                            2886
      jb=ix(k)                                                             2886
      je=ix(k+1)-1                                                         2886
      ak=a(k)                                                              2887
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2889 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2890
      if(au .gt. 0.0)goto 18011                                            2890
      a(k)=0.0                                                             2890
      goto 18021                                                           2891
18011 continue                                                             2892
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           2893
18021 continue                                                             2894
18001 continue                                                             2894
      if(a(k).eq.ak)goto 17991                                             2895
      if(mm(k) .ne. 0)goto 18041                                           2895
      nin=nin+1                                                            2895
      if(nin.gt.nx)goto 17992                                              2896
      mm(k)=nin                                                            2896
      m(nin)=k                                                             2897
18041 continue                                                             2898
      d=a(k)-ak                                                            2898
      dlx=max(dlx,v(k)*d**2)                                               2898
      dv=d/xs(k)                                                           2899
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2900
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2901
      uu=uu-dv*xb(k)                                                       2901
      tt=tt-dv*xm(k)                                                       2902
17991 continue                                                             2903
17992 continue                                                             2903
      if(nin.gt.nx)goto 17982                                              2903
      d=tt/ww-uu                                                           2904
      az=az+d                                                              2904
      dlx=max(dlx,ww*d**2)                                                 2904
      uu=uu+d                                                              2905
      if(dlx.lt.shr)goto 17982                                             2905
      if(nlp .le. maxit)goto 18061                                         2905
      jerr=-ilm                                                            2905
      return                                                               2905
18061 continue                                                             2906
18070 continue                                                             2906
18071 continue                                                             2906
      nlp=nlp+1                                                            2906
      dlx=0.0                                                              2907
18080 do 18081 l=1,nin                                                     2907
      k=m(l)                                                               2908
      jb=ix(k)                                                             2908
      je=ix(k+1)-1                                                         2908
      ak=a(k)                                                              2909
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2911 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2912
      if(au .gt. 0.0)goto 18101                                            2912
      a(k)=0.0                                                             2912
      goto 18111                                                           2913
18101 continue                                                             2914
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           2915
18111 continue                                                             2916
18091 continue                                                             2916
      if(a(k).eq.ak)goto 18081                                             2916
      d=a(k)-ak                                                            2916
      dlx=max(dlx,v(k)*d**2)                                               2917
      dv=d/xs(k)                                                           2917
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2918
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2919
      uu=uu-dv*xb(k)                                                       2919
      tt=tt-dv*xm(k)                                                       2920
18081 continue                                                             2921
18082 continue                                                             2921
      d=tt/ww-uu                                                           2921
      az=az+d                                                              2921
      dlx=max(dlx,ww*d**2)                                                 2921
      uu=uu+d                                                              2922
      if(dlx.lt.shr)goto 18072                                             2922
      if(nlp .le. maxit)goto 18131                                         2922
      jerr=-ilm                                                            2922
      return                                                               2922
18131 continue                                                             2923
      goto 18071                                                           2924
18072 continue                                                             2924
      goto 17981                                                           2925
17982 continue                                                             2925
      if(nin.gt.nx)goto 17962                                              2926
      euu=exp(sign(min(abs(uu),fmax),uu))                                  2927
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                2927
      ww=sum(w)                                                            2928
      wr=qy-w*(1.0-uu)                                                     2928
      tt=sum(wr)                                                           2929
      if(ww*(az-az0)**2 .ge. shr)goto 18151                                2929
      kx=0                                                                 2930
18160 do 18161 j=1,nin                                                     2930
      k=m(j)                                                               2931
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 18161                            2931
      kx=1                                                                 2931
      goto 18162                                                           2932
18161 continue                                                             2933
18162 continue                                                             2933
      if(kx .ne. 0)goto 18181                                              2934
18190 do 18191 j=1,ni                                                      2934
      if(ixx(j).eq.1)goto 18191                                            2934
      if(ju(j).eq.0)goto 18191                                             2935
      jb=ix(j)                                                             2935
      je=ix(j+1)-1                                                         2936
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2937
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   2939 
     *)-xb(j)*tt)/xs(j)
      if(ga(j) .le. al1*vp(j))goto 18211                                   2939
      ixx(j)=1                                                             2939
      kx=1                                                                 2939
18211 continue                                                             2940
18191 continue                                                             2941
18192 continue                                                             2941
      if(kx.eq.1) go to 10740                                              2942
      goto 17962                                                           2943
18181 continue                                                             2944
18151 continue                                                             2945
      goto 17961                                                           2946
17962 continue                                                             2946
      if(nin .le. nx)goto 18231                                            2946
      jerr=-10000-ilm                                                      2946
      goto 17882                                                           2946
18231 continue                                                             2947
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2947
      kin(ilm)=nin                                                         2948
      a0(ilm)=az                                                           2948
      alm(ilm)=al                                                          2948
      lmu=ilm                                                              2949
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        2950
      if(ilm.lt.mnl)goto 17881                                             2950
      if(flmin.ge.1.0)goto 17881                                           2951
      me=0                                                                 2951
18240 do 18241 j=1,nin                                                     2951
      if(ca(j,ilm).ne.0.0) me=me+1                                         2951
18241 continue                                                             2951
18242 continue                                                             2951
      if(me.gt.ne)goto 17882                                               2952
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17882              2953
      if(dev(ilm).gt.devmax)goto 17882                                     2954
17881 continue                                                             2955
17882 continue                                                             2955
      g=t+uu                                                               2956
11910 continue                                                             2956
      deallocate(t,w,wr,v,a,qy,xm,as,mm,ga,ixx)                            2957
      return                                                               2958
      end                                                                  2959
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       2960
      real x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)           2961
      integer ix(*),jx(*)                                                  2962
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 18261                                     2965
      jerr=8888                                                            2965
      return                                                               2965
18261 continue                                                             2966
      allocate(w(1:no),stat=jerr)                                          2967
      allocate(f(1:no),stat=ierr)                                          2967
      jerr=jerr+ierr                                                       2968
      if(jerr.ne.0) return                                                 2969
      w=max(0.0,q)                                                         2969
      sw=sum(w)                                                            2969
      if(sw .gt. 0.0)goto 18281                                            2969
      jerr=9999                                                            2969
      go to 11910                                                          2969
18281 continue                                                             2970
      yb=dot_product(w,y)/sw                                               2970
      fmax=log(huge(y(1))*0.1)                                             2971
18290 do 18291 lam=1,nlam                                                  2971
      f=a0(lam)                                                            2972
18300 do 18301 j=1,ni                                                      2972
      if(a(j,lam).eq.0.0)goto 18301                                        2972
      jb=ix(j)                                                             2972
      je=ix(j+1)-1                                                         2973
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          2974
18301 continue                                                             2975
18302 continue                                                             2975
      f=f+g                                                                2976
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2977
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2978
18291 continue                                                             2979
18292 continue                                                             2979
11910 continue                                                             2979
      deallocate(w,f)                                                      2980
      return                                                               2981
      end                                                                  2982
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   2983 
     *jerr)
      real x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(nlam)          2984
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 2985
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 18321                                     2988
      jerr=8888                                                            2988
      return                                                               2988
18321 continue                                                             2989
      allocate(w(1:no),stat=jerr)                                          2990
      allocate(f(1:no),stat=ierr)                                          2990
      jerr=jerr+ierr                                                       2991
      if(jerr.ne.0) return                                                 2992
      w=max(0.0,q)                                                         2992
      sw=sum(w)                                                            2992
      if(sw .gt. 0.0)goto 18341                                            2992
      jerr=9999                                                            2992
      go to 11910                                                          2992
18341 continue                                                             2993
      yb=dot_product(w,y)/sw                                               2993
      fmax=log(huge(y(1))*0.1)                                             2994
18350 do 18351 lam=1,nlam                                                  2994
      f=a0(lam)                                                            2995
18360 do 18361 k=1,nin(lam)                                                2995
      j=ia(k)                                                              2995
      jb=ix(j)                                                             2995
      je=ix(j+1)-1                                                         2996
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         2997
18361 continue                                                             2998
18362 continue                                                             2998
      f=f+g                                                                2999
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   3000
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                3001
18351 continue                                                             3002
18352 continue                                                             3002
11910 continue                                                             3002
      deallocate(w,f)                                                      3003
      return                                                               3004
      end                                                                  3005
      subroutine multelnet  (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flm   3008 
     *in,ulam,thr,isd,jsd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no,nr),w(no),vp(ni),ca(nx,nr,nlam)                   3009
      real ulam(nlam),a0(nr,nlam),rsq(nlam),alm(nlam),cl(2,ni)             3010
      integer jd(*),ia(nx),nin(nlam)                                       3011
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 18381                                    3014
      jerr=10000                                                           3014
      return                                                               3014
18381 continue                                                             3015
      allocate(vq(1:ni),stat=jerr)                                         3015
      if(jerr.ne.0) return                                                 3016
      vq=max(0.0,vp)                                                       3016
      vq=vq*ni/sum(vq)                                                     3017
      call multelnetn(parm,no,ni,nr,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam   3019 
     *,thr,isd,  jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      deallocate(vq)                                                       3020
      return                                                               3021
      end                                                                  3022
      subroutine multelnetn (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flm   3024 
     *in,ulam,thr,  isd,jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no,nr),w(no),ulam(nlam),cl(2,ni)              3025
      real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                  3026
      integer jd(*),ia(nx),nin(nlam)                                       3027
      real, dimension (:), allocatable :: xm,xs,xv,ym,ys                        
      integer, dimension (:), allocatable :: ju                                 
      real, dimension (:,:,:), allocatable :: clt                               
      allocate(clt(1:2,1:nr,1:ni),stat=jerr);                                   
      allocate(xm(1:ni),stat=ierr)                                         3033
      jerr=jerr+ierr                                                       3034
      allocate(xs(1:ni),stat=ierr)                                         3034
      jerr=jerr+ierr                                                       3035
      allocate(ym(1:nr),stat=ierr)                                         3035
      jerr=jerr+ierr                                                       3036
      allocate(ys(1:nr),stat=ierr)                                         3036
      jerr=jerr+ierr                                                       3037
      allocate(ju(1:ni),stat=ierr)                                         3037
      jerr=jerr+ierr                                                       3038
      allocate(xv(1:ni),stat=ierr)                                         3038
      jerr=jerr+ierr                                                       3039
      if(jerr.ne.0) return                                                 3040
      call chkvars(no,ni,x,ju)                                             3041
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3042
      if(maxval(ju) .gt. 0)goto 18401                                      3042
      jerr=7777                                                            3042
      return                                                               3042
18401 continue                                                             3043
      call multstandard1(no,ni,nr,x,y,w,isd,jsd,ju,xm,xs,ym,ys,xv,ys0,je   3044 
     *rr)
      if(jerr.ne.0) return                                                 3045
18410 do 18411 j=1,ni                                                      3045
18420 do 18421 k=1,nr                                                      3045
18430 do 18431 i=1,2                                                       3045
      clt(i,k,j)=cl(i,j)                                                   3045
18431 continue                                                             3045
18432 continue                                                             3045
18421 continue                                                             3045
18422 continue                                                             3045
18411 continue                                                             3046
18412 continue                                                             3046
      if(isd .le. 0)goto 18451                                             3046
18460 do 18461 j=1,ni                                                      3046
18470 do 18471 k=1,nr                                                      3046
18480 do 18481 i=1,2                                                       3046
      clt(i,k,j)=clt(i,k,j)*xs(j)                                          3046
18481 continue                                                             3046
18482 continue                                                             3046
18471 continue                                                             3046
18472 continue                                                             3046
18461 continue                                                             3046
18462 continue                                                             3046
18451 continue                                                             3047
      if(jsd .le. 0)goto 18501                                             3047
18510 do 18511 j=1,ni                                                      3047
18520 do 18521 k=1,nr                                                      3047
18530 do 18531 i=1,2                                                       3047
      clt(i,k,j)=clt(i,k,j)/ys(k)                                          3047
18531 continue                                                             3047
18532 continue                                                             3047
18521 continue                                                             3047
18522 continue                                                             3047
18511 continue                                                             3047
18512 continue                                                             3047
18501 continue                                                             3048
      call multelnet2(parm,ni,nr,ju,vp,clt,y,no,ne,nx,x,nlam,flmin,ulam,   3050 
     *thr,maxit,xv,  ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 3051
18540 do 18541 k=1,lmu                                                     3051
      nk=nin(k)                                                            3052
18550 do 18551 j=1,nr                                                      3053
18560 do 18561 l=1,nk                                                      3053
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  3053
18561 continue                                                             3054
18562 continue                                                             3054
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 3055
18551 continue                                                             3056
18552 continue                                                             3056
18541 continue                                                             3057
18542 continue                                                             3057
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                    3058
      return                                                               3059
      end                                                                  3060
      subroutine multstandard1 (no,ni,nr,x,y,w,isd,jsd,ju,xm,xs,ym,ys,xv   3061 
     *,ys0,jerr)
      real x(no,ni),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),ys(nr)      3062
      integer ju(ni)                                                       3063
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                          3066
      if(jerr.ne.0) return                                                 3067
      w=w/sum(w)                                                           3067
      v=sqrt(w)                                                            3068
18570 do 18571 j=1,ni                                                      3068
      if(ju(j).eq.0)goto 18571                                             3069
      xm(j)=dot_product(w,x(:,j))                                          3069
      x(:,j)=v*(x(:,j)-xm(j))                                              3070
      xv(j)=dot_product(x(:,j),x(:,j))                                     3070
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       3071
18571 continue                                                             3072
18572 continue                                                             3072
      if(isd .ne. 0)goto 18591                                             3072
      xs=1.0                                                               3072
      goto 18601                                                           3073
18591 continue                                                             3073
18610 do 18611 j=1,ni                                                      3073
      if(ju(j).eq.0)goto 18611                                             3073
      x(:,j)=x(:,j)/xs(j)                                                  3073
18611 continue                                                             3074
18612 continue                                                             3074
      xv=1.0                                                               3075
18601 continue                                                             3076
18581 continue                                                             3076
      ys0=0.0                                                              3077
18620 do 18621 j=1,nr                                                      3078
      ym(j)=dot_product(w,y(:,j))                                          3078
      y(:,j)=v*(y(:,j)-ym(j))                                              3079
      z=dot_product(y(:,j),y(:,j))                                         3080
      if(jsd .le. 0)goto 18641                                             3080
      ys(j)=sqrt(z)                                                        3080
      y(:,j)=y(:,j)/ys(j)                                                  3080
      goto 18651                                                           3081
18641 continue                                                             3081
      ys0=ys0+z                                                            3081
18651 continue                                                             3082
18631 continue                                                             3082
18621 continue                                                             3083
18622 continue                                                             3083
      if(jsd .ne. 0)goto 18671                                             3083
      ys=1.0                                                               3083
      goto 18681                                                           3083
18671 continue                                                             3083
      ys0=nr                                                               3083
18681 continue                                                             3084
18661 continue                                                             3084
      deallocate(v)                                                        3085
      return                                                               3086
      end                                                                  3087
      subroutine multelnet2(beta,ni,nr,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,   3089 
     *ulam,thri,  maxit,xv,ys0,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      real vp(ni),y(no,nr),x(no,ni),ulam(nlam),ao(nx,nr,nlam)              3090
      real rsqo(nlam),almo(nlam),xv(ni),cl(2,nr,ni)                        3091
      integer ju(ni),ia(nx),kin(nlam)                                      3092
      real, dimension (:), allocatable :: g,gk,del,gj                           
      integer, dimension (:), allocatable :: mm,ix,isc                          
      real, dimension (:,:), allocatable :: a                                   
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               3099
      allocate(gj(1:nr),stat=ierr)                                         3099
      jerr=jerr+ierr                                                       3100
      allocate(gk(1:nr),stat=ierr)                                         3100
      jerr=jerr+ierr                                                       3101
      allocate(del(1:nr),stat=ierr)                                        3101
      jerr=jerr+ierr                                                       3102
      allocate(mm(1:ni),stat=ierr)                                         3102
      jerr=jerr+ierr                                                       3103
      allocate(g(1:ni),stat=ierr)                                          3103
      jerr=jerr+ierr                                                       3104
      allocate(ix(1:ni),stat=ierr)                                         3104
      jerr=jerr+ierr                                                       3105
      allocate(isc(1:nr),stat=ierr)                                        3105
      jerr=jerr+ierr                                                       3106
      if(jerr.ne.0) return                                                 3107
      bta=beta                                                             3107
      omb=1.0-bta                                                          3107
      ix=0                                                                 3107
      thr=thri*ys0/nr                                                      3108
      if(flmin .ge. 1.0)goto 18701                                         3108
      eqs=max(eps,flmin)                                                   3108
      alf=eqs**(1.0/(nlam-1))                                              3108
18701 continue                                                             3109
      rsq=ys0                                                              3109
      a=0.0                                                                3109
      mm=0                                                                 3109
      nlp=0                                                                3109
      nin=nlp                                                              3109
      iz=0                                                                 3109
      mnl=min(mnlam,nlam)                                                  3109
      alm=0.0                                                              3110
18710 do 18711 j=1,ni                                                      3110
      if(ju(j).eq.0)goto 18711                                             3110
      g(j)=0.0                                                             3111
18720 do 18721 k=1,nr                                                      3111
      g(j)=g(j)+dot_product(y(:,k),x(:,j))**2                              3111
18721 continue                                                             3112
18722 continue                                                             3112
      g(j)=sqrt(g(j))                                                      3113
18711 continue                                                             3114
18712 continue                                                             3114
18730 do 18731 m=1,nlam                                                    3115
      if(flmin .lt. 1.0)goto 18751                                         3115
      alm=ulam(m)                                                          3115
      goto 18741                                                           3116
18751 if(m .le. 2)goto 18761                                               3116
      alm=alm*alf                                                          3116
      goto 18741                                                           3117
18761 if(m .ne. 1)goto 18771                                               3117
      alm=big                                                              3117
      goto 18781                                                           3118
18771 continue                                                             3118
      alm0=0.0                                                             3119
18790 do 18791 j=1,ni                                                      3119
      if(ju(j).eq.0)goto 18791                                             3120
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           3121
18791 continue                                                             3122
18792 continue                                                             3122
      alm0=alm0/max(bta,1.0e-3)                                            3122
      alm=alf*alm0                                                         3123
18781 continue                                                             3124
18741 continue                                                             3124
      dem=alm*omb                                                          3124
      ab=alm*bta                                                           3124
      rsq0=rsq                                                             3124
      jz=1                                                                 3125
      tlam=bta*(2.0*alm-alm0)                                              3126
18800 do 18801 k=1,ni                                                      3126
      if(ix(k).eq.1)goto 18801                                             3126
      if(ju(k).eq.0)goto 18801                                             3127
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       3128
18801 continue                                                             3129
18802 continue                                                             3129
18810 continue                                                             3129
18811 continue                                                             3129
      if(iz*jz.ne.0) go to 10290                                           3130
10740 continue                                                             3130
      nlp=nlp+1                                                            3130
      dlx=0.0                                                              3131
18820 do 18821 k=1,ni                                                      3131
      if(ix(k).eq.0)goto 18821                                             3131
      gkn=0.0                                                              3132
18830 do 18831 j=1,nr                                                      3132
      gj(j)=dot_product(y(:,j),x(:,k))                                     3133
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3133
      gkn=gkn+gk(j)**2                                                     3135
18831 continue                                                             3135
18832 continue                                                             3135
      gkn=sqrt(gkn)                                                        3135
      u=1.0-ab*vp(k)/gkn                                                   3135
      del=a(:,k)                                                           3136
      if(u .gt. 0.0)goto 18851                                             3136
      a(:,k)=0.0                                                           3136
      goto 18861                                                           3137
18851 continue                                                             3137
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3138
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3140 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3141
18861 continue                                                             3142
18841 continue                                                             3142
      del=a(:,k)-del                                                       3142
      if(maxval(abs(del)).le.0.0)goto 18821                                3143
18870 do 18871 j=1,nr                                                      3143
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3144
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3144
      dlx=max(dlx,xv(k)*del(j)**2)                                         3145
18871 continue                                                             3146
18872 continue                                                             3146
      if(mm(k) .ne. 0)goto 18891                                           3146
      nin=nin+1                                                            3146
      if(nin.gt.nx)goto 18822                                              3147
      mm(k)=nin                                                            3147
      ia(nin)=k                                                            3148
18891 continue                                                             3149
18821 continue                                                             3150
18822 continue                                                             3150
      if(nin.gt.nx)goto 18812                                              3151
      if(dlx .ge. thr)goto 18911                                           3151
      ixx=0                                                                3152
18920 do 18921 k=1,ni                                                      3152
      if(ix(k).eq.1)goto 18921                                             3152
      if(ju(k).eq.0)goto 18921                                             3152
      g(k)=0.0                                                             3153
18930 do 18931 j=1,nr                                                      3153
      g(k)=g(k)+dot_product(y(:,j),x(:,k))**2                              3153
18931 continue                                                             3154
18932 continue                                                             3154
      g(k)=sqrt(g(k))                                                      3155
      if(g(k) .le. ab*vp(k))goto 18951                                     3155
      ix(k)=1                                                              3155
      ixx=1                                                                3155
18951 continue                                                             3156
18921 continue                                                             3157
18922 continue                                                             3157
      if(ixx.eq.1) go to 10740                                             3158
      goto 18812                                                           3159
18911 continue                                                             3160
      if(nlp .le. maxit)goto 18971                                         3160
      jerr=-m                                                              3160
      return                                                               3160
18971 continue                                                             3161
10290 continue                                                             3161
      iz=1                                                                 3162
18980 continue                                                             3162
18981 continue                                                             3162
      nlp=nlp+1                                                            3162
      dlx=0.0                                                              3163
18990 do 18991 l=1,nin                                                     3163
      k=ia(l)                                                              3163
      gkn=0.0                                                              3164
19000 do 19001 j=1,nr                                                      3164
      gj(j)=dot_product(y(:,j),x(:,k))                                     3165
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3165
      gkn=gkn+gk(j)**2                                                     3167
19001 continue                                                             3167
19002 continue                                                             3167
      gkn=sqrt(gkn)                                                        3167
      u=1.0-ab*vp(k)/gkn                                                   3167
      del=a(:,k)                                                           3168
      if(u .gt. 0.0)goto 19021                                             3168
      a(:,k)=0.0                                                           3168
      goto 19031                                                           3169
19021 continue                                                             3169
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3170
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3172 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3173
19031 continue                                                             3174
19011 continue                                                             3174
      del=a(:,k)-del                                                       3174
      if(maxval(abs(del)).le.0.0)goto 18991                                3175
19040 do 19041 j=1,nr                                                      3175
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3176
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3176
      dlx=max(dlx,xv(k)*del(j)**2)                                         3177
19041 continue                                                             3178
19042 continue                                                             3178
18991 continue                                                             3179
18992 continue                                                             3179
      if(dlx.lt.thr)goto 18982                                             3179
      if(nlp .le. maxit)goto 19061                                         3179
      jerr=-m                                                              3179
      return                                                               3179
19061 continue                                                             3180
      goto 18981                                                           3181
18982 continue                                                             3181
      jz=0                                                                 3182
      goto 18811                                                           3183
18812 continue                                                             3183
      if(nin .le. nx)goto 19081                                            3183
      jerr=-10000-m                                                        3183
      goto 18732                                                           3183
19081 continue                                                             3184
      if(nin .le. 0)goto 19101                                             3184
19110 do 19111 j=1,nr                                                      3184
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3184
19111 continue                                                             3184
19112 continue                                                             3184
19101 continue                                                             3185
      kin(m)=nin                                                           3186
      rsqo(m)=1.0-rsq/ys0                                                  3186
      almo(m)=alm                                                          3186
      lmu=m                                                                3187
      if(m.lt.mnl)goto 18731                                               3187
      if(flmin.ge.1.0)goto 18731                                           3188
      me=0                                                                 3188
19120 do 19121 j=1,nin                                                     3188
      if(ao(j,1,m).ne.0.0) me=me+1                                         3188
19121 continue                                                             3188
19122 continue                                                             3188
      if(me.gt.ne)goto 18732                                               3189
      if(rsq0-rsq.lt.sml*rsq)goto 18732                                    3189
      if(rsqo(m).gt.rsqmax)goto 18732                                      3190
18731 continue                                                             3191
18732 continue                                                             3191
      deallocate(a,mm,g,ix,del,gj,gk)                                      3192
      return                                                               3193
      end                                                                  3194
      subroutine chkbnds(nr,gk,gkn,xv,cl,al1,al2,a,isc,jerr)               3195
      real gk(nr),cl(2,nr),a(nr)                                           3195
      integer isc(nr)                                                      3196
      kerr=0                                                               3196
      al1p=1.0+al1/xv                                                      3196
      al2p=al2/xv                                                          3196
      isc=0                                                                3197
      gsq=gkn**2                                                           3197
      asq=dot_product(a,a)                                                 3197
      usq=0.0                                                              3198
19130 continue                                                             3198
19131 continue                                                             3198
      vmx=0.0                                                              3199
19140 do 19141 k=1,nr                                                      3199
      v=max(a(k)-cl(2,k),cl(1,k)-a(k))                                     3200
      if(v .le. vmx)goto 19161                                             3200
      vmx=v                                                                3200
      kn=k                                                                 3200
19161 continue                                                             3201
19141 continue                                                             3202
19142 continue                                                             3202
      if(vmx.le.0.0)goto 19132                                             3202
      if(isc(kn).ne.0)goto 19132                                           3203
      gsq=gsq-gk(kn)**2                                                    3203
      g=sqrt(gsq)/xv                                                       3204
      if(a(kn).lt.cl(1,kn)) u=cl(1,kn)                                     3204
      if(a(kn).gt.cl(2,kn)) u=cl(2,kn)                                     3205
      usq=usq+u**2                                                         3206
      if(usq .ne. 0.0)goto 19181                                           3206
      b=max(0.0,(g-al2p)/al1p)                                             3206
      goto 19191                                                           3207
19181 continue                                                             3207
      b0=sqrt(asq-a(kn)**2)                                                3208
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                     3208
      if(kerr.ne.0)goto 19132                                              3209
19191 continue                                                             3210
19171 continue                                                             3210
      asq=usq+b**2                                                         3210
      if(asq .gt. 0.0)goto 19211                                           3210
      a=0.0                                                                3210
      goto 19132                                                           3210
19211 continue                                                             3211
      a(kn)=u                                                              3211
      isc(kn)=1                                                            3211
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                     3212
19220 do 19221 j=1,nr                                                      3212
      if(isc(j).eq.0) a(j)=f*gk(j)                                         3212
19221 continue                                                             3213
19222 continue                                                             3213
      goto 19131                                                           3214
19132 continue                                                             3214
      if(kerr.ne.0) jerr=kerr                                              3215
      return                                                               3216
      end                                                                  3217
      subroutine chkbnds1(nr,gk,gkn,xv,cl1,cl2,al1,al2,a,isc,jerr)         3218
      real gk(nr),a(nr)                                                    3218
      integer isc(nr)                                                      3219
      kerr=0                                                               3219
      al1p=1.0+al1/xv                                                      3219
      al2p=al2/xv                                                          3219
      isc=0                                                                3220
      gsq=gkn**2                                                           3220
      asq=dot_product(a,a)                                                 3220
      usq=0.0                                                              3221
19230 continue                                                             3221
19231 continue                                                             3221
      vmx=0.0                                                              3222
19240 do 19241 k=1,nr                                                      3222
      v=max(a(k)-cl2,cl1-a(k))                                             3223
      if(v .le. vmx)goto 19261                                             3223
      vmx=v                                                                3223
      kn=k                                                                 3223
19261 continue                                                             3224
19241 continue                                                             3225
19242 continue                                                             3225
      if(vmx.le.0.0)goto 19232                                             3225
      if(isc(kn).ne.0)goto 19232                                           3226
      gsq=gsq-gk(kn)**2                                                    3226
      g=sqrt(gsq)/xv                                                       3227
      if(a(kn).lt.cl1) u=cl1                                               3227
      if(a(kn).gt.cl2) u=cl2                                               3228
      usq=usq+u**2                                                         3229
      if(usq .ne. 0.0)goto 19281                                           3229
      b=max(0.0,(g-al2p)/al1p)                                             3229
      goto 19291                                                           3230
19281 continue                                                             3230
      b0=sqrt(asq-a(kn)**2)                                                3231
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                     3231
      if(kerr.ne.0)goto 19232                                              3232
19291 continue                                                             3233
19271 continue                                                             3233
      asq=usq+b**2                                                         3233
      if(asq .gt. 0.0)goto 19311                                           3233
      a=0.0                                                                3233
      goto 19232                                                           3233
19311 continue                                                             3234
      a(kn)=u                                                              3234
      isc(kn)=1                                                            3234
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                     3235
19320 do 19321 j=1,nr                                                      3235
      if(isc(j).eq.0) a(j)=f*gk(j)                                         3235
19321 continue                                                             3236
19322 continue                                                             3236
      goto 19231                                                           3237
19232 continue                                                             3237
      if(kerr.ne.0) jerr=kerr                                              3238
      return                                                               3239
      end                                                                  3240
      function bnorm(b0,al1p,al2p,g,usq,jerr)                              3241
      data thr,mxit /1.0e-10,100/                                          3242
      b=b0                                                                 3242
      zsq=b**2+usq                                                         3242
      if(zsq .gt. 0.0)goto 19341                                           3242
      bnorm=0.0                                                            3242
      return                                                               3242
19341 continue                                                             3243
      z=sqrt(zsq)                                                          3243
      f=b*(al1p+al2p/z)-g                                                  3243
      jerr=0                                                               3244
19350 do 19351 it=1,mxit                                                   3244
      b=b-f/(al1p+al2p*usq/(z*zsq))                                        3245
      zsq=b**2+usq                                                         3245
      if(zsq .gt. 0.0)goto 19371                                           3245
      bnorm=0.0                                                            3245
      return                                                               3245
19371 continue                                                             3246
      z=sqrt(zsq)                                                          3246
      f=b*(al1p+al2p/z)-g                                                  3247
      if(abs(f).le.thr)goto 19352                                          3247
      if(b .gt. 0.0)goto 19391                                             3247
      b=0.0                                                                3247
      goto 19352                                                           3247
19391 continue                                                             3248
19351 continue                                                             3249
19352 continue                                                             3249
      bnorm=b                                                              3249
      if(it.ge.mxit) jerr=90000                                            3250
      return                                                               3251
      entry chg_bnorm(arg,irg)                                             3251
      thr=arg                                                              3251
      mxit=irg                                                             3251
      return                                                               3252
      entry get_bnorm(arg,irg)                                             3252
      arg=thr                                                              3252
      irg=mxit                                                             3252
      return                                                               3253
      end                                                                  3254
      subroutine multsolns(ni,nx,nr,lmu,a,ia,nin,b)                        3255
      real a(nx,nr,lmu),b(ni,nr,lmu)                                       3255
      integer ia(nx),nin(lmu)                                              3256
19400 do 19401 lam=1,lmu                                                   3256
      call multuncomp(ni,nr,nx,a(1,1,lam),ia,nin(lam),b(1,1,lam))          3256
19401 continue                                                             3257
19402 continue                                                             3257
      return                                                               3258
      end                                                                  3259
      subroutine multuncomp(ni,nr,nx,ca,ia,nin,a)                          3260
      real ca(nx,nr),a(ni,nr)                                              3260
      integer ia(nx)                                                       3261
      a=0.0                                                                3262
      if(nin .le. 0)goto 19421                                             3262
19430 do 19431 j=1,nr                                                      3262
      a(ia(1:nin),j)=ca(1:nin,j)                                           3262
19431 continue                                                             3262
19432 continue                                                             3262
19421 continue                                                             3263
      return                                                               3264
      end                                                                  3265
      subroutine multmodval(nx,nr,a0,ca,ia,nin,n,x,f)                      3266
      real a0(nr),ca(nx,nr),x(n,*),f(nr,n)                                 3266
      integer ia(nx)                                                       3267
19440 do 19441 i=1,n                                                       3267
      f(:,i)=a0                                                            3267
19441 continue                                                             3267
19442 continue                                                             3267
      if(nin.le.0) return                                                  3268
19450 do 19451 i=1,n                                                       3268
19460 do 19461 j=1,nr                                                      3268
      f(j,i)=f(j,i)+dot_product(ca(1:nin,j),x(i,ia(1:nin)))                3268
19461 continue                                                             3268
19462 continue                                                             3268
19451 continue                                                             3269
19452 continue                                                             3269
      return                                                               3270
      end                                                                  3271
      subroutine multspelnet  (parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,   3274 
     *nlam,flmin,ulam,thr,isd,jsd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,j
     *err)
      real x(*),y(no,nr),w(no),vp(ni),ulam(nlam),cl(2,ni)                  3275
      real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                  3276
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3277
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 19481                                    3280
      jerr=10000                                                           3280
      return                                                               3280
19481 continue                                                             3281
      allocate(vq(1:ni),stat=jerr)                                         3281
      if(jerr.ne.0) return                                                 3282
      vq=max(0.0,vp)                                                       3282
      vq=vq*ni/sum(vq)                                                     3283
      call multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,fl   3285 
     *min,  ulam,thr,isd,jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      deallocate(vq)                                                       3286
      return                                                               3287
      end                                                                  3288
      subroutine multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,n   3290 
     *lam,flmin,  ulam,thr,isd,jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,je
     *rr)
      real x(*),vp(ni),y(no,nr),w(no),ulam(nlam),cl(2,ni)                  3291
      real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                  3292
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3293
      real, dimension (:), allocatable :: xm,xs,xv,ym,ys                        
      integer, dimension (:), allocatable :: ju                                 
      real, dimension (:,:,:), allocatable :: clt                               
      allocate(clt(1:2,1:nr,1:ni),stat=jerr)                                    
      allocate(xm(1:ni),stat=ierr)                                         3299
      jerr=jerr+ierr                                                       3300
      allocate(xs(1:ni),stat=ierr)                                         3300
      jerr=jerr+ierr                                                       3301
      allocate(ym(1:nr),stat=ierr)                                         3301
      jerr=jerr+ierr                                                       3302
      allocate(ys(1:nr),stat=ierr)                                         3302
      jerr=jerr+ierr                                                       3303
      allocate(ju(1:ni),stat=ierr)                                         3303
      jerr=jerr+ierr                                                       3304
      allocate(xv(1:ni),stat=ierr)                                         3304
      jerr=jerr+ierr                                                       3305
      if(jerr.ne.0) return                                                 3306
      call spchkvars(no,ni,x,ix,ju)                                        3307
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3308
      if(maxval(ju) .gt. 0)goto 19501                                      3308
      jerr=7777                                                            3308
      return                                                               3308
19501 continue                                                             3309
      call multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,xm,xs,ym,ys,x   3310 
     *v,ys0,jerr)
      if(jerr.ne.0) return                                                 3311
19510 do 19511 j=1,ni                                                      3311
19520 do 19521 k=1,nr                                                      3311
19530 do 19531 i=1,2                                                       3311
      clt(i,k,j)=cl(i,j)                                                   3311
19531 continue                                                             3311
19532 continue                                                             3311
19521 continue                                                             3311
19522 continue                                                             3311
19511 continue                                                             3312
19512 continue                                                             3312
      if(isd .le. 0)goto 19551                                             3312
19560 do 19561 j=1,ni                                                      3312
19570 do 19571 k=1,nr                                                      3312
19580 do 19581 i=1,2                                                       3312
      clt(i,k,j)=clt(i,k,j)*xs(j)                                          3312
19581 continue                                                             3312
19582 continue                                                             3312
19571 continue                                                             3312
19572 continue                                                             3312
19561 continue                                                             3312
19562 continue                                                             3312
19551 continue                                                             3313
      if(jsd .le. 0)goto 19601                                             3313
19610 do 19611 j=1,ni                                                      3313
19620 do 19621 k=1,nr                                                      3313
19630 do 19631 i=1,2                                                       3313
      clt(i,k,j)=clt(i,k,j)/ys(k)                                          3313
19631 continue                                                             3313
19632 continue                                                             3313
19621 continue                                                             3313
19622 continue                                                             3313
19611 continue                                                             3313
19612 continue                                                             3313
19601 continue                                                             3314
      call multspelnet2(parm,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,clt,nlam,f   3316 
     *lmin,  ulam,thr,maxit,xm,xs,xv,ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 3317
19640 do 19641 k=1,lmu                                                     3317
      nk=nin(k)                                                            3318
19650 do 19651 j=1,nr                                                      3319
19660 do 19661 l=1,nk                                                      3319
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  3319
19661 continue                                                             3320
19662 continue                                                             3320
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 3321
19651 continue                                                             3322
19652 continue                                                             3322
19641 continue                                                             3323
19642 continue                                                             3323
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                    3324
      return                                                               3325
      end                                                                  3326
      subroutine multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,  xm,xs   3328 
     *,ym,ys,xv,ys0,jerr)
      real x(*),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),ys(nr)          3329
      integer ix(*),jx(*),ju(ni)                                           3330
      w=w/sum(w)                                                           3331
19670 do 19671 j=1,ni                                                      3331
      if(ju(j).eq.0)goto 19671                                             3332
      jb=ix(j)                                                             3332
      je=ix(j+1)-1                                                         3332
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3333
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 3334
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       3335
19671 continue                                                             3336
19672 continue                                                             3336
      if(isd .ne. 0)goto 19691                                             3336
      xs=1.0                                                               3336
      goto 19701                                                           3336
19691 continue                                                             3336
      xv=1.0                                                               3336
19701 continue                                                             3337
19681 continue                                                             3337
      ys0=0.0                                                              3338
19710 do 19711 j=1,nr                                                      3339
      ym(j)=dot_product(w,y(:,j))                                          3339
      y(:,j)=y(:,j)-ym(j)                                                  3340
      z=dot_product(w,y(:,j)**2)                                           3341
      if(jsd .le. 0)goto 19731                                             3341
      ys(j)=sqrt(z)                                                        3341
      y(:,j)=y(:,j)/ys(j)                                                  3341
      goto 19741                                                           3342
19731 continue                                                             3342
      ys0=ys0+z                                                            3342
19741 continue                                                             3343
19721 continue                                                             3343
19711 continue                                                             3344
19712 continue                                                             3344
      if(jsd .ne. 0)goto 19761                                             3344
      ys=1.0                                                               3344
      goto 19771                                                           3344
19761 continue                                                             3344
      ys0=nr                                                               3344
19771 continue                                                             3345
19751 continue                                                             3345
      return                                                               3346
      end                                                                  3347
      subroutine multspelnet2(beta,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,n   3349 
     *lam,flmin,  ulam,thri,maxit,xm,xs,xv,ys0,lmu,ao,ia,kin,rsqo,almo,n
     *lp,jerr)
      real y(no,nr),w(no),x(*),vp(ni),ulam(nlam),cl(2,nr,ni)               3350
      real ao(nx,nr,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)       3351
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          3352
      real, dimension (:), allocatable :: g,gj,gk,del,o                         
      integer, dimension (:), allocatable :: mm,iy,isc                          
      real, dimension (:,:), allocatable :: a                                   
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               3359
      allocate(mm(1:ni),stat=ierr)                                         3359
      jerr=jerr+ierr                                                       3360
      allocate(g(1:ni),stat=ierr)                                          3360
      jerr=jerr+ierr                                                       3361
      allocate(gj(1:nr),stat=ierr)                                         3361
      jerr=jerr+ierr                                                       3362
      allocate(gk(1:nr),stat=ierr)                                         3362
      jerr=jerr+ierr                                                       3363
      allocate(del(1:nr),stat=ierr)                                        3363
      jerr=jerr+ierr                                                       3364
      allocate(o(1:nr),stat=ierr)                                          3364
      jerr=jerr+ierr                                                       3365
      allocate(iy(1:ni),stat=ierr)                                         3365
      jerr=jerr+ierr                                                       3366
      allocate(isc(1:nr),stat=ierr)                                        3366
      jerr=jerr+ierr                                                       3367
      if(jerr.ne.0) return                                                 3368
      bta=beta                                                             3368
      omb=1.0-bta                                                          3368
      alm=0.0                                                              3368
      iy=0                                                                 3368
      thr=thri*ys0/nr                                                      3369
      if(flmin .ge. 1.0)goto 19791                                         3369
      eqs=max(eps,flmin)                                                   3369
      alf=eqs**(1.0/(nlam-1))                                              3369
19791 continue                                                             3370
      rsq=ys0                                                              3370
      a=0.0                                                                3370
      mm=0                                                                 3370
      o=0.0                                                                3370
      nlp=0                                                                3370
      nin=nlp                                                              3370
      iz=0                                                                 3370
      mnl=min(mnlam,nlam)                                                  3371
19800 do 19801 j=1,ni                                                      3371
      if(ju(j).eq.0)goto 19801                                             3371
      jb=ix(j)                                                             3371
      je=ix(j+1)-1                                                         3371
      g(j)=0.0                                                             3372
19810 do 19811 k=1,nr                                                      3373
      g(j)=g(j)+(dot_product(y(jx(jb:je),k),w(jx(jb:je))*x(jb:je))/xs(j)   3374 
     *)**2
19811 continue                                                             3375
19812 continue                                                             3375
      g(j)=sqrt(g(j))                                                      3376
19801 continue                                                             3377
19802 continue                                                             3377
19820 do 19821 m=1,nlam                                                    3377
      alm0=alm                                                             3378
      if(flmin .lt. 1.0)goto 19841                                         3378
      alm=ulam(m)                                                          3378
      goto 19831                                                           3379
19841 if(m .le. 2)goto 19851                                               3379
      alm=alm*alf                                                          3379
      goto 19831                                                           3380
19851 if(m .ne. 1)goto 19861                                               3380
      alm=big                                                              3380
      goto 19871                                                           3381
19861 continue                                                             3381
      alm0=0.0                                                             3382
19880 do 19881 j=1,ni                                                      3382
      if(ju(j).eq.0)goto 19881                                             3383
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           3384
19881 continue                                                             3385
19882 continue                                                             3385
      alm0=alm0/max(bta,1.0e-3)                                            3385
      alm=alf*alm0                                                         3386
19871 continue                                                             3387
19831 continue                                                             3387
      dem=alm*omb                                                          3387
      ab=alm*bta                                                           3387
      rsq0=rsq                                                             3387
      jz=1                                                                 3388
      tlam=bta*(2.0*alm-alm0)                                              3389
19890 do 19891 k=1,ni                                                      3389
      if(iy(k).eq.1)goto 19891                                             3389
      if(ju(k).eq.0)goto 19891                                             3390
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       3391
19891 continue                                                             3392
19892 continue                                                             3392
19900 continue                                                             3392
19901 continue                                                             3392
      if(iz*jz.ne.0) go to 10290                                           3393
10740 continue                                                             3393
      nlp=nlp+1                                                            3393
      dlx=0.0                                                              3394
19910 do 19911 k=1,ni                                                      3394
      if(iy(k).eq.0)goto 19911                                             3394
      jb=ix(k)                                                             3394
      je=ix(k+1)-1                                                         3394
      gkn=0.0                                                              3395
19920 do 19921 j=1,nr                                                      3396
      gj(j)=dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(k)   3397
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3397
      gkn=gkn+gk(j)**2                                                     3398
19921 continue                                                             3399
19922 continue                                                             3399
      gkn=sqrt(gkn)                                                        3399
      u=1.0-ab*vp(k)/gkn                                                   3399
      del=a(:,k)                                                           3400
      if(u .gt. 0.0)goto 19941                                             3400
      a(:,k)=0.0                                                           3400
      goto 19951                                                           3401
19941 continue                                                             3401
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3402
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3404 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3405
19951 continue                                                             3406
19931 continue                                                             3406
      del=a(:,k)-del                                                       3406
      if(maxval(abs(del)).le.0.0)goto 19911                                3407
      if(mm(k) .ne. 0)goto 19971                                           3407
      nin=nin+1                                                            3407
      if(nin.gt.nx)goto 19912                                              3408
      mm(k)=nin                                                            3408
      ia(nin)=k                                                            3409
19971 continue                                                             3410
19980 do 19981 j=1,nr                                                      3410
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3411
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3412
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3412
      dlx=max(xv(k)*del(j)**2,dlx)                                         3413
19981 continue                                                             3414
19982 continue                                                             3414
19911 continue                                                             3415
19912 continue                                                             3415
      if(nin.gt.nx)goto 19902                                              3416
      if(dlx .ge. thr)goto 20001                                           3416
      ixx=0                                                                3417
20010 do 20011 j=1,ni                                                      3417
      if(iy(j).eq.1)goto 20011                                             3417
      if(ju(j).eq.0)goto 20011                                             3418
      jb=ix(j)                                                             3418
      je=ix(j+1)-1                                                         3418
      g(j)=0.0                                                             3419
20020 do 20021 k=1,nr                                                      3419
      g(j)=g(j)+  (dot_product(y(jx(jb:je),k)+o(k),w(jx(jb:je))*x(jb:je)   3421 
     *)/xs(j))**2
20021 continue                                                             3422
20022 continue                                                             3422
      g(j)=sqrt(g(j))                                                      3423
      if(g(j) .le. ab*vp(j))goto 20041                                     3423
      iy(j)=1                                                              3423
      ixx=1                                                                3423
20041 continue                                                             3424
20011 continue                                                             3425
20012 continue                                                             3425
      if(ixx.eq.1) go to 10740                                             3426
      goto 19902                                                           3427
20001 continue                                                             3428
      if(nlp .le. maxit)goto 20061                                         3428
      jerr=-m                                                              3428
      return                                                               3428
20061 continue                                                             3429
10290 continue                                                             3429
      iz=1                                                                 3430
20070 continue                                                             3430
20071 continue                                                             3430
      nlp=nlp+1                                                            3430
      dlx=0.0                                                              3431
20080 do 20081 l=1,nin                                                     3431
      k=ia(l)                                                              3431
      jb=ix(k)                                                             3431
      je=ix(k+1)-1                                                         3431
      gkn=0.0                                                              3432
20090 do 20091 j=1,nr                                                      3432
      gj(j)=  dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(   3434 
     *k)
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3434
      gkn=gkn+gk(j)**2                                                     3435
20091 continue                                                             3436
20092 continue                                                             3436
      gkn=sqrt(gkn)                                                        3436
      u=1.0-ab*vp(k)/gkn                                                   3436
      del=a(:,k)                                                           3437
      if(u .gt. 0.0)goto 20111                                             3437
      a(:,k)=0.0                                                           3437
      goto 20121                                                           3438
20111 continue                                                             3438
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3439
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3441 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3442
20121 continue                                                             3443
20101 continue                                                             3443
      del=a(:,k)-del                                                       3443
      if(maxval(abs(del)).le.0.0)goto 20081                                3444
20130 do 20131 j=1,nr                                                      3444
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3445
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3446
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3446
      dlx=max(xv(k)*del(j)**2,dlx)                                         3447
20131 continue                                                             3448
20132 continue                                                             3448
20081 continue                                                             3449
20082 continue                                                             3449
      if(dlx.lt.thr)goto 20072                                             3449
      if(nlp .le. maxit)goto 20151                                         3449
      jerr=-m                                                              3449
      return                                                               3449
20151 continue                                                             3450
      goto 20071                                                           3451
20072 continue                                                             3451
      jz=0                                                                 3452
      goto 19901                                                           3453
19902 continue                                                             3453
      if(nin .le. nx)goto 20171                                            3453
      jerr=-10000-m                                                        3453
      goto 19822                                                           3453
20171 continue                                                             3454
      if(nin .le. 0)goto 20191                                             3454
20200 do 20201 j=1,nr                                                      3454
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3454
20201 continue                                                             3454
20202 continue                                                             3454
20191 continue                                                             3455
      kin(m)=nin                                                           3456
      rsqo(m)=1.0-rsq/ys0                                                  3456
      almo(m)=alm                                                          3456
      lmu=m                                                                3457
      if(m.lt.mnl)goto 19821                                               3457
      if(flmin.ge.1.0)goto 19821                                           3458
      me=0                                                                 3458
20210 do 20211 j=1,nin                                                     3458
      if(ao(j,1,m).ne.0.0) me=me+1                                         3458
20211 continue                                                             3458
20212 continue                                                             3458
      if(me.gt.ne)goto 19822                                               3459
      if(rsq0-rsq.lt.sml*rsq)goto 19822                                    3459
      if(rsqo(m).gt.rsqmax)goto 19822                                      3460
19821 continue                                                             3461
19822 continue                                                             3461
      deallocate(a,mm,g,iy,gj,gk,del,o)                                    3462
      return                                                               3463
      end                                                                  3464
      subroutine multlognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,f   3466 
     *lmin,ulam,  shri,isd,maxit,xv,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr
     *)
      real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),cl(2,ni)     3467
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),xv(ni)            3468
      integer ju(ni),m(nx),kin(nlam)                                       3469
      real, dimension (:,:), allocatable :: q,r,b,bs                            
      real, dimension (:), allocatable :: sxp,sxpl,ga,gk,del                    
      integer, dimension (:), allocatable :: mm,is,ixx,isc                      
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no,1:nc),stat=ierr); jerr=jerr+ierr;                         
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               3478
      exmn=-exmx                                                           3479
      allocate(mm(1:ni),stat=ierr)                                         3479
      jerr=jerr+ierr                                                       3480
      allocate(is(1:max(nc,ni)),stat=ierr)                                 3480
      jerr=jerr+ierr                                                       3481
      allocate(sxp(1:no),stat=ierr)                                        3481
      jerr=jerr+ierr                                                       3482
      allocate(sxpl(1:no),stat=ierr)                                       3482
      jerr=jerr+ierr                                                       3483
      allocate(ga(1:ni),stat=ierr)                                         3483
      jerr=jerr+ierr                                                       3484
      allocate(ixx(1:ni),stat=ierr)                                        3484
      jerr=jerr+ierr                                                       3485
      allocate(gk(1:nc),stat=ierr)                                         3485
      jerr=jerr+ierr                                                       3486
      allocate(del(1:nc),stat=ierr)                                        3486
      jerr=jerr+ierr                                                       3487
      allocate(isc(1:nc),stat=ierr)                                        3487
      jerr=jerr+ierr                                                       3488
      if(jerr.ne.0) return                                                 3489
      pmax=1.0-pmin                                                        3489
      emin=pmin/pmax                                                       3489
      emax=1.0/emin                                                        3490
      bta=parm                                                             3490
      omb=1.0-bta                                                          3490
      dev1=0.0                                                             3490
      dev0=0.0                                                             3491
20220 do 20221 ic=1,nc                                                     3491
      q0=dot_product(w,y(:,ic))                                            3492
      if(q0 .gt. pmin)goto 20241                                           3492
      jerr =8000+ic                                                        3492
      return                                                               3492
20241 continue                                                             3493
      if(q0 .lt. pmax)goto 20261                                           3493
      jerr =9000+ic                                                        3493
      return                                                               3493
20261 continue                                                             3494
      b(0,ic)=log(q0)                                                      3494
      dev1=dev1-q0*b(0,ic)                                                 3494
      b(1:ni,ic)=0.0                                                       3495
20221 continue                                                             3496
20222 continue                                                             3496
      ixx=0                                                                3496
      al=0.0                                                               3497
      if(nonzero(no*nc,g) .ne. 0)goto 20281                                3498
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         3498
      sxp=0.0                                                              3499
20290 do 20291 ic=1,nc                                                     3499
      q(:,ic)=exp(b(0,ic))                                                 3499
      sxp=sxp+q(:,ic)                                                      3499
20291 continue                                                             3500
20292 continue                                                             3500
      goto 20301                                                           3501
20281 continue                                                             3501
20310 do 20311 i=1,no                                                      3501
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3501
20311 continue                                                             3501
20312 continue                                                             3501
      sxp=0.0                                                              3502
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 3502
      if(jerr.ne.0) return                                                 3503
      dev1=0.0                                                             3504
20320 do 20321 ic=1,nc                                                     3504
      q(:,ic)=b(0,ic)+g(:,ic)                                              3505
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             3506
      q(:,ic)=exp(q(:,ic))                                                 3506
      sxp=sxp+q(:,ic)                                                      3507
20321 continue                                                             3508
20322 continue                                                             3508
      sxpl=w*log(sxp)                                                      3508
20330 do 20331 ic=1,nc                                                     3508
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  3508
20331 continue                                                             3509
20332 continue                                                             3509
20301 continue                                                             3510
20271 continue                                                             3510
20340 do 20341 ic=1,nc                                                     3510
20350 do 20351 i=1,no                                                      3510
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               3510
20351 continue                                                             3510
20352 continue                                                             3510
20341 continue                                                             3511
20342 continue                                                             3511
      dev0=dev0+dev1                                                       3512
      if(flmin .ge. 1.0)goto 20371                                         3512
      eqs=max(eps,flmin)                                                   3512
      alf=eqs**(1.0/(nlam-1))                                              3512
20371 continue                                                             3513
      m=0                                                                  3513
      mm=0                                                                 3513
      nin=0                                                                3513
      nlp=0                                                                3513
      mnl=min(mnlam,nlam)                                                  3513
      bs=0.0                                                               3513
      shr=shri*dev0                                                        3514
      ga=0.0                                                               3515
20380 do 20381 ic=1,nc                                                     3515
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3516
20390 do 20391 j=1,ni                                                      3516
      if(ju(j).ne.0) ga(j)=ga(j)+dot_product(r(:,ic),x(:,j))**2            3516
20391 continue                                                             3517
20392 continue                                                             3517
20381 continue                                                             3518
20382 continue                                                             3518
      ga=sqrt(ga)                                                          3519
20400 do 20401 ilm=1,nlam                                                  3519
      al0=al                                                               3520
      if(flmin .lt. 1.0)goto 20421                                         3520
      al=ulam(ilm)                                                         3520
      goto 20411                                                           3521
20421 if(ilm .le. 2)goto 20431                                             3521
      al=al*alf                                                            3521
      goto 20411                                                           3522
20431 if(ilm .ne. 1)goto 20441                                             3522
      al=big                                                               3522
      goto 20451                                                           3523
20441 continue                                                             3523
      al0=0.0                                                              3524
20460 do 20461 j=1,ni                                                      3524
      if(ju(j).eq.0)goto 20461                                             3524
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3524
20461 continue                                                             3525
20462 continue                                                             3525
      al0=al0/max(bta,1.0e-3)                                              3525
      al=alf*al0                                                           3526
20451 continue                                                             3527
20411 continue                                                             3527
      al2=al*omb                                                           3527
      al1=al*bta                                                           3527
      tlam=bta*(2.0*al-al0)                                                3528
20470 do 20471 k=1,ni                                                      3528
      if(ixx(k).eq.1)goto 20471                                            3528
      if(ju(k).eq.0)goto 20471                                             3529
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     3530
20471 continue                                                             3531
20472 continue                                                             3531
10740 continue                                                             3532
20480 continue                                                             3532
20481 continue                                                             3532
      ix=0                                                                 3532
      jx=ix                                                                3532
      kx=jx                                                                3532
      t=0.0                                                                3533
20490 do 20491 ic=1,nc                                                     3533
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       3533
20491 continue                                                             3534
20492 continue                                                             3534
      if(t .ge. eps)goto 20511                                             3534
      kx=1                                                                 3534
      goto 20482                                                           3534
20511 continue                                                             3534
      t=2.0*t                                                              3534
      alt=al1/t                                                            3534
      al2t=al2/t                                                           3535
20520 do 20521 ic=1,nc                                                     3536
      bs(0,ic)=b(0,ic)                                                     3536
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          3537
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    3537
      d=sum(r(:,ic))                                                       3538
      b(0,ic)=b(0,ic)+d                                                    3538
      r(:,ic)=r(:,ic)-d*w                                                  3538
      dlx=max(dlx,d**2)                                                    3540
20521 continue                                                             3541
20522 continue                                                             3541
20530 continue                                                             3541
20531 continue                                                             3541
      nlp=nlp+nc                                                           3541
      dlx=0.0                                                              3542
20540 do 20541 k=1,ni                                                      3542
      if(ixx(k).eq.0)goto 20541                                            3542
      gkn=0.0                                                              3543
20550 do 20551 ic=1,nc                                                     3543
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     3544
      gkn=gkn+gk(ic)**2                                                    3545
20551 continue                                                             3546
20552 continue                                                             3546
      gkn=sqrt(gkn)                                                        3546
      u=1.0-alt*vp(k)/gkn                                                  3546
      del=b(k,:)                                                           3547
      if(u .gt. 0.0)goto 20571                                             3547
      b(k,:)=0.0                                                           3547
      goto 20581                                                           3548
20571 continue                                                             3548
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     3549
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(   3551 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 3552
20581 continue                                                             3553
20561 continue                                                             3553
      del=b(k,:)-del                                                       3553
      if(maxval(abs(del)).le.0.0)goto 20541                                3554
20590 do 20591 ic=1,nc                                                     3554
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3555
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     3556
20591 continue                                                             3557
20592 continue                                                             3557
      if(mm(k) .ne. 0)goto 20611                                           3557
      nin=nin+1                                                            3558
      if(nin .le. nx)goto 20631                                            3558
      jx=1                                                                 3558
      goto 20542                                                           3558
20631 continue                                                             3559
      mm(k)=nin                                                            3559
      m(nin)=k                                                             3560
20611 continue                                                             3561
20541 continue                                                             3562
20542 continue                                                             3562
      if(jx.gt.0)goto 20532                                                3562
      if(dlx.lt.shr)goto 20532                                             3563
      if(nlp .le. maxit)goto 20651                                         3563
      jerr=-ilm                                                            3563
      return                                                               3563
20651 continue                                                             3564
20660 continue                                                             3564
20661 continue                                                             3564
      nlp=nlp+nc                                                           3564
      dlx=0.0                                                              3565
20670 do 20671 l=1,nin                                                     3565
      k=m(l)                                                               3565
      gkn=0.0                                                              3566
20680 do 20681 ic=1,nc                                                     3566
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     3567
      gkn=gkn+gk(ic)**2                                                    3568
20681 continue                                                             3569
20682 continue                                                             3569
      gkn=sqrt(gkn)                                                        3569
      u=1.0-alt*vp(k)/gkn                                                  3569
      del=b(k,:)                                                           3570
      if(u .gt. 0.0)goto 20701                                             3570
      b(k,:)=0.0                                                           3570
      goto 20711                                                           3571
20701 continue                                                             3571
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     3572
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(   3574 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 3575
20711 continue                                                             3576
20691 continue                                                             3576
      del=b(k,:)-del                                                       3576
      if(maxval(abs(del)).le.0.0)goto 20671                                3577
20720 do 20721 ic=1,nc                                                     3577
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3578
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     3579
20721 continue                                                             3580
20722 continue                                                             3580
20671 continue                                                             3581
20672 continue                                                             3581
      if(dlx.lt.shr)goto 20662                                             3581
      if(nlp .le. maxit)goto 20741                                         3581
      jerr=-ilm                                                            3581
      return                                                               3581
20741 continue                                                             3583
      goto 20661                                                           3584
20662 continue                                                             3584
      goto 20531                                                           3585
20532 continue                                                             3585
      if(jx.gt.0)goto 20482                                                3586
20750 do 20751 ic=1,nc                                                     3587
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                                3588
      if(ix .ne. 0)goto 20771                                              3589
20780 do 20781 j=1,nin                                                     3589
      k=m(j)                                                               3590
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 20801                   3590
      ix=1                                                                 3590
      goto 20782                                                           3590
20801 continue                                                             3592
20781 continue                                                             3593
20782 continue                                                             3593
20771 continue                                                             3594
20810 do 20811 i=1,no                                                      3594
      fi=b(0,ic)+g(i,ic)                                                   3596
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         3597
      fi=min(max(exmn,fi),exmx)                                            3597
      sxp(i)=sxp(i)-q(i,ic)                                                3598
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    3599
      sxp(i)=sxp(i)+q(i,ic)                                                3600
20811 continue                                                             3601
20812 continue                                                             3601
20751 continue                                                             3602
20752 continue                                                             3602
      s=-sum(b(0,:))/nc                                                    3602
      b(0,:)=b(0,:)+s                                                      3603
      if(jx.gt.0)goto 20482                                                3604
      if(ix .ne. 0)goto 20831                                              3605
20840 do 20841 k=1,ni                                                      3605
      if(ixx(k).eq.1)goto 20841                                            3605
      if(ju(k).eq.0)goto 20841                                             3605
      ga(k)=0.0                                                            3605
20841 continue                                                             3606
20842 continue                                                             3606
20850 do 20851 ic=1,nc                                                     3606
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3607
20860 do 20861 k=1,ni                                                      3607
      if(ixx(k).eq.1)goto 20861                                            3607
      if(ju(k).eq.0)goto 20861                                             3608
      ga(k)=ga(k)+dot_product(r(:,ic),x(:,k))**2                           3609
20861 continue                                                             3610
20862 continue                                                             3610
20851 continue                                                             3611
20852 continue                                                             3611
      ga=sqrt(ga)                                                          3612
20870 do 20871 k=1,ni                                                      3612
      if(ixx(k).eq.1)goto 20871                                            3612
      if(ju(k).eq.0)goto 20871                                             3613
      if(ga(k) .le. al1*vp(k))goto 20891                                   3613
      ixx(k)=1                                                             3613
      ix=1                                                                 3613
20891 continue                                                             3614
20871 continue                                                             3615
20872 continue                                                             3615
      if(ix.eq.1) go to 10740                                              3616
      goto 20482                                                           3617
20831 continue                                                             3618
      goto 20481                                                           3619
20482 continue                                                             3619
      if(kx .le. 0)goto 20911                                              3619
      jerr=-20000-ilm                                                      3619
      goto 20402                                                           3619
20911 continue                                                             3620
      if(jx .le. 0)goto 20931                                              3620
      jerr=-10000-ilm                                                      3620
      goto 20402                                                           3620
20931 continue                                                             3620
      devi=0.0                                                             3621
20940 do 20941 ic=1,nc                                                     3622
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          3622
      a0(ic,ilm)=b(0,ic)                                                   3623
20950 do 20951 i=1,no                                                      3623
      if(y(i,ic).le.0.0)goto 20951                                         3624
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           3625
20951 continue                                                             3626
20952 continue                                                             3626
20941 continue                                                             3627
20942 continue                                                             3627
      kin(ilm)=nin                                                         3627
      alm(ilm)=al                                                          3627
      lmu=ilm                                                              3628
      dev(ilm)=(dev1-devi)/dev0                                            3629
      if(ilm.lt.mnl)goto 20401                                             3629
      if(flmin.ge.1.0)goto 20401                                           3630
      me=0                                                                 3630
20960 do 20961 j=1,nin                                                     3630
      if(a(j,1,ilm).ne.0.0) me=me+1                                        3630
20961 continue                                                             3630
20962 continue                                                             3630
      if(me.gt.ne)goto 20402                                               3631
      if(dev(ilm).gt.devmax)goto 20402                                     3631
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 20402                             3632
20401 continue                                                             3633
20402 continue                                                             3633
      g=log(q)                                                             3633
20970 do 20971 i=1,no                                                      3633
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3633
20971 continue                                                             3634
20972 continue                                                             3634
      deallocate(sxp,b,bs,r,q,mm,is,ga,ixx,gk,del,sxpl)                    3635
      return                                                               3636
      end                                                                  3637
      subroutine multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,   3639 
     *nx,nlam,  flmin,ulam,shri,isd,maxit,xv,xb,xs,lmu,a0,a,m,kin,dev0,d
     *ev,alm,nlp,jerr)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni),   3640 
     *xv(ni)
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(2,ni)          3641
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           3642
      real, dimension (:,:), allocatable :: q,r,b,bs                            
      real, dimension (:), allocatable :: sxp,sxpl,ga,gk,del,sc,svr             
      integer, dimension (:), allocatable :: mm,is,iy,isc                       
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               3651
      exmn=-exmx                                                           3652
      allocate(mm(1:ni),stat=ierr)                                         3652
      jerr=jerr+ierr                                                       3653
      allocate(ga(1:ni),stat=ierr)                                         3653
      jerr=jerr+ierr                                                       3654
      allocate(gk(1:nc),stat=ierr)                                         3654
      jerr=jerr+ierr                                                       3655
      allocate(del(1:nc),stat=ierr)                                        3655
      jerr=jerr+ierr                                                       3656
      allocate(iy(1:ni),stat=ierr)                                         3656
      jerr=jerr+ierr                                                       3657
      allocate(is(1:max(nc,ni)),stat=ierr)                                 3657
      jerr=jerr+ierr                                                       3658
      allocate(sxp(1:no),stat=ierr)                                        3658
      jerr=jerr+ierr                                                       3659
      allocate(sxpl(1:no),stat=ierr)                                       3659
      jerr=jerr+ierr                                                       3660
      allocate(svr(1:nc),stat=ierr)                                        3660
      jerr=jerr+ierr                                                       3661
      allocate(sc(1:no),stat=ierr)                                         3661
      jerr=jerr+ierr                                                       3662
      allocate(isc(1:nc),stat=ierr)                                        3662
      jerr=jerr+ierr                                                       3663
      if(jerr.ne.0) return                                                 3664
      pmax=1.0-pmin                                                        3664
      emin=pmin/pmax                                                       3664
      emax=1.0/emin                                                        3665
      bta=parm                                                             3665
      omb=1.0-bta                                                          3665
      dev1=0.0                                                             3665
      dev0=0.0                                                             3666
20980 do 20981 ic=1,nc                                                     3666
      q0=dot_product(w,y(:,ic))                                            3667
      if(q0 .gt. pmin)goto 21001                                           3667
      jerr =8000+ic                                                        3667
      return                                                               3667
21001 continue                                                             3668
      if(q0 .lt. pmax)goto 21021                                           3668
      jerr =9000+ic                                                        3668
      return                                                               3668
21021 continue                                                             3669
      b(1:ni,ic)=0.0                                                       3669
      b(0,ic)=log(q0)                                                      3669
      dev1=dev1-q0*b(0,ic)                                                 3670
20981 continue                                                             3671
20982 continue                                                             3671
      iy=0                                                                 3671
      al=0.0                                                               3672
      if(nonzero(no*nc,g) .ne. 0)goto 21041                                3673
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         3673
      sxp=0.0                                                              3674
21050 do 21051 ic=1,nc                                                     3674
      q(:,ic)=exp(b(0,ic))                                                 3674
      sxp=sxp+q(:,ic)                                                      3674
21051 continue                                                             3675
21052 continue                                                             3675
      goto 21061                                                           3676
21041 continue                                                             3676
21070 do 21071 i=1,no                                                      3676
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3676
21071 continue                                                             3676
21072 continue                                                             3676
      sxp=0.0                                                              3677
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 3677
      if(jerr.ne.0) return                                                 3678
      dev1=0.0                                                             3679
21080 do 21081 ic=1,nc                                                     3679
      q(:,ic)=b(0,ic)+g(:,ic)                                              3680
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             3681
      q(:,ic)=exp(q(:,ic))                                                 3681
      sxp=sxp+q(:,ic)                                                      3682
21081 continue                                                             3683
21082 continue                                                             3683
      sxpl=w*log(sxp)                                                      3683
21090 do 21091 ic=1,nc                                                     3683
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  3683
21091 continue                                                             3684
21092 continue                                                             3684
21061 continue                                                             3685
21031 continue                                                             3685
21100 do 21101 ic=1,nc                                                     3685
21110 do 21111 i=1,no                                                      3685
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               3685
21111 continue                                                             3685
21112 continue                                                             3685
21101 continue                                                             3686
21102 continue                                                             3686
      dev0=dev0+dev1                                                       3687
      if(flmin .ge. 1.0)goto 21131                                         3687
      eqs=max(eps,flmin)                                                   3687
      alf=eqs**(1.0/(nlam-1))                                              3687
21131 continue                                                             3688
      m=0                                                                  3688
      mm=0                                                                 3688
      nin=0                                                                3688
      nlp=0                                                                3688
      mnl=min(mnlam,nlam)                                                  3688
      bs=0.0                                                               3689
      shr=shri*dev0                                                        3689
      ga=0.0                                                               3690
21140 do 21141 ic=1,nc                                                     3690
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3690
      svr(ic)=sum(r(:,ic))                                                 3691
21150 do 21151 j=1,ni                                                      3691
      if(ju(j).eq.0)goto 21151                                             3692
      jb=ix(j)                                                             3692
      je=ix(j+1)-1                                                         3693
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             3694
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            3695
21151 continue                                                             3696
21152 continue                                                             3696
21141 continue                                                             3697
21142 continue                                                             3697
      ga=sqrt(ga)                                                          3698
21160 do 21161 ilm=1,nlam                                                  3698
      al0=al                                                               3699
      if(flmin .lt. 1.0)goto 21181                                         3699
      al=ulam(ilm)                                                         3699
      goto 21171                                                           3700
21181 if(ilm .le. 2)goto 21191                                             3700
      al=al*alf                                                            3700
      goto 21171                                                           3701
21191 if(ilm .ne. 1)goto 21201                                             3701
      al=big                                                               3701
      goto 21211                                                           3702
21201 continue                                                             3702
      al0=0.0                                                              3703
21220 do 21221 j=1,ni                                                      3703
      if(ju(j).eq.0)goto 21221                                             3703
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3703
21221 continue                                                             3704
21222 continue                                                             3704
      al0=al0/max(bta,1.0e-3)                                              3704
      al=alf*al0                                                           3705
21211 continue                                                             3706
21171 continue                                                             3706
      al2=al*omb                                                           3706
      al1=al*bta                                                           3706
      tlam=bta*(2.0*al-al0)                                                3707
21230 do 21231 k=1,ni                                                      3707
      if(iy(k).eq.1)goto 21231                                             3707
      if(ju(k).eq.0)goto 21231                                             3708
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      3709
21231 continue                                                             3710
21232 continue                                                             3710
10740 continue                                                             3711
21240 continue                                                             3711
21241 continue                                                             3711
      ixx=0                                                                3711
      jxx=ixx                                                              3711
      kxx=jxx                                                              3711
      t=0.0                                                                3712
21250 do 21251 ic=1,nc                                                     3712
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       3712
21251 continue                                                             3713
21252 continue                                                             3713
      if(t .ge. eps)goto 21271                                             3713
      kxx=1                                                                3713
      goto 21242                                                           3713
21271 continue                                                             3713
      t=2.0*t                                                              3713
      alt=al1/t                                                            3713
      al2t=al2/t                                                           3714
21280 do 21281 ic=1,nc                                                     3714
      bs(0,ic)=b(0,ic)                                                     3714
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          3715
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    3715
      svr(ic)=sum(r(:,ic))                                                 3716
      b(0,ic)=b(0,ic)+svr(ic)                                              3716
      r(:,ic)=r(:,ic)-svr(ic)*w                                            3717
      dlx=max(dlx,svr(ic)**2)                                              3718
21281 continue                                                             3719
21282 continue                                                             3719
21290 continue                                                             3719
21291 continue                                                             3719
      nlp=nlp+nc                                                           3719
      dlx=0.0                                                              3720
21300 do 21301 k=1,ni                                                      3720
      if(iy(k).eq.0)goto 21301                                             3721
      jb=ix(k)                                                             3721
      je=ix(k+1)-1                                                         3721
      del=b(k,:)                                                           3721
      gkn=0.0                                                              3722
21310 do 21311 ic=1,nc                                                     3723
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))-svr(ic)*xb(k))/xs(k)        3724
      gk(ic)=u+del(ic)*xv(k)                                               3724
      gkn=gkn+gk(ic)**2                                                    3725
21311 continue                                                             3726
21312 continue                                                             3726
      gkn=sqrt(gkn)                                                        3726
      u=1.0-alt*vp(k)/gkn                                                  3727
      if(u .gt. 0.0)goto 21331                                             3727
      b(k,:)=0.0                                                           3727
      goto 21341                                                           3728
21331 continue                                                             3729
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     3730
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(   3732 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 3733
21341 continue                                                             3734
21321 continue                                                             3734
      del=b(k,:)-del                                                       3734
      if(maxval(abs(del)).le.0.0)goto 21301                                3735
21350 do 21351 ic=1,nc                                                     3735
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3736
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   3738 
     *b(k))/xs(k)
21351 continue                                                             3739
21352 continue                                                             3739
      if(mm(k) .ne. 0)goto 21371                                           3739
      nin=nin+1                                                            3740
      if(nin .le. nx)goto 21391                                            3740
      jxx=1                                                                3740
      goto 21302                                                           3740
21391 continue                                                             3741
      mm(k)=nin                                                            3741
      m(nin)=k                                                             3742
21371 continue                                                             3743
21301 continue                                                             3744
21302 continue                                                             3744
      if(jxx.gt.0)goto 21292                                               3745
      if(dlx.lt.shr)goto 21292                                             3745
      if(nlp .le. maxit)goto 21411                                         3745
      jerr=-ilm                                                            3745
      return                                                               3745
21411 continue                                                             3746
21420 continue                                                             3746
21421 continue                                                             3746
      nlp=nlp+nc                                                           3746
      dlx=0.0                                                              3747
21430 do 21431 l=1,nin                                                     3747
      k=m(l)                                                               3747
      jb=ix(k)                                                             3747
      je=ix(k+1)-1                                                         3747
      del=b(k,:)                                                           3747
      gkn=0.0                                                              3748
21440 do 21441 ic=1,nc                                                     3749
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))  -svr(ic)*xb(k))/xs(k)      3751
      gk(ic)=u+del(ic)*xv(k)                                               3751
      gkn=gkn+gk(ic)**2                                                    3752
21441 continue                                                             3753
21442 continue                                                             3753
      gkn=sqrt(gkn)                                                        3753
      u=1.0-alt*vp(k)/gkn                                                  3754
      if(u .gt. 0.0)goto 21461                                             3754
      b(k,:)=0.0                                                           3754
      goto 21471                                                           3755
21461 continue                                                             3756
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     3757
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(   3759 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 3760
21471 continue                                                             3761
21451 continue                                                             3761
      del=b(k,:)-del                                                       3761
      if(maxval(abs(del)).le.0.0)goto 21431                                3762
21480 do 21481 ic=1,nc                                                     3762
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3763
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   3765 
     *b(k))/xs(k)
21481 continue                                                             3766
21482 continue                                                             3766
21431 continue                                                             3767
21432 continue                                                             3767
      if(dlx.lt.shr)goto 21422                                             3767
      if(nlp .le. maxit)goto 21501                                         3767
      jerr=-ilm                                                            3767
      return                                                               3767
21501 continue                                                             3769
      goto 21421                                                           3770
21422 continue                                                             3770
      goto 21291                                                           3771
21292 continue                                                             3771
      if(jxx.gt.0)goto 21242                                               3772
21510 do 21511 ic=1,nc                                                     3773
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                               3774
      if(ixx .ne. 0)goto 21531                                             3775
21540 do 21541 j=1,nin                                                     3775
      k=m(j)                                                               3776
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 21561                   3776
      ixx=1                                                                3776
      goto 21542                                                           3776
21561 continue                                                             3778
21541 continue                                                             3779
21542 continue                                                             3779
21531 continue                                                             3780
      sc=b(0,ic)+g(:,ic)                                                   3780
      b0=0.0                                                               3781
21570 do 21571 j=1,nin                                                     3781
      l=m(j)                                                               3781
      jb=ix(l)                                                             3781
      je=ix(l+1)-1                                                         3782
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   3783
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            3784
21571 continue                                                             3785
21572 continue                                                             3785
      sc=min(max(exmn,sc+b0),exmx)                                         3786
      sxp=sxp-q(:,ic)                                                      3787
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          3788
      sxp=sxp+q(:,ic)                                                      3789
21511 continue                                                             3790
21512 continue                                                             3790
      s=sum(b(0,:))/nc                                                     3790
      b(0,:)=b(0,:)-s                                                      3791
      if(jxx.gt.0)goto 21242                                               3792
      if(ixx .ne. 0)goto 21591                                             3793
21600 do 21601 j=1,ni                                                      3793
      if(iy(j).eq.1)goto 21601                                             3793
      if(ju(j).eq.0)goto 21601                                             3793
      ga(j)=0.0                                                            3793
21601 continue                                                             3794
21602 continue                                                             3794
21610 do 21611 ic=1,nc                                                     3794
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3795
21620 do 21621 j=1,ni                                                      3795
      if(iy(j).eq.1)goto 21621                                             3795
      if(ju(j).eq.0)goto 21621                                             3796
      jb=ix(j)                                                             3796
      je=ix(j+1)-1                                                         3797
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             3798
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            3799
21621 continue                                                             3800
21622 continue                                                             3800
21611 continue                                                             3801
21612 continue                                                             3801
      ga=sqrt(ga)                                                          3802
21630 do 21631 k=1,ni                                                      3802
      if(iy(k).eq.1)goto 21631                                             3802
      if(ju(k).eq.0)goto 21631                                             3803
      if(ga(k) .le. al1*vp(k))goto 21651                                   3803
      iy(k)=1                                                              3803
      ixx=1                                                                3803
21651 continue                                                             3804
21631 continue                                                             3805
21632 continue                                                             3805
      if(ixx.eq.1) go to 10740                                             3806
      goto 21242                                                           3807
21591 continue                                                             3808
      goto 21241                                                           3809
21242 continue                                                             3809
      if(kxx .le. 0)goto 21671                                             3809
      jerr=-20000-ilm                                                      3809
      goto 21162                                                           3809
21671 continue                                                             3810
      if(jxx .le. 0)goto 21691                                             3810
      jerr=-10000-ilm                                                      3810
      goto 21162                                                           3810
21691 continue                                                             3810
      devi=0.0                                                             3811
21700 do 21701 ic=1,nc                                                     3812
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          3812
      a0(ic,ilm)=b(0,ic)                                                   3813
21710 do 21711 i=1,no                                                      3813
      if(y(i,ic).le.0.0)goto 21711                                         3814
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           3815
21711 continue                                                             3816
21712 continue                                                             3816
21701 continue                                                             3817
21702 continue                                                             3817
      kin(ilm)=nin                                                         3817
      alm(ilm)=al                                                          3817
      lmu=ilm                                                              3818
      dev(ilm)=(dev1-devi)/dev0                                            3819
      if(ilm.lt.mnl)goto 21161                                             3819
      if(flmin.ge.1.0)goto 21161                                           3820
      me=0                                                                 3820
21720 do 21721 j=1,nin                                                     3820
      if(a(j,1,ilm).ne.0.0) me=me+1                                        3820
21721 continue                                                             3820
21722 continue                                                             3820
      if(me.gt.ne)goto 21162                                               3821
      if(dev(ilm).gt.devmax)goto 21162                                     3821
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 21162                             3822
21161 continue                                                             3823
21162 continue                                                             3823
      g=log(q)                                                             3823
21730 do 21731 i=1,no                                                      3823
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3823
21731 continue                                                             3824
21732 continue                                                             3824
      deallocate(sxp,b,bs,r,q,mm,is,sc,ga,iy,gk,del,sxpl)                  3825
      return                                                               3826
      end                                                                  3827
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
