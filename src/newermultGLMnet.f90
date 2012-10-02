c
c                          newGLMnet (9/20/12)
c
c
c                 Elastic net with squared-error loss
c
c dense predictor matrix:
c
c call elnet(ka,parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,maxit,
c            lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call spelnet(ka,parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,
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
c   isd = standarization flag:
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
c call multelnet(parm,no,ni,nr,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,
c                jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call multspelnet(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,
c             isd,jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse row format
c
c other inputs:
c
c   nr = number of response variables
c   y(no,nr) = response data matrix (overwritten)
c   all other inputs same as elnet/spelnet above
c
c output:
c
c   a0(nr,lmu) = intercept values for each solution
c   ca(nx,nr,lmu) = compressed coefficient values for each solution
c   all other outputs same as elnet/spelnet above
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
c call lognet (parm,no,ni,nc,x,y,o,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,
c              maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call splognet (parm,no,ni,nc,x,ix,jx,y,o,jd,vp,ne,nx,nlam,flmin,
c             ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse row format
c
c
c other inputs:
c
c   parm,no,ni,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,maxit
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
c call fishnet (parm,no,ni,x,y,o,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,
c               isd,maxit,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c sparse predictor matrix:
c
c call spfishnet (parm,no,ni,x,ix,jx,y,o,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,
c               isd,maxit,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
c
c    x, ix, jx = predictor data matrix in compressed sparse row format
c
c other inputs:
c
c   y(no) = observation response counts
c   o(no) = observation off-sets
c   parm,no,ni,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,maxit
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
c call coxnet (parm,no,ni,x,y,d,o,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,
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
c   parm,no,ni,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,maxit
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
      subroutine elnet  (ka,parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam    719 
     *,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam)                          720
      real ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)                          721
      integer jd(*),ia(nx),nin(nlam)                                        722
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10021                                     725
      jerr=10000                                                            725
      return                                                                725
10021 continue                                                              726
      allocate(vq(1:ni),stat=jerr)                                          726
      if(jerr.ne.0) return                                                  727
      vq=max(0.0,vp)                                                        727
      vq=vq*ni/sum(vq)                                                      728
      if(ka .ne. 1)goto 10041                                               729
      call elnetu  (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd    732 
     *,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                            733
10041 continue                                                              734
      call elnetn (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd,    737 
     *maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                              738
10031 continue                                                              738
      deallocate(vq)                                                        739
      return                                                                740
      end                                                                   741
      subroutine elnetu  (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,t    744 
     *hr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ulam(nlam)                           745
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         746
      integer jd(*),ia(nx),nin(nlam)                                        747
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           752
      allocate(xm(1:ni),stat=ierr)                                          752
      jerr=jerr+ierr                                                        753
      allocate(xs(1:ni),stat=ierr)                                          753
      jerr=jerr+ierr                                                        754
      allocate(ju(1:ni),stat=ierr)                                          754
      jerr=jerr+ierr                                                        755
      allocate(xv(1:ni),stat=ierr)                                          755
      jerr=jerr+ierr                                                        756
      allocate(vlam(1:nlam),stat=ierr)                                      756
      jerr=jerr+ierr                                                        757
      if(jerr.ne.0) return                                                  758
      call chkvars(no,ni,x,ju)                                              759
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  760
      if(maxval(ju) .gt. 0)goto 10071                                       760
      jerr=7777                                                             760
      return                                                                760
10071 continue                                                              761
      call standard(no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)               762
      if(jerr.ne.0) return                                                  763
      if(flmin.ge.1.0) vlam=ulam/ys                                         764
      call elnet1(parm,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxit,x    766 
     *v,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  767
10080 do 10081 k=1,lmu                                                      767
      alm(k)=ys*alm(k)                                                      767
      nk=nin(k)                                                             768
10090 do 10091 l=1,nk                                                       768
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          768
10091 continue                                                              769
10092 continue                                                              769
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         770
10081 continue                                                              771
10082 continue                                                              771
      deallocate(xm,xs,g,ju,xv,vlam)                                        772
      return                                                                773
      end                                                                   774
      subroutine standard (no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)        775
      real x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                  775
      integer ju(ni)                                                        776
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           779
      if(jerr.ne.0) return                                                  780
      w=w/sum(w)                                                            780
      v=sqrt(w)                                                             781
10100 do 10101 j=1,ni                                                       781
      if(ju(j).eq.0)goto 10101                                              782
      xm(j)=dot_product(w,x(:,j))                                           782
      x(:,j)=v*(x(:,j)-xm(j))                                               783
      xv(j)=dot_product(x(:,j),x(:,j))                                      783
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        784
10101 continue                                                              785
10102 continue                                                              785
      if(isd .ne. 0)goto 10121                                              785
      xs=1.0                                                                785
      goto 10131                                                            786
10121 continue                                                              787
10140 do 10141 j=1,ni                                                       787
      if(ju(j).eq.0)goto 10141                                              787
      x(:,j)=x(:,j)/xs(j)                                                   787
10141 continue                                                              788
10142 continue                                                              788
      xv=1.0                                                                789
10131 continue                                                              790
10111 continue                                                              790
      ym=dot_product(w,y)                                                   790
      y=v*(y-ym)                                                            790
      ys=sqrt(dot_product(y,y))                                             790
      y=y/ys                                                                790
      g=0.0                                                                 791
10150 do 10151 j=1,ni                                                       791
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                             791
10151 continue                                                              792
10152 continue                                                              792
      deallocate(v)                                                         793
      return                                                                794
      end                                                                   795
      subroutine elnet1 (beta,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,ulam,thr,    797 
     *maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    798 
     *9)
      real vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    799 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       800
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                           806
      jerr=jerr+ierr                                                        807
      allocate(mm(1:ni),stat=ierr)                                          807
      jerr=jerr+ierr                                                        808
      allocate(da(1:ni),stat=ierr)                                          808
      jerr=jerr+ierr                                                        809
      if(jerr.ne.0) return                                                  810
      bta=beta                                                              810
      omb=1.0-bta                                                           811
      if(flmin .ge. 1.0)goto 10171                                          811
      eqs=max(eps,flmin)                                                    811
      alf=eqs**(1.0/(nlam-1))                                               811
10171 continue                                                              812
      rsq=0.0                                                               812
      a=0.0                                                                 812
      mm=0                                                                  812
      nlp=0                                                                 812
      nin=nlp                                                               812
      iz=0                                                                  812
      mnl=min(mnlam,nlam)                                                   813
10180 do 10181 m=1,nlam                                                     814
      if(flmin .lt. 1.0)goto 10201                                          814
      alm=ulam(m)                                                           814
      goto 10191                                                            815
10201 if(m .le. 2)goto 10211                                                815
      alm=alm*alf                                                           815
      goto 10191                                                            816
10211 if(m .ne. 1)goto 10221                                                816
      alm=big                                                               816
      goto 10231                                                            817
10221 continue                                                              817
      alm=0.0                                                               818
10240 do 10241 j=1,ni                                                       818
      if(ju(j).eq.0)goto 10241                                              818
      if(vp(j).le.0.0)goto 10241                                            819
      alm=max(alm,abs(g(j))/vp(j))                                          820
10241 continue                                                              821
10242 continue                                                              821
      alm=alf*alm/max(bta,1.0e-3)                                           822
10231 continue                                                              823
10191 continue                                                              823
      dem=alm*omb                                                           823
      ab=alm*bta                                                            823
      rsq0=rsq                                                              823
      jz=1                                                                  824
10250 continue                                                              824
10251 continue                                                              824
      if(iz*jz.ne.0) go to 10260                                            824
      nlp=nlp+1                                                             824
      dlx=0.0                                                               825
10270 do 10271 k=1,ni                                                       825
      if(ju(k).eq.0)goto 10271                                              826
      ak=a(k)                                                               826
      u=g(k)+ak*xv(k)                                                       826
      v=abs(u)-vp(k)*ab                                                     826
      a(k)=0.0                                                              827
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         828
      if(a(k).eq.ak)goto 10271                                              829
      if(mm(k) .ne. 0)goto 10291                                            829
      nin=nin+1                                                             829
      if(nin.gt.nx)goto 10272                                               830
10300 do 10301 j=1,ni                                                       830
      if(ju(j).eq.0)goto 10301                                              831
      if(mm(j) .eq. 0)goto 10321                                            831
      c(j,nin)=c(k,mm(j))                                                   831
      goto 10301                                                            831
10321 continue                                                              832
      if(j .ne. k)goto 10341                                                832
      c(j,nin)=xv(j)                                                        832
      goto 10301                                                            832
10341 continue                                                              833
      c(j,nin)=dot_product(x(:,j),x(:,k))                                   834
10301 continue                                                              835
10302 continue                                                              835
      mm(k)=nin                                                             835
      ia(nin)=k                                                             836
10291 continue                                                              837
      del=a(k)-ak                                                           837
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      838
      dlx=max(xv(k)*del**2,dlx)                                             839
10350 do 10351 j=1,ni                                                       839
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                               839
10351 continue                                                              840
10352 continue                                                              840
10271 continue                                                              841
10272 continue                                                              841
      if(dlx.lt.thr)goto 10252                                              841
      if(nin.gt.nx)goto 10252                                               842
      if(nlp .le. maxit)goto 10371                                          842
      jerr=-m                                                               842
      return                                                                842
10371 continue                                                              843
10260 continue                                                              843
      iz=1                                                                  843
      da(1:nin)=a(ia(1:nin))                                                844
10380 continue                                                              844
10381 continue                                                              844
      nlp=nlp+1                                                             844
      dlx=0.0                                                               845
10390 do 10391 l=1,nin                                                      845
      k=ia(l)                                                               845
      ak=a(k)                                                               845
      u=g(k)+ak*xv(k)                                                       845
      v=abs(u)-vp(k)*ab                                                     846
      a(k)=0.0                                                              847
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         848
      if(a(k).eq.ak)goto 10391                                              849
      del=a(k)-ak                                                           849
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      850
      dlx=max(xv(k)*del**2,dlx)                                             851
10400 do 10401 j=1,nin                                                      851
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                  851
10401 continue                                                              852
10402 continue                                                              852
10391 continue                                                              853
10392 continue                                                              853
      if(dlx.lt.thr)goto 10382                                              853
      if(nlp .le. maxit)goto 10421                                          853
      jerr=-m                                                               853
      return                                                                853
10421 continue                                                              854
      goto 10381                                                            855
10382 continue                                                              855
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                      856
10430 do 10431 j=1,ni                                                       856
      if(mm(j).ne.0)goto 10431                                              857
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))            858
10431 continue                                                              859
10432 continue                                                              859
      jz=0                                                                  860
      goto 10251                                                            861
10252 continue                                                              861
      if(nin .le. nx)goto 10451                                             861
      jerr=-10000-m                                                         861
      goto 10182                                                            861
10451 continue                                                              862
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 862
      kin(m)=nin                                                            863
      rsqo(m)=rsq                                                           863
      almo(m)=alm                                                           863
      lmu=m                                                                 864
      if(m.lt.mnl)goto 10181                                                864
      if(flmin.ge.1.0)goto 10181                                            865
      me=0                                                                  865
10460 do 10461 j=1,nin                                                      865
      if(ao(j,m).ne.0.0) me=me+1                                            865
10461 continue                                                              865
10462 continue                                                              865
      if(me.gt.ne)goto 10182                                                866
      if(rsq-rsq0.lt.sml*rsq)goto 10182                                     866
      if(rsq.gt.rsqmax)goto 10182                                           867
10181 continue                                                              868
10182 continue                                                              868
      deallocate(a,mm,c,da)                                                 869
      return                                                                870
      end                                                                   871
      subroutine elnetn (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,th    873 
     *r,isd,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no),w(no),ulam(nlam)                           874
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         875
      integer jd(*),ia(nx),nin(nlam)                                        876
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          881
      allocate(xs(1:ni),stat=ierr)                                          881
      jerr=jerr+ierr                                                        882
      allocate(ju(1:ni),stat=ierr)                                          882
      jerr=jerr+ierr                                                        883
      allocate(xv(1:ni),stat=ierr)                                          883
      jerr=jerr+ierr                                                        884
      allocate(vlam(1:nlam),stat=ierr)                                      884
      jerr=jerr+ierr                                                        885
      if(jerr.ne.0) return                                                  886
      call chkvars(no,ni,x,ju)                                              887
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  888
      if(maxval(ju) .gt. 0)goto 10481                                       888
      jerr=7777                                                             888
      return                                                                888
10481 continue                                                              889
      call standard1(no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)                890
      if(jerr.ne.0) return                                                  891
      if(flmin.ge.1.0) vlam=ulam/ys                                         892
      call elnet2(parm,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,vlam,thr,maxit,x    894 
     *v,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  895
10490 do 10491 k=1,lmu                                                      895
      alm(k)=ys*alm(k)                                                      895
      nk=nin(k)                                                             896
10500 do 10501 l=1,nk                                                       896
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          896
10501 continue                                                              897
10502 continue                                                              897
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         898
10491 continue                                                              899
10492 continue                                                              899
      deallocate(xm,xs,ju,xv,vlam)                                          900
      return                                                                901
      end                                                                   902
      subroutine standard1 (no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)         903
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)                        903
      integer ju(ni)                                                        904
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           907
      if(jerr.ne.0) return                                                  908
      w=w/sum(w)                                                            908
      v=sqrt(w)                                                             909
10510 do 10511 j=1,ni                                                       909
      if(ju(j).eq.0)goto 10511                                              910
      xm(j)=dot_product(w,x(:,j))                                           910
      x(:,j)=v*(x(:,j)-xm(j))                                               911
      xv(j)=dot_product(x(:,j),x(:,j))                                      911
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        912
10511 continue                                                              913
10512 continue                                                              913
      if(isd .ne. 0)goto 10531                                              913
      xs=1.0                                                                913
      goto 10541                                                            914
10531 continue                                                              914
10550 do 10551 j=1,ni                                                       914
      if(ju(j).eq.0)goto 10551                                              914
      x(:,j)=x(:,j)/xs(j)                                                   914
10551 continue                                                              915
10552 continue                                                              915
      xv=1.0                                                                916
10541 continue                                                              917
10521 continue                                                              917
      ym=dot_product(w,y)                                                   917
      y=v*(y-ym)                                                            917
      ys=sqrt(dot_product(y,y))                                             917
      y=y/ys                                                                918
      deallocate(v)                                                         919
      return                                                                920
      end                                                                   921
      subroutine elnet2(beta,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,ulam,thr,m    923 
     *axit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    924 
     *9)
      real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    925 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       926
      real, dimension (:), allocatable :: a,g                                   
      integer, dimension (:), allocatable :: mm,ix                              
      allocate(a(1:ni),stat=jerr)                                           931
      allocate(mm(1:ni),stat=ierr)                                          931
      jerr=jerr+ierr                                                        932
      allocate(g(1:ni),stat=ierr)                                           932
      jerr=jerr+ierr                                                        933
      allocate(ix(1:ni),stat=ierr)                                          933
      jerr=jerr+ierr                                                        934
      if(jerr.ne.0) return                                                  935
      bta=beta                                                              935
      omb=1.0-bta                                                           935
      ix=0                                                                  936
      if(flmin .ge. 1.0)goto 10571                                          936
      eqs=max(eps,flmin)                                                    936
      alf=eqs**(1.0/(nlam-1))                                               936
10571 continue                                                              937
      rsq=0.0                                                               937
      a=0.0                                                                 937
      mm=0                                                                  937
      nlp=0                                                                 937
      nin=nlp                                                               937
      iz=0                                                                  937
      mnl=min(mnlam,nlam)                                                   937
      alm=0.0                                                               938
10580 do 10581 j=1,ni                                                       938
      if(ju(j).eq.0)goto 10581                                              938
      g(j)=abs(dot_product(y,x(:,j)))                                       938
10581 continue                                                              939
10582 continue                                                              939
10590 do 10591 m=1,nlam                                                     939
      alm0=alm                                                              940
      if(flmin .lt. 1.0)goto 10611                                          940
      alm=ulam(m)                                                           940
      goto 10601                                                            941
10611 if(m .le. 2)goto 10621                                                941
      alm=alm*alf                                                           941
      goto 10601                                                            942
10621 if(m .ne. 1)goto 10631                                                942
      alm=big                                                               942
      goto 10641                                                            943
10631 continue                                                              943
      alm0=0.0                                                              944
10650 do 10651 j=1,ni                                                       944
      if(ju(j).eq.0)goto 10651                                              944
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                            944
10651 continue                                                              945
10652 continue                                                              945
      alm0=alm0/max(bta,1.0e-3)                                             945
      alm=alf*alm0                                                          946
10641 continue                                                              947
10601 continue                                                              947
      dem=alm*omb                                                           947
      ab=alm*bta                                                            947
      rsq0=rsq                                                              947
      jz=1                                                                  948
      tlam=bta*(2.0*alm-alm0)                                               949
10660 do 10661 k=1,ni                                                       949
      if(ix(k).eq.1)goto 10661                                              949
      if(ju(k).eq.0)goto 10661                                              950
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                        951
10661 continue                                                              952
10662 continue                                                              952
10670 continue                                                              952
10671 continue                                                              952
      if(iz*jz.ne.0) go to 10260                                            953
10680 continue                                                              953
      nlp=nlp+1                                                             953
      dlx=0.0                                                               954
10690 do 10691 k=1,ni                                                       954
      if(ix(k).eq.0)goto 10691                                              954
      gk=dot_product(y,x(:,k))                                              955
      ak=a(k)                                                               955
      u=gk+ak*xv(k)                                                         955
      v=abs(u)-vp(k)*ab                                                     955
      a(k)=0.0                                                              956
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         957
      if(a(k).eq.ak)goto 10691                                              958
      if(mm(k) .ne. 0)goto 10711                                            958
      nin=nin+1                                                             958
      if(nin.gt.nx)goto 10692                                               959
      mm(k)=nin                                                             959
      ia(nin)=k                                                             960
10711 continue                                                              961
      del=a(k)-ak                                                           961
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        962
      y=y-del*x(:,k)                                                        962
      dlx=max(xv(k)*del**2,dlx)                                             963
10691 continue                                                              964
10692 continue                                                              964
      if(nin.gt.nx)goto 10672                                               965
      if(dlx .ge. thr)goto 10731                                            965
      ixx=0                                                                 966
10740 do 10741 k=1,ni                                                       966
      if(ix(k).eq.1)goto 10741                                              966
      if(ju(k).eq.0)goto 10741                                              967
      g(k)=abs(dot_product(y,x(:,k)))                                       968
      if(g(k) .le. ab*vp(k))goto 10761                                      968
      ix(k)=1                                                               968
      ixx=1                                                                 968
10761 continue                                                              969
10741 continue                                                              970
10742 continue                                                              970
      if(ixx.eq.1) go to 10680                                              971
      goto 10672                                                            972
10731 continue                                                              973
      if(nlp .le. maxit)goto 10781                                          973
      jerr=-m                                                               973
      return                                                                973
10781 continue                                                              974
10260 continue                                                              974
      iz=1                                                                  975
10790 continue                                                              975
10791 continue                                                              975
      nlp=nlp+1                                                             975
      dlx=0.0                                                               976
10800 do 10801 l=1,nin                                                      976
      k=ia(l)                                                               976
      gk=dot_product(y,x(:,k))                                              977
      ak=a(k)                                                               977
      u=gk+ak*xv(k)                                                         977
      v=abs(u)-vp(k)*ab                                                     977
      a(k)=0.0                                                              978
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         979
      if(a(k).eq.ak)goto 10801                                              980
      del=a(k)-ak                                                           980
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        981
      y=y-del*x(:,k)                                                        981
      dlx=max(xv(k)*del**2,dlx)                                             982
10801 continue                                                              983
10802 continue                                                              983
      if(dlx.lt.thr)goto 10792                                              983
      if(nlp .le. maxit)goto 10821                                          983
      jerr=-m                                                               983
      return                                                                983
10821 continue                                                              984
      goto 10791                                                            985
10792 continue                                                              985
      jz=0                                                                  986
      goto 10671                                                            987
10672 continue                                                              987
      if(nin .le. nx)goto 10841                                             987
      jerr=-10000-m                                                         987
      goto 10592                                                            987
10841 continue                                                              988
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 988
      kin(m)=nin                                                            989
      rsqo(m)=rsq                                                           989
      almo(m)=alm                                                           989
      lmu=m                                                                 990
      if(m.lt.mnl)goto 10591                                                990
      if(flmin.ge.1.0)goto 10591                                            991
      me=0                                                                  991
10850 do 10851 j=1,nin                                                      991
      if(ao(j,m).ne.0.0) me=me+1                                            991
10851 continue                                                              991
10852 continue                                                              991
      if(me.gt.ne)goto 10592                                                992
      if(rsq-rsq0.lt.sml*rsq)goto 10592                                     992
      if(rsq.gt.rsqmax)goto 10592                                           993
10591 continue                                                              994
10592 continue                                                              994
      deallocate(a,mm,g,ix)                                                 995
      return                                                                996
      end                                                                   997
      subroutine chkvars(no,ni,x,ju)                                        998
      real x(no,ni)                                                         998
      integer ju(ni)                                                        999
10860 do 10861 j=1,ni                                                       999
      ju(j)=0                                                               999
      t=x(1,j)                                                             1000
10870 do 10871 i=2,no                                                      1000
      if(x(i,j).eq.t)goto 10871                                            1000
      ju(j)=1                                                              1000
      goto 10872                                                           1000
10871 continue                                                             1001
10872 continue                                                             1001
10861 continue                                                             1002
10862 continue                                                             1002
      return                                                               1003
      end                                                                  1004
      subroutine uncomp(ni,ca,ia,nin,a)                                    1005
      real ca(*),a(ni)                                                     1005
      integer ia(*)                                                        1006
      a=0.0                                                                1006
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                  1007
      return                                                               1008
      end                                                                  1009
      subroutine modval(a0,ca,ia,nin,n,x,f)                                1010
      real ca(nin),x(n,*),f(n)                                             1010
      integer ia(nin)                                                      1011
      f=a0                                                                 1011
      if(nin.le.0) return                                                  1012
10880 do 10881 i=1,n                                                       1012
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      1012
10881 continue                                                             1013
10882 continue                                                             1013
      return                                                               1014
      end                                                                  1015
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,fl   1018 
     *min,ulam,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                              1019
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1020
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1021
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10901                                    1024
      jerr=10000                                                           1024
      return                                                               1024
10901 continue                                                             1025
      allocate(vq(1:ni),stat=jerr)                                         1025
      if(jerr.ne.0) return                                                 1026
      vq=max(0.0,vp)                                                       1026
      vq=vq*ni/sum(vq)                                                     1027
      if(ka .ne. 1)goto 10921                                              1028
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam   1031 
     *,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10931                                                           1032
10921 continue                                                             1033
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam,   1036 
     *thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10931 continue                                                             1037
10911 continue                                                             1037
      deallocate(vq)                                                       1038
      return                                                               1039
      end                                                                  1040
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmi   1043 
     *n,ulam,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                              1044
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1045
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1046
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                          1051
      allocate(xm(1:ni),stat=ierr)                                         1051
      jerr=jerr+ierr                                                       1052
      allocate(xs(1:ni),stat=ierr)                                         1052
      jerr=jerr+ierr                                                       1053
      allocate(ju(1:ni),stat=ierr)                                         1053
      jerr=jerr+ierr                                                       1054
      allocate(xv(1:ni),stat=ierr)                                         1054
      jerr=jerr+ierr                                                       1055
      allocate(vlam(1:nlam),stat=ierr)                                     1055
      jerr=jerr+ierr                                                       1056
      if(jerr.ne.0) return                                                 1057
      call spchkvars(no,ni,x,ix,ju)                                        1058
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1059
      if(maxval(ju) .gt. 0)goto 10951                                      1059
      jerr=7777                                                            1059
      return                                                               1059
10951 continue                                                             1060
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,jerr)      1061
      if(jerr.ne.0) return                                                 1062
      if(flmin.ge.1.0) vlam=ulam/ys                                        1063
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t   1065 
     *hr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1066
10960 do 10961 k=1,lmu                                                     1066
      alm(k)=ys*alm(k)                                                     1066
      nk=nin(k)                                                            1067
10970 do 10971 l=1,nk                                                      1067
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1067
10971 continue                                                             1068
10972 continue                                                             1068
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                        1069
10961 continue                                                             1070
10962 continue                                                             1070
      deallocate(xm,xs,g,ju,xv,vlam)                                       1071
      return                                                               1072
      end                                                                  1073
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,j   1074 
     *err)
      real x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                     1074
      integer ix(*),jx(*),ju(ni)                                           1075
      w=w/sum(w)                                                           1076
10980 do 10981 j=1,ni                                                      1076
      if(ju(j).eq.0)goto 10981                                             1077
      jb=ix(j)                                                             1077
      je=ix(j+1)-1                                                         1077
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1078
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1079
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1080
10981 continue                                                             1081
10982 continue                                                             1081
      if(isd .ne. 0)goto 11001                                             1081
      xs=1.0                                                               1081
      goto 11011                                                           1081
11001 continue                                                             1081
      xv=1.0                                                               1081
11011 continue                                                             1082
10991 continue                                                             1082
      ym=dot_product(w,y)                                                  1082
      y=y-ym                                                               1082
      ys=sqrt(dot_product(w,y**2))                                         1082
      y=y/ys                                                               1082
      g=0.0                                                                1083
11020 do 11021 j=1,ni                                                      1083
      if(ju(j).eq.0)goto 11021                                             1083
      jb=ix(j)                                                             1083
      je=ix(j+1)-1                                                         1084
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)           1085
11021 continue                                                             1086
11022 continue                                                             1086
      return                                                               1087
      end                                                                  1088
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,   1090 
     *ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99   1091 
     *9)
      real g(ni),vp(ni),x(*),ulam(nlam),w(no)                              1092
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)          1093
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1094
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                          1100
      jerr=jerr+ierr                                                       1101
      allocate(mm(1:ni),stat=ierr)                                         1101
      jerr=jerr+ierr                                                       1102
      allocate(da(1:ni),stat=ierr)                                         1102
      jerr=jerr+ierr                                                       1103
      if(jerr.ne.0) return                                                 1104
      bta=beta                                                             1104
      omb=1.0-bta                                                          1105
      if(flmin .ge. 1.0)goto 11041                                         1105
      eqs=max(eps,flmin)                                                   1105
      alf=eqs**(1.0/(nlam-1))                                              1105
11041 continue                                                             1106
      rsq=0.0                                                              1106
      a=0.0                                                                1106
      mm=0                                                                 1106
      nlp=0                                                                1106
      nin=nlp                                                              1106
      iz=0                                                                 1106
      mnl=min(mnlam,nlam)                                                  1107
11050 do 11051 m=1,nlam                                                    1108
      if(flmin .lt. 1.0)goto 11071                                         1108
      alm=ulam(m)                                                          1108
      goto 11061                                                           1109
11071 if(m .le. 2)goto 11081                                               1109
      alm=alm*alf                                                          1109
      goto 11061                                                           1110
11081 if(m .ne. 1)goto 11091                                               1110
      alm=big                                                              1110
      goto 11101                                                           1111
11091 continue                                                             1111
      alm=0.0                                                              1112
11110 do 11111 j=1,ni                                                      1112
      if(ju(j).eq.0)goto 11111                                             1112
      if(vp(j).le.0.0)goto 11111                                           1113
      alm=max(alm,abs(g(j))/vp(j))                                         1114
11111 continue                                                             1115
11112 continue                                                             1115
      alm=alf*alm/max(bta,1.0e-3)                                          1116
11101 continue                                                             1117
11061 continue                                                             1117
      dem=alm*omb                                                          1117
      ab=alm*bta                                                           1117
      rsq0=rsq                                                             1117
      jz=1                                                                 1118
11120 continue                                                             1118
11121 continue                                                             1118
      if(iz*jz.ne.0) go to 10260                                           1118
      nlp=nlp+1                                                            1118
      dlx=0.0                                                              1119
11130 do 11131 k=1,ni                                                      1119
      if(ju(k).eq.0)goto 11131                                             1120
      ak=a(k)                                                              1120
      u=g(k)+ak*xv(k)                                                      1120
      v=abs(u)-vp(k)*ab                                                    1120
      a(k)=0.0                                                             1121
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1122
      if(a(k).eq.ak)goto 11131                                             1123
      if(mm(k) .ne. 0)goto 11151                                           1123
      nin=nin+1                                                            1123
      if(nin.gt.nx)goto 11132                                              1124
11160 do 11161 j=1,ni                                                      1124
      if(ju(j).eq.0)goto 11161                                             1125
      if(mm(j) .eq. 0)goto 11181                                           1125
      c(j,nin)=c(k,mm(j))                                                  1125
      goto 11161                                                           1125
11181 continue                                                             1126
      if(j .ne. k)goto 11201                                               1126
      c(j,nin)=xv(j)                                                       1126
      goto 11161                                                           1126
11201 continue                                                             1127
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))       1129
11161 continue                                                             1130
11162 continue                                                             1130
      mm(k)=nin                                                            1130
      ia(nin)=k                                                            1131
11151 continue                                                             1132
      del=a(k)-ak                                                          1132
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1133
      dlx=max(xv(k)*del**2,dlx)                                            1134
11210 do 11211 j=1,ni                                                      1134
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                              1134
11211 continue                                                             1135
11212 continue                                                             1135
11131 continue                                                             1136
11132 continue                                                             1136
      if(dlx.lt.thr)goto 11122                                             1136
      if(nin.gt.nx)goto 11122                                              1137
      if(nlp .le. maxit)goto 11231                                         1137
      jerr=-m                                                              1137
      return                                                               1137
11231 continue                                                             1138
10260 continue                                                             1138
      iz=1                                                                 1138
      da(1:nin)=a(ia(1:nin))                                               1139
11240 continue                                                             1139
11241 continue                                                             1139
      nlp=nlp+1                                                            1139
      dlx=0.0                                                              1140
11250 do 11251 l=1,nin                                                     1140
      k=ia(l)                                                              1141
      ak=a(k)                                                              1141
      u=g(k)+ak*xv(k)                                                      1141
      v=abs(u)-vp(k)*ab                                                    1141
      a(k)=0.0                                                             1142
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1143
      if(a(k).eq.ak)goto 11251                                             1144
      del=a(k)-ak                                                          1144
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1145
      dlx=max(xv(k)*del**2,dlx)                                            1146
11260 do 11261 j=1,nin                                                     1146
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                 1146
11261 continue                                                             1147
11262 continue                                                             1147
11251 continue                                                             1148
11252 continue                                                             1148
      if(dlx.lt.thr)goto 11242                                             1148
      if(nlp .le. maxit)goto 11281                                         1148
      jerr=-m                                                              1148
      return                                                               1148
11281 continue                                                             1149
      goto 11241                                                           1150
11242 continue                                                             1150
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                     1151
11290 do 11291 j=1,ni                                                      1151
      if(mm(j).ne.0)goto 11291                                             1152
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))           1153
11291 continue                                                             1154
11292 continue                                                             1154
      jz=0                                                                 1155
      goto 11121                                                           1156
11122 continue                                                             1156
      if(nin .le. nx)goto 11311                                            1156
      jerr=-10000-m                                                        1156
      goto 11052                                                           1156
11311 continue                                                             1157
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1157
      kin(m)=nin                                                           1158
      rsqo(m)=rsq                                                          1158
      almo(m)=alm                                                          1158
      lmu=m                                                                1159
      if(m.lt.mnl)goto 11051                                               1159
      if(flmin.ge.1.0)goto 11051                                           1160
      me=0                                                                 1160
11320 do 11321 j=1,nin                                                     1160
      if(ao(j,m).ne.0.0) me=me+1                                           1160
11321 continue                                                             1160
11322 continue                                                             1160
      if(me.gt.ne)goto 11052                                               1161
      if(rsq-rsq0.lt.sml*rsq)goto 11052                                    1161
      if(rsq.gt.rsqmax)goto 11052                                          1162
11051 continue                                                             1163
11052 continue                                                             1163
      deallocate(a,mm,c,da)                                                1164
      return                                                               1165
      end                                                                  1166
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmin,   1168 
     *ulam,  thr,isd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),vp(ni),y(no),w(no),ulam(nlam)                              1169
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1170
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1171
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         1176
      allocate(xs(1:ni),stat=ierr)                                         1176
      jerr=jerr+ierr                                                       1177
      allocate(ju(1:ni),stat=ierr)                                         1177
      jerr=jerr+ierr                                                       1178
      allocate(xv(1:ni),stat=ierr)                                         1178
      jerr=jerr+ierr                                                       1179
      allocate(vlam(1:nlam),stat=ierr)                                     1179
      jerr=jerr+ierr                                                       1180
      if(jerr.ne.0) return                                                 1181
      call spchkvars(no,ni,x,ix,ju)                                        1182
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1183
      if(maxval(ju) .gt. 0)goto 11341                                      1183
      jerr=7777                                                            1183
      return                                                               1183
11341 continue                                                             1184
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,jerr)       1185
      if(jerr.ne.0) return                                                 1186
      if(flmin.ge.1.0) vlam=ulam/ys                                        1187
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t   1189 
     *hr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1190
11350 do 11351 k=1,lmu                                                     1190
      alm(k)=ys*alm(k)                                                     1190
      nk=nin(k)                                                            1191
11360 do 11361 l=1,nk                                                      1191
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1191
11361 continue                                                             1192
11362 continue                                                             1192
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                        1193
11351 continue                                                             1194
11352 continue                                                             1194
      deallocate(xm,xs,ju,xv,vlam)                                         1195
      return                                                               1196
      end                                                                  1197
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,je   1198 
     *rr)
      real x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)                           1198
      integer ix(*),jx(*),ju(ni)                                           1199
      w=w/sum(w)                                                           1200
11370 do 11371 j=1,ni                                                      1200
      if(ju(j).eq.0)goto 11371                                             1201
      jb=ix(j)                                                             1201
      je=ix(j+1)-1                                                         1201
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1202
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1203
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1204
11371 continue                                                             1205
11372 continue                                                             1205
      if(isd .ne. 0)goto 11391                                             1205
      xs=1.0                                                               1205
      goto 11401                                                           1205
11391 continue                                                             1205
      xv=1.0                                                               1205
11401 continue                                                             1206
11381 continue                                                             1206
      ym=dot_product(w,y)                                                  1206
      y=y-ym                                                               1206
      ys=sqrt(dot_product(w,y**2))                                         1206
      y=y/ys                                                               1207
      return                                                               1208
      end                                                                  1209
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,   1211 
     *ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99   1212 
     *9)
      real y(no),w(no),x(*),vp(ni),ulam(nlam)                              1213
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)          1214
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1215
      real, dimension (:), allocatable :: a,g                                   
      integer, dimension (:), allocatable :: mm,iy                              
      allocate(a(1:ni),stat=jerr)                                          1220
      allocate(mm(1:ni),stat=ierr)                                         1220
      jerr=jerr+ierr                                                       1221
      allocate(g(1:ni),stat=ierr)                                          1221
      jerr=jerr+ierr                                                       1222
      allocate(iy(1:ni),stat=ierr)                                         1222
      jerr=jerr+ierr                                                       1223
      if(jerr.ne.0) return                                                 1224
      bta=beta                                                             1224
      omb=1.0-bta                                                          1224
      alm=0.0                                                              1224
      iy=0                                                                 1225
      if(flmin .ge. 1.0)goto 11421                                         1225
      eqs=max(eps,flmin)                                                   1225
      alf=eqs**(1.0/(nlam-1))                                              1225
11421 continue                                                             1226
      rsq=0.0                                                              1226
      a=0.0                                                                1226
      mm=0                                                                 1226
      o=0.0                                                                1226
      nlp=0                                                                1226
      nin=nlp                                                              1226
      iz=0                                                                 1226
      mnl=min(mnlam,nlam)                                                  1227
11430 do 11431 j=1,ni                                                      1227
      if(ju(j).eq.0)goto 11431                                             1228
      jb=ix(j)                                                             1228
      je=ix(j+1)-1                                                         1229
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1230
11431 continue                                                             1231
11432 continue                                                             1231
11440 do 11441 m=1,nlam                                                    1231
      alm0=alm                                                             1232
      if(flmin .lt. 1.0)goto 11461                                         1232
      alm=ulam(m)                                                          1232
      goto 11451                                                           1233
11461 if(m .le. 2)goto 11471                                               1233
      alm=alm*alf                                                          1233
      goto 11451                                                           1234
11471 if(m .ne. 1)goto 11481                                               1234
      alm=big                                                              1234
      goto 11491                                                           1235
11481 continue                                                             1235
      alm0=0.0                                                             1236
11500 do 11501 j=1,ni                                                      1236
      if(ju(j).eq.0)goto 11501                                             1236
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           1236
11501 continue                                                             1237
11502 continue                                                             1237
      alm0=alm0/max(bta,1.0e-3)                                            1237
      alm=alf*alm0                                                         1238
11491 continue                                                             1239
11451 continue                                                             1239
      dem=alm*omb                                                          1239
      ab=alm*bta                                                           1239
      rsq0=rsq                                                             1239
      jz=1                                                                 1240
      tlam=bta*(2.0*alm-alm0)                                              1241
11510 do 11511 k=1,ni                                                      1241
      if(iy(k).eq.1)goto 11511                                             1241
      if(ju(k).eq.0)goto 11511                                             1242
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       1243
11511 continue                                                             1244
11512 continue                                                             1244
11520 continue                                                             1244
11521 continue                                                             1244
      if(iz*jz.ne.0) go to 10260                                           1245
10680 continue                                                             1245
      nlp=nlp+1                                                            1245
      dlx=0.0                                                              1246
11530 do 11531 k=1,ni                                                      1246
      if(iy(k).eq.0)goto 11531                                             1246
      jb=ix(k)                                                             1246
      je=ix(k+1)-1                                                         1247
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1248
      ak=a(k)                                                              1248
      u=gk+ak*xv(k)                                                        1248
      v=abs(u)-vp(k)*ab                                                    1248
      a(k)=0.0                                                             1249
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1250
      if(a(k).eq.ak)goto 11531                                             1251
      if(mm(k) .ne. 0)goto 11551                                           1251
      nin=nin+1                                                            1251
      if(nin.gt.nx)goto 11532                                              1252
      mm(k)=nin                                                            1252
      ia(nin)=k                                                            1253
11551 continue                                                             1254
      del=a(k)-ak                                                          1254
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1255
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1256
      o=o+del*xm(k)/xs(k)                                                  1256
      dlx=max(xv(k)*del**2,dlx)                                            1257
11531 continue                                                             1258
11532 continue                                                             1258
      if(nin.gt.nx)goto 11522                                              1259
      if(dlx .ge. thr)goto 11571                                           1259
      ixx=0                                                                1260
11580 do 11581 j=1,ni                                                      1260
      if(iy(j).eq.1)goto 11581                                             1260
      if(ju(j).eq.0)goto 11581                                             1261
      jb=ix(j)                                                             1261
      je=ix(j+1)-1                                                         1262
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1263
      if(g(j) .le. ab*vp(j))goto 11601                                     1263
      iy(j)=1                                                              1263
      ixx=1                                                                1263
11601 continue                                                             1264
11581 continue                                                             1265
11582 continue                                                             1265
      if(ixx.eq.1) go to 10680                                             1266
      goto 11522                                                           1267
11571 continue                                                             1268
      if(nlp .le. maxit)goto 11621                                         1268
      jerr=-m                                                              1268
      return                                                               1268
11621 continue                                                             1269
10260 continue                                                             1269
      iz=1                                                                 1270
11630 continue                                                             1270
11631 continue                                                             1270
      nlp=nlp+1                                                            1270
      dlx=0.0                                                              1271
11640 do 11641 l=1,nin                                                     1271
      k=ia(l)                                                              1271
      jb=ix(k)                                                             1271
      je=ix(k+1)-1                                                         1272
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1273
      ak=a(k)                                                              1273
      u=gk+ak*xv(k)                                                        1273
      v=abs(u)-vp(k)*ab                                                    1273
      a(k)=0.0                                                             1274
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1275
      if(a(k).eq.ak)goto 11641                                             1276
      del=a(k)-ak                                                          1276
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1277
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1278
      o=o+del*xm(k)/xs(k)                                                  1278
      dlx=max(xv(k)*del**2,dlx)                                            1279
11641 continue                                                             1280
11642 continue                                                             1280
      if(dlx.lt.thr)goto 11632                                             1280
      if(nlp .le. maxit)goto 11661                                         1280
      jerr=-m                                                              1280
      return                                                               1280
11661 continue                                                             1281
      goto 11631                                                           1282
11632 continue                                                             1282
      jz=0                                                                 1283
      goto 11521                                                           1284
11522 continue                                                             1284
      if(nin .le. nx)goto 11681                                            1284
      jerr=-10000-m                                                        1284
      goto 11442                                                           1284
11681 continue                                                             1285
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1285
      kin(m)=nin                                                           1286
      rsqo(m)=rsq                                                          1286
      almo(m)=alm                                                          1286
      lmu=m                                                                1287
      if(m.lt.mnl)goto 11441                                               1287
      if(flmin.ge.1.0)goto 11441                                           1288
      me=0                                                                 1288
11690 do 11691 j=1,nin                                                     1288
      if(ao(j,m).ne.0.0) me=me+1                                           1288
11691 continue                                                             1288
11692 continue                                                             1288
      if(me.gt.ne)goto 11442                                               1289
      if(rsq-rsq0.lt.sml*rsq)goto 11442                                    1289
      if(rsq.gt.rsqmax)goto 11442                                          1290
11441 continue                                                             1291
11442 continue                                                             1291
      deallocate(a,mm,g,iy)                                                1292
      return                                                               1293
      end                                                                  1294
      subroutine spchkvars(no,ni,x,ix,ju)                                  1295
      real x(*)                                                            1295
      integer ix(*),ju(ni)                                                 1296
11700 do 11701 j=1,ni                                                      1296
      ju(j)=0                                                              1296
      jb=ix(j)                                                             1296
      nj=ix(j+1)-jb                                                        1296
      if(nj.eq.0)goto 11701                                                1297
      je=ix(j+1)-1                                                         1298
      if(nj .ge. no)goto 11721                                             1298
11730 do 11731 i=jb,je                                                     1298
      if(x(i).eq.0.0)goto 11731                                            1298
      ju(j)=1                                                              1298
      goto 11732                                                           1298
11731 continue                                                             1298
11732 continue                                                             1298
      goto 11741                                                           1299
11721 continue                                                             1299
      t=x(jb)                                                              1299
11750 do 11751 i=jb+1,je                                                   1299
      if(x(i).eq.t)goto 11751                                              1299
      ju(j)=1                                                              1299
      goto 11752                                                           1299
11751 continue                                                             1299
11752 continue                                                             1299
11741 continue                                                             1300
11711 continue                                                             1300
11701 continue                                                             1301
11702 continue                                                             1301
      return                                                               1302
      end                                                                  1303
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                         1304
      real ca(*),x(*),f(n)                                                 1304
      integer ia(*),ix(*),jx(*)                                            1305
      f=a0                                                                 1306
11760 do 11761 j=1,nin                                                     1306
      k=ia(j)                                                              1306
      kb=ix(k)                                                             1306
      ke=ix(k+1)-1                                                         1307
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                             1308
11761 continue                                                             1309
11762 continue                                                             1309
      return                                                               1310
      end                                                                  1311
      function row_prod(i,j,ia,ja,ra,w)                                    1312
      integer ia(*),ja(*)                                                  1312
      real ra(*),w(*)                                                      1313
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(   1315 
     *i),ia(j+1)-ia(j),w)
      return                                                               1316
      end                                                                  1317
      function dot(x,y,mx,my,nx,ny,w)                                      1318
      real x(*),y(*),w(*)                                                  1318
      integer mx(*),my(*)                                                  1319
      i=1                                                                  1319
      j=i                                                                  1319
      s=0.0                                                                1320
11770 continue                                                             1320
11771 continue                                                             1320
11780 continue                                                             1321
11781 if(mx(i).ge.my(j))goto 11782                                         1321
      i=i+1                                                                1321
      if(i.gt.nx) go to 11790                                              1321
      goto 11781                                                           1322
11782 continue                                                             1322
      if(mx(i).eq.my(j)) go to 11800                                       1323
11810 continue                                                             1323
11811 if(my(j).ge.mx(i))goto 11812                                         1323
      j=j+1                                                                1323
      if(j.gt.ny) go to 11790                                              1323
      goto 11811                                                           1324
11812 continue                                                             1324
      if(mx(i).eq.my(j)) go to 11800                                       1324
      goto 11771                                                           1325
11800 continue                                                             1325
      s=s+w(mx(i))*x(i)*y(j)                                               1326
      i=i+1                                                                1326
      if(i.gt.nx)goto 11772                                                1326
      j=j+1                                                                1326
      if(j.gt.ny)goto 11772                                                1327
      goto 11771                                                           1328
11772 continue                                                             1328
11790 continue                                                             1328
      dot=s                                                                1329
      return                                                               1330
      end                                                                  1331
      subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,ne,nx,nlam,flmin,ulam   1333 
     *,thr,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)             1334
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1335
      integer jd(*),ia(nx),nin(nlam)                                       1336
      real, dimension (:), allocatable :: xm,xs,ww,vq,xv                        
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 11831                                    1340
      jerr=10000                                                           1340
      return                                                               1340
11831 continue                                                             1341
      allocate(ww(1:no),stat=jerr)                                         1342
      allocate(ju(1:ni),stat=ierr)                                         1342
      jerr=jerr+ierr                                                       1343
      allocate(vq(1:ni),stat=ierr)                                         1343
      jerr=jerr+ierr                                                       1344
      allocate(xm(1:ni),stat=ierr)                                         1344
      jerr=jerr+ierr                                                       1345
      if(kopt .ne. 2)goto 11851                                            1345
      allocate(xv(1:ni),stat=ierr)                                         1345
      jerr=jerr+ierr                                                       1345
11851 continue                                                             1346
      if(isd .le. 0)goto 11871                                             1346
      allocate(xs(1:ni),stat=ierr)                                         1346
      jerr=jerr+ierr                                                       1346
11871 continue                                                             1347
      if(jerr.ne.0) return                                                 1348
      call chkvars(no,ni,x,ju)                                             1349
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1350
      if(maxval(ju) .gt. 0)goto 11891                                      1350
      jerr=7777                                                            1350
      return                                                               1350
11891 continue                                                             1351
      vq=max(0.0,vp)                                                       1351
      vq=vq*ni/sum(vq)                                                     1352
11900 do 11901 i=1,no                                                      1352
      ww(i)=sum(y(i,:))                                                    1352
      y(i,:)=y(i,:)/ww(i)                                                  1352
11901 continue                                                             1352
11902 continue                                                             1352
      sw=sum(ww)                                                           1352
      ww=ww/sw                                                             1353
      if(nc .ne. 1)goto 11921                                              1353
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             1354
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,ne,nx,nlam,flmin   1356 
     *,ulam,  thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 11911                                                           1357
11921 if(kopt .ne. 2)goto 11931                                            1357
      call multlstandard1(no,ni,x,ww,ju,isd,xm,xs,xv)                      1358
      call multlognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ula   1360 
     *m,thr,  isd,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 11941                                                           1361
11931 continue                                                             1361
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             1362
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,th   1364 
     *r,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
11941 continue                                                             1365
11911 continue                                                             1365
      if(jerr.gt.0) return                                                 1365
      dev0=2.0*sw*dev0                                                     1366
11950 do 11951 k=1,lmu                                                     1366
      nk=nin(k)                                                            1367
11960 do 11961 ic=1,nc                                                     1367
      if(isd .le. 0)goto 11981                                             1367
11990 do 11991 l=1,nk                                                      1367
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1367
11991 continue                                                             1367
11992 continue                                                             1367
11981 continue                                                             1368
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1369
11961 continue                                                             1370
11962 continue                                                             1370
11951 continue                                                             1371
11952 continue                                                             1371
      deallocate(ww,ju,vq,xm)                                              1371
      if(isd.gt.0) deallocate(xs)                                          1372
      if(kopt.eq.2) deallocate(xv)                                         1373
      return                                                               1374
      end                                                                  1375
      subroutine lstandard1 (no,ni,x,w,ju,isd,xm,xs)                       1376
      real x(no,ni),w(no),xm(ni),xs(ni)                                    1376
      integer ju(ni)                                                       1377
12000 do 12001 j=1,ni                                                      1377
      if(ju(j).eq.0)goto 12001                                             1378
      xm(j)=dot_product(w,x(:,j))                                          1378
      x(:,j)=x(:,j)-xm(j)                                                  1379
      if(isd .le. 0)goto 12021                                             1379
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 1379
      x(:,j)=x(:,j)/xs(j)                                                  1379
12021 continue                                                             1380
12001 continue                                                             1381
12002 continue                                                             1381
      return                                                               1382
      end                                                                  1383
      subroutine multlstandard1 (no,ni,x,w,ju,isd,xm,xs,xv)                1384
      real x(no,ni),w(no),xm(ni),xs(ni),xv(ni)                             1384
      integer ju(ni)                                                       1385
12030 do 12031 j=1,ni                                                      1385
      if(ju(j).eq.0)goto 12031                                             1386
      xm(j)=dot_product(w,x(:,j))                                          1386
      x(:,j)=x(:,j)-xm(j)                                                  1387
      xv(j)=dot_product(w,x(:,j)**2)                                       1388
      if(isd .le. 0)goto 12051                                             1388
      xs(j)=sqrt(xv(j))                                                    1388
      x(:,j)=x(:,j)/xs(j)                                                  1388
      xv(j)=1.0                                                            1388
12051 continue                                                             1389
12031 continue                                                             1390
12032 continue                                                             1390
      return                                                               1391
      end                                                                  1392
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ulam   1394 
     *,shri,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1396 
     *5, devmax=0.999)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    1397
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1398
      integer ju(ni),m(nx),kin(nlam)                                       1399
      real, dimension (:), allocatable :: b,bs,v,r,xv,q,ga                      
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(b(0:ni),stat=jerr)                                          1404
      allocate(xv(1:ni),stat=ierr)                                         1404
      jerr=jerr+ierr                                                       1405
      allocate(ga(1:ni),stat=ierr)                                         1405
      jerr=jerr+ierr                                                       1406
      allocate(bs(0:ni),stat=ierr)                                         1406
      jerr=jerr+ierr                                                       1407
      allocate(mm(1:ni),stat=ierr)                                         1407
      jerr=jerr+ierr                                                       1408
      allocate(ixx(1:ni),stat=ierr)                                        1408
      jerr=jerr+ierr                                                       1409
      allocate(r(1:no),stat=ierr)                                          1409
      jerr=jerr+ierr                                                       1410
      allocate(v(1:no),stat=ierr)                                          1410
      jerr=jerr+ierr                                                       1411
      allocate(q(1:no),stat=ierr)                                          1411
      jerr=jerr+ierr                                                       1412
      if(jerr.ne.0) return                                                 1413
      fmax=log(1.0/pmin-1.0)                                               1413
      fmin=-fmax                                                           1413
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1414
      bta=parm                                                             1414
      omb=1.0-bta                                                          1415
      q0=dot_product(w,y)                                                  1415
      if(q0 .gt. pmin)goto 12071                                           1415
      jerr=8001                                                            1415
      return                                                               1415
12071 continue                                                             1416
      if(q0 .lt. 1.0-pmin)goto 12091                                       1416
      jerr=9001                                                            1416
      return                                                               1416
12091 continue                                                             1417
      ixx=0                                                                1417
      al=0.0                                                               1417
      bz=log(q0/(1.0-q0))                                                  1418
      if(nonzero(no,g) .ne. 0)goto 12111                                   1418
      vi=q0*(1.0-q0)                                                       1418
      b(0)=bz                                                              1418
      v=vi*w                                                               1419
      r=w*(y-q0)                                                           1419
      q=q0                                                                 1419
      xmz=vi                                                               1419
      dev1=-(bz*q0+log(1.0-q0))                                            1420
      goto 12121                                                           1421
12111 continue                                                             1421
      b(0)=azero(no,y,g,w,jerr)                                            1421
      if(jerr.ne.0) return                                                 1422
      q=1.0/(1.0+exp(-b(0)-g))                                             1422
      v=w*q*(1.0-q)                                                        1422
      r=w*(y-q)                                                            1422
      xmz=sum(v)                                                           1423
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1424
12121 continue                                                             1425
12101 continue                                                             1425
      if(kopt .le. 0)goto 12141                                            1426
      if(isd .le. 0)goto 12161                                             1426
      xv=0.25                                                              1426
      goto 12171                                                           1427
12161 continue                                                             1427
12180 do 12181 j=1,ni                                                      1427
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1427
12181 continue                                                             1427
12182 continue                                                             1427
12171 continue                                                             1428
12151 continue                                                             1428
12141 continue                                                             1429
      dev0=dev1                                                            1430
12190 do 12191 i=1,no                                                      1430
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1431
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1432
12191 continue                                                             1433
12192 continue                                                             1433
      if(flmin .ge. 1.0)goto 12211                                         1433
      eqs=max(eps,flmin)                                                   1433
      alf=eqs**(1.0/(nlam-1))                                              1433
12211 continue                                                             1434
      m=0                                                                  1434
      mm=0                                                                 1434
      nlp=0                                                                1434
      nin=nlp                                                              1434
      mnl=min(mnlam,nlam)                                                  1434
      bs=0.0                                                               1434
      b(1:ni)=0.0                                                          1435
      shr=shri*dev0                                                        1436
12220 do 12221 j=1,ni                                                      1436
      if(ju(j).eq.0)goto 12221                                             1436
      ga(j)=abs(dot_product(r,x(:,j)))                                     1436
12221 continue                                                             1437
12222 continue                                                             1437
12230 do 12231 ilm=1,nlam                                                  1437
      al0=al                                                               1438
      if(flmin .lt. 1.0)goto 12251                                         1438
      al=ulam(ilm)                                                         1438
      goto 12241                                                           1439
12251 if(ilm .le. 2)goto 12261                                             1439
      al=al*alf                                                            1439
      goto 12241                                                           1440
12261 if(ilm .ne. 1)goto 12271                                             1440
      al=big                                                               1440
      goto 12281                                                           1441
12271 continue                                                             1441
      al0=0.0                                                              1442
12290 do 12291 j=1,ni                                                      1442
      if(ju(j).eq.0)goto 12291                                             1442
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1442
12291 continue                                                             1443
12292 continue                                                             1443
      al0=al0/max(bta,1.0e-3)                                              1443
      al=alf*al0                                                           1444
12281 continue                                                             1445
12241 continue                                                             1445
      al2=al*omb                                                           1445
      al1=al*bta                                                           1445
      tlam=bta*(2.0*al-al0)                                                1446
12300 do 12301 k=1,ni                                                      1446
      if(ixx(k).eq.1)goto 12301                                            1446
      if(ju(k).eq.0)goto 12301                                             1447
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1448
12301 continue                                                             1449
12302 continue                                                             1449
10680 continue                                                             1450
12310 continue                                                             1450
12311 continue                                                             1450
      bs(0)=b(0)                                                           1450
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1451
      if(kopt .ne. 0)goto 12331                                            1452
12340 do 12341 j=1,ni                                                      1452
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                       1452
12341 continue                                                             1453
12342 continue                                                             1453
12331 continue                                                             1454
12350 continue                                                             1454
12351 continue                                                             1454
      nlp=nlp+1                                                            1454
      dlx=0.0                                                              1455
12360 do 12361 k=1,ni                                                      1455
      if(ixx(k).eq.0)goto 12361                                            1456
      bk=b(k)                                                              1456
      gk=dot_product(r,x(:,k))                                             1457
      u=gk+xv(k)*b(k)                                                      1457
      au=abs(u)-vp(k)*al1                                                  1458
      if(au .gt. 0.0)goto 12381                                            1458
      b(k)=0.0                                                             1458
      goto 12391                                                           1459
12381 continue                                                             1459
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1459
12391 continue                                                             1460
12371 continue                                                             1460
      d=b(k)-bk                                                            1460
      if(abs(d).le.0.0)goto 12361                                          1460
      dlx=max(dlx,xv(k)*d**2)                                              1461
      r=r-d*v*x(:,k)                                                       1462
      if(mm(k) .ne. 0)goto 12411                                           1462
      nin=nin+1                                                            1462
      if(nin.gt.nx)goto 12362                                              1463
      mm(k)=nin                                                            1463
      m(nin)=k                                                             1464
12411 continue                                                             1465
12361 continue                                                             1466
12362 continue                                                             1466
      if(nin.gt.nx)goto 12352                                              1467
      d=sum(r)/xmz                                                         1468
      if(d .eq. 0.0)goto 12431                                             1468
      b(0)=b(0)+d                                                          1468
      dlx=max(dlx,xmz*d**2)                                                1468
      r=r-d*v                                                              1468
12431 continue                                                             1469
      if(dlx.lt.shr)goto 12352                                             1469
      if(nlp .le. maxit)goto 12451                                         1469
      jerr=-ilm                                                            1469
      return                                                               1469
12451 continue                                                             1470
12460 continue                                                             1470
12461 continue                                                             1470
      nlp=nlp+1                                                            1470
      dlx=0.0                                                              1471
12470 do 12471 l=1,nin                                                     1471
      k=m(l)                                                               1471
      bk=b(k)                                                              1472
      gk=dot_product(r,x(:,k))                                             1473
      u=gk+xv(k)*b(k)                                                      1473
      au=abs(u)-vp(k)*al1                                                  1474
      if(au .gt. 0.0)goto 12491                                            1474
      b(k)=0.0                                                             1474
      goto 12501                                                           1475
12491 continue                                                             1475
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1475
12501 continue                                                             1476
12481 continue                                                             1476
      d=b(k)-bk                                                            1476
      if(abs(d).le.0.0)goto 12471                                          1476
      dlx=max(dlx,xv(k)*d**2)                                              1477
      r=r-d*v*x(:,k)                                                       1478
12471 continue                                                             1479
12472 continue                                                             1479
      d=sum(r)/xmz                                                         1480
      if(d .eq. 0.0)goto 12521                                             1480
      b(0)=b(0)+d                                                          1480
      dlx=max(dlx,xmz*d**2)                                                1480
      r=r-d*v                                                              1480
12521 continue                                                             1481
      if(dlx.lt.shr)goto 12462                                             1481
      if(nlp .le. maxit)goto 12541                                         1481
      jerr=-ilm                                                            1481
      return                                                               1481
12541 continue                                                             1482
      goto 12461                                                           1483
12462 continue                                                             1483
      goto 12351                                                           1484
12352 continue                                                             1484
      if(nin.gt.nx)goto 12312                                              1485
12550 do 12551 i=1,no                                                      1485
      fi=b(0)+g(i)                                                         1486
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1487
      if(fi .ge. fmin)goto 12571                                           1487
      q(i)=0.0                                                             1487
      goto 12561                                                           1487
12571 if(fi .le. fmax)goto 12581                                           1487
      q(i)=1.0                                                             1487
      goto 12591                                                           1488
12581 continue                                                             1488
      q(i)=1.0/(1.0+exp(-fi))                                              1488
12591 continue                                                             1489
12561 continue                                                             1489
12551 continue                                                             1490
12552 continue                                                             1490
      v=w*q*(1.0-q)                                                        1490
      xmz=sum(v)                                                           1490
      if(xmz.le.vmin)goto 12312                                            1490
      r=w*(y-q)                                                            1491
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 12611                           1491
      ix=0                                                                 1492
12620 do 12621 j=1,nin                                                     1492
      k=m(j)                                                               1493
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 12621                           1493
      ix=1                                                                 1493
      goto 12622                                                           1494
12621 continue                                                             1495
12622 continue                                                             1495
      if(ix .ne. 0)goto 12641                                              1496
12650 do 12651 k=1,ni                                                      1496
      if(ixx(k).eq.1)goto 12651                                            1496
      if(ju(k).eq.0)goto 12651                                             1497
      ga(k)=abs(dot_product(r,x(:,k)))                                     1498
      if(ga(k) .le. al1*vp(k))goto 12671                                   1498
      ixx(k)=1                                                             1498
      ix=1                                                                 1498
12671 continue                                                             1499
12651 continue                                                             1500
12652 continue                                                             1500
      if(ix.eq.1) go to 10680                                              1501
      goto 12312                                                           1502
12641 continue                                                             1503
12611 continue                                                             1504
      goto 12311                                                           1505
12312 continue                                                             1505
      if(nin .le. nx)goto 12691                                            1505
      jerr=-10000-ilm                                                      1505
      goto 12232                                                           1505
12691 continue                                                             1506
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1506
      kin(ilm)=nin                                                         1507
      a0(ilm)=b(0)                                                         1507
      alm(ilm)=al                                                          1507
      lmu=ilm                                                              1508
      devi=dev2(no,w,y,q,pmin)                                             1509
      dev(ilm)=(dev1-devi)/dev0                                            1509
      if(xmz.le.vmin)goto 12232                                            1510
      if(ilm.lt.mnl)goto 12231                                             1510
      if(flmin.ge.1.0)goto 12231                                           1511
      me=0                                                                 1511
12700 do 12701 j=1,nin                                                     1511
      if(a(j,ilm).ne.0.0) me=me+1                                          1511
12701 continue                                                             1511
12702 continue                                                             1511
      if(me.gt.ne)goto 12232                                               1512
      if(dev(ilm).gt.devmax)goto 12232                                     1512
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12232                             1513
12231 continue                                                             1514
12232 continue                                                             1514
      g=log(q/(1.0-q))                                                     1515
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                                  1516
      return                                                               1517
      end                                                                  1518
      function dev2(n,w,y,p,pmin)                                          1519
      real w(n),y(n),p(n)                                                  1520
      pmax=1.0-pmin                                                        1520
      s=0.0                                                                1521
12710 do 12711 i=1,n                                                       1521
      pi=min(max(pmin,p(i)),pmax)                                          1522
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1523
12711 continue                                                             1524
12712 continue                                                             1524
      dev2=s                                                               1525
      return                                                               1526
      end                                                                  1527
      function azero(n,y,g,q,jerr)                                         1528
      parameter(eps=1.0e-7)                                                1529
      real y(n),g(n),q(n)                                                  1530
      real, dimension (:), allocatable :: e,p,w                                 
      allocate(e(1:n),stat=jerr)                                           1534
      allocate(p(1:n),stat=ierr)                                           1534
      jerr=jerr+ierr                                                       1535
      allocate(w(1:n),stat=ierr)                                           1535
      jerr=jerr+ierr                                                       1536
      if(jerr.ne.0) return                                                 1537
      az=0.0                                                               1537
      e=exp(-g)                                                            1537
      qy=dot_product(q,y)                                                  1537
      p=1.0/(1.0+e)                                                        1538
12720 continue                                                             1538
12721 continue                                                             1538
      w=q*p*(1.0-p)                                                        1539
      d=(qy-dot_product(q,p))/sum(w)                                       1539
      az=az+d                                                              1539
      if(abs(d).lt.eps)goto 12722                                          1540
      ea0=exp(-az)                                                         1540
      p=1.0/(1.0+ea0*e)                                                    1541
      goto 12721                                                           1542
12722 continue                                                             1542
      azero=az                                                             1543
      deallocate(e,p,w)                                                    1544
      return                                                               1545
      end                                                                  1546
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ul   1548 
     *am,shri,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1550 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam)              1551
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1552
      integer ju(ni),m(nx),kin(nlam)                                       1553
      real, dimension (:,:), allocatable :: q                                   
      real, dimension (:), allocatable :: sxp,sxpl                              
      real, dimension (:), allocatable :: di,v,r,ga                             
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is,ixx                          
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no),stat=ierr)                                          1564
      jerr=jerr+ierr                                                       1565
      allocate(v(1:no),stat=ierr)                                          1565
      jerr=jerr+ierr                                                       1566
      allocate(mm(1:ni),stat=ierr)                                         1566
      jerr=jerr+ierr                                                       1567
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1567
      jerr=jerr+ierr                                                       1568
      allocate(sxp(1:no),stat=ierr)                                        1568
      jerr=jerr+ierr                                                       1569
      allocate(sxpl(1:no),stat=ierr)                                       1569
      jerr=jerr+ierr                                                       1570
      allocate(di(1:no),stat=ierr)                                         1570
      jerr=jerr+ierr                                                       1571
      allocate(ga(1:ni),stat=ierr)                                         1571
      jerr=jerr+ierr                                                       1572
      allocate(ixx(1:ni),stat=ierr)                                        1572
      jerr=jerr+ierr                                                       1573
      if(jerr.ne.0) return                                                 1574
      pmax=1.0-pmin                                                        1574
      emin=pmin/pmax                                                       1574
      emax=1.0/emin                                                        1575
      pfm=(1.0+pmin)*pmin                                                  1575
      pfx=(1.0-pmin)*pmax                                                  1575
      vmin=pfm*pmax                                                        1576
      bta=parm                                                             1576
      omb=1.0-bta                                                          1576
      dev1=0.0                                                             1576
      dev0=0.0                                                             1577
12730 do 12731 ic=1,nc                                                     1577
      q0=dot_product(w,y(:,ic))                                            1578
      if(q0 .gt. pmin)goto 12751                                           1578
      jerr =8000+ic                                                        1578
      return                                                               1578
12751 continue                                                             1579
      if(q0 .lt. 1.0-pmin)goto 12771                                       1579
      jerr =9000+ic                                                        1579
      return                                                               1579
12771 continue                                                             1580
      b(0,ic)=log(q0)                                                      1580
      dev1=dev1-q0*b(0,ic)                                                 1580
      b(1:ni,ic)=0.0                                                       1581
12731 continue                                                             1582
12732 continue                                                             1582
      ixx=0                                                                1582
      al=0.0                                                               1583
      if(nonzero(no*nc,g) .ne. 0)goto 12791                                1584
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1584
      sxp=0.0                                                              1585
12800 do 12801 ic=1,nc                                                     1585
      q(:,ic)=exp(b(0,ic))                                                 1585
      sxp=sxp+q(:,ic)                                                      1585
12801 continue                                                             1586
12802 continue                                                             1586
      goto 12811                                                           1587
12791 continue                                                             1587
12820 do 12821 i=1,no                                                      1587
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1587
12821 continue                                                             1587
12822 continue                                                             1587
      sxp=0.0                                                              1588
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1588
      if(jerr.ne.0) return                                                 1589
      dev1=0.0                                                             1590
12830 do 12831 ic=1,nc                                                     1590
      q(:,ic)=b(0,ic)+g(:,ic)                                              1591
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             1592
      q(:,ic)=exp(q(:,ic))                                                 1592
      sxp=sxp+q(:,ic)                                                      1593
12831 continue                                                             1594
12832 continue                                                             1594
      sxpl=w*log(sxp)                                                      1594
12840 do 12841 ic=1,nc                                                     1594
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  1594
12841 continue                                                             1595
12842 continue                                                             1595
12811 continue                                                             1596
12781 continue                                                             1596
12850 do 12851 ic=1,nc                                                     1596
12860 do 12861 i=1,no                                                      1596
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1596
12861 continue                                                             1596
12862 continue                                                             1596
12851 continue                                                             1597
12852 continue                                                             1597
      dev0=dev0+dev1                                                       1598
      if(kopt .le. 0)goto 12881                                            1599
      if(isd .le. 0)goto 12901                                             1599
      xv=0.25                                                              1599
      goto 12911                                                           1600
12901 continue                                                             1600
12920 do 12921 j=1,ni                                                      1600
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1600
12921 continue                                                             1600
12922 continue                                                             1600
12911 continue                                                             1601
12891 continue                                                             1601
12881 continue                                                             1602
      if(flmin .ge. 1.0)goto 12941                                         1602
      eqs=max(eps,flmin)                                                   1602
      alf=eqs**(1.0/(nlam-1))                                              1602
12941 continue                                                             1603
      m=0                                                                  1603
      mm=0                                                                 1603
      nin=0                                                                1603
      nlp=0                                                                1603
      mnl=min(mnlam,nlam)                                                  1603
      bs=0.0                                                               1603
      shr=shri*dev0                                                        1604
      ga=0.0                                                               1605
12950 do 12951 ic=1,nc                                                     1605
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1606
12960 do 12961 j=1,ni                                                      1606
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))           1606
12961 continue                                                             1607
12962 continue                                                             1607
12951 continue                                                             1608
12952 continue                                                             1608
12970 do 12971 ilm=1,nlam                                                  1608
      al0=al                                                               1609
      if(flmin .lt. 1.0)goto 12991                                         1609
      al=ulam(ilm)                                                         1609
      goto 12981                                                           1610
12991 if(ilm .le. 2)goto 13001                                             1610
      al=al*alf                                                            1610
      goto 12981                                                           1611
13001 if(ilm .ne. 1)goto 13011                                             1611
      al=big                                                               1611
      goto 13021                                                           1612
13011 continue                                                             1612
      al0=0.0                                                              1613
13030 do 13031 j=1,ni                                                      1613
      if(ju(j).eq.0)goto 13031                                             1613
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1613
13031 continue                                                             1614
13032 continue                                                             1614
      al0=al0/max(bta,1.0e-3)                                              1614
      al=alf*al0                                                           1615
13021 continue                                                             1616
12981 continue                                                             1616
      al2=al*omb                                                           1616
      al1=al*bta                                                           1616
      tlam=bta*(2.0*al-al0)                                                1617
13040 do 13041 k=1,ni                                                      1617
      if(ixx(k).eq.1)goto 13041                                            1617
      if(ju(k).eq.0)goto 13041                                             1618
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1619
13041 continue                                                             1620
13042 continue                                                             1620
10680 continue                                                             1621
13050 continue                                                             1621
13051 continue                                                             1621
      ix=0                                                                 1621
      jx=ix                                                                1621
      ig=0                                                                 1622
13060 do 13061 ic=1,nc                                                     1622
      bs(0,ic)=b(0,ic)                                                     1623
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1624
      xmz=0.0                                                              1625
13070 do 13071 i=1,no                                                      1625
      pic=q(i,ic)/sxp(i)                                                   1626
      if(pic .ge. pfm)goto 13091                                           1626
      pic=0.0                                                              1626
      v(i)=0.0                                                             1626
      goto 13081                                                           1627
13091 if(pic .le. pfx)goto 13101                                           1627
      pic=1.0                                                              1627
      v(i)=0.0                                                             1627
      goto 13111                                                           1628
13101 continue                                                             1628
      v(i)=w(i)*pic*(1.0-pic)                                              1628
      xmz=xmz+v(i)                                                         1628
13111 continue                                                             1629
13081 continue                                                             1629
      r(i)=w(i)*(y(i,ic)-pic)                                              1630
13071 continue                                                             1631
13072 continue                                                             1631
      if(xmz.le.vmin)goto 13061                                            1631
      ig=1                                                                 1632
      if(kopt .ne. 0)goto 13131                                            1633
13140 do 13141 j=1,ni                                                      1633
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                    1633
13141 continue                                                             1634
13142 continue                                                             1634
13131 continue                                                             1635
13150 continue                                                             1635
13151 continue                                                             1635
      nlp=nlp+1                                                            1635
      dlx=0.0                                                              1636
13160 do 13161 k=1,ni                                                      1636
      if(ixx(k).eq.0)goto 13161                                            1637
      bk=b(k,ic)                                                           1637
      gk=dot_product(r,x(:,k))                                             1638
      u=gk+xv(k,ic)*b(k,ic)                                                1638
      au=abs(u)-vp(k)*al1                                                  1639
      if(au .gt. 0.0)goto 13181                                            1639
      b(k,ic)=0.0                                                          1639
      goto 13191                                                           1640
13181 continue                                                             1640
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1640
13191 continue                                                             1641
13171 continue                                                             1641
      d=b(k,ic)-bk                                                         1641
      if(abs(d).le.0.0)goto 13161                                          1642
      dlx=max(dlx,xv(k,ic)*d**2)                                           1642
      r=r-d*v*x(:,k)                                                       1643
      if(mm(k) .ne. 0)goto 13211                                           1643
      nin=nin+1                                                            1644
      if(nin .le. nx)goto 13231                                            1644
      jx=1                                                                 1644
      goto 13162                                                           1644
13231 continue                                                             1645
      mm(k)=nin                                                            1645
      m(nin)=k                                                             1646
13211 continue                                                             1647
13161 continue                                                             1648
13162 continue                                                             1648
      if(jx.gt.0)goto 13152                                                1649
      d=sum(r)/xmz                                                         1650
      if(d .eq. 0.0)goto 13251                                             1650
      b(0,ic)=b(0,ic)+d                                                    1650
      dlx=max(dlx,xmz*d**2)                                                1650
      r=r-d*v                                                              1650
13251 continue                                                             1651
      if(dlx.lt.shr)goto 13152                                             1652
      if(nlp .le. maxit)goto 13271                                         1652
      jerr=-ilm                                                            1652
      return                                                               1652
13271 continue                                                             1653
13280 continue                                                             1653
13281 continue                                                             1653
      nlp=nlp+1                                                            1653
      dlx=0.0                                                              1654
13290 do 13291 l=1,nin                                                     1654
      k=m(l)                                                               1654
      bk=b(k,ic)                                                           1655
      gk=dot_product(r,x(:,k))                                             1656
      u=gk+xv(k,ic)*b(k,ic)                                                1656
      au=abs(u)-vp(k)*al1                                                  1657
      if(au .gt. 0.0)goto 13311                                            1657
      b(k,ic)=0.0                                                          1657
      goto 13321                                                           1658
13311 continue                                                             1658
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1658
13321 continue                                                             1659
13301 continue                                                             1659
      d=b(k,ic)-bk                                                         1659
      if(abs(d).le.0.0)goto 13291                                          1660
      dlx=max(dlx,xv(k,ic)*d**2)                                           1660
      r=r-d*v*x(:,k)                                                       1661
13291 continue                                                             1662
13292 continue                                                             1662
      d=sum(r)/xmz                                                         1663
      if(d .eq. 0.0)goto 13341                                             1663
      b(0,ic)=b(0,ic)+d                                                    1664
      dlx=max(dlx,xmz*d**2)                                                1664
      r=r-d*v                                                              1665
13341 continue                                                             1666
      if(dlx.lt.shr)goto 13282                                             1666
      if(nlp .le. maxit)goto 13361                                         1666
      jerr=-ilm                                                            1666
      return                                                               1666
13361 continue                                                             1667
      goto 13281                                                           1668
13282 continue                                                             1668
      goto 13151                                                           1669
13152 continue                                                             1669
      if(jx.gt.0)goto 13062                                                1670
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                            1671
      if(ix .ne. 0)goto 13381                                              1672
13390 do 13391 j=1,nin                                                     1672
      k=m(j)                                                               1673
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 13411                1673
      ix=1                                                                 1673
      goto 13392                                                           1673
13411 continue                                                             1674
13391 continue                                                             1675
13392 continue                                                             1675
13381 continue                                                             1676
13420 do 13421 i=1,no                                                      1676
      fi=b(0,ic)+g(i,ic)                                                   1678
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1679
      fi=min(max(exmn,fi),exmx)                                            1679
      sxp(i)=sxp(i)-q(i,ic)                                                1680
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    1681
      sxp(i)=sxp(i)+q(i,ic)                                                1682
13421 continue                                                             1683
13422 continue                                                             1683
13061 continue                                                             1684
13062 continue                                                             1684
      s=-sum(b(0,:))/nc                                                    1684
      b(0,:)=b(0,:)+s                                                      1684
      di=s                                                                 1685
13430 do 13431 j=1,nin                                                     1685
      l=m(j)                                                               1686
      if(vp(l) .gt. 0.0)goto 13451                                         1686
      s=sum(b(l,:))/nc                                                     1686
      goto 13461                                                           1687
13451 continue                                                             1687
      s=elc(parm,nc,b(l,:),is)                                             1687
13461 continue                                                             1688
13441 continue                                                             1688
      b(l,:)=b(l,:)-s                                                      1688
      di=di-s*x(:,l)                                                       1689
13431 continue                                                             1690
13432 continue                                                             1690
      di=exp(di)                                                           1690
      sxp=sxp*di                                                           1690
13470 do 13471 ic=1,nc                                                     1690
      q(:,ic)=q(:,ic)*di                                                   1690
13471 continue                                                             1691
13472 continue                                                             1691
      if(jx.gt.0)goto 13052                                                1691
      if(ig.eq.0)goto 13052                                                1692
      if(ix .ne. 0)goto 13491                                              1693
13500 do 13501 k=1,ni                                                      1693
      if(ixx(k).eq.1)goto 13501                                            1693
      if(ju(k).eq.0)goto 13501                                             1693
      ga(k)=0.0                                                            1693
13501 continue                                                             1694
13502 continue                                                             1694
13510 do 13511 ic=1,nc                                                     1694
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1695
13520 do 13521 k=1,ni                                                      1695
      if(ixx(k).eq.1)goto 13521                                            1695
      if(ju(k).eq.0)goto 13521                                             1696
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                          1697
13521 continue                                                             1698
13522 continue                                                             1698
13511 continue                                                             1699
13512 continue                                                             1699
13530 do 13531 k=1,ni                                                      1699
      if(ixx(k).eq.1)goto 13531                                            1699
      if(ju(k).eq.0)goto 13531                                             1700
      if(ga(k) .le. al1*vp(k))goto 13551                                   1700
      ixx(k)=1                                                             1700
      ix=1                                                                 1700
13551 continue                                                             1701
13531 continue                                                             1702
13532 continue                                                             1702
      if(ix.eq.1) go to 10680                                              1703
      goto 13052                                                           1704
13491 continue                                                             1705
      goto 13051                                                           1706
13052 continue                                                             1706
      if(jx .le. 0)goto 13571                                              1706
      jerr=-10000-ilm                                                      1706
      goto 12972                                                           1706
13571 continue                                                             1706
      devi=0.0                                                             1707
13580 do 13581 ic=1,nc                                                     1708
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1708
      a0(ic,ilm)=b(0,ic)                                                   1709
13590 do 13591 i=1,no                                                      1709
      if(y(i,ic).le.0.0)goto 13591                                         1710
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1711
13591 continue                                                             1712
13592 continue                                                             1712
13581 continue                                                             1713
13582 continue                                                             1713
      kin(ilm)=nin                                                         1713
      alm(ilm)=al                                                          1713
      lmu=ilm                                                              1714
      dev(ilm)=(dev1-devi)/dev0                                            1714
      if(ig.eq.0)goto 12972                                                1715
      if(ilm.lt.mnl)goto 12971                                             1715
      if(flmin.ge.1.0)goto 12971                                           1716
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 12972             1717
      if(dev(ilm).gt.devmax)goto 12972                                     1717
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12972                             1718
12971 continue                                                             1719
12972 continue                                                             1719
      g=log(q)                                                             1719
13600 do 13601 i=1,no                                                      1719
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1719
13601 continue                                                             1720
13602 continue                                                             1720
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                           1721
      return                                                               1722
      end                                                                  1723
      subroutine kazero(kk,n,y,g,q,az,jerr)                                1724
      parameter(eps=1.0e-7)                                                1725
      real y(n,kk),g(n,kk),q(n),az(kk)                                     1726
      real, dimension (:), allocatable :: s                                     
      real, dimension (:,:), allocatable :: e                                   
      allocate(e(1:n,1:kk),stat=jerr)                                           
      allocate(s(1:n),stat=ierr)                                           1731
      jerr=jerr+ierr                                                       1732
      if(jerr.ne.0) return                                                 1733
      az=0.0                                                               1733
      e=exp(g)                                                             1733
13610 do 13611 i=1,n                                                       1733
      s(i)=sum(e(i,:))                                                     1733
13611 continue                                                             1734
13612 continue                                                             1734
13620 continue                                                             1734
13621 continue                                                             1734
      dm=0.0                                                               1735
13630 do 13631 k=1,kk                                                      1735
      t=0.0                                                                1735
      u=t                                                                  1736
13640 do 13641 i=1,n                                                       1736
      pik=e(i,k)/s(i)                                                      1737
      t=t+q(i)*(y(i,k)-pik)                                                1737
      u=u+q(i)*pik*(1.0-pik)                                               1738
13641 continue                                                             1739
13642 continue                                                             1739
      d=t/u                                                                1739
      az(k)=az(k)+d                                                        1739
      ed=exp(d)                                                            1739
      dm=max(dm,abs(d))                                                    1740
13650 do 13651 i=1,n                                                       1740
      z=e(i,k)                                                             1740
      e(i,k)=z*ed                                                          1740
      s(i)=s(i)-z+e(i,k)                                                   1740
13651 continue                                                             1741
13652 continue                                                             1741
13631 continue                                                             1742
13632 continue                                                             1742
      if(dm.lt.eps)goto 13622                                              1742
      goto 13621                                                           1743
13622 continue                                                             1743
      az=az-sum(az)/kk                                                     1744
      deallocate(e,s)                                                      1745
      return                                                               1746
      end                                                                  1747
      function elc(parm,n,a,m)                                             1748
      real a(n)                                                            1748
      integer m(n)                                                         1749
      fn=n                                                                 1749
      am=sum(a)/fn                                                         1750
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 13671                       1750
      elc=am                                                               1750
      return                                                               1750
13671 continue                                                             1751
13680 do 13681 i=1,n                                                       1751
      m(i)=i                                                               1751
13681 continue                                                             1751
13682 continue                                                             1751
      call psort7(a,m,1,n)                                                 1752
      if(a(m(1)) .ne. a(m(n)))goto 13701                                   1752
      elc=a(1)                                                             1752
      return                                                               1752
13701 continue                                                             1753
      if(mod(n,2) .ne. 1)goto 13721                                        1753
      ad=a(m(n/2+1))                                                       1753
      goto 13731                                                           1754
13721 continue                                                             1754
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       1754
13731 continue                                                             1755
13711 continue                                                             1755
      if(parm .ne. 1.0)goto 13751                                          1755
      elc=ad                                                               1755
      return                                                               1755
13751 continue                                                             1756
      b1=min(am,ad)                                                        1756
      b2=max(am,ad)                                                        1756
      k2=1                                                                 1757
13760 continue                                                             1757
13761 if(a(m(k2)).gt.b1)goto 13762                                         1757
      k2=k2+1                                                              1757
      goto 13761                                                           1757
13762 continue                                                             1757
      k1=k2-1                                                              1758
13770 continue                                                             1758
13771 if(a(m(k2)).ge.b2)goto 13772                                         1758
      k2=k2+1                                                              1758
      goto 13771                                                           1759
13772 continue                                                             1759
      r=parm/((1.0-parm)*fn)                                               1759
      is=0                                                                 1759
      sm=n-2*(k1-1)                                                        1760
13780 do 13781 k=k1,k2-1                                                   1760
      sm=sm-2.0                                                            1760
      s=r*sm+am                                                            1761
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 13801                   1761
      is=k                                                                 1761
      goto 13782                                                           1761
13801 continue                                                             1762
13781 continue                                                             1763
13782 continue                                                             1763
      if(is .eq. 0)goto 13821                                              1763
      elc=s                                                                1763
      return                                                               1763
13821 continue                                                             1763
      r2=2.0*r                                                             1763
      s1=a(m(k1))                                                          1763
      am2=2.0*am                                                           1764
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    1764
      elc=s1                                                               1765
13830 do 13831 k=k1+1,k2                                                   1765
      s=a(m(k))                                                            1765
      if(s.eq.s1)goto 13831                                                1766
      c=r2*sum(abs(a-s))+s*(s-am2)                                         1767
      if(c .ge. cri)goto 13851                                             1767
      cri=c                                                                1767
      elc=s                                                                1767
13851 continue                                                             1767
      s1=s                                                                 1768
13831 continue                                                             1769
13832 continue                                                             1769
      return                                                               1770
      end                                                                  1771
      function nintot(ni,nx,nc,a,m,nin,is)                                 1772
      real a(nx,nc)                                                        1772
      integer m(nx),is(ni)                                                 1773
      is=0                                                                 1773
      nintot=0                                                             1774
13860 do 13861 ic=1,nc                                                     1774
13870 do 13871 j=1,nin                                                     1774
      k=m(j)                                                               1774
      if(is(k).ne.0)goto 13871                                             1775
      if(a(j,ic).eq.0.0)goto 13871                                         1775
      is(k)=k                                                              1775
      nintot=nintot+1                                                      1776
13871 continue                                                             1776
13872 continue                                                             1776
13861 continue                                                             1777
13862 continue                                                             1777
      return                                                               1778
      end                                                                  1779
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             1780
      real ca(nx,nc),a(ni,nc)                                              1780
      integer ia(nx)                                                       1781
      a=0.0                                                                1782
13880 do 13881 ic=1,nc                                                     1782
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            1782
13881 continue                                                             1783
13882 continue                                                             1783
      return                                                               1784
      end                                                                  1785
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      1786
      real a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                             1786
      integer ia(nx)                                                       1787
13890 do 13891 i=1,nt                                                      1787
13900 do 13901 ic=1,nc                                                     1787
      ans(ic,i)=a0(ic)                                                     1789
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   1790 
     *:nin)))
13901 continue                                                             1790
13902 continue                                                             1790
13891 continue                                                             1791
13892 continue                                                             1791
      return                                                               1792
      end                                                                  1793
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,ne,nx,nlam,fl   1795 
     *min,  ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      real x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)                 1796
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1797
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1798
      real, dimension (:), allocatable :: xm,xs,ww,vq,xv                        
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 13921                                    1802
      jerr=10000                                                           1802
      return                                                               1802
13921 continue                                                             1803
      allocate(ww(1:no),stat=jerr)                                         1804
      allocate(ju(1:ni),stat=ierr)                                         1804
      jerr=jerr+ierr                                                       1805
      allocate(vq(1:ni),stat=ierr)                                         1805
      jerr=jerr+ierr                                                       1806
      allocate(xm(1:ni),stat=ierr)                                         1806
      jerr=jerr+ierr                                                       1807
      allocate(xs(1:ni),stat=ierr)                                         1807
      jerr=jerr+ierr                                                       1808
      if(kopt .ne. 2)goto 13941                                            1808
      allocate(xv(1:ni),stat=ierr)                                         1808
      jerr=jerr+ierr                                                       1808
13941 continue                                                             1809
      if(jerr.ne.0) return                                                 1810
      call spchkvars(no,ni,x,ix,ju)                                        1811
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1812
      if(maxval(ju) .gt. 0)goto 13961                                      1812
      jerr=7777                                                            1812
      return                                                               1812
13961 continue                                                             1813
      vq=max(0.0,vp)                                                       1813
      vq=vq*ni/sum(vq)                                                     1814
13970 do 13971 i=1,no                                                      1814
      ww(i)=sum(y(i,:))                                                    1814
      y(i,:)=y(i,:)/ww(i)                                                  1814
13971 continue                                                             1814
13972 continue                                                             1814
      sw=sum(ww)                                                           1814
      ww=ww/sw                                                             1815
      if(nc .ne. 1)goto 13991                                              1815
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     1816
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,ne,nx,n   1818 
     *lam,flmin,  ulam,thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,de
     *v,alm,nlp,jerr)
      goto 13981                                                           1819
13991 if(kopt .ne. 2)goto 14001                                            1819
      call multsplstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs,xv)              1820
      call multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,   1822 
     *flmin,ulam, thr,isd,maxit,xv,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,n
     *lp,jerr)
      goto 14011                                                           1823
14001 continue                                                             1823
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     1824
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmi   1826 
     *n,ulam,  thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nl
     *p,jerr)
14011 continue                                                             1827
13981 continue                                                             1827
      if(jerr.gt.0) return                                                 1827
      dev0=2.0*sw*dev0                                                     1828
14020 do 14021 k=1,lmu                                                     1828
      nk=nin(k)                                                            1829
14030 do 14031 ic=1,nc                                                     1829
      if(isd .le. 0)goto 14051                                             1829
14060 do 14061 l=1,nk                                                      1829
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1829
14061 continue                                                             1829
14062 continue                                                             1829
14051 continue                                                             1830
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1831
14031 continue                                                             1832
14032 continue                                                             1832
14021 continue                                                             1833
14022 continue                                                             1833
      deallocate(ww,ju,vq,xm,xs)                                           1833
      if(kopt.eq.2) deallocate(xv)                                         1834
      return                                                               1835
      end                                                                  1836
      subroutine multsplstandard2(no,ni,x,ix,jx,w,ju,isd,xm,xs,xv)         1837
      real x(*),w(no),xm(ni),xs(ni),xv(ni)                                 1837
      integer ix(*),jx(*),ju(ni)                                           1838
14070 do 14071 j=1,ni                                                      1838
      if(ju(j).eq.0)goto 14071                                             1838
      jb=ix(j)                                                             1838
      je=ix(j+1)-1                                                         1839
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1840
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1841
      if(isd .le. 0)goto 14091                                             1841
      xs(j)=sqrt(xv(j))                                                    1841
      xv(j)=1.0                                                            1841
14091 continue                                                             1842
14071 continue                                                             1843
14072 continue                                                             1843
      if(isd.eq.0) xs=1.0                                                  1844
      return                                                               1845
      end                                                                  1846
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,xm,xs)                1847
      real x(*),w(no),xm(ni),xs(ni)                                        1847
      integer ix(*),jx(*),ju(ni)                                           1848
14100 do 14101 j=1,ni                                                      1848
      if(ju(j).eq.0)goto 14101                                             1848
      jb=ix(j)                                                             1848
      je=ix(j+1)-1                                                         1849
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1850
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   1851 
     *)**2)
14101 continue                                                             1852
14102 continue                                                             1852
      if(isd.eq.0) xs=1.0                                                  1853
      return                                                               1854
      end                                                                  1855
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam,   1857 
     *  flmin,ulam,shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,alm
     *,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1859 
     *5, devmax=0.999)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        1860
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1861
      real xb(ni),xs(ni)                                                   1861
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1862
      real, dimension (:), allocatable :: xm,b,bs,v,r,sc,xv,q,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(b(0:ni),stat=jerr)                                          1867
      allocate(xm(0:ni),stat=ierr)                                         1867
      jerr=jerr+ierr                                                       1868
      allocate(xv(1:ni),stat=ierr)                                         1868
      jerr=jerr+ierr                                                       1869
      allocate(bs(0:ni),stat=ierr)                                         1869
      jerr=jerr+ierr                                                       1870
      allocate(ga(1:ni),stat=ierr)                                         1870
      jerr=jerr+ierr                                                       1871
      allocate(mm(1:ni),stat=ierr)                                         1871
      jerr=jerr+ierr                                                       1872
      allocate(ixx(1:ni),stat=ierr)                                        1872
      jerr=jerr+ierr                                                       1873
      allocate(q(1:no),stat=ierr)                                          1873
      jerr=jerr+ierr                                                       1874
      allocate(r(1:no),stat=ierr)                                          1874
      jerr=jerr+ierr                                                       1875
      allocate(v(1:no),stat=ierr)                                          1875
      jerr=jerr+ierr                                                       1876
      allocate(sc(1:no),stat=ierr)                                         1876
      jerr=jerr+ierr                                                       1877
      if(jerr.ne.0) return                                                 1878
      fmax=log(1.0/pmin-1.0)                                               1878
      fmin=-fmax                                                           1878
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1879
      bta=parm                                                             1879
      omb=1.0-bta                                                          1880
      q0=dot_product(w,y)                                                  1880
      if(q0 .gt. pmin)goto 14121                                           1880
      jerr=8001                                                            1880
      return                                                               1880
14121 continue                                                             1881
      if(q0 .lt. 1.0-pmin)goto 14141                                       1881
      jerr=9001                                                            1881
      return                                                               1881
14141 continue                                                             1881
      bz=log(q0/(1.0-q0))                                                  1882
      if(nonzero(no,g) .ne. 0)goto 14161                                   1882
      vi=q0*(1.0-q0)                                                       1882
      b(0)=bz                                                              1882
      v=vi*w                                                               1883
      r=w*(y-q0)                                                           1883
      q=q0                                                                 1883
      xm(0)=vi                                                             1883
      dev1=-(bz*q0+log(1.0-q0))                                            1884
      goto 14171                                                           1885
14161 continue                                                             1885
      b(0)=azero(no,y,g,w,jerr)                                            1885
      if(jerr.ne.0) return                                                 1886
      q=1.0/(1.0+exp(-b(0)-g))                                             1886
      v=w*q*(1.0-q)                                                        1886
      r=w*(y-q)                                                            1886
      xm(0)=sum(v)                                                         1887
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1888
14171 continue                                                             1889
14151 continue                                                             1889
      if(kopt .le. 0)goto 14191                                            1890
      if(isd .le. 0)goto 14211                                             1890
      xv=0.25                                                              1890
      goto 14221                                                           1891
14211 continue                                                             1892
14230 do 14231 j=1,ni                                                      1892
      if(ju(j).eq.0)goto 14231                                             1892
      jb=ix(j)                                                             1892
      je=ix(j+1)-1                                                         1893
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          1894
14231 continue                                                             1895
14232 continue                                                             1895
14221 continue                                                             1896
14201 continue                                                             1896
14191 continue                                                             1897
      b(1:ni)=0.0                                                          1897
      dev0=dev1                                                            1898
14240 do 14241 i=1,no                                                      1898
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1899
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1900
14241 continue                                                             1901
14242 continue                                                             1901
      if(flmin .ge. 1.0)goto 14261                                         1901
      eqs=max(eps,flmin)                                                   1901
      alf=eqs**(1.0/(nlam-1))                                              1901
14261 continue                                                             1902
      m=0                                                                  1902
      mm=0                                                                 1902
      nin=0                                                                1902
      o=0.0                                                                1902
      svr=o                                                                1902
      mnl=min(mnlam,nlam)                                                  1902
      bs=0.0                                                               1902
      nlp=0                                                                1902
      nin=nlp                                                              1903
      shr=shri*dev0                                                        1903
      al=0.0                                                               1903
      ixx=0                                                                1904
14270 do 14271 j=1,ni                                                      1904
      if(ju(j).eq.0)goto 14271                                             1905
      jb=ix(j)                                                             1905
      je=ix(j+1)-1                                                         1905
      jn=ix(j+1)-ix(j)                                                     1906
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1907
      gj=dot_product(sc(1:jn),x(jb:je))                                    1908
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      1909
14271 continue                                                             1910
14272 continue                                                             1910
14280 do 14281 ilm=1,nlam                                                  1910
      al0=al                                                               1911
      if(flmin .lt. 1.0)goto 14301                                         1911
      al=ulam(ilm)                                                         1911
      goto 14291                                                           1912
14301 if(ilm .le. 2)goto 14311                                             1912
      al=al*alf                                                            1912
      goto 14291                                                           1913
14311 if(ilm .ne. 1)goto 14321                                             1913
      al=big                                                               1913
      goto 14331                                                           1914
14321 continue                                                             1914
      al0=0.0                                                              1915
14340 do 14341 j=1,ni                                                      1915
      if(ju(j).eq.0)goto 14341                                             1915
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1915
14341 continue                                                             1916
14342 continue                                                             1916
      al0=al0/max(bta,1.0e-3)                                              1916
      al=alf*al0                                                           1917
14331 continue                                                             1918
14291 continue                                                             1918
      al2=al*omb                                                           1918
      al1=al*bta                                                           1918
      tlam=bta*(2.0*al-al0)                                                1919
14350 do 14351 k=1,ni                                                      1919
      if(ixx(k).eq.1)goto 14351                                            1919
      if(ju(k).eq.0)goto 14351                                             1920
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1921
14351 continue                                                             1922
14352 continue                                                             1922
10680 continue                                                             1923
14360 continue                                                             1923
14361 continue                                                             1923
      bs(0)=b(0)                                                           1923
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1924
14370 do 14371 j=1,ni                                                      1924
      if(ixx(j).eq.0)goto 14371                                            1925
      jb=ix(j)                                                             1925
      je=ix(j+1)-1                                                         1925
      jn=ix(j+1)-ix(j)                                                     1926
      sc(1:jn)=v(jx(jb:je))                                                1927
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 1928
      if(kopt .ne. 0)goto 14391                                            1929
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              1930
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                1931
14391 continue                                                             1932
14371 continue                                                             1933
14372 continue                                                             1933
14400 continue                                                             1933
14401 continue                                                             1933
      nlp=nlp+1                                                            1933
      dlx=0.0                                                              1934
14410 do 14411 k=1,ni                                                      1934
      if(ixx(k).eq.0)goto 14411                                            1935
      jb=ix(k)                                                             1935
      je=ix(k+1)-1                                                         1935
      jn=ix(k+1)-ix(k)                                                     1935
      bk=b(k)                                                              1936
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1937
      gk=dot_product(sc(1:jn),x(jb:je))                                    1938
      gk=(gk-svr*xb(k))/xs(k)                                              1939
      u=gk+xv(k)*b(k)                                                      1939
      au=abs(u)-vp(k)*al1                                                  1940
      if(au .gt. 0.0)goto 14431                                            1940
      b(k)=0.0                                                             1940
      goto 14441                                                           1941
14431 continue                                                             1941
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1941
14441 continue                                                             1942
14421 continue                                                             1942
      d=b(k)-bk                                                            1942
      if(abs(d).le.0.0)goto 14411                                          1942
      dlx=max(dlx,xv(k)*d**2)                                              1943
      if(mm(k) .ne. 0)goto 14461                                           1943
      nin=nin+1                                                            1943
      if(nin.gt.nx)goto 14412                                              1944
      mm(k)=nin                                                            1944
      m(nin)=k                                                             1944
      sc(1:jn)=v(jx(jb:je))                                                1945
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 1946
14461 continue                                                             1947
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1948
      o=o+d*(xb(k)/xs(k))                                                  1949
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1950
14411 continue                                                             1951
14412 continue                                                             1951
      if(nin.gt.nx)goto 14402                                              1952
      d=svr/xm(0)                                                          1953
      if(d .eq. 0.0)goto 14481                                             1953
      b(0)=b(0)+d                                                          1953
      dlx=max(dlx,xm(0)*d**2)                                              1953
      r=r-d*v                                                              1953
14481 continue                                                             1954
      svr=svr-d*xm(0)                                                      1954
      if(dlx.lt.shr)goto 14402                                             1955
      if(nlp .le. maxit)goto 14501                                         1955
      jerr=-ilm                                                            1955
      return                                                               1955
14501 continue                                                             1956
14510 continue                                                             1956
14511 continue                                                             1956
      nlp=nlp+1                                                            1956
      dlx=0.0                                                              1957
14520 do 14521 l=1,nin                                                     1957
      k=m(l)                                                               1957
      jb=ix(k)                                                             1957
      je=ix(k+1)-1                                                         1958
      jn=ix(k+1)-ix(k)                                                     1958
      bk=b(k)                                                              1959
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1960
      gk=dot_product(sc(1:jn),x(jb:je))                                    1961
      gk=(gk-svr*xb(k))/xs(k)                                              1962
      u=gk+xv(k)*b(k)                                                      1962
      au=abs(u)-vp(k)*al1                                                  1963
      if(au .gt. 0.0)goto 14541                                            1963
      b(k)=0.0                                                             1963
      goto 14551                                                           1964
14541 continue                                                             1964
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1964
14551 continue                                                             1965
14531 continue                                                             1965
      d=b(k)-bk                                                            1965
      if(abs(d).le.0.0)goto 14521                                          1965
      dlx=max(dlx,xv(k)*d**2)                                              1966
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1967
      o=o+d*(xb(k)/xs(k))                                                  1968
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1969
14521 continue                                                             1970
14522 continue                                                             1970
      d=svr/xm(0)                                                          1971
      if(d .eq. 0.0)goto 14571                                             1971
      b(0)=b(0)+d                                                          1971
      dlx=max(dlx,xm(0)*d**2)                                              1971
      r=r-d*v                                                              1971
14571 continue                                                             1972
      svr=svr-d*xm(0)                                                      1972
      if(dlx.lt.shr)goto 14512                                             1973
      if(nlp .le. maxit)goto 14591                                         1973
      jerr=-ilm                                                            1973
      return                                                               1973
14591 continue                                                             1974
      goto 14511                                                           1975
14512 continue                                                             1975
      goto 14401                                                           1976
14402 continue                                                             1976
      if(nin.gt.nx)goto 14362                                              1977
      sc=b(0)                                                              1977
      b0=0.0                                                               1978
14600 do 14601 j=1,nin                                                     1978
      l=m(j)                                                               1978
      jb=ix(l)                                                             1978
      je=ix(l+1)-1                                                         1979
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      1980
      b0=b0-b(l)*xb(l)/xs(l)                                               1981
14601 continue                                                             1982
14602 continue                                                             1982
      sc=sc+b0                                                             1983
14610 do 14611 i=1,no                                                      1983
      fi=sc(i)+g(i)                                                        1984
      if(fi .ge. fmin)goto 14631                                           1984
      q(i)=0.0                                                             1984
      goto 14621                                                           1984
14631 if(fi .le. fmax)goto 14641                                           1984
      q(i)=1.0                                                             1984
      goto 14651                                                           1985
14641 continue                                                             1985
      q(i)=1.0/(1.0+exp(-fi))                                              1985
14651 continue                                                             1986
14621 continue                                                             1986
14611 continue                                                             1987
14612 continue                                                             1987
      v=w*q*(1.0-q)                                                        1987
      xm(0)=sum(v)                                                         1987
      if(xm(0).lt.vmin)goto 14362                                          1988
      r=w*(y-q)                                                            1988
      svr=sum(r)                                                           1988
      o=0.0                                                                1989
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 14671                         1989
      kx=0                                                                 1990
14680 do 14681 j=1,nin                                                     1990
      k=m(j)                                                               1991
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 14681                           1991
      kx=1                                                                 1991
      goto 14682                                                           1992
14681 continue                                                             1993
14682 continue                                                             1993
      if(kx .ne. 0)goto 14701                                              1994
14710 do 14711 j=1,ni                                                      1994
      if(ixx(j).eq.1)goto 14711                                            1994
      if(ju(j).eq.0)goto 14711                                             1995
      jb=ix(j)                                                             1995
      je=ix(j+1)-1                                                         1995
      jn=ix(j+1)-ix(j)                                                     1996
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1997
      gj=dot_product(sc(1:jn),x(jb:je))                                    1998
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      1999
      if(ga(j) .le. al1*vp(j))goto 14731                                   1999
      ixx(j)=1                                                             1999
      kx=1                                                                 1999
14731 continue                                                             2000
14711 continue                                                             2001
14712 continue                                                             2001
      if(kx.eq.1) go to 10680                                              2002
      goto 14362                                                           2003
14701 continue                                                             2004
14671 continue                                                             2005
      goto 14361                                                           2006
14362 continue                                                             2006
      if(nin .le. nx)goto 14751                                            2006
      jerr=-10000-ilm                                                      2006
      goto 14282                                                           2006
14751 continue                                                             2007
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                2007
      kin(ilm)=nin                                                         2008
      a0(ilm)=b(0)                                                         2008
      alm(ilm)=al                                                          2008
      lmu=ilm                                                              2009
      devi=dev2(no,w,y,q,pmin)                                             2010
      dev(ilm)=(dev1-devi)/dev0                                            2011
      if(ilm.lt.mnl)goto 14281                                             2011
      if(flmin.ge.1.0)goto 14281                                           2012
      me=0                                                                 2012
14760 do 14761 j=1,nin                                                     2012
      if(a(j,ilm).ne.0.0) me=me+1                                          2012
14761 continue                                                             2012
14762 continue                                                             2012
      if(me.gt.ne)goto 14282                                               2013
      if(dev(ilm).gt.devmax)goto 14282                                     2013
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 14282                             2014
      if(xm(0).lt.vmin)goto 14282                                          2015
14281 continue                                                             2016
14282 continue                                                             2016
      g=log(q/(1.0-q))                                                     2017
      deallocate(xm,b,bs,v,r,sc,xv,q,mm,ga,ixx)                            2018
      return                                                               2019
      end                                                                  2020
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam   2022 
     *,flmin,ulam,  shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,al
     *m,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   2024 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni)    2025
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   2026
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2027
      real, dimension (:,:), allocatable :: q                                   
      real, dimension (:), allocatable :: sxp,sxpl                              
      real, dimension (:), allocatable :: sc,xm,v,r,ga                          
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is,iy                           
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(xm(0:ni),stat=ierr)                                         2038
      jerr=jerr+ierr                                                       2039
      allocate(r(1:no),stat=ierr)                                          2039
      jerr=jerr+ierr                                                       2040
      allocate(v(1:no),stat=ierr)                                          2040
      jerr=jerr+ierr                                                       2041
      allocate(mm(1:ni),stat=ierr)                                         2041
      jerr=jerr+ierr                                                       2042
      allocate(ga(1:ni),stat=ierr)                                         2042
      jerr=jerr+ierr                                                       2043
      allocate(iy(1:ni),stat=ierr)                                         2043
      jerr=jerr+ierr                                                       2044
      allocate(is(1:max(nc,ni)),stat=ierr)                                 2044
      jerr=jerr+ierr                                                       2045
      allocate(sxp(1:no),stat=ierr)                                        2045
      jerr=jerr+ierr                                                       2046
      allocate(sxpl(1:no),stat=ierr)                                       2046
      jerr=jerr+ierr                                                       2047
      allocate(sc(1:no),stat=ierr)                                         2047
      jerr=jerr+ierr                                                       2048
      if(jerr.ne.0) return                                                 2049
      pmax=1.0-pmin                                                        2049
      emin=pmin/pmax                                                       2049
      emax=1.0/emin                                                        2050
      pfm=(1.0+pmin)*pmin                                                  2050
      pfx=(1.0-pmin)*pmax                                                  2050
      vmin=pfm*pmax                                                        2051
      bta=parm                                                             2051
      omb=1.0-bta                                                          2051
      dev1=0.0                                                             2051
      dev0=0.0                                                             2052
14770 do 14771 ic=1,nc                                                     2052
      q0=dot_product(w,y(:,ic))                                            2053
      if(q0 .gt. pmin)goto 14791                                           2053
      jerr =8000+ic                                                        2053
      return                                                               2053
14791 continue                                                             2054
      if(q0 .lt. 1.0-pmin)goto 14811                                       2054
      jerr =9000+ic                                                        2054
      return                                                               2054
14811 continue                                                             2055
      b(1:ni,ic)=0.0                                                       2055
      b(0,ic)=log(q0)                                                      2055
      dev1=dev1-q0*b(0,ic)                                                 2056
14771 continue                                                             2057
14772 continue                                                             2057
      iy=0                                                                 2057
      al=0.0                                                               2058
      if(nonzero(no*nc,g) .ne. 0)goto 14831                                2059
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         2059
      sxp=0.0                                                              2060
14840 do 14841 ic=1,nc                                                     2060
      q(:,ic)=exp(b(0,ic))                                                 2060
      sxp=sxp+q(:,ic)                                                      2060
14841 continue                                                             2061
14842 continue                                                             2061
      goto 14851                                                           2062
14831 continue                                                             2062
14860 do 14861 i=1,no                                                      2062
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2062
14861 continue                                                             2062
14862 continue                                                             2062
      sxp=0.0                                                              2063
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 2063
      if(jerr.ne.0) return                                                 2064
      dev1=0.0                                                             2065
14870 do 14871 ic=1,nc                                                     2065
      q(:,ic)=b(0,ic)+g(:,ic)                                              2066
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             2067
      q(:,ic)=exp(q(:,ic))                                                 2067
      sxp=sxp+q(:,ic)                                                      2068
14871 continue                                                             2069
14872 continue                                                             2069
      sxpl=w*log(sxp)                                                      2069
14880 do 14881 ic=1,nc                                                     2069
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  2069
14881 continue                                                             2070
14882 continue                                                             2070
14851 continue                                                             2071
14821 continue                                                             2071
14890 do 14891 ic=1,nc                                                     2071
14900 do 14901 i=1,no                                                      2071
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               2071
14901 continue                                                             2071
14902 continue                                                             2071
14891 continue                                                             2072
14892 continue                                                             2072
      dev0=dev0+dev1                                                       2073
      if(kopt .le. 0)goto 14921                                            2074
      if(isd .le. 0)goto 14941                                             2074
      xv=0.25                                                              2074
      goto 14951                                                           2075
14941 continue                                                             2076
14960 do 14961 j=1,ni                                                      2076
      if(ju(j).eq.0)goto 14961                                             2076
      jb=ix(j)                                                             2076
      je=ix(j+1)-1                                                         2077
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        2078
14961 continue                                                             2079
14962 continue                                                             2079
14951 continue                                                             2080
14931 continue                                                             2080
14921 continue                                                             2081
      if(flmin .ge. 1.0)goto 14981                                         2081
      eqs=max(eps,flmin)                                                   2081
      alf=eqs**(1.0/(nlam-1))                                              2081
14981 continue                                                             2082
      m=0                                                                  2082
      mm=0                                                                 2082
      nin=0                                                                2082
      nlp=0                                                                2082
      mnl=min(mnlam,nlam)                                                  2082
      bs=0.0                                                               2082
      svr=0.0                                                              2082
      o=0.0                                                                2083
      shr=shri*dev0                                                        2083
      ga=0.0                                                               2084
14990 do 14991 ic=1,nc                                                     2084
      v=q(:,ic)/sxp                                                        2084
      r=w*(y(:,ic)-v)                                                      2084
      v=w*v*(1.0-v)                                                        2085
15000 do 15001 j=1,ni                                                      2085
      if(ju(j).eq.0)goto 15001                                             2086
      jb=ix(j)                                                             2086
      je=ix(j+1)-1                                                         2086
      jn=ix(j+1)-ix(j)                                                     2087
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2088
      gj=dot_product(sc(1:jn),x(jb:je))                                    2089
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2090
15001 continue                                                             2091
15002 continue                                                             2091
14991 continue                                                             2092
14992 continue                                                             2092
15010 do 15011 ilm=1,nlam                                                  2092
      al0=al                                                               2093
      if(flmin .lt. 1.0)goto 15031                                         2093
      al=ulam(ilm)                                                         2093
      goto 15021                                                           2094
15031 if(ilm .le. 2)goto 15041                                             2094
      al=al*alf                                                            2094
      goto 15021                                                           2095
15041 if(ilm .ne. 1)goto 15051                                             2095
      al=big                                                               2095
      goto 15061                                                           2096
15051 continue                                                             2096
      al0=0.0                                                              2097
15070 do 15071 j=1,ni                                                      2097
      if(ju(j).eq.0)goto 15071                                             2097
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2097
15071 continue                                                             2098
15072 continue                                                             2098
      al0=al0/max(bta,1.0e-3)                                              2098
      al=alf*al0                                                           2099
15061 continue                                                             2100
15021 continue                                                             2100
      al2=al*omb                                                           2100
      al1=al*bta                                                           2100
      tlam=bta*(2.0*al-al0)                                                2101
15080 do 15081 k=1,ni                                                      2101
      if(iy(k).eq.1)goto 15081                                             2101
      if(ju(k).eq.0)goto 15081                                             2102
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      2103
15081 continue                                                             2104
15082 continue                                                             2104
10680 continue                                                             2105
15090 continue                                                             2105
15091 continue                                                             2105
      ixx=0                                                                2105
      jxx=ixx                                                              2105
      ig=0                                                                 2106
15100 do 15101 ic=1,nc                                                     2106
      bs(0,ic)=b(0,ic)                                                     2107
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          2108
      xm(0)=0.0                                                            2108
      svr=0.0                                                              2108
      o=0.0                                                                2109
15110 do 15111 i=1,no                                                      2109
      pic=q(i,ic)/sxp(i)                                                   2110
      if(pic .ge. pfm)goto 15131                                           2110
      pic=0.0                                                              2110
      v(i)=0.0                                                             2110
      goto 15121                                                           2111
15131 if(pic .le. pfx)goto 15141                                           2111
      pic=1.0                                                              2111
      v(i)=0.0                                                             2111
      goto 15151                                                           2112
15141 continue                                                             2112
      v(i)=w(i)*pic*(1.0-pic)                                              2112
      xm(0)=xm(0)+v(i)                                                     2112
15151 continue                                                             2113
15121 continue                                                             2113
      r(i)=w(i)*(y(i,ic)-pic)                                              2113
      svr=svr+r(i)                                                         2114
15111 continue                                                             2115
15112 continue                                                             2115
      if(xm(0).le.vmin)goto 15101                                          2115
      ig=1                                                                 2116
15160 do 15161 j=1,ni                                                      2116
      if(iy(j).eq.0)goto 15161                                             2117
      jb=ix(j)                                                             2117
      je=ix(j+1)-1                                                         2118
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             2119
      if(kopt .ne. 0)goto 15181                                            2120
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       2121
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          2122
15181 continue                                                             2123
15161 continue                                                             2124
15162 continue                                                             2124
15190 continue                                                             2124
15191 continue                                                             2124
      nlp=nlp+1                                                            2124
      dlx=0.0                                                              2125
15200 do 15201 k=1,ni                                                      2125
      if(iy(k).eq.0)goto 15201                                             2126
      jb=ix(k)                                                             2126
      je=ix(k+1)-1                                                         2126
      jn=ix(k+1)-ix(k)                                                     2126
      bk=b(k,ic)                                                           2127
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2128
      gk=dot_product(sc(1:jn),x(jb:je))                                    2129
      gk=(gk-svr*xb(k))/xs(k)                                              2130
      u=gk+xv(k,ic)*b(k,ic)                                                2130
      au=abs(u)-vp(k)*al1                                                  2131
      if(au .gt. 0.0)goto 15221                                            2131
      b(k,ic)=0.0                                                          2131
      goto 15231                                                           2132
15221 continue                                                             2132
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              2132
15231 continue                                                             2133
15211 continue                                                             2133
      d=b(k,ic)-bk                                                         2133
      if(abs(d).le.0.0)goto 15201                                          2134
      dlx=max(dlx,xv(k,ic)*d**2)                                           2135
      if(mm(k) .ne. 0)goto 15251                                           2135
      nin=nin+1                                                            2136
      if(nin .le. nx)goto 15271                                            2136
      jxx=1                                                                2136
      goto 15202                                                           2136
15271 continue                                                             2137
      mm(k)=nin                                                            2137
      m(nin)=k                                                             2138
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             2139
15251 continue                                                             2140
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2141
      o=o+d*(xb(k)/xs(k))                                                  2142
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2143
15201 continue                                                             2144
15202 continue                                                             2144
      if(jxx.gt.0)goto 15192                                               2145
      d=svr/xm(0)                                                          2146
      if(d .eq. 0.0)goto 15291                                             2146
      b(0,ic)=b(0,ic)+d                                                    2146
      dlx=max(dlx,xm(0)*d**2)                                              2147
      r=r-d*v                                                              2147
      svr=svr-d*xm(0)                                                      2148
15291 continue                                                             2149
      if(dlx.lt.shr)goto 15192                                             2149
      if(nlp .le. maxit)goto 15311                                         2149
      jerr=-ilm                                                            2149
      return                                                               2149
15311 continue                                                             2150
15320 continue                                                             2150
15321 continue                                                             2150
      nlp=nlp+1                                                            2150
      dlx=0.0                                                              2151
15330 do 15331 l=1,nin                                                     2151
      k=m(l)                                                               2151
      jb=ix(k)                                                             2151
      je=ix(k+1)-1                                                         2152
      jn=ix(k+1)-ix(k)                                                     2152
      bk=b(k,ic)                                                           2153
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2154
      gk=dot_product(sc(1:jn),x(jb:je))                                    2155
      gk=(gk-svr*xb(k))/xs(k)                                              2156
      u=gk+xv(k,ic)*b(k,ic)                                                2156
      au=abs(u)-vp(k)*al1                                                  2157
      if(au .gt. 0.0)goto 15351                                            2157
      b(k,ic)=0.0                                                          2157
      goto 15361                                                           2158
15351 continue                                                             2158
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              2158
15361 continue                                                             2159
15341 continue                                                             2159
      d=b(k,ic)-bk                                                         2159
      if(abs(d).le.0.0)goto 15331                                          2160
      dlx=max(dlx,xv(k,ic)*d**2)                                           2161
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2162
      o=o+d*(xb(k)/xs(k))                                                  2163
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2164
15331 continue                                                             2165
15332 continue                                                             2165
      d=svr/xm(0)                                                          2166
      if(d .eq. 0.0)goto 15381                                             2166
      b(0,ic)=b(0,ic)+d                                                    2166
      dlx=max(dlx,xm(0)*d**2)                                              2167
      r=r-d*v                                                              2167
      svr=svr-d*xm(0)                                                      2168
15381 continue                                                             2169
      if(dlx.lt.shr)goto 15322                                             2169
      if(nlp .le. maxit)goto 15401                                         2169
      jerr=-ilm                                                            2169
      return                                                               2169
15401 continue                                                             2170
      goto 15321                                                           2171
15322 continue                                                             2171
      goto 15191                                                           2172
15192 continue                                                             2172
      if(jxx.gt.0)goto 15102                                               2173
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                         2174
      if(ixx .ne. 0)goto 15421                                             2175
15430 do 15431 j=1,nin                                                     2175
      k=m(j)                                                               2176
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 15451                2176
      ixx=1                                                                2176
      goto 15432                                                           2176
15451 continue                                                             2177
15431 continue                                                             2178
15432 continue                                                             2178
15421 continue                                                             2179
      sc=b(0,ic)+g(:,ic)                                                   2179
      b0=0.0                                                               2180
15460 do 15461 j=1,nin                                                     2180
      l=m(j)                                                               2180
      jb=ix(l)                                                             2180
      je=ix(l+1)-1                                                         2181
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   2182
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            2183
15461 continue                                                             2184
15462 continue                                                             2184
      sc=min(max(exmn,sc+b0),exmx)                                         2185
      sxp=sxp-q(:,ic)                                                      2186
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          2187
      sxp=sxp+q(:,ic)                                                      2188
15101 continue                                                             2189
15102 continue                                                             2189
      s=-sum(b(0,:))/nc                                                    2189
      b(0,:)=b(0,:)+s                                                      2189
      sc=s                                                                 2189
      b0=0.0                                                               2190
15470 do 15471 j=1,nin                                                     2190
      l=m(j)                                                               2191
      if(vp(l) .gt. 0.0)goto 15491                                         2191
      s=sum(b(l,:))/nc                                                     2191
      goto 15501                                                           2192
15491 continue                                                             2192
      s=elc(parm,nc,b(l,:),is)                                             2192
15501 continue                                                             2193
15481 continue                                                             2193
      b(l,:)=b(l,:)-s                                                      2194
      jb=ix(l)                                                             2194
      je=ix(l+1)-1                                                         2195
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         2196
      b0=b0+s*xb(l)/xs(l)                                                  2197
15471 continue                                                             2198
15472 continue                                                             2198
      sc=sc+b0                                                             2198
      sc=exp(sc)                                                           2198
      sxp=sxp*sc                                                           2198
15510 do 15511 ic=1,nc                                                     2198
      q(:,ic)=q(:,ic)*sc                                                   2198
15511 continue                                                             2199
15512 continue                                                             2199
      if(jxx.gt.0)goto 15092                                               2199
      if(ig.eq.0)goto 15092                                                2200
      if(ixx .ne. 0)goto 15531                                             2201
15540 do 15541 j=1,ni                                                      2201
      if(iy(j).eq.1)goto 15541                                             2201
      if(ju(j).eq.0)goto 15541                                             2201
      ga(j)=0.0                                                            2201
15541 continue                                                             2202
15542 continue                                                             2202
15550 do 15551 ic=1,nc                                                     2202
      v=q(:,ic)/sxp                                                        2202
      r=w*(y(:,ic)-v)                                                      2202
      v=w*v*(1.0-v)                                                        2203
15560 do 15561 j=1,ni                                                      2203
      if(iy(j).eq.1)goto 15561                                             2203
      if(ju(j).eq.0)goto 15561                                             2204
      jb=ix(j)                                                             2204
      je=ix(j+1)-1                                                         2204
      jn=ix(j+1)-ix(j)                                                     2205
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2206
      gj=dot_product(sc(1:jn),x(jb:je))                                    2207
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2208
15561 continue                                                             2209
15562 continue                                                             2209
15551 continue                                                             2210
15552 continue                                                             2210
15570 do 15571 k=1,ni                                                      2210
      if(iy(k).eq.1)goto 15571                                             2210
      if(ju(k).eq.0)goto 15571                                             2211
      if(ga(k) .le. al1*vp(k))goto 15591                                   2211
      iy(k)=1                                                              2211
      ixx=1                                                                2211
15591 continue                                                             2212
15571 continue                                                             2213
15572 continue                                                             2213
      if(ixx.eq.1) go to 10680                                             2214
      goto 15092                                                           2215
15531 continue                                                             2216
      goto 15091                                                           2217
15092 continue                                                             2217
      if(jxx .le. 0)goto 15611                                             2217
      jerr=-10000-ilm                                                      2217
      goto 15012                                                           2217
15611 continue                                                             2217
      devi=0.0                                                             2218
15620 do 15621 ic=1,nc                                                     2219
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          2219
      a0(ic,ilm)=b(0,ic)                                                   2220
15630 do 15631 i=1,no                                                      2220
      if(y(i,ic).le.0.0)goto 15631                                         2221
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           2222
15631 continue                                                             2223
15632 continue                                                             2223
15621 continue                                                             2224
15622 continue                                                             2224
      kin(ilm)=nin                                                         2224
      alm(ilm)=al                                                          2224
      lmu=ilm                                                              2225
      dev(ilm)=(dev1-devi)/dev0                                            2225
      if(ig.eq.0)goto 15012                                                2226
      if(ilm.lt.mnl)goto 15011                                             2226
      if(flmin.ge.1.0)goto 15011                                           2227
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 15012             2228
      if(dev(ilm).gt.devmax)goto 15012                                     2228
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 15012                             2229
15011 continue                                                             2230
15012 continue                                                             2230
      g=log(q)                                                             2230
15640 do 15641 i=1,no                                                      2230
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2230
15641 continue                                                             2231
15642 continue                                                             2231
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy)                      2232
      return                                                               2233
      end                                                                  2234
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  2235
      real a0(nc),ca(nx,nc),x(*),f(nc,n)                                   2235
      integer ia(*),ix(*),jx(*)                                            2236
15650 do 15651 ic=1,nc                                                     2236
      f(ic,:)=a0(ic)                                                       2236
15651 continue                                                             2237
15652 continue                                                             2237
15660 do 15661 j=1,nin                                                     2237
      k=ia(j)                                                              2237
      kb=ix(k)                                                             2237
      ke=ix(k+1)-1                                                         2238
15670 do 15671 ic=1,nc                                                     2238
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    2238
15671 continue                                                             2239
15672 continue                                                             2239
15661 continue                                                             2240
15662 continue                                                             2240
      return                                                               2241
      end                                                                  2242
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,ne,nx,nlam,flmin,ula   2244 
     *m,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam)              2245
      real ca(nx,nlam),dev(nlam),alm(nlam)                                 2246
      integer jd(*),ia(nx),nin(nlam)                                       2247
      real, dimension (:), allocatable :: xs,ww,vq                              
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 15691                                    2251
      jerr=10000                                                           2251
      return                                                               2251
15691 continue                                                             2252
      allocate(ww(1:no),stat=jerr)                                         2253
      allocate(ju(1:ni),stat=ierr)                                         2253
      jerr=jerr+ierr                                                       2254
      allocate(vq(1:ni),stat=ierr)                                         2254
      jerr=jerr+ierr                                                       2255
      if(isd .le. 0)goto 15711                                             2255
      allocate(xs(1:ni),stat=ierr)                                         2255
      jerr=jerr+ierr                                                       2255
15711 continue                                                             2256
      if(jerr.ne.0) return                                                 2257
      call chkvars(no,ni,x,ju)                                             2258
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2259
      if(maxval(ju) .gt. 0)goto 15731                                      2259
      jerr=7777                                                            2259
      return                                                               2259
15731 continue                                                             2260
      vq=max(0.0,vp)                                                       2260
      vq=vq*ni/sum(vq)                                                     2261
      ww=max(0.0,w)                                                        2261
      sw=sum(ww)                                                           2262
      if(sw .gt. 0.0)goto 15751                                            2262
      jerr=9999                                                            2262
      return                                                               2262
15751 continue                                                             2262
      ww=ww/sw                                                             2263
      call cstandard(no,ni,x,ww,ju,isd,xs)                                 2264
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr   2266 
     *,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 2266
      dev0=2.0*sw*dev0                                                     2267
      if(isd .le. 0)goto 15771                                             2267
15780 do 15781 k=1,lmu                                                     2267
      nk=nin(k)                                                            2267
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   2267
15781 continue                                                             2267
15782 continue                                                             2267
15771 continue                                                             2268
      deallocate(ww,ju,vq)                                                 2268
      if(isd.gt.0) deallocate(xs)                                          2269
      return                                                               2270
      end                                                                  2271
      subroutine cstandard (no,ni,x,w,ju,isd,xs)                           2272
      real x(no,ni),w(no),xs(ni)                                           2272
      integer ju(ni)                                                       2273
15790 do 15791 j=1,ni                                                      2273
      if(ju(j).eq.0)goto 15791                                             2274
      xm=dot_product(w,x(:,j))                                             2274
      x(:,j)=x(:,j)-xm                                                     2275
      if(isd .le. 0)goto 15811                                             2275
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 2275
      x(:,j)=x(:,j)/xs(j)                                                  2275
15811 continue                                                             2276
15791 continue                                                             2277
15792 continue                                                             2277
      return                                                               2278
      end                                                                  2279
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,ne,nx,nlam,flmin,ula   2281 
     *m,cthri,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=0.001, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99)   2282
      real x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam)              2283
      real ao(nx,nlam),dev(nlam),alm(nlam)                                 2284
      integer ju(ni),m(nx),kin(nlam)                                       2285
      real, dimension (:), allocatable :: w,dk,v,xs,wr,a,as,f,dq                
      real, dimension (:), allocatable :: e,uu,ga                               
      integer, dimension (:), allocatable :: jp,kp,mm,ixx                       
      allocate(e(1:no),stat=jerr)                                          2291
      allocate(uu(1:no),stat=ierr)                                         2291
      jerr=jerr+ierr                                                       2292
      allocate(f(1:no),stat=ierr)                                          2292
      jerr=jerr+ierr                                                       2293
      allocate(w(1:no),stat=ierr)                                          2293
      jerr=jerr+ierr                                                       2294
      allocate(v(1:ni),stat=ierr)                                          2294
      jerr=jerr+ierr                                                       2295
      allocate(a(1:ni),stat=ierr)                                          2295
      jerr=jerr+ierr                                                       2296
      allocate(as(1:ni),stat=ierr)                                         2296
      jerr=jerr+ierr                                                       2297
      allocate(xs(1:ni),stat=ierr)                                         2297
      jerr=jerr+ierr                                                       2298
      allocate(ga(1:ni),stat=ierr)                                         2298
      jerr=jerr+ierr                                                       2299
      allocate(ixx(1:ni),stat=ierr)                                        2299
      jerr=jerr+ierr                                                       2300
      allocate(jp(1:no),stat=ierr)                                         2300
      jerr=jerr+ierr                                                       2301
      allocate(kp(1:no),stat=ierr)                                         2301
      jerr=jerr+ierr                                                       2302
      allocate(dk(1:no),stat=ierr)                                         2302
      jerr=jerr+ierr                                                       2303
      allocate(wr(1:no),stat=ierr)                                         2303
      jerr=jerr+ierr                                                       2304
      allocate(dq(1:no),stat=ierr)                                         2304
      jerr=jerr+ierr                                                       2305
      allocate(mm(1:ni),stat=ierr)                                         2305
      jerr=jerr+ierr                                                       2306
      if(jerr.ne.0)go to 11790                                             2307
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2308
      if(jerr.ne.0) go to 11790                                            2308
      alpha=parm                                                           2309
      oma=1.0-alpha                                                        2309
      nlm=0                                                                2309
      ixx=0                                                                2309
      al=0.0                                                               2310
      dq=d*q                                                               2310
      call died(no,nk,dq,kp,jp,dk)                                         2311
      a=0.0                                                                2311
      f(1)=0.0                                                             2311
      fmax=log(huge(f(1))*0.1)                                             2312
      if(nonzero(no,g) .eq. 0)goto 15831                                   2312
      f=g-dot_product(q,g)                                                 2313
      e=q*exp(sign(min(abs(f),fmax),f))                                    2314
      goto 15841                                                           2315
15831 continue                                                             2315
      f=0.0                                                                2315
      e=q                                                                  2315
15841 continue                                                             2316
15821 continue                                                             2316
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 2317
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         2317
      dev0=rr                                                              2318
15850 do 15851 i=1,no                                                      2318
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 15871                   2318
      w(i)=0.0                                                             2318
      wr(i)=w(i)                                                           2318
15871 continue                                                             2318
15851 continue                                                             2319
15852 continue                                                             2319
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2320
      if(jerr.ne.0) go to 11790                                            2321
      if(flmin .ge. 1.0)goto 15891                                         2321
      eqs=max(eps,flmin)                                                   2321
      alf=eqs**(1.0/(nlam-1))                                              2321
15891 continue                                                             2322
      m=0                                                                  2322
      mm=0                                                                 2322
      nlp=0                                                                2322
      nin=nlp                                                              2322
      mnl=min(mnlam,nlam)                                                  2322
      as=0.0                                                               2322
      cthr=cthri*dev0                                                      2323
15900 do 15901 j=1,ni                                                      2323
      if(ju(j).eq.0)goto 15901                                             2323
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2323
15901 continue                                                             2324
15902 continue                                                             2324
15910 do 15911 ilm=1,nlam                                                  2324
      al0=al                                                               2325
      if(flmin .lt. 1.0)goto 15931                                         2325
      al=ulam(ilm)                                                         2325
      goto 15921                                                           2326
15931 if(ilm .le. 2)goto 15941                                             2326
      al=al*alf                                                            2326
      goto 15921                                                           2327
15941 if(ilm .ne. 1)goto 15951                                             2327
      al=big                                                               2327
      goto 15961                                                           2328
15951 continue                                                             2328
      al0=0.0                                                              2329
15970 do 15971 j=1,ni                                                      2329
      if(ju(j).eq.0)goto 15971                                             2329
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2329
15971 continue                                                             2330
15972 continue                                                             2330
      al0=al0/max(parm,1.0e-3)                                             2330
      al=alf*al0                                                           2331
15961 continue                                                             2332
15921 continue                                                             2332
      sa=alpha*al                                                          2332
      omal=oma*al                                                          2332
      tlam=alpha*(2.0*al-al0)                                              2333
15980 do 15981 k=1,ni                                                      2333
      if(ixx(k).eq.1)goto 15981                                            2333
      if(ju(k).eq.0)goto 15981                                             2334
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2335
15981 continue                                                             2336
15982 continue                                                             2336
10680 continue                                                             2337
15990 continue                                                             2337
15991 continue                                                             2337
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2338
      call vars(no,ni,x,w,ixx,v)                                           2339
16000 continue                                                             2339
16001 continue                                                             2339
      nlp=nlp+1                                                            2339
      dli=0.0                                                              2340
16010 do 16011 j=1,ni                                                      2340
      if(ixx(j).eq.0)goto 16011                                            2341
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2342
      if(abs(u) .gt. vp(j)*sa)goto 16031                                   2342
      at=0.0                                                               2342
      goto 16041                                                           2343
16031 continue                                                             2343
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2343
16041 continue                                                             2344
16021 continue                                                             2344
      if(at .eq. a(j))goto 16061                                           2344
      del=at-a(j)                                                          2344
      a(j)=at                                                              2344
      dli=max(dli,v(j)*del**2)                                             2345
      wr=wr-del*w*x(:,j)                                                   2345
      f=f+del*x(:,j)                                                       2346
      if(mm(j) .ne. 0)goto 16081                                           2346
      nin=nin+1                                                            2346
      if(nin.gt.nx)goto 16012                                              2347
      mm(j)=nin                                                            2347
      m(nin)=j                                                             2348
16081 continue                                                             2349
16061 continue                                                             2350
16011 continue                                                             2351
16012 continue                                                             2351
      if(nin.gt.nx)goto 16002                                              2351
      if(dli.lt.cthr)goto 16002                                            2352
      if(nlp .le. maxit)goto 16101                                         2352
      jerr=-ilm                                                            2352
      return                                                               2352
16101 continue                                                             2353
16110 continue                                                             2353
16111 continue                                                             2353
      nlp=nlp+1                                                            2353
      dli=0.0                                                              2354
16120 do 16121 l=1,nin                                                     2354
      j=m(l)                                                               2355
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2356
      if(abs(u) .gt. vp(j)*sa)goto 16141                                   2356
      at=0.0                                                               2356
      goto 16151                                                           2357
16141 continue                                                             2357
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2357
16151 continue                                                             2358
16131 continue                                                             2358
      if(at .eq. a(j))goto 16171                                           2358
      del=at-a(j)                                                          2358
      a(j)=at                                                              2358
      dli=max(dli,v(j)*del**2)                                             2359
      wr=wr-del*w*x(:,j)                                                   2359
      f=f+del*x(:,j)                                                       2360
16171 continue                                                             2361
16121 continue                                                             2362
16122 continue                                                             2362
      if(dli.lt.cthr)goto 16112                                            2362
      if(nlp .le. maxit)goto 16191                                         2362
      jerr=-ilm                                                            2362
      return                                                               2362
16191 continue                                                             2363
      goto 16111                                                           2364
16112 continue                                                             2364
      goto 16001                                                           2365
16002 continue                                                             2365
      if(nin.gt.nx)goto 15992                                              2366
      e=q*exp(sign(min(abs(f),fmax),f))                                    2367
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2368
      if(jerr .eq. 0)goto 16211                                            2368
      jerr=jerr-ilm                                                        2368
      go to 11790                                                          2368
16211 continue                                                             2369
      ix=0                                                                 2370
16220 do 16221 j=1,nin                                                     2370
      k=m(j)                                                               2371
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 16221                           2371
      ix=1                                                                 2371
      goto 16222                                                           2371
16221 continue                                                             2372
16222 continue                                                             2372
      if(ix .ne. 0)goto 16241                                              2373
16250 do 16251 k=1,ni                                                      2373
      if(ixx(k).eq.1)goto 16251                                            2373
      if(ju(k).eq.0)goto 16251                                             2374
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2375
      if(ga(k) .le. sa*vp(k))goto 16271                                    2375
      ixx(k)=1                                                             2375
      ix=1                                                                 2375
16271 continue                                                             2376
16251 continue                                                             2377
16252 continue                                                             2377
      if(ix.eq.1) go to 10680                                              2378
      goto 15992                                                           2379
16241 continue                                                             2380
      goto 15991                                                           2381
15992 continue                                                             2381
      if(nin .le. nx)goto 16291                                            2381
      jerr=-10000-ilm                                                      2381
      goto 15912                                                           2381
16291 continue                                                             2382
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               2382
      kin(ilm)=nin                                                         2383
      alm(ilm)=al                                                          2383
      lmu=ilm                                                              2384
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   2385
      if(ilm.lt.mnl)goto 15911                                             2385
      if(flmin.ge.1.0)goto 15911                                           2386
      me=0                                                                 2386
16300 do 16301 j=1,nin                                                     2386
      if(ao(j,ilm).ne.0.0) me=me+1                                         2386
16301 continue                                                             2386
16302 continue                                                             2386
      if(me.gt.ne)goto 15912                                               2387
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 15912              2388
      if(dev(ilm).gt.devmax)goto 15912                                     2389
15911 continue                                                             2390
15912 continue                                                             2390
      g=f                                                                  2391
11790 continue                                                             2391
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)              2392
      return                                                               2393
      end                                                                  2394
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 2395
      real ca(nin),x(n,*),f(n)                                             2395
      integer ia(nin)                                                      2396
      f=0.0                                                                2396
      if(nin.le.0) return                                                  2397
16310 do 16311 i=1,n                                                       2397
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      2397
16311 continue                                                             2398
16312 continue                                                             2398
      return                                                               2399
      end                                                                  2400
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         2401
      real y(no),d(no),q(no)                                               2401
      integer jp(no),kp(*)                                                 2402
16320 do 16321 j=1,no                                                      2402
      jp(j)=j                                                              2402
16321 continue                                                             2402
16322 continue                                                             2402
      call psort7(y,jp,1,no)                                               2403
      nj=0                                                                 2403
16330 do 16331 j=1,no                                                      2403
      if(q(jp(j)).le.0.0)goto 16331                                        2403
      nj=nj+1                                                              2403
      jp(nj)=jp(j)                                                         2403
16331 continue                                                             2404
16332 continue                                                             2404
      if(nj .ne. 0)goto 16351                                              2404
      jerr=20000                                                           2404
      return                                                               2404
16351 continue                                                             2405
      j=1                                                                  2405
16360 continue                                                             2405
16361 if(d(jp(j)).gt.0.0)goto 16362                                        2405
      j=j+1                                                                2405
      if(j.gt.nj)goto 16362                                                2405
      goto 16361                                                           2406
16362 continue                                                             2406
      if(j .lt. nj-1)goto 16381                                            2406
      jerr=30000                                                           2406
      return                                                               2406
16381 continue                                                             2407
      j0=j-1                                                               2407
      nj=nj-j0                                                             2407
16390 do 16391 j=1,nj                                                      2407
      jp(j)=jp(j+j0)                                                       2407
16391 continue                                                             2408
16392 continue                                                             2408
      jerr=0                                                               2408
      nk=0                                                                 2408
      t0=y(jp(1))                                                          2408
      yk=t0                                                                2408
      j=2                                                                  2409
16400 continue                                                             2409
16401 continue                                                             2409
16410 continue                                                             2410
16411 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 16412                     2410
      j=j+1                                                                2410
      if(j.gt.nj)goto 16412                                                2410
      goto 16411                                                           2411
16412 continue                                                             2411
      nk=nk+1                                                              2411
      kp(nk)=j-1                                                           2411
      if(j.gt.nj)goto 16402                                                2412
      if(j .ne. nj)goto 16431                                              2412
      nk=nk+1                                                              2412
      kp(nk)=nj                                                            2412
      goto 16402                                                           2412
16431 continue                                                             2413
      yk=y(jp(j))                                                          2413
      j=j+1                                                                2414
      goto 16401                                                           2415
16402 continue                                                             2415
      return                                                               2416
      end                                                                  2417
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     2418
      real d(no),dk(nk),wr(no),w(no)                                       2419
      real e(no),u(no),b,c                                                 2419
      integer kp(nk),jp(no)                                                2420
      call usk(no,nk,kp,jp,e,u)                                            2421
      b=dk(1)/u(1)                                                         2421
      c=dk(1)/u(1)**2                                                      2421
      jerr=0                                                               2422
16440 do 16441 j=1,kp(1)                                                   2422
      i=jp(j)                                                              2423
      w(i)=e(i)*(b-e(i)*c)                                                 2423
      if(w(i) .gt. 0.0)goto 16461                                          2423
      jerr=-30000                                                          2423
      return                                                               2423
16461 continue                                                             2424
      wr(i)=d(i)-e(i)*b                                                    2425
16441 continue                                                             2426
16442 continue                                                             2426
16470 do 16471 k=2,nk                                                      2426
      j1=kp(k-1)+1                                                         2426
      j2=kp(k)                                                             2427
      b=b+dk(k)/u(k)                                                       2427
      c=c+dk(k)/u(k)**2                                                    2428
16480 do 16481 j=j1,j2                                                     2428
      i=jp(j)                                                              2429
      w(i)=e(i)*(b-e(i)*c)                                                 2429
      if(w(i) .gt. 0.0)goto 16501                                          2429
      jerr=-30000                                                          2429
      return                                                               2429
16501 continue                                                             2430
      wr(i)=d(i)-e(i)*b                                                    2431
16481 continue                                                             2432
16482 continue                                                             2432
16471 continue                                                             2433
16472 continue                                                             2433
      return                                                               2434
      end                                                                  2435
      subroutine vars(no,ni,x,w,ixx,v)                                     2436
      real x(no,ni),w(no),v(ni)                                            2436
      integer ixx(ni)                                                      2437
16510 do 16511 j=1,ni                                                      2437
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                        2437
16511 continue                                                             2438
16512 continue                                                             2438
      return                                                               2439
      end                                                                  2440
      subroutine died(no,nk,d,kp,jp,dk)                                    2441
      real d(no),dk(nk)                                                    2441
      integer kp(nk),jp(no)                                                2442
      dk(1)=sum(d(jp(1:kp(1))))                                            2443
16520 do 16521 k=2,nk                                                      2443
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  2443
16521 continue                                                             2444
16522 continue                                                             2444
      return                                                               2445
      end                                                                  2446
      subroutine usk(no,nk,kp,jp,e,u)                                      2447
      real e(no),u(nk),h                                                   2447
      integer kp(nk),jp(no)                                                2448
      h=0.0                                                                2449
16530 do 16531 k=nk,1,-1                                                   2449
      j2=kp(k)                                                             2450
      j1=1                                                                 2450
      if(k.gt.1) j1=kp(k-1)+1                                              2451
16540 do 16541 j=j2,j1,-1                                                  2451
      h=h+e(jp(j))                                                         2451
16541 continue                                                             2452
16542 continue                                                             2452
      u(k)=h                                                               2453
16531 continue                                                             2454
16532 continue                                                             2454
      return                                                               2455
      end                                                                  2456
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             2457
      real d(no),dk(nk),f(no)                                              2458
      integer kp(nk),jp(no)                                                2458
      real e(no),u(nk),s                                                   2459
      call usk(no,nk,kp,jp,e,u)                                            2459
      u=log(u)                                                             2460
      risk=dot_product(d,f)-dot_product(dk,u)                              2461
      return                                                               2462
      end                                                                  2463
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 2464
      real x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(nlam)          2465
      real, dimension (:), allocatable :: dk,f,xm,dq,q                          
      real, dimension (:), allocatable :: e,uu                                  
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          2471
      allocate(q(1:no),stat=ierr)                                          2471
      jerr=jerr+ierr                                                       2472
      allocate(uu(1:no),stat=ierr)                                         2472
      jerr=jerr+ierr                                                       2473
      allocate(f(1:no),stat=ierr)                                          2473
      jerr=jerr+ierr                                                       2474
      allocate(dk(1:no),stat=ierr)                                         2474
      jerr=jerr+ierr                                                       2475
      allocate(jp(1:no),stat=ierr)                                         2475
      jerr=jerr+ierr                                                       2476
      allocate(kp(1:no),stat=ierr)                                         2476
      jerr=jerr+ierr                                                       2477
      allocate(dq(1:no),stat=ierr)                                         2477
      jerr=jerr+ierr                                                       2478
      allocate(xm(1:ni),stat=ierr)                                         2478
      jerr=jerr+ierr                                                       2479
      if(jerr.ne.0) go to 11790                                            2480
      q=max(0.0,w)                                                         2480
      sw=sum(q)                                                            2481
      if(sw .gt. 0.0)goto 16561                                            2481
      jerr=9999                                                            2481
      go to 11790                                                          2481
16561 continue                                                             2482
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2483
      if(jerr.ne.0) go to 11790                                            2483
      fmax=log(huge(e(1))*0.1)                                             2484
      dq=d*q                                                               2484
      call died(no,nk,dq,kp,jp,dk)                                         2484
      gm=dot_product(q,g)/sw                                               2485
16570 do 16571 j=1,ni                                                      2485
      xm(j)=dot_product(q,x(:,j))/sw                                       2485
16571 continue                                                             2486
16572 continue                                                             2486
16580 do 16581 lam=1,nlam                                                  2487
16590 do 16591 i=1,no                                                      2487
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       2488
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        2489
16591 continue                                                             2490
16592 continue                                                             2490
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          2491
16581 continue                                                             2492
16582 continue                                                             2492
11790 continue                                                             2492
      deallocate(e,uu,dk,f,jp,kp,dq)                                       2493
      return                                                               2494
      end                                                                  2495
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,ne,nx,nlam,flmin,ulam   2497 
     *,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    2498
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2499
      integer jd(*),ia(nx),nin(nlam)                                       2500
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16611                                    2504
      jerr=10000                                                           2504
      return                                                               2504
16611 continue                                                             2505
      if(minval(y) .ge. 0.0)goto 16631                                     2505
      jerr=8888                                                            2505
      return                                                               2505
16631 continue                                                             2506
      allocate(ww(1:no),stat=jerr)                                         2507
      allocate(ju(1:ni),stat=ierr)                                         2507
      jerr=jerr+ierr                                                       2508
      allocate(vq(1:ni),stat=ierr)                                         2508
      jerr=jerr+ierr                                                       2509
      allocate(xm(1:ni),stat=ierr)                                         2509
      jerr=jerr+ierr                                                       2510
      if(isd .le. 0)goto 16651                                             2510
      allocate(xs(1:ni),stat=ierr)                                         2510
      jerr=jerr+ierr                                                       2510
16651 continue                                                             2511
      if(jerr.ne.0) return                                                 2512
      call chkvars(no,ni,x,ju)                                             2513
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2514
      if(maxval(ju) .gt. 0)goto 16671                                      2514
      jerr=7777                                                            2514
      go to 11790                                                          2514
16671 continue                                                             2515
      vq=max(0.0,vp)                                                       2515
      vq=vq*ni/sum(vq)                                                     2516
      ww=max(0.0,w)                                                        2516
      sw=sum(ww)                                                           2516
      if(sw .gt. 0.0)goto 16691                                            2516
      jerr=9999                                                            2516
      go to 11790                                                          2516
16691 continue                                                             2517
      ww=ww/sw                                                             2518
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             2519
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr,   2521 
     *  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11790                                            2521
      dev0=2.0*sw*dev0                                                     2522
16700 do 16701 k=1,lmu                                                     2522
      nk=nin(k)                                                            2523
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2524
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2525
16701 continue                                                             2526
16702 continue                                                             2526
11790 continue                                                             2526
      deallocate(ww,ju,vq,xm)                                              2526
      if(isd.gt.0) deallocate(xs)                                          2527
      return                                                               2528
      end                                                                  2529
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,ne,nx,nlam,flmin,ulam   2531 
     *,shri,  isd,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2532 
     *9)
      real x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)                    2533
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2534
      integer ju(ni),m(nx),kin(nlam)                                       2535
      real, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga                    
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(a(1:ni),stat=jerr)                                          2540
      allocate(as(1:ni),stat=ierr)                                         2540
      jerr=jerr+ierr                                                       2541
      allocate(t(1:no),stat=ierr)                                          2541
      jerr=jerr+ierr                                                       2542
      allocate(mm(1:ni),stat=ierr)                                         2542
      jerr=jerr+ierr                                                       2543
      allocate(ga(1:ni),stat=ierr)                                         2543
      jerr=jerr+ierr                                                       2544
      allocate(ixx(1:ni),stat=ierr)                                        2544
      jerr=jerr+ierr                                                       2545
      allocate(wr(1:no),stat=ierr)                                         2545
      jerr=jerr+ierr                                                       2546
      allocate(v(1:ni),stat=ierr)                                          2546
      jerr=jerr+ierr                                                       2547
      allocate(w(1:no),stat=ierr)                                          2547
      jerr=jerr+ierr                                                       2548
      allocate(f(1:no),stat=ierr)                                          2548
      jerr=jerr+ierr                                                       2549
      if(jerr.ne.0) return                                                 2550
      bta=parm                                                             2550
      omb=1.0-bta                                                          2551
      t=q*y                                                                2551
      yb=sum(t)                                                            2551
      fmax=log(huge(bta)*0.1)                                              2552
      if(nonzero(no,g) .ne. 0)goto 16721                                   2552
      w=q*yb                                                               2552
      az=log(yb)                                                           2552
      f=az                                                                 2552
      dv0=yb*(log(yb)-1.0)                                                 2552
      goto 16731                                                           2553
16721 continue                                                             2553
      w=q*exp(sign(min(abs(g),fmax),g))                                    2553
      v0=sum(w)                                                            2553
      eaz=yb/v0                                                            2554
      w=eaz*w                                                              2554
      az=log(eaz)                                                          2554
      f=az+g                                                               2555
      dv0=dot_product(t,g)-yb*(1.0-az)                                     2556
16731 continue                                                             2557
16711 continue                                                             2557
      a=0.0                                                                2557
      as=0.0                                                               2557
      wr=t-w                                                               2557
      v0=yb                                                                2557
      dvr=-yb                                                              2558
16740 do 16741 i=1,no                                                      2558
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               2558
16741 continue                                                             2558
16742 continue                                                             2558
      dvr=dvr-dv0                                                          2558
      dev0=dvr                                                             2559
      if(flmin .ge. 1.0)goto 16761                                         2559
      eqs=max(eps,flmin)                                                   2559
      alf=eqs**(1.0/(nlam-1))                                              2559
16761 continue                                                             2560
      m=0                                                                  2560
      mm=0                                                                 2560
      nlp=0                                                                2560
      nin=nlp                                                              2560
      mnl=min(mnlam,nlam)                                                  2560
      shr=shri*dev0                                                        2560
      ixx=0                                                                2560
      al=0.0                                                               2561
16770 do 16771 j=1,ni                                                      2561
      if(ju(j).eq.0)goto 16771                                             2561
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2561
16771 continue                                                             2562
16772 continue                                                             2562
16780 do 16781 ilm=1,nlam                                                  2562
      al0=al                                                               2563
      if(flmin .lt. 1.0)goto 16801                                         2563
      al=ulam(ilm)                                                         2563
      goto 16791                                                           2564
16801 if(ilm .le. 2)goto 16811                                             2564
      al=al*alf                                                            2564
      goto 16791                                                           2565
16811 if(ilm .ne. 1)goto 16821                                             2565
      al=big                                                               2565
      goto 16831                                                           2566
16821 continue                                                             2566
      al0=0.0                                                              2567
16840 do 16841 j=1,ni                                                      2567
      if(ju(j).eq.0)goto 16841                                             2567
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2567
16841 continue                                                             2568
16842 continue                                                             2568
      al0=al0/max(bta,1.0e-3)                                              2568
      al=alf*al0                                                           2569
16831 continue                                                             2570
16791 continue                                                             2570
      al2=al*omb                                                           2570
      al1=al*bta                                                           2570
      tlam=bta*(2.0*al-al0)                                                2571
16850 do 16851 k=1,ni                                                      2571
      if(ixx(k).eq.1)goto 16851                                            2571
      if(ju(k).eq.0)goto 16851                                             2572
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2573
16851 continue                                                             2574
16852 continue                                                             2574
10680 continue                                                             2575
16860 continue                                                             2575
16861 continue                                                             2575
      az0=az                                                               2576
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2577
16870 do 16871 j=1,ni                                                      2577
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                        2577
16871 continue                                                             2578
16872 continue                                                             2578
16880 continue                                                             2578
16881 continue                                                             2578
      nlp=nlp+1                                                            2578
      dlx=0.0                                                              2579
16890 do 16891 k=1,ni                                                      2579
      if(ixx(k).eq.0)goto 16891                                            2579
      ak=a(k)                                                              2580
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2580
      au=abs(u)-vp(k)*al1                                                  2581
      if(au .gt. 0.0)goto 16911                                            2581
      a(k)=0.0                                                             2581
      goto 16921                                                           2582
16911 continue                                                             2582
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2582
16921 continue                                                             2583
16901 continue                                                             2583
      if(a(k).eq.ak)goto 16891                                             2583
      d=a(k)-ak                                                            2583
      dlx=max(dlx,v(k)*d**2)                                               2584
      wr=wr-d*w*x(:,k)                                                     2584
      f=f+d*x(:,k)                                                         2585
      if(mm(k) .ne. 0)goto 16941                                           2585
      nin=nin+1                                                            2585
      if(nin.gt.nx)goto 16892                                              2586
      mm(k)=nin                                                            2586
      m(nin)=k                                                             2587
16941 continue                                                             2588
16891 continue                                                             2589
16892 continue                                                             2589
      if(nin.gt.nx)goto 16882                                              2589
      d=sum(wr)/v0                                                         2590
      az=az+d                                                              2590
      dlx=max(dlx,v0*d**2)                                                 2590
      wr=wr-d*w                                                            2590
      f=f+d                                                                2591
      if(dlx.lt.shr)goto 16882                                             2591
      if(nlp .le. maxit)goto 16961                                         2591
      jerr=-ilm                                                            2591
      return                                                               2591
16961 continue                                                             2592
16970 continue                                                             2592
16971 continue                                                             2592
      nlp=nlp+1                                                            2592
      dlx=0.0                                                              2593
16980 do 16981 l=1,nin                                                     2593
      k=m(l)                                                               2593
      ak=a(k)                                                              2594
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2594
      au=abs(u)-vp(k)*al1                                                  2595
      if(au .gt. 0.0)goto 17001                                            2595
      a(k)=0.0                                                             2595
      goto 17011                                                           2596
17001 continue                                                             2596
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2596
17011 continue                                                             2597
16991 continue                                                             2597
      if(a(k).eq.ak)goto 16981                                             2597
      d=a(k)-ak                                                            2597
      dlx=max(dlx,v(k)*d**2)                                               2598
      wr=wr-d*w*x(:,k)                                                     2598
      f=f+d*x(:,k)                                                         2600
16981 continue                                                             2600
16982 continue                                                             2600
      d=sum(wr)/v0                                                         2600
      az=az+d                                                              2600
      dlx=max(dlx,v0*d**2)                                                 2600
      wr=wr-d*w                                                            2600
      f=f+d                                                                2601
      if(dlx.lt.shr)goto 16972                                             2601
      if(nlp .le. maxit)goto 17031                                         2601
      jerr=-ilm                                                            2601
      return                                                               2601
17031 continue                                                             2602
      goto 16971                                                           2603
16972 continue                                                             2603
      goto 16881                                                           2604
16882 continue                                                             2604
      if(nin.gt.nx)goto 16862                                              2605
      w=q*exp(sign(min(abs(f),fmax),f))                                    2605
      v0=sum(w)                                                            2605
      wr=t-w                                                               2606
      if(v0*(az-az0)**2 .ge. shr)goto 17051                                2606
      ix=0                                                                 2607
17060 do 17061 j=1,nin                                                     2607
      k=m(j)                                                               2608
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 17061                            2608
      ix=1                                                                 2608
      goto 17062                                                           2609
17061 continue                                                             2610
17062 continue                                                             2610
      if(ix .ne. 0)goto 17081                                              2611
17090 do 17091 k=1,ni                                                      2611
      if(ixx(k).eq.1)goto 17091                                            2611
      if(ju(k).eq.0)goto 17091                                             2612
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2613
      if(ga(k) .le. al1*vp(k))goto 17111                                   2613
      ixx(k)=1                                                             2613
      ix=1                                                                 2613
17111 continue                                                             2614
17091 continue                                                             2615
17092 continue                                                             2615
      if(ix.eq.1) go to 10680                                              2616
      goto 16862                                                           2617
17081 continue                                                             2618
17051 continue                                                             2619
      goto 16861                                                           2620
16862 continue                                                             2620
      if(nin .le. nx)goto 17131                                            2620
      jerr=-10000-ilm                                                      2620
      goto 16782                                                           2620
17131 continue                                                             2621
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2621
      kin(ilm)=nin                                                         2622
      a0(ilm)=az                                                           2622
      alm(ilm)=al                                                          2622
      lmu=ilm                                                              2623
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               2624
      if(ilm.lt.mnl)goto 16781                                             2624
      if(flmin.ge.1.0)goto 16781                                           2625
      me=0                                                                 2625
17140 do 17141 j=1,nin                                                     2625
      if(ca(j,ilm).ne.0.0) me=me+1                                         2625
17141 continue                                                             2625
17142 continue                                                             2625
      if(me.gt.ne)goto 16782                                               2626
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 16782              2627
      if(dev(ilm).gt.devmax)goto 16782                                     2628
16781 continue                                                             2629
16782 continue                                                             2629
      g=f                                                                  2630
11790 continue                                                             2630
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                                2631
      return                                                               2632
      end                                                                  2633
      function nonzero(n,v)                                                2634
      real v(n)                                                            2635
      nonzero=0                                                            2635
17150 do 17151 i=1,n                                                       2635
      if(v(i) .eq. 0.0)goto 17171                                          2635
      nonzero=1                                                            2635
      return                                                               2635
17171 continue                                                             2635
17151 continue                                                             2636
17152 continue                                                             2636
      return                                                               2637
      end                                                                  2638
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               2639
      real a(nx,lmu),b(ni,lmu)                                             2639
      integer ia(nx),nin(lmu)                                              2640
17180 do 17181 lam=1,lmu                                                   2640
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        2640
17181 continue                                                             2641
17182 continue                                                             2641
      return                                                               2642
      end                                                                  2643
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           2644
      real a(nx,nc,lmu),b(ni,nc,lmu)                                       2644
      integer ia(nx),nin(lmu)                                              2645
17190 do 17191 lam=1,lmu                                                   2645
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             2645
17191 continue                                                             2646
17192 continue                                                             2646
      return                                                               2647
      end                                                                  2648
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               2649
      real x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)       2650
      real, dimension (:), allocatable :: w                                     
      if(minval(y) .ge. 0.0)goto 17211                                     2653
      jerr=8888                                                            2653
      return                                                               2653
17211 continue                                                             2654
      allocate(w(1:no),stat=jerr)                                          2654
      if(jerr.ne.0) return                                                 2655
      w=max(0.0,q)                                                         2655
      sw=sum(w)                                                            2655
      if(sw .gt. 0.0)goto 17231                                            2655
      jerr=9999                                                            2655
      go to 11790                                                          2655
17231 continue                                                             2656
      yb=dot_product(w,y)/sw                                               2656
      fmax=log(huge(y(1))*0.1)                                             2657
17240 do 17241 lam=1,nlam                                                  2657
      s=0.0                                                                2658
17250 do 17251 i=1,no                                                      2658
      if(w(i).le.0.0)goto 17251                                            2659
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          2660
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      2661
17251 continue                                                             2662
17252 continue                                                             2662
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2663
17241 continue                                                             2664
17242 continue                                                             2664
11790 continue                                                             2664
      deallocate(w)                                                        2665
      return                                                               2666
      end                                                                  2667
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,ne,nx,nlam,fl   2669 
     *min,ulam,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        2670
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2671
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           2672
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 17271                                    2676
      jerr=10000                                                           2676
      return                                                               2676
17271 continue                                                             2677
      if(minval(y) .ge. 0.0)goto 17291                                     2677
      jerr=8888                                                            2677
      return                                                               2677
17291 continue                                                             2678
      allocate(ww(1:no),stat=jerr)                                         2679
      allocate(ju(1:ni),stat=ierr)                                         2679
      jerr=jerr+ierr                                                       2680
      allocate(vq(1:ni),stat=ierr)                                         2680
      jerr=jerr+ierr                                                       2681
      allocate(xm(1:ni),stat=ierr)                                         2681
      jerr=jerr+ierr                                                       2682
      allocate(xs(1:ni),stat=ierr)                                         2682
      jerr=jerr+ierr                                                       2683
      if(jerr.ne.0) return                                                 2684
      call spchkvars(no,ni,x,ix,ju)                                        2685
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2686
      if(maxval(ju) .gt. 0)goto 17311                                      2686
      jerr=7777                                                            2686
      go to 11790                                                          2686
17311 continue                                                             2687
      vq=max(0.0,vp)                                                       2687
      vq=vq*ni/sum(vq)                                                     2688
      ww=max(0.0,w)                                                        2688
      sw=sum(ww)                                                           2688
      if(sw .gt. 0.0)goto 17331                                            2688
      jerr=9999                                                            2688
      go to 11790                                                          2688
17331 continue                                                             2689
      ww=ww/sw                                                             2690
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     2691
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmin,u   2693 
     *lam,thr,  isd,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11790                                            2693
      dev0=2.0*sw*dev0                                                     2694
17340 do 17341 k=1,lmu                                                     2694
      nk=nin(k)                                                            2695
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2696
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2697
17341 continue                                                             2698
17342 continue                                                             2698
11790 continue                                                             2698
      deallocate(ww,ju,vq,xm,xs)                                           2699
      return                                                               2700
      end                                                                  2701
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,ne,nx,nlam,fl   2703 
     *min,ulam,  shri,isd,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,j
     *err)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2704 
     *9)
      real x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),xs(ni)          2705
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2706
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2707
      real, dimension (:), allocatable :: qy,t,w,wr,v,a,as,xm,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(a(1:ni),stat=jerr)                                          2712
      allocate(as(1:ni),stat=ierr)                                         2712
      jerr=jerr+ierr                                                       2713
      allocate(t(1:no),stat=ierr)                                          2713
      jerr=jerr+ierr                                                       2714
      allocate(mm(1:ni),stat=ierr)                                         2714
      jerr=jerr+ierr                                                       2715
      allocate(ga(1:ni),stat=ierr)                                         2715
      jerr=jerr+ierr                                                       2716
      allocate(ixx(1:ni),stat=ierr)                                        2716
      jerr=jerr+ierr                                                       2717
      allocate(wr(1:no),stat=ierr)                                         2717
      jerr=jerr+ierr                                                       2718
      allocate(v(1:ni),stat=ierr)                                          2718
      jerr=jerr+ierr                                                       2719
      allocate(xm(1:ni),stat=ierr)                                         2719
      jerr=jerr+ierr                                                       2720
      allocate(w(1:no),stat=ierr)                                          2720
      jerr=jerr+ierr                                                       2721
      allocate(qy(1:no),stat=ierr)                                         2721
      jerr=jerr+ierr                                                       2722
      if(jerr.ne.0) return                                                 2723
      bta=parm                                                             2723
      omb=1.0-bta                                                          2723
      fmax=log(huge(bta)*0.1)                                              2724
      qy=q*y                                                               2724
      yb=sum(qy)                                                           2725
      if(nonzero(no,g) .ne. 0)goto 17361                                   2725
      w=q*yb                                                               2725
      az=log(yb)                                                           2725
      uu=az                                                                2726
      xm=yb*xb                                                             2726
      t=0.0                                                                2726
      dv0=yb*(log(yb)-1.0)                                                 2727
      goto 17371                                                           2728
17361 continue                                                             2728
      w=q*exp(sign(min(abs(g),fmax),g))                                    2728
      ww=sum(w)                                                            2728
      eaz=yb/ww                                                            2729
      w=eaz*w                                                              2729
      az=log(eaz)                                                          2729
      uu=az                                                                2729
      t=g                                                                  2729
      dv0=dot_product(qy,g)-yb*(1.0-az)                                    2730
17380 do 17381 j=1,ni                                                      2730
      if(ju(j).eq.0)goto 17381                                             2730
      jb=ix(j)                                                             2730
      je=ix(j+1)-1                                                         2731
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2732
17381 continue                                                             2733
17382 continue                                                             2733
17371 continue                                                             2734
17351 continue                                                             2734
      tt=yb*uu                                                             2734
      ww=yb                                                                2734
      wr=qy-q*(yb*(1.0-uu))                                                2734
      a=0.0                                                                2734
      as=0.0                                                               2735
      dvr=-yb                                                              2736
17390 do 17391 i=1,no                                                      2736
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             2736
17391 continue                                                             2736
17392 continue                                                             2736
      dvr=dvr-dv0                                                          2736
      dev0=dvr                                                             2737
      if(flmin .ge. 1.0)goto 17411                                         2737
      eqs=max(eps,flmin)                                                   2737
      alf=eqs**(1.0/(nlam-1))                                              2737
17411 continue                                                             2738
      m=0                                                                  2738
      mm=0                                                                 2738
      nlp=0                                                                2738
      nin=nlp                                                              2738
      mnl=min(mnlam,nlam)                                                  2738
      shr=shri*dev0                                                        2738
      al=0.0                                                               2738
      ixx=0                                                                2739
17420 do 17421 j=1,ni                                                      2739
      if(ju(j).eq.0)goto 17421                                             2740
      jb=ix(j)                                                             2740
      je=ix(j+1)-1                                                         2741
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   2743 
     *)-xb(j)*tt)/xs(j)
17421 continue                                                             2744
17422 continue                                                             2744
17430 do 17431 ilm=1,nlam                                                  2744
      al0=al                                                               2745
      if(flmin .lt. 1.0)goto 17451                                         2745
      al=ulam(ilm)                                                         2745
      goto 17441                                                           2746
17451 if(ilm .le. 2)goto 17461                                             2746
      al=al*alf                                                            2746
      goto 17441                                                           2747
17461 if(ilm .ne. 1)goto 17471                                             2747
      al=big                                                               2747
      goto 17481                                                           2748
17471 continue                                                             2748
      al0=0.0                                                              2749
17490 do 17491 j=1,ni                                                      2749
      if(ju(j).eq.0)goto 17491                                             2749
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2749
17491 continue                                                             2750
17492 continue                                                             2750
      al0=al0/max(bta,1.0e-3)                                              2750
      al=alf*al0                                                           2751
17481 continue                                                             2752
17441 continue                                                             2752
      al2=al*omb                                                           2752
      al1=al*bta                                                           2752
      tlam=bta*(2.0*al-al0)                                                2753
17500 do 17501 k=1,ni                                                      2753
      if(ixx(k).eq.1)goto 17501                                            2753
      if(ju(k).eq.0)goto 17501                                             2754
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2755
17501 continue                                                             2756
17502 continue                                                             2756
10680 continue                                                             2757
17510 continue                                                             2757
17511 continue                                                             2757
      az0=az                                                               2758
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2759
17520 do 17521 j=1,ni                                                      2759
      if(ixx(j).eq.0)goto 17521                                            2759
      jb=ix(j)                                                             2759
      je=ix(j+1)-1                                                         2760
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2761
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   2763 
     *b(j)**2)/xs(j)**2
17521 continue                                                             2764
17522 continue                                                             2764
17530 continue                                                             2764
17531 continue                                                             2764
      nlp=nlp+1                                                            2765
      dlx=0.0                                                              2766
17540 do 17541 k=1,ni                                                      2766
      if(ixx(k).eq.0)goto 17541                                            2766
      jb=ix(k)                                                             2766
      je=ix(k+1)-1                                                         2766
      ak=a(k)                                                              2767
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2769 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2770
      if(au .gt. 0.0)goto 17561                                            2770
      a(k)=0.0                                                             2770
      goto 17571                                                           2771
17561 continue                                                             2771
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2771
17571 continue                                                             2772
17551 continue                                                             2772
      if(a(k).eq.ak)goto 17541                                             2773
      if(mm(k) .ne. 0)goto 17591                                           2773
      nin=nin+1                                                            2773
      if(nin.gt.nx)goto 17542                                              2774
      mm(k)=nin                                                            2774
      m(nin)=k                                                             2775
17591 continue                                                             2776
      d=a(k)-ak                                                            2776
      dlx=max(dlx,v(k)*d**2)                                               2776
      dv=d/xs(k)                                                           2777
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2778
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2779
      uu=uu-dv*xb(k)                                                       2779
      tt=tt-dv*xm(k)                                                       2780
17541 continue                                                             2781
17542 continue                                                             2781
      if(nin.gt.nx)goto 17532                                              2781
      d=tt/ww-uu                                                           2782
      az=az+d                                                              2782
      dlx=max(dlx,ww*d**2)                                                 2782
      uu=uu+d                                                              2783
      if(dlx.lt.shr)goto 17532                                             2783
      if(nlp .le. maxit)goto 17611                                         2783
      jerr=-ilm                                                            2783
      return                                                               2783
17611 continue                                                             2784
17620 continue                                                             2784
17621 continue                                                             2784
      nlp=nlp+1                                                            2784
      dlx=0.0                                                              2785
17630 do 17631 l=1,nin                                                     2785
      k=m(l)                                                               2786
      jb=ix(k)                                                             2786
      je=ix(k+1)-1                                                         2786
      ak=a(k)                                                              2787
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2789 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2790
      if(au .gt. 0.0)goto 17651                                            2790
      a(k)=0.0                                                             2790
      goto 17661                                                           2791
17651 continue                                                             2791
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2791
17661 continue                                                             2792
17641 continue                                                             2792
      if(a(k).eq.ak)goto 17631                                             2792
      d=a(k)-ak                                                            2792
      dlx=max(dlx,v(k)*d**2)                                               2793
      dv=d/xs(k)                                                           2793
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2794
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2795
      uu=uu-dv*xb(k)                                                       2795
      tt=tt-dv*xm(k)                                                       2796
17631 continue                                                             2797
17632 continue                                                             2797
      d=tt/ww-uu                                                           2797
      az=az+d                                                              2797
      dlx=max(dlx,ww*d**2)                                                 2797
      uu=uu+d                                                              2798
      if(dlx.lt.shr)goto 17622                                             2798
      if(nlp .le. maxit)goto 17681                                         2798
      jerr=-ilm                                                            2798
      return                                                               2798
17681 continue                                                             2799
      goto 17621                                                           2800
17622 continue                                                             2800
      goto 17531                                                           2801
17532 continue                                                             2801
      if(nin.gt.nx)goto 17512                                              2802
      euu=exp(sign(min(abs(uu),fmax),uu))                                  2803
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                2803
      ww=sum(w)                                                            2804
      wr=qy-w*(1.0-uu)                                                     2804
      tt=sum(wr)                                                           2805
      if(ww*(az-az0)**2 .ge. shr)goto 17701                                2805
      kx=0                                                                 2806
17710 do 17711 j=1,nin                                                     2806
      k=m(j)                                                               2807
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 17711                            2807
      kx=1                                                                 2807
      goto 17712                                                           2808
17711 continue                                                             2809
17712 continue                                                             2809
      if(kx .ne. 0)goto 17731                                              2810
17740 do 17741 j=1,ni                                                      2810
      if(ixx(j).eq.1)goto 17741                                            2810
      if(ju(j).eq.0)goto 17741                                             2811
      jb=ix(j)                                                             2811
      je=ix(j+1)-1                                                         2812
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2813
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   2815 
     *)-xb(j)*tt)/xs(j)
      if(ga(j) .le. al1*vp(j))goto 17761                                   2815
      ixx(j)=1                                                             2815
      kx=1                                                                 2815
17761 continue                                                             2816
17741 continue                                                             2817
17742 continue                                                             2817
      if(kx.eq.1) go to 10680                                              2818
      goto 17512                                                           2819
17731 continue                                                             2820
17701 continue                                                             2821
      goto 17511                                                           2822
17512 continue                                                             2822
      if(nin .le. nx)goto 17781                                            2822
      jerr=-10000-ilm                                                      2822
      goto 17432                                                           2822
17781 continue                                                             2823
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2823
      kin(ilm)=nin                                                         2824
      a0(ilm)=az                                                           2824
      alm(ilm)=al                                                          2824
      lmu=ilm                                                              2825
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        2826
      if(ilm.lt.mnl)goto 17431                                             2826
      if(flmin.ge.1.0)goto 17431                                           2827
      me=0                                                                 2827
17790 do 17791 j=1,nin                                                     2827
      if(ca(j,ilm).ne.0.0) me=me+1                                         2827
17791 continue                                                             2827
17792 continue                                                             2827
      if(me.gt.ne)goto 17432                                               2828
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17432              2829
      if(dev(ilm).gt.devmax)goto 17432                                     2830
17431 continue                                                             2831
17432 continue                                                             2831
      g=t+uu                                                               2832
11790 continue                                                             2832
      deallocate(t,w,wr,v,a,qy,xm,as,mm,ga,ixx)                            2833
      return                                                               2834
      end                                                                  2835
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       2836
      real x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)           2837
      integer ix(*),jx(*)                                                  2838
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 17811                                     2841
      jerr=8888                                                            2841
      return                                                               2841
17811 continue                                                             2842
      allocate(w(1:no),stat=jerr)                                          2843
      allocate(f(1:no),stat=ierr)                                          2843
      jerr=jerr+ierr                                                       2844
      if(jerr.ne.0) return                                                 2845
      w=max(0.0,q)                                                         2845
      sw=sum(w)                                                            2845
      if(sw .gt. 0.0)goto 17831                                            2845
      jerr=9999                                                            2845
      go to 11790                                                          2845
17831 continue                                                             2846
      yb=dot_product(w,y)/sw                                               2846
      fmax=log(huge(y(1))*0.1)                                             2847
17840 do 17841 lam=1,nlam                                                  2847
      f=a0(lam)                                                            2848
17850 do 17851 j=1,ni                                                      2848
      if(a(j,lam).eq.0.0)goto 17851                                        2848
      jb=ix(j)                                                             2848
      je=ix(j+1)-1                                                         2849
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          2850
17851 continue                                                             2851
17852 continue                                                             2851
      f=f+g                                                                2852
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2853
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2854
17841 continue                                                             2855
17842 continue                                                             2855
11790 continue                                                             2855
      deallocate(w,f)                                                      2856
      return                                                               2857
      end                                                                  2858
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   2859 
     *jerr)
      real x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(nlam)          2860
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 2861
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 17871                                     2864
      jerr=8888                                                            2864
      return                                                               2864
17871 continue                                                             2865
      allocate(w(1:no),stat=jerr)                                          2866
      allocate(f(1:no),stat=ierr)                                          2866
      jerr=jerr+ierr                                                       2867
      if(jerr.ne.0) return                                                 2868
      w=max(0.0,q)                                                         2868
      sw=sum(w)                                                            2868
      if(sw .gt. 0.0)goto 17891                                            2868
      jerr=9999                                                            2868
      go to 11790                                                          2868
17891 continue                                                             2869
      yb=dot_product(w,y)/sw                                               2869
      fmax=log(huge(y(1))*0.1)                                             2870
17900 do 17901 lam=1,nlam                                                  2870
      f=a0(lam)                                                            2871
17910 do 17911 k=1,nin(lam)                                                2871
      j=ia(k)                                                              2871
      jb=ix(j)                                                             2871
      je=ix(j+1)-1                                                         2872
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         2873
17911 continue                                                             2874
17912 continue                                                             2874
      f=f+g                                                                2875
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2876
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2877
17901 continue                                                             2878
17902 continue                                                             2878
11790 continue                                                             2878
      deallocate(w,f)                                                      2879
      return                                                               2880
      end                                                                  2881
      subroutine multelnet  (parm,no,ni,nr,x,y,w,jd,vp,ne,nx,nlam,flmin,   2884 
     *ulam,thr,isd,jsd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no,nr),w(no),vp(ni),ca(nx,nr,nlam)                   2885
      real ulam(nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                      2886
      integer jd(*),ia(nx),nin(nlam)                                       2887
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 17931                                    2890
      jerr=10000                                                           2890
      return                                                               2890
17931 continue                                                             2891
      allocate(vq(1:ni),stat=jerr)                                         2891
      if(jerr.ne.0) return                                                 2892
      vq=max(0.0,vp)                                                       2892
      vq=vq*ni/sum(vq)                                                     2893
      call multelnetn(parm,no,ni,nr,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,th   2895 
     *r,isd,jsd,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      deallocate(vq)                                                       2896
      return                                                               2897
      end                                                                  2898
      subroutine multelnetn (parm,no,ni,nr,x,y,w,jd,vp,ne,nx,nlam,flmin,   2900 
     *ulam,thr,  isd,jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no,nr),w(no),ulam(nlam)                       2901
      real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                  2902
      integer jd(*),ia(nx),nin(nlam)                                       2903
      real, dimension (:), allocatable :: xm,xs,xv,ym,ys                        
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         2908
      allocate(xs(1:ni),stat=ierr)                                         2908
      jerr=jerr+ierr                                                       2909
      allocate(ym(1:nr),stat=ierr)                                         2909
      jerr=jerr+ierr                                                       2910
      allocate(ys(1:nr),stat=ierr)                                         2910
      jerr=jerr+ierr                                                       2911
      allocate(ju(1:ni),stat=ierr)                                         2911
      jerr=jerr+ierr                                                       2912
      allocate(xv(1:ni),stat=ierr)                                         2912
      jerr=jerr+ierr                                                       2913
      if(jerr.ne.0) return                                                 2914
      call chkvars(no,ni,x,ju)                                             2915
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2916
      if(maxval(ju) .gt. 0)goto 17951                                      2916
      jerr=7777                                                            2916
      return                                                               2916
17951 continue                                                             2917
      call multstandard1(no,ni,nr,x,y,w,isd,jsd,ju,xm,xs,ym,ys,xv,ys0,je   2918 
     *rr)
      if(jerr.ne.0) return                                                 2919
      call multelnet2(parm,ni,nr,ju,vp,y,no,ne,nx,x,nlam,flmin,ulam,thr,   2921 
     *maxit,xv,  ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 2922
17960 do 17961 k=1,lmu                                                     2922
      nk=nin(k)                                                            2923
17970 do 17971 j=1,nr                                                      2924
17980 do 17981 l=1,nk                                                      2924
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  2924
17981 continue                                                             2925
17982 continue                                                             2925
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 2926
17971 continue                                                             2927
17972 continue                                                             2927
17961 continue                                                             2928
17962 continue                                                             2928
      deallocate(xm,xs,ym,ys,ju,xv)                                        2929
      return                                                               2930
      end                                                                  2931
      subroutine multstandard1 (no,ni,nr,x,y,w,isd,jsd,ju,xm,xs,ym,ys,xv   2932 
     *,ys0,jerr)
      real x(no,ni),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),ys(nr)      2933
      integer ju(ni)                                                       2934
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                          2937
      if(jerr.ne.0) return                                                 2938
      w=w/sum(w)                                                           2938
      v=sqrt(w)                                                            2939
17990 do 17991 j=1,ni                                                      2939
      if(ju(j).eq.0)goto 17991                                             2940
      xm(j)=dot_product(w,x(:,j))                                          2940
      x(:,j)=v*(x(:,j)-xm(j))                                              2941
      xv(j)=dot_product(x(:,j),x(:,j))                                     2941
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       2942
17991 continue                                                             2943
17992 continue                                                             2943
      if(isd .ne. 0)goto 18011                                             2943
      xs=1.0                                                               2943
      goto 18021                                                           2944
18011 continue                                                             2944
18030 do 18031 j=1,ni                                                      2944
      if(ju(j).eq.0)goto 18031                                             2944
      x(:,j)=x(:,j)/xs(j)                                                  2944
18031 continue                                                             2945
18032 continue                                                             2945
      xv=1.0                                                               2946
18021 continue                                                             2947
18001 continue                                                             2947
      ys0=0.0                                                              2948
18040 do 18041 j=1,nr                                                      2949
      ym(j)=dot_product(w,y(:,j))                                          2949
      y(:,j)=v*(y(:,j)-ym(j))                                              2950
      z=dot_product(y(:,j),y(:,j))                                         2951
      if(jsd .le. 0)goto 18061                                             2951
      ys(j)=sqrt(z)                                                        2951
      y(:,j)=y(:,j)/ys(j)                                                  2951
      goto 18071                                                           2952
18061 continue                                                             2952
      ys0=ys0+z                                                            2952
18071 continue                                                             2953
18051 continue                                                             2953
18041 continue                                                             2954
18042 continue                                                             2954
      if(jsd .ne. 0)goto 18091                                             2954
      ys=1.0                                                               2954
      goto 18101                                                           2954
18091 continue                                                             2954
      ys0=nr                                                               2954
18101 continue                                                             2955
18081 continue                                                             2955
      deallocate(v)                                                        2956
      return                                                               2957
      end                                                                  2958
      subroutine multelnet2(beta,ni,nr,ju,vp,y,no,ne,nx,x,nlam,flmin,ula   2960 
     *m,thri,  maxit,xv,ys0,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99   2961 
     *9)
      real vp(ni),y(no,nr),x(no,ni),ulam(nlam),ao(nx,nr,nlam)              2962
      real rsqo(nlam),almo(nlam),xv(ni)                                    2963
      integer ju(ni),ia(nx),kin(nlam)                                      2964
      real, dimension (:), allocatable :: g,gk,del,gj                           
      integer, dimension (:), allocatable :: mm,ix                              
      real, dimension (:,:), allocatable :: a                                   
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      allocate(gj(1:nr),stat=ierr)                                         2970
      jerr=jerr+ierr                                                       2971
      allocate(gk(1:nr),stat=ierr)                                         2971
      jerr=jerr+ierr                                                       2972
      allocate(del(1:nr),stat=ierr)                                        2972
      jerr=jerr+ierr                                                       2973
      allocate(mm(1:ni),stat=ierr)                                         2973
      jerr=jerr+ierr                                                       2974
      allocate(g(1:ni),stat=ierr)                                          2974
      jerr=jerr+ierr                                                       2975
      allocate(ix(1:ni),stat=ierr)                                         2975
      jerr=jerr+ierr                                                       2976
      if(jerr.ne.0) return                                                 2977
      bta=beta                                                             2977
      omb=1.0-bta                                                          2977
      ix=0                                                                 2977
      thr=thri*ys0/nr                                                      2978
      if(flmin .ge. 1.0)goto 18121                                         2978
      eqs=max(eps,flmin)                                                   2978
      alf=eqs**(1.0/(nlam-1))                                              2978
18121 continue                                                             2979
      rsq=ys0                                                              2979
      a=0.0                                                                2979
      mm=0                                                                 2979
      nlp=0                                                                2979
      nin=nlp                                                              2979
      iz=0                                                                 2979
      mnl=min(mnlam,nlam)                                                  2979
      alm=0.0                                                              2980
18130 do 18131 j=1,ni                                                      2980
      if(ju(j).eq.0)goto 18131                                             2980
      g(j)=0.0                                                             2981
18140 do 18141 k=1,nr                                                      2981
      g(j)=g(j)+dot_product(y(:,k),x(:,j))**2                              2981
18141 continue                                                             2982
18142 continue                                                             2982
      g(j)=sqrt(g(j))                                                      2983
18131 continue                                                             2984
18132 continue                                                             2984
18150 do 18151 m=1,nlam                                                    2984
      alm0=alm                                                             2985
      if(flmin .lt. 1.0)goto 18171                                         2985
      alm=ulam(m)                                                          2985
      goto 18161                                                           2986
18171 if(m .le. 2)goto 18181                                               2986
      alm=alm*alf                                                          2986
      goto 18161                                                           2987
18181 if(m .ne. 1)goto 18191                                               2987
      alm=big                                                              2987
      goto 18201                                                           2988
18191 continue                                                             2988
      alm0=0.0                                                             2989
18210 do 18211 j=1,ni                                                      2989
      if(ju(j).eq.0)goto 18211                                             2990
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           2991
18211 continue                                                             2992
18212 continue                                                             2992
      alm0=alm0/max(bta,1.0e-3)                                            2992
      alm=alf*alm0                                                         2993
18201 continue                                                             2994
18161 continue                                                             2994
      dem=alm*omb                                                          2994
      ab=alm*bta                                                           2994
      rsq0=rsq                                                             2994
      jz=1                                                                 2995
      tlam=bta*(2.0*alm-alm0)                                              2996
18220 do 18221 k=1,ni                                                      2996
      if(ix(k).eq.1)goto 18221                                             2996
      if(ju(k).eq.0)goto 18221                                             2997
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       2998
18221 continue                                                             2999
18222 continue                                                             2999
18230 continue                                                             2999
18231 continue                                                             2999
      if(iz*jz.ne.0) go to 10260                                           3000
10680 continue                                                             3000
      nlp=nlp+1                                                            3000
      dlx=0.0                                                              3001
18240 do 18241 k=1,ni                                                      3001
      if(ix(k).eq.0)goto 18241                                             3001
      gkn=0.0                                                              3002
18250 do 18251 j=1,nr                                                      3002
      gj(j)=dot_product(y(:,j),x(:,k))                                     3003
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3003
      gkn=gkn+gk(j)**2                                                     3005
18251 continue                                                             3005
18252 continue                                                             3005
      gkn=sqrt(gkn)                                                        3005
      u=1.0-ab*vp(k)/gkn                                                   3005
      del=a(:,k)                                                           3006
      if(u .gt. 0.0)goto 18271                                             3006
      a(:,k)=0.0                                                           3006
      goto 18281                                                           3007
18271 continue                                                             3007
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3007
18281 continue                                                             3008
18261 continue                                                             3008
      del=a(:,k)-del                                                       3008
      if(maxval(abs(del)).le.0.0)goto 18241                                3009
18290 do 18291 j=1,nr                                                      3009
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3010
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3010
      dlx=max(dlx,xv(k)*del(j)**2)                                         3011
18291 continue                                                             3012
18292 continue                                                             3012
      if(mm(k) .ne. 0)goto 18311                                           3012
      nin=nin+1                                                            3012
      if(nin.gt.nx)goto 18242                                              3013
      mm(k)=nin                                                            3013
      ia(nin)=k                                                            3014
18311 continue                                                             3015
18241 continue                                                             3016
18242 continue                                                             3016
      if(nin.gt.nx)goto 18232                                              3017
      if(dlx .ge. thr)goto 18331                                           3017
      ixx=0                                                                3018
18340 do 18341 k=1,ni                                                      3018
      if(ix(k).eq.1)goto 18341                                             3018
      if(ju(k).eq.0)goto 18341                                             3018
      g(k)=0.0                                                             3019
18350 do 18351 j=1,nr                                                      3019
      g(k)=g(k)+dot_product(y(:,j),x(:,k))**2                              3019
18351 continue                                                             3020
18352 continue                                                             3020
      g(k)=sqrt(g(k))                                                      3021
      if(g(k) .le. ab*vp(k))goto 18371                                     3021
      ix(k)=1                                                              3021
      ixx=1                                                                3021
18371 continue                                                             3022
18341 continue                                                             3023
18342 continue                                                             3023
      if(ixx.eq.1) go to 10680                                             3024
      goto 18232                                                           3025
18331 continue                                                             3026
      if(nlp .le. maxit)goto 18391                                         3026
      jerr=-m                                                              3026
      return                                                               3026
18391 continue                                                             3027
10260 continue                                                             3027
      iz=1                                                                 3028
18400 continue                                                             3028
18401 continue                                                             3028
      nlp=nlp+1                                                            3028
      dlx=0.0                                                              3029
18410 do 18411 l=1,nin                                                     3029
      k=ia(l)                                                              3029
      gkn=0.0                                                              3030
18420 do 18421 j=1,nr                                                      3030
      gj(j)=dot_product(y(:,j),x(:,k))                                     3031
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3031
      gkn=gkn+gk(j)**2                                                     3033
18421 continue                                                             3033
18422 continue                                                             3033
      gkn=sqrt(gkn)                                                        3033
      u=1.0-ab*vp(k)/gkn                                                   3033
      del=a(:,k)                                                           3034
      if(u .gt. 0.0)goto 18441                                             3034
      a(:,k)=0.0                                                           3034
      goto 18451                                                           3035
18441 continue                                                             3035
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3035
18451 continue                                                             3036
18431 continue                                                             3036
      del=a(:,k)-del                                                       3036
      if(maxval(abs(del)).le.0.0)goto 18411                                3037
18460 do 18461 j=1,nr                                                      3037
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3038
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3038
      dlx=max(dlx,xv(k)*del(j)**2)                                         3039
18461 continue                                                             3040
18462 continue                                                             3040
18411 continue                                                             3041
18412 continue                                                             3041
      if(dlx.lt.thr)goto 18402                                             3041
      if(nlp .le. maxit)goto 18481                                         3041
      jerr=-m                                                              3041
      return                                                               3041
18481 continue                                                             3042
      goto 18401                                                           3043
18402 continue                                                             3043
      jz=0                                                                 3044
      goto 18231                                                           3045
18232 continue                                                             3045
      if(nin .le. nx)goto 18501                                            3045
      jerr=-10000-m                                                        3045
      goto 18152                                                           3045
18501 continue                                                             3046
      if(nin .le. 0)goto 18521                                             3046
18530 do 18531 j=1,nr                                                      3046
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3046
18531 continue                                                             3046
18532 continue                                                             3046
18521 continue                                                             3047
      kin(m)=nin                                                           3048
      rsqo(m)=1.0-rsq/ys0                                                  3048
      almo(m)=alm                                                          3048
      lmu=m                                                                3049
      if(m.lt.mnl)goto 18151                                               3049
      if(flmin.ge.1.0)goto 18151                                           3050
      me=0                                                                 3050
18540 do 18541 j=1,nin                                                     3050
      if(ao(j,1,m).ne.0.0) me=me+1                                         3050
18541 continue                                                             3050
18542 continue                                                             3050
      if(me.gt.ne)goto 18152                                               3051
      if(rsq0-rsq.lt.sml*rsq)goto 18152                                    3051
      if(rsqo(m).gt.rsqmax)goto 18152                                      3052
18151 continue                                                             3053
18152 continue                                                             3053
      deallocate(a,mm,g,ix,del,gj,gk)                                      3054
      return                                                               3055
      end                                                                  3056
      subroutine multsolns(ni,nx,nr,lmu,a,ia,nin,b)                        3057
      real a(nx,nr,lmu),b(ni,nr,lmu)                                       3057
      integer ia(nx),nin(lmu)                                              3058
18550 do 18551 lam=1,lmu                                                   3058
      call multuncomp(ni,nr,nx,a(1,1,lam),ia,nin(lam),b(1,1,lam))          3058
18551 continue                                                             3059
18552 continue                                                             3059
      return                                                               3060
      end                                                                  3061
      subroutine multuncomp(ni,nr,nx,ca,ia,nin,a)                          3062
      real ca(nx,nr),a(ni,nr)                                              3062
      integer ia(nx)                                                       3063
      a=0.0                                                                3064
      if(nin .le. 0)goto 18571                                             3064
18580 do 18581 j=1,nr                                                      3064
      a(ia(1:nin),j)=ca(1:nin,j)                                           3064
18581 continue                                                             3064
18582 continue                                                             3064
18571 continue                                                             3065
      return                                                               3066
      end                                                                  3067
      subroutine multmodval(nx,nr,a0,ca,ia,nin,n,x,f)                      3068
      real a0(nr),ca(nx,nr),x(n,*),f(nr,n)                                 3068
      integer ia(nx)                                                       3069
18590 do 18591 i=1,n                                                       3069
      f(:,i)=a0                                                            3069
18591 continue                                                             3069
18592 continue                                                             3069
      if(nin.le.0) return                                                  3070
18600 do 18601 i=1,n                                                       3070
18610 do 18611 j=1,nr                                                      3070
      f(j,i)=f(j,i)+dot_product(ca(1:nin,j),x(i,ia(1:nin)))                3070
18611 continue                                                             3070
18612 continue                                                             3070
18601 continue                                                             3071
18602 continue                                                             3071
      return                                                               3072
      end                                                                  3073
      subroutine multspelnet  (parm,no,ni,nr,x,ix,jx,y,w,jd,vp,ne,nx,nla   3076 
     *m,flmin,ulam,thr,isd,jsd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      real x(*),y(no,nr),w(no),vp(ni),ulam(nlam)                           3077
      real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                  3078
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3079
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 18631                                    3082
      jerr=10000                                                           3082
      return                                                               3082
18631 continue                                                             3083
      allocate(vq(1:ni),stat=jerr)                                         3083
      if(jerr.ne.0) return                                                 3084
      vq=max(0.0,vp)                                                       3084
      vq=vq*ni/sum(vq)                                                     3085
      call multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin   3087 
     *,  ulam,thr,isd,jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      deallocate(vq)                                                       3088
      return                                                               3089
      end                                                                  3090
      subroutine multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,ne,nx,nlam   3092 
     *,flmin,ulam,  thr,isd,jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),vp(ni),y(no,nr),w(no),ulam(nlam)                           3093
      real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                  3094
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3095
      real, dimension (:), allocatable :: xm,xs,xv,ym,ys                        
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         3100
      allocate(xs(1:ni),stat=ierr)                                         3100
      jerr=jerr+ierr                                                       3101
      allocate(ym(1:nr),stat=ierr)                                         3101
      jerr=jerr+ierr                                                       3102
      allocate(ys(1:nr),stat=ierr)                                         3102
      jerr=jerr+ierr                                                       3103
      allocate(ju(1:ni),stat=ierr)                                         3103
      jerr=jerr+ierr                                                       3104
      allocate(xv(1:ni),stat=ierr)                                         3104
      jerr=jerr+ierr                                                       3105
      if(jerr.ne.0) return                                                 3106
      call spchkvars(no,ni,x,ix,ju)                                        3107
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3108
      if(maxval(ju) .gt. 0)goto 18651                                      3108
      jerr=7777                                                            3108
      return                                                               3108
18651 continue                                                             3109
      call multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,xm,xs,ym,ys,x   3110 
     *v,ys0,jerr)
      if(jerr.ne.0) return                                                 3111
      call multspelnet2(parm,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin   3113 
     *,  ulam,thr,maxit,xm,xs,xv,ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 3114
18660 do 18661 k=1,lmu                                                     3114
      nk=nin(k)                                                            3115
18670 do 18671 j=1,nr                                                      3116
18680 do 18681 l=1,nk                                                      3116
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  3116
18681 continue                                                             3117
18682 continue                                                             3117
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 3118
18671 continue                                                             3119
18672 continue                                                             3119
18661 continue                                                             3120
18662 continue                                                             3120
      deallocate(xm,xs,ym,ys,ju,xv)                                        3121
      return                                                               3122
      end                                                                  3123
      subroutine multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,  xm,xs   3125 
     *,ym,ys,xv,ys0,jerr)
      real x(*),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),ys(nr)          3126
      integer ix(*),jx(*),ju(ni)                                           3127
      w=w/sum(w)                                                           3128
18690 do 18691 j=1,ni                                                      3128
      if(ju(j).eq.0)goto 18691                                             3129
      jb=ix(j)                                                             3129
      je=ix(j+1)-1                                                         3129
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3130
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 3131
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       3132
18691 continue                                                             3133
18692 continue                                                             3133
      if(isd .ne. 0)goto 18711                                             3133
      xs=1.0                                                               3133
      goto 18721                                                           3133
18711 continue                                                             3133
      xv=1.0                                                               3133
18721 continue                                                             3134
18701 continue                                                             3134
      ys0=0.0                                                              3135
18730 do 18731 j=1,nr                                                      3136
      ym(j)=dot_product(w,y(:,j))                                          3136
      y(:,j)=y(:,j)-ym(j)                                                  3137
      z=dot_product(w,y(:,j)**2)                                           3138
      if(jsd .le. 0)goto 18751                                             3138
      ys(j)=sqrt(z)                                                        3138
      y(:,j)=y(:,j)/ys(j)                                                  3138
      goto 18761                                                           3139
18751 continue                                                             3139
      ys0=ys0+z                                                            3139
18761 continue                                                             3140
18741 continue                                                             3140
18731 continue                                                             3141
18732 continue                                                             3141
      if(jsd .ne. 0)goto 18781                                             3141
      ys=1.0                                                               3141
      goto 18791                                                           3141
18781 continue                                                             3141
      ys0=nr                                                               3141
18791 continue                                                             3142
18771 continue                                                             3142
      return                                                               3143
      end                                                                  3144
      subroutine multspelnet2(beta,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam   3146 
     *,flmin,ulam,  thri,maxit,xm,xs,xv,ys0,lmu,ao,ia,kin,rsqo,almo,nlp,
     *jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99   3147 
     *9)
      real y(no,nr),w(no),x(*),vp(ni),ulam(nlam)                           3148
      real ao(nx,nr,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)       3149
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          3150
      real, dimension (:), allocatable :: g,gj,gk,del,o                         
      integer, dimension (:), allocatable :: mm,iy                              
      real, dimension (:,:), allocatable :: a                                   
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      allocate(mm(1:ni),stat=ierr)                                         3156
      jerr=jerr+ierr                                                       3157
      allocate(g(1:ni),stat=ierr)                                          3157
      jerr=jerr+ierr                                                       3158
      allocate(gj(1:nr),stat=ierr)                                         3158
      jerr=jerr+ierr                                                       3159
      allocate(gk(1:nr),stat=ierr)                                         3159
      jerr=jerr+ierr                                                       3160
      allocate(del(1:nr),stat=ierr)                                        3160
      jerr=jerr+ierr                                                       3161
      allocate(o(1:nr),stat=ierr)                                          3161
      jerr=jerr+ierr                                                       3162
      allocate(iy(1:ni),stat=ierr)                                         3162
      jerr=jerr+ierr                                                       3163
      if(jerr.ne.0) return                                                 3164
      bta=beta                                                             3164
      omb=1.0-bta                                                          3164
      alm=0.0                                                              3164
      iy=0                                                                 3164
      thr=thri*ys0/nr                                                      3165
      if(flmin .ge. 1.0)goto 18811                                         3165
      eqs=max(eps,flmin)                                                   3165
      alf=eqs**(1.0/(nlam-1))                                              3165
18811 continue                                                             3166
      rsq=ys0                                                              3166
      a=0.0                                                                3166
      mm=0                                                                 3166
      o=0.0                                                                3166
      nlp=0                                                                3166
      nin=nlp                                                              3166
      iz=0                                                                 3166
      mnl=min(mnlam,nlam)                                                  3167
18820 do 18821 j=1,ni                                                      3167
      if(ju(j).eq.0)goto 18821                                             3167
      jb=ix(j)                                                             3167
      je=ix(j+1)-1                                                         3167
      g(j)=0.0                                                             3168
18830 do 18831 k=1,nr                                                      3169
      g(j)=g(j)+(dot_product(y(jx(jb:je),k),w(jx(jb:je))*x(jb:je))/xs(j)   3170 
     *)**2
18831 continue                                                             3171
18832 continue                                                             3171
      g(j)=sqrt(g(j))                                                      3172
18821 continue                                                             3173
18822 continue                                                             3173
18840 do 18841 m=1,nlam                                                    3173
      alm0=alm                                                             3174
      if(flmin .lt. 1.0)goto 18861                                         3174
      alm=ulam(m)                                                          3174
      goto 18851                                                           3175
18861 if(m .le. 2)goto 18871                                               3175
      alm=alm*alf                                                          3175
      goto 18851                                                           3176
18871 if(m .ne. 1)goto 18881                                               3176
      alm=big                                                              3176
      goto 18891                                                           3177
18881 continue                                                             3177
      alm0=0.0                                                             3178
18900 do 18901 j=1,ni                                                      3178
      if(ju(j).eq.0)goto 18901                                             3179
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           3180
18901 continue                                                             3181
18902 continue                                                             3181
      alm0=alm0/max(bta,1.0e-3)                                            3181
      alm=alf*alm0                                                         3182
18891 continue                                                             3183
18851 continue                                                             3183
      dem=alm*omb                                                          3183
      ab=alm*bta                                                           3183
      rsq0=rsq                                                             3183
      jz=1                                                                 3184
      tlam=bta*(2.0*alm-alm0)                                              3185
18910 do 18911 k=1,ni                                                      3185
      if(iy(k).eq.1)goto 18911                                             3185
      if(ju(k).eq.0)goto 18911                                             3186
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       3187
18911 continue                                                             3188
18912 continue                                                             3188
18920 continue                                                             3188
18921 continue                                                             3188
      if(iz*jz.ne.0) go to 10260                                           3189
10680 continue                                                             3189
      nlp=nlp+1                                                            3189
      dlx=0.0                                                              3190
18930 do 18931 k=1,ni                                                      3190
      if(iy(k).eq.0)goto 18931                                             3190
      jb=ix(k)                                                             3190
      je=ix(k+1)-1                                                         3190
      gkn=0.0                                                              3191
18940 do 18941 j=1,nr                                                      3192
      gj(j)=dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(k)   3193
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3193
      gkn=gkn+gk(j)**2                                                     3194
18941 continue                                                             3195
18942 continue                                                             3195
      gkn=sqrt(gkn)                                                        3195
      u=1.0-ab*vp(k)/gkn                                                   3195
      del=a(:,k)                                                           3196
      if(u .gt. 0.0)goto 18961                                             3196
      a(:,k)=0.0                                                           3196
      goto 18971                                                           3197
18961 continue                                                             3197
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3197
18971 continue                                                             3198
18951 continue                                                             3198
      del=a(:,k)-del                                                       3198
      if(maxval(abs(del)).le.0.0)goto 18931                                3199
      if(mm(k) .ne. 0)goto 18991                                           3199
      nin=nin+1                                                            3199
      if(nin.gt.nx)goto 18932                                              3200
      mm(k)=nin                                                            3200
      ia(nin)=k                                                            3201
18991 continue                                                             3202
19000 do 19001 j=1,nr                                                      3202
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3203
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3204
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3204
      dlx=max(xv(k)*del(j)**2,dlx)                                         3205
19001 continue                                                             3206
19002 continue                                                             3206
18931 continue                                                             3207
18932 continue                                                             3207
      if(nin.gt.nx)goto 18922                                              3208
      if(dlx .ge. thr)goto 19021                                           3208
      ixx=0                                                                3209
19030 do 19031 j=1,ni                                                      3209
      if(iy(j).eq.1)goto 19031                                             3209
      if(ju(j).eq.0)goto 19031                                             3210
      jb=ix(j)                                                             3210
      je=ix(j+1)-1                                                         3210
      g(j)=0.0                                                             3211
19040 do 19041 k=1,nr                                                      3211
      g(j)=g(j)+  (dot_product(y(jx(jb:je),k)+o(k),w(jx(jb:je))*x(jb:je)   3213 
     *)/xs(j))**2
19041 continue                                                             3214
19042 continue                                                             3214
      g(j)=sqrt(g(j))                                                      3215
      if(g(j) .le. ab*vp(j))goto 19061                                     3215
      iy(j)=1                                                              3215
      ixx=1                                                                3215
19061 continue                                                             3216
19031 continue                                                             3217
19032 continue                                                             3217
      if(ixx.eq.1) go to 10680                                             3218
      goto 18922                                                           3219
19021 continue                                                             3220
      if(nlp .le. maxit)goto 19081                                         3220
      jerr=-m                                                              3220
      return                                                               3220
19081 continue                                                             3221
10260 continue                                                             3221
      iz=1                                                                 3222
19090 continue                                                             3222
19091 continue                                                             3222
      nlp=nlp+1                                                            3222
      dlx=0.0                                                              3223
19100 do 19101 l=1,nin                                                     3223
      k=ia(l)                                                              3223
      jb=ix(k)                                                             3223
      je=ix(k+1)-1                                                         3223
      gkn=0.0                                                              3224
19110 do 19111 j=1,nr                                                      3224
      gj(j)=  dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(   3226 
     *k)
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3226
      gkn=gkn+gk(j)**2                                                     3227
19111 continue                                                             3228
19112 continue                                                             3228
      gkn=sqrt(gkn)                                                        3228
      u=1.0-ab*vp(k)/gkn                                                   3228
      del=a(:,k)                                                           3229
      if(u .gt. 0.0)goto 19131                                             3229
      a(:,k)=0.0                                                           3229
      goto 19141                                                           3230
19131 continue                                                             3230
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3230
19141 continue                                                             3231
19121 continue                                                             3231
      del=a(:,k)-del                                                       3231
      if(maxval(abs(del)).le.0.0)goto 19101                                3232
19150 do 19151 j=1,nr                                                      3232
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3233
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3234
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3234
      dlx=max(xv(k)*del(j)**2,dlx)                                         3235
19151 continue                                                             3236
19152 continue                                                             3236
19101 continue                                                             3237
19102 continue                                                             3237
      if(dlx.lt.thr)goto 19092                                             3237
      if(nlp .le. maxit)goto 19171                                         3237
      jerr=-m                                                              3237
      return                                                               3237
19171 continue                                                             3238
      goto 19091                                                           3239
19092 continue                                                             3239
      jz=0                                                                 3240
      goto 18921                                                           3241
18922 continue                                                             3241
      if(nin .le. nx)goto 19191                                            3241
      jerr=-10000-m                                                        3241
      goto 18842                                                           3241
19191 continue                                                             3242
      if(nin .le. 0)goto 19211                                             3242
19220 do 19221 j=1,nr                                                      3242
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3242
19221 continue                                                             3242
19222 continue                                                             3242
19211 continue                                                             3243
      kin(m)=nin                                                           3244
      rsqo(m)=1.0-rsq/ys0                                                  3244
      almo(m)=alm                                                          3244
      lmu=m                                                                3245
      if(m.lt.mnl)goto 18841                                               3245
      if(flmin.ge.1.0)goto 18841                                           3246
      me=0                                                                 3246
19230 do 19231 j=1,nin                                                     3246
      if(ao(j,1,m).ne.0.0) me=me+1                                         3246
19231 continue                                                             3246
19232 continue                                                             3246
      if(me.gt.ne)goto 18842                                               3247
      if(rsq0-rsq.lt.sml*rsq)goto 18842                                    3247
      if(rsqo(m).gt.rsqmax)goto 18842                                      3248
18841 continue                                                             3249
18842 continue                                                             3249
      deallocate(a,mm,g,iy,gj,gk,del,o)                                    3250
      return                                                               3251
      end                                                                  3252
      subroutine multlognetn(parm,no,ni,nc,x,y,g,w,ju,vp,ne,nx,nlam,flmi   3254 
     *n,ulam,shri,  isd,maxit,xv,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   3256 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam)              3257
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),xv(ni)            3258
      integer ju(ni),m(nx),kin(nlam)                                       3259
      real, dimension (:,:), allocatable :: q,r,b,bs                            
      real, dimension (:), allocatable :: sxp,sxpl,ga,gk,del                    
      integer, dimension (:), allocatable :: mm,is,ixx                          
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no,1:nc),stat=ierr); jerr=jerr+ierr;                         
      allocate(mm(1:ni),stat=ierr)                                         3268
      jerr=jerr+ierr                                                       3269
      allocate(is(1:max(nc,ni)),stat=ierr)                                 3269
      jerr=jerr+ierr                                                       3270
      allocate(sxp(1:no),stat=ierr)                                        3270
      jerr=jerr+ierr                                                       3271
      allocate(sxpl(1:no),stat=ierr)                                       3271
      jerr=jerr+ierr                                                       3272
      allocate(ga(1:ni),stat=ierr)                                         3272
      jerr=jerr+ierr                                                       3273
      allocate(ixx(1:ni),stat=ierr)                                        3273
      jerr=jerr+ierr                                                       3274
      allocate(gk(1:nc),stat=ierr)                                         3274
      jerr=jerr+ierr                                                       3275
      allocate(del(1:nc),stat=ierr)                                        3275
      jerr=jerr+ierr                                                       3276
      if(jerr.ne.0) return                                                 3277
      pmax=1.0-pmin                                                        3277
      emin=pmin/pmax                                                       3277
      emax=1.0/emin                                                        3278
      bta=parm                                                             3278
      omb=1.0-bta                                                          3278
      dev1=0.0                                                             3278
      dev0=0.0                                                             3279
19240 do 19241 ic=1,nc                                                     3279
      q0=dot_product(w,y(:,ic))                                            3280
      if(q0 .gt. pmin)goto 19261                                           3280
      jerr =8000+ic                                                        3280
      return                                                               3280
19261 continue                                                             3281
      if(q0 .lt. pmax)goto 19281                                           3281
      jerr =9000+ic                                                        3281
      return                                                               3281
19281 continue                                                             3282
      b(0,ic)=log(q0)                                                      3282
      dev1=dev1-q0*b(0,ic)                                                 3282
      b(1:ni,ic)=0.0                                                       3283
19241 continue                                                             3284
19242 continue                                                             3284
      ixx=0                                                                3284
      al=0.0                                                               3285
      if(nonzero(no*nc,g) .ne. 0)goto 19301                                3286
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         3286
      sxp=0.0                                                              3287
19310 do 19311 ic=1,nc                                                     3287
      q(:,ic)=exp(b(0,ic))                                                 3287
      sxp=sxp+q(:,ic)                                                      3287
19311 continue                                                             3288
19312 continue                                                             3288
      goto 19321                                                           3289
19301 continue                                                             3289
19330 do 19331 i=1,no                                                      3289
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3289
19331 continue                                                             3289
19332 continue                                                             3289
      sxp=0.0                                                              3290
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 3290
      if(jerr.ne.0) return                                                 3291
      dev1=0.0                                                             3292
19340 do 19341 ic=1,nc                                                     3292
      q(:,ic)=b(0,ic)+g(:,ic)                                              3293
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             3294
      q(:,ic)=exp(q(:,ic))                                                 3294
      sxp=sxp+q(:,ic)                                                      3295
19341 continue                                                             3296
19342 continue                                                             3296
      sxpl=w*log(sxp)                                                      3296
19350 do 19351 ic=1,nc                                                     3296
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  3296
19351 continue                                                             3297
19352 continue                                                             3297
19321 continue                                                             3298
19291 continue                                                             3298
19360 do 19361 ic=1,nc                                                     3298
19370 do 19371 i=1,no                                                      3298
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               3298
19371 continue                                                             3298
19372 continue                                                             3298
19361 continue                                                             3299
19362 continue                                                             3299
      dev0=dev0+dev1                                                       3300
      if(flmin .ge. 1.0)goto 19391                                         3300
      eqs=max(eps,flmin)                                                   3300
      alf=eqs**(1.0/(nlam-1))                                              3300
19391 continue                                                             3301
      m=0                                                                  3301
      mm=0                                                                 3301
      nin=0                                                                3301
      nlp=0                                                                3301
      mnl=min(mnlam,nlam)                                                  3301
      bs=0.0                                                               3301
      shr=shri*dev0                                                        3302
      ga=0.0                                                               3303
19400 do 19401 ic=1,nc                                                     3303
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3304
19410 do 19411 j=1,ni                                                      3304
      if(ju(j).ne.0) ga(j)=ga(j)+dot_product(r(:,ic),x(:,j))**2            3304
19411 continue                                                             3305
19412 continue                                                             3305
19401 continue                                                             3306
19402 continue                                                             3306
      ga=sqrt(ga)                                                          3307
19420 do 19421 ilm=1,nlam                                                  3307
      al0=al                                                               3308
      if(flmin .lt. 1.0)goto 19441                                         3308
      al=ulam(ilm)                                                         3308
      goto 19431                                                           3309
19441 if(ilm .le. 2)goto 19451                                             3309
      al=al*alf                                                            3309
      goto 19431                                                           3310
19451 if(ilm .ne. 1)goto 19461                                             3310
      al=big                                                               3310
      goto 19471                                                           3311
19461 continue                                                             3311
      al0=0.0                                                              3312
19480 do 19481 j=1,ni                                                      3312
      if(ju(j).eq.0)goto 19481                                             3312
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3312
19481 continue                                                             3313
19482 continue                                                             3313
      al0=al0/max(bta,1.0e-3)                                              3313
      al=alf*al0                                                           3314
19471 continue                                                             3315
19431 continue                                                             3315
      al2=al*omb                                                           3315
      al1=al*bta                                                           3315
      tlam=bta*(2.0*al-al0)                                                3316
19490 do 19491 k=1,ni                                                      3316
      if(ixx(k).eq.1)goto 19491                                            3316
      if(ju(k).eq.0)goto 19491                                             3317
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     3318
19491 continue                                                             3319
19492 continue                                                             3319
10680 continue                                                             3320
19500 continue                                                             3320
19501 continue                                                             3320
      ix=0                                                                 3320
      jx=ix                                                                3320
      kx=jx                                                                3320
      t=0.0                                                                3321
19510 do 19511 ic=1,nc                                                     3321
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       3321
19511 continue                                                             3322
19512 continue                                                             3322
      if(t .ge. eps)goto 19531                                             3322
      kx=1                                                                 3322
      goto 19502                                                           3322
19531 continue                                                             3322
      t=2.0*t                                                              3322
      alt=al1/t                                                            3322
      al2t=al2/t                                                           3323
19540 do 19541 ic=1,nc                                                     3324
      bs(0,ic)=b(0,ic)                                                     3324
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          3325
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    3325
      d=sum(r(:,ic))                                                       3326
      b(0,ic)=b(0,ic)+d                                                    3326
      r(:,ic)=r(:,ic)-d*w                                                  3326
      dlx=max(dlx,d**2)                                                    3328
19541 continue                                                             3329
19542 continue                                                             3329
19550 continue                                                             3329
19551 continue                                                             3329
      nlp=nlp+nc                                                           3329
      dlx=0.0                                                              3330
19560 do 19561 k=1,ni                                                      3330
      if(ixx(k).eq.0)goto 19561                                            3330
      gkn=0.0                                                              3331
19570 do 19571 ic=1,nc                                                     3331
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     3332
      gkn=gkn+gk(ic)**2                                                    3333
19571 continue                                                             3334
19572 continue                                                             3334
      gkn=sqrt(gkn)                                                        3334
      u=1.0-alt*vp(k)/gkn                                                  3334
      del=b(k,:)                                                           3335
      if(u .gt. 0.0)goto 19591                                             3335
      b(k,:)=0.0                                                           3335
      goto 19601                                                           3336
19591 continue                                                             3336
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     3336
19601 continue                                                             3337
19581 continue                                                             3337
      del=b(k,:)-del                                                       3337
      if(maxval(abs(del)).le.0.0)goto 19561                                3338
19610 do 19611 ic=1,nc                                                     3338
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3339
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     3340
19611 continue                                                             3341
19612 continue                                                             3341
      if(mm(k) .ne. 0)goto 19631                                           3341
      nin=nin+1                                                            3342
      if(nin .le. nx)goto 19651                                            3342
      jx=1                                                                 3342
      goto 19562                                                           3342
19651 continue                                                             3343
      mm(k)=nin                                                            3343
      m(nin)=k                                                             3344
19631 continue                                                             3345
19561 continue                                                             3346
19562 continue                                                             3346
      if(jx.gt.0)goto 19552                                                3346
      if(dlx.lt.shr)goto 19552                                             3347
      if(nlp .le. maxit)goto 19671                                         3347
      jerr=-ilm                                                            3347
      return                                                               3347
19671 continue                                                             3348
19680 continue                                                             3348
19681 continue                                                             3348
      nlp=nlp+nc                                                           3348
      dlx=0.0                                                              3349
19690 do 19691 l=1,nin                                                     3349
      k=m(l)                                                               3349
      gkn=0.0                                                              3350
19700 do 19701 ic=1,nc                                                     3350
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     3351
      gkn=gkn+gk(ic)**2                                                    3352
19701 continue                                                             3353
19702 continue                                                             3353
      gkn=sqrt(gkn)                                                        3353
      u=1.0-alt*vp(k)/gkn                                                  3353
      del=b(k,:)                                                           3354
      if(u .gt. 0.0)goto 19721                                             3354
      b(k,:)=0.0                                                           3354
      goto 19731                                                           3355
19721 continue                                                             3355
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     3355
19731 continue                                                             3356
19711 continue                                                             3356
      del=b(k,:)-del                                                       3356
      if(maxval(abs(del)).le.0.0)goto 19691                                3357
19740 do 19741 ic=1,nc                                                     3357
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3358
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     3359
19741 continue                                                             3360
19742 continue                                                             3360
19691 continue                                                             3361
19692 continue                                                             3361
      if(dlx.lt.shr)goto 19682                                             3361
      if(nlp .le. maxit)goto 19761                                         3361
      jerr=-ilm                                                            3361
      return                                                               3361
19761 continue                                                             3363
      goto 19681                                                           3364
19682 continue                                                             3364
      goto 19551                                                           3365
19552 continue                                                             3365
      if(jx.gt.0)goto 19502                                                3366
19770 do 19771 ic=1,nc                                                     3367
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                                3368
      if(ix .ne. 0)goto 19791                                              3369
19800 do 19801 j=1,nin                                                     3369
      k=m(j)                                                               3370
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 19821                   3370
      ix=1                                                                 3370
      goto 19802                                                           3370
19821 continue                                                             3372
19801 continue                                                             3373
19802 continue                                                             3373
19791 continue                                                             3374
19830 do 19831 i=1,no                                                      3374
      fi=b(0,ic)+g(i,ic)                                                   3376
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         3377
      fi=min(max(exmn,fi),exmx)                                            3377
      sxp(i)=sxp(i)-q(i,ic)                                                3378
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    3379
      sxp(i)=sxp(i)+q(i,ic)                                                3380
19831 continue                                                             3381
19832 continue                                                             3381
19771 continue                                                             3382
19772 continue                                                             3382
      s=-sum(b(0,:))/nc                                                    3382
      b(0,:)=b(0,:)+s                                                      3383
      if(jx.gt.0)goto 19502                                                3384
      if(ix .ne. 0)goto 19851                                              3385
19860 do 19861 k=1,ni                                                      3385
      if(ixx(k).eq.1)goto 19861                                            3385
      if(ju(k).eq.0)goto 19861                                             3385
      ga(k)=0.0                                                            3385
19861 continue                                                             3386
19862 continue                                                             3386
19870 do 19871 ic=1,nc                                                     3386
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3387
19880 do 19881 k=1,ni                                                      3387
      if(ixx(k).eq.1)goto 19881                                            3387
      if(ju(k).eq.0)goto 19881                                             3388
      ga(k)=ga(k)+dot_product(r(:,ic),x(:,k))**2                           3389
19881 continue                                                             3390
19882 continue                                                             3390
19871 continue                                                             3391
19872 continue                                                             3391
      ga=sqrt(ga)                                                          3392
19890 do 19891 k=1,ni                                                      3392
      if(ixx(k).eq.1)goto 19891                                            3392
      if(ju(k).eq.0)goto 19891                                             3393
      if(ga(k) .le. al1*vp(k))goto 19911                                   3393
      ixx(k)=1                                                             3393
      ix=1                                                                 3393
19911 continue                                                             3394
19891 continue                                                             3395
19892 continue                                                             3395
      if(ix.eq.1) go to 10680                                              3396
      goto 19502                                                           3397
19851 continue                                                             3398
      goto 19501                                                           3399
19502 continue                                                             3399
      if(kx .le. 0)goto 19931                                              3399
      jerr=-20000-ilm                                                      3399
      goto 19422                                                           3399
19931 continue                                                             3400
      if(jx .le. 0)goto 19951                                              3400
      jerr=-10000-ilm                                                      3400
      goto 19422                                                           3400
19951 continue                                                             3400
      devi=0.0                                                             3401
19960 do 19961 ic=1,nc                                                     3402
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          3402
      a0(ic,ilm)=b(0,ic)                                                   3403
19970 do 19971 i=1,no                                                      3403
      if(y(i,ic).le.0.0)goto 19971                                         3404
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           3405
19971 continue                                                             3406
19972 continue                                                             3406
19961 continue                                                             3407
19962 continue                                                             3407
      kin(ilm)=nin                                                         3407
      alm(ilm)=al                                                          3407
      lmu=ilm                                                              3408
      dev(ilm)=(dev1-devi)/dev0                                            3409
      if(ilm.lt.mnl)goto 19421                                             3409
      if(flmin.ge.1.0)goto 19421                                           3410
      me=0                                                                 3410
19980 do 19981 j=1,nin                                                     3410
      if(a(j,1,ilm).ne.0.0) me=me+1                                        3410
19981 continue                                                             3410
19982 continue                                                             3410
      if(me.gt.ne)goto 19422                                               3411
      if(dev(ilm).gt.devmax)goto 19422                                     3411
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 19422                             3412
19421 continue                                                             3413
19422 continue                                                             3413
      g=log(q)                                                             3413
19990 do 19991 i=1,no                                                      3413
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3413
19991 continue                                                             3414
19992 continue                                                             3414
      deallocate(sxp,b,bs,r,q,mm,is,ga,ixx,gk,del,sxpl)                    3415
      return                                                               3416
      end                                                                  3417
      subroutine multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,ne,nx,   3419 
     *nlam,flmin,  ulam,shri,isd,maxit,xv,xb,xs,lmu,a0,a,m,kin,dev0,dev,
     *alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   3421 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni),   3422 
     *xv(ni)
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   3423
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           3424
      real, dimension (:,:), allocatable :: q,r,b,bs                            
      real, dimension (:), allocatable :: sxp,sxpl,ga,gk,del,sc,svr             
      integer, dimension (:), allocatable :: mm,is,iy                           
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(mm(1:ni),stat=ierr)                                         3433
      jerr=jerr+ierr                                                       3434
      allocate(ga(1:ni),stat=ierr)                                         3434
      jerr=jerr+ierr                                                       3435
      allocate(gk(1:nc),stat=ierr)                                         3435
      jerr=jerr+ierr                                                       3436
      allocate(del(1:nc),stat=ierr)                                        3436
      jerr=jerr+ierr                                                       3437
      allocate(iy(1:ni),stat=ierr)                                         3437
      jerr=jerr+ierr                                                       3438
      allocate(is(1:max(nc,ni)),stat=ierr)                                 3438
      jerr=jerr+ierr                                                       3439
      allocate(sxp(1:no),stat=ierr)                                        3439
      jerr=jerr+ierr                                                       3440
      allocate(sxpl(1:no),stat=ierr)                                       3440
      jerr=jerr+ierr                                                       3441
      allocate(svr(1:nc),stat=ierr)                                        3441
      jerr=jerr+ierr                                                       3442
      allocate(sc(1:no),stat=ierr)                                         3442
      jerr=jerr+ierr                                                       3443
      if(jerr.ne.0) return                                                 3444
      pmax=1.0-pmin                                                        3444
      emin=pmin/pmax                                                       3444
      emax=1.0/emin                                                        3445
      bta=parm                                                             3445
      omb=1.0-bta                                                          3445
      dev1=0.0                                                             3445
      dev0=0.0                                                             3446
20000 do 20001 ic=1,nc                                                     3446
      q0=dot_product(w,y(:,ic))                                            3447
      if(q0 .gt. pmin)goto 20021                                           3447
      jerr =8000+ic                                                        3447
      return                                                               3447
20021 continue                                                             3448
      if(q0 .lt. pmax)goto 20041                                           3448
      jerr =9000+ic                                                        3448
      return                                                               3448
20041 continue                                                             3449
      b(1:ni,ic)=0.0                                                       3449
      b(0,ic)=log(q0)                                                      3449
      dev1=dev1-q0*b(0,ic)                                                 3450
20001 continue                                                             3451
20002 continue                                                             3451
      iy=0                                                                 3451
      al=0.0                                                               3452
      if(nonzero(no*nc,g) .ne. 0)goto 20061                                3453
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         3453
      sxp=0.0                                                              3454
20070 do 20071 ic=1,nc                                                     3454
      q(:,ic)=exp(b(0,ic))                                                 3454
      sxp=sxp+q(:,ic)                                                      3454
20071 continue                                                             3455
20072 continue                                                             3455
      goto 20081                                                           3456
20061 continue                                                             3456
20090 do 20091 i=1,no                                                      3456
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3456
20091 continue                                                             3456
20092 continue                                                             3456
      sxp=0.0                                                              3457
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 3457
      if(jerr.ne.0) return                                                 3458
      dev1=0.0                                                             3459
20100 do 20101 ic=1,nc                                                     3459
      q(:,ic)=b(0,ic)+g(:,ic)                                              3460
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             3461
      q(:,ic)=exp(q(:,ic))                                                 3461
      sxp=sxp+q(:,ic)                                                      3462
20101 continue                                                             3463
20102 continue                                                             3463
      sxpl=w*log(sxp)                                                      3463
20110 do 20111 ic=1,nc                                                     3463
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  3463
20111 continue                                                             3464
20112 continue                                                             3464
20081 continue                                                             3465
20051 continue                                                             3465
20120 do 20121 ic=1,nc                                                     3465
20130 do 20131 i=1,no                                                      3465
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               3465
20131 continue                                                             3465
20132 continue                                                             3465
20121 continue                                                             3466
20122 continue                                                             3466
      dev0=dev0+dev1                                                       3467
      if(flmin .ge. 1.0)goto 20151                                         3467
      eqs=max(eps,flmin)                                                   3467
      alf=eqs**(1.0/(nlam-1))                                              3467
20151 continue                                                             3468
      m=0                                                                  3468
      mm=0                                                                 3468
      nin=0                                                                3468
      nlp=0                                                                3468
      mnl=min(mnlam,nlam)                                                  3468
      bs=0.0                                                               3469
      shr=shri*dev0                                                        3469
      ga=0.0                                                               3470
20160 do 20161 ic=1,nc                                                     3470
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3470
      svr(ic)=sum(r(:,ic))                                                 3471
20170 do 20171 j=1,ni                                                      3471
      if(ju(j).eq.0)goto 20171                                             3472
      jb=ix(j)                                                             3472
      je=ix(j+1)-1                                                         3473
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             3474
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            3475
20171 continue                                                             3476
20172 continue                                                             3476
20161 continue                                                             3477
20162 continue                                                             3477
      ga=sqrt(ga)                                                          3478
20180 do 20181 ilm=1,nlam                                                  3478
      al0=al                                                               3479
      if(flmin .lt. 1.0)goto 20201                                         3479
      al=ulam(ilm)                                                         3479
      goto 20191                                                           3480
20201 if(ilm .le. 2)goto 20211                                             3480
      al=al*alf                                                            3480
      goto 20191                                                           3481
20211 if(ilm .ne. 1)goto 20221                                             3481
      al=big                                                               3481
      goto 20231                                                           3482
20221 continue                                                             3482
      al0=0.0                                                              3483
20240 do 20241 j=1,ni                                                      3483
      if(ju(j).eq.0)goto 20241                                             3483
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3483
20241 continue                                                             3484
20242 continue                                                             3484
      al0=al0/max(bta,1.0e-3)                                              3484
      al=alf*al0                                                           3485
20231 continue                                                             3486
20191 continue                                                             3486
      al2=al*omb                                                           3486
      al1=al*bta                                                           3486
      tlam=bta*(2.0*al-al0)                                                3487
20250 do 20251 k=1,ni                                                      3487
      if(iy(k).eq.1)goto 20251                                             3487
      if(ju(k).eq.0)goto 20251                                             3488
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      3489
20251 continue                                                             3490
20252 continue                                                             3490
10680 continue                                                             3491
20260 continue                                                             3491
20261 continue                                                             3491
      ixx=0                                                                3491
      jxx=ixx                                                              3491
      kxx=jxx                                                              3491
      t=0.0                                                                3492
20270 do 20271 ic=1,nc                                                     3492
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       3492
20271 continue                                                             3493
20272 continue                                                             3493
      if(t .ge. eps)goto 20291                                             3493
      kxx=1                                                                3493
      goto 20262                                                           3493
20291 continue                                                             3493
      t=2.0*t                                                              3493
      alt=al1/t                                                            3493
      al2t=al2/t                                                           3494
20300 do 20301 ic=1,nc                                                     3494
      bs(0,ic)=b(0,ic)                                                     3494
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          3495
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    3495
      svr(ic)=sum(r(:,ic))                                                 3496
      b(0,ic)=b(0,ic)+svr(ic)                                              3496
      r(:,ic)=r(:,ic)-svr(ic)*w                                            3497
      dlx=max(dlx,svr(ic)**2)                                              3498
20301 continue                                                             3499
20302 continue                                                             3499
20310 continue                                                             3499
20311 continue                                                             3499
      nlp=nlp+nc                                                           3499
      dlx=0.0                                                              3500
20320 do 20321 k=1,ni                                                      3500
      if(iy(k).eq.0)goto 20321                                             3501
      jb=ix(k)                                                             3501
      je=ix(k+1)-1                                                         3501
      del=b(k,:)                                                           3501
      gkn=0.0                                                              3502
20330 do 20331 ic=1,nc                                                     3503
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))-svr(ic)*xb(k))/xs(k)        3504
      gk(ic)=u+del(ic)*xv(k)                                               3504
      gkn=gkn+gk(ic)**2                                                    3505
20331 continue                                                             3506
20332 continue                                                             3506
      gkn=sqrt(gkn)                                                        3506
      u=1.0-alt*vp(k)/gkn                                                  3507
      if(u .gt. 0.0)goto 20351                                             3507
      b(k,:)=0.0                                                           3507
      goto 20361                                                           3508
20351 continue                                                             3508
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     3508
20361 continue                                                             3509
20341 continue                                                             3509
      del=b(k,:)-del                                                       3509
      if(maxval(abs(del)).le.0.0)goto 20321                                3510
20370 do 20371 ic=1,nc                                                     3510
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3511
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   3513 
     *b(k))/xs(k)
20371 continue                                                             3514
20372 continue                                                             3514
      if(mm(k) .ne. 0)goto 20391                                           3514
      nin=nin+1                                                            3515
      if(nin .le. nx)goto 20411                                            3515
      jxx=1                                                                3515
      goto 20322                                                           3515
20411 continue                                                             3516
      mm(k)=nin                                                            3516
      m(nin)=k                                                             3517
20391 continue                                                             3518
20321 continue                                                             3519
20322 continue                                                             3519
      if(jxx.gt.0)goto 20312                                               3520
      if(dlx.lt.shr)goto 20312                                             3520
      if(nlp .le. maxit)goto 20431                                         3520
      jerr=-ilm                                                            3520
      return                                                               3520
20431 continue                                                             3521
20440 continue                                                             3521
20441 continue                                                             3521
      nlp=nlp+nc                                                           3521
      dlx=0.0                                                              3522
20450 do 20451 l=1,nin                                                     3522
      k=m(l)                                                               3522
      jb=ix(k)                                                             3522
      je=ix(k+1)-1                                                         3522
      del=b(k,:)                                                           3522
      gkn=0.0                                                              3523
20460 do 20461 ic=1,nc                                                     3524
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))  -svr(ic)*xb(k))/xs(k)      3526
      gk(ic)=u+del(ic)*xv(k)                                               3526
      gkn=gkn+gk(ic)**2                                                    3527
20461 continue                                                             3528
20462 continue                                                             3528
      gkn=sqrt(gkn)                                                        3528
      u=1.0-alt*vp(k)/gkn                                                  3529
      if(u .gt. 0.0)goto 20481                                             3529
      b(k,:)=0.0                                                           3529
      goto 20491                                                           3530
20481 continue                                                             3530
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     3530
20491 continue                                                             3531
20471 continue                                                             3531
      del=b(k,:)-del                                                       3531
      if(maxval(abs(del)).le.0.0)goto 20451                                3532
20500 do 20501 ic=1,nc                                                     3532
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3533
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   3535 
     *b(k))/xs(k)
20501 continue                                                             3536
20502 continue                                                             3536
20451 continue                                                             3537
20452 continue                                                             3537
      if(dlx.lt.shr)goto 20442                                             3537
      if(nlp .le. maxit)goto 20521                                         3537
      jerr=-ilm                                                            3537
      return                                                               3537
20521 continue                                                             3539
      goto 20441                                                           3540
20442 continue                                                             3540
      goto 20311                                                           3541
20312 continue                                                             3541
      if(jxx.gt.0)goto 20262                                               3542
20530 do 20531 ic=1,nc                                                     3543
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                               3544
      if(ixx .ne. 0)goto 20551                                             3545
20560 do 20561 j=1,nin                                                     3545
      k=m(j)                                                               3546
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 20581                   3546
      ixx=1                                                                3546
      goto 20562                                                           3546
20581 continue                                                             3548
20561 continue                                                             3549
20562 continue                                                             3549
20551 continue                                                             3550
      sc=b(0,ic)+g(:,ic)                                                   3550
      b0=0.0                                                               3551
20590 do 20591 j=1,nin                                                     3551
      l=m(j)                                                               3551
      jb=ix(l)                                                             3551
      je=ix(l+1)-1                                                         3552
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   3553
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            3554
20591 continue                                                             3555
20592 continue                                                             3555
      sc=min(max(exmn,sc+b0),exmx)                                         3556
      sxp=sxp-q(:,ic)                                                      3557
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          3558
      sxp=sxp+q(:,ic)                                                      3559
20531 continue                                                             3560
20532 continue                                                             3560
      s=sum(b(0,:))/nc                                                     3560
      b(0,:)=b(0,:)-s                                                      3561
      if(jxx.gt.0)goto 20262                                               3562
      if(ixx .ne. 0)goto 20611                                             3563
20620 do 20621 j=1,ni                                                      3563
      if(iy(j).eq.1)goto 20621                                             3563
      if(ju(j).eq.0)goto 20621                                             3563
      ga(j)=0.0                                                            3563
20621 continue                                                             3564
20622 continue                                                             3564
20630 do 20631 ic=1,nc                                                     3564
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3565
20640 do 20641 j=1,ni                                                      3565
      if(iy(j).eq.1)goto 20641                                             3565
      if(ju(j).eq.0)goto 20641                                             3566
      jb=ix(j)                                                             3566
      je=ix(j+1)-1                                                         3567
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             3568
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            3569
20641 continue                                                             3570
20642 continue                                                             3570
20631 continue                                                             3571
20632 continue                                                             3571
      ga=sqrt(ga)                                                          3572
20650 do 20651 k=1,ni                                                      3572
      if(iy(k).eq.1)goto 20651                                             3572
      if(ju(k).eq.0)goto 20651                                             3573
      if(ga(k) .le. al1*vp(k))goto 20671                                   3573
      iy(k)=1                                                              3573
      ixx=1                                                                3573
20671 continue                                                             3574
20651 continue                                                             3575
20652 continue                                                             3575
      if(ixx.eq.1) go to 10680                                             3576
      goto 20262                                                           3577
20611 continue                                                             3578
      goto 20261                                                           3579
20262 continue                                                             3579
      if(kxx .le. 0)goto 20691                                             3579
      jerr=-20000-ilm                                                      3579
      goto 20182                                                           3579
20691 continue                                                             3580
      if(jxx .le. 0)goto 20711                                             3580
      jerr=-10000-ilm                                                      3580
      goto 20182                                                           3580
20711 continue                                                             3580
      devi=0.0                                                             3581
20720 do 20721 ic=1,nc                                                     3582
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          3582
      a0(ic,ilm)=b(0,ic)                                                   3583
20730 do 20731 i=1,no                                                      3583
      if(y(i,ic).le.0.0)goto 20731                                         3584
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           3585
20731 continue                                                             3586
20732 continue                                                             3586
20721 continue                                                             3587
20722 continue                                                             3587
      kin(ilm)=nin                                                         3587
      alm(ilm)=al                                                          3587
      lmu=ilm                                                              3588
      dev(ilm)=(dev1-devi)/dev0                                            3589
      if(ilm.lt.mnl)goto 20181                                             3589
      if(flmin.ge.1.0)goto 20181                                           3590
      me=0                                                                 3590
20740 do 20741 j=1,nin                                                     3590
      if(a(j,1,ilm).ne.0.0) me=me+1                                        3590
20741 continue                                                             3590
20742 continue                                                             3590
      if(me.gt.ne)goto 20182                                               3591
      if(dev(ilm).gt.devmax)goto 20182                                     3591
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 20182                             3592
20181 continue                                                             3593
20182 continue                                                             3593
      g=log(q)                                                             3593
20750 do 20751 i=1,no                                                      3593
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3593
20751 continue                                                             3594
20752 continue                                                             3594
      deallocate(sxp,b,bs,r,q,mm,is,sc,ga,iy,gk,del,sxpl)                  3595
      return                                                               3596
      end                                                                  3597
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
