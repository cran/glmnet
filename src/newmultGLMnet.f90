c
c                          newGLMnet (6/9/12)
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
      subroutine elnet  (ka,parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam    718 
     *,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam)                          719
      real ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)                          720
      integer jd(*),ia(nx),nin(nlam)                                        721
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10021                                     724
      jerr=10000                                                            724
      return                                                                724
10021 continue                                                              725
      allocate(vq(1:ni),stat=jerr)                                          725
      if(jerr.ne.0) return                                                  726
      vq=max(0.0,vp)                                                        726
      vq=vq*ni/sum(vq)                                                      727
      if(ka .ne. 1)goto 10041                                               728
      call elnetu  (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd    731 
     *,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                            732
10041 continue                                                              733
      call elnetn (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd,    736 
     *maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                              737
10031 continue                                                              737
      deallocate(vq)                                                        738
      return                                                                739
      end                                                                   740
      subroutine elnetu  (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,t    743 
     *hr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ulam(nlam)                           744
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         745
      integer jd(*),ia(nx),nin(nlam)                                        746
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           751
      allocate(xm(1:ni),stat=ierr)                                          751
      jerr=jerr+ierr                                                        752
      allocate(xs(1:ni),stat=ierr)                                          752
      jerr=jerr+ierr                                                        753
      allocate(ju(1:ni),stat=ierr)                                          753
      jerr=jerr+ierr                                                        754
      allocate(xv(1:ni),stat=ierr)                                          754
      jerr=jerr+ierr                                                        755
      allocate(vlam(1:nlam),stat=ierr)                                      755
      jerr=jerr+ierr                                                        756
      if(jerr.ne.0) return                                                  757
      call chkvars(no,ni,x,ju)                                              758
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  759
      if(maxval(ju) .gt. 0)goto 10071                                       759
      jerr=7777                                                             759
      return                                                                759
10071 continue                                                              760
      call standard(no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)               761
      if(jerr.ne.0) return                                                  762
      if(flmin.ge.1.0) vlam=ulam/ys                                         763
      call elnet1(parm,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxit,x    765 
     *v,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  766
10080 do 10081 k=1,lmu                                                      766
      alm(k)=ys*alm(k)                                                      766
      nk=nin(k)                                                             767
10090 do 10091 l=1,nk                                                       767
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          767
10091 continue                                                              768
10092 continue                                                              768
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         769
10081 continue                                                              770
10082 continue                                                              770
      deallocate(xm,xs,g,ju,xv,vlam)                                        771
      return                                                                772
      end                                                                   773
      subroutine standard (no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)        774
      real x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                  774
      integer ju(ni)                                                        775
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           778
      if(jerr.ne.0) return                                                  779
      w=w/sum(w)                                                            779
      v=sqrt(w)                                                             780
10100 do 10101 j=1,ni                                                       780
      if(ju(j).eq.0)goto 10101                                              781
      xm(j)=dot_product(w,x(:,j))                                           781
      x(:,j)=v*(x(:,j)-xm(j))                                               782
      xv(j)=dot_product(x(:,j),x(:,j))                                      782
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        783
10101 continue                                                              784
10102 continue                                                              784
      if(isd .ne. 0)goto 10121                                              784
      xs=1.0                                                                784
      goto 10131                                                            785
10121 continue                                                              786
10140 do 10141 j=1,ni                                                       786
      if(ju(j).eq.0)goto 10141                                              786
      x(:,j)=x(:,j)/xs(j)                                                   786
10141 continue                                                              787
10142 continue                                                              787
      xv=1.0                                                                788
10131 continue                                                              789
10111 continue                                                              789
      ym=dot_product(w,y)                                                   789
      y=v*(y-ym)                                                            789
      ys=sqrt(dot_product(y,y))                                             789
      y=y/ys                                                                789
      g=0.0                                                                 790
10150 do 10151 j=1,ni                                                       790
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                             790
10151 continue                                                              791
10152 continue                                                              791
      deallocate(v)                                                         792
      return                                                                793
      end                                                                   794
      subroutine elnet1 (beta,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,ulam,thr,    796 
     *maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    797 
     *9)
      real vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    798 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       799
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                           805
      jerr=jerr+ierr                                                        806
      allocate(mm(1:ni),stat=ierr)                                          806
      jerr=jerr+ierr                                                        807
      allocate(da(1:ni),stat=ierr)                                          807
      jerr=jerr+ierr                                                        808
      if(jerr.ne.0) return                                                  809
      bta=beta                                                              809
      omb=1.0-bta                                                           810
      if(flmin .ge. 1.0)goto 10171                                          810
      eqs=max(eps,flmin)                                                    810
      alf=eqs**(1.0/(nlam-1))                                               810
10171 continue                                                              811
      rsq=0.0                                                               811
      a=0.0                                                                 811
      mm=0                                                                  811
      nlp=0                                                                 811
      nin=nlp                                                               811
      iz=0                                                                  811
      mnl=min(mnlam,nlam)                                                   812
10180 do 10181 m=1,nlam                                                     813
      if(flmin .lt. 1.0)goto 10201                                          813
      alm=ulam(m)                                                           813
      goto 10191                                                            814
10201 if(m .le. 2)goto 10211                                                814
      alm=alm*alf                                                           814
      goto 10191                                                            815
10211 if(m .ne. 1)goto 10221                                                815
      alm=big                                                               815
      goto 10231                                                            816
10221 continue                                                              816
      alm=0.0                                                               817
10240 do 10241 j=1,ni                                                       817
      if(ju(j).eq.0)goto 10241                                              817
      if(vp(j).le.0.0)goto 10241                                            818
      alm=max(alm,abs(g(j))/vp(j))                                          819
10241 continue                                                              820
10242 continue                                                              820
      alm=alf*alm/max(bta,1.0e-3)                                           821
10231 continue                                                              822
10191 continue                                                              822
      dem=alm*omb                                                           822
      ab=alm*bta                                                            822
      rsq0=rsq                                                              822
      jz=1                                                                  823
10250 continue                                                              823
10251 continue                                                              823
      if(iz*jz.ne.0) go to 10260                                            823
      nlp=nlp+1                                                             823
      dlx=0.0                                                               824
10270 do 10271 k=1,ni                                                       824
      if(ju(k).eq.0)goto 10271                                              825
      ak=a(k)                                                               825
      u=g(k)+ak*xv(k)                                                       825
      v=abs(u)-vp(k)*ab                                                     825
      a(k)=0.0                                                              826
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         827
      if(a(k).eq.ak)goto 10271                                              828
      if(mm(k) .ne. 0)goto 10291                                            828
      nin=nin+1                                                             828
      if(nin.gt.nx)goto 10272                                               829
10300 do 10301 j=1,ni                                                       829
      if(ju(j).eq.0)goto 10301                                              830
      if(mm(j) .eq. 0)goto 10321                                            830
      c(j,nin)=c(k,mm(j))                                                   830
      goto 10301                                                            830
10321 continue                                                              831
      if(j .ne. k)goto 10341                                                831
      c(j,nin)=xv(j)                                                        831
      goto 10301                                                            831
10341 continue                                                              832
      c(j,nin)=dot_product(x(:,j),x(:,k))                                   833
10301 continue                                                              834
10302 continue                                                              834
      mm(k)=nin                                                             834
      ia(nin)=k                                                             835
10291 continue                                                              836
      del=a(k)-ak                                                           836
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      837
      dlx=max(xv(k)*del**2,dlx)                                             838
10350 do 10351 j=1,ni                                                       838
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                               838
10351 continue                                                              839
10352 continue                                                              839
10271 continue                                                              840
10272 continue                                                              840
      if(dlx.lt.thr)goto 10252                                              840
      if(nin.gt.nx)goto 10252                                               841
      if(nlp .le. maxit)goto 10371                                          841
      jerr=-m                                                               841
      return                                                                841
10371 continue                                                              842
10260 continue                                                              842
      iz=1                                                                  842
      da(1:nin)=a(ia(1:nin))                                                843
10380 continue                                                              843
10381 continue                                                              843
      nlp=nlp+1                                                             843
      dlx=0.0                                                               844
10390 do 10391 l=1,nin                                                      844
      k=ia(l)                                                               844
      ak=a(k)                                                               844
      u=g(k)+ak*xv(k)                                                       844
      v=abs(u)-vp(k)*ab                                                     845
      a(k)=0.0                                                              846
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         847
      if(a(k).eq.ak)goto 10391                                              848
      del=a(k)-ak                                                           848
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      849
      dlx=max(xv(k)*del**2,dlx)                                             850
10400 do 10401 j=1,nin                                                      850
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                  850
10401 continue                                                              851
10402 continue                                                              851
10391 continue                                                              852
10392 continue                                                              852
      if(dlx.lt.thr)goto 10382                                              852
      if(nlp .le. maxit)goto 10421                                          852
      jerr=-m                                                               852
      return                                                                852
10421 continue                                                              853
      goto 10381                                                            854
10382 continue                                                              854
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                      855
10430 do 10431 j=1,ni                                                       855
      if(mm(j).ne.0)goto 10431                                              856
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))            857
10431 continue                                                              858
10432 continue                                                              858
      jz=0                                                                  859
      goto 10251                                                            860
10252 continue                                                              860
      if(nin .le. nx)goto 10451                                             860
      jerr=-10000-m                                                         860
      goto 10182                                                            860
10451 continue                                                              861
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 861
      kin(m)=nin                                                            862
      rsqo(m)=rsq                                                           862
      almo(m)=alm                                                           862
      lmu=m                                                                 863
      if(m.lt.mnl)goto 10181                                                863
      if(flmin.ge.1.0)goto 10181                                            864
      me=0                                                                  864
10460 do 10461 j=1,nin                                                      864
      if(ao(j,m).ne.0.0) me=me+1                                            864
10461 continue                                                              864
10462 continue                                                              864
      if(me.gt.ne)goto 10182                                                865
      if(rsq-rsq0.lt.sml*rsq)goto 10182                                     865
      if(rsq.gt.rsqmax)goto 10182                                           866
10181 continue                                                              867
10182 continue                                                              867
      deallocate(a,mm,c,da)                                                 868
      return                                                                869
      end                                                                   870
      subroutine elnetn (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,th    872 
     *r,isd,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no),w(no),ulam(nlam)                           873
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         874
      integer jd(*),ia(nx),nin(nlam)                                        875
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          880
      allocate(xs(1:ni),stat=ierr)                                          880
      jerr=jerr+ierr                                                        881
      allocate(ju(1:ni),stat=ierr)                                          881
      jerr=jerr+ierr                                                        882
      allocate(xv(1:ni),stat=ierr)                                          882
      jerr=jerr+ierr                                                        883
      allocate(vlam(1:nlam),stat=ierr)                                      883
      jerr=jerr+ierr                                                        884
      if(jerr.ne.0) return                                                  885
      call chkvars(no,ni,x,ju)                                              886
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  887
      if(maxval(ju) .gt. 0)goto 10481                                       887
      jerr=7777                                                             887
      return                                                                887
10481 continue                                                              888
      call standard1(no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)                889
      if(jerr.ne.0) return                                                  890
      if(flmin.ge.1.0) vlam=ulam/ys                                         891
      call elnet2(parm,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,vlam,thr,maxit,x    893 
     *v,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  894
10490 do 10491 k=1,lmu                                                      894
      alm(k)=ys*alm(k)                                                      894
      nk=nin(k)                                                             895
10500 do 10501 l=1,nk                                                       895
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          895
10501 continue                                                              896
10502 continue                                                              896
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         897
10491 continue                                                              898
10492 continue                                                              898
      deallocate(xm,xs,ju,xv,vlam)                                          899
      return                                                                900
      end                                                                   901
      subroutine standard1 (no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)         902
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)                        902
      integer ju(ni)                                                        903
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           906
      if(jerr.ne.0) return                                                  907
      w=w/sum(w)                                                            907
      v=sqrt(w)                                                             908
10510 do 10511 j=1,ni                                                       908
      if(ju(j).eq.0)goto 10511                                              909
      xm(j)=dot_product(w,x(:,j))                                           909
      x(:,j)=v*(x(:,j)-xm(j))                                               910
      xv(j)=dot_product(x(:,j),x(:,j))                                      910
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        911
10511 continue                                                              912
10512 continue                                                              912
      if(isd .ne. 0)goto 10531                                              912
      xs=1.0                                                                912
      goto 10541                                                            913
10531 continue                                                              913
10550 do 10551 j=1,ni                                                       913
      if(ju(j).eq.0)goto 10551                                              913
      x(:,j)=x(:,j)/xs(j)                                                   913
10551 continue                                                              914
10552 continue                                                              914
      xv=1.0                                                                915
10541 continue                                                              916
10521 continue                                                              916
      ym=dot_product(w,y)                                                   916
      y=v*(y-ym)                                                            916
      ys=sqrt(dot_product(y,y))                                             916
      y=y/ys                                                                917
      deallocate(v)                                                         918
      return                                                                919
      end                                                                   920
      subroutine elnet2(beta,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,ulam,thr,m    922 
     *axit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    923 
     *9)
      real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    924 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       925
      real, dimension (:), allocatable :: a,g                                   
      integer, dimension (:), allocatable :: mm,ix                              
      allocate(a(1:ni),stat=jerr)                                           930
      allocate(mm(1:ni),stat=ierr)                                          930
      jerr=jerr+ierr                                                        931
      allocate(g(1:ni),stat=ierr)                                           931
      jerr=jerr+ierr                                                        932
      allocate(ix(1:ni),stat=ierr)                                          932
      jerr=jerr+ierr                                                        933
      if(jerr.ne.0) return                                                  934
      bta=beta                                                              934
      omb=1.0-bta                                                           934
      ix=0                                                                  935
      if(flmin .ge. 1.0)goto 10571                                          935
      eqs=max(eps,flmin)                                                    935
      alf=eqs**(1.0/(nlam-1))                                               935
10571 continue                                                              936
      rsq=0.0                                                               936
      a=0.0                                                                 936
      mm=0                                                                  936
      nlp=0                                                                 936
      nin=nlp                                                               936
      iz=0                                                                  936
      mnl=min(mnlam,nlam)                                                   936
      alm=0.0                                                               937
10580 do 10581 j=1,ni                                                       937
      if(ju(j).eq.0)goto 10581                                              937
      g(j)=abs(dot_product(y,x(:,j)))                                       937
10581 continue                                                              938
10582 continue                                                              938
10590 do 10591 m=1,nlam                                                     938
      alm0=alm                                                              939
      if(flmin .lt. 1.0)goto 10611                                          939
      alm=ulam(m)                                                           939
      goto 10601                                                            940
10611 if(m .le. 2)goto 10621                                                940
      alm=alm*alf                                                           940
      goto 10601                                                            941
10621 if(m .ne. 1)goto 10631                                                941
      alm=big                                                               941
      goto 10641                                                            942
10631 continue                                                              942
      alm0=0.0                                                              943
10650 do 10651 j=1,ni                                                       943
      if(ju(j).eq.0)goto 10651                                              943
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                            943
10651 continue                                                              944
10652 continue                                                              944
      alm0=alm0/max(bta,1.0e-3)                                             944
      alm=alf*alm0                                                          945
10641 continue                                                              946
10601 continue                                                              946
      dem=alm*omb                                                           946
      ab=alm*bta                                                            946
      rsq0=rsq                                                              946
      jz=1                                                                  947
      tlam=bta*(2.0*alm-alm0)                                               948
10660 do 10661 k=1,ni                                                       948
      if(ix(k).eq.1)goto 10661                                              948
      if(ju(k).eq.0)goto 10661                                              949
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                        950
10661 continue                                                              951
10662 continue                                                              951
10670 continue                                                              951
10671 continue                                                              951
      if(iz*jz.ne.0) go to 10260                                            952
10680 continue                                                              952
      nlp=nlp+1                                                             952
      dlx=0.0                                                               953
10690 do 10691 k=1,ni                                                       953
      if(ix(k).eq.0)goto 10691                                              953
      gk=dot_product(y,x(:,k))                                              954
      ak=a(k)                                                               954
      u=gk+ak*xv(k)                                                         954
      v=abs(u)-vp(k)*ab                                                     954
      a(k)=0.0                                                              955
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         956
      if(a(k).eq.ak)goto 10691                                              957
      if(mm(k) .ne. 0)goto 10711                                            957
      nin=nin+1                                                             957
      if(nin.gt.nx)goto 10692                                               958
      mm(k)=nin                                                             958
      ia(nin)=k                                                             959
10711 continue                                                              960
      del=a(k)-ak                                                           960
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        961
      y=y-del*x(:,k)                                                        961
      dlx=max(xv(k)*del**2,dlx)                                             962
10691 continue                                                              963
10692 continue                                                              963
      if(nin.gt.nx)goto 10672                                               964
      if(dlx .ge. thr)goto 10731                                            964
      ixx=0                                                                 965
10740 do 10741 k=1,ni                                                       965
      if(ix(k).eq.1)goto 10741                                              965
      if(ju(k).eq.0)goto 10741                                              966
      g(k)=abs(dot_product(y,x(:,k)))                                       967
      if(g(k) .le. ab*vp(k))goto 10761                                      967
      ix(k)=1                                                               967
      ixx=1                                                                 967
10761 continue                                                              968
10741 continue                                                              969
10742 continue                                                              969
      if(ixx.eq.1) go to 10680                                              970
      goto 10672                                                            971
10731 continue                                                              972
      if(nlp .le. maxit)goto 10781                                          972
      jerr=-m                                                               972
      return                                                                972
10781 continue                                                              973
10260 continue                                                              973
      iz=1                                                                  974
10790 continue                                                              974
10791 continue                                                              974
      nlp=nlp+1                                                             974
      dlx=0.0                                                               975
10800 do 10801 l=1,nin                                                      975
      k=ia(l)                                                               975
      gk=dot_product(y,x(:,k))                                              976
      ak=a(k)                                                               976
      u=gk+ak*xv(k)                                                         976
      v=abs(u)-vp(k)*ab                                                     976
      a(k)=0.0                                                              977
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         978
      if(a(k).eq.ak)goto 10801                                              979
      del=a(k)-ak                                                           979
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        980
      y=y-del*x(:,k)                                                        980
      dlx=max(xv(k)*del**2,dlx)                                             981
10801 continue                                                              982
10802 continue                                                              982
      if(dlx.lt.thr)goto 10792                                              982
      if(nlp .le. maxit)goto 10821                                          982
      jerr=-m                                                               982
      return                                                                982
10821 continue                                                              983
      goto 10791                                                            984
10792 continue                                                              984
      jz=0                                                                  985
      goto 10671                                                            986
10672 continue                                                              986
      if(nin .le. nx)goto 10841                                             986
      jerr=-10000-m                                                         986
      goto 10592                                                            986
10841 continue                                                              987
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 987
      kin(m)=nin                                                            988
      rsqo(m)=rsq                                                           988
      almo(m)=alm                                                           988
      lmu=m                                                                 989
      if(m.lt.mnl)goto 10591                                                989
      if(flmin.ge.1.0)goto 10591                                            990
      me=0                                                                  990
10850 do 10851 j=1,nin                                                      990
      if(ao(j,m).ne.0.0) me=me+1                                            990
10851 continue                                                              990
10852 continue                                                              990
      if(me.gt.ne)goto 10592                                                991
      if(rsq-rsq0.lt.sml*rsq)goto 10592                                     991
      if(rsq.gt.rsqmax)goto 10592                                           992
10591 continue                                                              993
10592 continue                                                              993
      deallocate(a,mm,g,ix)                                                 994
      return                                                                995
      end                                                                   996
      subroutine chkvars(no,ni,x,ju)                                        997
      real x(no,ni)                                                         997
      integer ju(ni)                                                        998
10860 do 10861 j=1,ni                                                       998
      ju(j)=0                                                               998
      t=x(1,j)                                                              999
10870 do 10871 i=2,no                                                       999
      if(x(i,j).eq.t)goto 10871                                             999
      ju(j)=1                                                               999
      goto 10872                                                            999
10871 continue                                                             1000
10872 continue                                                             1000
10861 continue                                                             1001
10862 continue                                                             1001
      return                                                               1002
      end                                                                  1003
      subroutine uncomp(ni,ca,ia,nin,a)                                    1004
      real ca(*),a(ni)                                                     1004
      integer ia(*)                                                        1005
      a=0.0                                                                1005
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                  1006
      return                                                               1007
      end                                                                  1008
      subroutine modval(a0,ca,ia,nin,n,x,f)                                1009
      real ca(nin),x(n,*),f(n)                                             1009
      integer ia(nin)                                                      1010
      f=a0                                                                 1010
      if(nin.le.0) return                                                  1011
10880 do 10881 i=1,n                                                       1011
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      1011
10881 continue                                                             1012
10882 continue                                                             1012
      return                                                               1013
      end                                                                  1014
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,fl   1017 
     *min,ulam,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                              1018
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1019
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1020
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10901                                    1023
      jerr=10000                                                           1023
      return                                                               1023
10901 continue                                                             1024
      allocate(vq(1:ni),stat=jerr)                                         1024
      if(jerr.ne.0) return                                                 1025
      vq=max(0.0,vp)                                                       1025
      vq=vq*ni/sum(vq)                                                     1026
      if(ka .ne. 1)goto 10921                                              1027
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam   1030 
     *,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10931                                                           1031
10921 continue                                                             1032
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam,   1035 
     *thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10931 continue                                                             1036
10911 continue                                                             1036
      deallocate(vq)                                                       1037
      return                                                               1038
      end                                                                  1039
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmi   1042 
     *n,ulam,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                              1043
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1044
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1045
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                          1050
      allocate(xm(1:ni),stat=ierr)                                         1050
      jerr=jerr+ierr                                                       1051
      allocate(xs(1:ni),stat=ierr)                                         1051
      jerr=jerr+ierr                                                       1052
      allocate(ju(1:ni),stat=ierr)                                         1052
      jerr=jerr+ierr                                                       1053
      allocate(xv(1:ni),stat=ierr)                                         1053
      jerr=jerr+ierr                                                       1054
      allocate(vlam(1:nlam),stat=ierr)                                     1054
      jerr=jerr+ierr                                                       1055
      if(jerr.ne.0) return                                                 1056
      call spchkvars(no,ni,x,ix,ju)                                        1057
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1058
      if(maxval(ju) .gt. 0)goto 10951                                      1058
      jerr=7777                                                            1058
      return                                                               1058
10951 continue                                                             1059
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,jerr)      1060
      if(jerr.ne.0) return                                                 1061
      if(flmin.ge.1.0) vlam=ulam/ys                                        1062
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t   1064 
     *hr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1065
10960 do 10961 k=1,lmu                                                     1065
      alm(k)=ys*alm(k)                                                     1065
      nk=nin(k)                                                            1066
10970 do 10971 l=1,nk                                                      1066
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1066
10971 continue                                                             1067
10972 continue                                                             1067
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                        1068
10961 continue                                                             1069
10962 continue                                                             1069
      deallocate(xm,xs,g,ju,xv,vlam)                                       1070
      return                                                               1071
      end                                                                  1072
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,j   1073 
     *err)
      real x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                     1073
      integer ix(*),jx(*),ju(ni)                                           1074
      w=w/sum(w)                                                           1075
10980 do 10981 j=1,ni                                                      1075
      if(ju(j).eq.0)goto 10981                                             1076
      jb=ix(j)                                                             1076
      je=ix(j+1)-1                                                         1076
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1077
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1078
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1079
10981 continue                                                             1080
10982 continue                                                             1080
      if(isd .ne. 0)goto 11001                                             1080
      xs=1.0                                                               1080
      goto 11011                                                           1080
11001 continue                                                             1080
      xv=1.0                                                               1080
11011 continue                                                             1081
10991 continue                                                             1081
      ym=dot_product(w,y)                                                  1081
      y=y-ym                                                               1081
      ys=sqrt(dot_product(w,y**2))                                         1081
      y=y/ys                                                               1081
      g=0.0                                                                1082
11020 do 11021 j=1,ni                                                      1082
      if(ju(j).eq.0)goto 11021                                             1082
      jb=ix(j)                                                             1082
      je=ix(j+1)-1                                                         1083
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)           1084
11021 continue                                                             1085
11022 continue                                                             1085
      return                                                               1086
      end                                                                  1087
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,   1089 
     *ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99   1090 
     *9)
      real g(ni),vp(ni),x(*),ulam(nlam),w(no)                              1091
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)          1092
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1093
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                          1099
      jerr=jerr+ierr                                                       1100
      allocate(mm(1:ni),stat=ierr)                                         1100
      jerr=jerr+ierr                                                       1101
      allocate(da(1:ni),stat=ierr)                                         1101
      jerr=jerr+ierr                                                       1102
      if(jerr.ne.0) return                                                 1103
      bta=beta                                                             1103
      omb=1.0-bta                                                          1104
      if(flmin .ge. 1.0)goto 11041                                         1104
      eqs=max(eps,flmin)                                                   1104
      alf=eqs**(1.0/(nlam-1))                                              1104
11041 continue                                                             1105
      rsq=0.0                                                              1105
      a=0.0                                                                1105
      mm=0                                                                 1105
      nlp=0                                                                1105
      nin=nlp                                                              1105
      iz=0                                                                 1105
      mnl=min(mnlam,nlam)                                                  1106
11050 do 11051 m=1,nlam                                                    1107
      if(flmin .lt. 1.0)goto 11071                                         1107
      alm=ulam(m)                                                          1107
      goto 11061                                                           1108
11071 if(m .le. 2)goto 11081                                               1108
      alm=alm*alf                                                          1108
      goto 11061                                                           1109
11081 if(m .ne. 1)goto 11091                                               1109
      alm=big                                                              1109
      goto 11101                                                           1110
11091 continue                                                             1110
      alm=0.0                                                              1111
11110 do 11111 j=1,ni                                                      1111
      if(ju(j).eq.0)goto 11111                                             1111
      if(vp(j).le.0.0)goto 11111                                           1112
      alm=max(alm,abs(g(j))/vp(j))                                         1113
11111 continue                                                             1114
11112 continue                                                             1114
      alm=alf*alm/max(bta,1.0e-3)                                          1115
11101 continue                                                             1116
11061 continue                                                             1116
      dem=alm*omb                                                          1116
      ab=alm*bta                                                           1116
      rsq0=rsq                                                             1116
      jz=1                                                                 1117
11120 continue                                                             1117
11121 continue                                                             1117
      if(iz*jz.ne.0) go to 10260                                           1117
      nlp=nlp+1                                                            1117
      dlx=0.0                                                              1118
11130 do 11131 k=1,ni                                                      1118
      if(ju(k).eq.0)goto 11131                                             1119
      ak=a(k)                                                              1119
      u=g(k)+ak*xv(k)                                                      1119
      v=abs(u)-vp(k)*ab                                                    1119
      a(k)=0.0                                                             1120
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1121
      if(a(k).eq.ak)goto 11131                                             1122
      if(mm(k) .ne. 0)goto 11151                                           1122
      nin=nin+1                                                            1122
      if(nin.gt.nx)goto 11132                                              1123
11160 do 11161 j=1,ni                                                      1123
      if(ju(j).eq.0)goto 11161                                             1124
      if(mm(j) .eq. 0)goto 11181                                           1124
      c(j,nin)=c(k,mm(j))                                                  1124
      goto 11161                                                           1124
11181 continue                                                             1125
      if(j .ne. k)goto 11201                                               1125
      c(j,nin)=xv(j)                                                       1125
      goto 11161                                                           1125
11201 continue                                                             1126
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))       1128
11161 continue                                                             1129
11162 continue                                                             1129
      mm(k)=nin                                                            1129
      ia(nin)=k                                                            1130
11151 continue                                                             1131
      del=a(k)-ak                                                          1131
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1132
      dlx=max(xv(k)*del**2,dlx)                                            1133
11210 do 11211 j=1,ni                                                      1133
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                              1133
11211 continue                                                             1134
11212 continue                                                             1134
11131 continue                                                             1135
11132 continue                                                             1135
      if(dlx.lt.thr)goto 11122                                             1135
      if(nin.gt.nx)goto 11122                                              1136
      if(nlp .le. maxit)goto 11231                                         1136
      jerr=-m                                                              1136
      return                                                               1136
11231 continue                                                             1137
10260 continue                                                             1137
      iz=1                                                                 1137
      da(1:nin)=a(ia(1:nin))                                               1138
11240 continue                                                             1138
11241 continue                                                             1138
      nlp=nlp+1                                                            1138
      dlx=0.0                                                              1139
11250 do 11251 l=1,nin                                                     1139
      k=ia(l)                                                              1140
      ak=a(k)                                                              1140
      u=g(k)+ak*xv(k)                                                      1140
      v=abs(u)-vp(k)*ab                                                    1140
      a(k)=0.0                                                             1141
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1142
      if(a(k).eq.ak)goto 11251                                             1143
      del=a(k)-ak                                                          1143
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1144
      dlx=max(xv(k)*del**2,dlx)                                            1145
11260 do 11261 j=1,nin                                                     1145
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                 1145
11261 continue                                                             1146
11262 continue                                                             1146
11251 continue                                                             1147
11252 continue                                                             1147
      if(dlx.lt.thr)goto 11242                                             1147
      if(nlp .le. maxit)goto 11281                                         1147
      jerr=-m                                                              1147
      return                                                               1147
11281 continue                                                             1148
      goto 11241                                                           1149
11242 continue                                                             1149
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                     1150
11290 do 11291 j=1,ni                                                      1150
      if(mm(j).ne.0)goto 11291                                             1151
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))           1152
11291 continue                                                             1153
11292 continue                                                             1153
      jz=0                                                                 1154
      goto 11121                                                           1155
11122 continue                                                             1155
      if(nin .le. nx)goto 11311                                            1155
      jerr=-10000-m                                                        1155
      goto 11052                                                           1155
11311 continue                                                             1156
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1156
      kin(m)=nin                                                           1157
      rsqo(m)=rsq                                                          1157
      almo(m)=alm                                                          1157
      lmu=m                                                                1158
      if(m.lt.mnl)goto 11051                                               1158
      if(flmin.ge.1.0)goto 11051                                           1159
      me=0                                                                 1159
11320 do 11321 j=1,nin                                                     1159
      if(ao(j,m).ne.0.0) me=me+1                                           1159
11321 continue                                                             1159
11322 continue                                                             1159
      if(me.gt.ne)goto 11052                                               1160
      if(rsq-rsq0.lt.sml*rsq)goto 11052                                    1160
      if(rsq.gt.rsqmax)goto 11052                                          1161
11051 continue                                                             1162
11052 continue                                                             1162
      deallocate(a,mm,c,da)                                                1163
      return                                                               1164
      end                                                                  1165
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmin,   1167 
     *ulam,  thr,isd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),vp(ni),y(no),w(no),ulam(nlam)                              1168
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1169
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1170
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         1175
      allocate(xs(1:ni),stat=ierr)                                         1175
      jerr=jerr+ierr                                                       1176
      allocate(ju(1:ni),stat=ierr)                                         1176
      jerr=jerr+ierr                                                       1177
      allocate(xv(1:ni),stat=ierr)                                         1177
      jerr=jerr+ierr                                                       1178
      allocate(vlam(1:nlam),stat=ierr)                                     1178
      jerr=jerr+ierr                                                       1179
      if(jerr.ne.0) return                                                 1180
      call spchkvars(no,ni,x,ix,ju)                                        1181
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1182
      if(maxval(ju) .gt. 0)goto 11341                                      1182
      jerr=7777                                                            1182
      return                                                               1182
11341 continue                                                             1183
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,jerr)       1184
      if(jerr.ne.0) return                                                 1185
      if(flmin.ge.1.0) vlam=ulam/ys                                        1186
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t   1188 
     *hr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1189
11350 do 11351 k=1,lmu                                                     1189
      alm(k)=ys*alm(k)                                                     1189
      nk=nin(k)                                                            1190
11360 do 11361 l=1,nk                                                      1190
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1190
11361 continue                                                             1191
11362 continue                                                             1191
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                        1192
11351 continue                                                             1193
11352 continue                                                             1193
      deallocate(xm,xs,ju,xv,vlam)                                         1194
      return                                                               1195
      end                                                                  1196
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,je   1197 
     *rr)
      real x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)                           1197
      integer ix(*),jx(*),ju(ni)                                           1198
      w=w/sum(w)                                                           1199
11370 do 11371 j=1,ni                                                      1199
      if(ju(j).eq.0)goto 11371                                             1200
      jb=ix(j)                                                             1200
      je=ix(j+1)-1                                                         1200
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1201
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1202
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1203
11371 continue                                                             1204
11372 continue                                                             1204
      if(isd .ne. 0)goto 11391                                             1204
      xs=1.0                                                               1204
      goto 11401                                                           1204
11391 continue                                                             1204
      xv=1.0                                                               1204
11401 continue                                                             1205
11381 continue                                                             1205
      ym=dot_product(w,y)                                                  1205
      y=y-ym                                                               1205
      ys=sqrt(dot_product(w,y**2))                                         1205
      y=y/ys                                                               1206
      return                                                               1207
      end                                                                  1208
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,   1210 
     *ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99   1211 
     *9)
      real y(no),w(no),x(*),vp(ni),ulam(nlam)                              1212
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)          1213
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1214
      real, dimension (:), allocatable :: a,g                                   
      integer, dimension (:), allocatable :: mm,iy                              
      allocate(a(1:ni),stat=jerr)                                          1219
      allocate(mm(1:ni),stat=ierr)                                         1219
      jerr=jerr+ierr                                                       1220
      allocate(g(1:ni),stat=ierr)                                          1220
      jerr=jerr+ierr                                                       1221
      allocate(iy(1:ni),stat=ierr)                                         1221
      jerr=jerr+ierr                                                       1222
      if(jerr.ne.0) return                                                 1223
      bta=beta                                                             1223
      omb=1.0-bta                                                          1223
      alm=0.0                                                              1223
      iy=0                                                                 1224
      if(flmin .ge. 1.0)goto 11421                                         1224
      eqs=max(eps,flmin)                                                   1224
      alf=eqs**(1.0/(nlam-1))                                              1224
11421 continue                                                             1225
      rsq=0.0                                                              1225
      a=0.0                                                                1225
      mm=0                                                                 1225
      o=0.0                                                                1225
      nlp=0                                                                1225
      nin=nlp                                                              1225
      iz=0                                                                 1225
      mnl=min(mnlam,nlam)                                                  1226
11430 do 11431 j=1,ni                                                      1226
      if(ju(j).eq.0)goto 11431                                             1227
      jb=ix(j)                                                             1227
      je=ix(j+1)-1                                                         1228
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1229
11431 continue                                                             1230
11432 continue                                                             1230
11440 do 11441 m=1,nlam                                                    1230
      alm0=alm                                                             1231
      if(flmin .lt. 1.0)goto 11461                                         1231
      alm=ulam(m)                                                          1231
      goto 11451                                                           1232
11461 if(m .le. 2)goto 11471                                               1232
      alm=alm*alf                                                          1232
      goto 11451                                                           1233
11471 if(m .ne. 1)goto 11481                                               1233
      alm=big                                                              1233
      goto 11491                                                           1234
11481 continue                                                             1234
      alm0=0.0                                                             1235
11500 do 11501 j=1,ni                                                      1235
      if(ju(j).eq.0)goto 11501                                             1235
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           1235
11501 continue                                                             1236
11502 continue                                                             1236
      alm0=alm0/max(bta,1.0e-3)                                            1236
      alm=alf*alm0                                                         1237
11491 continue                                                             1238
11451 continue                                                             1238
      dem=alm*omb                                                          1238
      ab=alm*bta                                                           1238
      rsq0=rsq                                                             1238
      jz=1                                                                 1239
      tlam=bta*(2.0*alm-alm0)                                              1240
11510 do 11511 k=1,ni                                                      1240
      if(iy(k).eq.1)goto 11511                                             1240
      if(ju(k).eq.0)goto 11511                                             1241
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       1242
11511 continue                                                             1243
11512 continue                                                             1243
11520 continue                                                             1243
11521 continue                                                             1243
      if(iz*jz.ne.0) go to 10260                                           1244
10680 continue                                                             1244
      nlp=nlp+1                                                            1244
      dlx=0.0                                                              1245
11530 do 11531 k=1,ni                                                      1245
      if(iy(k).eq.0)goto 11531                                             1245
      jb=ix(k)                                                             1245
      je=ix(k+1)-1                                                         1246
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1247
      ak=a(k)                                                              1247
      u=gk+ak*xv(k)                                                        1247
      v=abs(u)-vp(k)*ab                                                    1247
      a(k)=0.0                                                             1248
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1249
      if(a(k).eq.ak)goto 11531                                             1250
      if(mm(k) .ne. 0)goto 11551                                           1250
      nin=nin+1                                                            1250
      if(nin.gt.nx)goto 11532                                              1251
      mm(k)=nin                                                            1251
      ia(nin)=k                                                            1252
11551 continue                                                             1253
      del=a(k)-ak                                                          1253
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1254
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1255
      o=o+del*xm(k)/xs(k)                                                  1255
      dlx=max(xv(k)*del**2,dlx)                                            1256
11531 continue                                                             1257
11532 continue                                                             1257
      if(nin.gt.nx)goto 11522                                              1258
      if(dlx .ge. thr)goto 11571                                           1258
      ixx=0                                                                1259
11580 do 11581 j=1,ni                                                      1259
      if(iy(j).eq.1)goto 11581                                             1259
      if(ju(j).eq.0)goto 11581                                             1260
      jb=ix(j)                                                             1260
      je=ix(j+1)-1                                                         1261
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1262
      if(g(j) .le. ab*vp(j))goto 11601                                     1262
      iy(j)=1                                                              1262
      ixx=1                                                                1262
11601 continue                                                             1263
11581 continue                                                             1264
11582 continue                                                             1264
      if(ixx.eq.1) go to 10680                                             1265
      goto 11522                                                           1266
11571 continue                                                             1267
      if(nlp .le. maxit)goto 11621                                         1267
      jerr=-m                                                              1267
      return                                                               1267
11621 continue                                                             1268
10260 continue                                                             1268
      iz=1                                                                 1269
11630 continue                                                             1269
11631 continue                                                             1269
      nlp=nlp+1                                                            1269
      dlx=0.0                                                              1270
11640 do 11641 l=1,nin                                                     1270
      k=ia(l)                                                              1270
      jb=ix(k)                                                             1270
      je=ix(k+1)-1                                                         1271
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1272
      ak=a(k)                                                              1272
      u=gk+ak*xv(k)                                                        1272
      v=abs(u)-vp(k)*ab                                                    1272
      a(k)=0.0                                                             1273
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1274
      if(a(k).eq.ak)goto 11641                                             1275
      del=a(k)-ak                                                          1275
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1276
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1277
      o=o+del*xm(k)/xs(k)                                                  1277
      dlx=max(xv(k)*del**2,dlx)                                            1278
11641 continue                                                             1279
11642 continue                                                             1279
      if(dlx.lt.thr)goto 11632                                             1279
      if(nlp .le. maxit)goto 11661                                         1279
      jerr=-m                                                              1279
      return                                                               1279
11661 continue                                                             1280
      goto 11631                                                           1281
11632 continue                                                             1281
      jz=0                                                                 1282
      goto 11521                                                           1283
11522 continue                                                             1283
      if(nin .le. nx)goto 11681                                            1283
      jerr=-10000-m                                                        1283
      goto 11442                                                           1283
11681 continue                                                             1284
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1284
      kin(m)=nin                                                           1285
      rsqo(m)=rsq                                                          1285
      almo(m)=alm                                                          1285
      lmu=m                                                                1286
      if(m.lt.mnl)goto 11441                                               1286
      if(flmin.ge.1.0)goto 11441                                           1287
      me=0                                                                 1287
11690 do 11691 j=1,nin                                                     1287
      if(ao(j,m).ne.0.0) me=me+1                                           1287
11691 continue                                                             1287
11692 continue                                                             1287
      if(me.gt.ne)goto 11442                                               1288
      if(rsq-rsq0.lt.sml*rsq)goto 11442                                    1288
      if(rsq.gt.rsqmax)goto 11442                                          1289
11441 continue                                                             1290
11442 continue                                                             1290
      deallocate(a,mm,g,iy)                                                1291
      return                                                               1292
      end                                                                  1293
      subroutine spchkvars(no,ni,x,ix,ju)                                  1294
      real x(*)                                                            1294
      integer ix(*),ju(ni)                                                 1295
11700 do 11701 j=1,ni                                                      1295
      ju(j)=0                                                              1295
      jb=ix(j)                                                             1295
      nj=ix(j+1)-jb                                                        1295
      if(nj.eq.0)goto 11701                                                1296
      je=ix(j+1)-1                                                         1297
      if(nj .ge. no)goto 11721                                             1297
11730 do 11731 i=jb,je                                                     1297
      if(x(i).eq.0.0)goto 11731                                            1297
      ju(j)=1                                                              1297
      goto 11732                                                           1297
11731 continue                                                             1297
11732 continue                                                             1297
      goto 11741                                                           1298
11721 continue                                                             1298
      t=x(jb)                                                              1298
11750 do 11751 i=jb+1,je                                                   1298
      if(x(i).eq.t)goto 11751                                              1298
      ju(j)=1                                                              1298
      goto 11752                                                           1298
11751 continue                                                             1298
11752 continue                                                             1298
11741 continue                                                             1299
11711 continue                                                             1299
11701 continue                                                             1300
11702 continue                                                             1300
      return                                                               1301
      end                                                                  1302
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                         1303
      real ca(*),x(*),f(n)                                                 1303
      integer ia(*),ix(*),jx(*)                                            1304
      f=a0                                                                 1305
11760 do 11761 j=1,nin                                                     1305
      k=ia(j)                                                              1305
      kb=ix(k)                                                             1305
      ke=ix(k+1)-1                                                         1306
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                             1307
11761 continue                                                             1308
11762 continue                                                             1308
      return                                                               1309
      end                                                                  1310
      function row_prod(i,j,ia,ja,ra,w)                                    1311
      integer ia(*),ja(*)                                                  1311
      real ra(*),w(*)                                                      1312
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(   1314 
     *i),ia(j+1)-ia(j),w)
      return                                                               1315
      end                                                                  1316
      function dot(x,y,mx,my,nx,ny,w)                                      1317
      real x(*),y(*),w(*)                                                  1317
      integer mx(*),my(*)                                                  1318
      i=1                                                                  1318
      j=i                                                                  1318
      s=0.0                                                                1319
11770 continue                                                             1319
11771 continue                                                             1319
11780 continue                                                             1320
11781 if(mx(i).ge.my(j))goto 11782                                         1320
      i=i+1                                                                1320
      if(i.gt.nx) go to 11790                                              1320
      goto 11781                                                           1321
11782 continue                                                             1321
      if(mx(i).eq.my(j)) go to 11800                                       1322
11810 continue                                                             1322
11811 if(my(j).ge.mx(i))goto 11812                                         1322
      j=j+1                                                                1322
      if(j.gt.ny) go to 11790                                              1322
      goto 11811                                                           1323
11812 continue                                                             1323
      if(mx(i).eq.my(j)) go to 11800                                       1323
      goto 11771                                                           1324
11800 continue                                                             1324
      s=s+w(mx(i))*x(i)*y(j)                                               1325
      i=i+1                                                                1325
      if(i.gt.nx)goto 11772                                                1325
      j=j+1                                                                1325
      if(j.gt.ny)goto 11772                                                1326
      goto 11771                                                           1327
11772 continue                                                             1327
11790 continue                                                             1327
      dot=s                                                                1328
      return                                                               1329
      end                                                                  1330
      subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,ne,nx,nlam,flmin,ulam   1332 
     *,thr,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)             1333
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1334
      integer jd(*),ia(nx),nin(nlam)                                       1335
      real, dimension (:), allocatable :: xm,xs,ww,vq,xv                        
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 11831                                    1339
      jerr=10000                                                           1339
      return                                                               1339
11831 continue                                                             1340
      allocate(ww(1:no),stat=jerr)                                         1341
      allocate(ju(1:ni),stat=ierr)                                         1341
      jerr=jerr+ierr                                                       1342
      allocate(vq(1:ni),stat=ierr)                                         1342
      jerr=jerr+ierr                                                       1343
      allocate(xm(1:ni),stat=ierr)                                         1343
      jerr=jerr+ierr                                                       1344
      if(kopt .ne. 2)goto 11851                                            1344
      allocate(xv(1:ni),stat=ierr)                                         1344
      jerr=jerr+ierr                                                       1344
11851 continue                                                             1345
      if(isd .le. 0)goto 11871                                             1345
      allocate(xs(1:ni),stat=ierr)                                         1345
      jerr=jerr+ierr                                                       1345
11871 continue                                                             1346
      if(jerr.ne.0) return                                                 1347
      call chkvars(no,ni,x,ju)                                             1348
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1349
      if(maxval(ju) .gt. 0)goto 11891                                      1349
      jerr=7777                                                            1349
      return                                                               1349
11891 continue                                                             1350
      vq=max(0.0,vp)                                                       1350
      vq=vq*ni/sum(vq)                                                     1351
11900 do 11901 i=1,no                                                      1351
      ww(i)=sum(y(i,:))                                                    1351
      y(i,:)=y(i,:)/ww(i)                                                  1351
11901 continue                                                             1351
11902 continue                                                             1351
      sw=sum(ww)                                                           1351
      ww=ww/sw                                                             1352
      if(nc .ne. 1)goto 11921                                              1352
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             1353
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,ne,nx,nlam,flmin   1355 
     *,ulam,  thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 11911                                                           1356
11921 if(kopt .ne. 2)goto 11931                                            1356
      call multlstandard1(no,ni,x,ww,ju,isd,xm,xs,xv)                      1357
      call multlognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ula   1359 
     *m,thr,  isd,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 11941                                                           1360
11931 continue                                                             1360
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             1361
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,th   1363 
     *r,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
11941 continue                                                             1364
11911 continue                                                             1364
      if(jerr.gt.0) return                                                 1364
      dev0=2.0*sw*dev0                                                     1365
11950 do 11951 k=1,lmu                                                     1365
      nk=nin(k)                                                            1366
11960 do 11961 ic=1,nc                                                     1366
      if(isd .le. 0)goto 11981                                             1366
11990 do 11991 l=1,nk                                                      1366
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1366
11991 continue                                                             1366
11992 continue                                                             1366
11981 continue                                                             1367
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1368
11961 continue                                                             1369
11962 continue                                                             1369
11951 continue                                                             1370
11952 continue                                                             1370
      deallocate(ww,ju,vq,xm)                                              1370
      if(isd.gt.0) deallocate(xs)                                          1371
      if(kopt.eq.2) deallocate(xv)                                         1372
      return                                                               1373
      end                                                                  1374
      subroutine lstandard1 (no,ni,x,w,ju,isd,xm,xs)                       1375
      real x(no,ni),w(no),xm(ni),xs(ni)                                    1375
      integer ju(ni)                                                       1376
12000 do 12001 j=1,ni                                                      1376
      if(ju(j).eq.0)goto 12001                                             1377
      xm(j)=dot_product(w,x(:,j))                                          1377
      x(:,j)=x(:,j)-xm(j)                                                  1378
      if(isd .le. 0)goto 12021                                             1378
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 1378
      x(:,j)=x(:,j)/xs(j)                                                  1378
12021 continue                                                             1379
12001 continue                                                             1380
12002 continue                                                             1380
      return                                                               1381
      end                                                                  1382
      subroutine multlstandard1 (no,ni,x,w,ju,isd,xm,xs,xv)                1383
      real x(no,ni),w(no),xm(ni),xs(ni),xv(ni)                             1383
      integer ju(ni)                                                       1384
12030 do 12031 j=1,ni                                                      1384
      if(ju(j).eq.0)goto 12031                                             1385
      xm(j)=dot_product(w,x(:,j))                                          1385
      x(:,j)=x(:,j)-xm(j)                                                  1386
      xv(j)=dot_product(w,x(:,j)**2)                                       1387
      if(isd .le. 0)goto 12051                                             1387
      xs(j)=sqrt(xv(j))                                                    1387
      x(:,j)=x(:,j)/xs(j)                                                  1387
      xv(j)=1.0                                                            1387
12051 continue                                                             1388
12031 continue                                                             1389
12032 continue                                                             1389
      return                                                               1390
      end                                                                  1391
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ulam   1393 
     *,shri,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1395 
     *5, devmax=0.999)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    1396
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1397
      integer ju(ni),m(nx),kin(nlam)                                       1398
      real, dimension (:), allocatable :: b,bs,v,r,xv,q,ga                      
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(b(0:ni),stat=jerr)                                          1403
      allocate(xv(1:ni),stat=ierr)                                         1403
      jerr=jerr+ierr                                                       1404
      allocate(ga(1:ni),stat=ierr)                                         1404
      jerr=jerr+ierr                                                       1405
      allocate(bs(0:ni),stat=ierr)                                         1405
      jerr=jerr+ierr                                                       1406
      allocate(mm(1:ni),stat=ierr)                                         1406
      jerr=jerr+ierr                                                       1407
      allocate(ixx(1:ni),stat=ierr)                                        1407
      jerr=jerr+ierr                                                       1408
      allocate(r(1:no),stat=ierr)                                          1408
      jerr=jerr+ierr                                                       1409
      allocate(v(1:no),stat=ierr)                                          1409
      jerr=jerr+ierr                                                       1410
      allocate(q(1:no),stat=ierr)                                          1410
      jerr=jerr+ierr                                                       1411
      if(jerr.ne.0) return                                                 1412
      fmax=log(1.0/pmin-1.0)                                               1412
      fmin=-fmax                                                           1412
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1413
      bta=parm                                                             1413
      omb=1.0-bta                                                          1414
      q0=dot_product(w,y)                                                  1414
      if(q0 .gt. pmin)goto 12071                                           1414
      jerr=8001                                                            1414
      return                                                               1414
12071 continue                                                             1415
      if(q0 .lt. 1.0-pmin)goto 12091                                       1415
      jerr=9001                                                            1415
      return                                                               1415
12091 continue                                                             1416
      ixx=0                                                                1416
      al=0.0                                                               1416
      bz=log(q0/(1.0-q0))                                                  1417
      if(nonzero(no,g) .ne. 0)goto 12111                                   1417
      vi=q0*(1.0-q0)                                                       1417
      b(0)=bz                                                              1417
      v=vi*w                                                               1418
      r=w*(y-q0)                                                           1418
      q=q0                                                                 1418
      xmz=vi                                                               1418
      dev1=-(bz*q0+log(1.0-q0))                                            1419
      goto 12121                                                           1420
12111 continue                                                             1420
      b(0)=azero(no,y,g,w,jerr)                                            1420
      if(jerr.ne.0) return                                                 1421
      q=1.0/(1.0+exp(-b(0)-g))                                             1421
      v=w*q*(1.0-q)                                                        1421
      r=w*(y-q)                                                            1421
      xmz=sum(v)                                                           1422
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1423
12121 continue                                                             1424
12101 continue                                                             1424
      if(kopt .le. 0)goto 12141                                            1425
      if(isd .le. 0)goto 12161                                             1425
      xv=0.25                                                              1425
      goto 12171                                                           1426
12161 continue                                                             1426
12180 do 12181 j=1,ni                                                      1426
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1426
12181 continue                                                             1426
12182 continue                                                             1426
12171 continue                                                             1427
12151 continue                                                             1427
12141 continue                                                             1428
      dev0=dev1                                                            1429
12190 do 12191 i=1,no                                                      1429
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1430
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1431
12191 continue                                                             1432
12192 continue                                                             1432
      if(flmin .ge. 1.0)goto 12211                                         1432
      eqs=max(eps,flmin)                                                   1432
      alf=eqs**(1.0/(nlam-1))                                              1432
12211 continue                                                             1433
      m=0                                                                  1433
      mm=0                                                                 1433
      nlp=0                                                                1433
      nin=nlp                                                              1433
      mnl=min(mnlam,nlam)                                                  1433
      bs=0.0                                                               1433
      b(1:ni)=0.0                                                          1434
      shr=shri*dev0                                                        1435
12220 do 12221 j=1,ni                                                      1435
      if(ju(j).eq.0)goto 12221                                             1435
      ga(j)=abs(dot_product(r,x(:,j)))                                     1435
12221 continue                                                             1436
12222 continue                                                             1436
12230 do 12231 ilm=1,nlam                                                  1436
      al0=al                                                               1437
      if(flmin .lt. 1.0)goto 12251                                         1437
      al=ulam(ilm)                                                         1437
      goto 12241                                                           1438
12251 if(ilm .le. 2)goto 12261                                             1438
      al=al*alf                                                            1438
      goto 12241                                                           1439
12261 if(ilm .ne. 1)goto 12271                                             1439
      al=big                                                               1439
      goto 12281                                                           1440
12271 continue                                                             1440
      al0=0.0                                                              1441
12290 do 12291 j=1,ni                                                      1441
      if(ju(j).eq.0)goto 12291                                             1441
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1441
12291 continue                                                             1442
12292 continue                                                             1442
      al0=al0/max(bta,1.0e-3)                                              1442
      al=alf*al0                                                           1443
12281 continue                                                             1444
12241 continue                                                             1444
      al2=al*omb                                                           1444
      al1=al*bta                                                           1444
      tlam=bta*(2.0*al-al0)                                                1445
12300 do 12301 k=1,ni                                                      1445
      if(ixx(k).eq.1)goto 12301                                            1445
      if(ju(k).eq.0)goto 12301                                             1446
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1447
12301 continue                                                             1448
12302 continue                                                             1448
10680 continue                                                             1449
12310 continue                                                             1449
12311 continue                                                             1449
      bs(0)=b(0)                                                           1449
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1450
      if(kopt .ne. 0)goto 12331                                            1451
12340 do 12341 j=1,ni                                                      1451
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                       1451
12341 continue                                                             1452
12342 continue                                                             1452
12331 continue                                                             1453
12350 continue                                                             1453
12351 continue                                                             1453
      nlp=nlp+1                                                            1453
      dlx=0.0                                                              1454
12360 do 12361 k=1,ni                                                      1454
      if(ixx(k).eq.0)goto 12361                                            1455
      bk=b(k)                                                              1455
      gk=dot_product(r,x(:,k))                                             1456
      u=gk+xv(k)*b(k)                                                      1456
      au=abs(u)-vp(k)*al1                                                  1457
      if(au .gt. 0.0)goto 12381                                            1457
      b(k)=0.0                                                             1457
      goto 12391                                                           1458
12381 continue                                                             1458
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1458
12391 continue                                                             1459
12371 continue                                                             1459
      d=b(k)-bk                                                            1459
      if(abs(d).le.0.0)goto 12361                                          1459
      dlx=max(dlx,xv(k)*d**2)                                              1460
      r=r-d*v*x(:,k)                                                       1461
      if(mm(k) .ne. 0)goto 12411                                           1461
      nin=nin+1                                                            1461
      if(nin.gt.nx)goto 12362                                              1462
      mm(k)=nin                                                            1462
      m(nin)=k                                                             1463
12411 continue                                                             1464
12361 continue                                                             1465
12362 continue                                                             1465
      if(nin.gt.nx)goto 12352                                              1466
      d=sum(r)/xmz                                                         1467
      if(d .eq. 0.0)goto 12431                                             1467
      b(0)=b(0)+d                                                          1467
      dlx=max(dlx,xmz*d**2)                                                1467
      r=r-d*v                                                              1467
12431 continue                                                             1468
      if(dlx.lt.shr)goto 12352                                             1468
      if(nlp .le. maxit)goto 12451                                         1468
      jerr=-ilm                                                            1468
      return                                                               1468
12451 continue                                                             1469
12460 continue                                                             1469
12461 continue                                                             1469
      nlp=nlp+1                                                            1469
      dlx=0.0                                                              1470
12470 do 12471 l=1,nin                                                     1470
      k=m(l)                                                               1470
      bk=b(k)                                                              1471
      gk=dot_product(r,x(:,k))                                             1472
      u=gk+xv(k)*b(k)                                                      1472
      au=abs(u)-vp(k)*al1                                                  1473
      if(au .gt. 0.0)goto 12491                                            1473
      b(k)=0.0                                                             1473
      goto 12501                                                           1474
12491 continue                                                             1474
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1474
12501 continue                                                             1475
12481 continue                                                             1475
      d=b(k)-bk                                                            1475
      if(abs(d).le.0.0)goto 12471                                          1475
      dlx=max(dlx,xv(k)*d**2)                                              1476
      r=r-d*v*x(:,k)                                                       1477
12471 continue                                                             1478
12472 continue                                                             1478
      d=sum(r)/xmz                                                         1479
      if(d .eq. 0.0)goto 12521                                             1479
      b(0)=b(0)+d                                                          1479
      dlx=max(dlx,xmz*d**2)                                                1479
      r=r-d*v                                                              1479
12521 continue                                                             1480
      if(dlx.lt.shr)goto 12462                                             1480
      if(nlp .le. maxit)goto 12541                                         1480
      jerr=-ilm                                                            1480
      return                                                               1480
12541 continue                                                             1481
      goto 12461                                                           1482
12462 continue                                                             1482
      goto 12351                                                           1483
12352 continue                                                             1483
      if(nin.gt.nx)goto 12312                                              1484
12550 do 12551 i=1,no                                                      1484
      fi=b(0)+g(i)                                                         1485
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1486
      if(fi .ge. fmin)goto 12571                                           1486
      q(i)=0.0                                                             1486
      goto 12561                                                           1486
12571 if(fi .le. fmax)goto 12581                                           1486
      q(i)=1.0                                                             1486
      goto 12591                                                           1487
12581 continue                                                             1487
      q(i)=1.0/(1.0+exp(-fi))                                              1487
12591 continue                                                             1488
12561 continue                                                             1488
12551 continue                                                             1489
12552 continue                                                             1489
      v=w*q*(1.0-q)                                                        1489
      xmz=sum(v)                                                           1489
      if(xmz.le.vmin)goto 12312                                            1489
      r=w*(y-q)                                                            1490
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 12611                           1490
      ix=0                                                                 1491
12620 do 12621 j=1,nin                                                     1491
      k=m(j)                                                               1492
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 12621                           1492
      ix=1                                                                 1492
      goto 12622                                                           1493
12621 continue                                                             1494
12622 continue                                                             1494
      if(ix .ne. 0)goto 12641                                              1495
12650 do 12651 k=1,ni                                                      1495
      if(ixx(k).eq.1)goto 12651                                            1495
      if(ju(k).eq.0)goto 12651                                             1496
      ga(k)=abs(dot_product(r,x(:,k)))                                     1497
      if(ga(k) .le. al1*vp(k))goto 12671                                   1497
      ixx(k)=1                                                             1497
      ix=1                                                                 1497
12671 continue                                                             1498
12651 continue                                                             1499
12652 continue                                                             1499
      if(ix.eq.1) go to 10680                                              1500
      goto 12312                                                           1501
12641 continue                                                             1502
12611 continue                                                             1503
      goto 12311                                                           1504
12312 continue                                                             1504
      if(nin .le. nx)goto 12691                                            1504
      jerr=-10000-ilm                                                      1504
      goto 12232                                                           1504
12691 continue                                                             1505
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1505
      kin(ilm)=nin                                                         1506
      a0(ilm)=b(0)                                                         1506
      alm(ilm)=al                                                          1506
      lmu=ilm                                                              1507
      devi=dev2(no,w,y,q,pmin)                                             1508
      dev(ilm)=(dev1-devi)/dev0                                            1508
      if(xmz.le.vmin)goto 12232                                            1509
      if(ilm.lt.mnl)goto 12231                                             1509
      if(flmin.ge.1.0)goto 12231                                           1510
      me=0                                                                 1510
12700 do 12701 j=1,nin                                                     1510
      if(a(j,ilm).ne.0.0) me=me+1                                          1510
12701 continue                                                             1510
12702 continue                                                             1510
      if(me.gt.ne)goto 12232                                               1511
      if(dev(ilm).gt.devmax)goto 12232                                     1511
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12232                             1512
12231 continue                                                             1513
12232 continue                                                             1513
      g=log(q/(1.0-q))                                                     1514
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                                  1515
      return                                                               1516
      end                                                                  1517
      function dev2(n,w,y,p,pmin)                                          1518
      real w(n),y(n),p(n)                                                  1519
      pmax=1.0-pmin                                                        1519
      s=0.0                                                                1520
12710 do 12711 i=1,n                                                       1520
      pi=min(max(pmin,p(i)),pmax)                                          1521
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1522
12711 continue                                                             1523
12712 continue                                                             1523
      dev2=s                                                               1524
      return                                                               1525
      end                                                                  1526
      function azero(n,y,g,q,jerr)                                         1527
      parameter(eps=1.0e-7)                                                1528
      real y(n),g(n),q(n)                                                  1529
      real, dimension (:), allocatable :: e,p,w                                 
      allocate(e(1:n),stat=jerr)                                           1533
      allocate(p(1:n),stat=ierr)                                           1533
      jerr=jerr+ierr                                                       1534
      allocate(w(1:n),stat=ierr)                                           1534
      jerr=jerr+ierr                                                       1535
      if(jerr.ne.0) return                                                 1536
      az=0.0                                                               1536
      e=exp(-g)                                                            1536
      qy=dot_product(q,y)                                                  1536
      p=1.0/(1.0+e)                                                        1537
12720 continue                                                             1537
12721 continue                                                             1537
      w=q*p*(1.0-p)                                                        1538
      d=(qy-dot_product(q,p))/sum(w)                                       1538
      az=az+d                                                              1538
      if(abs(d).lt.eps)goto 12722                                          1539
      ea0=exp(-az)                                                         1539
      p=1.0/(1.0+ea0*e)                                                    1540
      goto 12721                                                           1541
12722 continue                                                             1541
      azero=az                                                             1542
      deallocate(e,p,w)                                                    1543
      return                                                               1544
      end                                                                  1545
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ul   1547 
     *am,shri,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1549 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam)              1550
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1551
      integer ju(ni),m(nx),kin(nlam)                                       1552
      real, dimension (:,:), allocatable :: q                                   
      real, dimension (:), allocatable :: sxp,sxpl                              
      real, dimension (:), allocatable :: di,v,r,ga                             
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is,ixx                          
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no),stat=ierr)                                          1563
      jerr=jerr+ierr                                                       1564
      allocate(v(1:no),stat=ierr)                                          1564
      jerr=jerr+ierr                                                       1565
      allocate(mm(1:ni),stat=ierr)                                         1565
      jerr=jerr+ierr                                                       1566
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1566
      jerr=jerr+ierr                                                       1567
      allocate(sxp(1:no),stat=ierr)                                        1567
      jerr=jerr+ierr                                                       1568
      allocate(sxpl(1:no),stat=ierr)                                       1568
      jerr=jerr+ierr                                                       1569
      allocate(di(1:no),stat=ierr)                                         1569
      jerr=jerr+ierr                                                       1570
      allocate(ga(1:ni),stat=ierr)                                         1570
      jerr=jerr+ierr                                                       1571
      allocate(ixx(1:ni),stat=ierr)                                        1571
      jerr=jerr+ierr                                                       1572
      if(jerr.ne.0) return                                                 1573
      pmax=1.0-pmin                                                        1573
      emin=pmin/pmax                                                       1573
      emax=1.0/emin                                                        1574
      pfm=(1.0+pmin)*pmin                                                  1574
      pfx=(1.0-pmin)*pmax                                                  1574
      vmin=pfm*pmax                                                        1575
      bta=parm                                                             1575
      omb=1.0-bta                                                          1575
      dev1=0.0                                                             1575
      dev0=0.0                                                             1576
12730 do 12731 ic=1,nc                                                     1576
      q0=dot_product(w,y(:,ic))                                            1577
      if(q0 .gt. pmin)goto 12751                                           1577
      jerr =8000+ic                                                        1577
      return                                                               1577
12751 continue                                                             1578
      if(q0 .lt. 1.0-pmin)goto 12771                                       1578
      jerr =9000+ic                                                        1578
      return                                                               1578
12771 continue                                                             1579
      b(0,ic)=log(q0)                                                      1579
      dev1=dev1-q0*b(0,ic)                                                 1579
      b(1:ni,ic)=0.0                                                       1580
12731 continue                                                             1581
12732 continue                                                             1581
      ixx=0                                                                1581
      al=0.0                                                               1582
      if(nonzero(no*nc,g) .ne. 0)goto 12791                                1583
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1583
      sxp=0.0                                                              1584
12800 do 12801 ic=1,nc                                                     1584
      q(:,ic)=exp(b(0,ic))                                                 1584
      sxp=sxp+q(:,ic)                                                      1584
12801 continue                                                             1585
12802 continue                                                             1585
      goto 12811                                                           1586
12791 continue                                                             1586
12820 do 12821 i=1,no                                                      1586
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1586
12821 continue                                                             1586
12822 continue                                                             1586
      sxp=0.0                                                              1587
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1587
      if(jerr.ne.0) return                                                 1588
      dev1=0.0                                                             1589
12830 do 12831 ic=1,nc                                                     1589
      q(:,ic)=b(0,ic)+g(:,ic)                                              1590
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             1591
      q(:,ic)=exp(q(:,ic))                                                 1591
      sxp=sxp+q(:,ic)                                                      1592
12831 continue                                                             1593
12832 continue                                                             1593
      sxpl=w*log(sxp)                                                      1593
12840 do 12841 ic=1,nc                                                     1593
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  1593
12841 continue                                                             1594
12842 continue                                                             1594
12811 continue                                                             1595
12781 continue                                                             1595
12850 do 12851 ic=1,nc                                                     1595
12860 do 12861 i=1,no                                                      1595
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1595
12861 continue                                                             1595
12862 continue                                                             1595
12851 continue                                                             1596
12852 continue                                                             1596
      dev0=dev0+dev1                                                       1597
      if(kopt .le. 0)goto 12881                                            1598
      if(isd .le. 0)goto 12901                                             1598
      xv=0.25                                                              1598
      goto 12911                                                           1599
12901 continue                                                             1599
12920 do 12921 j=1,ni                                                      1599
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1599
12921 continue                                                             1599
12922 continue                                                             1599
12911 continue                                                             1600
12891 continue                                                             1600
12881 continue                                                             1601
      if(flmin .ge. 1.0)goto 12941                                         1601
      eqs=max(eps,flmin)                                                   1601
      alf=eqs**(1.0/(nlam-1))                                              1601
12941 continue                                                             1602
      m=0                                                                  1602
      mm=0                                                                 1602
      nin=0                                                                1602
      nlp=0                                                                1602
      mnl=min(mnlam,nlam)                                                  1602
      bs=0.0                                                               1602
      shr=shri*dev0                                                        1603
      ga=0.0                                                               1604
12950 do 12951 ic=1,nc                                                     1604
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1605
12960 do 12961 j=1,ni                                                      1605
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))           1605
12961 continue                                                             1606
12962 continue                                                             1606
12951 continue                                                             1607
12952 continue                                                             1607
12970 do 12971 ilm=1,nlam                                                  1607
      al0=al                                                               1608
      if(flmin .lt. 1.0)goto 12991                                         1608
      al=ulam(ilm)                                                         1608
      goto 12981                                                           1609
12991 if(ilm .le. 2)goto 13001                                             1609
      al=al*alf                                                            1609
      goto 12981                                                           1610
13001 if(ilm .ne. 1)goto 13011                                             1610
      al=big                                                               1610
      goto 13021                                                           1611
13011 continue                                                             1611
      al0=0.0                                                              1612
13030 do 13031 j=1,ni                                                      1612
      if(ju(j).eq.0)goto 13031                                             1612
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1612
13031 continue                                                             1613
13032 continue                                                             1613
      al0=al0/max(bta,1.0e-3)                                              1613
      al=alf*al0                                                           1614
13021 continue                                                             1615
12981 continue                                                             1615
      al2=al*omb                                                           1615
      al1=al*bta                                                           1615
      tlam=bta*(2.0*al-al0)                                                1616
13040 do 13041 k=1,ni                                                      1616
      if(ixx(k).eq.1)goto 13041                                            1616
      if(ju(k).eq.0)goto 13041                                             1617
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1618
13041 continue                                                             1619
13042 continue                                                             1619
10680 continue                                                             1620
13050 continue                                                             1620
13051 continue                                                             1620
      ix=0                                                                 1620
      jx=ix                                                                1620
      ig=0                                                                 1621
13060 do 13061 ic=1,nc                                                     1621
      bs(0,ic)=b(0,ic)                                                     1622
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1623
      xmz=0.0                                                              1624
13070 do 13071 i=1,no                                                      1624
      pic=q(i,ic)/sxp(i)                                                   1625
      if(pic .ge. pfm)goto 13091                                           1625
      pic=0.0                                                              1625
      v(i)=0.0                                                             1625
      goto 13081                                                           1626
13091 if(pic .le. pfx)goto 13101                                           1626
      pic=1.0                                                              1626
      v(i)=0.0                                                             1626
      goto 13111                                                           1627
13101 continue                                                             1627
      v(i)=w(i)*pic*(1.0-pic)                                              1627
      xmz=xmz+v(i)                                                         1627
13111 continue                                                             1628
13081 continue                                                             1628
      r(i)=w(i)*(y(i,ic)-pic)                                              1629
13071 continue                                                             1630
13072 continue                                                             1630
      if(xmz.le.vmin)goto 13061                                            1630
      ig=1                                                                 1631
      if(kopt .ne. 0)goto 13131                                            1632
13140 do 13141 j=1,ni                                                      1632
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                    1632
13141 continue                                                             1633
13142 continue                                                             1633
13131 continue                                                             1634
13150 continue                                                             1634
13151 continue                                                             1634
      nlp=nlp+1                                                            1634
      dlx=0.0                                                              1635
13160 do 13161 k=1,ni                                                      1635
      if(ixx(k).eq.0)goto 13161                                            1636
      bk=b(k,ic)                                                           1636
      gk=dot_product(r,x(:,k))                                             1637
      u=gk+xv(k,ic)*b(k,ic)                                                1637
      au=abs(u)-vp(k)*al1                                                  1638
      if(au .gt. 0.0)goto 13181                                            1638
      b(k,ic)=0.0                                                          1638
      goto 13191                                                           1639
13181 continue                                                             1639
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1639
13191 continue                                                             1640
13171 continue                                                             1640
      d=b(k,ic)-bk                                                         1640
      if(abs(d).le.0.0)goto 13161                                          1641
      dlx=max(dlx,xv(k,ic)*d**2)                                           1641
      r=r-d*v*x(:,k)                                                       1642
      if(mm(k) .ne. 0)goto 13211                                           1642
      nin=nin+1                                                            1643
      if(nin .le. nx)goto 13231                                            1643
      jx=1                                                                 1643
      goto 13162                                                           1643
13231 continue                                                             1644
      mm(k)=nin                                                            1644
      m(nin)=k                                                             1645
13211 continue                                                             1646
13161 continue                                                             1647
13162 continue                                                             1647
      if(jx.gt.0)goto 13152                                                1648
      d=sum(r)/xmz                                                         1649
      if(d .eq. 0.0)goto 13251                                             1649
      b(0,ic)=b(0,ic)+d                                                    1649
      dlx=max(dlx,xmz*d**2)                                                1649
      r=r-d*v                                                              1649
13251 continue                                                             1650
      if(dlx.lt.shr)goto 13152                                             1651
      if(nlp .le. maxit)goto 13271                                         1651
      jerr=-ilm                                                            1651
      return                                                               1651
13271 continue                                                             1652
13280 continue                                                             1652
13281 continue                                                             1652
      nlp=nlp+1                                                            1652
      dlx=0.0                                                              1653
13290 do 13291 l=1,nin                                                     1653
      k=m(l)                                                               1653
      bk=b(k,ic)                                                           1654
      gk=dot_product(r,x(:,k))                                             1655
      u=gk+xv(k,ic)*b(k,ic)                                                1655
      au=abs(u)-vp(k)*al1                                                  1656
      if(au .gt. 0.0)goto 13311                                            1656
      b(k,ic)=0.0                                                          1656
      goto 13321                                                           1657
13311 continue                                                             1657
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1657
13321 continue                                                             1658
13301 continue                                                             1658
      d=b(k,ic)-bk                                                         1658
      if(abs(d).le.0.0)goto 13291                                          1659
      dlx=max(dlx,xv(k,ic)*d**2)                                           1659
      r=r-d*v*x(:,k)                                                       1660
13291 continue                                                             1661
13292 continue                                                             1661
      d=sum(r)/xmz                                                         1662
      if(d .eq. 0.0)goto 13341                                             1662
      b(0,ic)=b(0,ic)+d                                                    1663
      dlx=max(dlx,xmz*d**2)                                                1663
      r=r-d*v                                                              1664
13341 continue                                                             1665
      if(dlx.lt.shr)goto 13282                                             1665
      if(nlp .le. maxit)goto 13361                                         1665
      jerr=-ilm                                                            1665
      return                                                               1665
13361 continue                                                             1666
      goto 13281                                                           1667
13282 continue                                                             1667
      goto 13151                                                           1668
13152 continue                                                             1668
      if(jx.gt.0)goto 13062                                                1669
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                            1670
      if(ix .ne. 0)goto 13381                                              1671
13390 do 13391 j=1,nin                                                     1671
      k=m(j)                                                               1672
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 13411                1672
      ix=1                                                                 1672
      goto 13392                                                           1672
13411 continue                                                             1673
13391 continue                                                             1674
13392 continue                                                             1674
13381 continue                                                             1675
13420 do 13421 i=1,no                                                      1675
      fi=b(0,ic)+g(i,ic)                                                   1677
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1678
      fi=min(max(exmn,fi),exmx)                                            1678
      sxp(i)=sxp(i)-q(i,ic)                                                1679
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    1680
      sxp(i)=sxp(i)+q(i,ic)                                                1681
13421 continue                                                             1682
13422 continue                                                             1682
13061 continue                                                             1683
13062 continue                                                             1683
      s=-sum(b(0,:))/nc                                                    1683
      b(0,:)=b(0,:)+s                                                      1683
      di=s                                                                 1684
13430 do 13431 j=1,nin                                                     1684
      l=m(j)                                                               1685
      if(vp(l) .gt. 0.0)goto 13451                                         1685
      s=sum(b(l,:))/nc                                                     1685
      goto 13461                                                           1686
13451 continue                                                             1686
      s=elc(parm,nc,b(l,:),is)                                             1686
13461 continue                                                             1687
13441 continue                                                             1687
      b(l,:)=b(l,:)-s                                                      1687
      di=di-s*x(:,l)                                                       1688
13431 continue                                                             1689
13432 continue                                                             1689
      di=exp(di)                                                           1689
      sxp=sxp*di                                                           1689
13470 do 13471 ic=1,nc                                                     1689
      q(:,ic)=q(:,ic)*di                                                   1689
13471 continue                                                             1690
13472 continue                                                             1690
      if(jx.gt.0)goto 13052                                                1690
      if(ig.eq.0)goto 13052                                                1691
      if(ix .ne. 0)goto 13491                                              1692
13500 do 13501 k=1,ni                                                      1692
      if(ixx(k).eq.1)goto 13501                                            1692
      if(ju(k).eq.0)goto 13501                                             1692
      ga(k)=0.0                                                            1692
13501 continue                                                             1693
13502 continue                                                             1693
13510 do 13511 ic=1,nc                                                     1693
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1694
13520 do 13521 k=1,ni                                                      1694
      if(ixx(k).eq.1)goto 13521                                            1694
      if(ju(k).eq.0)goto 13521                                             1695
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                          1696
13521 continue                                                             1697
13522 continue                                                             1697
13511 continue                                                             1698
13512 continue                                                             1698
13530 do 13531 k=1,ni                                                      1698
      if(ixx(k).eq.1)goto 13531                                            1698
      if(ju(k).eq.0)goto 13531                                             1699
      if(ga(k) .le. al1*vp(k))goto 13551                                   1699
      ixx(k)=1                                                             1699
      ix=1                                                                 1699
13551 continue                                                             1700
13531 continue                                                             1701
13532 continue                                                             1701
      if(ix.eq.1) go to 10680                                              1702
      goto 13052                                                           1703
13491 continue                                                             1704
      goto 13051                                                           1705
13052 continue                                                             1705
      if(jx .le. 0)goto 13571                                              1705
      jerr=-10000-ilm                                                      1705
      goto 12972                                                           1705
13571 continue                                                             1705
      devi=0.0                                                             1706
13580 do 13581 ic=1,nc                                                     1707
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1707
      a0(ic,ilm)=b(0,ic)                                                   1708
13590 do 13591 i=1,no                                                      1708
      if(y(i,ic).le.0.0)goto 13591                                         1709
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1710
13591 continue                                                             1711
13592 continue                                                             1711
13581 continue                                                             1712
13582 continue                                                             1712
      kin(ilm)=nin                                                         1712
      alm(ilm)=al                                                          1712
      lmu=ilm                                                              1713
      dev(ilm)=(dev1-devi)/dev0                                            1713
      if(ig.eq.0)goto 12972                                                1714
      if(ilm.lt.mnl)goto 12971                                             1714
      if(flmin.ge.1.0)goto 12971                                           1715
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 12972             1716
      if(dev(ilm).gt.devmax)goto 12972                                     1716
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12972                             1717
12971 continue                                                             1718
12972 continue                                                             1718
      g=log(q)                                                             1718
13600 do 13601 i=1,no                                                      1718
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1718
13601 continue                                                             1719
13602 continue                                                             1719
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                           1720
      return                                                               1721
      end                                                                  1722
      subroutine kazero(kk,n,y,g,q,az,jerr)                                1723
      parameter(eps=1.0e-7)                                                1724
      real y(n,kk),g(n,kk),q(n),az(kk)                                     1725
      real, dimension (:), allocatable :: s                                     
      real, dimension (:,:), allocatable :: e                                   
      allocate(e(1:n,1:kk),stat=jerr)                                           
      allocate(s(1:n),stat=ierr)                                           1730
      jerr=jerr+ierr                                                       1731
      if(jerr.ne.0) return                                                 1732
      az=0.0                                                               1732
      e=exp(g)                                                             1732
13610 do 13611 i=1,n                                                       1732
      s(i)=sum(e(i,:))                                                     1732
13611 continue                                                             1733
13612 continue                                                             1733
13620 continue                                                             1733
13621 continue                                                             1733
      dm=0.0                                                               1734
13630 do 13631 k=1,kk                                                      1734
      t=0.0                                                                1734
      u=t                                                                  1735
13640 do 13641 i=1,n                                                       1735
      pik=e(i,k)/s(i)                                                      1736
      t=t+q(i)*(y(i,k)-pik)                                                1736
      u=u+q(i)*pik*(1.0-pik)                                               1737
13641 continue                                                             1738
13642 continue                                                             1738
      d=t/u                                                                1738
      az(k)=az(k)+d                                                        1738
      ed=exp(d)                                                            1738
      dm=max(dm,abs(d))                                                    1739
13650 do 13651 i=1,n                                                       1739
      z=e(i,k)                                                             1739
      e(i,k)=z*ed                                                          1739
      s(i)=s(i)-z+e(i,k)                                                   1739
13651 continue                                                             1740
13652 continue                                                             1740
13631 continue                                                             1741
13632 continue                                                             1741
      if(dm.lt.eps)goto 13622                                              1741
      goto 13621                                                           1742
13622 continue                                                             1742
      az=az-sum(az)/kk                                                     1743
      deallocate(e,s)                                                      1744
      return                                                               1745
      end                                                                  1746
      function elc(parm,n,a,m)                                             1747
      real a(n)                                                            1747
      integer m(n)                                                         1748
      fn=n                                                                 1748
      am=sum(a)/fn                                                         1749
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 13671                       1749
      elc=am                                                               1749
      return                                                               1749
13671 continue                                                             1750
13680 do 13681 i=1,n                                                       1750
      m(i)=i                                                               1750
13681 continue                                                             1750
13682 continue                                                             1750
      call psort7(a,m,1,n)                                                 1751
      if(a(m(1)) .ne. a(m(n)))goto 13701                                   1751
      elc=a(1)                                                             1751
      return                                                               1751
13701 continue                                                             1752
      if(mod(n,2) .ne. 1)goto 13721                                        1752
      ad=a(m(n/2+1))                                                       1752
      goto 13731                                                           1753
13721 continue                                                             1753
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       1753
13731 continue                                                             1754
13711 continue                                                             1754
      if(parm .ne. 1.0)goto 13751                                          1754
      elc=ad                                                               1754
      return                                                               1754
13751 continue                                                             1755
      b1=min(am,ad)                                                        1755
      b2=max(am,ad)                                                        1755
      k2=1                                                                 1756
13760 continue                                                             1756
13761 if(a(m(k2)).gt.b1)goto 13762                                         1756
      k2=k2+1                                                              1756
      goto 13761                                                           1756
13762 continue                                                             1756
      k1=k2-1                                                              1757
13770 continue                                                             1757
13771 if(a(m(k2)).ge.b2)goto 13772                                         1757
      k2=k2+1                                                              1757
      goto 13771                                                           1758
13772 continue                                                             1758
      r=parm/((1.0-parm)*fn)                                               1758
      is=0                                                                 1758
      sm=n-2*(k1-1)                                                        1759
13780 do 13781 k=k1,k2-1                                                   1759
      sm=sm-2.0                                                            1759
      s=r*sm+am                                                            1760
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 13801                   1760
      is=k                                                                 1760
      goto 13782                                                           1760
13801 continue                                                             1761
13781 continue                                                             1762
13782 continue                                                             1762
      if(is .eq. 0)goto 13821                                              1762
      elc=s                                                                1762
      return                                                               1762
13821 continue                                                             1762
      r2=2.0*r                                                             1762
      s1=a(m(k1))                                                          1762
      am2=2.0*am                                                           1763
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    1763
      elc=s1                                                               1764
13830 do 13831 k=k1+1,k2                                                   1764
      s=a(m(k))                                                            1764
      if(s.eq.s1)goto 13831                                                1765
      c=r2*sum(abs(a-s))+s*(s-am2)                                         1766
      if(c .ge. cri)goto 13851                                             1766
      cri=c                                                                1766
      elc=s                                                                1766
13851 continue                                                             1766
      s1=s                                                                 1767
13831 continue                                                             1768
13832 continue                                                             1768
      return                                                               1769
      end                                                                  1770
      function nintot(ni,nx,nc,a,m,nin,is)                                 1771
      real a(nx,nc)                                                        1771
      integer m(nx),is(ni)                                                 1772
      is=0                                                                 1772
      nintot=0                                                             1773
13860 do 13861 ic=1,nc                                                     1773
13870 do 13871 j=1,nin                                                     1773
      k=m(j)                                                               1773
      if(is(k).ne.0)goto 13871                                             1774
      if(a(j,ic).eq.0.0)goto 13871                                         1774
      is(k)=k                                                              1774
      nintot=nintot+1                                                      1775
13871 continue                                                             1775
13872 continue                                                             1775
13861 continue                                                             1776
13862 continue                                                             1776
      return                                                               1777
      end                                                                  1778
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             1779
      real ca(nx,nc),a(ni,nc)                                              1779
      integer ia(nx)                                                       1780
      a=0.0                                                                1781
13880 do 13881 ic=1,nc                                                     1781
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            1781
13881 continue                                                             1782
13882 continue                                                             1782
      return                                                               1783
      end                                                                  1784
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      1785
      real a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                             1785
      integer ia(nx)                                                       1786
13890 do 13891 i=1,nt                                                      1786
13900 do 13901 ic=1,nc                                                     1786
      ans(ic,i)=a0(ic)                                                     1788
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   1789 
     *:nin)))
13901 continue                                                             1789
13902 continue                                                             1789
13891 continue                                                             1790
13892 continue                                                             1790
      return                                                               1791
      end                                                                  1792
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,ne,nx,nlam,fl   1794 
     *min,  ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      real x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)                 1795
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1796
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1797
      real, dimension (:), allocatable :: xm,xs,ww,vq,xv                        
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 13921                                    1801
      jerr=10000                                                           1801
      return                                                               1801
13921 continue                                                             1802
      allocate(ww(1:no),stat=jerr)                                         1803
      allocate(ju(1:ni),stat=ierr)                                         1803
      jerr=jerr+ierr                                                       1804
      allocate(vq(1:ni),stat=ierr)                                         1804
      jerr=jerr+ierr                                                       1805
      allocate(xm(1:ni),stat=ierr)                                         1805
      jerr=jerr+ierr                                                       1806
      allocate(xs(1:ni),stat=ierr)                                         1806
      jerr=jerr+ierr                                                       1807
      if(kopt .ne. 2)goto 13941                                            1807
      allocate(xv(1:ni),stat=ierr)                                         1807
      jerr=jerr+ierr                                                       1807
13941 continue                                                             1808
      if(jerr.ne.0) return                                                 1809
      call spchkvars(no,ni,x,ix,ju)                                        1810
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1811
      if(maxval(ju) .gt. 0)goto 13961                                      1811
      jerr=7777                                                            1811
      return                                                               1811
13961 continue                                                             1812
      vq=max(0.0,vp)                                                       1812
      vq=vq*ni/sum(vq)                                                     1813
13970 do 13971 i=1,no                                                      1813
      ww(i)=sum(y(i,:))                                                    1813
      y(i,:)=y(i,:)/ww(i)                                                  1813
13971 continue                                                             1813
13972 continue                                                             1813
      sw=sum(ww)                                                           1813
      ww=ww/sw                                                             1814
      if(nc .ne. 1)goto 13991                                              1814
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     1815
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,ne,nx,n   1817 
     *lam,flmin,  ulam,thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,de
     *v,alm,nlp,jerr)
      goto 13981                                                           1818
13991 if(kopt .ne. 2)goto 14001                                            1818
      call multsplstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs,xv)              1819
      call multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,   1821 
     *flmin,ulam, thr,isd,maxit,xv,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,n
     *lp,jerr)
      goto 14011                                                           1822
14001 continue                                                             1822
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     1823
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmi   1825 
     *n,ulam,  thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nl
     *p,jerr)
14011 continue                                                             1826
13981 continue                                                             1826
      if(jerr.gt.0) return                                                 1826
      dev0=2.0*sw*dev0                                                     1827
14020 do 14021 k=1,lmu                                                     1827
      nk=nin(k)                                                            1828
14030 do 14031 ic=1,nc                                                     1828
      if(isd .le. 0)goto 14051                                             1828
14060 do 14061 l=1,nk                                                      1828
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1828
14061 continue                                                             1828
14062 continue                                                             1828
14051 continue                                                             1829
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1830
14031 continue                                                             1831
14032 continue                                                             1831
14021 continue                                                             1832
14022 continue                                                             1832
      deallocate(ww,ju,vq,xm,xs)                                           1832
      if(kopt.eq.2) deallocate(xv)                                         1833
      return                                                               1834
      end                                                                  1835
      subroutine multsplstandard2(no,ni,x,ix,jx,w,ju,isd,xm,xs,xv)         1836
      real x(*),w(no),xm(ni),xs(ni),xv(ni)                                 1836
      integer ix(*),jx(*),ju(ni)                                           1837
14070 do 14071 j=1,ni                                                      1837
      if(ju(j).eq.0)goto 14071                                             1837
      jb=ix(j)                                                             1837
      je=ix(j+1)-1                                                         1838
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1839
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1840
      if(isd .le. 0)goto 14091                                             1840
      xs(j)=sqrt(xv(j))                                                    1840
      xv(j)=1.0                                                            1840
14091 continue                                                             1841
14071 continue                                                             1842
14072 continue                                                             1842
      if(isd.eq.0) xs=1.0                                                  1843
      return                                                               1844
      end                                                                  1845
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,xm,xs)                1846
      real x(*),w(no),xm(ni),xs(ni)                                        1846
      integer ix(*),jx(*),ju(ni)                                           1847
14100 do 14101 j=1,ni                                                      1847
      if(ju(j).eq.0)goto 14101                                             1847
      jb=ix(j)                                                             1847
      je=ix(j+1)-1                                                         1848
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1849
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   1850 
     *)**2)
14101 continue                                                             1851
14102 continue                                                             1851
      if(isd.eq.0) xs=1.0                                                  1852
      return                                                               1853
      end                                                                  1854
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam,   1856 
     *  flmin,ulam,shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,alm
     *,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1858 
     *5, devmax=0.999)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        1859
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1860
      real xb(ni),xs(ni)                                                   1860
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1861
      real, dimension (:), allocatable :: xm,b,bs,v,r,sc,xv,q,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(b(0:ni),stat=jerr)                                          1866
      allocate(xm(0:ni),stat=ierr)                                         1866
      jerr=jerr+ierr                                                       1867
      allocate(xv(1:ni),stat=ierr)                                         1867
      jerr=jerr+ierr                                                       1868
      allocate(bs(0:ni),stat=ierr)                                         1868
      jerr=jerr+ierr                                                       1869
      allocate(ga(1:ni),stat=ierr)                                         1869
      jerr=jerr+ierr                                                       1870
      allocate(mm(1:ni),stat=ierr)                                         1870
      jerr=jerr+ierr                                                       1871
      allocate(ixx(1:ni),stat=ierr)                                        1871
      jerr=jerr+ierr                                                       1872
      allocate(q(1:no),stat=ierr)                                          1872
      jerr=jerr+ierr                                                       1873
      allocate(r(1:no),stat=ierr)                                          1873
      jerr=jerr+ierr                                                       1874
      allocate(v(1:no),stat=ierr)                                          1874
      jerr=jerr+ierr                                                       1875
      allocate(sc(1:no),stat=ierr)                                         1875
      jerr=jerr+ierr                                                       1876
      if(jerr.ne.0) return                                                 1877
      fmax=log(1.0/pmin-1.0)                                               1877
      fmin=-fmax                                                           1877
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1878
      bta=parm                                                             1878
      omb=1.0-bta                                                          1879
      q0=dot_product(w,y)                                                  1879
      if(q0 .gt. pmin)goto 14121                                           1879
      jerr=8001                                                            1879
      return                                                               1879
14121 continue                                                             1880
      if(q0 .lt. 1.0-pmin)goto 14141                                       1880
      jerr=9001                                                            1880
      return                                                               1880
14141 continue                                                             1880
      bz=log(q0/(1.0-q0))                                                  1881
      if(nonzero(no,g) .ne. 0)goto 14161                                   1881
      vi=q0*(1.0-q0)                                                       1881
      b(0)=bz                                                              1881
      v=vi*w                                                               1882
      r=w*(y-q0)                                                           1882
      q=q0                                                                 1882
      xm(0)=vi                                                             1882
      dev1=-(bz*q0+log(1.0-q0))                                            1883
      goto 14171                                                           1884
14161 continue                                                             1884
      b(0)=azero(no,y,g,w,jerr)                                            1884
      if(jerr.ne.0) return                                                 1885
      q=1.0/(1.0+exp(-b(0)-g))                                             1885
      v=w*q*(1.0-q)                                                        1885
      r=w*(y-q)                                                            1885
      xm(0)=sum(v)                                                         1886
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1887
14171 continue                                                             1888
14151 continue                                                             1888
      if(kopt .le. 0)goto 14191                                            1889
      if(isd .le. 0)goto 14211                                             1889
      xv=0.25                                                              1889
      goto 14221                                                           1890
14211 continue                                                             1891
14230 do 14231 j=1,ni                                                      1891
      if(ju(j).eq.0)goto 14231                                             1891
      jb=ix(j)                                                             1891
      je=ix(j+1)-1                                                         1892
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          1893
14231 continue                                                             1894
14232 continue                                                             1894
14221 continue                                                             1895
14201 continue                                                             1895
14191 continue                                                             1896
      b(1:ni)=0.0                                                          1896
      dev0=dev1                                                            1897
14240 do 14241 i=1,no                                                      1897
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1898
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1899
14241 continue                                                             1900
14242 continue                                                             1900
      if(flmin .ge. 1.0)goto 14261                                         1900
      eqs=max(eps,flmin)                                                   1900
      alf=eqs**(1.0/(nlam-1))                                              1900
14261 continue                                                             1901
      m=0                                                                  1901
      mm=0                                                                 1901
      nin=0                                                                1901
      o=0.0                                                                1901
      svr=o                                                                1901
      mnl=min(mnlam,nlam)                                                  1901
      bs=0.0                                                               1901
      nlp=0                                                                1901
      nin=nlp                                                              1902
      shr=shri*dev0                                                        1902
      al=0.0                                                               1902
      ixx=0                                                                1903
14270 do 14271 j=1,ni                                                      1903
      if(ju(j).eq.0)goto 14271                                             1904
      jb=ix(j)                                                             1904
      je=ix(j+1)-1                                                         1904
      jn=ix(j+1)-ix(j)                                                     1905
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1906
      gj=dot_product(sc(1:jn),x(jb:je))                                    1907
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      1908
14271 continue                                                             1909
14272 continue                                                             1909
14280 do 14281 ilm=1,nlam                                                  1909
      al0=al                                                               1910
      if(flmin .lt. 1.0)goto 14301                                         1910
      al=ulam(ilm)                                                         1910
      goto 14291                                                           1911
14301 if(ilm .le. 2)goto 14311                                             1911
      al=al*alf                                                            1911
      goto 14291                                                           1912
14311 if(ilm .ne. 1)goto 14321                                             1912
      al=big                                                               1912
      goto 14331                                                           1913
14321 continue                                                             1913
      al0=0.0                                                              1914
14340 do 14341 j=1,ni                                                      1914
      if(ju(j).eq.0)goto 14341                                             1914
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1914
14341 continue                                                             1915
14342 continue                                                             1915
      al0=al0/max(bta,1.0e-3)                                              1915
      al=alf*al0                                                           1916
14331 continue                                                             1917
14291 continue                                                             1917
      al2=al*omb                                                           1917
      al1=al*bta                                                           1917
      tlam=bta*(2.0*al-al0)                                                1918
14350 do 14351 k=1,ni                                                      1918
      if(ixx(k).eq.1)goto 14351                                            1918
      if(ju(k).eq.0)goto 14351                                             1919
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1920
14351 continue                                                             1921
14352 continue                                                             1921
10680 continue                                                             1922
14360 continue                                                             1922
14361 continue                                                             1922
      bs(0)=b(0)                                                           1922
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1923
14370 do 14371 j=1,ni                                                      1923
      if(ixx(j).eq.0)goto 14371                                            1924
      jb=ix(j)                                                             1924
      je=ix(j+1)-1                                                         1924
      jn=ix(j+1)-ix(j)                                                     1925
      sc(1:jn)=v(jx(jb:je))                                                1926
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 1927
      if(kopt .ne. 0)goto 14391                                            1928
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              1929
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                1930
14391 continue                                                             1931
14371 continue                                                             1932
14372 continue                                                             1932
14400 continue                                                             1932
14401 continue                                                             1932
      nlp=nlp+1                                                            1932
      dlx=0.0                                                              1933
14410 do 14411 k=1,ni                                                      1933
      if(ixx(k).eq.0)goto 14411                                            1934
      jb=ix(k)                                                             1934
      je=ix(k+1)-1                                                         1934
      jn=ix(k+1)-ix(k)                                                     1934
      bk=b(k)                                                              1935
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1936
      gk=dot_product(sc(1:jn),x(jb:je))                                    1937
      gk=(gk-svr*xb(k))/xs(k)                                              1938
      u=gk+xv(k)*b(k)                                                      1938
      au=abs(u)-vp(k)*al1                                                  1939
      if(au .gt. 0.0)goto 14431                                            1939
      b(k)=0.0                                                             1939
      goto 14441                                                           1940
14431 continue                                                             1940
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1940
14441 continue                                                             1941
14421 continue                                                             1941
      d=b(k)-bk                                                            1941
      if(abs(d).le.0.0)goto 14411                                          1941
      dlx=max(dlx,xv(k)*d**2)                                              1942
      if(mm(k) .ne. 0)goto 14461                                           1942
      nin=nin+1                                                            1942
      if(nin.gt.nx)goto 14412                                              1943
      mm(k)=nin                                                            1943
      m(nin)=k                                                             1943
      sc(1:jn)=v(jx(jb:je))                                                1944
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 1945
14461 continue                                                             1946
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1947
      o=o+d*(xb(k)/xs(k))                                                  1948
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1949
14411 continue                                                             1950
14412 continue                                                             1950
      if(nin.gt.nx)goto 14402                                              1951
      d=svr/xm(0)                                                          1952
      if(d .eq. 0.0)goto 14481                                             1952
      b(0)=b(0)+d                                                          1952
      dlx=max(dlx,xm(0)*d**2)                                              1952
      r=r-d*v                                                              1952
14481 continue                                                             1953
      svr=svr-d*xm(0)                                                      1953
      if(dlx.lt.shr)goto 14402                                             1954
      if(nlp .le. maxit)goto 14501                                         1954
      jerr=-ilm                                                            1954
      return                                                               1954
14501 continue                                                             1955
14510 continue                                                             1955
14511 continue                                                             1955
      nlp=nlp+1                                                            1955
      dlx=0.0                                                              1956
14520 do 14521 l=1,nin                                                     1956
      k=m(l)                                                               1956
      jb=ix(k)                                                             1956
      je=ix(k+1)-1                                                         1957
      jn=ix(k+1)-ix(k)                                                     1957
      bk=b(k)                                                              1958
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1959
      gk=dot_product(sc(1:jn),x(jb:je))                                    1960
      gk=(gk-svr*xb(k))/xs(k)                                              1961
      u=gk+xv(k)*b(k)                                                      1961
      au=abs(u)-vp(k)*al1                                                  1962
      if(au .gt. 0.0)goto 14541                                            1962
      b(k)=0.0                                                             1962
      goto 14551                                                           1963
14541 continue                                                             1963
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1963
14551 continue                                                             1964
14531 continue                                                             1964
      d=b(k)-bk                                                            1964
      if(abs(d).le.0.0)goto 14521                                          1964
      dlx=max(dlx,xv(k)*d**2)                                              1965
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1966
      o=o+d*(xb(k)/xs(k))                                                  1967
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1968
14521 continue                                                             1969
14522 continue                                                             1969
      d=svr/xm(0)                                                          1970
      if(d .eq. 0.0)goto 14571                                             1970
      b(0)=b(0)+d                                                          1970
      dlx=max(dlx,xm(0)*d**2)                                              1970
      r=r-d*v                                                              1970
14571 continue                                                             1971
      svr=svr-d*xm(0)                                                      1971
      if(dlx.lt.shr)goto 14512                                             1972
      if(nlp .le. maxit)goto 14591                                         1972
      jerr=-ilm                                                            1972
      return                                                               1972
14591 continue                                                             1973
      goto 14511                                                           1974
14512 continue                                                             1974
      goto 14401                                                           1975
14402 continue                                                             1975
      if(nin.gt.nx)goto 14362                                              1976
      sc=b(0)                                                              1976
      b0=0.0                                                               1977
14600 do 14601 j=1,nin                                                     1977
      l=m(j)                                                               1977
      jb=ix(l)                                                             1977
      je=ix(l+1)-1                                                         1978
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      1979
      b0=b0-b(l)*xb(l)/xs(l)                                               1980
14601 continue                                                             1981
14602 continue                                                             1981
      sc=sc+b0                                                             1982
14610 do 14611 i=1,no                                                      1982
      fi=sc(i)+g(i)                                                        1983
      if(fi .ge. fmin)goto 14631                                           1983
      q(i)=0.0                                                             1983
      goto 14621                                                           1983
14631 if(fi .le. fmax)goto 14641                                           1983
      q(i)=1.0                                                             1983
      goto 14651                                                           1984
14641 continue                                                             1984
      q(i)=1.0/(1.0+exp(-fi))                                              1984
14651 continue                                                             1985
14621 continue                                                             1985
14611 continue                                                             1986
14612 continue                                                             1986
      v=w*q*(1.0-q)                                                        1986
      xm(0)=sum(v)                                                         1986
      if(xm(0).lt.vmin)goto 14362                                          1987
      r=w*(y-q)                                                            1987
      svr=sum(r)                                                           1987
      o=0.0                                                                1988
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 14671                         1988
      kx=0                                                                 1989
14680 do 14681 j=1,nin                                                     1989
      k=m(j)                                                               1990
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 14681                           1990
      kx=1                                                                 1990
      goto 14682                                                           1991
14681 continue                                                             1992
14682 continue                                                             1992
      if(kx .ne. 0)goto 14701                                              1993
14710 do 14711 j=1,ni                                                      1993
      if(ixx(j).eq.1)goto 14711                                            1993
      if(ju(j).eq.0)goto 14711                                             1994
      jb=ix(j)                                                             1994
      je=ix(j+1)-1                                                         1994
      jn=ix(j+1)-ix(j)                                                     1995
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1996
      gj=dot_product(sc(1:jn),x(jb:je))                                    1997
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      1998
      if(ga(j) .le. al1*vp(j))goto 14731                                   1998
      ixx(j)=1                                                             1998
      kx=1                                                                 1998
14731 continue                                                             1999
14711 continue                                                             2000
14712 continue                                                             2000
      if(kx.eq.1) go to 10680                                              2001
      goto 14362                                                           2002
14701 continue                                                             2003
14671 continue                                                             2004
      goto 14361                                                           2005
14362 continue                                                             2005
      if(nin .le. nx)goto 14751                                            2005
      jerr=-10000-ilm                                                      2005
      goto 14282                                                           2005
14751 continue                                                             2006
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                2006
      kin(ilm)=nin                                                         2007
      a0(ilm)=b(0)                                                         2007
      alm(ilm)=al                                                          2007
      lmu=ilm                                                              2008
      devi=dev2(no,w,y,q,pmin)                                             2009
      dev(ilm)=(dev1-devi)/dev0                                            2010
      if(ilm.lt.mnl)goto 14281                                             2010
      if(flmin.ge.1.0)goto 14281                                           2011
      me=0                                                                 2011
14760 do 14761 j=1,nin                                                     2011
      if(a(j,ilm).ne.0.0) me=me+1                                          2011
14761 continue                                                             2011
14762 continue                                                             2011
      if(me.gt.ne)goto 14282                                               2012
      if(dev(ilm).gt.devmax)goto 14282                                     2012
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 14282                             2013
      if(xm(0).lt.vmin)goto 14282                                          2014
14281 continue                                                             2015
14282 continue                                                             2015
      g=log(q/(1.0-q))                                                     2016
      deallocate(xm,b,bs,v,r,sc,xv,q,mm,ga,ixx)                            2017
      return                                                               2018
      end                                                                  2019
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam   2021 
     *,flmin,ulam,  shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,al
     *m,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   2023 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni)    2024
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   2025
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2026
      real, dimension (:,:), allocatable :: q                                   
      real, dimension (:), allocatable :: sxp,sxpl                              
      real, dimension (:), allocatable :: sc,xm,v,r,ga                          
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is,iy                           
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(xm(0:ni),stat=ierr)                                         2037
      jerr=jerr+ierr                                                       2038
      allocate(r(1:no),stat=ierr)                                          2038
      jerr=jerr+ierr                                                       2039
      allocate(v(1:no),stat=ierr)                                          2039
      jerr=jerr+ierr                                                       2040
      allocate(mm(1:ni),stat=ierr)                                         2040
      jerr=jerr+ierr                                                       2041
      allocate(ga(1:ni),stat=ierr)                                         2041
      jerr=jerr+ierr                                                       2042
      allocate(iy(1:ni),stat=ierr)                                         2042
      jerr=jerr+ierr                                                       2043
      allocate(is(1:max(nc,ni)),stat=ierr)                                 2043
      jerr=jerr+ierr                                                       2044
      allocate(sxp(1:no),stat=ierr)                                        2044
      jerr=jerr+ierr                                                       2045
      allocate(sxpl(1:no),stat=ierr)                                       2045
      jerr=jerr+ierr                                                       2046
      allocate(sc(1:no),stat=ierr)                                         2046
      jerr=jerr+ierr                                                       2047
      if(jerr.ne.0) return                                                 2048
      pmax=1.0-pmin                                                        2048
      emin=pmin/pmax                                                       2048
      emax=1.0/emin                                                        2049
      pfm=(1.0+pmin)*pmin                                                  2049
      pfx=(1.0-pmin)*pmax                                                  2049
      vmin=pfm*pmax                                                        2050
      bta=parm                                                             2050
      omb=1.0-bta                                                          2050
      dev1=0.0                                                             2050
      dev0=0.0                                                             2051
14770 do 14771 ic=1,nc                                                     2051
      q0=dot_product(w,y(:,ic))                                            2052
      if(q0 .gt. pmin)goto 14791                                           2052
      jerr =8000+ic                                                        2052
      return                                                               2052
14791 continue                                                             2053
      if(q0 .lt. 1.0-pmin)goto 14811                                       2053
      jerr =9000+ic                                                        2053
      return                                                               2053
14811 continue                                                             2054
      b(1:ni,ic)=0.0                                                       2054
      b(0,ic)=log(q0)                                                      2054
      dev1=dev1-q0*b(0,ic)                                                 2055
14771 continue                                                             2056
14772 continue                                                             2056
      iy=0                                                                 2056
      al=0.0                                                               2057
      if(nonzero(no*nc,g) .ne. 0)goto 14831                                2058
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         2058
      sxp=0.0                                                              2059
14840 do 14841 ic=1,nc                                                     2059
      q(:,ic)=exp(b(0,ic))                                                 2059
      sxp=sxp+q(:,ic)                                                      2059
14841 continue                                                             2060
14842 continue                                                             2060
      goto 14851                                                           2061
14831 continue                                                             2061
14860 do 14861 i=1,no                                                      2061
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2061
14861 continue                                                             2061
14862 continue                                                             2061
      sxp=0.0                                                              2062
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 2062
      if(jerr.ne.0) return                                                 2063
      dev1=0.0                                                             2064
14870 do 14871 ic=1,nc                                                     2064
      q(:,ic)=b(0,ic)+g(:,ic)                                              2065
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             2066
      q(:,ic)=exp(q(:,ic))                                                 2066
      sxp=sxp+q(:,ic)                                                      2067
14871 continue                                                             2068
14872 continue                                                             2068
      sxpl=w*log(sxp)                                                      2068
14880 do 14881 ic=1,nc                                                     2068
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  2068
14881 continue                                                             2069
14882 continue                                                             2069
14851 continue                                                             2070
14821 continue                                                             2070
14890 do 14891 ic=1,nc                                                     2070
14900 do 14901 i=1,no                                                      2070
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               2070
14901 continue                                                             2070
14902 continue                                                             2070
14891 continue                                                             2071
14892 continue                                                             2071
      dev0=dev0+dev1                                                       2072
      if(kopt .le. 0)goto 14921                                            2073
      if(isd .le. 0)goto 14941                                             2073
      xv=0.25                                                              2073
      goto 14951                                                           2074
14941 continue                                                             2075
14960 do 14961 j=1,ni                                                      2075
      if(ju(j).eq.0)goto 14961                                             2075
      jb=ix(j)                                                             2075
      je=ix(j+1)-1                                                         2076
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        2077
14961 continue                                                             2078
14962 continue                                                             2078
14951 continue                                                             2079
14931 continue                                                             2079
14921 continue                                                             2080
      if(flmin .ge. 1.0)goto 14981                                         2080
      eqs=max(eps,flmin)                                                   2080
      alf=eqs**(1.0/(nlam-1))                                              2080
14981 continue                                                             2081
      m=0                                                                  2081
      mm=0                                                                 2081
      nin=0                                                                2081
      nlp=0                                                                2081
      mnl=min(mnlam,nlam)                                                  2081
      bs=0.0                                                               2081
      svr=0.0                                                              2081
      o=0.0                                                                2082
      shr=shri*dev0                                                        2082
      ga=0.0                                                               2083
14990 do 14991 ic=1,nc                                                     2083
      v=q(:,ic)/sxp                                                        2083
      r=w*(y(:,ic)-v)                                                      2083
      v=w*v*(1.0-v)                                                        2084
15000 do 15001 j=1,ni                                                      2084
      if(ju(j).eq.0)goto 15001                                             2085
      jb=ix(j)                                                             2085
      je=ix(j+1)-1                                                         2085
      jn=ix(j+1)-ix(j)                                                     2086
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2087
      gj=dot_product(sc(1:jn),x(jb:je))                                    2088
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2089
15001 continue                                                             2090
15002 continue                                                             2090
14991 continue                                                             2091
14992 continue                                                             2091
15010 do 15011 ilm=1,nlam                                                  2091
      al0=al                                                               2092
      if(flmin .lt. 1.0)goto 15031                                         2092
      al=ulam(ilm)                                                         2092
      goto 15021                                                           2093
15031 if(ilm .le. 2)goto 15041                                             2093
      al=al*alf                                                            2093
      goto 15021                                                           2094
15041 if(ilm .ne. 1)goto 15051                                             2094
      al=big                                                               2094
      goto 15061                                                           2095
15051 continue                                                             2095
      al0=0.0                                                              2096
15070 do 15071 j=1,ni                                                      2096
      if(ju(j).eq.0)goto 15071                                             2096
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2096
15071 continue                                                             2097
15072 continue                                                             2097
      al0=al0/max(bta,1.0e-3)                                              2097
      al=alf*al0                                                           2098
15061 continue                                                             2099
15021 continue                                                             2099
      al2=al*omb                                                           2099
      al1=al*bta                                                           2099
      tlam=bta*(2.0*al-al0)                                                2100
15080 do 15081 k=1,ni                                                      2100
      if(iy(k).eq.1)goto 15081                                             2100
      if(ju(k).eq.0)goto 15081                                             2101
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      2102
15081 continue                                                             2103
15082 continue                                                             2103
10680 continue                                                             2104
15090 continue                                                             2104
15091 continue                                                             2104
      ixx=0                                                                2104
      jxx=ixx                                                              2104
      ig=0                                                                 2105
15100 do 15101 ic=1,nc                                                     2105
      bs(0,ic)=b(0,ic)                                                     2106
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          2107
      xm(0)=0.0                                                            2107
      svr=0.0                                                              2107
      o=0.0                                                                2108
15110 do 15111 i=1,no                                                      2108
      pic=q(i,ic)/sxp(i)                                                   2109
      if(pic .ge. pfm)goto 15131                                           2109
      pic=0.0                                                              2109
      v(i)=0.0                                                             2109
      goto 15121                                                           2110
15131 if(pic .le. pfx)goto 15141                                           2110
      pic=1.0                                                              2110
      v(i)=0.0                                                             2110
      goto 15151                                                           2111
15141 continue                                                             2111
      v(i)=w(i)*pic*(1.0-pic)                                              2111
      xm(0)=xm(0)+v(i)                                                     2111
15151 continue                                                             2112
15121 continue                                                             2112
      r(i)=w(i)*(y(i,ic)-pic)                                              2112
      svr=svr+r(i)                                                         2113
15111 continue                                                             2114
15112 continue                                                             2114
      if(xm(0).le.vmin)goto 15101                                          2114
      ig=1                                                                 2115
15160 do 15161 j=1,ni                                                      2115
      if(iy(j).eq.0)goto 15161                                             2116
      jb=ix(j)                                                             2116
      je=ix(j+1)-1                                                         2117
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             2118
      if(kopt .ne. 0)goto 15181                                            2119
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       2120
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          2121
15181 continue                                                             2122
15161 continue                                                             2123
15162 continue                                                             2123
15190 continue                                                             2123
15191 continue                                                             2123
      nlp=nlp+1                                                            2123
      dlx=0.0                                                              2124
15200 do 15201 k=1,ni                                                      2124
      if(iy(k).eq.0)goto 15201                                             2125
      jb=ix(k)                                                             2125
      je=ix(k+1)-1                                                         2125
      jn=ix(k+1)-ix(k)                                                     2125
      bk=b(k,ic)                                                           2126
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2127
      gk=dot_product(sc(1:jn),x(jb:je))                                    2128
      gk=(gk-svr*xb(k))/xs(k)                                              2129
      u=gk+xv(k,ic)*b(k,ic)                                                2129
      au=abs(u)-vp(k)*al1                                                  2130
      if(au .gt. 0.0)goto 15221                                            2130
      b(k,ic)=0.0                                                          2130
      goto 15231                                                           2131
15221 continue                                                             2131
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              2131
15231 continue                                                             2132
15211 continue                                                             2132
      d=b(k,ic)-bk                                                         2132
      if(abs(d).le.0.0)goto 15201                                          2133
      dlx=max(dlx,xv(k,ic)*d**2)                                           2134
      if(mm(k) .ne. 0)goto 15251                                           2134
      nin=nin+1                                                            2135
      if(nin .le. nx)goto 15271                                            2135
      jxx=1                                                                2135
      goto 15202                                                           2135
15271 continue                                                             2136
      mm(k)=nin                                                            2136
      m(nin)=k                                                             2137
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             2138
15251 continue                                                             2139
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2140
      o=o+d*(xb(k)/xs(k))                                                  2141
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2142
15201 continue                                                             2143
15202 continue                                                             2143
      if(jxx.gt.0)goto 15192                                               2144
      d=svr/xm(0)                                                          2145
      if(d .eq. 0.0)goto 15291                                             2145
      b(0,ic)=b(0,ic)+d                                                    2145
      dlx=max(dlx,xm(0)*d**2)                                              2146
      r=r-d*v                                                              2146
      svr=svr-d*xm(0)                                                      2147
15291 continue                                                             2148
      if(dlx.lt.shr)goto 15192                                             2148
      if(nlp .le. maxit)goto 15311                                         2148
      jerr=-ilm                                                            2148
      return                                                               2148
15311 continue                                                             2149
15320 continue                                                             2149
15321 continue                                                             2149
      nlp=nlp+1                                                            2149
      dlx=0.0                                                              2150
15330 do 15331 l=1,nin                                                     2150
      k=m(l)                                                               2150
      jb=ix(k)                                                             2150
      je=ix(k+1)-1                                                         2151
      jn=ix(k+1)-ix(k)                                                     2151
      bk=b(k,ic)                                                           2152
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2153
      gk=dot_product(sc(1:jn),x(jb:je))                                    2154
      gk=(gk-svr*xb(k))/xs(k)                                              2155
      u=gk+xv(k,ic)*b(k,ic)                                                2155
      au=abs(u)-vp(k)*al1                                                  2156
      if(au .gt. 0.0)goto 15351                                            2156
      b(k,ic)=0.0                                                          2156
      goto 15361                                                           2157
15351 continue                                                             2157
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              2157
15361 continue                                                             2158
15341 continue                                                             2158
      d=b(k,ic)-bk                                                         2158
      if(abs(d).le.0.0)goto 15331                                          2159
      dlx=max(dlx,xv(k,ic)*d**2)                                           2160
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2161
      o=o+d*(xb(k)/xs(k))                                                  2162
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2163
15331 continue                                                             2164
15332 continue                                                             2164
      d=svr/xm(0)                                                          2165
      if(d .eq. 0.0)goto 15381                                             2165
      b(0,ic)=b(0,ic)+d                                                    2165
      dlx=max(dlx,xm(0)*d**2)                                              2166
      r=r-d*v                                                              2166
      svr=svr-d*xm(0)                                                      2167
15381 continue                                                             2168
      if(dlx.lt.shr)goto 15322                                             2168
      if(nlp .le. maxit)goto 15401                                         2168
      jerr=-ilm                                                            2168
      return                                                               2168
15401 continue                                                             2169
      goto 15321                                                           2170
15322 continue                                                             2170
      goto 15191                                                           2171
15192 continue                                                             2171
      if(jxx.gt.0)goto 15102                                               2172
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                         2173
      if(ixx .ne. 0)goto 15421                                             2174
15430 do 15431 j=1,nin                                                     2174
      k=m(j)                                                               2175
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 15451                2175
      ixx=1                                                                2175
      goto 15432                                                           2175
15451 continue                                                             2176
15431 continue                                                             2177
15432 continue                                                             2177
15421 continue                                                             2178
      sc=b(0,ic)+g(:,ic)                                                   2178
      b0=0.0                                                               2179
15460 do 15461 j=1,nin                                                     2179
      l=m(j)                                                               2179
      jb=ix(l)                                                             2179
      je=ix(l+1)-1                                                         2180
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   2181
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            2182
15461 continue                                                             2183
15462 continue                                                             2183
      sc=min(max(exmn,sc+b0),exmx)                                         2184
      sxp=sxp-q(:,ic)                                                      2185
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          2186
      sxp=sxp+q(:,ic)                                                      2187
15101 continue                                                             2188
15102 continue                                                             2188
      s=-sum(b(0,:))/nc                                                    2188
      b(0,:)=b(0,:)+s                                                      2188
      sc=s                                                                 2188
      b0=0.0                                                               2189
15470 do 15471 j=1,nin                                                     2189
      l=m(j)                                                               2190
      if(vp(l) .gt. 0.0)goto 15491                                         2190
      s=sum(b(l,:))/nc                                                     2190
      goto 15501                                                           2191
15491 continue                                                             2191
      s=elc(parm,nc,b(l,:),is)                                             2191
15501 continue                                                             2192
15481 continue                                                             2192
      b(l,:)=b(l,:)-s                                                      2193
      jb=ix(l)                                                             2193
      je=ix(l+1)-1                                                         2194
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         2195
      b0=b0+s*xb(l)/xs(l)                                                  2196
15471 continue                                                             2197
15472 continue                                                             2197
      sc=sc+b0                                                             2197
      sc=exp(sc)                                                           2197
      sxp=sxp*sc                                                           2197
15510 do 15511 ic=1,nc                                                     2197
      q(:,ic)=q(:,ic)*sc                                                   2197
15511 continue                                                             2198
15512 continue                                                             2198
      if(jxx.gt.0)goto 15092                                               2198
      if(ig.eq.0)goto 15092                                                2199
      if(ixx .ne. 0)goto 15531                                             2200
15540 do 15541 j=1,ni                                                      2200
      if(iy(j).eq.1)goto 15541                                             2200
      if(ju(j).eq.0)goto 15541                                             2200
      ga(j)=0.0                                                            2200
15541 continue                                                             2201
15542 continue                                                             2201
15550 do 15551 ic=1,nc                                                     2201
      v=q(:,ic)/sxp                                                        2201
      r=w*(y(:,ic)-v)                                                      2201
      v=w*v*(1.0-v)                                                        2202
15560 do 15561 j=1,ni                                                      2202
      if(iy(j).eq.1)goto 15561                                             2202
      if(ju(j).eq.0)goto 15561                                             2203
      jb=ix(j)                                                             2203
      je=ix(j+1)-1                                                         2203
      jn=ix(j+1)-ix(j)                                                     2204
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2205
      gj=dot_product(sc(1:jn),x(jb:je))                                    2206
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2207
15561 continue                                                             2208
15562 continue                                                             2208
15551 continue                                                             2209
15552 continue                                                             2209
15570 do 15571 k=1,ni                                                      2209
      if(iy(k).eq.1)goto 15571                                             2209
      if(ju(k).eq.0)goto 15571                                             2210
      if(ga(k) .le. al1*vp(k))goto 15591                                   2210
      iy(k)=1                                                              2210
      ixx=1                                                                2210
15591 continue                                                             2211
15571 continue                                                             2212
15572 continue                                                             2212
      if(ixx.eq.1) go to 10680                                             2213
      goto 15092                                                           2214
15531 continue                                                             2215
      goto 15091                                                           2216
15092 continue                                                             2216
      if(jxx .le. 0)goto 15611                                             2216
      jerr=-10000-ilm                                                      2216
      goto 15012                                                           2216
15611 continue                                                             2216
      devi=0.0                                                             2217
15620 do 15621 ic=1,nc                                                     2218
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          2218
      a0(ic,ilm)=b(0,ic)                                                   2219
15630 do 15631 i=1,no                                                      2219
      if(y(i,ic).le.0.0)goto 15631                                         2220
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           2221
15631 continue                                                             2222
15632 continue                                                             2222
15621 continue                                                             2223
15622 continue                                                             2223
      kin(ilm)=nin                                                         2223
      alm(ilm)=al                                                          2223
      lmu=ilm                                                              2224
      dev(ilm)=(dev1-devi)/dev0                                            2224
      if(ig.eq.0)goto 15012                                                2225
      if(ilm.lt.mnl)goto 15011                                             2225
      if(flmin.ge.1.0)goto 15011                                           2226
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 15012             2227
      if(dev(ilm).gt.devmax)goto 15012                                     2227
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 15012                             2228
15011 continue                                                             2229
15012 continue                                                             2229
      g=log(q)                                                             2229
15640 do 15641 i=1,no                                                      2229
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2229
15641 continue                                                             2230
15642 continue                                                             2230
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy)                      2231
      return                                                               2232
      end                                                                  2233
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  2234
      real a0(nc),ca(nx,nc),x(*),f(nc,n)                                   2234
      integer ia(*),ix(*),jx(*)                                            2235
15650 do 15651 ic=1,nc                                                     2235
      f(ic,:)=a0(ic)                                                       2235
15651 continue                                                             2236
15652 continue                                                             2236
15660 do 15661 j=1,nin                                                     2236
      k=ia(j)                                                              2236
      kb=ix(k)                                                             2236
      ke=ix(k+1)-1                                                         2237
15670 do 15671 ic=1,nc                                                     2237
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    2237
15671 continue                                                             2238
15672 continue                                                             2238
15661 continue                                                             2239
15662 continue                                                             2239
      return                                                               2240
      end                                                                  2241
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,ne,nx,nlam,flmin,ula   2243 
     *m,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam)              2244
      real ca(nx,nlam),dev(nlam),alm(nlam)                                 2245
      integer jd(*),ia(nx),nin(nlam)                                       2246
      real, dimension (:), allocatable :: xs,ww,vq                              
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 15691                                    2250
      jerr=10000                                                           2250
      return                                                               2250
15691 continue                                                             2251
      allocate(ww(1:no),stat=jerr)                                         2252
      allocate(ju(1:ni),stat=ierr)                                         2252
      jerr=jerr+ierr                                                       2253
      allocate(vq(1:ni),stat=ierr)                                         2253
      jerr=jerr+ierr                                                       2254
      if(isd .le. 0)goto 15711                                             2254
      allocate(xs(1:ni),stat=ierr)                                         2254
      jerr=jerr+ierr                                                       2254
15711 continue                                                             2255
      if(jerr.ne.0) return                                                 2256
      call chkvars(no,ni,x,ju)                                             2257
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2258
      if(maxval(ju) .gt. 0)goto 15731                                      2258
      jerr=7777                                                            2258
      return                                                               2258
15731 continue                                                             2259
      vq=max(0.0,vp)                                                       2259
      vq=vq*ni/sum(vq)                                                     2260
      ww=max(0.0,w)                                                        2260
      sw=sum(ww)                                                           2261
      if(sw .gt. 0.0)goto 15751                                            2261
      jerr=9999                                                            2261
      return                                                               2261
15751 continue                                                             2261
      ww=ww/sw                                                             2262
      call cstandard(no,ni,x,ww,ju,isd,xs)                                 2263
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr   2265 
     *,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 2265
      dev0=2.0*sw*dev0                                                     2266
      if(isd .le. 0)goto 15771                                             2266
15780 do 15781 k=1,lmu                                                     2266
      nk=nin(k)                                                            2266
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   2266
15781 continue                                                             2266
15782 continue                                                             2266
15771 continue                                                             2267
      deallocate(ww,ju,vq)                                                 2267
      if(isd.gt.0) deallocate(xs)                                          2268
      return                                                               2269
      end                                                                  2270
      subroutine cstandard (no,ni,x,w,ju,isd,xs)                           2271
      real x(no,ni),w(no),xs(ni)                                           2271
      integer ju(ni)                                                       2272
15790 do 15791 j=1,ni                                                      2272
      if(ju(j).eq.0)goto 15791                                             2273
      xm=dot_product(w,x(:,j))                                             2273
      x(:,j)=x(:,j)-xm                                                     2274
      if(isd .le. 0)goto 15811                                             2274
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 2274
      x(:,j)=x(:,j)/xs(j)                                                  2274
15811 continue                                                             2275
15791 continue                                                             2276
15792 continue                                                             2276
      return                                                               2277
      end                                                                  2278
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,ne,nx,nlam,flmin,ula   2280 
     *m,cthri,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=0.001, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99)   2281
      real x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam)              2282
      real ao(nx,nlam),dev(nlam),alm(nlam)                                 2283
      integer ju(ni),m(nx),kin(nlam)                                       2284
      real, dimension (:), allocatable :: w,dk,v,xs,wr,a,as,f,dq                
      real, dimension (:), allocatable :: e,uu,ga                               
      integer, dimension (:), allocatable :: jp,kp,mm,ixx                       
      allocate(e(1:no),stat=jerr)                                          2290
      allocate(uu(1:no),stat=ierr)                                         2290
      jerr=jerr+ierr                                                       2291
      allocate(f(1:no),stat=ierr)                                          2291
      jerr=jerr+ierr                                                       2292
      allocate(w(1:no),stat=ierr)                                          2292
      jerr=jerr+ierr                                                       2293
      allocate(v(1:ni),stat=ierr)                                          2293
      jerr=jerr+ierr                                                       2294
      allocate(a(1:ni),stat=ierr)                                          2294
      jerr=jerr+ierr                                                       2295
      allocate(as(1:ni),stat=ierr)                                         2295
      jerr=jerr+ierr                                                       2296
      allocate(xs(1:ni),stat=ierr)                                         2296
      jerr=jerr+ierr                                                       2297
      allocate(ga(1:ni),stat=ierr)                                         2297
      jerr=jerr+ierr                                                       2298
      allocate(ixx(1:ni),stat=ierr)                                        2298
      jerr=jerr+ierr                                                       2299
      allocate(jp(1:no),stat=ierr)                                         2299
      jerr=jerr+ierr                                                       2300
      allocate(kp(1:no),stat=ierr)                                         2300
      jerr=jerr+ierr                                                       2301
      allocate(dk(1:no),stat=ierr)                                         2301
      jerr=jerr+ierr                                                       2302
      allocate(wr(1:no),stat=ierr)                                         2302
      jerr=jerr+ierr                                                       2303
      allocate(dq(1:no),stat=ierr)                                         2303
      jerr=jerr+ierr                                                       2304
      allocate(mm(1:ni),stat=ierr)                                         2304
      jerr=jerr+ierr                                                       2305
      if(jerr.ne.0)go to 11790                                             2306
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2307
      if(jerr.ne.0) go to 11790                                            2307
      alpha=parm                                                           2308
      oma=1.0-alpha                                                        2308
      nlm=0                                                                2308
      ixx=0                                                                2308
      al=0.0                                                               2309
      dq=d*q                                                               2309
      call died(no,nk,dq,kp,jp,dk)                                         2310
      a=0.0                                                                2310
      f(1)=0.0                                                             2310
      fmax=log(huge(f(1))*0.1)                                             2311
      if(nonzero(no,g) .eq. 0)goto 15831                                   2311
      f=g-dot_product(q,g)                                                 2312
      e=q*exp(sign(min(abs(f),fmax),f))                                    2313
      goto 15841                                                           2314
15831 continue                                                             2314
      f=0.0                                                                2314
      e=q                                                                  2314
15841 continue                                                             2315
15821 continue                                                             2315
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 2316
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         2316
      dev0=rr                                                              2317
15850 do 15851 i=1,no                                                      2317
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 15871                   2317
      w(i)=0.0                                                             2317
      wr(i)=w(i)                                                           2317
15871 continue                                                             2317
15851 continue                                                             2318
15852 continue                                                             2318
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2319
      if(jerr.ne.0) go to 11790                                            2320
      if(flmin .ge. 1.0)goto 15891                                         2320
      eqs=max(eps,flmin)                                                   2320
      alf=eqs**(1.0/(nlam-1))                                              2320
15891 continue                                                             2321
      m=0                                                                  2321
      mm=0                                                                 2321
      nlp=0                                                                2321
      nin=nlp                                                              2321
      mnl=min(mnlam,nlam)                                                  2321
      as=0.0                                                               2321
      cthr=cthri*dev0                                                      2322
15900 do 15901 j=1,ni                                                      2322
      if(ju(j).eq.0)goto 15901                                             2322
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2322
15901 continue                                                             2323
15902 continue                                                             2323
15910 do 15911 ilm=1,nlam                                                  2323
      al0=al                                                               2324
      if(flmin .lt. 1.0)goto 15931                                         2324
      al=ulam(ilm)                                                         2324
      goto 15921                                                           2325
15931 if(ilm .le. 2)goto 15941                                             2325
      al=al*alf                                                            2325
      goto 15921                                                           2326
15941 if(ilm .ne. 1)goto 15951                                             2326
      al=big                                                               2326
      goto 15961                                                           2327
15951 continue                                                             2327
      al0=0.0                                                              2328
15970 do 15971 j=1,ni                                                      2328
      if(ju(j).eq.0)goto 15971                                             2328
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2328
15971 continue                                                             2329
15972 continue                                                             2329
      al0=al0/max(parm,1.0e-3)                                             2329
      al=alf*al0                                                           2330
15961 continue                                                             2331
15921 continue                                                             2331
      sa=alpha*al                                                          2331
      omal=oma*al                                                          2331
      tlam=alpha*(2.0*al-al0)                                              2332
15980 do 15981 k=1,ni                                                      2332
      if(ixx(k).eq.1)goto 15981                                            2332
      if(ju(k).eq.0)goto 15981                                             2333
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2334
15981 continue                                                             2335
15982 continue                                                             2335
10680 continue                                                             2336
15990 continue                                                             2336
15991 continue                                                             2336
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2337
      call vars(no,ni,x,w,ixx,v)                                           2338
16000 continue                                                             2338
16001 continue                                                             2338
      nlp=nlp+1                                                            2338
      dli=0.0                                                              2339
16010 do 16011 j=1,ni                                                      2339
      if(ixx(j).eq.0)goto 16011                                            2340
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2341
      if(abs(u) .gt. vp(j)*sa)goto 16031                                   2341
      at=0.0                                                               2341
      goto 16041                                                           2342
16031 continue                                                             2342
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2342
16041 continue                                                             2343
16021 continue                                                             2343
      if(at .eq. a(j))goto 16061                                           2343
      del=at-a(j)                                                          2343
      a(j)=at                                                              2343
      dli=max(dli,v(j)*del**2)                                             2344
      wr=wr-del*w*x(:,j)                                                   2344
      f=f+del*x(:,j)                                                       2345
      if(mm(j) .ne. 0)goto 16081                                           2345
      nin=nin+1                                                            2345
      if(nin.gt.nx)goto 16012                                              2346
      mm(j)=nin                                                            2346
      m(nin)=j                                                             2347
16081 continue                                                             2348
16061 continue                                                             2349
16011 continue                                                             2350
16012 continue                                                             2350
      if(nin.gt.nx)goto 16002                                              2350
      if(dli.lt.cthr)goto 16002                                            2351
      if(nlp .le. maxit)goto 16101                                         2351
      jerr=-ilm                                                            2351
      return                                                               2351
16101 continue                                                             2352
16110 continue                                                             2352
16111 continue                                                             2352
      nlp=nlp+1                                                            2352
      dli=0.0                                                              2353
16120 do 16121 l=1,nin                                                     2353
      j=m(l)                                                               2354
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2355
      if(abs(u) .gt. vp(j)*sa)goto 16141                                   2355
      at=0.0                                                               2355
      goto 16151                                                           2356
16141 continue                                                             2356
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2356
16151 continue                                                             2357
16131 continue                                                             2357
      if(at .eq. a(j))goto 16171                                           2357
      del=at-a(j)                                                          2357
      a(j)=at                                                              2357
      dli=max(dli,v(j)*del**2)                                             2358
      wr=wr-del*w*x(:,j)                                                   2358
      f=f+del*x(:,j)                                                       2359
16171 continue                                                             2360
16121 continue                                                             2361
16122 continue                                                             2361
      if(dli.lt.cthr)goto 16112                                            2361
      if(nlp .le. maxit)goto 16191                                         2361
      jerr=-ilm                                                            2361
      return                                                               2361
16191 continue                                                             2362
      goto 16111                                                           2363
16112 continue                                                             2363
      goto 16001                                                           2364
16002 continue                                                             2364
      if(nin.gt.nx)goto 15992                                              2365
      e=q*exp(sign(min(abs(f),fmax),f))                                    2366
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2367
      if(jerr .eq. 0)goto 16211                                            2367
      jerr=jerr-ilm                                                        2367
      go to 11790                                                          2367
16211 continue                                                             2368
      ix=0                                                                 2369
16220 do 16221 j=1,nin                                                     2369
      k=m(j)                                                               2370
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 16221                           2370
      ix=1                                                                 2370
      goto 16222                                                           2370
16221 continue                                                             2371
16222 continue                                                             2371
      if(ix .ne. 0)goto 16241                                              2372
16250 do 16251 k=1,ni                                                      2372
      if(ixx(k).eq.1)goto 16251                                            2372
      if(ju(k).eq.0)goto 16251                                             2373
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2374
      if(ga(k) .le. sa*vp(k))goto 16271                                    2374
      ixx(k)=1                                                             2374
      ix=1                                                                 2374
16271 continue                                                             2375
16251 continue                                                             2376
16252 continue                                                             2376
      if(ix.eq.1) go to 10680                                              2377
      goto 15992                                                           2378
16241 continue                                                             2379
      goto 15991                                                           2380
15992 continue                                                             2380
      if(nin .le. nx)goto 16291                                            2380
      jerr=-10000-ilm                                                      2380
      goto 15912                                                           2380
16291 continue                                                             2381
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               2381
      kin(ilm)=nin                                                         2382
      alm(ilm)=al                                                          2382
      lmu=ilm                                                              2383
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   2384
      if(ilm.lt.mnl)goto 15911                                             2384
      if(flmin.ge.1.0)goto 15911                                           2385
      me=0                                                                 2385
16300 do 16301 j=1,nin                                                     2385
      if(ao(j,ilm).ne.0.0) me=me+1                                         2385
16301 continue                                                             2385
16302 continue                                                             2385
      if(me.gt.ne)goto 15912                                               2386
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 15912              2387
      if(dev(ilm).gt.devmax)goto 15912                                     2388
15911 continue                                                             2389
15912 continue                                                             2389
      g=f                                                                  2390
11790 continue                                                             2390
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)              2391
      return                                                               2392
      end                                                                  2393
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 2394
      real ca(nin),x(n,*),f(n)                                             2394
      integer ia(nin)                                                      2395
      f=0.0                                                                2395
      if(nin.le.0) return                                                  2396
16310 do 16311 i=1,n                                                       2396
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      2396
16311 continue                                                             2397
16312 continue                                                             2397
      return                                                               2398
      end                                                                  2399
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         2400
      real y(no),d(no),q(no)                                               2400
      integer jp(no),kp(*)                                                 2401
16320 do 16321 j=1,no                                                      2401
      jp(j)=j                                                              2401
16321 continue                                                             2401
16322 continue                                                             2401
      call psort7(y,jp,1,no)                                               2402
      nj=0                                                                 2402
16330 do 16331 j=1,no                                                      2402
      if(q(jp(j)).le.0.0)goto 16331                                        2402
      nj=nj+1                                                              2402
      jp(nj)=jp(j)                                                         2402
16331 continue                                                             2403
16332 continue                                                             2403
      if(nj .ne. 0)goto 16351                                              2403
      jerr=20000                                                           2403
      return                                                               2403
16351 continue                                                             2404
      j=1                                                                  2404
16360 continue                                                             2404
16361 if(d(jp(j)).gt.0.0)goto 16362                                        2404
      j=j+1                                                                2404
      if(j.gt.nj)goto 16362                                                2404
      goto 16361                                                           2405
16362 continue                                                             2405
      if(j .lt. nj-1)goto 16381                                            2405
      jerr=30000                                                           2405
      return                                                               2405
16381 continue                                                             2406
      j0=j-1                                                               2406
      nj=nj-j0                                                             2406
16390 do 16391 j=1,nj                                                      2406
      jp(j)=jp(j+j0)                                                       2406
16391 continue                                                             2407
16392 continue                                                             2407
      jerr=0                                                               2407
      nk=0                                                                 2407
      t0=y(jp(1))                                                          2407
      yk=t0                                                                2407
      j=2                                                                  2408
16400 continue                                                             2408
16401 continue                                                             2408
16410 continue                                                             2409
16411 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 16412                     2409
      j=j+1                                                                2409
      if(j.gt.nj)goto 16412                                                2409
      goto 16411                                                           2410
16412 continue                                                             2410
      nk=nk+1                                                              2410
      kp(nk)=j-1                                                           2410
      if(j.gt.nj)goto 16402                                                2411
      if(j .ne. nj)goto 16431                                              2411
      nk=nk+1                                                              2411
      kp(nk)=nj                                                            2411
      goto 16402                                                           2411
16431 continue                                                             2412
      yk=y(jp(j))                                                          2412
      j=j+1                                                                2413
      goto 16401                                                           2414
16402 continue                                                             2414
      return                                                               2415
      end                                                                  2416
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     2417
      real d(no),dk(nk),wr(no),w(no)                                       2418
      real e(no),u(no),b,c                                                 2418
      integer kp(nk),jp(no)                                                2419
      call usk(no,nk,kp,jp,e,u)                                            2420
      b=dk(1)/u(1)                                                         2420
      c=dk(1)/u(1)**2                                                      2420
      jerr=0                                                               2421
16440 do 16441 j=1,kp(1)                                                   2421
      i=jp(j)                                                              2422
      w(i)=e(i)*(b-e(i)*c)                                                 2422
      if(w(i) .gt. 0.0)goto 16461                                          2422
      jerr=-3                                                              2422
      return                                                               2422
16461 continue                                                             2423
      wr(i)=d(i)-e(i)*b                                                    2424
16441 continue                                                             2425
16442 continue                                                             2425
16470 do 16471 k=2,nk                                                      2425
      j1=kp(k-1)+1                                                         2425
      j2=kp(k)                                                             2426
      b=b+dk(k)/u(k)                                                       2426
      c=c+dk(k)/u(k)**2                                                    2427
16480 do 16481 j=j1,j2                                                     2427
      i=jp(j)                                                              2428
      w(i)=e(i)*(b-e(i)*c)                                                 2428
      if(w(i) .gt. 0.0)goto 16501                                          2428
      jerr=-30000                                                          2428
      return                                                               2428
16501 continue                                                             2429
      wr(i)=d(i)-e(i)*b                                                    2430
16481 continue                                                             2431
16482 continue                                                             2431
16471 continue                                                             2432
16472 continue                                                             2432
      return                                                               2433
      end                                                                  2434
      subroutine vars(no,ni,x,w,ixx,v)                                     2435
      real x(no,ni),w(no),v(ni)                                            2435
      integer ixx(ni)                                                      2436
16510 do 16511 j=1,ni                                                      2436
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                        2436
16511 continue                                                             2437
16512 continue                                                             2437
      return                                                               2438
      end                                                                  2439
      subroutine died(no,nk,d,kp,jp,dk)                                    2440
      real d(no),dk(nk)                                                    2440
      integer kp(nk),jp(no)                                                2441
      dk(1)=sum(d(jp(1:kp(1))))                                            2442
16520 do 16521 k=2,nk                                                      2442
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  2442
16521 continue                                                             2443
16522 continue                                                             2443
      return                                                               2444
      end                                                                  2445
      subroutine usk(no,nk,kp,jp,e,u)                                      2446
      real e(no),u(nk),h                                                   2446
      integer kp(nk),jp(no)                                                2447
      h=0.0                                                                2448
16530 do 16531 k=nk,1,-1                                                   2448
      j2=kp(k)                                                             2449
      j1=1                                                                 2449
      if(k.gt.1) j1=kp(k-1)+1                                              2450
16540 do 16541 j=j2,j1,-1                                                  2450
      h=h+e(jp(j))                                                         2450
16541 continue                                                             2451
16542 continue                                                             2451
      u(k)=h                                                               2452
16531 continue                                                             2453
16532 continue                                                             2453
      return                                                               2454
      end                                                                  2455
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             2456
      real d(no),dk(nk),f(no)                                              2457
      integer kp(nk),jp(no)                                                2457
      real e(no),u(nk),s                                                   2458
      call usk(no,nk,kp,jp,e,u)                                            2458
      u=log(u)                                                             2459
      risk=dot_product(d,f)-dot_product(dk,u)                              2460
      return                                                               2461
      end                                                                  2462
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 2463
      real x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(nlam)          2464
      real, dimension (:), allocatable :: dk,f,xm,dq,q                          
      real, dimension (:), allocatable :: e,uu                                  
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          2470
      allocate(q(1:no),stat=ierr)                                          2470
      jerr=jerr+ierr                                                       2471
      allocate(uu(1:no),stat=ierr)                                         2471
      jerr=jerr+ierr                                                       2472
      allocate(f(1:no),stat=ierr)                                          2472
      jerr=jerr+ierr                                                       2473
      allocate(dk(1:no),stat=ierr)                                         2473
      jerr=jerr+ierr                                                       2474
      allocate(jp(1:no),stat=ierr)                                         2474
      jerr=jerr+ierr                                                       2475
      allocate(kp(1:no),stat=ierr)                                         2475
      jerr=jerr+ierr                                                       2476
      allocate(dq(1:no),stat=ierr)                                         2476
      jerr=jerr+ierr                                                       2477
      allocate(xm(1:ni),stat=ierr)                                         2477
      jerr=jerr+ierr                                                       2478
      if(jerr.ne.0) go to 11790                                            2479
      q=max(0.0,w)                                                         2479
      sw=sum(q)                                                            2480
      if(sw .gt. 0.0)goto 16561                                            2480
      jerr=9999                                                            2480
      go to 11790                                                          2480
16561 continue                                                             2481
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2482
      if(jerr.ne.0) go to 11790                                            2482
      fmax=log(huge(e(1))*0.1)                                             2483
      dq=d*q                                                               2483
      call died(no,nk,dq,kp,jp,dk)                                         2483
      gm=dot_product(q,g)/sw                                               2484
16570 do 16571 j=1,ni                                                      2484
      xm(j)=dot_product(q,x(:,j))/sw                                       2484
16571 continue                                                             2485
16572 continue                                                             2485
16580 do 16581 lam=1,nlam                                                  2486
16590 do 16591 i=1,no                                                      2486
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       2487
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        2488
16591 continue                                                             2489
16592 continue                                                             2489
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          2490
16581 continue                                                             2491
16582 continue                                                             2491
11790 continue                                                             2491
      deallocate(e,uu,dk,f,jp,kp,dq)                                       2492
      return                                                               2493
      end                                                                  2494
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,ne,nx,nlam,flmin,ulam   2496 
     *,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    2497
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2498
      integer jd(*),ia(nx),nin(nlam)                                       2499
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16611                                    2503
      jerr=10000                                                           2503
      return                                                               2503
16611 continue                                                             2504
      if(minval(y) .ge. 0.0)goto 16631                                     2504
      jerr=8888                                                            2504
      return                                                               2504
16631 continue                                                             2505
      allocate(ww(1:no),stat=jerr)                                         2506
      allocate(ju(1:ni),stat=ierr)                                         2506
      jerr=jerr+ierr                                                       2507
      allocate(vq(1:ni),stat=ierr)                                         2507
      jerr=jerr+ierr                                                       2508
      allocate(xm(1:ni),stat=ierr)                                         2508
      jerr=jerr+ierr                                                       2509
      if(isd .le. 0)goto 16651                                             2509
      allocate(xs(1:ni),stat=ierr)                                         2509
      jerr=jerr+ierr                                                       2509
16651 continue                                                             2510
      if(jerr.ne.0) return                                                 2511
      call chkvars(no,ni,x,ju)                                             2512
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2513
      if(maxval(ju) .gt. 0)goto 16671                                      2513
      jerr=7777                                                            2513
      go to 11790                                                          2513
16671 continue                                                             2514
      vq=max(0.0,vp)                                                       2514
      vq=vq*ni/sum(vq)                                                     2515
      ww=max(0.0,w)                                                        2515
      sw=sum(ww)                                                           2515
      if(sw .gt. 0.0)goto 16691                                            2515
      jerr=9999                                                            2515
      go to 11790                                                          2515
16691 continue                                                             2516
      ww=ww/sw                                                             2517
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             2518
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr,   2520 
     *  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11790                                            2520
      dev0=2.0*sw*dev0                                                     2521
16700 do 16701 k=1,lmu                                                     2521
      nk=nin(k)                                                            2522
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2523
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2524
16701 continue                                                             2525
16702 continue                                                             2525
11790 continue                                                             2525
      deallocate(ww,ju,vq,xm)                                              2525
      if(isd.gt.0) deallocate(xs)                                          2526
      return                                                               2527
      end                                                                  2528
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,ne,nx,nlam,flmin,ulam   2530 
     *,shri,  isd,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2531 
     *9)
      real x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)                    2532
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2533
      integer ju(ni),m(nx),kin(nlam)                                       2534
      real, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga                    
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(a(1:ni),stat=jerr)                                          2539
      allocate(as(1:ni),stat=ierr)                                         2539
      jerr=jerr+ierr                                                       2540
      allocate(t(1:no),stat=ierr)                                          2540
      jerr=jerr+ierr                                                       2541
      allocate(mm(1:ni),stat=ierr)                                         2541
      jerr=jerr+ierr                                                       2542
      allocate(ga(1:ni),stat=ierr)                                         2542
      jerr=jerr+ierr                                                       2543
      allocate(ixx(1:ni),stat=ierr)                                        2543
      jerr=jerr+ierr                                                       2544
      allocate(wr(1:no),stat=ierr)                                         2544
      jerr=jerr+ierr                                                       2545
      allocate(v(1:ni),stat=ierr)                                          2545
      jerr=jerr+ierr                                                       2546
      allocate(w(1:no),stat=ierr)                                          2546
      jerr=jerr+ierr                                                       2547
      allocate(f(1:no),stat=ierr)                                          2547
      jerr=jerr+ierr                                                       2548
      if(jerr.ne.0) return                                                 2549
      bta=parm                                                             2549
      omb=1.0-bta                                                          2550
      t=q*y                                                                2550
      yb=sum(t)                                                            2550
      fmax=log(huge(bta)*0.1)                                              2551
      if(nonzero(no,g) .ne. 0)goto 16721                                   2551
      w=q*yb                                                               2551
      az=log(yb)                                                           2551
      f=az                                                                 2551
      dv0=yb*(log(yb)-1.0)                                                 2551
      goto 16731                                                           2552
16721 continue                                                             2552
      w=q*exp(sign(min(abs(g),fmax),g))                                    2552
      v0=sum(w)                                                            2552
      eaz=yb/v0                                                            2553
      w=eaz*w                                                              2553
      az=log(eaz)                                                          2553
      f=az+g                                                               2554
      dv0=dot_product(t,g)-yb*(1.0-az)                                     2555
16731 continue                                                             2556
16711 continue                                                             2556
      a=0.0                                                                2556
      as=0.0                                                               2556
      wr=t-w                                                               2556
      v0=yb                                                                2556
      dvr=-yb                                                              2557
16740 do 16741 i=1,no                                                      2557
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               2557
16741 continue                                                             2557
16742 continue                                                             2557
      dvr=dvr-dv0                                                          2557
      dev0=dvr                                                             2558
      if(flmin .ge. 1.0)goto 16761                                         2558
      eqs=max(eps,flmin)                                                   2558
      alf=eqs**(1.0/(nlam-1))                                              2558
16761 continue                                                             2559
      m=0                                                                  2559
      mm=0                                                                 2559
      nlp=0                                                                2559
      nin=nlp                                                              2559
      mnl=min(mnlam,nlam)                                                  2559
      shr=shri*dev0                                                        2559
      ixx=0                                                                2559
      al=0.0                                                               2560
16770 do 16771 j=1,ni                                                      2560
      if(ju(j).eq.0)goto 16771                                             2560
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2560
16771 continue                                                             2561
16772 continue                                                             2561
16780 do 16781 ilm=1,nlam                                                  2561
      al0=al                                                               2562
      if(flmin .lt. 1.0)goto 16801                                         2562
      al=ulam(ilm)                                                         2562
      goto 16791                                                           2563
16801 if(ilm .le. 2)goto 16811                                             2563
      al=al*alf                                                            2563
      goto 16791                                                           2564
16811 if(ilm .ne. 1)goto 16821                                             2564
      al=big                                                               2564
      goto 16831                                                           2565
16821 continue                                                             2565
      al0=0.0                                                              2566
16840 do 16841 j=1,ni                                                      2566
      if(ju(j).eq.0)goto 16841                                             2566
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2566
16841 continue                                                             2567
16842 continue                                                             2567
      al0=al0/max(bta,1.0e-3)                                              2567
      al=alf*al0                                                           2568
16831 continue                                                             2569
16791 continue                                                             2569
      al2=al*omb                                                           2569
      al1=al*bta                                                           2569
      tlam=bta*(2.0*al-al0)                                                2570
16850 do 16851 k=1,ni                                                      2570
      if(ixx(k).eq.1)goto 16851                                            2570
      if(ju(k).eq.0)goto 16851                                             2571
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2572
16851 continue                                                             2573
16852 continue                                                             2573
10680 continue                                                             2574
16860 continue                                                             2574
16861 continue                                                             2574
      az0=az                                                               2575
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2576
16870 do 16871 j=1,ni                                                      2576
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                        2576
16871 continue                                                             2577
16872 continue                                                             2577
16880 continue                                                             2577
16881 continue                                                             2577
      nlp=nlp+1                                                            2577
      dlx=0.0                                                              2578
16890 do 16891 k=1,ni                                                      2578
      if(ixx(k).eq.0)goto 16891                                            2578
      ak=a(k)                                                              2579
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2579
      au=abs(u)-vp(k)*al1                                                  2580
      if(au .gt. 0.0)goto 16911                                            2580
      a(k)=0.0                                                             2580
      goto 16921                                                           2581
16911 continue                                                             2581
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2581
16921 continue                                                             2582
16901 continue                                                             2582
      if(a(k).eq.ak)goto 16891                                             2582
      d=a(k)-ak                                                            2582
      dlx=max(dlx,v(k)*d**2)                                               2583
      wr=wr-d*w*x(:,k)                                                     2583
      f=f+d*x(:,k)                                                         2584
      if(mm(k) .ne. 0)goto 16941                                           2584
      nin=nin+1                                                            2584
      if(nin.gt.nx)goto 16892                                              2585
      mm(k)=nin                                                            2585
      m(nin)=k                                                             2586
16941 continue                                                             2587
16891 continue                                                             2588
16892 continue                                                             2588
      if(nin.gt.nx)goto 16882                                              2588
      d=sum(wr)/v0                                                         2589
      az=az+d                                                              2589
      dlx=max(dlx,v0*d**2)                                                 2589
      wr=wr-d*w                                                            2589
      f=f+d                                                                2590
      if(dlx.lt.shr)goto 16882                                             2590
      if(nlp .le. maxit)goto 16961                                         2590
      jerr=-ilm                                                            2590
      return                                                               2590
16961 continue                                                             2591
16970 continue                                                             2591
16971 continue                                                             2591
      nlp=nlp+1                                                            2591
      dlx=0.0                                                              2592
16980 do 16981 l=1,nin                                                     2592
      k=m(l)                                                               2592
      ak=a(k)                                                              2593
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2593
      au=abs(u)-vp(k)*al1                                                  2594
      if(au .gt. 0.0)goto 17001                                            2594
      a(k)=0.0                                                             2594
      goto 17011                                                           2595
17001 continue                                                             2595
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2595
17011 continue                                                             2596
16991 continue                                                             2596
      if(a(k).eq.ak)goto 16981                                             2596
      d=a(k)-ak                                                            2596
      dlx=max(dlx,v(k)*d**2)                                               2597
      wr=wr-d*w*x(:,k)                                                     2597
      f=f+d*x(:,k)                                                         2599
16981 continue                                                             2599
16982 continue                                                             2599
      d=sum(wr)/v0                                                         2599
      az=az+d                                                              2599
      dlx=max(dlx,v0*d**2)                                                 2599
      wr=wr-d*w                                                            2599
      f=f+d                                                                2600
      if(dlx.lt.shr)goto 16972                                             2600
      if(nlp .le. maxit)goto 17031                                         2600
      jerr=-ilm                                                            2600
      return                                                               2600
17031 continue                                                             2601
      goto 16971                                                           2602
16972 continue                                                             2602
      goto 16881                                                           2603
16882 continue                                                             2603
      if(nin.gt.nx)goto 16862                                              2604
      w=q*exp(sign(min(abs(f),fmax),f))                                    2604
      v0=sum(w)                                                            2604
      wr=t-w                                                               2605
      if(v0*(az-az0)**2 .ge. shr)goto 17051                                2605
      ix=0                                                                 2606
17060 do 17061 j=1,nin                                                     2606
      k=m(j)                                                               2607
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 17061                            2607
      ix=1                                                                 2607
      goto 17062                                                           2608
17061 continue                                                             2609
17062 continue                                                             2609
      if(ix .ne. 0)goto 17081                                              2610
17090 do 17091 k=1,ni                                                      2610
      if(ixx(k).eq.1)goto 17091                                            2610
      if(ju(k).eq.0)goto 17091                                             2611
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2612
      if(ga(k) .le. al1*vp(k))goto 17111                                   2612
      ixx(k)=1                                                             2612
      ix=1                                                                 2612
17111 continue                                                             2613
17091 continue                                                             2614
17092 continue                                                             2614
      if(ix.eq.1) go to 10680                                              2615
      goto 16862                                                           2616
17081 continue                                                             2617
17051 continue                                                             2618
      goto 16861                                                           2619
16862 continue                                                             2619
      if(nin .le. nx)goto 17131                                            2619
      jerr=-10000-ilm                                                      2619
      goto 16782                                                           2619
17131 continue                                                             2620
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2620
      kin(ilm)=nin                                                         2621
      a0(ilm)=az                                                           2621
      alm(ilm)=al                                                          2621
      lmu=ilm                                                              2622
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               2623
      if(ilm.lt.mnl)goto 16781                                             2623
      if(flmin.ge.1.0)goto 16781                                           2624
      me=0                                                                 2624
17140 do 17141 j=1,nin                                                     2624
      if(ca(j,ilm).ne.0.0) me=me+1                                         2624
17141 continue                                                             2624
17142 continue                                                             2624
      if(me.gt.ne)goto 16782                                               2625
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 16782              2626
      if(dev(ilm).gt.devmax)goto 16782                                     2627
16781 continue                                                             2628
16782 continue                                                             2628
      g=f                                                                  2629
11790 continue                                                             2629
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                                2630
      return                                                               2631
      end                                                                  2632
      function nonzero(n,v)                                                2633
      real v(n)                                                            2634
      nonzero=0                                                            2634
17150 do 17151 i=1,n                                                       2634
      if(v(i) .eq. 0.0)goto 17171                                          2634
      nonzero=1                                                            2634
      return                                                               2634
17171 continue                                                             2634
17151 continue                                                             2635
17152 continue                                                             2635
      return                                                               2636
      end                                                                  2637
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               2638
      real a(nx,lmu),b(ni,lmu)                                             2638
      integer ia(nx),nin(lmu)                                              2639
17180 do 17181 lam=1,lmu                                                   2639
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        2639
17181 continue                                                             2640
17182 continue                                                             2640
      return                                                               2641
      end                                                                  2642
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           2643
      real a(nx,nc,lmu),b(ni,nc,lmu)                                       2643
      integer ia(nx),nin(lmu)                                              2644
17190 do 17191 lam=1,lmu                                                   2644
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             2644
17191 continue                                                             2645
17192 continue                                                             2645
      return                                                               2646
      end                                                                  2647
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               2648
      real x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)       2649
      real, dimension (:), allocatable :: w                                     
      if(minval(y) .ge. 0.0)goto 17211                                     2652
      jerr=8888                                                            2652
      return                                                               2652
17211 continue                                                             2653
      allocate(w(1:no),stat=jerr)                                          2653
      if(jerr.ne.0) return                                                 2654
      w=max(0.0,q)                                                         2654
      sw=sum(w)                                                            2654
      if(sw .gt. 0.0)goto 17231                                            2654
      jerr=9999                                                            2654
      go to 11790                                                          2654
17231 continue                                                             2655
      yb=dot_product(w,y)/sw                                               2655
      fmax=log(huge(y(1))*0.1)                                             2656
17240 do 17241 lam=1,nlam                                                  2656
      s=0.0                                                                2657
17250 do 17251 i=1,no                                                      2657
      if(w(i).le.0.0)goto 17251                                            2658
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          2659
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      2660
17251 continue                                                             2661
17252 continue                                                             2661
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2662
17241 continue                                                             2663
17242 continue                                                             2663
11790 continue                                                             2663
      deallocate(w)                                                        2664
      return                                                               2665
      end                                                                  2666
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,ne,nx,nlam,fl   2668 
     *min,ulam,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        2669
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2670
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           2671
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 17271                                    2675
      jerr=10000                                                           2675
      return                                                               2675
17271 continue                                                             2676
      if(minval(y) .ge. 0.0)goto 17291                                     2676
      jerr=8888                                                            2676
      return                                                               2676
17291 continue                                                             2677
      allocate(ww(1:no),stat=jerr)                                         2678
      allocate(ju(1:ni),stat=ierr)                                         2678
      jerr=jerr+ierr                                                       2679
      allocate(vq(1:ni),stat=ierr)                                         2679
      jerr=jerr+ierr                                                       2680
      allocate(xm(1:ni),stat=ierr)                                         2680
      jerr=jerr+ierr                                                       2681
      allocate(xs(1:ni),stat=ierr)                                         2681
      jerr=jerr+ierr                                                       2682
      if(jerr.ne.0) return                                                 2683
      call spchkvars(no,ni,x,ix,ju)                                        2684
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2685
      if(maxval(ju) .gt. 0)goto 17311                                      2685
      jerr=7777                                                            2685
      go to 11790                                                          2685
17311 continue                                                             2686
      vq=max(0.0,vp)                                                       2686
      vq=vq*ni/sum(vq)                                                     2687
      ww=max(0.0,w)                                                        2687
      sw=sum(ww)                                                           2687
      if(sw .gt. 0.0)goto 17331                                            2687
      jerr=9999                                                            2687
      go to 11790                                                          2687
17331 continue                                                             2688
      ww=ww/sw                                                             2689
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     2690
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmin,u   2692 
     *lam,thr,  isd,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11790                                            2692
      dev0=2.0*sw*dev0                                                     2693
17340 do 17341 k=1,lmu                                                     2693
      nk=nin(k)                                                            2694
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2695
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2696
17341 continue                                                             2697
17342 continue                                                             2697
11790 continue                                                             2697
      deallocate(ww,ju,vq,xm,xs)                                           2698
      return                                                               2699
      end                                                                  2700
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,ne,nx,nlam,fl   2702 
     *min,ulam,  shri,isd,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,j
     *err)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2703 
     *9)
      real x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),xs(ni)          2704
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2705
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2706
      real, dimension (:), allocatable :: qy,t,w,wr,v,a,as,xm,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(a(1:ni),stat=jerr)                                          2711
      allocate(as(1:ni),stat=ierr)                                         2711
      jerr=jerr+ierr                                                       2712
      allocate(t(1:no),stat=ierr)                                          2712
      jerr=jerr+ierr                                                       2713
      allocate(mm(1:ni),stat=ierr)                                         2713
      jerr=jerr+ierr                                                       2714
      allocate(ga(1:ni),stat=ierr)                                         2714
      jerr=jerr+ierr                                                       2715
      allocate(ixx(1:ni),stat=ierr)                                        2715
      jerr=jerr+ierr                                                       2716
      allocate(wr(1:no),stat=ierr)                                         2716
      jerr=jerr+ierr                                                       2717
      allocate(v(1:ni),stat=ierr)                                          2717
      jerr=jerr+ierr                                                       2718
      allocate(xm(1:ni),stat=ierr)                                         2718
      jerr=jerr+ierr                                                       2719
      allocate(w(1:no),stat=ierr)                                          2719
      jerr=jerr+ierr                                                       2720
      allocate(qy(1:no),stat=ierr)                                         2720
      jerr=jerr+ierr                                                       2721
      if(jerr.ne.0) return                                                 2722
      bta=parm                                                             2722
      omb=1.0-bta                                                          2722
      fmax=log(huge(bta)*0.1)                                              2723
      qy=q*y                                                               2723
      yb=sum(qy)                                                           2724
      if(nonzero(no,g) .ne. 0)goto 17361                                   2724
      w=q*yb                                                               2724
      az=log(yb)                                                           2724
      uu=az                                                                2725
      xm=yb*xb                                                             2725
      t=0.0                                                                2725
      dv0=yb*(log(yb)-1.0)                                                 2726
      goto 17371                                                           2727
17361 continue                                                             2727
      w=q*exp(sign(min(abs(g),fmax),g))                                    2727
      ww=sum(w)                                                            2727
      eaz=yb/ww                                                            2728
      w=eaz*w                                                              2728
      az=log(eaz)                                                          2728
      uu=az                                                                2728
      t=g                                                                  2728
      dv0=dot_product(qy,g)-yb*(1.0-az)                                    2729
17380 do 17381 j=1,ni                                                      2729
      if(ju(j).eq.0)goto 17381                                             2729
      jb=ix(j)                                                             2729
      je=ix(j+1)-1                                                         2730
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2731
17381 continue                                                             2732
17382 continue                                                             2732
17371 continue                                                             2733
17351 continue                                                             2733
      tt=yb*uu                                                             2733
      ww=yb                                                                2733
      wr=qy-q*(yb*(1.0-uu))                                                2733
      a=0.0                                                                2733
      as=0.0                                                               2734
      dvr=-yb                                                              2735
17390 do 17391 i=1,no                                                      2735
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             2735
17391 continue                                                             2735
17392 continue                                                             2735
      dvr=dvr-dv0                                                          2735
      dev0=dvr                                                             2736
      if(flmin .ge. 1.0)goto 17411                                         2736
      eqs=max(eps,flmin)                                                   2736
      alf=eqs**(1.0/(nlam-1))                                              2736
17411 continue                                                             2737
      m=0                                                                  2737
      mm=0                                                                 2737
      nlp=0                                                                2737
      nin=nlp                                                              2737
      mnl=min(mnlam,nlam)                                                  2737
      shr=shri*dev0                                                        2737
      al=0.0                                                               2737
      ixx=0                                                                2738
17420 do 17421 j=1,ni                                                      2738
      if(ju(j).eq.0)goto 17421                                             2739
      jb=ix(j)                                                             2739
      je=ix(j+1)-1                                                         2740
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   2742 
     *)-xb(j)*tt)/xs(j)
17421 continue                                                             2743
17422 continue                                                             2743
17430 do 17431 ilm=1,nlam                                                  2743
      al0=al                                                               2744
      if(flmin .lt. 1.0)goto 17451                                         2744
      al=ulam(ilm)                                                         2744
      goto 17441                                                           2745
17451 if(ilm .le. 2)goto 17461                                             2745
      al=al*alf                                                            2745
      goto 17441                                                           2746
17461 if(ilm .ne. 1)goto 17471                                             2746
      al=big                                                               2746
      goto 17481                                                           2747
17471 continue                                                             2747
      al0=0.0                                                              2748
17490 do 17491 j=1,ni                                                      2748
      if(ju(j).eq.0)goto 17491                                             2748
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2748
17491 continue                                                             2749
17492 continue                                                             2749
      al0=al0/max(bta,1.0e-3)                                              2749
      al=alf*al0                                                           2750
17481 continue                                                             2751
17441 continue                                                             2751
      al2=al*omb                                                           2751
      al1=al*bta                                                           2751
      tlam=bta*(2.0*al-al0)                                                2752
17500 do 17501 k=1,ni                                                      2752
      if(ixx(k).eq.1)goto 17501                                            2752
      if(ju(k).eq.0)goto 17501                                             2753
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2754
17501 continue                                                             2755
17502 continue                                                             2755
10680 continue                                                             2756
17510 continue                                                             2756
17511 continue                                                             2756
      az0=az                                                               2757
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2758
17520 do 17521 j=1,ni                                                      2758
      if(ixx(j).eq.0)goto 17521                                            2758
      jb=ix(j)                                                             2758
      je=ix(j+1)-1                                                         2759
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2760
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   2762 
     *b(j)**2)/xs(j)**2
17521 continue                                                             2763
17522 continue                                                             2763
17530 continue                                                             2763
17531 continue                                                             2763
      nlp=nlp+1                                                            2764
      dlx=0.0                                                              2765
17540 do 17541 k=1,ni                                                      2765
      if(ixx(k).eq.0)goto 17541                                            2765
      jb=ix(k)                                                             2765
      je=ix(k+1)-1                                                         2765
      ak=a(k)                                                              2766
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2768 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2769
      if(au .gt. 0.0)goto 17561                                            2769
      a(k)=0.0                                                             2769
      goto 17571                                                           2770
17561 continue                                                             2770
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2770
17571 continue                                                             2771
17551 continue                                                             2771
      if(a(k).eq.ak)goto 17541                                             2772
      if(mm(k) .ne. 0)goto 17591                                           2772
      nin=nin+1                                                            2772
      if(nin.gt.nx)goto 17542                                              2773
      mm(k)=nin                                                            2773
      m(nin)=k                                                             2774
17591 continue                                                             2775
      d=a(k)-ak                                                            2775
      dlx=max(dlx,v(k)*d**2)                                               2775
      dv=d/xs(k)                                                           2776
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2777
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2778
      uu=uu-dv*xb(k)                                                       2778
      tt=tt-dv*xm(k)                                                       2779
17541 continue                                                             2780
17542 continue                                                             2780
      if(nin.gt.nx)goto 17532                                              2780
      d=tt/ww-uu                                                           2781
      az=az+d                                                              2781
      dlx=max(dlx,ww*d**2)                                                 2781
      uu=uu+d                                                              2782
      if(dlx.lt.shr)goto 17532                                             2782
      if(nlp .le. maxit)goto 17611                                         2782
      jerr=-ilm                                                            2782
      return                                                               2782
17611 continue                                                             2783
17620 continue                                                             2783
17621 continue                                                             2783
      nlp=nlp+1                                                            2783
      dlx=0.0                                                              2784
17630 do 17631 l=1,nin                                                     2784
      k=m(l)                                                               2785
      jb=ix(k)                                                             2785
      je=ix(k+1)-1                                                         2785
      ak=a(k)                                                              2786
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2788 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2789
      if(au .gt. 0.0)goto 17651                                            2789
      a(k)=0.0                                                             2789
      goto 17661                                                           2790
17651 continue                                                             2790
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2790
17661 continue                                                             2791
17641 continue                                                             2791
      if(a(k).eq.ak)goto 17631                                             2791
      d=a(k)-ak                                                            2791
      dlx=max(dlx,v(k)*d**2)                                               2792
      dv=d/xs(k)                                                           2792
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2793
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2794
      uu=uu-dv*xb(k)                                                       2794
      tt=tt-dv*xm(k)                                                       2795
17631 continue                                                             2796
17632 continue                                                             2796
      d=tt/ww-uu                                                           2796
      az=az+d                                                              2796
      dlx=max(dlx,ww*d**2)                                                 2796
      uu=uu+d                                                              2797
      if(dlx.lt.shr)goto 17622                                             2797
      if(nlp .le. maxit)goto 17681                                         2797
      jerr=-ilm                                                            2797
      return                                                               2797
17681 continue                                                             2798
      goto 17621                                                           2799
17622 continue                                                             2799
      goto 17531                                                           2800
17532 continue                                                             2800
      if(nin.gt.nx)goto 17512                                              2801
      euu=exp(sign(min(abs(uu),fmax),uu))                                  2802
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                2802
      ww=sum(w)                                                            2803
      wr=qy-w*(1.0-uu)                                                     2803
      tt=sum(wr)                                                           2804
      if(ww*(az-az0)**2 .ge. shr)goto 17701                                2804
      kx=0                                                                 2805
17710 do 17711 j=1,nin                                                     2805
      k=m(j)                                                               2806
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 17711                            2806
      kx=1                                                                 2806
      goto 17712                                                           2807
17711 continue                                                             2808
17712 continue                                                             2808
      if(kx .ne. 0)goto 17731                                              2809
17740 do 17741 j=1,ni                                                      2809
      if(ixx(j).eq.1)goto 17741                                            2809
      if(ju(j).eq.0)goto 17741                                             2810
      jb=ix(j)                                                             2810
      je=ix(j+1)-1                                                         2811
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2812
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   2814 
     *)-xb(j)*tt)/xs(j)
      if(ga(j) .le. al1*vp(j))goto 17761                                   2814
      ixx(j)=1                                                             2814
      kx=1                                                                 2814
17761 continue                                                             2815
17741 continue                                                             2816
17742 continue                                                             2816
      if(kx.eq.1) go to 10680                                              2817
      goto 17512                                                           2818
17731 continue                                                             2819
17701 continue                                                             2820
      goto 17511                                                           2821
17512 continue                                                             2821
      if(nin .le. nx)goto 17781                                            2821
      jerr=-10000-ilm                                                      2821
      goto 17432                                                           2821
17781 continue                                                             2822
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2822
      kin(ilm)=nin                                                         2823
      a0(ilm)=az                                                           2823
      alm(ilm)=al                                                          2823
      lmu=ilm                                                              2824
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        2825
      if(ilm.lt.mnl)goto 17431                                             2825
      if(flmin.ge.1.0)goto 17431                                           2826
      me=0                                                                 2826
17790 do 17791 j=1,nin                                                     2826
      if(ca(j,ilm).ne.0.0) me=me+1                                         2826
17791 continue                                                             2826
17792 continue                                                             2826
      if(me.gt.ne)goto 17432                                               2827
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17432              2828
      if(dev(ilm).gt.devmax)goto 17432                                     2829
17431 continue                                                             2830
17432 continue                                                             2830
      g=t+uu                                                               2831
11790 continue                                                             2831
      deallocate(t,w,wr,v,a,qy,xm,as,mm,ga,ixx)                            2832
      return                                                               2833
      end                                                                  2834
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       2835
      real x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)           2836
      integer ix(*),jx(*)                                                  2837
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 17811                                     2840
      jerr=8888                                                            2840
      return                                                               2840
17811 continue                                                             2841
      allocate(w(1:no),stat=jerr)                                          2842
      allocate(f(1:no),stat=ierr)                                          2842
      jerr=jerr+ierr                                                       2843
      if(jerr.ne.0) return                                                 2844
      w=max(0.0,q)                                                         2844
      sw=sum(w)                                                            2844
      if(sw .gt. 0.0)goto 17831                                            2844
      jerr=9999                                                            2844
      go to 11790                                                          2844
17831 continue                                                             2845
      yb=dot_product(w,y)/sw                                               2845
      fmax=log(huge(y(1))*0.1)                                             2846
17840 do 17841 lam=1,nlam                                                  2846
      f=a0(lam)                                                            2847
17850 do 17851 j=1,ni                                                      2847
      if(a(j,lam).eq.0.0)goto 17851                                        2847
      jb=ix(j)                                                             2847
      je=ix(j+1)-1                                                         2848
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          2849
17851 continue                                                             2850
17852 continue                                                             2850
      f=f+g                                                                2851
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2852
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2853
17841 continue                                                             2854
17842 continue                                                             2854
11790 continue                                                             2854
      deallocate(w,f)                                                      2855
      return                                                               2856
      end                                                                  2857
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   2858 
     *jerr)
      real x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(nlam)          2859
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 2860
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 17871                                     2863
      jerr=8888                                                            2863
      return                                                               2863
17871 continue                                                             2864
      allocate(w(1:no),stat=jerr)                                          2865
      allocate(f(1:no),stat=ierr)                                          2865
      jerr=jerr+ierr                                                       2866
      if(jerr.ne.0) return                                                 2867
      w=max(0.0,q)                                                         2867
      sw=sum(w)                                                            2867
      if(sw .gt. 0.0)goto 17891                                            2867
      jerr=9999                                                            2867
      go to 11790                                                          2867
17891 continue                                                             2868
      yb=dot_product(w,y)/sw                                               2868
      fmax=log(huge(y(1))*0.1)                                             2869
17900 do 17901 lam=1,nlam                                                  2869
      f=a0(lam)                                                            2870
17910 do 17911 k=1,nin(lam)                                                2870
      j=ia(k)                                                              2870
      jb=ix(j)                                                             2870
      je=ix(j+1)-1                                                         2871
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         2872
17911 continue                                                             2873
17912 continue                                                             2873
      f=f+g                                                                2874
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2875
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2876
17901 continue                                                             2877
17902 continue                                                             2877
11790 continue                                                             2877
      deallocate(w,f)                                                      2878
      return                                                               2879
      end                                                                  2880
      subroutine multelnet  (parm,no,ni,nr,x,y,w,jd,vp,ne,nx,nlam,flmin,   2883 
     *ulam,thr,isd,jsd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no,nr),w(no),vp(ni),ca(nx,nr,nlam)                   2884
      real ulam(nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                      2885
      integer jd(*),ia(nx),nin(nlam)                                       2886
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 17931                                    2889
      jerr=10000                                                           2889
      return                                                               2889
17931 continue                                                             2890
      allocate(vq(1:ni),stat=jerr)                                         2890
      if(jerr.ne.0) return                                                 2891
      vq=max(0.0,vp)                                                       2891
      vq=vq*ni/sum(vq)                                                     2892
      call multelnetn(parm,no,ni,nr,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,th   2894 
     *r,isd,jsd,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      deallocate(vq)                                                       2895
      return                                                               2896
      end                                                                  2897
      subroutine multelnetn (parm,no,ni,nr,x,y,w,jd,vp,ne,nx,nlam,flmin,   2899 
     *ulam,thr,  isd,jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no,nr),w(no),ulam(nlam)                       2900
      real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                  2901
      integer jd(*),ia(nx),nin(nlam)                                       2902
      real, dimension (:), allocatable :: xm,xs,xv,ym,ys                        
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         2907
      allocate(xs(1:ni),stat=ierr)                                         2907
      jerr=jerr+ierr                                                       2908
      allocate(ym(1:nr),stat=ierr)                                         2908
      jerr=jerr+ierr                                                       2909
      allocate(ys(1:nr),stat=ierr)                                         2909
      jerr=jerr+ierr                                                       2910
      allocate(ju(1:ni),stat=ierr)                                         2910
      jerr=jerr+ierr                                                       2911
      allocate(xv(1:ni),stat=ierr)                                         2911
      jerr=jerr+ierr                                                       2912
      if(jerr.ne.0) return                                                 2913
      call chkvars(no,ni,x,ju)                                             2914
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2915
      if(maxval(ju) .gt. 0)goto 17951                                      2915
      jerr=7777                                                            2915
      return                                                               2915
17951 continue                                                             2916
      call multstandard1(no,ni,nr,x,y,w,isd,jsd,ju,xm,xs,ym,ys,xv,ys0,je   2917 
     *rr)
      if(jerr.ne.0) return                                                 2918
      call multelnet2(parm,ni,nr,ju,vp,y,no,ne,nx,x,nlam,flmin,ulam,thr,   2920 
     *maxit,xv,  ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 2921
17960 do 17961 k=1,lmu                                                     2921
      nk=nin(k)                                                            2922
17970 do 17971 j=1,nr                                                      2923
17980 do 17981 l=1,nk                                                      2923
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  2923
17981 continue                                                             2924
17982 continue                                                             2924
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 2925
17971 continue                                                             2926
17972 continue                                                             2926
17961 continue                                                             2927
17962 continue                                                             2927
      deallocate(xm,xs,ym,ys,ju,xv)                                        2928
      return                                                               2929
      end                                                                  2930
      subroutine multstandard1 (no,ni,nr,x,y,w,isd,jsd,ju,xm,xs,ym,ys,xv   2931 
     *,ys0,jerr)
      real x(no,ni),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),ys(nr)      2932
      integer ju(ni)                                                       2933
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                          2936
      if(jerr.ne.0) return                                                 2937
      w=w/sum(w)                                                           2937
      v=sqrt(w)                                                            2938
17990 do 17991 j=1,ni                                                      2938
      if(ju(j).eq.0)goto 17991                                             2939
      xm(j)=dot_product(w,x(:,j))                                          2939
      x(:,j)=v*(x(:,j)-xm(j))                                              2940
      xv(j)=dot_product(x(:,j),x(:,j))                                     2940
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       2941
17991 continue                                                             2942
17992 continue                                                             2942
      if(isd .ne. 0)goto 18011                                             2942
      xs=1.0                                                               2942
      goto 18021                                                           2943
18011 continue                                                             2943
18030 do 18031 j=1,ni                                                      2943
      if(ju(j).eq.0)goto 18031                                             2943
      x(:,j)=x(:,j)/xs(j)                                                  2943
18031 continue                                                             2944
18032 continue                                                             2944
      xv=1.0                                                               2945
18021 continue                                                             2946
18001 continue                                                             2946
      ys0=0.0                                                              2947
18040 do 18041 j=1,nr                                                      2948
      ym(j)=dot_product(w,y(:,j))                                          2948
      y(:,j)=v*(y(:,j)-ym(j))                                              2949
      z=dot_product(y(:,j),y(:,j))                                         2950
      if(jsd .le. 0)goto 18061                                             2950
      ys(j)=sqrt(z)                                                        2950
      y(:,j)=y(:,j)/ys(j)                                                  2950
      goto 18071                                                           2951
18061 continue                                                             2951
      ys0=ys0+z                                                            2951
18071 continue                                                             2952
18051 continue                                                             2952
18041 continue                                                             2953
18042 continue                                                             2953
      if(jsd .ne. 0)goto 18091                                             2953
      ys=1.0                                                               2953
      goto 18101                                                           2953
18091 continue                                                             2953
      ys0=nr                                                               2953
18101 continue                                                             2954
18081 continue                                                             2954
      deallocate(v)                                                        2955
      return                                                               2956
      end                                                                  2957
      subroutine multelnet2(beta,ni,nr,ju,vp,y,no,ne,nx,x,nlam,flmin,ula   2959 
     *m,thri,  maxit,xv,ys0,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99   2960 
     *9)
      real vp(ni),y(no,nr),x(no,ni),ulam(nlam),ao(nx,nr,nlam)              2961
      real rsqo(nlam),almo(nlam),xv(ni)                                    2962
      integer ju(ni),ia(nx),kin(nlam)                                      2963
      real, dimension (:), allocatable :: g,gk,del,gj                           
      integer, dimension (:), allocatable :: mm,ix                              
      real, dimension (:,:), allocatable :: a                                   
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      allocate(gj(1:nr),stat=ierr)                                         2969
      jerr=jerr+ierr                                                       2970
      allocate(gk(1:nr),stat=ierr)                                         2970
      jerr=jerr+ierr                                                       2971
      allocate(del(1:nr),stat=ierr)                                        2971
      jerr=jerr+ierr                                                       2972
      allocate(mm(1:ni),stat=ierr)                                         2972
      jerr=jerr+ierr                                                       2973
      allocate(g(1:ni),stat=ierr)                                          2973
      jerr=jerr+ierr                                                       2974
      allocate(ix(1:ni),stat=ierr)                                         2974
      jerr=jerr+ierr                                                       2975
      if(jerr.ne.0) return                                                 2976
      bta=beta                                                             2976
      omb=1.0-bta                                                          2976
      ix=0                                                                 2976
      thr=thri*ys0/nr                                                      2977
      if(flmin .ge. 1.0)goto 18121                                         2977
      eqs=max(eps,flmin)                                                   2977
      alf=eqs**(1.0/(nlam-1))                                              2977
18121 continue                                                             2978
      rsq=ys0                                                              2978
      a=0.0                                                                2978
      mm=0                                                                 2978
      nlp=0                                                                2978
      nin=nlp                                                              2978
      iz=0                                                                 2978
      mnl=min(mnlam,nlam)                                                  2978
      alm=0.0                                                              2979
18130 do 18131 j=1,ni                                                      2979
      if(ju(j).eq.0)goto 18131                                             2979
      g(j)=0.0                                                             2980
18140 do 18141 k=1,nr                                                      2980
      g(j)=g(j)+dot_product(y(:,k),x(:,j))**2                              2980
18141 continue                                                             2981
18142 continue                                                             2981
      g(j)=sqrt(g(j))                                                      2982
18131 continue                                                             2983
18132 continue                                                             2983
18150 do 18151 m=1,nlam                                                    2983
      alm0=alm                                                             2984
      if(flmin .lt. 1.0)goto 18171                                         2984
      alm=ulam(m)                                                          2984
      goto 18161                                                           2985
18171 if(m .le. 2)goto 18181                                               2985
      alm=alm*alf                                                          2985
      goto 18161                                                           2986
18181 if(m .ne. 1)goto 18191                                               2986
      alm=big                                                              2986
      goto 18201                                                           2987
18191 continue                                                             2987
      alm0=0.0                                                             2988
18210 do 18211 j=1,ni                                                      2988
      if(ju(j).eq.0)goto 18211                                             2989
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           2990
18211 continue                                                             2991
18212 continue                                                             2991
      alm0=alm0/max(bta,1.0e-3)                                            2991
      alm=alf*alm0                                                         2992
18201 continue                                                             2993
18161 continue                                                             2993
      dem=alm*omb                                                          2993
      ab=alm*bta                                                           2993
      rsq0=rsq                                                             2993
      jz=1                                                                 2994
      tlam=bta*(2.0*alm-alm0)                                              2995
18220 do 18221 k=1,ni                                                      2995
      if(ix(k).eq.1)goto 18221                                             2995
      if(ju(k).eq.0)goto 18221                                             2996
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       2997
18221 continue                                                             2998
18222 continue                                                             2998
18230 continue                                                             2998
18231 continue                                                             2998
      if(iz*jz.ne.0) go to 10260                                           2999
10680 continue                                                             2999
      nlp=nlp+1                                                            2999
      dlx=0.0                                                              3000
18240 do 18241 k=1,ni                                                      3000
      if(ix(k).eq.0)goto 18241                                             3000
      gkn=0.0                                                              3001
18250 do 18251 j=1,nr                                                      3001
      gj(j)=dot_product(y(:,j),x(:,k))                                     3002
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3002
      gkn=gkn+gk(j)**2                                                     3004
18251 continue                                                             3004
18252 continue                                                             3004
      gkn=sqrt(gkn)                                                        3004
      u=1.0-ab*vp(k)/gkn                                                   3004
      del=a(:,k)                                                           3005
      if(u .gt. 0.0)goto 18271                                             3005
      a(:,k)=0.0                                                           3005
      goto 18281                                                           3006
18271 continue                                                             3006
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3006
18281 continue                                                             3007
18261 continue                                                             3007
      del=a(:,k)-del                                                       3007
      if(maxval(abs(del)).le.0.0)goto 18241                                3008
18290 do 18291 j=1,nr                                                      3008
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3009
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3009
      dlx=max(dlx,xv(k)*del(j)**2)                                         3010
18291 continue                                                             3011
18292 continue                                                             3011
      if(mm(k) .ne. 0)goto 18311                                           3011
      nin=nin+1                                                            3011
      if(nin.gt.nx)goto 18242                                              3012
      mm(k)=nin                                                            3012
      ia(nin)=k                                                            3013
18311 continue                                                             3014
18241 continue                                                             3015
18242 continue                                                             3015
      if(nin.gt.nx)goto 18232                                              3016
      if(dlx .ge. thr)goto 18331                                           3016
      ixx=0                                                                3017
18340 do 18341 k=1,ni                                                      3017
      if(ix(k).eq.1)goto 18341                                             3017
      if(ju(k).eq.0)goto 18341                                             3017
      g(k)=0.0                                                             3018
18350 do 18351 j=1,nr                                                      3018
      g(k)=g(k)+dot_product(y(:,j),x(:,k))**2                              3018
18351 continue                                                             3019
18352 continue                                                             3019
      g(k)=sqrt(g(k))                                                      3020
      if(g(k) .le. ab*vp(k))goto 18371                                     3020
      ix(k)=1                                                              3020
      ixx=1                                                                3020
18371 continue                                                             3021
18341 continue                                                             3022
18342 continue                                                             3022
      if(ixx.eq.1) go to 10680                                             3023
      goto 18232                                                           3024
18331 continue                                                             3025
      if(nlp .le. maxit)goto 18391                                         3025
      jerr=-m                                                              3025
      return                                                               3025
18391 continue                                                             3026
10260 continue                                                             3026
      iz=1                                                                 3027
18400 continue                                                             3027
18401 continue                                                             3027
      nlp=nlp+1                                                            3027
      dlx=0.0                                                              3028
18410 do 18411 l=1,nin                                                     3028
      k=ia(l)                                                              3028
      gkn=0.0                                                              3029
18420 do 18421 j=1,nr                                                      3029
      gj(j)=dot_product(y(:,j),x(:,k))                                     3030
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3030
      gkn=gkn+gk(j)**2                                                     3032
18421 continue                                                             3032
18422 continue                                                             3032
      gkn=sqrt(gkn)                                                        3032
      u=1.0-ab*vp(k)/gkn                                                   3032
      del=a(:,k)                                                           3033
      if(u .gt. 0.0)goto 18441                                             3033
      a(:,k)=0.0                                                           3033
      goto 18451                                                           3034
18441 continue                                                             3034
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3034
18451 continue                                                             3035
18431 continue                                                             3035
      del=a(:,k)-del                                                       3035
      if(maxval(abs(del)).le.0.0)goto 18411                                3036
18460 do 18461 j=1,nr                                                      3036
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3037
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3037
      dlx=max(dlx,xv(k)*del(j)**2)                                         3038
18461 continue                                                             3039
18462 continue                                                             3039
18411 continue                                                             3040
18412 continue                                                             3040
      if(dlx.lt.thr)goto 18402                                             3040
      if(nlp .le. maxit)goto 18481                                         3040
      jerr=-m                                                              3040
      return                                                               3040
18481 continue                                                             3041
      goto 18401                                                           3042
18402 continue                                                             3042
      jz=0                                                                 3043
      goto 18231                                                           3044
18232 continue                                                             3044
      if(nin .le. nx)goto 18501                                            3044
      jerr=-10000-m                                                        3044
      goto 18152                                                           3044
18501 continue                                                             3045
      if(nin .le. 0)goto 18521                                             3045
18530 do 18531 j=1,nr                                                      3045
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3045
18531 continue                                                             3045
18532 continue                                                             3045
18521 continue                                                             3046
      kin(m)=nin                                                           3047
      rsqo(m)=1.0-rsq/ys0                                                  3047
      almo(m)=alm                                                          3047
      lmu=m                                                                3048
      if(m.lt.mnl)goto 18151                                               3048
      if(flmin.ge.1.0)goto 18151                                           3049
      me=0                                                                 3049
18540 do 18541 j=1,nin                                                     3049
      if(ao(j,1,m).ne.0.0) me=me+1                                         3049
18541 continue                                                             3049
18542 continue                                                             3049
      if(me.gt.ne)goto 18152                                               3050
      if(rsq0-rsq.lt.sml*rsq)goto 18152                                    3050
      if(rsqo(m).gt.rsqmax)goto 18152                                      3051
18151 continue                                                             3052
18152 continue                                                             3052
      deallocate(a,mm,g,ix,del,gj,gk)                                      3053
      return                                                               3054
      end                                                                  3055
      subroutine multsolns(ni,nx,nr,lmu,a,ia,nin,b)                        3056
      real a(nx,nr,lmu),b(ni,nr,lmu)                                       3056
      integer ia(nx),nin(lmu)                                              3057
18550 do 18551 lam=1,lmu                                                   3057
      call multuncomp(ni,nr,nx,a(1,1,lam),ia,nin(lam),b(1,1,lam))          3057
18551 continue                                                             3058
18552 continue                                                             3058
      return                                                               3059
      end                                                                  3060
      subroutine multuncomp(ni,nr,nx,ca,ia,nin,a)                          3061
      real ca(nx,nr),a(ni,nr)                                              3061
      integer ia(nx)                                                       3062
      a=0.0                                                                3063
      if(nin .le. 0)goto 18571                                             3063
18580 do 18581 j=1,nr                                                      3063
      a(ia(1:nin),j)=ca(1:nin,j)                                           3063
18581 continue                                                             3063
18582 continue                                                             3063
18571 continue                                                             3064
      return                                                               3065
      end                                                                  3066
      subroutine multmodval(nx,nr,a0,ca,ia,nin,n,x,f)                      3067
      real a0(nr),ca(nx,nr),x(n,*),f(nr,n)                                 3067
      integer ia(nx)                                                       3068
18590 do 18591 i=1,n                                                       3068
      f(:,i)=a0                                                            3068
18591 continue                                                             3068
18592 continue                                                             3068
      if(nin.le.0) return                                                  3069
18600 do 18601 i=1,n                                                       3069
18610 do 18611 j=1,nr                                                      3069
      f(j,i)=f(j,i)+dot_product(ca(1:nin,j),x(i,ia(1:nin)))                3069
18611 continue                                                             3069
18612 continue                                                             3069
18601 continue                                                             3070
18602 continue                                                             3070
      return                                                               3071
      end                                                                  3072
      subroutine multspelnet  (parm,no,ni,nr,x,ix,jx,y,w,jd,vp,ne,nx,nla   3075 
     *m,flmin,ulam,thr,isd,jsd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      real x(*),y(no,nr),w(no),vp(ni),ulam(nlam)                           3076
      real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                  3077
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3078
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 18631                                    3081
      jerr=10000                                                           3081
      return                                                               3081
18631 continue                                                             3082
      allocate(vq(1:ni),stat=jerr)                                         3082
      if(jerr.ne.0) return                                                 3083
      vq=max(0.0,vp)                                                       3083
      vq=vq*ni/sum(vq)                                                     3084
      call multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin   3086 
     *,  ulam,thr,isd,jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      deallocate(vq)                                                       3087
      return                                                               3088
      end                                                                  3089
      subroutine multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,ne,nx,nlam   3091 
     *,flmin,ulam,  thr,isd,jsd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),vp(ni),y(no,nr),w(no),ulam(nlam)                           3092
      real ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)                  3093
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3094
      real, dimension (:), allocatable :: xm,xs,xv,ym,ys                        
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         3099
      allocate(xs(1:ni),stat=ierr)                                         3099
      jerr=jerr+ierr                                                       3100
      allocate(ym(1:nr),stat=ierr)                                         3100
      jerr=jerr+ierr                                                       3101
      allocate(ys(1:nr),stat=ierr)                                         3101
      jerr=jerr+ierr                                                       3102
      allocate(ju(1:ni),stat=ierr)                                         3102
      jerr=jerr+ierr                                                       3103
      allocate(xv(1:ni),stat=ierr)                                         3103
      jerr=jerr+ierr                                                       3104
      if(jerr.ne.0) return                                                 3105
      call spchkvars(no,ni,x,ix,ju)                                        3106
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3107
      if(maxval(ju) .gt. 0)goto 18651                                      3107
      jerr=7777                                                            3107
      return                                                               3107
18651 continue                                                             3108
      call multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,xm,xs,ym,ys,x   3109 
     *v,ys0,jerr)
      if(jerr.ne.0) return                                                 3110
      call multspelnet2(parm,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin   3112 
     *,  ulam,thr,maxit,xm,xs,xv,ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 3113
18660 do 18661 k=1,lmu                                                     3113
      nk=nin(k)                                                            3114
18670 do 18671 j=1,nr                                                      3115
18680 do 18681 l=1,nk                                                      3115
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  3115
18681 continue                                                             3116
18682 continue                                                             3116
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 3117
18671 continue                                                             3118
18672 continue                                                             3118
18661 continue                                                             3119
18662 continue                                                             3119
      deallocate(xm,xs,ym,ys,ju,xv)                                        3120
      return                                                               3121
      end                                                                  3122
      subroutine multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,  xm,xs   3124 
     *,ym,ys,xv,ys0,jerr)
      real x(*),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),ys(nr)          3125
      integer ix(*),jx(*),ju(ni)                                           3126
      w=w/sum(w)                                                           3127
18690 do 18691 j=1,ni                                                      3127
      if(ju(j).eq.0)goto 18691                                             3128
      jb=ix(j)                                                             3128
      je=ix(j+1)-1                                                         3128
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3129
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 3130
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       3131
18691 continue                                                             3132
18692 continue                                                             3132
      if(isd .ne. 0)goto 18711                                             3132
      xs=1.0                                                               3132
      goto 18721                                                           3132
18711 continue                                                             3132
      xv=1.0                                                               3132
18721 continue                                                             3133
18701 continue                                                             3133
      ys0=0.0                                                              3134
18730 do 18731 j=1,nr                                                      3135
      ym(j)=dot_product(w,y(:,j))                                          3135
      y(:,j)=y(:,j)-ym(j)                                                  3136
      z=dot_product(w,y(:,j)**2)                                           3137
      if(jsd .le. 0)goto 18751                                             3137
      ys(j)=sqrt(z)                                                        3137
      y(:,j)=y(:,j)/ys(j)                                                  3137
      goto 18761                                                           3138
18751 continue                                                             3138
      ys0=ys0+z                                                            3138
18761 continue                                                             3139
18741 continue                                                             3139
18731 continue                                                             3140
18732 continue                                                             3140
      if(jsd .ne. 0)goto 18781                                             3140
      ys=1.0                                                               3140
      goto 18791                                                           3140
18781 continue                                                             3140
      ys0=nr                                                               3140
18791 continue                                                             3141
18771 continue                                                             3141
      return                                                               3142
      end                                                                  3143
      subroutine multspelnet2(beta,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam   3145 
     *,flmin,ulam,  thri,maxit,xm,xs,xv,ys0,lmu,ao,ia,kin,rsqo,almo,nlp,
     *jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99   3146 
     *9)
      real y(no,nr),w(no),x(*),vp(ni),ulam(nlam)                           3147
      real ao(nx,nr,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)       3148
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          3149
      real, dimension (:), allocatable :: g,gj,gk,del,o                         
      integer, dimension (:), allocatable :: mm,iy                              
      real, dimension (:,:), allocatable :: a                                   
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      allocate(mm(1:ni),stat=ierr)                                         3155
      jerr=jerr+ierr                                                       3156
      allocate(g(1:ni),stat=ierr)                                          3156
      jerr=jerr+ierr                                                       3157
      allocate(gj(1:nr),stat=ierr)                                         3157
      jerr=jerr+ierr                                                       3158
      allocate(gk(1:nr),stat=ierr)                                         3158
      jerr=jerr+ierr                                                       3159
      allocate(del(1:nr),stat=ierr)                                        3159
      jerr=jerr+ierr                                                       3160
      allocate(o(1:nr),stat=ierr)                                          3160
      jerr=jerr+ierr                                                       3161
      allocate(iy(1:ni),stat=ierr)                                         3161
      jerr=jerr+ierr                                                       3162
      if(jerr.ne.0) return                                                 3163
      bta=beta                                                             3163
      omb=1.0-bta                                                          3163
      alm=0.0                                                              3163
      iy=0                                                                 3163
      thr=thri*ys0/nr                                                      3164
      if(flmin .ge. 1.0)goto 18811                                         3164
      eqs=max(eps,flmin)                                                   3164
      alf=eqs**(1.0/(nlam-1))                                              3164
18811 continue                                                             3165
      rsq=ys0                                                              3165
      a=0.0                                                                3165
      mm=0                                                                 3165
      o=0.0                                                                3165
      nlp=0                                                                3165
      nin=nlp                                                              3165
      iz=0                                                                 3165
      mnl=min(mnlam,nlam)                                                  3166
18820 do 18821 j=1,ni                                                      3166
      if(ju(j).eq.0)goto 18821                                             3166
      jb=ix(j)                                                             3166
      je=ix(j+1)-1                                                         3166
      g(j)=0.0                                                             3167
18830 do 18831 k=1,nr                                                      3168
      g(j)=g(j)+(dot_product(y(jx(jb:je),k),w(jx(jb:je))*x(jb:je))/xs(j)   3169 
     *)**2
18831 continue                                                             3170
18832 continue                                                             3170
      g(j)=sqrt(g(j))                                                      3171
18821 continue                                                             3172
18822 continue                                                             3172
18840 do 18841 m=1,nlam                                                    3172
      alm0=alm                                                             3173
      if(flmin .lt. 1.0)goto 18861                                         3173
      alm=ulam(m)                                                          3173
      goto 18851                                                           3174
18861 if(m .le. 2)goto 18871                                               3174
      alm=alm*alf                                                          3174
      goto 18851                                                           3175
18871 if(m .ne. 1)goto 18881                                               3175
      alm=big                                                              3175
      goto 18891                                                           3176
18881 continue                                                             3176
      alm0=0.0                                                             3177
18900 do 18901 j=1,ni                                                      3177
      if(ju(j).eq.0)goto 18901                                             3178
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           3179
18901 continue                                                             3180
18902 continue                                                             3180
      alm0=alm0/max(bta,1.0e-3)                                            3180
      alm=alf*alm0                                                         3181
18891 continue                                                             3182
18851 continue                                                             3182
      dem=alm*omb                                                          3182
      ab=alm*bta                                                           3182
      rsq0=rsq                                                             3182
      jz=1                                                                 3183
      tlam=bta*(2.0*alm-alm0)                                              3184
18910 do 18911 k=1,ni                                                      3184
      if(iy(k).eq.1)goto 18911                                             3184
      if(ju(k).eq.0)goto 18911                                             3185
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       3186
18911 continue                                                             3187
18912 continue                                                             3187
18920 continue                                                             3187
18921 continue                                                             3187
      if(iz*jz.ne.0) go to 10260                                           3188
10680 continue                                                             3188
      nlp=nlp+1                                                            3188
      dlx=0.0                                                              3189
18930 do 18931 k=1,ni                                                      3189
      if(iy(k).eq.0)goto 18931                                             3189
      jb=ix(k)                                                             3189
      je=ix(k+1)-1                                                         3189
      gkn=0.0                                                              3190
18940 do 18941 j=1,nr                                                      3191
      gj(j)=dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(k)   3192
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3192
      gkn=gkn+gk(j)**2                                                     3193
18941 continue                                                             3194
18942 continue                                                             3194
      gkn=sqrt(gkn)                                                        3194
      u=1.0-ab*vp(k)/gkn                                                   3194
      del=a(:,k)                                                           3195
      if(u .gt. 0.0)goto 18961                                             3195
      a(:,k)=0.0                                                           3195
      goto 18971                                                           3196
18961 continue                                                             3196
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3196
18971 continue                                                             3197
18951 continue                                                             3197
      del=a(:,k)-del                                                       3197
      if(maxval(abs(del)).le.0.0)goto 18931                                3198
      if(mm(k) .ne. 0)goto 18991                                           3198
      nin=nin+1                                                            3198
      if(nin.gt.nx)goto 18932                                              3199
      mm(k)=nin                                                            3199
      ia(nin)=k                                                            3200
18991 continue                                                             3201
19000 do 19001 j=1,nr                                                      3201
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3202
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3203
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3203
      dlx=max(xv(k)*del(j)**2,dlx)                                         3204
19001 continue                                                             3205
19002 continue                                                             3205
18931 continue                                                             3206
18932 continue                                                             3206
      if(nin.gt.nx)goto 18922                                              3207
      if(dlx .ge. thr)goto 19021                                           3207
      ixx=0                                                                3208
19030 do 19031 j=1,ni                                                      3208
      if(iy(j).eq.1)goto 19031                                             3208
      if(ju(j).eq.0)goto 19031                                             3209
      jb=ix(j)                                                             3209
      je=ix(j+1)-1                                                         3209
      g(j)=0.0                                                             3210
19040 do 19041 k=1,nr                                                      3210
      g(j)=g(j)+  (dot_product(y(jx(jb:je),k)+o(k),w(jx(jb:je))*x(jb:je)   3212 
     *)/xs(j))**2
19041 continue                                                             3213
19042 continue                                                             3213
      g(j)=sqrt(g(j))                                                      3214
      if(g(j) .le. ab*vp(j))goto 19061                                     3214
      iy(j)=1                                                              3214
      ixx=1                                                                3214
19061 continue                                                             3215
19031 continue                                                             3216
19032 continue                                                             3216
      if(ixx.eq.1) go to 10680                                             3217
      goto 18922                                                           3218
19021 continue                                                             3219
      if(nlp .le. maxit)goto 19081                                         3219
      jerr=-m                                                              3219
      return                                                               3219
19081 continue                                                             3220
10260 continue                                                             3220
      iz=1                                                                 3221
19090 continue                                                             3221
19091 continue                                                             3221
      nlp=nlp+1                                                            3221
      dlx=0.0                                                              3222
19100 do 19101 l=1,nin                                                     3222
      k=ia(l)                                                              3222
      jb=ix(k)                                                             3222
      je=ix(k+1)-1                                                         3222
      gkn=0.0                                                              3223
19110 do 19111 j=1,nr                                                      3223
      gj(j)=  dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(   3225 
     *k)
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3225
      gkn=gkn+gk(j)**2                                                     3226
19111 continue                                                             3227
19112 continue                                                             3227
      gkn=sqrt(gkn)                                                        3227
      u=1.0-ab*vp(k)/gkn                                                   3227
      del=a(:,k)                                                           3228
      if(u .gt. 0.0)goto 19131                                             3228
      a(:,k)=0.0                                                           3228
      goto 19141                                                           3229
19131 continue                                                             3229
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3229
19141 continue                                                             3230
19121 continue                                                             3230
      del=a(:,k)-del                                                       3230
      if(maxval(abs(del)).le.0.0)goto 19101                                3231
19150 do 19151 j=1,nr                                                      3231
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3232
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3233
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3233
      dlx=max(xv(k)*del(j)**2,dlx)                                         3234
19151 continue                                                             3235
19152 continue                                                             3235
19101 continue                                                             3236
19102 continue                                                             3236
      if(dlx.lt.thr)goto 19092                                             3236
      if(nlp .le. maxit)goto 19171                                         3236
      jerr=-m                                                              3236
      return                                                               3236
19171 continue                                                             3237
      goto 19091                                                           3238
19092 continue                                                             3238
      jz=0                                                                 3239
      goto 18921                                                           3240
18922 continue                                                             3240
      if(nin .le. nx)goto 19191                                            3240
      jerr=-10000-m                                                        3240
      goto 18842                                                           3240
19191 continue                                                             3241
      if(nin .le. 0)goto 19211                                             3241
19220 do 19221 j=1,nr                                                      3241
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3241
19221 continue                                                             3241
19222 continue                                                             3241
19211 continue                                                             3242
      kin(m)=nin                                                           3243
      rsqo(m)=1.0-rsq/ys0                                                  3243
      almo(m)=alm                                                          3243
      lmu=m                                                                3244
      if(m.lt.mnl)goto 18841                                               3244
      if(flmin.ge.1.0)goto 18841                                           3245
      me=0                                                                 3245
19230 do 19231 j=1,nin                                                     3245
      if(ao(j,1,m).ne.0.0) me=me+1                                         3245
19231 continue                                                             3245
19232 continue                                                             3245
      if(me.gt.ne)goto 18842                                               3246
      if(rsq0-rsq.lt.sml*rsq)goto 18842                                    3246
      if(rsqo(m).gt.rsqmax)goto 18842                                      3247
18841 continue                                                             3248
18842 continue                                                             3248
      deallocate(a,mm,g,iy,gj,gk,del,o)                                    3249
      return                                                               3250
      end                                                                  3251
      subroutine multlognetn(parm,no,ni,nc,x,y,g,w,ju,vp,ne,nx,nlam,flmi   3253 
     *n,ulam,shri,  isd,maxit,xv,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   3255 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam)              3256
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),xv(ni)            3257
      integer ju(ni),m(nx),kin(nlam)                                       3258
      real, dimension (:,:), allocatable :: q,r,b,bs                            
      real, dimension (:), allocatable :: sxp,sxpl,ga,gk,del                    
      integer, dimension (:), allocatable :: mm,is,ixx                          
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no,1:nc),stat=ierr); jerr=jerr+ierr;                         
      allocate(mm(1:ni),stat=ierr)                                         3267
      jerr=jerr+ierr                                                       3268
      allocate(is(1:max(nc,ni)),stat=ierr)                                 3268
      jerr=jerr+ierr                                                       3269
      allocate(sxp(1:no),stat=ierr)                                        3269
      jerr=jerr+ierr                                                       3270
      allocate(sxpl(1:no),stat=ierr)                                       3270
      jerr=jerr+ierr                                                       3271
      allocate(ga(1:ni),stat=ierr)                                         3271
      jerr=jerr+ierr                                                       3272
      allocate(ixx(1:ni),stat=ierr)                                        3272
      jerr=jerr+ierr                                                       3273
      allocate(gk(1:nc),stat=ierr)                                         3273
      jerr=jerr+ierr                                                       3274
      allocate(del(1:nc),stat=ierr)                                        3274
      jerr=jerr+ierr                                                       3275
      if(jerr.ne.0) return                                                 3276
      pmax=1.0-pmin                                                        3276
      emin=pmin/pmax                                                       3276
      emax=1.0/emin                                                        3277
      bta=parm                                                             3277
      omb=1.0-bta                                                          3277
      dev1=0.0                                                             3277
      dev0=0.0                                                             3278
19240 do 19241 ic=1,nc                                                     3278
      q0=dot_product(w,y(:,ic))                                            3279
      if(q0 .gt. pmin)goto 19261                                           3279
      jerr =8000+ic                                                        3279
      return                                                               3279
19261 continue                                                             3280
      if(q0 .lt. pmax)goto 19281                                           3280
      jerr =9000+ic                                                        3280
      return                                                               3280
19281 continue                                                             3281
      b(0,ic)=log(q0)                                                      3281
      dev1=dev1-q0*b(0,ic)                                                 3281
      b(1:ni,ic)=0.0                                                       3282
19241 continue                                                             3283
19242 continue                                                             3283
      ixx=0                                                                3283
      al=0.0                                                               3284
      if(nonzero(no*nc,g) .ne. 0)goto 19301                                3285
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         3285
      sxp=0.0                                                              3286
19310 do 19311 ic=1,nc                                                     3286
      q(:,ic)=exp(b(0,ic))                                                 3286
      sxp=sxp+q(:,ic)                                                      3286
19311 continue                                                             3287
19312 continue                                                             3287
      goto 19321                                                           3288
19301 continue                                                             3288
19330 do 19331 i=1,no                                                      3288
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3288
19331 continue                                                             3288
19332 continue                                                             3288
      sxp=0.0                                                              3289
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 3289
      if(jerr.ne.0) return                                                 3290
      dev1=0.0                                                             3291
19340 do 19341 ic=1,nc                                                     3291
      q(:,ic)=b(0,ic)+g(:,ic)                                              3292
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             3293
      q(:,ic)=exp(q(:,ic))                                                 3293
      sxp=sxp+q(:,ic)                                                      3294
19341 continue                                                             3295
19342 continue                                                             3295
      sxpl=w*log(sxp)                                                      3295
19350 do 19351 ic=1,nc                                                     3295
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  3295
19351 continue                                                             3296
19352 continue                                                             3296
19321 continue                                                             3297
19291 continue                                                             3297
19360 do 19361 ic=1,nc                                                     3297
19370 do 19371 i=1,no                                                      3297
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               3297
19371 continue                                                             3297
19372 continue                                                             3297
19361 continue                                                             3298
19362 continue                                                             3298
      dev0=dev0+dev1                                                       3299
      if(flmin .ge. 1.0)goto 19391                                         3299
      eqs=max(eps,flmin)                                                   3299
      alf=eqs**(1.0/(nlam-1))                                              3299
19391 continue                                                             3300
      m=0                                                                  3300
      mm=0                                                                 3300
      nin=0                                                                3300
      nlp=0                                                                3300
      mnl=min(mnlam,nlam)                                                  3300
      bs=0.0                                                               3300
      shr=shri*dev0                                                        3301
      ga=0.0                                                               3302
19400 do 19401 ic=1,nc                                                     3302
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3303
19410 do 19411 j=1,ni                                                      3303
      if(ju(j).ne.0) ga(j)=ga(j)+dot_product(r(:,ic),x(:,j))**2            3303
19411 continue                                                             3304
19412 continue                                                             3304
19401 continue                                                             3305
19402 continue                                                             3305
      ga=sqrt(ga)                                                          3306
19420 do 19421 ilm=1,nlam                                                  3306
      al0=al                                                               3307
      if(flmin .lt. 1.0)goto 19441                                         3307
      al=ulam(ilm)                                                         3307
      goto 19431                                                           3308
19441 if(ilm .le. 2)goto 19451                                             3308
      al=al*alf                                                            3308
      goto 19431                                                           3309
19451 if(ilm .ne. 1)goto 19461                                             3309
      al=big                                                               3309
      goto 19471                                                           3310
19461 continue                                                             3310
      al0=0.0                                                              3311
19480 do 19481 j=1,ni                                                      3311
      if(ju(j).eq.0)goto 19481                                             3311
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3311
19481 continue                                                             3312
19482 continue                                                             3312
      al0=al0/max(bta,1.0e-3)                                              3312
      al=alf*al0                                                           3313
19471 continue                                                             3314
19431 continue                                                             3314
      al2=al*omb                                                           3314
      al1=al*bta                                                           3314
      tlam=bta*(2.0*al-al0)                                                3315
19490 do 19491 k=1,ni                                                      3315
      if(ixx(k).eq.1)goto 19491                                            3315
      if(ju(k).eq.0)goto 19491                                             3316
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     3317
19491 continue                                                             3318
19492 continue                                                             3318
10680 continue                                                             3319
19500 continue                                                             3319
19501 continue                                                             3319
      ix=0                                                                 3319
      jx=ix                                                                3320
19510 do 19511 ic=1,nc                                                     3321
      bs(0,ic)=b(0,ic)                                                     3321
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          3322
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3322
      d=sum(r(:,ic))                                                       3323
      b(0,ic)=b(0,ic)+d                                                    3323
      r(:,ic)=r(:,ic)-d*w                                                  3323
      dlx=max(dlx,d**2)                                                    3325
19511 continue                                                             3326
19512 continue                                                             3326
19520 continue                                                             3326
19521 continue                                                             3326
      nlp=nlp+nc                                                           3326
      dlx=0.0                                                              3327
19530 do 19531 k=1,ni                                                      3327
      if(ixx(k).eq.0)goto 19531                                            3327
      gkn=0.0                                                              3328
19540 do 19541 ic=1,nc                                                     3328
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     3329
      gkn=gkn+gk(ic)**2                                                    3330
19541 continue                                                             3331
19542 continue                                                             3331
      gkn=sqrt(gkn)                                                        3331
      u=1.0-al1*vp(k)/gkn                                                  3331
      del=b(k,:)                                                           3332
      if(u .gt. 0.0)goto 19561                                             3332
      b(k,:)=0.0                                                           3332
      goto 19571                                                           3333
19561 continue                                                             3333
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2))                                      3333
19571 continue                                                             3334
19551 continue                                                             3334
      del=b(k,:)-del                                                       3334
      if(maxval(abs(del)).le.0.0)goto 19531                                3335
19580 do 19581 ic=1,nc                                                     3335
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3336
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     3337
19581 continue                                                             3338
19582 continue                                                             3338
      if(mm(k) .ne. 0)goto 19601                                           3338
      nin=nin+1                                                            3339
      if(nin .le. nx)goto 19621                                            3339
      jx=1                                                                 3339
      goto 19532                                                           3339
19621 continue                                                             3340
      mm(k)=nin                                                            3340
      m(nin)=k                                                             3341
19601 continue                                                             3342
19531 continue                                                             3343
19532 continue                                                             3343
      if(jx.gt.0)goto 19522                                                3343
      if(dlx.lt.shr)goto 19522                                             3344
      if(nlp .le. maxit)goto 19641                                         3344
      jerr=-ilm                                                            3344
      return                                                               3344
19641 continue                                                             3345
19650 continue                                                             3345
19651 continue                                                             3345
      nlp=nlp+nc                                                           3345
      dlx=0.0                                                              3346
19660 do 19661 l=1,nin                                                     3346
      k=m(l)                                                               3346
      gkn=0.0                                                              3347
19670 do 19671 ic=1,nc                                                     3347
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     3348
      gkn=gkn+gk(ic)**2                                                    3349
19671 continue                                                             3350
19672 continue                                                             3350
      gkn=sqrt(gkn)                                                        3350
      u=1.0-al1*vp(k)/gkn                                                  3350
      del=b(k,:)                                                           3351
      if(u .gt. 0.0)goto 19691                                             3351
      b(k,:)=0.0                                                           3351
      goto 19701                                                           3352
19691 continue                                                             3352
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2))                                      3352
19701 continue                                                             3353
19681 continue                                                             3353
      del=b(k,:)-del                                                       3353
      if(maxval(abs(del)).le.0.0)goto 19661                                3354
19710 do 19711 ic=1,nc                                                     3354
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3355
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     3356
19711 continue                                                             3357
19712 continue                                                             3357
19661 continue                                                             3358
19662 continue                                                             3358
      if(dlx.lt.shr)goto 19652                                             3358
      if(nlp .le. maxit)goto 19731                                         3358
      jerr=-ilm                                                            3358
      return                                                               3358
19731 continue                                                             3360
      goto 19651                                                           3361
19652 continue                                                             3361
      goto 19521                                                           3362
19522 continue                                                             3362
      if(jx.gt.0)goto 19502                                                3363
19740 do 19741 ic=1,nc                                                     3364
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                                3365
      if(ix .ne. 0)goto 19761                                              3366
19770 do 19771 j=1,nin                                                     3366
      k=m(j)                                                               3367
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 19791                   3367
      ix=1                                                                 3367
      goto 19772                                                           3367
19791 continue                                                             3369
19771 continue                                                             3370
19772 continue                                                             3370
19761 continue                                                             3371
19800 do 19801 i=1,no                                                      3371
      fi=b(0,ic)+g(i,ic)                                                   3373
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         3374
      fi=min(max(exmn,fi),exmx)                                            3374
      sxp(i)=sxp(i)-q(i,ic)                                                3375
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    3376
      sxp(i)=sxp(i)+q(i,ic)                                                3377
19801 continue                                                             3378
19802 continue                                                             3378
19741 continue                                                             3379
19742 continue                                                             3379
      s=-sum(b(0,:))/nc                                                    3379
      b(0,:)=b(0,:)+s                                                      3380
      if(jx.gt.0)goto 19502                                                3381
      if(ix .ne. 0)goto 19821                                              3382
19830 do 19831 k=1,ni                                                      3382
      if(ixx(k).eq.1)goto 19831                                            3382
      if(ju(k).eq.0)goto 19831                                             3382
      ga(k)=0.0                                                            3382
19831 continue                                                             3383
19832 continue                                                             3383
19840 do 19841 ic=1,nc                                                     3383
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3384
19850 do 19851 k=1,ni                                                      3384
      if(ixx(k).eq.1)goto 19851                                            3384
      if(ju(k).eq.0)goto 19851                                             3385
      ga(k)=ga(k)+dot_product(r(:,ic),x(:,k))**2                           3386
19851 continue                                                             3387
19852 continue                                                             3387
19841 continue                                                             3388
19842 continue                                                             3388
      ga=sqrt(ga)                                                          3389
19860 do 19861 k=1,ni                                                      3389
      if(ixx(k).eq.1)goto 19861                                            3389
      if(ju(k).eq.0)goto 19861                                             3390
      if(ga(k) .le. al1*vp(k))goto 19881                                   3390
      ixx(k)=1                                                             3390
      ix=1                                                                 3390
19881 continue                                                             3391
19861 continue                                                             3392
19862 continue                                                             3392
      if(ix.eq.1) go to 10680                                              3393
      goto 19502                                                           3394
19821 continue                                                             3395
      goto 19501                                                           3396
19502 continue                                                             3396
      if(jx .le. 0)goto 19901                                              3396
      jerr=-10000-ilm                                                      3396
      goto 19422                                                           3396
19901 continue                                                             3396
      devi=0.0                                                             3397
19910 do 19911 ic=1,nc                                                     3398
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          3398
      a0(ic,ilm)=b(0,ic)                                                   3399
19920 do 19921 i=1,no                                                      3399
      if(y(i,ic).le.0.0)goto 19921                                         3400
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           3401
19921 continue                                                             3402
19922 continue                                                             3402
19911 continue                                                             3403
19912 continue                                                             3403
      kin(ilm)=nin                                                         3403
      alm(ilm)=al                                                          3403
      lmu=ilm                                                              3404
      dev(ilm)=(dev1-devi)/dev0                                            3405
      if(ilm.lt.mnl)goto 19421                                             3405
      if(flmin.ge.1.0)goto 19421                                           3406
      me=0                                                                 3406
19930 do 19931 j=1,nin                                                     3406
      if(a(j,1,ilm).ne.0.0) me=me+1                                        3406
19931 continue                                                             3406
19932 continue                                                             3406
      if(me.gt.ne)goto 19422                                               3407
      if(dev(ilm).gt.devmax)goto 19422                                     3407
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 19422                             3408
19421 continue                                                             3409
19422 continue                                                             3409
      g=log(q)                                                             3409
19940 do 19941 i=1,no                                                      3409
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3409
19941 continue                                                             3410
19942 continue                                                             3410
      deallocate(sxp,b,bs,r,q,mm,is,ga,ixx,gk,del,sxpl)                    3411
      return                                                               3412
      end                                                                  3413
      subroutine multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,ne,nx,   3415 
     *nlam,flmin,  ulam,shri,isd,maxit,xv,xb,xs,lmu,a0,a,m,kin,dev0,dev,
     *alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   3417 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni),   3418 
     *xv(ni)
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   3419
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           3420
      real, dimension (:,:), allocatable :: q,r,b,bs                            
      real, dimension (:), allocatable :: sxp,sxpl,ga,gk,del,sc,svr             
      integer, dimension (:), allocatable :: mm,is,iy                           
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(mm(1:ni),stat=ierr)                                         3429
      jerr=jerr+ierr                                                       3430
      allocate(ga(1:ni),stat=ierr)                                         3430
      jerr=jerr+ierr                                                       3431
      allocate(gk(1:nc),stat=ierr)                                         3431
      jerr=jerr+ierr                                                       3432
      allocate(del(1:nc),stat=ierr)                                        3432
      jerr=jerr+ierr                                                       3433
      allocate(iy(1:ni),stat=ierr)                                         3433
      jerr=jerr+ierr                                                       3434
      allocate(is(1:max(nc,ni)),stat=ierr)                                 3434
      jerr=jerr+ierr                                                       3435
      allocate(sxp(1:no),stat=ierr)                                        3435
      jerr=jerr+ierr                                                       3436
      allocate(sxpl(1:no),stat=ierr)                                       3436
      jerr=jerr+ierr                                                       3437
      allocate(svr(1:nc),stat=ierr)                                        3437
      jerr=jerr+ierr                                                       3438
      allocate(sc(1:no),stat=ierr)                                         3438
      jerr=jerr+ierr                                                       3439
      if(jerr.ne.0) return                                                 3440
      pmax=1.0-pmin                                                        3440
      emin=pmin/pmax                                                       3440
      emax=1.0/emin                                                        3441
      bta=parm                                                             3441
      omb=1.0-bta                                                          3441
      dev1=0.0                                                             3441
      dev0=0.0                                                             3442
19950 do 19951 ic=1,nc                                                     3442
      q0=dot_product(w,y(:,ic))                                            3443
      if(q0 .gt. pmin)goto 19971                                           3443
      jerr =8000+ic                                                        3443
      return                                                               3443
19971 continue                                                             3444
      if(q0 .lt. pmax)goto 19991                                           3444
      jerr =9000+ic                                                        3444
      return                                                               3444
19991 continue                                                             3445
      b(1:ni,ic)=0.0                                                       3445
      b(0,ic)=log(q0)                                                      3445
      dev1=dev1-q0*b(0,ic)                                                 3446
19951 continue                                                             3447
19952 continue                                                             3447
      iy=0                                                                 3447
      al=0.0                                                               3448
      if(nonzero(no*nc,g) .ne. 0)goto 20011                                3449
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         3449
      sxp=0.0                                                              3450
20020 do 20021 ic=1,nc                                                     3450
      q(:,ic)=exp(b(0,ic))                                                 3450
      sxp=sxp+q(:,ic)                                                      3450
20021 continue                                                             3451
20022 continue                                                             3451
      goto 20031                                                           3452
20011 continue                                                             3452
20040 do 20041 i=1,no                                                      3452
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3452
20041 continue                                                             3452
20042 continue                                                             3452
      sxp=0.0                                                              3453
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 3453
      if(jerr.ne.0) return                                                 3454
      dev1=0.0                                                             3455
20050 do 20051 ic=1,nc                                                     3455
      q(:,ic)=b(0,ic)+g(:,ic)                                              3456
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             3457
      q(:,ic)=exp(q(:,ic))                                                 3457
      sxp=sxp+q(:,ic)                                                      3458
20051 continue                                                             3459
20052 continue                                                             3459
      sxpl=w*log(sxp)                                                      3459
20060 do 20061 ic=1,nc                                                     3459
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  3459
20061 continue                                                             3460
20062 continue                                                             3460
20031 continue                                                             3461
20001 continue                                                             3461
20070 do 20071 ic=1,nc                                                     3461
20080 do 20081 i=1,no                                                      3461
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               3461
20081 continue                                                             3461
20082 continue                                                             3461
20071 continue                                                             3462
20072 continue                                                             3462
      dev0=dev0+dev1                                                       3463
      if(flmin .ge. 1.0)goto 20101                                         3463
      eqs=max(eps,flmin)                                                   3463
      alf=eqs**(1.0/(nlam-1))                                              3463
20101 continue                                                             3464
      m=0                                                                  3464
      mm=0                                                                 3464
      nin=0                                                                3464
      nlp=0                                                                3464
      mnl=min(mnlam,nlam)                                                  3464
      bs=0.0                                                               3465
      shr=shri*dev0                                                        3465
      ga=0.0                                                               3466
20110 do 20111 ic=1,nc                                                     3466
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3466
      svr(ic)=sum(r(:,ic))                                                 3467
20120 do 20121 j=1,ni                                                      3467
      if(ju(j).eq.0)goto 20121                                             3468
      jb=ix(j)                                                             3468
      je=ix(j+1)-1                                                         3469
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             3470
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            3471
20121 continue                                                             3472
20122 continue                                                             3472
20111 continue                                                             3473
20112 continue                                                             3473
      ga=sqrt(ga)                                                          3474
20130 do 20131 ilm=1,nlam                                                  3474
      al0=al                                                               3475
      if(flmin .lt. 1.0)goto 20151                                         3475
      al=ulam(ilm)                                                         3475
      goto 20141                                                           3476
20151 if(ilm .le. 2)goto 20161                                             3476
      al=al*alf                                                            3476
      goto 20141                                                           3477
20161 if(ilm .ne. 1)goto 20171                                             3477
      al=big                                                               3477
      goto 20181                                                           3478
20171 continue                                                             3478
      al0=0.0                                                              3479
20190 do 20191 j=1,ni                                                      3479
      if(ju(j).eq.0)goto 20191                                             3479
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3479
20191 continue                                                             3480
20192 continue                                                             3480
      al0=al0/max(bta,1.0e-3)                                              3480
      al=alf*al0                                                           3481
20181 continue                                                             3482
20141 continue                                                             3482
      al2=al*omb                                                           3482
      al1=al*bta                                                           3482
      tlam=bta*(2.0*al-al0)                                                3483
20200 do 20201 k=1,ni                                                      3483
      if(iy(k).eq.1)goto 20201                                             3483
      if(ju(k).eq.0)goto 20201                                             3484
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      3485
20201 continue                                                             3486
20202 continue                                                             3486
10680 continue                                                             3487
20210 continue                                                             3487
20211 continue                                                             3487
      ixx=0                                                                3487
      jxx=ixx                                                              3488
20220 do 20221 ic=1,nc                                                     3488
      bs(0,ic)=b(0,ic)                                                     3488
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          3489
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3489
      svr(ic)=sum(r(:,ic))                                                 3490
      b(0,ic)=b(0,ic)+svr(ic)                                              3490
      r(:,ic)=r(:,ic)-svr(ic)*w                                            3491
      dlx=max(dlx,svr(ic)**2)                                              3492
20221 continue                                                             3493
20222 continue                                                             3493
20230 continue                                                             3493
20231 continue                                                             3493
      nlp=nlp+nc                                                           3493
      dlx=0.0                                                              3494
20240 do 20241 k=1,ni                                                      3494
      if(iy(k).eq.0)goto 20241                                             3495
      jb=ix(k)                                                             3495
      je=ix(k+1)-1                                                         3495
      del=b(k,:)                                                           3495
      gkn=0.0                                                              3496
20250 do 20251 ic=1,nc                                                     3497
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))-svr(ic)*xb(k))/xs(k)        3498
      gk(ic)=u+del(ic)*xv(k)                                               3498
      gkn=gkn+gk(ic)**2                                                    3499
20251 continue                                                             3500
20252 continue                                                             3500
      gkn=sqrt(gkn)                                                        3500
      u=1.0-al1*vp(k)/gkn                                                  3501
      if(u .gt. 0.0)goto 20271                                             3501
      b(k,:)=0.0                                                           3501
      goto 20281                                                           3502
20271 continue                                                             3502
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2))                                      3502
20281 continue                                                             3503
20261 continue                                                             3503
      del=b(k,:)-del                                                       3503
      if(maxval(abs(del)).le.0.0)goto 20241                                3504
20290 do 20291 ic=1,nc                                                     3504
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3505
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   3507 
     *b(k))/xs(k)
20291 continue                                                             3508
20292 continue                                                             3508
      if(mm(k) .ne. 0)goto 20311                                           3508
      nin=nin+1                                                            3509
      if(nin .le. nx)goto 20331                                            3509
      jxx=1                                                                3509
      goto 20242                                                           3509
20331 continue                                                             3510
      mm(k)=nin                                                            3510
      m(nin)=k                                                             3511
20311 continue                                                             3512
20241 continue                                                             3513
20242 continue                                                             3513
      if(jxx.gt.0)goto 20232                                               3514
      if(dlx.lt.shr)goto 20232                                             3514
      if(nlp .le. maxit)goto 20351                                         3514
      jerr=-ilm                                                            3514
      return                                                               3514
20351 continue                                                             3515
20360 continue                                                             3515
20361 continue                                                             3515
      nlp=nlp+nc                                                           3515
      dlx=0.0                                                              3516
20370 do 20371 l=1,nin                                                     3516
      k=m(l)                                                               3516
      jb=ix(k)                                                             3516
      je=ix(k+1)-1                                                         3516
      del=b(k,:)                                                           3516
      gkn=0.0                                                              3517
20380 do 20381 ic=1,nc                                                     3518
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))  -svr(ic)*xb(k))/xs(k)      3520
      gk(ic)=u+del(ic)*xv(k)                                               3520
      gkn=gkn+gk(ic)**2                                                    3521
20381 continue                                                             3522
20382 continue                                                             3522
      gkn=sqrt(gkn)                                                        3522
      u=1.0-al1*vp(k)/gkn                                                  3523
      if(u .gt. 0.0)goto 20401                                             3523
      b(k,:)=0.0                                                           3523
      goto 20411                                                           3524
20401 continue                                                             3524
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2))                                      3524
20411 continue                                                             3525
20391 continue                                                             3525
      del=b(k,:)-del                                                       3525
      if(maxval(abs(del)).le.0.0)goto 20371                                3526
20420 do 20421 ic=1,nc                                                     3526
      dlx=max(dlx,xv(k)*del(ic)**2)                                        3527
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   3529 
     *b(k))/xs(k)
20421 continue                                                             3530
20422 continue                                                             3530
20371 continue                                                             3531
20372 continue                                                             3531
      if(dlx.lt.shr)goto 20362                                             3531
      if(nlp .le. maxit)goto 20441                                         3531
      jerr=-ilm                                                            3531
      return                                                               3531
20441 continue                                                             3533
      goto 20361                                                           3534
20362 continue                                                             3534
      goto 20231                                                           3535
20232 continue                                                             3535
      if(jxx.gt.0)goto 20212                                               3536
20450 do 20451 ic=1,nc                                                     3537
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                               3538
      if(ixx .ne. 0)goto 20471                                             3539
20480 do 20481 j=1,nin                                                     3539
      k=m(j)                                                               3540
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 20501                   3540
      ixx=1                                                                3540
      goto 20482                                                           3540
20501 continue                                                             3542
20481 continue                                                             3543
20482 continue                                                             3543
20471 continue                                                             3544
      sc=b(0,ic)+g(:,ic)                                                   3544
      b0=0.0                                                               3545
20510 do 20511 j=1,nin                                                     3545
      l=m(j)                                                               3545
      jb=ix(l)                                                             3545
      je=ix(l+1)-1                                                         3546
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   3547
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            3548
20511 continue                                                             3549
20512 continue                                                             3549
      sc=min(max(exmn,sc+b0),exmx)                                         3550
      sxp=sxp-q(:,ic)                                                      3551
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          3552
      sxp=sxp+q(:,ic)                                                      3553
20451 continue                                                             3554
20452 continue                                                             3554
      s=sum(b(0,:))/nc                                                     3554
      b(0,:)=b(0,:)-s                                                      3555
      if(jxx.gt.0)goto 20212                                               3556
      if(ixx .ne. 0)goto 20531                                             3557
20540 do 20541 j=1,ni                                                      3557
      if(iy(j).eq.1)goto 20541                                             3557
      if(ju(j).eq.0)goto 20541                                             3557
      ga(j)=0.0                                                            3557
20541 continue                                                             3558
20542 continue                                                             3558
20550 do 20551 ic=1,nc                                                     3558
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      3559
20560 do 20561 j=1,ni                                                      3559
      if(iy(j).eq.1)goto 20561                                             3559
      if(ju(j).eq.0)goto 20561                                             3560
      jb=ix(j)                                                             3560
      je=ix(j+1)-1                                                         3561
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             3562
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            3563
20561 continue                                                             3564
20562 continue                                                             3564
20551 continue                                                             3565
20552 continue                                                             3565
      ga=sqrt(ga)                                                          3566
20570 do 20571 k=1,ni                                                      3566
      if(iy(k).eq.1)goto 20571                                             3566
      if(ju(k).eq.0)goto 20571                                             3567
      if(ga(k) .le. al1*vp(k))goto 20591                                   3567
      iy(k)=1                                                              3567
      ixx=1                                                                3567
20591 continue                                                             3568
20571 continue                                                             3569
20572 continue                                                             3569
      if(ixx.eq.1) go to 10680                                             3570
      goto 20212                                                           3571
20531 continue                                                             3572
      goto 20211                                                           3573
20212 continue                                                             3573
      if(jxx .le. 0)goto 20611                                             3573
      jerr=-10000-ilm                                                      3573
      goto 20132                                                           3573
20611 continue                                                             3573
      devi=0.0                                                             3574
20620 do 20621 ic=1,nc                                                     3575
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          3575
      a0(ic,ilm)=b(0,ic)                                                   3576
20630 do 20631 i=1,no                                                      3576
      if(y(i,ic).le.0.0)goto 20631                                         3577
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           3578
20631 continue                                                             3579
20632 continue                                                             3579
20621 continue                                                             3580
20622 continue                                                             3580
      kin(ilm)=nin                                                         3580
      alm(ilm)=al                                                          3580
      lmu=ilm                                                              3581
      dev(ilm)=(dev1-devi)/dev0                                            3582
      if(ilm.lt.mnl)goto 20131                                             3582
      if(flmin.ge.1.0)goto 20131                                           3583
      me=0                                                                 3583
20640 do 20641 j=1,nin                                                     3583
      if(a(j,1,ilm).ne.0.0) me=me+1                                        3583
20641 continue                                                             3583
20642 continue                                                             3583
      if(me.gt.ne)goto 20132                                               3584
      if(dev(ilm).gt.devmax)goto 20132                                     3584
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 20132                             3585
20131 continue                                                             3586
20132 continue                                                             3586
      g=log(q)                                                             3586
20650 do 20651 i=1,no                                                      3586
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3586
20651 continue                                                             3587
20652 continue                                                             3587
      deallocate(sxp,b,bs,r,q,mm,is,sc,ga,iy,gk,del,sxpl)                  3588
      return                                                               3589
      end                                                                  3590
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
