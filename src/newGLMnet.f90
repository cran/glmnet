c
c
c                          newGLMnet (2/17/12)
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
      subroutine elnet  (ka,parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam    609 
     *,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam)                          610
      real ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)                          611
      integer jd(*),ia(nx),nin(nlam)                                        612
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10021                                     615
      jerr=10000                                                            615
      return                                                                615
10021 continue                                                              616
      allocate(vq(1:ni),stat=jerr)                                          616
      if(jerr.ne.0) return                                                  617
      vq=max(0.0,vp)                                                        617
      vq=vq*ni/sum(vq)                                                      618
      if(ka .ne. 1)goto 10041                                               619
      call elnetu  (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd    622 
     *,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                            623
10041 continue                                                              624
      call elnetn (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd,    627 
     *maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                              628
10031 continue                                                              628
      deallocate(vq)                                                        629
      return                                                                630
      end                                                                   631
      subroutine elnetu  (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,t    634 
     *hr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ulam(nlam)                           635
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         636
      integer jd(*),ia(nx),nin(nlam)                                        637
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           642
      allocate(xm(1:ni),stat=ierr)                                          642
      jerr=jerr+ierr                                                        643
      allocate(xs(1:ni),stat=ierr)                                          643
      jerr=jerr+ierr                                                        644
      allocate(ju(1:ni),stat=ierr)                                          644
      jerr=jerr+ierr                                                        645
      allocate(xv(1:ni),stat=ierr)                                          645
      jerr=jerr+ierr                                                        646
      allocate(vlam(1:nlam),stat=ierr)                                      646
      jerr=jerr+ierr                                                        647
      if(jerr.ne.0) return                                                  648
      call chkvars(no,ni,x,ju)                                              649
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  650
      if(maxval(ju) .gt. 0)goto 10071                                       650
      jerr=7777                                                             650
      return                                                                650
10071 continue                                                              651
      call standard(no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)               652
      if(jerr.ne.0) return                                                  653
      if(flmin.ge.1.0) vlam=ulam/ys                                         654
      call elnet1(parm,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxit,x    656 
     *v,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  657
10080 do 10081 k=1,lmu                                                      657
      alm(k)=ys*alm(k)                                                      657
      nk=nin(k)                                                             658
10090 do 10091 l=1,nk                                                       658
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          658
10091 continue                                                              659
10092 continue                                                              659
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         660
10081 continue                                                              661
10082 continue                                                              661
      deallocate(xm,xs,g,ju,xv,vlam)                                        662
      return                                                                663
      end                                                                   664
      subroutine standard (no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)        665
      real x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                  665
      integer ju(ni)                                                        666
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           669
      if(jerr.ne.0) return                                                  670
      w=w/sum(w)                                                            670
      v=sqrt(w)                                                             671
10100 do 10101 j=1,ni                                                       671
      if(ju(j).eq.0)goto 10101                                              672
      xm(j)=dot_product(w,x(:,j))                                           672
      x(:,j)=v*(x(:,j)-xm(j))                                               673
      xv(j)=dot_product(x(:,j),x(:,j))                                      673
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        674
10101 continue                                                              675
10102 continue                                                              675
      if(isd .ne. 0)goto 10121                                              675
      xs=1.0                                                                675
      goto 10131                                                            676
10121 continue                                                              677
10140 do 10141 j=1,ni                                                       677
      if(ju(j).eq.0)goto 10141                                              677
      x(:,j)=x(:,j)/xs(j)                                                   677
10141 continue                                                              678
10142 continue                                                              678
      xv=1.0                                                                679
10131 continue                                                              680
10111 continue                                                              680
      ym=dot_product(w,y)                                                   680
      y=v*(y-ym)                                                            680
      ys=sqrt(dot_product(y,y))                                             680
      y=y/ys                                                                680
      g=0.0                                                                 681
10150 do 10151 j=1,ni                                                       681
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                             681
10151 continue                                                              682
10152 continue                                                              682
      deallocate(v)                                                         683
      return                                                                684
      end                                                                   685
      subroutine elnet1 (beta,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,ulam,thr,    687 
     *maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    688 
     *9)
      real vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    689 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       690
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                           696
      jerr=jerr+ierr                                                        697
      allocate(mm(1:ni),stat=ierr)                                          697
      jerr=jerr+ierr                                                        698
      allocate(da(1:ni),stat=ierr)                                          698
      jerr=jerr+ierr                                                        699
      if(jerr.ne.0) return                                                  700
      bta=beta                                                              700
      omb=1.0-bta                                                           701
      if(flmin .ge. 1.0)goto 10171                                          701
      eqs=max(eps,flmin)                                                    701
      alf=eqs**(1.0/(nlam-1))                                               701
10171 continue                                                              702
      rsq=0.0                                                               702
      a=0.0                                                                 702
      mm=0                                                                  702
      nlp=0                                                                 702
      nin=nlp                                                               702
      iz=0                                                                  702
      mnl=min(mnlam,nlam)                                                   703
10180 do 10181 m=1,nlam                                                     704
      if(flmin .lt. 1.0)goto 10201                                          704
      alm=ulam(m)                                                           704
      goto 10191                                                            705
10201 if(m .le. 2)goto 10211                                                705
      alm=alm*alf                                                           705
      goto 10191                                                            706
10211 if(m .ne. 1)goto 10221                                                706
      alm=big                                                               706
      goto 10231                                                            707
10221 continue                                                              707
      alm=0.0                                                               708
10240 do 10241 j=1,ni                                                       708
      if(ju(j).eq.0)goto 10241                                              708
      if(vp(j).le.0.0)goto 10241                                            709
      alm=max(alm,abs(g(j))/vp(j))                                          710
10241 continue                                                              711
10242 continue                                                              711
      alm=alf*alm/max(bta,1.0e-3)                                           712
10231 continue                                                              713
10191 continue                                                              713
      dem=alm*omb                                                           713
      ab=alm*bta                                                            713
      rsq0=rsq                                                              713
      jz=1                                                                  714
10250 continue                                                              714
10251 continue                                                              714
      if(iz*jz.ne.0) go to 10260                                            714
      nlp=nlp+1                                                             714
      dlx=0.0                                                               715
10270 do 10271 k=1,ni                                                       715
      if(ju(k).eq.0)goto 10271                                              716
      ak=a(k)                                                               716
      u=g(k)+ak*xv(k)                                                       716
      v=abs(u)-vp(k)*ab                                                     716
      a(k)=0.0                                                              717
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         718
      if(a(k).eq.ak)goto 10271                                              719
      if(mm(k) .ne. 0)goto 10291                                            719
      nin=nin+1                                                             719
      if(nin.gt.nx)goto 10272                                               720
10300 do 10301 j=1,ni                                                       720
      if(ju(j).eq.0)goto 10301                                              721
      if(mm(j) .eq. 0)goto 10321                                            721
      c(j,nin)=c(k,mm(j))                                                   721
      goto 10301                                                            721
10321 continue                                                              722
      if(j .ne. k)goto 10341                                                722
      c(j,nin)=xv(j)                                                        722
      goto 10301                                                            722
10341 continue                                                              723
      c(j,nin)=dot_product(x(:,j),x(:,k))                                   724
10301 continue                                                              725
10302 continue                                                              725
      mm(k)=nin                                                             725
      ia(nin)=k                                                             726
10291 continue                                                              727
      del=a(k)-ak                                                           727
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      728
      dlx=max(xv(k)*del**2,dlx)                                             729
10350 do 10351 j=1,ni                                                       729
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                               729
10351 continue                                                              730
10352 continue                                                              730
10271 continue                                                              731
10272 continue                                                              731
      if(dlx.lt.thr)goto 10252                                              731
      if(nin.gt.nx)goto 10252                                               732
      if(nlp .le. maxit)goto 10371                                          732
      jerr=-m                                                               732
      return                                                                732
10371 continue                                                              733
10260 continue                                                              733
      iz=1                                                                  733
      da(1:nin)=a(ia(1:nin))                                                734
10380 continue                                                              734
10381 continue                                                              734
      nlp=nlp+1                                                             734
      dlx=0.0                                                               735
10390 do 10391 l=1,nin                                                      735
      k=ia(l)                                                               735
      ak=a(k)                                                               735
      u=g(k)+ak*xv(k)                                                       735
      v=abs(u)-vp(k)*ab                                                     736
      a(k)=0.0                                                              737
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         738
      if(a(k).eq.ak)goto 10391                                              739
      del=a(k)-ak                                                           739
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      740
      dlx=max(xv(k)*del**2,dlx)                                             741
10400 do 10401 j=1,nin                                                      741
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                  741
10401 continue                                                              742
10402 continue                                                              742
10391 continue                                                              743
10392 continue                                                              743
      if(dlx.lt.thr)goto 10382                                              743
      if(nlp .le. maxit)goto 10421                                          743
      jerr=-m                                                               743
      return                                                                743
10421 continue                                                              744
      goto 10381                                                            745
10382 continue                                                              745
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                      746
10430 do 10431 j=1,ni                                                       746
      if(mm(j).ne.0)goto 10431                                              747
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))            748
10431 continue                                                              749
10432 continue                                                              749
      jz=0                                                                  750
      goto 10251                                                            751
10252 continue                                                              751
      if(nin .le. nx)goto 10451                                             751
      jerr=-10000-m                                                         751
      goto 10182                                                            751
10451 continue                                                              752
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 752
      kin(m)=nin                                                            753
      rsqo(m)=rsq                                                           753
      almo(m)=alm                                                           753
      lmu=m                                                                 754
      if(m.lt.mnl)goto 10181                                                754
      if(flmin.ge.1.0)goto 10181                                            755
      me=0                                                                  755
10460 do 10461 j=1,nin                                                      755
      if(ao(j,m).ne.0.0) me=me+1                                            755
10461 continue                                                              755
10462 continue                                                              755
      if(me.gt.ne)goto 10182                                                756
      if(rsq-rsq0.lt.sml*rsq)goto 10182                                     756
      if(rsq.gt.rsqmax)goto 10182                                           757
10181 continue                                                              758
10182 continue                                                              758
      deallocate(a,mm,c,da)                                                 759
      return                                                                760
      end                                                                   761
      subroutine elnetn (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,th    763 
     *r,isd,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no),w(no),ulam(nlam)                           764
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         765
      integer jd(*),ia(nx),nin(nlam)                                        766
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          771
      allocate(xs(1:ni),stat=ierr)                                          771
      jerr=jerr+ierr                                                        772
      allocate(ju(1:ni),stat=ierr)                                          772
      jerr=jerr+ierr                                                        773
      allocate(xv(1:ni),stat=ierr)                                          773
      jerr=jerr+ierr                                                        774
      allocate(vlam(1:nlam),stat=ierr)                                      774
      jerr=jerr+ierr                                                        775
      if(jerr.ne.0) return                                                  776
      call chkvars(no,ni,x,ju)                                              777
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  778
      if(maxval(ju) .gt. 0)goto 10481                                       778
      jerr=7777                                                             778
      return                                                                778
10481 continue                                                              779
      call standard1(no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)                780
      if(jerr.ne.0) return                                                  781
      if(flmin.ge.1.0) vlam=ulam/ys                                         782
      call elnet2(parm,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,vlam,thr,maxit,x    784 
     *v,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  785
10490 do 10491 k=1,lmu                                                      785
      alm(k)=ys*alm(k)                                                      785
      nk=nin(k)                                                             786
10500 do 10501 l=1,nk                                                       786
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          786
10501 continue                                                              787
10502 continue                                                              787
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         788
10491 continue                                                              789
10492 continue                                                              789
      deallocate(xm,xs,ju,xv,vlam)                                          790
      return                                                                791
      end                                                                   792
      subroutine standard1 (no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)         793
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)                        793
      integer ju(ni)                                                        794
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           797
      if(jerr.ne.0) return                                                  798
      w=w/sum(w)                                                            798
      v=sqrt(w)                                                             799
10510 do 10511 j=1,ni                                                       799
      if(ju(j).eq.0)goto 10511                                              800
      xm(j)=dot_product(w,x(:,j))                                           800
      x(:,j)=v*(x(:,j)-xm(j))                                               801
      xv(j)=dot_product(x(:,j),x(:,j))                                      801
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        802
10511 continue                                                              803
10512 continue                                                              803
      if(isd .ne. 0)goto 10531                                              803
      xs=1.0                                                                803
      goto 10541                                                            804
10531 continue                                                              804
10550 do 10551 j=1,ni                                                       804
      if(ju(j).eq.0)goto 10551                                              804
      x(:,j)=x(:,j)/xs(j)                                                   804
10551 continue                                                              805
10552 continue                                                              805
      xv=1.0                                                                806
10541 continue                                                              807
10521 continue                                                              807
      ym=dot_product(w,y)                                                   807
      y=v*(y-ym)                                                            807
      ys=sqrt(dot_product(y,y))                                             807
      y=y/ys                                                                808
      deallocate(v)                                                         809
      return                                                                810
      end                                                                   811
      subroutine elnet2(beta,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,ulam,thr,m    813 
     *axit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    814 
     *9)
      real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    815 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       816
      real, dimension (:), allocatable :: a,g                                   
      integer, dimension (:), allocatable :: mm,ix                              
      allocate(a(1:ni),stat=jerr)                                           821
      allocate(mm(1:ni),stat=ierr)                                          821
      jerr=jerr+ierr                                                        822
      allocate(g(1:ni),stat=ierr)                                           822
      jerr=jerr+ierr                                                        823
      allocate(ix(1:ni),stat=ierr)                                          823
      jerr=jerr+ierr                                                        824
      if(jerr.ne.0) return                                                  825
      bta=beta                                                              825
      omb=1.0-bta                                                           825
      ix=0                                                                  826
      if(flmin .ge. 1.0)goto 10571                                          826
      eqs=max(eps,flmin)                                                    826
      alf=eqs**(1.0/(nlam-1))                                               826
10571 continue                                                              827
      rsq=0.0                                                               827
      a=0.0                                                                 827
      mm=0                                                                  827
      nlp=0                                                                 827
      nin=nlp                                                               827
      iz=0                                                                  827
      mnl=min(mnlam,nlam)                                                   827
      alm=0.0                                                               828
10580 do 10581 j=1,ni                                                       828
      if(ju(j).eq.0)goto 10581                                              828
      g(j)=abs(dot_product(y,x(:,j)))                                       828
10581 continue                                                              829
10582 continue                                                              829
10590 do 10591 m=1,nlam                                                     829
      alm0=alm                                                              830
      if(flmin .lt. 1.0)goto 10611                                          830
      alm=ulam(m)                                                           830
      goto 10601                                                            831
10611 if(m .le. 2)goto 10621                                                831
      alm=alm*alf                                                           831
      goto 10601                                                            832
10621 if(m .ne. 1)goto 10631                                                832
      alm=big                                                               832
      goto 10641                                                            833
10631 continue                                                              833
      alm0=0.0                                                              834
10650 do 10651 j=1,ni                                                       834
      if(ju(j).eq.0)goto 10651                                              834
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                            834
10651 continue                                                              835
10652 continue                                                              835
      alm0=alm0/max(bta,1.0e-3)                                             835
      alm=alf*alm0                                                          836
10641 continue                                                              837
10601 continue                                                              837
      dem=alm*omb                                                           837
      ab=alm*bta                                                            837
      rsq0=rsq                                                              837
      jz=1                                                                  838
      tlam=bta*(2.0*alm-alm0)                                               839
10660 do 10661 k=1,ni                                                       839
      if(ix(k).eq.1)goto 10661                                              839
      if(ju(k).eq.0)goto 10661                                              840
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                        841
10661 continue                                                              842
10662 continue                                                              842
10670 continue                                                              842
10671 continue                                                              842
      if(iz*jz.ne.0) go to 10260                                            843
10680 continue                                                              843
      nlp=nlp+1                                                             843
      dlx=0.0                                                               844
10690 do 10691 k=1,ni                                                       844
      if(ix(k).eq.0)goto 10691                                              844
      gk=dot_product(y,x(:,k))                                              845
      ak=a(k)                                                               845
      u=gk+ak*xv(k)                                                         845
      v=abs(u)-vp(k)*ab                                                     845
      a(k)=0.0                                                              846
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         847
      if(a(k).eq.ak)goto 10691                                              848
      if(mm(k) .ne. 0)goto 10711                                            848
      nin=nin+1                                                             848
      if(nin.gt.nx)goto 10692                                               849
      mm(k)=nin                                                             849
      ia(nin)=k                                                             850
10711 continue                                                              851
      del=a(k)-ak                                                           851
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        852
      y=y-del*x(:,k)                                                        852
      dlx=max(xv(k)*del**2,dlx)                                             853
10691 continue                                                              854
10692 continue                                                              854
      if(nin.gt.nx)goto 10672                                               855
      if(dlx .ge. thr)goto 10731                                            855
      ixx=0                                                                 856
10740 do 10741 k=1,ni                                                       856
      if(ix(k).eq.1)goto 10741                                              856
      if(ju(k).eq.0)goto 10741                                              857
      g(k)=abs(dot_product(y,x(:,k)))                                       858
      if(g(k) .le. ab*vp(k))goto 10761                                      858
      ix(k)=1                                                               858
      ixx=1                                                                 858
10761 continue                                                              859
10741 continue                                                              860
10742 continue                                                              860
      if(ixx.eq.1) go to 10680                                              861
      goto 10672                                                            862
10731 continue                                                              863
      if(nlp .le. maxit)goto 10781                                          863
      jerr=-m                                                               863
      return                                                                863
10781 continue                                                              864
10260 continue                                                              864
      iz=1                                                                  865
10790 continue                                                              865
10791 continue                                                              865
      nlp=nlp+1                                                             865
      dlx=0.0                                                               866
10800 do 10801 l=1,nin                                                      866
      k=ia(l)                                                               866
      gk=dot_product(y,x(:,k))                                              867
      ak=a(k)                                                               867
      u=gk+ak*xv(k)                                                         867
      v=abs(u)-vp(k)*ab                                                     867
      a(k)=0.0                                                              868
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         869
      if(a(k).eq.ak)goto 10801                                              870
      del=a(k)-ak                                                           870
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        871
      y=y-del*x(:,k)                                                        871
      dlx=max(xv(k)*del**2,dlx)                                             872
10801 continue                                                              873
10802 continue                                                              873
      if(dlx.lt.thr)goto 10792                                              873
      if(nlp .le. maxit)goto 10821                                          873
      jerr=-m                                                               873
      return                                                                873
10821 continue                                                              874
      goto 10791                                                            875
10792 continue                                                              875
      jz=0                                                                  876
      goto 10671                                                            877
10672 continue                                                              877
      if(nin .le. nx)goto 10841                                             877
      jerr=-10000-m                                                         877
      goto 10592                                                            877
10841 continue                                                              878
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 878
      kin(m)=nin                                                            879
      rsqo(m)=rsq                                                           879
      almo(m)=alm                                                           879
      lmu=m                                                                 880
      if(m.lt.mnl)goto 10591                                                880
      if(flmin.ge.1.0)goto 10591                                            881
      me=0                                                                  881
10850 do 10851 j=1,nin                                                      881
      if(ao(j,m).ne.0.0) me=me+1                                            881
10851 continue                                                              881
10852 continue                                                              881
      if(me.gt.ne)goto 10592                                                882
      if(rsq-rsq0.lt.sml*rsq)goto 10592                                     882
      if(rsq.gt.rsqmax)goto 10592                                           883
10591 continue                                                              884
10592 continue                                                              884
      deallocate(a,mm,g,ix)                                                 885
      return                                                                886
      end                                                                   887
      subroutine chkvars(no,ni,x,ju)                                        888
      real x(no,ni)                                                         888
      integer ju(ni)                                                        889
10860 do 10861 j=1,ni                                                       889
      ju(j)=0                                                               889
      t=x(1,j)                                                              890
10870 do 10871 i=2,no                                                       890
      if(x(i,j).eq.t)goto 10871                                             890
      ju(j)=1                                                               890
      goto 10872                                                            890
10871 continue                                                              891
10872 continue                                                              891
10861 continue                                                              892
10862 continue                                                              892
      return                                                                893
      end                                                                   894
      subroutine uncomp(ni,ca,ia,nin,a)                                     895
      real ca(*),a(ni)                                                      895
      integer ia(*)                                                         896
      a=0.0                                                                 896
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                   897
      return                                                                898
      end                                                                   899
      subroutine modval(a0,ca,ia,nin,n,x,f)                                 900
      real ca(nin),x(n,*),f(n)                                              900
      integer ia(nin)                                                       901
      f=a0                                                                  901
      if(nin.le.0) return                                                   902
10880 do 10881 i=1,n                                                        902
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                       902
10881 continue                                                              903
10882 continue                                                              903
      return                                                                904
      end                                                                   905
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,fl    908 
     *min,ulam,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                               909
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         910
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                            911
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10901                                     914
      jerr=10000                                                            914
      return                                                                914
10901 continue                                                              915
      allocate(vq(1:ni),stat=jerr)                                          915
      if(jerr.ne.0) return                                                  916
      vq=max(0.0,vp)                                                        916
      vq=vq*ni/sum(vq)                                                      917
      if(ka .ne. 1)goto 10921                                               918
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam    921 
     *,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10931                                                            922
10921 continue                                                              923
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam,    926 
     *thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10931 continue                                                              927
10911 continue                                                              927
      deallocate(vq)                                                        928
      return                                                                929
      end                                                                   930
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmi    933 
     *n,ulam,thr,isd,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                               934
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         935
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                            936
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           941
      allocate(xm(1:ni),stat=ierr)                                          941
      jerr=jerr+ierr                                                        942
      allocate(xs(1:ni),stat=ierr)                                          942
      jerr=jerr+ierr                                                        943
      allocate(ju(1:ni),stat=ierr)                                          943
      jerr=jerr+ierr                                                        944
      allocate(xv(1:ni),stat=ierr)                                          944
      jerr=jerr+ierr                                                        945
      allocate(vlam(1:nlam),stat=ierr)                                      945
      jerr=jerr+ierr                                                        946
      if(jerr.ne.0) return                                                  947
      call spchkvars(no,ni,x,ix,ju)                                         948
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  949
      if(maxval(ju) .gt. 0)goto 10951                                       949
      jerr=7777                                                             949
      return                                                                949
10951 continue                                                              950
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,jerr)       951
      if(jerr.ne.0) return                                                  952
      if(flmin.ge.1.0) vlam=ulam/ys                                         953
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t    955 
     *hr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  956
10960 do 10961 k=1,lmu                                                      956
      alm(k)=ys*alm(k)                                                      956
      nk=nin(k)                                                             957
10970 do 10971 l=1,nk                                                       957
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          957
10971 continue                                                              958
10972 continue                                                              958
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         959
10961 continue                                                              960
10962 continue                                                              960
      deallocate(xm,xs,g,ju,xv,vlam)                                        961
      return                                                                962
      end                                                                   963
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,j    964 
     *err)
      real x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                      964
      integer ix(*),jx(*),ju(ni)                                            965
      w=w/sum(w)                                                            966
10980 do 10981 j=1,ni                                                       966
      if(ju(j).eq.0)goto 10981                                              967
      jb=ix(j)                                                              967
      je=ix(j+1)-1                                                          967
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                              968
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                  969
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        970
10981 continue                                                              971
10982 continue                                                              971
      if(isd .ne. 0)goto 11001                                              971
      xs=1.0                                                                971
      goto 11011                                                            971
11001 continue                                                              971
      xv=1.0                                                                971
11011 continue                                                              972
10991 continue                                                              972
      ym=dot_product(w,y)                                                   972
      y=y-ym                                                                972
      ys=sqrt(dot_product(w,y**2))                                          972
      y=y/ys                                                                972
      g=0.0                                                                 973
11020 do 11021 j=1,ni                                                       973
      if(ju(j).eq.0)goto 11021                                              973
      jb=ix(j)                                                              973
      je=ix(j+1)-1                                                          974
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)            975
11021 continue                                                              976
11022 continue                                                              976
      return                                                                977
      end                                                                   978
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,    980 
     *ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    981 
     *9)
      real g(ni),vp(ni),x(*),ulam(nlam),w(no)                               982
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)           983
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                           984
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                           990
      jerr=jerr+ierr                                                        991
      allocate(mm(1:ni),stat=ierr)                                          991
      jerr=jerr+ierr                                                        992
      allocate(da(1:ni),stat=ierr)                                          992
      jerr=jerr+ierr                                                        993
      if(jerr.ne.0) return                                                  994
      bta=beta                                                              994
      omb=1.0-bta                                                           995
      if(flmin .ge. 1.0)goto 11041                                          995
      eqs=max(eps,flmin)                                                    995
      alf=eqs**(1.0/(nlam-1))                                               995
11041 continue                                                              996
      rsq=0.0                                                               996
      a=0.0                                                                 996
      mm=0                                                                  996
      nlp=0                                                                 996
      nin=nlp                                                               996
      iz=0                                                                  996
      mnl=min(mnlam,nlam)                                                   997
11050 do 11051 m=1,nlam                                                     998
      if(flmin .lt. 1.0)goto 11071                                          998
      alm=ulam(m)                                                           998
      goto 11061                                                            999
11071 if(m .le. 2)goto 11081                                                999
      alm=alm*alf                                                           999
      goto 11061                                                           1000
11081 if(m .ne. 1)goto 11091                                               1000
      alm=big                                                              1000
      goto 11101                                                           1001
11091 continue                                                             1001
      alm=0.0                                                              1002
11110 do 11111 j=1,ni                                                      1002
      if(ju(j).eq.0)goto 11111                                             1002
      if(vp(j).le.0.0)goto 11111                                           1003
      alm=max(alm,abs(g(j))/vp(j))                                         1004
11111 continue                                                             1005
11112 continue                                                             1005
      alm=alf*alm/max(bta,1.0e-3)                                          1006
11101 continue                                                             1007
11061 continue                                                             1007
      dem=alm*omb                                                          1007
      ab=alm*bta                                                           1007
      rsq0=rsq                                                             1007
      jz=1                                                                 1008
11120 continue                                                             1008
11121 continue                                                             1008
      if(iz*jz.ne.0) go to 10260                                           1008
      nlp=nlp+1                                                            1008
      dlx=0.0                                                              1009
11130 do 11131 k=1,ni                                                      1009
      if(ju(k).eq.0)goto 11131                                             1010
      ak=a(k)                                                              1010
      u=g(k)+ak*xv(k)                                                      1010
      v=abs(u)-vp(k)*ab                                                    1010
      a(k)=0.0                                                             1011
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1012
      if(a(k).eq.ak)goto 11131                                             1013
      if(mm(k) .ne. 0)goto 11151                                           1013
      nin=nin+1                                                            1013
      if(nin.gt.nx)goto 11132                                              1014
11160 do 11161 j=1,ni                                                      1014
      if(ju(j).eq.0)goto 11161                                             1015
      if(mm(j) .eq. 0)goto 11181                                           1015
      c(j,nin)=c(k,mm(j))                                                  1015
      goto 11161                                                           1015
11181 continue                                                             1016
      if(j .ne. k)goto 11201                                               1016
      c(j,nin)=xv(j)                                                       1016
      goto 11161                                                           1016
11201 continue                                                             1017
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))       1019
11161 continue                                                             1020
11162 continue                                                             1020
      mm(k)=nin                                                            1020
      ia(nin)=k                                                            1021
11151 continue                                                             1022
      del=a(k)-ak                                                          1022
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1023
      dlx=max(xv(k)*del**2,dlx)                                            1024
11210 do 11211 j=1,ni                                                      1024
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                              1024
11211 continue                                                             1025
11212 continue                                                             1025
11131 continue                                                             1026
11132 continue                                                             1026
      if(dlx.lt.thr)goto 11122                                             1026
      if(nin.gt.nx)goto 11122                                              1027
      if(nlp .le. maxit)goto 11231                                         1027
      jerr=-m                                                              1027
      return                                                               1027
11231 continue                                                             1028
10260 continue                                                             1028
      iz=1                                                                 1028
      da(1:nin)=a(ia(1:nin))                                               1029
11240 continue                                                             1029
11241 continue                                                             1029
      nlp=nlp+1                                                            1029
      dlx=0.0                                                              1030
11250 do 11251 l=1,nin                                                     1030
      k=ia(l)                                                              1031
      ak=a(k)                                                              1031
      u=g(k)+ak*xv(k)                                                      1031
      v=abs(u)-vp(k)*ab                                                    1031
      a(k)=0.0                                                             1032
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1033
      if(a(k).eq.ak)goto 11251                                             1034
      del=a(k)-ak                                                          1034
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1035
      dlx=max(xv(k)*del**2,dlx)                                            1036
11260 do 11261 j=1,nin                                                     1036
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                 1036
11261 continue                                                             1037
11262 continue                                                             1037
11251 continue                                                             1038
11252 continue                                                             1038
      if(dlx.lt.thr)goto 11242                                             1038
      if(nlp .le. maxit)goto 11281                                         1038
      jerr=-m                                                              1038
      return                                                               1038
11281 continue                                                             1039
      goto 11241                                                           1040
11242 continue                                                             1040
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                     1041
11290 do 11291 j=1,ni                                                      1041
      if(mm(j).ne.0)goto 11291                                             1042
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))           1043
11291 continue                                                             1044
11292 continue                                                             1044
      jz=0                                                                 1045
      goto 11121                                                           1046
11122 continue                                                             1046
      if(nin .le. nx)goto 11311                                            1046
      jerr=-10000-m                                                        1046
      goto 11052                                                           1046
11311 continue                                                             1047
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1047
      kin(m)=nin                                                           1048
      rsqo(m)=rsq                                                          1048
      almo(m)=alm                                                          1048
      lmu=m                                                                1049
      if(m.lt.mnl)goto 11051                                               1049
      if(flmin.ge.1.0)goto 11051                                           1050
      me=0                                                                 1050
11320 do 11321 j=1,nin                                                     1050
      if(ao(j,m).ne.0.0) me=me+1                                           1050
11321 continue                                                             1050
11322 continue                                                             1050
      if(me.gt.ne)goto 11052                                               1051
      if(rsq-rsq0.lt.sml*rsq)goto 11052                                    1051
      if(rsq.gt.rsqmax)goto 11052                                          1052
11051 continue                                                             1053
11052 continue                                                             1053
      deallocate(a,mm,c,da)                                                1054
      return                                                               1055
      end                                                                  1056
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmin,   1058 
     *ulam,  thr,isd,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),vp(ni),y(no),w(no),ulam(nlam)                              1059
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1060
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1061
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         1066
      allocate(xs(1:ni),stat=ierr)                                         1066
      jerr=jerr+ierr                                                       1067
      allocate(ju(1:ni),stat=ierr)                                         1067
      jerr=jerr+ierr                                                       1068
      allocate(xv(1:ni),stat=ierr)                                         1068
      jerr=jerr+ierr                                                       1069
      allocate(vlam(1:nlam),stat=ierr)                                     1069
      jerr=jerr+ierr                                                       1070
      if(jerr.ne.0) return                                                 1071
      call spchkvars(no,ni,x,ix,ju)                                        1072
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1073
      if(maxval(ju) .gt. 0)goto 11341                                      1073
      jerr=7777                                                            1073
      return                                                               1073
11341 continue                                                             1074
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,jerr)       1075
      if(jerr.ne.0) return                                                 1076
      if(flmin.ge.1.0) vlam=ulam/ys                                        1077
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t   1079 
     *hr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1080
11350 do 11351 k=1,lmu                                                     1080
      alm(k)=ys*alm(k)                                                     1080
      nk=nin(k)                                                            1081
11360 do 11361 l=1,nk                                                      1081
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1081
11361 continue                                                             1082
11362 continue                                                             1082
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                        1083
11351 continue                                                             1084
11352 continue                                                             1084
      deallocate(xm,xs,ju,xv,vlam)                                         1085
      return                                                               1086
      end                                                                  1087
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,je   1088 
     *rr)
      real x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)                           1088
      integer ix(*),jx(*),ju(ni)                                           1089
      w=w/sum(w)                                                           1090
11370 do 11371 j=1,ni                                                      1090
      if(ju(j).eq.0)goto 11371                                             1091
      jb=ix(j)                                                             1091
      je=ix(j+1)-1                                                         1091
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1092
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1093
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1094
11371 continue                                                             1095
11372 continue                                                             1095
      if(isd .ne. 0)goto 11391                                             1095
      xs=1.0                                                               1095
      goto 11401                                                           1095
11391 continue                                                             1095
      xv=1.0                                                               1095
11401 continue                                                             1096
11381 continue                                                             1096
      ym=dot_product(w,y)                                                  1096
      y=y-ym                                                               1096
      ys=sqrt(dot_product(w,y**2))                                         1096
      y=y/ys                                                               1097
      return                                                               1098
      end                                                                  1099
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,   1101 
     *ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99   1102 
     *9)
      real y(no),w(no),x(*),vp(ni),ulam(nlam)                              1103
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)          1104
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1105
      real, dimension (:), allocatable :: a,g                                   
      integer, dimension (:), allocatable :: mm,iy                              
      allocate(a(1:ni),stat=jerr)                                          1110
      allocate(mm(1:ni),stat=ierr)                                         1110
      jerr=jerr+ierr                                                       1111
      allocate(g(1:ni),stat=ierr)                                          1111
      jerr=jerr+ierr                                                       1112
      allocate(iy(1:ni),stat=ierr)                                         1112
      jerr=jerr+ierr                                                       1113
      if(jerr.ne.0) return                                                 1114
      bta=beta                                                             1114
      omb=1.0-bta                                                          1114
      alm=0.0                                                              1114
      iy=0                                                                 1115
      if(flmin .ge. 1.0)goto 11421                                         1115
      eqs=max(eps,flmin)                                                   1115
      alf=eqs**(1.0/(nlam-1))                                              1115
11421 continue                                                             1116
      rsq=0.0                                                              1116
      a=0.0                                                                1116
      mm=0                                                                 1116
      o=0.0                                                                1116
      nlp=0                                                                1116
      nin=nlp                                                              1116
      iz=0                                                                 1116
      mnl=min(mnlam,nlam)                                                  1117
11430 do 11431 j=1,ni                                                      1117
      if(ju(j).eq.0)goto 11431                                             1118
      jb=ix(j)                                                             1118
      je=ix(j+1)-1                                                         1119
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1120
11431 continue                                                             1121
11432 continue                                                             1121
11440 do 11441 m=1,nlam                                                    1121
      alm0=alm                                                             1122
      if(flmin .lt. 1.0)goto 11461                                         1122
      alm=ulam(m)                                                          1122
      goto 11451                                                           1123
11461 if(m .le. 2)goto 11471                                               1123
      alm=alm*alf                                                          1123
      goto 11451                                                           1124
11471 if(m .ne. 1)goto 11481                                               1124
      alm=big                                                              1124
      goto 11491                                                           1125
11481 continue                                                             1125
      alm0=0.0                                                             1126
11500 do 11501 j=1,ni                                                      1126
      if(ju(j).eq.0)goto 11501                                             1126
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           1126
11501 continue                                                             1127
11502 continue                                                             1127
      alm0=alm0/max(bta,1.0e-3)                                            1127
      alm=alf*alm0                                                         1128
11491 continue                                                             1129
11451 continue                                                             1129
      dem=alm*omb                                                          1129
      ab=alm*bta                                                           1129
      rsq0=rsq                                                             1129
      jz=1                                                                 1130
      tlam=bta*(2.0*alm-alm0)                                              1131
11510 do 11511 k=1,ni                                                      1131
      if(iy(k).eq.1)goto 11511                                             1131
      if(ju(k).eq.0)goto 11511                                             1132
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       1133
11511 continue                                                             1134
11512 continue                                                             1134
11520 continue                                                             1134
11521 continue                                                             1134
      if(iz*jz.ne.0) go to 10260                                           1135
10680 continue                                                             1135
      nlp=nlp+1                                                            1135
      dlx=0.0                                                              1136
11530 do 11531 k=1,ni                                                      1136
      if(iy(k).eq.0)goto 11531                                             1136
      jb=ix(k)                                                             1136
      je=ix(k+1)-1                                                         1137
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1138
      ak=a(k)                                                              1138
      u=gk+ak*xv(k)                                                        1138
      v=abs(u)-vp(k)*ab                                                    1138
      a(k)=0.0                                                             1139
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1140
      if(a(k).eq.ak)goto 11531                                             1141
      if(mm(k) .ne. 0)goto 11551                                           1141
      nin=nin+1                                                            1141
      if(nin.gt.nx)goto 11532                                              1142
      mm(k)=nin                                                            1142
      ia(nin)=k                                                            1143
11551 continue                                                             1144
      del=a(k)-ak                                                          1144
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1145
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1146
      o=o+del*xm(k)/xs(k)                                                  1146
      dlx=max(xv(k)*del**2,dlx)                                            1147
11531 continue                                                             1148
11532 continue                                                             1148
      if(nin.gt.nx)goto 11522                                              1149
      if(dlx .ge. thr)goto 11571                                           1149
      ixx=0                                                                1150
11580 do 11581 j=1,ni                                                      1150
      if(iy(j).eq.1)goto 11581                                             1150
      if(ju(j).eq.0)goto 11581                                             1151
      jb=ix(j)                                                             1151
      je=ix(j+1)-1                                                         1152
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1153
      if(g(j) .le. ab*vp(j))goto 11601                                     1153
      iy(j)=1                                                              1153
      ixx=1                                                                1153
11601 continue                                                             1154
11581 continue                                                             1155
11582 continue                                                             1155
      if(ixx.eq.1) go to 10680                                             1156
      goto 11522                                                           1157
11571 continue                                                             1158
      if(nlp .le. maxit)goto 11621                                         1158
      jerr=-m                                                              1158
      return                                                               1158
11621 continue                                                             1159
10260 continue                                                             1159
      iz=1                                                                 1160
11630 continue                                                             1160
11631 continue                                                             1160
      nlp=nlp+1                                                            1160
      dlx=0.0                                                              1161
11640 do 11641 l=1,nin                                                     1161
      k=ia(l)                                                              1161
      jb=ix(k)                                                             1161
      je=ix(k+1)-1                                                         1162
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1163
      ak=a(k)                                                              1163
      u=gk+ak*xv(k)                                                        1163
      v=abs(u)-vp(k)*ab                                                    1163
      a(k)=0.0                                                             1164
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1165
      if(a(k).eq.ak)goto 11641                                             1166
      del=a(k)-ak                                                          1166
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1167
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1168
      o=o+del*xm(k)/xs(k)                                                  1168
      dlx=max(xv(k)*del**2,dlx)                                            1169
11641 continue                                                             1170
11642 continue                                                             1170
      if(dlx.lt.thr)goto 11632                                             1170
      if(nlp .le. maxit)goto 11661                                         1170
      jerr=-m                                                              1170
      return                                                               1170
11661 continue                                                             1171
      goto 11631                                                           1172
11632 continue                                                             1172
      jz=0                                                                 1173
      goto 11521                                                           1174
11522 continue                                                             1174
      if(nin .le. nx)goto 11681                                            1174
      jerr=-10000-m                                                        1174
      goto 11442                                                           1174
11681 continue                                                             1175
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1175
      kin(m)=nin                                                           1176
      rsqo(m)=rsq                                                          1176
      almo(m)=alm                                                          1176
      lmu=m                                                                1177
      if(m.lt.mnl)goto 11441                                               1177
      if(flmin.ge.1.0)goto 11441                                           1178
      me=0                                                                 1178
11690 do 11691 j=1,nin                                                     1178
      if(ao(j,m).ne.0.0) me=me+1                                           1178
11691 continue                                                             1178
11692 continue                                                             1178
      if(me.gt.ne)goto 11442                                               1179
      if(rsq-rsq0.lt.sml*rsq)goto 11442                                    1179
      if(rsq.gt.rsqmax)goto 11442                                          1180
11441 continue                                                             1181
11442 continue                                                             1181
      deallocate(a,mm,g,iy)                                                1182
      return                                                               1183
      end                                                                  1184
      subroutine spchkvars(no,ni,x,ix,ju)                                  1185
      real x(*)                                                            1185
      integer ix(*),ju(ni)                                                 1186
11700 do 11701 j=1,ni                                                      1186
      ju(j)=0                                                              1186
      jb=ix(j)                                                             1186
      nj=ix(j+1)-jb                                                        1186
      if(nj.eq.0)goto 11701                                                1187
      je=ix(j+1)-1                                                         1188
      if(nj .ge. no)goto 11721                                             1188
11730 do 11731 i=jb,je                                                     1188
      if(x(i).eq.0.0)goto 11731                                            1188
      ju(j)=1                                                              1188
      goto 11732                                                           1188
11731 continue                                                             1188
11732 continue                                                             1188
      goto 11741                                                           1189
11721 continue                                                             1189
      t=x(jb)                                                              1189
11750 do 11751 i=jb+1,je                                                   1189
      if(x(i).eq.t)goto 11751                                              1189
      ju(j)=1                                                              1189
      goto 11752                                                           1189
11751 continue                                                             1189
11752 continue                                                             1189
11741 continue                                                             1190
11711 continue                                                             1190
11701 continue                                                             1191
11702 continue                                                             1191
      return                                                               1192
      end                                                                  1193
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                         1194
      real ca(*),x(*),f(n)                                                 1194
      integer ia(*),ix(*),jx(*)                                            1195
      f=a0                                                                 1196
11760 do 11761 j=1,nin                                                     1196
      k=ia(j)                                                              1196
      kb=ix(k)                                                             1196
      ke=ix(k+1)-1                                                         1197
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                             1198
11761 continue                                                             1199
11762 continue                                                             1199
      return                                                               1200
      end                                                                  1201
      function row_prod(i,j,ia,ja,ra,w)                                    1202
      integer ia(*),ja(*)                                                  1202
      real ra(*),w(*)                                                      1203
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(   1205 
     *i),ia(j+1)-ia(j),w)
      return                                                               1206
      end                                                                  1207
      function dot(x,y,mx,my,nx,ny,w)                                      1208
      real x(*),y(*),w(*)                                                  1208
      integer mx(*),my(*)                                                  1209
      i=1                                                                  1209
      j=i                                                                  1209
      s=0.0                                                                1210
11770 continue                                                             1210
11771 continue                                                             1210
11780 continue                                                             1211
11781 if(mx(i).ge.my(j))goto 11782                                         1211
      i=i+1                                                                1211
      if(i.gt.nx) go to 11790                                              1211
      goto 11781                                                           1212
11782 continue                                                             1212
      if(mx(i).eq.my(j)) go to 11800                                       1213
11810 continue                                                             1213
11811 if(my(j).ge.mx(i))goto 11812                                         1213
      j=j+1                                                                1213
      if(j.gt.ny) go to 11790                                              1213
      goto 11811                                                           1214
11812 continue                                                             1214
      if(mx(i).eq.my(j)) go to 11800                                       1214
      goto 11771                                                           1215
11800 continue                                                             1215
      s=s+w(mx(i))*x(i)*y(j)                                               1216
      i=i+1                                                                1216
      if(i.gt.nx)goto 11772                                                1216
      j=j+1                                                                1216
      if(j.gt.ny)goto 11772                                                1217
      goto 11771                                                           1218
11772 continue                                                             1218
11790 continue                                                             1218
      dot=s                                                                1219
      return                                                               1220
      end                                                                  1221
      subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,ne,nx,nlam,flmin,ulam   1223 
     *,thr,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)             1224
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1225
      integer jd(*),ia(nx),nin(nlam)                                       1226
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 11831                                    1230
      jerr=10000                                                           1230
      return                                                               1230
11831 continue                                                             1231
      allocate(ww(1:no),stat=jerr)                                         1232
      allocate(ju(1:ni),stat=ierr)                                         1232
      jerr=jerr+ierr                                                       1233
      allocate(vq(1:ni),stat=ierr)                                         1233
      jerr=jerr+ierr                                                       1234
      allocate(xm(1:ni),stat=ierr)                                         1234
      jerr=jerr+ierr                                                       1235
      if(isd .le. 0)goto 11851                                             1235
      allocate(xs(1:ni),stat=ierr)                                         1235
      jerr=jerr+ierr                                                       1235
11851 continue                                                             1236
      if(jerr.ne.0) return                                                 1237
      call chkvars(no,ni,x,ju)                                             1238
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1239
      if(maxval(ju) .gt. 0)goto 11871                                      1239
      jerr=7777                                                            1239
      return                                                               1239
11871 continue                                                             1240
      vq=max(0.0,vp)                                                       1240
      vq=vq*ni/sum(vq)                                                     1241
11880 do 11881 i=1,no                                                      1241
      ww(i)=sum(y(i,:))                                                    1241
      y(i,:)=y(i,:)/ww(i)                                                  1241
11881 continue                                                             1241
11882 continue                                                             1241
      sw=sum(ww)                                                           1241
      ww=ww/sw                                                             1242
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             1243
      if(nc .ne. 1)goto 11901                                              1244
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,ne,nx,nlam,flmin   1246 
     *,ulam,  thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 11911                                                           1247
11901 continue                                                             1248
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,th   1250 
     *r,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
11911 continue                                                             1251
11891 continue                                                             1251
      if(jerr.gt.0) return                                                 1251
      dev0=2.0*sw*dev0                                                     1252
11920 do 11921 k=1,lmu                                                     1252
      nk=nin(k)                                                            1253
11930 do 11931 ic=1,nc                                                     1253
      if(isd .le. 0)goto 11951                                             1253
11960 do 11961 l=1,nk                                                      1253
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1253
11961 continue                                                             1253
11962 continue                                                             1253
11951 continue                                                             1254
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1255
11931 continue                                                             1256
11932 continue                                                             1256
11921 continue                                                             1257
11922 continue                                                             1257
      deallocate(ww,ju,vq,xm)                                              1257
      if(isd.gt.0) deallocate(xs)                                          1258
      return                                                               1259
      end                                                                  1260
      subroutine lstandard1 (no,ni,x,w,ju,isd,xm,xs)                       1261
      real x(no,ni),w(no),xm(ni),xs(ni)                                    1261
      integer ju(ni)                                                       1262
11970 do 11971 j=1,ni                                                      1262
      if(ju(j).eq.0)goto 11971                                             1263
      xm(j)=dot_product(w,x(:,j))                                          1263
      x(:,j)=x(:,j)-xm(j)                                                  1264
      if(isd .le. 0)goto 11991                                             1264
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 1264
      x(:,j)=x(:,j)/xs(j)                                                  1264
11991 continue                                                             1265
11971 continue                                                             1266
11972 continue                                                             1266
      return                                                               1267
      end                                                                  1268
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ulam   1270 
     *,shri,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1272 
     *5, devmax=0.999)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    1273
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1274
      integer ju(ni),m(nx),kin(nlam)                                       1275
      real, dimension (:), allocatable :: b,bs,v,r,xv,q,ga                      
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(b(0:ni),stat=jerr)                                          1280
      allocate(xv(1:ni),stat=ierr)                                         1280
      jerr=jerr+ierr                                                       1281
      allocate(ga(1:ni),stat=ierr)                                         1281
      jerr=jerr+ierr                                                       1282
      allocate(bs(0:ni),stat=ierr)                                         1282
      jerr=jerr+ierr                                                       1283
      allocate(mm(1:ni),stat=ierr)                                         1283
      jerr=jerr+ierr                                                       1284
      allocate(ixx(1:ni),stat=ierr)                                        1284
      jerr=jerr+ierr                                                       1285
      allocate(r(1:no),stat=ierr)                                          1285
      jerr=jerr+ierr                                                       1286
      allocate(v(1:no),stat=ierr)                                          1286
      jerr=jerr+ierr                                                       1287
      allocate(q(1:no),stat=ierr)                                          1287
      jerr=jerr+ierr                                                       1288
      if(jerr.ne.0) return                                                 1289
      fmax=log(1.0/pmin-1.0)                                               1289
      fmin=-fmax                                                           1289
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1290
      bta=parm                                                             1290
      omb=1.0-bta                                                          1291
      q0=dot_product(w,y)                                                  1291
      if(q0 .gt. pmin)goto 12011                                           1291
      jerr=8001                                                            1291
      return                                                               1291
12011 continue                                                             1292
      if(q0 .lt. 1.0-pmin)goto 12031                                       1292
      jerr=9001                                                            1292
      return                                                               1292
12031 continue                                                             1293
      ixx=0                                                                1293
      al=0.0                                                               1293
      bz=log(q0/(1.0-q0))                                                  1294
      if(nonzero(no,g) .ne. 0)goto 12051                                   1294
      vi=q0*(1.0-q0)                                                       1294
      b(0)=bz                                                              1294
      v=vi*w                                                               1295
      r=w*(y-q0)                                                           1295
      q=q0                                                                 1295
      xmz=vi                                                               1295
      dev1=-(bz*q0+log(1.0-q0))                                            1296
      goto 12061                                                           1297
12051 continue                                                             1297
      b(0)=azero(no,y,g,w,jerr)                                            1297
      if(jerr.ne.0) return                                                 1298
      q=1.0/(1.0+exp(-b(0)-g))                                             1298
      v=w*q*(1.0-q)                                                        1298
      r=w*(y-q)                                                            1298
      xmz=sum(v)                                                           1299
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1300
12061 continue                                                             1301
12041 continue                                                             1301
      if(kopt .le. 0)goto 12081                                            1302
      if(isd .le. 0)goto 12101                                             1302
      xv=0.25                                                              1302
      goto 12111                                                           1303
12101 continue                                                             1303
12120 do 12121 j=1,ni                                                      1303
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1303
12121 continue                                                             1303
12122 continue                                                             1303
12111 continue                                                             1304
12091 continue                                                             1304
12081 continue                                                             1305
      dev0=dev1                                                            1306
12130 do 12131 i=1,no                                                      1306
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1307
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1308
12131 continue                                                             1309
12132 continue                                                             1309
      if(flmin .ge. 1.0)goto 12151                                         1309
      eqs=max(eps,flmin)                                                   1309
      alf=eqs**(1.0/(nlam-1))                                              1309
12151 continue                                                             1310
      m=0                                                                  1310
      mm=0                                                                 1310
      nlp=0                                                                1310
      nin=nlp                                                              1310
      mnl=min(mnlam,nlam)                                                  1310
      bs=0.0                                                               1310
      b(1:ni)=0.0                                                          1311
      shr=shri*dev0                                                        1312
12160 do 12161 j=1,ni                                                      1312
      if(ju(j).eq.0)goto 12161                                             1312
      ga(j)=abs(dot_product(r,x(:,j)))                                     1312
12161 continue                                                             1313
12162 continue                                                             1313
12170 do 12171 ilm=1,nlam                                                  1313
      al0=al                                                               1314
      if(flmin .lt. 1.0)goto 12191                                         1314
      al=ulam(ilm)                                                         1314
      goto 12181                                                           1315
12191 if(ilm .le. 2)goto 12201                                             1315
      al=al*alf                                                            1315
      goto 12181                                                           1316
12201 if(ilm .ne. 1)goto 12211                                             1316
      al=big                                                               1316
      goto 12221                                                           1317
12211 continue                                                             1317
      al0=0.0                                                              1318
12230 do 12231 j=1,ni                                                      1318
      if(ju(j).eq.0)goto 12231                                             1318
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1318
12231 continue                                                             1319
12232 continue                                                             1319
      al0=al0/max(bta,1.0e-3)                                              1319
      al=alf*al0                                                           1320
12221 continue                                                             1321
12181 continue                                                             1321
      al2=al*omb                                                           1321
      al1=al*bta                                                           1321
      tlam=bta*(2.0*al-al0)                                                1322
12240 do 12241 k=1,ni                                                      1322
      if(ixx(k).eq.1)goto 12241                                            1322
      if(ju(k).eq.0)goto 12241                                             1323
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1324
12241 continue                                                             1325
12242 continue                                                             1325
10680 continue                                                             1326
12250 continue                                                             1326
12251 continue                                                             1326
      bs(0)=b(0)                                                           1326
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1327
      if(kopt .ne. 0)goto 12271                                            1328
12280 do 12281 j=1,ni                                                      1328
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                       1328
12281 continue                                                             1329
12282 continue                                                             1329
12271 continue                                                             1330
12290 continue                                                             1330
12291 continue                                                             1330
      nlp=nlp+1                                                            1330
      dlx=0.0                                                              1331
12300 do 12301 k=1,ni                                                      1331
      if(ixx(k).eq.0)goto 12301                                            1332
      bk=b(k)                                                              1332
      gk=dot_product(r,x(:,k))                                             1333
      u=gk+xv(k)*b(k)                                                      1333
      au=abs(u)-vp(k)*al1                                                  1334
      if(au .gt. 0.0)goto 12321                                            1334
      b(k)=0.0                                                             1334
      goto 12331                                                           1335
12321 continue                                                             1335
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1335
12331 continue                                                             1336
12311 continue                                                             1336
      d=b(k)-bk                                                            1336
      if(abs(d).le.0.0)goto 12301                                          1336
      dlx=max(dlx,xv(k)*d**2)                                              1337
      r=r-d*v*x(:,k)                                                       1338
      if(mm(k) .ne. 0)goto 12351                                           1338
      nin=nin+1                                                            1338
      if(nin.gt.nx)goto 12302                                              1339
      mm(k)=nin                                                            1339
      m(nin)=k                                                             1340
12351 continue                                                             1341
12301 continue                                                             1342
12302 continue                                                             1342
      if(nin.gt.nx)goto 12292                                              1343
      d=sum(r)/xmz                                                         1344
      if(d .eq. 0.0)goto 12371                                             1344
      b(0)=b(0)+d                                                          1344
      dlx=max(dlx,xmz*d**2)                                                1344
      r=r-d*v                                                              1344
12371 continue                                                             1345
      if(dlx.lt.shr)goto 12292                                             1345
      if(nlp .le. maxit)goto 12391                                         1345
      jerr=-ilm                                                            1345
      return                                                               1345
12391 continue                                                             1346
12400 continue                                                             1346
12401 continue                                                             1346
      nlp=nlp+1                                                            1346
      dlx=0.0                                                              1347
12410 do 12411 l=1,nin                                                     1347
      k=m(l)                                                               1347
      bk=b(k)                                                              1348
      gk=dot_product(r,x(:,k))                                             1349
      u=gk+xv(k)*b(k)                                                      1349
      au=abs(u)-vp(k)*al1                                                  1350
      if(au .gt. 0.0)goto 12431                                            1350
      b(k)=0.0                                                             1350
      goto 12441                                                           1351
12431 continue                                                             1351
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1351
12441 continue                                                             1352
12421 continue                                                             1352
      d=b(k)-bk                                                            1352
      if(abs(d).le.0.0)goto 12411                                          1352
      dlx=max(dlx,xv(k)*d**2)                                              1353
      r=r-d*v*x(:,k)                                                       1354
12411 continue                                                             1355
12412 continue                                                             1355
      d=sum(r)/xmz                                                         1356
      if(d .eq. 0.0)goto 12461                                             1356
      b(0)=b(0)+d                                                          1356
      dlx=max(dlx,xmz*d**2)                                                1356
      r=r-d*v                                                              1356
12461 continue                                                             1357
      if(dlx.lt.shr)goto 12402                                             1357
      if(nlp .le. maxit)goto 12481                                         1357
      jerr=-ilm                                                            1357
      return                                                               1357
12481 continue                                                             1358
      goto 12401                                                           1359
12402 continue                                                             1359
      goto 12291                                                           1360
12292 continue                                                             1360
      if(nin.gt.nx)goto 12252                                              1361
12490 do 12491 i=1,no                                                      1361
      fi=b(0)+g(i)                                                         1362
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1363
      if(fi .ge. fmin)goto 12511                                           1363
      q(i)=0.0                                                             1363
      goto 12501                                                           1363
12511 if(fi .le. fmax)goto 12521                                           1363
      q(i)=1.0                                                             1363
      goto 12531                                                           1364
12521 continue                                                             1364
      q(i)=1.0/(1.0+exp(-fi))                                              1364
12531 continue                                                             1365
12501 continue                                                             1365
12491 continue                                                             1366
12492 continue                                                             1366
      v=w*q*(1.0-q)                                                        1366
      xmz=sum(v)                                                           1366
      if(xmz.le.vmin)goto 12252                                            1366
      r=w*(y-q)                                                            1367
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 12551                           1367
      ix=0                                                                 1368
12560 do 12561 j=1,nin                                                     1368
      k=m(j)                                                               1369
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 12561                           1369
      ix=1                                                                 1369
      goto 12562                                                           1370
12561 continue                                                             1371
12562 continue                                                             1371
      if(ix .ne. 0)goto 12581                                              1372
12590 do 12591 k=1,ni                                                      1372
      if(ixx(k).eq.1)goto 12591                                            1372
      if(ju(k).eq.0)goto 12591                                             1373
      ga(k)=abs(dot_product(r,x(:,k)))                                     1374
      if(ga(k) .le. al1*vp(k))goto 12611                                   1374
      ixx(k)=1                                                             1374
      ix=1                                                                 1374
12611 continue                                                             1375
12591 continue                                                             1376
12592 continue                                                             1376
      if(ix.eq.1) go to 10680                                              1377
      goto 12252                                                           1378
12581 continue                                                             1379
12551 continue                                                             1380
      goto 12251                                                           1381
12252 continue                                                             1381
      if(nin .le. nx)goto 12631                                            1381
      jerr=-10000-ilm                                                      1381
      goto 12172                                                           1381
12631 continue                                                             1382
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1382
      kin(ilm)=nin                                                         1383
      a0(ilm)=b(0)                                                         1383
      alm(ilm)=al                                                          1383
      lmu=ilm                                                              1384
      devi=dev2(no,w,y,q,pmin)                                             1385
      dev(ilm)=(dev1-devi)/dev0                                            1385
      if(xmz.le.vmin)goto 12172                                            1386
      if(ilm.lt.mnl)goto 12171                                             1386
      if(flmin.ge.1.0)goto 12171                                           1387
      me=0                                                                 1387
12640 do 12641 j=1,nin                                                     1387
      if(a(j,ilm).ne.0.0) me=me+1                                          1387
12641 continue                                                             1387
12642 continue                                                             1387
      if(me.gt.ne)goto 12172                                               1388
      if(dev(ilm).gt.devmax)goto 12172                                     1388
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12172                             1389
12171 continue                                                             1390
12172 continue                                                             1390
      g=log(q/(1.0-q))                                                     1391
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                                  1392
      return                                                               1393
      end                                                                  1394
      function dev2(n,w,y,p,pmin)                                          1395
      real w(n),y(n),p(n)                                                  1396
      pmax=1.0-pmin                                                        1396
      s=0.0                                                                1397
12650 do 12651 i=1,n                                                       1397
      pi=min(max(pmin,p(i)),pmax)                                          1398
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1399
12651 continue                                                             1400
12652 continue                                                             1400
      dev2=s                                                               1401
      return                                                               1402
      end                                                                  1403
      function azero(n,y,g,q,jerr)                                         1404
      parameter(eps=1.0e-7)                                                1405
      real y(n),g(n),q(n)                                                  1406
      real, dimension (:), allocatable :: e,p,w                                 
      allocate(e(1:n),stat=jerr)                                           1410
      allocate(p(1:n),stat=ierr)                                           1410
      jerr=jerr+ierr                                                       1411
      allocate(w(1:n),stat=ierr)                                           1411
      jerr=jerr+ierr                                                       1412
      if(jerr.ne.0) return                                                 1413
      az=0.0                                                               1413
      e=exp(-g)                                                            1413
      qy=dot_product(q,y)                                                  1413
      p=1.0/(1.0+e)                                                        1414
12660 continue                                                             1414
12661 continue                                                             1414
      w=q*p*(1.0-p)                                                        1415
      d=(qy-dot_product(q,p))/sum(w)                                       1415
      az=az+d                                                              1415
      if(abs(d).lt.eps)goto 12662                                          1416
      ea0=exp(-az)                                                         1416
      p=1.0/(1.0+ea0*e)                                                    1417
      goto 12661                                                           1418
12662 continue                                                             1418
      azero=az                                                             1419
      deallocate(e,p,w)                                                    1420
      return                                                               1421
      end                                                                  1422
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ul   1424 
     *am,shri,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1426 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam)              1427
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1428
      integer ju(ni),m(nx),kin(nlam)                                       1429
      real, dimension (:,:), allocatable :: q                                   
      real, dimension (:), allocatable :: sxp,sxpl                              
      real, dimension (:), allocatable :: di,v,r,ga                             
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is,ixx                          
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no),stat=ierr)                                          1440
      jerr=jerr+ierr                                                       1441
      allocate(v(1:no),stat=ierr)                                          1441
      jerr=jerr+ierr                                                       1442
      allocate(mm(1:ni),stat=ierr)                                         1442
      jerr=jerr+ierr                                                       1443
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1443
      jerr=jerr+ierr                                                       1444
      allocate(sxp(1:no),stat=ierr)                                        1444
      jerr=jerr+ierr                                                       1445
      allocate(sxpl(1:no),stat=ierr)                                       1445
      jerr=jerr+ierr                                                       1446
      allocate(di(1:no),stat=ierr)                                         1446
      jerr=jerr+ierr                                                       1447
      allocate(ga(1:ni),stat=ierr)                                         1447
      jerr=jerr+ierr                                                       1448
      allocate(ixx(1:ni),stat=ierr)                                        1448
      jerr=jerr+ierr                                                       1449
      if(jerr.ne.0) return                                                 1450
      pmax=1.0-pmin                                                        1450
      emin=pmin/pmax                                                       1450
      emax=1.0/emin                                                        1451
      pfm=(1.0+pmin)*pmin                                                  1451
      pfx=(1.0-pmin)*pmax                                                  1451
      vmin=pfm*pmax                                                        1452
      bta=parm                                                             1452
      omb=1.0-bta                                                          1452
      dev1=0.0                                                             1452
      dev0=0.0                                                             1453
12670 do 12671 ic=1,nc                                                     1453
      q0=dot_product(w,y(:,ic))                                            1454
      if(q0 .gt. pmin)goto 12691                                           1454
      jerr =8000+ic                                                        1454
      return                                                               1454
12691 continue                                                             1455
      if(q0 .lt. 1.0-pmin)goto 12711                                       1455
      jerr =9000+ic                                                        1455
      return                                                               1455
12711 continue                                                             1456
      b(0,ic)=log(q0)                                                      1456
      dev1=dev1-q0*b(0,ic)                                                 1456
      b(1:ni,ic)=0.0                                                       1457
12671 continue                                                             1458
12672 continue                                                             1458
      ixx=0                                                                1458
      al=0.0                                                               1459
      if(nonzero(no*nc,g) .ne. 0)goto 12731                                1460
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1460
      sxp=0.0                                                              1461
12740 do 12741 ic=1,nc                                                     1461
      q(:,ic)=exp(b(0,ic))                                                 1461
      sxp=sxp+q(:,ic)                                                      1461
12741 continue                                                             1462
12742 continue                                                             1462
      goto 12751                                                           1463
12731 continue                                                             1463
12760 do 12761 i=1,no                                                      1463
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1463
12761 continue                                                             1463
12762 continue                                                             1463
      sxp=0.0                                                              1464
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1464
      if(jerr.ne.0) return                                                 1465
      dev1=0.0                                                             1466
12770 do 12771 ic=1,nc                                                     1466
      q(:,ic)=b(0,ic)+g(:,ic)                                              1467
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             1468
      q(:,ic)=exp(q(:,ic))                                                 1468
      sxp=sxp+q(:,ic)                                                      1469
12771 continue                                                             1470
12772 continue                                                             1470
      sxpl=w*log(sxp)                                                      1470
12780 do 12781 ic=1,nc                                                     1470
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  1470
12781 continue                                                             1471
12782 continue                                                             1471
12751 continue                                                             1472
12721 continue                                                             1472
12790 do 12791 ic=1,nc                                                     1472
12800 do 12801 i=1,no                                                      1472
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1472
12801 continue                                                             1472
12802 continue                                                             1472
12791 continue                                                             1473
12792 continue                                                             1473
      dev0=dev0+dev1                                                       1474
      if(kopt .le. 0)goto 12821                                            1475
      if(isd .le. 0)goto 12841                                             1475
      xv=0.25                                                              1475
      goto 12851                                                           1476
12841 continue                                                             1476
12860 do 12861 j=1,ni                                                      1476
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1476
12861 continue                                                             1476
12862 continue                                                             1476
12851 continue                                                             1477
12831 continue                                                             1477
12821 continue                                                             1478
      if(flmin .ge. 1.0)goto 12881                                         1478
      eqs=max(eps,flmin)                                                   1478
      alf=eqs**(1.0/(nlam-1))                                              1478
12881 continue                                                             1479
      m=0                                                                  1479
      mm=0                                                                 1479
      nin=0                                                                1479
      nlp=0                                                                1479
      mnl=min(mnlam,nlam)                                                  1479
      bs=0.0                                                               1479
      shr=shri*dev0                                                        1480
      ga=0.0                                                               1481
12890 do 12891 ic=1,nc                                                     1481
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1482
12900 do 12901 j=1,ni                                                      1482
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))           1482
12901 continue                                                             1483
12902 continue                                                             1483
12891 continue                                                             1484
12892 continue                                                             1484
12910 do 12911 ilm=1,nlam                                                  1484
      al0=al                                                               1485
      if(flmin .lt. 1.0)goto 12931                                         1485
      al=ulam(ilm)                                                         1485
      goto 12921                                                           1486
12931 if(ilm .le. 2)goto 12941                                             1486
      al=al*alf                                                            1486
      goto 12921                                                           1487
12941 if(ilm .ne. 1)goto 12951                                             1487
      al=big                                                               1487
      goto 12961                                                           1488
12951 continue                                                             1488
      al0=0.0                                                              1489
12970 do 12971 j=1,ni                                                      1489
      if(ju(j).eq.0)goto 12971                                             1489
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1489
12971 continue                                                             1490
12972 continue                                                             1490
      al0=al0/max(bta,1.0e-3)                                              1490
      al=alf*al0                                                           1491
12961 continue                                                             1492
12921 continue                                                             1492
      al2=al*omb                                                           1492
      al1=al*bta                                                           1492
      tlam=bta*(2.0*al-al0)                                                1493
12980 do 12981 k=1,ni                                                      1493
      if(ixx(k).eq.1)goto 12981                                            1493
      if(ju(k).eq.0)goto 12981                                             1494
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1495
12981 continue                                                             1496
12982 continue                                                             1496
10680 continue                                                             1497
12990 continue                                                             1497
12991 continue                                                             1497
      ix=0                                                                 1497
      jx=ix                                                                1497
      ig=0                                                                 1498
13000 do 13001 ic=1,nc                                                     1498
      bs(0,ic)=b(0,ic)                                                     1499
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1500
      xmz=0.0                                                              1501
13010 do 13011 i=1,no                                                      1501
      pic=q(i,ic)/sxp(i)                                                   1502
      if(pic .ge. pfm)goto 13031                                           1502
      pic=0.0                                                              1502
      v(i)=0.0                                                             1502
      goto 13021                                                           1503
13031 if(pic .le. pfx)goto 13041                                           1503
      pic=1.0                                                              1503
      v(i)=0.0                                                             1503
      goto 13051                                                           1504
13041 continue                                                             1504
      v(i)=w(i)*pic*(1.0-pic)                                              1504
      xmz=xmz+v(i)                                                         1504
13051 continue                                                             1505
13021 continue                                                             1505
      r(i)=w(i)*(y(i,ic)-pic)                                              1506
13011 continue                                                             1507
13012 continue                                                             1507
      if(xmz.le.vmin)goto 13001                                            1507
      ig=1                                                                 1508
      if(kopt .ne. 0)goto 13071                                            1509
13080 do 13081 j=1,ni                                                      1509
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                    1509
13081 continue                                                             1510
13082 continue                                                             1510
13071 continue                                                             1511
13090 continue                                                             1511
13091 continue                                                             1511
      nlp=nlp+1                                                            1511
      dlx=0.0                                                              1512
13100 do 13101 k=1,ni                                                      1512
      if(ixx(k).eq.0)goto 13101                                            1513
      bk=b(k,ic)                                                           1513
      gk=dot_product(r,x(:,k))                                             1514
      u=gk+xv(k,ic)*b(k,ic)                                                1514
      au=abs(u)-vp(k)*al1                                                  1515
      if(au .gt. 0.0)goto 13121                                            1515
      b(k,ic)=0.0                                                          1515
      goto 13131                                                           1516
13121 continue                                                             1516
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1516
13131 continue                                                             1517
13111 continue                                                             1517
      d=b(k,ic)-bk                                                         1517
      if(abs(d).le.0.0)goto 13101                                          1518
      dlx=max(dlx,xv(k,ic)*d**2)                                           1518
      r=r-d*v*x(:,k)                                                       1519
      if(mm(k) .ne. 0)goto 13151                                           1519
      nin=nin+1                                                            1520
      if(nin .le. nx)goto 13171                                            1520
      jx=1                                                                 1520
      goto 13102                                                           1520
13171 continue                                                             1521
      mm(k)=nin                                                            1521
      m(nin)=k                                                             1522
13151 continue                                                             1523
13101 continue                                                             1524
13102 continue                                                             1524
      if(jx.gt.0)goto 13092                                                1525
      d=sum(r)/xmz                                                         1526
      if(d .eq. 0.0)goto 13191                                             1526
      b(0,ic)=b(0,ic)+d                                                    1526
      dlx=max(dlx,xmz*d**2)                                                1526
      r=r-d*v                                                              1526
13191 continue                                                             1527
      if(dlx.lt.shr)goto 13092                                             1528
      if(nlp .le. maxit)goto 13211                                         1528
      jerr=-ilm                                                            1528
      return                                                               1528
13211 continue                                                             1529
13220 continue                                                             1529
13221 continue                                                             1529
      nlp=nlp+1                                                            1529
      dlx=0.0                                                              1530
13230 do 13231 l=1,nin                                                     1530
      k=m(l)                                                               1530
      bk=b(k,ic)                                                           1531
      gk=dot_product(r,x(:,k))                                             1532
      u=gk+xv(k,ic)*b(k,ic)                                                1532
      au=abs(u)-vp(k)*al1                                                  1533
      if(au .gt. 0.0)goto 13251                                            1533
      b(k,ic)=0.0                                                          1533
      goto 13261                                                           1534
13251 continue                                                             1534
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1534
13261 continue                                                             1535
13241 continue                                                             1535
      d=b(k,ic)-bk                                                         1535
      if(abs(d).le.0.0)goto 13231                                          1536
      dlx=max(dlx,xv(k,ic)*d**2)                                           1536
      r=r-d*v*x(:,k)                                                       1537
13231 continue                                                             1538
13232 continue                                                             1538
      d=sum(r)/xmz                                                         1539
      if(d .eq. 0.0)goto 13281                                             1539
      b(0,ic)=b(0,ic)+d                                                    1540
      dlx=max(dlx,xmz*d**2)                                                1540
      r=r-d*v                                                              1541
13281 continue                                                             1542
      if(dlx.lt.shr)goto 13222                                             1542
      if(nlp .le. maxit)goto 13301                                         1542
      jerr=-ilm                                                            1542
      return                                                               1542
13301 continue                                                             1543
      goto 13221                                                           1544
13222 continue                                                             1544
      goto 13091                                                           1545
13092 continue                                                             1545
      if(jx.gt.0)goto 13002                                                1546
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                            1547
      if(ix .ne. 0)goto 13321                                              1548
13330 do 13331 j=1,nin                                                     1548
      k=m(j)                                                               1549
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 13351                1549
      ix=1                                                                 1549
      goto 13332                                                           1549
13351 continue                                                             1550
13331 continue                                                             1551
13332 continue                                                             1551
13321 continue                                                             1552
13360 do 13361 i=1,no                                                      1552
      fi=b(0,ic)+g(i,ic)                                                   1554
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1555
      fi=min(max(exmn,fi),exmx)                                            1555
      sxp(i)=sxp(i)-q(i,ic)                                                1556
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    1557
      sxp(i)=sxp(i)+q(i,ic)                                                1558
13361 continue                                                             1559
13362 continue                                                             1559
13001 continue                                                             1560
13002 continue                                                             1560
      s=-sum(b(0,:))/nc                                                    1560
      b(0,:)=b(0,:)+s                                                      1560
      di=s                                                                 1561
13370 do 13371 j=1,nin                                                     1561
      l=m(j)                                                               1562
      if(vp(l) .gt. 0.0)goto 13391                                         1562
      s=sum(b(l,:))/nc                                                     1562
      goto 13401                                                           1563
13391 continue                                                             1563
      s=elc(parm,nc,b(l,:),is)                                             1563
13401 continue                                                             1564
13381 continue                                                             1564
      b(l,:)=b(l,:)-s                                                      1564
      di=di-s*x(:,l)                                                       1565
13371 continue                                                             1566
13372 continue                                                             1566
      di=exp(di)                                                           1566
      sxp=sxp*di                                                           1566
13410 do 13411 ic=1,nc                                                     1566
      q(:,ic)=q(:,ic)*di                                                   1566
13411 continue                                                             1567
13412 continue                                                             1567
      if(jx.gt.0)goto 12992                                                1567
      if(ig.eq.0)goto 12992                                                1568
      if(ix .ne. 0)goto 13431                                              1569
13440 do 13441 k=1,ni                                                      1569
      if(ixx(k).eq.1)goto 13441                                            1569
      if(ju(k).eq.0)goto 13441                                             1569
      ga(k)=0.0                                                            1569
13441 continue                                                             1570
13442 continue                                                             1570
13450 do 13451 ic=1,nc                                                     1570
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1571
13460 do 13461 k=1,ni                                                      1571
      if(ixx(k).eq.1)goto 13461                                            1571
      if(ju(k).eq.0)goto 13461                                             1572
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                          1573
13461 continue                                                             1574
13462 continue                                                             1574
13451 continue                                                             1575
13452 continue                                                             1575
13470 do 13471 k=1,ni                                                      1575
      if(ixx(k).eq.1)goto 13471                                            1575
      if(ju(k).eq.0)goto 13471                                             1576
      if(ga(k) .le. al1*vp(k))goto 13491                                   1576
      ixx(k)=1                                                             1576
      ix=1                                                                 1576
13491 continue                                                             1577
13471 continue                                                             1578
13472 continue                                                             1578
      if(ix.eq.1) go to 10680                                              1579
      goto 12992                                                           1580
13431 continue                                                             1581
      goto 12991                                                           1582
12992 continue                                                             1582
      if(jx .le. 0)goto 13511                                              1582
      jerr=-10000-ilm                                                      1582
      goto 12912                                                           1582
13511 continue                                                             1582
      devi=0.0                                                             1583
13520 do 13521 ic=1,nc                                                     1584
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1584
      a0(ic,ilm)=b(0,ic)                                                   1585
13530 do 13531 i=1,no                                                      1585
      if(y(i,ic).le.0.0)goto 13531                                         1586
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1587
13531 continue                                                             1588
13532 continue                                                             1588
13521 continue                                                             1589
13522 continue                                                             1589
      kin(ilm)=nin                                                         1589
      alm(ilm)=al                                                          1589
      lmu=ilm                                                              1590
      dev(ilm)=(dev1-devi)/dev0                                            1590
      if(ig.eq.0)goto 12912                                                1591
      if(ilm.lt.mnl)goto 12911                                             1591
      if(flmin.ge.1.0)goto 12911                                           1592
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 12912             1593
      if(dev(ilm).gt.devmax)goto 12912                                     1593
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12912                             1594
12911 continue                                                             1595
12912 continue                                                             1595
      g=log(q)                                                             1595
13540 do 13541 i=1,no                                                      1595
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1595
13541 continue                                                             1596
13542 continue                                                             1596
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                           1597
      return                                                               1598
      end                                                                  1599
      subroutine kazero(kk,n,y,g,q,az,jerr)                                1600
      parameter(eps=1.0e-7)                                                1601
      real y(n,kk),g(n,kk),q(n),az(kk)                                     1602
      real, dimension (:), allocatable :: s                                     
      real, dimension (:,:), allocatable :: e                                   
      allocate(e(1:n,1:kk),stat=jerr)                                           
      allocate(s(1:n),stat=ierr)                                           1607
      jerr=jerr+ierr                                                       1608
      if(jerr.ne.0) return                                                 1609
      az=0.0                                                               1609
      e=exp(g)                                                             1609
13550 do 13551 i=1,n                                                       1609
      s(i)=sum(e(i,:))                                                     1609
13551 continue                                                             1610
13552 continue                                                             1610
13560 continue                                                             1610
13561 continue                                                             1610
      dm=0.0                                                               1611
13570 do 13571 k=1,kk                                                      1611
      t=0.0                                                                1611
      u=t                                                                  1612
13580 do 13581 i=1,n                                                       1612
      pik=e(i,k)/s(i)                                                      1613
      t=t+q(i)*(y(i,k)-pik)                                                1613
      u=u+q(i)*pik*(1.0-pik)                                               1614
13581 continue                                                             1615
13582 continue                                                             1615
      d=t/u                                                                1615
      az(k)=az(k)+d                                                        1615
      ed=exp(d)                                                            1615
      dm=max(dm,abs(d))                                                    1616
13590 do 13591 i=1,n                                                       1616
      z=e(i,k)                                                             1616
      e(i,k)=z*ed                                                          1616
      s(i)=s(i)-z+e(i,k)                                                   1616
13591 continue                                                             1617
13592 continue                                                             1617
13571 continue                                                             1618
13572 continue                                                             1618
      if(dm.lt.eps)goto 13562                                              1618
      goto 13561                                                           1619
13562 continue                                                             1619
      az=az-sum(az)/kk                                                     1620
      deallocate(e,s)                                                      1621
      return                                                               1622
      end                                                                  1623
      function elc(parm,n,a,m)                                             1624
      real a(n)                                                            1624
      integer m(n)                                                         1625
      fn=n                                                                 1625
      am=sum(a)/fn                                                         1626
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 13611                       1626
      elc=am                                                               1626
      return                                                               1626
13611 continue                                                             1627
13620 do 13621 i=1,n                                                       1627
      m(i)=i                                                               1627
13621 continue                                                             1627
13622 continue                                                             1627
      call psort7(a,m,1,n)                                                 1628
      if(a(m(1)) .ne. a(m(n)))goto 13641                                   1628
      elc=a(1)                                                             1628
      return                                                               1628
13641 continue                                                             1629
      if(mod(n,2) .ne. 1)goto 13661                                        1629
      ad=a(m(n/2+1))                                                       1629
      goto 13671                                                           1630
13661 continue                                                             1630
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       1630
13671 continue                                                             1631
13651 continue                                                             1631
      if(parm .ne. 1.0)goto 13691                                          1631
      elc=ad                                                               1631
      return                                                               1631
13691 continue                                                             1632
      b1=min(am,ad)                                                        1632
      b2=max(am,ad)                                                        1632
      k2=1                                                                 1633
13700 continue                                                             1633
13701 if(a(m(k2)).gt.b1)goto 13702                                         1633
      k2=k2+1                                                              1633
      goto 13701                                                           1633
13702 continue                                                             1633
      k1=k2-1                                                              1634
13710 continue                                                             1634
13711 if(a(m(k2)).ge.b2)goto 13712                                         1634
      k2=k2+1                                                              1634
      goto 13711                                                           1635
13712 continue                                                             1635
      r=parm/((1.0-parm)*fn)                                               1635
      is=0                                                                 1635
      sm=n-2*(k1-1)                                                        1636
13720 do 13721 k=k1,k2-1                                                   1636
      sm=sm-2.0                                                            1636
      s=r*sm+am                                                            1637
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 13741                   1637
      is=k                                                                 1637
      goto 13722                                                           1637
13741 continue                                                             1638
13721 continue                                                             1639
13722 continue                                                             1639
      if(is .eq. 0)goto 13761                                              1639
      elc=s                                                                1639
      return                                                               1639
13761 continue                                                             1639
      r2=2.0*r                                                             1639
      s1=a(m(k1))                                                          1639
      am2=2.0*am                                                           1640
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    1640
      elc=s1                                                               1641
13770 do 13771 k=k1+1,k2                                                   1641
      s=a(m(k))                                                            1641
      if(s.eq.s1)goto 13771                                                1642
      c=r2*sum(abs(a-s))+s*(s-am2)                                         1643
      if(c .ge. cri)goto 13791                                             1643
      cri=c                                                                1643
      elc=s                                                                1643
13791 continue                                                             1643
      s1=s                                                                 1644
13771 continue                                                             1645
13772 continue                                                             1645
      return                                                               1646
      end                                                                  1647
      function nintot(ni,nx,nc,a,m,nin,is)                                 1648
      real a(nx,nc)                                                        1648
      integer m(nx),is(ni)                                                 1649
      is=0                                                                 1649
      nintot=0                                                             1650
13800 do 13801 ic=1,nc                                                     1650
13810 do 13811 j=1,nin                                                     1650
      k=m(j)                                                               1650
      if(is(k).ne.0)goto 13811                                             1651
      if(a(j,ic).eq.0.0)goto 13811                                         1651
      is(k)=k                                                              1651
      nintot=nintot+1                                                      1652
13811 continue                                                             1652
13812 continue                                                             1652
13801 continue                                                             1653
13802 continue                                                             1653
      return                                                               1654
      end                                                                  1655
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             1656
      real ca(nx,nc),a(ni,nc)                                              1656
      integer ia(nx)                                                       1657
      a=0.0                                                                1658
13820 do 13821 ic=1,nc                                                     1658
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            1658
13821 continue                                                             1659
13822 continue                                                             1659
      return                                                               1660
      end                                                                  1661
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      1662
      real a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                             1662
      integer ia(nx)                                                       1663
13830 do 13831 i=1,nt                                                      1663
13840 do 13841 ic=1,nc                                                     1663
      ans(ic,i)=a0(ic)                                                     1665
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   1666 
     *:nin)))
13841 continue                                                             1666
13842 continue                                                             1666
13831 continue                                                             1667
13832 continue                                                             1667
      return                                                               1668
      end                                                                  1669
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,ne,nx,nlam,fl   1671 
     *min,  ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      real x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)                 1672
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1673
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1674
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 13861                                    1678
      jerr=10000                                                           1678
      return                                                               1678
13861 continue                                                             1679
      allocate(ww(1:no),stat=jerr)                                         1680
      allocate(ju(1:ni),stat=ierr)                                         1680
      jerr=jerr+ierr                                                       1681
      allocate(vq(1:ni),stat=ierr)                                         1681
      jerr=jerr+ierr                                                       1682
      allocate(xm(1:ni),stat=ierr)                                         1682
      jerr=jerr+ierr                                                       1683
      allocate(xs(1:ni),stat=ierr)                                         1683
      jerr=jerr+ierr                                                       1684
      if(jerr.ne.0) return                                                 1685
      call spchkvars(no,ni,x,ix,ju)                                        1686
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1687
      if(maxval(ju) .gt. 0)goto 13881                                      1687
      jerr=7777                                                            1687
      return                                                               1687
13881 continue                                                             1688
      vq=max(0.0,vp)                                                       1688
      vq=vq*ni/sum(vq)                                                     1689
13890 do 13891 i=1,no                                                      1689
      ww(i)=sum(y(i,:))                                                    1689
      y(i,:)=y(i,:)/ww(i)                                                  1689
13891 continue                                                             1689
13892 continue                                                             1689
      sw=sum(ww)                                                           1689
      ww=ww/sw                                                             1690
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     1691
      if(nc .ne. 1)goto 13911                                              1692
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,ne,nx,n   1694 
     *lam,flmin,  ulam,thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,de
     *v,alm,nlp,jerr)
      goto 13921                                                           1695
13911 continue                                                             1696
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmi   1698 
     *n,ulam,  thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nl
     *p,jerr)
13921 continue                                                             1699
13901 continue                                                             1699
      if(jerr.gt.0) return                                                 1699
      dev0=2.0*sw*dev0                                                     1700
13930 do 13931 k=1,lmu                                                     1700
      nk=nin(k)                                                            1701
13940 do 13941 ic=1,nc                                                     1701
      if(isd .le. 0)goto 13961                                             1701
13970 do 13971 l=1,nk                                                      1701
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1701
13971 continue                                                             1701
13972 continue                                                             1701
13961 continue                                                             1702
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1703
13941 continue                                                             1704
13942 continue                                                             1704
13931 continue                                                             1705
13932 continue                                                             1705
      deallocate(ww,ju,vq,xm,xs)                                           1706
      return                                                               1707
      end                                                                  1708
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,xm,xs)                1709
      real x(*),w(no),xm(ni),xs(ni)                                        1709
      integer ix(*),jx(*),ju(ni)                                           1710
13980 do 13981 j=1,ni                                                      1710
      if(ju(j).eq.0)goto 13981                                             1710
      jb=ix(j)                                                             1710
      je=ix(j+1)-1                                                         1711
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1712
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   1713 
     *)**2)
13981 continue                                                             1714
13982 continue                                                             1714
      if(isd.eq.0) xs=1.0                                                  1715
      return                                                               1716
      end                                                                  1717
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam,   1719 
     *  flmin,ulam,shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,alm
     *,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1721 
     *5, devmax=0.999)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        1722
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1723
      real xb(ni),xs(ni)                                                   1723
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1724
      real, dimension (:), allocatable :: xm,b,bs,v,r,sc,xv,q,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(b(0:ni),stat=jerr)                                          1729
      allocate(xm(0:ni),stat=ierr)                                         1729
      jerr=jerr+ierr                                                       1730
      allocate(xv(1:ni),stat=ierr)                                         1730
      jerr=jerr+ierr                                                       1731
      allocate(bs(0:ni),stat=ierr)                                         1731
      jerr=jerr+ierr                                                       1732
      allocate(ga(1:ni),stat=ierr)                                         1732
      jerr=jerr+ierr                                                       1733
      allocate(mm(1:ni),stat=ierr)                                         1733
      jerr=jerr+ierr                                                       1734
      allocate(ixx(1:ni),stat=ierr)                                        1734
      jerr=jerr+ierr                                                       1735
      allocate(q(1:no),stat=ierr)                                          1735
      jerr=jerr+ierr                                                       1736
      allocate(r(1:no),stat=ierr)                                          1736
      jerr=jerr+ierr                                                       1737
      allocate(v(1:no),stat=ierr)                                          1737
      jerr=jerr+ierr                                                       1738
      allocate(sc(1:no),stat=ierr)                                         1738
      jerr=jerr+ierr                                                       1739
      if(jerr.ne.0) return                                                 1740
      fmax=log(1.0/pmin-1.0)                                               1740
      fmin=-fmax                                                           1740
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1741
      bta=parm                                                             1741
      omb=1.0-bta                                                          1742
      q0=dot_product(w,y)                                                  1742
      if(q0 .gt. pmin)goto 14001                                           1742
      jerr=8001                                                            1742
      return                                                               1742
14001 continue                                                             1743
      if(q0 .lt. 1.0-pmin)goto 14021                                       1743
      jerr=9001                                                            1743
      return                                                               1743
14021 continue                                                             1743
      bz=log(q0/(1.0-q0))                                                  1744
      if(nonzero(no,g) .ne. 0)goto 14041                                   1744
      vi=q0*(1.0-q0)                                                       1744
      b(0)=bz                                                              1744
      v=vi*w                                                               1745
      r=w*(y-q0)                                                           1745
      q=q0                                                                 1745
      xm(0)=vi                                                             1745
      dev1=-(bz*q0+log(1.0-q0))                                            1746
      goto 14051                                                           1747
14041 continue                                                             1747
      b(0)=azero(no,y,g,w,jerr)                                            1747
      if(jerr.ne.0) return                                                 1748
      q=1.0/(1.0+exp(-b(0)-g))                                             1748
      v=w*q*(1.0-q)                                                        1748
      r=w*(y-q)                                                            1748
      xm(0)=sum(v)                                                         1749
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1750
14051 continue                                                             1751
14031 continue                                                             1751
      if(kopt .le. 0)goto 14071                                            1752
      if(isd .le. 0)goto 14091                                             1752
      xv=0.25                                                              1752
      goto 14101                                                           1753
14091 continue                                                             1754
14110 do 14111 j=1,ni                                                      1754
      if(ju(j).eq.0)goto 14111                                             1754
      jb=ix(j)                                                             1754
      je=ix(j+1)-1                                                         1755
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          1756
14111 continue                                                             1757
14112 continue                                                             1757
14101 continue                                                             1758
14081 continue                                                             1758
14071 continue                                                             1759
      b(1:ni)=0.0                                                          1759
      dev0=dev1                                                            1760
14120 do 14121 i=1,no                                                      1760
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1761
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1762
14121 continue                                                             1763
14122 continue                                                             1763
      if(flmin .ge. 1.0)goto 14141                                         1763
      eqs=max(eps,flmin)                                                   1763
      alf=eqs**(1.0/(nlam-1))                                              1763
14141 continue                                                             1764
      m=0                                                                  1764
      mm=0                                                                 1764
      nin=0                                                                1764
      o=0.0                                                                1764
      svr=o                                                                1764
      mnl=min(mnlam,nlam)                                                  1764
      bs=0.0                                                               1764
      nlp=0                                                                1764
      nin=nlp                                                              1765
      shr=shri*dev0                                                        1765
      al=0.0                                                               1765
      ixx=0                                                                1766
14150 do 14151 j=1,ni                                                      1766
      if(ju(j).eq.0)goto 14151                                             1767
      jb=ix(j)                                                             1767
      je=ix(j+1)-1                                                         1767
      jn=ix(j+1)-ix(j)                                                     1768
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1769
      gj=dot_product(sc(1:jn),x(jb:je))                                    1770
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      1771
14151 continue                                                             1772
14152 continue                                                             1772
14160 do 14161 ilm=1,nlam                                                  1772
      al0=al                                                               1773
      if(flmin .lt. 1.0)goto 14181                                         1773
      al=ulam(ilm)                                                         1773
      goto 14171                                                           1774
14181 if(ilm .le. 2)goto 14191                                             1774
      al=al*alf                                                            1774
      goto 14171                                                           1775
14191 if(ilm .ne. 1)goto 14201                                             1775
      al=big                                                               1775
      goto 14211                                                           1776
14201 continue                                                             1776
      al0=0.0                                                              1777
14220 do 14221 j=1,ni                                                      1777
      if(ju(j).eq.0)goto 14221                                             1777
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1777
14221 continue                                                             1778
14222 continue                                                             1778
      al0=al0/max(bta,1.0e-3)                                              1778
      al=alf*al0                                                           1779
14211 continue                                                             1780
14171 continue                                                             1780
      al2=al*omb                                                           1780
      al1=al*bta                                                           1780
      tlam=bta*(2.0*al-al0)                                                1781
14230 do 14231 k=1,ni                                                      1781
      if(ixx(k).eq.1)goto 14231                                            1781
      if(ju(k).eq.0)goto 14231                                             1782
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1783
14231 continue                                                             1784
14232 continue                                                             1784
10680 continue                                                             1785
14240 continue                                                             1785
14241 continue                                                             1785
      bs(0)=b(0)                                                           1785
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1786
14250 do 14251 j=1,ni                                                      1786
      if(ixx(j).eq.0)goto 14251                                            1787
      jb=ix(j)                                                             1787
      je=ix(j+1)-1                                                         1787
      jn=ix(j+1)-ix(j)                                                     1788
      sc(1:jn)=v(jx(jb:je))                                                1789
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 1790
      if(kopt .ne. 0)goto 14271                                            1791
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              1792
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                1793
14271 continue                                                             1794
14251 continue                                                             1795
14252 continue                                                             1795
14280 continue                                                             1795
14281 continue                                                             1795
      nlp=nlp+1                                                            1795
      dlx=0.0                                                              1796
14290 do 14291 k=1,ni                                                      1796
      if(ixx(k).eq.0)goto 14291                                            1797
      jb=ix(k)                                                             1797
      je=ix(k+1)-1                                                         1797
      jn=ix(k+1)-ix(k)                                                     1797
      bk=b(k)                                                              1798
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1799
      gk=dot_product(sc(1:jn),x(jb:je))                                    1800
      gk=(gk-svr*xb(k))/xs(k)                                              1801
      u=gk+xv(k)*b(k)                                                      1801
      au=abs(u)-vp(k)*al1                                                  1802
      if(au .gt. 0.0)goto 14311                                            1802
      b(k)=0.0                                                             1802
      goto 14321                                                           1803
14311 continue                                                             1803
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1803
14321 continue                                                             1804
14301 continue                                                             1804
      d=b(k)-bk                                                            1804
      if(abs(d).le.0.0)goto 14291                                          1804
      dlx=max(dlx,xv(k)*d**2)                                              1805
      if(mm(k) .ne. 0)goto 14341                                           1805
      nin=nin+1                                                            1805
      if(nin.gt.nx)goto 14292                                              1806
      mm(k)=nin                                                            1806
      m(nin)=k                                                             1806
      sc(1:jn)=v(jx(jb:je))                                                1807
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 1808
14341 continue                                                             1809
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1810
      o=o+d*(xb(k)/xs(k))                                                  1811
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1812
14291 continue                                                             1813
14292 continue                                                             1813
      if(nin.gt.nx)goto 14282                                              1814
      d=svr/xm(0)                                                          1815
      if(d .eq. 0.0)goto 14361                                             1815
      b(0)=b(0)+d                                                          1815
      dlx=max(dlx,xm(0)*d**2)                                              1815
      r=r-d*v                                                              1815
14361 continue                                                             1816
      svr=svr-d*xm(0)                                                      1816
      if(dlx.lt.shr)goto 14282                                             1817
      if(nlp .le. maxit)goto 14381                                         1817
      jerr=-ilm                                                            1817
      return                                                               1817
14381 continue                                                             1818
14390 continue                                                             1818
14391 continue                                                             1818
      nlp=nlp+1                                                            1818
      dlx=0.0                                                              1819
14400 do 14401 l=1,nin                                                     1819
      k=m(l)                                                               1819
      jb=ix(k)                                                             1819
      je=ix(k+1)-1                                                         1820
      jn=ix(k+1)-ix(k)                                                     1820
      bk=b(k)                                                              1821
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1822
      gk=dot_product(sc(1:jn),x(jb:je))                                    1823
      gk=(gk-svr*xb(k))/xs(k)                                              1824
      u=gk+xv(k)*b(k)                                                      1824
      au=abs(u)-vp(k)*al1                                                  1825
      if(au .gt. 0.0)goto 14421                                            1825
      b(k)=0.0                                                             1825
      goto 14431                                                           1826
14421 continue                                                             1826
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1826
14431 continue                                                             1827
14411 continue                                                             1827
      d=b(k)-bk                                                            1827
      if(abs(d).le.0.0)goto 14401                                          1827
      dlx=max(dlx,xv(k)*d**2)                                              1828
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1829
      o=o+d*(xb(k)/xs(k))                                                  1830
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1831
14401 continue                                                             1832
14402 continue                                                             1832
      d=svr/xm(0)                                                          1833
      if(d .eq. 0.0)goto 14451                                             1833
      b(0)=b(0)+d                                                          1833
      dlx=max(dlx,xm(0)*d**2)                                              1833
      r=r-d*v                                                              1833
14451 continue                                                             1834
      svr=svr-d*xm(0)                                                      1834
      if(dlx.lt.shr)goto 14392                                             1835
      if(nlp .le. maxit)goto 14471                                         1835
      jerr=-ilm                                                            1835
      return                                                               1835
14471 continue                                                             1836
      goto 14391                                                           1837
14392 continue                                                             1837
      goto 14281                                                           1838
14282 continue                                                             1838
      if(nin.gt.nx)goto 14242                                              1839
      sc=b(0)                                                              1839
      b0=0.0                                                               1840
14480 do 14481 j=1,nin                                                     1840
      l=m(j)                                                               1840
      jb=ix(l)                                                             1840
      je=ix(l+1)-1                                                         1841
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      1842
      b0=b0-b(l)*xb(l)/xs(l)                                               1843
14481 continue                                                             1844
14482 continue                                                             1844
      sc=sc+b0                                                             1845
14490 do 14491 i=1,no                                                      1845
      fi=sc(i)+g(i)                                                        1846
      if(fi .ge. fmin)goto 14511                                           1846
      q(i)=0.0                                                             1846
      goto 14501                                                           1846
14511 if(fi .le. fmax)goto 14521                                           1846
      q(i)=1.0                                                             1846
      goto 14531                                                           1847
14521 continue                                                             1847
      q(i)=1.0/(1.0+exp(-fi))                                              1847
14531 continue                                                             1848
14501 continue                                                             1848
14491 continue                                                             1849
14492 continue                                                             1849
      v=w*q*(1.0-q)                                                        1849
      xm(0)=sum(v)                                                         1849
      if(xm(0).lt.vmin)goto 14242                                          1850
      r=w*(y-q)                                                            1850
      svr=sum(r)                                                           1850
      o=0.0                                                                1851
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 14551                         1851
      kx=0                                                                 1852
14560 do 14561 j=1,nin                                                     1852
      k=m(j)                                                               1853
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 14561                           1853
      kx=1                                                                 1853
      goto 14562                                                           1854
14561 continue                                                             1855
14562 continue                                                             1855
      if(kx .ne. 0)goto 14581                                              1856
14590 do 14591 j=1,ni                                                      1856
      if(ixx(j).eq.1)goto 14591                                            1856
      if(ju(j).eq.0)goto 14591                                             1857
      jb=ix(j)                                                             1857
      je=ix(j+1)-1                                                         1857
      jn=ix(j+1)-ix(j)                                                     1858
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1859
      gj=dot_product(sc(1:jn),x(jb:je))                                    1860
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      1861
      if(ga(j) .le. al1*vp(j))goto 14611                                   1861
      ixx(j)=1                                                             1861
      kx=1                                                                 1861
14611 continue                                                             1862
14591 continue                                                             1863
14592 continue                                                             1863
      if(kx.eq.1) go to 10680                                              1864
      goto 14242                                                           1865
14581 continue                                                             1866
14551 continue                                                             1867
      goto 14241                                                           1868
14242 continue                                                             1868
      if(nin .le. nx)goto 14631                                            1868
      jerr=-10000-ilm                                                      1868
      goto 14162                                                           1868
14631 continue                                                             1869
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1869
      kin(ilm)=nin                                                         1870
      a0(ilm)=b(0)                                                         1870
      alm(ilm)=al                                                          1870
      lmu=ilm                                                              1871
      devi=dev2(no,w,y,q,pmin)                                             1872
      dev(ilm)=(dev1-devi)/dev0                                            1873
      if(ilm.lt.mnl)goto 14161                                             1873
      if(flmin.ge.1.0)goto 14161                                           1874
      me=0                                                                 1874
14640 do 14641 j=1,nin                                                     1874
      if(a(j,ilm).ne.0.0) me=me+1                                          1874
14641 continue                                                             1874
14642 continue                                                             1874
      if(me.gt.ne)goto 14162                                               1875
      if(dev(ilm).gt.devmax)goto 14162                                     1875
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 14162                             1876
      if(xm(0).lt.vmin)goto 14162                                          1877
14161 continue                                                             1878
14162 continue                                                             1878
      g=log(q/(1.0-q))                                                     1879
      deallocate(xm,b,bs,v,r,sc,xv,q,mm,ga,ixx)                            1880
      return                                                               1881
      end                                                                  1882
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam   1884 
     *,flmin,ulam,  shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,al
     *m,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1886 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni)    1887
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1888
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1889
      real, dimension (:,:), allocatable :: q                                   
      real, dimension (:), allocatable :: sxp,sxpl                              
      real, dimension (:), allocatable :: sc,xm,v,r,ga                          
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is,iy                           
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(xm(0:ni),stat=ierr)                                         1900
      jerr=jerr+ierr                                                       1901
      allocate(r(1:no),stat=ierr)                                          1901
      jerr=jerr+ierr                                                       1902
      allocate(v(1:no),stat=ierr)                                          1902
      jerr=jerr+ierr                                                       1903
      allocate(mm(1:ni),stat=ierr)                                         1903
      jerr=jerr+ierr                                                       1904
      allocate(ga(1:ni),stat=ierr)                                         1904
      jerr=jerr+ierr                                                       1905
      allocate(iy(1:ni),stat=ierr)                                         1905
      jerr=jerr+ierr                                                       1906
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1906
      jerr=jerr+ierr                                                       1907
      allocate(sxp(1:no),stat=ierr)                                        1907
      jerr=jerr+ierr                                                       1908
      allocate(sxpl(1:no),stat=ierr)                                       1908
      jerr=jerr+ierr                                                       1909
      allocate(sc(1:no),stat=ierr)                                         1909
      jerr=jerr+ierr                                                       1910
      if(jerr.ne.0) return                                                 1911
      pmax=1.0-pmin                                                        1911
      emin=pmin/pmax                                                       1911
      emax=1.0/emin                                                        1912
      pfm=(1.0+pmin)*pmin                                                  1912
      pfx=(1.0-pmin)*pmax                                                  1912
      vmin=pfm*pmax                                                        1913
      bta=parm                                                             1913
      omb=1.0-bta                                                          1913
      dev1=0.0                                                             1913
      dev0=0.0                                                             1914
14650 do 14651 ic=1,nc                                                     1914
      q0=dot_product(w,y(:,ic))                                            1915
      if(q0 .gt. pmin)goto 14671                                           1915
      jerr =8000+ic                                                        1915
      return                                                               1915
14671 continue                                                             1916
      if(q0 .lt. 1.0-pmin)goto 14691                                       1916
      jerr =9000+ic                                                        1916
      return                                                               1916
14691 continue                                                             1917
      b(1:ni,ic)=0.0                                                       1917
      b(0,ic)=log(q0)                                                      1917
      dev1=dev1-q0*b(0,ic)                                                 1918
14651 continue                                                             1919
14652 continue                                                             1919
      iy=0                                                                 1919
      al=0.0                                                               1920
      if(nonzero(no*nc,g) .ne. 0)goto 14711                                1921
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1921
      sxp=0.0                                                              1922
14720 do 14721 ic=1,nc                                                     1922
      q(:,ic)=exp(b(0,ic))                                                 1922
      sxp=sxp+q(:,ic)                                                      1922
14721 continue                                                             1923
14722 continue                                                             1923
      goto 14731                                                           1924
14711 continue                                                             1924
14740 do 14741 i=1,no                                                      1924
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1924
14741 continue                                                             1924
14742 continue                                                             1924
      sxp=0.0                                                              1925
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1925
      if(jerr.ne.0) return                                                 1926
      dev1=0.0                                                             1927
14750 do 14751 ic=1,nc                                                     1927
      q(:,ic)=b(0,ic)+g(:,ic)                                              1928
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             1929
      q(:,ic)=exp(q(:,ic))                                                 1929
      sxp=sxp+q(:,ic)                                                      1930
14751 continue                                                             1931
14752 continue                                                             1931
      sxpl=w*log(sxp)                                                      1931
14760 do 14761 ic=1,nc                                                     1931
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  1931
14761 continue                                                             1932
14762 continue                                                             1932
14731 continue                                                             1933
14701 continue                                                             1933
14770 do 14771 ic=1,nc                                                     1933
14780 do 14781 i=1,no                                                      1933
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1933
14781 continue                                                             1933
14782 continue                                                             1933
14771 continue                                                             1934
14772 continue                                                             1934
      dev0=dev0+dev1                                                       1935
      if(kopt .le. 0)goto 14801                                            1936
      if(isd .le. 0)goto 14821                                             1936
      xv=0.25                                                              1936
      goto 14831                                                           1937
14821 continue                                                             1938
14840 do 14841 j=1,ni                                                      1938
      if(ju(j).eq.0)goto 14841                                             1938
      jb=ix(j)                                                             1938
      je=ix(j+1)-1                                                         1939
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        1940
14841 continue                                                             1941
14842 continue                                                             1941
14831 continue                                                             1942
14811 continue                                                             1942
14801 continue                                                             1943
      if(flmin .ge. 1.0)goto 14861                                         1943
      eqs=max(eps,flmin)                                                   1943
      alf=eqs**(1.0/(nlam-1))                                              1943
14861 continue                                                             1944
      m=0                                                                  1944
      mm=0                                                                 1944
      nin=0                                                                1944
      nlp=0                                                                1944
      mnl=min(mnlam,nlam)                                                  1944
      bs=0.0                                                               1944
      svr=0.0                                                              1944
      o=0.0                                                                1945
      shr=shri*dev0                                                        1945
      ga=0.0                                                               1946
14870 do 14871 ic=1,nc                                                     1946
      v=q(:,ic)/sxp                                                        1946
      r=w*(y(:,ic)-v)                                                      1946
      v=w*v*(1.0-v)                                                        1947
14880 do 14881 j=1,ni                                                      1947
      if(ju(j).eq.0)goto 14881                                             1948
      jb=ix(j)                                                             1948
      je=ix(j+1)-1                                                         1948
      jn=ix(j+1)-ix(j)                                                     1949
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1950
      gj=dot_product(sc(1:jn),x(jb:je))                                    1951
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             1952
14881 continue                                                             1953
14882 continue                                                             1953
14871 continue                                                             1954
14872 continue                                                             1954
14890 do 14891 ilm=1,nlam                                                  1954
      al0=al                                                               1955
      if(flmin .lt. 1.0)goto 14911                                         1955
      al=ulam(ilm)                                                         1955
      goto 14901                                                           1956
14911 if(ilm .le. 2)goto 14921                                             1956
      al=al*alf                                                            1956
      goto 14901                                                           1957
14921 if(ilm .ne. 1)goto 14931                                             1957
      al=big                                                               1957
      goto 14941                                                           1958
14931 continue                                                             1958
      al0=0.0                                                              1959
14950 do 14951 j=1,ni                                                      1959
      if(ju(j).eq.0)goto 14951                                             1959
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1959
14951 continue                                                             1960
14952 continue                                                             1960
      al0=al0/max(bta,1.0e-3)                                              1960
      al=alf*al0                                                           1961
14941 continue                                                             1962
14901 continue                                                             1962
      al2=al*omb                                                           1962
      al1=al*bta                                                           1962
      tlam=bta*(2.0*al-al0)                                                1963
14960 do 14961 k=1,ni                                                      1963
      if(iy(k).eq.1)goto 14961                                             1963
      if(ju(k).eq.0)goto 14961                                             1964
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      1965
14961 continue                                                             1966
14962 continue                                                             1966
10680 continue                                                             1967
14970 continue                                                             1967
14971 continue                                                             1967
      ixx=0                                                                1967
      jxx=ixx                                                              1967
      ig=0                                                                 1968
14980 do 14981 ic=1,nc                                                     1968
      bs(0,ic)=b(0,ic)                                                     1969
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1970
      xm(0)=0.0                                                            1970
      svr=0.0                                                              1970
      o=0.0                                                                1971
14990 do 14991 i=1,no                                                      1971
      pic=q(i,ic)/sxp(i)                                                   1972
      if(pic .ge. pfm)goto 15011                                           1972
      pic=0.0                                                              1972
      v(i)=0.0                                                             1972
      goto 15001                                                           1973
15011 if(pic .le. pfx)goto 15021                                           1973
      pic=1.0                                                              1973
      v(i)=0.0                                                             1973
      goto 15031                                                           1974
15021 continue                                                             1974
      v(i)=w(i)*pic*(1.0-pic)                                              1974
      xm(0)=xm(0)+v(i)                                                     1974
15031 continue                                                             1975
15001 continue                                                             1975
      r(i)=w(i)*(y(i,ic)-pic)                                              1975
      svr=svr+r(i)                                                         1976
14991 continue                                                             1977
14992 continue                                                             1977
      if(xm(0).le.vmin)goto 14981                                          1977
      ig=1                                                                 1978
15040 do 15041 j=1,ni                                                      1978
      if(iy(j).eq.0)goto 15041                                             1979
      jb=ix(j)                                                             1979
      je=ix(j+1)-1                                                         1980
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             1981
      if(kopt .ne. 0)goto 15061                                            1982
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       1983
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          1984
15061 continue                                                             1985
15041 continue                                                             1986
15042 continue                                                             1986
15070 continue                                                             1986
15071 continue                                                             1986
      nlp=nlp+1                                                            1986
      dlx=0.0                                                              1987
15080 do 15081 k=1,ni                                                      1987
      if(iy(k).eq.0)goto 15081                                             1988
      jb=ix(k)                                                             1988
      je=ix(k+1)-1                                                         1988
      jn=ix(k+1)-ix(k)                                                     1988
      bk=b(k,ic)                                                           1989
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1990
      gk=dot_product(sc(1:jn),x(jb:je))                                    1991
      gk=(gk-svr*xb(k))/xs(k)                                              1992
      u=gk+xv(k,ic)*b(k,ic)                                                1992
      au=abs(u)-vp(k)*al1                                                  1993
      if(au .gt. 0.0)goto 15101                                            1993
      b(k,ic)=0.0                                                          1993
      goto 15111                                                           1994
15101 continue                                                             1994
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1994
15111 continue                                                             1995
15091 continue                                                             1995
      d=b(k,ic)-bk                                                         1995
      if(abs(d).le.0.0)goto 15081                                          1996
      dlx=max(dlx,xv(k,ic)*d**2)                                           1997
      if(mm(k) .ne. 0)goto 15131                                           1997
      nin=nin+1                                                            1998
      if(nin .le. nx)goto 15151                                            1998
      jxx=1                                                                1998
      goto 15082                                                           1998
15151 continue                                                             1999
      mm(k)=nin                                                            1999
      m(nin)=k                                                             2000
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             2001
15131 continue                                                             2002
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2003
      o=o+d*(xb(k)/xs(k))                                                  2004
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2005
15081 continue                                                             2006
15082 continue                                                             2006
      if(jxx.gt.0)goto 15072                                               2007
      d=svr/xm(0)                                                          2008
      if(d .eq. 0.0)goto 15171                                             2008
      b(0,ic)=b(0,ic)+d                                                    2008
      dlx=max(dlx,xm(0)*d**2)                                              2009
      r=r-d*v                                                              2009
      svr=svr-d*xm(0)                                                      2010
15171 continue                                                             2011
      if(dlx.lt.shr)goto 15072                                             2011
      if(nlp .le. maxit)goto 15191                                         2011
      jerr=-ilm                                                            2011
      return                                                               2011
15191 continue                                                             2012
15200 continue                                                             2012
15201 continue                                                             2012
      nlp=nlp+1                                                            2012
      dlx=0.0                                                              2013
15210 do 15211 l=1,nin                                                     2013
      k=m(l)                                                               2013
      jb=ix(k)                                                             2013
      je=ix(k+1)-1                                                         2014
      jn=ix(k+1)-ix(k)                                                     2014
      bk=b(k,ic)                                                           2015
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2016
      gk=dot_product(sc(1:jn),x(jb:je))                                    2017
      gk=(gk-svr*xb(k))/xs(k)                                              2018
      u=gk+xv(k,ic)*b(k,ic)                                                2018
      au=abs(u)-vp(k)*al1                                                  2019
      if(au .gt. 0.0)goto 15231                                            2019
      b(k,ic)=0.0                                                          2019
      goto 15241                                                           2020
15231 continue                                                             2020
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              2020
15241 continue                                                             2021
15221 continue                                                             2021
      d=b(k,ic)-bk                                                         2021
      if(abs(d).le.0.0)goto 15211                                          2022
      dlx=max(dlx,xv(k,ic)*d**2)                                           2023
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2024
      o=o+d*(xb(k)/xs(k))                                                  2025
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2026
15211 continue                                                             2027
15212 continue                                                             2027
      d=svr/xm(0)                                                          2028
      if(d .eq. 0.0)goto 15261                                             2028
      b(0,ic)=b(0,ic)+d                                                    2028
      dlx=max(dlx,xm(0)*d**2)                                              2029
      r=r-d*v                                                              2029
      svr=svr-d*xm(0)                                                      2030
15261 continue                                                             2031
      if(dlx.lt.shr)goto 15202                                             2031
      if(nlp .le. maxit)goto 15281                                         2031
      jerr=-ilm                                                            2031
      return                                                               2031
15281 continue                                                             2032
      goto 15201                                                           2033
15202 continue                                                             2033
      goto 15071                                                           2034
15072 continue                                                             2034
      if(jxx.gt.0)goto 14982                                               2035
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                         2036
      if(ixx .ne. 0)goto 15301                                             2037
15310 do 15311 j=1,nin                                                     2037
      k=m(j)                                                               2038
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 15331                2038
      ixx=1                                                                2038
      goto 15312                                                           2038
15331 continue                                                             2039
15311 continue                                                             2040
15312 continue                                                             2040
15301 continue                                                             2041
      sc=b(0,ic)+g(:,ic)                                                   2041
      b0=0.0                                                               2042
15340 do 15341 j=1,nin                                                     2042
      l=m(j)                                                               2042
      jb=ix(l)                                                             2042
      je=ix(l+1)-1                                                         2043
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   2044
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            2045
15341 continue                                                             2046
15342 continue                                                             2046
      sc=min(max(exmn,sc+b0),exmx)                                         2047
      sxp=sxp-q(:,ic)                                                      2048
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          2049
      sxp=sxp+q(:,ic)                                                      2050
14981 continue                                                             2051
14982 continue                                                             2051
      s=-sum(b(0,:))/nc                                                    2051
      b(0,:)=b(0,:)+s                                                      2051
      sc=s                                                                 2051
      b0=0.0                                                               2052
15350 do 15351 j=1,nin                                                     2052
      l=m(j)                                                               2053
      if(vp(l) .gt. 0.0)goto 15371                                         2053
      s=sum(b(l,:))/nc                                                     2053
      goto 15381                                                           2054
15371 continue                                                             2054
      s=elc(parm,nc,b(l,:),is)                                             2054
15381 continue                                                             2055
15361 continue                                                             2055
      b(l,:)=b(l,:)-s                                                      2056
      jb=ix(l)                                                             2056
      je=ix(l+1)-1                                                         2057
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         2058
      b0=b0+s*xb(l)/xs(l)                                                  2059
15351 continue                                                             2060
15352 continue                                                             2060
      sc=sc+b0                                                             2060
      sc=exp(sc)                                                           2060
      sxp=sxp*sc                                                           2060
15390 do 15391 ic=1,nc                                                     2060
      q(:,ic)=q(:,ic)*sc                                                   2060
15391 continue                                                             2061
15392 continue                                                             2061
      if(jxx.gt.0)goto 14972                                               2061
      if(ig.eq.0)goto 14972                                                2062
      if(ixx .ne. 0)goto 15411                                             2063
15420 do 15421 j=1,ni                                                      2063
      if(iy(j).eq.1)goto 15421                                             2063
      if(ju(j).eq.0)goto 15421                                             2063
      ga(j)=0.0                                                            2063
15421 continue                                                             2064
15422 continue                                                             2064
15430 do 15431 ic=1,nc                                                     2064
      v=q(:,ic)/sxp                                                        2064
      r=w*(y(:,ic)-v)                                                      2064
      v=w*v*(1.0-v)                                                        2065
15440 do 15441 j=1,ni                                                      2065
      if(iy(j).eq.1)goto 15441                                             2065
      if(ju(j).eq.0)goto 15441                                             2066
      jb=ix(j)                                                             2066
      je=ix(j+1)-1                                                         2066
      jn=ix(j+1)-ix(j)                                                     2067
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2068
      gj=dot_product(sc(1:jn),x(jb:je))                                    2069
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2070
15441 continue                                                             2071
15442 continue                                                             2071
15431 continue                                                             2072
15432 continue                                                             2072
15450 do 15451 k=1,ni                                                      2072
      if(iy(k).eq.1)goto 15451                                             2072
      if(ju(k).eq.0)goto 15451                                             2073
      if(ga(k) .le. al1*vp(k))goto 15471                                   2073
      iy(k)=1                                                              2073
      ixx=1                                                                2073
15471 continue                                                             2074
15451 continue                                                             2075
15452 continue                                                             2075
      if(ixx.eq.1) go to 10680                                             2076
      goto 14972                                                           2077
15411 continue                                                             2078
      goto 14971                                                           2079
14972 continue                                                             2079
      if(jxx .le. 0)goto 15491                                             2079
      jerr=-10000-ilm                                                      2079
      goto 14892                                                           2079
15491 continue                                                             2079
      devi=0.0                                                             2080
15500 do 15501 ic=1,nc                                                     2081
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          2081
      a0(ic,ilm)=b(0,ic)                                                   2082
15510 do 15511 i=1,no                                                      2082
      if(y(i,ic).le.0.0)goto 15511                                         2083
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           2084
15511 continue                                                             2085
15512 continue                                                             2085
15501 continue                                                             2086
15502 continue                                                             2086
      kin(ilm)=nin                                                         2086
      alm(ilm)=al                                                          2086
      lmu=ilm                                                              2087
      dev(ilm)=(dev1-devi)/dev0                                            2087
      if(ig.eq.0)goto 14892                                                2088
      if(ilm.lt.mnl)goto 14891                                             2088
      if(flmin.ge.1.0)goto 14891                                           2089
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 14892             2090
      if(dev(ilm).gt.devmax)goto 14892                                     2090
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 14892                             2091
14891 continue                                                             2092
14892 continue                                                             2092
      g=log(q)                                                             2092
15520 do 15521 i=1,no                                                      2092
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2092
15521 continue                                                             2093
15522 continue                                                             2093
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy)                      2094
      return                                                               2095
      end                                                                  2096
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  2097
      real a0(nc),ca(nx,nc),x(*),f(nc,n)                                   2097
      integer ia(*),ix(*),jx(*)                                            2098
15530 do 15531 ic=1,nc                                                     2098
      f(ic,:)=a0(ic)                                                       2098
15531 continue                                                             2099
15532 continue                                                             2099
15540 do 15541 j=1,nin                                                     2099
      k=ia(j)                                                              2099
      kb=ix(k)                                                             2099
      ke=ix(k+1)-1                                                         2100
15550 do 15551 ic=1,nc                                                     2100
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    2100
15551 continue                                                             2101
15552 continue                                                             2101
15541 continue                                                             2102
15542 continue                                                             2102
      return                                                               2103
      end                                                                  2104
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,ne,nx,nlam,flmin,ula   2106 
     *m,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam)              2107
      real ca(nx,nlam),dev(nlam),alm(nlam)                                 2108
      integer jd(*),ia(nx),nin(nlam)                                       2109
      real, dimension (:), allocatable :: xs,ww,vq                              
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 15571                                    2113
      jerr=10000                                                           2113
      return                                                               2113
15571 continue                                                             2114
      allocate(ww(1:no),stat=jerr)                                         2115
      allocate(ju(1:ni),stat=ierr)                                         2115
      jerr=jerr+ierr                                                       2116
      allocate(vq(1:ni),stat=ierr)                                         2116
      jerr=jerr+ierr                                                       2117
      if(isd .le. 0)goto 15591                                             2117
      allocate(xs(1:ni),stat=ierr)                                         2117
      jerr=jerr+ierr                                                       2117
15591 continue                                                             2118
      if(jerr.ne.0) return                                                 2119
      call chkvars(no,ni,x,ju)                                             2120
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2121
      if(maxval(ju) .gt. 0)goto 15611                                      2121
      jerr=7777                                                            2121
      return                                                               2121
15611 continue                                                             2122
      vq=max(0.0,vp)                                                       2122
      vq=vq*ni/sum(vq)                                                     2123
      ww=max(0.0,w)                                                        2123
      sw=sum(ww)                                                           2124
      if(sw .gt. 0.0)goto 15631                                            2124
      jerr=9999                                                            2124
      return                                                               2124
15631 continue                                                             2124
      ww=ww/sw                                                             2125
      call cstandard(no,ni,x,ww,ju,isd,xs)                                 2126
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr   2128 
     *,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 2128
      dev0=2.0*sw*dev0                                                     2129
      if(isd .le. 0)goto 15651                                             2129
15660 do 15661 k=1,lmu                                                     2129
      nk=nin(k)                                                            2129
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   2129
15661 continue                                                             2129
15662 continue                                                             2129
15651 continue                                                             2130
      deallocate(ww,ju,vq)                                                 2130
      if(isd.gt.0) deallocate(xs)                                          2131
      return                                                               2132
      end                                                                  2133
      subroutine cstandard (no,ni,x,w,ju,isd,xs)                           2134
      real x(no,ni),w(no),xs(ni)                                           2134
      integer ju(ni)                                                       2135
15670 do 15671 j=1,ni                                                      2135
      if(ju(j).eq.0)goto 15671                                             2136
      xm=dot_product(w,x(:,j))                                             2136
      x(:,j)=x(:,j)-xm                                                     2137
      if(isd .le. 0)goto 15691                                             2137
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 2137
      x(:,j)=x(:,j)/xs(j)                                                  2137
15691 continue                                                             2138
15671 continue                                                             2139
15672 continue                                                             2139
      return                                                               2140
      end                                                                  2141
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,ne,nx,nlam,flmin,ula   2143 
     *m,cthri,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=0.001, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99)   2144
      real x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam)              2145
      real ao(nx,nlam),dev(nlam),alm(nlam)                                 2146
      integer ju(ni),m(nx),kin(nlam)                                       2147
      real, dimension (:), allocatable :: w,dk,v,xs,wr,a,as,f,dq                
      real, dimension (:), allocatable :: e,uu,ga                               
      integer, dimension (:), allocatable :: jp,kp,mm,ixx                       
      allocate(e(1:no),stat=jerr)                                          2153
      allocate(uu(1:no),stat=ierr)                                         2153
      jerr=jerr+ierr                                                       2154
      allocate(f(1:no),stat=ierr)                                          2154
      jerr=jerr+ierr                                                       2155
      allocate(w(1:no),stat=ierr)                                          2155
      jerr=jerr+ierr                                                       2156
      allocate(v(1:ni),stat=ierr)                                          2156
      jerr=jerr+ierr                                                       2157
      allocate(a(1:ni),stat=ierr)                                          2157
      jerr=jerr+ierr                                                       2158
      allocate(as(1:ni),stat=ierr)                                         2158
      jerr=jerr+ierr                                                       2159
      allocate(xs(1:ni),stat=ierr)                                         2159
      jerr=jerr+ierr                                                       2160
      allocate(ga(1:ni),stat=ierr)                                         2160
      jerr=jerr+ierr                                                       2161
      allocate(ixx(1:ni),stat=ierr)                                        2161
      jerr=jerr+ierr                                                       2162
      allocate(jp(1:no),stat=ierr)                                         2162
      jerr=jerr+ierr                                                       2163
      allocate(kp(1:no),stat=ierr)                                         2163
      jerr=jerr+ierr                                                       2164
      allocate(dk(1:no),stat=ierr)                                         2164
      jerr=jerr+ierr                                                       2165
      allocate(wr(1:no),stat=ierr)                                         2165
      jerr=jerr+ierr                                                       2166
      allocate(dq(1:no),stat=ierr)                                         2166
      jerr=jerr+ierr                                                       2167
      allocate(mm(1:ni),stat=ierr)                                         2167
      jerr=jerr+ierr                                                       2168
      if(jerr.ne.0)go to 11790                                             2169
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2170
      if(jerr.ne.0) go to 11790                                            2170
      alpha=parm                                                           2171
      oma=1.0-alpha                                                        2171
      nlm=0                                                                2171
      ixx=0                                                                2171
      al=0.0                                                               2172
      dq=d*q                                                               2172
      call died(no,nk,dq,kp,jp,dk)                                         2173
      a=0.0                                                                2173
      f(1)=0.0                                                             2173
      fmax=log(huge(f(1))*0.1)                                             2174
      if(nonzero(no,g) .eq. 0)goto 15711                                   2174
      f=g-dot_product(q,g)                                                 2175
      e=q*exp(sign(min(abs(f),fmax),f))                                    2176
      goto 15721                                                           2177
15711 continue                                                             2177
      f=0.0                                                                2177
      e=q                                                                  2177
15721 continue                                                             2178
15701 continue                                                             2178
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 2179
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         2179
      dev0=rr                                                              2180
15730 do 15731 i=1,no                                                      2180
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 15751                   2180
      w(i)=0.0                                                             2180
      wr(i)=w(i)                                                           2180
15751 continue                                                             2180
15731 continue                                                             2181
15732 continue                                                             2181
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2182
      if(jerr.ne.0) go to 11790                                            2183
      if(flmin .ge. 1.0)goto 15771                                         2183
      eqs=max(eps,flmin)                                                   2183
      alf=eqs**(1.0/(nlam-1))                                              2183
15771 continue                                                             2184
      m=0                                                                  2184
      mm=0                                                                 2184
      nlp=0                                                                2184
      nin=nlp                                                              2184
      mnl=min(mnlam,nlam)                                                  2184
      as=0.0                                                               2184
      cthr=cthri*dev0                                                      2185
15780 do 15781 j=1,ni                                                      2185
      if(ju(j).eq.0)goto 15781                                             2185
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2185
15781 continue                                                             2186
15782 continue                                                             2186
15790 do 15791 ilm=1,nlam                                                  2186
      al0=al                                                               2187
      if(flmin .lt. 1.0)goto 15811                                         2187
      al=ulam(ilm)                                                         2187
      goto 15801                                                           2188
15811 if(ilm .le. 2)goto 15821                                             2188
      al=al*alf                                                            2188
      goto 15801                                                           2189
15821 if(ilm .ne. 1)goto 15831                                             2189
      al=big                                                               2189
      goto 15841                                                           2190
15831 continue                                                             2190
      al0=0.0                                                              2191
15850 do 15851 j=1,ni                                                      2191
      if(ju(j).eq.0)goto 15851                                             2191
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2191
15851 continue                                                             2192
15852 continue                                                             2192
      al0=al0/max(parm,1.0e-3)                                             2192
      al=alf*al0                                                           2193
15841 continue                                                             2194
15801 continue                                                             2194
      sa=alpha*al                                                          2194
      omal=oma*al                                                          2194
      tlam=alpha*(2.0*al-al0)                                              2195
15860 do 15861 k=1,ni                                                      2195
      if(ixx(k).eq.1)goto 15861                                            2195
      if(ju(k).eq.0)goto 15861                                             2196
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2197
15861 continue                                                             2198
15862 continue                                                             2198
10680 continue                                                             2199
15870 continue                                                             2199
15871 continue                                                             2199
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2200
      call vars(no,ni,x,w,ixx,v)                                           2201
15880 continue                                                             2201
15881 continue                                                             2201
      nlp=nlp+1                                                            2201
      dli=0.0                                                              2202
15890 do 15891 j=1,ni                                                      2202
      if(ixx(j).eq.0)goto 15891                                            2203
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2204
      if(abs(u) .gt. vp(j)*sa)goto 15911                                   2204
      at=0.0                                                               2204
      goto 15921                                                           2205
15911 continue                                                             2205
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2205
15921 continue                                                             2206
15901 continue                                                             2206
      if(at .eq. a(j))goto 15941                                           2206
      del=at-a(j)                                                          2206
      a(j)=at                                                              2206
      dli=max(dli,v(j)*del**2)                                             2207
      wr=wr-del*w*x(:,j)                                                   2207
      f=f+del*x(:,j)                                                       2208
      if(mm(j) .ne. 0)goto 15961                                           2208
      nin=nin+1                                                            2208
      if(nin.gt.nx)goto 15892                                              2209
      mm(j)=nin                                                            2209
      m(nin)=j                                                             2210
15961 continue                                                             2211
15941 continue                                                             2212
15891 continue                                                             2213
15892 continue                                                             2213
      if(nin.gt.nx)goto 15882                                              2213
      if(dli.lt.cthr)goto 15882                                            2214
      if(nlp .le. maxit)goto 15981                                         2214
      jerr=-ilm                                                            2214
      return                                                               2214
15981 continue                                                             2215
15990 continue                                                             2215
15991 continue                                                             2215
      nlp=nlp+1                                                            2215
      dli=0.0                                                              2216
16000 do 16001 l=1,nin                                                     2216
      j=m(l)                                                               2217
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2218
      if(abs(u) .gt. vp(j)*sa)goto 16021                                   2218
      at=0.0                                                               2218
      goto 16031                                                           2219
16021 continue                                                             2219
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2219
16031 continue                                                             2220
16011 continue                                                             2220
      if(at .eq. a(j))goto 16051                                           2220
      del=at-a(j)                                                          2220
      a(j)=at                                                              2220
      dli=max(dli,v(j)*del**2)                                             2221
      wr=wr-del*w*x(:,j)                                                   2221
      f=f+del*x(:,j)                                                       2222
16051 continue                                                             2223
16001 continue                                                             2224
16002 continue                                                             2224
      if(dli.lt.cthr)goto 15992                                            2224
      if(nlp .le. maxit)goto 16071                                         2224
      jerr=-ilm                                                            2224
      return                                                               2224
16071 continue                                                             2225
      goto 15991                                                           2226
15992 continue                                                             2226
      goto 15881                                                           2227
15882 continue                                                             2227
      if(nin.gt.nx)goto 15872                                              2228
      e=q*exp(sign(min(abs(f),fmax),f))                                    2229
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2230
      if(jerr .eq. 0)goto 16091                                            2230
      jerr=jerr-ilm                                                        2230
      go to 11790                                                          2230
16091 continue                                                             2231
      ix=0                                                                 2232
16100 do 16101 j=1,nin                                                     2232
      k=m(j)                                                               2233
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 16101                           2233
      ix=1                                                                 2233
      goto 16102                                                           2233
16101 continue                                                             2234
16102 continue                                                             2234
      if(ix .ne. 0)goto 16121                                              2235
16130 do 16131 k=1,ni                                                      2235
      if(ixx(k).eq.1)goto 16131                                            2235
      if(ju(k).eq.0)goto 16131                                             2236
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2237
      if(ga(k) .le. sa*vp(k))goto 16151                                    2237
      ixx(k)=1                                                             2237
      ix=1                                                                 2237
16151 continue                                                             2238
16131 continue                                                             2239
16132 continue                                                             2239
      if(ix.eq.1) go to 10680                                              2240
      goto 15872                                                           2241
16121 continue                                                             2242
      goto 15871                                                           2243
15872 continue                                                             2243
      if(nin .le. nx)goto 16171                                            2243
      jerr=-10000-ilm                                                      2243
      goto 15792                                                           2243
16171 continue                                                             2244
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               2244
      kin(ilm)=nin                                                         2245
      alm(ilm)=al                                                          2245
      lmu=ilm                                                              2246
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   2247
      if(ilm.lt.mnl)goto 15791                                             2247
      if(flmin.ge.1.0)goto 15791                                           2248
      me=0                                                                 2248
16180 do 16181 j=1,nin                                                     2248
      if(ao(j,ilm).ne.0.0) me=me+1                                         2248
16181 continue                                                             2248
16182 continue                                                             2248
      if(me.gt.ne)goto 15792                                               2249
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 15792              2250
      if(dev(ilm).gt.devmax)goto 15792                                     2251
15791 continue                                                             2252
15792 continue                                                             2252
      g=f                                                                  2253
11790 continue                                                             2253
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)              2254
      return                                                               2255
      end                                                                  2256
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 2257
      real ca(nin),x(n,*),f(n)                                             2257
      integer ia(nin)                                                      2258
      f=0.0                                                                2258
      if(nin.le.0) return                                                  2259
16190 do 16191 i=1,n                                                       2259
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      2259
16191 continue                                                             2260
16192 continue                                                             2260
      return                                                               2261
      end                                                                  2262
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         2263
      real y(no),d(no),q(no)                                               2263
      integer jp(no),kp(*)                                                 2264
16200 do 16201 j=1,no                                                      2264
      jp(j)=j                                                              2264
16201 continue                                                             2264
16202 continue                                                             2264
      call psort7(y,jp,1,no)                                               2265
      nj=0                                                                 2265
16210 do 16211 j=1,no                                                      2265
      if(q(jp(j)).le.0.0)goto 16211                                        2265
      nj=nj+1                                                              2265
      jp(nj)=jp(j)                                                         2265
16211 continue                                                             2266
16212 continue                                                             2266
      if(nj .ne. 0)goto 16231                                              2266
      jerr=20000                                                           2266
      return                                                               2266
16231 continue                                                             2267
      j=1                                                                  2267
16240 continue                                                             2267
16241 if(d(jp(j)).gt.0.0)goto 16242                                        2267
      j=j+1                                                                2267
      if(j.gt.nj)goto 16242                                                2267
      goto 16241                                                           2268
16242 continue                                                             2268
      if(j .lt. nj-1)goto 16261                                            2268
      jerr=30000                                                           2268
      return                                                               2268
16261 continue                                                             2269
      j0=j-1                                                               2269
      nj=nj-j0                                                             2269
16270 do 16271 j=1,nj                                                      2269
      jp(j)=jp(j+j0)                                                       2269
16271 continue                                                             2270
16272 continue                                                             2270
      jerr=0                                                               2270
      nk=0                                                                 2270
      t0=y(jp(1))                                                          2270
      yk=t0                                                                2270
      j=2                                                                  2271
16280 continue                                                             2271
16281 continue                                                             2271
16290 continue                                                             2272
16291 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 16292                     2272
      j=j+1                                                                2272
      if(j.gt.nj)goto 16292                                                2272
      goto 16291                                                           2273
16292 continue                                                             2273
      nk=nk+1                                                              2273
      kp(nk)=j-1                                                           2273
      if(j.gt.nj)goto 16282                                                2274
      if(j .ne. nj)goto 16311                                              2274
      nk=nk+1                                                              2274
      kp(nk)=nj                                                            2274
      goto 16282                                                           2274
16311 continue                                                             2275
      yk=y(jp(j))                                                          2275
      j=j+1                                                                2276
      goto 16281                                                           2277
16282 continue                                                             2277
      return                                                               2278
      end                                                                  2279
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     2280
      real d(no),dk(nk),wr(no),w(no)                                       2281
      real e(no),u(no),b,c                                                 2281
      integer kp(nk),jp(no)                                                2282
      call usk(no,nk,kp,jp,e,u)                                            2283
      b=dk(1)/u(1)                                                         2283
      c=dk(1)/u(1)**2                                                      2283
      jerr=0                                                               2284
16320 do 16321 j=1,kp(1)                                                   2284
      i=jp(j)                                                              2285
      w(i)=e(i)*(b-e(i)*c)                                                 2285
      if(w(i) .gt. 0.0)goto 16341                                          2285
      jerr=-3                                                              2285
      return                                                               2285
16341 continue                                                             2286
      wr(i)=d(i)-e(i)*b                                                    2287
16321 continue                                                             2288
16322 continue                                                             2288
16350 do 16351 k=2,nk                                                      2288
      j1=kp(k-1)+1                                                         2288
      j2=kp(k)                                                             2289
      b=b+dk(k)/u(k)                                                       2289
      c=c+dk(k)/u(k)**2                                                    2290
16360 do 16361 j=j1,j2                                                     2290
      i=jp(j)                                                              2291
      w(i)=e(i)*(b-e(i)*c)                                                 2291
      if(w(i) .gt. 0.0)goto 16381                                          2291
      jerr=-30000                                                          2291
      return                                                               2291
16381 continue                                                             2292
      wr(i)=d(i)-e(i)*b                                                    2293
16361 continue                                                             2294
16362 continue                                                             2294
16351 continue                                                             2295
16352 continue                                                             2295
      return                                                               2296
      end                                                                  2297
      subroutine vars(no,ni,x,w,ixx,v)                                     2298
      real x(no,ni),w(no),v(ni)                                            2298
      integer ixx(ni)                                                      2299
16390 do 16391 j=1,ni                                                      2299
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                        2299
16391 continue                                                             2300
16392 continue                                                             2300
      return                                                               2301
      end                                                                  2302
      subroutine died(no,nk,d,kp,jp,dk)                                    2303
      real d(no),dk(nk)                                                    2303
      integer kp(nk),jp(no)                                                2304
      dk(1)=sum(d(jp(1:kp(1))))                                            2305
16400 do 16401 k=2,nk                                                      2305
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  2305
16401 continue                                                             2306
16402 continue                                                             2306
      return                                                               2307
      end                                                                  2308
      subroutine usk(no,nk,kp,jp,e,u)                                      2309
      real e(no),u(nk),h                                                   2309
      integer kp(nk),jp(no)                                                2310
      h=0.0                                                                2311
16410 do 16411 k=nk,1,-1                                                   2311
      j2=kp(k)                                                             2312
      j1=1                                                                 2312
      if(k.gt.1) j1=kp(k-1)+1                                              2313
16420 do 16421 j=j2,j1,-1                                                  2313
      h=h+e(jp(j))                                                         2313
16421 continue                                                             2314
16422 continue                                                             2314
      u(k)=h                                                               2315
16411 continue                                                             2316
16412 continue                                                             2316
      return                                                               2317
      end                                                                  2318
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             2319
      real d(no),dk(nk),f(no)                                              2320
      integer kp(nk),jp(no)                                                2320
      real e(no),u(nk),s                                                   2321
      call usk(no,nk,kp,jp,e,u)                                            2321
      u=log(u)                                                             2322
      risk=dot_product(d,f)-dot_product(dk,u)                              2323
      return                                                               2324
      end                                                                  2325
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 2326
      real x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(nlam)          2327
      real, dimension (:), allocatable :: dk,f,xm,dq,q                          
      real, dimension (:), allocatable :: e,uu                                  
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          2333
      allocate(q(1:no),stat=ierr)                                          2333
      jerr=jerr+ierr                                                       2334
      allocate(uu(1:no),stat=ierr)                                         2334
      jerr=jerr+ierr                                                       2335
      allocate(f(1:no),stat=ierr)                                          2335
      jerr=jerr+ierr                                                       2336
      allocate(dk(1:no),stat=ierr)                                         2336
      jerr=jerr+ierr                                                       2337
      allocate(jp(1:no),stat=ierr)                                         2337
      jerr=jerr+ierr                                                       2338
      allocate(kp(1:no),stat=ierr)                                         2338
      jerr=jerr+ierr                                                       2339
      allocate(dq(1:no),stat=ierr)                                         2339
      jerr=jerr+ierr                                                       2340
      allocate(xm(1:ni),stat=ierr)                                         2340
      jerr=jerr+ierr                                                       2341
      if(jerr.ne.0) go to 11790                                            2342
      q=max(0.0,w)                                                         2342
      sw=sum(q)                                                            2343
      if(sw .gt. 0.0)goto 16441                                            2343
      jerr=9999                                                            2343
      go to 11790                                                          2343
16441 continue                                                             2344
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2345
      if(jerr.ne.0) go to 11790                                            2345
      fmax=log(huge(e(1))*0.1)                                             2346
      dq=d*q                                                               2346
      call died(no,nk,dq,kp,jp,dk)                                         2346
      gm=dot_product(q,g)/sw                                               2347
16450 do 16451 j=1,ni                                                      2347
      xm(j)=dot_product(q,x(:,j))/sw                                       2347
16451 continue                                                             2348
16452 continue                                                             2348
16460 do 16461 lam=1,nlam                                                  2349
16470 do 16471 i=1,no                                                      2349
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       2350
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        2351
16471 continue                                                             2352
16472 continue                                                             2352
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          2353
16461 continue                                                             2354
16462 continue                                                             2354
11790 continue                                                             2354
      deallocate(e,uu,dk,f,jp,kp,dq)                                       2355
      return                                                               2356
      end                                                                  2357
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,ne,nx,nlam,flmin,ulam   2359 
     *,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    2360
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2361
      integer jd(*),ia(nx),nin(nlam)                                       2362
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16491                                    2366
      jerr=10000                                                           2366
      return                                                               2366
16491 continue                                                             2367
      if(minval(y) .ge. 0.0)goto 16511                                     2367
      jerr=8888                                                            2367
      return                                                               2367
16511 continue                                                             2368
      allocate(ww(1:no),stat=jerr)                                         2369
      allocate(ju(1:ni),stat=ierr)                                         2369
      jerr=jerr+ierr                                                       2370
      allocate(vq(1:ni),stat=ierr)                                         2370
      jerr=jerr+ierr                                                       2371
      allocate(xm(1:ni),stat=ierr)                                         2371
      jerr=jerr+ierr                                                       2372
      if(isd .le. 0)goto 16531                                             2372
      allocate(xs(1:ni),stat=ierr)                                         2372
      jerr=jerr+ierr                                                       2372
16531 continue                                                             2373
      if(jerr.ne.0) return                                                 2374
      call chkvars(no,ni,x,ju)                                             2375
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2376
      if(maxval(ju) .gt. 0)goto 16551                                      2376
      jerr=7777                                                            2376
      go to 11790                                                          2376
16551 continue                                                             2377
      vq=max(0.0,vp)                                                       2377
      vq=vq*ni/sum(vq)                                                     2378
      ww=max(0.0,w)                                                        2378
      sw=sum(ww)                                                           2378
      if(sw .gt. 0.0)goto 16571                                            2378
      jerr=9999                                                            2378
      go to 11790                                                          2378
16571 continue                                                             2379
      ww=ww/sw                                                             2380
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             2381
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr,   2383 
     *  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11790                                            2383
      dev0=2.0*sw*dev0                                                     2384
16580 do 16581 k=1,lmu                                                     2384
      nk=nin(k)                                                            2385
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2386
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2387
16581 continue                                                             2388
16582 continue                                                             2388
11790 continue                                                             2388
      deallocate(ww,ju,vq,xm)                                              2388
      if(isd.gt.0) deallocate(xs)                                          2389
      return                                                               2390
      end                                                                  2391
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,ne,nx,nlam,flmin,ulam   2393 
     *,shri,  isd,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2394 
     *9)
      real x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)                    2395
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2396
      integer ju(ni),m(nx),kin(nlam)                                       2397
      real, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga                    
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(a(1:ni),stat=jerr)                                          2402
      allocate(as(1:ni),stat=ierr)                                         2402
      jerr=jerr+ierr                                                       2403
      allocate(t(1:no),stat=ierr)                                          2403
      jerr=jerr+ierr                                                       2404
      allocate(mm(1:ni),stat=ierr)                                         2404
      jerr=jerr+ierr                                                       2405
      allocate(ga(1:ni),stat=ierr)                                         2405
      jerr=jerr+ierr                                                       2406
      allocate(ixx(1:ni),stat=ierr)                                        2406
      jerr=jerr+ierr                                                       2407
      allocate(wr(1:no),stat=ierr)                                         2407
      jerr=jerr+ierr                                                       2408
      allocate(v(1:ni),stat=ierr)                                          2408
      jerr=jerr+ierr                                                       2409
      allocate(w(1:no),stat=ierr)                                          2409
      jerr=jerr+ierr                                                       2410
      allocate(f(1:no),stat=ierr)                                          2410
      jerr=jerr+ierr                                                       2411
      if(jerr.ne.0) return                                                 2412
      bta=parm                                                             2412
      omb=1.0-bta                                                          2413
      t=q*y                                                                2413
      yb=sum(t)                                                            2413
      fmax=log(huge(bta)*0.1)                                              2414
      if(nonzero(no,g) .ne. 0)goto 16601                                   2414
      w=q*yb                                                               2414
      az=log(yb)                                                           2414
      f=az                                                                 2414
      dv0=yb*(log(yb)-1.0)                                                 2414
      goto 16611                                                           2415
16601 continue                                                             2415
      w=q*exp(sign(min(abs(g),fmax),g))                                    2415
      v0=sum(w)                                                            2415
      eaz=yb/v0                                                            2416
      w=eaz*w                                                              2416
      az=log(eaz)                                                          2416
      f=az+g                                                               2417
      dv0=dot_product(t,g)-yb*(1.0-az)                                     2418
16611 continue                                                             2419
16591 continue                                                             2419
      a=0.0                                                                2419
      as=0.0                                                               2419
      wr=t-w                                                               2419
      v0=yb                                                                2419
      dvr=-yb                                                              2420
16620 do 16621 i=1,no                                                      2420
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               2420
16621 continue                                                             2420
16622 continue                                                             2420
      dvr=dvr-dv0                                                          2420
      dev0=dvr                                                             2421
      if(flmin .ge. 1.0)goto 16641                                         2421
      eqs=max(eps,flmin)                                                   2421
      alf=eqs**(1.0/(nlam-1))                                              2421
16641 continue                                                             2422
      m=0                                                                  2422
      mm=0                                                                 2422
      nlp=0                                                                2422
      nin=nlp                                                              2422
      mnl=min(mnlam,nlam)                                                  2422
      shr=shri*dev0                                                        2422
      ixx=0                                                                2422
      al=0.0                                                               2423
16650 do 16651 j=1,ni                                                      2423
      if(ju(j).eq.0)goto 16651                                             2423
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2423
16651 continue                                                             2424
16652 continue                                                             2424
16660 do 16661 ilm=1,nlam                                                  2424
      al0=al                                                               2425
      if(flmin .lt. 1.0)goto 16681                                         2425
      al=ulam(ilm)                                                         2425
      goto 16671                                                           2426
16681 if(ilm .le. 2)goto 16691                                             2426
      al=al*alf                                                            2426
      goto 16671                                                           2427
16691 if(ilm .ne. 1)goto 16701                                             2427
      al=big                                                               2427
      goto 16711                                                           2428
16701 continue                                                             2428
      al0=0.0                                                              2429
16720 do 16721 j=1,ni                                                      2429
      if(ju(j).eq.0)goto 16721                                             2429
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2429
16721 continue                                                             2430
16722 continue                                                             2430
      al0=al0/max(bta,1.0e-3)                                              2430
      al=alf*al0                                                           2431
16711 continue                                                             2432
16671 continue                                                             2432
      al2=al*omb                                                           2432
      al1=al*bta                                                           2432
      tlam=bta*(2.0*al-al0)                                                2433
16730 do 16731 k=1,ni                                                      2433
      if(ixx(k).eq.1)goto 16731                                            2433
      if(ju(k).eq.0)goto 16731                                             2434
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2435
16731 continue                                                             2436
16732 continue                                                             2436
10680 continue                                                             2437
16740 continue                                                             2437
16741 continue                                                             2437
      az0=az                                                               2438
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2439
16750 do 16751 j=1,ni                                                      2439
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                        2439
16751 continue                                                             2440
16752 continue                                                             2440
16760 continue                                                             2440
16761 continue                                                             2440
      nlp=nlp+1                                                            2440
      dlx=0.0                                                              2441
16770 do 16771 k=1,ni                                                      2441
      if(ixx(k).eq.0)goto 16771                                            2441
      ak=a(k)                                                              2442
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2442
      au=abs(u)-vp(k)*al1                                                  2443
      if(au .gt. 0.0)goto 16791                                            2443
      a(k)=0.0                                                             2443
      goto 16801                                                           2444
16791 continue                                                             2444
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2444
16801 continue                                                             2445
16781 continue                                                             2445
      if(a(k).eq.ak)goto 16771                                             2445
      d=a(k)-ak                                                            2445
      dlx=max(dlx,v(k)*d**2)                                               2446
      wr=wr-d*w*x(:,k)                                                     2446
      f=f+d*x(:,k)                                                         2447
      if(mm(k) .ne. 0)goto 16821                                           2447
      nin=nin+1                                                            2447
      if(nin.gt.nx)goto 16772                                              2448
      mm(k)=nin                                                            2448
      m(nin)=k                                                             2449
16821 continue                                                             2450
16771 continue                                                             2451
16772 continue                                                             2451
      if(nin.gt.nx)goto 16762                                              2451
      d=sum(wr)/v0                                                         2452
      az=az+d                                                              2452
      dlx=max(dlx,v0*d**2)                                                 2452
      wr=wr-d*w                                                            2452
      f=f+d                                                                2453
      if(dlx.lt.shr)goto 16762                                             2453
      if(nlp .le. maxit)goto 16841                                         2453
      jerr=-ilm                                                            2453
      return                                                               2453
16841 continue                                                             2454
16850 continue                                                             2454
16851 continue                                                             2454
      nlp=nlp+1                                                            2454
      dlx=0.0                                                              2455
16860 do 16861 l=1,nin                                                     2455
      k=m(l)                                                               2455
      ak=a(k)                                                              2456
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2456
      au=abs(u)-vp(k)*al1                                                  2457
      if(au .gt. 0.0)goto 16881                                            2457
      a(k)=0.0                                                             2457
      goto 16891                                                           2458
16881 continue                                                             2458
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2458
16891 continue                                                             2459
16871 continue                                                             2459
      if(a(k).eq.ak)goto 16861                                             2459
      d=a(k)-ak                                                            2459
      dlx=max(dlx,v(k)*d**2)                                               2460
      wr=wr-d*w*x(:,k)                                                     2460
      f=f+d*x(:,k)                                                         2462
16861 continue                                                             2462
16862 continue                                                             2462
      d=sum(wr)/v0                                                         2462
      az=az+d                                                              2462
      dlx=max(dlx,v0*d**2)                                                 2462
      wr=wr-d*w                                                            2462
      f=f+d                                                                2463
      if(dlx.lt.shr)goto 16852                                             2463
      if(nlp .le. maxit)goto 16911                                         2463
      jerr=-ilm                                                            2463
      return                                                               2463
16911 continue                                                             2464
      goto 16851                                                           2465
16852 continue                                                             2465
      goto 16761                                                           2466
16762 continue                                                             2466
      if(nin.gt.nx)goto 16742                                              2467
      w=q*exp(sign(min(abs(f),fmax),f))                                    2467
      v0=sum(w)                                                            2467
      wr=t-w                                                               2468
      if(v0*(az-az0)**2 .ge. shr)goto 16931                                2468
      ix=0                                                                 2469
16940 do 16941 j=1,nin                                                     2469
      k=m(j)                                                               2470
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 16941                            2470
      ix=1                                                                 2470
      goto 16942                                                           2471
16941 continue                                                             2472
16942 continue                                                             2472
      if(ix .ne. 0)goto 16961                                              2473
16970 do 16971 k=1,ni                                                      2473
      if(ixx(k).eq.1)goto 16971                                            2473
      if(ju(k).eq.0)goto 16971                                             2474
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2475
      if(ga(k) .le. al1*vp(k))goto 16991                                   2475
      ixx(k)=1                                                             2475
      ix=1                                                                 2475
16991 continue                                                             2476
16971 continue                                                             2477
16972 continue                                                             2477
      if(ix.eq.1) go to 10680                                              2478
      goto 16742                                                           2479
16961 continue                                                             2480
16931 continue                                                             2481
      goto 16741                                                           2482
16742 continue                                                             2482
      if(nin .le. nx)goto 17011                                            2482
      jerr=-10000-ilm                                                      2482
      goto 16662                                                           2482
17011 continue                                                             2483
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2483
      kin(ilm)=nin                                                         2484
      a0(ilm)=az                                                           2484
      alm(ilm)=al                                                          2484
      lmu=ilm                                                              2485
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               2486
      if(ilm.lt.mnl)goto 16661                                             2486
      if(flmin.ge.1.0)goto 16661                                           2487
      me=0                                                                 2487
17020 do 17021 j=1,nin                                                     2487
      if(ca(j,ilm).ne.0.0) me=me+1                                         2487
17021 continue                                                             2487
17022 continue                                                             2487
      if(me.gt.ne)goto 16662                                               2488
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 16662              2489
      if(dev(ilm).gt.devmax)goto 16662                                     2490
16661 continue                                                             2491
16662 continue                                                             2491
      g=f                                                                  2492
11790 continue                                                             2492
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                                2493
      return                                                               2494
      end                                                                  2495
      function nonzero(n,v)                                                2496
      real v(n)                                                            2497
      nonzero=0                                                            2497
17030 do 17031 i=1,n                                                       2497
      if(v(i) .eq. 0.0)goto 17051                                          2497
      nonzero=1                                                            2497
      return                                                               2497
17051 continue                                                             2497
17031 continue                                                             2498
17032 continue                                                             2498
      return                                                               2499
      end                                                                  2500
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               2501
      real a(nx,lmu),b(ni,lmu)                                             2501
      integer ia(nx),nin(lmu)                                              2502
17060 do 17061 lam=1,lmu                                                   2502
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        2502
17061 continue                                                             2503
17062 continue                                                             2503
      return                                                               2504
      end                                                                  2505
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           2506
      real a(nx,nc,lmu),b(ni,nc,lmu)                                       2506
      integer ia(nx),nin(lmu)                                              2507
17070 do 17071 lam=1,lmu                                                   2507
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             2507
17071 continue                                                             2508
17072 continue                                                             2508
      return                                                               2509
      end                                                                  2510
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               2511
      real x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)       2512
      real, dimension (:), allocatable :: w                                     
      if(minval(y) .ge. 0.0)goto 17091                                     2515
      jerr=8888                                                            2515
      return                                                               2515
17091 continue                                                             2516
      allocate(w(1:no),stat=jerr)                                          2516
      if(jerr.ne.0) return                                                 2517
      w=max(0.0,q)                                                         2517
      sw=sum(w)                                                            2517
      if(sw .gt. 0.0)goto 17111                                            2517
      jerr=9999                                                            2517
      go to 11790                                                          2517
17111 continue                                                             2518
      yb=dot_product(w,y)/sw                                               2518
      fmax=log(huge(y(1))*0.1)                                             2519
17120 do 17121 lam=1,nlam                                                  2519
      s=0.0                                                                2520
17130 do 17131 i=1,no                                                      2520
      if(w(i).le.0.0)goto 17131                                            2521
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          2522
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      2523
17131 continue                                                             2524
17132 continue                                                             2524
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2525
17121 continue                                                             2526
17122 continue                                                             2526
11790 continue                                                             2526
      deallocate(w)                                                        2527
      return                                                               2528
      end                                                                  2529
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,ne,nx,nlam,fl   2531 
     *min,ulam,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        2532
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2533
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           2534
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 17151                                    2538
      jerr=10000                                                           2538
      return                                                               2538
17151 continue                                                             2539
      if(minval(y) .ge. 0.0)goto 17171                                     2539
      jerr=8888                                                            2539
      return                                                               2539
17171 continue                                                             2540
      allocate(ww(1:no),stat=jerr)                                         2541
      allocate(ju(1:ni),stat=ierr)                                         2541
      jerr=jerr+ierr                                                       2542
      allocate(vq(1:ni),stat=ierr)                                         2542
      jerr=jerr+ierr                                                       2543
      allocate(xm(1:ni),stat=ierr)                                         2543
      jerr=jerr+ierr                                                       2544
      allocate(xs(1:ni),stat=ierr)                                         2544
      jerr=jerr+ierr                                                       2545
      if(jerr.ne.0) return                                                 2546
      call spchkvars(no,ni,x,ix,ju)                                        2547
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2548
      if(maxval(ju) .gt. 0)goto 17191                                      2548
      jerr=7777                                                            2548
      go to 11790                                                          2548
17191 continue                                                             2549
      vq=max(0.0,vp)                                                       2549
      vq=vq*ni/sum(vq)                                                     2550
      ww=max(0.0,w)                                                        2550
      sw=sum(ww)                                                           2550
      if(sw .gt. 0.0)goto 17211                                            2550
      jerr=9999                                                            2550
      go to 11790                                                          2550
17211 continue                                                             2551
      ww=ww/sw                                                             2552
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     2553
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmin,u   2555 
     *lam,thr,  isd,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11790                                            2555
      dev0=2.0*sw*dev0                                                     2556
17220 do 17221 k=1,lmu                                                     2556
      nk=nin(k)                                                            2557
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2558
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2559
17221 continue                                                             2560
17222 continue                                                             2560
11790 continue                                                             2560
      deallocate(ww,ju,vq,xm,xs)                                           2561
      return                                                               2562
      end                                                                  2563
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,ne,nx,nlam,fl   2565 
     *min,ulam,  shri,isd,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,j
     *err)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2566 
     *9)
      real x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),xs(ni)          2567
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2568
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2569
      real, dimension (:), allocatable :: qy,t,w,wr,v,a,as,xm,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(a(1:ni),stat=jerr)                                          2574
      allocate(as(1:ni),stat=ierr)                                         2574
      jerr=jerr+ierr                                                       2575
      allocate(t(1:no),stat=ierr)                                          2575
      jerr=jerr+ierr                                                       2576
      allocate(mm(1:ni),stat=ierr)                                         2576
      jerr=jerr+ierr                                                       2577
      allocate(ga(1:ni),stat=ierr)                                         2577
      jerr=jerr+ierr                                                       2578
      allocate(ixx(1:ni),stat=ierr)                                        2578
      jerr=jerr+ierr                                                       2579
      allocate(wr(1:no),stat=ierr)                                         2579
      jerr=jerr+ierr                                                       2580
      allocate(v(1:ni),stat=ierr)                                          2580
      jerr=jerr+ierr                                                       2581
      allocate(xm(1:ni),stat=ierr)                                         2581
      jerr=jerr+ierr                                                       2582
      allocate(w(1:no),stat=ierr)                                          2582
      jerr=jerr+ierr                                                       2583
      allocate(qy(1:no),stat=ierr)                                         2583
      jerr=jerr+ierr                                                       2584
      if(jerr.ne.0) return                                                 2585
      bta=parm                                                             2585
      omb=1.0-bta                                                          2585
      fmax=log(huge(bta)*0.1)                                              2586
      qy=q*y                                                               2586
      yb=sum(qy)                                                           2587
      if(nonzero(no,g) .ne. 0)goto 17241                                   2587
      w=q*yb                                                               2587
      az=log(yb)                                                           2587
      uu=az                                                                2588
      xm=yb*xb                                                             2588
      t=0.0                                                                2588
      dv0=yb*(log(yb)-1.0)                                                 2589
      goto 17251                                                           2590
17241 continue                                                             2590
      w=q*exp(sign(min(abs(g),fmax),g))                                    2590
      ww=sum(w)                                                            2590
      eaz=yb/ww                                                            2591
      w=eaz*w                                                              2591
      az=log(eaz)                                                          2591
      uu=az                                                                2591
      t=g                                                                  2591
      dv0=dot_product(qy,g)-yb*(1.0-az)                                    2592
17260 do 17261 j=1,ni                                                      2592
      if(ju(j).eq.0)goto 17261                                             2592
      jb=ix(j)                                                             2592
      je=ix(j+1)-1                                                         2593
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2594
17261 continue                                                             2595
17262 continue                                                             2595
17251 continue                                                             2596
17231 continue                                                             2596
      tt=yb*uu                                                             2596
      ww=yb                                                                2596
      wr=qy-q*(yb*(1.0-uu))                                                2596
      a=0.0                                                                2596
      as=0.0                                                               2597
      dvr=-yb                                                              2598
17270 do 17271 i=1,no                                                      2598
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             2598
17271 continue                                                             2598
17272 continue                                                             2598
      dvr=dvr-dv0                                                          2598
      dev0=dvr                                                             2599
      if(flmin .ge. 1.0)goto 17291                                         2599
      eqs=max(eps,flmin)                                                   2599
      alf=eqs**(1.0/(nlam-1))                                              2599
17291 continue                                                             2600
      m=0                                                                  2600
      mm=0                                                                 2600
      nlp=0                                                                2600
      nin=nlp                                                              2600
      mnl=min(mnlam,nlam)                                                  2600
      shr=shri*dev0                                                        2600
      al=0.0                                                               2600
      ixx=0                                                                2601
17300 do 17301 j=1,ni                                                      2601
      if(ju(j).eq.0)goto 17301                                             2602
      jb=ix(j)                                                             2602
      je=ix(j+1)-1                                                         2603
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   2605 
     *)-xb(j)*tt)/xs(j)
17301 continue                                                             2606
17302 continue                                                             2606
17310 do 17311 ilm=1,nlam                                                  2606
      al0=al                                                               2607
      if(flmin .lt. 1.0)goto 17331                                         2607
      al=ulam(ilm)                                                         2607
      goto 17321                                                           2608
17331 if(ilm .le. 2)goto 17341                                             2608
      al=al*alf                                                            2608
      goto 17321                                                           2609
17341 if(ilm .ne. 1)goto 17351                                             2609
      al=big                                                               2609
      goto 17361                                                           2610
17351 continue                                                             2610
      al0=0.0                                                              2611
17370 do 17371 j=1,ni                                                      2611
      if(ju(j).eq.0)goto 17371                                             2611
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2611
17371 continue                                                             2612
17372 continue                                                             2612
      al0=al0/max(bta,1.0e-3)                                              2612
      al=alf*al0                                                           2613
17361 continue                                                             2614
17321 continue                                                             2614
      al2=al*omb                                                           2614
      al1=al*bta                                                           2614
      tlam=bta*(2.0*al-al0)                                                2615
17380 do 17381 k=1,ni                                                      2615
      if(ixx(k).eq.1)goto 17381                                            2615
      if(ju(k).eq.0)goto 17381                                             2616
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2617
17381 continue                                                             2618
17382 continue                                                             2618
10680 continue                                                             2619
17390 continue                                                             2619
17391 continue                                                             2619
      az0=az                                                               2620
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2621
17400 do 17401 j=1,ni                                                      2621
      if(ixx(j).eq.0)goto 17401                                            2621
      jb=ix(j)                                                             2621
      je=ix(j+1)-1                                                         2622
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2623
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   2625 
     *b(j)**2)/xs(j)**2
17401 continue                                                             2626
17402 continue                                                             2626
17410 continue                                                             2626
17411 continue                                                             2626
      nlp=nlp+1                                                            2627
      dlx=0.0                                                              2628
17420 do 17421 k=1,ni                                                      2628
      if(ixx(k).eq.0)goto 17421                                            2628
      jb=ix(k)                                                             2628
      je=ix(k+1)-1                                                         2628
      ak=a(k)                                                              2629
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2631 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2632
      if(au .gt. 0.0)goto 17441                                            2632
      a(k)=0.0                                                             2632
      goto 17451                                                           2633
17441 continue                                                             2633
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2633
17451 continue                                                             2634
17431 continue                                                             2634
      if(a(k).eq.ak)goto 17421                                             2635
      if(mm(k) .ne. 0)goto 17471                                           2635
      nin=nin+1                                                            2635
      if(nin.gt.nx)goto 17422                                              2636
      mm(k)=nin                                                            2636
      m(nin)=k                                                             2637
17471 continue                                                             2638
      d=a(k)-ak                                                            2638
      dlx=max(dlx,v(k)*d**2)                                               2638
      dv=d/xs(k)                                                           2639
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2640
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2641
      uu=uu-dv*xb(k)                                                       2641
      tt=tt-dv*xm(k)                                                       2642
17421 continue                                                             2643
17422 continue                                                             2643
      if(nin.gt.nx)goto 17412                                              2643
      d=tt/ww-uu                                                           2644
      az=az+d                                                              2644
      dlx=max(dlx,ww*d**2)                                                 2644
      uu=uu+d                                                              2645
      if(dlx.lt.shr)goto 17412                                             2645
      if(nlp .le. maxit)goto 17491                                         2645
      jerr=-ilm                                                            2645
      return                                                               2645
17491 continue                                                             2646
17500 continue                                                             2646
17501 continue                                                             2646
      nlp=nlp+1                                                            2646
      dlx=0.0                                                              2647
17510 do 17511 l=1,nin                                                     2647
      k=m(l)                                                               2648
      jb=ix(k)                                                             2648
      je=ix(k+1)-1                                                         2648
      ak=a(k)                                                              2649
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2651 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2652
      if(au .gt. 0.0)goto 17531                                            2652
      a(k)=0.0                                                             2652
      goto 17541                                                           2653
17531 continue                                                             2653
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2653
17541 continue                                                             2654
17521 continue                                                             2654
      if(a(k).eq.ak)goto 17511                                             2654
      d=a(k)-ak                                                            2654
      dlx=max(dlx,v(k)*d**2)                                               2655
      dv=d/xs(k)                                                           2655
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2656
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2657
      uu=uu-dv*xb(k)                                                       2657
      tt=tt-dv*xm(k)                                                       2658
17511 continue                                                             2659
17512 continue                                                             2659
      d=tt/ww-uu                                                           2659
      az=az+d                                                              2659
      dlx=max(dlx,ww*d**2)                                                 2659
      uu=uu+d                                                              2660
      if(dlx.lt.shr)goto 17502                                             2660
      if(nlp .le. maxit)goto 17561                                         2660
      jerr=-ilm                                                            2660
      return                                                               2660
17561 continue                                                             2661
      goto 17501                                                           2662
17502 continue                                                             2662
      goto 17411                                                           2663
17412 continue                                                             2663
      if(nin.gt.nx)goto 17392                                              2664
      euu=exp(sign(min(abs(uu),fmax),uu))                                  2665
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                2665
      ww=sum(w)                                                            2666
      wr=qy-w*(1.0-uu)                                                     2666
      tt=sum(wr)                                                           2667
      if(ww*(az-az0)**2 .ge. shr)goto 17581                                2667
      kx=0                                                                 2668
17590 do 17591 j=1,nin                                                     2668
      k=m(j)                                                               2669
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 17591                            2669
      kx=1                                                                 2669
      goto 17592                                                           2670
17591 continue                                                             2671
17592 continue                                                             2671
      if(kx .ne. 0)goto 17611                                              2672
17620 do 17621 j=1,ni                                                      2672
      if(ixx(j).eq.1)goto 17621                                            2672
      if(ju(j).eq.0)goto 17621                                             2673
      jb=ix(j)                                                             2673
      je=ix(j+1)-1                                                         2674
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2675
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   2677 
     *)-xb(j)*tt)/xs(j)
      if(ga(j) .le. al1*vp(j))goto 17641                                   2677
      ixx(j)=1                                                             2677
      kx=1                                                                 2677
17641 continue                                                             2678
17621 continue                                                             2679
17622 continue                                                             2679
      if(kx.eq.1) go to 10680                                              2680
      goto 17392                                                           2681
17611 continue                                                             2682
17581 continue                                                             2683
      goto 17391                                                           2684
17392 continue                                                             2684
      if(nin .le. nx)goto 17661                                            2684
      jerr=-10000-ilm                                                      2684
      goto 17312                                                           2684
17661 continue                                                             2685
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2685
      kin(ilm)=nin                                                         2686
      a0(ilm)=az                                                           2686
      alm(ilm)=al                                                          2686
      lmu=ilm                                                              2687
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        2688
      if(ilm.lt.mnl)goto 17311                                             2688
      if(flmin.ge.1.0)goto 17311                                           2689
      me=0                                                                 2689
17670 do 17671 j=1,nin                                                     2689
      if(ca(j,ilm).ne.0.0) me=me+1                                         2689
17671 continue                                                             2689
17672 continue                                                             2689
      if(me.gt.ne)goto 17312                                               2690
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17312              2691
      if(dev(ilm).gt.devmax)goto 17312                                     2692
17311 continue                                                             2693
17312 continue                                                             2693
      g=t+uu                                                               2694
11790 continue                                                             2694
      deallocate(t,w,wr,v,a,qy,xm,as,mm,ga,ixx)                            2695
      return                                                               2696
      end                                                                  2697
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       2698
      real x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)           2699
      integer ix(*),jx(*)                                                  2700
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 17691                                     2703
      jerr=8888                                                            2703
      return                                                               2703
17691 continue                                                             2704
      allocate(w(1:no),stat=jerr)                                          2705
      allocate(f(1:no),stat=ierr)                                          2705
      jerr=jerr+ierr                                                       2706
      if(jerr.ne.0) return                                                 2707
      w=max(0.0,q)                                                         2707
      sw=sum(w)                                                            2707
      if(sw .gt. 0.0)goto 17711                                            2707
      jerr=9999                                                            2707
      go to 11790                                                          2707
17711 continue                                                             2708
      yb=dot_product(w,y)/sw                                               2708
      fmax=log(huge(y(1))*0.1)                                             2709
17720 do 17721 lam=1,nlam                                                  2709
      f=a0(lam)                                                            2710
17730 do 17731 j=1,ni                                                      2710
      if(a(j,lam).eq.0.0)goto 17731                                        2710
      jb=ix(j)                                                             2710
      je=ix(j+1)-1                                                         2711
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          2712
17731 continue                                                             2713
17732 continue                                                             2713
      f=f+g                                                                2714
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2715
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2716
17721 continue                                                             2717
17722 continue                                                             2717
11790 continue                                                             2717
      deallocate(w,f)                                                      2718
      return                                                               2719
      end                                                                  2720
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   2721 
     *jerr)
      real x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(nlam)          2722
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 2723
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 17751                                     2726
      jerr=8888                                                            2726
      return                                                               2726
17751 continue                                                             2727
      allocate(w(1:no),stat=jerr)                                          2728
      allocate(f(1:no),stat=ierr)                                          2728
      jerr=jerr+ierr                                                       2729
      if(jerr.ne.0) return                                                 2730
      w=max(0.0,q)                                                         2730
      sw=sum(w)                                                            2730
      if(sw .gt. 0.0)goto 17771                                            2730
      jerr=9999                                                            2730
      go to 11790                                                          2730
17771 continue                                                             2731
      yb=dot_product(w,y)/sw                                               2731
      fmax=log(huge(y(1))*0.1)                                             2732
17780 do 17781 lam=1,nlam                                                  2732
      f=a0(lam)                                                            2733
17790 do 17791 k=1,nin(lam)                                                2733
      j=ia(k)                                                              2733
      jb=ix(j)                                                             2733
      je=ix(j+1)-1                                                         2734
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         2735
17791 continue                                                             2736
17792 continue                                                             2736
      f=f+g                                                                2737
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2738
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2739
17781 continue                                                             2740
17782 continue                                                             2740
11790 continue                                                             2740
      deallocate(w,f)                                                      2741
      return                                                               2742
      end                                                                  2743
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
