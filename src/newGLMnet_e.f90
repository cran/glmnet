c
c                          newGLMnet (4/18/11)
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
      xmz=vi                                                               1296
      goto 12061                                                           1297
12051 continue                                                             1297
      b(0)=azero(no,y,g,w,jerr)                                            1297
      if(jerr.ne.0) return                                                 1298
      q=1.0/(1.0+exp(-b(0)-g))                                             1298
      v=w*q*(1.0-q)                                                        1298
      r=w*(y-q)                                                            1298
      xmz=sum(v)                                                           1299
12061 continue                                                             1300
12041 continue                                                             1300
      if(kopt .le. 0)goto 12081                                            1301
      if(isd .le. 0)goto 12101                                             1301
      xv=0.25                                                              1301
      goto 12111                                                           1302
12101 continue                                                             1302
12120 do 12121 j=1,ni                                                      1302
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1302
12121 continue                                                             1302
12122 continue                                                             1302
12111 continue                                                             1303
12091 continue                                                             1303
12081 continue                                                             1304
      dev1=-(bz*q0+log(1.0-q0))                                            1304
      dev0=dev1                                                            1305
12130 do 12131 i=1,no                                                      1305
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1306
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1307
12131 continue                                                             1308
12132 continue                                                             1308
      if(flmin .ge. 1.0)goto 12151                                         1308
      eqs=max(eps,flmin)                                                   1308
      alf=eqs**(1.0/(nlam-1))                                              1308
12151 continue                                                             1309
      m=0                                                                  1309
      mm=0                                                                 1309
      nlp=0                                                                1309
      nin=nlp                                                              1309
      mnl=min(mnlam,nlam)                                                  1309
      bs=0.0                                                               1309
      b(1:ni)=0.0                                                          1310
      shr=shri*dev0                                                        1311
12160 do 12161 j=1,ni                                                      1311
      if(ju(j).eq.0)goto 12161                                             1311
      ga(j)=abs(dot_product(r,x(:,j)))                                     1311
12161 continue                                                             1312
12162 continue                                                             1312
12170 do 12171 ilm=1,nlam                                                  1312
      al0=al                                                               1313
      if(flmin .lt. 1.0)goto 12191                                         1313
      al=ulam(ilm)                                                         1313
      goto 12181                                                           1314
12191 if(ilm .le. 2)goto 12201                                             1314
      al=al*alf                                                            1314
      goto 12181                                                           1315
12201 if(ilm .ne. 1)goto 12211                                             1315
      al=big                                                               1315
      goto 12221                                                           1316
12211 continue                                                             1316
      al0=0.0                                                              1317
12230 do 12231 j=1,ni                                                      1317
      if(ju(j).eq.0)goto 12231                                             1317
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1317
12231 continue                                                             1318
12232 continue                                                             1318
      al0=al0/max(bta,1.0e-3)                                              1318
      al=alf*al0                                                           1319
12221 continue                                                             1320
12181 continue                                                             1320
      al2=al*omb                                                           1320
      al1=al*bta                                                           1320
      tlam=bta*(2.0*al-al0)                                                1321
12240 do 12241 k=1,ni                                                      1321
      if(ixx(k).eq.1)goto 12241                                            1321
      if(ju(k).eq.0)goto 12241                                             1322
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1323
12241 continue                                                             1324
12242 continue                                                             1324
10680 continue                                                             1325
12250 continue                                                             1325
12251 continue                                                             1325
      bs(0)=b(0)                                                           1325
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1326
      if(kopt .ne. 0)goto 12271                                            1327
12280 do 12281 j=1,ni                                                      1327
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                       1327
12281 continue                                                             1328
12282 continue                                                             1328
12271 continue                                                             1329
12290 continue                                                             1329
12291 continue                                                             1329
      nlp=nlp+1                                                            1329
      dlx=0.0                                                              1330
12300 do 12301 k=1,ni                                                      1330
      if(ixx(k).eq.0)goto 12301                                            1331
      bk=b(k)                                                              1331
      gk=dot_product(r,x(:,k))                                             1332
      u=gk+xv(k)*b(k)                                                      1332
      au=abs(u)-vp(k)*al1                                                  1333
      if(au .gt. 0.0)goto 12321                                            1333
      b(k)=0.0                                                             1333
      goto 12331                                                           1334
12321 continue                                                             1334
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1334
12331 continue                                                             1335
12311 continue                                                             1335
      d=b(k)-bk                                                            1335
      if(abs(d).le.0.0)goto 12301                                          1335
      dlx=max(dlx,xv(k)*d**2)                                              1336
      r=r-d*v*x(:,k)                                                       1337
      if(mm(k) .ne. 0)goto 12351                                           1337
      nin=nin+1                                                            1337
      if(nin.gt.nx)goto 12302                                              1338
      mm(k)=nin                                                            1338
      m(nin)=k                                                             1339
12351 continue                                                             1340
12301 continue                                                             1341
12302 continue                                                             1341
      if(nin.gt.nx)goto 12292                                              1342
      d=sum(r)/xmz                                                         1343
      if(d .eq. 0.0)goto 12371                                             1343
      b(0)=b(0)+d                                                          1343
      dlx=max(dlx,xmz*d**2)                                                1343
      r=r-d*v                                                              1343
12371 continue                                                             1344
      if(dlx.lt.shr)goto 12292                                             1344
      if(nlp .le. maxit)goto 12391                                         1344
      jerr=-ilm                                                            1344
      return                                                               1344
12391 continue                                                             1345
12400 continue                                                             1345
12401 continue                                                             1345
      nlp=nlp+1                                                            1345
      dlx=0.0                                                              1346
12410 do 12411 l=1,nin                                                     1346
      k=m(l)                                                               1346
      bk=b(k)                                                              1347
      gk=dot_product(r,x(:,k))                                             1348
      u=gk+xv(k)*b(k)                                                      1348
      au=abs(u)-vp(k)*al1                                                  1349
      if(au .gt. 0.0)goto 12431                                            1349
      b(k)=0.0                                                             1349
      goto 12441                                                           1350
12431 continue                                                             1350
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1350
12441 continue                                                             1351
12421 continue                                                             1351
      d=b(k)-bk                                                            1351
      if(abs(d).le.0.0)goto 12411                                          1351
      dlx=max(dlx,xv(k)*d**2)                                              1352
      r=r-d*v*x(:,k)                                                       1353
12411 continue                                                             1354
12412 continue                                                             1354
      d=sum(r)/xmz                                                         1355
      if(d .eq. 0.0)goto 12461                                             1355
      b(0)=b(0)+d                                                          1355
      dlx=max(dlx,xmz*d**2)                                                1355
      r=r-d*v                                                              1355
12461 continue                                                             1356
      if(dlx.lt.shr)goto 12402                                             1356
      if(nlp .le. maxit)goto 12481                                         1356
      jerr=-ilm                                                            1356
      return                                                               1356
12481 continue                                                             1357
      goto 12401                                                           1358
12402 continue                                                             1358
      goto 12291                                                           1359
12292 continue                                                             1359
      if(nin.gt.nx)goto 12252                                              1360
12490 do 12491 i=1,no                                                      1360
      fi=b(0)+g(i)                                                         1361
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1362
      if(fi .ge. fmin)goto 12511                                           1362
      q(i)=0.0                                                             1362
      goto 12501                                                           1362
12511 if(fi .le. fmax)goto 12521                                           1362
      q(i)=1.0                                                             1362
      goto 12531                                                           1363
12521 continue                                                             1363
      q(i)=1.0/(1.0+exp(-fi))                                              1363
12531 continue                                                             1364
12501 continue                                                             1364
12491 continue                                                             1365
12492 continue                                                             1365
      v=w*q*(1.0-q)                                                        1365
      xmz=sum(v)                                                           1365
      if(xmz.le.vmin)goto 12252                                            1365
      r=w*(y-q)                                                            1366
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 12551                           1366
      ix=0                                                                 1367
12560 do 12561 j=1,nin                                                     1367
      k=m(j)                                                               1368
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 12561                           1368
      ix=1                                                                 1368
      goto 12562                                                           1369
12561 continue                                                             1370
12562 continue                                                             1370
      if(ix .ne. 0)goto 12581                                              1371
12590 do 12591 k=1,ni                                                      1371
      if(ixx(k).eq.1)goto 12591                                            1371
      if(ju(k).eq.0)goto 12591                                             1372
      ga(k)=abs(dot_product(r,x(:,k)))                                     1373
      if(ga(k) .le. al1*vp(k))goto 12611                                   1373
      ixx(k)=1                                                             1373
      ix=1                                                                 1373
12611 continue                                                             1374
12591 continue                                                             1375
12592 continue                                                             1375
      if(ix.eq.1) go to 10680                                              1376
      goto 12252                                                           1377
12581 continue                                                             1378
12551 continue                                                             1379
      goto 12251                                                           1380
12252 continue                                                             1380
      if(nin .le. nx)goto 12631                                            1380
      jerr=-10000-ilm                                                      1380
      goto 12172                                                           1380
12631 continue                                                             1381
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1381
      kin(ilm)=nin                                                         1382
      a0(ilm)=b(0)                                                         1382
      alm(ilm)=al                                                          1382
      lmu=ilm                                                              1383
      devi=dev2(no,w,y,q,pmin)                                             1384
      dev(ilm)=(dev1-devi)/dev0                                            1384
      if(xmz.le.vmin)goto 12172                                            1385
      if(ilm.lt.mnl)goto 12171                                             1385
      if(flmin.ge.1.0)goto 12171                                           1386
      me=0                                                                 1386
12640 do 12641 j=1,nin                                                     1386
      if(a(j,ilm).ne.0.0) me=me+1                                          1386
12641 continue                                                             1386
12642 continue                                                             1386
      if(me.gt.ne)goto 12172                                               1387
      if(dev(ilm).gt.devmax)goto 12172                                     1387
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12172                             1388
12171 continue                                                             1389
12172 continue                                                             1389
      g=log(q/(1.0-q))                                                     1390
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                                  1391
      return                                                               1392
      end                                                                  1393
      function dev2(n,w,y,p,pmin)                                          1394
      real w(n),y(n),p(n)                                                  1395
      pmax=1.0-pmin                                                        1395
      s=0.0                                                                1396
12650 do 12651 i=1,n                                                       1396
      pi=min(max(pmin,p(i)),pmax)                                          1397
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1398
12651 continue                                                             1399
12652 continue                                                             1399
      dev2=s                                                               1400
      return                                                               1401
      end                                                                  1402
      function azero(n,y,g,q,jerr)                                         1403
      parameter(eps=1.0e-7)                                                1404
      real y(n),g(n),q(n)                                                  1405
      real, dimension (:), allocatable :: e,p,w                                 
      allocate(e(1:n),stat=jerr)                                           1409
      allocate(p(1:n),stat=ierr)                                           1409
      jerr=jerr+ierr                                                       1410
      allocate(w(1:n),stat=ierr)                                           1410
      jerr=jerr+ierr                                                       1411
      if(jerr.ne.0) return                                                 1412
      az=0.0                                                               1412
      e=exp(-g)                                                            1412
      qy=dot_product(q,y)                                                  1412
      p=1.0/(1.0+e)                                                        1413
12660 continue                                                             1413
12661 continue                                                             1413
      w=q*p*(1.0-p)                                                        1414
      d=(qy-dot_product(q,p))/sum(w)                                       1414
      az=az+d                                                              1414
      if(abs(d).lt.eps)goto 12662                                          1415
      ea0=exp(-az)                                                         1415
      p=1.0/(1.0+ea0*e)                                                    1416
      goto 12661                                                           1417
12662 continue                                                             1417
      azero=az                                                             1418
      deallocate(e,p,w)                                                    1419
      return                                                               1420
      end                                                                  1421
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ul   1423 
     *am,shri,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1425 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam)              1426
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1427
      integer ju(ni),m(nx),kin(nlam)                                       1428
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp                       
      real, dimension (:), allocatable :: di,v,r,ga                             
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is,ixx                          
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no),stat=ierr)                                          1439
      jerr=jerr+ierr                                                       1440
      allocate(v(1:no),stat=ierr)                                          1440
      jerr=jerr+ierr                                                       1441
      allocate(mm(1:ni),stat=ierr)                                         1441
      jerr=jerr+ierr                                                       1442
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1442
      jerr=jerr+ierr                                                       1443
      allocate(sxp(1:no),stat=ierr)                                        1443
      jerr=jerr+ierr                                                       1444
      allocate(di(1:no),stat=ierr)                                         1444
      jerr=jerr+ierr                                                       1445
      allocate(ga(1:ni),stat=ierr)                                         1445
      jerr=jerr+ierr                                                       1446
      allocate(ixx(1:ni),stat=ierr)                                        1446
      jerr=jerr+ierr                                                       1447
      if(jerr.ne.0) return                                                 1448
      pmax=1.0-pmin                                                        1448
      emin=pmin/pmax                                                       1448
      emax=1.0/emin                                                        1449
      pfm=(1.0+pmin)*pmin                                                  1449
      pfx=(1.0-pmin)*pmax                                                  1449
      vmin=pfm*pmax                                                        1450
      bta=parm                                                             1450
      omb=1.0-bta                                                          1450
      dev1=0.0                                                             1450
      dev0=0.0                                                             1451
12670 do 12671 ic=1,nc                                                     1451
      q0=dot_product(w,y(:,ic))                                            1452
      if(q0 .gt. pmin)goto 12691                                           1452
      jerr =8000+ic                                                        1452
      return                                                               1452
12691 continue                                                             1453
      if(q0 .lt. 1.0-pmin)goto 12711                                       1453
      jerr =9000+ic                                                        1453
      return                                                               1453
12711 continue                                                             1454
      b(0,ic)=log(q0)                                                      1454
      dev1=dev1-q0*b(0,ic)                                                 1454
      b(1:ni,ic)=0.0                                                       1455
12720 do 12721 i=1,no                                                      1455
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1455
12721 continue                                                             1456
12722 continue                                                             1456
12671 continue                                                             1457
12672 continue                                                             1457
      dev0=dev0+dev1                                                       1457
      ixx=0                                                                1457
      al=0.0                                                               1458
      if(nonzero(no*nc,g) .ne. 0)goto 12741                                1459
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1459
      sxp=0.0                                                              1460
12750 do 12751 ic=1,nc                                                     1460
      q(:,ic)=exp(b(0,ic))                                                 1460
      sxp=sxp+q(:,ic)                                                      1460
12751 continue                                                             1461
12752 continue                                                             1461
      goto 12761                                                           1462
12741 continue                                                             1462
12770 do 12771 i=1,no                                                      1462
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1462
12771 continue                                                             1462
12772 continue                                                             1462
      sxp=0.0                                                              1463
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1463
      if(jerr.ne.0) return                                                 1464
12780 do 12781 ic=1,nc                                                     1464
      q(:,ic)=exp(b(0,ic)+g(:,ic))                                         1464
      sxp=sxp+q(:,ic)                                                      1464
12781 continue                                                             1465
12782 continue                                                             1465
12761 continue                                                             1466
12731 continue                                                             1466
      if(kopt .le. 0)goto 12801                                            1467
      if(isd .le. 0)goto 12821                                             1467
      xv=0.25                                                              1467
      goto 12831                                                           1468
12821 continue                                                             1468
12840 do 12841 j=1,ni                                                      1468
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1468
12841 continue                                                             1468
12842 continue                                                             1468
12831 continue                                                             1469
12811 continue                                                             1469
12801 continue                                                             1470
      if(flmin .ge. 1.0)goto 12861                                         1470
      eqs=max(eps,flmin)                                                   1470
      alf=eqs**(1.0/(nlam-1))                                              1470
12861 continue                                                             1471
      m=0                                                                  1471
      mm=0                                                                 1471
      nin=0                                                                1471
      nlp=0                                                                1471
      mnl=min(mnlam,nlam)                                                  1471
      bs=0.0                                                               1471
      shr=shri*dev0                                                        1472
      ga=0.0                                                               1473
12870 do 12871 ic=1,nc                                                     1473
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1474
12880 do 12881 j=1,ni                                                      1474
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))           1474
12881 continue                                                             1475
12882 continue                                                             1475
12871 continue                                                             1476
12872 continue                                                             1476
12890 do 12891 ilm=1,nlam                                                  1476
      al0=al                                                               1477
      if(flmin .lt. 1.0)goto 12911                                         1477
      al=ulam(ilm)                                                         1477
      goto 12901                                                           1478
12911 if(ilm .le. 2)goto 12921                                             1478
      al=al*alf                                                            1478
      goto 12901                                                           1479
12921 if(ilm .ne. 1)goto 12931                                             1479
      al=big                                                               1479
      goto 12941                                                           1480
12931 continue                                                             1480
      al0=0.0                                                              1481
12950 do 12951 j=1,ni                                                      1481
      if(ju(j).eq.0)goto 12951                                             1481
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1481
12951 continue                                                             1482
12952 continue                                                             1482
      al0=al0/max(bta,1.0e-3)                                              1482
      al=alf*al0                                                           1483
12941 continue                                                             1484
12901 continue                                                             1484
      al2=al*omb                                                           1484
      al1=al*bta                                                           1484
      tlam=bta*(2.0*al-al0)                                                1485
12960 do 12961 k=1,ni                                                      1485
      if(ixx(k).eq.1)goto 12961                                            1485
      if(ju(k).eq.0)goto 12961                                             1486
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1487
12961 continue                                                             1488
12962 continue                                                             1488
10680 continue                                                             1489
12970 continue                                                             1489
12971 continue                                                             1489
      ix=0                                                                 1489
      jx=ix                                                                1489
      ig=0                                                                 1490
12980 do 12981 ic=1,nc                                                     1490
      bs(0,ic)=b(0,ic)                                                     1491
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1492
      xmz=0.0                                                              1493
12990 do 12991 i=1,no                                                      1493
      pic=q(i,ic)/sxp(i)                                                   1494
      if(pic .ge. pfm)goto 13011                                           1494
      pic=0.0                                                              1494
      v(i)=0.0                                                             1494
      goto 13001                                                           1495
13011 if(pic .le. pfx)goto 13021                                           1495
      pic=1.0                                                              1495
      v(i)=0.0                                                             1495
      goto 13031                                                           1496
13021 continue                                                             1496
      v(i)=w(i)*pic*(1.0-pic)                                              1496
      xmz=xmz+v(i)                                                         1496
13031 continue                                                             1497
13001 continue                                                             1497
      r(i)=w(i)*(y(i,ic)-pic)                                              1498
12991 continue                                                             1499
12992 continue                                                             1499
      if(xmz.le.vmin)goto 12981                                            1499
      ig=1                                                                 1500
      if(kopt .ne. 0)goto 13051                                            1501
13060 do 13061 j=1,ni                                                      1501
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                    1501
13061 continue                                                             1502
13062 continue                                                             1502
13051 continue                                                             1503
13070 continue                                                             1503
13071 continue                                                             1503
      nlp=nlp+1                                                            1503
      dlx=0.0                                                              1504
13080 do 13081 k=1,ni                                                      1504
      if(ixx(k).eq.0)goto 13081                                            1505
      bk=b(k,ic)                                                           1505
      gk=dot_product(r,x(:,k))                                             1506
      u=gk+xv(k,ic)*b(k,ic)                                                1506
      au=abs(u)-vp(k)*al1                                                  1507
      if(au .gt. 0.0)goto 13101                                            1507
      b(k,ic)=0.0                                                          1507
      goto 13111                                                           1508
13101 continue                                                             1508
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1508
13111 continue                                                             1509
13091 continue                                                             1509
      d=b(k,ic)-bk                                                         1509
      if(abs(d).le.0.0)goto 13081                                          1510
      dlx=max(dlx,xv(k,ic)*d**2)                                           1510
      r=r-d*v*x(:,k)                                                       1511
      if(mm(k) .ne. 0)goto 13131                                           1511
      nin=nin+1                                                            1512
      if(nin .le. nx)goto 13151                                            1512
      jx=1                                                                 1512
      goto 13082                                                           1512
13151 continue                                                             1513
      mm(k)=nin                                                            1513
      m(nin)=k                                                             1514
13131 continue                                                             1515
13081 continue                                                             1516
13082 continue                                                             1516
      if(jx.gt.0)goto 13072                                                1517
      d=sum(r)/xmz                                                         1518
      if(d .eq. 0.0)goto 13171                                             1518
      b(0,ic)=b(0,ic)+d                                                    1518
      dlx=max(dlx,xmz*d**2)                                                1518
      r=r-d*v                                                              1518
13171 continue                                                             1519
      if(dlx.lt.shr)goto 13072                                             1520
      if(nlp .le. maxit)goto 13191                                         1520
      jerr=-ilm                                                            1520
      return                                                               1520
13191 continue                                                             1521
13200 continue                                                             1521
13201 continue                                                             1521
      nlp=nlp+1                                                            1521
      dlx=0.0                                                              1522
13210 do 13211 l=1,nin                                                     1522
      k=m(l)                                                               1522
      bk=b(k,ic)                                                           1523
      gk=dot_product(r,x(:,k))                                             1524
      u=gk+xv(k,ic)*b(k,ic)                                                1524
      au=abs(u)-vp(k)*al1                                                  1525
      if(au .gt. 0.0)goto 13231                                            1525
      b(k,ic)=0.0                                                          1525
      goto 13241                                                           1526
13231 continue                                                             1526
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1526
13241 continue                                                             1527
13221 continue                                                             1527
      d=b(k,ic)-bk                                                         1527
      if(abs(d).le.0.0)goto 13211                                          1528
      dlx=max(dlx,xv(k,ic)*d**2)                                           1528
      r=r-d*v*x(:,k)                                                       1529
13211 continue                                                             1530
13212 continue                                                             1530
      d=sum(r)/xmz                                                         1531
      if(d .eq. 0.0)goto 13261                                             1531
      b(0,ic)=b(0,ic)+d                                                    1532
      dlx=max(dlx,xmz*d**2)                                                1532
      r=r-d*v                                                              1533
13261 continue                                                             1534
      if(dlx.lt.shr)goto 13202                                             1534
      if(nlp .le. maxit)goto 13281                                         1534
      jerr=-ilm                                                            1534
      return                                                               1534
13281 continue                                                             1535
      goto 13201                                                           1536
13202 continue                                                             1536
      goto 13071                                                           1537
13072 continue                                                             1537
      if(jx.gt.0)goto 12982                                                1538
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                            1539
      if(ix .ne. 0)goto 13301                                              1540
13310 do 13311 j=1,nin                                                     1540
      k=m(j)                                                               1541
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 13331                1541
      ix=1                                                                 1541
      goto 13312                                                           1541
13331 continue                                                             1542
13311 continue                                                             1543
13312 continue                                                             1543
13301 continue                                                             1544
13340 do 13341 i=1,no                                                      1544
      fi=b(0,ic)+g(i,ic)                                                   1546
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1547
      fi=min(max(exmn,fi),exmx)                                            1547
      sxp(i)=sxp(i)-q(i,ic)                                                1548
      q(i,ic)=min(max(emin*sxp(i),exp(dble(fi))),emax*sxp(i))              1549
      sxp(i)=sxp(i)+q(i,ic)                                                1550
13341 continue                                                             1551
13342 continue                                                             1551
12981 continue                                                             1552
12982 continue                                                             1552
      s=-sum(b(0,:))/nc                                                    1552
      b(0,:)=b(0,:)+s                                                      1552
      di=s                                                                 1553
13350 do 13351 j=1,nin                                                     1553
      l=m(j)                                                               1554
      if(vp(l) .gt. 0.0)goto 13371                                         1554
      s=sum(b(l,:))/nc                                                     1554
      goto 13381                                                           1555
13371 continue                                                             1555
      s=elc(parm,nc,b(l,:),is)                                             1555
13381 continue                                                             1556
13361 continue                                                             1556
      b(l,:)=b(l,:)-s                                                      1556
      di=di-s*x(:,l)                                                       1557
13351 continue                                                             1558
13352 continue                                                             1558
      di=exp(di)                                                           1558
      sxp=sxp*di                                                           1558
13390 do 13391 ic=1,nc                                                     1558
      q(:,ic)=q(:,ic)*di                                                   1558
13391 continue                                                             1559
13392 continue                                                             1559
      if(jx.gt.0)goto 12972                                                1559
      if(ig.eq.0)goto 12972                                                1560
      if(ix .ne. 0)goto 13411                                              1560
      ga=0.0                                                               1561
13420 do 13421 ic=1,nc                                                     1561
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1562
13430 do 13431 k=1,ni                                                      1562
      if(ixx(k).eq.1)goto 13431                                            1562
      if(ju(k).eq.0)goto 13431                                             1563
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                          1564
13431 continue                                                             1565
13432 continue                                                             1565
13421 continue                                                             1566
13422 continue                                                             1566
13440 do 13441 k=1,ni                                                      1566
      if(ixx(k).eq.1)goto 13441                                            1566
      if(ju(k).eq.0)goto 13441                                             1567
      if(ga(k) .le. al1*vp(k))goto 13461                                   1567
      ixx(k)=1                                                             1567
      ix=1                                                                 1567
13461 continue                                                             1568
13441 continue                                                             1569
13442 continue                                                             1569
      if(ix.eq.1) go to 10680                                              1570
      goto 12972                                                           1571
13411 continue                                                             1572
      goto 12971                                                           1573
12972 continue                                                             1573
      if(jx .le. 0)goto 13481                                              1573
      jerr=-10000-ilm                                                      1573
      goto 12892                                                           1573
13481 continue                                                             1573
      devi=0.0                                                             1574
13490 do 13491 ic=1,nc                                                     1575
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1575
      a0(ic,ilm)=b(0,ic)                                                   1576
13500 do 13501 i=1,no                                                      1576
      if(y(i,ic).le.0.0)goto 13501                                         1577
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1578
13501 continue                                                             1579
13502 continue                                                             1579
13491 continue                                                             1580
13492 continue                                                             1580
      kin(ilm)=nin                                                         1580
      alm(ilm)=al                                                          1580
      lmu=ilm                                                              1581
      dev(ilm)=(dev1-devi)/dev0                                            1581
      if(ig.eq.0)goto 12892                                                1582
      if(ilm.lt.mnl)goto 12891                                             1582
      if(flmin.ge.1.0)goto 12891                                           1583
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 12892             1584
      if(dev(ilm).gt.devmax)goto 12892                                     1584
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12892                             1585
12891 continue                                                             1586
12892 continue                                                             1586
      g=log(q)                                                             1586
13510 do 13511 i=1,no                                                      1586
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1586
13511 continue                                                             1587
13512 continue                                                             1587
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                           1588
      return                                                               1589
      end                                                                  1590
      subroutine kazero(kk,n,y,g,q,az,jerr)                                1591
      parameter(eps=1.0e-7)                                                1592
      real y(n,kk),g(n,kk),q(n),az(kk)                                     1593
      real, dimension (:), allocatable :: s                                     
      real, dimension (:,:), allocatable :: e                                   
      allocate(e(1:n,1:kk),stat=jerr)                                           
      allocate(s(1:n),stat=ierr)                                           1598
      jerr=jerr+ierr                                                       1599
      if(jerr.ne.0) return                                                 1600
      az=0.0                                                               1600
      e=exp(g)                                                             1600
13520 do 13521 i=1,n                                                       1600
      s(i)=sum(e(i,:))                                                     1600
13521 continue                                                             1601
13522 continue                                                             1601
13530 continue                                                             1601
13531 continue                                                             1601
      dm=0.0                                                               1602
13540 do 13541 k=1,kk                                                      1602
      t=0.0                                                                1602
      u=t                                                                  1603
13550 do 13551 i=1,n                                                       1603
      pik=e(i,k)/s(i)                                                      1604
      t=t+q(i)*(y(i,k)-pik)                                                1604
      u=u+q(i)*pik*(1.0-pik)                                               1605
13551 continue                                                             1606
13552 continue                                                             1606
      d=t/u                                                                1606
      az(k)=az(k)+d                                                        1606
      ed=exp(d)                                                            1606
      dm=max(dm,abs(d))                                                    1607
13560 do 13561 i=1,n                                                       1607
      z=e(i,k)                                                             1607
      e(i,k)=z*ed                                                          1607
      s(i)=s(i)-z+e(i,k)                                                   1607
13561 continue                                                             1608
13562 continue                                                             1608
13541 continue                                                             1609
13542 continue                                                             1609
      if(dm.lt.eps)goto 13532                                              1609
      goto 13531                                                           1610
13532 continue                                                             1610
      az=az-sum(az)/kk                                                     1611
      deallocate(e,s)                                                      1612
      return                                                               1613
      end                                                                  1614
      function elc(parm,n,a,m)                                             1615
      real a(n)                                                            1615
      integer m(n)                                                         1616
      fn=n                                                                 1616
      am=sum(a)/fn                                                         1617
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 13581                       1617
      elc=am                                                               1617
      return                                                               1617
13581 continue                                                             1618
13590 do 13591 i=1,n                                                       1618
      m(i)=i                                                               1618
13591 continue                                                             1618
13592 continue                                                             1618
      call psort7(a,m,1,n)                                                 1619
      if(a(m(1)) .ne. a(m(n)))goto 13611                                   1619
      elc=a(1)                                                             1619
      return                                                               1619
13611 continue                                                             1620
      if(mod(n,2) .ne. 1)goto 13631                                        1620
      ad=a(m(n/2+1))                                                       1620
      goto 13641                                                           1621
13631 continue                                                             1621
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       1621
13641 continue                                                             1622
13621 continue                                                             1622
      if(parm .ne. 1.0)goto 13661                                          1622
      elc=ad                                                               1622
      return                                                               1622
13661 continue                                                             1623
      b1=min(am,ad)                                                        1623
      b2=max(am,ad)                                                        1623
      k2=1                                                                 1624
13670 continue                                                             1624
13671 if(a(m(k2)).gt.b1)goto 13672                                         1624
      k2=k2+1                                                              1624
      goto 13671                                                           1624
13672 continue                                                             1624
      k1=k2-1                                                              1625
13680 continue                                                             1625
13681 if(a(m(k2)).ge.b2)goto 13682                                         1625
      k2=k2+1                                                              1625
      goto 13681                                                           1626
13682 continue                                                             1626
      r=parm/((1.0-parm)*fn)                                               1626
      is=0                                                                 1626
      sm=n-2*(k1-1)                                                        1627
13690 do 13691 k=k1,k2-1                                                   1627
      sm=sm-2.0                                                            1627
      s=r*sm+am                                                            1628
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 13711                   1628
      is=k                                                                 1628
      goto 13692                                                           1628
13711 continue                                                             1629
13691 continue                                                             1630
13692 continue                                                             1630
      if(is .eq. 0)goto 13731                                              1630
      elc=s                                                                1630
      return                                                               1630
13731 continue                                                             1630
      r2=2.0*r                                                             1630
      s1=a(m(k1))                                                          1630
      am2=2.0*am                                                           1631
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    1631
      elc=s1                                                               1632
13740 do 13741 k=k1+1,k2                                                   1632
      s=a(m(k))                                                            1632
      if(s.eq.s1)goto 13741                                                1633
      c=r2*sum(abs(a-s))+s*(s-am2)                                         1634
      if(c .ge. cri)goto 13761                                             1634
      cri=c                                                                1634
      elc=s                                                                1634
13761 continue                                                             1634
      s1=s                                                                 1635
13741 continue                                                             1636
13742 continue                                                             1636
      return                                                               1637
      end                                                                  1638
      function nintot(ni,nx,nc,a,m,nin,is)                                 1639
      real a(nx,nc)                                                        1639
      integer m(nx),is(ni)                                                 1640
      is=0                                                                 1640
      nintot=0                                                             1641
13770 do 13771 ic=1,nc                                                     1641
13780 do 13781 j=1,nin                                                     1641
      k=m(j)                                                               1641
      if(is(k).ne.0)goto 13781                                             1642
      if(a(j,ic).eq.0.0)goto 13781                                         1642
      is(k)=k                                                              1642
      nintot=nintot+1                                                      1643
13781 continue                                                             1643
13782 continue                                                             1643
13771 continue                                                             1644
13772 continue                                                             1644
      return                                                               1645
      end                                                                  1646
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             1647
      real ca(nx,nc),a(ni,nc)                                              1647
      integer ia(nx)                                                       1648
      a=0.0                                                                1649
13790 do 13791 ic=1,nc                                                     1649
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            1649
13791 continue                                                             1650
13792 continue                                                             1650
      return                                                               1651
      end                                                                  1652
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      1653
      real a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                             1653
      integer ia(nx)                                                       1654
13800 do 13801 i=1,nt                                                      1654
13810 do 13811 ic=1,nc                                                     1654
      ans(ic,i)=a0(ic)                                                     1656
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   1657 
     *:nin)))
13811 continue                                                             1657
13812 continue                                                             1657
13801 continue                                                             1658
13802 continue                                                             1658
      return                                                               1659
      end                                                                  1660
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,ne,nx,nlam,fl   1662 
     *min,  ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      real x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)                 1663
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1664
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1665
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 13831                                    1669
      jerr=10000                                                           1669
      return                                                               1669
13831 continue                                                             1670
      allocate(ww(1:no),stat=jerr)                                         1671
      allocate(ju(1:ni),stat=ierr)                                         1671
      jerr=jerr+ierr                                                       1672
      allocate(vq(1:ni),stat=ierr)                                         1672
      jerr=jerr+ierr                                                       1673
      allocate(xm(1:ni),stat=ierr)                                         1673
      jerr=jerr+ierr                                                       1674
      allocate(xs(1:ni),stat=ierr)                                         1674
      jerr=jerr+ierr                                                       1675
      if(jerr.ne.0) return                                                 1676
      call spchkvars(no,ni,x,ix,ju)                                        1677
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1678
      if(maxval(ju) .gt. 0)goto 13851                                      1678
      jerr=7777                                                            1678
      return                                                               1678
13851 continue                                                             1679
      vq=max(0.0,vp)                                                       1679
      vq=vq*ni/sum(vq)                                                     1680
13860 do 13861 i=1,no                                                      1680
      ww(i)=sum(y(i,:))                                                    1680
      y(i,:)=y(i,:)/ww(i)                                                  1680
13861 continue                                                             1680
13862 continue                                                             1680
      sw=sum(ww)                                                           1680
      ww=ww/sw                                                             1681
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     1682
      if(nc .ne. 1)goto 13881                                              1683
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,ne,nx,n   1685 
     *lam,flmin,  ulam,thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,de
     *v,alm,nlp,jerr)
      goto 13891                                                           1686
13881 continue                                                             1687
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmi   1689 
     *n,ulam,  thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nl
     *p,jerr)
13891 continue                                                             1690
13871 continue                                                             1690
      if(jerr.gt.0) return                                                 1690
      dev0=2.0*sw*dev0                                                     1691
13900 do 13901 k=1,lmu                                                     1691
      nk=nin(k)                                                            1692
13910 do 13911 ic=1,nc                                                     1692
      if(isd .le. 0)goto 13931                                             1692
13940 do 13941 l=1,nk                                                      1692
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1692
13941 continue                                                             1692
13942 continue                                                             1692
13931 continue                                                             1693
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1694
13911 continue                                                             1695
13912 continue                                                             1695
13901 continue                                                             1696
13902 continue                                                             1696
      deallocate(ww,ju,vq,xm,xs)                                           1697
      return                                                               1698
      end                                                                  1699
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,xm,xs)                1700
      real x(*),w(no),xm(ni),xs(ni)                                        1700
      integer ix(*),jx(*),ju(ni)                                           1701
13950 do 13951 j=1,ni                                                      1701
      if(ju(j).eq.0)goto 13951                                             1701
      jb=ix(j)                                                             1701
      je=ix(j+1)-1                                                         1702
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1703
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   1704 
     *)**2)
13951 continue                                                             1705
13952 continue                                                             1705
      if(isd.eq.0) xs=1.0                                                  1706
      return                                                               1707
      end                                                                  1708
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam,   1710 
     *  flmin,ulam,shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,alm
     *,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1712 
     *5, devmax=0.999)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        1713
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1714
      real xb(ni),xs(ni)                                                   1714
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1715
      real, dimension (:), allocatable :: xm,b,bs,v,r,sc,xv,q,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(b(0:ni),stat=jerr)                                          1720
      allocate(xm(0:ni),stat=ierr)                                         1720
      jerr=jerr+ierr                                                       1721
      allocate(xv(1:ni),stat=ierr)                                         1721
      jerr=jerr+ierr                                                       1722
      allocate(bs(0:ni),stat=ierr)                                         1722
      jerr=jerr+ierr                                                       1723
      allocate(ga(1:ni),stat=ierr)                                         1723
      jerr=jerr+ierr                                                       1724
      allocate(mm(1:ni),stat=ierr)                                         1724
      jerr=jerr+ierr                                                       1725
      allocate(ixx(1:ni),stat=ierr)                                        1725
      jerr=jerr+ierr                                                       1726
      allocate(q(1:no),stat=ierr)                                          1726
      jerr=jerr+ierr                                                       1727
      allocate(r(1:no),stat=ierr)                                          1727
      jerr=jerr+ierr                                                       1728
      allocate(v(1:no),stat=ierr)                                          1728
      jerr=jerr+ierr                                                       1729
      allocate(sc(1:no),stat=ierr)                                         1729
      jerr=jerr+ierr                                                       1730
      if(jerr.ne.0) return                                                 1731
      fmax=log(1.0/pmin-1.0)                                               1731
      fmin=-fmax                                                           1731
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1732
      bta=parm                                                             1732
      omb=1.0-bta                                                          1733
      q0=dot_product(w,y)                                                  1733
      if(q0 .gt. pmin)goto 13971                                           1733
      jerr=8001                                                            1733
      return                                                               1733
13971 continue                                                             1734
      if(q0 .lt. 1.0-pmin)goto 13991                                       1734
      jerr=9001                                                            1734
      return                                                               1734
13991 continue                                                             1734
      bz=log(q0/(1.0-q0))                                                  1735
      if(nonzero(no,g) .ne. 0)goto 14011                                   1735
      vi=q0*(1.0-q0)                                                       1735
      b(0)=bz                                                              1735
      v=vi*w                                                               1736
      r=w*(y-q0)                                                           1736
      q=q0                                                                 1736
      xm(0)=vi                                                             1737
      goto 14021                                                           1738
14011 continue                                                             1738
      b(0)=azero(no,y,g,w,jerr)                                            1738
      if(jerr.ne.0) return                                                 1739
      q=1.0/(1.0+exp(-b(0)-g))                                             1739
      v=w*q*(1.0-q)                                                        1739
      r=w*(y-q)                                                            1739
      xm(0)=sum(v)                                                         1740
14021 continue                                                             1741
14001 continue                                                             1741
      if(kopt .le. 0)goto 14041                                            1742
      if(isd .le. 0)goto 14061                                             1742
      xv=0.25                                                              1742
      goto 14071                                                           1743
14061 continue                                                             1744
14080 do 14081 j=1,ni                                                      1744
      if(ju(j).eq.0)goto 14081                                             1744
      jb=ix(j)                                                             1744
      je=ix(j+1)-1                                                         1745
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          1746
14081 continue                                                             1747
14082 continue                                                             1747
14071 continue                                                             1748
14051 continue                                                             1748
14041 continue                                                             1749
      b(1:ni)=0.0                                                          1749
      dev1=-(bz*q0+log(1.0-q0))                                            1749
      dev0=dev1                                                            1750
14090 do 14091 i=1,no                                                      1750
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1751
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1752
14091 continue                                                             1753
14092 continue                                                             1753
      if(flmin .ge. 1.0)goto 14111                                         1753
      eqs=max(eps,flmin)                                                   1753
      alf=eqs**(1.0/(nlam-1))                                              1753
14111 continue                                                             1754
      m=0                                                                  1754
      mm=0                                                                 1754
      nin=0                                                                1754
      o=0.0                                                                1754
      svr=o                                                                1754
      mnl=min(mnlam,nlam)                                                  1754
      bs=0.0                                                               1754
      nlp=0                                                                1754
      nin=nlp                                                              1755
      shr=shri*dev0                                                        1755
      al=0.0                                                               1755
      ixx=0                                                                1756
14120 do 14121 j=1,ni                                                      1756
      if(ju(j).eq.0)goto 14121                                             1757
      jb=ix(j)                                                             1757
      je=ix(j+1)-1                                                         1757
      jn=ix(j+1)-ix(j)                                                     1758
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1759
      gj=dot_product(sc(1:jn),x(jb:je))                                    1760
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      1761
14121 continue                                                             1762
14122 continue                                                             1762
14130 do 14131 ilm=1,nlam                                                  1762
      al0=al                                                               1763
      if(flmin .lt. 1.0)goto 14151                                         1763
      al=ulam(ilm)                                                         1763
      goto 14141                                                           1764
14151 if(ilm .le. 2)goto 14161                                             1764
      al=al*alf                                                            1764
      goto 14141                                                           1765
14161 if(ilm .ne. 1)goto 14171                                             1765
      al=big                                                               1765
      goto 14181                                                           1766
14171 continue                                                             1766
      al0=0.0                                                              1767
14190 do 14191 j=1,ni                                                      1767
      if(ju(j).eq.0)goto 14191                                             1767
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1767
14191 continue                                                             1768
14192 continue                                                             1768
      al0=al0/max(bta,1.0e-3)                                              1768
      al=alf*al0                                                           1769
14181 continue                                                             1770
14141 continue                                                             1770
      al2=al*omb                                                           1770
      al1=al*bta                                                           1770
      tlam=bta*(2.0*al-al0)                                                1771
14200 do 14201 k=1,ni                                                      1771
      if(ixx(k).eq.1)goto 14201                                            1771
      if(ju(k).eq.0)goto 14201                                             1772
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1773
14201 continue                                                             1774
14202 continue                                                             1774
10680 continue                                                             1775
14210 continue                                                             1775
14211 continue                                                             1775
      bs(0)=b(0)                                                           1775
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1776
14220 do 14221 j=1,ni                                                      1776
      if(ixx(j).eq.0)goto 14221                                            1777
      jb=ix(j)                                                             1777
      je=ix(j+1)-1                                                         1777
      jn=ix(j+1)-ix(j)                                                     1778
      sc(1:jn)=v(jx(jb:je))                                                1779
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 1780
      if(kopt .ne. 0)goto 14241                                            1781
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              1782
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                1783
14241 continue                                                             1784
14221 continue                                                             1785
14222 continue                                                             1785
14250 continue                                                             1785
14251 continue                                                             1785
      nlp=nlp+1                                                            1785
      dlx=0.0                                                              1786
14260 do 14261 k=1,ni                                                      1786
      if(ixx(k).eq.0)goto 14261                                            1787
      jb=ix(k)                                                             1787
      je=ix(k+1)-1                                                         1787
      jn=ix(k+1)-ix(k)                                                     1787
      bk=b(k)                                                              1788
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1789
      gk=dot_product(sc(1:jn),x(jb:je))                                    1790
      gk=(gk-svr*xb(k))/xs(k)                                              1791
      u=gk+xv(k)*b(k)                                                      1791
      au=abs(u)-vp(k)*al1                                                  1792
      if(au .gt. 0.0)goto 14281                                            1792
      b(k)=0.0                                                             1792
      goto 14291                                                           1793
14281 continue                                                             1793
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1793
14291 continue                                                             1794
14271 continue                                                             1794
      d=b(k)-bk                                                            1794
      if(abs(d).le.0.0)goto 14261                                          1794
      dlx=max(dlx,xv(k)*d**2)                                              1795
      if(mm(k) .ne. 0)goto 14311                                           1795
      nin=nin+1                                                            1795
      if(nin.gt.nx)goto 14262                                              1796
      mm(k)=nin                                                            1796
      m(nin)=k                                                             1796
      sc(1:jn)=v(jx(jb:je))                                                1797
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 1798
14311 continue                                                             1799
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1800
      o=o+d*(xb(k)/xs(k))                                                  1801
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1802
14261 continue                                                             1803
14262 continue                                                             1803
      if(nin.gt.nx)goto 14252                                              1804
      d=svr/xm(0)                                                          1805
      if(d .eq. 0.0)goto 14331                                             1805
      b(0)=b(0)+d                                                          1805
      dlx=max(dlx,xm(0)*d**2)                                              1805
      r=r-d*v                                                              1805
14331 continue                                                             1806
      svr=svr-d*xm(0)                                                      1806
      if(dlx.lt.shr)goto 14252                                             1807
      if(nlp .le. maxit)goto 14351                                         1807
      jerr=-ilm                                                            1807
      return                                                               1807
14351 continue                                                             1808
14360 continue                                                             1808
14361 continue                                                             1808
      nlp=nlp+1                                                            1808
      dlx=0.0                                                              1809
14370 do 14371 l=1,nin                                                     1809
      k=m(l)                                                               1809
      jb=ix(k)                                                             1809
      je=ix(k+1)-1                                                         1810
      jn=ix(k+1)-ix(k)                                                     1810
      bk=b(k)                                                              1811
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1812
      gk=dot_product(sc(1:jn),x(jb:je))                                    1813
      gk=(gk-svr*xb(k))/xs(k)                                              1814
      u=gk+xv(k)*b(k)                                                      1814
      au=abs(u)-vp(k)*al1                                                  1815
      if(au .gt. 0.0)goto 14391                                            1815
      b(k)=0.0                                                             1815
      goto 14401                                                           1816
14391 continue                                                             1816
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1816
14401 continue                                                             1817
14381 continue                                                             1817
      d=b(k)-bk                                                            1817
      if(abs(d).le.0.0)goto 14371                                          1817
      dlx=max(dlx,xv(k)*d**2)                                              1818
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1819
      o=o+d*(xb(k)/xs(k))                                                  1820
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1821
14371 continue                                                             1822
14372 continue                                                             1822
      d=svr/xm(0)                                                          1823
      if(d .eq. 0.0)goto 14421                                             1823
      b(0)=b(0)+d                                                          1823
      dlx=max(dlx,xm(0)*d**2)                                              1823
      r=r-d*v                                                              1823
14421 continue                                                             1824
      svr=svr-d*xm(0)                                                      1824
      if(dlx.lt.shr)goto 14362                                             1825
      if(nlp .le. maxit)goto 14441                                         1825
      jerr=-ilm                                                            1825
      return                                                               1825
14441 continue                                                             1826
      goto 14361                                                           1827
14362 continue                                                             1827
      goto 14251                                                           1828
14252 continue                                                             1828
      if(nin.gt.nx)goto 14212                                              1829
      sc=b(0)                                                              1829
      b0=0.0                                                               1830
14450 do 14451 j=1,nin                                                     1830
      l=m(j)                                                               1830
      jb=ix(l)                                                             1830
      je=ix(l+1)-1                                                         1831
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      1832
      b0=b0-b(l)*xb(l)/xs(l)                                               1833
14451 continue                                                             1834
14452 continue                                                             1834
      sc=sc+b0                                                             1835
14460 do 14461 i=1,no                                                      1835
      fi=sc(i)+g(i)                                                        1836
      if(fi .ge. fmin)goto 14481                                           1836
      q(i)=0.0                                                             1836
      goto 14471                                                           1836
14481 if(fi .le. fmax)goto 14491                                           1836
      q(i)=1.0                                                             1836
      goto 14501                                                           1837
14491 continue                                                             1837
      q(i)=1.0/(1.0+exp(-fi))                                              1837
14501 continue                                                             1838
14471 continue                                                             1838
14461 continue                                                             1839
14462 continue                                                             1839
      v=w*q*(1.0-q)                                                        1839
      xm(0)=sum(v)                                                         1839
      if(xm(0).lt.vmin)goto 14212                                          1840
      r=w*(y-q)                                                            1840
      svr=sum(r)                                                           1840
      o=0.0                                                                1841
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 14521                         1841
      kx=0                                                                 1842
14530 do 14531 j=1,nin                                                     1842
      k=m(j)                                                               1843
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 14531                           1843
      kx=1                                                                 1843
      goto 14532                                                           1844
14531 continue                                                             1845
14532 continue                                                             1845
      if(kx .ne. 0)goto 14551                                              1846
14560 do 14561 j=1,ni                                                      1846
      if(ixx(j).eq.1)goto 14561                                            1846
      if(ju(j).eq.0)goto 14561                                             1847
      jb=ix(j)                                                             1847
      je=ix(j+1)-1                                                         1847
      jn=ix(j+1)-ix(j)                                                     1848
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1849
      gj=dot_product(sc(1:jn),x(jb:je))                                    1850
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      1851
      if(ga(j) .le. al1*vp(j))goto 14581                                   1851
      ixx(j)=1                                                             1851
      kx=1                                                                 1851
14581 continue                                                             1852
14561 continue                                                             1853
14562 continue                                                             1853
      if(kx.eq.1) go to 10680                                              1854
      goto 14212                                                           1855
14551 continue                                                             1856
14521 continue                                                             1857
      goto 14211                                                           1858
14212 continue                                                             1858
      if(nin .le. nx)goto 14601                                            1858
      jerr=-10000-ilm                                                      1858
      goto 14132                                                           1858
14601 continue                                                             1859
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1859
      kin(ilm)=nin                                                         1860
      a0(ilm)=b(0)                                                         1860
      alm(ilm)=al                                                          1860
      lmu=ilm                                                              1861
      devi=dev2(no,w,y,q,pmin)                                             1862
      dev(ilm)=(dev1-devi)/dev0                                            1863
      if(ilm.lt.mnl)goto 14131                                             1863
      if(flmin.ge.1.0)goto 14131                                           1864
      me=0                                                                 1864
14610 do 14611 j=1,nin                                                     1864
      if(a(j,ilm).ne.0.0) me=me+1                                          1864
14611 continue                                                             1864
14612 continue                                                             1864
      if(me.gt.ne)goto 14132                                               1865
      if(dev(ilm).gt.devmax)goto 14132                                     1865
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 14132                             1866
      if(xm(0).lt.vmin)goto 14132                                          1867
14131 continue                                                             1868
14132 continue                                                             1868
      g=log(q/(1.0-q))                                                     1869
      deallocate(xm,b,bs,v,r,sc,xv,q,mm,ga,ixx)                            1870
      return                                                               1871
      end                                                                  1872
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam   1874 
     *,flmin,ulam,  shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,al
     *m,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1876 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni)    1877
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1878
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1879
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp                       
      real, dimension (:), allocatable :: sc,xm,v,r,ga                          
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is,iy                           
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(xm(0:ni),stat=ierr)                                         1890
      jerr=jerr+ierr                                                       1891
      allocate(r(1:no),stat=ierr)                                          1891
      jerr=jerr+ierr                                                       1892
      allocate(v(1:no),stat=ierr)                                          1892
      jerr=jerr+ierr                                                       1893
      allocate(mm(1:ni),stat=ierr)                                         1893
      jerr=jerr+ierr                                                       1894
      allocate(ga(1:ni),stat=ierr)                                         1894
      jerr=jerr+ierr                                                       1895
      allocate(iy(1:ni),stat=ierr)                                         1895
      jerr=jerr+ierr                                                       1896
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1896
      jerr=jerr+ierr                                                       1897
      allocate(sxp(1:no),stat=ierr)                                        1897
      jerr=jerr+ierr                                                       1898
      allocate(sc(1:no),stat=ierr)                                         1898
      jerr=jerr+ierr                                                       1899
      if(jerr.ne.0) return                                                 1900
      pmax=1.0-pmin                                                        1900
      emin=pmin/pmax                                                       1900
      emax=1.0/emin                                                        1901
      pfm=(1.0+pmin)*pmin                                                  1901
      pfx=(1.0-pmin)*pmax                                                  1901
      vmin=pfm*pmax                                                        1902
      bta=parm                                                             1902
      omb=1.0-bta                                                          1902
      dev1=0.0                                                             1902
      dev0=0.0                                                             1903
14620 do 14621 ic=1,nc                                                     1903
      q0=dot_product(w,y(:,ic))                                            1904
      if(q0 .gt. pmin)goto 14641                                           1904
      jerr =8000+ic                                                        1904
      return                                                               1904
14641 continue                                                             1905
      if(q0 .lt. 1.0-pmin)goto 14661                                       1905
      jerr =9000+ic                                                        1905
      return                                                               1905
14661 continue                                                             1906
      b(1:ni,ic)=0.0                                                       1906
      b(0,ic)=log(q0)                                                      1906
      dev1=dev1-q0*b(0,ic)                                                 1907
14670 do 14671 i=1,no                                                      1907
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1907
14671 continue                                                             1908
14672 continue                                                             1908
14621 continue                                                             1909
14622 continue                                                             1909
      dev0=dev0+dev1                                                       1909
      iy=0                                                                 1909
      al=0.0                                                               1910
      if(nonzero(no*nc,g) .ne. 0)goto 14691                                1911
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1911
      sxp=0.0                                                              1912
14700 do 14701 ic=1,nc                                                     1912
      q(:,ic)=exp(b(0,ic))                                                 1912
      sxp=sxp+q(:,ic)                                                      1912
14701 continue                                                             1913
14702 continue                                                             1913
      goto 14711                                                           1914
14691 continue                                                             1914
14720 do 14721 i=1,no                                                      1914
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1914
14721 continue                                                             1914
14722 continue                                                             1914
      sxp=0.0                                                              1915
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1915
      if(jerr.ne.0) return                                                 1916
14730 do 14731 ic=1,nc                                                     1916
      q(:,ic)=exp(b(0,ic)+g(:,ic))                                         1916
      sxp=sxp+q(:,ic)                                                      1916
14731 continue                                                             1917
14732 continue                                                             1917
14711 continue                                                             1918
14681 continue                                                             1918
      if(kopt .le. 0)goto 14751                                            1919
      if(isd .le. 0)goto 14771                                             1919
      xv=0.25                                                              1919
      goto 14781                                                           1920
14771 continue                                                             1921
14790 do 14791 j=1,ni                                                      1921
      if(ju(j).eq.0)goto 14791                                             1921
      jb=ix(j)                                                             1921
      je=ix(j+1)-1                                                         1922
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        1923
14791 continue                                                             1924
14792 continue                                                             1924
14781 continue                                                             1925
14761 continue                                                             1925
14751 continue                                                             1926
      if(flmin .ge. 1.0)goto 14811                                         1926
      eqs=max(eps,flmin)                                                   1926
      alf=eqs**(1.0/(nlam-1))                                              1926
14811 continue                                                             1927
      m=0                                                                  1927
      mm=0                                                                 1927
      nin=0                                                                1927
      nlp=0                                                                1927
      mnl=min(mnlam,nlam)                                                  1927
      bs=0.0                                                               1927
      svr=0.0                                                              1927
      o=0.0                                                                1928
      shr=shri*dev0                                                        1928
      ga=0.0                                                               1929
14820 do 14821 ic=1,nc                                                     1929
      v=q(:,ic)/sxp                                                        1929
      r=w*(y(:,ic)-v)                                                      1929
      v=w*v*(1.0-v)                                                        1930
14830 do 14831 j=1,ni                                                      1930
      if(ju(j).eq.0)goto 14831                                             1931
      jb=ix(j)                                                             1931
      je=ix(j+1)-1                                                         1931
      jn=ix(j+1)-ix(j)                                                     1932
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1933
      gj=dot_product(sc(1:jn),x(jb:je))                                    1934
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             1935
14831 continue                                                             1936
14832 continue                                                             1936
14821 continue                                                             1937
14822 continue                                                             1937
14840 do 14841 ilm=1,nlam                                                  1937
      al0=al                                                               1938
      if(flmin .lt. 1.0)goto 14861                                         1938
      al=ulam(ilm)                                                         1938
      goto 14851                                                           1939
14861 if(ilm .le. 2)goto 14871                                             1939
      al=al*alf                                                            1939
      goto 14851                                                           1940
14871 if(ilm .ne. 1)goto 14881                                             1940
      al=big                                                               1940
      goto 14891                                                           1941
14881 continue                                                             1941
      al0=0.0                                                              1942
14900 do 14901 j=1,ni                                                      1942
      if(ju(j).eq.0)goto 14901                                             1942
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1942
14901 continue                                                             1943
14902 continue                                                             1943
      al0=al0/max(bta,1.0e-3)                                              1943
      al=alf*al0                                                           1944
14891 continue                                                             1945
14851 continue                                                             1945
      al2=al*omb                                                           1945
      al1=al*bta                                                           1945
      tlam=bta*(2.0*al-al0)                                                1946
14910 do 14911 k=1,ni                                                      1946
      if(iy(k).eq.1)goto 14911                                             1946
      if(ju(k).eq.0)goto 14911                                             1947
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      1948
14911 continue                                                             1949
14912 continue                                                             1949
10680 continue                                                             1950
14920 continue                                                             1950
14921 continue                                                             1950
      ixx=0                                                                1950
      jxx=ixx                                                              1950
      ig=0                                                                 1951
14930 do 14931 ic=1,nc                                                     1951
      bs(0,ic)=b(0,ic)                                                     1952
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1953
      xm(0)=0.0                                                            1953
      svr=0.0                                                              1953
      o=0.0                                                                1954
14940 do 14941 i=1,no                                                      1954
      pic=q(i,ic)/sxp(i)                                                   1955
      if(pic .ge. pfm)goto 14961                                           1955
      pic=0.0                                                              1955
      v(i)=0.0                                                             1955
      goto 14951                                                           1956
14961 if(pic .le. pfx)goto 14971                                           1956
      pic=1.0                                                              1956
      v(i)=0.0                                                             1956
      goto 14981                                                           1957
14971 continue                                                             1957
      v(i)=w(i)*pic*(1.0-pic)                                              1957
      xm(0)=xm(0)+v(i)                                                     1957
14981 continue                                                             1958
14951 continue                                                             1958
      r(i)=w(i)*(y(i,ic)-pic)                                              1958
      svr=svr+r(i)                                                         1959
14941 continue                                                             1960
14942 continue                                                             1960
      if(xm(0).le.vmin)goto 14931                                          1960
      ig=1                                                                 1961
14990 do 14991 j=1,ni                                                      1961
      if(iy(j).eq.0)goto 14991                                             1962
      jb=ix(j)                                                             1962
      je=ix(j+1)-1                                                         1963
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             1964
      if(kopt .ne. 0)goto 15011                                            1965
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       1966
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          1967
15011 continue                                                             1968
14991 continue                                                             1969
14992 continue                                                             1969
15020 continue                                                             1969
15021 continue                                                             1969
      nlp=nlp+1                                                            1969
      dlx=0.0                                                              1970
15030 do 15031 k=1,ni                                                      1970
      if(iy(k).eq.0)goto 15031                                             1971
      jb=ix(k)                                                             1971
      je=ix(k+1)-1                                                         1971
      jn=ix(k+1)-ix(k)                                                     1971
      bk=b(k,ic)                                                           1972
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1973
      gk=dot_product(sc(1:jn),x(jb:je))                                    1974
      gk=(gk-svr*xb(k))/xs(k)                                              1975
      u=gk+xv(k,ic)*b(k,ic)                                                1975
      au=abs(u)-vp(k)*al1                                                  1976
      if(au .gt. 0.0)goto 15051                                            1976
      b(k,ic)=0.0                                                          1976
      goto 15061                                                           1977
15051 continue                                                             1977
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1977
15061 continue                                                             1978
15041 continue                                                             1978
      d=b(k,ic)-bk                                                         1978
      if(abs(d).le.0.0)goto 15031                                          1979
      dlx=max(dlx,xv(k,ic)*d**2)                                           1980
      if(mm(k) .ne. 0)goto 15081                                           1980
      nin=nin+1                                                            1981
      if(nin .le. nx)goto 15101                                            1981
      jxx=1                                                                1981
      goto 15032                                                           1981
15101 continue                                                             1982
      mm(k)=nin                                                            1982
      m(nin)=k                                                             1983
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             1984
15081 continue                                                             1985
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1986
      o=o+d*(xb(k)/xs(k))                                                  1987
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1988
15031 continue                                                             1989
15032 continue                                                             1989
      if(jxx.gt.0)goto 15022                                               1990
      d=svr/xm(0)                                                          1991
      if(d .eq. 0.0)goto 15121                                             1991
      b(0,ic)=b(0,ic)+d                                                    1991
      dlx=max(dlx,xm(0)*d**2)                                              1992
      r=r-d*v                                                              1992
      svr=svr-d*xm(0)                                                      1993
15121 continue                                                             1994
      if(dlx.lt.shr)goto 15022                                             1994
      if(nlp .le. maxit)goto 15141                                         1994
      jerr=-ilm                                                            1994
      return                                                               1994
15141 continue                                                             1995
15150 continue                                                             1995
15151 continue                                                             1995
      nlp=nlp+1                                                            1995
      dlx=0.0                                                              1996
15160 do 15161 l=1,nin                                                     1996
      k=m(l)                                                               1996
      jb=ix(k)                                                             1996
      je=ix(k+1)-1                                                         1997
      jn=ix(k+1)-ix(k)                                                     1997
      bk=b(k,ic)                                                           1998
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1999
      gk=dot_product(sc(1:jn),x(jb:je))                                    2000
      gk=(gk-svr*xb(k))/xs(k)                                              2001
      u=gk+xv(k,ic)*b(k,ic)                                                2001
      au=abs(u)-vp(k)*al1                                                  2002
      if(au .gt. 0.0)goto 15181                                            2002
      b(k,ic)=0.0                                                          2002
      goto 15191                                                           2003
15181 continue                                                             2003
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              2003
15191 continue                                                             2004
15171 continue                                                             2004
      d=b(k,ic)-bk                                                         2004
      if(abs(d).le.0.0)goto 15161                                          2005
      dlx=max(dlx,xv(k,ic)*d**2)                                           2006
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2007
      o=o+d*(xb(k)/xs(k))                                                  2008
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2009
15161 continue                                                             2010
15162 continue                                                             2010
      d=svr/xm(0)                                                          2011
      if(d .eq. 0.0)goto 15211                                             2011
      b(0,ic)=b(0,ic)+d                                                    2011
      dlx=max(dlx,xm(0)*d**2)                                              2012
      r=r-d*v                                                              2012
      svr=svr-d*xm(0)                                                      2013
15211 continue                                                             2014
      if(dlx.lt.shr)goto 15152                                             2014
      if(nlp .le. maxit)goto 15231                                         2014
      jerr=-ilm                                                            2014
      return                                                               2014
15231 continue                                                             2015
      goto 15151                                                           2016
15152 continue                                                             2016
      goto 15021                                                           2017
15022 continue                                                             2017
      if(jxx.gt.0)goto 14932                                               2018
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                         2019
      if(ixx .ne. 0)goto 15251                                             2020
15260 do 15261 j=1,nin                                                     2020
      k=m(j)                                                               2021
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 15281                2021
      ixx=1                                                                2021
      goto 15262                                                           2021
15281 continue                                                             2022
15261 continue                                                             2023
15262 continue                                                             2023
15251 continue                                                             2024
      sc=b(0,ic)+g(:,ic)                                                   2024
      b0=0.0                                                               2025
15290 do 15291 j=1,nin                                                     2025
      l=m(j)                                                               2025
      jb=ix(l)                                                             2025
      je=ix(l+1)-1                                                         2026
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   2027
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            2028
15291 continue                                                             2029
15292 continue                                                             2029
      sc=min(max(exmn,sc+b0),exmx)                                         2030
      sxp=sxp-q(:,ic)                                                      2031
      q(:,ic)=min(max(emin*sxp,exp(dble(sc))),emax*sxp)                    2032
      sxp=sxp+q(:,ic)                                                      2033
14931 continue                                                             2034
14932 continue                                                             2034
      s=-sum(b(0,:))/nc                                                    2034
      b(0,:)=b(0,:)+s                                                      2034
      sc=s                                                                 2034
      b0=0.0                                                               2035
15300 do 15301 j=1,nin                                                     2035
      l=m(j)                                                               2036
      if(vp(l) .gt. 0.0)goto 15321                                         2036
      s=sum(b(l,:))/nc                                                     2036
      goto 15331                                                           2037
15321 continue                                                             2037
      s=elc(parm,nc,b(l,:),is)                                             2037
15331 continue                                                             2038
15311 continue                                                             2038
      b(l,:)=b(l,:)-s                                                      2039
      jb=ix(l)                                                             2039
      je=ix(l+1)-1                                                         2040
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         2041
      b0=b0+s*xb(l)/xs(l)                                                  2042
15301 continue                                                             2043
15302 continue                                                             2043
      sc=sc+b0                                                             2043
      sc=exp(sc)                                                           2043
      sxp=sxp*sc                                                           2043
15340 do 15341 ic=1,nc                                                     2043
      q(:,ic)=q(:,ic)*sc                                                   2043
15341 continue                                                             2044
15342 continue                                                             2044
      if(jxx.gt.0)goto 14922                                               2044
      if(ig.eq.0)goto 14922                                                2045
      if(ixx .ne. 0)goto 15361                                             2045
      ga=0.0                                                               2046
15370 do 15371 ic=1,nc                                                     2046
      v=q(:,ic)/sxp                                                        2046
      r=w*(y(:,ic)-v)                                                      2046
      v=w*v*(1.0-v)                                                        2047
15380 do 15381 j=1,ni                                                      2047
      if(iy(j).eq.1)goto 15381                                             2047
      if(ju(j).eq.0)goto 15381                                             2048
      jb=ix(j)                                                             2048
      je=ix(j+1)-1                                                         2048
      jn=ix(j+1)-ix(j)                                                     2049
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2050
      gj=dot_product(sc(1:jn),x(jb:je))                                    2051
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2052
15381 continue                                                             2053
15382 continue                                                             2053
15371 continue                                                             2054
15372 continue                                                             2054
15390 do 15391 k=1,ni                                                      2054
      if(iy(k).eq.1)goto 15391                                             2054
      if(ju(k).eq.0)goto 15391                                             2055
      if(ga(k) .le. al1*vp(k))goto 15411                                   2055
      iy(k)=1                                                              2055
      ixx=1                                                                2055
15411 continue                                                             2056
15391 continue                                                             2057
15392 continue                                                             2057
      if(ixx.eq.1) go to 10680                                             2058
      goto 14922                                                           2059
15361 continue                                                             2060
      goto 14921                                                           2061
14922 continue                                                             2061
      if(jxx .le. 0)goto 15431                                             2061
      jerr=-10000-ilm                                                      2061
      goto 14842                                                           2061
15431 continue                                                             2061
      devi=0.0                                                             2062
15440 do 15441 ic=1,nc                                                     2063
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          2063
      a0(ic,ilm)=b(0,ic)                                                   2064
15450 do 15451 i=1,no                                                      2064
      if(y(i,ic).le.0.0)goto 15451                                         2065
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           2066
15451 continue                                                             2067
15452 continue                                                             2067
15441 continue                                                             2068
15442 continue                                                             2068
      kin(ilm)=nin                                                         2068
      alm(ilm)=al                                                          2068
      lmu=ilm                                                              2069
      dev(ilm)=(dev1-devi)/dev0                                            2069
      if(ig.eq.0)goto 14842                                                2070
      if(ilm.lt.mnl)goto 14841                                             2070
      if(flmin.ge.1.0)goto 14841                                           2071
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 14842             2072
      if(dev(ilm).gt.devmax)goto 14842                                     2072
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 14842                             2073
14841 continue                                                             2074
14842 continue                                                             2074
      g=log(q)                                                             2074
15460 do 15461 i=1,no                                                      2074
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2074
15461 continue                                                             2075
15462 continue                                                             2075
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy)                      2076
      return                                                               2077
      end                                                                  2078
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  2079
      real a0(nc),ca(nx,nc),x(*),f(nc,n)                                   2079
      integer ia(*),ix(*),jx(*)                                            2080
15470 do 15471 ic=1,nc                                                     2080
      f(ic,:)=a0(ic)                                                       2080
15471 continue                                                             2081
15472 continue                                                             2081
15480 do 15481 j=1,nin                                                     2081
      k=ia(j)                                                              2081
      kb=ix(k)                                                             2081
      ke=ix(k+1)-1                                                         2082
15490 do 15491 ic=1,nc                                                     2082
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    2082
15491 continue                                                             2083
15492 continue                                                             2083
15481 continue                                                             2084
15482 continue                                                             2084
      return                                                               2085
      end                                                                  2086
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,ne,nx,nlam,flmin,ula   2088 
     *m,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam)              2089
      real ca(nx,nlam),dev(nlam),alm(nlam)                                 2090
      integer jd(*),ia(nx),nin(nlam)                                       2091
      real, dimension (:), allocatable :: xs,ww,vq                              
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 15511                                    2095
      jerr=10000                                                           2095
      return                                                               2095
15511 continue                                                             2096
      allocate(ww(1:no),stat=jerr)                                         2097
      allocate(ju(1:ni),stat=ierr)                                         2097
      jerr=jerr+ierr                                                       2098
      allocate(vq(1:ni),stat=ierr)                                         2098
      jerr=jerr+ierr                                                       2099
      if(isd .le. 0)goto 15531                                             2099
      allocate(xs(1:ni),stat=ierr)                                         2099
      jerr=jerr+ierr                                                       2099
15531 continue                                                             2100
      if(jerr.ne.0) return                                                 2101
      call chkvars(no,ni,x,ju)                                             2102
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2103
      if(maxval(ju) .gt. 0)goto 15551                                      2103
      jerr=7777                                                            2103
      return                                                               2103
15551 continue                                                             2104
      vq=max(0.0,vp)                                                       2104
      vq=vq*ni/sum(vq)                                                     2105
      ww=max(0.0,w)                                                        2105
      sw=sum(ww)                                                           2106
      if(sw .gt. 0.0)goto 15571                                            2106
      jerr=9999                                                            2106
      return                                                               2106
15571 continue                                                             2106
      ww=ww/sw                                                             2107
      call cstandard(no,ni,x,ww,ju,isd,xs)                                 2108
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr   2110 
     *,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 2110
      dev0=2.0*sw*dev0                                                     2111
      if(isd .le. 0)goto 15591                                             2111
15600 do 15601 k=1,lmu                                                     2111
      nk=nin(k)                                                            2111
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   2111
15601 continue                                                             2111
15602 continue                                                             2111
15591 continue                                                             2112
      deallocate(ww,ju,vq)                                                 2112
      if(isd.gt.0) deallocate(xs)                                          2113
      return                                                               2114
      end                                                                  2115
      subroutine cstandard (no,ni,x,w,ju,isd,xs)                           2116
      real x(no,ni),w(no),xs(ni)                                           2116
      integer ju(ni)                                                       2117
15610 do 15611 j=1,ni                                                      2117
      if(ju(j).eq.0)goto 15611                                             2118
      xm=dot_product(w,x(:,j))                                             2118
      x(:,j)=x(:,j)-xm                                                     2119
      if(isd .le. 0)goto 15631                                             2119
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 2119
      x(:,j)=x(:,j)/xs(j)                                                  2119
15631 continue                                                             2120
15611 continue                                                             2121
15612 continue                                                             2121
      return                                                               2122
      end                                                                  2123
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,ne,nx,nlam,flmin,ula   2125 
     *m,cthri,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=0.001, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99)   2126
      real x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam)              2127
      real ao(nx,nlam),dev(nlam),alm(nlam)                                 2128
      integer ju(ni),m(nx),kin(nlam)                                       2129
      real, dimension (:), allocatable :: w,dk,v,xs,wr,a,as,f,dq                
      real, dimension (:), allocatable :: e,uu,ga                               
      integer, dimension (:), allocatable :: jp,kp,mm,ixx                       
      allocate(e(1:no),stat=jerr)                                          2135
      allocate(uu(1:no),stat=ierr)                                         2135
      jerr=jerr+ierr                                                       2136
      allocate(f(1:no),stat=ierr)                                          2136
      jerr=jerr+ierr                                                       2137
      allocate(w(1:no),stat=ierr)                                          2137
      jerr=jerr+ierr                                                       2138
      allocate(v(1:ni),stat=ierr)                                          2138
      jerr=jerr+ierr                                                       2139
      allocate(a(1:ni),stat=ierr)                                          2139
      jerr=jerr+ierr                                                       2140
      allocate(as(1:ni),stat=ierr)                                         2140
      jerr=jerr+ierr                                                       2141
      allocate(xs(1:ni),stat=ierr)                                         2141
      jerr=jerr+ierr                                                       2142
      allocate(ga(1:ni),stat=ierr)                                         2142
      jerr=jerr+ierr                                                       2143
      allocate(ixx(1:ni),stat=ierr)                                        2143
      jerr=jerr+ierr                                                       2144
      allocate(jp(1:no),stat=ierr)                                         2144
      jerr=jerr+ierr                                                       2145
      allocate(kp(1:no),stat=ierr)                                         2145
      jerr=jerr+ierr                                                       2146
      allocate(dk(1:no),stat=ierr)                                         2146
      jerr=jerr+ierr                                                       2147
      allocate(wr(1:no),stat=ierr)                                         2147
      jerr=jerr+ierr                                                       2148
      allocate(dq(1:no),stat=ierr)                                         2148
      jerr=jerr+ierr                                                       2149
      allocate(mm(1:ni),stat=ierr)                                         2149
      jerr=jerr+ierr                                                       2150
      if(jerr.ne.0)go to 11790                                             2151
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2152
      if(jerr.ne.0) go to 11790                                            2152
      alpha=parm                                                           2153
      oma=1.0-alpha                                                        2153
      nlm=0                                                                2153
      ixx=0                                                                2153
      al=0.0                                                               2154
      dq=d*q                                                               2154
      call died(no,nk,dq,kp,jp,dk)                                         2155
      a=0.0                                                                2155
      f=0.0                                                                2155
      e=q                                                                  2155
      fmax=log(huge(f(1))*0.1)                                             2156
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 2157
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         2157
      dev0=rr                                                              2158
15640 do 15641 i=1,no                                                      2158
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 15661                   2158
      w(i)=0.0                                                             2158
      wr(i)=w(i)                                                           2158
15661 continue                                                             2158
15641 continue                                                             2159
15642 continue                                                             2159
      if(nonzero(no,g) .eq. 0)goto 15681                                   2159
      f=g-dot_product(q,g)                                                 2160
      e=q*exp(sign(min(abs(f),fmax),f))                                    2161
15681 continue                                                             2162
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2163
      if(jerr.ne.0) go to 11790                                            2164
      if(flmin .ge. 1.0)goto 15701                                         2164
      eqs=max(eps,flmin)                                                   2164
      alf=eqs**(1.0/(nlam-1))                                              2164
15701 continue                                                             2165
      m=0                                                                  2165
      mm=0                                                                 2165
      nlp=0                                                                2165
      nin=nlp                                                              2165
      mnl=min(mnlam,nlam)                                                  2165
      as=0.0                                                               2165
      cthr=cthri*dev0                                                      2166
15710 do 15711 j=1,ni                                                      2166
      if(ju(j).eq.0)goto 15711                                             2166
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2166
15711 continue                                                             2167
15712 continue                                                             2167
15720 do 15721 ilm=1,nlam                                                  2167
      al0=al                                                               2168
      if(flmin .lt. 1.0)goto 15741                                         2168
      al=ulam(ilm)                                                         2168
      goto 15731                                                           2169
15741 if(ilm .le. 2)goto 15751                                             2169
      al=al*alf                                                            2169
      goto 15731                                                           2170
15751 if(ilm .ne. 1)goto 15761                                             2170
      al=big                                                               2170
      goto 15771                                                           2171
15761 continue                                                             2171
      al0=0.0                                                              2172
15780 do 15781 j=1,ni                                                      2172
      if(ju(j).eq.0)goto 15781                                             2172
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2172
15781 continue                                                             2173
15782 continue                                                             2173
      al0=al0/max(parm,1.0e-3)                                             2173
      al=alf*al0                                                           2174
15771 continue                                                             2175
15731 continue                                                             2175
      sa=alpha*al                                                          2175
      omal=oma*al                                                          2175
      tlam=alpha*(2.0*al-al0)                                              2176
15790 do 15791 k=1,ni                                                      2176
      if(ixx(k).eq.1)goto 15791                                            2176
      if(ju(k).eq.0)goto 15791                                             2177
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2178
15791 continue                                                             2179
15792 continue                                                             2179
10680 continue                                                             2180
15800 continue                                                             2180
15801 continue                                                             2180
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2181
      call vars(no,ni,x,w,ixx,v)                                           2182
15810 continue                                                             2182
15811 continue                                                             2182
      nlp=nlp+1                                                            2182
      dli=0.0                                                              2183
15820 do 15821 j=1,ni                                                      2183
      if(ixx(j).eq.0)goto 15821                                            2184
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2185
      if(abs(u) .gt. vp(j)*sa)goto 15841                                   2185
      at=0.0                                                               2185
      goto 15851                                                           2186
15841 continue                                                             2186
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2186
15851 continue                                                             2187
15831 continue                                                             2187
      if(at .eq. a(j))goto 15871                                           2187
      del=at-a(j)                                                          2187
      a(j)=at                                                              2187
      dli=max(dli,v(j)*del**2)                                             2188
      wr=wr-del*w*x(:,j)                                                   2188
      f=f+del*x(:,j)                                                       2189
      if(mm(j) .ne. 0)goto 15891                                           2189
      nin=nin+1                                                            2189
      if(nin.gt.nx)goto 15822                                              2190
      mm(j)=nin                                                            2190
      m(nin)=j                                                             2191
15891 continue                                                             2192
15871 continue                                                             2193
15821 continue                                                             2194
15822 continue                                                             2194
      if(nin.gt.nx)goto 15812                                              2194
      if(dli.lt.cthr)goto 15812                                            2195
      if(nlp .le. maxit)goto 15911                                         2195
      jerr=-ilm                                                            2195
      return                                                               2195
15911 continue                                                             2196
15920 continue                                                             2196
15921 continue                                                             2196
      nlp=nlp+1                                                            2196
      dli=0.0                                                              2197
15930 do 15931 l=1,nin                                                     2197
      j=m(l)                                                               2198
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2199
      if(abs(u) .gt. vp(j)*sa)goto 15951                                   2199
      at=0.0                                                               2199
      goto 15961                                                           2200
15951 continue                                                             2200
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2200
15961 continue                                                             2201
15941 continue                                                             2201
      if(at .eq. a(j))goto 15981                                           2201
      del=at-a(j)                                                          2201
      a(j)=at                                                              2201
      dli=max(dli,v(j)*del**2)                                             2202
      wr=wr-del*w*x(:,j)                                                   2202
      f=f+del*x(:,j)                                                       2203
15981 continue                                                             2204
15931 continue                                                             2205
15932 continue                                                             2205
      if(dli.lt.cthr)goto 15922                                            2205
      if(nlp .le. maxit)goto 16001                                         2205
      jerr=-ilm                                                            2205
      return                                                               2205
16001 continue                                                             2206
      goto 15921                                                           2207
15922 continue                                                             2207
      goto 15811                                                           2208
15812 continue                                                             2208
      if(nin.gt.nx)goto 15802                                              2209
      e=q*exp(sign(min(abs(f),fmax),f))                                    2210
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2211
      if(jerr .eq. 0)goto 16021                                            2211
      jerr=jerr-ilm                                                        2211
      go to 11790                                                          2211
16021 continue                                                             2212
      ix=0                                                                 2213
16030 do 16031 j=1,nin                                                     2213
      k=m(j)                                                               2214
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 16031                           2214
      ix=1                                                                 2214
      goto 16032                                                           2214
16031 continue                                                             2215
16032 continue                                                             2215
      if(ix .ne. 0)goto 16051                                              2216
16060 do 16061 k=1,ni                                                      2216
      if(ixx(k).eq.1)goto 16061                                            2216
      if(ju(k).eq.0)goto 16061                                             2217
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2218
      if(ga(k) .le. sa*vp(k))goto 16081                                    2218
      ixx(k)=1                                                             2218
      ix=1                                                                 2218
16081 continue                                                             2219
16061 continue                                                             2220
16062 continue                                                             2220
      if(ix.eq.1) go to 10680                                              2221
      goto 15802                                                           2222
16051 continue                                                             2223
      goto 15801                                                           2224
15802 continue                                                             2224
      if(nin .le. nx)goto 16101                                            2224
      jerr=-10000-ilm                                                      2224
      goto 15722                                                           2224
16101 continue                                                             2225
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               2225
      kin(ilm)=nin                                                         2226
      alm(ilm)=al                                                          2226
      lmu=ilm                                                              2227
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   2228
      if(ilm.lt.mnl)goto 15721                                             2228
      if(flmin.ge.1.0)goto 15721                                           2229
      me=0                                                                 2229
16110 do 16111 j=1,nin                                                     2229
      if(ao(j,ilm).ne.0.0) me=me+1                                         2229
16111 continue                                                             2229
16112 continue                                                             2229
      if(me.gt.ne)goto 15722                                               2230
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 15722              2231
      if(dev(ilm).gt.devmax)goto 15722                                     2232
15721 continue                                                             2233
15722 continue                                                             2233
      g=f                                                                  2234
11790 continue                                                             2234
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)              2235
      return                                                               2236
      end                                                                  2237
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 2238
      real ca(nin),x(n,*),f(n)                                             2238
      integer ia(nin)                                                      2239
      f=0.0                                                                2239
      if(nin.le.0) return                                                  2240
16120 do 16121 i=1,n                                                       2240
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      2240
16121 continue                                                             2241
16122 continue                                                             2241
      return                                                               2242
      end                                                                  2243
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         2244
      real y(no),d(no),q(no)                                               2244
      integer jp(no),kp(*)                                                 2245
16130 do 16131 j=1,no                                                      2245
      jp(j)=j                                                              2245
16131 continue                                                             2245
16132 continue                                                             2245
      call psort7(y,jp,1,no)                                               2246
      nj=0                                                                 2246
16140 do 16141 j=1,no                                                      2246
      if(q(jp(j)).le.0.0)goto 16141                                        2246
      nj=nj+1                                                              2246
      jp(nj)=jp(j)                                                         2246
16141 continue                                                             2247
16142 continue                                                             2247
      if(nj .ne. 0)goto 16161                                              2247
      jerr=20000                                                           2247
      return                                                               2247
16161 continue                                                             2248
      j=1                                                                  2248
16170 continue                                                             2248
16171 if(d(jp(j)).gt.0.0)goto 16172                                        2248
      j=j+1                                                                2248
      if(j.gt.nj)goto 16172                                                2248
      goto 16171                                                           2249
16172 continue                                                             2249
      if(j .lt. nj-1)goto 16191                                            2249
      jerr=30000                                                           2249
      return                                                               2249
16191 continue                                                             2250
      j0=j-1                                                               2250
      nj=nj-j0                                                             2250
16200 do 16201 j=1,nj                                                      2250
      jp(j)=jp(j+j0)                                                       2250
16201 continue                                                             2251
16202 continue                                                             2251
      jerr=0                                                               2251
      nk=0                                                                 2251
      t0=y(jp(1))                                                          2251
      yk=t0                                                                2251
      j=2                                                                  2252
16210 continue                                                             2252
16211 continue                                                             2252
16220 continue                                                             2253
16221 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 16222                     2253
      j=j+1                                                                2253
      if(j.gt.nj)goto 16222                                                2253
      goto 16221                                                           2254
16222 continue                                                             2254
      nk=nk+1                                                              2254
      kp(nk)=j-1                                                           2254
      if(j.gt.nj)goto 16212                                                2255
      if(j .ne. nj)goto 16241                                              2255
      nk=nk+1                                                              2255
      kp(nk)=nj                                                            2255
      goto 16212                                                           2255
16241 continue                                                             2256
      yk=y(jp(j))                                                          2256
      j=j+1                                                                2257
      goto 16211                                                           2258
16212 continue                                                             2258
      return                                                               2259
      end                                                                  2260
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     2261
      real d(no),dk(nk),wr(no),w(no)                                       2262
      real e(no),u(no),b,c                                                 2262
      integer kp(nk),jp(no)                                                2263
      call usk(no,nk,kp,jp,e,u)                                            2264
      b=dk(1)/u(1)                                                         2264
      c=dk(1)/u(1)**2                                                      2264
      jerr=0                                                               2265
16250 do 16251 j=1,kp(1)                                                   2265
      i=jp(j)                                                              2266
      w(i)=e(i)*(b-e(i)*c)                                                 2266
      if(w(i) .gt. 0.0)goto 16271                                          2266
      jerr=-3                                                              2266
      return                                                               2266
16271 continue                                                             2267
      wr(i)=d(i)-e(i)*b                                                    2268
16251 continue                                                             2269
16252 continue                                                             2269
16280 do 16281 k=2,nk                                                      2269
      j1=kp(k-1)+1                                                         2269
      j2=kp(k)                                                             2270
      b=b+dk(k)/u(k)                                                       2270
      c=c+dk(k)/u(k)**2                                                    2271
16290 do 16291 j=j1,j2                                                     2271
      i=jp(j)                                                              2272
      w(i)=e(i)*(b-e(i)*c)                                                 2272
      if(w(i) .gt. 0.0)goto 16311                                          2272
      jerr=-30000                                                          2272
      return                                                               2272
16311 continue                                                             2273
      wr(i)=d(i)-e(i)*b                                                    2274
16291 continue                                                             2275
16292 continue                                                             2275
16281 continue                                                             2276
16282 continue                                                             2276
      return                                                               2277
      end                                                                  2278
      subroutine vars(no,ni,x,w,ixx,v)                                     2279
      real x(no,ni),w(no),v(ni)                                            2279
      integer ixx(ni)                                                      2280
16320 do 16321 j=1,ni                                                      2280
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                        2280
16321 continue                                                             2281
16322 continue                                                             2281
      return                                                               2282
      end                                                                  2283
      subroutine died(no,nk,d,kp,jp,dk)                                    2284
      real d(no),dk(nk)                                                    2284
      integer kp(nk),jp(no)                                                2285
      dk(1)=sum(d(jp(1:kp(1))))                                            2286
16330 do 16331 k=2,nk                                                      2286
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  2286
16331 continue                                                             2287
16332 continue                                                             2287
      return                                                               2288
      end                                                                  2289
      subroutine usk(no,nk,kp,jp,e,u)                                      2290
      real e(no),u(nk),h                                                   2290
      integer kp(nk),jp(no)                                                2291
      h=0.0                                                                2292
16340 do 16341 k=nk,1,-1                                                   2292
      j2=kp(k)                                                             2293
      j1=1                                                                 2293
      if(k.gt.1) j1=kp(k-1)+1                                              2294
16350 do 16351 j=j2,j1,-1                                                  2294
      h=h+e(jp(j))                                                         2294
16351 continue                                                             2295
16352 continue                                                             2295
      u(k)=h                                                               2296
16341 continue                                                             2297
16342 continue                                                             2297
      return                                                               2298
      end                                                                  2299
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             2300
      real d(no),dk(nk),f(no)                                              2301
      integer kp(nk),jp(no)                                                2301
      real e(no),u(nk),s                                                   2302
      call usk(no,nk,kp,jp,e,u)                                            2302
      u=log(u)                                                             2303
      risk=dot_product(d,f)-dot_product(dk,u)                              2304
      return                                                               2305
      end                                                                  2306
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 2307
      real x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(nlam)          2308
      real, dimension (:), allocatable :: dk,f,xm,dq,q                          
      real, dimension (:), allocatable :: e,uu                                  
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          2314
      allocate(q(1:no),stat=ierr)                                          2314
      jerr=jerr+ierr                                                       2315
      allocate(uu(1:no),stat=ierr)                                         2315
      jerr=jerr+ierr                                                       2316
      allocate(f(1:no),stat=ierr)                                          2316
      jerr=jerr+ierr                                                       2317
      allocate(dk(1:no),stat=ierr)                                         2317
      jerr=jerr+ierr                                                       2318
      allocate(jp(1:no),stat=ierr)                                         2318
      jerr=jerr+ierr                                                       2319
      allocate(kp(1:no),stat=ierr)                                         2319
      jerr=jerr+ierr                                                       2320
      allocate(dq(1:no),stat=ierr)                                         2320
      jerr=jerr+ierr                                                       2321
      allocate(xm(1:ni),stat=ierr)                                         2321
      jerr=jerr+ierr                                                       2322
      if(jerr.ne.0) go to 11790                                            2323
      q=max(0.0,w)                                                         2323
      sw=sum(q)                                                            2324
      if(sw .gt. 0.0)goto 16371                                            2324
      jerr=9999                                                            2324
      go to 11790                                                          2324
16371 continue                                                             2325
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2326
      if(jerr.ne.0) go to 11790                                            2326
      fmax=log(huge(e(1))*0.1)                                             2327
      dq=d*q                                                               2327
      call died(no,nk,dq,kp,jp,dk)                                         2327
      gm=dot_product(q,g)/sw                                               2328
16380 do 16381 j=1,ni                                                      2328
      xm(j)=dot_product(q,x(:,j))/sw                                       2328
16381 continue                                                             2329
16382 continue                                                             2329
16390 do 16391 lam=1,nlam                                                  2330
16400 do 16401 i=1,no                                                      2330
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       2331
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        2332
16401 continue                                                             2333
16402 continue                                                             2333
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          2334
16391 continue                                                             2335
16392 continue                                                             2335
11790 continue                                                             2335
      deallocate(e,uu,dk,f,jp,kp,dq)                                       2336
      return                                                               2337
      end                                                                  2338
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,ne,nx,nlam,flmin,ulam   2340 
     *,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    2341
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2342
      integer jd(*),ia(nx),nin(nlam)                                       2343
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16421                                    2347
      jerr=10000                                                           2347
      return                                                               2347
16421 continue                                                             2348
      if(minval(y) .ge. 0.0)goto 16441                                     2348
      jerr=8888                                                            2348
      return                                                               2348
16441 continue                                                             2349
      allocate(ww(1:no),stat=jerr)                                         2350
      allocate(ju(1:ni),stat=ierr)                                         2350
      jerr=jerr+ierr                                                       2351
      allocate(vq(1:ni),stat=ierr)                                         2351
      jerr=jerr+ierr                                                       2352
      allocate(xm(1:ni),stat=ierr)                                         2352
      jerr=jerr+ierr                                                       2353
      if(isd .le. 0)goto 16461                                             2353
      allocate(xs(1:ni),stat=ierr)                                         2353
      jerr=jerr+ierr                                                       2353
16461 continue                                                             2354
      if(jerr.ne.0) return                                                 2355
      call chkvars(no,ni,x,ju)                                             2356
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2357
      if(maxval(ju) .gt. 0)goto 16481                                      2357
      jerr=7777                                                            2357
      go to 11790                                                          2357
16481 continue                                                             2358
      vq=max(0.0,vp)                                                       2358
      vq=vq*ni/sum(vq)                                                     2359
      ww=max(0.0,w)                                                        2359
      sw=sum(ww)                                                           2359
      if(sw .gt. 0.0)goto 16501                                            2359
      jerr=9999                                                            2359
      go to 11790                                                          2359
16501 continue                                                             2360
      ww=ww/sw                                                             2361
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             2362
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr,   2364 
     *  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11790                                            2364
      dev0=2.0*sw*dev0                                                     2365
16510 do 16511 k=1,lmu                                                     2365
      nk=nin(k)                                                            2366
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2367
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2368
16511 continue                                                             2369
16512 continue                                                             2369
11790 continue                                                             2369
      deallocate(ww,ju,vq,xm)                                              2369
      if(isd.gt.0) deallocate(xs)                                          2370
      return                                                               2371
      end                                                                  2372
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,ne,nx,nlam,flmin,ulam   2374 
     *,shri,  isd,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2375 
     *9)
      real x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)                    2376
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2377
      integer ju(ni),m(nx),kin(nlam)                                       2378
      real, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga                    
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(a(1:ni),stat=jerr)                                          2383
      allocate(as(1:ni),stat=ierr)                                         2383
      jerr=jerr+ierr                                                       2384
      allocate(t(1:no),stat=ierr)                                          2384
      jerr=jerr+ierr                                                       2385
      allocate(mm(1:ni),stat=ierr)                                         2385
      jerr=jerr+ierr                                                       2386
      allocate(ga(1:ni),stat=ierr)                                         2386
      jerr=jerr+ierr                                                       2387
      allocate(ixx(1:ni),stat=ierr)                                        2387
      jerr=jerr+ierr                                                       2388
      allocate(wr(1:no),stat=ierr)                                         2388
      jerr=jerr+ierr                                                       2389
      allocate(v(1:ni),stat=ierr)                                          2389
      jerr=jerr+ierr                                                       2390
      allocate(w(1:no),stat=ierr)                                          2390
      jerr=jerr+ierr                                                       2391
      allocate(f(1:no),stat=ierr)                                          2391
      jerr=jerr+ierr                                                       2392
      if(jerr.ne.0) return                                                 2393
      bta=parm                                                             2393
      omb=1.0-bta                                                          2394
      t=q*y                                                                2394
      yb=sum(t)                                                            2394
      fmax=log(huge(bta)*0.1)                                              2395
      if(nonzero(no,g) .ne. 0)goto 16531                                   2395
      w=q*yb                                                               2395
      az=log(yb)                                                           2395
      f=az                                                                 2395
      goto 16541                                                           2396
16531 continue                                                             2396
      w=q*exp(sign(min(abs(g),fmax),g))                                    2396
      v0=sum(w)                                                            2396
      eaz=yb/v0                                                            2397
      w=eaz*w                                                              2397
      az=log(eaz)                                                          2397
      f=az+g                                                               2398
16541 continue                                                             2399
16521 continue                                                             2399
      a=0.0                                                                2399
      as=0.0                                                               2399
      wr=t-w                                                               2399
      v0=yb                                                                2399
      dv0=yb*(log(yb)-1.0)                                                 2399
      dvr=-yb                                                              2400
16550 do 16551 i=1,no                                                      2400
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               2400
16551 continue                                                             2400
16552 continue                                                             2400
      dvr=dvr-dv0                                                          2400
      dev0=dvr                                                             2401
      if(flmin .ge. 1.0)goto 16571                                         2401
      eqs=max(eps,flmin)                                                   2401
      alf=eqs**(1.0/(nlam-1))                                              2401
16571 continue                                                             2402
      m=0                                                                  2402
      mm=0                                                                 2402
      nlp=0                                                                2402
      nin=nlp                                                              2402
      mnl=min(mnlam,nlam)                                                  2402
      shr=shri*dev0                                                        2402
      ixx=0                                                                2402
      al=0.0                                                               2403
16580 do 16581 j=1,ni                                                      2403
      if(ju(j).eq.0)goto 16581                                             2403
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2403
16581 continue                                                             2404
16582 continue                                                             2404
16590 do 16591 ilm=1,nlam                                                  2404
      al0=al                                                               2405
      if(flmin .lt. 1.0)goto 16611                                         2405
      al=ulam(ilm)                                                         2405
      goto 16601                                                           2406
16611 if(ilm .le. 2)goto 16621                                             2406
      al=al*alf                                                            2406
      goto 16601                                                           2407
16621 if(ilm .ne. 1)goto 16631                                             2407
      al=big                                                               2407
      goto 16641                                                           2408
16631 continue                                                             2408
      al0=0.0                                                              2409
16650 do 16651 j=1,ni                                                      2409
      if(ju(j).eq.0)goto 16651                                             2409
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2409
16651 continue                                                             2410
16652 continue                                                             2410
      al0=al0/max(bta,1.0e-3)                                              2410
      al=alf*al0                                                           2411
16641 continue                                                             2412
16601 continue                                                             2412
      al2=al*omb                                                           2412
      al1=al*bta                                                           2412
      tlam=bta*(2.0*al-al0)                                                2413
16660 do 16661 k=1,ni                                                      2413
      if(ixx(k).eq.1)goto 16661                                            2413
      if(ju(k).eq.0)goto 16661                                             2414
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2415
16661 continue                                                             2416
16662 continue                                                             2416
10680 continue                                                             2417
16670 continue                                                             2417
16671 continue                                                             2417
      az0=az                                                               2418
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2419
16680 do 16681 j=1,ni                                                      2419
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                        2419
16681 continue                                                             2420
16682 continue                                                             2420
16690 continue                                                             2420
16691 continue                                                             2420
      nlp=nlp+1                                                            2420
      dlx=0.0                                                              2421
16700 do 16701 k=1,ni                                                      2421
      if(ixx(k).eq.0)goto 16701                                            2421
      ak=a(k)                                                              2422
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2422
      au=abs(u)-vp(k)*al1                                                  2423
      if(au .gt. 0.0)goto 16721                                            2423
      a(k)=0.0                                                             2423
      goto 16731                                                           2424
16721 continue                                                             2424
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2424
16731 continue                                                             2425
16711 continue                                                             2425
      if(a(k).eq.ak)goto 16701                                             2425
      d=a(k)-ak                                                            2425
      dlx=max(dlx,v(k)*d**2)                                               2426
      wr=wr-d*w*x(:,k)                                                     2426
      f=f+d*x(:,k)                                                         2427
      if(mm(k) .ne. 0)goto 16751                                           2427
      nin=nin+1                                                            2427
      if(nin.gt.nx)goto 16702                                              2428
      mm(k)=nin                                                            2428
      m(nin)=k                                                             2429
16751 continue                                                             2430
16701 continue                                                             2431
16702 continue                                                             2431
      if(nin.gt.nx)goto 16692                                              2431
      d=sum(wr)/v0                                                         2432
      az=az+d                                                              2432
      dlx=max(dlx,v0*d**2)                                                 2432
      wr=wr-d*w                                                            2432
      f=f+d                                                                2433
      if(dlx.lt.shr)goto 16692                                             2433
      if(nlp .le. maxit)goto 16771                                         2433
      jerr=-ilm                                                            2433
      return                                                               2433
16771 continue                                                             2434
16780 continue                                                             2434
16781 continue                                                             2434
      nlp=nlp+1                                                            2434
      dlx=0.0                                                              2435
16790 do 16791 l=1,nin                                                     2435
      k=m(l)                                                               2435
      ak=a(k)                                                              2436
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2436
      au=abs(u)-vp(k)*al1                                                  2437
      if(au .gt. 0.0)goto 16811                                            2437
      a(k)=0.0                                                             2437
      goto 16821                                                           2438
16811 continue                                                             2438
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2438
16821 continue                                                             2439
16801 continue                                                             2439
      if(a(k).eq.ak)goto 16791                                             2439
      d=a(k)-ak                                                            2439
      dlx=max(dlx,v(k)*d**2)                                               2440
      wr=wr-d*w*x(:,k)                                                     2440
      f=f+d*x(:,k)                                                         2442
16791 continue                                                             2442
16792 continue                                                             2442
      d=sum(wr)/v0                                                         2442
      az=az+d                                                              2442
      dlx=max(dlx,v0*d**2)                                                 2442
      wr=wr-d*w                                                            2442
      f=f+d                                                                2443
      if(dlx.lt.shr)goto 16782                                             2443
      if(nlp .le. maxit)goto 16841                                         2443
      jerr=-ilm                                                            2443
      return                                                               2443
16841 continue                                                             2444
      goto 16781                                                           2445
16782 continue                                                             2445
      goto 16691                                                           2446
16692 continue                                                             2446
      if(nin.gt.nx)goto 16672                                              2447
      w=q*exp(sign(min(abs(f),fmax),f))                                    2447
      v0=sum(w)                                                            2447
      wr=t-w                                                               2448
      if(v0*(az-az0)**2 .ge. shr)goto 16861                                2448
      ix=0                                                                 2449
16870 do 16871 j=1,nin                                                     2449
      k=m(j)                                                               2450
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 16871                            2450
      ix=1                                                                 2450
      goto 16872                                                           2451
16871 continue                                                             2452
16872 continue                                                             2452
      if(ix .ne. 0)goto 16891                                              2453
16900 do 16901 k=1,ni                                                      2453
      if(ixx(k).eq.1)goto 16901                                            2453
      if(ju(k).eq.0)goto 16901                                             2454
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2455
      if(ga(k) .le. al1*vp(k))goto 16921                                   2455
      ixx(k)=1                                                             2455
      ix=1                                                                 2455
16921 continue                                                             2456
16901 continue                                                             2457
16902 continue                                                             2457
      if(ix.eq.1) go to 10680                                              2458
      goto 16672                                                           2459
16891 continue                                                             2460
16861 continue                                                             2461
      goto 16671                                                           2462
16672 continue                                                             2462
      if(nin .le. nx)goto 16941                                            2462
      jerr=-10000-ilm                                                      2462
      goto 16592                                                           2462
16941 continue                                                             2463
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2463
      kin(ilm)=nin                                                         2464
      a0(ilm)=az                                                           2464
      alm(ilm)=al                                                          2464
      lmu=ilm                                                              2465
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               2466
      if(ilm.lt.mnl)goto 16591                                             2466
      if(flmin.ge.1.0)goto 16591                                           2467
      me=0                                                                 2467
16950 do 16951 j=1,nin                                                     2467
      if(ca(j,ilm).ne.0.0) me=me+1                                         2467
16951 continue                                                             2467
16952 continue                                                             2467
      if(me.gt.ne)goto 16592                                               2468
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 16592              2469
      if(dev(ilm).gt.devmax)goto 16592                                     2470
16591 continue                                                             2471
16592 continue                                                             2471
      g=f                                                                  2472
11790 continue                                                             2472
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                                2473
      return                                                               2474
      end                                                                  2475
      function nonzero(n,v)                                                2476
      real v(n)                                                            2477
      nonzero=0                                                            2477
16960 do 16961 i=1,n                                                       2477
      if(v(i) .eq. 0.0)goto 16981                                          2477
      nonzero=1                                                            2477
      return                                                               2477
16981 continue                                                             2477
16961 continue                                                             2478
16962 continue                                                             2478
      return                                                               2479
      end                                                                  2480
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               2481
      real a(nx,lmu),b(ni,lmu)                                             2481
      integer ia(nx),nin(lmu)                                              2482
16990 do 16991 lam=1,lmu                                                   2482
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        2482
16991 continue                                                             2483
16992 continue                                                             2483
      return                                                               2484
      end                                                                  2485
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           2486
      real a(nx,nc,lmu),b(ni,nc,lmu)                                       2486
      integer ia(nx),nin(lmu)                                              2487
17000 do 17001 lam=1,lmu                                                   2487
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             2487
17001 continue                                                             2488
17002 continue                                                             2488
      return                                                               2489
      end                                                                  2490
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               2491
      real x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)       2492
      real, dimension (:), allocatable :: w                                     
      if(minval(y) .ge. 0.0)goto 17021                                     2495
      jerr=8888                                                            2495
      return                                                               2495
17021 continue                                                             2496
      allocate(w(1:no),stat=jerr)                                          2496
      if(jerr.ne.0) return                                                 2497
      w=max(0.0,q)                                                         2497
      sw=sum(w)                                                            2497
      if(sw .gt. 0.0)goto 17041                                            2497
      jerr=9999                                                            2497
      go to 11790                                                          2497
17041 continue                                                             2498
      yb=dot_product(w,y)/sw                                               2498
      fmax=log(huge(y(1))*0.1)                                             2499
17050 do 17051 lam=1,nlam                                                  2499
      s=0.0                                                                2500
17060 do 17061 i=1,no                                                      2500
      if(w(i).le.0.0)goto 17061                                            2501
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          2502
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      2503
17061 continue                                                             2504
17062 continue                                                             2504
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2505
17051 continue                                                             2506
17052 continue                                                             2506
11790 continue                                                             2506
      deallocate(w)                                                        2507
      return                                                               2508
      end                                                                  2509
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,ne,nx,nlam,fl   2511 
     *min,ulam,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        2512
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2513
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           2514
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 17081                                    2518
      jerr=10000                                                           2518
      return                                                               2518
17081 continue                                                             2519
      if(minval(y) .ge. 0.0)goto 17101                                     2519
      jerr=8888                                                            2519
      return                                                               2519
17101 continue                                                             2520
      allocate(ww(1:no),stat=jerr)                                         2521
      allocate(ju(1:ni),stat=ierr)                                         2521
      jerr=jerr+ierr                                                       2522
      allocate(vq(1:ni),stat=ierr)                                         2522
      jerr=jerr+ierr                                                       2523
      allocate(xm(1:ni),stat=ierr)                                         2523
      jerr=jerr+ierr                                                       2524
      allocate(xs(1:ni),stat=ierr)                                         2524
      jerr=jerr+ierr                                                       2525
      if(jerr.ne.0) return                                                 2526
      call spchkvars(no,ni,x,ix,ju)                                        2527
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2528
      if(maxval(ju) .gt. 0)goto 17121                                      2528
      jerr=7777                                                            2528
      go to 11790                                                          2528
17121 continue                                                             2529
      vq=max(0.0,vp)                                                       2529
      vq=vq*ni/sum(vq)                                                     2530
      ww=max(0.0,w)                                                        2530
      sw=sum(ww)                                                           2530
      if(sw .gt. 0.0)goto 17141                                            2530
      jerr=9999                                                            2530
      go to 11790                                                          2530
17141 continue                                                             2531
      ww=ww/sw                                                             2532
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     2533
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmin,u   2535 
     *lam,thr,  isd,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11790                                            2535
      dev0=2.0*sw*dev0                                                     2536
17150 do 17151 k=1,lmu                                                     2536
      nk=nin(k)                                                            2537
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2538
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2539
17151 continue                                                             2540
17152 continue                                                             2540
11790 continue                                                             2540
      deallocate(ww,ju,vq,xm,xs)                                           2541
      return                                                               2542
      end                                                                  2543
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,ne,nx,nlam,fl   2545 
     *min,ulam,  shri,isd,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,j
     *err)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2546 
     *9)
      real x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),xs(ni)          2547
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2548
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2549
      real, dimension (:), allocatable :: qy,t,w,wr,v,a,as,xm,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(a(1:ni),stat=jerr)                                          2554
      allocate(as(1:ni),stat=ierr)                                         2554
      jerr=jerr+ierr                                                       2555
      allocate(t(1:no),stat=ierr)                                          2555
      jerr=jerr+ierr                                                       2556
      allocate(mm(1:ni),stat=ierr)                                         2556
      jerr=jerr+ierr                                                       2557
      allocate(ga(1:ni),stat=ierr)                                         2557
      jerr=jerr+ierr                                                       2558
      allocate(ixx(1:ni),stat=ierr)                                        2558
      jerr=jerr+ierr                                                       2559
      allocate(wr(1:no),stat=ierr)                                         2559
      jerr=jerr+ierr                                                       2560
      allocate(v(1:ni),stat=ierr)                                          2560
      jerr=jerr+ierr                                                       2561
      allocate(xm(1:ni),stat=ierr)                                         2561
      jerr=jerr+ierr                                                       2562
      allocate(w(1:no),stat=ierr)                                          2562
      jerr=jerr+ierr                                                       2563
      allocate(qy(1:no),stat=ierr)                                         2563
      jerr=jerr+ierr                                                       2564
      if(jerr.ne.0) return                                                 2565
      bta=parm                                                             2565
      omb=1.0-bta                                                          2565
      fmax=log(huge(bta)*0.1)                                              2566
      qy=q*y                                                               2566
      yb=sum(qy)                                                           2567
      if(nonzero(no,g) .ne. 0)goto 17171                                   2567
      w=q*yb                                                               2567
      az=log(yb)                                                           2567
      uu=az                                                                2568
      xm=yb*xb                                                             2568
      t=0.0                                                                2569
      goto 17181                                                           2570
17171 continue                                                             2570
      w=q*exp(sign(min(abs(g),fmax),g))                                    2570
      ww=sum(w)                                                            2570
      eaz=yb/ww                                                            2571
      w=eaz*w                                                              2571
      az=log(eaz)                                                          2571
      uu=az                                                                2571
      t=g                                                                  2572
17190 do 17191 j=1,ni                                                      2572
      if(ju(j).eq.0)goto 17191                                             2572
      jb=ix(j)                                                             2572
      je=ix(j+1)-1                                                         2573
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2574
17191 continue                                                             2575
17192 continue                                                             2575
17181 continue                                                             2576
17161 continue                                                             2576
      tt=yb*uu                                                             2576
      ww=yb                                                                2576
      wr=qy-q*(yb*(1.0-uu))                                                2576
      a=0.0                                                                2576
      as=0.0                                                               2577
      dv0=yb*(log(yb)-1.0)                                                 2577
      dvr=-yb                                                              2578
17200 do 17201 i=1,no                                                      2578
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             2578
17201 continue                                                             2578
17202 continue                                                             2578
      dvr=dvr-dv0                                                          2578
      dev0=dvr                                                             2579
      if(flmin .ge. 1.0)goto 17221                                         2579
      eqs=max(eps,flmin)                                                   2579
      alf=eqs**(1.0/(nlam-1))                                              2579
17221 continue                                                             2580
      m=0                                                                  2580
      mm=0                                                                 2580
      nlp=0                                                                2580
      nin=nlp                                                              2580
      mnl=min(mnlam,nlam)                                                  2580
      shr=shri*dev0                                                        2580
      al=0.0                                                               2580
      ixx=0                                                                2581
17230 do 17231 j=1,ni                                                      2581
      if(ju(j).eq.0)goto 17231                                             2582
      jb=ix(j)                                                             2582
      je=ix(j+1)-1                                                         2583
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   2585 
     *)-xb(j)*tt)/xs(j)
17231 continue                                                             2586
17232 continue                                                             2586
17240 do 17241 ilm=1,nlam                                                  2586
      al0=al                                                               2587
      if(flmin .lt. 1.0)goto 17261                                         2587
      al=ulam(ilm)                                                         2587
      goto 17251                                                           2588
17261 if(ilm .le. 2)goto 17271                                             2588
      al=al*alf                                                            2588
      goto 17251                                                           2589
17271 if(ilm .ne. 1)goto 17281                                             2589
      al=big                                                               2589
      goto 17291                                                           2590
17281 continue                                                             2590
      al0=0.0                                                              2591
17300 do 17301 j=1,ni                                                      2591
      if(ju(j).eq.0)goto 17301                                             2591
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2591
17301 continue                                                             2592
17302 continue                                                             2592
      al0=al0/max(bta,1.0e-3)                                              2592
      al=alf*al0                                                           2593
17291 continue                                                             2594
17251 continue                                                             2594
      al2=al*omb                                                           2594
      al1=al*bta                                                           2594
      tlam=bta*(2.0*al-al0)                                                2595
17310 do 17311 k=1,ni                                                      2595
      if(ixx(k).eq.1)goto 17311                                            2595
      if(ju(k).eq.0)goto 17311                                             2596
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2597
17311 continue                                                             2598
17312 continue                                                             2598
10680 continue                                                             2599
17320 continue                                                             2599
17321 continue                                                             2599
      az0=az                                                               2600
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2601
17330 do 17331 j=1,ni                                                      2601
      if(ixx(j).eq.0)goto 17331                                            2601
      jb=ix(j)                                                             2601
      je=ix(j+1)-1                                                         2602
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2603
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   2605 
     *b(j)**2)/xs(j)**2
17331 continue                                                             2606
17332 continue                                                             2606
17340 continue                                                             2606
17341 continue                                                             2606
      nlp=nlp+1                                                            2607
      dlx=0.0                                                              2608
17350 do 17351 k=1,ni                                                      2608
      if(ixx(k).eq.0)goto 17351                                            2608
      jb=ix(k)                                                             2608
      je=ix(k+1)-1                                                         2608
      ak=a(k)                                                              2609
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2611 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2612
      if(au .gt. 0.0)goto 17371                                            2612
      a(k)=0.0                                                             2612
      goto 17381                                                           2613
17371 continue                                                             2613
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2613
17381 continue                                                             2614
17361 continue                                                             2614
      if(a(k).eq.ak)goto 17351                                             2615
      if(mm(k) .ne. 0)goto 17401                                           2615
      nin=nin+1                                                            2615
      if(nin.gt.nx)goto 17352                                              2616
      mm(k)=nin                                                            2616
      m(nin)=k                                                             2617
17401 continue                                                             2618
      d=a(k)-ak                                                            2618
      dlx=max(dlx,v(k)*d**2)                                               2618
      dv=d/xs(k)                                                           2619
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2620
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2621
      uu=uu-dv*xb(k)                                                       2621
      tt=tt-dv*xm(k)                                                       2622
17351 continue                                                             2623
17352 continue                                                             2623
      if(nin.gt.nx)goto 17342                                              2623
      d=tt/ww-uu                                                           2624
      az=az+d                                                              2624
      dlx=max(dlx,ww*d**2)                                                 2624
      uu=uu+d                                                              2625
      if(dlx.lt.shr)goto 17342                                             2625
      if(nlp .le. maxit)goto 17421                                         2625
      jerr=-ilm                                                            2625
      return                                                               2625
17421 continue                                                             2626
17430 continue                                                             2626
17431 continue                                                             2626
      nlp=nlp+1                                                            2626
      dlx=0.0                                                              2627
17440 do 17441 l=1,nin                                                     2627
      k=m(l)                                                               2628
      jb=ix(k)                                                             2628
      je=ix(k+1)-1                                                         2628
      ak=a(k)                                                              2629
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2631 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2632
      if(au .gt. 0.0)goto 17461                                            2632
      a(k)=0.0                                                             2632
      goto 17471                                                           2633
17461 continue                                                             2633
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2633
17471 continue                                                             2634
17451 continue                                                             2634
      if(a(k).eq.ak)goto 17441                                             2634
      d=a(k)-ak                                                            2634
      dlx=max(dlx,v(k)*d**2)                                               2635
      dv=d/xs(k)                                                           2635
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2636
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2637
      uu=uu-dv*xb(k)                                                       2637
      tt=tt-dv*xm(k)                                                       2638
17441 continue                                                             2639
17442 continue                                                             2639
      d=tt/ww-uu                                                           2639
      az=az+d                                                              2639
      dlx=max(dlx,ww*d**2)                                                 2639
      uu=uu+d                                                              2640
      if(dlx.lt.shr)goto 17432                                             2640
      if(nlp .le. maxit)goto 17491                                         2640
      jerr=-ilm                                                            2640
      return                                                               2640
17491 continue                                                             2641
      goto 17431                                                           2642
17432 continue                                                             2642
      goto 17341                                                           2643
17342 continue                                                             2643
      if(nin.gt.nx)goto 17322                                              2644
      euu=exp(sign(min(abs(uu),fmax),uu))                                  2645
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                2645
      ww=sum(w)                                                            2646
      wr=qy-w*(1.0-uu)                                                     2646
      tt=sum(wr)                                                           2647
      if(ww*(az-az0)**2 .ge. shr)goto 17511                                2647
      kx=0                                                                 2648
17520 do 17521 j=1,nin                                                     2648
      k=m(j)                                                               2649
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 17521                            2649
      kx=1                                                                 2649
      goto 17522                                                           2650
17521 continue                                                             2651
17522 continue                                                             2651
      if(kx .ne. 0)goto 17541                                              2652
17550 do 17551 j=1,ni                                                      2652
      if(ixx(j).eq.1)goto 17551                                            2652
      if(ju(j).eq.0)goto 17551                                             2653
      jb=ix(j)                                                             2653
      je=ix(j+1)-1                                                         2654
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2655
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   2657 
     *)-xb(j)*tt)/xs(j)
      if(ga(j) .le. al1*vp(j))goto 17571                                   2657
      ixx(j)=1                                                             2657
      kx=1                                                                 2657
17571 continue                                                             2658
17551 continue                                                             2659
17552 continue                                                             2659
      if(kx.eq.1) go to 10680                                              2660
      goto 17322                                                           2661
17541 continue                                                             2662
17511 continue                                                             2663
      goto 17321                                                           2664
17322 continue                                                             2664
      if(nin .le. nx)goto 17591                                            2664
      jerr=-10000-ilm                                                      2664
      goto 17242                                                           2664
17591 continue                                                             2665
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2665
      kin(ilm)=nin                                                         2666
      a0(ilm)=az                                                           2666
      alm(ilm)=al                                                          2666
      lmu=ilm                                                              2667
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        2668
      if(ilm.lt.mnl)goto 17241                                             2668
      if(flmin.ge.1.0)goto 17241                                           2669
      me=0                                                                 2669
17600 do 17601 j=1,nin                                                     2669
      if(ca(j,ilm).ne.0.0) me=me+1                                         2669
17601 continue                                                             2669
17602 continue                                                             2669
      if(me.gt.ne)goto 17242                                               2670
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17242              2671
      if(dev(ilm).gt.devmax)goto 17242                                     2672
17241 continue                                                             2673
17242 continue                                                             2673
      g=t+uu                                                               2674
11790 continue                                                             2674
      deallocate(t,w,wr,v,a,qy,xm,as,mm,ga,ixx)                            2675
      return                                                               2676
      end                                                                  2677
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       2678
      real x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)           2679
      integer ix(*),jx(*)                                                  2680
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 17621                                     2683
      jerr=8888                                                            2683
      return                                                               2683
17621 continue                                                             2684
      allocate(w(1:no),stat=jerr)                                          2685
      allocate(f(1:no),stat=ierr)                                          2685
      jerr=jerr+ierr                                                       2686
      if(jerr.ne.0) return                                                 2687
      w=max(0.0,q)                                                         2687
      sw=sum(w)                                                            2687
      if(sw .gt. 0.0)goto 17641                                            2687
      jerr=9999                                                            2687
      go to 11790                                                          2687
17641 continue                                                             2688
      yb=dot_product(w,y)/sw                                               2688
      fmax=log(huge(y(1))*0.1)                                             2689
17650 do 17651 lam=1,nlam                                                  2689
      f=a0(lam)                                                            2690
17660 do 17661 j=1,ni                                                      2690
      if(a(j,lam).eq.0.0)goto 17661                                        2690
      jb=ix(j)                                                             2690
      je=ix(j+1)-1                                                         2691
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          2692
17661 continue                                                             2693
17662 continue                                                             2693
      f=f+g                                                                2694
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2695
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2696
17651 continue                                                             2697
17652 continue                                                             2697
11790 continue                                                             2697
      deallocate(w,f)                                                      2698
      return                                                               2699
      end                                                                  2700
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   2701 
     *jerr)
      real x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(nlam)          2702
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 2703
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 17681                                     2706
      jerr=8888                                                            2706
      return                                                               2706
17681 continue                                                             2707
      allocate(w(1:no),stat=jerr)                                          2708
      allocate(f(1:no),stat=ierr)                                          2708
      jerr=jerr+ierr                                                       2709
      if(jerr.ne.0) return                                                 2710
      w=max(0.0,q)                                                         2710
      sw=sum(w)                                                            2710
      if(sw .gt. 0.0)goto 17701                                            2710
      jerr=9999                                                            2710
      go to 11790                                                          2710
17701 continue                                                             2711
      yb=dot_product(w,y)/sw                                               2711
      fmax=log(huge(y(1))*0.1)                                             2712
17710 do 17711 lam=1,nlam                                                  2712
      f=a0(lam)                                                            2713
17720 do 17721 k=1,nin(lam)                                                2713
      j=ia(k)                                                              2713
      jb=ix(j)                                                             2713
      je=ix(j+1)-1                                                         2714
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         2715
17721 continue                                                             2716
17722 continue                                                             2716
      f=f+g                                                                2717
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2718
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2719
17711 continue                                                             2720
17712 continue                                                             2720
11790 continue                                                             2720
      deallocate(w,f)                                                      2721
      return                                                               2722
      end                                                                  2723
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
