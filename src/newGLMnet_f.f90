c
c                          newGLMnet (6/1/11)
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
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp,sxpl                  
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
      q(i,ic)=min(max(emin*sxp(i),exp(dble(fi))),emax*sxp(i))              1557
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
      if(ix .ne. 0)goto 13431                                              1568
      ga=0.0                                                               1569
13440 do 13441 ic=1,nc                                                     1569
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1570
13450 do 13451 k=1,ni                                                      1570
      if(ixx(k).eq.1)goto 13451                                            1570
      if(ju(k).eq.0)goto 13451                                             1571
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                          1572
13451 continue                                                             1573
13452 continue                                                             1573
13441 continue                                                             1574
13442 continue                                                             1574
13460 do 13461 k=1,ni                                                      1574
      if(ixx(k).eq.1)goto 13461                                            1574
      if(ju(k).eq.0)goto 13461                                             1575
      if(ga(k) .le. al1*vp(k))goto 13481                                   1575
      ixx(k)=1                                                             1575
      ix=1                                                                 1575
13481 continue                                                             1576
13461 continue                                                             1577
13462 continue                                                             1577
      if(ix.eq.1) go to 10680                                              1578
      goto 12992                                                           1579
13431 continue                                                             1580
      goto 12991                                                           1581
12992 continue                                                             1581
      if(jx .le. 0)goto 13501                                              1581
      jerr=-10000-ilm                                                      1581
      goto 12912                                                           1581
13501 continue                                                             1581
      devi=0.0                                                             1582
13510 do 13511 ic=1,nc                                                     1583
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1583
      a0(ic,ilm)=b(0,ic)                                                   1584
13520 do 13521 i=1,no                                                      1584
      if(y(i,ic).le.0.0)goto 13521                                         1585
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1586
13521 continue                                                             1587
13522 continue                                                             1587
13511 continue                                                             1588
13512 continue                                                             1588
      kin(ilm)=nin                                                         1588
      alm(ilm)=al                                                          1588
      lmu=ilm                                                              1589
      dev(ilm)=(dev1-devi)/dev0                                            1589
      if(ig.eq.0)goto 12912                                                1590
      if(ilm.lt.mnl)goto 12911                                             1590
      if(flmin.ge.1.0)goto 12911                                           1591
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 12912             1592
      if(dev(ilm).gt.devmax)goto 12912                                     1592
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12912                             1593
12911 continue                                                             1594
12912 continue                                                             1594
      g=log(q)                                                             1594
13530 do 13531 i=1,no                                                      1594
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1594
13531 continue                                                             1595
13532 continue                                                             1595
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                           1596
      return                                                               1597
      end                                                                  1598
      subroutine kazero(kk,n,y,g,q,az,jerr)                                1599
      parameter(eps=1.0e-7)                                                1600
      real y(n,kk),g(n,kk),q(n),az(kk)                                     1601
      real, dimension (:), allocatable :: s                                     
      real, dimension (:,:), allocatable :: e                                   
      allocate(e(1:n,1:kk),stat=jerr)                                           
      allocate(s(1:n),stat=ierr)                                           1606
      jerr=jerr+ierr                                                       1607
      if(jerr.ne.0) return                                                 1608
      az=0.0                                                               1608
      e=exp(g)                                                             1608
13540 do 13541 i=1,n                                                       1608
      s(i)=sum(e(i,:))                                                     1608
13541 continue                                                             1609
13542 continue                                                             1609
13550 continue                                                             1609
13551 continue                                                             1609
      dm=0.0                                                               1610
13560 do 13561 k=1,kk                                                      1610
      t=0.0                                                                1610
      u=t                                                                  1611
13570 do 13571 i=1,n                                                       1611
      pik=e(i,k)/s(i)                                                      1612
      t=t+q(i)*(y(i,k)-pik)                                                1612
      u=u+q(i)*pik*(1.0-pik)                                               1613
13571 continue                                                             1614
13572 continue                                                             1614
      d=t/u                                                                1614
      az(k)=az(k)+d                                                        1614
      ed=exp(d)                                                            1614
      dm=max(dm,abs(d))                                                    1615
13580 do 13581 i=1,n                                                       1615
      z=e(i,k)                                                             1615
      e(i,k)=z*ed                                                          1615
      s(i)=s(i)-z+e(i,k)                                                   1615
13581 continue                                                             1616
13582 continue                                                             1616
13561 continue                                                             1617
13562 continue                                                             1617
      if(dm.lt.eps)goto 13552                                              1617
      goto 13551                                                           1618
13552 continue                                                             1618
      az=az-sum(az)/kk                                                     1619
      deallocate(e,s)                                                      1620
      return                                                               1621
      end                                                                  1622
      function elc(parm,n,a,m)                                             1623
      real a(n)                                                            1623
      integer m(n)                                                         1624
      fn=n                                                                 1624
      am=sum(a)/fn                                                         1625
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 13601                       1625
      elc=am                                                               1625
      return                                                               1625
13601 continue                                                             1626
13610 do 13611 i=1,n                                                       1626
      m(i)=i                                                               1626
13611 continue                                                             1626
13612 continue                                                             1626
      call psort7(a,m,1,n)                                                 1627
      if(a(m(1)) .ne. a(m(n)))goto 13631                                   1627
      elc=a(1)                                                             1627
      return                                                               1627
13631 continue                                                             1628
      if(mod(n,2) .ne. 1)goto 13651                                        1628
      ad=a(m(n/2+1))                                                       1628
      goto 13661                                                           1629
13651 continue                                                             1629
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       1629
13661 continue                                                             1630
13641 continue                                                             1630
      if(parm .ne. 1.0)goto 13681                                          1630
      elc=ad                                                               1630
      return                                                               1630
13681 continue                                                             1631
      b1=min(am,ad)                                                        1631
      b2=max(am,ad)                                                        1631
      k2=1                                                                 1632
13690 continue                                                             1632
13691 if(a(m(k2)).gt.b1)goto 13692                                         1632
      k2=k2+1                                                              1632
      goto 13691                                                           1632
13692 continue                                                             1632
      k1=k2-1                                                              1633
13700 continue                                                             1633
13701 if(a(m(k2)).ge.b2)goto 13702                                         1633
      k2=k2+1                                                              1633
      goto 13701                                                           1634
13702 continue                                                             1634
      r=parm/((1.0-parm)*fn)                                               1634
      is=0                                                                 1634
      sm=n-2*(k1-1)                                                        1635
13710 do 13711 k=k1,k2-1                                                   1635
      sm=sm-2.0                                                            1635
      s=r*sm+am                                                            1636
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 13731                   1636
      is=k                                                                 1636
      goto 13712                                                           1636
13731 continue                                                             1637
13711 continue                                                             1638
13712 continue                                                             1638
      if(is .eq. 0)goto 13751                                              1638
      elc=s                                                                1638
      return                                                               1638
13751 continue                                                             1638
      r2=2.0*r                                                             1638
      s1=a(m(k1))                                                          1638
      am2=2.0*am                                                           1639
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    1639
      elc=s1                                                               1640
13760 do 13761 k=k1+1,k2                                                   1640
      s=a(m(k))                                                            1640
      if(s.eq.s1)goto 13761                                                1641
      c=r2*sum(abs(a-s))+s*(s-am2)                                         1642
      if(c .ge. cri)goto 13781                                             1642
      cri=c                                                                1642
      elc=s                                                                1642
13781 continue                                                             1642
      s1=s                                                                 1643
13761 continue                                                             1644
13762 continue                                                             1644
      return                                                               1645
      end                                                                  1646
      function nintot(ni,nx,nc,a,m,nin,is)                                 1647
      real a(nx,nc)                                                        1647
      integer m(nx),is(ni)                                                 1648
      is=0                                                                 1648
      nintot=0                                                             1649
13790 do 13791 ic=1,nc                                                     1649
13800 do 13801 j=1,nin                                                     1649
      k=m(j)                                                               1649
      if(is(k).ne.0)goto 13801                                             1650
      if(a(j,ic).eq.0.0)goto 13801                                         1650
      is(k)=k                                                              1650
      nintot=nintot+1                                                      1651
13801 continue                                                             1651
13802 continue                                                             1651
13791 continue                                                             1652
13792 continue                                                             1652
      return                                                               1653
      end                                                                  1654
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             1655
      real ca(nx,nc),a(ni,nc)                                              1655
      integer ia(nx)                                                       1656
      a=0.0                                                                1657
13810 do 13811 ic=1,nc                                                     1657
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            1657
13811 continue                                                             1658
13812 continue                                                             1658
      return                                                               1659
      end                                                                  1660
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      1661
      real a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                             1661
      integer ia(nx)                                                       1662
13820 do 13821 i=1,nt                                                      1662
13830 do 13831 ic=1,nc                                                     1662
      ans(ic,i)=a0(ic)                                                     1664
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   1665 
     *:nin)))
13831 continue                                                             1665
13832 continue                                                             1665
13821 continue                                                             1666
13822 continue                                                             1666
      return                                                               1667
      end                                                                  1668
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,ne,nx,nlam,fl   1670 
     *min,  ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      real x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)                 1671
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1672
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1673
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 13851                                    1677
      jerr=10000                                                           1677
      return                                                               1677
13851 continue                                                             1678
      allocate(ww(1:no),stat=jerr)                                         1679
      allocate(ju(1:ni),stat=ierr)                                         1679
      jerr=jerr+ierr                                                       1680
      allocate(vq(1:ni),stat=ierr)                                         1680
      jerr=jerr+ierr                                                       1681
      allocate(xm(1:ni),stat=ierr)                                         1681
      jerr=jerr+ierr                                                       1682
      allocate(xs(1:ni),stat=ierr)                                         1682
      jerr=jerr+ierr                                                       1683
      if(jerr.ne.0) return                                                 1684
      call spchkvars(no,ni,x,ix,ju)                                        1685
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1686
      if(maxval(ju) .gt. 0)goto 13871                                      1686
      jerr=7777                                                            1686
      return                                                               1686
13871 continue                                                             1687
      vq=max(0.0,vp)                                                       1687
      vq=vq*ni/sum(vq)                                                     1688
13880 do 13881 i=1,no                                                      1688
      ww(i)=sum(y(i,:))                                                    1688
      y(i,:)=y(i,:)/ww(i)                                                  1688
13881 continue                                                             1688
13882 continue                                                             1688
      sw=sum(ww)                                                           1688
      ww=ww/sw                                                             1689
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     1690
      if(nc .ne. 1)goto 13901                                              1691
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,ne,nx,n   1693 
     *lam,flmin,  ulam,thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,de
     *v,alm,nlp,jerr)
      goto 13911                                                           1694
13901 continue                                                             1695
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmi   1697 
     *n,ulam,  thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nl
     *p,jerr)
13911 continue                                                             1698
13891 continue                                                             1698
      if(jerr.gt.0) return                                                 1698
      dev0=2.0*sw*dev0                                                     1699
13920 do 13921 k=1,lmu                                                     1699
      nk=nin(k)                                                            1700
13930 do 13931 ic=1,nc                                                     1700
      if(isd .le. 0)goto 13951                                             1700
13960 do 13961 l=1,nk                                                      1700
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1700
13961 continue                                                             1700
13962 continue                                                             1700
13951 continue                                                             1701
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1702
13931 continue                                                             1703
13932 continue                                                             1703
13921 continue                                                             1704
13922 continue                                                             1704
      deallocate(ww,ju,vq,xm,xs)                                           1705
      return                                                               1706
      end                                                                  1707
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,xm,xs)                1708
      real x(*),w(no),xm(ni),xs(ni)                                        1708
      integer ix(*),jx(*),ju(ni)                                           1709
13970 do 13971 j=1,ni                                                      1709
      if(ju(j).eq.0)goto 13971                                             1709
      jb=ix(j)                                                             1709
      je=ix(j+1)-1                                                         1710
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1711
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   1712 
     *)**2)
13971 continue                                                             1713
13972 continue                                                             1713
      if(isd.eq.0) xs=1.0                                                  1714
      return                                                               1715
      end                                                                  1716
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam,   1718 
     *  flmin,ulam,shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,alm
     *,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1720 
     *5, devmax=0.999)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        1721
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1722
      real xb(ni),xs(ni)                                                   1722
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1723
      real, dimension (:), allocatable :: xm,b,bs,v,r,sc,xv,q,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(b(0:ni),stat=jerr)                                          1728
      allocate(xm(0:ni),stat=ierr)                                         1728
      jerr=jerr+ierr                                                       1729
      allocate(xv(1:ni),stat=ierr)                                         1729
      jerr=jerr+ierr                                                       1730
      allocate(bs(0:ni),stat=ierr)                                         1730
      jerr=jerr+ierr                                                       1731
      allocate(ga(1:ni),stat=ierr)                                         1731
      jerr=jerr+ierr                                                       1732
      allocate(mm(1:ni),stat=ierr)                                         1732
      jerr=jerr+ierr                                                       1733
      allocate(ixx(1:ni),stat=ierr)                                        1733
      jerr=jerr+ierr                                                       1734
      allocate(q(1:no),stat=ierr)                                          1734
      jerr=jerr+ierr                                                       1735
      allocate(r(1:no),stat=ierr)                                          1735
      jerr=jerr+ierr                                                       1736
      allocate(v(1:no),stat=ierr)                                          1736
      jerr=jerr+ierr                                                       1737
      allocate(sc(1:no),stat=ierr)                                         1737
      jerr=jerr+ierr                                                       1738
      if(jerr.ne.0) return                                                 1739
      fmax=log(1.0/pmin-1.0)                                               1739
      fmin=-fmax                                                           1739
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1740
      bta=parm                                                             1740
      omb=1.0-bta                                                          1741
      q0=dot_product(w,y)                                                  1741
      if(q0 .gt. pmin)goto 13991                                           1741
      jerr=8001                                                            1741
      return                                                               1741
13991 continue                                                             1742
      if(q0 .lt. 1.0-pmin)goto 14011                                       1742
      jerr=9001                                                            1742
      return                                                               1742
14011 continue                                                             1742
      bz=log(q0/(1.0-q0))                                                  1743
      if(nonzero(no,g) .ne. 0)goto 14031                                   1743
      vi=q0*(1.0-q0)                                                       1743
      b(0)=bz                                                              1743
      v=vi*w                                                               1744
      r=w*(y-q0)                                                           1744
      q=q0                                                                 1744
      xm(0)=vi                                                             1744
      dev1=-(bz*q0+log(1.0-q0))                                            1745
      goto 14041                                                           1746
14031 continue                                                             1746
      b(0)=azero(no,y,g,w,jerr)                                            1746
      if(jerr.ne.0) return                                                 1747
      q=1.0/(1.0+exp(-b(0)-g))                                             1747
      v=w*q*(1.0-q)                                                        1747
      r=w*(y-q)                                                            1747
      xm(0)=sum(v)                                                         1748
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1749
14041 continue                                                             1750
14021 continue                                                             1750
      if(kopt .le. 0)goto 14061                                            1751
      if(isd .le. 0)goto 14081                                             1751
      xv=0.25                                                              1751
      goto 14091                                                           1752
14081 continue                                                             1753
14100 do 14101 j=1,ni                                                      1753
      if(ju(j).eq.0)goto 14101                                             1753
      jb=ix(j)                                                             1753
      je=ix(j+1)-1                                                         1754
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          1755
14101 continue                                                             1756
14102 continue                                                             1756
14091 continue                                                             1757
14071 continue                                                             1757
14061 continue                                                             1758
      b(1:ni)=0.0                                                          1758
      dev0=dev1                                                            1759
14110 do 14111 i=1,no                                                      1759
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1760
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1761
14111 continue                                                             1762
14112 continue                                                             1762
      if(flmin .ge. 1.0)goto 14131                                         1762
      eqs=max(eps,flmin)                                                   1762
      alf=eqs**(1.0/(nlam-1))                                              1762
14131 continue                                                             1763
      m=0                                                                  1763
      mm=0                                                                 1763
      nin=0                                                                1763
      o=0.0                                                                1763
      svr=o                                                                1763
      mnl=min(mnlam,nlam)                                                  1763
      bs=0.0                                                               1763
      nlp=0                                                                1763
      nin=nlp                                                              1764
      shr=shri*dev0                                                        1764
      al=0.0                                                               1764
      ixx=0                                                                1765
14140 do 14141 j=1,ni                                                      1765
      if(ju(j).eq.0)goto 14141                                             1766
      jb=ix(j)                                                             1766
      je=ix(j+1)-1                                                         1766
      jn=ix(j+1)-ix(j)                                                     1767
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1768
      gj=dot_product(sc(1:jn),x(jb:je))                                    1769
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      1770
14141 continue                                                             1771
14142 continue                                                             1771
14150 do 14151 ilm=1,nlam                                                  1771
      al0=al                                                               1772
      if(flmin .lt. 1.0)goto 14171                                         1772
      al=ulam(ilm)                                                         1772
      goto 14161                                                           1773
14171 if(ilm .le. 2)goto 14181                                             1773
      al=al*alf                                                            1773
      goto 14161                                                           1774
14181 if(ilm .ne. 1)goto 14191                                             1774
      al=big                                                               1774
      goto 14201                                                           1775
14191 continue                                                             1775
      al0=0.0                                                              1776
14210 do 14211 j=1,ni                                                      1776
      if(ju(j).eq.0)goto 14211                                             1776
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1776
14211 continue                                                             1777
14212 continue                                                             1777
      al0=al0/max(bta,1.0e-3)                                              1777
      al=alf*al0                                                           1778
14201 continue                                                             1779
14161 continue                                                             1779
      al2=al*omb                                                           1779
      al1=al*bta                                                           1779
      tlam=bta*(2.0*al-al0)                                                1780
14220 do 14221 k=1,ni                                                      1780
      if(ixx(k).eq.1)goto 14221                                            1780
      if(ju(k).eq.0)goto 14221                                             1781
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1782
14221 continue                                                             1783
14222 continue                                                             1783
10680 continue                                                             1784
14230 continue                                                             1784
14231 continue                                                             1784
      bs(0)=b(0)                                                           1784
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1785
14240 do 14241 j=1,ni                                                      1785
      if(ixx(j).eq.0)goto 14241                                            1786
      jb=ix(j)                                                             1786
      je=ix(j+1)-1                                                         1786
      jn=ix(j+1)-ix(j)                                                     1787
      sc(1:jn)=v(jx(jb:je))                                                1788
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 1789
      if(kopt .ne. 0)goto 14261                                            1790
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              1791
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                1792
14261 continue                                                             1793
14241 continue                                                             1794
14242 continue                                                             1794
14270 continue                                                             1794
14271 continue                                                             1794
      nlp=nlp+1                                                            1794
      dlx=0.0                                                              1795
14280 do 14281 k=1,ni                                                      1795
      if(ixx(k).eq.0)goto 14281                                            1796
      jb=ix(k)                                                             1796
      je=ix(k+1)-1                                                         1796
      jn=ix(k+1)-ix(k)                                                     1796
      bk=b(k)                                                              1797
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1798
      gk=dot_product(sc(1:jn),x(jb:je))                                    1799
      gk=(gk-svr*xb(k))/xs(k)                                              1800
      u=gk+xv(k)*b(k)                                                      1800
      au=abs(u)-vp(k)*al1                                                  1801
      if(au .gt. 0.0)goto 14301                                            1801
      b(k)=0.0                                                             1801
      goto 14311                                                           1802
14301 continue                                                             1802
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1802
14311 continue                                                             1803
14291 continue                                                             1803
      d=b(k)-bk                                                            1803
      if(abs(d).le.0.0)goto 14281                                          1803
      dlx=max(dlx,xv(k)*d**2)                                              1804
      if(mm(k) .ne. 0)goto 14331                                           1804
      nin=nin+1                                                            1804
      if(nin.gt.nx)goto 14282                                              1805
      mm(k)=nin                                                            1805
      m(nin)=k                                                             1805
      sc(1:jn)=v(jx(jb:je))                                                1806
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 1807
14331 continue                                                             1808
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1809
      o=o+d*(xb(k)/xs(k))                                                  1810
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1811
14281 continue                                                             1812
14282 continue                                                             1812
      if(nin.gt.nx)goto 14272                                              1813
      d=svr/xm(0)                                                          1814
      if(d .eq. 0.0)goto 14351                                             1814
      b(0)=b(0)+d                                                          1814
      dlx=max(dlx,xm(0)*d**2)                                              1814
      r=r-d*v                                                              1814
14351 continue                                                             1815
      svr=svr-d*xm(0)                                                      1815
      if(dlx.lt.shr)goto 14272                                             1816
      if(nlp .le. maxit)goto 14371                                         1816
      jerr=-ilm                                                            1816
      return                                                               1816
14371 continue                                                             1817
14380 continue                                                             1817
14381 continue                                                             1817
      nlp=nlp+1                                                            1817
      dlx=0.0                                                              1818
14390 do 14391 l=1,nin                                                     1818
      k=m(l)                                                               1818
      jb=ix(k)                                                             1818
      je=ix(k+1)-1                                                         1819
      jn=ix(k+1)-ix(k)                                                     1819
      bk=b(k)                                                              1820
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1821
      gk=dot_product(sc(1:jn),x(jb:je))                                    1822
      gk=(gk-svr*xb(k))/xs(k)                                              1823
      u=gk+xv(k)*b(k)                                                      1823
      au=abs(u)-vp(k)*al1                                                  1824
      if(au .gt. 0.0)goto 14411                                            1824
      b(k)=0.0                                                             1824
      goto 14421                                                           1825
14411 continue                                                             1825
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1825
14421 continue                                                             1826
14401 continue                                                             1826
      d=b(k)-bk                                                            1826
      if(abs(d).le.0.0)goto 14391                                          1826
      dlx=max(dlx,xv(k)*d**2)                                              1827
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1828
      o=o+d*(xb(k)/xs(k))                                                  1829
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1830
14391 continue                                                             1831
14392 continue                                                             1831
      d=svr/xm(0)                                                          1832
      if(d .eq. 0.0)goto 14441                                             1832
      b(0)=b(0)+d                                                          1832
      dlx=max(dlx,xm(0)*d**2)                                              1832
      r=r-d*v                                                              1832
14441 continue                                                             1833
      svr=svr-d*xm(0)                                                      1833
      if(dlx.lt.shr)goto 14382                                             1834
      if(nlp .le. maxit)goto 14461                                         1834
      jerr=-ilm                                                            1834
      return                                                               1834
14461 continue                                                             1835
      goto 14381                                                           1836
14382 continue                                                             1836
      goto 14271                                                           1837
14272 continue                                                             1837
      if(nin.gt.nx)goto 14232                                              1838
      sc=b(0)                                                              1838
      b0=0.0                                                               1839
14470 do 14471 j=1,nin                                                     1839
      l=m(j)                                                               1839
      jb=ix(l)                                                             1839
      je=ix(l+1)-1                                                         1840
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      1841
      b0=b0-b(l)*xb(l)/xs(l)                                               1842
14471 continue                                                             1843
14472 continue                                                             1843
      sc=sc+b0                                                             1844
14480 do 14481 i=1,no                                                      1844
      fi=sc(i)+g(i)                                                        1845
      if(fi .ge. fmin)goto 14501                                           1845
      q(i)=0.0                                                             1845
      goto 14491                                                           1845
14501 if(fi .le. fmax)goto 14511                                           1845
      q(i)=1.0                                                             1845
      goto 14521                                                           1846
14511 continue                                                             1846
      q(i)=1.0/(1.0+exp(-fi))                                              1846
14521 continue                                                             1847
14491 continue                                                             1847
14481 continue                                                             1848
14482 continue                                                             1848
      v=w*q*(1.0-q)                                                        1848
      xm(0)=sum(v)                                                         1848
      if(xm(0).lt.vmin)goto 14232                                          1849
      r=w*(y-q)                                                            1849
      svr=sum(r)                                                           1849
      o=0.0                                                                1850
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 14541                         1850
      kx=0                                                                 1851
14550 do 14551 j=1,nin                                                     1851
      k=m(j)                                                               1852
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 14551                           1852
      kx=1                                                                 1852
      goto 14552                                                           1853
14551 continue                                                             1854
14552 continue                                                             1854
      if(kx .ne. 0)goto 14571                                              1855
14580 do 14581 j=1,ni                                                      1855
      if(ixx(j).eq.1)goto 14581                                            1855
      if(ju(j).eq.0)goto 14581                                             1856
      jb=ix(j)                                                             1856
      je=ix(j+1)-1                                                         1856
      jn=ix(j+1)-ix(j)                                                     1857
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1858
      gj=dot_product(sc(1:jn),x(jb:je))                                    1859
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      1860
      if(ga(j) .le. al1*vp(j))goto 14601                                   1860
      ixx(j)=1                                                             1860
      kx=1                                                                 1860
14601 continue                                                             1861
14581 continue                                                             1862
14582 continue                                                             1862
      if(kx.eq.1) go to 10680                                              1863
      goto 14232                                                           1864
14571 continue                                                             1865
14541 continue                                                             1866
      goto 14231                                                           1867
14232 continue                                                             1867
      if(nin .le. nx)goto 14621                                            1867
      jerr=-10000-ilm                                                      1867
      goto 14152                                                           1867
14621 continue                                                             1868
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1868
      kin(ilm)=nin                                                         1869
      a0(ilm)=b(0)                                                         1869
      alm(ilm)=al                                                          1869
      lmu=ilm                                                              1870
      devi=dev2(no,w,y,q,pmin)                                             1871
      dev(ilm)=(dev1-devi)/dev0                                            1872
      if(ilm.lt.mnl)goto 14151                                             1872
      if(flmin.ge.1.0)goto 14151                                           1873
      me=0                                                                 1873
14630 do 14631 j=1,nin                                                     1873
      if(a(j,ilm).ne.0.0) me=me+1                                          1873
14631 continue                                                             1873
14632 continue                                                             1873
      if(me.gt.ne)goto 14152                                               1874
      if(dev(ilm).gt.devmax)goto 14152                                     1874
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 14152                             1875
      if(xm(0).lt.vmin)goto 14152                                          1876
14151 continue                                                             1877
14152 continue                                                             1877
      g=log(q/(1.0-q))                                                     1878
      deallocate(xm,b,bs,v,r,sc,xv,q,mm,ga,ixx)                            1879
      return                                                               1880
      end                                                                  1881
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam   1883 
     *,flmin,ulam,  shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,al
     *m,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1885 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni)    1886
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1887
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1888
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp,sxpl                  
      real, dimension (:), allocatable :: sc,xm,v,r,ga                          
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is,iy                           
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(xm(0:ni),stat=ierr)                                         1899
      jerr=jerr+ierr                                                       1900
      allocate(r(1:no),stat=ierr)                                          1900
      jerr=jerr+ierr                                                       1901
      allocate(v(1:no),stat=ierr)                                          1901
      jerr=jerr+ierr                                                       1902
      allocate(mm(1:ni),stat=ierr)                                         1902
      jerr=jerr+ierr                                                       1903
      allocate(ga(1:ni),stat=ierr)                                         1903
      jerr=jerr+ierr                                                       1904
      allocate(iy(1:ni),stat=ierr)                                         1904
      jerr=jerr+ierr                                                       1905
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1905
      jerr=jerr+ierr                                                       1906
      allocate(sxp(1:no),stat=ierr)                                        1906
      jerr=jerr+ierr                                                       1907
      allocate(sxpl(1:no),stat=ierr)                                       1907
      jerr=jerr+ierr                                                       1908
      allocate(sc(1:no),stat=ierr)                                         1908
      jerr=jerr+ierr                                                       1909
      if(jerr.ne.0) return                                                 1910
      pmax=1.0-pmin                                                        1910
      emin=pmin/pmax                                                       1910
      emax=1.0/emin                                                        1911
      pfm=(1.0+pmin)*pmin                                                  1911
      pfx=(1.0-pmin)*pmax                                                  1911
      vmin=pfm*pmax                                                        1912
      bta=parm                                                             1912
      omb=1.0-bta                                                          1912
      dev1=0.0                                                             1912
      dev0=0.0                                                             1913
14640 do 14641 ic=1,nc                                                     1913
      q0=dot_product(w,y(:,ic))                                            1914
      if(q0 .gt. pmin)goto 14661                                           1914
      jerr =8000+ic                                                        1914
      return                                                               1914
14661 continue                                                             1915
      if(q0 .lt. 1.0-pmin)goto 14681                                       1915
      jerr =9000+ic                                                        1915
      return                                                               1915
14681 continue                                                             1916
      b(1:ni,ic)=0.0                                                       1916
      b(0,ic)=log(q0)                                                      1916
      dev1=dev1-q0*b(0,ic)                                                 1917
14641 continue                                                             1918
14642 continue                                                             1918
      iy=0                                                                 1918
      al=0.0                                                               1919
      if(nonzero(no*nc,g) .ne. 0)goto 14701                                1920
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1920
      sxp=0.0                                                              1921
14710 do 14711 ic=1,nc                                                     1921
      q(:,ic)=exp(b(0,ic))                                                 1921
      sxp=sxp+q(:,ic)                                                      1921
14711 continue                                                             1922
14712 continue                                                             1922
      goto 14721                                                           1923
14701 continue                                                             1923
14730 do 14731 i=1,no                                                      1923
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1923
14731 continue                                                             1923
14732 continue                                                             1923
      sxp=0.0                                                              1924
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1924
      if(jerr.ne.0) return                                                 1925
      dev1=0.0                                                             1926
14740 do 14741 ic=1,nc                                                     1926
      q(:,ic)=b(0,ic)+g(:,ic)                                              1927
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             1928
      q(:,ic)=exp(q(:,ic))                                                 1928
      sxp=sxp+q(:,ic)                                                      1929
14741 continue                                                             1930
14742 continue                                                             1930
      sxpl=w*log(sxp)                                                      1930
14750 do 14751 ic=1,nc                                                     1930
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  1930
14751 continue                                                             1931
14752 continue                                                             1931
14721 continue                                                             1932
14691 continue                                                             1932
14760 do 14761 ic=1,nc                                                     1932
14770 do 14771 i=1,no                                                      1932
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1932
14771 continue                                                             1932
14772 continue                                                             1932
14761 continue                                                             1933
14762 continue                                                             1933
      dev0=dev0+dev1                                                       1934
      if(kopt .le. 0)goto 14791                                            1935
      if(isd .le. 0)goto 14811                                             1935
      xv=0.25                                                              1935
      goto 14821                                                           1936
14811 continue                                                             1937
14830 do 14831 j=1,ni                                                      1937
      if(ju(j).eq.0)goto 14831                                             1937
      jb=ix(j)                                                             1937
      je=ix(j+1)-1                                                         1938
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        1939
14831 continue                                                             1940
14832 continue                                                             1940
14821 continue                                                             1941
14801 continue                                                             1941
14791 continue                                                             1942
      if(flmin .ge. 1.0)goto 14851                                         1942
      eqs=max(eps,flmin)                                                   1942
      alf=eqs**(1.0/(nlam-1))                                              1942
14851 continue                                                             1943
      m=0                                                                  1943
      mm=0                                                                 1943
      nin=0                                                                1943
      nlp=0                                                                1943
      mnl=min(mnlam,nlam)                                                  1943
      bs=0.0                                                               1943
      svr=0.0                                                              1943
      o=0.0                                                                1944
      shr=shri*dev0                                                        1944
      ga=0.0                                                               1945
14860 do 14861 ic=1,nc                                                     1945
      v=q(:,ic)/sxp                                                        1945
      r=w*(y(:,ic)-v)                                                      1945
      v=w*v*(1.0-v)                                                        1946
14870 do 14871 j=1,ni                                                      1946
      if(ju(j).eq.0)goto 14871                                             1947
      jb=ix(j)                                                             1947
      je=ix(j+1)-1                                                         1947
      jn=ix(j+1)-ix(j)                                                     1948
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1949
      gj=dot_product(sc(1:jn),x(jb:je))                                    1950
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             1951
14871 continue                                                             1952
14872 continue                                                             1952
14861 continue                                                             1953
14862 continue                                                             1953
14880 do 14881 ilm=1,nlam                                                  1953
      al0=al                                                               1954
      if(flmin .lt. 1.0)goto 14901                                         1954
      al=ulam(ilm)                                                         1954
      goto 14891                                                           1955
14901 if(ilm .le. 2)goto 14911                                             1955
      al=al*alf                                                            1955
      goto 14891                                                           1956
14911 if(ilm .ne. 1)goto 14921                                             1956
      al=big                                                               1956
      goto 14931                                                           1957
14921 continue                                                             1957
      al0=0.0                                                              1958
14940 do 14941 j=1,ni                                                      1958
      if(ju(j).eq.0)goto 14941                                             1958
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1958
14941 continue                                                             1959
14942 continue                                                             1959
      al0=al0/max(bta,1.0e-3)                                              1959
      al=alf*al0                                                           1960
14931 continue                                                             1961
14891 continue                                                             1961
      al2=al*omb                                                           1961
      al1=al*bta                                                           1961
      tlam=bta*(2.0*al-al0)                                                1962
14950 do 14951 k=1,ni                                                      1962
      if(iy(k).eq.1)goto 14951                                             1962
      if(ju(k).eq.0)goto 14951                                             1963
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      1964
14951 continue                                                             1965
14952 continue                                                             1965
10680 continue                                                             1966
14960 continue                                                             1966
14961 continue                                                             1966
      ixx=0                                                                1966
      jxx=ixx                                                              1966
      ig=0                                                                 1967
14970 do 14971 ic=1,nc                                                     1967
      bs(0,ic)=b(0,ic)                                                     1968
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1969
      xm(0)=0.0                                                            1969
      svr=0.0                                                              1969
      o=0.0                                                                1970
14980 do 14981 i=1,no                                                      1970
      pic=q(i,ic)/sxp(i)                                                   1971
      if(pic .ge. pfm)goto 15001                                           1971
      pic=0.0                                                              1971
      v(i)=0.0                                                             1971
      goto 14991                                                           1972
15001 if(pic .le. pfx)goto 15011                                           1972
      pic=1.0                                                              1972
      v(i)=0.0                                                             1972
      goto 15021                                                           1973
15011 continue                                                             1973
      v(i)=w(i)*pic*(1.0-pic)                                              1973
      xm(0)=xm(0)+v(i)                                                     1973
15021 continue                                                             1974
14991 continue                                                             1974
      r(i)=w(i)*(y(i,ic)-pic)                                              1974
      svr=svr+r(i)                                                         1975
14981 continue                                                             1976
14982 continue                                                             1976
      if(xm(0).le.vmin)goto 14971                                          1976
      ig=1                                                                 1977
15030 do 15031 j=1,ni                                                      1977
      if(iy(j).eq.0)goto 15031                                             1978
      jb=ix(j)                                                             1978
      je=ix(j+1)-1                                                         1979
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             1980
      if(kopt .ne. 0)goto 15051                                            1981
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       1982
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          1983
15051 continue                                                             1984
15031 continue                                                             1985
15032 continue                                                             1985
15060 continue                                                             1985
15061 continue                                                             1985
      nlp=nlp+1                                                            1985
      dlx=0.0                                                              1986
15070 do 15071 k=1,ni                                                      1986
      if(iy(k).eq.0)goto 15071                                             1987
      jb=ix(k)                                                             1987
      je=ix(k+1)-1                                                         1987
      jn=ix(k+1)-ix(k)                                                     1987
      bk=b(k,ic)                                                           1988
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1989
      gk=dot_product(sc(1:jn),x(jb:je))                                    1990
      gk=(gk-svr*xb(k))/xs(k)                                              1991
      u=gk+xv(k,ic)*b(k,ic)                                                1991
      au=abs(u)-vp(k)*al1                                                  1992
      if(au .gt. 0.0)goto 15091                                            1992
      b(k,ic)=0.0                                                          1992
      goto 15101                                                           1993
15091 continue                                                             1993
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1993
15101 continue                                                             1994
15081 continue                                                             1994
      d=b(k,ic)-bk                                                         1994
      if(abs(d).le.0.0)goto 15071                                          1995
      dlx=max(dlx,xv(k,ic)*d**2)                                           1996
      if(mm(k) .ne. 0)goto 15121                                           1996
      nin=nin+1                                                            1997
      if(nin .le. nx)goto 15141                                            1997
      jxx=1                                                                1997
      goto 15072                                                           1997
15141 continue                                                             1998
      mm(k)=nin                                                            1998
      m(nin)=k                                                             1999
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             2000
15121 continue                                                             2001
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2002
      o=o+d*(xb(k)/xs(k))                                                  2003
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2004
15071 continue                                                             2005
15072 continue                                                             2005
      if(jxx.gt.0)goto 15062                                               2006
      d=svr/xm(0)                                                          2007
      if(d .eq. 0.0)goto 15161                                             2007
      b(0,ic)=b(0,ic)+d                                                    2007
      dlx=max(dlx,xm(0)*d**2)                                              2008
      r=r-d*v                                                              2008
      svr=svr-d*xm(0)                                                      2009
15161 continue                                                             2010
      if(dlx.lt.shr)goto 15062                                             2010
      if(nlp .le. maxit)goto 15181                                         2010
      jerr=-ilm                                                            2010
      return                                                               2010
15181 continue                                                             2011
15190 continue                                                             2011
15191 continue                                                             2011
      nlp=nlp+1                                                            2011
      dlx=0.0                                                              2012
15200 do 15201 l=1,nin                                                     2012
      k=m(l)                                                               2012
      jb=ix(k)                                                             2012
      je=ix(k+1)-1                                                         2013
      jn=ix(k+1)-ix(k)                                                     2013
      bk=b(k,ic)                                                           2014
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2015
      gk=dot_product(sc(1:jn),x(jb:je))                                    2016
      gk=(gk-svr*xb(k))/xs(k)                                              2017
      u=gk+xv(k,ic)*b(k,ic)                                                2017
      au=abs(u)-vp(k)*al1                                                  2018
      if(au .gt. 0.0)goto 15221                                            2018
      b(k,ic)=0.0                                                          2018
      goto 15231                                                           2019
15221 continue                                                             2019
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              2019
15231 continue                                                             2020
15211 continue                                                             2020
      d=b(k,ic)-bk                                                         2020
      if(abs(d).le.0.0)goto 15201                                          2021
      dlx=max(dlx,xv(k,ic)*d**2)                                           2022
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2023
      o=o+d*(xb(k)/xs(k))                                                  2024
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2025
15201 continue                                                             2026
15202 continue                                                             2026
      d=svr/xm(0)                                                          2027
      if(d .eq. 0.0)goto 15251                                             2027
      b(0,ic)=b(0,ic)+d                                                    2027
      dlx=max(dlx,xm(0)*d**2)                                              2028
      r=r-d*v                                                              2028
      svr=svr-d*xm(0)                                                      2029
15251 continue                                                             2030
      if(dlx.lt.shr)goto 15192                                             2030
      if(nlp .le. maxit)goto 15271                                         2030
      jerr=-ilm                                                            2030
      return                                                               2030
15271 continue                                                             2031
      goto 15191                                                           2032
15192 continue                                                             2032
      goto 15061                                                           2033
15062 continue                                                             2033
      if(jxx.gt.0)goto 14972                                               2034
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                         2035
      if(ixx .ne. 0)goto 15291                                             2036
15300 do 15301 j=1,nin                                                     2036
      k=m(j)                                                               2037
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 15321                2037
      ixx=1                                                                2037
      goto 15302                                                           2037
15321 continue                                                             2038
15301 continue                                                             2039
15302 continue                                                             2039
15291 continue                                                             2040
      sc=b(0,ic)+g(:,ic)                                                   2040
      b0=0.0                                                               2041
15330 do 15331 j=1,nin                                                     2041
      l=m(j)                                                               2041
      jb=ix(l)                                                             2041
      je=ix(l+1)-1                                                         2042
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   2043
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            2044
15331 continue                                                             2045
15332 continue                                                             2045
      sc=min(max(exmn,sc+b0),exmx)                                         2046
      sxp=sxp-q(:,ic)                                                      2047
      q(:,ic)=min(max(emin*sxp,exp(dble(sc))),emax*sxp)                    2048
      sxp=sxp+q(:,ic)                                                      2049
14971 continue                                                             2050
14972 continue                                                             2050
      s=-sum(b(0,:))/nc                                                    2050
      b(0,:)=b(0,:)+s                                                      2050
      sc=s                                                                 2050
      b0=0.0                                                               2051
15340 do 15341 j=1,nin                                                     2051
      l=m(j)                                                               2052
      if(vp(l) .gt. 0.0)goto 15361                                         2052
      s=sum(b(l,:))/nc                                                     2052
      goto 15371                                                           2053
15361 continue                                                             2053
      s=elc(parm,nc,b(l,:),is)                                             2053
15371 continue                                                             2054
15351 continue                                                             2054
      b(l,:)=b(l,:)-s                                                      2055
      jb=ix(l)                                                             2055
      je=ix(l+1)-1                                                         2056
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         2057
      b0=b0+s*xb(l)/xs(l)                                                  2058
15341 continue                                                             2059
15342 continue                                                             2059
      sc=sc+b0                                                             2059
      sc=exp(sc)                                                           2059
      sxp=sxp*sc                                                           2059
15380 do 15381 ic=1,nc                                                     2059
      q(:,ic)=q(:,ic)*sc                                                   2059
15381 continue                                                             2060
15382 continue                                                             2060
      if(jxx.gt.0)goto 14962                                               2060
      if(ig.eq.0)goto 14962                                                2061
      if(ixx .ne. 0)goto 15401                                             2061
      ga=0.0                                                               2062
15410 do 15411 ic=1,nc                                                     2062
      v=q(:,ic)/sxp                                                        2062
      r=w*(y(:,ic)-v)                                                      2062
      v=w*v*(1.0-v)                                                        2063
15420 do 15421 j=1,ni                                                      2063
      if(iy(j).eq.1)goto 15421                                             2063
      if(ju(j).eq.0)goto 15421                                             2064
      jb=ix(j)                                                             2064
      je=ix(j+1)-1                                                         2064
      jn=ix(j+1)-ix(j)                                                     2065
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2066
      gj=dot_product(sc(1:jn),x(jb:je))                                    2067
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2068
15421 continue                                                             2069
15422 continue                                                             2069
15411 continue                                                             2070
15412 continue                                                             2070
15430 do 15431 k=1,ni                                                      2070
      if(iy(k).eq.1)goto 15431                                             2070
      if(ju(k).eq.0)goto 15431                                             2071
      if(ga(k) .le. al1*vp(k))goto 15451                                   2071
      iy(k)=1                                                              2071
      ixx=1                                                                2071
15451 continue                                                             2072
15431 continue                                                             2073
15432 continue                                                             2073
      if(ixx.eq.1) go to 10680                                             2074
      goto 14962                                                           2075
15401 continue                                                             2076
      goto 14961                                                           2077
14962 continue                                                             2077
      if(jxx .le. 0)goto 15471                                             2077
      jerr=-10000-ilm                                                      2077
      goto 14882                                                           2077
15471 continue                                                             2077
      devi=0.0                                                             2078
15480 do 15481 ic=1,nc                                                     2079
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          2079
      a0(ic,ilm)=b(0,ic)                                                   2080
15490 do 15491 i=1,no                                                      2080
      if(y(i,ic).le.0.0)goto 15491                                         2081
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           2082
15491 continue                                                             2083
15492 continue                                                             2083
15481 continue                                                             2084
15482 continue                                                             2084
      kin(ilm)=nin                                                         2084
      alm(ilm)=al                                                          2084
      lmu=ilm                                                              2085
      dev(ilm)=(dev1-devi)/dev0                                            2085
      if(ig.eq.0)goto 14882                                                2086
      if(ilm.lt.mnl)goto 14881                                             2086
      if(flmin.ge.1.0)goto 14881                                           2087
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 14882             2088
      if(dev(ilm).gt.devmax)goto 14882                                     2088
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 14882                             2089
14881 continue                                                             2090
14882 continue                                                             2090
      g=log(q)                                                             2090
15500 do 15501 i=1,no                                                      2090
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2090
15501 continue                                                             2091
15502 continue                                                             2091
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy)                      2092
      return                                                               2093
      end                                                                  2094
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  2095
      real a0(nc),ca(nx,nc),x(*),f(nc,n)                                   2095
      integer ia(*),ix(*),jx(*)                                            2096
15510 do 15511 ic=1,nc                                                     2096
      f(ic,:)=a0(ic)                                                       2096
15511 continue                                                             2097
15512 continue                                                             2097
15520 do 15521 j=1,nin                                                     2097
      k=ia(j)                                                              2097
      kb=ix(k)                                                             2097
      ke=ix(k+1)-1                                                         2098
15530 do 15531 ic=1,nc                                                     2098
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    2098
15531 continue                                                             2099
15532 continue                                                             2099
15521 continue                                                             2100
15522 continue                                                             2100
      return                                                               2101
      end                                                                  2102
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,ne,nx,nlam,flmin,ula   2104 
     *m,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam)              2105
      real ca(nx,nlam),dev(nlam),alm(nlam)                                 2106
      integer jd(*),ia(nx),nin(nlam)                                       2107
      real, dimension (:), allocatable :: xs,ww,vq                              
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 15551                                    2111
      jerr=10000                                                           2111
      return                                                               2111
15551 continue                                                             2112
      allocate(ww(1:no),stat=jerr)                                         2113
      allocate(ju(1:ni),stat=ierr)                                         2113
      jerr=jerr+ierr                                                       2114
      allocate(vq(1:ni),stat=ierr)                                         2114
      jerr=jerr+ierr                                                       2115
      if(isd .le. 0)goto 15571                                             2115
      allocate(xs(1:ni),stat=ierr)                                         2115
      jerr=jerr+ierr                                                       2115
15571 continue                                                             2116
      if(jerr.ne.0) return                                                 2117
      call chkvars(no,ni,x,ju)                                             2118
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2119
      if(maxval(ju) .gt. 0)goto 15591                                      2119
      jerr=7777                                                            2119
      return                                                               2119
15591 continue                                                             2120
      vq=max(0.0,vp)                                                       2120
      vq=vq*ni/sum(vq)                                                     2121
      ww=max(0.0,w)                                                        2121
      sw=sum(ww)                                                           2122
      if(sw .gt. 0.0)goto 15611                                            2122
      jerr=9999                                                            2122
      return                                                               2122
15611 continue                                                             2122
      ww=ww/sw                                                             2123
      call cstandard(no,ni,x,ww,ju,isd,xs)                                 2124
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr   2126 
     *,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 2126
      dev0=2.0*sw*dev0                                                     2127
      if(isd .le. 0)goto 15631                                             2127
15640 do 15641 k=1,lmu                                                     2127
      nk=nin(k)                                                            2127
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   2127
15641 continue                                                             2127
15642 continue                                                             2127
15631 continue                                                             2128
      deallocate(ww,ju,vq)                                                 2128
      if(isd.gt.0) deallocate(xs)                                          2129
      return                                                               2130
      end                                                                  2131
      subroutine cstandard (no,ni,x,w,ju,isd,xs)                           2132
      real x(no,ni),w(no),xs(ni)                                           2132
      integer ju(ni)                                                       2133
15650 do 15651 j=1,ni                                                      2133
      if(ju(j).eq.0)goto 15651                                             2134
      xm=dot_product(w,x(:,j))                                             2134
      x(:,j)=x(:,j)-xm                                                     2135
      if(isd .le. 0)goto 15671                                             2135
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 2135
      x(:,j)=x(:,j)/xs(j)                                                  2135
15671 continue                                                             2136
15651 continue                                                             2137
15652 continue                                                             2137
      return                                                               2138
      end                                                                  2139
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,ne,nx,nlam,flmin,ula   2141 
     *m,cthri,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=0.001, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99)   2142
      real x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam)              2143
      real ao(nx,nlam),dev(nlam),alm(nlam)                                 2144
      integer ju(ni),m(nx),kin(nlam)                                       2145
      real, dimension (:), allocatable :: w,dk,v,xs,wr,a,as,f,dq                
      real, dimension (:), allocatable :: e,uu,ga                               
      integer, dimension (:), allocatable :: jp,kp,mm,ixx                       
      allocate(e(1:no),stat=jerr)                                          2151
      allocate(uu(1:no),stat=ierr)                                         2151
      jerr=jerr+ierr                                                       2152
      allocate(f(1:no),stat=ierr)                                          2152
      jerr=jerr+ierr                                                       2153
      allocate(w(1:no),stat=ierr)                                          2153
      jerr=jerr+ierr                                                       2154
      allocate(v(1:ni),stat=ierr)                                          2154
      jerr=jerr+ierr                                                       2155
      allocate(a(1:ni),stat=ierr)                                          2155
      jerr=jerr+ierr                                                       2156
      allocate(as(1:ni),stat=ierr)                                         2156
      jerr=jerr+ierr                                                       2157
      allocate(xs(1:ni),stat=ierr)                                         2157
      jerr=jerr+ierr                                                       2158
      allocate(ga(1:ni),stat=ierr)                                         2158
      jerr=jerr+ierr                                                       2159
      allocate(ixx(1:ni),stat=ierr)                                        2159
      jerr=jerr+ierr                                                       2160
      allocate(jp(1:no),stat=ierr)                                         2160
      jerr=jerr+ierr                                                       2161
      allocate(kp(1:no),stat=ierr)                                         2161
      jerr=jerr+ierr                                                       2162
      allocate(dk(1:no),stat=ierr)                                         2162
      jerr=jerr+ierr                                                       2163
      allocate(wr(1:no),stat=ierr)                                         2163
      jerr=jerr+ierr                                                       2164
      allocate(dq(1:no),stat=ierr)                                         2164
      jerr=jerr+ierr                                                       2165
      allocate(mm(1:ni),stat=ierr)                                         2165
      jerr=jerr+ierr                                                       2166
      if(jerr.ne.0)go to 11790                                             2167
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2168
      if(jerr.ne.0) go to 11790                                            2168
      alpha=parm                                                           2169
      oma=1.0-alpha                                                        2169
      nlm=0                                                                2169
      ixx=0                                                                2169
      al=0.0                                                               2170
      dq=d*q                                                               2170
      call died(no,nk,dq,kp,jp,dk)                                         2171
      a=0.0                                                                2171
      f(1)=0.0                                                             2171
      fmax=log(huge(f(1))*0.1)                                             2172
      if(nonzero(no,g) .eq. 0)goto 15691                                   2172
      f=g-dot_product(q,g)                                                 2173
      e=q*exp(sign(min(abs(f),fmax),f))                                    2174
      goto 15701                                                           2175
15691 continue                                                             2175
      f=0.0                                                                2175
      e=q                                                                  2175
15701 continue                                                             2176
15681 continue                                                             2176
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 2177
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         2177
      dev0=rr                                                              2178
15710 do 15711 i=1,no                                                      2178
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 15731                   2178
      w(i)=0.0                                                             2178
      wr(i)=w(i)                                                           2178
15731 continue                                                             2178
15711 continue                                                             2179
15712 continue                                                             2179
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2180
      if(jerr.ne.0) go to 11790                                            2181
      if(flmin .ge. 1.0)goto 15751                                         2181
      eqs=max(eps,flmin)                                                   2181
      alf=eqs**(1.0/(nlam-1))                                              2181
15751 continue                                                             2182
      m=0                                                                  2182
      mm=0                                                                 2182
      nlp=0                                                                2182
      nin=nlp                                                              2182
      mnl=min(mnlam,nlam)                                                  2182
      as=0.0                                                               2182
      cthr=cthri*dev0                                                      2183
15760 do 15761 j=1,ni                                                      2183
      if(ju(j).eq.0)goto 15761                                             2183
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2183
15761 continue                                                             2184
15762 continue                                                             2184
15770 do 15771 ilm=1,nlam                                                  2184
      al0=al                                                               2185
      if(flmin .lt. 1.0)goto 15791                                         2185
      al=ulam(ilm)                                                         2185
      goto 15781                                                           2186
15791 if(ilm .le. 2)goto 15801                                             2186
      al=al*alf                                                            2186
      goto 15781                                                           2187
15801 if(ilm .ne. 1)goto 15811                                             2187
      al=big                                                               2187
      goto 15821                                                           2188
15811 continue                                                             2188
      al0=0.0                                                              2189
15830 do 15831 j=1,ni                                                      2189
      if(ju(j).eq.0)goto 15831                                             2189
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2189
15831 continue                                                             2190
15832 continue                                                             2190
      al0=al0/max(parm,1.0e-3)                                             2190
      al=alf*al0                                                           2191
15821 continue                                                             2192
15781 continue                                                             2192
      sa=alpha*al                                                          2192
      omal=oma*al                                                          2192
      tlam=alpha*(2.0*al-al0)                                              2193
15840 do 15841 k=1,ni                                                      2193
      if(ixx(k).eq.1)goto 15841                                            2193
      if(ju(k).eq.0)goto 15841                                             2194
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2195
15841 continue                                                             2196
15842 continue                                                             2196
10680 continue                                                             2197
15850 continue                                                             2197
15851 continue                                                             2197
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2198
      call vars(no,ni,x,w,ixx,v)                                           2199
15860 continue                                                             2199
15861 continue                                                             2199
      nlp=nlp+1                                                            2199
      dli=0.0                                                              2200
15870 do 15871 j=1,ni                                                      2200
      if(ixx(j).eq.0)goto 15871                                            2201
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2202
      if(abs(u) .gt. vp(j)*sa)goto 15891                                   2202
      at=0.0                                                               2202
      goto 15901                                                           2203
15891 continue                                                             2203
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2203
15901 continue                                                             2204
15881 continue                                                             2204
      if(at .eq. a(j))goto 15921                                           2204
      del=at-a(j)                                                          2204
      a(j)=at                                                              2204
      dli=max(dli,v(j)*del**2)                                             2205
      wr=wr-del*w*x(:,j)                                                   2205
      f=f+del*x(:,j)                                                       2206
      if(mm(j) .ne. 0)goto 15941                                           2206
      nin=nin+1                                                            2206
      if(nin.gt.nx)goto 15872                                              2207
      mm(j)=nin                                                            2207
      m(nin)=j                                                             2208
15941 continue                                                             2209
15921 continue                                                             2210
15871 continue                                                             2211
15872 continue                                                             2211
      if(nin.gt.nx)goto 15862                                              2211
      if(dli.lt.cthr)goto 15862                                            2212
      if(nlp .le. maxit)goto 15961                                         2212
      jerr=-ilm                                                            2212
      return                                                               2212
15961 continue                                                             2213
15970 continue                                                             2213
15971 continue                                                             2213
      nlp=nlp+1                                                            2213
      dli=0.0                                                              2214
15980 do 15981 l=1,nin                                                     2214
      j=m(l)                                                               2215
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2216
      if(abs(u) .gt. vp(j)*sa)goto 16001                                   2216
      at=0.0                                                               2216
      goto 16011                                                           2217
16001 continue                                                             2217
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2217
16011 continue                                                             2218
15991 continue                                                             2218
      if(at .eq. a(j))goto 16031                                           2218
      del=at-a(j)                                                          2218
      a(j)=at                                                              2218
      dli=max(dli,v(j)*del**2)                                             2219
      wr=wr-del*w*x(:,j)                                                   2219
      f=f+del*x(:,j)                                                       2220
16031 continue                                                             2221
15981 continue                                                             2222
15982 continue                                                             2222
      if(dli.lt.cthr)goto 15972                                            2222
      if(nlp .le. maxit)goto 16051                                         2222
      jerr=-ilm                                                            2222
      return                                                               2222
16051 continue                                                             2223
      goto 15971                                                           2224
15972 continue                                                             2224
      goto 15861                                                           2225
15862 continue                                                             2225
      if(nin.gt.nx)goto 15852                                              2226
      e=q*exp(sign(min(abs(f),fmax),f))                                    2227
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2228
      if(jerr .eq. 0)goto 16071                                            2228
      jerr=jerr-ilm                                                        2228
      go to 11790                                                          2228
16071 continue                                                             2229
      ix=0                                                                 2230
16080 do 16081 j=1,nin                                                     2230
      k=m(j)                                                               2231
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 16081                           2231
      ix=1                                                                 2231
      goto 16082                                                           2231
16081 continue                                                             2232
16082 continue                                                             2232
      if(ix .ne. 0)goto 16101                                              2233
16110 do 16111 k=1,ni                                                      2233
      if(ixx(k).eq.1)goto 16111                                            2233
      if(ju(k).eq.0)goto 16111                                             2234
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2235
      if(ga(k) .le. sa*vp(k))goto 16131                                    2235
      ixx(k)=1                                                             2235
      ix=1                                                                 2235
16131 continue                                                             2236
16111 continue                                                             2237
16112 continue                                                             2237
      if(ix.eq.1) go to 10680                                              2238
      goto 15852                                                           2239
16101 continue                                                             2240
      goto 15851                                                           2241
15852 continue                                                             2241
      if(nin .le. nx)goto 16151                                            2241
      jerr=-10000-ilm                                                      2241
      goto 15772                                                           2241
16151 continue                                                             2242
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               2242
      kin(ilm)=nin                                                         2243
      alm(ilm)=al                                                          2243
      lmu=ilm                                                              2244
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   2245
      if(ilm.lt.mnl)goto 15771                                             2245
      if(flmin.ge.1.0)goto 15771                                           2246
      me=0                                                                 2246
16160 do 16161 j=1,nin                                                     2246
      if(ao(j,ilm).ne.0.0) me=me+1                                         2246
16161 continue                                                             2246
16162 continue                                                             2246
      if(me.gt.ne)goto 15772                                               2247
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 15772              2248
      if(dev(ilm).gt.devmax)goto 15772                                     2249
15771 continue                                                             2250
15772 continue                                                             2250
      g=f                                                                  2251
11790 continue                                                             2251
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)              2252
      return                                                               2253
      end                                                                  2254
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 2255
      real ca(nin),x(n,*),f(n)                                             2255
      integer ia(nin)                                                      2256
      f=0.0                                                                2256
      if(nin.le.0) return                                                  2257
16170 do 16171 i=1,n                                                       2257
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      2257
16171 continue                                                             2258
16172 continue                                                             2258
      return                                                               2259
      end                                                                  2260
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         2261
      real y(no),d(no),q(no)                                               2261
      integer jp(no),kp(*)                                                 2262
16180 do 16181 j=1,no                                                      2262
      jp(j)=j                                                              2262
16181 continue                                                             2262
16182 continue                                                             2262
      call psort7(y,jp,1,no)                                               2263
      nj=0                                                                 2263
16190 do 16191 j=1,no                                                      2263
      if(q(jp(j)).le.0.0)goto 16191                                        2263
      nj=nj+1                                                              2263
      jp(nj)=jp(j)                                                         2263
16191 continue                                                             2264
16192 continue                                                             2264
      if(nj .ne. 0)goto 16211                                              2264
      jerr=20000                                                           2264
      return                                                               2264
16211 continue                                                             2265
      j=1                                                                  2265
16220 continue                                                             2265
16221 if(d(jp(j)).gt.0.0)goto 16222                                        2265
      j=j+1                                                                2265
      if(j.gt.nj)goto 16222                                                2265
      goto 16221                                                           2266
16222 continue                                                             2266
      if(j .lt. nj-1)goto 16241                                            2266
      jerr=30000                                                           2266
      return                                                               2266
16241 continue                                                             2267
      j0=j-1                                                               2267
      nj=nj-j0                                                             2267
16250 do 16251 j=1,nj                                                      2267
      jp(j)=jp(j+j0)                                                       2267
16251 continue                                                             2268
16252 continue                                                             2268
      jerr=0                                                               2268
      nk=0                                                                 2268
      t0=y(jp(1))                                                          2268
      yk=t0                                                                2268
      j=2                                                                  2269
16260 continue                                                             2269
16261 continue                                                             2269
16270 continue                                                             2270
16271 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 16272                     2270
      j=j+1                                                                2270
      if(j.gt.nj)goto 16272                                                2270
      goto 16271                                                           2271
16272 continue                                                             2271
      nk=nk+1                                                              2271
      kp(nk)=j-1                                                           2271
      if(j.gt.nj)goto 16262                                                2272
      if(j .ne. nj)goto 16291                                              2272
      nk=nk+1                                                              2272
      kp(nk)=nj                                                            2272
      goto 16262                                                           2272
16291 continue                                                             2273
      yk=y(jp(j))                                                          2273
      j=j+1                                                                2274
      goto 16261                                                           2275
16262 continue                                                             2275
      return                                                               2276
      end                                                                  2277
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     2278
      real d(no),dk(nk),wr(no),w(no)                                       2279
      real e(no),u(no),b,c                                                 2279
      integer kp(nk),jp(no)                                                2280
      call usk(no,nk,kp,jp,e,u)                                            2281
      b=dk(1)/u(1)                                                         2281
      c=dk(1)/u(1)**2                                                      2281
      jerr=0                                                               2282
16300 do 16301 j=1,kp(1)                                                   2282
      i=jp(j)                                                              2283
      w(i)=e(i)*(b-e(i)*c)                                                 2283
      if(w(i) .gt. 0.0)goto 16321                                          2283
      jerr=-3                                                              2283
      return                                                               2283
16321 continue                                                             2284
      wr(i)=d(i)-e(i)*b                                                    2285
16301 continue                                                             2286
16302 continue                                                             2286
16330 do 16331 k=2,nk                                                      2286
      j1=kp(k-1)+1                                                         2286
      j2=kp(k)                                                             2287
      b=b+dk(k)/u(k)                                                       2287
      c=c+dk(k)/u(k)**2                                                    2288
16340 do 16341 j=j1,j2                                                     2288
      i=jp(j)                                                              2289
      w(i)=e(i)*(b-e(i)*c)                                                 2289
      if(w(i) .gt. 0.0)goto 16361                                          2289
      jerr=-30000                                                          2289
      return                                                               2289
16361 continue                                                             2290
      wr(i)=d(i)-e(i)*b                                                    2291
16341 continue                                                             2292
16342 continue                                                             2292
16331 continue                                                             2293
16332 continue                                                             2293
      return                                                               2294
      end                                                                  2295
      subroutine vars(no,ni,x,w,ixx,v)                                     2296
      real x(no,ni),w(no),v(ni)                                            2296
      integer ixx(ni)                                                      2297
16370 do 16371 j=1,ni                                                      2297
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                        2297
16371 continue                                                             2298
16372 continue                                                             2298
      return                                                               2299
      end                                                                  2300
      subroutine died(no,nk,d,kp,jp,dk)                                    2301
      real d(no),dk(nk)                                                    2301
      integer kp(nk),jp(no)                                                2302
      dk(1)=sum(d(jp(1:kp(1))))                                            2303
16380 do 16381 k=2,nk                                                      2303
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  2303
16381 continue                                                             2304
16382 continue                                                             2304
      return                                                               2305
      end                                                                  2306
      subroutine usk(no,nk,kp,jp,e,u)                                      2307
      real e(no),u(nk),h                                                   2307
      integer kp(nk),jp(no)                                                2308
      h=0.0                                                                2309
16390 do 16391 k=nk,1,-1                                                   2309
      j2=kp(k)                                                             2310
      j1=1                                                                 2310
      if(k.gt.1) j1=kp(k-1)+1                                              2311
16400 do 16401 j=j2,j1,-1                                                  2311
      h=h+e(jp(j))                                                         2311
16401 continue                                                             2312
16402 continue                                                             2312
      u(k)=h                                                               2313
16391 continue                                                             2314
16392 continue                                                             2314
      return                                                               2315
      end                                                                  2316
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             2317
      real d(no),dk(nk),f(no)                                              2318
      integer kp(nk),jp(no)                                                2318
      real e(no),u(nk),s                                                   2319
      call usk(no,nk,kp,jp,e,u)                                            2319
      u=log(u)                                                             2320
      risk=dot_product(d,f)-dot_product(dk,u)                              2321
      return                                                               2322
      end                                                                  2323
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 2324
      real x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(nlam)          2325
      real, dimension (:), allocatable :: dk,f,xm,dq,q                          
      real, dimension (:), allocatable :: e,uu                                  
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          2331
      allocate(q(1:no),stat=ierr)                                          2331
      jerr=jerr+ierr                                                       2332
      allocate(uu(1:no),stat=ierr)                                         2332
      jerr=jerr+ierr                                                       2333
      allocate(f(1:no),stat=ierr)                                          2333
      jerr=jerr+ierr                                                       2334
      allocate(dk(1:no),stat=ierr)                                         2334
      jerr=jerr+ierr                                                       2335
      allocate(jp(1:no),stat=ierr)                                         2335
      jerr=jerr+ierr                                                       2336
      allocate(kp(1:no),stat=ierr)                                         2336
      jerr=jerr+ierr                                                       2337
      allocate(dq(1:no),stat=ierr)                                         2337
      jerr=jerr+ierr                                                       2338
      allocate(xm(1:ni),stat=ierr)                                         2338
      jerr=jerr+ierr                                                       2339
      if(jerr.ne.0) go to 11790                                            2340
      q=max(0.0,w)                                                         2340
      sw=sum(q)                                                            2341
      if(sw .gt. 0.0)goto 16421                                            2341
      jerr=9999                                                            2341
      go to 11790                                                          2341
16421 continue                                                             2342
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2343
      if(jerr.ne.0) go to 11790                                            2343
      fmax=log(huge(e(1))*0.1)                                             2344
      dq=d*q                                                               2344
      call died(no,nk,dq,kp,jp,dk)                                         2344
      gm=dot_product(q,g)/sw                                               2345
16430 do 16431 j=1,ni                                                      2345
      xm(j)=dot_product(q,x(:,j))/sw                                       2345
16431 continue                                                             2346
16432 continue                                                             2346
16440 do 16441 lam=1,nlam                                                  2347
16450 do 16451 i=1,no                                                      2347
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       2348
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        2349
16451 continue                                                             2350
16452 continue                                                             2350
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          2351
16441 continue                                                             2352
16442 continue                                                             2352
11790 continue                                                             2352
      deallocate(e,uu,dk,f,jp,kp,dq)                                       2353
      return                                                               2354
      end                                                                  2355
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,ne,nx,nlam,flmin,ulam   2357 
     *,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    2358
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2359
      integer jd(*),ia(nx),nin(nlam)                                       2360
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16471                                    2364
      jerr=10000                                                           2364
      return                                                               2364
16471 continue                                                             2365
      if(minval(y) .ge. 0.0)goto 16491                                     2365
      jerr=8888                                                            2365
      return                                                               2365
16491 continue                                                             2366
      allocate(ww(1:no),stat=jerr)                                         2367
      allocate(ju(1:ni),stat=ierr)                                         2367
      jerr=jerr+ierr                                                       2368
      allocate(vq(1:ni),stat=ierr)                                         2368
      jerr=jerr+ierr                                                       2369
      allocate(xm(1:ni),stat=ierr)                                         2369
      jerr=jerr+ierr                                                       2370
      if(isd .le. 0)goto 16511                                             2370
      allocate(xs(1:ni),stat=ierr)                                         2370
      jerr=jerr+ierr                                                       2370
16511 continue                                                             2371
      if(jerr.ne.0) return                                                 2372
      call chkvars(no,ni,x,ju)                                             2373
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2374
      if(maxval(ju) .gt. 0)goto 16531                                      2374
      jerr=7777                                                            2374
      go to 11790                                                          2374
16531 continue                                                             2375
      vq=max(0.0,vp)                                                       2375
      vq=vq*ni/sum(vq)                                                     2376
      ww=max(0.0,w)                                                        2376
      sw=sum(ww)                                                           2376
      if(sw .gt. 0.0)goto 16551                                            2376
      jerr=9999                                                            2376
      go to 11790                                                          2376
16551 continue                                                             2377
      ww=ww/sw                                                             2378
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             2379
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr,   2381 
     *  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11790                                            2381
      dev0=2.0*sw*dev0                                                     2382
16560 do 16561 k=1,lmu                                                     2382
      nk=nin(k)                                                            2383
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2384
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2385
16561 continue                                                             2386
16562 continue                                                             2386
11790 continue                                                             2386
      deallocate(ww,ju,vq,xm)                                              2386
      if(isd.gt.0) deallocate(xs)                                          2387
      return                                                               2388
      end                                                                  2389
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,ne,nx,nlam,flmin,ulam   2391 
     *,shri,  isd,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2392 
     *9)
      real x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)                    2393
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2394
      integer ju(ni),m(nx),kin(nlam)                                       2395
      real, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga                    
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(a(1:ni),stat=jerr)                                          2400
      allocate(as(1:ni),stat=ierr)                                         2400
      jerr=jerr+ierr                                                       2401
      allocate(t(1:no),stat=ierr)                                          2401
      jerr=jerr+ierr                                                       2402
      allocate(mm(1:ni),stat=ierr)                                         2402
      jerr=jerr+ierr                                                       2403
      allocate(ga(1:ni),stat=ierr)                                         2403
      jerr=jerr+ierr                                                       2404
      allocate(ixx(1:ni),stat=ierr)                                        2404
      jerr=jerr+ierr                                                       2405
      allocate(wr(1:no),stat=ierr)                                         2405
      jerr=jerr+ierr                                                       2406
      allocate(v(1:ni),stat=ierr)                                          2406
      jerr=jerr+ierr                                                       2407
      allocate(w(1:no),stat=ierr)                                          2407
      jerr=jerr+ierr                                                       2408
      allocate(f(1:no),stat=ierr)                                          2408
      jerr=jerr+ierr                                                       2409
      if(jerr.ne.0) return                                                 2410
      bta=parm                                                             2410
      omb=1.0-bta                                                          2411
      t=q*y                                                                2411
      yb=sum(t)                                                            2411
      fmax=log(huge(bta)*0.1)                                              2412
      if(nonzero(no,g) .ne. 0)goto 16581                                   2412
      w=q*yb                                                               2412
      az=log(yb)                                                           2412
      f=az                                                                 2412
      dv0=yb*(log(yb)-1.0)                                                 2412
      goto 16591                                                           2413
16581 continue                                                             2413
      w=q*exp(sign(min(abs(g),fmax),g))                                    2413
      v0=sum(w)                                                            2413
      eaz=yb/v0                                                            2414
      w=eaz*w                                                              2414
      az=log(eaz)                                                          2414
      f=az+g                                                               2415
      dv0=dot_product(t,g)-yb*(1.0-az)                                     2416
16591 continue                                                             2417
16571 continue                                                             2417
      a=0.0                                                                2417
      as=0.0                                                               2417
      wr=t-w                                                               2417
      v0=yb                                                                2417
      dvr=-yb                                                              2418
16600 do 16601 i=1,no                                                      2418
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               2418
16601 continue                                                             2418
16602 continue                                                             2418
      dvr=dvr-dv0                                                          2418
      dev0=dvr                                                             2419
      if(flmin .ge. 1.0)goto 16621                                         2419
      eqs=max(eps,flmin)                                                   2419
      alf=eqs**(1.0/(nlam-1))                                              2419
16621 continue                                                             2420
      m=0                                                                  2420
      mm=0                                                                 2420
      nlp=0                                                                2420
      nin=nlp                                                              2420
      mnl=min(mnlam,nlam)                                                  2420
      shr=shri*dev0                                                        2420
      ixx=0                                                                2420
      al=0.0                                                               2421
16630 do 16631 j=1,ni                                                      2421
      if(ju(j).eq.0)goto 16631                                             2421
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2421
16631 continue                                                             2422
16632 continue                                                             2422
16640 do 16641 ilm=1,nlam                                                  2422
      al0=al                                                               2423
      if(flmin .lt. 1.0)goto 16661                                         2423
      al=ulam(ilm)                                                         2423
      goto 16651                                                           2424
16661 if(ilm .le. 2)goto 16671                                             2424
      al=al*alf                                                            2424
      goto 16651                                                           2425
16671 if(ilm .ne. 1)goto 16681                                             2425
      al=big                                                               2425
      goto 16691                                                           2426
16681 continue                                                             2426
      al0=0.0                                                              2427
16700 do 16701 j=1,ni                                                      2427
      if(ju(j).eq.0)goto 16701                                             2427
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2427
16701 continue                                                             2428
16702 continue                                                             2428
      al0=al0/max(bta,1.0e-3)                                              2428
      al=alf*al0                                                           2429
16691 continue                                                             2430
16651 continue                                                             2430
      al2=al*omb                                                           2430
      al1=al*bta                                                           2430
      tlam=bta*(2.0*al-al0)                                                2431
16710 do 16711 k=1,ni                                                      2431
      if(ixx(k).eq.1)goto 16711                                            2431
      if(ju(k).eq.0)goto 16711                                             2432
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2433
16711 continue                                                             2434
16712 continue                                                             2434
10680 continue                                                             2435
16720 continue                                                             2435
16721 continue                                                             2435
      az0=az                                                               2436
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2437
16730 do 16731 j=1,ni                                                      2437
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                        2437
16731 continue                                                             2438
16732 continue                                                             2438
16740 continue                                                             2438
16741 continue                                                             2438
      nlp=nlp+1                                                            2438
      dlx=0.0                                                              2439
16750 do 16751 k=1,ni                                                      2439
      if(ixx(k).eq.0)goto 16751                                            2439
      ak=a(k)                                                              2440
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2440
      au=abs(u)-vp(k)*al1                                                  2441
      if(au .gt. 0.0)goto 16771                                            2441
      a(k)=0.0                                                             2441
      goto 16781                                                           2442
16771 continue                                                             2442
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2442
16781 continue                                                             2443
16761 continue                                                             2443
      if(a(k).eq.ak)goto 16751                                             2443
      d=a(k)-ak                                                            2443
      dlx=max(dlx,v(k)*d**2)                                               2444
      wr=wr-d*w*x(:,k)                                                     2444
      f=f+d*x(:,k)                                                         2445
      if(mm(k) .ne. 0)goto 16801                                           2445
      nin=nin+1                                                            2445
      if(nin.gt.nx)goto 16752                                              2446
      mm(k)=nin                                                            2446
      m(nin)=k                                                             2447
16801 continue                                                             2448
16751 continue                                                             2449
16752 continue                                                             2449
      if(nin.gt.nx)goto 16742                                              2449
      d=sum(wr)/v0                                                         2450
      az=az+d                                                              2450
      dlx=max(dlx,v0*d**2)                                                 2450
      wr=wr-d*w                                                            2450
      f=f+d                                                                2451
      if(dlx.lt.shr)goto 16742                                             2451
      if(nlp .le. maxit)goto 16821                                         2451
      jerr=-ilm                                                            2451
      return                                                               2451
16821 continue                                                             2452
16830 continue                                                             2452
16831 continue                                                             2452
      nlp=nlp+1                                                            2452
      dlx=0.0                                                              2453
16840 do 16841 l=1,nin                                                     2453
      k=m(l)                                                               2453
      ak=a(k)                                                              2454
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2454
      au=abs(u)-vp(k)*al1                                                  2455
      if(au .gt. 0.0)goto 16861                                            2455
      a(k)=0.0                                                             2455
      goto 16871                                                           2456
16861 continue                                                             2456
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2456
16871 continue                                                             2457
16851 continue                                                             2457
      if(a(k).eq.ak)goto 16841                                             2457
      d=a(k)-ak                                                            2457
      dlx=max(dlx,v(k)*d**2)                                               2458
      wr=wr-d*w*x(:,k)                                                     2458
      f=f+d*x(:,k)                                                         2460
16841 continue                                                             2460
16842 continue                                                             2460
      d=sum(wr)/v0                                                         2460
      az=az+d                                                              2460
      dlx=max(dlx,v0*d**2)                                                 2460
      wr=wr-d*w                                                            2460
      f=f+d                                                                2461
      if(dlx.lt.shr)goto 16832                                             2461
      if(nlp .le. maxit)goto 16891                                         2461
      jerr=-ilm                                                            2461
      return                                                               2461
16891 continue                                                             2462
      goto 16831                                                           2463
16832 continue                                                             2463
      goto 16741                                                           2464
16742 continue                                                             2464
      if(nin.gt.nx)goto 16722                                              2465
      w=q*exp(sign(min(abs(f),fmax),f))                                    2465
      v0=sum(w)                                                            2465
      wr=t-w                                                               2466
      if(v0*(az-az0)**2 .ge. shr)goto 16911                                2466
      ix=0                                                                 2467
16920 do 16921 j=1,nin                                                     2467
      k=m(j)                                                               2468
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 16921                            2468
      ix=1                                                                 2468
      goto 16922                                                           2469
16921 continue                                                             2470
16922 continue                                                             2470
      if(ix .ne. 0)goto 16941                                              2471
16950 do 16951 k=1,ni                                                      2471
      if(ixx(k).eq.1)goto 16951                                            2471
      if(ju(k).eq.0)goto 16951                                             2472
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2473
      if(ga(k) .le. al1*vp(k))goto 16971                                   2473
      ixx(k)=1                                                             2473
      ix=1                                                                 2473
16971 continue                                                             2474
16951 continue                                                             2475
16952 continue                                                             2475
      if(ix.eq.1) go to 10680                                              2476
      goto 16722                                                           2477
16941 continue                                                             2478
16911 continue                                                             2479
      goto 16721                                                           2480
16722 continue                                                             2480
      if(nin .le. nx)goto 16991                                            2480
      jerr=-10000-ilm                                                      2480
      goto 16642                                                           2480
16991 continue                                                             2481
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2481
      kin(ilm)=nin                                                         2482
      a0(ilm)=az                                                           2482
      alm(ilm)=al                                                          2482
      lmu=ilm                                                              2483
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               2484
      if(ilm.lt.mnl)goto 16641                                             2484
      if(flmin.ge.1.0)goto 16641                                           2485
      me=0                                                                 2485
17000 do 17001 j=1,nin                                                     2485
      if(ca(j,ilm).ne.0.0) me=me+1                                         2485
17001 continue                                                             2485
17002 continue                                                             2485
      if(me.gt.ne)goto 16642                                               2486
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 16642              2487
      if(dev(ilm).gt.devmax)goto 16642                                     2488
16641 continue                                                             2489
16642 continue                                                             2489
      g=f                                                                  2490
11790 continue                                                             2490
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                                2491
      return                                                               2492
      end                                                                  2493
      function nonzero(n,v)                                                2494
      real v(n)                                                            2495
      nonzero=0                                                            2495
17010 do 17011 i=1,n                                                       2495
      if(v(i) .eq. 0.0)goto 17031                                          2495
      nonzero=1                                                            2495
      return                                                               2495
17031 continue                                                             2495
17011 continue                                                             2496
17012 continue                                                             2496
      return                                                               2497
      end                                                                  2498
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               2499
      real a(nx,lmu),b(ni,lmu)                                             2499
      integer ia(nx),nin(lmu)                                              2500
17040 do 17041 lam=1,lmu                                                   2500
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        2500
17041 continue                                                             2501
17042 continue                                                             2501
      return                                                               2502
      end                                                                  2503
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           2504
      real a(nx,nc,lmu),b(ni,nc,lmu)                                       2504
      integer ia(nx),nin(lmu)                                              2505
17050 do 17051 lam=1,lmu                                                   2505
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             2505
17051 continue                                                             2506
17052 continue                                                             2506
      return                                                               2507
      end                                                                  2508
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               2509
      real x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)       2510
      real, dimension (:), allocatable :: w                                     
      if(minval(y) .ge. 0.0)goto 17071                                     2513
      jerr=8888                                                            2513
      return                                                               2513
17071 continue                                                             2514
      allocate(w(1:no),stat=jerr)                                          2514
      if(jerr.ne.0) return                                                 2515
      w=max(0.0,q)                                                         2515
      sw=sum(w)                                                            2515
      if(sw .gt. 0.0)goto 17091                                            2515
      jerr=9999                                                            2515
      go to 11790                                                          2515
17091 continue                                                             2516
      yb=dot_product(w,y)/sw                                               2516
      fmax=log(huge(y(1))*0.1)                                             2517
17100 do 17101 lam=1,nlam                                                  2517
      s=0.0                                                                2518
17110 do 17111 i=1,no                                                      2518
      if(w(i).le.0.0)goto 17111                                            2519
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          2520
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      2521
17111 continue                                                             2522
17112 continue                                                             2522
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2523
17101 continue                                                             2524
17102 continue                                                             2524
11790 continue                                                             2524
      deallocate(w)                                                        2525
      return                                                               2526
      end                                                                  2527
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,ne,nx,nlam,fl   2529 
     *min,ulam,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        2530
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2531
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           2532
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 17131                                    2536
      jerr=10000                                                           2536
      return                                                               2536
17131 continue                                                             2537
      if(minval(y) .ge. 0.0)goto 17151                                     2537
      jerr=8888                                                            2537
      return                                                               2537
17151 continue                                                             2538
      allocate(ww(1:no),stat=jerr)                                         2539
      allocate(ju(1:ni),stat=ierr)                                         2539
      jerr=jerr+ierr                                                       2540
      allocate(vq(1:ni),stat=ierr)                                         2540
      jerr=jerr+ierr                                                       2541
      allocate(xm(1:ni),stat=ierr)                                         2541
      jerr=jerr+ierr                                                       2542
      allocate(xs(1:ni),stat=ierr)                                         2542
      jerr=jerr+ierr                                                       2543
      if(jerr.ne.0) return                                                 2544
      call spchkvars(no,ni,x,ix,ju)                                        2545
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2546
      if(maxval(ju) .gt. 0)goto 17171                                      2546
      jerr=7777                                                            2546
      go to 11790                                                          2546
17171 continue                                                             2547
      vq=max(0.0,vp)                                                       2547
      vq=vq*ni/sum(vq)                                                     2548
      ww=max(0.0,w)                                                        2548
      sw=sum(ww)                                                           2548
      if(sw .gt. 0.0)goto 17191                                            2548
      jerr=9999                                                            2548
      go to 11790                                                          2548
17191 continue                                                             2549
      ww=ww/sw                                                             2550
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     2551
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmin,u   2553 
     *lam,thr,  isd,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11790                                            2553
      dev0=2.0*sw*dev0                                                     2554
17200 do 17201 k=1,lmu                                                     2554
      nk=nin(k)                                                            2555
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2556
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2557
17201 continue                                                             2558
17202 continue                                                             2558
11790 continue                                                             2558
      deallocate(ww,ju,vq,xm,xs)                                           2559
      return                                                               2560
      end                                                                  2561
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,ne,nx,nlam,fl   2563 
     *min,ulam,  shri,isd,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,j
     *err)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2564 
     *9)
      real x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),xs(ni)          2565
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2566
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2567
      real, dimension (:), allocatable :: qy,t,w,wr,v,a,as,xm,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      allocate(a(1:ni),stat=jerr)                                          2572
      allocate(as(1:ni),stat=ierr)                                         2572
      jerr=jerr+ierr                                                       2573
      allocate(t(1:no),stat=ierr)                                          2573
      jerr=jerr+ierr                                                       2574
      allocate(mm(1:ni),stat=ierr)                                         2574
      jerr=jerr+ierr                                                       2575
      allocate(ga(1:ni),stat=ierr)                                         2575
      jerr=jerr+ierr                                                       2576
      allocate(ixx(1:ni),stat=ierr)                                        2576
      jerr=jerr+ierr                                                       2577
      allocate(wr(1:no),stat=ierr)                                         2577
      jerr=jerr+ierr                                                       2578
      allocate(v(1:ni),stat=ierr)                                          2578
      jerr=jerr+ierr                                                       2579
      allocate(xm(1:ni),stat=ierr)                                         2579
      jerr=jerr+ierr                                                       2580
      allocate(w(1:no),stat=ierr)                                          2580
      jerr=jerr+ierr                                                       2581
      allocate(qy(1:no),stat=ierr)                                         2581
      jerr=jerr+ierr                                                       2582
      if(jerr.ne.0) return                                                 2583
      bta=parm                                                             2583
      omb=1.0-bta                                                          2583
      fmax=log(huge(bta)*0.1)                                              2584
      qy=q*y                                                               2584
      yb=sum(qy)                                                           2585
      if(nonzero(no,g) .ne. 0)goto 17221                                   2585
      w=q*yb                                                               2585
      az=log(yb)                                                           2585
      uu=az                                                                2586
      xm=yb*xb                                                             2586
      t=0.0                                                                2586
      dv0=yb*(log(yb)-1.0)                                                 2587
      goto 17231                                                           2588
17221 continue                                                             2588
      w=q*exp(sign(min(abs(g),fmax),g))                                    2588
      ww=sum(w)                                                            2588
      eaz=yb/ww                                                            2589
      w=eaz*w                                                              2589
      az=log(eaz)                                                          2589
      uu=az                                                                2589
      t=g                                                                  2589
      dv0=dot_product(qy,g)-yb*(1.0-az)                                    2590
17240 do 17241 j=1,ni                                                      2590
      if(ju(j).eq.0)goto 17241                                             2590
      jb=ix(j)                                                             2590
      je=ix(j+1)-1                                                         2591
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2592
17241 continue                                                             2593
17242 continue                                                             2593
17231 continue                                                             2594
17211 continue                                                             2594
      tt=yb*uu                                                             2594
      ww=yb                                                                2594
      wr=qy-q*(yb*(1.0-uu))                                                2594
      a=0.0                                                                2594
      as=0.0                                                               2595
      dvr=-yb                                                              2596
17250 do 17251 i=1,no                                                      2596
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             2596
17251 continue                                                             2596
17252 continue                                                             2596
      dvr=dvr-dv0                                                          2596
      dev0=dvr                                                             2597
      if(flmin .ge. 1.0)goto 17271                                         2597
      eqs=max(eps,flmin)                                                   2597
      alf=eqs**(1.0/(nlam-1))                                              2597
17271 continue                                                             2598
      m=0                                                                  2598
      mm=0                                                                 2598
      nlp=0                                                                2598
      nin=nlp                                                              2598
      mnl=min(mnlam,nlam)                                                  2598
      shr=shri*dev0                                                        2598
      al=0.0                                                               2598
      ixx=0                                                                2599
17280 do 17281 j=1,ni                                                      2599
      if(ju(j).eq.0)goto 17281                                             2600
      jb=ix(j)                                                             2600
      je=ix(j+1)-1                                                         2601
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   2603 
     *)-xb(j)*tt)/xs(j)
17281 continue                                                             2604
17282 continue                                                             2604
17290 do 17291 ilm=1,nlam                                                  2604
      al0=al                                                               2605
      if(flmin .lt. 1.0)goto 17311                                         2605
      al=ulam(ilm)                                                         2605
      goto 17301                                                           2606
17311 if(ilm .le. 2)goto 17321                                             2606
      al=al*alf                                                            2606
      goto 17301                                                           2607
17321 if(ilm .ne. 1)goto 17331                                             2607
      al=big                                                               2607
      goto 17341                                                           2608
17331 continue                                                             2608
      al0=0.0                                                              2609
17350 do 17351 j=1,ni                                                      2609
      if(ju(j).eq.0)goto 17351                                             2609
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2609
17351 continue                                                             2610
17352 continue                                                             2610
      al0=al0/max(bta,1.0e-3)                                              2610
      al=alf*al0                                                           2611
17341 continue                                                             2612
17301 continue                                                             2612
      al2=al*omb                                                           2612
      al1=al*bta                                                           2612
      tlam=bta*(2.0*al-al0)                                                2613
17360 do 17361 k=1,ni                                                      2613
      if(ixx(k).eq.1)goto 17361                                            2613
      if(ju(k).eq.0)goto 17361                                             2614
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2615
17361 continue                                                             2616
17362 continue                                                             2616
10680 continue                                                             2617
17370 continue                                                             2617
17371 continue                                                             2617
      az0=az                                                               2618
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2619
17380 do 17381 j=1,ni                                                      2619
      if(ixx(j).eq.0)goto 17381                                            2619
      jb=ix(j)                                                             2619
      je=ix(j+1)-1                                                         2620
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2621
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   2623 
     *b(j)**2)/xs(j)**2
17381 continue                                                             2624
17382 continue                                                             2624
17390 continue                                                             2624
17391 continue                                                             2624
      nlp=nlp+1                                                            2625
      dlx=0.0                                                              2626
17400 do 17401 k=1,ni                                                      2626
      if(ixx(k).eq.0)goto 17401                                            2626
      jb=ix(k)                                                             2626
      je=ix(k+1)-1                                                         2626
      ak=a(k)                                                              2627
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2629 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2630
      if(au .gt. 0.0)goto 17421                                            2630
      a(k)=0.0                                                             2630
      goto 17431                                                           2631
17421 continue                                                             2631
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2631
17431 continue                                                             2632
17411 continue                                                             2632
      if(a(k).eq.ak)goto 17401                                             2633
      if(mm(k) .ne. 0)goto 17451                                           2633
      nin=nin+1                                                            2633
      if(nin.gt.nx)goto 17402                                              2634
      mm(k)=nin                                                            2634
      m(nin)=k                                                             2635
17451 continue                                                             2636
      d=a(k)-ak                                                            2636
      dlx=max(dlx,v(k)*d**2)                                               2636
      dv=d/xs(k)                                                           2637
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2638
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2639
      uu=uu-dv*xb(k)                                                       2639
      tt=tt-dv*xm(k)                                                       2640
17401 continue                                                             2641
17402 continue                                                             2641
      if(nin.gt.nx)goto 17392                                              2641
      d=tt/ww-uu                                                           2642
      az=az+d                                                              2642
      dlx=max(dlx,ww*d**2)                                                 2642
      uu=uu+d                                                              2643
      if(dlx.lt.shr)goto 17392                                             2643
      if(nlp .le. maxit)goto 17471                                         2643
      jerr=-ilm                                                            2643
      return                                                               2643
17471 continue                                                             2644
17480 continue                                                             2644
17481 continue                                                             2644
      nlp=nlp+1                                                            2644
      dlx=0.0                                                              2645
17490 do 17491 l=1,nin                                                     2645
      k=m(l)                                                               2646
      jb=ix(k)                                                             2646
      je=ix(k+1)-1                                                         2646
      ak=a(k)                                                              2647
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2649 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2650
      if(au .gt. 0.0)goto 17511                                            2650
      a(k)=0.0                                                             2650
      goto 17521                                                           2651
17511 continue                                                             2651
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2651
17521 continue                                                             2652
17501 continue                                                             2652
      if(a(k).eq.ak)goto 17491                                             2652
      d=a(k)-ak                                                            2652
      dlx=max(dlx,v(k)*d**2)                                               2653
      dv=d/xs(k)                                                           2653
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2654
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2655
      uu=uu-dv*xb(k)                                                       2655
      tt=tt-dv*xm(k)                                                       2656
17491 continue                                                             2657
17492 continue                                                             2657
      d=tt/ww-uu                                                           2657
      az=az+d                                                              2657
      dlx=max(dlx,ww*d**2)                                                 2657
      uu=uu+d                                                              2658
      if(dlx.lt.shr)goto 17482                                             2658
      if(nlp .le. maxit)goto 17541                                         2658
      jerr=-ilm                                                            2658
      return                                                               2658
17541 continue                                                             2659
      goto 17481                                                           2660
17482 continue                                                             2660
      goto 17391                                                           2661
17392 continue                                                             2661
      if(nin.gt.nx)goto 17372                                              2662
      euu=exp(sign(min(abs(uu),fmax),uu))                                  2663
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                2663
      ww=sum(w)                                                            2664
      wr=qy-w*(1.0-uu)                                                     2664
      tt=sum(wr)                                                           2665
      if(ww*(az-az0)**2 .ge. shr)goto 17561                                2665
      kx=0                                                                 2666
17570 do 17571 j=1,nin                                                     2666
      k=m(j)                                                               2667
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 17571                            2667
      kx=1                                                                 2667
      goto 17572                                                           2668
17571 continue                                                             2669
17572 continue                                                             2669
      if(kx .ne. 0)goto 17591                                              2670
17600 do 17601 j=1,ni                                                      2670
      if(ixx(j).eq.1)goto 17601                                            2670
      if(ju(j).eq.0)goto 17601                                             2671
      jb=ix(j)                                                             2671
      je=ix(j+1)-1                                                         2672
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2673
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   2675 
     *)-xb(j)*tt)/xs(j)
      if(ga(j) .le. al1*vp(j))goto 17621                                   2675
      ixx(j)=1                                                             2675
      kx=1                                                                 2675
17621 continue                                                             2676
17601 continue                                                             2677
17602 continue                                                             2677
      if(kx.eq.1) go to 10680                                              2678
      goto 17372                                                           2679
17591 continue                                                             2680
17561 continue                                                             2681
      goto 17371                                                           2682
17372 continue                                                             2682
      if(nin .le. nx)goto 17641                                            2682
      jerr=-10000-ilm                                                      2682
      goto 17292                                                           2682
17641 continue                                                             2683
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2683
      kin(ilm)=nin                                                         2684
      a0(ilm)=az                                                           2684
      alm(ilm)=al                                                          2684
      lmu=ilm                                                              2685
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        2686
      if(ilm.lt.mnl)goto 17291                                             2686
      if(flmin.ge.1.0)goto 17291                                           2687
      me=0                                                                 2687
17650 do 17651 j=1,nin                                                     2687
      if(ca(j,ilm).ne.0.0) me=me+1                                         2687
17651 continue                                                             2687
17652 continue                                                             2687
      if(me.gt.ne)goto 17292                                               2688
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17292              2689
      if(dev(ilm).gt.devmax)goto 17292                                     2690
17291 continue                                                             2691
17292 continue                                                             2691
      g=t+uu                                                               2692
11790 continue                                                             2692
      deallocate(t,w,wr,v,a,qy,xm,as,mm,ga,ixx)                            2693
      return                                                               2694
      end                                                                  2695
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       2696
      real x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)           2697
      integer ix(*),jx(*)                                                  2698
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 17671                                     2701
      jerr=8888                                                            2701
      return                                                               2701
17671 continue                                                             2702
      allocate(w(1:no),stat=jerr)                                          2703
      allocate(f(1:no),stat=ierr)                                          2703
      jerr=jerr+ierr                                                       2704
      if(jerr.ne.0) return                                                 2705
      w=max(0.0,q)                                                         2705
      sw=sum(w)                                                            2705
      if(sw .gt. 0.0)goto 17691                                            2705
      jerr=9999                                                            2705
      go to 11790                                                          2705
17691 continue                                                             2706
      yb=dot_product(w,y)/sw                                               2706
      fmax=log(huge(y(1))*0.1)                                             2707
17700 do 17701 lam=1,nlam                                                  2707
      f=a0(lam)                                                            2708
17710 do 17711 j=1,ni                                                      2708
      if(a(j,lam).eq.0.0)goto 17711                                        2708
      jb=ix(j)                                                             2708
      je=ix(j+1)-1                                                         2709
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          2710
17711 continue                                                             2711
17712 continue                                                             2711
      f=f+g                                                                2712
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2713
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2714
17701 continue                                                             2715
17702 continue                                                             2715
11790 continue                                                             2715
      deallocate(w,f)                                                      2716
      return                                                               2717
      end                                                                  2718
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   2719 
     *jerr)
      real x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(nlam)          2720
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 2721
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 17731                                     2724
      jerr=8888                                                            2724
      return                                                               2724
17731 continue                                                             2725
      allocate(w(1:no),stat=jerr)                                          2726
      allocate(f(1:no),stat=ierr)                                          2726
      jerr=jerr+ierr                                                       2727
      if(jerr.ne.0) return                                                 2728
      w=max(0.0,q)                                                         2728
      sw=sum(w)                                                            2728
      if(sw .gt. 0.0)goto 17751                                            2728
      jerr=9999                                                            2728
      go to 11790                                                          2728
17751 continue                                                             2729
      yb=dot_product(w,y)/sw                                               2729
      fmax=log(huge(y(1))*0.1)                                             2730
17760 do 17761 lam=1,nlam                                                  2730
      f=a0(lam)                                                            2731
17770 do 17771 k=1,nin(lam)                                                2731
      j=ia(k)                                                              2731
      jb=ix(j)                                                             2731
      je=ix(j+1)-1                                                         2732
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         2733
17771 continue                                                             2734
17772 continue                                                             2734
      f=f+g                                                                2735
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2736
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2737
17761 continue                                                             2738
17762 continue                                                             2738
11790 continue                                                             2738
      deallocate(w,f)                                                      2739
      return                                                               2740
      end                                                                  2741
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
