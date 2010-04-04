c
c                          newGLMnet (3/16/10)
c
c
c                 Elastic net with squared-error loss
c
c dense predictor matrix:
c
c call elnet(ka,parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd,
c            lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x(no,ni) = predictor data matrix flat file (overwritten)
c
c
c sparse predictor matrix:
c
c call spelnet(ka,parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,
c             isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c   x, ix, jx = predictor data matrix in compressed sparse row format
c
c
c other inputs:
c
c   ka = algorithm flag
c      ka=1 => covariance updating algorithm
c      ka=2 => naive algorithm
c   parm = family member index (0 <= parm <= 1)
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
c      iterations stop when the maximum standardized coefficient
c      change from the previous iteration is less than thr
c      (suggested value, thr=1.0e-4)
c   isd = standarization flag:
c      isd = 0 => regression on original predictor variables
c      isd = 1 => regression on standardized predictor variables
c      Note: output solutions always reference original
c            variables locations and scales.
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
c   nlp = total passes over the data summed over all lamda values
c   jerr = error flag:
c      jerr  = 0 => no error
c      jerr != 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 10000 => maxval(vp) <= 0.0
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
c   parm,no,ni,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd = same as elnet above.
c
c   nc = number of classes (distinct outcome values)
c        nc=1 => binomial two-class logistic regression
c            (all output references class 1)
c   y(no,max(2,nc)) = number of each class at each design point(overwritten)
c   o(no,nc) = observation off-sets for each class
c   maxit = maximum number of iterations allowed for any lamda value
c           (suggested value, maxit = 500)
c   kopt = optimization flag
c      kopt = 0 => Newton-Raphson
c      kpot = 1 => modified Newton-Raphson (recommended)
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
C      jerr < 0 => non fatal error - partial output returned
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations. Solutions for
c            larger lamdas returned
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value. Solutions for
c            larger lamdas returned
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
c   parm,no,ni,w,jd,vp,ne,nx,nlam,flmin,ulam,thr,isd = same as elnet above
c   maxit = same as lognet above
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
C      jerr < 0 => non fatal error - partial output returned:
c                  solutions for larger (1:lmu) lamdas returned.
c         jerr = -1 => inner loop convergence not reached
c            after maxit (see above) iterations.
c         jerr = -2 => outer loop convergence not reached
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
C      jerr < 0 => non fatal error - partial output:
c                  solutions for larger (1:lmu) lamdas returned.
c         jerr = -1 => inner loop convergence not reached
c            after maxit (see above) iterations.
c         jerr = -2 => outer loop convergence not reached
c            after maxit (see above) iterations.
c         jerr = -3 => numerical error.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
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
      subroutine elnet  (ka,parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam    605 
     *,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam)                          606
      real ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)                          607
      integer jd(*),ia(nx),nin(nlam)                                        608
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10021                                     611
      jerr=10000                                                            611
      return                                                                611
10021 continue                                                              612
      allocate(vq(1:ni),stat=jerr)                                          612
      if(jerr.ne.0) return                                                  613
      vq=max(0.0,vp)                                                        613
      vq=vq*ni/sum(vq)                                                      614
      if(ka .ne. 1)goto 10041                                               615
      call elnetu  (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd    618 
     *,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                            619
10041 continue                                                              620
      call elnetn (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd,    623 
     *  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                              624
10031 continue                                                              624
      deallocate(vq)                                                        625
      return                                                                626
      end                                                                   627
      subroutine elnetu  (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,t    630 
     *hr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ulam(nlam)                           631
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         632
      integer jd(*),ia(nx),nin(nlam)                                        633
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           638
      allocate(xm(1:ni),stat=ierr)                                          638
      jerr=jerr+ierr                                                        639
      allocate(xs(1:ni),stat=ierr)                                          639
      jerr=jerr+ierr                                                        640
      allocate(ju(1:ni),stat=ierr)                                          640
      jerr=jerr+ierr                                                        641
      allocate(xv(1:ni),stat=ierr)                                          641
      jerr=jerr+ierr                                                        642
      allocate(vlam(1:nlam),stat=ierr)                                      642
      jerr=jerr+ierr                                                        643
      if(jerr.ne.0) return                                                  644
      call chkvars(no,ni,x,ju)                                              645
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  646
      if(maxval(ju) .gt. 0)goto 10071                                       646
      jerr=7777                                                             646
      return                                                                646
10071 continue                                                              647
      call standard(no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)               648
      if(jerr.ne.0) return                                                  649
      if(flmin.ge.1.0) vlam=ulam/ys                                         650
      call elnet1(parm,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,vlam,thr,xv,  lm    652 
     *u,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                  653
10080 do 10081 k=1,lmu                                                      653
      alm(k)=ys*alm(k)                                                      653
      nk=nin(k)                                                             654
10090 do 10091 l=1,nk                                                       654
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          654
10091 continue                                                              655
10092 continue                                                              655
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         656
10081 continue                                                              657
10082 continue                                                              657
      deallocate(xm,xs,g,ju,xv,vlam)                                        658
      return                                                                659
      end                                                                   660
      subroutine standard (no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)        661
      real x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                  661
      integer ju(ni)                                                        662
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           665
      if(jerr.ne.0) return                                                  666
      w=w/sum(w)                                                            666
      v=sqrt(w)                                                             667
10100 do 10101 j=1,ni                                                       667
      if(ju(j).eq.0)goto 10101                                              668
      xm(j)=dot_product(w,x(:,j))                                           668
      x(:,j)=v*(x(:,j)-xm(j))                                               669
      xv(j)=dot_product(x(:,j),x(:,j))                                      670
10101 continue                                                              671
10102 continue                                                              671
      if(isd .ne. 0)goto 10121                                              671
      xs=1.0                                                                671
      goto 10131                                                            672
10121 continue                                                              673
10140 do 10141 j=1,ni                                                       673
      if(ju(j).eq.0)goto 10141                                              673
      xs(j)=sqrt(xv(j))                                                     673
      x(:,j)=x(:,j)/xs(j)                                                   673
10141 continue                                                              674
10142 continue                                                              674
      xv=1.0                                                                675
10131 continue                                                              676
10111 continue                                                              676
      ym=dot_product(w,y)                                                   676
      y=v*(y-ym)                                                            676
      ys=sqrt(dot_product(y,y))                                             676
      y=y/ys                                                                676
      g=0.0                                                                 677
10150 do 10151 j=1,ni                                                       677
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                             677
10151 continue                                                              678
10152 continue                                                              678
      deallocate(v)                                                         679
      return                                                                680
      end                                                                   681
      subroutine elnet1 (beta,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,ulam,thr,    683 
     *xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    684 
     *9)
      real vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    685 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       686
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                           692
      jerr=jerr+ierr                                                        693
      allocate(mm(1:ni),stat=ierr)                                          693
      jerr=jerr+ierr                                                        694
      allocate(da(1:ni),stat=ierr)                                          694
      jerr=jerr+ierr                                                        695
      if(jerr.ne.0) return                                                  696
      bta=max(beta,1.0e-3)                                                  696
      omb=1.0-bta                                                           697
      if(flmin .ge. 1.0)goto 10171                                          697
      eqs=max(eps,flmin)                                                    697
      alf=eqs**(1.0/(nlam-1))                                               697
10171 continue                                                              698
      rsq=0.0                                                               698
      a=0.0                                                                 698
      mm=0                                                                  698
      nlp=0                                                                 698
      nin=nlp                                                               698
      iz=0                                                                  698
      mnl=min(mnlam,nlam)                                                   699
10180 do 10181 m=1,nlam                                                     700
      if(flmin .lt. 1.0)goto 10201                                          700
      alm=ulam(m)                                                           700
      goto 10191                                                            701
10201 if(m .le. 2)goto 10211                                                701
      alm=alm*alf                                                           701
      goto 10191                                                            702
10211 if(m .ne. 1)goto 10221                                                702
      alm=big                                                               702
      goto 10231                                                            703
10221 continue                                                              703
      alm=0.0                                                               704
10240 do 10241 j=1,ni                                                       704
      if(ju(j).eq.0)goto 10241                                              704
      if(vp(j).le.0.0)goto 10241                                            705
      alm=max(alm,abs(g(j))/vp(j))                                          706
10241 continue                                                              707
10242 continue                                                              707
      alm=alf*alm/bta                                                       708
10231 continue                                                              709
10191 continue                                                              709
      dem=alm*omb                                                           709
      ab=alm*bta                                                            709
      rsq0=rsq                                                              709
      jz=1                                                                  710
10250 continue                                                              710
10251 continue                                                              710
      if(iz*jz.ne.0) go to 10260                                            710
      nlp=nlp+1                                                             710
      dlx=0.0                                                               711
10270 do 10271 k=1,ni                                                       711
      if(ju(k).eq.0)goto 10271                                              712
      ak=a(k)                                                               712
      u=g(k)+ak*xv(k)                                                       712
      v=abs(u)-vp(k)*ab                                                     712
      a(k)=0.0                                                              713
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         714
      if(a(k).eq.ak)goto 10271                                              715
      if(mm(k) .ne. 0)goto 10291                                            715
      nin=nin+1                                                             715
      if(nin.gt.nx)goto 10272                                               716
10300 do 10301 j=1,ni                                                       716
      if(ju(j).eq.0)goto 10301                                              717
      if(mm(j) .eq. 0)goto 10321                                            717
      c(j,nin)=c(k,mm(j))                                                   717
      goto 10301                                                            717
10321 continue                                                              718
      if(j .ne. k)goto 10341                                                718
      c(j,nin)=xv(j)                                                        718
      goto 10301                                                            718
10341 continue                                                              719
      c(j,nin)=dot_product(x(:,j),x(:,k))                                   720
10301 continue                                                              721
10302 continue                                                              721
      mm(k)=nin                                                             721
      ia(nin)=k                                                             722
10291 continue                                                              723
      del=a(k)-ak                                                           723
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      724
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                     725
10350 do 10351 j=1,ni                                                       725
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                               725
10351 continue                                                              726
10352 continue                                                              726
10271 continue                                                              727
10272 continue                                                              727
      if(dlx.lt.thr)goto 10252                                              727
      if(nin.gt.nx)goto 10252                                               728
10260 continue                                                              728
      iz=1                                                                  728
      da(1:nin)=a(ia(1:nin))                                                729
10360 continue                                                              729
10361 continue                                                              729
      nlp=nlp+1                                                             729
      dlx=0.0                                                               730
10370 do 10371 l=1,nin                                                      730
      k=ia(l)                                                               730
      ak=a(k)                                                               730
      u=g(k)+ak*xv(k)                                                       730
      v=abs(u)-vp(k)*ab                                                     731
      a(k)=0.0                                                              732
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         733
      if(a(k).eq.ak)goto 10371                                              734
      del=a(k)-ak                                                           734
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      735
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                     736
10380 do 10381 j=1,nin                                                      736
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                  736
10381 continue                                                              737
10382 continue                                                              737
10371 continue                                                              738
10372 continue                                                              738
      if(dlx.lt.thr)goto 10362                                              738
      goto 10361                                                            739
10362 continue                                                              739
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                      740
10390 do 10391 j=1,ni                                                       740
      if(mm(j).ne.0)goto 10391                                              741
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))            742
10391 continue                                                              743
10392 continue                                                              743
      jz=0                                                                  744
      goto 10251                                                            745
10252 continue                                                              745
      if(nin.gt.nx)goto 10182                                               746
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 746
      kin(m)=nin                                                            747
      rsqo(m)=rsq                                                           747
      almo(m)=alm                                                           747
      lmu=m                                                                 748
      if(m.lt.mnl)goto 10181                                                748
      if(flmin.ge.1.0)goto 10181                                            749
      me=0                                                                  749
10400 do 10401 j=1,nin                                                      749
      if(ao(j,m).ne.0.0) me=me+1                                            749
10401 continue                                                              749
10402 continue                                                              749
      if(me.gt.ne)goto 10182                                                750
      if(rsq-rsq0.lt.sml*rsq)goto 10182                                     750
      if(rsq.gt.rsqmax)goto 10182                                           751
10181 continue                                                              752
10182 continue                                                              752
      deallocate(a,mm,c,da)                                                 753
      return                                                                754
      end                                                                   755
      subroutine elnetn (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,th    757 
     *r,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no),w(no),ulam(nlam)                           758
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         759
      integer jd(*),ia(nx),nin(nlam)                                        760
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          765
      allocate(xs(1:ni),stat=ierr)                                          765
      jerr=jerr+ierr                                                        766
      allocate(ju(1:ni),stat=ierr)                                          766
      jerr=jerr+ierr                                                        767
      allocate(xv(1:ni),stat=ierr)                                          767
      jerr=jerr+ierr                                                        768
      allocate(vlam(1:nlam),stat=ierr)                                      768
      jerr=jerr+ierr                                                        769
      if(jerr.ne.0) return                                                  770
      call chkvars(no,ni,x,ju)                                              771
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  772
      if(maxval(ju) .gt. 0)goto 10421                                       772
      jerr=7777                                                             772
      return                                                                772
10421 continue                                                              773
      call standard1(no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)                774
      if(jerr.ne.0) return                                                  775
      if(flmin.ge.1.0) vlam=ulam/ys                                         776
      call elnet2(parm,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,vlam,thr,xv,  lm    778 
     *u,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                  779
10430 do 10431 k=1,lmu                                                      779
      alm(k)=ys*alm(k)                                                      779
      nk=nin(k)                                                             780
10440 do 10441 l=1,nk                                                       780
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          780
10441 continue                                                              781
10442 continue                                                              781
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         782
10431 continue                                                              783
10432 continue                                                              783
      deallocate(xm,xs,ju,xv,vlam)                                          784
      return                                                                785
      end                                                                   786
      subroutine standard1 (no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)         787
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)                        787
      integer ju(ni)                                                        788
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           791
      if(jerr.ne.0) return                                                  792
      w=w/sum(w)                                                            792
      v=sqrt(w)                                                             793
10450 do 10451 j=1,ni                                                       793
      if(ju(j).eq.0)goto 10451                                              794
      xm(j)=dot_product(w,x(:,j))                                           794
      x(:,j)=v*(x(:,j)-xm(j))                                               795
      xv(j)=dot_product(x(:,j),x(:,j))                                      796
10451 continue                                                              797
10452 continue                                                              797
      if(isd .ne. 0)goto 10471                                              797
      xs=1.0                                                                797
      goto 10481                                                            798
10471 continue                                                              798
10490 do 10491 j=1,ni                                                       798
      if(ju(j).eq.0)goto 10491                                              798
      xs(j)=sqrt(xv(j))                                                     798
      x(:,j)=x(:,j)/xs(j)                                                   798
10491 continue                                                              799
10492 continue                                                              799
      xv=1.0                                                                800
10481 continue                                                              801
10461 continue                                                              801
      ym=dot_product(w,y)                                                   801
      y=v*(y-ym)                                                            801
      ys=sqrt(dot_product(y,y))                                             801
      y=y/ys                                                                802
      deallocate(v)                                                         803
      return                                                                804
      end                                                                   805
      subroutine elnet2(beta,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,ulam,thr,x    807 
     *v,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    808 
     *9)
      real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    809 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       810
      real, dimension (:), allocatable :: a                                     
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                           815
      allocate(mm(1:ni),stat=ierr)                                          815
      jerr=jerr+ierr                                                        816
      if(jerr.ne.0) return                                                  817
      bta=max(beta,1.0e-3)                                                  817
      omb=1.0-bta                                                           818
      if(flmin .ge. 1.0)goto 10511                                          818
      eqs=max(eps,flmin)                                                    818
      alf=eqs**(1.0/(nlam-1))                                               818
10511 continue                                                              819
      rsq=0.0                                                               819
      a=0.0                                                                 819
      mm=0                                                                  819
      nlp=0                                                                 819
      nin=nlp                                                               819
      iz=0                                                                  819
      mnl=min(mnlam,nlam)                                                   820
10520 do 10521 m=1,nlam                                                     821
      if(flmin .lt. 1.0)goto 10541                                          821
      alm=ulam(m)                                                           821
      goto 10531                                                            822
10541 if(m .le. 2)goto 10551                                                822
      alm=alm*alf                                                           822
      goto 10531                                                            823
10551 if(m .ne. 1)goto 10561                                                823
      alm=big                                                               823
      goto 10571                                                            824
10561 continue                                                              824
      alm=0.0                                                               825
10580 do 10581 j=1,ni                                                       825
      if(ju(j).eq.0)goto 10581                                              825
      if(vp(j).le.0.0)goto 10581                                            826
      alm=max(alm,abs(dot_product(y,x(:,j)))/vp(j))                         827
10581 continue                                                              828
10582 continue                                                              828
      alm=alf*alm/bta                                                       829
10571 continue                                                              830
10531 continue                                                              830
      dem=alm*omb                                                           830
      ab=alm*bta                                                            830
      rsq0=rsq                                                              830
      jz=1                                                                  831
10590 continue                                                              831
10591 continue                                                              831
      if(iz*jz.ne.0) go to 10260                                            831
      nlp=nlp+1                                                             831
      dlx=0.0                                                               832
10600 do 10601 k=1,ni                                                       832
      if(ju(k).eq.0)goto 10601                                              832
      gk=dot_product(y,x(:,k))                                              833
      ak=a(k)                                                               833
      u=gk+ak*xv(k)                                                         833
      v=abs(u)-vp(k)*ab                                                     833
      a(k)=0.0                                                              834
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         835
      if(a(k).eq.ak)goto 10601                                              836
      if(mm(k) .ne. 0)goto 10621                                            836
      nin=nin+1                                                             836
      if(nin.gt.nx)goto 10602                                               837
      mm(k)=nin                                                             837
      ia(nin)=k                                                             838
10621 continue                                                              839
      del=a(k)-ak                                                           839
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        840
      y=y-del*x(:,k)                                                        840
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                     841
10601 continue                                                              842
10602 continue                                                              842
      if(dlx.lt.thr)goto 10592                                              842
      if(nin.gt.nx)goto 10592                                               843
10260 continue                                                              843
      iz=1                                                                  844
10630 continue                                                              844
10631 continue                                                              844
      nlp=nlp+1                                                             844
      dlx=0.0                                                               845
10640 do 10641 l=1,nin                                                      845
      k=ia(l)                                                               845
      gk=dot_product(y,x(:,k))                                              846
      ak=a(k)                                                               846
      u=gk+ak*xv(k)                                                         846
      v=abs(u)-vp(k)*ab                                                     846
      a(k)=0.0                                                              847
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         848
      if(a(k).eq.ak)goto 10641                                              849
      del=a(k)-ak                                                           849
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        850
      y=y-del*x(:,k)                                                        850
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                     851
10641 continue                                                              852
10642 continue                                                              852
      if(dlx.lt.thr)goto 10632                                              852
      goto 10631                                                            853
10632 continue                                                              853
      jz=0                                                                  854
      goto 10591                                                            855
10592 continue                                                              855
      if(nin.gt.nx)goto 10522                                               856
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 856
      kin(m)=nin                                                            857
      rsqo(m)=rsq                                                           857
      almo(m)=alm                                                           857
      lmu=m                                                                 858
      if(m.lt.mnl)goto 10521                                                858
      if(flmin.ge.1.0)goto 10521                                            859
      me=0                                                                  859
10650 do 10651 j=1,nin                                                      859
      if(ao(j,m).ne.0.0) me=me+1                                            859
10651 continue                                                              859
10652 continue                                                              859
      if(me.gt.ne)goto 10522                                                860
      if(rsq-rsq0.lt.sml*rsq)goto 10522                                     860
      if(rsq.gt.rsqmax)goto 10522                                           861
10521 continue                                                              862
10522 continue                                                              862
      deallocate(a,mm)                                                      863
      return                                                                864
      end                                                                   865
      subroutine chkvars(no,ni,x,ju)                                        866
      real x(no,ni)                                                         866
      integer ju(ni)                                                        867
10660 do 10661 j=1,ni                                                       867
      ju(j)=0                                                               867
      t=x(1,j)                                                              868
10670 do 10671 i=2,no                                                       868
      if(x(i,j).eq.t)goto 10671                                             868
      ju(j)=1                                                               868
      goto 10672                                                            868
10671 continue                                                              869
10672 continue                                                              869
10661 continue                                                              870
10662 continue                                                              870
      return                                                                871
      end                                                                   872
      subroutine uncomp(ni,ca,ia,nin,a)                                     873
      real ca(*),a(ni)                                                      873
      integer ia(*)                                                         874
      a=0.0                                                                 874
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                   875
      return                                                                876
      end                                                                   877
      subroutine modval(a0,ca,ia,nin,n,x,f)                                 878
      real ca(nin),x(n,*),f(n)                                              878
      integer ia(nin)                                                       879
      f=a0                                                                  879
      if(nin.le.0) return                                                   880
10680 do 10681 i=1,n                                                        880
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                       880
10681 continue                                                              881
10682 continue                                                              881
      return                                                                882
      end                                                                   883
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,fl    886 
     *min,ulam,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                               887
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         888
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                            889
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10701                                     892
      jerr=10000                                                            892
      return                                                                892
10701 continue                                                              893
      allocate(vq(1:ni),stat=jerr)                                          893
      if(jerr.ne.0) return                                                  894
      vq=max(0.0,vp)                                                        894
      vq=vq*ni/sum(vq)                                                      895
      if(ka .ne. 1)goto 10721                                               896
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam    899 
     *,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10731                                                            900
10721 continue                                                              901
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam,    904 
     *thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10731 continue                                                              905
10711 continue                                                              905
      deallocate(vq)                                                        906
      return                                                                907
      end                                                                   908
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmi    911 
     *n,ulam,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                               912
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         913
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                            914
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           919
      allocate(xm(1:ni),stat=ierr)                                          919
      jerr=jerr+ierr                                                        920
      allocate(xs(1:ni),stat=ierr)                                          920
      jerr=jerr+ierr                                                        921
      allocate(ju(1:ni),stat=ierr)                                          921
      jerr=jerr+ierr                                                        922
      allocate(xv(1:ni),stat=ierr)                                          922
      jerr=jerr+ierr                                                        923
      allocate(vlam(1:nlam),stat=ierr)                                      923
      jerr=jerr+ierr                                                        924
      if(jerr.ne.0) return                                                  925
      call spchkvars(no,ni,x,ix,ju)                                         926
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  927
      if(maxval(ju) .gt. 0)goto 10751                                       927
      jerr=7777                                                             927
      return                                                                927
10751 continue                                                              928
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,jerr)       929
      if(jerr.ne.0) return                                                  930
      if(flmin.ge.1.0) vlam=ulam/ys                                         931
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t    933 
     *hr,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                  934
10760 do 10761 k=1,lmu                                                      934
      alm(k)=ys*alm(k)                                                      934
      nk=nin(k)                                                             935
10770 do 10771 l=1,nk                                                       935
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          935
10771 continue                                                              936
10772 continue                                                              936
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         937
10761 continue                                                              938
10762 continue                                                              938
      deallocate(xm,xs,g,ju,xv,vlam)                                        939
      return                                                                940
      end                                                                   941
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,j    942 
     *err)
      real x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                      942
      integer ix(*),jx(*),ju(ni)                                            943
      w=w/sum(w)                                                            944
10780 do 10781 j=1,ni                                                       944
      if(ju(j).eq.0)goto 10781                                              945
      jb=ix(j)                                                              945
      je=ix(j+1)-1                                                          945
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                              946
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                  947
10781 continue                                                              948
10782 continue                                                              948
      if(isd .ne. 0)goto 10801                                              948
      xs=1.0                                                                948
      goto 10811                                                            949
10801 continue                                                              949
10820 do 10821 j=1,ni                                                       949
      if(ju(j).ne.0) xs(j)=sqrt(xv(j))                                      949
10821 continue                                                              949
10822 continue                                                              949
      xv=1.0                                                                949
10811 continue                                                              950
10791 continue                                                              950
      ym=dot_product(w,y)                                                   950
      y=y-ym                                                                950
      ys=sqrt(dot_product(w,y**2))                                          950
      y=y/ys                                                                950
      g=0.0                                                                 951
10830 do 10831 j=1,ni                                                       951
      if(ju(j).eq.0)goto 10831                                              951
      jb=ix(j)                                                              951
      je=ix(j+1)-1                                                          952
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)            953
10831 continue                                                              954
10832 continue                                                              954
      return                                                                955
      end                                                                   956
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,    958 
     *ulam,  thr,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    959 
     *9)
      real g(ni),vp(ni),x(*),ulam(nlam),w(no)                               960
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)           961
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                           962
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                           968
      jerr=jerr+ierr                                                        969
      allocate(mm(1:ni),stat=ierr)                                          969
      jerr=jerr+ierr                                                        970
      allocate(da(1:ni),stat=ierr)                                          970
      jerr=jerr+ierr                                                        971
      if(jerr.ne.0) return                                                  972
      bta=max(beta,1.0e-3)                                                  972
      omb=1.0-bta                                                           973
      if(flmin .ge. 1.0)goto 10851                                          973
      eqs=max(eps,flmin)                                                    973
      alf=eqs**(1.0/(nlam-1))                                               973
10851 continue                                                              974
      rsq=0.0                                                               974
      a=0.0                                                                 974
      mm=0                                                                  974
      nlp=0                                                                 974
      nin=nlp                                                               974
      iz=0                                                                  974
      mnl=min(mnlam,nlam)                                                   975
10860 do 10861 m=1,nlam                                                     976
      if(flmin .lt. 1.0)goto 10881                                          976
      alm=ulam(m)                                                           976
      goto 10871                                                            977
10881 if(m .le. 2)goto 10891                                                977
      alm=alm*alf                                                           977
      goto 10871                                                            978
10891 if(m .ne. 1)goto 10901                                                978
      alm=big                                                               978
      goto 10911                                                            979
10901 continue                                                              979
      alm=0.0                                                               980
10920 do 10921 j=1,ni                                                       980
      if(ju(j).eq.0)goto 10921                                              980
      if(vp(j).le.0.0)goto 10921                                            981
      alm=max(alm,abs(g(j))/vp(j))                                          982
10921 continue                                                              983
10922 continue                                                              983
      alm=alf*alm/bta                                                       984
10911 continue                                                              985
10871 continue                                                              985
      dem=alm*omb                                                           985
      ab=alm*bta                                                            985
      rsq0=rsq                                                              985
      jz=1                                                                  986
10930 continue                                                              986
10931 continue                                                              986
      if(iz*jz.ne.0) go to 10260                                            986
      nlp=nlp+1                                                             986
      dlx=0.0                                                               987
10940 do 10941 k=1,ni                                                       987
      if(ju(k).eq.0)goto 10941                                              988
      ak=a(k)                                                               988
      u=g(k)+ak*xv(k)                                                       988
      v=abs(u)-vp(k)*ab                                                     988
      a(k)=0.0                                                              989
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         990
      if(a(k).eq.ak)goto 10941                                              991
      if(mm(k) .ne. 0)goto 10961                                            991
      nin=nin+1                                                             991
      if(nin.gt.nx)goto 10942                                               992
10970 do 10971 j=1,ni                                                       992
      if(ju(j).eq.0)goto 10971                                              993
      if(mm(j) .eq. 0)goto 10991                                            993
      c(j,nin)=c(k,mm(j))                                                   993
      goto 10971                                                            993
10991 continue                                                              994
      if(j .ne. k)goto 11011                                                994
      c(j,nin)=xv(j)                                                        994
      goto 10971                                                            994
11011 continue                                                              995
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))        997
10971 continue                                                              998
10972 continue                                                              998
      mm(k)=nin                                                             998
      ia(nin)=k                                                             999
10961 continue                                                             1000
      del=a(k)-ak                                                          1000
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1001
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                    1002
11020 do 11021 j=1,ni                                                      1002
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                              1002
11021 continue                                                             1003
11022 continue                                                             1003
10941 continue                                                             1004
10942 continue                                                             1004
      if(dlx.lt.thr)goto 10932                                             1004
      if(nin.gt.nx)goto 10932                                              1005
10260 continue                                                             1005
      iz=1                                                                 1005
      da(1:nin)=a(ia(1:nin))                                               1006
11030 continue                                                             1006
11031 continue                                                             1006
      nlp=nlp+1                                                            1006
      dlx=0.0                                                              1007
11040 do 11041 l=1,nin                                                     1007
      k=ia(l)                                                              1008
      ak=a(k)                                                              1008
      u=g(k)+ak*xv(k)                                                      1008
      v=abs(u)-vp(k)*ab                                                    1008
      a(k)=0.0                                                             1009
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1010
      if(a(k).eq.ak)goto 11041                                             1011
      del=a(k)-ak                                                          1011
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1012
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                    1013
11050 do 11051 j=1,nin                                                     1013
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                 1013
11051 continue                                                             1014
11052 continue                                                             1014
11041 continue                                                             1015
11042 continue                                                             1015
      if(dlx.lt.thr)goto 11032                                             1015
      goto 11031                                                           1016
11032 continue                                                             1016
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                     1017
11060 do 11061 j=1,ni                                                      1017
      if(mm(j).ne.0)goto 11061                                             1018
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))           1019
11061 continue                                                             1020
11062 continue                                                             1020
      jz=0                                                                 1021
      goto 10931                                                           1022
10932 continue                                                             1022
      if(nin.gt.nx)goto 10862                                              1023
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1023
      kin(m)=nin                                                           1024
      rsqo(m)=rsq                                                          1024
      almo(m)=alm                                                          1024
      lmu=m                                                                1025
      if(m.lt.mnl)goto 10861                                               1025
      if(flmin.ge.1.0)goto 10861                                           1026
      me=0                                                                 1026
11070 do 11071 j=1,nin                                                     1026
      if(ao(j,m).ne.0.0) me=me+1                                           1026
11071 continue                                                             1026
11072 continue                                                             1026
      if(me.gt.ne)goto 10862                                               1027
      if(rsq-rsq0.lt.sml*rsq)goto 10862                                    1027
      if(rsq.gt.rsqmax)goto 10862                                          1028
10861 continue                                                             1029
10862 continue                                                             1029
      deallocate(a,mm,c,da)                                                1030
      return                                                               1031
      end                                                                  1032
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmin,   1034 
     *ulam,  thr,isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),vp(ni),y(no),w(no),ulam(nlam)                              1035
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1036
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1037
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         1042
      allocate(xs(1:ni),stat=ierr)                                         1042
      jerr=jerr+ierr                                                       1043
      allocate(ju(1:ni),stat=ierr)                                         1043
      jerr=jerr+ierr                                                       1044
      allocate(xv(1:ni),stat=ierr)                                         1044
      jerr=jerr+ierr                                                       1045
      allocate(vlam(1:nlam),stat=ierr)                                     1045
      jerr=jerr+ierr                                                       1046
      if(jerr.ne.0) return                                                 1047
      call spchkvars(no,ni,x,ix,ju)                                        1048
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1049
      if(maxval(ju) .gt. 0)goto 11091                                      1049
      jerr=7777                                                            1049
      return                                                               1049
11091 continue                                                             1050
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,jerr)       1051
      if(jerr.ne.0) return                                                 1052
      if(flmin.ge.1.0) vlam=ulam/ys                                        1053
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t   1055 
     *hr,xm,xs,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                 1056
11100 do 11101 k=1,lmu                                                     1056
      alm(k)=ys*alm(k)                                                     1056
      nk=nin(k)                                                            1057
11110 do 11111 l=1,nk                                                      1057
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1057
11111 continue                                                             1058
11112 continue                                                             1058
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                        1059
11101 continue                                                             1060
11102 continue                                                             1060
      deallocate(xm,xs,ju,xv,vlam)                                         1061
      return                                                               1062
      end                                                                  1063
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,je   1064 
     *rr)
      real x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)                           1064
      integer ix(*),jx(*),ju(ni)                                           1065
      w=w/sum(w)                                                           1066
11120 do 11121 j=1,ni                                                      1066
      if(ju(j).eq.0)goto 11121                                             1067
      jb=ix(j)                                                             1067
      je=ix(j+1)-1                                                         1067
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1068
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1069
11121 continue                                                             1070
11122 continue                                                             1070
      if(isd .ne. 0)goto 11141                                             1070
      xs=1.0                                                               1070
      goto 11151                                                           1071
11141 continue                                                             1071
11160 do 11161 j=1,ni                                                      1071
      if(ju(j).ne.0) xs(j)=sqrt(xv(j))                                     1071
11161 continue                                                             1071
11162 continue                                                             1071
      xv=1.0                                                               1071
11151 continue                                                             1072
11131 continue                                                             1072
      ym=dot_product(w,y)                                                  1072
      y=y-ym                                                               1072
      ys=sqrt(dot_product(w,y**2))                                         1072
      y=y/ys                                                               1073
      return                                                               1074
      end                                                                  1075
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,   1077 
     *ulam,  thr,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99   1078 
     *9)
      real y(no),w(no),x(*),vp(ni),ulam(nlam)                              1079
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)          1080
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1081
      real, dimension (:), allocatable :: a                                     
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                          1086
      allocate(mm(1:ni),stat=ierr)                                         1086
      jerr=jerr+ierr                                                       1087
      if(jerr.ne.0) return                                                 1088
      bta=max(beta,1.0e-3)                                                 1088
      omb=1.0-bta                                                          1089
      if(flmin .ge. 1.0)goto 11181                                         1089
      eqs=max(eps,flmin)                                                   1089
      alf=eqs**(1.0/(nlam-1))                                              1089
11181 continue                                                             1090
      rsq=0.0                                                              1090
      a=0.0                                                                1090
      mm=0                                                                 1090
      o=0.0                                                                1090
      nlp=0                                                                1090
      nin=nlp                                                              1090
      iz=0                                                                 1090
      mnl=min(mnlam,nlam)                                                  1091
11190 do 11191 m=1,nlam                                                    1092
      if(flmin .lt. 1.0)goto 11211                                         1092
      alm=ulam(m)                                                          1092
      goto 11201                                                           1093
11211 if(m .le. 2)goto 11221                                               1093
      alm=alm*alf                                                          1093
      goto 11201                                                           1094
11221 if(m .ne. 1)goto 11231                                               1094
      alm=big                                                              1094
      goto 11241                                                           1095
11231 continue                                                             1095
      alm=0.0                                                              1096
11250 do 11251 j=1,ni                                                      1096
      if(ju(j).eq.0)goto 11251                                             1096
      if(vp(j).le.0.0)goto 11251                                           1097
      jb=ix(j)                                                             1097
      je=ix(j+1)-1                                                         1098
      alm=max(alm,abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))    1100 
     * /(vp(j)*xs(j))))
11251 continue                                                             1101
11252 continue                                                             1101
      alm=alf*alm/bta                                                      1102
11241 continue                                                             1103
11201 continue                                                             1103
      dem=alm*omb                                                          1103
      ab=alm*bta                                                           1103
      rsq0=rsq                                                             1103
      jz=1                                                                 1104
11260 continue                                                             1104
11261 continue                                                             1104
      if(iz*jz.ne.0) go to 10260                                           1104
      nlp=nlp+1                                                            1104
      dlx=0.0                                                              1105
11270 do 11271 k=1,ni                                                      1105
      if(ju(k).eq.0)goto 11271                                             1105
      jb=ix(k)                                                             1105
      je=ix(k+1)-1                                                         1106
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1107
      ak=a(k)                                                              1107
      u=gk+ak*xv(k)                                                        1107
      v=abs(u)-vp(k)*ab                                                    1107
      a(k)=0.0                                                             1108
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1109
      if(a(k).eq.ak)goto 11271                                             1110
      if(mm(k) .ne. 0)goto 11291                                           1110
      nin=nin+1                                                            1110
      if(nin.gt.nx)goto 11272                                              1111
      mm(k)=nin                                                            1111
      ia(nin)=k                                                            1112
11291 continue                                                             1113
      del=a(k)-ak                                                          1113
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1114
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1115
      o=o+del*xm(k)/xs(k)                                                  1115
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                    1116
11271 continue                                                             1117
11272 continue                                                             1117
      if(dlx.lt.thr)goto 11262                                             1117
      if(nin.gt.nx)goto 11262                                              1118
10260 continue                                                             1118
      iz=1                                                                 1119
11300 continue                                                             1119
11301 continue                                                             1119
      nlp=nlp+1                                                            1119
      dlx=0.0                                                              1120
11310 do 11311 l=1,nin                                                     1120
      k=ia(l)                                                              1120
      jb=ix(k)                                                             1120
      je=ix(k+1)-1                                                         1121
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1122
      ak=a(k)                                                              1122
      u=gk+ak*xv(k)                                                        1122
      v=abs(u)-vp(k)*ab                                                    1122
      a(k)=0.0                                                             1123
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1124
      if(a(k).eq.ak)goto 11311                                             1125
      del=a(k)-ak                                                          1125
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1126
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1127
      o=o+del*xm(k)/xs(k)                                                  1127
      dlx=max(abs(del)/sqrt(xv(k)),dlx)                                    1128
11311 continue                                                             1129
11312 continue                                                             1129
      if(dlx.lt.thr)goto 11302                                             1129
      goto 11301                                                           1130
11302 continue                                                             1130
      jz=0                                                                 1131
      goto 11261                                                           1132
11262 continue                                                             1132
      if(nin.gt.nx)goto 11192                                              1133
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1133
      kin(m)=nin                                                           1134
      rsqo(m)=rsq                                                          1134
      almo(m)=alm                                                          1134
      lmu=m                                                                1135
      if(m.lt.mnl)goto 11191                                               1135
      if(flmin.ge.1.0)goto 11191                                           1136
      me=0                                                                 1136
11320 do 11321 j=1,nin                                                     1136
      if(ao(j,m).ne.0.0) me=me+1                                           1136
11321 continue                                                             1136
11322 continue                                                             1136
      if(me.gt.ne)goto 11192                                               1137
      if(rsq-rsq0.lt.sml*rsq)goto 11192                                    1137
      if(rsq.gt.rsqmax)goto 11192                                          1138
11191 continue                                                             1139
11192 continue                                                             1139
      deallocate(a,mm)                                                     1140
      return                                                               1141
      end                                                                  1142
      subroutine spchkvars(no,ni,x,ix,ju)                                  1143
      real x(*)                                                            1143
      integer ix(*),ju(ni)                                                 1144
11330 do 11331 j=1,ni                                                      1144
      ju(j)=0                                                              1144
      jb=ix(j)                                                             1144
      nj=ix(j+1)-jb                                                        1144
      if(nj.eq.0)goto 11331                                                1145
      je=ix(j+1)-1                                                         1146
      if(nj .ge. no)goto 11351                                             1146
11360 do 11361 i=jb,je                                                     1146
      if(x(i).eq.0.0)goto 11361                                            1146
      ju(j)=1                                                              1146
      goto 11362                                                           1146
11361 continue                                                             1146
11362 continue                                                             1146
      goto 11371                                                           1147
11351 continue                                                             1147
      t=x(jb)                                                              1147
11380 do 11381 i=jb+1,je                                                   1147
      if(x(i).eq.t)goto 11381                                              1147
      ju(j)=1                                                              1147
      goto 11382                                                           1147
11381 continue                                                             1147
11382 continue                                                             1147
11371 continue                                                             1148
11341 continue                                                             1148
11331 continue                                                             1149
11332 continue                                                             1149
      return                                                               1150
      end                                                                  1151
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                         1152
      real ca(*),x(*),f(n)                                                 1152
      integer ia(*),ix(*),jx(*)                                            1153
      f=a0                                                                 1154
11390 do 11391 j=1,nin                                                     1154
      k=ia(j)                                                              1154
      kb=ix(k)                                                             1154
      ke=ix(k+1)-1                                                         1155
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                             1156
11391 continue                                                             1157
11392 continue                                                             1157
      return                                                               1158
      end                                                                  1159
      function row_prod(i,j,ia,ja,ra,w)                                    1160
      integer ia(*),ja(*)                                                  1160
      real ra(*),w(*)                                                      1161
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(   1163 
     *i),ia(j+1)-ia(j),w)
      return                                                               1164
      end                                                                  1165
      function dot(x,y,mx,my,nx,ny,w)                                      1166
      real x(*),y(*),w(*)                                                  1166
      integer mx(*),my(*)                                                  1167
      i=1                                                                  1167
      j=i                                                                  1167
      s=0.0                                                                1168
11400 continue                                                             1168
11401 continue                                                             1168
11410 continue                                                             1169
11411 if(mx(i).ge.my(j))goto 11412                                         1169
      i=i+1                                                                1169
      if(i.gt.nx) go to 11420                                              1169
      goto 11411                                                           1170
11412 continue                                                             1170
      if(mx(i).eq.my(j)) go to 11430                                       1171
11440 continue                                                             1171
11441 if(my(j).ge.mx(i))goto 11442                                         1171
      j=j+1                                                                1171
      if(j.gt.ny) go to 11420                                              1171
      goto 11441                                                           1172
11442 continue                                                             1172
      if(mx(i).eq.my(j)) go to 11430                                       1172
      goto 11401                                                           1173
11430 continue                                                             1173
      s=s+w(mx(i))*x(i)*y(j)                                               1174
      i=i+1                                                                1174
      if(i.gt.nx)goto 11402                                                1174
      j=j+1                                                                1174
      if(j.gt.ny)goto 11402                                                1175
      goto 11401                                                           1176
11402 continue                                                             1176
11420 continue                                                             1176
      dot=s                                                                1177
      return                                                               1178
      end                                                                  1179
      subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,ne,nx,nlam,flmin,ulam   1181 
     *,thr,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)             1182
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1183
      integer jd(*),ia(nx),nin(nlam)                                       1184
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 11461                                    1188
      jerr=10000                                                           1188
      return                                                               1188
11461 continue                                                             1189
      allocate(ww(1:no),stat=jerr)                                         1190
      allocate(ju(1:ni),stat=ierr)                                         1190
      jerr=jerr+ierr                                                       1191
      allocate(vq(1:ni),stat=ierr)                                         1191
      jerr=jerr+ierr                                                       1192
      allocate(xm(1:ni),stat=ierr)                                         1192
      jerr=jerr+ierr                                                       1193
      if(isd .le. 0)goto 11481                                             1193
      allocate(xs(1:ni),stat=ierr)                                         1193
      jerr=jerr+ierr                                                       1193
11481 continue                                                             1194
      if(jerr.ne.0) return                                                 1195
      call chkvars(no,ni,x,ju)                                             1196
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1197
      if(maxval(ju) .gt. 0)goto 11501                                      1197
      jerr=7777                                                            1197
      return                                                               1197
11501 continue                                                             1198
      vq=max(0.0,vp)                                                       1198
      vq=vq*ni/sum(vq)                                                     1199
11510 do 11511 i=1,no                                                      1199
      ww(i)=sum(y(i,:))                                                    1199
      y(i,:)=y(i,:)/ww(i)                                                  1199
11511 continue                                                             1199
11512 continue                                                             1199
      sw=sum(ww)                                                           1199
      ww=ww/sw                                                             1200
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             1201
      if(nc .ne. 1)goto 11531                                              1202
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,ne,nx,nlam,flmin   1204 
     *,ulam,  thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 11541                                                           1205
11531 continue                                                             1206
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,th   1208 
     *r,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
11541 continue                                                             1209
11521 continue                                                             1209
      if(jerr.gt.0) return                                                 1209
      dev0=2.0*sw*dev0                                                     1210
11550 do 11551 k=1,lmu                                                     1210
      nk=nin(k)                                                            1211
11560 do 11561 ic=1,nc                                                     1211
      if(isd .le. 0)goto 11581                                             1211
11590 do 11591 l=1,nk                                                      1211
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1211
11591 continue                                                             1211
11592 continue                                                             1211
11581 continue                                                             1212
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1213
11561 continue                                                             1214
11562 continue                                                             1214
11551 continue                                                             1215
11552 continue                                                             1215
      deallocate(ww,ju,vq,xm)                                              1215
      if(isd.gt.0) deallocate(xs)                                          1216
      return                                                               1217
      end                                                                  1218
      subroutine lstandard1 (no,ni,x,w,ju,isd,xm,xs)                       1219
      real x(no,ni),w(no),xm(ni),xs(ni)                                    1219
      integer ju(ni)                                                       1220
11600 do 11601 j=1,ni                                                      1220
      if(ju(j).eq.0)goto 11601                                             1221
      xm(j)=dot_product(w,x(:,j))                                          1221
      x(1:no,j)=x(1:no,j)-xm(j)                                            1222
      if(isd .le. 0)goto 11621                                             1223
      xs(j)=sqrt(dot_product(w*x(:,j),x(:,j)))                             1224
      x(1:no,j)=x(1:no,j)/xs(j)                                            1225
11621 continue                                                             1226
11601 continue                                                             1227
11602 continue                                                             1227
      return                                                               1228
      end                                                                  1229
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ulam   1231 
     *,shr,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1233 
     *5, devmax=0.999)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    1234
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1235
      integer ju(ni),m(nx),kin(nlam)                                       1236
      real, dimension (:), allocatable :: b,bs,v,r,xv,q                         
      integer, dimension (:), allocatable :: mm                                 
      allocate(b(0:ni),stat=jerr)                                          1241
      allocate(xv(1:ni),stat=ierr)                                         1241
      jerr=jerr+ierr                                                       1242
      allocate(bs(0:ni),stat=ierr)                                         1242
      jerr=jerr+ierr                                                       1243
      allocate(mm(1:ni),stat=ierr)                                         1243
      jerr=jerr+ierr                                                       1244
      allocate(r(1:no),stat=ierr)                                          1244
      jerr=jerr+ierr                                                       1245
      allocate(v(1:no),stat=ierr)                                          1245
      jerr=jerr+ierr                                                       1246
      allocate(q(1:no),stat=ierr)                                          1246
      jerr=jerr+ierr                                                       1247
      if(jerr.ne.0) return                                                 1248
      fmax=log(1.0/pmin-1.0)                                               1248
      fmin=-fmax                                                           1248
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1249
      bta=max(parm,1.0e-3)                                                 1249
      omb=1.0-bta                                                          1250
      q0=dot_product(w,y)                                                  1250
      if(q0 .gt. pmin)goto 11641                                           1250
      jerr=8001                                                            1250
      return                                                               1250
11641 continue                                                             1251
      if(q0 .lt. 1.0-pmin)goto 11661                                       1251
      jerr=9001                                                            1251
      return                                                               1251
11661 continue                                                             1252
      bz=log(q0/(1.0-q0))                                                  1253
      if(nonzero(no,g) .ne. 0)goto 11681                                   1253
      vi=q0*(1.0-q0)                                                       1253
      b(0)=bz                                                              1253
      v=vi*w                                                               1254
      r=w*(y-q0)                                                           1254
      q=q0                                                                 1254
      xmz=vi                                                               1255
      goto 11691                                                           1256
11681 continue                                                             1256
      b(0)=azero(no,y,g,w,jerr)                                            1256
      if(jerr.ne.0) return                                                 1257
      q=1.0/(1.0+exp(-b(0)-g))                                             1257
      v=w*q*(1.0-q)                                                        1257
      r=w*(y-q)                                                            1257
      xmz=sum(v)                                                           1258
11691 continue                                                             1259
11671 continue                                                             1259
      if(isd .le. 0)goto 11711                                             1259
      xv=0.25                                                              1259
      goto 11721                                                           1260
11711 continue                                                             1260
11730 do 11731 j=1,ni                                                      1260
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1260
11731 continue                                                             1260
11732 continue                                                             1260
11721 continue                                                             1261
11701 continue                                                             1261
      dev1=-(bz*q0+log(1.0-q0))                                            1261
      dev0=dev1                                                            1262
11740 do 11741 i=1,no                                                      1262
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1263
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1264
11741 continue                                                             1265
11742 continue                                                             1265
      if(flmin .ge. 1.0)goto 11761                                         1265
      eqs=max(eps,flmin)                                                   1265
      alf=eqs**(1.0/(nlam-1))                                              1265
11761 continue                                                             1266
      m=0                                                                  1266
      mm=0                                                                 1266
      nlp=0                                                                1266
      nin=nlp                                                              1266
      mnl=min(mnlam,nlam)                                                  1266
      bs=0.0                                                               1266
      b(1:ni)=0.0                                                          1267
11770 do 11771 ilm=1,nlam                                                  1268
      if(flmin .lt. 1.0)goto 11791                                         1268
      al=ulam(ilm)                                                         1268
      goto 11781                                                           1269
11791 if(ilm .le. 2)goto 11801                                             1269
      al=al*alf                                                            1269
      goto 11781                                                           1270
11801 if(ilm .ne. 1)goto 11811                                             1270
      al=big                                                               1270
      goto 11821                                                           1271
11811 continue                                                             1271
      al=0.0                                                               1272
11830 do 11831 j=1,ni                                                      1272
      if(ju(j).eq.0)goto 11831                                             1272
      if(vp(j).le.0.0)goto 11831                                           1273
      al=max(al,abs(dot_product(r,x(:,j)))/vp(j))                          1274
11831 continue                                                             1275
11832 continue                                                             1275
      al=alf*al/bta                                                        1276
11821 continue                                                             1277
11781 continue                                                             1277
      al2=al*omb                                                           1277
      al1=al*bta                                                           1277
      nit=0                                                                1278
11840 continue                                                             1278
11841 continue                                                             1278
      bs(0)=b(0)                                                           1278
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1279
11850 continue                                                             1279
11851 continue                                                             1279
      nlp=nlp+1                                                            1279
      dlx=0.0                                                              1280
11860 do 11861 k=1,ni                                                      1280
      if(ju(k).eq.0)goto 11861                                             1281
      bk=b(k)                                                              1281
      gk=dot_product(r,x(:,k))                                             1282
      u=gk+xv(k)*b(k)                                                      1282
      au=abs(u)-vp(k)*al1                                                  1283
      if(au .gt. 0.0)goto 11881                                            1283
      b(k)=0.0                                                             1283
      goto 11891                                                           1284
11881 continue                                                             1284
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1284
11891 continue                                                             1285
11871 continue                                                             1285
      d=b(k)-bk                                                            1285
      if(abs(d).le.0.0)goto 11861                                          1285
      dlx=max(dlx,abs(d))                                                  1286
      r=r-d*v*x(:,k)                                                       1287
      if(mm(k) .ne. 0)goto 11911                                           1287
      nin=nin+1                                                            1287
      if(nin.gt.nx)goto 11862                                              1288
      mm(k)=nin                                                            1288
      m(nin)=k                                                             1289
11911 continue                                                             1290
11861 continue                                                             1291
11862 continue                                                             1291
      if(nin.gt.nx)goto 11852                                              1292
      d=sum(r)/xmz                                                         1293
      if(d .eq. 0.0)goto 11931                                             1293
      b(0)=b(0)+d                                                          1293
      dlx=max(dlx,abs(d))                                                  1293
      r=r-d*v                                                              1293
11931 continue                                                             1294
      if(dlx.lt.shr)goto 11852                                             1295
11940 continue                                                             1295
11941 continue                                                             1295
      nlp=nlp+1                                                            1295
      dlx=0.0                                                              1296
11950 do 11951 l=1,nin                                                     1296
      k=m(l)                                                               1296
      bk=b(k)                                                              1296
      gk=dot_product(r,x(:,k))                                             1297
      u=gk+xv(k)*b(k)                                                      1297
      au=abs(u)-vp(k)*al1                                                  1298
      if(au .gt. 0.0)goto 11971                                            1298
      b(k)=0.0                                                             1298
      goto 11981                                                           1299
11971 continue                                                             1299
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1299
11981 continue                                                             1300
11961 continue                                                             1300
      d=b(k)-bk                                                            1300
      if(abs(d).le.0.0)goto 11951                                          1300
      dlx=max(dlx,abs(d))                                                  1301
      r=r-d*v*x(:,k)                                                       1302
11951 continue                                                             1303
11952 continue                                                             1303
      d=sum(r)/xmz                                                         1304
      if(d .eq. 0.0)goto 12001                                             1304
      b(0)=b(0)+d                                                          1304
      dlx=max(dlx,abs(d))                                                  1304
      r=r-d*v                                                              1304
12001 continue                                                             1306
      if(dlx.lt.shr)goto 11942                                             1306
      goto 11941                                                           1307
11942 continue                                                             1307
      goto 11851                                                           1308
11852 continue                                                             1308
      if(nin.gt.nx)goto 11842                                              1309
12010 do 12011 i=1,no                                                      1309
      fi=b(0)+g(i)                                                         1310
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1311
      if(fi .ge. fmin)goto 12031                                           1311
      q(i)=0.0                                                             1311
      goto 12021                                                           1311
12031 if(fi .le. fmax)goto 12041                                           1311
      q(i)=1.0                                                             1311
      goto 12051                                                           1312
12041 continue                                                             1312
      q(i)=1.0/(1.0+exp(-fi))                                              1312
12051 continue                                                             1313
12021 continue                                                             1313
12011 continue                                                             1314
12012 continue                                                             1314
      v=w*q*(1.0-q)                                                        1314
      xmz=sum(v)                                                           1314
      if(xmz.le.vmin)goto 11842                                            1315
      if(abs(b(0)-bs(0)) .ge. shr)goto 12071                               1315
      ix=0                                                                 1316
12080 do 12081 j=1,nin                                                     1316
      if(abs(b(m(j))-bs(m(j))).lt.shr)goto 12081                           1316
      ix=1                                                                 1316
      goto 12082                                                           1316
12081 continue                                                             1317
12082 continue                                                             1317
      if(ix.eq.0)goto 11842                                                1318
12071 continue                                                             1319
      r=w*(y-q)                                                            1320
      if(kopt .ne. 0)goto 12101                                            1321
12110 do 12111 j=1,nin                                                     1321
      xv(m(j))=dot_product(v,x(:,m(j))**2)                                 1321
12111 continue                                                             1322
12112 continue                                                             1322
12101 continue                                                             1323
      nit=nit+1                                                            1323
      if(nit .le. maxit)goto 12131                                         1323
      jerr=-ilm                                                            1323
      return                                                               1323
12131 continue                                                             1324
      goto 11841                                                           1325
11842 continue                                                             1325
      if(nin .le. nx)goto 12151                                            1325
      jerr=-10000-ilm                                                      1325
      goto 11772                                                           1325
12151 continue                                                             1326
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1326
      kin(ilm)=nin                                                         1327
      a0(ilm)=b(0)                                                         1327
      alm(ilm)=al                                                          1327
      lmu=ilm                                                              1328
      devi=dev2(no,w,y,q,pmin)                                             1329
      dev(ilm)=(dev1-devi)/dev0                                            1329
      if(xmz.le.vmin)goto 11772                                            1330
      if(ilm.lt.mnl)goto 11771                                             1330
      if(flmin.ge.1.0)goto 11771                                           1331
      me=0                                                                 1331
12160 do 12161 j=1,nin                                                     1331
      if(a(j,ilm).ne.0.0) me=me+1                                          1331
12161 continue                                                             1331
12162 continue                                                             1331
      if(me.gt.ne)goto 11772                                               1332
      if(dev(ilm).gt.devmax)goto 11772                                     1332
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 11772                             1333
11771 continue                                                             1334
11772 continue                                                             1334
      g=log(q/(1.0-q))                                                     1335
      deallocate(b,bs,v,r,xv,q,mm)                                         1336
      return                                                               1337
      end                                                                  1338
      function dev2(n,w,y,p,pmin)                                          1339
      real w(n),y(n),p(n)                                                  1340
      pmax=1.0-pmin                                                        1340
      s=0.0                                                                1341
12170 do 12171 i=1,n                                                       1341
      pi=min(max(pmin,p(i)),pmax)                                          1342
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1343
12171 continue                                                             1344
12172 continue                                                             1344
      dev2=s                                                               1345
      return                                                               1346
      end                                                                  1347
      function azero(n,y,g,q,jerr)                                         1348
      parameter(eps=1.0e-7)                                                1349
      real y(n),g(n),q(n)                                                  1350
      real, dimension (:), allocatable :: e,p,w                                 
      allocate(e(1:n),stat=jerr)                                           1354
      allocate(p(1:n),stat=ierr)                                           1354
      jerr=jerr+ierr                                                       1355
      allocate(w(1:n),stat=ierr)                                           1355
      jerr=jerr+ierr                                                       1356
      if(jerr.ne.0) return                                                 1357
      az=0.0                                                               1357
      e=exp(-g)                                                            1357
      qy=dot_product(q,y)                                                  1357
      p=1.0/(1.0+e)                                                        1358
12180 continue                                                             1358
12181 continue                                                             1358
      w=q*p*(1.0-p)                                                        1359
      d=(qy-dot_product(q,p))/sum(w)                                       1359
      az=az+d                                                              1359
      if(abs(d).lt.eps)goto 12182                                          1360
      ea0=exp(-az)                                                         1360
      p=1.0/(1.0+ea0*e)                                                    1361
      goto 12181                                                           1362
12182 continue                                                             1362
      azero=az                                                             1363
      deallocate(e,p,w)                                                    1364
      return                                                               1365
      end                                                                  1366
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ul   1368 
     *am,shr,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1370 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam)              1371
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1372
      integer ju(ni),m(nx),kin(nlam)                                       1373
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp                       
      real, dimension (:), allocatable :: di,v,r                                
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is                              
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no),stat=ierr)                                          1384
      jerr=jerr+ierr                                                       1385
      allocate(v(1:no),stat=ierr)                                          1385
      jerr=jerr+ierr                                                       1386
      allocate(mm(1:ni),stat=ierr)                                         1386
      jerr=jerr+ierr                                                       1387
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1387
      jerr=jerr+ierr                                                       1388
      allocate(sxp(1:no),stat=ierr)                                        1388
      jerr=jerr+ierr                                                       1389
      allocate(di(1:no),stat=ierr)                                         1389
      jerr=jerr+ierr                                                       1390
      if(jerr.ne.0) return                                                 1391
      pmax=1.0-pmin                                                        1391
      emin=pmin/pmax                                                       1391
      emax=1.0/emin                                                        1392
      pfm=(1.0+pmin)*pmin                                                  1392
      pfx=(1.0-pmin)*pmax                                                  1392
      vmin=pfm*pmax                                                        1393
      bta=max(parm,1.0e-3)                                                 1393
      omb=1.0-bta                                                          1393
      dev1=0.0                                                             1393
      dev0=0.0                                                             1394
12190 do 12191 ic=1,nc                                                     1394
      q0=dot_product(w,y(:,ic))                                            1395
      if(q0 .gt. pmin)goto 12211                                           1395
      jerr =8000+ic                                                        1395
      return                                                               1395
12211 continue                                                             1396
      if(q0 .lt. 1.0-pmin)goto 12231                                       1396
      jerr =9000+ic                                                        1396
      return                                                               1396
12231 continue                                                             1397
      b(0,ic)=log(q0)                                                      1397
      dev1=dev1-q0*b(0,ic)                                                 1397
      b(1:ni,ic)=0.0                                                       1398
12240 do 12241 i=1,no                                                      1398
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1398
12241 continue                                                             1399
12242 continue                                                             1399
12191 continue                                                             1400
12192 continue                                                             1400
      dev0=dev0+dev1                                                       1401
      if(nonzero(no*nc,g) .ne. 0)goto 12261                                1402
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1402
      sxp=0.0                                                              1403
12270 do 12271 ic=1,nc                                                     1403
      q(:,ic)=exp(b(0,ic))                                                 1403
      sxp=sxp+q(:,ic)                                                      1403
12271 continue                                                             1404
12272 continue                                                             1404
      goto 12281                                                           1405
12261 continue                                                             1405
12290 do 12291 i=1,no                                                      1405
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1405
12291 continue                                                             1405
12292 continue                                                             1405
      sxp=0.0                                                              1406
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1406
      if(jerr.ne.0) return                                                 1407
12300 do 12301 ic=1,nc                                                     1407
      q(:,ic)=exp(b(0,ic)+g(:,ic))                                         1407
      sxp=sxp+q(:,ic)                                                      1407
12301 continue                                                             1408
12302 continue                                                             1408
12281 continue                                                             1409
12251 continue                                                             1409
      if(isd .le. 0)goto 12321                                             1409
      xv=0.25                                                              1409
      goto 12331                                                           1410
12321 continue                                                             1410
12340 do 12341 j=1,ni                                                      1410
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1410
12341 continue                                                             1410
12342 continue                                                             1410
12331 continue                                                             1411
12311 continue                                                             1411
      if(flmin .ge. 1.0)goto 12361                                         1411
      eqs=max(eps,flmin)                                                   1411
      alf=eqs**(1.0/(nlam-1))                                              1411
12361 continue                                                             1412
      m=0                                                                  1412
      mm=0                                                                 1412
      nin=0                                                                1412
      nlp=0                                                                1412
      mnl=min(mnlam,nlam)                                                  1412
      bs=0.0                                                               1413
12370 do 12371 ilm=1,nlam                                                  1414
      if(flmin .lt. 1.0)goto 12391                                         1414
      al=ulam(ilm)                                                         1414
      goto 12381                                                           1415
12391 if(ilm .le. 2)goto 12401                                             1415
      al=al*alf                                                            1415
      goto 12381                                                           1416
12401 if(ilm .ne. 1)goto 12411                                             1416
      al=big                                                               1416
      goto 12421                                                           1417
12411 continue                                                             1417
      al=0.0                                                               1418
12430 do 12431 ic=1,nc                                                     1418
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1419
12440 do 12441 j=1,ni                                                      1419
      if(ju(j).eq.0)goto 12441                                             1419
      if(vp(j).le.0.0)goto 12441                                           1420
      al=max(al,abs(dot_product(r,x(:,j)))/vp(j))                          1421
12441 continue                                                             1422
12442 continue                                                             1422
12431 continue                                                             1423
12432 continue                                                             1423
      al=alf*al/bta                                                        1424
12421 continue                                                             1425
12381 continue                                                             1425
      al2=al*omb                                                           1425
      al1=al*bta                                                           1425
      nit=0                                                                1426
12450 continue                                                             1426
12451 continue                                                             1426
      ix=0                                                                 1426
      jx=ix                                                                1426
      ig=0                                                                 1427
12460 do 12461 ic=1,nc                                                     1427
      bs(0,ic)=b(0,ic)                                                     1428
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1429
      xmz=0.0                                                              1430
12470 do 12471 i=1,no                                                      1430
      pic=q(i,ic)/sxp(i)                                                   1431
      if(pic .ge. pfm)goto 12491                                           1431
      pic=0.0                                                              1431
      v(i)=0.0                                                             1431
      goto 12481                                                           1432
12491 if(pic .le. pfx)goto 12501                                           1432
      pic=1.0                                                              1432
      v(i)=0.0                                                             1432
      goto 12511                                                           1433
12501 continue                                                             1433
      v(i)=w(i)*pic*(1.0-pic)                                              1433
      xmz=xmz+v(i)                                                         1433
12511 continue                                                             1434
12481 continue                                                             1434
      r(i)=w(i)*(y(i,ic)-pic)                                              1435
12471 continue                                                             1436
12472 continue                                                             1436
      if(xmz.le.vmin)goto 12461                                            1436
      ig=1                                                                 1437
      if(kopt .ne. 0)goto 12531                                            1438
12540 do 12541 j=1,nin                                                     1438
      xv(m(j),ic)=dot_product(v,x(:,m(j))**2)                              1438
12541 continue                                                             1439
12542 continue                                                             1439
12531 continue                                                             1440
12550 continue                                                             1440
12551 continue                                                             1440
      nlp=nlp+1                                                            1440
      dlx=0.0                                                              1441
12560 do 12561 k=1,ni                                                      1441
      if(ju(k).eq.0)goto 12561                                             1442
      bk=b(k,ic)                                                           1442
      gk=dot_product(r,x(:,k))                                             1443
      u=gk+xv(k,ic)*b(k,ic)                                                1443
      au=abs(u)-vp(k)*al1                                                  1444
      if(au .gt. 0.0)goto 12581                                            1444
      b(k,ic)=0.0                                                          1444
      goto 12591                                                           1445
12581 continue                                                             1445
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1445
12591 continue                                                             1446
12571 continue                                                             1446
      d=b(k,ic)-bk                                                         1446
      if(abs(d).le.0.0)goto 12561                                          1446
      dlx=max(dlx,abs(d))                                                  1447
      r=r-d*v*x(:,k)                                                       1448
      if(mm(k) .ne. 0)goto 12611                                           1448
      nin=nin+1                                                            1449
      if(nin .le. nx)goto 12631                                            1449
      jx=1                                                                 1449
      goto 12562                                                           1449
12631 continue                                                             1450
      mm(k)=nin                                                            1450
      m(nin)=k                                                             1451
12611 continue                                                             1452
12561 continue                                                             1453
12562 continue                                                             1453
      if(jx.gt.0)goto 12552                                                1454
      d=sum(r)/xmz                                                         1455
      if(d .eq. 0.0)goto 12651                                             1455
      b(0,ic)=b(0,ic)+d                                                    1455
      dlx=max(dlx,abs(d))                                                  1455
      r=r-d*v                                                              1455
12651 continue                                                             1456
      if(dlx.lt.shr)goto 12552                                             1457
12660 continue                                                             1457
12661 continue                                                             1457
      nlp=nlp+1                                                            1457
      dlx=0.0                                                              1458
12670 do 12671 l=1,nin                                                     1458
      k=m(l)                                                               1458
      bk=b(k,ic)                                                           1459
      gk=dot_product(r,x(:,k))                                             1460
      u=gk+xv(k,ic)*b(k,ic)                                                1460
      au=abs(u)-vp(k)*al1                                                  1461
      if(au .gt. 0.0)goto 12691                                            1461
      b(k,ic)=0.0                                                          1461
      goto 12701                                                           1462
12691 continue                                                             1462
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1462
12701 continue                                                             1463
12681 continue                                                             1463
      d=b(k,ic)-bk                                                         1463
      if(abs(d).le.0.0)goto 12671                                          1464
      dlx=max(dlx,abs(d))                                                  1464
      r=r-d*v*x(:,k)                                                       1465
12671 continue                                                             1466
12672 continue                                                             1466
      d=sum(r)/xmz                                                         1467
      if(d .eq. 0.0)goto 12721                                             1467
      b(0,ic)=b(0,ic)+d                                                    1468
      dlx=max(dlx,abs(d))                                                  1468
      r=r-d*v                                                              1469
12721 continue                                                             1470
      if(dlx.lt.shr)goto 12662                                             1470
      goto 12661                                                           1471
12662 continue                                                             1471
      goto 12551                                                           1472
12552 continue                                                             1472
      if(jx.gt.0)goto 12462                                                1473
      if(abs(b(0,ic)-bs(0,ic)).gt.shr) ix=1                                1474
      if(ix .ne. 0)goto 12741                                              1475
12750 do 12751 j=1,nin                                                     1476
      if(abs(b(m(j),ic)-bs(m(j),ic)) .le. shr)goto 12771                   1476
      ix=1                                                                 1476
      goto 12752                                                           1476
12771 continue                                                             1477
12751 continue                                                             1478
12752 continue                                                             1478
12741 continue                                                             1479
12780 do 12781 i=1,no                                                      1479
      fi=b(0,ic)+g(i,ic)                                                   1481
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1482
      fi=min(max(exmn,fi),exmx)                                            1482
      sxp(i)=sxp(i)-q(i,ic)                                                1483
      q(i,ic)=min(max(emin*sxp(i),exp(dble(fi))),emax*sxp(i))              1484
      sxp(i)=sxp(i)+q(i,ic)                                                1485
12781 continue                                                             1486
12782 continue                                                             1486
12461 continue                                                             1487
12462 continue                                                             1487
      s=-sum(b(0,:))/nc                                                    1487
      b(0,:)=b(0,:)+s                                                      1487
      di=s                                                                 1488
12790 do 12791 j=1,nin                                                     1488
      l=m(j)                                                               1489
      if(vp(l) .gt. 0.0)goto 12811                                         1489
      s=sum(b(l,:))/nc                                                     1489
      goto 12821                                                           1490
12811 continue                                                             1490
      s=elc(parm,nc,b(l,:),is)                                             1490
12821 continue                                                             1491
12801 continue                                                             1491
      b(l,:)=b(l,:)-s                                                      1491
      di=di-s*x(:,l)                                                       1492
12791 continue                                                             1493
12792 continue                                                             1493
      di=exp(di)                                                           1493
      sxp=sxp*di                                                           1493
12830 do 12831 ic=1,nc                                                     1493
      q(:,ic)=q(:,ic)*di                                                   1493
12831 continue                                                             1494
12832 continue                                                             1494
      if(jx.gt.0)goto 12452                                                1494
      if(ix.eq.0)goto 12452                                                1494
      if(ig.eq.0)goto 12452                                                1495
      nit=nit+1                                                            1495
      if(nit .le. maxit)goto 12851                                         1495
      jerr=-ilm                                                            1495
      return                                                               1495
12851 continue                                                             1496
      goto 12451                                                           1497
12452 continue                                                             1497
      if(jx .le. 0)goto 12871                                              1497
      jerr=-10000-ilm                                                      1497
      goto 12372                                                           1497
12871 continue                                                             1497
      devi=0.0                                                             1498
12880 do 12881 ic=1,nc                                                     1499
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1499
      a0(ic,ilm)=b(0,ic)                                                   1500
12890 do 12891 i=1,no                                                      1500
      if(y(i,ic).le.0.0)goto 12891                                         1501
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1502
12891 continue                                                             1503
12892 continue                                                             1503
12881 continue                                                             1504
12882 continue                                                             1504
      kin(ilm)=nin                                                         1504
      alm(ilm)=al                                                          1504
      lmu=ilm                                                              1505
      dev(ilm)=(dev1-devi)/dev0                                            1505
      if(ig.eq.0)goto 12372                                                1506
      if(ilm.lt.mnl)goto 12371                                             1506
      if(flmin.ge.1.0)goto 12371                                           1507
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 12372             1508
      if(dev(ilm).gt.devmax)goto 12372                                     1508
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12372                             1509
12371 continue                                                             1510
12372 continue                                                             1510
      g=log(q)                                                             1510
12900 do 12901 i=1,no                                                      1510
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1510
12901 continue                                                             1511
12902 continue                                                             1511
      deallocate(sxp,b,bs,v,r,xv,q,mm,is)                                  1512
      return                                                               1513
      end                                                                  1514
      subroutine kazero(kk,n,y,g,q,az,jerr)                                1515
      parameter(eps=1.0e-7)                                                1516
      real y(n,kk),g(n,kk),q(n),az(kk)                                     1517
      real, dimension (:), allocatable :: s                                     
      real, dimension (:,:), allocatable :: e                                   
      allocate(e(1:n,1:kk),stat=jerr)                                           
      allocate(s(1:n),stat=ierr)                                           1522
      jerr=jerr+ierr                                                       1523
      if(jerr.ne.0) return                                                 1524
      az=0.0                                                               1524
      e=exp(g)                                                             1524
12910 do 12911 i=1,n                                                       1524
      s(i)=sum(e(i,:))                                                     1524
12911 continue                                                             1525
12912 continue                                                             1525
12920 continue                                                             1525
12921 continue                                                             1525
      dm=0.0                                                               1526
12930 do 12931 k=1,kk                                                      1526
      t=0.0                                                                1526
      u=t                                                                  1527
12940 do 12941 i=1,n                                                       1527
      pik=e(i,k)/s(i)                                                      1528
      t=t+q(i)*(y(i,k)-pik)                                                1528
      u=u+q(i)*pik*(1.0-pik)                                               1529
12941 continue                                                             1530
12942 continue                                                             1530
      d=t/u                                                                1530
      az(k)=az(k)+d                                                        1530
      ed=exp(d)                                                            1530
      dm=max(dm,abs(d))                                                    1531
12950 do 12951 i=1,n                                                       1531
      z=e(i,k)                                                             1531
      e(i,k)=z*ed                                                          1531
      s(i)=s(i)-z+e(i,k)                                                   1531
12951 continue                                                             1532
12952 continue                                                             1532
12931 continue                                                             1533
12932 continue                                                             1533
      if(dm.lt.eps)goto 12922                                              1533
      goto 12921                                                           1534
12922 continue                                                             1534
      az=az-sum(az)/kk                                                     1535
      deallocate(e,s)                                                      1536
      return                                                               1537
      end                                                                  1538
      function elc(parm,n,a,m)                                             1539
      real a(n)                                                            1539
      integer m(n)                                                         1540
      fn=n                                                                 1540
      am=sum(a)/fn                                                         1541
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 12971                       1541
      elc=am                                                               1541
      return                                                               1541
12971 continue                                                             1542
12980 do 12981 i=1,n                                                       1542
      m(i)=i                                                               1542
12981 continue                                                             1542
12982 continue                                                             1542
      call psort7(a,m,1,n)                                                 1543
      if(a(m(1)) .ne. a(m(n)))goto 13001                                   1543
      elc=a(1)                                                             1543
      return                                                               1543
13001 continue                                                             1544
      if(mod(n,2) .ne. 1)goto 13021                                        1544
      ad=a(m(n/2+1))                                                       1544
      goto 13031                                                           1545
13021 continue                                                             1545
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       1545
13031 continue                                                             1546
13011 continue                                                             1546
      if(parm .ne. 1.0)goto 13051                                          1546
      elc=ad                                                               1546
      return                                                               1546
13051 continue                                                             1547
      b1=min(am,ad)                                                        1547
      b2=max(am,ad)                                                        1547
      k2=1                                                                 1548
13060 continue                                                             1548
13061 if(a(m(k2)).gt.b1)goto 13062                                         1548
      k2=k2+1                                                              1548
      goto 13061                                                           1548
13062 continue                                                             1548
      k1=k2-1                                                              1549
13070 continue                                                             1549
13071 if(a(m(k2)).ge.b2)goto 13072                                         1549
      k2=k2+1                                                              1549
      goto 13071                                                           1550
13072 continue                                                             1550
      r=parm/((1.0-parm)*fn)                                               1550
      is=0                                                                 1550
      sm=n-2*(k1-1)                                                        1551
13080 do 13081 k=k1,k2-1                                                   1551
      sm=sm-2.0                                                            1551
      s=r*sm+am                                                            1552
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 13101                   1552
      is=k                                                                 1552
      goto 13082                                                           1552
13101 continue                                                             1553
13081 continue                                                             1554
13082 continue                                                             1554
      if(is .eq. 0)goto 13121                                              1554
      elc=s                                                                1554
      return                                                               1554
13121 continue                                                             1554
      r2=2.0*r                                                             1554
      s1=a(m(k1))                                                          1554
      am2=2.0*am                                                           1555
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    1555
      elc=s1                                                               1556
13130 do 13131 k=k1+1,k2                                                   1556
      s=a(m(k))                                                            1556
      if(s.eq.s1)goto 13131                                                1557
      c=r2*sum(abs(a-s))+s*(s-am2)                                         1558
      if(c .ge. cri)goto 13151                                             1558
      cri=c                                                                1558
      elc=s                                                                1558
13151 continue                                                             1558
      s1=s                                                                 1559
13131 continue                                                             1560
13132 continue                                                             1560
      return                                                               1561
      end                                                                  1562
      function nintot(ni,nx,nc,a,m,nin,is)                                 1563
      real a(nx,nc)                                                        1563
      integer m(nx),is(ni)                                                 1564
      is=0                                                                 1564
      nintot=0                                                             1565
13160 do 13161 ic=1,nc                                                     1565
13170 do 13171 j=1,nin                                                     1565
      k=m(j)                                                               1565
      if(is(k).ne.0)goto 13171                                             1566
      if(a(j,ic).eq.0.0)goto 13171                                         1566
      is(k)=k                                                              1566
      nintot=nintot+1                                                      1567
13171 continue                                                             1567
13172 continue                                                             1567
13161 continue                                                             1568
13162 continue                                                             1568
      return                                                               1569
      end                                                                  1570
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             1571
      real ca(nx,nc),a(ni,nc)                                              1571
      integer ia(nx)                                                       1572
      a=0.0                                                                1573
13180 do 13181 ic=1,nc                                                     1573
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            1573
13181 continue                                                             1574
13182 continue                                                             1574
      return                                                               1575
      end                                                                  1576
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      1577
      real a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                             1577
      integer ia(nx)                                                       1578
13190 do 13191 i=1,nt                                                      1578
13200 do 13201 ic=1,nc                                                     1578
      ans(ic,i)=a0(ic)                                                     1580
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   1581 
     *:nin)))
13201 continue                                                             1581
13202 continue                                                             1581
13191 continue                                                             1582
13192 continue                                                             1582
      return                                                               1583
      end                                                                  1584
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,ne,nx,nlam,fl   1586 
     *min,  ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      real x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)                 1587
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1588
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1589
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 13221                                    1593
      jerr=10000                                                           1593
      return                                                               1593
13221 continue                                                             1594
      allocate(ww(1:no),stat=jerr)                                         1595
      allocate(ju(1:ni),stat=ierr)                                         1595
      jerr=jerr+ierr                                                       1596
      allocate(vq(1:ni),stat=ierr)                                         1596
      jerr=jerr+ierr                                                       1597
      allocate(xm(1:ni),stat=ierr)                                         1597
      jerr=jerr+ierr                                                       1598
      allocate(xs(1:ni),stat=ierr)                                         1598
      jerr=jerr+ierr                                                       1599
      if(jerr.ne.0) return                                                 1600
      call spchkvars(no,ni,x,ix,ju)                                        1601
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1602
      if(maxval(ju) .gt. 0)goto 13241                                      1602
      jerr=7777                                                            1602
      return                                                               1602
13241 continue                                                             1603
      vq=max(0.0,vp)                                                       1603
      vq=vq*ni/sum(vq)                                                     1604
13250 do 13251 i=1,no                                                      1604
      ww(i)=sum(y(i,:))                                                    1604
      y(i,:)=y(i,:)/ww(i)                                                  1604
13251 continue                                                             1604
13252 continue                                                             1604
      sw=sum(ww)                                                           1604
      ww=ww/sw                                                             1605
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     1606
      if(nc .ne. 1)goto 13271                                              1607
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,ne,nx,n   1609 
     *lam,flmin,  ulam,thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,de
     *v,alm,nlp,jerr)
      goto 13281                                                           1610
13271 continue                                                             1611
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmi   1613 
     *n,ulam,thr, isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp
     *,jerr)
13281 continue                                                             1614
13261 continue                                                             1614
      if(jerr.gt.0) return                                                 1614
      dev0=2.0*sw*dev0                                                     1615
13290 do 13291 k=1,lmu                                                     1615
      nk=nin(k)                                                            1616
13300 do 13301 ic=1,nc                                                     1616
      if(isd .le. 0)goto 13321                                             1616
13330 do 13331 l=1,nk                                                      1616
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1616
13331 continue                                                             1616
13332 continue                                                             1616
13321 continue                                                             1617
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1618
13301 continue                                                             1619
13302 continue                                                             1619
13291 continue                                                             1620
13292 continue                                                             1620
      deallocate(ww,ju,vq,xm,xs)                                           1621
      return                                                               1622
      end                                                                  1623
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,xm,xs)                1624
      real x(*),w(no),xm(ni),xs(ni)                                        1624
      integer ix(*),jx(*),ju(ni)                                           1625
13340 do 13341 j=1,ni                                                      1625
      if(ju(j).eq.0)goto 13341                                             1625
      jb=ix(j)                                                             1625
      je=ix(j+1)-1                                                         1626
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1627
      if(isd.gt.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   1628 
     *)**2)
13341 continue                                                             1629
13342 continue                                                             1629
      if(isd.eq.0) xs=1.0                                                  1630
      return                                                               1631
      end                                                                  1632
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam,   1634 
     *  flmin,ulam,shr,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,alm,
     *nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1636 
     *5, devmax=0.999)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        1637
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1638
      real xb(ni),xs(ni)                                                   1638
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1639
      real, dimension (:), allocatable :: xm,b,bs,v,r,sc,xv,q                   
      integer, dimension (:), allocatable :: mm                                 
      allocate(b(0:ni),stat=jerr)                                          1644
      allocate(xm(0:ni),stat=ierr)                                         1644
      jerr=jerr+ierr                                                       1645
      allocate(xv(1:ni),stat=ierr)                                         1645
      jerr=jerr+ierr                                                       1646
      allocate(bs(0:ni),stat=ierr)                                         1646
      jerr=jerr+ierr                                                       1647
      allocate(mm(1:ni),stat=ierr)                                         1647
      jerr=jerr+ierr                                                       1648
      allocate(q(1:no),stat=ierr)                                          1648
      jerr=jerr+ierr                                                       1649
      allocate(r(1:no),stat=ierr)                                          1649
      jerr=jerr+ierr                                                       1650
      allocate(v(1:no),stat=ierr)                                          1650
      jerr=jerr+ierr                                                       1651
      allocate(sc(1:no),stat=ierr)                                         1651
      jerr=jerr+ierr                                                       1652
      if(jerr.ne.0) return                                                 1653
      fmax=log(1.0/pmin-1.0)                                               1653
      fmin=-fmax                                                           1653
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1654
      bta=max(parm,1.0e-3)                                                 1654
      omb=1.0-bta                                                          1655
      q0=dot_product(w,y)                                                  1655
      if(q0 .gt. pmin)goto 13361                                           1655
      jerr=8001                                                            1655
      return                                                               1655
13361 continue                                                             1656
      if(q0 .lt. 1.0-pmin)goto 13381                                       1656
      jerr=9001                                                            1656
      return                                                               1656
13381 continue                                                             1656
      bz=log(q0/(1.0-q0))                                                  1657
      if(nonzero(no,g) .ne. 0)goto 13401                                   1657
      vi=q0*(1.0-q0)                                                       1657
      b(0)=bz                                                              1657
      v=vi*w                                                               1658
      r=w*(y-q0)                                                           1658
      q=q0                                                                 1658
      xm(0)=vi                                                             1659
      goto 13411                                                           1660
13401 continue                                                             1660
      b(0)=azero(no,y,g,w,jerr)                                            1660
      if(jerr.ne.0) return                                                 1661
      q=1.0/(1.0+exp(-b(0)-g))                                             1661
      v=w*q*(1.0-q)                                                        1661
      r=w*(y-q)                                                            1661
      xm(0)=sum(v)                                                         1662
13411 continue                                                             1663
13391 continue                                                             1663
      if(isd .le. 0)goto 13431                                             1663
      xv=0.25                                                              1663
      goto 13441                                                           1664
13431 continue                                                             1665
13450 do 13451 j=1,ni                                                      1665
      if(ju(j).eq.0)goto 13451                                             1665
      jb=ix(j)                                                             1665
      je=ix(j+1)-1                                                         1666
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          1667
13451 continue                                                             1668
13452 continue                                                             1668
13441 continue                                                             1669
13421 continue                                                             1669
      b(1:ni)=0.0                                                          1669
      dev1=-(bz*q0+log(1.0-q0))                                            1669
      dev0=dev1                                                            1670
13460 do 13461 i=1,no                                                      1670
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1671
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1672
13461 continue                                                             1673
13462 continue                                                             1673
      if(flmin .ge. 1.0)goto 13481                                         1673
      eqs=max(eps,flmin)                                                   1673
      alf=eqs**(1.0/(nlam-1))                                              1673
13481 continue                                                             1674
      m=0                                                                  1674
      mm=0                                                                 1674
      nin=0                                                                1674
      o=0.0                                                                1674
      svr=o                                                                1674
      mnl=min(mnlam,nlam)                                                  1674
      bs=0.0                                                               1674
      nlp=0                                                                1674
      nin=nlp                                                              1675
13490 do 13491 ilm=1,nlam                                                  1676
      if(flmin .lt. 1.0)goto 13511                                         1676
      al=ulam(ilm)                                                         1676
      goto 13501                                                           1677
13511 if(ilm .le. 2)goto 13521                                             1677
      al=al*alf                                                            1677
      goto 13501                                                           1678
13521 if(ilm .ne. 1)goto 13531                                             1678
      al=big                                                               1678
      goto 13541                                                           1679
13531 continue                                                             1679
      al=0.0                                                               1680
13550 do 13551 j=1,ni                                                      1680
      if(ju(j).eq.0)goto 13551                                             1680
      if(vp(j).le.0.0)goto 13551                                           1681
      jb=ix(j)                                                             1681
      je=ix(j+1)-1                                                         1681
      jn=ix(j+1)-ix(j)                                                     1682
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1683
      gj=dot_product(sc(1:jn),x(jb:je))                                    1684
      gj=(gj-svr*xb(j))/xs(j)                                              1685
      al=max(al,abs(gj)/vp(j))                                             1686
13551 continue                                                             1687
13552 continue                                                             1687
      al=alf*al/bta                                                        1688
13541 continue                                                             1689
13501 continue                                                             1689
      al2=al*omb                                                           1689
      al1=al*bta                                                           1689
      nit=0                                                                1690
13560 continue                                                             1690
13561 continue                                                             1690
      bs(0)=b(0)                                                           1690
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1691
13570 continue                                                             1691
13571 continue                                                             1691
      nlp=nlp+1                                                            1691
      dlx=0.0                                                              1692
13580 do 13581 k=1,ni                                                      1692
      if(ju(k).eq.0)goto 13581                                             1693
      jb=ix(k)                                                             1693
      je=ix(k+1)-1                                                         1693
      jn=ix(k+1)-ix(k)                                                     1693
      bk=b(k)                                                              1694
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1695
      gk=dot_product(sc(1:jn),x(jb:je))                                    1696
      gk=(gk-svr*xb(k))/xs(k)                                              1697
      u=gk+xv(k)*b(k)                                                      1697
      au=abs(u)-vp(k)*al1                                                  1698
      if(au .gt. 0.0)goto 13601                                            1698
      b(k)=0.0                                                             1698
      goto 13611                                                           1699
13601 continue                                                             1699
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1699
13611 continue                                                             1700
13591 continue                                                             1700
      d=b(k)-bk                                                            1700
      if(abs(d).le.0.0)goto 13581                                          1700
      dlx=max(dlx,abs(d))                                                  1701
      if(mm(k) .ne. 0)goto 13631                                           1701
      nin=nin+1                                                            1701
      if(nin.gt.nx)goto 13582                                              1702
      mm(k)=nin                                                            1702
      m(nin)=k                                                             1702
      sc(1:jn)=v(jx(jb:je))                                                1703
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 1704
13631 continue                                                             1705
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1706
      o=o+d*(xb(k)/xs(k))                                                  1707
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1708
13581 continue                                                             1709
13582 continue                                                             1709
      if(nin.gt.nx)goto 13572                                              1710
      d=svr/xm(0)                                                          1711
      if(d .eq. 0.0)goto 13651                                             1711
      b(0)=b(0)+d                                                          1711
      dlx=max(dlx,abs(d))                                                  1711
      r=r-d*v                                                              1711
13651 continue                                                             1712
      svr=svr-d*xm(0)                                                      1712
      if(dlx.lt.shr)goto 13572                                             1713
13660 continue                                                             1713
13661 continue                                                             1713
      nlp=nlp+1                                                            1713
      dlx=0.0                                                              1714
13670 do 13671 l=1,nin                                                     1714
      k=m(l)                                                               1714
      jb=ix(k)                                                             1714
      je=ix(k+1)-1                                                         1715
      jn=ix(k+1)-ix(k)                                                     1715
      bk=b(k)                                                              1716
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1717
      gk=dot_product(sc(1:jn),x(jb:je))                                    1718
      gk=(gk-svr*xb(k))/xs(k)                                              1719
      u=gk+xv(k)*b(k)                                                      1719
      au=abs(u)-vp(k)*al1                                                  1720
      if(au .gt. 0.0)goto 13691                                            1720
      b(k)=0.0                                                             1720
      goto 13701                                                           1721
13691 continue                                                             1721
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1721
13701 continue                                                             1722
13681 continue                                                             1722
      d=b(k)-bk                                                            1722
      if(abs(d).le.0.0)goto 13671                                          1722
      dlx=max(dlx,abs(d))                                                  1723
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1724
      o=o+d*(xb(k)/xs(k))                                                  1725
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1726
13671 continue                                                             1727
13672 continue                                                             1727
      d=svr/xm(0)                                                          1728
      if(d .eq. 0.0)goto 13721                                             1728
      b(0)=b(0)+d                                                          1728
      dlx=max(dlx,abs(d))                                                  1728
      r=r-d*v                                                              1728
13721 continue                                                             1729
      svr=svr-d*xm(0)                                                      1730
      if(dlx.lt.shr)goto 13662                                             1730
      goto 13661                                                           1731
13662 continue                                                             1731
      goto 13571                                                           1732
13572 continue                                                             1732
      if(nin.gt.nx)goto 13562                                              1733
      sc=b(0)                                                              1733
      b0=0.0                                                               1734
13730 do 13731 j=1,nin                                                     1734
      l=m(j)                                                               1734
      jb=ix(l)                                                             1734
      je=ix(l+1)-1                                                         1735
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      1736
      b0=b0-b(l)*xb(l)/xs(l)                                               1737
13731 continue                                                             1738
13732 continue                                                             1738
      sc=sc+b0                                                             1739
13740 do 13741 i=1,no                                                      1739
      fi=sc(i)+g(i)                                                        1740
      if(fi .ge. fmin)goto 13761                                           1740
      q(i)=0.0                                                             1740
      goto 13751                                                           1740
13761 if(fi .le. fmax)goto 13771                                           1740
      q(i)=1.0                                                             1740
      goto 13781                                                           1741
13771 continue                                                             1741
      q(i)=1.0/(1.0+exp(-fi))                                              1741
13781 continue                                                             1742
13751 continue                                                             1742
13741 continue                                                             1743
13742 continue                                                             1743
      v=w*q*(1.0-q)                                                        1743
      xm(0)=sum(v)                                                         1743
      if(xm(0).lt.vmin)goto 13562                                          1744
      if(abs(b(0)-bs(0)) .ge. shr)goto 13801                               1744
      kx=0                                                                 1745
13810 do 13811 j=1,nin                                                     1745
      if(abs(b(m(j))-bs(m(j))).lt.shr)goto 13811                           1745
      kx=1                                                                 1745
      goto 13812                                                           1745
13811 continue                                                             1746
13812 continue                                                             1746
      if(kx.eq.0)goto 13562                                                1747
13801 continue                                                             1748
      r=w*(y-q)                                                            1748
      svr=sum(r)                                                           1748
      o=0.0                                                                1749
13820 do 13821 l=1,nin                                                     1749
      j=m(l)                                                               1750
      jb=ix(j)                                                             1750
      je=ix(j+1)-1                                                         1750
      jn=ix(j+1)-ix(j)                                                     1751
      sc(1:jn)=v(jx(jb:je))                                                1752
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 1753
      if(kopt .ne. 0)goto 13841                                            1754
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              1755
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                1757
13841 continue                                                             1758
13821 continue                                                             1759
13822 continue                                                             1759
      nit=nit+1                                                            1759
      if(nit .le. maxit)goto 13861                                         1759
      jerr=-ilm                                                            1759
      return                                                               1759
13861 continue                                                             1760
      goto 13561                                                           1761
13562 continue                                                             1761
      if(nin .le. nx)goto 13881                                            1761
      jerr=-10000-ilm                                                      1761
      goto 13492                                                           1761
13881 continue                                                             1762
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1762
      kin(ilm)=nin                                                         1763
      a0(ilm)=b(0)                                                         1763
      alm(ilm)=al                                                          1763
      lmu=ilm                                                              1764
      devi=dev2(no,w,y,q,pmin)                                             1765
      dev(ilm)=(dev1-devi)/dev0                                            1766
      if(ilm.lt.mnl)goto 13491                                             1766
      if(flmin.ge.1.0)goto 13491                                           1767
      me=0                                                                 1767
13890 do 13891 j=1,nin                                                     1767
      if(a(j,ilm).ne.0.0) me=me+1                                          1767
13891 continue                                                             1767
13892 continue                                                             1767
      if(me.gt.ne)goto 13492                                               1768
      if(dev(ilm).gt.devmax)goto 13492                                     1768
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13492                             1769
      if(xm(0).lt.vmin)goto 13492                                          1770
13491 continue                                                             1771
13492 continue                                                             1771
      g=log(q/(1.0-q))                                                     1772
      deallocate(xm,b,bs,v,r,sc,xv,q,mm)                                   1773
      return                                                               1774
      end                                                                  1775
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam   1777 
     *,flmin,ulam,  shr,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,alm
     *,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1779 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni)    1780
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1781
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1782
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp                       
      real, dimension (:), allocatable :: sc,xm,v,r                             
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is                              
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(xm(0:ni),stat=ierr)                                         1793
      jerr=jerr+ierr                                                       1794
      allocate(r(1:no),stat=ierr)                                          1794
      jerr=jerr+ierr                                                       1795
      allocate(v(1:no),stat=ierr)                                          1795
      jerr=jerr+ierr                                                       1796
      allocate(mm(1:ni),stat=ierr)                                         1796
      jerr=jerr+ierr                                                       1797
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1797
      jerr=jerr+ierr                                                       1798
      allocate(sxp(1:no),stat=ierr)                                        1798
      jerr=jerr+ierr                                                       1799
      allocate(sc(1:no),stat=ierr)                                         1799
      jerr=jerr+ierr                                                       1800
      if(jerr.ne.0) return                                                 1801
      pmax=1.0-pmin                                                        1801
      emin=pmin/pmax                                                       1801
      emax=1.0/emin                                                        1802
      pfm=(1.0+pmin)*pmin                                                  1802
      pfx=(1.0-pmin)*pmax                                                  1802
      vmin=pfm*pmax                                                        1803
      bta=max(parm,1.0e-3)                                                 1803
      omb=1.0-bta                                                          1803
      dev1=0.0                                                             1803
      dev0=0.0                                                             1804
13900 do 13901 ic=1,nc                                                     1804
      q0=dot_product(w,y(:,ic))                                            1805
      if(q0 .gt. pmin)goto 13921                                           1805
      jerr =8000+ic                                                        1805
      return                                                               1805
13921 continue                                                             1806
      if(q0 .lt. 1.0-pmin)goto 13941                                       1806
      jerr =9000+ic                                                        1806
      return                                                               1806
13941 continue                                                             1807
      b(1:ni,ic)=0.0                                                       1807
      b(0,ic)=log(q0)                                                      1807
      dev1=dev1-q0*b(0,ic)                                                 1808
13950 do 13951 i=1,no                                                      1808
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1808
13951 continue                                                             1809
13952 continue                                                             1809
13901 continue                                                             1810
13902 continue                                                             1810
      dev0=dev0+dev1                                                       1811
      if(nonzero(no*nc,g) .ne. 0)goto 13971                                1812
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1812
      sxp=0.0                                                              1813
13980 do 13981 ic=1,nc                                                     1813
      q(:,ic)=exp(b(0,ic))                                                 1813
      sxp=sxp+q(:,ic)                                                      1813
13981 continue                                                             1814
13982 continue                                                             1814
      goto 13991                                                           1815
13971 continue                                                             1815
14000 do 14001 i=1,no                                                      1815
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1815
14001 continue                                                             1815
14002 continue                                                             1815
      sxp=0.0                                                              1816
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1816
      if(jerr.ne.0) return                                                 1817
14010 do 14011 ic=1,nc                                                     1817
      q(:,ic)=exp(b(0,ic)+g(:,ic))                                         1817
      sxp=sxp+q(:,ic)                                                      1817
14011 continue                                                             1818
14012 continue                                                             1818
13991 continue                                                             1819
13961 continue                                                             1819
      if(isd .le. 0)goto 14031                                             1819
      xv=0.25                                                              1819
      goto 14041                                                           1820
14031 continue                                                             1821
14050 do 14051 j=1,ni                                                      1821
      if(ju(j).eq.0)goto 14051                                             1821
      jb=ix(j)                                                             1821
      je=ix(j+1)-1                                                         1822
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        1823
14051 continue                                                             1824
14052 continue                                                             1824
14041 continue                                                             1825
14021 continue                                                             1825
      if(flmin .ge. 1.0)goto 14071                                         1825
      eqs=max(eps,flmin)                                                   1825
      alf=eqs**(1.0/(nlam-1))                                              1825
14071 continue                                                             1826
      m=0                                                                  1826
      mm=0                                                                 1826
      nin=0                                                                1826
      nlp=0                                                                1826
      mnl=min(mnlam,nlam)                                                  1826
      bs=0.0                                                               1826
      svr=0.0                                                              1826
      o=0.0                                                                1827
14080 do 14081 ilm=1,nlam                                                  1828
      if(flmin .lt. 1.0)goto 14101                                         1828
      al=ulam(ilm)                                                         1828
      goto 14091                                                           1829
14101 if(ilm .le. 2)goto 14111                                             1829
      al=al*alf                                                            1829
      goto 14091                                                           1830
14111 if(ilm .ne. 1)goto 14121                                             1830
      al=big                                                               1830
      goto 14131                                                           1831
14121 continue                                                             1831
      al=0.0                                                               1832
14140 do 14141 ic=1,nc                                                     1832
      v=q(:,ic)/sxp                                                        1832
      r=w*(y(:,ic)-v)                                                      1832
      v=w*v*(1.0-v)                                                        1833
14150 do 14151 j=1,ni                                                      1833
      if(ju(j).eq.0)goto 14151                                             1833
      if(vp(j).le.0.0)goto 14151                                           1834
      jb=ix(j)                                                             1834
      je=ix(j+1)-1                                                         1834
      jn=ix(j+1)-ix(j)                                                     1835
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1836
      gj=dot_product(sc(1:jn),x(jb:je))                                    1837
      gj=(gj-svr*xb(j))/xs(j)                                              1838
      al=max(al,abs(gj)/vp(j))                                             1839
14151 continue                                                             1840
14152 continue                                                             1840
14141 continue                                                             1841
14142 continue                                                             1841
      al=alf*al/bta                                                        1842
14131 continue                                                             1843
14091 continue                                                             1843
      al2=al*omb                                                           1843
      al1=al*bta                                                           1843
      nit=0                                                                1844
14160 continue                                                             1844
14161 continue                                                             1844
      ixx=0                                                                1844
      jxx=ixx                                                              1844
      ig=0                                                                 1845
14170 do 14171 ic=1,nc                                                     1845
      bs(0,ic)=b(0,ic)                                                     1846
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1847
      xm(0)=0.0                                                            1847
      svr=0.0                                                              1847
      o=0.0                                                                1848
14180 do 14181 i=1,no                                                      1848
      pic=q(i,ic)/sxp(i)                                                   1849
      if(pic .ge. pfm)goto 14201                                           1849
      pic=0.0                                                              1849
      v(i)=0.0                                                             1849
      goto 14191                                                           1850
14201 if(pic .le. pfx)goto 14211                                           1850
      pic=1.0                                                              1850
      v(i)=0.0                                                             1850
      goto 14221                                                           1851
14211 continue                                                             1851
      v(i)=w(i)*pic*(1.0-pic)                                              1851
      xm(0)=xm(0)+v(i)                                                     1851
14221 continue                                                             1852
14191 continue                                                             1852
      r(i)=w(i)*(y(i,ic)-pic)                                              1852
      svr=svr+r(i)                                                         1853
14181 continue                                                             1854
14182 continue                                                             1854
      if(xm(0).le.vmin)goto 14171                                          1854
      ig=1                                                                 1855
14230 do 14231 l=1,nin                                                     1855
      j=m(l)                                                               1856
      jb=ix(j)                                                             1856
      je=ix(j+1)-1                                                         1857
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             1858
      if(kopt .ne. 0)goto 14251                                            1859
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       1860
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          1861
14251 continue                                                             1862
14231 continue                                                             1863
14232 continue                                                             1863
14260 continue                                                             1863
14261 continue                                                             1863
      nlp=nlp+1                                                            1863
      dlx=0.0                                                              1864
14270 do 14271 k=1,ni                                                      1864
      if(ju(k).eq.0)goto 14271                                             1865
      jb=ix(k)                                                             1865
      je=ix(k+1)-1                                                         1865
      jn=ix(k+1)-ix(k)                                                     1865
      bk=b(k,ic)                                                           1866
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1867
      gk=dot_product(sc(1:jn),x(jb:je))                                    1868
      gk=(gk-svr*xb(k))/xs(k)                                              1869
      u=gk+xv(k,ic)*b(k,ic)                                                1869
      au=abs(u)-vp(k)*al1                                                  1870
      if(au .gt. 0.0)goto 14291                                            1870
      b(k,ic)=0.0                                                          1870
      goto 14301                                                           1871
14291 continue                                                             1871
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1871
14301 continue                                                             1872
14281 continue                                                             1872
      d=b(k,ic)-bk                                                         1872
      if(abs(d).le.0.0)goto 14271                                          1872
      dlx=max(dlx,abs(d))                                                  1873
      if(mm(k) .ne. 0)goto 14321                                           1873
      nin=nin+1                                                            1874
      if(nin .le. nx)goto 14341                                            1874
      jxx=1                                                                1874
      goto 14272                                                           1874
14341 continue                                                             1875
      mm(k)=nin                                                            1875
      m(nin)=k                                                             1876
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             1877
14321 continue                                                             1878
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1879
      o=o+d*(xb(k)/xs(k))                                                  1880
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1881
14271 continue                                                             1882
14272 continue                                                             1882
      if(jxx.gt.0)goto 14262                                               1883
      d=svr/xm(0)                                                          1884
      if(d .eq. 0.0)goto 14361                                             1884
      b(0,ic)=b(0,ic)+d                                                    1884
      dlx=max(dlx,abs(d))                                                  1885
      r=r-d*v                                                              1885
      svr=svr-d*xm(0)                                                      1886
14361 continue                                                             1887
      if(dlx.lt.shr)goto 14262                                             1888
14370 continue                                                             1888
14371 continue                                                             1888
      nlp=nlp+1                                                            1888
      dlx=0.0                                                              1889
14380 do 14381 l=1,nin                                                     1889
      k=m(l)                                                               1889
      jb=ix(k)                                                             1889
      je=ix(k+1)-1                                                         1890
      jn=ix(k+1)-ix(k)                                                     1890
      bk=b(k,ic)                                                           1891
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1892
      gk=dot_product(sc(1:jn),x(jb:je))                                    1893
      gk=(gk-svr*xb(k))/xs(k)                                              1894
      u=gk+xv(k,ic)*b(k,ic)                                                1894
      au=abs(u)-vp(k)*al1                                                  1895
      if(au .gt. 0.0)goto 14401                                            1895
      b(k,ic)=0.0                                                          1895
      goto 14411                                                           1896
14401 continue                                                             1896
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1896
14411 continue                                                             1897
14391 continue                                                             1897
      d=b(k,ic)-bk                                                         1897
      if(abs(d).le.0.0)goto 14381                                          1898
      dlx=max(dlx,abs(d))                                                  1899
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1900
      o=o+d*(xb(k)/xs(k))                                                  1901
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1902
14381 continue                                                             1903
14382 continue                                                             1903
      d=svr/xm(0)                                                          1904
      if(d .eq. 0.0)goto 14431                                             1904
      b(0,ic)=b(0,ic)+d                                                    1904
      dlx=max(dlx,abs(d))                                                  1905
      r=r-d*v                                                              1905
      svr=svr-d*xm(0)                                                      1906
14431 continue                                                             1907
      if(dlx.lt.shr)goto 14372                                             1907
      goto 14371                                                           1908
14372 continue                                                             1908
      goto 14261                                                           1909
14262 continue                                                             1909
      if(jxx.gt.0)goto 14172                                               1910
      if(abs(b(0,ic)-bs(0,ic)).gt.shr) ixx=1                               1911
      if(ixx .ne. 0)goto 14451                                             1912
14460 do 14461 j=1,nin                                                     1913
      if(abs(b(m(j),ic)-bs(m(j),ic)) .le. shr)goto 14481                   1913
      ixx=1                                                                1913
      goto 14462                                                           1913
14481 continue                                                             1914
14461 continue                                                             1915
14462 continue                                                             1915
14451 continue                                                             1916
      sc=b(0,ic)+g(:,ic)                                                   1916
      b0=0.0                                                               1917
14490 do 14491 j=1,nin                                                     1917
      l=m(j)                                                               1917
      jb=ix(l)                                                             1917
      je=ix(l+1)-1                                                         1918
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   1919
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            1920
14491 continue                                                             1921
14492 continue                                                             1921
      sc=min(max(exmn,sc+b0),exmx)                                         1922
      sxp=sxp-q(:,ic)                                                      1923
      q(:,ic)=min(max(emin*sxp,exp(dble(sc))),emax*sxp)                    1924
      sxp=sxp+q(:,ic)                                                      1925
14171 continue                                                             1926
14172 continue                                                             1926
      s=-sum(b(0,:))/nc                                                    1926
      b(0,:)=b(0,:)+s                                                      1926
      sc=s                                                                 1926
      b0=0.0                                                               1927
14500 do 14501 j=1,nin                                                     1927
      l=m(j)                                                               1928
      if(vp(l) .gt. 0.0)goto 14521                                         1928
      s=sum(b(l,:))/nc                                                     1928
      goto 14531                                                           1929
14521 continue                                                             1929
      s=elc(parm,nc,b(l,:),is)                                             1929
14531 continue                                                             1930
14511 continue                                                             1930
      b(l,:)=b(l,:)-s                                                      1931
      jb=ix(l)                                                             1931
      je=ix(l+1)-1                                                         1932
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         1933
      b0=b0+s*xb(l)/xs(l)                                                  1934
14501 continue                                                             1935
14502 continue                                                             1935
      sc=sc+b0                                                             1935
      sc=exp(sc)                                                           1935
      sxp=sxp*sc                                                           1935
14540 do 14541 ic=1,nc                                                     1935
      q(:,ic)=q(:,ic)*sc                                                   1935
14541 continue                                                             1936
14542 continue                                                             1936
      if(jxx.gt.0)goto 14162                                               1936
      if(ixx.eq.0)goto 14162                                               1936
      if(ig.eq.0)goto 14162                                                1937
      nit=nit+1                                                            1937
      if(nit .le. maxit)goto 14561                                         1937
      jerr=-ilm                                                            1937
      return                                                               1937
14561 continue                                                             1938
      goto 14161                                                           1939
14162 continue                                                             1939
      if(jxx .le. 0)goto 14581                                             1939
      jerr=-10000-ilm                                                      1939
      goto 14082                                                           1939
14581 continue                                                             1939
      devi=0.0                                                             1940
14590 do 14591 ic=1,nc                                                     1941
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1941
      a0(ic,ilm)=b(0,ic)                                                   1942
14600 do 14601 i=1,no                                                      1942
      if(y(i,ic).le.0.0)goto 14601                                         1943
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1944
14601 continue                                                             1945
14602 continue                                                             1945
14591 continue                                                             1946
14592 continue                                                             1946
      kin(ilm)=nin                                                         1946
      alm(ilm)=al                                                          1946
      lmu=ilm                                                              1947
      dev(ilm)=(dev1-devi)/dev0                                            1947
      if(ig.eq.0)goto 14082                                                1948
      if(ilm.lt.mnl)goto 14081                                             1948
      if(flmin.ge.1.0)goto 14081                                           1949
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 14082             1950
      if(dev(ilm).gt.devmax)goto 14082                                     1950
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 14082                             1951
14081 continue                                                             1952
14082 continue                                                             1952
      g=log(q)                                                             1952
14610 do 14611 i=1,no                                                      1952
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1952
14611 continue                                                             1953
14612 continue                                                             1953
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc)                            1954
      return                                                               1955
      end                                                                  1956
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  1957
      real a0(nc),ca(nx,nc),x(*),f(nc,n)                                   1957
      integer ia(*),ix(*),jx(*)                                            1958
14620 do 14621 ic=1,nc                                                     1958
      f(ic,:)=a0(ic)                                                       1958
14621 continue                                                             1959
14622 continue                                                             1959
14630 do 14631 j=1,nin                                                     1959
      k=ia(j)                                                              1959
      kb=ix(k)                                                             1959
      ke=ix(k+1)-1                                                         1960
14640 do 14641 ic=1,nc                                                     1960
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    1960
14641 continue                                                             1961
14642 continue                                                             1961
14631 continue                                                             1962
14632 continue                                                             1962
      return                                                               1963
      end                                                                  1964
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,ne,nx,nlam,flmin,ula   1966 
     *m,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam)              1967
      real ca(nx,nlam),dev(nlam),alm(nlam)                                 1968
      integer jd(*),ia(nx),nin(nlam)                                       1969
      real, dimension (:), allocatable :: xs,ww,vq                              
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 14661                                    1973
      jerr=10000                                                           1973
      return                                                               1973
14661 continue                                                             1974
      allocate(ww(1:no),stat=jerr)                                         1975
      allocate(ju(1:ni),stat=ierr)                                         1975
      jerr=jerr+ierr                                                       1976
      allocate(vq(1:ni),stat=ierr)                                         1976
      jerr=jerr+ierr                                                       1977
      if(isd .le. 0)goto 14681                                             1977
      allocate(xs(1:ni),stat=ierr)                                         1977
      jerr=jerr+ierr                                                       1977
14681 continue                                                             1978
      if(jerr.ne.0) return                                                 1979
      call chkvars(no,ni,x,ju)                                             1980
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1981
      if(maxval(ju) .gt. 0)goto 14701                                      1981
      jerr=7777                                                            1981
      return                                                               1981
14701 continue                                                             1982
      vq=max(0.0,vp)                                                       1982
      vq=vq*ni/sum(vq)                                                     1983
      ww=max(0.0,w)                                                        1983
      sw=sum(ww)                                                           1984
      if(sw .gt. 0.0)goto 14721                                            1984
      jerr=9999                                                            1984
      return                                                               1984
14721 continue                                                             1984
      ww=ww/sw                                                             1985
      call cstandard(no,ni,x,ww,ju,isd,xs)                                 1986
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr   1988 
     *,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1988
      dev0=2.0*sw*dev0                                                     1989
      if(isd .le. 0)goto 14741                                             1989
14750 do 14751 k=1,lmu                                                     1989
      nk=nin(k)                                                            1989
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   1989
14751 continue                                                             1989
14752 continue                                                             1989
14741 continue                                                             1990
      deallocate(ww,ju,vq)                                                 1990
      if(isd.gt.0) deallocate(xs)                                          1991
      return                                                               1992
      end                                                                  1993
      subroutine cstandard (no,ni,x,w,ju,isd,xs)                           1994
      real x(no,ni),w(no),xs(ni)                                           1994
      integer ju(ni)                                                       1995
14760 do 14761 j=1,ni                                                      1995
      if(ju(j).eq.0)goto 14761                                             1996
      xm=dot_product(w,x(:,j))                                             1996
      x(:,j)=x(:,j)-xm                                                     1997
      if(isd .le. 0)goto 14781                                             1998
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 1998
      x(:,j)=x(:,j)/xs(j)                                                  1999
14781 continue                                                             2000
14761 continue                                                             2001
14762 continue                                                             2001
      return                                                               2002
      end                                                                  2003
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,ne,nx,nlam,flmin,ula   2005 
     *m,cthr,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=0.001, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99)   2006
      real x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam)              2007
      real ao(nx,nlam),dev(nlam),alm(nlam)                                 2008
      integer ju(ni),m(nx),kin(nlam)                                       2009
      real, dimension (:), allocatable :: w,dk,v,xs,wr,a,as,f,dq                
      real, dimension (:), allocatable :: e,uu                                  
      integer, dimension (:), allocatable :: jp,kp,mm                           
      allocate(e(1:no),stat=jerr)                                          2015
      allocate(uu(1:no),stat=ierr)                                         2015
      jerr=jerr+ierr                                                       2016
      allocate(f(1:no),stat=ierr)                                          2016
      jerr=jerr+ierr                                                       2017
      allocate(w(1:no),stat=ierr)                                          2017
      jerr=jerr+ierr                                                       2018
      allocate(v(1:ni),stat=ierr)                                          2018
      jerr=jerr+ierr                                                       2019
      allocate(a(1:ni),stat=ierr)                                          2019
      jerr=jerr+ierr                                                       2020
      allocate(as(1:ni),stat=ierr)                                         2020
      jerr=jerr+ierr                                                       2021
      allocate(xs(1:ni),stat=ierr)                                         2021
      jerr=jerr+ierr                                                       2022
      allocate(jp(1:no),stat=ierr)                                         2022
      jerr=jerr+ierr                                                       2023
      allocate(kp(1:no),stat=ierr)                                         2023
      jerr=jerr+ierr                                                       2024
      allocate(dk(1:no),stat=ierr)                                         2024
      jerr=jerr+ierr                                                       2025
      allocate(wr(1:no),stat=ierr)                                         2025
      jerr=jerr+ierr                                                       2026
      allocate(dq(1:no),stat=ierr)                                         2026
      jerr=jerr+ierr                                                       2027
      allocate(mm(1:ni),stat=ierr)                                         2027
      jerr=jerr+ierr                                                       2028
      if(jerr.ne.0)go to 11420                                             2029
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2030
      if(jerr.ne.0) go to 11420                                            2030
      alpha=max(parm,1.0e-3)                                               2031
      oma=1.0-alpha                                                        2031
      nlm=0                                                                2032
      dq=d*q                                                               2032
      call died(no,nk,dq,kp,jp,dk)                                         2033
      a=0.0                                                                2033
      f=0.0                                                                2033
      e=q                                                                  2033
      fmax=log(huge(f(1))*0.1)                                             2034
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 2035
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         2035
      dev0=rr                                                              2036
14790 do 14791 i=1,no                                                      2036
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 14811                   2036
      w(i)=0.0                                                             2036
      wr(i)=w(i)                                                           2036
14811 continue                                                             2036
14791 continue                                                             2037
14792 continue                                                             2037
      if(nonzero(no,g) .eq. 0)goto 14831                                   2037
      f=g-dot_product(q,g)                                                 2038
      e=q*exp(sign(min(abs(f),fmax),f))                                    2039
14831 continue                                                             2040
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2041
      if(jerr.ne.0) go to 11420                                            2042
      call vars(no,ni,x,w,v)                                               2043
      if(flmin .ge. 1.0)goto 14851                                         2043
      eqs=max(eps,flmin)                                                   2043
      alf=eqs**(1.0/(nlam-1))                                              2043
14851 continue                                                             2044
      m=0                                                                  2044
      mm=0                                                                 2044
      nlp=0                                                                2044
      nin=nlp                                                              2044
      mnl=min(mnlam,nlam)                                                  2044
      as=0.0                                                               2045
14860 do 14861 ilm=1,nlam                                                  2046
      if(flmin .lt. 1.0)goto 14881                                         2046
      al=ulam(ilm)                                                         2046
      goto 14871                                                           2047
14881 if(ilm .le. 2)goto 14891                                             2047
      al=al*alf                                                            2047
      goto 14871                                                           2048
14891 if(ilm .ne. 1)goto 14901                                             2048
      al=big                                                               2048
      goto 14911                                                           2049
14901 continue                                                             2049
      al=0.0                                                               2050
14920 do 14921 j=1,ni                                                      2050
      if(ju(j).eq.0)goto 14921                                             2050
      if(vp(j).le.0.0)goto 14921                                           2051
      al=max(al,abs(dot_product(wr,x(:,j)))/vp(j))                         2052
14921 continue                                                             2053
14922 continue                                                             2053
      al=alf*al/alpha                                                      2054
14911 continue                                                             2055
14871 continue                                                             2055
      sa=alpha*al                                                          2055
      omal=oma*al                                                          2055
      nit=0                                                                2056
14930 do 14931 ito=1,maxit                                                 2056
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2057
14940 do 14941 iti=1,maxit                                                 2057
      nlp=nlp+1                                                            2057
      dli=0.0                                                              2058
14950 do 14951 j=1,ni                                                      2058
      if(ju(j).eq.0)goto 14951                                             2059
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2060
      if(abs(u) .gt. vp(j)*sa)goto 14971                                   2060
      at=0.0                                                               2060
      goto 14981                                                           2061
14971 continue                                                             2061
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2061
14981 continue                                                             2062
14961 continue                                                             2062
      if(at .eq. a(j))goto 15001                                           2062
      del=at-a(j)                                                          2062
      a(j)=at                                                              2062
      dli=max(dli,abs(del))                                                2063
      wr=wr-del*w*x(:,j)                                                   2063
      f=f+del*x(:,j)                                                       2064
      if(mm(j) .ne. 0)goto 15021                                           2064
      nin=nin+1                                                            2064
      if(nin.gt.nx)goto 14952                                              2065
      mm(j)=nin                                                            2065
      m(nin)=j                                                             2066
15021 continue                                                             2067
15001 continue                                                             2068
14951 continue                                                             2069
14952 continue                                                             2069
      if(nin.gt.nx)goto 14942                                              2069
      if(dli.lt.cthr)goto 14942                                            2070
15030 do 15031 ita=1,maxit                                                 2070
      nlp=nlp+1                                                            2070
      dli=0.0                                                              2071
15040 do 15041 l=1,nin                                                     2071
      j=m(l)                                                               2072
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2073
      if(abs(u) .gt. vp(j)*sa)goto 15061                                   2073
      at=0.0                                                               2073
      goto 15071                                                           2074
15061 continue                                                             2074
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2074
15071 continue                                                             2075
15051 continue                                                             2075
      if(at .eq. a(j))goto 15091                                           2075
      del=at-a(j)                                                          2075
      a(j)=at                                                              2075
      dli=max(dli,abs(del))                                                2076
      wr=wr-del*w*x(:,j)                                                   2076
      f=f+del*x(:,j)                                                       2077
15091 continue                                                             2078
15041 continue                                                             2079
15042 continue                                                             2079
      if(dli.lt.cthr)goto 15032                                            2080
15031 continue                                                             2081
15032 continue                                                             2081
      if(dli .lt. cthr)goto 15111                                          2081
      jerr=-1                                                              2081
      go to 11420                                                          2081
15111 continue                                                             2082
14941 continue                                                             2083
14942 continue                                                             2083
      if(nin.gt.nx)goto 14932                                              2083
      if(dli .lt. cthr)goto 15131                                          2083
      jerr=-1                                                              2083
      go to 11420                                                          2083
15131 continue                                                             2084
      e=q*exp(sign(min(abs(f),fmax),f))                                    2085
      ix=0                                                                 2085
15140 do 15141 j=1,nin                                                     2085
      if(abs(a(m(j))-as(m(j))).lt.cthr)goto 15141                          2085
      ix=1                                                                 2085
      goto 15142                                                           2085
15141 continue                                                             2086
15142 continue                                                             2086
      if(ix.eq.0)goto 14932                                                2087
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2088
      if(jerr.ne.0) go to 11420                                            2089
      call vars(no,ni,x,w,v)                                               2090
14931 continue                                                             2091
14932 continue                                                             2091
      if(ix .eq. 0)goto 15161                                              2091
      jerr=-2                                                              2091
      goto 14862                                                           2091
15161 continue                                                             2092
      if(nin .le. nx)goto 15181                                            2092
      jerr=-10000-ilm                                                      2092
      goto 14862                                                           2092
15181 continue                                                             2093
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               2093
      kin(ilm)=nin                                                         2094
      alm(ilm)=al                                                          2094
      lmu=ilm                                                              2095
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   2096
      if(ilm.lt.mnl)goto 14861                                             2096
      if(flmin.ge.1.0)goto 14861                                           2097
      me=0                                                                 2097
15190 do 15191 j=1,nin                                                     2097
      if(ao(j,ilm).ne.0.0) me=me+1                                         2097
15191 continue                                                             2097
15192 continue                                                             2097
      if(me.gt.ne)goto 14862                                               2098
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 14862              2099
      if(dev(ilm).gt.devmax)goto 14862                                     2100
14861 continue                                                             2101
14862 continue                                                             2101
      g=f                                                                  2102
11420 continue                                                             2102
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm)                     2103
      return                                                               2104
      end                                                                  2105
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 2106
      real ca(nin),x(n,*),f(n)                                             2106
      integer ia(nin)                                                      2107
      f=0.0                                                                2107
      if(nin.le.0) return                                                  2108
15200 do 15201 i=1,n                                                       2108
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      2108
15201 continue                                                             2109
15202 continue                                                             2109
      return                                                               2110
      end                                                                  2111
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         2112
      real y(no),d(no),q(no)                                               2112
      integer jp(no),kp(*)                                                 2113
15210 do 15211 j=1,no                                                      2113
      jp(j)=j                                                              2113
15211 continue                                                             2113
15212 continue                                                             2113
      call psort7(y,jp,1,no)                                               2114
      nj=0                                                                 2114
15220 do 15221 j=1,no                                                      2114
      if(q(jp(j)).le.0.0)goto 15221                                        2114
      nj=nj+1                                                              2114
      jp(nj)=jp(j)                                                         2114
15221 continue                                                             2115
15222 continue                                                             2115
      if(nj .ne. 0)goto 15241                                              2115
      jerr=-4                                                              2115
      return                                                               2115
15241 continue                                                             2116
      j=1                                                                  2116
15250 continue                                                             2116
15251 if(d(jp(j)).gt.0.0)goto 15252                                        2116
      j=j+1                                                                2116
      if(j.gt.nj)goto 15252                                                2116
      goto 15251                                                           2117
15252 continue                                                             2117
      if(j .lt. nj-1)goto 15271                                            2117
      jerr=-5                                                              2117
      return                                                               2117
15271 continue                                                             2118
      j0=j-1                                                               2118
      nj=nj-j0                                                             2118
15280 do 15281 j=1,nj                                                      2118
      jp(j)=jp(j+j0)                                                       2118
15281 continue                                                             2119
15282 continue                                                             2119
      jerr=0                                                               2119
      nk=0                                                                 2119
      t0=y(jp(1))                                                          2119
      yk=t0                                                                2119
      j=2                                                                  2120
15290 continue                                                             2120
15291 continue                                                             2120
15300 continue                                                             2121
15301 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 15302                     2121
      j=j+1                                                                2121
      if(j.gt.nj)goto 15302                                                2121
      goto 15301                                                           2122
15302 continue                                                             2122
      nk=nk+1                                                              2122
      kp(nk)=j-1                                                           2122
      if(j.gt.nj)goto 15292                                                2123
      if(j .ne. nj)goto 15321                                              2123
      nk=nk+1                                                              2123
      kp(nk)=nj                                                            2123
      goto 15292                                                           2123
15321 continue                                                             2124
      yk=y(jp(j))                                                          2124
      j=j+1                                                                2125
      goto 15291                                                           2126
15292 continue                                                             2126
      return                                                               2127
      end                                                                  2128
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     2129
      real d(no),dk(nk),wr(no),w(no)                                       2130
      real e(no),u(no),b,c                                                 2130
      integer kp(nk),jp(no)                                                2131
      call usk(no,nk,kp,jp,e,u)                                            2132
      b=dk(1)/u(1)                                                         2132
      c=dk(1)/u(1)**2                                                      2132
      jerr=0                                                               2133
15330 do 15331 j=1,kp(1)                                                   2133
      i=jp(j)                                                              2134
      w(i)=e(i)*(b-e(i)*c)                                                 2134
      if(w(i) .gt. 0.0)goto 15351                                          2134
      jerr=-3                                                              2134
      return                                                               2134
15351 continue                                                             2135
      wr(i)=d(i)-e(i)*b                                                    2136
15331 continue                                                             2137
15332 continue                                                             2137
15360 do 15361 k=2,nk                                                      2137
      j1=kp(k-1)+1                                                         2137
      j2=kp(k)                                                             2138
      b=b+dk(k)/u(k)                                                       2138
      c=c+dk(k)/u(k)**2                                                    2139
15370 do 15371 j=j1,j2                                                     2139
      i=jp(j)                                                              2140
      w(i)=e(i)*(b-e(i)*c)                                                 2140
      if(w(i) .gt. 0.0)goto 15391                                          2140
      jerr=-3                                                              2140
      return                                                               2140
15391 continue                                                             2141
      wr(i)=d(i)-e(i)*b                                                    2142
15371 continue                                                             2143
15372 continue                                                             2143
15361 continue                                                             2144
15362 continue                                                             2144
      return                                                               2145
      end                                                                  2146
      subroutine vars(no,ni,x,w,v)                                         2147
      real x(no,ni),w(no),v(ni)                                            2148
15400 do 15401 j=1,ni                                                      2148
      v(j)=dot_product(w,x(:,j)**2)                                        2148
15401 continue                                                             2149
15402 continue                                                             2149
      return                                                               2150
      end                                                                  2151
      subroutine died(no,nk,d,kp,jp,dk)                                    2152
      real d(no),dk(nk)                                                    2152
      integer kp(nk),jp(no)                                                2153
      dk(1)=sum(d(jp(1:kp(1))))                                            2154
15410 do 15411 k=2,nk                                                      2154
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  2154
15411 continue                                                             2155
15412 continue                                                             2155
      return                                                               2156
      end                                                                  2157
      subroutine usk(no,nk,kp,jp,e,u)                                      2158
      real e(no),u(nk),h                                                   2158
      integer kp(nk),jp(no)                                                2159
      h=0.0                                                                2160
15420 do 15421 k=nk,1,-1                                                   2160
      j2=kp(k)                                                             2161
      j1=1                                                                 2161
      if(k.gt.1) j1=kp(k-1)+1                                              2162
15430 do 15431 j=j2,j1,-1                                                  2162
      h=h+e(jp(j))                                                         2162
15431 continue                                                             2163
15432 continue                                                             2163
      u(k)=h                                                               2164
15421 continue                                                             2165
15422 continue                                                             2165
      return                                                               2166
      end                                                                  2167
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             2168
      real d(no),dk(nk),f(no)                                              2169
      integer kp(nk),jp(no)                                                2169
      real e(no),u(nk),s                                                   2170
      call usk(no,nk,kp,jp,e,u)                                            2170
      u=log(u)                                                             2171
      risk=dot_product(d,f)-dot_product(dk,u)                              2172
      return                                                               2173
      end                                                                  2174
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 2175
      real x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(nlam)          2176
      real, dimension (:), allocatable :: dk,f,xm,dq,q                          
      real, dimension (:), allocatable :: e,uu                                  
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          2182
      allocate(q(1:no),stat=ierr)                                          2182
      jerr=jerr+ierr                                                       2183
      allocate(uu(1:no),stat=ierr)                                         2183
      jerr=jerr+ierr                                                       2184
      allocate(f(1:no),stat=ierr)                                          2184
      jerr=jerr+ierr                                                       2185
      allocate(dk(1:no),stat=ierr)                                         2185
      jerr=jerr+ierr                                                       2186
      allocate(jp(1:no),stat=ierr)                                         2186
      jerr=jerr+ierr                                                       2187
      allocate(kp(1:no),stat=ierr)                                         2187
      jerr=jerr+ierr                                                       2188
      allocate(dq(1:no),stat=ierr)                                         2188
      jerr=jerr+ierr                                                       2189
      allocate(xm(1:ni),stat=ierr)                                         2189
      jerr=jerr+ierr                                                       2190
      if(jerr.ne.0) go to 11420                                            2191
      q=max(0.0,w)                                                         2191
      sw=sum(q)                                                            2192
      if(sw .gt. 0.0)goto 15451                                            2192
      jerr=9999                                                            2192
      go to 11420                                                          2192
15451 continue                                                             2193
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2194
      if(jerr.ne.0) go to 11420                                            2194
      fmax=log(huge(e(1))*0.1)                                             2195
      dq=d*q                                                               2195
      call died(no,nk,dq,kp,jp,dk)                                         2195
      gm=dot_product(q,g)/sw                                               2196
15460 do 15461 j=1,ni                                                      2196
      xm(j)=dot_product(q,x(:,j))/sw                                       2196
15461 continue                                                             2197
15462 continue                                                             2197
15470 do 15471 lam=1,nlam                                                  2198
15480 do 15481 i=1,no                                                      2198
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       2199
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        2200
15481 continue                                                             2201
15482 continue                                                             2201
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          2202
15471 continue                                                             2203
15472 continue                                                             2203
11420 continue                                                             2203
      deallocate(e,uu,dk,f,jp,kp,dq)                                       2204
      return                                                               2205
      end                                                                  2206
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,ne,nx,nlam,flmin,ulam   2208 
     *,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    2209
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2210
      integer jd(*),ia(nx),nin(nlam)                                       2211
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 15501                                    2215
      jerr=10000                                                           2215
      return                                                               2215
15501 continue                                                             2216
      if(minval(y) .ge. 0.0)goto 15521                                     2216
      jerr=8888                                                            2216
      return                                                               2216
15521 continue                                                             2217
      allocate(ww(1:no),stat=jerr)                                         2218
      allocate(ju(1:ni),stat=ierr)                                         2218
      jerr=jerr+ierr                                                       2219
      allocate(vq(1:ni),stat=ierr)                                         2219
      jerr=jerr+ierr                                                       2220
      allocate(xm(1:ni),stat=ierr)                                         2220
      jerr=jerr+ierr                                                       2221
      if(isd .le. 0)goto 15541                                             2221
      allocate(xs(1:ni),stat=ierr)                                         2221
      jerr=jerr+ierr                                                       2221
15541 continue                                                             2222
      if(jerr.ne.0) return                                                 2223
      call chkvars(no,ni,x,ju)                                             2224
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2225
      if(maxval(ju) .gt. 0)goto 15561                                      2225
      jerr=7777                                                            2225
      go to 11420                                                          2225
15561 continue                                                             2226
      vq=max(0.0,vp)                                                       2226
      vq=vq*ni/sum(vq)                                                     2227
      ww=max(0.0,w)                                                        2227
      sw=sum(ww)                                                           2227
      if(sw .gt. 0.0)goto 15581                                            2227
      jerr=9999                                                            2227
      go to 11420                                                          2227
15581 continue                                                             2228
      ww=ww/sw                                                             2229
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             2230
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr,   2232 
     *  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11420                                            2232
      dev0=2.0*sw*dev0                                                     2233
15590 do 15591 k=1,lmu                                                     2233
      nk=nin(k)                                                            2234
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2235
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2236
15591 continue                                                             2237
15592 continue                                                             2237
11420 continue                                                             2237
      deallocate(ww,ju,vq,xm)                                              2237
      if(isd.gt.0) deallocate(xs)                                          2238
      return                                                               2239
      end                                                                  2240
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,ne,nx,nlam,flmin,ulam   2242 
     *,shr,  isd,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2243 
     *9)
      real x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)                    2244
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2245
      integer ju(ni),m(nx),kin(nlam)                                       2246
      real, dimension (:), allocatable :: t,w,wr,v,a,f,as                       
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                          2251
      allocate(as(1:ni),stat=ierr)                                         2251
      jerr=jerr+ierr                                                       2252
      allocate(t(1:no),stat=ierr)                                          2252
      jerr=jerr+ierr                                                       2253
      allocate(mm(1:ni),stat=ierr)                                         2253
      jerr=jerr+ierr                                                       2254
      allocate(wr(1:no),stat=ierr)                                         2254
      jerr=jerr+ierr                                                       2255
      allocate(v(1:ni),stat=ierr)                                          2255
      jerr=jerr+ierr                                                       2256
      allocate(w(1:no),stat=ierr)                                          2256
      jerr=jerr+ierr                                                       2257
      allocate(f(1:no),stat=ierr)                                          2257
      jerr=jerr+ierr                                                       2258
      if(jerr.ne.0) return                                                 2259
      bta=max(parm,1.0e-3)                                                 2259
      omb=1.0-bta                                                          2260
      t=q*y                                                                2260
      yb=sum(t)                                                            2260
      fmax=log(huge(bta)*0.1)                                              2261
      if(nonzero(no,g) .ne. 0)goto 15611                                   2261
      w=q*yb                                                               2261
      az=log(yb)                                                           2261
      f=az                                                                 2261
      goto 15621                                                           2262
15611 continue                                                             2262
      w=q*exp(sign(min(abs(g),fmax),g))                                    2262
      v0=sum(w)                                                            2262
      eaz=yb/v0                                                            2263
      w=eaz*w                                                              2263
      az=log(eaz)                                                          2263
      f=az+g                                                               2264
15621 continue                                                             2265
15601 continue                                                             2265
      a=0.0                                                                2265
      wr=t-w                                                               2265
      v0=yb                                                                2265
      dv0=yb*(log(yb)-1.0)                                                 2265
      dvr=-yb                                                              2266
15630 do 15631 i=1,no                                                      2266
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               2266
15631 continue                                                             2266
15632 continue                                                             2266
      dvr=dvr-dv0                                                          2266
      dev0=dvr                                                             2267
15640 do 15641 j=1,ni                                                      2267
      if(ju(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                         2267
15641 continue                                                             2268
15642 continue                                                             2268
      if(flmin .ge. 1.0)goto 15661                                         2268
      eqs=max(eps,flmin)                                                   2268
      alf=eqs**(1.0/(nlam-1))                                              2268
15661 continue                                                             2269
      m=0                                                                  2269
      mm=0                                                                 2269
      nlp=0                                                                2269
      nin=nlp                                                              2269
      mnl=min(mnlam,nlam)                                                  2270
15670 do 15671 ilm=1,nlam                                                  2271
      if(flmin .lt. 1.0)goto 15691                                         2271
      al=ulam(ilm)                                                         2271
      goto 15681                                                           2272
15691 if(ilm .le. 2)goto 15701                                             2272
      al=al*alf                                                            2272
      goto 15681                                                           2273
15701 if(ilm .ne. 1)goto 15711                                             2273
      al=big                                                               2273
      goto 15721                                                           2274
15711 continue                                                             2274
      al=0.0                                                               2275
15730 do 15731 j=1,ni                                                      2275
      if(ju(j).eq.0)goto 15731                                             2275
      if(vp(j).le.0.0)goto 15731                                           2276
      al=max(al,abs(dot_product(wr,x(:,j)))/vp(j))                         2277
15731 continue                                                             2278
15732 continue                                                             2278
      al=alf*al/bta                                                        2279
15721 continue                                                             2280
15681 continue                                                             2280
      al2=al*omb                                                           2280
      al1=al*bta                                                           2280
      nit=0                                                                2281
15740 continue                                                             2281
15741 continue                                                             2281
      nit=nit+1                                                            2281
      if(nit .le. maxit)goto 15761                                         2281
      jerr=-2                                                              2281
      go to 11420                                                          2281
15761 continue                                                             2282
      az0=az                                                               2282
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2282
      nit1=0                                                               2283
15770 continue                                                             2283
15771 continue                                                             2283
      nit1=nit1+1                                                          2283
      if(nit1 .le. maxit)goto 15791                                        2283
      jerr=-1                                                              2283
      go to 11420                                                          2283
15791 continue                                                             2284
      nlp=nlp+1                                                            2284
      dlx=0.0                                                              2285
15800 do 15801 k=1,ni                                                      2285
      if(ju(k).eq.0)goto 15801                                             2285
      ak=a(k)                                                              2286
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2286
      au=abs(u)-vp(k)*al1                                                  2287
      if(au .gt. 0.0)goto 15821                                            2287
      a(k)=0.0                                                             2287
      goto 15831                                                           2288
15821 continue                                                             2288
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2288
15831 continue                                                             2289
15811 continue                                                             2289
      if(a(k).eq.ak)goto 15801                                             2289
      d=a(k)-ak                                                            2289
      dlx=max(dlx,abs(d))                                                  2290
      wr=wr-d*w*x(:,k)                                                     2290
      f=f+d*x(:,k)                                                         2291
      if(mm(k) .ne. 0)goto 15851                                           2291
      nin=nin+1                                                            2291
      if(nin.gt.nx)goto 15802                                              2292
      mm(k)=nin                                                            2292
      m(nin)=k                                                             2293
15851 continue                                                             2294
15801 continue                                                             2295
15802 continue                                                             2295
      if(nin.gt.nx)goto 15772                                              2295
      d=sum(wr)/v0                                                         2296
      az=az+d                                                              2296
      dlx=max(dlx,abs(d))                                                  2296
      wr=wr-d*w                                                            2296
      f=f+d                                                                2297
      if(dlx.lt.shr)goto 15772                                             2297
      nit2=0                                                               2298
15860 continue                                                             2298
15861 continue                                                             2298
      nit2=nit2+1                                                          2298
      if(nit2 .le. maxit)goto 15881                                        2298
      jerr=-1                                                              2298
      go to 11420                                                          2298
15881 continue                                                             2299
      nlp=nlp+1                                                            2299
      dlx=0.0                                                              2300
15890 do 15891 l=1,nin                                                     2300
      k=m(l)                                                               2300
      ak=a(k)                                                              2301
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2301
      au=abs(u)-vp(k)*al1                                                  2302
      if(au .gt. 0.0)goto 15911                                            2302
      a(k)=0.0                                                             2302
      goto 15921                                                           2303
15911 continue                                                             2303
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2303
15921 continue                                                             2304
15901 continue                                                             2304
      if(a(k).eq.ak)goto 15891                                             2304
      d=a(k)-ak                                                            2304
      dlx=max(dlx,abs(d))                                                  2305
      wr=wr-d*w*x(:,k)                                                     2305
      f=f+d*x(:,k)                                                         2307
15891 continue                                                             2307
15892 continue                                                             2307
      d=sum(wr)/v0                                                         2307
      az=az+d                                                              2307
      dlx=max(dlx,abs(d))                                                  2307
      wr=wr-d*w                                                            2307
      f=f+d                                                                2308
      if(dlx.lt.shr)goto 15862                                             2308
      goto 15861                                                           2309
15862 continue                                                             2309
      goto 15771                                                           2310
15772 continue                                                             2310
      if(nin.gt.nx)goto 15742                                              2311
      w=q*exp(sign(min(abs(f),fmax),f))                                    2311
      v0=sum(w)                                                            2312
      if(abs(az-az0) .ge. shr)goto 15941                                   2312
      ix=0                                                                 2313
15950 do 15951 j=1,nin                                                     2313
      if(abs(a(m(j))-as(m(j))).lt.shr)goto 15951                           2313
      ix=1                                                                 2313
      goto 15952                                                           2313
15951 continue                                                             2314
15952 continue                                                             2314
      if(ix.eq.0)goto 15742                                                2315
15941 continue                                                             2316
      wr=t-w                                                               2316
15960 do 15961 j=1,ni                                                      2316
      if(ju(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                         2316
15961 continue                                                             2317
15962 continue                                                             2317
      goto 15741                                                           2318
15742 continue                                                             2318
      if(nin .le. nx)goto 15981                                            2318
      jerr=-10000-ilm                                                      2318
      goto 15672                                                           2318
15981 continue                                                             2319
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2319
      kin(ilm)=nin                                                         2320
      a0(ilm)=az                                                           2320
      alm(ilm)=al                                                          2320
      lmu=ilm                                                              2321
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               2322
      if(ilm.lt.mnl)goto 15671                                             2322
      if(flmin.ge.1.0)goto 15671                                           2323
      me=0                                                                 2323
15990 do 15991 j=1,nin                                                     2323
      if(ca(j,ilm).ne.0.0) me=me+1                                         2323
15991 continue                                                             2323
15992 continue                                                             2323
      if(me.gt.ne)goto 15672                                               2324
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 15672              2325
      if(dev(ilm).gt.devmax)goto 15672                                     2326
15671 continue                                                             2327
15672 continue                                                             2327
      g=f                                                                  2328
11420 continue                                                             2328
      deallocate(t,w,wr,v,a,f,as,mm)                                       2329
      return                                                               2330
      end                                                                  2331
      function nonzero(n,v)                                                2332
      real v(n)                                                            2333
      nonzero=0                                                            2333
16000 do 16001 i=1,n                                                       2333
      if(v(i) .eq. 0.0)goto 16021                                          2333
      nonzero=1                                                            2333
      return                                                               2333
16021 continue                                                             2333
16001 continue                                                             2334
16002 continue                                                             2334
      return                                                               2335
      end                                                                  2336
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               2337
      real a(nx,lmu),b(ni,lmu)                                             2337
      integer ia(nx),nin(lmu)                                              2338
16030 do 16031 lam=1,lmu                                                   2338
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        2338
16031 continue                                                             2339
16032 continue                                                             2339
      return                                                               2340
      end                                                                  2341
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           2342
      real a(nx,nc,lmu),b(ni,nc,lmu)                                       2342
      integer ia(nx),nin(lmu)                                              2343
16040 do 16041 lam=1,lmu                                                   2343
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             2343
16041 continue                                                             2344
16042 continue                                                             2344
      return                                                               2345
      end                                                                  2346
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               2347
      real x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)       2348
      real, dimension (:), allocatable :: w                                     
      if(minval(y) .ge. 0.0)goto 16061                                     2351
      jerr=8888                                                            2351
      return                                                               2351
16061 continue                                                             2352
      allocate(w(1:no),stat=jerr)                                          2352
      if(jerr.ne.0) return                                                 2353
      w=max(0.0,q)                                                         2353
      sw=sum(w)                                                            2353
      if(sw .gt. 0.0)goto 16081                                            2353
      jerr=9999                                                            2353
      go to 11420                                                          2353
16081 continue                                                             2354
      yb=dot_product(w,y)/sw                                               2354
      fmax=log(huge(y(1))*0.1)                                             2355
16090 do 16091 lam=1,nlam                                                  2355
      s=0.0                                                                2356
16100 do 16101 i=1,no                                                      2356
      if(w(i).le.0.0)goto 16101                                            2357
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          2358
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      2359
16101 continue                                                             2360
16102 continue                                                             2360
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2361
16091 continue                                                             2362
16092 continue                                                             2362
11420 continue                                                             2362
      deallocate(w)                                                        2363
      return                                                               2364
      end                                                                  2365
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,ne,nx,nlam,fl   2367 
     *min,ulam,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        2368
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2369
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           2370
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16121                                    2374
      jerr=10000                                                           2374
      return                                                               2374
16121 continue                                                             2375
      if(minval(y) .ge. 0.0)goto 16141                                     2375
      jerr=8888                                                            2375
      return                                                               2375
16141 continue                                                             2376
      allocate(ww(1:no),stat=jerr)                                         2377
      allocate(ju(1:ni),stat=ierr)                                         2377
      jerr=jerr+ierr                                                       2378
      allocate(vq(1:ni),stat=ierr)                                         2378
      jerr=jerr+ierr                                                       2379
      allocate(xm(1:ni),stat=ierr)                                         2379
      jerr=jerr+ierr                                                       2380
      allocate(xs(1:ni),stat=ierr)                                         2380
      jerr=jerr+ierr                                                       2381
      if(jerr.ne.0) return                                                 2382
      call spchkvars(no,ni,x,ix,ju)                                        2383
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2384
      if(maxval(ju) .gt. 0)goto 16161                                      2384
      jerr=7777                                                            2384
      go to 11420                                                          2384
16161 continue                                                             2385
      vq=max(0.0,vp)                                                       2385
      vq=vq*ni/sum(vq)                                                     2386
      ww=max(0.0,w)                                                        2386
      sw=sum(ww)                                                           2386
      if(sw .gt. 0.0)goto 16181                                            2386
      jerr=9999                                                            2386
      go to 11420                                                          2386
16181 continue                                                             2387
      ww=ww/sw                                                             2388
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     2389
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmin,u   2391 
     *lam,thr,  isd,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11420                                            2391
      dev0=2.0*sw*dev0                                                     2392
16190 do 16191 k=1,lmu                                                     2392
      nk=nin(k)                                                            2393
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2394
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2395
16191 continue                                                             2396
16192 continue                                                             2396
11420 continue                                                             2396
      deallocate(ww,ju,vq,xm,xs)                                           2397
      return                                                               2398
      end                                                                  2399
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,ne,nx,nlam,fl   2401 
     *min,ulam,shr,  isd,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,je
     *rr)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2402 
     *9)
      real x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),xs(ni)          2403
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2404
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2405
      real, dimension (:), allocatable :: qy,t,w,wr,v,a,as,xm                   
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                          2410
      allocate(as(1:ni),stat=ierr)                                         2410
      jerr=jerr+ierr                                                       2411
      allocate(t(1:no),stat=ierr)                                          2411
      jerr=jerr+ierr                                                       2412
      allocate(mm(1:ni),stat=ierr)                                         2412
      jerr=jerr+ierr                                                       2413
      allocate(wr(1:no),stat=ierr)                                         2413
      jerr=jerr+ierr                                                       2414
      allocate(v(1:ni),stat=ierr)                                          2414
      jerr=jerr+ierr                                                       2415
      allocate(xm(1:ni),stat=ierr)                                         2415
      jerr=jerr+ierr                                                       2416
      allocate(w(1:no),stat=ierr)                                          2416
      jerr=jerr+ierr                                                       2417
      allocate(qy(1:no),stat=ierr)                                         2417
      jerr=jerr+ierr                                                       2418
      if(jerr.ne.0) return                                                 2419
      bta=max(parm,1.0e-3)                                                 2419
      omb=1.0-bta                                                          2419
      fmax=log(huge(bta)*0.1)                                              2420
      qy=q*y                                                               2420
      yb=sum(qy)                                                           2421
      if(nonzero(no,g) .ne. 0)goto 16211                                   2421
      w=q*yb                                                               2421
      az=log(yb)                                                           2421
      uu=az                                                                2422
      xm=yb*xb                                                             2422
      v=yb                                                                 2422
      t=0.0                                                                2423
      goto 16221                                                           2424
16211 continue                                                             2424
      w=q*exp(sign(min(abs(g),fmax),g))                                    2424
      ww=sum(w)                                                            2424
      eaz=yb/ww                                                            2425
      w=eaz*w                                                              2425
      az=log(eaz)                                                          2425
      uu=az                                                                2425
      t=g                                                                  2426
16230 do 16231 j=1,ni                                                      2426
      if(ju(j).eq.0)goto 16231                                             2426
      jb=ix(j)                                                             2426
      je=ix(j+1)-1                                                         2427
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2428
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   2430 
     *b(j)**2)/xs(j)**2
16231 continue                                                             2431
16232 continue                                                             2431
16221 continue                                                             2432
16201 continue                                                             2432
      tt=yb*uu                                                             2432
      ww=yb                                                                2432
      wr=qy-q*(yb*(1.0-uu))                                                2432
      a=0.0                                                                2433
      dv0=yb*(log(yb)-1.0)                                                 2433
      dvr=-yb                                                              2434
16240 do 16241 i=1,no                                                      2434
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             2434
16241 continue                                                             2434
16242 continue                                                             2434
      dvr=dvr-dv0                                                          2434
      dev0=dvr                                                             2435
      if(flmin .ge. 1.0)goto 16261                                         2435
      eqs=max(eps,flmin)                                                   2435
      alf=eqs**(1.0/(nlam-1))                                              2435
16261 continue                                                             2436
      m=0                                                                  2436
      mm=0                                                                 2436
      nlp=0                                                                2436
      nin=nlp                                                              2436
      mnl=min(mnlam,nlam)                                                  2437
16270 do 16271 ilm=1,nlam                                                  2438
      if(flmin .lt. 1.0)goto 16291                                         2438
      al=ulam(ilm)                                                         2438
      goto 16281                                                           2439
16291 if(ilm .le. 2)goto 16301                                             2439
      al=al*alf                                                            2439
      goto 16281                                                           2440
16301 if(ilm .ne. 1)goto 16311                                             2440
      al=big                                                               2440
      goto 16321                                                           2441
16311 continue                                                             2441
      al=0.0                                                               2442
16330 do 16331 j=1,ni                                                      2442
      if(ju(j).eq.0)goto 16331                                             2442
      if(vp(j).le.0.0)goto 16331                                           2443
      jb=ix(j)                                                             2443
      je=ix(j+1)-1                                                         2444
      al=max(al,  abs(dot_product(qy(jx(jb:je)),x(jb:je))-xm(j))/(xs(j)*   2446 
     *vp(j)))
16331 continue                                                             2447
16332 continue                                                             2447
      al=alf*al/bta                                                        2448
16321 continue                                                             2449
16281 continue                                                             2449
      al2=al*omb                                                           2449
      al1=al*bta                                                           2449
      nit=0                                                                2450
16340 continue                                                             2450
16341 continue                                                             2450
      nit=nit+1                                                            2450
      if(nit .le. maxit)goto 16361                                         2450
      jerr=-2                                                              2450
      go to 11420                                                          2450
16361 continue                                                             2451
      az0=az                                                               2451
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2451
      nit1=0                                                               2452
16370 continue                                                             2452
16371 continue                                                             2452
      nit1=nit1+1                                                          2452
      if(nit1 .le. maxit)goto 16391                                        2452
      jerr=-1                                                              2452
      go to 11420                                                          2452
16391 continue                                                             2453
      nlp=nlp+1                                                            2453
      dlx=0.0                                                              2454
16400 do 16401 k=1,ni                                                      2454
      if(ju(k).eq.0)goto 16401                                             2454
      jb=ix(k)                                                             2454
      je=ix(k+1)-1                                                         2454
      ak=a(k)                                                              2455
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2457 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2458
      if(au .gt. 0.0)goto 16421                                            2458
      a(k)=0.0                                                             2458
      goto 16431                                                           2459
16421 continue                                                             2459
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2459
16431 continue                                                             2460
16411 continue                                                             2460
      if(a(k).eq.ak)goto 16401                                             2461
      if(mm(k) .ne. 0)goto 16451                                           2461
      nin=nin+1                                                            2461
      if(nin.gt.nx)goto 16402                                              2462
      mm(k)=nin                                                            2462
      m(nin)=k                                                             2463
16451 continue                                                             2464
      d=a(k)-ak                                                            2464
      dlx=max(dlx,abs(d))                                                  2464
      dv=d/xs(k)                                                           2465
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2466
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2467
      uu=uu-dv*xb(k)                                                       2467
      tt=tt-dv*xm(k)                                                       2468
16401 continue                                                             2469
16402 continue                                                             2469
      if(nin.gt.nx)goto 16372                                              2469
      d=tt/ww-uu                                                           2470
      az=az+d                                                              2470
      dlx=max(dlx,abs(d))                                                  2470
      uu=uu+d                                                              2471
      if(dlx.lt.shr)goto 16372                                             2471
      nit2=0                                                               2472
16460 continue                                                             2472
16461 continue                                                             2472
      nit2=nit2+1                                                          2472
      if(nit2 .le. maxit)goto 16481                                        2472
      jerr=-1                                                              2472
      go to 11420                                                          2472
16481 continue                                                             2473
      nlp=nlp+1                                                            2473
      dlx=0.0                                                              2474
16490 do 16491 l=1,nin                                                     2474
      k=m(l)                                                               2474
      jb=ix(k)                                                             2474
      je=ix(k+1)-1                                                         2474
      ak=a(k)                                                              2475
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2477 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2478
      if(au .gt. 0.0)goto 16511                                            2478
      a(k)=0.0                                                             2478
      goto 16521                                                           2479
16511 continue                                                             2479
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2479
16521 continue                                                             2480
16501 continue                                                             2480
      if(a(k).eq.ak)goto 16491                                             2480
      d=a(k)-ak                                                            2480
      dlx=max(dlx,abs(d))                                                  2480
      dv=d/xs(k)                                                           2481
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2482
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2483
      uu=uu-dv*xb(k)                                                       2483
      tt=tt-dv*xm(k)                                                       2484
16491 continue                                                             2485
16492 continue                                                             2485
      d=tt/ww-uu                                                           2485
      az=az+d                                                              2485
      dlx=max(dlx,abs(d))                                                  2485
      uu=uu+d                                                              2486
      if(dlx.lt.shr)goto 16462                                             2486
      goto 16461                                                           2487
16462 continue                                                             2487
      goto 16371                                                           2488
16372 continue                                                             2488
      if(nin.gt.nx)goto 16342                                              2489
      euu=exp(sign(min(abs(uu),fmax),uu))                                  2490
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                2490
      ww=sum(w)                                                            2491
      if(abs(az-az0) .ge. shr)goto 16541                                   2491
      ixx=0                                                                2492
16550 do 16551 j=1,nin                                                     2492
      if(abs(a(m(j))-as(m(j))).lt.shr)goto 16551                           2492
      ixx=1                                                                2492
      goto 16552                                                           2492
16551 continue                                                             2493
16552 continue                                                             2493
      if(ixx.eq.0)goto 16342                                               2494
16541 continue                                                             2495
      wr=qy-w*(1.0-uu)                                                     2495
      tt=sum(wr)                                                           2496
16560 do 16561 j=1,ni                                                      2496
      if(ju(j).eq.0)goto 16561                                             2496
      jb=ix(j)                                                             2496
      je=ix(j+1)-1                                                         2497
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2498
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   2500 
     *b(j)**2)/xs(j)**2
16561 continue                                                             2501
16562 continue                                                             2501
      goto 16341                                                           2502
16342 continue                                                             2502
      if(nin .le. nx)goto 16581                                            2502
      jerr=-10000-ilm                                                      2502
      goto 16272                                                           2502
16581 continue                                                             2503
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2503
      kin(ilm)=nin                                                         2504
      a0(ilm)=az                                                           2504
      alm(ilm)=al                                                          2504
      lmu=ilm                                                              2505
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        2506
      if(ilm.lt.mnl)goto 16271                                             2506
      if(flmin.ge.1.0)goto 16271                                           2507
      me=0                                                                 2507
16590 do 16591 j=1,nin                                                     2507
      if(ca(j,ilm).ne.0.0) me=me+1                                         2507
16591 continue                                                             2507
16592 continue                                                             2507
      if(me.gt.ne)goto 16272                                               2508
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 16272              2509
      if(dev(ilm).gt.devmax)goto 16272                                     2510
16271 continue                                                             2511
16272 continue                                                             2511
      g=t+uu                                                               2512
11420 continue                                                             2512
      deallocate(t,w,wr,v,a,qy,xm,as,mm)                                   2513
      return                                                               2514
      end                                                                  2515
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       2516
      real x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)           2517
      integer ix(*),jx(*)                                                  2518
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 16611                                     2521
      jerr=8888                                                            2521
      return                                                               2521
16611 continue                                                             2522
      allocate(w(1:no),stat=jerr)                                          2523
      allocate(f(1:no),stat=ierr)                                          2523
      jerr=jerr+ierr                                                       2524
      if(jerr.ne.0) return                                                 2525
      w=max(0.0,q)                                                         2525
      sw=sum(w)                                                            2525
      if(sw .gt. 0.0)goto 16631                                            2525
      jerr=9999                                                            2525
      go to 11420                                                          2525
16631 continue                                                             2526
      yb=dot_product(w,y)/sw                                               2526
      fmax=log(huge(y(1))*0.1)                                             2527
16640 do 16641 lam=1,nlam                                                  2527
      f=a0(lam)                                                            2528
16650 do 16651 j=1,ni                                                      2528
      if(a(j,lam).eq.0.0)goto 16651                                        2528
      jb=ix(j)                                                             2528
      je=ix(j+1)-1                                                         2529
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          2530
16651 continue                                                             2531
16652 continue                                                             2531
      f=f+g                                                                2532
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2533
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2534
16641 continue                                                             2535
16642 continue                                                             2535
11420 continue                                                             2535
      deallocate(w,f)                                                      2536
      return                                                               2537
      end                                                                  2538
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   2539 
     *jerr)
      real x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(nlam)          2540
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 2541
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 16671                                     2544
      jerr=8888                                                            2544
      return                                                               2544
16671 continue                                                             2545
      allocate(w(1:no),stat=jerr)                                          2546
      allocate(f(1:no),stat=ierr)                                          2546
      jerr=jerr+ierr                                                       2547
      if(jerr.ne.0) return                                                 2548
      w=max(0.0,q)                                                         2548
      sw=sum(w)                                                            2548
      if(sw .gt. 0.0)goto 16691                                            2548
      jerr=9999                                                            2548
      go to 11420                                                          2548
16691 continue                                                             2549
      yb=dot_product(w,y)/sw                                               2549
      fmax=log(huge(y(1))*0.1)                                             2550
16700 do 16701 lam=1,nlam                                                  2550
      f=a0(lam)                                                            2551
16710 do 16711 k=1,nin(lam)                                                2551
      j=ia(k)                                                              2551
      jb=ix(j)                                                             2551
      je=ix(j+1)-1                                                         2552
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         2553
16711 continue                                                             2554
16712 continue                                                             2554
      f=f+g                                                                2555
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2556
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2557
16701 continue                                                             2558
16702 continue                                                             2558
11420 continue                                                             2558
      deallocate(w,f)                                                      2559
      return                                                               2560
      end                                                                  2561
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
