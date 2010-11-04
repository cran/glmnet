c
c                          newGLMnet (10/15/10)
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
      call standard(no,ni,x,y,w,isd,ju,thr,sthr,g,xm,xs,ym,ys,xv,jerr)      648
      if(jerr.ne.0) return                                                  649
      if(flmin.ge.1.0) vlam=ulam/ys                                         650
      call elnet1(parm,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,vlam,sthr,xv,  l    652 
     *mu,ca,ia,nin,rsq,alm,nlp,jerr)
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
      subroutine standard (no,ni,x,y,w,isd,ju,thr,sthr,g,xm,xs,ym,ys,xv,    661 
     *jerr)
      real x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                  661
      integer ju(ni)                                                        662
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           665
      if(jerr.ne.0) return                                                  666
      w=w/sum(w)                                                            666
      v=sqrt(w)                                                             666
      ssd=0.0                                                               666
      nv=0                                                                  667
10100 do 10101 j=1,ni                                                       667
      if(ju(j).eq.0)goto 10101                                              668
      xm(j)=dot_product(w,x(:,j))                                           668
      x(:,j)=v*(x(:,j)-xm(j))                                               669
      xv(j)=dot_product(x(:,j),x(:,j))                                      669
      xs(j)=sqrt(xv(j))                                                     670
      ssd=ssd+xs(j)                                                         670
      nv=nv+1                                                               671
10101 continue                                                              672
10102 continue                                                              672
      if(isd .ne. 0)goto 10121                                              672
      xs=1.0                                                                672
      sthr=thr*ssd/nv                                                       672
      goto 10131                                                            673
10121 continue                                                              674
10140 do 10141 j=1,ni                                                       674
      if(ju(j).eq.0)goto 10141                                              674
      x(:,j)=x(:,j)/xs(j)                                                   674
10141 continue                                                              675
10142 continue                                                              675
      xv=1.0                                                                675
      sthr=thr                                                              676
10131 continue                                                              677
10111 continue                                                              677
      ym=dot_product(w,y)                                                   677
      y=v*(y-ym)                                                            677
      ys=sqrt(dot_product(y,y))                                             677
      y=y/ys                                                                677
      g=0.0                                                                 678
10150 do 10151 j=1,ni                                                       678
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                             678
10151 continue                                                              679
10152 continue                                                              679
      deallocate(v)                                                         680
      return                                                                681
      end                                                                   682
      subroutine elnet1 (beta,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,ulam,thr,    684 
     *xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    685 
     *9)
      real vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    686 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       687
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                           693
      jerr=jerr+ierr                                                        694
      allocate(mm(1:ni),stat=ierr)                                          694
      jerr=jerr+ierr                                                        695
      allocate(da(1:ni),stat=ierr)                                          695
      jerr=jerr+ierr                                                        696
      if(jerr.ne.0) return                                                  697
      bta=max(beta,1.0e-3)                                                  697
      omb=1.0-bta                                                           698
      if(flmin .ge. 1.0)goto 10171                                          698
      eqs=max(eps,flmin)                                                    698
      alf=eqs**(1.0/(nlam-1))                                               698
10171 continue                                                              699
      rsq=0.0                                                               699
      a=0.0                                                                 699
      mm=0                                                                  699
      nlp=0                                                                 699
      nin=nlp                                                               699
      iz=0                                                                  699
      mnl=min(mnlam,nlam)                                                   700
10180 do 10181 m=1,nlam                                                     701
      if(flmin .lt. 1.0)goto 10201                                          701
      alm=ulam(m)                                                           701
      goto 10191                                                            702
10201 if(m .le. 2)goto 10211                                                702
      alm=alm*alf                                                           702
      goto 10191                                                            703
10211 if(m .ne. 1)goto 10221                                                703
      alm=big                                                               703
      goto 10231                                                            704
10221 continue                                                              704
      alm=0.0                                                               705
10240 do 10241 j=1,ni                                                       705
      if(ju(j).eq.0)goto 10241                                              705
      if(vp(j).le.0.0)goto 10241                                            706
      alm=max(alm,abs(g(j))/vp(j))                                          707
10241 continue                                                              708
10242 continue                                                              708
      alm=alf*alm/bta                                                       709
10231 continue                                                              710
10191 continue                                                              710
      dem=alm*omb                                                           710
      ab=alm*bta                                                            710
      rsq0=rsq                                                              710
      jz=1                                                                  711
10250 continue                                                              711
10251 continue                                                              711
      if(iz*jz.ne.0) go to 10260                                            711
      nlp=nlp+1                                                             711
      dlx=0.0                                                               712
10270 do 10271 k=1,ni                                                       712
      if(ju(k).eq.0)goto 10271                                              713
      ak=a(k)                                                               713
      u=g(k)+ak*xv(k)                                                       713
      v=abs(u)-vp(k)*ab                                                     713
      a(k)=0.0                                                              714
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         715
      if(a(k).eq.ak)goto 10271                                              716
      if(mm(k) .ne. 0)goto 10291                                            716
      nin=nin+1                                                             716
      if(nin.gt.nx)goto 10272                                               717
10300 do 10301 j=1,ni                                                       717
      if(ju(j).eq.0)goto 10301                                              718
      if(mm(j) .eq. 0)goto 10321                                            718
      c(j,nin)=c(k,mm(j))                                                   718
      goto 10301                                                            718
10321 continue                                                              719
      if(j .ne. k)goto 10341                                                719
      c(j,nin)=xv(j)                                                        719
      goto 10301                                                            719
10341 continue                                                              720
      c(j,nin)=dot_product(x(:,j),x(:,k))                                   721
10301 continue                                                              722
10302 continue                                                              722
      mm(k)=nin                                                             722
      ia(nin)=k                                                             723
10291 continue                                                              724
      del=a(k)-ak                                                           724
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      725
      dlx=max(xv(k)*abs(del),dlx)                                           726
10350 do 10351 j=1,ni                                                       726
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                               726
10351 continue                                                              727
10352 continue                                                              727
10271 continue                                                              728
10272 continue                                                              728
      if(dlx.lt.thr)goto 10252                                              728
      if(nin.gt.nx)goto 10252                                               729
10260 continue                                                              729
      iz=1                                                                  729
      da(1:nin)=a(ia(1:nin))                                                730
10360 continue                                                              730
10361 continue                                                              730
      nlp=nlp+1                                                             730
      dlx=0.0                                                               731
10370 do 10371 l=1,nin                                                      731
      k=ia(l)                                                               731
      ak=a(k)                                                               731
      u=g(k)+ak*xv(k)                                                       731
      v=abs(u)-vp(k)*ab                                                     732
      a(k)=0.0                                                              733
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         734
      if(a(k).eq.ak)goto 10371                                              735
      del=a(k)-ak                                                           735
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      736
      dlx=max(xv(k)*abs(del),dlx)                                           737
10380 do 10381 j=1,nin                                                      737
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                  737
10381 continue                                                              738
10382 continue                                                              738
10371 continue                                                              739
10372 continue                                                              739
      if(dlx.lt.thr)goto 10362                                              739
      goto 10361                                                            740
10362 continue                                                              740
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                      741
10390 do 10391 j=1,ni                                                       741
      if(mm(j).ne.0)goto 10391                                              742
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))            743
10391 continue                                                              744
10392 continue                                                              744
      jz=0                                                                  745
      goto 10251                                                            746
10252 continue                                                              746
      if(nin.gt.nx)goto 10182                                               747
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 747
      kin(m)=nin                                                            748
      rsqo(m)=rsq                                                           748
      almo(m)=alm                                                           748
      lmu=m                                                                 749
      if(m.lt.mnl)goto 10181                                                749
      if(flmin.ge.1.0)goto 10181                                            750
      me=0                                                                  750
10400 do 10401 j=1,nin                                                      750
      if(ao(j,m).ne.0.0) me=me+1                                            750
10401 continue                                                              750
10402 continue                                                              750
      if(me.gt.ne)goto 10182                                                751
      if(rsq-rsq0.lt.sml*rsq)goto 10182                                     751
      if(rsq.gt.rsqmax)goto 10182                                           752
10181 continue                                                              753
10182 continue                                                              753
      deallocate(a,mm,c,da)                                                 754
      return                                                                755
      end                                                                   756
      subroutine elnetn (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,th    758 
     *r,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no),w(no),ulam(nlam)                           759
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         760
      integer jd(*),ia(nx),nin(nlam)                                        761
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          766
      allocate(xs(1:ni),stat=ierr)                                          766
      jerr=jerr+ierr                                                        767
      allocate(ju(1:ni),stat=ierr)                                          767
      jerr=jerr+ierr                                                        768
      allocate(xv(1:ni),stat=ierr)                                          768
      jerr=jerr+ierr                                                        769
      allocate(vlam(1:nlam),stat=ierr)                                      769
      jerr=jerr+ierr                                                        770
      if(jerr.ne.0) return                                                  771
      call chkvars(no,ni,x,ju)                                              772
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  773
      if(maxval(ju) .gt. 0)goto 10421                                       773
      jerr=7777                                                             773
      return                                                                773
10421 continue                                                              774
      call standard1(no,ni,x,y,w,isd,ju,thr,sthr,xm,xs,ym,ys,xv,jerr)       775
      if(jerr.ne.0) return                                                  776
      if(flmin.ge.1.0) vlam=ulam/ys                                         777
      call elnet2(parm,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,vlam,sthr,xv,  l    779 
     *mu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                  780
10430 do 10431 k=1,lmu                                                      780
      alm(k)=ys*alm(k)                                                      780
      nk=nin(k)                                                             781
10440 do 10441 l=1,nk                                                       781
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          781
10441 continue                                                              782
10442 continue                                                              782
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         783
10431 continue                                                              784
10432 continue                                                              784
      deallocate(xm,xs,ju,xv,vlam)                                          785
      return                                                                786
      end                                                                   787
      subroutine standard1 (no,ni,x,y,w,isd,ju,thr,sthr,xm,xs,ym,ys,xv,j    788 
     *err)
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)                        788
      integer ju(ni)                                                        789
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           792
      if(jerr.ne.0) return                                                  793
      w=w/sum(w)                                                            793
      v=sqrt(w)                                                             793
      ssd=0.0                                                               793
      nv=0                                                                  794
10450 do 10451 j=1,ni                                                       794
      if(ju(j).eq.0)goto 10451                                              795
      xm(j)=dot_product(w,x(:,j))                                           795
      x(:,j)=v*(x(:,j)-xm(j))                                               796
      xv(j)=dot_product(x(:,j),x(:,j))                                      796
      xs(j)=sqrt(xv(j))                                                     797
      ssd=ssd+xs(j)                                                         797
      nv=nv+1                                                               798
10451 continue                                                              799
10452 continue                                                              799
      if(isd .ne. 0)goto 10471                                              799
      xs=1.0                                                                799
      sthr=thr*ssd/nv                                                       799
      goto 10481                                                            800
10471 continue                                                              800
10490 do 10491 j=1,ni                                                       800
      if(ju(j).eq.0)goto 10491                                              800
      x(:,j)=x(:,j)/xs(j)                                                   800
10491 continue                                                              801
10492 continue                                                              801
      xv=1.0                                                                801
      sthr=thr                                                              802
10481 continue                                                              803
10461 continue                                                              803
      ym=dot_product(w,y)                                                   803
      y=v*(y-ym)                                                            803
      ys=sqrt(dot_product(y,y))                                             803
      y=y/ys                                                                804
      deallocate(v)                                                         805
      return                                                                806
      end                                                                   807
      subroutine elnet2(beta,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,ulam,thr,x    809 
     *v,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    810 
     *9)
      real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    811 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       812
      real, dimension (:), allocatable :: a                                     
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                           817
      allocate(mm(1:ni),stat=ierr)                                          817
      jerr=jerr+ierr                                                        818
      if(jerr.ne.0) return                                                  819
      bta=max(beta,1.0e-3)                                                  819
      omb=1.0-bta                                                           820
      if(flmin .ge. 1.0)goto 10511                                          820
      eqs=max(eps,flmin)                                                    820
      alf=eqs**(1.0/(nlam-1))                                               820
10511 continue                                                              821
      rsq=0.0                                                               821
      a=0.0                                                                 821
      mm=0                                                                  821
      nlp=0                                                                 821
      nin=nlp                                                               821
      iz=0                                                                  821
      mnl=min(mnlam,nlam)                                                   822
10520 do 10521 m=1,nlam                                                     823
      if(flmin .lt. 1.0)goto 10541                                          823
      alm=ulam(m)                                                           823
      goto 10531                                                            824
10541 if(m .le. 2)goto 10551                                                824
      alm=alm*alf                                                           824
      goto 10531                                                            825
10551 if(m .ne. 1)goto 10561                                                825
      alm=big                                                               825
      goto 10571                                                            826
10561 continue                                                              826
      alm=0.0                                                               827
10580 do 10581 j=1,ni                                                       827
      if(ju(j).eq.0)goto 10581                                              827
      if(vp(j).le.0.0)goto 10581                                            828
      alm=max(alm,abs(dot_product(y,x(:,j)))/vp(j))                         829
10581 continue                                                              830
10582 continue                                                              830
      alm=alf*alm/bta                                                       831
10571 continue                                                              832
10531 continue                                                              832
      dem=alm*omb                                                           832
      ab=alm*bta                                                            832
      rsq0=rsq                                                              832
      jz=1                                                                  833
10590 continue                                                              833
10591 continue                                                              833
      if(iz*jz.ne.0) go to 10260                                            833
      nlp=nlp+1                                                             833
      dlx=0.0                                                               834
10600 do 10601 k=1,ni                                                       834
      if(ju(k).eq.0)goto 10601                                              834
      gk=dot_product(y,x(:,k))                                              835
      ak=a(k)                                                               835
      u=gk+ak*xv(k)                                                         835
      v=abs(u)-vp(k)*ab                                                     835
      a(k)=0.0                                                              836
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         837
      if(a(k).eq.ak)goto 10601                                              838
      if(mm(k) .ne. 0)goto 10621                                            838
      nin=nin+1                                                             838
      if(nin.gt.nx)goto 10602                                               839
      mm(k)=nin                                                             839
      ia(nin)=k                                                             840
10621 continue                                                              841
      del=a(k)-ak                                                           841
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        842
      y=y-del*x(:,k)                                                        842
      dlx=max(xv(k)*abs(del),dlx)                                           843
10601 continue                                                              844
10602 continue                                                              844
      if(dlx.lt.thr)goto 10592                                              844
      if(nin.gt.nx)goto 10592                                               845
10260 continue                                                              845
      iz=1                                                                  846
10630 continue                                                              846
10631 continue                                                              846
      nlp=nlp+1                                                             846
      dlx=0.0                                                               847
10640 do 10641 l=1,nin                                                      847
      k=ia(l)                                                               847
      gk=dot_product(y,x(:,k))                                              848
      ak=a(k)                                                               848
      u=gk+ak*xv(k)                                                         848
      v=abs(u)-vp(k)*ab                                                     848
      a(k)=0.0                                                              849
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         850
      if(a(k).eq.ak)goto 10641                                              851
      del=a(k)-ak                                                           851
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        852
      y=y-del*x(:,k)                                                        852
      dlx=max(xv(k)*abs(del),dlx)                                           853
10641 continue                                                              854
10642 continue                                                              854
      if(dlx.lt.thr)goto 10632                                              854
      goto 10631                                                            855
10632 continue                                                              855
      jz=0                                                                  856
      goto 10591                                                            857
10592 continue                                                              857
      if(nin.gt.nx)goto 10522                                               858
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 858
      kin(m)=nin                                                            859
      rsqo(m)=rsq                                                           859
      almo(m)=alm                                                           859
      lmu=m                                                                 860
      if(m.lt.mnl)goto 10521                                                860
      if(flmin.ge.1.0)goto 10521                                            861
      me=0                                                                  861
10650 do 10651 j=1,nin                                                      861
      if(ao(j,m).ne.0.0) me=me+1                                            861
10651 continue                                                              861
10652 continue                                                              861
      if(me.gt.ne)goto 10522                                                862
      if(rsq-rsq0.lt.sml*rsq)goto 10522                                     862
      if(rsq.gt.rsqmax)goto 10522                                           863
10521 continue                                                              864
10522 continue                                                              864
      deallocate(a,mm)                                                      865
      return                                                                866
      end                                                                   867
      subroutine chkvars(no,ni,x,ju)                                        868
      real x(no,ni)                                                         868
      integer ju(ni)                                                        869
10660 do 10661 j=1,ni                                                       869
      ju(j)=0                                                               869
      t=x(1,j)                                                              870
10670 do 10671 i=2,no                                                       870
      if(x(i,j).eq.t)goto 10671                                             870
      ju(j)=1                                                               870
      goto 10672                                                            870
10671 continue                                                              871
10672 continue                                                              871
10661 continue                                                              872
10662 continue                                                              872
      return                                                                873
      end                                                                   874
      subroutine uncomp(ni,ca,ia,nin,a)                                     875
      real ca(*),a(ni)                                                      875
      integer ia(*)                                                         876
      a=0.0                                                                 876
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                   877
      return                                                                878
      end                                                                   879
      subroutine modval(a0,ca,ia,nin,n,x,f)                                 880
      real ca(nin),x(n,*),f(n)                                              880
      integer ia(nin)                                                       881
      f=a0                                                                  881
      if(nin.le.0) return                                                   882
10680 do 10681 i=1,n                                                        882
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                       882
10681 continue                                                              883
10682 continue                                                              883
      return                                                                884
      end                                                                   885
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,fl    888 
     *min,ulam,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                               889
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         890
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                            891
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10701                                     894
      jerr=10000                                                            894
      return                                                                894
10701 continue                                                              895
      allocate(vq(1:ni),stat=jerr)                                          895
      if(jerr.ne.0) return                                                  896
      vq=max(0.0,vp)                                                        896
      vq=vq*ni/sum(vq)                                                      897
      if(ka .ne. 1)goto 10721                                               898
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam    901 
     *,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10731                                                            902
10721 continue                                                              903
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam,    906 
     *thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10731 continue                                                              907
10711 continue                                                              907
      deallocate(vq)                                                        908
      return                                                                909
      end                                                                   910
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmi    913 
     *n,ulam,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                               914
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         915
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                            916
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           921
      allocate(xm(1:ni),stat=ierr)                                          921
      jerr=jerr+ierr                                                        922
      allocate(xs(1:ni),stat=ierr)                                          922
      jerr=jerr+ierr                                                        923
      allocate(ju(1:ni),stat=ierr)                                          923
      jerr=jerr+ierr                                                        924
      allocate(xv(1:ni),stat=ierr)                                          924
      jerr=jerr+ierr                                                        925
      allocate(vlam(1:nlam),stat=ierr)                                      925
      jerr=jerr+ierr                                                        926
      if(jerr.ne.0) return                                                  927
      call spchkvars(no,ni,x,ix,ju)                                         928
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  929
      if(maxval(ju) .gt. 0)goto 10751                                       929
      jerr=7777                                                             929
      return                                                                929
10751 continue                                                              930
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,thr,sthr,g,xm,xs,ym,ys,xv    931 
     *,jerr)
      if(jerr.ne.0) return                                                  932
      if(flmin.ge.1.0) vlam=ulam/ys                                         933
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,s    935 
     *thr,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                  936
10760 do 10761 k=1,lmu                                                      936
      alm(k)=ys*alm(k)                                                      936
      nk=nin(k)                                                             937
10770 do 10771 l=1,nk                                                       937
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          937
10771 continue                                                              938
10772 continue                                                              938
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         939
10761 continue                                                              940
10762 continue                                                              940
      deallocate(xm,xs,g,ju,xv,vlam)                                        941
      return                                                                942
      end                                                                   943
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,thr,sthr,g,xm,xs,y    944 
     *m,ys,xv,jerr)
      real x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                      944
      integer ix(*),jx(*),ju(ni)                                            945
      w=w/sum(w)                                                            945
      ssd=0.0                                                               945
      nv=0                                                                  946
10780 do 10781 j=1,ni                                                       946
      if(ju(j).eq.0)goto 10781                                              947
      jb=ix(j)                                                              947
      je=ix(j+1)-1                                                          947
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                              948
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                  948
      xs(j)=sqrt(xv(j))                                                     949
      ssd=ssd+xs(j)                                                         949
      nv=nv+1                                                               950
10781 continue                                                              951
10782 continue                                                              951
      if(isd .ne. 0)goto 10801                                              951
      xs=1.0                                                                951
      sthr=thr*ssd/nv                                                       951
      goto 10811                                                            952
10801 continue                                                              952
      xv=1.0                                                                952
      sthr=thr                                                              952
10811 continue                                                              953
10791 continue                                                              953
      ym=dot_product(w,y)                                                   953
      y=y-ym                                                                953
      ys=sqrt(dot_product(w,y**2))                                          953
      y=y/ys                                                                953
      g=0.0                                                                 954
10820 do 10821 j=1,ni                                                       954
      if(ju(j).eq.0)goto 10821                                              954
      jb=ix(j)                                                              954
      je=ix(j+1)-1                                                          955
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)            956
10821 continue                                                              957
10822 continue                                                              957
      return                                                                958
      end                                                                   959
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,    961 
     *ulam,  thr,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    962 
     *9)
      real g(ni),vp(ni),x(*),ulam(nlam),w(no)                               963
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)           964
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                           965
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                           971
      jerr=jerr+ierr                                                        972
      allocate(mm(1:ni),stat=ierr)                                          972
      jerr=jerr+ierr                                                        973
      allocate(da(1:ni),stat=ierr)                                          973
      jerr=jerr+ierr                                                        974
      if(jerr.ne.0) return                                                  975
      bta=max(beta,1.0e-3)                                                  975
      omb=1.0-bta                                                           976
      if(flmin .ge. 1.0)goto 10841                                          976
      eqs=max(eps,flmin)                                                    976
      alf=eqs**(1.0/(nlam-1))                                               976
10841 continue                                                              977
      rsq=0.0                                                               977
      a=0.0                                                                 977
      mm=0                                                                  977
      nlp=0                                                                 977
      nin=nlp                                                               977
      iz=0                                                                  977
      mnl=min(mnlam,nlam)                                                   978
10850 do 10851 m=1,nlam                                                     979
      if(flmin .lt. 1.0)goto 10871                                          979
      alm=ulam(m)                                                           979
      goto 10861                                                            980
10871 if(m .le. 2)goto 10881                                                980
      alm=alm*alf                                                           980
      goto 10861                                                            981
10881 if(m .ne. 1)goto 10891                                                981
      alm=big                                                               981
      goto 10901                                                            982
10891 continue                                                              982
      alm=0.0                                                               983
10910 do 10911 j=1,ni                                                       983
      if(ju(j).eq.0)goto 10911                                              983
      if(vp(j).le.0.0)goto 10911                                            984
      alm=max(alm,abs(g(j))/vp(j))                                          985
10911 continue                                                              986
10912 continue                                                              986
      alm=alf*alm/bta                                                       987
10901 continue                                                              988
10861 continue                                                              988
      dem=alm*omb                                                           988
      ab=alm*bta                                                            988
      rsq0=rsq                                                              988
      jz=1                                                                  989
10920 continue                                                              989
10921 continue                                                              989
      if(iz*jz.ne.0) go to 10260                                            989
      nlp=nlp+1                                                             989
      dlx=0.0                                                               990
10930 do 10931 k=1,ni                                                       990
      if(ju(k).eq.0)goto 10931                                              991
      ak=a(k)                                                               991
      u=g(k)+ak*xv(k)                                                       991
      v=abs(u)-vp(k)*ab                                                     991
      a(k)=0.0                                                              992
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         993
      if(a(k).eq.ak)goto 10931                                              994
      if(mm(k) .ne. 0)goto 10951                                            994
      nin=nin+1                                                             994
      if(nin.gt.nx)goto 10932                                               995
10960 do 10961 j=1,ni                                                       995
      if(ju(j).eq.0)goto 10961                                              996
      if(mm(j) .eq. 0)goto 10981                                            996
      c(j,nin)=c(k,mm(j))                                                   996
      goto 10961                                                            996
10981 continue                                                              997
      if(j .ne. k)goto 11001                                                997
      c(j,nin)=xv(j)                                                        997
      goto 10961                                                            997
11001 continue                                                              998
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))       1000
10961 continue                                                             1001
10962 continue                                                             1001
      mm(k)=nin                                                            1001
      ia(nin)=k                                                            1002
10951 continue                                                             1003
      del=a(k)-ak                                                          1003
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1004
      dlx=max(xv(k)*abs(del),dlx)                                          1005
11010 do 11011 j=1,ni                                                      1005
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                              1005
11011 continue                                                             1006
11012 continue                                                             1006
10931 continue                                                             1007
10932 continue                                                             1007
      if(dlx.lt.thr)goto 10922                                             1007
      if(nin.gt.nx)goto 10922                                              1008
10260 continue                                                             1008
      iz=1                                                                 1008
      da(1:nin)=a(ia(1:nin))                                               1009
11020 continue                                                             1009
11021 continue                                                             1009
      nlp=nlp+1                                                            1009
      dlx=0.0                                                              1010
11030 do 11031 l=1,nin                                                     1010
      k=ia(l)                                                              1011
      ak=a(k)                                                              1011
      u=g(k)+ak*xv(k)                                                      1011
      v=abs(u)-vp(k)*ab                                                    1011
      a(k)=0.0                                                             1012
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1013
      if(a(k).eq.ak)goto 11031                                             1014
      del=a(k)-ak                                                          1014
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1015
      dlx=max(xv(k)*abs(del),dlx)                                          1016
11040 do 11041 j=1,nin                                                     1016
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                 1016
11041 continue                                                             1017
11042 continue                                                             1017
11031 continue                                                             1018
11032 continue                                                             1018
      if(dlx.lt.thr)goto 11022                                             1018
      goto 11021                                                           1019
11022 continue                                                             1019
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                     1020
11050 do 11051 j=1,ni                                                      1020
      if(mm(j).ne.0)goto 11051                                             1021
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))           1022
11051 continue                                                             1023
11052 continue                                                             1023
      jz=0                                                                 1024
      goto 10921                                                           1025
10922 continue                                                             1025
      if(nin.gt.nx)goto 10852                                              1026
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1026
      kin(m)=nin                                                           1027
      rsqo(m)=rsq                                                          1027
      almo(m)=alm                                                          1027
      lmu=m                                                                1028
      if(m.lt.mnl)goto 10851                                               1028
      if(flmin.ge.1.0)goto 10851                                           1029
      me=0                                                                 1029
11060 do 11061 j=1,nin                                                     1029
      if(ao(j,m).ne.0.0) me=me+1                                           1029
11061 continue                                                             1029
11062 continue                                                             1029
      if(me.gt.ne)goto 10852                                               1030
      if(rsq-rsq0.lt.sml*rsq)goto 10852                                    1030
      if(rsq.gt.rsqmax)goto 10852                                          1031
10851 continue                                                             1032
10852 continue                                                             1032
      deallocate(a,mm,c,da)                                                1033
      return                                                               1034
      end                                                                  1035
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmin,   1037 
     *ulam,  thr,isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),vp(ni),y(no),w(no),ulam(nlam)                              1038
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1039
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1040
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         1045
      allocate(xs(1:ni),stat=ierr)                                         1045
      jerr=jerr+ierr                                                       1046
      allocate(ju(1:ni),stat=ierr)                                         1046
      jerr=jerr+ierr                                                       1047
      allocate(xv(1:ni),stat=ierr)                                         1047
      jerr=jerr+ierr                                                       1048
      allocate(vlam(1:nlam),stat=ierr)                                     1048
      jerr=jerr+ierr                                                       1049
      if(jerr.ne.0) return                                                 1050
      call spchkvars(no,ni,x,ix,ju)                                        1051
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1052
      if(maxval(ju) .gt. 0)goto 11081                                      1052
      jerr=7777                                                            1052
      return                                                               1052
11081 continue                                                             1053
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,thr,sthr,xm,xs,ym,ys,xv,   1054 
     *jerr)
      if(jerr.ne.0) return                                                 1055
      if(flmin.ge.1.0) vlam=ulam/ys                                        1056
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,s   1058 
     *thr,xm,xs,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                 1059
11090 do 11091 k=1,lmu                                                     1059
      alm(k)=ys*alm(k)                                                     1059
      nk=nin(k)                                                            1060
11100 do 11101 l=1,nk                                                      1060
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1060
11101 continue                                                             1061
11102 continue                                                             1061
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                        1062
11091 continue                                                             1063
11092 continue                                                             1063
      deallocate(xm,xs,ju,xv,vlam)                                         1064
      return                                                               1065
      end                                                                  1066
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,thr,sthr,xm,xs,ym   1067 
     *,ys,xv,jerr)
      real x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)                           1067
      integer ix(*),jx(*),ju(ni)                                           1068
      w=w/sum(w)                                                           1068
      ssd=0.0                                                              1068
      nv=0                                                                 1069
11110 do 11111 j=1,ni                                                      1069
      if(ju(j).eq.0)goto 11111                                             1070
      jb=ix(j)                                                             1070
      je=ix(j+1)-1                                                         1070
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1071
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1072
      xs(j)=sqrt(xv(j))                                                    1072
      ssd=ssd+xs(j)                                                        1072
      nv=nv+1                                                              1073
11111 continue                                                             1074
11112 continue                                                             1074
      if(isd .ne. 0)goto 11131                                             1074
      xs=1.0                                                               1074
      sthr=thr*ssd/nv                                                      1074
      goto 11141                                                           1075
11131 continue                                                             1075
      xv=1.0                                                               1075
      sthr=thr                                                             1075
11141 continue                                                             1076
11121 continue                                                             1076
      ym=dot_product(w,y)                                                  1076
      y=y-ym                                                               1076
      ys=sqrt(dot_product(w,y**2))                                         1076
      y=y/ys                                                               1077
      return                                                               1078
      end                                                                  1079
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,   1081 
     *ulam,  thr,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99   1082 
     *9)
      real y(no),w(no),x(*),vp(ni),ulam(nlam)                              1083
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)          1084
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1085
      real, dimension (:), allocatable :: a                                     
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                          1090
      allocate(mm(1:ni),stat=ierr)                                         1090
      jerr=jerr+ierr                                                       1091
      if(jerr.ne.0) return                                                 1092
      bta=max(beta,1.0e-3)                                                 1092
      omb=1.0-bta                                                          1093
      if(flmin .ge. 1.0)goto 11161                                         1093
      eqs=max(eps,flmin)                                                   1093
      alf=eqs**(1.0/(nlam-1))                                              1093
11161 continue                                                             1094
      rsq=0.0                                                              1094
      a=0.0                                                                1094
      mm=0                                                                 1094
      o=0.0                                                                1094
      nlp=0                                                                1094
      nin=nlp                                                              1094
      iz=0                                                                 1094
      mnl=min(mnlam,nlam)                                                  1095
11170 do 11171 m=1,nlam                                                    1096
      if(flmin .lt. 1.0)goto 11191                                         1096
      alm=ulam(m)                                                          1096
      goto 11181                                                           1097
11191 if(m .le. 2)goto 11201                                               1097
      alm=alm*alf                                                          1097
      goto 11181                                                           1098
11201 if(m .ne. 1)goto 11211                                               1098
      alm=big                                                              1098
      goto 11221                                                           1099
11211 continue                                                             1099
      alm=0.0                                                              1100
11230 do 11231 j=1,ni                                                      1100
      if(ju(j).eq.0)goto 11231                                             1100
      if(vp(j).le.0.0)goto 11231                                           1101
      jb=ix(j)                                                             1101
      je=ix(j+1)-1                                                         1102
      alm=max(alm,abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))    1104 
     * /(vp(j)*xs(j))))
11231 continue                                                             1105
11232 continue                                                             1105
      alm=alf*alm/bta                                                      1106
11221 continue                                                             1107
11181 continue                                                             1107
      dem=alm*omb                                                          1107
      ab=alm*bta                                                           1107
      rsq0=rsq                                                             1107
      jz=1                                                                 1108
11240 continue                                                             1108
11241 continue                                                             1108
      if(iz*jz.ne.0) go to 10260                                           1108
      nlp=nlp+1                                                            1108
      dlx=0.0                                                              1109
11250 do 11251 k=1,ni                                                      1109
      if(ju(k).eq.0)goto 11251                                             1109
      jb=ix(k)                                                             1109
      je=ix(k+1)-1                                                         1110
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1111
      ak=a(k)                                                              1111
      u=gk+ak*xv(k)                                                        1111
      v=abs(u)-vp(k)*ab                                                    1111
      a(k)=0.0                                                             1112
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1113
      if(a(k).eq.ak)goto 11251                                             1114
      if(mm(k) .ne. 0)goto 11271                                           1114
      nin=nin+1                                                            1114
      if(nin.gt.nx)goto 11252                                              1115
      mm(k)=nin                                                            1115
      ia(nin)=k                                                            1116
11271 continue                                                             1117
      del=a(k)-ak                                                          1117
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1118
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1119
      o=o+del*xm(k)/xs(k)                                                  1119
      dlx=max(xv(k)*abs(del),dlx)                                          1120
11251 continue                                                             1121
11252 continue                                                             1121
      if(dlx.lt.thr)goto 11242                                             1121
      if(nin.gt.nx)goto 11242                                              1122
10260 continue                                                             1122
      iz=1                                                                 1123
11280 continue                                                             1123
11281 continue                                                             1123
      nlp=nlp+1                                                            1123
      dlx=0.0                                                              1124
11290 do 11291 l=1,nin                                                     1124
      k=ia(l)                                                              1124
      jb=ix(k)                                                             1124
      je=ix(k+1)-1                                                         1125
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1126
      ak=a(k)                                                              1126
      u=gk+ak*xv(k)                                                        1126
      v=abs(u)-vp(k)*ab                                                    1126
      a(k)=0.0                                                             1127
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1128
      if(a(k).eq.ak)goto 11291                                             1129
      del=a(k)-ak                                                          1129
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1130
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1131
      o=o+del*xm(k)/xs(k)                                                  1131
      dlx=max(xv(k)*abs(del),dlx)                                          1132
11291 continue                                                             1133
11292 continue                                                             1133
      if(dlx.lt.thr)goto 11282                                             1133
      goto 11281                                                           1134
11282 continue                                                             1134
      jz=0                                                                 1135
      goto 11241                                                           1136
11242 continue                                                             1136
      if(nin.gt.nx)goto 11172                                              1137
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1137
      kin(m)=nin                                                           1138
      rsqo(m)=rsq                                                          1138
      almo(m)=alm                                                          1138
      lmu=m                                                                1139
      if(m.lt.mnl)goto 11171                                               1139
      if(flmin.ge.1.0)goto 11171                                           1140
      me=0                                                                 1140
11300 do 11301 j=1,nin                                                     1140
      if(ao(j,m).ne.0.0) me=me+1                                           1140
11301 continue                                                             1140
11302 continue                                                             1140
      if(me.gt.ne)goto 11172                                               1141
      if(rsq-rsq0.lt.sml*rsq)goto 11172                                    1141
      if(rsq.gt.rsqmax)goto 11172                                          1142
11171 continue                                                             1143
11172 continue                                                             1143
      deallocate(a,mm)                                                     1144
      return                                                               1145
      end                                                                  1146
      subroutine spchkvars(no,ni,x,ix,ju)                                  1147
      real x(*)                                                            1147
      integer ix(*),ju(ni)                                                 1148
11310 do 11311 j=1,ni                                                      1148
      ju(j)=0                                                              1148
      jb=ix(j)                                                             1148
      nj=ix(j+1)-jb                                                        1148
      if(nj.eq.0)goto 11311                                                1149
      je=ix(j+1)-1                                                         1150
      if(nj .ge. no)goto 11331                                             1150
11340 do 11341 i=jb,je                                                     1150
      if(x(i).eq.0.0)goto 11341                                            1150
      ju(j)=1                                                              1150
      goto 11342                                                           1150
11341 continue                                                             1150
11342 continue                                                             1150
      goto 11351                                                           1151
11331 continue                                                             1151
      t=x(jb)                                                              1151
11360 do 11361 i=jb+1,je                                                   1151
      if(x(i).eq.t)goto 11361                                              1151
      ju(j)=1                                                              1151
      goto 11362                                                           1151
11361 continue                                                             1151
11362 continue                                                             1151
11351 continue                                                             1152
11321 continue                                                             1152
11311 continue                                                             1153
11312 continue                                                             1153
      return                                                               1154
      end                                                                  1155
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                         1156
      real ca(*),x(*),f(n)                                                 1156
      integer ia(*),ix(*),jx(*)                                            1157
      f=a0                                                                 1158
11370 do 11371 j=1,nin                                                     1158
      k=ia(j)                                                              1158
      kb=ix(k)                                                             1158
      ke=ix(k+1)-1                                                         1159
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                             1160
11371 continue                                                             1161
11372 continue                                                             1161
      return                                                               1162
      end                                                                  1163
      function row_prod(i,j,ia,ja,ra,w)                                    1164
      integer ia(*),ja(*)                                                  1164
      real ra(*),w(*)                                                      1165
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(   1167 
     *i),ia(j+1)-ia(j),w)
      return                                                               1168
      end                                                                  1169
      function dot(x,y,mx,my,nx,ny,w)                                      1170
      real x(*),y(*),w(*)                                                  1170
      integer mx(*),my(*)                                                  1171
      i=1                                                                  1171
      j=i                                                                  1171
      s=0.0                                                                1172
11380 continue                                                             1172
11381 continue                                                             1172
11390 continue                                                             1173
11391 if(mx(i).ge.my(j))goto 11392                                         1173
      i=i+1                                                                1173
      if(i.gt.nx) go to 11400                                              1173
      goto 11391                                                           1174
11392 continue                                                             1174
      if(mx(i).eq.my(j)) go to 11410                                       1175
11420 continue                                                             1175
11421 if(my(j).ge.mx(i))goto 11422                                         1175
      j=j+1                                                                1175
      if(j.gt.ny) go to 11400                                              1175
      goto 11421                                                           1176
11422 continue                                                             1176
      if(mx(i).eq.my(j)) go to 11410                                       1176
      goto 11381                                                           1177
11410 continue                                                             1177
      s=s+w(mx(i))*x(i)*y(j)                                               1178
      i=i+1                                                                1178
      if(i.gt.nx)goto 11382                                                1178
      j=j+1                                                                1178
      if(j.gt.ny)goto 11382                                                1179
      goto 11381                                                           1180
11382 continue                                                             1180
11400 continue                                                             1180
      dot=s                                                                1181
      return                                                               1182
      end                                                                  1183
      subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,ne,nx,nlam,flmin,ulam   1185 
     *,thr,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)             1186
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1187
      integer jd(*),ia(nx),nin(nlam)                                       1188
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 11441                                    1192
      jerr=10000                                                           1192
      return                                                               1192
11441 continue                                                             1193
      allocate(ww(1:no),stat=jerr)                                         1194
      allocate(ju(1:ni),stat=ierr)                                         1194
      jerr=jerr+ierr                                                       1195
      allocate(vq(1:ni),stat=ierr)                                         1195
      jerr=jerr+ierr                                                       1196
      allocate(xm(1:ni),stat=ierr)                                         1196
      jerr=jerr+ierr                                                       1197
      if(isd .le. 0)goto 11461                                             1197
      allocate(xs(1:ni),stat=ierr)                                         1197
      jerr=jerr+ierr                                                       1197
11461 continue                                                             1198
      if(jerr.ne.0) return                                                 1199
      call chkvars(no,ni,x,ju)                                             1200
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1201
      if(maxval(ju) .gt. 0)goto 11481                                      1201
      jerr=7777                                                            1201
      return                                                               1201
11481 continue                                                             1202
      vq=max(0.0,vp)                                                       1202
      vq=vq*ni/sum(vq)                                                     1203
11490 do 11491 i=1,no                                                      1203
      ww(i)=sum(y(i,:))                                                    1203
      y(i,:)=y(i,:)/ww(i)                                                  1203
11491 continue                                                             1203
11492 continue                                                             1203
      sw=sum(ww)                                                           1203
      ww=ww/sw                                                             1204
      call lstandard1(no,ni,x,ww,ju,isd,thr,sthr,xm,xs)                    1205
      if(nc .ne. 1)goto 11511                                              1206
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,ne,nx,nlam,flmin   1208 
     *,ulam,  sthr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr
     *)
      goto 11521                                                           1209
11511 continue                                                             1210
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,st   1212 
     *hr,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
11521 continue                                                             1213
11501 continue                                                             1213
      if(jerr.gt.0) return                                                 1213
      dev0=2.0*sw*dev0                                                     1214
11530 do 11531 k=1,lmu                                                     1214
      nk=nin(k)                                                            1215
11540 do 11541 ic=1,nc                                                     1215
      if(isd .le. 0)goto 11561                                             1215
11570 do 11571 l=1,nk                                                      1215
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1215
11571 continue                                                             1215
11572 continue                                                             1215
11561 continue                                                             1216
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1217
11541 continue                                                             1218
11542 continue                                                             1218
11531 continue                                                             1219
11532 continue                                                             1219
      deallocate(ww,ju,vq,xm)                                              1219
      if(isd.gt.0) deallocate(xs)                                          1220
      return                                                               1221
      end                                                                  1222
      subroutine lstandard1 (no,ni,x,w,ju,isd,thr,sthr,xm,xs)              1223
      real x(no,ni),w(no),xm(ni),xs(ni)                                    1223
      integer ju(ni)                                                       1224
      ssd=0.0                                                              1224
      nv=0                                                                 1225
11580 do 11581 j=1,ni                                                      1225
      if(ju(j).eq.0)goto 11581                                             1226
      xm(j)=dot_product(w,x(:,j))                                          1226
      x(:,j)=x(:,j)-xm(j)                                                  1227
      xsj=sqrt(dot_product(w,x(:,j)**2))                                   1227
      ssd=ssd+xsj                                                          1227
      nv=nv+1                                                              1228
      if(isd .le. 0)goto 11601                                             1228
      xs(j)=xsj                                                            1228
      x(:,j)=x(:,j)/xsj                                                    1228
11601 continue                                                             1229
11581 continue                                                             1230
11582 continue                                                             1230
      sthr=thr                                                             1230
      if(isd.eq.0) sthr=sthr*ssd/nv                                        1231
      return                                                               1232
      end                                                                  1233
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ulam   1235 
     *,shr,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1237 
     *5, devmax=0.999)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    1238
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1239
      integer ju(ni),m(nx),kin(nlam)                                       1240
      real, dimension (:), allocatable :: b,bs,v,r,xv,q                         
      integer, dimension (:), allocatable :: mm                                 
      allocate(b(0:ni),stat=jerr)                                          1245
      allocate(xv(1:ni),stat=ierr)                                         1245
      jerr=jerr+ierr                                                       1246
      allocate(bs(0:ni),stat=ierr)                                         1246
      jerr=jerr+ierr                                                       1247
      allocate(mm(1:ni),stat=ierr)                                         1247
      jerr=jerr+ierr                                                       1248
      allocate(r(1:no),stat=ierr)                                          1248
      jerr=jerr+ierr                                                       1249
      allocate(v(1:no),stat=ierr)                                          1249
      jerr=jerr+ierr                                                       1250
      allocate(q(1:no),stat=ierr)                                          1250
      jerr=jerr+ierr                                                       1251
      if(jerr.ne.0) return                                                 1252
      fmax=log(1.0/pmin-1.0)                                               1252
      fmin=-fmax                                                           1252
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1253
      bta=max(parm,1.0e-3)                                                 1253
      omb=1.0-bta                                                          1254
      q0=dot_product(w,y)                                                  1254
      if(q0 .gt. pmin)goto 11621                                           1254
      jerr=8001                                                            1254
      return                                                               1254
11621 continue                                                             1255
      if(q0 .lt. 1.0-pmin)goto 11641                                       1255
      jerr=9001                                                            1255
      return                                                               1255
11641 continue                                                             1256
      bz=log(q0/(1.0-q0))                                                  1257
      if(nonzero(no,g) .ne. 0)goto 11661                                   1257
      vi=q0*(1.0-q0)                                                       1257
      b(0)=bz                                                              1257
      v=vi*w                                                               1258
      r=w*(y-q0)                                                           1258
      q=q0                                                                 1258
      xmz=vi                                                               1259
      goto 11671                                                           1260
11661 continue                                                             1260
      b(0)=azero(no,y,g,w,jerr)                                            1260
      if(jerr.ne.0) return                                                 1261
      q=1.0/(1.0+exp(-b(0)-g))                                             1261
      v=w*q*(1.0-q)                                                        1261
      r=w*(y-q)                                                            1261
      xmz=sum(v)                                                           1262
11671 continue                                                             1263
11651 continue                                                             1263
      if(isd .le. 0)goto 11691                                             1263
      xv=0.25                                                              1263
      goto 11701                                                           1264
11691 continue                                                             1264
11710 do 11711 j=1,ni                                                      1264
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1264
11711 continue                                                             1264
11712 continue                                                             1264
11701 continue                                                             1265
11681 continue                                                             1265
      dev1=-(bz*q0+log(1.0-q0))                                            1265
      dev0=dev1                                                            1266
11720 do 11721 i=1,no                                                      1266
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1267
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1268
11721 continue                                                             1269
11722 continue                                                             1269
      if(flmin .ge. 1.0)goto 11741                                         1269
      eqs=max(eps,flmin)                                                   1269
      alf=eqs**(1.0/(nlam-1))                                              1269
11741 continue                                                             1270
      m=0                                                                  1270
      mm=0                                                                 1270
      nlp=0                                                                1270
      nin=nlp                                                              1270
      mnl=min(mnlam,nlam)                                                  1270
      bs=0.0                                                               1270
      b(1:ni)=0.0                                                          1271
11750 do 11751 ilm=1,nlam                                                  1272
      if(flmin .lt. 1.0)goto 11771                                         1272
      al=ulam(ilm)                                                         1272
      goto 11761                                                           1273
11771 if(ilm .le. 2)goto 11781                                             1273
      al=al*alf                                                            1273
      goto 11761                                                           1274
11781 if(ilm .ne. 1)goto 11791                                             1274
      al=big                                                               1274
      goto 11801                                                           1275
11791 continue                                                             1275
      al=0.0                                                               1276
11810 do 11811 j=1,ni                                                      1276
      if(ju(j).eq.0)goto 11811                                             1276
      if(vp(j).le.0.0)goto 11811                                           1277
      al=max(al,abs(dot_product(r,x(:,j)))/vp(j))                          1278
11811 continue                                                             1279
11812 continue                                                             1279
      al=alf*al/bta                                                        1280
11801 continue                                                             1281
11761 continue                                                             1281
      al2=al*omb                                                           1281
      al1=al*bta                                                           1281
      nit=0                                                                1282
11820 continue                                                             1282
11821 continue                                                             1282
      bs(0)=b(0)                                                           1282
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1283
11830 continue                                                             1283
11831 continue                                                             1283
      nlp=nlp+1                                                            1283
      dlx=0.0                                                              1284
11840 do 11841 k=1,ni                                                      1284
      if(ju(k).eq.0)goto 11841                                             1285
      bk=b(k)                                                              1285
      gk=dot_product(r,x(:,k))                                             1286
      u=gk+xv(k)*b(k)                                                      1286
      au=abs(u)-vp(k)*al1                                                  1287
      if(au .gt. 0.0)goto 11861                                            1287
      b(k)=0.0                                                             1287
      goto 11871                                                           1288
11861 continue                                                             1288
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1288
11871 continue                                                             1289
11851 continue                                                             1289
      d=b(k)-bk                                                            1289
      if(abs(d).le.0.0)goto 11841                                          1289
      dlx=max(dlx,xv(k)*abs(d))                                            1290
      r=r-d*v*x(:,k)                                                       1291
      if(mm(k) .ne. 0)goto 11891                                           1291
      nin=nin+1                                                            1291
      if(nin.gt.nx)goto 11842                                              1292
      mm(k)=nin                                                            1292
      m(nin)=k                                                             1293
11891 continue                                                             1294
11841 continue                                                             1295
11842 continue                                                             1295
      if(nin.gt.nx)goto 11832                                              1296
      d=sum(r)/xmz                                                         1297
      if(d .eq. 0.0)goto 11911                                             1297
      b(0)=b(0)+d                                                          1297
      dlx=max(dlx,xmz*abs(d))                                              1297
      r=r-d*v                                                              1297
11911 continue                                                             1298
      if(dlx.lt.shr)goto 11832                                             1299
11920 continue                                                             1299
11921 continue                                                             1299
      nlp=nlp+1                                                            1299
      dlx=0.0                                                              1300
11930 do 11931 l=1,nin                                                     1300
      k=m(l)                                                               1300
      bk=b(k)                                                              1300
      gk=dot_product(r,x(:,k))                                             1301
      u=gk+xv(k)*b(k)                                                      1301
      au=abs(u)-vp(k)*al1                                                  1302
      if(au .gt. 0.0)goto 11951                                            1302
      b(k)=0.0                                                             1302
      goto 11961                                                           1303
11951 continue                                                             1303
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1303
11961 continue                                                             1304
11941 continue                                                             1304
      d=b(k)-bk                                                            1304
      if(abs(d).le.0.0)goto 11931                                          1304
      dlx=max(dlx,xv(k)*abs(d))                                            1305
      r=r-d*v*x(:,k)                                                       1306
11931 continue                                                             1307
11932 continue                                                             1307
      d=sum(r)/xmz                                                         1308
      if(d .eq. 0.0)goto 11981                                             1308
      b(0)=b(0)+d                                                          1308
      dlx=max(dlx,xmz*abs(d))                                              1308
      r=r-d*v                                                              1308
11981 continue                                                             1310
      if(dlx.lt.shr)goto 11922                                             1310
      goto 11921                                                           1311
11922 continue                                                             1311
      goto 11831                                                           1312
11832 continue                                                             1312
      if(nin.gt.nx)goto 11822                                              1313
11990 do 11991 i=1,no                                                      1313
      fi=b(0)+g(i)                                                         1314
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1315
      if(fi .ge. fmin)goto 12011                                           1315
      q(i)=0.0                                                             1315
      goto 12001                                                           1315
12011 if(fi .le. fmax)goto 12021                                           1315
      q(i)=1.0                                                             1315
      goto 12031                                                           1316
12021 continue                                                             1316
      q(i)=1.0/(1.0+exp(-fi))                                              1316
12031 continue                                                             1317
12001 continue                                                             1317
11991 continue                                                             1318
11992 continue                                                             1318
      v=w*q*(1.0-q)                                                        1318
      xmz=sum(v)                                                           1318
      if(xmz.le.vmin)goto 11822                                            1318
      r=w*(y-q)                                                            1319
      if(kopt .ne. 0)goto 12051                                            1320
12060 do 12061 j=1,nin                                                     1320
      xv(m(j))=dot_product(v,x(:,m(j))**2)                                 1320
12061 continue                                                             1321
12062 continue                                                             1321
12051 continue                                                             1322
      if(xmz*abs(b(0)-bs(0)) .ge. shr)goto 12081                           1322
      ix=0                                                                 1323
12090 do 12091 j=1,nin                                                     1323
      if(xv(m(j))*abs(b(m(j))-bs(m(j))).lt.shr)goto 12091                  1323
      ix=1                                                                 1323
      goto 12092                                                           1323
12091 continue                                                             1324
12092 continue                                                             1324
      if(ix.eq.0)goto 11822                                                1325
12081 continue                                                             1326
      nit=nit+1                                                            1326
      if(nit .le. maxit)goto 12111                                         1326
      jerr=-ilm                                                            1326
      return                                                               1326
12111 continue                                                             1327
      goto 11821                                                           1328
11822 continue                                                             1328
      if(nin .le. nx)goto 12131                                            1328
      jerr=-10000-ilm                                                      1328
      goto 11752                                                           1328
12131 continue                                                             1329
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1329
      kin(ilm)=nin                                                         1330
      a0(ilm)=b(0)                                                         1330
      alm(ilm)=al                                                          1330
      lmu=ilm                                                              1331
      devi=dev2(no,w,y,q,pmin)                                             1332
      dev(ilm)=(dev1-devi)/dev0                                            1332
      if(xmz.le.vmin)goto 11752                                            1333
      if(ilm.lt.mnl)goto 11751                                             1333
      if(flmin.ge.1.0)goto 11751                                           1334
      me=0                                                                 1334
12140 do 12141 j=1,nin                                                     1334
      if(a(j,ilm).ne.0.0) me=me+1                                          1334
12141 continue                                                             1334
12142 continue                                                             1334
      if(me.gt.ne)goto 11752                                               1335
      if(dev(ilm).gt.devmax)goto 11752                                     1335
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 11752                             1336
11751 continue                                                             1337
11752 continue                                                             1337
      g=log(q/(1.0-q))                                                     1338
      deallocate(b,bs,v,r,xv,q,mm)                                         1339
      return                                                               1340
      end                                                                  1341
      function dev2(n,w,y,p,pmin)                                          1342
      real w(n),y(n),p(n)                                                  1343
      pmax=1.0-pmin                                                        1343
      s=0.0                                                                1344
12150 do 12151 i=1,n                                                       1344
      pi=min(max(pmin,p(i)),pmax)                                          1345
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1346
12151 continue                                                             1347
12152 continue                                                             1347
      dev2=s                                                               1348
      return                                                               1349
      end                                                                  1350
      function azero(n,y,g,q,jerr)                                         1351
      parameter(eps=1.0e-7)                                                1352
      real y(n),g(n),q(n)                                                  1353
      real, dimension (:), allocatable :: e,p,w                                 
      allocate(e(1:n),stat=jerr)                                           1357
      allocate(p(1:n),stat=ierr)                                           1357
      jerr=jerr+ierr                                                       1358
      allocate(w(1:n),stat=ierr)                                           1358
      jerr=jerr+ierr                                                       1359
      if(jerr.ne.0) return                                                 1360
      az=0.0                                                               1360
      e=exp(-g)                                                            1360
      qy=dot_product(q,y)                                                  1360
      p=1.0/(1.0+e)                                                        1361
12160 continue                                                             1361
12161 continue                                                             1361
      w=q*p*(1.0-p)                                                        1362
      d=(qy-dot_product(q,p))/sum(w)                                       1362
      az=az+d                                                              1362
      if(abs(d).lt.eps)goto 12162                                          1363
      ea0=exp(-az)                                                         1363
      p=1.0/(1.0+ea0*e)                                                    1364
      goto 12161                                                           1365
12162 continue                                                             1365
      azero=az                                                             1366
      deallocate(e,p,w)                                                    1367
      return                                                               1368
      end                                                                  1369
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ul   1371 
     *am,shr,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1373 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam)              1374
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1375
      integer ju(ni),m(nx),kin(nlam)                                       1376
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp                       
      real, dimension (:), allocatable :: di,v,r                                
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is                              
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no),stat=ierr)                                          1387
      jerr=jerr+ierr                                                       1388
      allocate(v(1:no),stat=ierr)                                          1388
      jerr=jerr+ierr                                                       1389
      allocate(mm(1:ni),stat=ierr)                                         1389
      jerr=jerr+ierr                                                       1390
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1390
      jerr=jerr+ierr                                                       1391
      allocate(sxp(1:no),stat=ierr)                                        1391
      jerr=jerr+ierr                                                       1392
      allocate(di(1:no),stat=ierr)                                         1392
      jerr=jerr+ierr                                                       1393
      if(jerr.ne.0) return                                                 1394
      pmax=1.0-pmin                                                        1394
      emin=pmin/pmax                                                       1394
      emax=1.0/emin                                                        1395
      pfm=(1.0+pmin)*pmin                                                  1395
      pfx=(1.0-pmin)*pmax                                                  1395
      vmin=pfm*pmax                                                        1396
      bta=max(parm,1.0e-3)                                                 1396
      omb=1.0-bta                                                          1396
      dev1=0.0                                                             1396
      dev0=0.0                                                             1397
12170 do 12171 ic=1,nc                                                     1397
      q0=dot_product(w,y(:,ic))                                            1398
      if(q0 .gt. pmin)goto 12191                                           1398
      jerr =8000+ic                                                        1398
      return                                                               1398
12191 continue                                                             1399
      if(q0 .lt. 1.0-pmin)goto 12211                                       1399
      jerr =9000+ic                                                        1399
      return                                                               1399
12211 continue                                                             1400
      b(0,ic)=log(q0)                                                      1400
      dev1=dev1-q0*b(0,ic)                                                 1400
      b(1:ni,ic)=0.0                                                       1401
12220 do 12221 i=1,no                                                      1401
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1401
12221 continue                                                             1402
12222 continue                                                             1402
12171 continue                                                             1403
12172 continue                                                             1403
      dev0=dev0+dev1                                                       1404
      if(nonzero(no*nc,g) .ne. 0)goto 12241                                1405
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1405
      sxp=0.0                                                              1406
12250 do 12251 ic=1,nc                                                     1406
      q(:,ic)=exp(b(0,ic))                                                 1406
      sxp=sxp+q(:,ic)                                                      1406
12251 continue                                                             1407
12252 continue                                                             1407
      goto 12261                                                           1408
12241 continue                                                             1408
12270 do 12271 i=1,no                                                      1408
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1408
12271 continue                                                             1408
12272 continue                                                             1408
      sxp=0.0                                                              1409
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1409
      if(jerr.ne.0) return                                                 1410
12280 do 12281 ic=1,nc                                                     1410
      q(:,ic)=exp(b(0,ic)+g(:,ic))                                         1410
      sxp=sxp+q(:,ic)                                                      1410
12281 continue                                                             1411
12282 continue                                                             1411
12261 continue                                                             1412
12231 continue                                                             1412
      if(isd .le. 0)goto 12301                                             1412
      xv=0.25                                                              1412
      goto 12311                                                           1413
12301 continue                                                             1413
12320 do 12321 j=1,ni                                                      1413
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1413
12321 continue                                                             1413
12322 continue                                                             1413
12311 continue                                                             1414
12291 continue                                                             1414
      if(flmin .ge. 1.0)goto 12341                                         1414
      eqs=max(eps,flmin)                                                   1414
      alf=eqs**(1.0/(nlam-1))                                              1414
12341 continue                                                             1415
      m=0                                                                  1415
      mm=0                                                                 1415
      nin=0                                                                1415
      nlp=0                                                                1415
      mnl=min(mnlam,nlam)                                                  1415
      bs=0.0                                                               1416
12350 do 12351 ilm=1,nlam                                                  1417
      if(flmin .lt. 1.0)goto 12371                                         1417
      al=ulam(ilm)                                                         1417
      goto 12361                                                           1418
12371 if(ilm .le. 2)goto 12381                                             1418
      al=al*alf                                                            1418
      goto 12361                                                           1419
12381 if(ilm .ne. 1)goto 12391                                             1419
      al=big                                                               1419
      goto 12401                                                           1420
12391 continue                                                             1420
      al=0.0                                                               1421
12410 do 12411 ic=1,nc                                                     1421
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1422
12420 do 12421 j=1,ni                                                      1422
      if(ju(j).eq.0)goto 12421                                             1422
      if(vp(j).le.0.0)goto 12421                                           1423
      al=max(al,abs(dot_product(r,x(:,j)))/vp(j))                          1424
12421 continue                                                             1425
12422 continue                                                             1425
12411 continue                                                             1426
12412 continue                                                             1426
      al=alf*al/bta                                                        1427
12401 continue                                                             1428
12361 continue                                                             1428
      al2=al*omb                                                           1428
      al1=al*bta                                                           1428
      nit=0                                                                1429
12430 continue                                                             1429
12431 continue                                                             1429
      ix=0                                                                 1429
      jx=ix                                                                1429
      ig=0                                                                 1430
12440 do 12441 ic=1,nc                                                     1430
      bs(0,ic)=b(0,ic)                                                     1431
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1432
      xmz=0.0                                                              1433
12450 do 12451 i=1,no                                                      1433
      pic=q(i,ic)/sxp(i)                                                   1434
      if(pic .ge. pfm)goto 12471                                           1434
      pic=0.0                                                              1434
      v(i)=0.0                                                             1434
      goto 12461                                                           1435
12471 if(pic .le. pfx)goto 12481                                           1435
      pic=1.0                                                              1435
      v(i)=0.0                                                             1435
      goto 12491                                                           1436
12481 continue                                                             1436
      v(i)=w(i)*pic*(1.0-pic)                                              1436
      xmz=xmz+v(i)                                                         1436
12491 continue                                                             1437
12461 continue                                                             1437
      r(i)=w(i)*(y(i,ic)-pic)                                              1438
12451 continue                                                             1439
12452 continue                                                             1439
      if(xmz.le.vmin)goto 12441                                            1439
      ig=1                                                                 1440
      if(kopt .ne. 0)goto 12511                                            1441
12520 do 12521 j=1,nin                                                     1441
      xv(m(j),ic)=dot_product(v,x(:,m(j))**2)                              1441
12521 continue                                                             1442
12522 continue                                                             1442
12511 continue                                                             1443
12530 continue                                                             1443
12531 continue                                                             1443
      nlp=nlp+1                                                            1443
      dlx=0.0                                                              1444
12540 do 12541 k=1,ni                                                      1444
      if(ju(k).eq.0)goto 12541                                             1445
      bk=b(k,ic)                                                           1445
      gk=dot_product(r,x(:,k))                                             1446
      u=gk+xv(k,ic)*b(k,ic)                                                1446
      au=abs(u)-vp(k)*al1                                                  1447
      if(au .gt. 0.0)goto 12561                                            1447
      b(k,ic)=0.0                                                          1447
      goto 12571                                                           1448
12561 continue                                                             1448
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1448
12571 continue                                                             1449
12551 continue                                                             1449
      d=b(k,ic)-bk                                                         1449
      if(abs(d).le.0.0)goto 12541                                          1450
      dlx=max(dlx,xv(k,ic)*abs(d))                                         1450
      r=r-d*v*x(:,k)                                                       1451
      if(mm(k) .ne. 0)goto 12591                                           1451
      nin=nin+1                                                            1452
      if(nin .le. nx)goto 12611                                            1452
      jx=1                                                                 1452
      goto 12542                                                           1452
12611 continue                                                             1453
      mm(k)=nin                                                            1453
      m(nin)=k                                                             1454
12591 continue                                                             1455
12541 continue                                                             1456
12542 continue                                                             1456
      if(jx.gt.0)goto 12532                                                1457
      d=sum(r)/xmz                                                         1458
      if(d .eq. 0.0)goto 12631                                             1458
      b(0,ic)=b(0,ic)+d                                                    1458
      dlx=max(dlx,xmz*abs(d))                                              1458
      r=r-d*v                                                              1458
12631 continue                                                             1459
      if(dlx.lt.shr)goto 12532                                             1460
12640 continue                                                             1460
12641 continue                                                             1460
      nlp=nlp+1                                                            1460
      dlx=0.0                                                              1461
12650 do 12651 l=1,nin                                                     1461
      k=m(l)                                                               1461
      bk=b(k,ic)                                                           1462
      gk=dot_product(r,x(:,k))                                             1463
      u=gk+xv(k,ic)*b(k,ic)                                                1463
      au=abs(u)-vp(k)*al1                                                  1464
      if(au .gt. 0.0)goto 12671                                            1464
      b(k,ic)=0.0                                                          1464
      goto 12681                                                           1465
12671 continue                                                             1465
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1465
12681 continue                                                             1466
12661 continue                                                             1466
      d=b(k,ic)-bk                                                         1466
      if(abs(d).le.0.0)goto 12651                                          1467
      dlx=max(dlx,xv(k,ic)*abs(d))                                         1467
      r=r-d*v*x(:,k)                                                       1468
12651 continue                                                             1469
12652 continue                                                             1469
      d=sum(r)/xmz                                                         1470
      if(d .eq. 0.0)goto 12701                                             1470
      b(0,ic)=b(0,ic)+d                                                    1471
      dlx=max(dlx,xmz*abs(d))                                              1471
      r=r-d*v                                                              1472
12701 continue                                                             1473
      if(dlx.lt.shr)goto 12642                                             1473
      goto 12641                                                           1474
12642 continue                                                             1474
      goto 12531                                                           1475
12532 continue                                                             1475
      if(jx.gt.0)goto 12442                                                1476
      if(xmz*abs(b(0,ic)-bs(0,ic)).gt.shr) ix=1                            1477
      if(ix .ne. 0)goto 12721                                              1478
12730 do 12731 j=1,nin                                                     1479
      if(xv(m(j),ic)*abs(b(m(j),ic)-bs(m(j),ic)) .le. shr)goto 12751       1479
      ix=1                                                                 1479
      goto 12732                                                           1479
12751 continue                                                             1480
12731 continue                                                             1481
12732 continue                                                             1481
12721 continue                                                             1482
12760 do 12761 i=1,no                                                      1482
      fi=b(0,ic)+g(i,ic)                                                   1484
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1485
      fi=min(max(exmn,fi),exmx)                                            1485
      sxp(i)=sxp(i)-q(i,ic)                                                1486
      q(i,ic)=min(max(emin*sxp(i),exp(dble(fi))),emax*sxp(i))              1487
      sxp(i)=sxp(i)+q(i,ic)                                                1488
12761 continue                                                             1489
12762 continue                                                             1489
12441 continue                                                             1490
12442 continue                                                             1490
      s=-sum(b(0,:))/nc                                                    1490
      b(0,:)=b(0,:)+s                                                      1490
      di=s                                                                 1491
12770 do 12771 j=1,nin                                                     1491
      l=m(j)                                                               1492
      if(vp(l) .gt. 0.0)goto 12791                                         1492
      s=sum(b(l,:))/nc                                                     1492
      goto 12801                                                           1493
12791 continue                                                             1493
      s=elc(parm,nc,b(l,:),is)                                             1493
12801 continue                                                             1494
12781 continue                                                             1494
      b(l,:)=b(l,:)-s                                                      1494
      di=di-s*x(:,l)                                                       1495
12771 continue                                                             1496
12772 continue                                                             1496
      di=exp(di)                                                           1496
      sxp=sxp*di                                                           1496
12810 do 12811 ic=1,nc                                                     1496
      q(:,ic)=q(:,ic)*di                                                   1496
12811 continue                                                             1497
12812 continue                                                             1497
      if(jx.gt.0)goto 12432                                                1497
      if(ix.eq.0)goto 12432                                                1497
      if(ig.eq.0)goto 12432                                                1498
      nit=nit+1                                                            1498
      if(nit .le. maxit)goto 12831                                         1498
      jerr=-ilm                                                            1498
      return                                                               1498
12831 continue                                                             1499
      goto 12431                                                           1500
12432 continue                                                             1500
      if(jx .le. 0)goto 12851                                              1500
      jerr=-10000-ilm                                                      1500
      goto 12352                                                           1500
12851 continue                                                             1500
      devi=0.0                                                             1501
12860 do 12861 ic=1,nc                                                     1502
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1502
      a0(ic,ilm)=b(0,ic)                                                   1503
12870 do 12871 i=1,no                                                      1503
      if(y(i,ic).le.0.0)goto 12871                                         1504
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1505
12871 continue                                                             1506
12872 continue                                                             1506
12861 continue                                                             1507
12862 continue                                                             1507
      kin(ilm)=nin                                                         1507
      alm(ilm)=al                                                          1507
      lmu=ilm                                                              1508
      dev(ilm)=(dev1-devi)/dev0                                            1508
      if(ig.eq.0)goto 12352                                                1509
      if(ilm.lt.mnl)goto 12351                                             1509
      if(flmin.ge.1.0)goto 12351                                           1510
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 12352             1511
      if(dev(ilm).gt.devmax)goto 12352                                     1511
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12352                             1512
12351 continue                                                             1513
12352 continue                                                             1513
      g=log(q)                                                             1513
12880 do 12881 i=1,no                                                      1513
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1513
12881 continue                                                             1514
12882 continue                                                             1514
      deallocate(sxp,b,bs,v,r,xv,q,mm,is)                                  1515
      return                                                               1516
      end                                                                  1517
      subroutine kazero(kk,n,y,g,q,az,jerr)                                1518
      parameter(eps=1.0e-7)                                                1519
      real y(n,kk),g(n,kk),q(n),az(kk)                                     1520
      real, dimension (:), allocatable :: s                                     
      real, dimension (:,:), allocatable :: e                                   
      allocate(e(1:n,1:kk),stat=jerr)                                           
      allocate(s(1:n),stat=ierr)                                           1525
      jerr=jerr+ierr                                                       1526
      if(jerr.ne.0) return                                                 1527
      az=0.0                                                               1527
      e=exp(g)                                                             1527
12890 do 12891 i=1,n                                                       1527
      s(i)=sum(e(i,:))                                                     1527
12891 continue                                                             1528
12892 continue                                                             1528
12900 continue                                                             1528
12901 continue                                                             1528
      dm=0.0                                                               1529
12910 do 12911 k=1,kk                                                      1529
      t=0.0                                                                1529
      u=t                                                                  1530
12920 do 12921 i=1,n                                                       1530
      pik=e(i,k)/s(i)                                                      1531
      t=t+q(i)*(y(i,k)-pik)                                                1531
      u=u+q(i)*pik*(1.0-pik)                                               1532
12921 continue                                                             1533
12922 continue                                                             1533
      d=t/u                                                                1533
      az(k)=az(k)+d                                                        1533
      ed=exp(d)                                                            1533
      dm=max(dm,abs(d))                                                    1534
12930 do 12931 i=1,n                                                       1534
      z=e(i,k)                                                             1534
      e(i,k)=z*ed                                                          1534
      s(i)=s(i)-z+e(i,k)                                                   1534
12931 continue                                                             1535
12932 continue                                                             1535
12911 continue                                                             1536
12912 continue                                                             1536
      if(dm.lt.eps)goto 12902                                              1536
      goto 12901                                                           1537
12902 continue                                                             1537
      az=az-sum(az)/kk                                                     1538
      deallocate(e,s)                                                      1539
      return                                                               1540
      end                                                                  1541
      function elc(parm,n,a,m)                                             1542
      real a(n)                                                            1542
      integer m(n)                                                         1543
      fn=n                                                                 1543
      am=sum(a)/fn                                                         1544
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 12951                       1544
      elc=am                                                               1544
      return                                                               1544
12951 continue                                                             1545
12960 do 12961 i=1,n                                                       1545
      m(i)=i                                                               1545
12961 continue                                                             1545
12962 continue                                                             1545
      call psort7(a,m,1,n)                                                 1546
      if(a(m(1)) .ne. a(m(n)))goto 12981                                   1546
      elc=a(1)                                                             1546
      return                                                               1546
12981 continue                                                             1547
      if(mod(n,2) .ne. 1)goto 13001                                        1547
      ad=a(m(n/2+1))                                                       1547
      goto 13011                                                           1548
13001 continue                                                             1548
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       1548
13011 continue                                                             1549
12991 continue                                                             1549
      if(parm .ne. 1.0)goto 13031                                          1549
      elc=ad                                                               1549
      return                                                               1549
13031 continue                                                             1550
      b1=min(am,ad)                                                        1550
      b2=max(am,ad)                                                        1550
      k2=1                                                                 1551
13040 continue                                                             1551
13041 if(a(m(k2)).gt.b1)goto 13042                                         1551
      k2=k2+1                                                              1551
      goto 13041                                                           1551
13042 continue                                                             1551
      k1=k2-1                                                              1552
13050 continue                                                             1552
13051 if(a(m(k2)).ge.b2)goto 13052                                         1552
      k2=k2+1                                                              1552
      goto 13051                                                           1553
13052 continue                                                             1553
      r=parm/((1.0-parm)*fn)                                               1553
      is=0                                                                 1553
      sm=n-2*(k1-1)                                                        1554
13060 do 13061 k=k1,k2-1                                                   1554
      sm=sm-2.0                                                            1554
      s=r*sm+am                                                            1555
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 13081                   1555
      is=k                                                                 1555
      goto 13062                                                           1555
13081 continue                                                             1556
13061 continue                                                             1557
13062 continue                                                             1557
      if(is .eq. 0)goto 13101                                              1557
      elc=s                                                                1557
      return                                                               1557
13101 continue                                                             1557
      r2=2.0*r                                                             1557
      s1=a(m(k1))                                                          1557
      am2=2.0*am                                                           1558
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    1558
      elc=s1                                                               1559
13110 do 13111 k=k1+1,k2                                                   1559
      s=a(m(k))                                                            1559
      if(s.eq.s1)goto 13111                                                1560
      c=r2*sum(abs(a-s))+s*(s-am2)                                         1561
      if(c .ge. cri)goto 13131                                             1561
      cri=c                                                                1561
      elc=s                                                                1561
13131 continue                                                             1561
      s1=s                                                                 1562
13111 continue                                                             1563
13112 continue                                                             1563
      return                                                               1564
      end                                                                  1565
      function nintot(ni,nx,nc,a,m,nin,is)                                 1566
      real a(nx,nc)                                                        1566
      integer m(nx),is(ni)                                                 1567
      is=0                                                                 1567
      nintot=0                                                             1568
13140 do 13141 ic=1,nc                                                     1568
13150 do 13151 j=1,nin                                                     1568
      k=m(j)                                                               1568
      if(is(k).ne.0)goto 13151                                             1569
      if(a(j,ic).eq.0.0)goto 13151                                         1569
      is(k)=k                                                              1569
      nintot=nintot+1                                                      1570
13151 continue                                                             1570
13152 continue                                                             1570
13141 continue                                                             1571
13142 continue                                                             1571
      return                                                               1572
      end                                                                  1573
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             1574
      real ca(nx,nc),a(ni,nc)                                              1574
      integer ia(nx)                                                       1575
      a=0.0                                                                1576
13160 do 13161 ic=1,nc                                                     1576
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            1576
13161 continue                                                             1577
13162 continue                                                             1577
      return                                                               1578
      end                                                                  1579
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      1580
      real a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                             1580
      integer ia(nx)                                                       1581
13170 do 13171 i=1,nt                                                      1581
13180 do 13181 ic=1,nc                                                     1581
      ans(ic,i)=a0(ic)                                                     1583
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   1584 
     *:nin)))
13181 continue                                                             1584
13182 continue                                                             1584
13171 continue                                                             1585
13172 continue                                                             1585
      return                                                               1586
      end                                                                  1587
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,ne,nx,nlam,fl   1589 
     *min,  ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      real x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)                 1590
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1591
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1592
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 13201                                    1596
      jerr=10000                                                           1596
      return                                                               1596
13201 continue                                                             1597
      allocate(ww(1:no),stat=jerr)                                         1598
      allocate(ju(1:ni),stat=ierr)                                         1598
      jerr=jerr+ierr                                                       1599
      allocate(vq(1:ni),stat=ierr)                                         1599
      jerr=jerr+ierr                                                       1600
      allocate(xm(1:ni),stat=ierr)                                         1600
      jerr=jerr+ierr                                                       1601
      allocate(xs(1:ni),stat=ierr)                                         1601
      jerr=jerr+ierr                                                       1602
      if(jerr.ne.0) return                                                 1603
      call spchkvars(no,ni,x,ix,ju)                                        1604
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1605
      if(maxval(ju) .gt. 0)goto 13221                                      1605
      jerr=7777                                                            1605
      return                                                               1605
13221 continue                                                             1606
      vq=max(0.0,vp)                                                       1606
      vq=vq*ni/sum(vq)                                                     1607
13230 do 13231 i=1,no                                                      1607
      ww(i)=sum(y(i,:))                                                    1607
      y(i,:)=y(i,:)/ww(i)                                                  1607
13231 continue                                                             1607
13232 continue                                                             1607
      sw=sum(ww)                                                           1607
      ww=ww/sw                                                             1608
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,thr,sthr,xm,xs)            1609
      if(nc .ne. 1)goto 13251                                              1610
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,ne,nx,n   1612 
     *lam,flmin,  ulam,sthr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,d
     *ev,alm,nlp,jerr)
      goto 13261                                                           1613
13251 continue                                                             1614
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmi   1616 
     *n,ulam,  sthr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,n
     *lp,jerr)
13261 continue                                                             1617
13241 continue                                                             1617
      if(jerr.gt.0) return                                                 1617
      dev0=2.0*sw*dev0                                                     1618
13270 do 13271 k=1,lmu                                                     1618
      nk=nin(k)                                                            1619
13280 do 13281 ic=1,nc                                                     1619
      if(isd .le. 0)goto 13301                                             1619
13310 do 13311 l=1,nk                                                      1619
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1619
13311 continue                                                             1619
13312 continue                                                             1619
13301 continue                                                             1620
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1621
13281 continue                                                             1622
13282 continue                                                             1622
13271 continue                                                             1623
13272 continue                                                             1623
      deallocate(ww,ju,vq,xm,xs)                                           1624
      return                                                               1625
      end                                                                  1626
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,thr,sthr,xm,xs)       1627
      real x(*),w(no),xm(ni),xs(ni)                                        1627
      integer ix(*),jx(*),ju(ni)                                           1628
      ssd=0.0                                                              1628
      nv=0                                                                 1629
13320 do 13321 j=1,ni                                                      1629
      if(ju(j).eq.0)goto 13321                                             1629
      jb=ix(j)                                                             1629
      je=ix(j+1)-1                                                         1630
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1631
      xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2)           1632
      ssd=ssd+xs(j)                                                        1632
      nv=nv+1                                                              1633
13321 continue                                                             1634
13322 continue                                                             1634
      if(isd .ne. 0)goto 13341                                             1634
      xs=1.0                                                               1634
      sthr=thr*ssd/nv                                                      1634
      goto 13351                                                           1634
13341 continue                                                             1634
      sthr=thr                                                             1634
13351 continue                                                             1635
13331 continue                                                             1635
      return                                                               1636
      end                                                                  1637
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam,   1639 
     *  flmin,ulam,shr,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,alm,
     *nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1641 
     *5, devmax=0.999)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        1642
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1643
      real xb(ni),xs(ni)                                                   1643
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1644
      real, dimension (:), allocatable :: xm,b,bs,v,r,sc,xv,q                   
      integer, dimension (:), allocatable :: mm                                 
      allocate(b(0:ni),stat=jerr)                                          1649
      allocate(xm(0:ni),stat=ierr)                                         1649
      jerr=jerr+ierr                                                       1650
      allocate(xv(1:ni),stat=ierr)                                         1650
      jerr=jerr+ierr                                                       1651
      allocate(bs(0:ni),stat=ierr)                                         1651
      jerr=jerr+ierr                                                       1652
      allocate(mm(1:ni),stat=ierr)                                         1652
      jerr=jerr+ierr                                                       1653
      allocate(q(1:no),stat=ierr)                                          1653
      jerr=jerr+ierr                                                       1654
      allocate(r(1:no),stat=ierr)                                          1654
      jerr=jerr+ierr                                                       1655
      allocate(v(1:no),stat=ierr)                                          1655
      jerr=jerr+ierr                                                       1656
      allocate(sc(1:no),stat=ierr)                                         1656
      jerr=jerr+ierr                                                       1657
      if(jerr.ne.0) return                                                 1658
      fmax=log(1.0/pmin-1.0)                                               1658
      fmin=-fmax                                                           1658
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1659
      bta=max(parm,1.0e-3)                                                 1659
      omb=1.0-bta                                                          1660
      q0=dot_product(w,y)                                                  1660
      if(q0 .gt. pmin)goto 13371                                           1660
      jerr=8001                                                            1660
      return                                                               1660
13371 continue                                                             1661
      if(q0 .lt. 1.0-pmin)goto 13391                                       1661
      jerr=9001                                                            1661
      return                                                               1661
13391 continue                                                             1661
      bz=log(q0/(1.0-q0))                                                  1662
      if(nonzero(no,g) .ne. 0)goto 13411                                   1662
      vi=q0*(1.0-q0)                                                       1662
      b(0)=bz                                                              1662
      v=vi*w                                                               1663
      r=w*(y-q0)                                                           1663
      q=q0                                                                 1663
      xm(0)=vi                                                             1664
      goto 13421                                                           1665
13411 continue                                                             1665
      b(0)=azero(no,y,g,w,jerr)                                            1665
      if(jerr.ne.0) return                                                 1666
      q=1.0/(1.0+exp(-b(0)-g))                                             1666
      v=w*q*(1.0-q)                                                        1666
      r=w*(y-q)                                                            1666
      xm(0)=sum(v)                                                         1667
13421 continue                                                             1668
13401 continue                                                             1668
      if(isd .le. 0)goto 13441                                             1668
      xv=0.25                                                              1668
      goto 13451                                                           1669
13441 continue                                                             1670
13460 do 13461 j=1,ni                                                      1670
      if(ju(j).eq.0)goto 13461                                             1670
      jb=ix(j)                                                             1670
      je=ix(j+1)-1                                                         1671
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          1672
13461 continue                                                             1673
13462 continue                                                             1673
13451 continue                                                             1674
13431 continue                                                             1674
      b(1:ni)=0.0                                                          1674
      dev1=-(bz*q0+log(1.0-q0))                                            1674
      dev0=dev1                                                            1675
13470 do 13471 i=1,no                                                      1675
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1676
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1677
13471 continue                                                             1678
13472 continue                                                             1678
      if(flmin .ge. 1.0)goto 13491                                         1678
      eqs=max(eps,flmin)                                                   1678
      alf=eqs**(1.0/(nlam-1))                                              1678
13491 continue                                                             1679
      m=0                                                                  1679
      mm=0                                                                 1679
      nin=0                                                                1679
      o=0.0                                                                1679
      svr=o                                                                1679
      mnl=min(mnlam,nlam)                                                  1679
      bs=0.0                                                               1679
      nlp=0                                                                1679
      nin=nlp                                                              1680
13500 do 13501 ilm=1,nlam                                                  1681
      if(flmin .lt. 1.0)goto 13521                                         1681
      al=ulam(ilm)                                                         1681
      goto 13511                                                           1682
13521 if(ilm .le. 2)goto 13531                                             1682
      al=al*alf                                                            1682
      goto 13511                                                           1683
13531 if(ilm .ne. 1)goto 13541                                             1683
      al=big                                                               1683
      goto 13551                                                           1684
13541 continue                                                             1684
      al=0.0                                                               1685
13560 do 13561 j=1,ni                                                      1685
      if(ju(j).eq.0)goto 13561                                             1685
      if(vp(j).le.0.0)goto 13561                                           1686
      jb=ix(j)                                                             1686
      je=ix(j+1)-1                                                         1686
      jn=ix(j+1)-ix(j)                                                     1687
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1688
      gj=dot_product(sc(1:jn),x(jb:je))                                    1689
      gj=(gj-svr*xb(j))/xs(j)                                              1690
      al=max(al,abs(gj)/vp(j))                                             1691
13561 continue                                                             1692
13562 continue                                                             1692
      al=alf*al/bta                                                        1693
13551 continue                                                             1694
13511 continue                                                             1694
      al2=al*omb                                                           1694
      al1=al*bta                                                           1694
      nit=0                                                                1695
13570 continue                                                             1695
13571 continue                                                             1695
      bs(0)=b(0)                                                           1695
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1696
13580 continue                                                             1696
13581 continue                                                             1696
      nlp=nlp+1                                                            1696
      dlx=0.0                                                              1697
13590 do 13591 k=1,ni                                                      1697
      if(ju(k).eq.0)goto 13591                                             1698
      jb=ix(k)                                                             1698
      je=ix(k+1)-1                                                         1698
      jn=ix(k+1)-ix(k)                                                     1698
      bk=b(k)                                                              1699
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1700
      gk=dot_product(sc(1:jn),x(jb:je))                                    1701
      gk=(gk-svr*xb(k))/xs(k)                                              1702
      u=gk+xv(k)*b(k)                                                      1702
      au=abs(u)-vp(k)*al1                                                  1703
      if(au .gt. 0.0)goto 13611                                            1703
      b(k)=0.0                                                             1703
      goto 13621                                                           1704
13611 continue                                                             1704
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1704
13621 continue                                                             1705
13601 continue                                                             1705
      d=b(k)-bk                                                            1705
      if(abs(d).le.0.0)goto 13591                                          1705
      dlx=max(dlx,xv(k)*abs(d))                                            1706
      if(mm(k) .ne. 0)goto 13641                                           1706
      nin=nin+1                                                            1706
      if(nin.gt.nx)goto 13592                                              1707
      mm(k)=nin                                                            1707
      m(nin)=k                                                             1707
      sc(1:jn)=v(jx(jb:je))                                                1708
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 1709
13641 continue                                                             1710
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1711
      o=o+d*(xb(k)/xs(k))                                                  1712
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1713
13591 continue                                                             1714
13592 continue                                                             1714
      if(nin.gt.nx)goto 13582                                              1715
      d=svr/xm(0)                                                          1716
      if(d .eq. 0.0)goto 13661                                             1716
      b(0)=b(0)+d                                                          1716
      dlx=max(dlx,abs(svr))                                                1716
      r=r-d*v                                                              1716
13661 continue                                                             1717
      svr=svr-d*xm(0)                                                      1717
      if(dlx.lt.shr)goto 13582                                             1718
13670 continue                                                             1718
13671 continue                                                             1718
      nlp=nlp+1                                                            1718
      dlx=0.0                                                              1719
13680 do 13681 l=1,nin                                                     1719
      k=m(l)                                                               1719
      jb=ix(k)                                                             1719
      je=ix(k+1)-1                                                         1720
      jn=ix(k+1)-ix(k)                                                     1720
      bk=b(k)                                                              1721
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1722
      gk=dot_product(sc(1:jn),x(jb:je))                                    1723
      gk=(gk-svr*xb(k))/xs(k)                                              1724
      u=gk+xv(k)*b(k)                                                      1724
      au=abs(u)-vp(k)*al1                                                  1725
      if(au .gt. 0.0)goto 13701                                            1725
      b(k)=0.0                                                             1725
      goto 13711                                                           1726
13701 continue                                                             1726
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1726
13711 continue                                                             1727
13691 continue                                                             1727
      d=b(k)-bk                                                            1727
      if(abs(d).le.0.0)goto 13681                                          1727
      dlx=max(dlx,xv(k)*abs(d))                                            1728
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1729
      o=o+d*(xb(k)/xs(k))                                                  1730
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1731
13681 continue                                                             1732
13682 continue                                                             1732
      d=svr/xm(0)                                                          1733
      if(d .eq. 0.0)goto 13731                                             1733
      b(0)=b(0)+d                                                          1733
      dlx=max(dlx,abs(svr))                                                1733
      r=r-d*v                                                              1733
13731 continue                                                             1734
      svr=svr-d*xm(0)                                                      1735
      if(dlx.lt.shr)goto 13672                                             1735
      goto 13671                                                           1736
13672 continue                                                             1736
      goto 13581                                                           1737
13582 continue                                                             1737
      if(nin.gt.nx)goto 13572                                              1738
      sc=b(0)                                                              1738
      b0=0.0                                                               1739
13740 do 13741 j=1,nin                                                     1739
      l=m(j)                                                               1739
      jb=ix(l)                                                             1739
      je=ix(l+1)-1                                                         1740
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      1741
      b0=b0-b(l)*xb(l)/xs(l)                                               1742
13741 continue                                                             1743
13742 continue                                                             1743
      sc=sc+b0                                                             1744
13750 do 13751 i=1,no                                                      1744
      fi=sc(i)+g(i)                                                        1745
      if(fi .ge. fmin)goto 13771                                           1745
      q(i)=0.0                                                             1745
      goto 13761                                                           1745
13771 if(fi .le. fmax)goto 13781                                           1745
      q(i)=1.0                                                             1745
      goto 13791                                                           1746
13781 continue                                                             1746
      q(i)=1.0/(1.0+exp(-fi))                                              1746
13791 continue                                                             1747
13761 continue                                                             1747
13751 continue                                                             1748
13752 continue                                                             1748
      v=w*q*(1.0-q)                                                        1748
      xm(0)=sum(v)                                                         1748
      if(xm(0).lt.vmin)goto 13572                                          1749
      r=w*(y-q)                                                            1749
      svr=sum(r)                                                           1749
      o=0.0                                                                1750
13800 do 13801 l=1,nin                                                     1750
      j=m(l)                                                               1751
      jb=ix(j)                                                             1751
      je=ix(j+1)-1                                                         1751
      jn=ix(j+1)-ix(j)                                                     1752
      sc(1:jn)=v(jx(jb:je))                                                1753
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 1754
      if(kopt .ne. 0)goto 13821                                            1755
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              1756
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                1757
13821 continue                                                             1758
13801 continue                                                             1759
13802 continue                                                             1759
      if(xm(0)*abs(b(0)-bs(0)) .ge. shr)goto 13841                         1759
      kx=0                                                                 1760
13850 do 13851 j=1,nin                                                     1760
      if(xv(m(j))*abs(b(m(j))-bs(m(j))).lt.shr)goto 13851                  1760
      kx=1                                                                 1760
      goto 13852                                                           1760
13851 continue                                                             1761
13852 continue                                                             1761
      if(kx.eq.0)goto 13572                                                1762
13841 continue                                                             1763
      nit=nit+1                                                            1763
      if(nit .le. maxit)goto 13871                                         1763
      jerr=-ilm                                                            1763
      return                                                               1763
13871 continue                                                             1764
      goto 13571                                                           1765
13572 continue                                                             1765
      if(nin .le. nx)goto 13891                                            1765
      jerr=-10000-ilm                                                      1765
      goto 13502                                                           1765
13891 continue                                                             1766
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1766
      kin(ilm)=nin                                                         1767
      a0(ilm)=b(0)                                                         1767
      alm(ilm)=al                                                          1767
      lmu=ilm                                                              1768
      devi=dev2(no,w,y,q,pmin)                                             1769
      dev(ilm)=(dev1-devi)/dev0                                            1770
      if(ilm.lt.mnl)goto 13501                                             1770
      if(flmin.ge.1.0)goto 13501                                           1771
      me=0                                                                 1771
13900 do 13901 j=1,nin                                                     1771
      if(a(j,ilm).ne.0.0) me=me+1                                          1771
13901 continue                                                             1771
13902 continue                                                             1771
      if(me.gt.ne)goto 13502                                               1772
      if(dev(ilm).gt.devmax)goto 13502                                     1772
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13502                             1773
      if(xm(0).lt.vmin)goto 13502                                          1774
13501 continue                                                             1775
13502 continue                                                             1775
      g=log(q/(1.0-q))                                                     1776
      deallocate(xm,b,bs,v,r,sc,xv,q,mm)                                   1777
      return                                                               1778
      end                                                                  1779
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam   1781 
     *,flmin,ulam,  shr,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,alm
     *,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1783 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni)    1784
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1785
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1786
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp                       
      real, dimension (:), allocatable :: sc,xm,v,r                             
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is                              
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(xm(0:ni),stat=ierr)                                         1797
      jerr=jerr+ierr                                                       1798
      allocate(r(1:no),stat=ierr)                                          1798
      jerr=jerr+ierr                                                       1799
      allocate(v(1:no),stat=ierr)                                          1799
      jerr=jerr+ierr                                                       1800
      allocate(mm(1:ni),stat=ierr)                                         1800
      jerr=jerr+ierr                                                       1801
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1801
      jerr=jerr+ierr                                                       1802
      allocate(sxp(1:no),stat=ierr)                                        1802
      jerr=jerr+ierr                                                       1803
      allocate(sc(1:no),stat=ierr)                                         1803
      jerr=jerr+ierr                                                       1804
      if(jerr.ne.0) return                                                 1805
      pmax=1.0-pmin                                                        1805
      emin=pmin/pmax                                                       1805
      emax=1.0/emin                                                        1806
      pfm=(1.0+pmin)*pmin                                                  1806
      pfx=(1.0-pmin)*pmax                                                  1806
      vmin=pfm*pmax                                                        1807
      bta=max(parm,1.0e-3)                                                 1807
      omb=1.0-bta                                                          1807
      dev1=0.0                                                             1807
      dev0=0.0                                                             1808
13910 do 13911 ic=1,nc                                                     1808
      q0=dot_product(w,y(:,ic))                                            1809
      if(q0 .gt. pmin)goto 13931                                           1809
      jerr =8000+ic                                                        1809
      return                                                               1809
13931 continue                                                             1810
      if(q0 .lt. 1.0-pmin)goto 13951                                       1810
      jerr =9000+ic                                                        1810
      return                                                               1810
13951 continue                                                             1811
      b(1:ni,ic)=0.0                                                       1811
      b(0,ic)=log(q0)                                                      1811
      dev1=dev1-q0*b(0,ic)                                                 1812
13960 do 13961 i=1,no                                                      1812
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1812
13961 continue                                                             1813
13962 continue                                                             1813
13911 continue                                                             1814
13912 continue                                                             1814
      dev0=dev0+dev1                                                       1815
      if(nonzero(no*nc,g) .ne. 0)goto 13981                                1816
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1816
      sxp=0.0                                                              1817
13990 do 13991 ic=1,nc                                                     1817
      q(:,ic)=exp(b(0,ic))                                                 1817
      sxp=sxp+q(:,ic)                                                      1817
13991 continue                                                             1818
13992 continue                                                             1818
      goto 14001                                                           1819
13981 continue                                                             1819
14010 do 14011 i=1,no                                                      1819
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1819
14011 continue                                                             1819
14012 continue                                                             1819
      sxp=0.0                                                              1820
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1820
      if(jerr.ne.0) return                                                 1821
14020 do 14021 ic=1,nc                                                     1821
      q(:,ic)=exp(b(0,ic)+g(:,ic))                                         1821
      sxp=sxp+q(:,ic)                                                      1821
14021 continue                                                             1822
14022 continue                                                             1822
14001 continue                                                             1823
13971 continue                                                             1823
      if(isd .le. 0)goto 14041                                             1823
      xv=0.25                                                              1823
      goto 14051                                                           1824
14041 continue                                                             1825
14060 do 14061 j=1,ni                                                      1825
      if(ju(j).eq.0)goto 14061                                             1825
      jb=ix(j)                                                             1825
      je=ix(j+1)-1                                                         1826
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        1827
14061 continue                                                             1828
14062 continue                                                             1828
14051 continue                                                             1829
14031 continue                                                             1829
      if(flmin .ge. 1.0)goto 14081                                         1829
      eqs=max(eps,flmin)                                                   1829
      alf=eqs**(1.0/(nlam-1))                                              1829
14081 continue                                                             1830
      m=0                                                                  1830
      mm=0                                                                 1830
      nin=0                                                                1830
      nlp=0                                                                1830
      mnl=min(mnlam,nlam)                                                  1830
      bs=0.0                                                               1830
      svr=0.0                                                              1830
      o=0.0                                                                1831
14090 do 14091 ilm=1,nlam                                                  1832
      if(flmin .lt. 1.0)goto 14111                                         1832
      al=ulam(ilm)                                                         1832
      goto 14101                                                           1833
14111 if(ilm .le. 2)goto 14121                                             1833
      al=al*alf                                                            1833
      goto 14101                                                           1834
14121 if(ilm .ne. 1)goto 14131                                             1834
      al=big                                                               1834
      goto 14141                                                           1835
14131 continue                                                             1835
      al=0.0                                                               1836
14150 do 14151 ic=1,nc                                                     1836
      v=q(:,ic)/sxp                                                        1836
      r=w*(y(:,ic)-v)                                                      1836
      v=w*v*(1.0-v)                                                        1837
14160 do 14161 j=1,ni                                                      1837
      if(ju(j).eq.0)goto 14161                                             1837
      if(vp(j).le.0.0)goto 14161                                           1838
      jb=ix(j)                                                             1838
      je=ix(j+1)-1                                                         1838
      jn=ix(j+1)-ix(j)                                                     1839
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1840
      gj=dot_product(sc(1:jn),x(jb:je))                                    1841
      gj=(gj-svr*xb(j))/xs(j)                                              1842
      al=max(al,abs(gj)/vp(j))                                             1843
14161 continue                                                             1844
14162 continue                                                             1844
14151 continue                                                             1845
14152 continue                                                             1845
      al=alf*al/bta                                                        1846
14141 continue                                                             1847
14101 continue                                                             1847
      al2=al*omb                                                           1847
      al1=al*bta                                                           1847
      nit=0                                                                1848
14170 continue                                                             1848
14171 continue                                                             1848
      ixx=0                                                                1848
      jxx=ixx                                                              1848
      ig=0                                                                 1849
14180 do 14181 ic=1,nc                                                     1849
      bs(0,ic)=b(0,ic)                                                     1850
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1851
      xm(0)=0.0                                                            1851
      svr=0.0                                                              1851
      o=0.0                                                                1852
14190 do 14191 i=1,no                                                      1852
      pic=q(i,ic)/sxp(i)                                                   1853
      if(pic .ge. pfm)goto 14211                                           1853
      pic=0.0                                                              1853
      v(i)=0.0                                                             1853
      goto 14201                                                           1854
14211 if(pic .le. pfx)goto 14221                                           1854
      pic=1.0                                                              1854
      v(i)=0.0                                                             1854
      goto 14231                                                           1855
14221 continue                                                             1855
      v(i)=w(i)*pic*(1.0-pic)                                              1855
      xm(0)=xm(0)+v(i)                                                     1855
14231 continue                                                             1856
14201 continue                                                             1856
      r(i)=w(i)*(y(i,ic)-pic)                                              1856
      svr=svr+r(i)                                                         1857
14191 continue                                                             1858
14192 continue                                                             1858
      if(xm(0).le.vmin)goto 14181                                          1858
      ig=1                                                                 1859
14240 do 14241 l=1,nin                                                     1859
      j=m(l)                                                               1860
      jb=ix(j)                                                             1860
      je=ix(j+1)-1                                                         1861
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             1862
      if(kopt .ne. 0)goto 14261                                            1863
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       1864
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          1865
14261 continue                                                             1866
14241 continue                                                             1867
14242 continue                                                             1867
14270 continue                                                             1867
14271 continue                                                             1867
      nlp=nlp+1                                                            1867
      dlx=0.0                                                              1868
14280 do 14281 k=1,ni                                                      1868
      if(ju(k).eq.0)goto 14281                                             1869
      jb=ix(k)                                                             1869
      je=ix(k+1)-1                                                         1869
      jn=ix(k+1)-ix(k)                                                     1869
      bk=b(k,ic)                                                           1870
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1871
      gk=dot_product(sc(1:jn),x(jb:je))                                    1872
      gk=(gk-svr*xb(k))/xs(k)                                              1873
      u=gk+xv(k,ic)*b(k,ic)                                                1873
      au=abs(u)-vp(k)*al1                                                  1874
      if(au .gt. 0.0)goto 14301                                            1874
      b(k,ic)=0.0                                                          1874
      goto 14311                                                           1875
14301 continue                                                             1875
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1875
14311 continue                                                             1876
14291 continue                                                             1876
      d=b(k,ic)-bk                                                         1876
      if(abs(d).le.0.0)goto 14281                                          1877
      dlx=max(dlx,xv(k,ic)*abs(d))                                         1878
      if(mm(k) .ne. 0)goto 14331                                           1878
      nin=nin+1                                                            1879
      if(nin .le. nx)goto 14351                                            1879
      jxx=1                                                                1879
      goto 14282                                                           1879
14351 continue                                                             1880
      mm(k)=nin                                                            1880
      m(nin)=k                                                             1881
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             1882
14331 continue                                                             1883
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1884
      o=o+d*(xb(k)/xs(k))                                                  1885
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1886
14281 continue                                                             1887
14282 continue                                                             1887
      if(jxx.gt.0)goto 14272                                               1888
      d=svr/xm(0)                                                          1889
      if(d .eq. 0.0)goto 14371                                             1889
      b(0,ic)=b(0,ic)+d                                                    1889
      dlx=max(dlx,abs(svr))                                                1890
      r=r-d*v                                                              1890
      svr=svr-d*xm(0)                                                      1891
14371 continue                                                             1892
      if(dlx.lt.shr)goto 14272                                             1893
14380 continue                                                             1893
14381 continue                                                             1893
      nlp=nlp+1                                                            1893
      dlx=0.0                                                              1894
14390 do 14391 l=1,nin                                                     1894
      k=m(l)                                                               1894
      jb=ix(k)                                                             1894
      je=ix(k+1)-1                                                         1895
      jn=ix(k+1)-ix(k)                                                     1895
      bk=b(k,ic)                                                           1896
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1897
      gk=dot_product(sc(1:jn),x(jb:je))                                    1898
      gk=(gk-svr*xb(k))/xs(k)                                              1899
      u=gk+xv(k,ic)*b(k,ic)                                                1899
      au=abs(u)-vp(k)*al1                                                  1900
      if(au .gt. 0.0)goto 14411                                            1900
      b(k,ic)=0.0                                                          1900
      goto 14421                                                           1901
14411 continue                                                             1901
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1901
14421 continue                                                             1902
14401 continue                                                             1902
      d=b(k,ic)-bk                                                         1902
      if(abs(d).le.0.0)goto 14391                                          1903
      dlx=max(dlx,xv(k,ic)*abs(d))                                         1904
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1905
      o=o+d*(xb(k)/xs(k))                                                  1906
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1907
14391 continue                                                             1908
14392 continue                                                             1908
      d=svr/xm(0)                                                          1909
      if(d .eq. 0.0)goto 14441                                             1909
      b(0,ic)=b(0,ic)+d                                                    1909
      dlx=max(dlx,abs(svr))                                                1910
      r=r-d*v                                                              1910
      svr=svr-d*xm(0)                                                      1911
14441 continue                                                             1912
      if(dlx.lt.shr)goto 14382                                             1912
      goto 14381                                                           1913
14382 continue                                                             1913
      goto 14271                                                           1914
14272 continue                                                             1914
      if(jxx.gt.0)goto 14182                                               1915
      if(xm(0)*abs(b(0,ic)-bs(0,ic)).gt.shr) ixx=1                         1916
      if(ixx .ne. 0)goto 14461                                             1917
14470 do 14471 j=1,nin                                                     1918
      if(xv(m(j),ic)*abs(b(m(j),ic)-bs(m(j),ic)) .le. shr)goto 14491       1918
      ixx=1                                                                1918
      goto 14472                                                           1918
14491 continue                                                             1919
14471 continue                                                             1920
14472 continue                                                             1920
14461 continue                                                             1921
      sc=b(0,ic)+g(:,ic)                                                   1921
      b0=0.0                                                               1922
14500 do 14501 j=1,nin                                                     1922
      l=m(j)                                                               1922
      jb=ix(l)                                                             1922
      je=ix(l+1)-1                                                         1923
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   1924
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            1925
14501 continue                                                             1926
14502 continue                                                             1926
      sc=min(max(exmn,sc+b0),exmx)                                         1927
      sxp=sxp-q(:,ic)                                                      1928
      q(:,ic)=min(max(emin*sxp,exp(dble(sc))),emax*sxp)                    1929
      sxp=sxp+q(:,ic)                                                      1930
14181 continue                                                             1931
14182 continue                                                             1931
      s=-sum(b(0,:))/nc                                                    1931
      b(0,:)=b(0,:)+s                                                      1931
      sc=s                                                                 1931
      b0=0.0                                                               1932
14510 do 14511 j=1,nin                                                     1932
      l=m(j)                                                               1933
      if(vp(l) .gt. 0.0)goto 14531                                         1933
      s=sum(b(l,:))/nc                                                     1933
      goto 14541                                                           1934
14531 continue                                                             1934
      s=elc(parm,nc,b(l,:),is)                                             1934
14541 continue                                                             1935
14521 continue                                                             1935
      b(l,:)=b(l,:)-s                                                      1936
      jb=ix(l)                                                             1936
      je=ix(l+1)-1                                                         1937
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         1938
      b0=b0+s*xb(l)/xs(l)                                                  1939
14511 continue                                                             1940
14512 continue                                                             1940
      sc=sc+b0                                                             1940
      sc=exp(sc)                                                           1940
      sxp=sxp*sc                                                           1940
14550 do 14551 ic=1,nc                                                     1940
      q(:,ic)=q(:,ic)*sc                                                   1940
14551 continue                                                             1941
14552 continue                                                             1941
      if(jxx.gt.0)goto 14172                                               1941
      if(ixx.eq.0)goto 14172                                               1941
      if(ig.eq.0)goto 14172                                                1942
      nit=nit+1                                                            1942
      if(nit .le. maxit)goto 14571                                         1942
      jerr=-ilm                                                            1942
      return                                                               1942
14571 continue                                                             1943
      goto 14171                                                           1944
14172 continue                                                             1944
      if(jxx .le. 0)goto 14591                                             1944
      jerr=-10000-ilm                                                      1944
      goto 14092                                                           1944
14591 continue                                                             1944
      devi=0.0                                                             1945
14600 do 14601 ic=1,nc                                                     1946
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1946
      a0(ic,ilm)=b(0,ic)                                                   1947
14610 do 14611 i=1,no                                                      1947
      if(y(i,ic).le.0.0)goto 14611                                         1948
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1949
14611 continue                                                             1950
14612 continue                                                             1950
14601 continue                                                             1951
14602 continue                                                             1951
      kin(ilm)=nin                                                         1951
      alm(ilm)=al                                                          1951
      lmu=ilm                                                              1952
      dev(ilm)=(dev1-devi)/dev0                                            1952
      if(ig.eq.0)goto 14092                                                1953
      if(ilm.lt.mnl)goto 14091                                             1953
      if(flmin.ge.1.0)goto 14091                                           1954
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 14092             1955
      if(dev(ilm).gt.devmax)goto 14092                                     1955
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 14092                             1956
14091 continue                                                             1957
14092 continue                                                             1957
      g=log(q)                                                             1957
14620 do 14621 i=1,no                                                      1957
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1957
14621 continue                                                             1958
14622 continue                                                             1958
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc)                            1959
      return                                                               1960
      end                                                                  1961
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  1962
      real a0(nc),ca(nx,nc),x(*),f(nc,n)                                   1962
      integer ia(*),ix(*),jx(*)                                            1963
14630 do 14631 ic=1,nc                                                     1963
      f(ic,:)=a0(ic)                                                       1963
14631 continue                                                             1964
14632 continue                                                             1964
14640 do 14641 j=1,nin                                                     1964
      k=ia(j)                                                              1964
      kb=ix(k)                                                             1964
      ke=ix(k+1)-1                                                         1965
14650 do 14651 ic=1,nc                                                     1965
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    1965
14651 continue                                                             1966
14652 continue                                                             1966
14641 continue                                                             1967
14642 continue                                                             1967
      return                                                               1968
      end                                                                  1969
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,ne,nx,nlam,flmin,ula   1971 
     *m,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam)              1972
      real ca(nx,nlam),dev(nlam),alm(nlam)                                 1973
      integer jd(*),ia(nx),nin(nlam)                                       1974
      real, dimension (:), allocatable :: xs,ww,vq                              
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 14671                                    1978
      jerr=10000                                                           1978
      return                                                               1978
14671 continue                                                             1979
      allocate(ww(1:no),stat=jerr)                                         1980
      allocate(ju(1:ni),stat=ierr)                                         1980
      jerr=jerr+ierr                                                       1981
      allocate(vq(1:ni),stat=ierr)                                         1981
      jerr=jerr+ierr                                                       1982
      if(isd .le. 0)goto 14691                                             1982
      allocate(xs(1:ni),stat=ierr)                                         1982
      jerr=jerr+ierr                                                       1982
14691 continue                                                             1983
      if(jerr.ne.0) return                                                 1984
      call chkvars(no,ni,x,ju)                                             1985
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1986
      if(maxval(ju) .gt. 0)goto 14711                                      1986
      jerr=7777                                                            1986
      return                                                               1986
14711 continue                                                             1987
      vq=max(0.0,vp)                                                       1987
      vq=vq*ni/sum(vq)                                                     1988
      ww=max(0.0,w)                                                        1988
      sw=sum(ww)                                                           1989
      if(sw .gt. 0.0)goto 14731                                            1989
      jerr=9999                                                            1989
      return                                                               1989
14731 continue                                                             1989
      ww=ww/sw                                                             1990
      call cstandard(no,ni,x,ww,ju,isd,thr,sthr,xs)                        1991
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,sth   1993 
     *r,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1993
      dev0=2.0*sw*dev0                                                     1994
      if(isd .le. 0)goto 14751                                             1994
14760 do 14761 k=1,lmu                                                     1994
      nk=nin(k)                                                            1994
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   1994
14761 continue                                                             1994
14762 continue                                                             1994
14751 continue                                                             1995
      deallocate(ww,ju,vq)                                                 1995
      if(isd.gt.0) deallocate(xs)                                          1996
      return                                                               1997
      end                                                                  1998
      subroutine cstandard (no,ni,x,w,ju,isd,thr,sthr,xs)                  1999
      real x(no,ni),w(no),xs(ni)                                           1999
      integer ju(ni)                                                       2000
      ssd=0.0                                                              2000
      nv=0                                                                 2001
14770 do 14771 j=1,ni                                                      2001
      if(ju(j).eq.0)goto 14771                                             2002
      xm=dot_product(w,x(:,j))                                             2002
      x(:,j)=x(:,j)-xm                                                     2003
      xsj=sqrt(dot_product(w,x(:,j)**2))                                   2004
      ssd=ssd+xsj                                                          2004
      nv=nv+1                                                              2005
      if(isd .le. 0)goto 14791                                             2005
      xs(j)=xsj                                                            2005
      x(:,j)=x(:,j)/xsj                                                    2005
14791 continue                                                             2006
14771 continue                                                             2007
14772 continue                                                             2007
      sthr=thr                                                             2007
      if(isd.eq.0) sthr=sthr*ssd/nv                                        2008
      return                                                               2009
      end                                                                  2010
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,ne,nx,nlam,flmin,ula   2012 
     *m,cthr,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=0.001, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99)   2013
      real x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam)              2014
      real ao(nx,nlam),dev(nlam),alm(nlam)                                 2015
      integer ju(ni),m(nx),kin(nlam)                                       2016
      real, dimension (:), allocatable :: w,dk,v,xs,wr,a,as,f,dq                
      real, dimension (:), allocatable :: e,uu                                  
      integer, dimension (:), allocatable :: jp,kp,mm                           
      allocate(e(1:no),stat=jerr)                                          2022
      allocate(uu(1:no),stat=ierr)                                         2022
      jerr=jerr+ierr                                                       2023
      allocate(f(1:no),stat=ierr)                                          2023
      jerr=jerr+ierr                                                       2024
      allocate(w(1:no),stat=ierr)                                          2024
      jerr=jerr+ierr                                                       2025
      allocate(v(1:ni),stat=ierr)                                          2025
      jerr=jerr+ierr                                                       2026
      allocate(a(1:ni),stat=ierr)                                          2026
      jerr=jerr+ierr                                                       2027
      allocate(as(1:ni),stat=ierr)                                         2027
      jerr=jerr+ierr                                                       2028
      allocate(xs(1:ni),stat=ierr)                                         2028
      jerr=jerr+ierr                                                       2029
      allocate(jp(1:no),stat=ierr)                                         2029
      jerr=jerr+ierr                                                       2030
      allocate(kp(1:no),stat=ierr)                                         2030
      jerr=jerr+ierr                                                       2031
      allocate(dk(1:no),stat=ierr)                                         2031
      jerr=jerr+ierr                                                       2032
      allocate(wr(1:no),stat=ierr)                                         2032
      jerr=jerr+ierr                                                       2033
      allocate(dq(1:no),stat=ierr)                                         2033
      jerr=jerr+ierr                                                       2034
      allocate(mm(1:ni),stat=ierr)                                         2034
      jerr=jerr+ierr                                                       2035
      if(jerr.ne.0)go to 11400                                             2036
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2037
      if(jerr.ne.0) go to 11400                                            2037
      alpha=max(parm,1.0e-3)                                               2038
      oma=1.0-alpha                                                        2038
      nlm=0                                                                2039
      dq=d*q                                                               2039
      call died(no,nk,dq,kp,jp,dk)                                         2040
      a=0.0                                                                2040
      f=0.0                                                                2040
      e=q                                                                  2040
      fmax=log(huge(f(1))*0.1)                                             2041
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 2042
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         2042
      dev0=rr                                                              2043
14800 do 14801 i=1,no                                                      2043
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 14821                   2043
      w(i)=0.0                                                             2043
      wr(i)=w(i)                                                           2043
14821 continue                                                             2043
14801 continue                                                             2044
14802 continue                                                             2044
      if(nonzero(no,g) .eq. 0)goto 14841                                   2044
      f=g-dot_product(q,g)                                                 2045
      e=q*exp(sign(min(abs(f),fmax),f))                                    2046
14841 continue                                                             2047
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2048
      if(jerr.ne.0) go to 11400                                            2049
      call vars(no,ni,x,w,v)                                               2050
      if(flmin .ge. 1.0)goto 14861                                         2050
      eqs=max(eps,flmin)                                                   2050
      alf=eqs**(1.0/(nlam-1))                                              2050
14861 continue                                                             2051
      m=0                                                                  2051
      mm=0                                                                 2051
      nlp=0                                                                2051
      nin=nlp                                                              2051
      mnl=min(mnlam,nlam)                                                  2051
      as=0.0                                                               2052
14870 do 14871 ilm=1,nlam                                                  2053
      if(flmin .lt. 1.0)goto 14891                                         2053
      al=ulam(ilm)                                                         2053
      goto 14881                                                           2054
14891 if(ilm .le. 2)goto 14901                                             2054
      al=al*alf                                                            2054
      goto 14881                                                           2055
14901 if(ilm .ne. 1)goto 14911                                             2055
      al=big                                                               2055
      goto 14921                                                           2056
14911 continue                                                             2056
      al=0.0                                                               2057
14930 do 14931 j=1,ni                                                      2057
      if(ju(j).eq.0)goto 14931                                             2057
      if(vp(j).le.0.0)goto 14931                                           2058
      al=max(al,abs(dot_product(wr,x(:,j)))/vp(j))                         2059
14931 continue                                                             2060
14932 continue                                                             2060
      al=alf*al/alpha                                                      2061
14921 continue                                                             2062
14881 continue                                                             2062
      sa=alpha*al                                                          2062
      omal=oma*al                                                          2062
      nit=0                                                                2063
14940 do 14941 ito=1,maxit                                                 2063
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2064
14950 do 14951 iti=1,maxit                                                 2064
      nlp=nlp+1                                                            2064
      dli=0.0                                                              2065
14960 do 14961 j=1,ni                                                      2065
      if(ju(j).eq.0)goto 14961                                             2066
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2067
      if(abs(u) .gt. vp(j)*sa)goto 14981                                   2067
      at=0.0                                                               2067
      goto 14991                                                           2068
14981 continue                                                             2068
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2068
14991 continue                                                             2069
14971 continue                                                             2069
      if(at .eq. a(j))goto 15011                                           2069
      del=at-a(j)                                                          2069
      a(j)=at                                                              2069
      dli=max(dli,v(j)*abs(del))                                           2070
      wr=wr-del*w*x(:,j)                                                   2070
      f=f+del*x(:,j)                                                       2071
      if(mm(j) .ne. 0)goto 15031                                           2071
      nin=nin+1                                                            2071
      if(nin.gt.nx)goto 14962                                              2072
      mm(j)=nin                                                            2072
      m(nin)=j                                                             2073
15031 continue                                                             2074
15011 continue                                                             2075
14961 continue                                                             2076
14962 continue                                                             2076
      if(nin.gt.nx)goto 14952                                              2076
      if(dli.lt.cthr)goto 14952                                            2077
15040 do 15041 ita=1,maxit                                                 2077
      nlp=nlp+1                                                            2077
      dli=0.0                                                              2078
15050 do 15051 l=1,nin                                                     2078
      j=m(l)                                                               2079
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2080
      if(abs(u) .gt. vp(j)*sa)goto 15071                                   2080
      at=0.0                                                               2080
      goto 15081                                                           2081
15071 continue                                                             2081
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2081
15081 continue                                                             2082
15061 continue                                                             2082
      if(at .eq. a(j))goto 15101                                           2082
      del=at-a(j)                                                          2082
      a(j)=at                                                              2082
      dli=max(dli,v(j)*abs(del))                                           2083
      wr=wr-del*w*x(:,j)                                                   2083
      f=f+del*x(:,j)                                                       2084
15101 continue                                                             2085
15051 continue                                                             2086
15052 continue                                                             2086
      if(dli.lt.cthr)goto 15042                                            2087
15041 continue                                                             2088
15042 continue                                                             2088
      if(dli .lt. cthr)goto 15121                                          2088
      jerr=-1                                                              2088
      go to 11400                                                          2088
15121 continue                                                             2089
14951 continue                                                             2090
14952 continue                                                             2090
      if(nin.gt.nx)goto 14942                                              2090
      if(dli .lt. cthr)goto 15141                                          2090
      jerr=-1                                                              2090
      go to 11400                                                          2090
15141 continue                                                             2091
      e=q*exp(sign(min(abs(f),fmax),f))                                    2092
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2093
      if(jerr.ne.0) go to 11400                                            2094
      call vars(no,ni,x,w,v)                                               2095
      ix=0                                                                 2096
15150 do 15151 j=1,nin                                                     2096
      if(v(m(j))*abs(a(m(j))-as(m(j))).lt.cthr)goto 15151                  2096
      ix=1                                                                 2096
      goto 15152                                                           2096
15151 continue                                                             2097
15152 continue                                                             2097
      if(ix.eq.0)goto 14942                                                2098
14941 continue                                                             2099
14942 continue                                                             2099
      if(ix .eq. 0)goto 15171                                              2099
      jerr=-2                                                              2099
      goto 14872                                                           2099
15171 continue                                                             2100
      if(nin .le. nx)goto 15191                                            2100
      jerr=-10000-ilm                                                      2100
      goto 14872                                                           2100
15191 continue                                                             2101
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               2101
      kin(ilm)=nin                                                         2102
      alm(ilm)=al                                                          2102
      lmu=ilm                                                              2103
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   2104
      if(ilm.lt.mnl)goto 14871                                             2104
      if(flmin.ge.1.0)goto 14871                                           2105
      me=0                                                                 2105
15200 do 15201 j=1,nin                                                     2105
      if(ao(j,ilm).ne.0.0) me=me+1                                         2105
15201 continue                                                             2105
15202 continue                                                             2105
      if(me.gt.ne)goto 14872                                               2106
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 14872              2107
      if(dev(ilm).gt.devmax)goto 14872                                     2108
14871 continue                                                             2109
14872 continue                                                             2109
      g=f                                                                  2110
11400 continue                                                             2110
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm)                     2111
      return                                                               2112
      end                                                                  2113
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 2114
      real ca(nin),x(n,*),f(n)                                             2114
      integer ia(nin)                                                      2115
      f=0.0                                                                2115
      if(nin.le.0) return                                                  2116
15210 do 15211 i=1,n                                                       2116
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      2116
15211 continue                                                             2117
15212 continue                                                             2117
      return                                                               2118
      end                                                                  2119
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         2120
      real y(no),d(no),q(no)                                               2120
      integer jp(no),kp(*)                                                 2121
15220 do 15221 j=1,no                                                      2121
      jp(j)=j                                                              2121
15221 continue                                                             2121
15222 continue                                                             2121
      call psort7(y,jp,1,no)                                               2122
      nj=0                                                                 2122
15230 do 15231 j=1,no                                                      2122
      if(q(jp(j)).le.0.0)goto 15231                                        2122
      nj=nj+1                                                              2122
      jp(nj)=jp(j)                                                         2122
15231 continue                                                             2123
15232 continue                                                             2123
      if(nj .ne. 0)goto 15251                                              2123
      jerr=-4                                                              2123
      return                                                               2123
15251 continue                                                             2124
      j=1                                                                  2124
15260 continue                                                             2124
15261 if(d(jp(j)).gt.0.0)goto 15262                                        2124
      j=j+1                                                                2124
      if(j.gt.nj)goto 15262                                                2124
      goto 15261                                                           2125
15262 continue                                                             2125
      if(j .lt. nj-1)goto 15281                                            2125
      jerr=-5                                                              2125
      return                                                               2125
15281 continue                                                             2126
      j0=j-1                                                               2126
      nj=nj-j0                                                             2126
15290 do 15291 j=1,nj                                                      2126
      jp(j)=jp(j+j0)                                                       2126
15291 continue                                                             2127
15292 continue                                                             2127
      jerr=0                                                               2127
      nk=0                                                                 2127
      t0=y(jp(1))                                                          2127
      yk=t0                                                                2127
      j=2                                                                  2128
15300 continue                                                             2128
15301 continue                                                             2128
15310 continue                                                             2129
15311 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 15312                     2129
      j=j+1                                                                2129
      if(j.gt.nj)goto 15312                                                2129
      goto 15311                                                           2130
15312 continue                                                             2130
      nk=nk+1                                                              2130
      kp(nk)=j-1                                                           2130
      if(j.gt.nj)goto 15302                                                2131
      if(j .ne. nj)goto 15331                                              2131
      nk=nk+1                                                              2131
      kp(nk)=nj                                                            2131
      goto 15302                                                           2131
15331 continue                                                             2132
      yk=y(jp(j))                                                          2132
      j=j+1                                                                2133
      goto 15301                                                           2134
15302 continue                                                             2134
      return                                                               2135
      end                                                                  2136
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     2137
      real d(no),dk(nk),wr(no),w(no)                                       2138
      real e(no),u(no),b,c                                                 2138
      integer kp(nk),jp(no)                                                2139
      call usk(no,nk,kp,jp,e,u)                                            2140
      b=dk(1)/u(1)                                                         2140
      c=dk(1)/u(1)**2                                                      2140
      jerr=0                                                               2141
15340 do 15341 j=1,kp(1)                                                   2141
      i=jp(j)                                                              2142
      w(i)=e(i)*(b-e(i)*c)                                                 2142
      if(w(i) .gt. 0.0)goto 15361                                          2142
      jerr=-3                                                              2142
      return                                                               2142
15361 continue                                                             2143
      wr(i)=d(i)-e(i)*b                                                    2144
15341 continue                                                             2145
15342 continue                                                             2145
15370 do 15371 k=2,nk                                                      2145
      j1=kp(k-1)+1                                                         2145
      j2=kp(k)                                                             2146
      b=b+dk(k)/u(k)                                                       2146
      c=c+dk(k)/u(k)**2                                                    2147
15380 do 15381 j=j1,j2                                                     2147
      i=jp(j)                                                              2148
      w(i)=e(i)*(b-e(i)*c)                                                 2148
      if(w(i) .gt. 0.0)goto 15401                                          2148
      jerr=-3                                                              2148
      return                                                               2148
15401 continue                                                             2149
      wr(i)=d(i)-e(i)*b                                                    2150
15381 continue                                                             2151
15382 continue                                                             2151
15371 continue                                                             2152
15372 continue                                                             2152
      return                                                               2153
      end                                                                  2154
      subroutine vars(no,ni,x,w,v)                                         2155
      real x(no,ni),w(no),v(ni)                                            2156
15410 do 15411 j=1,ni                                                      2156
      v(j)=dot_product(w,x(:,j)**2)                                        2156
15411 continue                                                             2157
15412 continue                                                             2157
      return                                                               2158
      end                                                                  2159
      subroutine died(no,nk,d,kp,jp,dk)                                    2160
      real d(no),dk(nk)                                                    2160
      integer kp(nk),jp(no)                                                2161
      dk(1)=sum(d(jp(1:kp(1))))                                            2162
15420 do 15421 k=2,nk                                                      2162
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  2162
15421 continue                                                             2163
15422 continue                                                             2163
      return                                                               2164
      end                                                                  2165
      subroutine usk(no,nk,kp,jp,e,u)                                      2166
      real e(no),u(nk),h                                                   2166
      integer kp(nk),jp(no)                                                2167
      h=0.0                                                                2168
15430 do 15431 k=nk,1,-1                                                   2168
      j2=kp(k)                                                             2169
      j1=1                                                                 2169
      if(k.gt.1) j1=kp(k-1)+1                                              2170
15440 do 15441 j=j2,j1,-1                                                  2170
      h=h+e(jp(j))                                                         2170
15441 continue                                                             2171
15442 continue                                                             2171
      u(k)=h                                                               2172
15431 continue                                                             2173
15432 continue                                                             2173
      return                                                               2174
      end                                                                  2175
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             2176
      real d(no),dk(nk),f(no)                                              2177
      integer kp(nk),jp(no)                                                2177
      real e(no),u(nk),s                                                   2178
      call usk(no,nk,kp,jp,e,u)                                            2178
      u=log(u)                                                             2179
      risk=dot_product(d,f)-dot_product(dk,u)                              2180
      return                                                               2181
      end                                                                  2182
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 2183
      real x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(nlam)          2184
      real, dimension (:), allocatable :: dk,f,xm,dq,q                          
      real, dimension (:), allocatable :: e,uu                                  
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          2190
      allocate(q(1:no),stat=ierr)                                          2190
      jerr=jerr+ierr                                                       2191
      allocate(uu(1:no),stat=ierr)                                         2191
      jerr=jerr+ierr                                                       2192
      allocate(f(1:no),stat=ierr)                                          2192
      jerr=jerr+ierr                                                       2193
      allocate(dk(1:no),stat=ierr)                                         2193
      jerr=jerr+ierr                                                       2194
      allocate(jp(1:no),stat=ierr)                                         2194
      jerr=jerr+ierr                                                       2195
      allocate(kp(1:no),stat=ierr)                                         2195
      jerr=jerr+ierr                                                       2196
      allocate(dq(1:no),stat=ierr)                                         2196
      jerr=jerr+ierr                                                       2197
      allocate(xm(1:ni),stat=ierr)                                         2197
      jerr=jerr+ierr                                                       2198
      if(jerr.ne.0) go to 11400                                            2199
      q=max(0.0,w)                                                         2199
      sw=sum(q)                                                            2200
      if(sw .gt. 0.0)goto 15461                                            2200
      jerr=9999                                                            2200
      go to 11400                                                          2200
15461 continue                                                             2201
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2202
      if(jerr.ne.0) go to 11400                                            2202
      fmax=log(huge(e(1))*0.1)                                             2203
      dq=d*q                                                               2203
      call died(no,nk,dq,kp,jp,dk)                                         2203
      gm=dot_product(q,g)/sw                                               2204
15470 do 15471 j=1,ni                                                      2204
      xm(j)=dot_product(q,x(:,j))/sw                                       2204
15471 continue                                                             2205
15472 continue                                                             2205
15480 do 15481 lam=1,nlam                                                  2206
15490 do 15491 i=1,no                                                      2206
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       2207
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        2208
15491 continue                                                             2209
15492 continue                                                             2209
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          2210
15481 continue                                                             2211
15482 continue                                                             2211
11400 continue                                                             2211
      deallocate(e,uu,dk,f,jp,kp,dq)                                       2212
      return                                                               2213
      end                                                                  2214
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,ne,nx,nlam,flmin,ulam   2216 
     *,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    2217
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2218
      integer jd(*),ia(nx),nin(nlam)                                       2219
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 15511                                    2223
      jerr=10000                                                           2223
      return                                                               2223
15511 continue                                                             2224
      if(minval(y) .ge. 0.0)goto 15531                                     2224
      jerr=8888                                                            2224
      return                                                               2224
15531 continue                                                             2225
      allocate(ww(1:no),stat=jerr)                                         2226
      allocate(ju(1:ni),stat=ierr)                                         2226
      jerr=jerr+ierr                                                       2227
      allocate(vq(1:ni),stat=ierr)                                         2227
      jerr=jerr+ierr                                                       2228
      allocate(xm(1:ni),stat=ierr)                                         2228
      jerr=jerr+ierr                                                       2229
      if(isd .le. 0)goto 15551                                             2229
      allocate(xs(1:ni),stat=ierr)                                         2229
      jerr=jerr+ierr                                                       2229
15551 continue                                                             2230
      if(jerr.ne.0) return                                                 2231
      call chkvars(no,ni,x,ju)                                             2232
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2233
      if(maxval(ju) .gt. 0)goto 15571                                      2233
      jerr=7777                                                            2233
      go to 11400                                                          2233
15571 continue                                                             2234
      vq=max(0.0,vp)                                                       2234
      vq=vq*ni/sum(vq)                                                     2235
      ww=max(0.0,w)                                                        2235
      sw=sum(ww)                                                           2235
      if(sw .gt. 0.0)goto 15591                                            2235
      jerr=9999                                                            2235
      go to 11400                                                          2235
15591 continue                                                             2236
      ww=ww/sw                                                             2237
      call lstandard1(no,ni,x,ww,ju,isd,thr,sthr,xm,xs)                    2238
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,sthr   2240 
     *,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11400                                            2240
      dev0=2.0*sw*dev0                                                     2241
15600 do 15601 k=1,lmu                                                     2241
      nk=nin(k)                                                            2242
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2243
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2244
15601 continue                                                             2245
15602 continue                                                             2245
11400 continue                                                             2245
      deallocate(ww,ju,vq,xm)                                              2245
      if(isd.gt.0) deallocate(xs)                                          2246
      return                                                               2247
      end                                                                  2248
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,ne,nx,nlam,flmin,ulam   2250 
     *,shr,  isd,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2251 
     *9)
      real x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)                    2252
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2253
      integer ju(ni),m(nx),kin(nlam)                                       2254
      real, dimension (:), allocatable :: t,w,wr,v,a,f,as                       
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                          2259
      allocate(as(1:ni),stat=ierr)                                         2259
      jerr=jerr+ierr                                                       2260
      allocate(t(1:no),stat=ierr)                                          2260
      jerr=jerr+ierr                                                       2261
      allocate(mm(1:ni),stat=ierr)                                         2261
      jerr=jerr+ierr                                                       2262
      allocate(wr(1:no),stat=ierr)                                         2262
      jerr=jerr+ierr                                                       2263
      allocate(v(1:ni),stat=ierr)                                          2263
      jerr=jerr+ierr                                                       2264
      allocate(w(1:no),stat=ierr)                                          2264
      jerr=jerr+ierr                                                       2265
      allocate(f(1:no),stat=ierr)                                          2265
      jerr=jerr+ierr                                                       2266
      if(jerr.ne.0) return                                                 2267
      bta=max(parm,1.0e-3)                                                 2267
      omb=1.0-bta                                                          2268
      t=q*y                                                                2268
      yb=sum(t)                                                            2268
      fmax=log(huge(bta)*0.1)                                              2269
      if(nonzero(no,g) .ne. 0)goto 15621                                   2269
      w=q*yb                                                               2269
      az=log(yb)                                                           2269
      f=az                                                                 2269
      goto 15631                                                           2270
15621 continue                                                             2270
      w=q*exp(sign(min(abs(g),fmax),g))                                    2270
      v0=sum(w)                                                            2270
      eaz=yb/v0                                                            2271
      w=eaz*w                                                              2271
      az=log(eaz)                                                          2271
      f=az+g                                                               2272
15631 continue                                                             2273
15611 continue                                                             2273
      a=0.0                                                                2273
      wr=t-w                                                               2273
      v0=yb                                                                2273
      dv0=yb*(log(yb)-1.0)                                                 2273
      dvr=-yb                                                              2274
15640 do 15641 i=1,no                                                      2274
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               2274
15641 continue                                                             2274
15642 continue                                                             2274
      dvr=dvr-dv0                                                          2274
      dev0=dvr                                                             2275
15650 do 15651 j=1,ni                                                      2275
      if(ju(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                         2275
15651 continue                                                             2276
15652 continue                                                             2276
      if(flmin .ge. 1.0)goto 15671                                         2276
      eqs=max(eps,flmin)                                                   2276
      alf=eqs**(1.0/(nlam-1))                                              2276
15671 continue                                                             2277
      m=0                                                                  2277
      mm=0                                                                 2277
      nlp=0                                                                2277
      nin=nlp                                                              2277
      mnl=min(mnlam,nlam)                                                  2278
15680 do 15681 ilm=1,nlam                                                  2279
      if(flmin .lt. 1.0)goto 15701                                         2279
      al=ulam(ilm)                                                         2279
      goto 15691                                                           2280
15701 if(ilm .le. 2)goto 15711                                             2280
      al=al*alf                                                            2280
      goto 15691                                                           2281
15711 if(ilm .ne. 1)goto 15721                                             2281
      al=big                                                               2281
      goto 15731                                                           2282
15721 continue                                                             2282
      al=0.0                                                               2283
15740 do 15741 j=1,ni                                                      2283
      if(ju(j).eq.0)goto 15741                                             2283
      if(vp(j).le.0.0)goto 15741                                           2284
      al=max(al,abs(dot_product(wr,x(:,j)))/vp(j))                         2285
15741 continue                                                             2286
15742 continue                                                             2286
      al=alf*al/bta                                                        2287
15731 continue                                                             2288
15691 continue                                                             2288
      al2=al*omb                                                           2288
      al1=al*bta                                                           2288
      nit=0                                                                2289
15750 continue                                                             2289
15751 continue                                                             2289
      nit=nit+1                                                            2289
      if(nit .le. maxit)goto 15771                                         2289
      jerr=-2                                                              2289
      go to 11400                                                          2289
15771 continue                                                             2290
      az0=az                                                               2290
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2290
      nit1=0                                                               2291
15780 continue                                                             2291
15781 continue                                                             2291
      nit1=nit1+1                                                          2291
      if(nit1 .le. maxit)goto 15801                                        2291
      jerr=-1                                                              2291
      go to 11400                                                          2291
15801 continue                                                             2292
      nlp=nlp+1                                                            2292
      dlx=0.0                                                              2293
15810 do 15811 k=1,ni                                                      2293
      if(ju(k).eq.0)goto 15811                                             2293
      ak=a(k)                                                              2294
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2294
      au=abs(u)-vp(k)*al1                                                  2295
      if(au .gt. 0.0)goto 15831                                            2295
      a(k)=0.0                                                             2295
      goto 15841                                                           2296
15831 continue                                                             2296
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2296
15841 continue                                                             2297
15821 continue                                                             2297
      if(a(k).eq.ak)goto 15811                                             2297
      d=a(k)-ak                                                            2297
      dlx=max(dlx,v(k)*abs(d))                                             2298
      wr=wr-d*w*x(:,k)                                                     2298
      f=f+d*x(:,k)                                                         2299
      if(mm(k) .ne. 0)goto 15861                                           2299
      nin=nin+1                                                            2299
      if(nin.gt.nx)goto 15812                                              2300
      mm(k)=nin                                                            2300
      m(nin)=k                                                             2301
15861 continue                                                             2302
15811 continue                                                             2303
15812 continue                                                             2303
      if(nin.gt.nx)goto 15782                                              2303
      d=sum(wr)/v0                                                         2304
      az=az+d                                                              2304
      dlx=max(dlx,v0*abs(d))                                               2304
      wr=wr-d*w                                                            2304
      f=f+d                                                                2305
      if(dlx.lt.shr)goto 15782                                             2305
      nit2=0                                                               2306
15870 continue                                                             2306
15871 continue                                                             2306
      nit2=nit2+1                                                          2306
      if(nit2 .le. maxit)goto 15891                                        2306
      jerr=-1                                                              2306
      go to 11400                                                          2306
15891 continue                                                             2307
      nlp=nlp+1                                                            2307
      dlx=0.0                                                              2308
15900 do 15901 l=1,nin                                                     2308
      k=m(l)                                                               2308
      ak=a(k)                                                              2309
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2309
      au=abs(u)-vp(k)*al1                                                  2310
      if(au .gt. 0.0)goto 15921                                            2310
      a(k)=0.0                                                             2310
      goto 15931                                                           2311
15921 continue                                                             2311
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2311
15931 continue                                                             2312
15911 continue                                                             2312
      if(a(k).eq.ak)goto 15901                                             2312
      d=a(k)-ak                                                            2312
      dlx=max(dlx,v(k)*abs(d))                                             2313
      wr=wr-d*w*x(:,k)                                                     2313
      f=f+d*x(:,k)                                                         2315
15901 continue                                                             2315
15902 continue                                                             2315
      d=sum(wr)/v0                                                         2315
      az=az+d                                                              2315
      dlx=max(dlx,v0*abs(d))                                               2315
      wr=wr-d*w                                                            2315
      f=f+d                                                                2316
      if(dlx.lt.shr)goto 15872                                             2316
      goto 15871                                                           2317
15872 continue                                                             2317
      goto 15781                                                           2318
15782 continue                                                             2318
      if(nin.gt.nx)goto 15752                                              2319
      w=q*exp(sign(min(abs(f),fmax),f))                                    2319
      v0=sum(w)                                                            2320
      wr=t-w                                                               2320
15940 do 15941 j=1,ni                                                      2320
      if(ju(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                         2320
15941 continue                                                             2321
15942 continue                                                             2321
      if(v0*abs(az-az0) .ge. shr)goto 15961                                2321
      ix=0                                                                 2322
15970 do 15971 j=1,nin                                                     2322
      if(v(m(j))*abs(a(m(j))-as(m(j))).lt.shr)goto 15971                   2322
      ix=1                                                                 2322
      goto 15972                                                           2322
15971 continue                                                             2323
15972 continue                                                             2323
      if(ix.eq.0)goto 15752                                                2324
15961 continue                                                             2325
      goto 15751                                                           2326
15752 continue                                                             2326
      if(nin .le. nx)goto 15991                                            2326
      jerr=-10000-ilm                                                      2326
      goto 15682                                                           2326
15991 continue                                                             2327
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2327
      kin(ilm)=nin                                                         2328
      a0(ilm)=az                                                           2328
      alm(ilm)=al                                                          2328
      lmu=ilm                                                              2329
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               2330
      if(ilm.lt.mnl)goto 15681                                             2330
      if(flmin.ge.1.0)goto 15681                                           2331
      me=0                                                                 2331
16000 do 16001 j=1,nin                                                     2331
      if(ca(j,ilm).ne.0.0) me=me+1                                         2331
16001 continue                                                             2331
16002 continue                                                             2331
      if(me.gt.ne)goto 15682                                               2332
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 15682              2333
      if(dev(ilm).gt.devmax)goto 15682                                     2334
15681 continue                                                             2335
15682 continue                                                             2335
      g=f                                                                  2336
11400 continue                                                             2336
      deallocate(t,w,wr,v,a,f,as,mm)                                       2337
      return                                                               2338
      end                                                                  2339
      function nonzero(n,v)                                                2340
      real v(n)                                                            2341
      nonzero=0                                                            2341
16010 do 16011 i=1,n                                                       2341
      if(v(i) .eq. 0.0)goto 16031                                          2341
      nonzero=1                                                            2341
      return                                                               2341
16031 continue                                                             2341
16011 continue                                                             2342
16012 continue                                                             2342
      return                                                               2343
      end                                                                  2344
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               2345
      real a(nx,lmu),b(ni,lmu)                                             2345
      integer ia(nx),nin(lmu)                                              2346
16040 do 16041 lam=1,lmu                                                   2346
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        2346
16041 continue                                                             2347
16042 continue                                                             2347
      return                                                               2348
      end                                                                  2349
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           2350
      real a(nx,nc,lmu),b(ni,nc,lmu)                                       2350
      integer ia(nx),nin(lmu)                                              2351
16050 do 16051 lam=1,lmu                                                   2351
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             2351
16051 continue                                                             2352
16052 continue                                                             2352
      return                                                               2353
      end                                                                  2354
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               2355
      real x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)       2356
      real, dimension (:), allocatable :: w                                     
      if(minval(y) .ge. 0.0)goto 16071                                     2359
      jerr=8888                                                            2359
      return                                                               2359
16071 continue                                                             2360
      allocate(w(1:no),stat=jerr)                                          2360
      if(jerr.ne.0) return                                                 2361
      w=max(0.0,q)                                                         2361
      sw=sum(w)                                                            2361
      if(sw .gt. 0.0)goto 16091                                            2361
      jerr=9999                                                            2361
      go to 11400                                                          2361
16091 continue                                                             2362
      yb=dot_product(w,y)/sw                                               2362
      fmax=log(huge(y(1))*0.1)                                             2363
16100 do 16101 lam=1,nlam                                                  2363
      s=0.0                                                                2364
16110 do 16111 i=1,no                                                      2364
      if(w(i).le.0.0)goto 16111                                            2365
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          2366
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      2367
16111 continue                                                             2368
16112 continue                                                             2368
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2369
16101 continue                                                             2370
16102 continue                                                             2370
11400 continue                                                             2370
      deallocate(w)                                                        2371
      return                                                               2372
      end                                                                  2373
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,ne,nx,nlam,fl   2375 
     *min,ulam,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        2376
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2377
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           2378
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16131                                    2382
      jerr=10000                                                           2382
      return                                                               2382
16131 continue                                                             2383
      if(minval(y) .ge. 0.0)goto 16151                                     2383
      jerr=8888                                                            2383
      return                                                               2383
16151 continue                                                             2384
      allocate(ww(1:no),stat=jerr)                                         2385
      allocate(ju(1:ni),stat=ierr)                                         2385
      jerr=jerr+ierr                                                       2386
      allocate(vq(1:ni),stat=ierr)                                         2386
      jerr=jerr+ierr                                                       2387
      allocate(xm(1:ni),stat=ierr)                                         2387
      jerr=jerr+ierr                                                       2388
      allocate(xs(1:ni),stat=ierr)                                         2388
      jerr=jerr+ierr                                                       2389
      if(jerr.ne.0) return                                                 2390
      call spchkvars(no,ni,x,ix,ju)                                        2391
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2392
      if(maxval(ju) .gt. 0)goto 16171                                      2392
      jerr=7777                                                            2392
      go to 11400                                                          2392
16171 continue                                                             2393
      vq=max(0.0,vp)                                                       2393
      vq=vq*ni/sum(vq)                                                     2394
      ww=max(0.0,w)                                                        2394
      sw=sum(ww)                                                           2394
      if(sw .gt. 0.0)goto 16191                                            2394
      jerr=9999                                                            2394
      go to 11400                                                          2394
16191 continue                                                             2395
      ww=ww/sw                                                             2396
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,thr,sthr,xm,xs)            2397
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmin,u   2399 
     *lam,sthr,  isd,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11400                                            2399
      dev0=2.0*sw*dev0                                                     2400
16200 do 16201 k=1,lmu                                                     2400
      nk=nin(k)                                                            2401
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2402
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2403
16201 continue                                                             2404
16202 continue                                                             2404
11400 continue                                                             2404
      deallocate(ww,ju,vq,xm,xs)                                           2405
      return                                                               2406
      end                                                                  2407
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,ne,nx,nlam,fl   2409 
     *min,ulam,shr,  isd,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,je
     *rr)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2410 
     *9)
      real x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),xs(ni)          2411
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2412
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2413
      real, dimension (:), allocatable :: qy,t,w,wr,v,a,as,xm                   
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                          2418
      allocate(as(1:ni),stat=ierr)                                         2418
      jerr=jerr+ierr                                                       2419
      allocate(t(1:no),stat=ierr)                                          2419
      jerr=jerr+ierr                                                       2420
      allocate(mm(1:ni),stat=ierr)                                         2420
      jerr=jerr+ierr                                                       2421
      allocate(wr(1:no),stat=ierr)                                         2421
      jerr=jerr+ierr                                                       2422
      allocate(v(1:ni),stat=ierr)                                          2422
      jerr=jerr+ierr                                                       2423
      allocate(xm(1:ni),stat=ierr)                                         2423
      jerr=jerr+ierr                                                       2424
      allocate(w(1:no),stat=ierr)                                          2424
      jerr=jerr+ierr                                                       2425
      allocate(qy(1:no),stat=ierr)                                         2425
      jerr=jerr+ierr                                                       2426
      if(jerr.ne.0) return                                                 2427
      bta=max(parm,1.0e-3)                                                 2427
      omb=1.0-bta                                                          2427
      fmax=log(huge(bta)*0.1)                                              2428
      qy=q*y                                                               2428
      yb=sum(qy)                                                           2429
      if(nonzero(no,g) .ne. 0)goto 16221                                   2429
      w=q*yb                                                               2429
      az=log(yb)                                                           2429
      uu=az                                                                2430
      xm=yb*xb                                                             2430
      v=yb                                                                 2430
      t=0.0                                                                2431
      goto 16231                                                           2432
16221 continue                                                             2432
      w=q*exp(sign(min(abs(g),fmax),g))                                    2432
      ww=sum(w)                                                            2432
      eaz=yb/ww                                                            2433
      w=eaz*w                                                              2433
      az=log(eaz)                                                          2433
      uu=az                                                                2433
      t=g                                                                  2434
16240 do 16241 j=1,ni                                                      2434
      if(ju(j).eq.0)goto 16241                                             2434
      jb=ix(j)                                                             2434
      je=ix(j+1)-1                                                         2435
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2436
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   2438 
     *b(j)**2)/xs(j)**2
16241 continue                                                             2439
16242 continue                                                             2439
16231 continue                                                             2440
16211 continue                                                             2440
      tt=yb*uu                                                             2440
      ww=yb                                                                2440
      wr=qy-q*(yb*(1.0-uu))                                                2440
      a=0.0                                                                2441
      dv0=yb*(log(yb)-1.0)                                                 2441
      dvr=-yb                                                              2442
16250 do 16251 i=1,no                                                      2442
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             2442
16251 continue                                                             2442
16252 continue                                                             2442
      dvr=dvr-dv0                                                          2442
      dev0=dvr                                                             2443
      if(flmin .ge. 1.0)goto 16271                                         2443
      eqs=max(eps,flmin)                                                   2443
      alf=eqs**(1.0/(nlam-1))                                              2443
16271 continue                                                             2444
      m=0                                                                  2444
      mm=0                                                                 2444
      nlp=0                                                                2444
      nin=nlp                                                              2444
      mnl=min(mnlam,nlam)                                                  2445
16280 do 16281 ilm=1,nlam                                                  2446
      if(flmin .lt. 1.0)goto 16301                                         2446
      al=ulam(ilm)                                                         2446
      goto 16291                                                           2447
16301 if(ilm .le. 2)goto 16311                                             2447
      al=al*alf                                                            2447
      goto 16291                                                           2448
16311 if(ilm .ne. 1)goto 16321                                             2448
      al=big                                                               2448
      goto 16331                                                           2449
16321 continue                                                             2449
      al=0.0                                                               2450
16340 do 16341 j=1,ni                                                      2450
      if(ju(j).eq.0)goto 16341                                             2450
      if(vp(j).le.0.0)goto 16341                                           2451
      jb=ix(j)                                                             2451
      je=ix(j+1)-1                                                         2452
      al=max(al,  abs(dot_product(qy(jx(jb:je)),x(jb:je))-xm(j))/(xs(j)*   2454 
     *vp(j)))
16341 continue                                                             2455
16342 continue                                                             2455
      al=alf*al/bta                                                        2456
16331 continue                                                             2457
16291 continue                                                             2457
      al2=al*omb                                                           2457
      al1=al*bta                                                           2457
      nit=0                                                                2458
16350 continue                                                             2458
16351 continue                                                             2458
      nit=nit+1                                                            2458
      if(nit .le. maxit)goto 16371                                         2458
      jerr=-2                                                              2458
      go to 11400                                                          2458
16371 continue                                                             2459
      az0=az                                                               2459
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2459
      nit1=0                                                               2460
16380 continue                                                             2460
16381 continue                                                             2460
      nit1=nit1+1                                                          2460
      if(nit1 .le. maxit)goto 16401                                        2460
      jerr=-1                                                              2460
      go to 11400                                                          2460
16401 continue                                                             2461
      nlp=nlp+1                                                            2461
      dlx=0.0                                                              2462
16410 do 16411 k=1,ni                                                      2462
      if(ju(k).eq.0)goto 16411                                             2462
      jb=ix(k)                                                             2462
      je=ix(k+1)-1                                                         2462
      ak=a(k)                                                              2463
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2465 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2466
      if(au .gt. 0.0)goto 16431                                            2466
      a(k)=0.0                                                             2466
      goto 16441                                                           2467
16431 continue                                                             2467
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2467
16441 continue                                                             2468
16421 continue                                                             2468
      if(a(k).eq.ak)goto 16411                                             2469
      if(mm(k) .ne. 0)goto 16461                                           2469
      nin=nin+1                                                            2469
      if(nin.gt.nx)goto 16412                                              2470
      mm(k)=nin                                                            2470
      m(nin)=k                                                             2471
16461 continue                                                             2472
      d=a(k)-ak                                                            2472
      dlx=max(dlx,v(k)*abs(d))                                             2472
      dv=d/xs(k)                                                           2473
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2474
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2475
      uu=uu-dv*xb(k)                                                       2475
      tt=tt-dv*xm(k)                                                       2476
16411 continue                                                             2477
16412 continue                                                             2477
      if(nin.gt.nx)goto 16382                                              2477
      d=tt/ww-uu                                                           2478
      az=az+d                                                              2478
      dlx=max(dlx,ww*abs(d))                                               2478
      uu=uu+d                                                              2479
      if(dlx.lt.shr)goto 16382                                             2479
      nit2=0                                                               2480
16470 continue                                                             2480
16471 continue                                                             2480
      nit2=nit2+1                                                          2480
      if(nit2 .le. maxit)goto 16491                                        2480
      jerr=-1                                                              2480
      go to 11400                                                          2480
16491 continue                                                             2481
      nlp=nlp+1                                                            2481
      dlx=0.0                                                              2482
16500 do 16501 l=1,nin                                                     2482
      k=m(l)                                                               2482
      jb=ix(k)                                                             2482
      je=ix(k+1)-1                                                         2482
      ak=a(k)                                                              2483
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2485 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2486
      if(au .gt. 0.0)goto 16521                                            2486
      a(k)=0.0                                                             2486
      goto 16531                                                           2487
16521 continue                                                             2487
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2487
16531 continue                                                             2488
16511 continue                                                             2488
      if(a(k).eq.ak)goto 16501                                             2488
      d=a(k)-ak                                                            2488
      dlx=max(dlx,v(k)*abs(d))                                             2489
      dv=d/xs(k)                                                           2489
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2490
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2491
      uu=uu-dv*xb(k)                                                       2491
      tt=tt-dv*xm(k)                                                       2492
16501 continue                                                             2493
16502 continue                                                             2493
      d=tt/ww-uu                                                           2493
      az=az+d                                                              2493
      dlx=max(dlx,ww*abs(d))                                               2493
      uu=uu+d                                                              2494
      if(dlx.lt.shr)goto 16472                                             2494
      goto 16471                                                           2495
16472 continue                                                             2495
      goto 16381                                                           2496
16382 continue                                                             2496
      if(nin.gt.nx)goto 16352                                              2497
      euu=exp(sign(min(abs(uu),fmax),uu))                                  2498
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                2498
      ww=sum(w)                                                            2499
      wr=qy-w*(1.0-uu)                                                     2499
      tt=sum(wr)                                                           2500
16540 do 16541 j=1,ni                                                      2500
      if(ju(j).eq.0)goto 16541                                             2500
      jb=ix(j)                                                             2500
      je=ix(j+1)-1                                                         2501
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2502
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   2504 
     *b(j)**2)/xs(j)**2
16541 continue                                                             2505
16542 continue                                                             2505
      if(ww*abs(az-az0) .ge. shr)goto 16561                                2505
      ixx=0                                                                2506
16570 do 16571 j=1,nin                                                     2506
      if(v(m(j))*abs(a(m(j))-as(m(j))).lt.shr)goto 16571                   2506
      ixx=1                                                                2506
      goto 16572                                                           2506
16571 continue                                                             2507
16572 continue                                                             2507
      if(ixx.eq.0)goto 16352                                               2508
16561 continue                                                             2509
      goto 16351                                                           2510
16352 continue                                                             2510
      if(nin .le. nx)goto 16591                                            2510
      jerr=-10000-ilm                                                      2510
      goto 16282                                                           2510
16591 continue                                                             2511
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2511
      kin(ilm)=nin                                                         2512
      a0(ilm)=az                                                           2512
      alm(ilm)=al                                                          2512
      lmu=ilm                                                              2513
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        2514
      if(ilm.lt.mnl)goto 16281                                             2514
      if(flmin.ge.1.0)goto 16281                                           2515
      me=0                                                                 2515
16600 do 16601 j=1,nin                                                     2515
      if(ca(j,ilm).ne.0.0) me=me+1                                         2515
16601 continue                                                             2515
16602 continue                                                             2515
      if(me.gt.ne)goto 16282                                               2516
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 16282              2517
      if(dev(ilm).gt.devmax)goto 16282                                     2518
16281 continue                                                             2519
16282 continue                                                             2519
      g=t+uu                                                               2520
11400 continue                                                             2520
      deallocate(t,w,wr,v,a,qy,xm,as,mm)                                   2521
      return                                                               2522
      end                                                                  2523
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       2524
      real x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)           2525
      integer ix(*),jx(*)                                                  2526
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 16621                                     2529
      jerr=8888                                                            2529
      return                                                               2529
16621 continue                                                             2530
      allocate(w(1:no),stat=jerr)                                          2531
      allocate(f(1:no),stat=ierr)                                          2531
      jerr=jerr+ierr                                                       2532
      if(jerr.ne.0) return                                                 2533
      w=max(0.0,q)                                                         2533
      sw=sum(w)                                                            2533
      if(sw .gt. 0.0)goto 16641                                            2533
      jerr=9999                                                            2533
      go to 11400                                                          2533
16641 continue                                                             2534
      yb=dot_product(w,y)/sw                                               2534
      fmax=log(huge(y(1))*0.1)                                             2535
16650 do 16651 lam=1,nlam                                                  2535
      f=a0(lam)                                                            2536
16660 do 16661 j=1,ni                                                      2536
      if(a(j,lam).eq.0.0)goto 16661                                        2536
      jb=ix(j)                                                             2536
      je=ix(j+1)-1                                                         2537
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          2538
16661 continue                                                             2539
16662 continue                                                             2539
      f=f+g                                                                2540
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2541
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2542
16651 continue                                                             2543
16652 continue                                                             2543
11400 continue                                                             2543
      deallocate(w,f)                                                      2544
      return                                                               2545
      end                                                                  2546
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   2547 
     *jerr)
      real x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(nlam)          2548
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 2549
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 16681                                     2552
      jerr=8888                                                            2552
      return                                                               2552
16681 continue                                                             2553
      allocate(w(1:no),stat=jerr)                                          2554
      allocate(f(1:no),stat=ierr)                                          2554
      jerr=jerr+ierr                                                       2555
      if(jerr.ne.0) return                                                 2556
      w=max(0.0,q)                                                         2556
      sw=sum(w)                                                            2556
      if(sw .gt. 0.0)goto 16701                                            2556
      jerr=9999                                                            2556
      go to 11400                                                          2556
16701 continue                                                             2557
      yb=dot_product(w,y)/sw                                               2557
      fmax=log(huge(y(1))*0.1)                                             2558
16710 do 16711 lam=1,nlam                                                  2558
      f=a0(lam)                                                            2559
16720 do 16721 k=1,nin(lam)                                                2559
      j=ia(k)                                                              2559
      jb=ix(j)                                                             2559
      je=ix(j+1)-1                                                         2560
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         2561
16721 continue                                                             2562
16722 continue                                                             2562
      f=f+g                                                                2563
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2564
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2565
16711 continue                                                             2566
16712 continue                                                             2566
11400 continue                                                             2566
      deallocate(w,f)                                                      2567
      return                                                               2568
      end                                                                  2569
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
