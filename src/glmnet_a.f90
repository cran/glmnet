c
c                          newGLMnet (1/26/11)
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
      subroutine elnet  (ka,parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam    606 
     *,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam)                          607
      real ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)                          608
      integer jd(*),ia(nx),nin(nlam)                                        609
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10021                                     612
      jerr=10000                                                            612
      return                                                                612
10021 continue                                                              613
      allocate(vq(1:ni),stat=jerr)                                          613
      if(jerr.ne.0) return                                                  614
      vq=max(0.0,vp)                                                        614
      vq=vq*ni/sum(vq)                                                      615
      if(ka .ne. 1)goto 10041                                               616
      call elnetu  (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd    619 
     *,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                            620
10041 continue                                                              621
      call elnetn (parm,no,ni,x,y,w,jd,vq,ne,nx,nlam,flmin,ulam,thr,isd,    624 
     *  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                              625
10031 continue                                                              625
      deallocate(vq)                                                        626
      return                                                                627
      end                                                                   628
      subroutine elnetu  (parm,no,ni,x,y,w,jd,vp,ne,nx,nlam,flmin,ulam,t    631 
     *hr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ulam(nlam)                           632
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         633
      integer jd(*),ia(nx),nin(nlam)                                        634
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           639
      allocate(xm(1:ni),stat=ierr)                                          639
      jerr=jerr+ierr                                                        640
      allocate(xs(1:ni),stat=ierr)                                          640
      jerr=jerr+ierr                                                        641
      allocate(ju(1:ni),stat=ierr)                                          641
      jerr=jerr+ierr                                                        642
      allocate(xv(1:ni),stat=ierr)                                          642
      jerr=jerr+ierr                                                        643
      allocate(vlam(1:nlam),stat=ierr)                                      643
      jerr=jerr+ierr                                                        644
      if(jerr.ne.0) return                                                  645
      call chkvars(no,ni,x,ju)                                              646
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  647
      if(maxval(ju) .gt. 0)goto 10071                                       647
      jerr=7777                                                             647
      return                                                                647
10071 continue                                                              648
      call standard(no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)               649
      if(jerr.ne.0) return                                                  650
      if(flmin.ge.1.0) vlam=ulam/ys                                         651
      call elnet1(parm,ni,ju,vp,g,no,ne,nx,x,nlam,flmin,vlam,thr,xv,  lm    653 
     *u,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                  654
10080 do 10081 k=1,lmu                                                      654
      alm(k)=ys*alm(k)                                                      654
      nk=nin(k)                                                             655
10090 do 10091 l=1,nk                                                       655
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          655
10091 continue                                                              656
10092 continue                                                              656
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         657
10081 continue                                                              658
10082 continue                                                              658
      deallocate(xm,xs,g,ju,xv,vlam)                                        659
      return                                                                660
      end                                                                   661
      subroutine standard (no,ni,x,y,w,isd,ju,g,xm,xs,ym,ys,xv,jerr)        662
      real x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                  662
      integer ju(ni)                                                        663
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           666
      if(jerr.ne.0) return                                                  667
      w=w/sum(w)                                                            667
      v=sqrt(w)                                                             668
10100 do 10101 j=1,ni                                                       668
      if(ju(j).eq.0)goto 10101                                              669
      xm(j)=dot_product(w,x(:,j))                                           669
      x(:,j)=v*(x(:,j)-xm(j))                                               670
      xv(j)=dot_product(x(:,j),x(:,j))                                      670
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        671
10101 continue                                                              672
10102 continue                                                              672
      if(isd .ne. 0)goto 10121                                              672
      xs=1.0                                                                672
      goto 10131                                                            673
10121 continue                                                              674
10140 do 10141 j=1,ni                                                       674
      if(ju(j).eq.0)goto 10141                                              674
      x(:,j)=x(:,j)/xs(j)                                                   674
10141 continue                                                              675
10142 continue                                                              675
      xv=1.0                                                                676
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
      dlx=max(xv(k)*del**2,dlx)                                             726
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
      dlx=max(xv(k)*del**2,dlx)                                             737
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
      call standard1(no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)                775
      if(jerr.ne.0) return                                                  776
      if(flmin.ge.1.0) vlam=ulam/ys                                         777
      call elnet2(parm,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,vlam,thr,xv,  lm    779 
     *u,ca,ia,nin,rsq,alm,nlp,jerr)
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
      subroutine standard1 (no,ni,x,y,w,isd,ju,xm,xs,ym,ys,xv,jerr)         788
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)                        788
      integer ju(ni)                                                        789
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           792
      if(jerr.ne.0) return                                                  793
      w=w/sum(w)                                                            793
      v=sqrt(w)                                                             794
10450 do 10451 j=1,ni                                                       794
      if(ju(j).eq.0)goto 10451                                              795
      xm(j)=dot_product(w,x(:,j))                                           795
      x(:,j)=v*(x(:,j)-xm(j))                                               796
      xv(j)=dot_product(x(:,j),x(:,j))                                      796
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        797
10451 continue                                                              798
10452 continue                                                              798
      if(isd .ne. 0)goto 10471                                              798
      xs=1.0                                                                798
      goto 10481                                                            799
10471 continue                                                              799
10490 do 10491 j=1,ni                                                       799
      if(ju(j).eq.0)goto 10491                                              799
      x(:,j)=x(:,j)/xs(j)                                                   799
10491 continue                                                              800
10492 continue                                                              800
      xv=1.0                                                                801
10481 continue                                                              802
10461 continue                                                              802
      ym=dot_product(w,y)                                                   802
      y=v*(y-ym)                                                            802
      ys=sqrt(dot_product(y,y))                                             802
      y=y/ys                                                                803
      deallocate(v)                                                         804
      return                                                                805
      end                                                                   806
      subroutine elnet2(beta,ni,ju,vp,y,no,ne,nx,x,nlam,flmin,ulam,thr,x    808 
     *v,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    809 
     *9)
      real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    810 
     *nlam),xv(ni)
      integer ju(ni),ia(nx),kin(nlam)                                       811
      real, dimension (:), allocatable :: a                                     
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                           816
      allocate(mm(1:ni),stat=ierr)                                          816
      jerr=jerr+ierr                                                        817
      if(jerr.ne.0) return                                                  818
      bta=max(beta,1.0e-3)                                                  818
      omb=1.0-bta                                                           819
      if(flmin .ge. 1.0)goto 10511                                          819
      eqs=max(eps,flmin)                                                    819
      alf=eqs**(1.0/(nlam-1))                                               819
10511 continue                                                              820
      rsq=0.0                                                               820
      a=0.0                                                                 820
      mm=0                                                                  820
      nlp=0                                                                 820
      nin=nlp                                                               820
      iz=0                                                                  820
      mnl=min(mnlam,nlam)                                                   821
10520 do 10521 m=1,nlam                                                     822
      if(flmin .lt. 1.0)goto 10541                                          822
      alm=ulam(m)                                                           822
      goto 10531                                                            823
10541 if(m .le. 2)goto 10551                                                823
      alm=alm*alf                                                           823
      goto 10531                                                            824
10551 if(m .ne. 1)goto 10561                                                824
      alm=big                                                               824
      goto 10571                                                            825
10561 continue                                                              825
      alm=0.0                                                               826
10580 do 10581 j=1,ni                                                       826
      if(ju(j).eq.0)goto 10581                                              826
      if(vp(j).le.0.0)goto 10581                                            827
      alm=max(alm,abs(dot_product(y,x(:,j)))/vp(j))                         828
10581 continue                                                              829
10582 continue                                                              829
      alm=alf*alm/bta                                                       830
10571 continue                                                              831
10531 continue                                                              831
      dem=alm*omb                                                           831
      ab=alm*bta                                                            831
      rsq0=rsq                                                              831
      jz=1                                                                  832
10590 continue                                                              832
10591 continue                                                              832
      if(iz*jz.ne.0) go to 10260                                            832
      nlp=nlp+1                                                             832
      dlx=0.0                                                               833
10600 do 10601 k=1,ni                                                       833
      if(ju(k).eq.0)goto 10601                                              833
      gk=dot_product(y,x(:,k))                                              834
      ak=a(k)                                                               834
      u=gk+ak*xv(k)                                                         834
      v=abs(u)-vp(k)*ab                                                     834
      a(k)=0.0                                                              835
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         836
      if(a(k).eq.ak)goto 10601                                              837
      if(mm(k) .ne. 0)goto 10621                                            837
      nin=nin+1                                                             837
      if(nin.gt.nx)goto 10602                                               838
      mm(k)=nin                                                             838
      ia(nin)=k                                                             839
10621 continue                                                              840
      del=a(k)-ak                                                           840
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        841
      y=y-del*x(:,k)                                                        841
      dlx=max(xv(k)*del**2,dlx)                                             842
10601 continue                                                              843
10602 continue                                                              843
      if(dlx.lt.thr)goto 10592                                              843
      if(nin.gt.nx)goto 10592                                               844
10260 continue                                                              844
      iz=1                                                                  845
10630 continue                                                              845
10631 continue                                                              845
      nlp=nlp+1                                                             845
      dlx=0.0                                                               846
10640 do 10641 l=1,nin                                                      846
      k=ia(l)                                                               846
      gk=dot_product(y,x(:,k))                                              847
      ak=a(k)                                                               847
      u=gk+ak*xv(k)                                                         847
      v=abs(u)-vp(k)*ab                                                     847
      a(k)=0.0                                                              848
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         849
      if(a(k).eq.ak)goto 10641                                              850
      del=a(k)-ak                                                           850
      rsq=rsq+del*(2.0*gk-del*xv(k))                                        851
      y=y-del*x(:,k)                                                        851
      dlx=max(xv(k)*del**2,dlx)                                             852
10641 continue                                                              853
10642 continue                                                              853
      if(dlx.lt.thr)goto 10632                                              853
      goto 10631                                                            854
10632 continue                                                              854
      jz=0                                                                  855
      goto 10591                                                            856
10592 continue                                                              856
      if(nin.gt.nx)goto 10522                                               857
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 857
      kin(m)=nin                                                            858
      rsqo(m)=rsq                                                           858
      almo(m)=alm                                                           858
      lmu=m                                                                 859
      if(m.lt.mnl)goto 10521                                                859
      if(flmin.ge.1.0)goto 10521                                            860
      me=0                                                                  860
10650 do 10651 j=1,nin                                                      860
      if(ao(j,m).ne.0.0) me=me+1                                            860
10651 continue                                                              860
10652 continue                                                              860
      if(me.gt.ne)goto 10522                                                861
      if(rsq-rsq0.lt.sml*rsq)goto 10522                                     861
      if(rsq.gt.rsqmax)goto 10522                                           862
10521 continue                                                              863
10522 continue                                                              863
      deallocate(a,mm)                                                      864
      return                                                                865
      end                                                                   866
      subroutine chkvars(no,ni,x,ju)                                        867
      real x(no,ni)                                                         867
      integer ju(ni)                                                        868
10660 do 10661 j=1,ni                                                       868
      ju(j)=0                                                               868
      t=x(1,j)                                                              869
10670 do 10671 i=2,no                                                       869
      if(x(i,j).eq.t)goto 10671                                             869
      ju(j)=1                                                               869
      goto 10672                                                            869
10671 continue                                                              870
10672 continue                                                              870
10661 continue                                                              871
10662 continue                                                              871
      return                                                                872
      end                                                                   873
      subroutine uncomp(ni,ca,ia,nin,a)                                     874
      real ca(*),a(ni)                                                      874
      integer ia(*)                                                         875
      a=0.0                                                                 875
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                   876
      return                                                                877
      end                                                                   878
      subroutine modval(a0,ca,ia,nin,n,x,f)                                 879
      real ca(nin),x(n,*),f(n)                                              879
      integer ia(nin)                                                       880
      f=a0                                                                  880
      if(nin.le.0) return                                                   881
10680 do 10681 i=1,n                                                        881
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                       881
10681 continue                                                              882
10682 continue                                                              882
      return                                                                883
      end                                                                   884
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,fl    887 
     *min,ulam,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                               888
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         889
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                            890
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10701                                     893
      jerr=10000                                                            893
      return                                                                893
10701 continue                                                              894
      allocate(vq(1:ni),stat=jerr)                                          894
      if(jerr.ne.0) return                                                  895
      vq=max(0.0,vp)                                                        895
      vq=vq*ni/sum(vq)                                                      896
      if(ka .ne. 1)goto 10721                                               897
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam    900 
     *,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10731                                                            901
10721 continue                                                              902
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,ne,nx,nlam,flmin,ulam,    905 
     *thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10731 continue                                                              906
10711 continue                                                              906
      deallocate(vq)                                                        907
      return                                                                908
      end                                                                   909
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmi    912 
     *n,ulam,thr,isd,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),y(no),w(no),vp(ni),ulam(nlam)                               913
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         914
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                            915
      real, dimension (:), allocatable :: xm,xs,g,xv,vlam                       
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           920
      allocate(xm(1:ni),stat=ierr)                                          920
      jerr=jerr+ierr                                                        921
      allocate(xs(1:ni),stat=ierr)                                          921
      jerr=jerr+ierr                                                        922
      allocate(ju(1:ni),stat=ierr)                                          922
      jerr=jerr+ierr                                                        923
      allocate(xv(1:ni),stat=ierr)                                          923
      jerr=jerr+ierr                                                        924
      allocate(vlam(1:nlam),stat=ierr)                                      924
      jerr=jerr+ierr                                                        925
      if(jerr.ne.0) return                                                  926
      call spchkvars(no,ni,x,ix,ju)                                         927
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  928
      if(maxval(ju) .gt. 0)goto 10751                                       928
      jerr=7777                                                             928
      return                                                                928
10751 continue                                                              929
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,jerr)       930
      if(jerr.ne.0) return                                                  931
      if(flmin.ge.1.0) vlam=ulam/ys                                         932
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t    934 
     *hr,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                  935
10760 do 10761 k=1,lmu                                                      935
      alm(k)=ys*alm(k)                                                      935
      nk=nin(k)                                                             936
10770 do 10771 l=1,nk                                                       936
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          936
10771 continue                                                              937
10772 continue                                                              937
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                         938
10761 continue                                                              939
10762 continue                                                              939
      deallocate(xm,xs,g,ju,xv,vlam)                                        940
      return                                                                941
      end                                                                   942
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,g,xm,xs,ym,ys,xv,j    943 
     *err)
      real x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)                      943
      integer ix(*),jx(*),ju(ni)                                            944
      w=w/sum(w)                                                            945
10780 do 10781 j=1,ni                                                       945
      if(ju(j).eq.0)goto 10781                                              946
      jb=ix(j)                                                              946
      je=ix(j+1)-1                                                          946
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                              947
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                  948
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        949
10781 continue                                                              950
10782 continue                                                              950
      if(isd .ne. 0)goto 10801                                              950
      xs=1.0                                                                950
      goto 10811                                                            950
10801 continue                                                              950
      xv=1.0                                                                950
10811 continue                                                              951
10791 continue                                                              951
      ym=dot_product(w,y)                                                   951
      y=y-ym                                                                951
      ys=sqrt(dot_product(w,y**2))                                          951
      y=y/ys                                                                951
      g=0.0                                                                 952
10820 do 10821 j=1,ni                                                       952
      if(ju(j).eq.0)goto 10821                                              952
      jb=ix(j)                                                              952
      je=ix(j+1)-1                                                          953
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)            954
10821 continue                                                              955
10822 continue                                                              955
      return                                                                956
      end                                                                   957
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,nlam,flmin,    959 
     *ulam,  thr,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99    960 
     *9)
      real g(ni),vp(ni),x(*),ulam(nlam),w(no)                               961
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)           962
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                           963
      real, dimension (:), allocatable :: a,da                                  
      integer, dimension (:), allocatable :: mm                                 
      real, dimension (:,:), allocatable :: c                                   
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      allocate(a(1:ni),stat=ierr)                                           969
      jerr=jerr+ierr                                                        970
      allocate(mm(1:ni),stat=ierr)                                          970
      jerr=jerr+ierr                                                        971
      allocate(da(1:ni),stat=ierr)                                          971
      jerr=jerr+ierr                                                        972
      if(jerr.ne.0) return                                                  973
      bta=max(beta,1.0e-3)                                                  973
      omb=1.0-bta                                                           974
      if(flmin .ge. 1.0)goto 10841                                          974
      eqs=max(eps,flmin)                                                    974
      alf=eqs**(1.0/(nlam-1))                                               974
10841 continue                                                              975
      rsq=0.0                                                               975
      a=0.0                                                                 975
      mm=0                                                                  975
      nlp=0                                                                 975
      nin=nlp                                                               975
      iz=0                                                                  975
      mnl=min(mnlam,nlam)                                                   976
10850 do 10851 m=1,nlam                                                     977
      if(flmin .lt. 1.0)goto 10871                                          977
      alm=ulam(m)                                                           977
      goto 10861                                                            978
10871 if(m .le. 2)goto 10881                                                978
      alm=alm*alf                                                           978
      goto 10861                                                            979
10881 if(m .ne. 1)goto 10891                                                979
      alm=big                                                               979
      goto 10901                                                            980
10891 continue                                                              980
      alm=0.0                                                               981
10910 do 10911 j=1,ni                                                       981
      if(ju(j).eq.0)goto 10911                                              981
      if(vp(j).le.0.0)goto 10911                                            982
      alm=max(alm,abs(g(j))/vp(j))                                          983
10911 continue                                                              984
10912 continue                                                              984
      alm=alf*alm/bta                                                       985
10901 continue                                                              986
10861 continue                                                              986
      dem=alm*omb                                                           986
      ab=alm*bta                                                            986
      rsq0=rsq                                                              986
      jz=1                                                                  987
10920 continue                                                              987
10921 continue                                                              987
      if(iz*jz.ne.0) go to 10260                                            987
      nlp=nlp+1                                                             987
      dlx=0.0                                                               988
10930 do 10931 k=1,ni                                                       988
      if(ju(k).eq.0)goto 10931                                              989
      ak=a(k)                                                               989
      u=g(k)+ak*xv(k)                                                       989
      v=abs(u)-vp(k)*ab                                                     989
      a(k)=0.0                                                              990
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                         991
      if(a(k).eq.ak)goto 10931                                              992
      if(mm(k) .ne. 0)goto 10951                                            992
      nin=nin+1                                                             992
      if(nin.gt.nx)goto 10932                                               993
10960 do 10961 j=1,ni                                                       993
      if(ju(j).eq.0)goto 10961                                              994
      if(mm(j) .eq. 0)goto 10981                                            994
      c(j,nin)=c(k,mm(j))                                                   994
      goto 10961                                                            994
10981 continue                                                              995
      if(j .ne. k)goto 11001                                                995
      c(j,nin)=xv(j)                                                        995
      goto 10961                                                            995
11001 continue                                                              996
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))        998
10961 continue                                                              999
10962 continue                                                              999
      mm(k)=nin                                                             999
      ia(nin)=k                                                            1000
10951 continue                                                             1001
      del=a(k)-ak                                                          1001
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1002
      dlx=max(xv(k)*del**2,dlx)                                            1003
11010 do 11011 j=1,ni                                                      1003
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                              1003
11011 continue                                                             1004
11012 continue                                                             1004
10931 continue                                                             1005
10932 continue                                                             1005
      if(dlx.lt.thr)goto 10922                                             1005
      if(nin.gt.nx)goto 10922                                              1006
10260 continue                                                             1006
      iz=1                                                                 1006
      da(1:nin)=a(ia(1:nin))                                               1007
11020 continue                                                             1007
11021 continue                                                             1007
      nlp=nlp+1                                                            1007
      dlx=0.0                                                              1008
11030 do 11031 l=1,nin                                                     1008
      k=ia(l)                                                              1009
      ak=a(k)                                                              1009
      u=g(k)+ak*xv(k)                                                      1009
      v=abs(u)-vp(k)*ab                                                    1009
      a(k)=0.0                                                             1010
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1011
      if(a(k).eq.ak)goto 11031                                             1012
      del=a(k)-ak                                                          1012
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1013
      dlx=max(xv(k)*del**2,dlx)                                            1014
11040 do 11041 j=1,nin                                                     1014
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                 1014
11041 continue                                                             1015
11042 continue                                                             1015
11031 continue                                                             1016
11032 continue                                                             1016
      if(dlx.lt.thr)goto 11022                                             1016
      goto 11021                                                           1017
11022 continue                                                             1017
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                     1018
11050 do 11051 j=1,ni                                                      1018
      if(mm(j).ne.0)goto 11051                                             1019
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))           1020
11051 continue                                                             1021
11052 continue                                                             1021
      jz=0                                                                 1022
      goto 10921                                                           1023
10922 continue                                                             1023
      if(nin.gt.nx)goto 10852                                              1024
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1024
      kin(m)=nin                                                           1025
      rsqo(m)=rsq                                                          1025
      almo(m)=alm                                                          1025
      lmu=m                                                                1026
      if(m.lt.mnl)goto 10851                                               1026
      if(flmin.ge.1.0)goto 10851                                           1027
      me=0                                                                 1027
11060 do 11061 j=1,nin                                                     1027
      if(ao(j,m).ne.0.0) me=me+1                                           1027
11061 continue                                                             1027
11062 continue                                                             1027
      if(me.gt.ne)goto 10852                                               1028
      if(rsq-rsq0.lt.sml*rsq)goto 10852                                    1028
      if(rsq.gt.rsqmax)goto 10852                                          1029
10851 continue                                                             1030
10852 continue                                                             1030
      deallocate(a,mm,c,da)                                                1031
      return                                                               1032
      end                                                                  1033
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,ne,nx,nlam,flmin,   1035 
     *ulam,  thr,isd,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(*),vp(ni),y(no),w(no),ulam(nlam)                              1036
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                        1037
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1038
      real, dimension (:), allocatable :: xm,xs,xv,vlam                         
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         1043
      allocate(xs(1:ni),stat=ierr)                                         1043
      jerr=jerr+ierr                                                       1044
      allocate(ju(1:ni),stat=ierr)                                         1044
      jerr=jerr+ierr                                                       1045
      allocate(xv(1:ni),stat=ierr)                                         1045
      jerr=jerr+ierr                                                       1046
      allocate(vlam(1:nlam),stat=ierr)                                     1046
      jerr=jerr+ierr                                                       1047
      if(jerr.ne.0) return                                                 1048
      call spchkvars(no,ni,x,ix,ju)                                        1049
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1050
      if(maxval(ju) .gt. 0)goto 11081                                      1050
      jerr=7777                                                            1050
      return                                                               1050
11081 continue                                                             1051
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,jerr)       1052
      if(jerr.ne.0) return                                                 1053
      if(flmin.ge.1.0) vlam=ulam/ys                                        1054
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,vlam,t   1056 
     *hr,xm,xs,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.ne.0) return                                                 1057
11090 do 11091 k=1,lmu                                                     1057
      alm(k)=ys*alm(k)                                                     1057
      nk=nin(k)                                                            1058
11100 do 11101 l=1,nk                                                      1058
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1058
11101 continue                                                             1059
11102 continue                                                             1059
      a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))                        1060
11091 continue                                                             1061
11092 continue                                                             1061
      deallocate(xm,xs,ju,xv,vlam)                                         1062
      return                                                               1063
      end                                                                  1064
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,xm,xs,ym,ys,xv,je   1065 
     *rr)
      real x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)                           1065
      integer ix(*),jx(*),ju(ni)                                           1066
      w=w/sum(w)                                                           1067
11110 do 11111 j=1,ni                                                      1067
      if(ju(j).eq.0)goto 11111                                             1068
      jb=ix(j)                                                             1068
      je=ix(j+1)-1                                                         1068
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1069
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1070
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1071
11111 continue                                                             1072
11112 continue                                                             1072
      if(isd .ne. 0)goto 11131                                             1072
      xs=1.0                                                               1072
      goto 11141                                                           1072
11131 continue                                                             1072
      xv=1.0                                                               1072
11141 continue                                                             1073
11121 continue                                                             1073
      ym=dot_product(w,y)                                                  1073
      y=y-ym                                                               1073
      ys=sqrt(dot_product(w,y**2))                                         1073
      y=y/ys                                                               1074
      return                                                               1075
      end                                                                  1076
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,nlam,flmin,   1078 
     *ulam,  thr,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, big=9.9e30, mnlam=5, rsqmax=0.99   1079 
     *9)
      real y(no),w(no),x(*),vp(ni),ulam(nlam)                              1080
      real ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),xv(ni)          1081
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1082
      real, dimension (:), allocatable :: a                                     
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                          1087
      allocate(mm(1:ni),stat=ierr)                                         1087
      jerr=jerr+ierr                                                       1088
      if(jerr.ne.0) return                                                 1089
      bta=max(beta,1.0e-3)                                                 1089
      omb=1.0-bta                                                          1090
      if(flmin .ge. 1.0)goto 11161                                         1090
      eqs=max(eps,flmin)                                                   1090
      alf=eqs**(1.0/(nlam-1))                                              1090
11161 continue                                                             1091
      rsq=0.0                                                              1091
      a=0.0                                                                1091
      mm=0                                                                 1091
      o=0.0                                                                1091
      nlp=0                                                                1091
      nin=nlp                                                              1091
      iz=0                                                                 1091
      mnl=min(mnlam,nlam)                                                  1092
11170 do 11171 m=1,nlam                                                    1093
      if(flmin .lt. 1.0)goto 11191                                         1093
      alm=ulam(m)                                                          1093
      goto 11181                                                           1094
11191 if(m .le. 2)goto 11201                                               1094
      alm=alm*alf                                                          1094
      goto 11181                                                           1095
11201 if(m .ne. 1)goto 11211                                               1095
      alm=big                                                              1095
      goto 11221                                                           1096
11211 continue                                                             1096
      alm=0.0                                                              1097
11230 do 11231 j=1,ni                                                      1097
      if(ju(j).eq.0)goto 11231                                             1097
      if(vp(j).le.0.0)goto 11231                                           1098
      jb=ix(j)                                                             1098
      je=ix(j+1)-1                                                         1099
      alm=max(alm,abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))    1101 
     * /(vp(j)*xs(j))))
11231 continue                                                             1102
11232 continue                                                             1102
      alm=alf*alm/bta                                                      1103
11221 continue                                                             1104
11181 continue                                                             1104
      dem=alm*omb                                                          1104
      ab=alm*bta                                                           1104
      rsq0=rsq                                                             1104
      jz=1                                                                 1105
11240 continue                                                             1105
11241 continue                                                             1105
      if(iz*jz.ne.0) go to 10260                                           1105
      nlp=nlp+1                                                            1105
      dlx=0.0                                                              1106
11250 do 11251 k=1,ni                                                      1106
      if(ju(k).eq.0)goto 11251                                             1106
      jb=ix(k)                                                             1106
      je=ix(k+1)-1                                                         1107
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1108
      ak=a(k)                                                              1108
      u=gk+ak*xv(k)                                                        1108
      v=abs(u)-vp(k)*ab                                                    1108
      a(k)=0.0                                                             1109
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1110
      if(a(k).eq.ak)goto 11251                                             1111
      if(mm(k) .ne. 0)goto 11271                                           1111
      nin=nin+1                                                            1111
      if(nin.gt.nx)goto 11252                                              1112
      mm(k)=nin                                                            1112
      ia(nin)=k                                                            1113
11271 continue                                                             1114
      del=a(k)-ak                                                          1114
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1115
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1116
      o=o+del*xm(k)/xs(k)                                                  1116
      dlx=max(xv(k)*del**2,dlx)                                            1117
11251 continue                                                             1118
11252 continue                                                             1118
      if(dlx.lt.thr)goto 11242                                             1118
      if(nin.gt.nx)goto 11242                                              1119
10260 continue                                                             1119
      iz=1                                                                 1120
11280 continue                                                             1120
11281 continue                                                             1120
      nlp=nlp+1                                                            1120
      dlx=0.0                                                              1121
11290 do 11291 l=1,nin                                                     1121
      k=ia(l)                                                              1121
      jb=ix(k)                                                             1121
      je=ix(k+1)-1                                                         1122
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1123
      ak=a(k)                                                              1123
      u=gk+ak*xv(k)                                                        1123
      v=abs(u)-vp(k)*ab                                                    1123
      a(k)=0.0                                                             1124
      if(v.gt.0.0) a(k)=sign(v,u)/(xv(k)+vp(k)*dem)                        1125
      if(a(k).eq.ak)goto 11291                                             1126
      del=a(k)-ak                                                          1126
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1127
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1128
      o=o+del*xm(k)/xs(k)                                                  1128
      dlx=max(xv(k)*del**2,dlx)                                            1129
11291 continue                                                             1130
11292 continue                                                             1130
      if(dlx.lt.thr)goto 11282                                             1130
      goto 11281                                                           1131
11282 continue                                                             1131
      jz=0                                                                 1132
      goto 11241                                                           1133
11242 continue                                                             1133
      if(nin.gt.nx)goto 11172                                              1134
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1134
      kin(m)=nin                                                           1135
      rsqo(m)=rsq                                                          1135
      almo(m)=alm                                                          1135
      lmu=m                                                                1136
      if(m.lt.mnl)goto 11171                                               1136
      if(flmin.ge.1.0)goto 11171                                           1137
      me=0                                                                 1137
11300 do 11301 j=1,nin                                                     1137
      if(ao(j,m).ne.0.0) me=me+1                                           1137
11301 continue                                                             1137
11302 continue                                                             1137
      if(me.gt.ne)goto 11172                                               1138
      if(rsq-rsq0.lt.sml*rsq)goto 11172                                    1138
      if(rsq.gt.rsqmax)goto 11172                                          1139
11171 continue                                                             1140
11172 continue                                                             1140
      deallocate(a,mm)                                                     1141
      return                                                               1142
      end                                                                  1143
      subroutine spchkvars(no,ni,x,ix,ju)                                  1144
      real x(*)                                                            1144
      integer ix(*),ju(ni)                                                 1145
11310 do 11311 j=1,ni                                                      1145
      ju(j)=0                                                              1145
      jb=ix(j)                                                             1145
      nj=ix(j+1)-jb                                                        1145
      if(nj.eq.0)goto 11311                                                1146
      je=ix(j+1)-1                                                         1147
      if(nj .ge. no)goto 11331                                             1147
11340 do 11341 i=jb,je                                                     1147
      if(x(i).eq.0.0)goto 11341                                            1147
      ju(j)=1                                                              1147
      goto 11342                                                           1147
11341 continue                                                             1147
11342 continue                                                             1147
      goto 11351                                                           1148
11331 continue                                                             1148
      t=x(jb)                                                              1148
11360 do 11361 i=jb+1,je                                                   1148
      if(x(i).eq.t)goto 11361                                              1148
      ju(j)=1                                                              1148
      goto 11362                                                           1148
11361 continue                                                             1148
11362 continue                                                             1148
11351 continue                                                             1149
11321 continue                                                             1149
11311 continue                                                             1150
11312 continue                                                             1150
      return                                                               1151
      end                                                                  1152
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                         1153
      real ca(*),x(*),f(n)                                                 1153
      integer ia(*),ix(*),jx(*)                                            1154
      f=a0                                                                 1155
11370 do 11371 j=1,nin                                                     1155
      k=ia(j)                                                              1155
      kb=ix(k)                                                             1155
      ke=ix(k+1)-1                                                         1156
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                             1157
11371 continue                                                             1158
11372 continue                                                             1158
      return                                                               1159
      end                                                                  1160
      function row_prod(i,j,ia,ja,ra,w)                                    1161
      integer ia(*),ja(*)                                                  1161
      real ra(*),w(*)                                                      1162
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(   1164 
     *i),ia(j+1)-ia(j),w)
      return                                                               1165
      end                                                                  1166
      function dot(x,y,mx,my,nx,ny,w)                                      1167
      real x(*),y(*),w(*)                                                  1167
      integer mx(*),my(*)                                                  1168
      i=1                                                                  1168
      j=i                                                                  1168
      s=0.0                                                                1169
11380 continue                                                             1169
11381 continue                                                             1169
11390 continue                                                             1170
11391 if(mx(i).ge.my(j))goto 11392                                         1170
      i=i+1                                                                1170
      if(i.gt.nx) go to 11400                                              1170
      goto 11391                                                           1171
11392 continue                                                             1171
      if(mx(i).eq.my(j)) go to 11410                                       1172
11420 continue                                                             1172
11421 if(my(j).ge.mx(i))goto 11422                                         1172
      j=j+1                                                                1172
      if(j.gt.ny) go to 11400                                              1172
      goto 11421                                                           1173
11422 continue                                                             1173
      if(mx(i).eq.my(j)) go to 11410                                       1173
      goto 11381                                                           1174
11410 continue                                                             1174
      s=s+w(mx(i))*x(i)*y(j)                                               1175
      i=i+1                                                                1175
      if(i.gt.nx)goto 11382                                                1175
      j=j+1                                                                1175
      if(j.gt.ny)goto 11382                                                1176
      goto 11381                                                           1177
11382 continue                                                             1177
11400 continue                                                             1177
      dot=s                                                                1178
      return                                                               1179
      end                                                                  1180
      subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,ne,nx,nlam,flmin,ulam   1182 
     *,thr,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)             1183
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1184
      integer jd(*),ia(nx),nin(nlam)                                       1185
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 11441                                    1189
      jerr=10000                                                           1189
      return                                                               1189
11441 continue                                                             1190
      allocate(ww(1:no),stat=jerr)                                         1191
      allocate(ju(1:ni),stat=ierr)                                         1191
      jerr=jerr+ierr                                                       1192
      allocate(vq(1:ni),stat=ierr)                                         1192
      jerr=jerr+ierr                                                       1193
      allocate(xm(1:ni),stat=ierr)                                         1193
      jerr=jerr+ierr                                                       1194
      if(isd .le. 0)goto 11461                                             1194
      allocate(xs(1:ni),stat=ierr)                                         1194
      jerr=jerr+ierr                                                       1194
11461 continue                                                             1195
      if(jerr.ne.0) return                                                 1196
      call chkvars(no,ni,x,ju)                                             1197
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1198
      if(maxval(ju) .gt. 0)goto 11481                                      1198
      jerr=7777                                                            1198
      return                                                               1198
11481 continue                                                             1199
      vq=max(0.0,vp)                                                       1199
      vq=vq*ni/sum(vq)                                                     1200
11490 do 11491 i=1,no                                                      1200
      ww(i)=sum(y(i,:))                                                    1200
      y(i,:)=y(i,:)/ww(i)                                                  1200
11491 continue                                                             1200
11492 continue                                                             1200
      sw=sum(ww)                                                           1200
      ww=ww/sw                                                             1201
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             1202
      if(nc .ne. 1)goto 11511                                              1203
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,ne,nx,nlam,flmin   1205 
     *,ulam,  thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 11521                                                           1206
11511 continue                                                             1207
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,th   1209 
     *r,  isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
11521 continue                                                             1210
11501 continue                                                             1210
      if(jerr.gt.0) return                                                 1210
      dev0=2.0*sw*dev0                                                     1211
11530 do 11531 k=1,lmu                                                     1211
      nk=nin(k)                                                            1212
11540 do 11541 ic=1,nc                                                     1212
      if(isd .le. 0)goto 11561                                             1212
11570 do 11571 l=1,nk                                                      1212
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1212
11571 continue                                                             1212
11572 continue                                                             1212
11561 continue                                                             1213
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1214
11541 continue                                                             1215
11542 continue                                                             1215
11531 continue                                                             1216
11532 continue                                                             1216
      deallocate(ww,ju,vq,xm)                                              1216
      if(isd.gt.0) deallocate(xs)                                          1217
      return                                                               1218
      end                                                                  1219
      subroutine lstandard1 (no,ni,x,w,ju,isd,xm,xs)                       1220
      real x(no,ni),w(no),xm(ni),xs(ni)                                    1220
      integer ju(ni)                                                       1221
11580 do 11581 j=1,ni                                                      1221
      if(ju(j).eq.0)goto 11581                                             1222
      xm(j)=dot_product(w,x(:,j))                                          1222
      x(:,j)=x(:,j)-xm(j)                                                  1223
      if(isd .le. 0)goto 11601                                             1223
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 1223
      x(:,j)=x(:,j)/xs(j)                                                  1223
11601 continue                                                             1224
11581 continue                                                             1225
11582 continue                                                             1225
      return                                                               1226
      end                                                                  1227
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ulam   1229 
     *,shri,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1231 
     *5, devmax=0.999)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    1232
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1233
      integer ju(ni),m(nx),kin(nlam)                                       1234
      real, dimension (:), allocatable :: b,bs,v,r,xv,q                         
      integer, dimension (:), allocatable :: mm                                 
      allocate(b(0:ni),stat=jerr)                                          1239
      allocate(xv(1:ni),stat=ierr)                                         1239
      jerr=jerr+ierr                                                       1240
      allocate(bs(0:ni),stat=ierr)                                         1240
      jerr=jerr+ierr                                                       1241
      allocate(mm(1:ni),stat=ierr)                                         1241
      jerr=jerr+ierr                                                       1242
      allocate(r(1:no),stat=ierr)                                          1242
      jerr=jerr+ierr                                                       1243
      allocate(v(1:no),stat=ierr)                                          1243
      jerr=jerr+ierr                                                       1244
      allocate(q(1:no),stat=ierr)                                          1244
      jerr=jerr+ierr                                                       1245
      if(jerr.ne.0) return                                                 1246
      fmax=log(1.0/pmin-1.0)                                               1246
      fmin=-fmax                                                           1246
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1247
      bta=max(parm,1.0e-3)                                                 1247
      omb=1.0-bta                                                          1248
      q0=dot_product(w,y)                                                  1248
      if(q0 .gt. pmin)goto 11621                                           1248
      jerr=8001                                                            1248
      return                                                               1248
11621 continue                                                             1249
      if(q0 .lt. 1.0-pmin)goto 11641                                       1249
      jerr=9001                                                            1249
      return                                                               1249
11641 continue                                                             1250
      bz=log(q0/(1.0-q0))                                                  1251
      if(nonzero(no,g) .ne. 0)goto 11661                                   1251
      vi=q0*(1.0-q0)                                                       1251
      b(0)=bz                                                              1251
      v=vi*w                                                               1252
      r=w*(y-q0)                                                           1252
      q=q0                                                                 1252
      xmz=vi                                                               1253
      goto 11671                                                           1254
11661 continue                                                             1254
      b(0)=azero(no,y,g,w,jerr)                                            1254
      if(jerr.ne.0) return                                                 1255
      q=1.0/(1.0+exp(-b(0)-g))                                             1255
      v=w*q*(1.0-q)                                                        1255
      r=w*(y-q)                                                            1255
      xmz=sum(v)                                                           1256
11671 continue                                                             1257
11651 continue                                                             1257
      if(isd .le. 0)goto 11691                                             1257
      xv=0.25                                                              1257
      goto 11701                                                           1258
11691 continue                                                             1258
11710 do 11711 j=1,ni                                                      1258
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1258
11711 continue                                                             1258
11712 continue                                                             1258
11701 continue                                                             1259
11681 continue                                                             1259
      dev1=-(bz*q0+log(1.0-q0))                                            1259
      dev0=dev1                                                            1260
11720 do 11721 i=1,no                                                      1260
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1261
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1262
11721 continue                                                             1263
11722 continue                                                             1263
      if(flmin .ge. 1.0)goto 11741                                         1263
      eqs=max(eps,flmin)                                                   1263
      alf=eqs**(1.0/(nlam-1))                                              1263
11741 continue                                                             1264
      m=0                                                                  1264
      mm=0                                                                 1264
      nlp=0                                                                1264
      nin=nlp                                                              1264
      mnl=min(mnlam,nlam)                                                  1264
      bs=0.0                                                               1264
      b(1:ni)=0.0                                                          1265
      shr=shri*dev0                                                        1266
11750 do 11751 ilm=1,nlam                                                  1267
      if(flmin .lt. 1.0)goto 11771                                         1267
      al=ulam(ilm)                                                         1267
      goto 11761                                                           1268
11771 if(ilm .le. 2)goto 11781                                             1268
      al=al*alf                                                            1268
      goto 11761                                                           1269
11781 if(ilm .ne. 1)goto 11791                                             1269
      al=big                                                               1269
      goto 11801                                                           1270
11791 continue                                                             1270
      al=0.0                                                               1271
11810 do 11811 j=1,ni                                                      1271
      if(ju(j).eq.0)goto 11811                                             1271
      if(vp(j).le.0.0)goto 11811                                           1272
      al=max(al,abs(dot_product(r,x(:,j)))/vp(j))                          1273
11811 continue                                                             1274
11812 continue                                                             1274
      al=alf*al/bta                                                        1275
11801 continue                                                             1276
11761 continue                                                             1276
      al2=al*omb                                                           1276
      al1=al*bta                                                           1276
      nit=0                                                                1277
11820 continue                                                             1277
11821 continue                                                             1277
      bs(0)=b(0)                                                           1277
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1278
11830 continue                                                             1278
11831 continue                                                             1278
      nlp=nlp+1                                                            1278
      dlx=0.0                                                              1279
11840 do 11841 k=1,ni                                                      1279
      if(ju(k).eq.0)goto 11841                                             1280
      bk=b(k)                                                              1280
      gk=dot_product(r,x(:,k))                                             1281
      u=gk+xv(k)*b(k)                                                      1281
      au=abs(u)-vp(k)*al1                                                  1282
      if(au .gt. 0.0)goto 11861                                            1282
      b(k)=0.0                                                             1282
      goto 11871                                                           1283
11861 continue                                                             1283
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1283
11871 continue                                                             1284
11851 continue                                                             1284
      d=b(k)-bk                                                            1284
      if(abs(d).le.0.0)goto 11841                                          1284
      dlx=max(dlx,xv(k)*d**2)                                              1285
      r=r-d*v*x(:,k)                                                       1286
      if(mm(k) .ne. 0)goto 11891                                           1286
      nin=nin+1                                                            1286
      if(nin.gt.nx)goto 11842                                              1287
      mm(k)=nin                                                            1287
      m(nin)=k                                                             1288
11891 continue                                                             1289
11841 continue                                                             1290
11842 continue                                                             1290
      if(nin.gt.nx)goto 11832                                              1291
      d=sum(r)/xmz                                                         1292
      if(d .eq. 0.0)goto 11911                                             1292
      b(0)=b(0)+d                                                          1292
      dlx=max(dlx,xmz*d**2)                                                1292
      r=r-d*v                                                              1292
11911 continue                                                             1293
      if(dlx.lt.shr)goto 11832                                             1294
11920 continue                                                             1294
11921 continue                                                             1294
      nlp=nlp+1                                                            1294
      dlx=0.0                                                              1295
11930 do 11931 l=1,nin                                                     1295
      k=m(l)                                                               1295
      bk=b(k)                                                              1295
      gk=dot_product(r,x(:,k))                                             1296
      u=gk+xv(k)*b(k)                                                      1296
      au=abs(u)-vp(k)*al1                                                  1297
      if(au .gt. 0.0)goto 11951                                            1297
      b(k)=0.0                                                             1297
      goto 11961                                                           1298
11951 continue                                                             1298
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1298
11961 continue                                                             1299
11941 continue                                                             1299
      d=b(k)-bk                                                            1299
      if(abs(d).le.0.0)goto 11931                                          1299
      dlx=max(dlx,xv(k)*d**2)                                              1300
      r=r-d*v*x(:,k)                                                       1301
11931 continue                                                             1302
11932 continue                                                             1302
      d=sum(r)/xmz                                                         1303
      if(d .eq. 0.0)goto 11981                                             1303
      b(0)=b(0)+d                                                          1303
      dlx=max(dlx,xmz*d**2)                                                1303
      r=r-d*v                                                              1303
11981 continue                                                             1305
      if(dlx.lt.shr)goto 11922                                             1305
      goto 11921                                                           1306
11922 continue                                                             1306
      goto 11831                                                           1307
11832 continue                                                             1307
      if(nin.gt.nx)goto 11822                                              1308
11990 do 11991 i=1,no                                                      1308
      fi=b(0)+g(i)                                                         1309
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1310
      if(fi .ge. fmin)goto 12011                                           1310
      q(i)=0.0                                                             1310
      goto 12001                                                           1310
12011 if(fi .le. fmax)goto 12021                                           1310
      q(i)=1.0                                                             1310
      goto 12031                                                           1311
12021 continue                                                             1311
      q(i)=1.0/(1.0+exp(-fi))                                              1311
12031 continue                                                             1312
12001 continue                                                             1312
11991 continue                                                             1313
11992 continue                                                             1313
      v=w*q*(1.0-q)                                                        1313
      xmz=sum(v)                                                           1313
      if(xmz.le.vmin)goto 11822                                            1313
      r=w*(y-q)                                                            1314
      if(kopt .ne. 0)goto 12051                                            1315
12060 do 12061 j=1,nin                                                     1315
      xv(m(j))=dot_product(v,x(:,m(j))**2)                                 1315
12061 continue                                                             1316
12062 continue                                                             1316
12051 continue                                                             1317
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 12081                           1317
      ix=0                                                                 1318
12090 do 12091 j=1,nin                                                     1318
      if(xv(m(j))*(b(m(j))-bs(m(j)))**2.lt.shr)goto 12091                  1318
      ix=1                                                                 1318
      goto 12092                                                           1318
12091 continue                                                             1319
12092 continue                                                             1319
      if(ix.eq.0)goto 11822                                                1320
12081 continue                                                             1321
      nit=nit+1                                                            1321
      if(nit .le. maxit)goto 12111                                         1321
      jerr=-ilm                                                            1321
      return                                                               1321
12111 continue                                                             1322
      goto 11821                                                           1323
11822 continue                                                             1323
      if(nin .le. nx)goto 12131                                            1323
      jerr=-10000-ilm                                                      1323
      goto 11752                                                           1323
12131 continue                                                             1324
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1324
      kin(ilm)=nin                                                         1325
      a0(ilm)=b(0)                                                         1325
      alm(ilm)=al                                                          1325
      lmu=ilm                                                              1326
      devi=dev2(no,w,y,q,pmin)                                             1327
      dev(ilm)=(dev1-devi)/dev0                                            1327
      if(xmz.le.vmin)goto 11752                                            1328
      if(ilm.lt.mnl)goto 11751                                             1328
      if(flmin.ge.1.0)goto 11751                                           1329
      me=0                                                                 1329
12140 do 12141 j=1,nin                                                     1329
      if(a(j,ilm).ne.0.0) me=me+1                                          1329
12141 continue                                                             1329
12142 continue                                                             1329
      if(me.gt.ne)goto 11752                                               1330
      if(dev(ilm).gt.devmax)goto 11752                                     1330
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 11752                             1331
11751 continue                                                             1332
11752 continue                                                             1332
      g=log(q/(1.0-q))                                                     1333
      deallocate(b,bs,v,r,xv,q,mm)                                         1334
      return                                                               1335
      end                                                                  1336
      function dev2(n,w,y,p,pmin)                                          1337
      real w(n),y(n),p(n)                                                  1338
      pmax=1.0-pmin                                                        1338
      s=0.0                                                                1339
12150 do 12151 i=1,n                                                       1339
      pi=min(max(pmin,p(i)),pmax)                                          1340
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1341
12151 continue                                                             1342
12152 continue                                                             1342
      dev2=s                                                               1343
      return                                                               1344
      end                                                                  1345
      function azero(n,y,g,q,jerr)                                         1346
      parameter(eps=1.0e-7)                                                1347
      real y(n),g(n),q(n)                                                  1348
      real, dimension (:), allocatable :: e,p,w                                 
      allocate(e(1:n),stat=jerr)                                           1352
      allocate(p(1:n),stat=ierr)                                           1352
      jerr=jerr+ierr                                                       1353
      allocate(w(1:n),stat=ierr)                                           1353
      jerr=jerr+ierr                                                       1354
      if(jerr.ne.0) return                                                 1355
      az=0.0                                                               1355
      e=exp(-g)                                                            1355
      qy=dot_product(q,y)                                                  1355
      p=1.0/(1.0+e)                                                        1356
12160 continue                                                             1356
12161 continue                                                             1356
      w=q*p*(1.0-p)                                                        1357
      d=(qy-dot_product(q,p))/sum(w)                                       1357
      az=az+d                                                              1357
      if(abs(d).lt.eps)goto 12162                                          1358
      ea0=exp(-az)                                                         1358
      p=1.0/(1.0+ea0*e)                                                    1359
      goto 12161                                                           1360
12162 continue                                                             1360
      azero=az                                                             1361
      deallocate(e,p,w)                                                    1362
      return                                                               1363
      end                                                                  1364
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,ne,nx,nlam,flmin,ul   1366 
     *am,shri,  isd,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1368 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam)              1369
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1370
      integer ju(ni),m(nx),kin(nlam)                                       1371
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp                       
      real, dimension (:), allocatable :: di,v,r                                
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is                              
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(r(1:no),stat=ierr)                                          1382
      jerr=jerr+ierr                                                       1383
      allocate(v(1:no),stat=ierr)                                          1383
      jerr=jerr+ierr                                                       1384
      allocate(mm(1:ni),stat=ierr)                                         1384
      jerr=jerr+ierr                                                       1385
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1385
      jerr=jerr+ierr                                                       1386
      allocate(sxp(1:no),stat=ierr)                                        1386
      jerr=jerr+ierr                                                       1387
      allocate(di(1:no),stat=ierr)                                         1387
      jerr=jerr+ierr                                                       1388
      if(jerr.ne.0) return                                                 1389
      pmax=1.0-pmin                                                        1389
      emin=pmin/pmax                                                       1389
      emax=1.0/emin                                                        1390
      pfm=(1.0+pmin)*pmin                                                  1390
      pfx=(1.0-pmin)*pmax                                                  1390
      vmin=pfm*pmax                                                        1391
      bta=max(parm,1.0e-3)                                                 1391
      omb=1.0-bta                                                          1391
      dev1=0.0                                                             1391
      dev0=0.0                                                             1392
12170 do 12171 ic=1,nc                                                     1392
      q0=dot_product(w,y(:,ic))                                            1393
      if(q0 .gt. pmin)goto 12191                                           1393
      jerr =8000+ic                                                        1393
      return                                                               1393
12191 continue                                                             1394
      if(q0 .lt. 1.0-pmin)goto 12211                                       1394
      jerr =9000+ic                                                        1394
      return                                                               1394
12211 continue                                                             1395
      b(0,ic)=log(q0)                                                      1395
      dev1=dev1-q0*b(0,ic)                                                 1395
      b(1:ni,ic)=0.0                                                       1396
12220 do 12221 i=1,no                                                      1396
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1396
12221 continue                                                             1397
12222 continue                                                             1397
12171 continue                                                             1398
12172 continue                                                             1398
      dev0=dev0+dev1                                                       1399
      if(nonzero(no*nc,g) .ne. 0)goto 12241                                1400
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1400
      sxp=0.0                                                              1401
12250 do 12251 ic=1,nc                                                     1401
      q(:,ic)=exp(b(0,ic))                                                 1401
      sxp=sxp+q(:,ic)                                                      1401
12251 continue                                                             1402
12252 continue                                                             1402
      goto 12261                                                           1403
12241 continue                                                             1403
12270 do 12271 i=1,no                                                      1403
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1403
12271 continue                                                             1403
12272 continue                                                             1403
      sxp=0.0                                                              1404
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1404
      if(jerr.ne.0) return                                                 1405
12280 do 12281 ic=1,nc                                                     1405
      q(:,ic)=exp(b(0,ic)+g(:,ic))                                         1405
      sxp=sxp+q(:,ic)                                                      1405
12281 continue                                                             1406
12282 continue                                                             1406
12261 continue                                                             1407
12231 continue                                                             1407
      if(isd .le. 0)goto 12301                                             1407
      xv=0.25                                                              1407
      goto 12311                                                           1408
12301 continue                                                             1408
12320 do 12321 j=1,ni                                                      1408
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1408
12321 continue                                                             1408
12322 continue                                                             1408
12311 continue                                                             1409
12291 continue                                                             1409
      if(flmin .ge. 1.0)goto 12341                                         1409
      eqs=max(eps,flmin)                                                   1409
      alf=eqs**(1.0/(nlam-1))                                              1409
12341 continue                                                             1410
      m=0                                                                  1410
      mm=0                                                                 1410
      nin=0                                                                1410
      nlp=0                                                                1410
      mnl=min(mnlam,nlam)                                                  1410
      bs=0.0                                                               1410
      shr=shri*dev0                                                        1411
12350 do 12351 ilm=1,nlam                                                  1412
      if(flmin .lt. 1.0)goto 12371                                         1412
      al=ulam(ilm)                                                         1412
      goto 12361                                                           1413
12371 if(ilm .le. 2)goto 12381                                             1413
      al=al*alf                                                            1413
      goto 12361                                                           1414
12381 if(ilm .ne. 1)goto 12391                                             1414
      al=big                                                               1414
      goto 12401                                                           1415
12391 continue                                                             1415
      al=0.0                                                               1416
12410 do 12411 ic=1,nc                                                     1416
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1417
12420 do 12421 j=1,ni                                                      1417
      if(ju(j).eq.0)goto 12421                                             1417
      if(vp(j).le.0.0)goto 12421                                           1418
      al=max(al,abs(dot_product(r,x(:,j)))/vp(j))                          1419
12421 continue                                                             1420
12422 continue                                                             1420
12411 continue                                                             1421
12412 continue                                                             1421
      al=alf*al/bta                                                        1422
12401 continue                                                             1423
12361 continue                                                             1423
      al2=al*omb                                                           1423
      al1=al*bta                                                           1423
      nit=0                                                                1424
12430 continue                                                             1424
12431 continue                                                             1424
      ix=0                                                                 1424
      jx=ix                                                                1424
      ig=0                                                                 1425
12440 do 12441 ic=1,nc                                                     1425
      bs(0,ic)=b(0,ic)                                                     1426
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1427
      xmz=0.0                                                              1428
12450 do 12451 i=1,no                                                      1428
      pic=q(i,ic)/sxp(i)                                                   1429
      if(pic .ge. pfm)goto 12471                                           1429
      pic=0.0                                                              1429
      v(i)=0.0                                                             1429
      goto 12461                                                           1430
12471 if(pic .le. pfx)goto 12481                                           1430
      pic=1.0                                                              1430
      v(i)=0.0                                                             1430
      goto 12491                                                           1431
12481 continue                                                             1431
      v(i)=w(i)*pic*(1.0-pic)                                              1431
      xmz=xmz+v(i)                                                         1431
12491 continue                                                             1432
12461 continue                                                             1432
      r(i)=w(i)*(y(i,ic)-pic)                                              1433
12451 continue                                                             1434
12452 continue                                                             1434
      if(xmz.le.vmin)goto 12441                                            1434
      ig=1                                                                 1435
      if(kopt .ne. 0)goto 12511                                            1436
12520 do 12521 j=1,nin                                                     1436
      xv(m(j),ic)=dot_product(v,x(:,m(j))**2)                              1436
12521 continue                                                             1437
12522 continue                                                             1437
12511 continue                                                             1438
12530 continue                                                             1438
12531 continue                                                             1438
      nlp=nlp+1                                                            1438
      dlx=0.0                                                              1439
12540 do 12541 k=1,ni                                                      1439
      if(ju(k).eq.0)goto 12541                                             1440
      bk=b(k,ic)                                                           1440
      gk=dot_product(r,x(:,k))                                             1441
      u=gk+xv(k,ic)*b(k,ic)                                                1441
      au=abs(u)-vp(k)*al1                                                  1442
      if(au .gt. 0.0)goto 12561                                            1442
      b(k,ic)=0.0                                                          1442
      goto 12571                                                           1443
12561 continue                                                             1443
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1443
12571 continue                                                             1444
12551 continue                                                             1444
      d=b(k,ic)-bk                                                         1444
      if(abs(d).le.0.0)goto 12541                                          1445
      dlx=max(dlx,xv(k,ic)*d**2)                                           1445
      r=r-d*v*x(:,k)                                                       1446
      if(mm(k) .ne. 0)goto 12591                                           1446
      nin=nin+1                                                            1447
      if(nin .le. nx)goto 12611                                            1447
      jx=1                                                                 1447
      goto 12542                                                           1447
12611 continue                                                             1448
      mm(k)=nin                                                            1448
      m(nin)=k                                                             1449
12591 continue                                                             1450
12541 continue                                                             1451
12542 continue                                                             1451
      if(jx.gt.0)goto 12532                                                1452
      d=sum(r)/xmz                                                         1453
      if(d .eq. 0.0)goto 12631                                             1453
      b(0,ic)=b(0,ic)+d                                                    1453
      dlx=max(dlx,xmz*d**2)                                                1453
      r=r-d*v                                                              1453
12631 continue                                                             1454
      if(dlx.lt.shr)goto 12532                                             1455
12640 continue                                                             1455
12641 continue                                                             1455
      nlp=nlp+1                                                            1455
      dlx=0.0                                                              1456
12650 do 12651 l=1,nin                                                     1456
      k=m(l)                                                               1456
      bk=b(k,ic)                                                           1457
      gk=dot_product(r,x(:,k))                                             1458
      u=gk+xv(k,ic)*b(k,ic)                                                1458
      au=abs(u)-vp(k)*al1                                                  1459
      if(au .gt. 0.0)goto 12671                                            1459
      b(k,ic)=0.0                                                          1459
      goto 12681                                                           1460
12671 continue                                                             1460
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1460
12681 continue                                                             1461
12661 continue                                                             1461
      d=b(k,ic)-bk                                                         1461
      if(abs(d).le.0.0)goto 12651                                          1462
      dlx=max(dlx,xv(k,ic)*d**2)                                           1462
      r=r-d*v*x(:,k)                                                       1463
12651 continue                                                             1464
12652 continue                                                             1464
      d=sum(r)/xmz                                                         1465
      if(d .eq. 0.0)goto 12701                                             1465
      b(0,ic)=b(0,ic)+d                                                    1466
      dlx=max(dlx,xmz*d**2)                                                1466
      r=r-d*v                                                              1467
12701 continue                                                             1468
      if(dlx.lt.shr)goto 12642                                             1468
      goto 12641                                                           1469
12642 continue                                                             1469
      goto 12531                                                           1470
12532 continue                                                             1470
      if(jx.gt.0)goto 12442                                                1471
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                            1472
      if(ix .ne. 0)goto 12721                                              1473
12730 do 12731 j=1,nin                                                     1474
      if(xv(m(j),ic)*(b(m(j),ic)-bs(m(j),ic))**2 .le. shr)goto 12751       1474
      ix=1                                                                 1474
      goto 12732                                                           1474
12751 continue                                                             1475
12731 continue                                                             1476
12732 continue                                                             1476
12721 continue                                                             1477
12760 do 12761 i=1,no                                                      1477
      fi=b(0,ic)+g(i,ic)                                                   1479
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1480
      fi=min(max(exmn,fi),exmx)                                            1480
      sxp(i)=sxp(i)-q(i,ic)                                                1481
      q(i,ic)=min(max(emin*sxp(i),exp(dble(fi))),emax*sxp(i))              1482
      sxp(i)=sxp(i)+q(i,ic)                                                1483
12761 continue                                                             1484
12762 continue                                                             1484
12441 continue                                                             1485
12442 continue                                                             1485
      s=-sum(b(0,:))/nc                                                    1485
      b(0,:)=b(0,:)+s                                                      1485
      di=s                                                                 1486
12770 do 12771 j=1,nin                                                     1486
      l=m(j)                                                               1487
      if(vp(l) .gt. 0.0)goto 12791                                         1487
      s=sum(b(l,:))/nc                                                     1487
      goto 12801                                                           1488
12791 continue                                                             1488
      s=elc(parm,nc,b(l,:),is)                                             1488
12801 continue                                                             1489
12781 continue                                                             1489
      b(l,:)=b(l,:)-s                                                      1489
      di=di-s*x(:,l)                                                       1490
12771 continue                                                             1491
12772 continue                                                             1491
      di=exp(di)                                                           1491
      sxp=sxp*di                                                           1491
12810 do 12811 ic=1,nc                                                     1491
      q(:,ic)=q(:,ic)*di                                                   1491
12811 continue                                                             1492
12812 continue                                                             1492
      if(jx.gt.0)goto 12432                                                1492
      if(ix.eq.0)goto 12432                                                1492
      if(ig.eq.0)goto 12432                                                1493
      nit=nit+1                                                            1493
      if(nit .le. maxit)goto 12831                                         1493
      jerr=-ilm                                                            1493
      return                                                               1493
12831 continue                                                             1494
      goto 12431                                                           1495
12432 continue                                                             1495
      if(jx .le. 0)goto 12851                                              1495
      jerr=-10000-ilm                                                      1495
      goto 12352                                                           1495
12851 continue                                                             1495
      devi=0.0                                                             1496
12860 do 12861 ic=1,nc                                                     1497
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1497
      a0(ic,ilm)=b(0,ic)                                                   1498
12870 do 12871 i=1,no                                                      1498
      if(y(i,ic).le.0.0)goto 12871                                         1499
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1500
12871 continue                                                             1501
12872 continue                                                             1501
12861 continue                                                             1502
12862 continue                                                             1502
      kin(ilm)=nin                                                         1502
      alm(ilm)=al                                                          1502
      lmu=ilm                                                              1503
      dev(ilm)=(dev1-devi)/dev0                                            1503
      if(ig.eq.0)goto 12352                                                1504
      if(ilm.lt.mnl)goto 12351                                             1504
      if(flmin.ge.1.0)goto 12351                                           1505
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 12352             1506
      if(dev(ilm).gt.devmax)goto 12352                                     1506
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12352                             1507
12351 continue                                                             1508
12352 continue                                                             1508
      g=log(q)                                                             1508
12880 do 12881 i=1,no                                                      1508
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1508
12881 continue                                                             1509
12882 continue                                                             1509
      deallocate(sxp,b,bs,v,r,xv,q,mm,is)                                  1510
      return                                                               1511
      end                                                                  1512
      subroutine kazero(kk,n,y,g,q,az,jerr)                                1513
      parameter(eps=1.0e-7)                                                1514
      real y(n,kk),g(n,kk),q(n),az(kk)                                     1515
      real, dimension (:), allocatable :: s                                     
      real, dimension (:,:), allocatable :: e                                   
      allocate(e(1:n,1:kk),stat=jerr)                                           
      allocate(s(1:n),stat=ierr)                                           1520
      jerr=jerr+ierr                                                       1521
      if(jerr.ne.0) return                                                 1522
      az=0.0                                                               1522
      e=exp(g)                                                             1522
12890 do 12891 i=1,n                                                       1522
      s(i)=sum(e(i,:))                                                     1522
12891 continue                                                             1523
12892 continue                                                             1523
12900 continue                                                             1523
12901 continue                                                             1523
      dm=0.0                                                               1524
12910 do 12911 k=1,kk                                                      1524
      t=0.0                                                                1524
      u=t                                                                  1525
12920 do 12921 i=1,n                                                       1525
      pik=e(i,k)/s(i)                                                      1526
      t=t+q(i)*(y(i,k)-pik)                                                1526
      u=u+q(i)*pik*(1.0-pik)                                               1527
12921 continue                                                             1528
12922 continue                                                             1528
      d=t/u                                                                1528
      az(k)=az(k)+d                                                        1528
      ed=exp(d)                                                            1528
      dm=max(dm,abs(d))                                                    1529
12930 do 12931 i=1,n                                                       1529
      z=e(i,k)                                                             1529
      e(i,k)=z*ed                                                          1529
      s(i)=s(i)-z+e(i,k)                                                   1529
12931 continue                                                             1530
12932 continue                                                             1530
12911 continue                                                             1531
12912 continue                                                             1531
      if(dm.lt.eps)goto 12902                                              1531
      goto 12901                                                           1532
12902 continue                                                             1532
      az=az-sum(az)/kk                                                     1533
      deallocate(e,s)                                                      1534
      return                                                               1535
      end                                                                  1536
      function elc(parm,n,a,m)                                             1537
      real a(n)                                                            1537
      integer m(n)                                                         1538
      fn=n                                                                 1538
      am=sum(a)/fn                                                         1539
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 12951                       1539
      elc=am                                                               1539
      return                                                               1539
12951 continue                                                             1540
12960 do 12961 i=1,n                                                       1540
      m(i)=i                                                               1540
12961 continue                                                             1540
12962 continue                                                             1540
      call psort7(a,m,1,n)                                                 1541
      if(a(m(1)) .ne. a(m(n)))goto 12981                                   1541
      elc=a(1)                                                             1541
      return                                                               1541
12981 continue                                                             1542
      if(mod(n,2) .ne. 1)goto 13001                                        1542
      ad=a(m(n/2+1))                                                       1542
      goto 13011                                                           1543
13001 continue                                                             1543
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       1543
13011 continue                                                             1544
12991 continue                                                             1544
      if(parm .ne. 1.0)goto 13031                                          1544
      elc=ad                                                               1544
      return                                                               1544
13031 continue                                                             1545
      b1=min(am,ad)                                                        1545
      b2=max(am,ad)                                                        1545
      k2=1                                                                 1546
13040 continue                                                             1546
13041 if(a(m(k2)).gt.b1)goto 13042                                         1546
      k2=k2+1                                                              1546
      goto 13041                                                           1546
13042 continue                                                             1546
      k1=k2-1                                                              1547
13050 continue                                                             1547
13051 if(a(m(k2)).ge.b2)goto 13052                                         1547
      k2=k2+1                                                              1547
      goto 13051                                                           1548
13052 continue                                                             1548
      r=parm/((1.0-parm)*fn)                                               1548
      is=0                                                                 1548
      sm=n-2*(k1-1)                                                        1549
13060 do 13061 k=k1,k2-1                                                   1549
      sm=sm-2.0                                                            1549
      s=r*sm+am                                                            1550
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 13081                   1550
      is=k                                                                 1550
      goto 13062                                                           1550
13081 continue                                                             1551
13061 continue                                                             1552
13062 continue                                                             1552
      if(is .eq. 0)goto 13101                                              1552
      elc=s                                                                1552
      return                                                               1552
13101 continue                                                             1552
      r2=2.0*r                                                             1552
      s1=a(m(k1))                                                          1552
      am2=2.0*am                                                           1553
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    1553
      elc=s1                                                               1554
13110 do 13111 k=k1+1,k2                                                   1554
      s=a(m(k))                                                            1554
      if(s.eq.s1)goto 13111                                                1555
      c=r2*sum(abs(a-s))+s*(s-am2)                                         1556
      if(c .ge. cri)goto 13131                                             1556
      cri=c                                                                1556
      elc=s                                                                1556
13131 continue                                                             1556
      s1=s                                                                 1557
13111 continue                                                             1558
13112 continue                                                             1558
      return                                                               1559
      end                                                                  1560
      function nintot(ni,nx,nc,a,m,nin,is)                                 1561
      real a(nx,nc)                                                        1561
      integer m(nx),is(ni)                                                 1562
      is=0                                                                 1562
      nintot=0                                                             1563
13140 do 13141 ic=1,nc                                                     1563
13150 do 13151 j=1,nin                                                     1563
      k=m(j)                                                               1563
      if(is(k).ne.0)goto 13151                                             1564
      if(a(j,ic).eq.0.0)goto 13151                                         1564
      is(k)=k                                                              1564
      nintot=nintot+1                                                      1565
13151 continue                                                             1565
13152 continue                                                             1565
13141 continue                                                             1566
13142 continue                                                             1566
      return                                                               1567
      end                                                                  1568
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             1569
      real ca(nx,nc),a(ni,nc)                                              1569
      integer ia(nx)                                                       1570
      a=0.0                                                                1571
13160 do 13161 ic=1,nc                                                     1571
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            1571
13161 continue                                                             1572
13162 continue                                                             1572
      return                                                               1573
      end                                                                  1574
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      1575
      real a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                             1575
      integer ia(nx)                                                       1576
13170 do 13171 i=1,nt                                                      1576
13180 do 13181 ic=1,nc                                                     1576
      ans(ic,i)=a0(ic)                                                     1578
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   1579 
     *:nin)))
13181 continue                                                             1579
13182 continue                                                             1579
13171 continue                                                             1580
13172 continue                                                             1580
      return                                                               1581
      end                                                                  1582
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,ne,nx,nlam,fl   1584 
     *min,  ulam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      real x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)                 1585
      real ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                  1586
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1587
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 13201                                    1591
      jerr=10000                                                           1591
      return                                                               1591
13201 continue                                                             1592
      allocate(ww(1:no),stat=jerr)                                         1593
      allocate(ju(1:ni),stat=ierr)                                         1593
      jerr=jerr+ierr                                                       1594
      allocate(vq(1:ni),stat=ierr)                                         1594
      jerr=jerr+ierr                                                       1595
      allocate(xm(1:ni),stat=ierr)                                         1595
      jerr=jerr+ierr                                                       1596
      allocate(xs(1:ni),stat=ierr)                                         1596
      jerr=jerr+ierr                                                       1597
      if(jerr.ne.0) return                                                 1598
      call spchkvars(no,ni,x,ix,ju)                                        1599
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1600
      if(maxval(ju) .gt. 0)goto 13221                                      1600
      jerr=7777                                                            1600
      return                                                               1600
13221 continue                                                             1601
      vq=max(0.0,vp)                                                       1601
      vq=vq*ni/sum(vq)                                                     1602
13230 do 13231 i=1,no                                                      1602
      ww(i)=sum(y(i,:))                                                    1602
      y(i,:)=y(i,:)/ww(i)                                                  1602
13231 continue                                                             1602
13232 continue                                                             1602
      sw=sum(ww)                                                           1602
      ww=ww/sw                                                             1603
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     1604
      if(nc .ne. 1)goto 13251                                              1605
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,ne,nx,n   1607 
     *lam,flmin,  ulam,thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,de
     *v,alm,nlp,jerr)
      goto 13261                                                           1608
13251 continue                                                             1609
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmi   1611 
     *n,ulam,  thr,isd,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nl
     *p,jerr)
13261 continue                                                             1612
13241 continue                                                             1612
      if(jerr.gt.0) return                                                 1612
      dev0=2.0*sw*dev0                                                     1613
13270 do 13271 k=1,lmu                                                     1613
      nk=nin(k)                                                            1614
13280 do 13281 ic=1,nc                                                     1614
      if(isd .le. 0)goto 13301                                             1614
13310 do 13311 l=1,nk                                                      1614
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1614
13311 continue                                                             1614
13312 continue                                                             1614
13301 continue                                                             1615
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1616
13281 continue                                                             1617
13282 continue                                                             1617
13271 continue                                                             1618
13272 continue                                                             1618
      deallocate(ww,ju,vq,xm,xs)                                           1619
      return                                                               1620
      end                                                                  1621
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,xm,xs)                1622
      real x(*),w(no),xm(ni),xs(ni)                                        1622
      integer ix(*),jx(*),ju(ni)                                           1623
13320 do 13321 j=1,ni                                                      1623
      if(ju(j).eq.0)goto 13321                                             1623
      jb=ix(j)                                                             1623
      je=ix(j+1)-1                                                         1624
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1625
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   1626 
     *)**2)
13321 continue                                                             1627
13322 continue                                                             1627
      if(isd.eq.0) xs=1.0                                                  1628
      return                                                               1629
      end                                                                  1630
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam,   1632 
     *  flmin,ulam,shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,alm
     *,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1634 
     *5, devmax=0.999)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        1635
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                         1636
      real xb(ni),xs(ni)                                                   1636
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1637
      real, dimension (:), allocatable :: xm,b,bs,v,r,sc,xv,q                   
      integer, dimension (:), allocatable :: mm                                 
      allocate(b(0:ni),stat=jerr)                                          1642
      allocate(xm(0:ni),stat=ierr)                                         1642
      jerr=jerr+ierr                                                       1643
      allocate(xv(1:ni),stat=ierr)                                         1643
      jerr=jerr+ierr                                                       1644
      allocate(bs(0:ni),stat=ierr)                                         1644
      jerr=jerr+ierr                                                       1645
      allocate(mm(1:ni),stat=ierr)                                         1645
      jerr=jerr+ierr                                                       1646
      allocate(q(1:no),stat=ierr)                                          1646
      jerr=jerr+ierr                                                       1647
      allocate(r(1:no),stat=ierr)                                          1647
      jerr=jerr+ierr                                                       1648
      allocate(v(1:no),stat=ierr)                                          1648
      jerr=jerr+ierr                                                       1649
      allocate(sc(1:no),stat=ierr)                                         1649
      jerr=jerr+ierr                                                       1650
      if(jerr.ne.0) return                                                 1651
      fmax=log(1.0/pmin-1.0)                                               1651
      fmin=-fmax                                                           1651
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1652
      bta=max(parm,1.0e-3)                                                 1652
      omb=1.0-bta                                                          1653
      q0=dot_product(w,y)                                                  1653
      if(q0 .gt. pmin)goto 13341                                           1653
      jerr=8001                                                            1653
      return                                                               1653
13341 continue                                                             1654
      if(q0 .lt. 1.0-pmin)goto 13361                                       1654
      jerr=9001                                                            1654
      return                                                               1654
13361 continue                                                             1654
      bz=log(q0/(1.0-q0))                                                  1655
      if(nonzero(no,g) .ne. 0)goto 13381                                   1655
      vi=q0*(1.0-q0)                                                       1655
      b(0)=bz                                                              1655
      v=vi*w                                                               1656
      r=w*(y-q0)                                                           1656
      q=q0                                                                 1656
      xm(0)=vi                                                             1657
      goto 13391                                                           1658
13381 continue                                                             1658
      b(0)=azero(no,y,g,w,jerr)                                            1658
      if(jerr.ne.0) return                                                 1659
      q=1.0/(1.0+exp(-b(0)-g))                                             1659
      v=w*q*(1.0-q)                                                        1659
      r=w*(y-q)                                                            1659
      xm(0)=sum(v)                                                         1660
13391 continue                                                             1661
13371 continue                                                             1661
      if(isd .le. 0)goto 13411                                             1661
      xv=0.25                                                              1661
      goto 13421                                                           1662
13411 continue                                                             1663
13430 do 13431 j=1,ni                                                      1663
      if(ju(j).eq.0)goto 13431                                             1663
      jb=ix(j)                                                             1663
      je=ix(j+1)-1                                                         1664
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          1665
13431 continue                                                             1666
13432 continue                                                             1666
13421 continue                                                             1667
13401 continue                                                             1667
      b(1:ni)=0.0                                                          1667
      dev1=-(bz*q0+log(1.0-q0))                                            1667
      dev0=dev1                                                            1668
13440 do 13441 i=1,no                                                      1668
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1669
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1670
13441 continue                                                             1671
13442 continue                                                             1671
      if(flmin .ge. 1.0)goto 13461                                         1671
      eqs=max(eps,flmin)                                                   1671
      alf=eqs**(1.0/(nlam-1))                                              1671
13461 continue                                                             1672
      m=0                                                                  1672
      mm=0                                                                 1672
      nin=0                                                                1672
      o=0.0                                                                1672
      svr=o                                                                1672
      mnl=min(mnlam,nlam)                                                  1672
      bs=0.0                                                               1672
      nlp=0                                                                1672
      nin=nlp                                                              1673
      shr=shri*dev0                                                        1674
13470 do 13471 ilm=1,nlam                                                  1675
      if(flmin .lt. 1.0)goto 13491                                         1675
      al=ulam(ilm)                                                         1675
      goto 13481                                                           1676
13491 if(ilm .le. 2)goto 13501                                             1676
      al=al*alf                                                            1676
      goto 13481                                                           1677
13501 if(ilm .ne. 1)goto 13511                                             1677
      al=big                                                               1677
      goto 13521                                                           1678
13511 continue                                                             1678
      al=0.0                                                               1679
13530 do 13531 j=1,ni                                                      1679
      if(ju(j).eq.0)goto 13531                                             1679
      if(vp(j).le.0.0)goto 13531                                           1680
      jb=ix(j)                                                             1680
      je=ix(j+1)-1                                                         1680
      jn=ix(j+1)-ix(j)                                                     1681
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1682
      gj=dot_product(sc(1:jn),x(jb:je))                                    1683
      gj=(gj-svr*xb(j))/xs(j)                                              1684
      al=max(al,abs(gj)/vp(j))                                             1685
13531 continue                                                             1686
13532 continue                                                             1686
      al=alf*al/bta                                                        1687
13521 continue                                                             1688
13481 continue                                                             1688
      al2=al*omb                                                           1688
      al1=al*bta                                                           1688
      nit=0                                                                1689
13540 continue                                                             1689
13541 continue                                                             1689
      bs(0)=b(0)                                                           1689
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1690
13550 continue                                                             1690
13551 continue                                                             1690
      nlp=nlp+1                                                            1690
      dlx=0.0                                                              1691
13560 do 13561 k=1,ni                                                      1691
      if(ju(k).eq.0)goto 13561                                             1692
      jb=ix(k)                                                             1692
      je=ix(k+1)-1                                                         1692
      jn=ix(k+1)-ix(k)                                                     1692
      bk=b(k)                                                              1693
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1694
      gk=dot_product(sc(1:jn),x(jb:je))                                    1695
      gk=(gk-svr*xb(k))/xs(k)                                              1696
      u=gk+xv(k)*b(k)                                                      1696
      au=abs(u)-vp(k)*al1                                                  1697
      if(au .gt. 0.0)goto 13581                                            1697
      b(k)=0.0                                                             1697
      goto 13591                                                           1698
13581 continue                                                             1698
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1698
13591 continue                                                             1699
13571 continue                                                             1699
      d=b(k)-bk                                                            1699
      if(abs(d).le.0.0)goto 13561                                          1699
      dlx=max(dlx,xv(k)*d**2)                                              1700
      if(mm(k) .ne. 0)goto 13611                                           1700
      nin=nin+1                                                            1700
      if(nin.gt.nx)goto 13562                                              1701
      mm(k)=nin                                                            1701
      m(nin)=k                                                             1701
      sc(1:jn)=v(jx(jb:je))                                                1702
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 1703
13611 continue                                                             1704
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1705
      o=o+d*(xb(k)/xs(k))                                                  1706
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1707
13561 continue                                                             1708
13562 continue                                                             1708
      if(nin.gt.nx)goto 13552                                              1709
      d=svr/xm(0)                                                          1710
      if(d .eq. 0.0)goto 13631                                             1710
      b(0)=b(0)+d                                                          1710
      dlx=max(dlx,xm(0)*d**2)                                              1710
      r=r-d*v                                                              1710
13631 continue                                                             1711
      svr=svr-d*xm(0)                                                      1711
      if(dlx.lt.shr)goto 13552                                             1712
13640 continue                                                             1712
13641 continue                                                             1712
      nlp=nlp+1                                                            1712
      dlx=0.0                                                              1713
13650 do 13651 l=1,nin                                                     1713
      k=m(l)                                                               1713
      jb=ix(k)                                                             1713
      je=ix(k+1)-1                                                         1714
      jn=ix(k+1)-ix(k)                                                     1714
      bk=b(k)                                                              1715
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 1716
      gk=dot_product(sc(1:jn),x(jb:je))                                    1717
      gk=(gk-svr*xb(k))/xs(k)                                              1718
      u=gk+xv(k)*b(k)                                                      1718
      au=abs(u)-vp(k)*al1                                                  1719
      if(au .gt. 0.0)goto 13671                                            1719
      b(k)=0.0                                                             1719
      goto 13681                                                           1720
13671 continue                                                             1720
      b(k)=sign(au,u)/(xv(k)+vp(k)*al2)                                    1720
13681 continue                                                             1721
13661 continue                                                             1721
      d=b(k)-bk                                                            1721
      if(abs(d).le.0.0)goto 13651                                          1721
      dlx=max(dlx,xv(k)*d**2)                                              1722
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1723
      o=o+d*(xb(k)/xs(k))                                                  1724
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1725
13651 continue                                                             1726
13652 continue                                                             1726
      d=svr/xm(0)                                                          1727
      if(d .eq. 0.0)goto 13701                                             1727
      b(0)=b(0)+d                                                          1727
      dlx=max(dlx,xm(0)*d**2)                                              1727
      r=r-d*v                                                              1727
13701 continue                                                             1728
      svr=svr-d*xm(0)                                                      1729
      if(dlx.lt.shr)goto 13642                                             1729
      goto 13641                                                           1730
13642 continue                                                             1730
      goto 13551                                                           1731
13552 continue                                                             1731
      if(nin.gt.nx)goto 13542                                              1732
      sc=b(0)                                                              1732
      b0=0.0                                                               1733
13710 do 13711 j=1,nin                                                     1733
      l=m(j)                                                               1733
      jb=ix(l)                                                             1733
      je=ix(l+1)-1                                                         1734
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      1735
      b0=b0-b(l)*xb(l)/xs(l)                                               1736
13711 continue                                                             1737
13712 continue                                                             1737
      sc=sc+b0                                                             1738
13720 do 13721 i=1,no                                                      1738
      fi=sc(i)+g(i)                                                        1739
      if(fi .ge. fmin)goto 13741                                           1739
      q(i)=0.0                                                             1739
      goto 13731                                                           1739
13741 if(fi .le. fmax)goto 13751                                           1739
      q(i)=1.0                                                             1739
      goto 13761                                                           1740
13751 continue                                                             1740
      q(i)=1.0/(1.0+exp(-fi))                                              1740
13761 continue                                                             1741
13731 continue                                                             1741
13721 continue                                                             1742
13722 continue                                                             1742
      v=w*q*(1.0-q)                                                        1742
      xm(0)=sum(v)                                                         1742
      if(xm(0).lt.vmin)goto 13542                                          1743
      r=w*(y-q)                                                            1743
      svr=sum(r)                                                           1743
      o=0.0                                                                1744
13770 do 13771 l=1,nin                                                     1744
      j=m(l)                                                               1745
      jb=ix(j)                                                             1745
      je=ix(j+1)-1                                                         1745
      jn=ix(j+1)-ix(j)                                                     1746
      sc(1:jn)=v(jx(jb:je))                                                1747
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 1748
      if(kopt .ne. 0)goto 13791                                            1749
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              1750
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                1751
13791 continue                                                             1752
13771 continue                                                             1753
13772 continue                                                             1753
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 13811                         1753
      kx=0                                                                 1754
13820 do 13821 j=1,nin                                                     1754
      if(xv(m(j))*(b(m(j))-bs(m(j)))**2.lt.shr)goto 13821                  1754
      kx=1                                                                 1754
      goto 13822                                                           1754
13821 continue                                                             1755
13822 continue                                                             1755
      if(kx.eq.0)goto 13542                                                1756
13811 continue                                                             1757
      nit=nit+1                                                            1757
      if(nit .le. maxit)goto 13841                                         1757
      jerr=-ilm                                                            1757
      return                                                               1757
13841 continue                                                             1758
      goto 13541                                                           1759
13542 continue                                                             1759
      if(nin .le. nx)goto 13861                                            1759
      jerr=-10000-ilm                                                      1759
      goto 13472                                                           1759
13861 continue                                                             1760
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1760
      kin(ilm)=nin                                                         1761
      a0(ilm)=b(0)                                                         1761
      alm(ilm)=al                                                          1761
      lmu=ilm                                                              1762
      devi=dev2(no,w,y,q,pmin)                                             1763
      dev(ilm)=(dev1-devi)/dev0                                            1764
      if(ilm.lt.mnl)goto 13471                                             1764
      if(flmin.ge.1.0)goto 13471                                           1765
      me=0                                                                 1765
13870 do 13871 j=1,nin                                                     1765
      if(a(j,ilm).ne.0.0) me=me+1                                          1765
13871 continue                                                             1765
13872 continue                                                             1765
      if(me.gt.ne)goto 13472                                               1766
      if(dev(ilm).gt.devmax)goto 13472                                     1766
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13472                             1767
      if(xm(0).lt.vmin)goto 13472                                          1768
13471 continue                                                             1769
13472 continue                                                             1769
      g=log(q/(1.0-q))                                                     1770
      deallocate(xm,b,bs,v,r,sc,xv,q,mm)                                   1771
      return                                                               1772
      end                                                                  1773
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,ne,nx,nlam   1775 
     *,flmin,ulam,  shri,isd,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev0,dev,al
     *m,nlp,jerr)
      parameter(sml=1.0e-5, eps=1.0e-6, pmin=1.0e-5,  big=9.9e30, mnlam=   1777 
     *5, devmax=0.999, exmx=250.0, exmn=-exmx)
      real x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb(ni),xs(ni)    1778
      real a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam)                   1779
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           1780
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp                       
      real, dimension (:), allocatable :: sc,xm,v,r                             
      real, dimension (:,:), allocatable :: b,bs,xv                             
      integer, dimension (:), allocatable :: mm,is                              
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      allocate(xv(1:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(bs(0:ni,1:nc),stat=ierr); jerr=jerr+ierr                         
      allocate(q(1:no,1:nc),stat=ierr); jerr=jerr+ierr                          
      allocate(xm(0:ni),stat=ierr)                                         1791
      jerr=jerr+ierr                                                       1792
      allocate(r(1:no),stat=ierr)                                          1792
      jerr=jerr+ierr                                                       1793
      allocate(v(1:no),stat=ierr)                                          1793
      jerr=jerr+ierr                                                       1794
      allocate(mm(1:ni),stat=ierr)                                         1794
      jerr=jerr+ierr                                                       1795
      allocate(is(1:max(nc,ni)),stat=ierr)                                 1795
      jerr=jerr+ierr                                                       1796
      allocate(sxp(1:no),stat=ierr)                                        1796
      jerr=jerr+ierr                                                       1797
      allocate(sc(1:no),stat=ierr)                                         1797
      jerr=jerr+ierr                                                       1798
      if(jerr.ne.0) return                                                 1799
      pmax=1.0-pmin                                                        1799
      emin=pmin/pmax                                                       1799
      emax=1.0/emin                                                        1800
      pfm=(1.0+pmin)*pmin                                                  1800
      pfx=(1.0-pmin)*pmax                                                  1800
      vmin=pfm*pmax                                                        1801
      bta=max(parm,1.0e-3)                                                 1801
      omb=1.0-bta                                                          1801
      dev1=0.0                                                             1801
      dev0=0.0                                                             1802
13880 do 13881 ic=1,nc                                                     1802
      q0=dot_product(w,y(:,ic))                                            1803
      if(q0 .gt. pmin)goto 13901                                           1803
      jerr =8000+ic                                                        1803
      return                                                               1803
13901 continue                                                             1804
      if(q0 .lt. 1.0-pmin)goto 13921                                       1804
      jerr =9000+ic                                                        1804
      return                                                               1804
13921 continue                                                             1805
      b(1:ni,ic)=0.0                                                       1805
      b(0,ic)=log(q0)                                                      1805
      dev1=dev1-q0*b(0,ic)                                                 1806
13930 do 13931 i=1,no                                                      1806
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1806
13931 continue                                                             1807
13932 continue                                                             1807
13881 continue                                                             1808
13882 continue                                                             1808
      dev0=dev0+dev1                                                       1809
      if(nonzero(no*nc,g) .ne. 0)goto 13951                                1810
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1810
      sxp=0.0                                                              1811
13960 do 13961 ic=1,nc                                                     1811
      q(:,ic)=exp(b(0,ic))                                                 1811
      sxp=sxp+q(:,ic)                                                      1811
13961 continue                                                             1812
13962 continue                                                             1812
      goto 13971                                                           1813
13951 continue                                                             1813
13980 do 13981 i=1,no                                                      1813
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1813
13981 continue                                                             1813
13982 continue                                                             1813
      sxp=0.0                                                              1814
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1814
      if(jerr.ne.0) return                                                 1815
13990 do 13991 ic=1,nc                                                     1815
      q(:,ic)=exp(b(0,ic)+g(:,ic))                                         1815
      sxp=sxp+q(:,ic)                                                      1815
13991 continue                                                             1816
13992 continue                                                             1816
13971 continue                                                             1817
13941 continue                                                             1817
      if(isd .le. 0)goto 14011                                             1817
      xv=0.25                                                              1817
      goto 14021                                                           1818
14011 continue                                                             1819
14030 do 14031 j=1,ni                                                      1819
      if(ju(j).eq.0)goto 14031                                             1819
      jb=ix(j)                                                             1819
      je=ix(j+1)-1                                                         1820
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        1821
14031 continue                                                             1822
14032 continue                                                             1822
14021 continue                                                             1823
14001 continue                                                             1823
      if(flmin .ge. 1.0)goto 14051                                         1823
      eqs=max(eps,flmin)                                                   1823
      alf=eqs**(1.0/(nlam-1))                                              1823
14051 continue                                                             1824
      m=0                                                                  1824
      mm=0                                                                 1824
      nin=0                                                                1824
      nlp=0                                                                1824
      mnl=min(mnlam,nlam)                                                  1824
      bs=0.0                                                               1824
      svr=0.0                                                              1824
      o=0.0                                                                1825
      shr=shri*dev0                                                        1826
14060 do 14061 ilm=1,nlam                                                  1827
      if(flmin .lt. 1.0)goto 14081                                         1827
      al=ulam(ilm)                                                         1827
      goto 14071                                                           1828
14081 if(ilm .le. 2)goto 14091                                             1828
      al=al*alf                                                            1828
      goto 14071                                                           1829
14091 if(ilm .ne. 1)goto 14101                                             1829
      al=big                                                               1829
      goto 14111                                                           1830
14101 continue                                                             1830
      al=0.0                                                               1831
14120 do 14121 ic=1,nc                                                     1831
      v=q(:,ic)/sxp                                                        1831
      r=w*(y(:,ic)-v)                                                      1831
      v=w*v*(1.0-v)                                                        1832
14130 do 14131 j=1,ni                                                      1832
      if(ju(j).eq.0)goto 14131                                             1832
      if(vp(j).le.0.0)goto 14131                                           1833
      jb=ix(j)                                                             1833
      je=ix(j+1)-1                                                         1833
      jn=ix(j+1)-ix(j)                                                     1834
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1835
      gj=dot_product(sc(1:jn),x(jb:je))                                    1836
      gj=(gj-svr*xb(j))/xs(j)                                              1837
      al=max(al,abs(gj)/vp(j))                                             1838
14131 continue                                                             1839
14132 continue                                                             1839
14121 continue                                                             1840
14122 continue                                                             1840
      al=alf*al/bta                                                        1841
14111 continue                                                             1842
14071 continue                                                             1842
      al2=al*omb                                                           1842
      al1=al*bta                                                           1842
      nit=0                                                                1843
14140 continue                                                             1843
14141 continue                                                             1843
      ixx=0                                                                1843
      jxx=ixx                                                              1843
      ig=0                                                                 1844
14150 do 14151 ic=1,nc                                                     1844
      bs(0,ic)=b(0,ic)                                                     1845
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1846
      xm(0)=0.0                                                            1846
      svr=0.0                                                              1846
      o=0.0                                                                1847
14160 do 14161 i=1,no                                                      1847
      pic=q(i,ic)/sxp(i)                                                   1848
      if(pic .ge. pfm)goto 14181                                           1848
      pic=0.0                                                              1848
      v(i)=0.0                                                             1848
      goto 14171                                                           1849
14181 if(pic .le. pfx)goto 14191                                           1849
      pic=1.0                                                              1849
      v(i)=0.0                                                             1849
      goto 14201                                                           1850
14191 continue                                                             1850
      v(i)=w(i)*pic*(1.0-pic)                                              1850
      xm(0)=xm(0)+v(i)                                                     1850
14201 continue                                                             1851
14171 continue                                                             1851
      r(i)=w(i)*(y(i,ic)-pic)                                              1851
      svr=svr+r(i)                                                         1852
14161 continue                                                             1853
14162 continue                                                             1853
      if(xm(0).le.vmin)goto 14151                                          1853
      ig=1                                                                 1854
14210 do 14211 l=1,nin                                                     1854
      j=m(l)                                                               1855
      jb=ix(j)                                                             1855
      je=ix(j+1)-1                                                         1856
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             1857
      if(kopt .ne. 0)goto 14231                                            1858
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       1859
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          1860
14231 continue                                                             1861
14211 continue                                                             1862
14212 continue                                                             1862
14240 continue                                                             1862
14241 continue                                                             1862
      nlp=nlp+1                                                            1862
      dlx=0.0                                                              1863
14250 do 14251 k=1,ni                                                      1863
      if(ju(k).eq.0)goto 14251                                             1864
      jb=ix(k)                                                             1864
      je=ix(k+1)-1                                                         1864
      jn=ix(k+1)-ix(k)                                                     1864
      bk=b(k,ic)                                                           1865
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 1866
      gk=dot_product(sc(1:jn),x(jb:je))                                    1867
      gk=(gk-svr*xb(k))/xs(k)                                              1868
      u=gk+xv(k,ic)*b(k,ic)                                                1868
      au=abs(u)-vp(k)*al1                                                  1869
      if(au .gt. 0.0)goto 14271                                            1869
      b(k,ic)=0.0                                                          1869
      goto 14281                                                           1870
14271 continue                                                             1870
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1870
14281 continue                                                             1871
14261 continue                                                             1871
      d=b(k,ic)-bk                                                         1871
      if(abs(d).le.0.0)goto 14251                                          1872
      dlx=max(dlx,xv(k,ic)*d**2)                                           1873
      if(mm(k) .ne. 0)goto 14301                                           1873
      nin=nin+1                                                            1874
      if(nin .le. nx)goto 14321                                            1874
      jxx=1                                                                1874
      goto 14252                                                           1874
14321 continue                                                             1875
      mm(k)=nin                                                            1875
      m(nin)=k                                                             1876
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             1877
14301 continue                                                             1878
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1879
      o=o+d*(xb(k)/xs(k))                                                  1880
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1881
14251 continue                                                             1882
14252 continue                                                             1882
      if(jxx.gt.0)goto 14242                                               1883
      d=svr/xm(0)                                                          1884
      if(d .eq. 0.0)goto 14341                                             1884
      b(0,ic)=b(0,ic)+d                                                    1884
      dlx=max(dlx,xm(0)*d**2)                                              1885
      r=r-d*v                                                              1885
      svr=svr-d*xm(0)                                                      1886
14341 continue                                                             1887
      if(dlx.lt.shr)goto 14242                                             1888
14350 continue                                                             1888
14351 continue                                                             1888
      nlp=nlp+1                                                            1888
      dlx=0.0                                                              1889
14360 do 14361 l=1,nin                                                     1889
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
      if(au .gt. 0.0)goto 14381                                            1895
      b(k,ic)=0.0                                                          1895
      goto 14391                                                           1896
14381 continue                                                             1896
      b(k,ic)=sign(au,u)/(xv(k,ic)+vp(k)*al2)                              1896
14391 continue                                                             1897
14371 continue                                                             1897
      d=b(k,ic)-bk                                                         1897
      if(abs(d).le.0.0)goto 14361                                          1898
      dlx=max(dlx,xv(k,ic)*d**2)                                           1899
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              1900
      o=o+d*(xb(k)/xs(k))                                                  1901
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  1902
14361 continue                                                             1903
14362 continue                                                             1903
      d=svr/xm(0)                                                          1904
      if(d .eq. 0.0)goto 14411                                             1904
      b(0,ic)=b(0,ic)+d                                                    1904
      dlx=max(dlx,xm(0)*d**2)                                              1905
      r=r-d*v                                                              1905
      svr=svr-d*xm(0)                                                      1906
14411 continue                                                             1907
      if(dlx.lt.shr)goto 14352                                             1907
      goto 14351                                                           1908
14352 continue                                                             1908
      goto 14241                                                           1909
14242 continue                                                             1909
      if(jxx.gt.0)goto 14152                                               1910
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                         1911
      if(ixx .ne. 0)goto 14431                                             1912
14440 do 14441 j=1,nin                                                     1913
      if(xv(m(j),ic)*(b(m(j),ic)-bs(m(j),ic))**2 .le. shr)goto 14461       1913
      ixx=1                                                                1913
      goto 14442                                                           1913
14461 continue                                                             1914
14441 continue                                                             1915
14442 continue                                                             1915
14431 continue                                                             1916
      sc=b(0,ic)+g(:,ic)                                                   1916
      b0=0.0                                                               1917
14470 do 14471 j=1,nin                                                     1917
      l=m(j)                                                               1917
      jb=ix(l)                                                             1917
      je=ix(l+1)-1                                                         1918
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   1919
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            1920
14471 continue                                                             1921
14472 continue                                                             1921
      sc=min(max(exmn,sc+b0),exmx)                                         1922
      sxp=sxp-q(:,ic)                                                      1923
      q(:,ic)=min(max(emin*sxp,exp(dble(sc))),emax*sxp)                    1924
      sxp=sxp+q(:,ic)                                                      1925
14151 continue                                                             1926
14152 continue                                                             1926
      s=-sum(b(0,:))/nc                                                    1926
      b(0,:)=b(0,:)+s                                                      1926
      sc=s                                                                 1926
      b0=0.0                                                               1927
14480 do 14481 j=1,nin                                                     1927
      l=m(j)                                                               1928
      if(vp(l) .gt. 0.0)goto 14501                                         1928
      s=sum(b(l,:))/nc                                                     1928
      goto 14511                                                           1929
14501 continue                                                             1929
      s=elc(parm,nc,b(l,:),is)                                             1929
14511 continue                                                             1930
14491 continue                                                             1930
      b(l,:)=b(l,:)-s                                                      1931
      jb=ix(l)                                                             1931
      je=ix(l+1)-1                                                         1932
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         1933
      b0=b0+s*xb(l)/xs(l)                                                  1934
14481 continue                                                             1935
14482 continue                                                             1935
      sc=sc+b0                                                             1935
      sc=exp(sc)                                                           1935
      sxp=sxp*sc                                                           1935
14520 do 14521 ic=1,nc                                                     1935
      q(:,ic)=q(:,ic)*sc                                                   1935
14521 continue                                                             1936
14522 continue                                                             1936
      if(jxx.gt.0)goto 14142                                               1936
      if(ixx.eq.0)goto 14142                                               1936
      if(ig.eq.0)goto 14142                                                1937
      nit=nit+1                                                            1937
      if(nit .le. maxit)goto 14541                                         1937
      jerr=-ilm                                                            1937
      return                                                               1937
14541 continue                                                             1938
      goto 14141                                                           1939
14142 continue                                                             1939
      if(jxx .le. 0)goto 14561                                             1939
      jerr=-10000-ilm                                                      1939
      goto 14062                                                           1939
14561 continue                                                             1939
      devi=0.0                                                             1940
14570 do 14571 ic=1,nc                                                     1941
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1941
      a0(ic,ilm)=b(0,ic)                                                   1942
14580 do 14581 i=1,no                                                      1942
      if(y(i,ic).le.0.0)goto 14581                                         1943
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1944
14581 continue                                                             1945
14582 continue                                                             1945
14571 continue                                                             1946
14572 continue                                                             1946
      kin(ilm)=nin                                                         1946
      alm(ilm)=al                                                          1946
      lmu=ilm                                                              1947
      dev(ilm)=(dev1-devi)/dev0                                            1947
      if(ig.eq.0)goto 14062                                                1948
      if(ilm.lt.mnl)goto 14061                                             1948
      if(flmin.ge.1.0)goto 14061                                           1949
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 14062             1950
      if(dev(ilm).gt.devmax)goto 14062                                     1950
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 14062                             1951
14061 continue                                                             1952
14062 continue                                                             1952
      g=log(q)                                                             1952
14590 do 14591 i=1,no                                                      1952
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1952
14591 continue                                                             1953
14592 continue                                                             1953
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc)                            1954
      return                                                               1955
      end                                                                  1956
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  1957
      real a0(nc),ca(nx,nc),x(*),f(nc,n)                                   1957
      integer ia(*),ix(*),jx(*)                                            1958
14600 do 14601 ic=1,nc                                                     1958
      f(ic,:)=a0(ic)                                                       1958
14601 continue                                                             1959
14602 continue                                                             1959
14610 do 14611 j=1,nin                                                     1959
      k=ia(j)                                                              1959
      kb=ix(k)                                                             1959
      ke=ix(k+1)-1                                                         1960
14620 do 14621 ic=1,nc                                                     1960
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    1960
14621 continue                                                             1961
14622 continue                                                             1961
14611 continue                                                             1962
14612 continue                                                             1962
      return                                                               1963
      end                                                                  1964
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,ne,nx,nlam,flmin,ula   1966 
     *m,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam)              1967
      real ca(nx,nlam),dev(nlam),alm(nlam)                                 1968
      integer jd(*),ia(nx),nin(nlam)                                       1969
      real, dimension (:), allocatable :: xs,ww,vq                              
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 14641                                    1973
      jerr=10000                                                           1973
      return                                                               1973
14641 continue                                                             1974
      allocate(ww(1:no),stat=jerr)                                         1975
      allocate(ju(1:ni),stat=ierr)                                         1975
      jerr=jerr+ierr                                                       1976
      allocate(vq(1:ni),stat=ierr)                                         1976
      jerr=jerr+ierr                                                       1977
      if(isd .le. 0)goto 14661                                             1977
      allocate(xs(1:ni),stat=ierr)                                         1977
      jerr=jerr+ierr                                                       1977
14661 continue                                                             1978
      if(jerr.ne.0) return                                                 1979
      call chkvars(no,ni,x,ju)                                             1980
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1981
      if(maxval(ju) .gt. 0)goto 14681                                      1981
      jerr=7777                                                            1981
      return                                                               1981
14681 continue                                                             1982
      vq=max(0.0,vp)                                                       1982
      vq=vq*ni/sum(vq)                                                     1983
      ww=max(0.0,w)                                                        1983
      sw=sum(ww)                                                           1984
      if(sw .gt. 0.0)goto 14701                                            1984
      jerr=9999                                                            1984
      return                                                               1984
14701 continue                                                             1984
      ww=ww/sw                                                             1985
      call cstandard(no,ni,x,ww,ju,isd,xs)                                 1986
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr   1988 
     *,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1988
      dev0=2.0*sw*dev0                                                     1989
      if(isd .le. 0)goto 14721                                             1989
14730 do 14731 k=1,lmu                                                     1989
      nk=nin(k)                                                            1989
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   1989
14731 continue                                                             1989
14732 continue                                                             1989
14721 continue                                                             1990
      deallocate(ww,ju,vq)                                                 1990
      if(isd.gt.0) deallocate(xs)                                          1991
      return                                                               1992
      end                                                                  1993
      subroutine cstandard (no,ni,x,w,ju,isd,xs)                           1994
      real x(no,ni),w(no),xs(ni)                                           1994
      integer ju(ni)                                                       1995
14740 do 14741 j=1,ni                                                      1995
      if(ju(j).eq.0)goto 14741                                             1996
      xm=dot_product(w,x(:,j))                                             1996
      x(:,j)=x(:,j)-xm                                                     1997
      if(isd .le. 0)goto 14761                                             1997
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 1997
      x(:,j)=x(:,j)/xs(j)                                                  1997
14761 continue                                                             1998
14741 continue                                                             1999
14742 continue                                                             1999
      return                                                               2000
      end                                                                  2001
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,ne,nx,nlam,flmin,ula   2003 
     *m,cthri,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=0.001, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99)   2004
      real x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam)              2005
      real ao(nx,nlam),dev(nlam),alm(nlam)                                 2006
      integer ju(ni),m(nx),kin(nlam)                                       2007
      real, dimension (:), allocatable :: w,dk,v,xs,wr,a,as,f,dq                
      real, dimension (:), allocatable :: e,uu                                  
      integer, dimension (:), allocatable :: jp,kp,mm                           
      allocate(e(1:no),stat=jerr)                                          2013
      allocate(uu(1:no),stat=ierr)                                         2013
      jerr=jerr+ierr                                                       2014
      allocate(f(1:no),stat=ierr)                                          2014
      jerr=jerr+ierr                                                       2015
      allocate(w(1:no),stat=ierr)                                          2015
      jerr=jerr+ierr                                                       2016
      allocate(v(1:ni),stat=ierr)                                          2016
      jerr=jerr+ierr                                                       2017
      allocate(a(1:ni),stat=ierr)                                          2017
      jerr=jerr+ierr                                                       2018
      allocate(as(1:ni),stat=ierr)                                         2018
      jerr=jerr+ierr                                                       2019
      allocate(xs(1:ni),stat=ierr)                                         2019
      jerr=jerr+ierr                                                       2020
      allocate(jp(1:no),stat=ierr)                                         2020
      jerr=jerr+ierr                                                       2021
      allocate(kp(1:no),stat=ierr)                                         2021
      jerr=jerr+ierr                                                       2022
      allocate(dk(1:no),stat=ierr)                                         2022
      jerr=jerr+ierr                                                       2023
      allocate(wr(1:no),stat=ierr)                                         2023
      jerr=jerr+ierr                                                       2024
      allocate(dq(1:no),stat=ierr)                                         2024
      jerr=jerr+ierr                                                       2025
      allocate(mm(1:ni),stat=ierr)                                         2025
      jerr=jerr+ierr                                                       2026
      if(jerr.ne.0)go to 11400                                             2027
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2028
      if(jerr.ne.0) go to 11400                                            2028
      alpha=max(parm,1.0e-3)                                               2029
      oma=1.0-alpha                                                        2029
      nlm=0                                                                2030
      dq=d*q                                                               2030
      call died(no,nk,dq,kp,jp,dk)                                         2031
      a=0.0                                                                2031
      f=0.0                                                                2031
      e=q                                                                  2031
      fmax=log(huge(f(1))*0.1)                                             2032
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 2033
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         2033
      dev0=rr                                                              2034
14770 do 14771 i=1,no                                                      2034
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 14791                   2034
      w(i)=0.0                                                             2034
      wr(i)=w(i)                                                           2034
14791 continue                                                             2034
14771 continue                                                             2035
14772 continue                                                             2035
      if(nonzero(no,g) .eq. 0)goto 14811                                   2035
      f=g-dot_product(q,g)                                                 2036
      e=q*exp(sign(min(abs(f),fmax),f))                                    2037
14811 continue                                                             2038
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2039
      if(jerr.ne.0) go to 11400                                            2040
      call vars(no,ni,x,w,v)                                               2041
      if(flmin .ge. 1.0)goto 14831                                         2041
      eqs=max(eps,flmin)                                                   2041
      alf=eqs**(1.0/(nlam-1))                                              2041
14831 continue                                                             2042
      m=0                                                                  2042
      mm=0                                                                 2042
      nlp=0                                                                2042
      nin=nlp                                                              2042
      mnl=min(mnlam,nlam)                                                  2042
      as=0.0                                                               2042
      cthr=cthri*dev0                                                      2043
14840 do 14841 ilm=1,nlam                                                  2044
      if(flmin .lt. 1.0)goto 14861                                         2044
      al=ulam(ilm)                                                         2044
      goto 14851                                                           2045
14861 if(ilm .le. 2)goto 14871                                             2045
      al=al*alf                                                            2045
      goto 14851                                                           2046
14871 if(ilm .ne. 1)goto 14881                                             2046
      al=big                                                               2046
      goto 14891                                                           2047
14881 continue                                                             2047
      al=0.0                                                               2048
14900 do 14901 j=1,ni                                                      2048
      if(ju(j).eq.0)goto 14901                                             2048
      if(vp(j).le.0.0)goto 14901                                           2049
      al=max(al,abs(dot_product(wr,x(:,j)))/vp(j))                         2050
14901 continue                                                             2051
14902 continue                                                             2051
      al=alf*al/alpha                                                      2052
14891 continue                                                             2053
14851 continue                                                             2053
      sa=alpha*al                                                          2053
      omal=oma*al                                                          2053
      nit=0                                                                2054
14910 do 14911 ito=1,maxit                                                 2054
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2055
14920 do 14921 iti=1,maxit                                                 2055
      nlp=nlp+1                                                            2055
      dli=0.0                                                              2056
14930 do 14931 j=1,ni                                                      2056
      if(ju(j).eq.0)goto 14931                                             2057
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2058
      if(abs(u) .gt. vp(j)*sa)goto 14951                                   2058
      at=0.0                                                               2058
      goto 14961                                                           2059
14951 continue                                                             2059
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2059
14961 continue                                                             2060
14941 continue                                                             2060
      if(at .eq. a(j))goto 14981                                           2060
      del=at-a(j)                                                          2060
      a(j)=at                                                              2060
      dli=max(dli,v(j)*del**2)                                             2061
      wr=wr-del*w*x(:,j)                                                   2061
      f=f+del*x(:,j)                                                       2062
      if(mm(j) .ne. 0)goto 15001                                           2062
      nin=nin+1                                                            2062
      if(nin.gt.nx)goto 14932                                              2063
      mm(j)=nin                                                            2063
      m(nin)=j                                                             2064
15001 continue                                                             2065
14981 continue                                                             2066
14931 continue                                                             2067
14932 continue                                                             2067
      if(nin.gt.nx)goto 14922                                              2067
      if(dli.lt.cthr)goto 14922                                            2068
15010 do 15011 ita=1,maxit                                                 2068
      nlp=nlp+1                                                            2068
      dli=0.0                                                              2069
15020 do 15021 l=1,nin                                                     2069
      j=m(l)                                                               2070
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2071
      if(abs(u) .gt. vp(j)*sa)goto 15041                                   2071
      at=0.0                                                               2071
      goto 15051                                                           2072
15041 continue                                                             2072
      at=sign(abs(u)-vp(j)*sa,u)/(v(j)+vp(j)*omal)                         2072
15051 continue                                                             2073
15031 continue                                                             2073
      if(at .eq. a(j))goto 15071                                           2073
      del=at-a(j)                                                          2073
      a(j)=at                                                              2073
      dli=max(dli,v(j)*del**2)                                             2074
      wr=wr-del*w*x(:,j)                                                   2074
      f=f+del*x(:,j)                                                       2075
15071 continue                                                             2076
15021 continue                                                             2077
15022 continue                                                             2077
      if(dli.lt.cthr)goto 15012                                            2078
15011 continue                                                             2079
15012 continue                                                             2079
      if(dli .lt. cthr)goto 15091                                          2079
      jerr=-1                                                              2079
      go to 11400                                                          2079
15091 continue                                                             2080
14921 continue                                                             2081
14922 continue                                                             2081
      if(nin.gt.nx)goto 14912                                              2081
      if(dli .lt. cthr)goto 15111                                          2081
      jerr=-1                                                              2081
      go to 11400                                                          2081
15111 continue                                                             2082
      e=q*exp(sign(min(abs(f),fmax),f))                                    2083
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2084
      if(jerr.ne.0) go to 11400                                            2085
      call vars(no,ni,x,w,v)                                               2086
      ix=0                                                                 2087
15120 do 15121 j=1,nin                                                     2087
      if(v(m(j))*(a(m(j))-as(m(j)))**2.lt.cthr)goto 15121                  2087
      ix=1                                                                 2087
      goto 15122                                                           2087
15121 continue                                                             2088
15122 continue                                                             2088
      if(ix.eq.0)goto 14912                                                2089
14911 continue                                                             2090
14912 continue                                                             2090
      if(ix .eq. 0)goto 15141                                              2090
      jerr=-2                                                              2090
      goto 14842                                                           2090
15141 continue                                                             2091
      if(nin .le. nx)goto 15161                                            2091
      jerr=-10000-ilm                                                      2091
      goto 14842                                                           2091
15161 continue                                                             2092
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               2092
      kin(ilm)=nin                                                         2093
      alm(ilm)=al                                                          2093
      lmu=ilm                                                              2094
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   2095
      if(ilm.lt.mnl)goto 14841                                             2095
      if(flmin.ge.1.0)goto 14841                                           2096
      me=0                                                                 2096
15170 do 15171 j=1,nin                                                     2096
      if(ao(j,ilm).ne.0.0) me=me+1                                         2096
15171 continue                                                             2096
15172 continue                                                             2096
      if(me.gt.ne)goto 14842                                               2097
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 14842              2098
      if(dev(ilm).gt.devmax)goto 14842                                     2099
14841 continue                                                             2100
14842 continue                                                             2100
      g=f                                                                  2101
11400 continue                                                             2101
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm)                     2102
      return                                                               2103
      end                                                                  2104
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 2105
      real ca(nin),x(n,*),f(n)                                             2105
      integer ia(nin)                                                      2106
      f=0.0                                                                2106
      if(nin.le.0) return                                                  2107
15180 do 15181 i=1,n                                                       2107
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      2107
15181 continue                                                             2108
15182 continue                                                             2108
      return                                                               2109
      end                                                                  2110
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         2111
      real y(no),d(no),q(no)                                               2111
      integer jp(no),kp(*)                                                 2112
15190 do 15191 j=1,no                                                      2112
      jp(j)=j                                                              2112
15191 continue                                                             2112
15192 continue                                                             2112
      call psort7(y,jp,1,no)                                               2113
      nj=0                                                                 2113
15200 do 15201 j=1,no                                                      2113
      if(q(jp(j)).le.0.0)goto 15201                                        2113
      nj=nj+1                                                              2113
      jp(nj)=jp(j)                                                         2113
15201 continue                                                             2114
15202 continue                                                             2114
      if(nj .ne. 0)goto 15221                                              2114
      jerr=-4                                                              2114
      return                                                               2114
15221 continue                                                             2115
      j=1                                                                  2115
15230 continue                                                             2115
15231 if(d(jp(j)).gt.0.0)goto 15232                                        2115
      j=j+1                                                                2115
      if(j.gt.nj)goto 15232                                                2115
      goto 15231                                                           2116
15232 continue                                                             2116
      if(j .lt. nj-1)goto 15251                                            2116
      jerr=-5                                                              2116
      return                                                               2116
15251 continue                                                             2117
      j0=j-1                                                               2117
      nj=nj-j0                                                             2117
15260 do 15261 j=1,nj                                                      2117
      jp(j)=jp(j+j0)                                                       2117
15261 continue                                                             2118
15262 continue                                                             2118
      jerr=0                                                               2118
      nk=0                                                                 2118
      t0=y(jp(1))                                                          2118
      yk=t0                                                                2118
      j=2                                                                  2119
15270 continue                                                             2119
15271 continue                                                             2119
15280 continue                                                             2120
15281 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 15282                     2120
      j=j+1                                                                2120
      if(j.gt.nj)goto 15282                                                2120
      goto 15281                                                           2121
15282 continue                                                             2121
      nk=nk+1                                                              2121
      kp(nk)=j-1                                                           2121
      if(j.gt.nj)goto 15272                                                2122
      if(j .ne. nj)goto 15301                                              2122
      nk=nk+1                                                              2122
      kp(nk)=nj                                                            2122
      goto 15272                                                           2122
15301 continue                                                             2123
      yk=y(jp(j))                                                          2123
      j=j+1                                                                2124
      goto 15271                                                           2125
15272 continue                                                             2125
      return                                                               2126
      end                                                                  2127
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     2128
      real d(no),dk(nk),wr(no),w(no)                                       2129
      real e(no),u(no),b,c                                                 2129
      integer kp(nk),jp(no)                                                2130
      call usk(no,nk,kp,jp,e,u)                                            2131
      b=dk(1)/u(1)                                                         2131
      c=dk(1)/u(1)**2                                                      2131
      jerr=0                                                               2132
15310 do 15311 j=1,kp(1)                                                   2132
      i=jp(j)                                                              2133
      w(i)=e(i)*(b-e(i)*c)                                                 2133
      if(w(i) .gt. 0.0)goto 15331                                          2133
      jerr=-3                                                              2133
      return                                                               2133
15331 continue                                                             2134
      wr(i)=d(i)-e(i)*b                                                    2135
15311 continue                                                             2136
15312 continue                                                             2136
15340 do 15341 k=2,nk                                                      2136
      j1=kp(k-1)+1                                                         2136
      j2=kp(k)                                                             2137
      b=b+dk(k)/u(k)                                                       2137
      c=c+dk(k)/u(k)**2                                                    2138
15350 do 15351 j=j1,j2                                                     2138
      i=jp(j)                                                              2139
      w(i)=e(i)*(b-e(i)*c)                                                 2139
      if(w(i) .gt. 0.0)goto 15371                                          2139
      jerr=-3                                                              2139
      return                                                               2139
15371 continue                                                             2140
      wr(i)=d(i)-e(i)*b                                                    2141
15351 continue                                                             2142
15352 continue                                                             2142
15341 continue                                                             2143
15342 continue                                                             2143
      return                                                               2144
      end                                                                  2145
      subroutine vars(no,ni,x,w,v)                                         2146
      real x(no,ni),w(no),v(ni)                                            2147
15380 do 15381 j=1,ni                                                      2147
      v(j)=dot_product(w,x(:,j)**2)                                        2147
15381 continue                                                             2148
15382 continue                                                             2148
      return                                                               2149
      end                                                                  2150
      subroutine died(no,nk,d,kp,jp,dk)                                    2151
      real d(no),dk(nk)                                                    2151
      integer kp(nk),jp(no)                                                2152
      dk(1)=sum(d(jp(1:kp(1))))                                            2153
15390 do 15391 k=2,nk                                                      2153
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  2153
15391 continue                                                             2154
15392 continue                                                             2154
      return                                                               2155
      end                                                                  2156
      subroutine usk(no,nk,kp,jp,e,u)                                      2157
      real e(no),u(nk),h                                                   2157
      integer kp(nk),jp(no)                                                2158
      h=0.0                                                                2159
15400 do 15401 k=nk,1,-1                                                   2159
      j2=kp(k)                                                             2160
      j1=1                                                                 2160
      if(k.gt.1) j1=kp(k-1)+1                                              2161
15410 do 15411 j=j2,j1,-1                                                  2161
      h=h+e(jp(j))                                                         2161
15411 continue                                                             2162
15412 continue                                                             2162
      u(k)=h                                                               2163
15401 continue                                                             2164
15402 continue                                                             2164
      return                                                               2165
      end                                                                  2166
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             2167
      real d(no),dk(nk),f(no)                                              2168
      integer kp(nk),jp(no)                                                2168
      real e(no),u(nk),s                                                   2169
      call usk(no,nk,kp,jp,e,u)                                            2169
      u=log(u)                                                             2170
      risk=dot_product(d,f)-dot_product(dk,u)                              2171
      return                                                               2172
      end                                                                  2173
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 2174
      real x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(nlam)          2175
      real, dimension (:), allocatable :: dk,f,xm,dq,q                          
      real, dimension (:), allocatable :: e,uu                                  
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          2181
      allocate(q(1:no),stat=ierr)                                          2181
      jerr=jerr+ierr                                                       2182
      allocate(uu(1:no),stat=ierr)                                         2182
      jerr=jerr+ierr                                                       2183
      allocate(f(1:no),stat=ierr)                                          2183
      jerr=jerr+ierr                                                       2184
      allocate(dk(1:no),stat=ierr)                                         2184
      jerr=jerr+ierr                                                       2185
      allocate(jp(1:no),stat=ierr)                                         2185
      jerr=jerr+ierr                                                       2186
      allocate(kp(1:no),stat=ierr)                                         2186
      jerr=jerr+ierr                                                       2187
      allocate(dq(1:no),stat=ierr)                                         2187
      jerr=jerr+ierr                                                       2188
      allocate(xm(1:ni),stat=ierr)                                         2188
      jerr=jerr+ierr                                                       2189
      if(jerr.ne.0) go to 11400                                            2190
      q=max(0.0,w)                                                         2190
      sw=sum(q)                                                            2191
      if(sw .gt. 0.0)goto 15431                                            2191
      jerr=9999                                                            2191
      go to 11400                                                          2191
15431 continue                                                             2192
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2193
      if(jerr.ne.0) go to 11400                                            2193
      fmax=log(huge(e(1))*0.1)                                             2194
      dq=d*q                                                               2194
      call died(no,nk,dq,kp,jp,dk)                                         2194
      gm=dot_product(q,g)/sw                                               2195
15440 do 15441 j=1,ni                                                      2195
      xm(j)=dot_product(q,x(:,j))/sw                                       2195
15441 continue                                                             2196
15442 continue                                                             2196
15450 do 15451 lam=1,nlam                                                  2197
15460 do 15461 i=1,no                                                      2197
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       2198
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        2199
15461 continue                                                             2200
15462 continue                                                             2200
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          2201
15451 continue                                                             2202
15452 continue                                                             2202
11400 continue                                                             2202
      deallocate(e,uu,dk,f,jp,kp,dq)                                       2203
      return                                                               2204
      end                                                                  2205
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,ne,nx,nlam,flmin,ulam   2207 
     *,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)                    2208
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2209
      integer jd(*),ia(nx),nin(nlam)                                       2210
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 15481                                    2214
      jerr=10000                                                           2214
      return                                                               2214
15481 continue                                                             2215
      if(minval(y) .ge. 0.0)goto 15501                                     2215
      jerr=8888                                                            2215
      return                                                               2215
15501 continue                                                             2216
      allocate(ww(1:no),stat=jerr)                                         2217
      allocate(ju(1:ni),stat=ierr)                                         2217
      jerr=jerr+ierr                                                       2218
      allocate(vq(1:ni),stat=ierr)                                         2218
      jerr=jerr+ierr                                                       2219
      allocate(xm(1:ni),stat=ierr)                                         2219
      jerr=jerr+ierr                                                       2220
      if(isd .le. 0)goto 15521                                             2220
      allocate(xs(1:ni),stat=ierr)                                         2220
      jerr=jerr+ierr                                                       2220
15521 continue                                                             2221
      if(jerr.ne.0) return                                                 2222
      call chkvars(no,ni,x,ju)                                             2223
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2224
      if(maxval(ju) .gt. 0)goto 15541                                      2224
      jerr=7777                                                            2224
      go to 11400                                                          2224
15541 continue                                                             2225
      vq=max(0.0,vp)                                                       2225
      vq=vq*ni/sum(vq)                                                     2226
      ww=max(0.0,w)                                                        2226
      sw=sum(ww)                                                           2226
      if(sw .gt. 0.0)goto 15561                                            2226
      jerr=9999                                                            2226
      go to 11400                                                          2226
15561 continue                                                             2227
      ww=ww/sw                                                             2228
      call lstandard1(no,ni,x,ww,ju,isd,xm,xs)                             2229
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,ne,nx,nlam,flmin,ulam,thr,   2231 
     *  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11400                                            2231
      dev0=2.0*sw*dev0                                                     2232
15570 do 15571 k=1,lmu                                                     2232
      nk=nin(k)                                                            2233
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2234
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2235
15571 continue                                                             2236
15572 continue                                                             2236
11400 continue                                                             2236
      deallocate(ww,ju,vq,xm)                                              2236
      if(isd.gt.0) deallocate(xs)                                          2237
      return                                                               2238
      end                                                                  2239
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,ne,nx,nlam,flmin,ulam   2241 
     *,shri,  isd,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2242 
     *9)
      real x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)                    2243
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2244
      integer ju(ni),m(nx),kin(nlam)                                       2245
      real, dimension (:), allocatable :: t,w,wr,v,a,f,as                       
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                          2250
      allocate(as(1:ni),stat=ierr)                                         2250
      jerr=jerr+ierr                                                       2251
      allocate(t(1:no),stat=ierr)                                          2251
      jerr=jerr+ierr                                                       2252
      allocate(mm(1:ni),stat=ierr)                                         2252
      jerr=jerr+ierr                                                       2253
      allocate(wr(1:no),stat=ierr)                                         2253
      jerr=jerr+ierr                                                       2254
      allocate(v(1:ni),stat=ierr)                                          2254
      jerr=jerr+ierr                                                       2255
      allocate(w(1:no),stat=ierr)                                          2255
      jerr=jerr+ierr                                                       2256
      allocate(f(1:no),stat=ierr)                                          2256
      jerr=jerr+ierr                                                       2257
      if(jerr.ne.0) return                                                 2258
      bta=max(parm,1.0e-3)                                                 2258
      omb=1.0-bta                                                          2259
      t=q*y                                                                2259
      yb=sum(t)                                                            2259
      fmax=log(huge(bta)*0.1)                                              2260
      if(nonzero(no,g) .ne. 0)goto 15591                                   2260
      w=q*yb                                                               2260
      az=log(yb)                                                           2260
      f=az                                                                 2260
      goto 15601                                                           2261
15591 continue                                                             2261
      w=q*exp(sign(min(abs(g),fmax),g))                                    2261
      v0=sum(w)                                                            2261
      eaz=yb/v0                                                            2262
      w=eaz*w                                                              2262
      az=log(eaz)                                                          2262
      f=az+g                                                               2263
15601 continue                                                             2264
15581 continue                                                             2264
      a=0.0                                                                2264
      wr=t-w                                                               2264
      v0=yb                                                                2264
      dv0=yb*(log(yb)-1.0)                                                 2264
      dvr=-yb                                                              2265
15610 do 15611 i=1,no                                                      2265
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               2265
15611 continue                                                             2265
15612 continue                                                             2265
      dvr=dvr-dv0                                                          2265
      dev0=dvr                                                             2266
15620 do 15621 j=1,ni                                                      2266
      if(ju(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                         2266
15621 continue                                                             2267
15622 continue                                                             2267
      if(flmin .ge. 1.0)goto 15641                                         2267
      eqs=max(eps,flmin)                                                   2267
      alf=eqs**(1.0/(nlam-1))                                              2267
15641 continue                                                             2268
      m=0                                                                  2268
      mm=0                                                                 2268
      nlp=0                                                                2268
      nin=nlp                                                              2268
      mnl=min(mnlam,nlam)                                                  2268
      shr=shri*dev0                                                        2269
15650 do 15651 ilm=1,nlam                                                  2270
      if(flmin .lt. 1.0)goto 15671                                         2270
      al=ulam(ilm)                                                         2270
      goto 15661                                                           2271
15671 if(ilm .le. 2)goto 15681                                             2271
      al=al*alf                                                            2271
      goto 15661                                                           2272
15681 if(ilm .ne. 1)goto 15691                                             2272
      al=big                                                               2272
      goto 15701                                                           2273
15691 continue                                                             2273
      al=0.0                                                               2274
15710 do 15711 j=1,ni                                                      2274
      if(ju(j).eq.0)goto 15711                                             2274
      if(vp(j).le.0.0)goto 15711                                           2275
      al=max(al,abs(dot_product(wr,x(:,j)))/vp(j))                         2276
15711 continue                                                             2277
15712 continue                                                             2277
      al=alf*al/bta                                                        2278
15701 continue                                                             2279
15661 continue                                                             2279
      al2=al*omb                                                           2279
      al1=al*bta                                                           2279
      nit=0                                                                2280
15720 continue                                                             2280
15721 continue                                                             2280
      nit=nit+1                                                            2280
      if(nit .le. maxit)goto 15741                                         2280
      jerr=-2                                                              2280
      go to 11400                                                          2280
15741 continue                                                             2281
      az0=az                                                               2281
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2281
      nit1=0                                                               2282
15750 continue                                                             2282
15751 continue                                                             2282
      nit1=nit1+1                                                          2282
      if(nit1 .le. maxit)goto 15771                                        2282
      jerr=-1                                                              2282
      go to 11400                                                          2282
15771 continue                                                             2283
      nlp=nlp+1                                                            2283
      dlx=0.0                                                              2284
15780 do 15781 k=1,ni                                                      2284
      if(ju(k).eq.0)goto 15781                                             2284
      ak=a(k)                                                              2285
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2285
      au=abs(u)-vp(k)*al1                                                  2286
      if(au .gt. 0.0)goto 15801                                            2286
      a(k)=0.0                                                             2286
      goto 15811                                                           2287
15801 continue                                                             2287
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2287
15811 continue                                                             2288
15791 continue                                                             2288
      if(a(k).eq.ak)goto 15781                                             2288
      d=a(k)-ak                                                            2288
      dlx=max(dlx,v(k)*d**2)                                               2289
      wr=wr-d*w*x(:,k)                                                     2289
      f=f+d*x(:,k)                                                         2290
      if(mm(k) .ne. 0)goto 15831                                           2290
      nin=nin+1                                                            2290
      if(nin.gt.nx)goto 15782                                              2291
      mm(k)=nin                                                            2291
      m(nin)=k                                                             2292
15831 continue                                                             2293
15781 continue                                                             2294
15782 continue                                                             2294
      if(nin.gt.nx)goto 15752                                              2294
      d=sum(wr)/v0                                                         2295
      az=az+d                                                              2295
      dlx=max(dlx,v0*d**2)                                                 2295
      wr=wr-d*w                                                            2295
      f=f+d                                                                2296
      if(dlx.lt.shr)goto 15752                                             2296
      nit2=0                                                               2297
15840 continue                                                             2297
15841 continue                                                             2297
      nit2=nit2+1                                                          2297
      if(nit2 .le. maxit)goto 15861                                        2297
      jerr=-1                                                              2297
      go to 11400                                                          2297
15861 continue                                                             2298
      nlp=nlp+1                                                            2298
      dlx=0.0                                                              2299
15870 do 15871 l=1,nin                                                     2299
      k=m(l)                                                               2299
      ak=a(k)                                                              2300
      u=dot_product(wr,x(:,k))+v(k)*ak                                     2300
      au=abs(u)-vp(k)*al1                                                  2301
      if(au .gt. 0.0)goto 15891                                            2301
      a(k)=0.0                                                             2301
      goto 15901                                                           2302
15891 continue                                                             2302
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2302
15901 continue                                                             2303
15881 continue                                                             2303
      if(a(k).eq.ak)goto 15871                                             2303
      d=a(k)-ak                                                            2303
      dlx=max(dlx,v(k)*d**2)                                               2304
      wr=wr-d*w*x(:,k)                                                     2304
      f=f+d*x(:,k)                                                         2306
15871 continue                                                             2306
15872 continue                                                             2306
      d=sum(wr)/v0                                                         2306
      az=az+d                                                              2306
      dlx=max(dlx,v0*d**2)                                                 2306
      wr=wr-d*w                                                            2306
      f=f+d                                                                2307
      if(dlx.lt.shr)goto 15842                                             2307
      goto 15841                                                           2308
15842 continue                                                             2308
      goto 15751                                                           2309
15752 continue                                                             2309
      if(nin.gt.nx)goto 15722                                              2310
      w=q*exp(sign(min(abs(f),fmax),f))                                    2310
      v0=sum(w)                                                            2311
      wr=t-w                                                               2311
15910 do 15911 j=1,ni                                                      2311
      if(ju(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                         2311
15911 continue                                                             2312
15912 continue                                                             2312
      if(v0*abs(az-az0) .ge. shr)goto 15931                                2312
      ix=0                                                                 2313
15940 do 15941 j=1,nin                                                     2313
      if(v(m(j))*abs(a(m(j))-as(m(j))).lt.shr)goto 15941                   2313
      ix=1                                                                 2313
      goto 15942                                                           2313
15941 continue                                                             2314
15942 continue                                                             2314
      if(ix.eq.0)goto 15722                                                2315
15931 continue                                                             2316
      goto 15721                                                           2317
15722 continue                                                             2317
      if(nin .le. nx)goto 15961                                            2317
      jerr=-10000-ilm                                                      2317
      goto 15652                                                           2317
15961 continue                                                             2318
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2318
      kin(ilm)=nin                                                         2319
      a0(ilm)=az                                                           2319
      alm(ilm)=al                                                          2319
      lmu=ilm                                                              2320
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               2321
      if(ilm.lt.mnl)goto 15651                                             2321
      if(flmin.ge.1.0)goto 15651                                           2322
      me=0                                                                 2322
15970 do 15971 j=1,nin                                                     2322
      if(ca(j,ilm).ne.0.0) me=me+1                                         2322
15971 continue                                                             2322
15972 continue                                                             2322
      if(me.gt.ne)goto 15652                                               2323
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 15652              2324
      if(dev(ilm).gt.devmax)goto 15652                                     2325
15651 continue                                                             2326
15652 continue                                                             2326
      g=f                                                                  2327
11400 continue                                                             2327
      deallocate(t,w,wr,v,a,f,as,mm)                                       2328
      return                                                               2329
      end                                                                  2330
      function nonzero(n,v)                                                2331
      real v(n)                                                            2332
      nonzero=0                                                            2332
15980 do 15981 i=1,n                                                       2332
      if(v(i) .eq. 0.0)goto 16001                                          2332
      nonzero=1                                                            2332
      return                                                               2332
16001 continue                                                             2332
15981 continue                                                             2333
15982 continue                                                             2333
      return                                                               2334
      end                                                                  2335
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               2336
      real a(nx,lmu),b(ni,lmu)                                             2336
      integer ia(nx),nin(lmu)                                              2337
16010 do 16011 lam=1,lmu                                                   2337
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        2337
16011 continue                                                             2338
16012 continue                                                             2338
      return                                                               2339
      end                                                                  2340
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           2341
      real a(nx,nc,lmu),b(ni,nc,lmu)                                       2341
      integer ia(nx),nin(lmu)                                              2342
16020 do 16021 lam=1,lmu                                                   2342
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             2342
16021 continue                                                             2343
16022 continue                                                             2343
      return                                                               2344
      end                                                                  2345
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               2346
      real x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)       2347
      real, dimension (:), allocatable :: w                                     
      if(minval(y) .ge. 0.0)goto 16041                                     2350
      jerr=8888                                                            2350
      return                                                               2350
16041 continue                                                             2351
      allocate(w(1:no),stat=jerr)                                          2351
      if(jerr.ne.0) return                                                 2352
      w=max(0.0,q)                                                         2352
      sw=sum(w)                                                            2352
      if(sw .gt. 0.0)goto 16061                                            2352
      jerr=9999                                                            2352
      go to 11400                                                          2352
16061 continue                                                             2353
      yb=dot_product(w,y)/sw                                               2353
      fmax=log(huge(y(1))*0.1)                                             2354
16070 do 16071 lam=1,nlam                                                  2354
      s=0.0                                                                2355
16080 do 16081 i=1,no                                                      2355
      if(w(i).le.0.0)goto 16081                                            2356
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          2357
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      2358
16081 continue                                                             2359
16082 continue                                                             2359
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2360
16071 continue                                                             2361
16072 continue                                                             2361
11400 continue                                                             2361
      deallocate(w)                                                        2362
      return                                                               2363
      end                                                                  2364
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,ne,nx,nlam,fl   2366 
     *min,ulam,thr,  isd,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      real x(*),y(no),g(no),w(no),vp(ni),ulam(nlam)                        2367
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2368
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           2369
      real, dimension (:), allocatable :: xm,xs,ww,vq                           
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16101                                    2373
      jerr=10000                                                           2373
      return                                                               2373
16101 continue                                                             2374
      if(minval(y) .ge. 0.0)goto 16121                                     2374
      jerr=8888                                                            2374
      return                                                               2374
16121 continue                                                             2375
      allocate(ww(1:no),stat=jerr)                                         2376
      allocate(ju(1:ni),stat=ierr)                                         2376
      jerr=jerr+ierr                                                       2377
      allocate(vq(1:ni),stat=ierr)                                         2377
      jerr=jerr+ierr                                                       2378
      allocate(xm(1:ni),stat=ierr)                                         2378
      jerr=jerr+ierr                                                       2379
      allocate(xs(1:ni),stat=ierr)                                         2379
      jerr=jerr+ierr                                                       2380
      if(jerr.ne.0) return                                                 2381
      call spchkvars(no,ni,x,ix,ju)                                        2382
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2383
      if(maxval(ju) .gt. 0)goto 16141                                      2383
      jerr=7777                                                            2383
      go to 11400                                                          2383
16141 continue                                                             2384
      vq=max(0.0,vp)                                                       2384
      vq=vq*ni/sum(vq)                                                     2385
      ww=max(0.0,w)                                                        2385
      sw=sum(ww)                                                           2385
      if(sw .gt. 0.0)goto 16161                                            2385
      jerr=9999                                                            2385
      go to 11400                                                          2385
16161 continue                                                             2386
      ww=ww/sw                                                             2387
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,xm,xs)                     2388
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,ne,nx,nlam,flmin,u   2390 
     *lam,thr,  isd,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 11400                                            2390
      dev0=2.0*sw*dev0                                                     2391
16170 do 16171 k=1,lmu                                                     2391
      nk=nin(k)                                                            2392
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2393
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2394
16171 continue                                                             2395
16172 continue                                                             2395
11400 continue                                                             2395
      deallocate(ww,ju,vq,xm,xs)                                           2396
      return                                                               2397
      end                                                                  2398
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,ne,nx,nlam,fl   2400 
     *min,ulam,  shri,isd,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,j
     *err)
      parameter(sml=1.0e-4, eps=1.0e-6, big=9.9e35, mnlam=5, devmax=0.99   2401 
     *9)
      real x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),xs(ni)          2402
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)                        2403
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2404
      real, dimension (:), allocatable :: qy,t,w,wr,v,a,as,xm                   
      integer, dimension (:), allocatable :: mm                                 
      allocate(a(1:ni),stat=jerr)                                          2409
      allocate(as(1:ni),stat=ierr)                                         2409
      jerr=jerr+ierr                                                       2410
      allocate(t(1:no),stat=ierr)                                          2410
      jerr=jerr+ierr                                                       2411
      allocate(mm(1:ni),stat=ierr)                                         2411
      jerr=jerr+ierr                                                       2412
      allocate(wr(1:no),stat=ierr)                                         2412
      jerr=jerr+ierr                                                       2413
      allocate(v(1:ni),stat=ierr)                                          2413
      jerr=jerr+ierr                                                       2414
      allocate(xm(1:ni),stat=ierr)                                         2414
      jerr=jerr+ierr                                                       2415
      allocate(w(1:no),stat=ierr)                                          2415
      jerr=jerr+ierr                                                       2416
      allocate(qy(1:no),stat=ierr)                                         2416
      jerr=jerr+ierr                                                       2417
      if(jerr.ne.0) return                                                 2418
      bta=max(parm,1.0e-3)                                                 2418
      omb=1.0-bta                                                          2418
      fmax=log(huge(bta)*0.1)                                              2419
      qy=q*y                                                               2419
      yb=sum(qy)                                                           2420
      if(nonzero(no,g) .ne. 0)goto 16191                                   2420
      w=q*yb                                                               2420
      az=log(yb)                                                           2420
      uu=az                                                                2421
      xm=yb*xb                                                             2421
      v=yb                                                                 2421
      t=0.0                                                                2422
      goto 16201                                                           2423
16191 continue                                                             2423
      w=q*exp(sign(min(abs(g),fmax),g))                                    2423
      ww=sum(w)                                                            2423
      eaz=yb/ww                                                            2424
      w=eaz*w                                                              2424
      az=log(eaz)                                                          2424
      uu=az                                                                2424
      t=g                                                                  2425
16210 do 16211 j=1,ni                                                      2425
      if(ju(j).eq.0)goto 16211                                             2425
      jb=ix(j)                                                             2425
      je=ix(j+1)-1                                                         2426
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2427
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   2429 
     *b(j)**2)/xs(j)**2
16211 continue                                                             2430
16212 continue                                                             2430
16201 continue                                                             2431
16181 continue                                                             2431
      tt=yb*uu                                                             2431
      ww=yb                                                                2431
      wr=qy-q*(yb*(1.0-uu))                                                2431
      a=0.0                                                                2432
      dv0=yb*(log(yb)-1.0)                                                 2432
      dvr=-yb                                                              2433
16220 do 16221 i=1,no                                                      2433
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             2433
16221 continue                                                             2433
16222 continue                                                             2433
      dvr=dvr-dv0                                                          2433
      dev0=dvr                                                             2434
      if(flmin .ge. 1.0)goto 16241                                         2434
      eqs=max(eps,flmin)                                                   2434
      alf=eqs**(1.0/(nlam-1))                                              2434
16241 continue                                                             2435
      m=0                                                                  2435
      mm=0                                                                 2435
      nlp=0                                                                2435
      nin=nlp                                                              2435
      mnl=min(mnlam,nlam)                                                  2435
      shr=shri*dev0                                                        2436
16250 do 16251 ilm=1,nlam                                                  2437
      if(flmin .lt. 1.0)goto 16271                                         2437
      al=ulam(ilm)                                                         2437
      goto 16261                                                           2438
16271 if(ilm .le. 2)goto 16281                                             2438
      al=al*alf                                                            2438
      goto 16261                                                           2439
16281 if(ilm .ne. 1)goto 16291                                             2439
      al=big                                                               2439
      goto 16301                                                           2440
16291 continue                                                             2440
      al=0.0                                                               2441
16310 do 16311 j=1,ni                                                      2441
      if(ju(j).eq.0)goto 16311                                             2441
      if(vp(j).le.0.0)goto 16311                                           2442
      jb=ix(j)                                                             2442
      je=ix(j+1)-1                                                         2443
      al=max(al,  abs(dot_product(qy(jx(jb:je)),x(jb:je))-xm(j))/(xs(j)*   2445 
     *vp(j)))
16311 continue                                                             2446
16312 continue                                                             2446
      al=alf*al/bta                                                        2447
16301 continue                                                             2448
16261 continue                                                             2448
      al2=al*omb                                                           2448
      al1=al*bta                                                           2448
      nit=0                                                                2449
16320 continue                                                             2449
16321 continue                                                             2449
      nit=nit+1                                                            2449
      if(nit .le. maxit)goto 16341                                         2449
      jerr=-2                                                              2449
      go to 11400                                                          2449
16341 continue                                                             2450
      az0=az                                                               2450
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2450
      nit1=0                                                               2451
16350 continue                                                             2451
16351 continue                                                             2451
      nit1=nit1+1                                                          2451
      if(nit1 .le. maxit)goto 16371                                        2451
      jerr=-1                                                              2451
      go to 11400                                                          2451
16371 continue                                                             2452
      nlp=nlp+1                                                            2452
      dlx=0.0                                                              2453
16380 do 16381 k=1,ni                                                      2453
      if(ju(k).eq.0)goto 16381                                             2453
      jb=ix(k)                                                             2453
      je=ix(k+1)-1                                                         2453
      ak=a(k)                                                              2454
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2456 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2457
      if(au .gt. 0.0)goto 16401                                            2457
      a(k)=0.0                                                             2457
      goto 16411                                                           2458
16401 continue                                                             2458
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2458
16411 continue                                                             2459
16391 continue                                                             2459
      if(a(k).eq.ak)goto 16381                                             2460
      if(mm(k) .ne. 0)goto 16431                                           2460
      nin=nin+1                                                            2460
      if(nin.gt.nx)goto 16382                                              2461
      mm(k)=nin                                                            2461
      m(nin)=k                                                             2462
16431 continue                                                             2463
      d=a(k)-ak                                                            2463
      dlx=max(dlx,v(k)*d**2)                                               2463
      dv=d/xs(k)                                                           2464
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2465
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2466
      uu=uu-dv*xb(k)                                                       2466
      tt=tt-dv*xm(k)                                                       2467
16381 continue                                                             2468
16382 continue                                                             2468
      if(nin.gt.nx)goto 16352                                              2468
      d=tt/ww-uu                                                           2469
      az=az+d                                                              2469
      dlx=max(dlx,ww*d**2)                                                 2469
      uu=uu+d                                                              2470
      if(dlx.lt.shr)goto 16352                                             2470
      nit2=0                                                               2471
16440 continue                                                             2471
16441 continue                                                             2471
      nit2=nit2+1                                                          2471
      if(nit2 .le. maxit)goto 16461                                        2471
      jerr=-1                                                              2471
      go to 11400                                                          2471
16461 continue                                                             2472
      nlp=nlp+1                                                            2472
      dlx=0.0                                                              2473
16470 do 16471 l=1,nin                                                     2473
      k=m(l)                                                               2473
      jb=ix(k)                                                             2473
      je=ix(k+1)-1                                                         2473
      ak=a(k)                                                              2474
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   2476 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  2477
      if(au .gt. 0.0)goto 16491                                            2477
      a(k)=0.0                                                             2477
      goto 16501                                                           2478
16491 continue                                                             2478
      a(k)=sign(au,u)/(v(k)+vp(k)*al2)                                     2478
16501 continue                                                             2479
16481 continue                                                             2479
      if(a(k).eq.ak)goto 16471                                             2479
      d=a(k)-ak                                                            2479
      dlx=max(dlx,v(k)*d**2)                                               2480
      dv=d/xs(k)                                                           2480
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 2481
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                2482
      uu=uu-dv*xb(k)                                                       2482
      tt=tt-dv*xm(k)                                                       2483
16471 continue                                                             2484
16472 continue                                                             2484
      d=tt/ww-uu                                                           2484
      az=az+d                                                              2484
      dlx=max(dlx,ww*d**2)                                                 2484
      uu=uu+d                                                              2485
      if(dlx.lt.shr)goto 16442                                             2485
      goto 16441                                                           2486
16442 continue                                                             2486
      goto 16351                                                           2487
16352 continue                                                             2487
      if(nin.gt.nx)goto 16322                                              2488
      euu=exp(sign(min(abs(uu),fmax),uu))                                  2489
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                2489
      ww=sum(w)                                                            2490
      wr=qy-w*(1.0-uu)                                                     2490
      tt=sum(wr)                                                           2491
16510 do 16511 j=1,ni                                                      2491
      if(ju(j).eq.0)goto 16511                                             2491
      jb=ix(j)                                                             2491
      je=ix(j+1)-1                                                         2492
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2493
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   2495 
     *b(j)**2)/xs(j)**2
16511 continue                                                             2496
16512 continue                                                             2496
      if(ww*(az-az0)**2 .ge. shr)goto 16531                                2496
      ixx=0                                                                2497
16540 do 16541 j=1,nin                                                     2497
      if(v(m(j))*(a(m(j))-as(m(j)))**2.lt.shr)goto 16541                   2497
      ixx=1                                                                2497
      goto 16542                                                           2497
16541 continue                                                             2498
16542 continue                                                             2498
      if(ixx.eq.0)goto 16322                                               2499
16531 continue                                                             2500
      goto 16321                                                           2501
16322 continue                                                             2501
      if(nin .le. nx)goto 16561                                            2501
      jerr=-10000-ilm                                                      2501
      goto 16252                                                           2501
16561 continue                                                             2502
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               2502
      kin(ilm)=nin                                                         2503
      a0(ilm)=az                                                           2503
      alm(ilm)=al                                                          2503
      lmu=ilm                                                              2504
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        2505
      if(ilm.lt.mnl)goto 16251                                             2505
      if(flmin.ge.1.0)goto 16251                                           2506
      me=0                                                                 2506
16570 do 16571 j=1,nin                                                     2506
      if(ca(j,ilm).ne.0.0) me=me+1                                         2506
16571 continue                                                             2506
16572 continue                                                             2506
      if(me.gt.ne)goto 16252                                               2507
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 16252              2508
      if(dev(ilm).gt.devmax)goto 16252                                     2509
16251 continue                                                             2510
16252 continue                                                             2510
      g=t+uu                                                               2511
11400 continue                                                             2511
      deallocate(t,w,wr,v,a,qy,xm,as,mm)                                   2512
      return                                                               2513
      end                                                                  2514
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       2515
      real x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(nlam)           2516
      integer ix(*),jx(*)                                                  2517
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 16591                                     2520
      jerr=8888                                                            2520
      return                                                               2520
16591 continue                                                             2521
      allocate(w(1:no),stat=jerr)                                          2522
      allocate(f(1:no),stat=ierr)                                          2522
      jerr=jerr+ierr                                                       2523
      if(jerr.ne.0) return                                                 2524
      w=max(0.0,q)                                                         2524
      sw=sum(w)                                                            2524
      if(sw .gt. 0.0)goto 16611                                            2524
      jerr=9999                                                            2524
      go to 11400                                                          2524
16611 continue                                                             2525
      yb=dot_product(w,y)/sw                                               2525
      fmax=log(huge(y(1))*0.1)                                             2526
16620 do 16621 lam=1,nlam                                                  2526
      f=a0(lam)                                                            2527
16630 do 16631 j=1,ni                                                      2527
      if(a(j,lam).eq.0.0)goto 16631                                        2527
      jb=ix(j)                                                             2527
      je=ix(j+1)-1                                                         2528
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          2529
16631 continue                                                             2530
16632 continue                                                             2530
      f=f+g                                                                2531
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2532
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2533
16621 continue                                                             2534
16622 continue                                                             2534
11400 continue                                                             2534
      deallocate(w,f)                                                      2535
      return                                                               2536
      end                                                                  2537
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   2538 
     *jerr)
      real x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(nlam)          2539
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 2540
      real, dimension (:), allocatable :: w,f                                   
      if(minval(y) .ge. 0.0)goto 16651                                     2543
      jerr=8888                                                            2543
      return                                                               2543
16651 continue                                                             2544
      allocate(w(1:no),stat=jerr)                                          2545
      allocate(f(1:no),stat=ierr)                                          2545
      jerr=jerr+ierr                                                       2546
      if(jerr.ne.0) return                                                 2547
      w=max(0.0,q)                                                         2547
      sw=sum(w)                                                            2547
      if(sw .gt. 0.0)goto 16671                                            2547
      jerr=9999                                                            2547
      go to 11400                                                          2547
16671 continue                                                             2548
      yb=dot_product(w,y)/sw                                               2548
      fmax=log(huge(y(1))*0.1)                                             2549
16680 do 16681 lam=1,nlam                                                  2549
      f=a0(lam)                                                            2550
16690 do 16691 k=1,nin(lam)                                                2550
      j=ia(k)                                                              2550
      jb=ix(j)                                                             2550
      je=ix(j+1)-1                                                         2551
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         2552
16691 continue                                                             2553
16692 continue                                                             2553
      f=f+g                                                                2554
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   2555
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                2556
16681 continue                                                             2557
16682 continue                                                             2557
11400 continue                                                             2557
      deallocate(w,f)                                                      2558
      return                                                               2559
      end                                                                  2560
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
