"new experimental WLS routine"
"alm0: previous lambda value"
"almc: current lambda value"
"alpha: alpha"
"m: current lambda iteration no."
"no: no of observations"
"ni: no of variables"
"x: x matrix"
"r: weighted residual! v * (y - yhat)"
"v: weights"
"intr: 1 if fitting intercept, 0 otherwise"
"ju: ju(k) = 1 if feature k is included for consideration in model"
"vp: relative penalties on features (sum to ni)"
"cl: coordinate limits"
"nx: limit for max no. of variables ever to be nonzero"
"thr: threshold for dlx"
"maxit: max no. of passes over the data for all lambda values"
"a, aint: warm start for coefs, intercept"
"g: abs(dot_product(r,x(:,j)))"
"ia: mapping nino to k"
"iy: ever-active set (for compatibility with sparse version)"
"iz: flag for loop. 0 for first lambda, 1 subsequently"
"mm: mapping k to nino"
"nino: no. of features that have ever been nonzero"
"nlp: no. of passes over the data"
"jerr: error code"
subroutine wls(alm0,almc,alpha,m,no,ni,x,r,v,intr,ju,vp,cl,nx,thr,maxit,
   a,aint,g,ia,iy,iz,mm,nino,rsqc,nlp,jerr);
implicit double precision(a-h,o-z);
double precision x(no,ni),r(no),a(ni),vp(ni),cl(2,ni);
double precision v(no),g(ni);
integer iy(ni),ia(nx),ju(ni),mm(ni);
%fortran
      double precision, dimension (:), allocatable :: xv
%mortran
allocate(xv(1:ni),stat=jerr); if(jerr.ne.0) return;

"compute g initialization"
<j=1,ni; if(ju(j).eq.0) next; g(j)=abs(dot_product(r,x(:,j)));>

"compute xv"
<j=1,ni; if(iy(j).gt.0) xv(j)=dot_product(v,x(:,j)**2);>

"compute xmz (for intercept later)"
xmz = sum(v);

"ab: lambda * alpha, dem: lambda * (1 - alpha)"
ab=almc*alpha; dem=almc*(1.0-alpha); 

"strong rules: iy(k) = 1 if we don't discard feature k"
tlam=alpha*(2.0*almc-alm0);
<k=1,ni; if(iy(k).eq.1) next; if(ju(k).eq.0) next;
  if(g(k).gt.tlam*vp(k)) <iy(k)=1; xv(k)=dot_product(v,x(:,k)**2);>
>

jz = 1;

loop < if(iz*jz.ne.0) go to :b:;

   :again: nlp=nlp+1; dlx=0.0;

   <k=1,ni;
      "if feature discarded by strong rules, skip over it"  
      if(iy(k).eq.0) next;

      "check if ST threshold for descent is met"
      "if it goes to 0, set a(k) = 0.0"
      "if not, set a(k) to the post-gradient descent value"
      "u is the kth partial residual"
      gk=dot_product(r,x(:,k));
      ak=a(k); 
      u=gk+ak*xv(k); 
      au=abs(u)-vp(k)*ab;
      if au.le.0.0 < a(k)=0.0;>
      else <
         a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*dem)));
      >

      "if the update didn't change the coefficient value, go to"
      "the next feature/variable"
      if(a(k).eq.ak) next;

      "if coef for feature k was previously 0, we now have a "
      "new non-zero coef. update nino, mm(k) and ia(nino)."
      if mm(k).eq.0 < nino=nino+1; 
         if(nino.gt.nx) exit;
         mm(k)=nino; ia(nino)=k;
      >

      "update residual r, rsqc, and dlx (diff exp from wls)"
      d=a(k)-ak;
      rsqc=rsqc+d*(2.0*gk-d*xv(k));
      r=r-d*v*x(:,k); 
      dlx=max(xv(k)*d**2,dlx);
   >

   "if we've gone over max no. of vars allowed to enter all"
   "models, leave the loop"
   if(nino.gt.nx) exit;

   "updating of intercept term"
   d=0.0; if(intr.ne.0) d=sum(r)/xmz;
   if d.ne.0.0 <
      aint=aint+d;
      rsqc=rsqc+d*(2.0*sum(r)-d*xmz);
      dlx=max(dlx,xmz*d**2);
      r=r-d*v;>

   "in wls, this leads to KKT checks. here, we exit"
   "the loop instead."
   if(dlx.lt.thr) < ixx=0;
        <k=1,ni; if(iy(k).eq.1) next; if(ju(k).eq.0) next;
            g(k)=abs(dot_product(r,x(:,k)));
            if g(k).gt.ab*vp(k) < 
               iy(k)=1; xv(k)=dot_product(v,x(:,k)**2);
               ixx=1;
            >
        >
        if(ixx.eq.1) go to :again:;
        exit;
    >

   "if we've gone over max iterations, return w error"
   if nlp.gt.maxit < jerr=-m; return;>

   "this is like the :b: loop in wls (M)"
   :b: iz = 1;
   loop < nlp=nlp+1; dlx=0.0;
      <l=1,nino; k=ia(l); 
         gk=dot_product(r,x(:,k));
         ak=a(k); u=gk+ak*xv(k); au=abs(u)-vp(k)*ab;
         if au.le.0.0 < a(k)=0.0;>
         else <
           a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*dem)));
         >
         if(a(k).eq.ak) next;
         d=a(k)-ak; 
         rsqc=rsqc+d*(2.0*gk-d*xv(k));
         r=r-d*v*x(:,k); dlx=max(xv(k)*d**2,dlx);
      >
      "updating of intercept term"
      d=0.0; if(intr.ne.0) d=sum(r)/xmz;
      if d.ne.0.0 <
        aint=aint+d;
        rsqc=rsqc+d*(2.0*sum(r)-d*xmz);
        dlx=max(dlx,xmz*d**2);
        r=r-d*v;>

      if(dlx.lt.thr) exit;
      if nlp.gt.maxit < jerr=-m; return;>
   >
   "set jz = 0 so that we have to go into :again: tag"
   "to check KKT conditions."
   jz=0;
>

deallocate(xv);
return;
end;

"new experimental WLS routine for sparse x"
"Differences in arguments from wls subroutine:"
"x: x@x instead of full matrix x"
"ix: x@p + 1"
"jx: x@i + 1"
"xm: vector of feature means"
"xs: vector of feature SDs"
"iy: ever-active set (previously ix)"
subroutine spwls(alm0,almc,alpha,m,no,ni,x,ix,jx,xm,xs,r,v,intr,ju,
    vp,cl,nx,thr,maxit,a,aint,g,ia,iy,iz,mm,nino,rsqc,nlp,jerr);
implicit double precision(a-h,o-z);
double precision x(*),xm(ni),xs(ni),r(no),a(ni),vp(ni),cl(2,ni);
double precision v(no),g(ni);
integer ix(*),jx(*),iy(ni),ia(nx),ju(ni),mm(ni);
%fortran
      double precision, dimension (:), allocatable :: xv
%mortran
allocate(xv(1:ni),stat=jerr); if(jerr.ne.0) return;

"compute xmz and sum of r (for intercept later)"
xmz = sum(v);
rsum = sum(r);

"compute g initialization"
<j=1,ni; if(ju(j).eq.0) next;
    jb=ix(j); je=ix(j+1)-1;
    g(j)=abs(dot_product(r(jx(jb:je)),x(jb:je))-rsum*xm(j))/xs(j);
>

"compute xv"
<j=1,ni; if(iy(j).gt.0) 
    <jb=ix(j); je=ix(j+1)-1;
        xv(j)=dot_product(v(jx(jb:je)),x(jb:je)**2);
        xv(j)=xv(j)-2*xm(j)*dot_product(v(jx(jb:je)),x(jb:je));
        xv(j)=(xv(j)+xmz*xm(j)**2)/xs(j)**2;
    >
>

"ab: lambda * alpha, dem: lambda * (1 - alpha)"
ab=almc*alpha; dem=almc*(1.0-alpha); 

"strong rules: ix(k) = 1 if we don't discard feature k"
tlam=alpha*(2.0*almc-alm0);
<k=1,ni; if(iy(k).eq.1) next; if(ju(k).eq.0) next;
  if(g(k).gt.tlam*vp(k)) <iy(k)=1; 
      jb=ix(k); je=ix(k+1)-1;
      xv(k)=dot_product(v(jx(jb:je)),x(jb:je)**2);
      xv(k)=xv(k)-2*xm(k)*dot_product(v(jx(jb:je)),x(jb:je));
      xv(k)=(xv(k)+xmz*xm(k)**2)/xs(k)**2;
  >
>

jz = 1;

loop < if(iz*jz.ne.0) go to :b:;

   :again: nlp=nlp+1; dlx=0.0;

   <k=1,ni;
      "if feature discarded by strong rules, skip over it"  
      if(iy(k).eq.0) next;

      "check if ST threshold for descent is met"
      "if it goes to 0, set a(k) = 0.0"
      "if not, set a(k) to the post-gradient descent value"
      "u is the kth partial residual"
      jb=ix(k); je=ix(k+1)-1;
      gk=(dot_product(r(jx(jb:je)),x(jb:je))-rsum*xm(k))/xs(k);
      ak=a(k); 
      u=gk+ak*xv(k); 
      au=abs(u)-vp(k)*ab;
      if au.le.0.0 < a(k)=0.0;>
      else <
         a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*dem)));
      >

      "if the update didn't change the coefficient value, go to"
      "the next feature/variable"
      if(a(k).eq.ak) next;

      "if coef for feature k was previously 0, we now have a "
      "new non-zero coef. update nino, mm(k) and ia(nino)."
      if mm(k).eq.0 < nino=nino+1; 
         if(nino.gt.nx) exit;
         mm(k)=nino; ia(nino)=k;
      >

      "update residual r, rsqc, and dlx (diff exp from wls)"
      d=a(k)-ak;
      rsqc=rsqc+d*(2.0*gk-d*xv(k));
      jb=ix(k); je=ix(k+1)-1;
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k);
      r=r+d*v*xm(k)/xs(k);
      rsum=sum(r);
      dlx=max(xv(k)*d**2,dlx);
   >

   "if we've gone over max no. of vars allowed to enter all"
   "models, leave the loop"
   if(nino.gt.nx) exit;

   "updating of intercept term"
   d=0.0; if(intr.ne.0) d=rsum/xmz;
   if d.ne.0.0 <
      aint=aint+d;
      rsqc=rsqc+d*(2.0*rsum-d*xmz);
      dlx=max(dlx,xmz*d**2);
      r=r-d*v;
      rsum=sum(r);
  >

   "in wls, this leads to KKT checks. here, we exit"
   "the loop instead."
   if(dlx.lt.thr) < ixx=0;
        <k=1,ni; if(iy(k).eq.1) next; if(ju(k).eq.0) next;
            jb=ix(k); je=ix(k+1)-1;
            g(k)=dot_product(r(jx(jb:je)),x(jb:je));
            g(k)=abs(g(k)-rsum*xm(k))/xs(k);
            if g(k).gt.ab*vp(k) < 
                iy(k)=1;
                xv(k)=dot_product(v(jx(jb:je)),x(jb:je)**2);
                vx=dot_product(v(jx(jb:je)),x(jb:je));
                xv(k)=xv(k)-2*xm(k)*vx;
                xv(k)=(xv(k)+xmz*xm(k)**2)/xs(k)**2;
                ixx=1;
            >
        >
        if(ixx.eq.1) go to :again:;
        exit;
    >

   "if we've gone over max iterations, return w error"
   if nlp.gt.maxit < jerr=-m; return;>

   "this is like the :b: loop in wls (M)"
   :b: iz = 1;
   loop < nlp=nlp+1; dlx=0.0;
      <l=1,nino; k=ia(l); 
         jb=ix(k); je=ix(k+1)-1;
         gk=(dot_product(r(jx(jb:je)),x(jb:je))-rsum*xm(k))/xs(k);
         ak=a(k); u=gk+ak*xv(k); au=abs(u)-vp(k)*ab;
         if au.le.0.0 < a(k)=0.0;>
         else <
           a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*dem)));
         >
         if(a(k).eq.ak) next;
         d=a(k)-ak; 
         rsqc=rsqc+d*(2.0*gk-d*xv(k));
         jb=ix(k); je=ix(k+1)-1;
         r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k);
         r=r+d*v*xm(k)/xs(k);
         rsum=sum(r);
         dlx=max(xv(k)*d**2,dlx);
      >
      "updating of intercept term"
      d=0.0; if(intr.ne.0) d=rsum/xmz;
      if d.ne.0.0 <
        aint=aint+d;
        rsqc=rsqc+d*(2.0*rsum-d*xmz);
        dlx=max(dlx,xmz*d**2);
        r=r-d*v;
        rsum=rsum-d*xmz;>

      if(dlx.lt.thr) exit;
      if nlp.gt.maxit < jerr=-m; return;>
   >
   "set jz = 0 so that we have to go into :again: tag"
   "to check KKT conditions."
   jz=0;
>

deallocate(xv);
return;
end;

subroutine get_int_parms2(epsnr,mxitnr);
implicit double precision(a-h,o-z);
data epsnr0,mxitnr0
  /1.0d-8,25/;
epsnr=epsnr0; mxitnr=mxitnr0; 
return;
entry chg_epsnr(arg); epsnr0=arg; return;
entry chg_mxitnr(irg); mxitnr0=irg; return;
end;

%%
