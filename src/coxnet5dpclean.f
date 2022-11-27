c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))             
      subroutine coxnet(parm,no,ni,x,y,d,g,w,jd,vp,cl,ne,nx,nlam,flmin,u
     *lam,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam
     *)
      double precision ca(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)         
      integer jd(*),ia(nx),nin(nlam)                                    
      double precision, dimension (:), allocatable :: xs,ww,vq          
      integer, dimension (:), allocatable :: ju                         
      if(maxval(vp) .gt. 0.0)goto 10021                                 
      jerr=10000                                                        
      return                                                            
10021 continue                                                          
      allocate(ww(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      if(isd .le. 0)goto 10041                                          
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
10041 continue                                                          
      call chkvars(no,ni,x,ju)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 10061                                   
      jerr=7777                                                         
      return                                                            
10061 continue                                                          
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      ww=max(0d0,w)                                                     
      sw=sum(ww)                                                        
      if(sw .gt. 0.0)goto 10081                                         
      jerr=9999                                                         
      return                                                            
10081 continue                                                          
      ww=ww/sw                                                          
      call cstandard(no,ni,x,ww,ju,isd,xs)                              
      if(isd .le. 0)goto 10101                                          
      do 10111 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
10111 continue                                                          
      continue                                                          
10101 continue                                                          
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,
     *thr,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                              
      dev0=2.0*sw*dev0                                                  
      if(isd .le. 0)goto 10131                                          
      do 10141 k=1,lmu                                                  
      nk=nin(k)                                                         
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                
10141 continue                                                          
      continue                                                          
10131 continue                                                          
      deallocate(ww,ju,vq)                                              
      if(isd.gt.0) deallocate(xs)                                       
      return                                                            
      end                                                               
      subroutine cstandard(no,ni,x,w,ju,isd,xs)                         
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),w(no),xs(ni)                            
      integer ju(ni)                                                    
      do 10151 j=1,ni                                                   
      if(ju(j).eq.0)goto 10151                                          
      xm=dot_product(w,x(:,j))                                          
      x(:,j)=x(:,j)-xm                                                  
      if(isd .le. 0)goto 10171                                          
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                              
      x(:,j)=x(:,j)/xs(j)                                               
10171 continue                                                          
10151 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,cl,ne,nx,nlam,flmin,
     *ulam,cthri,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam
     *)
      double precision ao(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)         
      integer ju(ni),m(nx),kin(nlam)                                    
      double precision, dimension (:), allocatable :: w,dk,v,xs,wr      
      double precision, dimension (:), allocatable :: a,as,f,dq         
      double precision, dimension (:), allocatable :: e,uu,ga           
      integer, dimension (:), allocatable :: jp,kp,mm,ixx               
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx,itrace)     
      isd = isd*1                                                       
      sml=sml*100.0                                                     
      devmax=devmax*0.99/0.999                                          
      allocate(e(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(uu(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(f(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(w(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(v(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(a(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(as(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ga(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ixx(1:ni),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(jp(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(kp(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(dk(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(wr(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(dq(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                            
      if(jerr.ne.0) go to 10180                                         
      alpha=parm                                                        
      oma=1.0-alpha                                                     
      nlm=0                                                             
      ixx=0                                                             
      al=0.0                                                            
      dq=d*q                                                            
      call died(no,nk,dq,kp,jp,dk)                                      
      a=0.0                                                             
      f(1)=0.0                                                          
      fmax=log(huge(f(1))*0.1)                                          
      if(nonzero(no,g) .eq. 0)goto 10201                                
      f=g-dot_product(q,g)                                              
      e=q*exp(sign(min(abs(f),fmax),f))                                 
      goto 10211                                                        
10201 continue                                                          
      f=0.0                                                             
      e=q                                                               
10211 continue                                                          
      continue                                                          
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                              
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                      
      dev0=rr                                                           
      do 10221 i=1,no                                                   
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 10241                
      w(i)=0.0                                                          
      wr(i)=w(i)                                                        
10241 continue                                                          
10221 continue                                                          
      continue                                                          
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                      
      if(jerr.ne.0) go to 10180                                         
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 10261                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
10261 continue                                                          
      m=0                                                               
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      mnl=min(mnlam,nlam)                                               
      as=0.0                                                            
      cthr=cthri*dev0                                                   
      do 10271 j=1,ni                                                   
      if(ju(j).eq.0)goto 10271                                          
      ga(j)=abs(dot_product(wr,x(:,j)))                                 
10271 continue                                                          
      continue                                                          
      do 10281 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 10301                                      
      al=ulam(ilm)                                                      
      goto 10291                                                        
10301 if(ilm .le. 2)goto 10311                                          
      al=al*alf                                                         
      goto 10291                                                        
10311 if(ilm .ne. 1)goto 10321                                          
      al=big                                                            
      goto 10331                                                        
10321 continue                                                          
      al0=0.0                                                           
      do 10341 j=1,ni                                                   
      if(ju(j).eq.0)goto 10341                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
10341 continue                                                          
      continue                                                          
      al0=al0/max(parm,1.0d-3)                                          
      al=alf*al0                                                        
10331 continue                                                          
10291 continue                                                          
      sa=alpha*al                                                       
      omal=oma*al                                                       
      tlam=alpha*(2.0*al-al0)                                           
      do 10351 k=1,ni                                                   
      if(ixx(k).eq.1)goto 10351                                         
      if(ju(k).eq.0)goto 10351                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
10351 continue                                                          
      continue                                                          
10360 continue                                                          
      continue                                                          
10371 continue                                                          
      if(nlp .le. maxit)goto 10391                                      
      jerr=-ilm                                                         
      return                                                            
10391 continue                                                          
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                             
      call vars(no,ni,x,w,ixx,v)                                        
      continue                                                          
10401 continue                                                          
      nlp=nlp+1                                                         
      dli=0.0                                                           
      do 10411 j=1,ni                                                   
      if(ixx(j).eq.0)goto 10411                                         
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                
      if(abs(u) .gt. vp(j)*sa)goto 10431                                
      at=0.0                                                            
      goto 10441                                                        
10431 continue                                                          
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o
     *mal)))
10441 continue                                                          
      continue                                                          
      if(at .eq. a(j))goto 10461                                        
      del=at-a(j)                                                       
      a(j)=at                                                           
      dli=max(dli,v(j)*del**2)                                          
      wr=wr-del*w*x(:,j)                                                
      f=f+del*x(:,j)                                                    
      if(mm(j) .ne. 0)goto 10481                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 10412                                           
      mm(j)=nin                                                         
      m(nin)=j                                                          
10481 continue                                                          
10461 continue                                                          
10411 continue                                                          
10412 continue                                                          
      if(nin.gt.nx)goto 10402                                           
      if(dli.lt.cthr)goto 10402                                         
      if(nlp .le. maxit)goto 10501                                      
      jerr=-ilm                                                         
      return                                                            
10501 continue                                                          
      continue                                                          
10511 continue                                                          
      nlp=nlp+1                                                         
      dli=0.0                                                           
      do 10521 l=1,nin                                                  
      j=m(l)                                                            
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                
      if(abs(u) .gt. vp(j)*sa)goto 10541                                
      at=0.0                                                            
      goto 10551                                                        
10541 continue                                                          
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o
     *mal)))
10551 continue                                                          
      continue                                                          
      if(at .eq. a(j))goto 10571                                        
      del=at-a(j)                                                       
      a(j)=at                                                           
      dli=max(dli,v(j)*del**2)                                          
      wr=wr-del*w*x(:,j)                                                
      f=f+del*x(:,j)                                                    
10571 continue                                                          
10521 continue                                                          
      continue                                                          
      if(dli.lt.cthr)goto 10512                                         
      if(nlp .le. maxit)goto 10591                                      
      jerr=-ilm                                                         
      return                                                            
10591 continue                                                          
      goto 10511                                                        
10512 continue                                                          
      goto 10401                                                        
10402 continue                                                          
      if(nin.gt.nx)goto 10372                                           
      e=q*exp(sign(min(abs(f),fmax),f))                                 
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                      
      if(jerr .eq. 0)goto 10611                                         
      jerr=jerr-ilm                                                     
      go to 10180                                                       
10611 continue                                                          
      ix=0                                                              
      do 10621 j=1,nin                                                  
      k=m(j)                                                            
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 10621                        
      ix=1                                                              
      goto 10622                                                        
10621 continue                                                          
10622 continue                                                          
      if(ix .ne. 0)goto 10641                                           
      do 10651 k=1,ni                                                   
      if(ixx(k).eq.1)goto 10651                                         
      if(ju(k).eq.0)goto 10651                                          
      ga(k)=abs(dot_product(wr,x(:,k)))                                 
      if(ga(k) .le. sa*vp(k))goto 10671                                 
      ixx(k)=1                                                          
      ix=1                                                              
10671 continue                                                          
10651 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10360                                           
      goto 10372                                                        
10641 continue                                                          
      goto 10371                                                        
10372 continue                                                          
      if(nin .le. nx)goto 10691                                         
      jerr=-10000-ilm                                                   
      goto 10282                                                        
10691 continue                                                          
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                            
      kin(ilm)=nin                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                
      if(ilm.lt.mnl)goto 10281                                          
      if(flmin.ge.1.0)goto 10281                                        
      me=0                                                              
      do 10701 j=1,nin                                                  
      if(ao(j,ilm).ne.0.0) me=me+1                                      
10701 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 10282                                            
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 10282           
      if(dev(ilm).gt.devmax)goto 10282                                  
10281 continue                                                          
10282 continue                                                          
      g=f                                                               
10180 continue                                                          
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)           
      return                                                            
      end                                                               
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                      
      implicit double precision(a-h,o-z)                                
      double precision y(no),d(no),q(no)                                
      integer jp(no),kp(*)                                              
      do 10711 j=1,no                                                   
      jp(j)=j                                                           
10711 continue                                                          
      continue                                                          
      call psort7(y,jp,1,no)                                            
      nj=0                                                              
      do 10721 j=1,no                                                   
      if(q(jp(j)).le.0.0)goto 10721                                     
      nj=nj+1                                                           
      jp(nj)=jp(j)                                                      
10721 continue                                                          
      continue                                                          
      if(nj .ne. 0)goto 10741                                           
      jerr=20000                                                        
      return                                                            
10741 continue                                                          
      j=1                                                               
      continue                                                          
10751 if(d(jp(j)).gt.0.0)goto 10752                                     
      j=j+1                                                             
      if(j.gt.nj)goto 10752                                             
      goto 10751                                                        
10752 continue                                                          
      if(j .lt. nj-1)goto 10771                                         
      jerr=30000                                                        
      return                                                            
10771 continue                                                          
      t0=y(jp(j))                                                       
      j0=j-1                                                            
      if(j0 .le. 0)goto 10791                                           
      continue                                                          
10801 if(y(jp(j0)).lt.t0)goto 10802                                     
      j0=j0-1                                                           
      if(j0.eq.0)goto 10802                                             
      goto 10801                                                        
10802 continue                                                          
      if(j0 .le. 0)goto 10821                                           
      nj=nj-j0                                                          
      do 10831 j=1,nj                                                   
      jp(j)=jp(j+j0)                                                    
10831 continue                                                          
      continue                                                          
10821 continue                                                          
10791 continue                                                          
      jerr=0                                                            
      nk=0                                                              
      yk=t0                                                             
      j=2                                                               
      continue                                                          
10841 continue                                                          
      continue                                                          
10851 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 10852                  
      j=j+1                                                             
      if(j.gt.nj)goto 10852                                             
      goto 10851                                                        
10852 continue                                                          
      nk=nk+1                                                           
      kp(nk)=j-1                                                        
      if(j.gt.nj)goto 10842                                             
      if(j .ne. nj)goto 10871                                           
      nk=nk+1                                                           
      kp(nk)=nj                                                         
      goto 10842                                                        
10871 continue                                                          
      yk=y(jp(j))                                                       
      j=j+1                                                             
      goto 10841                                                        
10842 continue                                                          
      return                                                            
      end                                                               
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                  
      implicit double precision(a-h,o-z)                                
      double precision d(no),dk(nk),wr(no),w(no)                        
      double precision e(no),u(no),b,c                                  
      integer kp(nk),jp(no)                                             
      call usk(no,nk,kp,jp,e,u)                                         
      b=dk(1)/u(1)                                                      
      c=dk(1)/u(1)**2                                                   
      jerr=0                                                            
      do 10881 j=1,kp(1)                                                
      i=jp(j)                                                           
      w(i)=e(i)*(b-e(i)*c)                                              
      if(w(i) .gt. 0.0)goto 10901                                       
      jerr=-30000                                                       
      return                                                            
10901 continue                                                          
      wr(i)=d(i)-e(i)*b                                                 
10881 continue                                                          
      continue                                                          
      do 10911 k=2,nk                                                   
      j1=kp(k-1)+1                                                      
      j2=kp(k)                                                          
      b=b+dk(k)/u(k)                                                    
      c=c+dk(k)/u(k)**2                                                 
      do 10921 j=j1,j2                                                  
      i=jp(j)                                                           
      w(i)=e(i)*(b-e(i)*c)                                              
      if(w(i) .gt. 0.0)goto 10941                                       
      jerr=-30000                                                       
      return                                                            
10941 continue                                                          
      wr(i)=d(i)-e(i)*b                                                 
10921 continue                                                          
      continue                                                          
10911 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine vars(no,ni,x,w,ixx,v)                                  
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),w(no),v(ni)                             
      integer ixx(ni)                                                   
      do 10951 j=1,ni                                                   
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                     
10951 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine died(no,nk,d,kp,jp,dk)                                 
      implicit double precision(a-h,o-z)                                
      double precision d(no),dk(nk)                                     
      integer kp(nk),jp(no)                                             
      dk(1)=sum(d(jp(1:kp(1))))                                         
      do 10961 k=2,nk                                                   
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                               
10961 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine usk(no,nk,kp,jp,e,u)                                   
      implicit double precision(a-h,o-z)                                
      double precision e(no),u(nk),h                                    
      integer kp(nk),jp(no)                                             
      h=0.0                                                             
      do 10971 k=nk,1,-1                                                
      j2=kp(k)                                                          
      j1=1                                                              
      if(k.gt.1) j1=kp(k-1)+1                                           
      do 10981 j=j2,j1,-1                                               
      h=h+e(jp(j))                                                      
10981 continue                                                          
      continue                                                          
      u(k)=h                                                            
10971 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                          
      implicit double precision(a-h,o-z)                                
      double precision d(no),dk(nk),f(no)                               
      integer kp(nk),jp(no)                                             
      double precision e(no),u(nk)                                      
      ni = ni*1                                                         
      call usk(no,nk,kp,jp,e,u)                                         
      u=log(u)                                                          
      risk=dot_product(d,f)-dot_product(dk,u)                           
      return                                                            
      end                                                               
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)              
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(
     *nlam)
      double precision, dimension (:), allocatable :: dk,f,xm,dq,q      
      double precision, dimension (:), allocatable :: e,uu              
      integer, dimension (:), allocatable :: jp,kp                      
      allocate(e(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(q(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(uu(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(f(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(dk(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(jp(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(kp(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(dq(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      q=max(0d0,w)                                                      
      sw=sum(q)                                                         
      if(sw .gt. 0.0)goto 11001                                         
      jerr=9999                                                         
      go to 10180                                                       
11001 continue                                                          
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                            
      if(jerr.ne.0) go to 10180                                         
      fmax=log(huge(e(1))*0.1)                                          
      dq=d*q                                                            
      call died(no,nk,dq,kp,jp,dk)                                      
      gm=dot_product(q,g)/sw                                            
      do 11011 j=1,ni                                                   
      xm(j)=dot_product(q,x(:,j))/sw                                    
11011 continue                                                          
      continue                                                          
      do 11021 lam=1,nlam                                               
      do 11031 i=1,no                                                   
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                    
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                     
11031 continue                                                          
      continue                                                          
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                       
11021 continue                                                          
      continue                                                          
10180 continue                                                          
      deallocate(e,uu,dk,f,jp,kp,dq)                                    
      return                                                            
      end                                                               
      subroutine chkvars(no,ni,x,ju)                                    
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni)                                         
      integer ju(ni)                                                    
      do 11041 j=1,ni                                                   
      ju(j)=0                                                           
      t=x(1,j)                                                          
      do 11051 i=2,no                                                   
      if(x(i,j).eq.t)goto 11051                                         
      ju(j)=1                                                           
      goto 11052                                                        
11051 continue                                                          
11052 continue                                                          
11041 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      function bnorm(b0,al1p,al2p,g,usq,jerr)                           
      implicit double precision(a-h,o-z)                                
      data thr,mxit /1.0d-10,100/                                       
      b=b0                                                              
      zsq=b**2+usq                                                      
      if(zsq .gt. 0.0)goto 11071                                        
      bnorm=0.0                                                         
      return                                                            
11071 continue                                                          
      z=sqrt(zsq)                                                       
      f=b*(al1p+al2p/z)-g                                               
      jerr=0                                                            
      do 11081 it=1,mxit                                                
      b=b-f/(al1p+al2p*usq/(z*zsq))                                     
      zsq=b**2+usq                                                      
      if(zsq .gt. 0.0)goto 11101                                        
      bnorm=0.0                                                         
      return                                                            
11101 continue                                                          
      z=sqrt(zsq)                                                       
      f=b*(al1p+al2p/z)-g                                               
      if(abs(f).le.thr)goto 11082                                       
      if(b .gt. 0.0)goto 11121                                          
      b=0.0                                                             
      goto 11082                                                        
11121 continue                                                          
11081 continue                                                          
11082 continue                                                          
      bnorm=b                                                           
      if(it.ge.mxit) jerr=90000                                         
      return                                                            
      entry chg_bnorm(arg,irg)                                          
      bnorm = 0.0                                                       
      thr=arg                                                           
      mxit=irg                                                          
      return                                                            
      entry get_bnorm(arg,irg)                                          
      bnorm = 0.0                                                       
      arg=thr                                                           
      irg=mxit                                                          
      return                                                            
      end                                                               
      subroutine get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx,itrace
     *)
      implicit double precision(a-h,o-z)                                
      data sml0,eps0,big0,mnlam0,rsqmax0,pmin0,exmx0,itrace0  /1.0d-5,1.
     *0d-6,9.9d35,5,0.999,1.0d-9,250.0,0/
      sml=sml0                                                          
      eps=eps0                                                          
      big=big0                                                          
      mnlam=mnlam0                                                      
      rsqmax=rsqmax0                                                    
      pmin=pmin0                                                        
      exmx=exmx0                                                        
      itrace=itrace0                                                    
      return                                                            
      entry chg_fract_dev(arg)                                          
      sml0=arg                                                          
      return                                                            
      entry chg_dev_max(arg)                                            
      rsqmax0=arg                                                       
      return                                                            
      entry chg_min_flmin(arg)                                          
      eps0=arg                                                          
      return                                                            
      entry chg_big(arg)                                                
      big0=arg                                                          
      return                                                            
      entry chg_min_lambdas(irg)                                        
      mnlam0=irg                                                        
      return                                                            
      entry chg_min_null_prob(arg)                                      
      pmin0=arg                                                         
      return                                                            
      entry chg_max_exp(arg)                                            
      exmx0=arg                                                         
      return                                                            
      entry chg_itrace(irg)                                             
      itrace0=irg                                                       
      return                                                            
      end                                                               
      subroutine get_int_parms2(epsnr,mxitnr)                           
      implicit double precision(a-h,o-z)                                
      data epsnr0,mxitnr0  /1.0d-6,25/                                  
      epsnr=epsnr0                                                      
      mxitnr=mxitnr0                                                    
      return                                                            
      entry chg_epsnr(arg)                                              
      epsnr0=arg                                                        
      return                                                            
      entry chg_mxitnr(irg)                                             
      mxitnr0=irg                                                       
      return                                                            
      end                                                               
      function nonzero(n,v)                                             
      implicit double precision(a-h,o-z)                                
      double precision v(n)                                             
      nonzero=0                                                         
      do 11131 i=1,n                                                    
      if(v(i) .eq. 0.0)goto 11151                                       
      nonzero=1                                                         
      return                                                            
11151 continue                                                          
11131 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine psort7(v,a,ii,jj)                                      
      implicit double precision(a-h,o-z)                                
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
      double precision v                                                
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
