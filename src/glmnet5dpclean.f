c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))             
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
      subroutine elnet(ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,  flmin,u
     *lam,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam),cl(2,ni) 
      double precision ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)          
      integer jd(*),ia(nx),nin(nlam)                                    
      double precision, dimension (:), allocatable :: vq;               
      if(maxval(vp) .gt. 0.0)goto 10021                                 
      jerr=10000                                                        
      return                                                            
10021 continue                                                          
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      if(ka .ne. 1)goto 10041                                           
      call elnetu  (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,
     *isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                        
10041 continue                                                          
      call elnetn (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,i
     *sd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                          
      continue                                                          
      deallocate(vq)                                                    
      return                                                            
      end                                                               
      subroutine elnetu(parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,  flmin,ula
     *m,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)  
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)         
      integer jd(*),ia(nx),nin(nlam)                                    
      double precision, dimension (:), allocatable :: xm,xs,g,xv,vlam   
      integer, dimension (:), allocatable :: ju                         
      allocate(g(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vlam(1:nlam),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      call chkvars(no,ni,x,ju)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 10071                                   
      jerr=7777                                                         
      return                                                            
10071 continue                                                          
      call standard(no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr)      
      if(jerr.ne.0) return                                              
      cl=cl/ys                                                          
      if(isd .le. 0)goto 10091                                          
      do 10101 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
10101 continue                                                          
      continue                                                          
10091 continue                                                          
      if(flmin.ge.1.0) vlam=ulam/ys                                     
      call elnet1(parm,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxi
     *t,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                              
      do 10111 k=1,lmu                                                  
      alm(k)=ys*alm(k)                                                  
      nk=nin(k)                                                         
      do 10121 l=1,nk                                                   
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                      
10121 continue                                                          
      continue                                                          
      a0(k)=0.0                                                         
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))       
10111 continue                                                          
      continue                                                          
      deallocate(xm,xs,g,ju,xv,vlam)                                    
      return                                                            
      end                                                               
      subroutine standard(no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)  
      integer ju(ni)                                                    
      double precision, dimension (:), allocatable :: v                 
      allocate(v(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      w=w/sum(w)                                                        
      v=sqrt(w)                                                         
      if(intr .ne. 0)goto 10141                                         
      ym=0.0                                                            
      y=v*y                                                             
      ys=sqrt(dot_product(y,y))                                         
      y=y/ys                                                            
      do 10151 j=1,ni                                                   
      if(ju(j).eq.0)goto 10151                                          
      xm(j)=0.0                                                         
      x(:,j)=v*x(:,j)                                                   
      xv(j)=dot_product(x(:,j),x(:,j))                                  
      if(isd .eq. 0)goto 10171                                          
      xbq=dot_product(v,x(:,j))**2                                      
      vc=xv(j)-xbq                                                      
      xs(j)=sqrt(vc)                                                    
      x(:,j)=x(:,j)/xs(j)                                               
      xv(j)=1.0+xbq/vc                                                  
      goto 10181                                                        
10171 continue                                                          
      xs(j)=1.0                                                         
10181 continue                                                          
      continue                                                          
10151 continue                                                          
      continue                                                          
      goto 10191                                                        
10141 continue                                                          
      do 10201 j=1,ni                                                   
      if(ju(j).eq.0)goto 10201                                          
      xm(j)=dot_product(w,x(:,j))                                       
      x(:,j)=v*(x(:,j)-xm(j))                                           
      xv(j)=dot_product(x(:,j),x(:,j))                                  
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                    
10201 continue                                                          
      continue                                                          
      if(isd .ne. 0)goto 10221                                          
      xs=1.0                                                            
      goto 10231                                                        
10221 continue                                                          
      do 10241 j=1,ni                                                   
      if(ju(j).eq.0)goto 10241                                          
      x(:,j)=x(:,j)/xs(j)                                               
10241 continue                                                          
      continue                                                          
      xv=1.0                                                            
10231 continue                                                          
      continue                                                          
      ym=dot_product(w,y)                                               
      y=v*(y-ym)                                                        
      ys=sqrt(dot_product(y,y))                                         
      y=y/ys                                                            
10191 continue                                                          
      continue                                                          
      g=0.0                                                             
      do 10251 j=1,ni                                                   
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                         
10251 continue                                                          
      continue                                                          
      deallocate(v)                                                     
      return                                                            
      end                                                               
      subroutine elnet1(beta,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,ulam,th
     *r,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam)     
      double precision rsqo(nlam),almo(nlam),xv(ni)                     
      double precision cl(2,ni)                                         
      integer ju(ni),ia(nx),kin(nlam)                                   
      double precision, dimension (:), allocatable :: a,da              
      integer, dimension (:), allocatable :: mm                         
      double precision, dimension (:,:), allocatable :: c               
      allocate(c(1:ni,1:nx),stat=jerr)                                  
      if(jerr.ne.0) return;                                             
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx,itrace)     
      allocate(a(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(da(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      bta=beta                                                          
      omb=1.0-bta                                                       
      alm=0.0                                                           
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 10271                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
10271 continue                                                          
      rsq=0.0                                                           
      a=0.0                                                             
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      iz=0                                                              
      mnl=min(mnlam,nlam)                                               
      do 10281 m=1,nlam                                                 
      if(itrace.ne.0) call setpb(m-1)                                   
      if(flmin .lt. 1.0)goto 10301                                      
      alm=ulam(m)                                                       
      goto 10291                                                        
10301 if(m .le. 2)goto 10311                                            
      alm=alm*alf                                                       
      goto 10291                                                        
10311 if(m .ne. 1)goto 10321                                            
      alm=big                                                           
      goto 10331                                                        
10321 continue                                                          
      alm=0.0                                                           
      do 10341 j=1,ni                                                   
      if(ju(j).eq.0)goto 10341                                          
      if(vp(j).le.0.0)goto 10341                                        
      alm=max(alm,abs(g(j))/vp(j))                                      
10341 continue                                                          
      continue                                                          
      alm=alf*alm/max(bta,1.0d-3)                                       
10331 continue                                                          
10291 continue                                                          
      dem=alm*omb                                                       
      ab=alm*bta                                                        
      rsq0=rsq                                                          
      jz=1                                                              
      continue                                                          
10351 continue                                                          
      if(iz*jz.ne.0) go to 10360                                        
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10371 k=1,ni                                                   
      if(ju(k).eq.0)goto 10371                                          
      ak=a(k)                                                           
      u=g(k)+ak*xv(k)                                                   
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 10371                                          
      if(mm(k) .ne. 0)goto 10391                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 10372                                           
      do 10401 j=1,ni                                                   
      if(ju(j).eq.0)goto 10401                                          
      if(mm(j) .eq. 0)goto 10421                                        
      c(j,nin)=c(k,mm(j))                                               
      goto 10401                                                        
10421 continue                                                          
      if(j .ne. k)goto 10441                                            
      c(j,nin)=xv(j)                                                    
      goto 10401                                                        
10441 continue                                                          
      c(j,nin)=dot_product(x(:,j),x(:,k))                               
10401 continue                                                          
      continue                                                          
      mm(k)=nin                                                         
      ia(nin)=k                                                         
10391 continue                                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                  
      dlx=max(xv(k)*del**2,dlx)                                         
      do 10451 j=1,ni                                                   
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                           
10451 continue                                                          
      continue                                                          
10371 continue                                                          
10372 continue                                                          
      if(dlx.lt.thr)goto 10352                                          
      if(nin.gt.nx)goto 10352                                           
      if(nlp .le. maxit)goto 10471                                      
      jerr=-m                                                           
      return                                                            
10471 continue                                                          
10360 continue                                                          
      iz=1                                                              
      da(1:nin)=a(ia(1:nin))                                            
      continue                                                          
10481 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10491 l=1,nin                                                  
      k=ia(l)                                                           
      ak=a(k)                                                           
      u=g(k)+ak*xv(k)                                                   
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 10491                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                  
      dlx=max(xv(k)*del**2,dlx)                                         
      do 10501 j=1,nin                                                  
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                              
10501 continue                                                          
      continue                                                          
10491 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 10482                                          
      if(nlp .le. maxit)goto 10521                                      
      jerr=-m                                                           
      return                                                            
10521 continue                                                          
      goto 10481                                                        
10482 continue                                                          
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                  
      do 10531 j=1,ni                                                   
      if(mm(j).ne.0)goto 10531                                          
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))        
10531 continue                                                          
      continue                                                          
      jz=0                                                              
      goto 10351                                                        
10352 continue                                                          
      if(nin .le. nx)goto 10551                                         
      jerr=-10000-m                                                     
      goto 10282                                                        
10551 continue                                                          
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                             
      kin(m)=nin                                                        
      rsqo(m)=rsq                                                       
      almo(m)=alm                                                       
      lmu=m                                                             
      if(m.lt.mnl)goto 10281                                            
      if(flmin.ge.1.0)goto 10281                                        
      me=0                                                              
      do 10561 j=1,nin                                                  
      if(ao(j,m).ne.0.0) me=me+1                                        
10561 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 10282                                            
      if(rsq-rsq0.lt.sml*rsq)goto 10282                                 
      if(rsq.gt.rsqmax)goto 10282                                       
10281 continue                                                          
10282 continue                                                          
      deallocate(a,mm,c,da)                                             
      return                                                            
      end                                                               
      subroutine elnetn(parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam,
     *thr,isd,  intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision vp(ni),x(no,ni),y(no),w(no),ulam(nlam),cl(2,ni)  
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)         
      integer jd(*),ia(nx),nin(nlam)                                    
      double precision, dimension (:), allocatable :: xm,xs,xv,vlam     
      integer, dimension (:), allocatable :: ju                         
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vlam(1:nlam),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      call chkvars(no,ni,x,ju)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 10581                                   
      jerr=7777                                                         
      return                                                            
10581 continue                                                          
      call standard1(no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)       
      if(jerr.ne.0) return                                              
      cl=cl/ys                                                          
      if(isd .le. 0)goto 10601                                          
      do 10611 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
10611 continue                                                          
      continue                                                          
10601 continue                                                          
      if(flmin.ge.1.0) vlam=ulam/ys                                     
      call elnet2(parm,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,vlam,thr,maxi
     *t,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                              
      do 10621 k=1,lmu                                                  
      alm(k)=ys*alm(k)                                                  
      nk=nin(k)                                                         
      do 10631 l=1,nk                                                   
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                      
10631 continue                                                          
      continue                                                          
      a0(k)=0.0                                                         
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))       
10621 continue                                                          
      continue                                                          
      deallocate(xm,xs,ju,xv,vlam)                                      
      return                                                            
      end                                                               
      subroutine standard1(no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr) 
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)        
      integer ju(ni)                                                    
      double precision, dimension (:), allocatable :: v                 
      allocate(v(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      w=w/sum(w)                                                        
      v=sqrt(w)                                                         
      if(intr .ne. 0)goto 10651                                         
      ym=0.0                                                            
      y=v*y                                                             
      ys=sqrt(dot_product(y,y))                                         
      y=y/ys                                                            
      do 10661 j=1,ni                                                   
      if(ju(j).eq.0)goto 10661                                          
      xm(j)=0.0                                                         
      x(:,j)=v*x(:,j)                                                   
      xv(j)=dot_product(x(:,j),x(:,j))                                  
      if(isd .eq. 0)goto 10681                                          
      xbq=dot_product(v,x(:,j))**2                                      
      vc=xv(j)-xbq                                                      
      xs(j)=sqrt(vc)                                                    
      x(:,j)=x(:,j)/xs(j)                                               
      xv(j)=1.0+xbq/vc                                                  
      goto 10691                                                        
10681 continue                                                          
      xs(j)=1.0                                                         
10691 continue                                                          
      continue                                                          
10661 continue                                                          
      continue                                                          
      go to 10700                                                       
10651 continue                                                          
      do 10711 j=1,ni                                                   
      if(ju(j).eq.0)goto 10711                                          
      xm(j)=dot_product(w,x(:,j))                                       
      x(:,j)=v*(x(:,j)-xm(j))                                           
      xv(j)=dot_product(x(:,j),x(:,j))                                  
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                    
10711 continue                                                          
      continue                                                          
      if(isd .ne. 0)goto 10731                                          
      xs=1.0                                                            
      goto 10741                                                        
10731 continue                                                          
      do 10751 j=1,ni                                                   
      if(ju(j).eq.0)goto 10751                                          
      x(:,j)=x(:,j)/xs(j)                                               
10751 continue                                                          
      continue                                                          
      xv=1.0                                                            
10741 continue                                                          
      continue                                                          
      ym=dot_product(w,y)                                               
      y=v*(y-ym)                                                        
      ys=sqrt(dot_product(y,y))                                         
      y=y/ys                                                            
10700 continue                                                          
      deallocate(v)                                                     
      return                                                            
      end                                                               
      subroutine elnet2(beta,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,ulam,th
     *r,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam)     
      double precision rsqo(nlam),almo(nlam),xv(ni)                     
      double precision cl(2,ni)                                         
      integer ju(ni),ia(nx),kin(nlam)                                   
      double precision, dimension (:), allocatable :: a,g               
      integer, dimension (:), allocatable :: mm,ix                      
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx,itrace)     
      allocate(a(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(g(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(ix(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      bta=beta                                                          
      omb=1.0-bta                                                       
      ix=0                                                              
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 10771                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
10771 continue                                                          
      rsq=0.0                                                           
      a=0.0                                                             
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      iz=0                                                              
      mnl=min(mnlam,nlam)                                               
      alm=0.0                                                           
      do 10781 j=1,ni                                                   
      if(ju(j).eq.0)goto 10781                                          
      g(j)=abs(dot_product(y,x(:,j)))                                   
10781 continue                                                          
      continue                                                          
      do 10791 m=1,nlam                                                 
      if(itrace.ne.0) call setpb(m-1)                                   
      alm0=alm                                                          
      if(flmin .lt. 1.0)goto 10811                                      
      alm=ulam(m)                                                       
      goto 10801                                                        
10811 if(m .le. 2)goto 10821                                            
      alm=alm*alf                                                       
      goto 10801                                                        
10821 if(m .ne. 1)goto 10831                                            
      alm=big                                                           
      goto 10841                                                        
10831 continue                                                          
      alm0=0.0                                                          
      do 10851 j=1,ni                                                   
      if(ju(j).eq.0)goto 10851                                          
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                        
10851 continue                                                          
      continue                                                          
      alm0=alm0/max(bta,1.0d-3)                                         
      alm=alf*alm0                                                      
10841 continue                                                          
10801 continue                                                          
      dem=alm*omb                                                       
      ab=alm*bta                                                        
      rsq0=rsq                                                          
      jz=1                                                              
      tlam=bta*(2.0*alm-alm0)                                           
      do 10861 k=1,ni                                                   
      if(ix(k).eq.1)goto 10861                                          
      if(ju(k).eq.0)goto 10861                                          
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                    
10861 continue                                                          
      continue                                                          
      continue                                                          
10871 continue                                                          
      if(iz*jz.ne.0) go to 10360                                        
10880 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10891 k=1,ni                                                   
      if(ix(k).eq.0)goto 10891                                          
      gk=dot_product(y,x(:,k))                                          
      ak=a(k)                                                           
      u=gk+ak*xv(k)                                                     
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 10891                                          
      if(mm(k) .ne. 0)goto 10911                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 10892                                           
      mm(k)=nin                                                         
      ia(nin)=k                                                         
10911 continue                                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*gk-del*xv(k))                                    
      y=y-del*x(:,k)                                                    
      dlx=max(xv(k)*del**2,dlx)                                         
10891 continue                                                          
10892 continue                                                          
      if(nin.gt.nx)goto 10872                                           
      if(dlx .ge. thr)goto 10931                                        
      ixx=0                                                             
      do 10941 k=1,ni                                                   
      if(ix(k).eq.1)goto 10941                                          
      if(ju(k).eq.0)goto 10941                                          
      g(k)=abs(dot_product(y,x(:,k)))                                   
      if(g(k) .le. ab*vp(k))goto 10961                                  
      ix(k)=1                                                           
      ixx=1                                                             
10961 continue                                                          
10941 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10880                                          
      goto 10872                                                        
10931 continue                                                          
      if(nlp .le. maxit)goto 10981                                      
      jerr=-m                                                           
      return                                                            
10981 continue                                                          
10360 continue                                                          
      iz=1                                                              
      continue                                                          
10991 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 11001 l=1,nin                                                  
      k=ia(l)                                                           
      gk=dot_product(y,x(:,k))                                          
      ak=a(k)                                                           
      u=gk+ak*xv(k)                                                     
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 11001                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*gk-del*xv(k))                                    
      y=y-del*x(:,k)                                                    
      dlx=max(xv(k)*del**2,dlx)                                         
11001 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 10992                                          
      if(nlp .le. maxit)goto 11021                                      
      jerr=-m                                                           
      return                                                            
11021 continue                                                          
      goto 10991                                                        
10992 continue                                                          
      jz=0                                                              
      goto 10871                                                        
10872 continue                                                          
      if(nin .le. nx)goto 11041                                         
      jerr=-10000-m                                                     
      goto 10792                                                        
11041 continue                                                          
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                             
      kin(m)=nin                                                        
      rsqo(m)=rsq                                                       
      almo(m)=alm                                                       
      lmu=m                                                             
      if(m.lt.mnl)goto 10791                                            
      if(flmin.ge.1.0)goto 10791                                        
      me=0                                                              
      do 11051 j=1,nin                                                  
      if(ao(j,m).ne.0.0) me=me+1                                        
11051 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 10792                                            
      if(rsq-rsq0.lt.sml*rsq)goto 10792                                 
      if(rsq.gt.rsqmax)goto 10792                                       
10791 continue                                                          
10792 continue                                                          
      deallocate(a,mm,g,ix)                                             
      return                                                            
      end                                                               
      subroutine chkvars(no,ni,x,ju)                                    
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni)                                         
      integer ju(ni)                                                    
      do 11061 j=1,ni                                                   
      ju(j)=0                                                           
      t=x(1,j)                                                          
      do 11071 i=2,no                                                   
      if(x(i,j).eq.t)goto 11071                                         
      ju(j)=1                                                           
      goto 11072                                                        
11071 continue                                                          
11072 continue                                                          
11061 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine uncomp(ni,ca,ia,nin,a)                                 
      implicit double precision(a-h,o-z)                                
      double precision ca(*),a(ni)                                      
      integer ia(*)                                                     
      a=0.0                                                             
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                               
      return                                                            
      end                                                               
      subroutine modval(a0,ca,ia,nin,n,x,f)                             
      implicit double precision(a-h,o-z)                                
      double precision ca(nin),x(n,*),f(n)                              
      integer ia(nin)                                                   
      f=a0                                                              
      if(nin.le.0) return                                               
      do 11081 i=1,n                                                    
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                   
11081 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine spelnet(ka,parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam, 
     * flmin,ulam,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)      
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)         
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                        
      double precision, dimension (:), allocatable :: vq;               
      if(maxval(vp) .gt. 0.0)goto 11101                                 
      jerr=10000                                                        
      return                                                            
11101 continue                                                          
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      if(ka .ne. 1)goto 11121                                           
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,u
     *lam,thr,isd,  intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 11131                                                        
11121 continue                                                          
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,ul
     *am,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
11131 continue                                                          
      continue                                                          
      deallocate(vq)                                                    
      return                                                            
      end                                                               
      subroutine spelnetu(parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,  f
     *lmin,ulam,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)      
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)         
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                        
      double precision, dimension (:), allocatable :: xm,xs,g,xv,vlam   
      integer, dimension (:), allocatable :: ju                         
      allocate(g(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vlam(1:nlam),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      call spchkvars(no,ni,x,ix,ju)                                     
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 11151                                   
      jerr=7777                                                         
      return                                                            
11151 continue                                                          
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys,xv,jer
     *r)
      if(jerr.ne.0) return                                              
      cl=cl/ys                                                          
      if(isd .le. 0)goto 11171                                          
      do 11181 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
11181 continue                                                          
      continue                                                          
11171 continue                                                          
      if(flmin.ge.1.0) vlam=ulam/ys                                     
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                              
      do 11191 k=1,lmu                                                  
      alm(k)=ys*alm(k)                                                  
      nk=nin(k)                                                         
      do 11201 l=1,nk                                                   
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                      
11201 continue                                                          
      continue                                                          
      a0(k)=0.0                                                         
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))       
11191 continue                                                          
      continue                                                          
      deallocate(xm,xs,g,ju,xv,vlam)                                    
      return                                                            
      end                                                               
      subroutine spstandard(no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys,
     *xv,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)      
      integer ix(*),jx(*),ju(ni)                                        
      jerr = jerr*1                                                     
      w=w/sum(w)                                                        
      if(intr .ne. 0)goto 11221                                         
      ym=0.0                                                            
      ys=sqrt(dot_product(w,y**2))                                      
      y=y/ys                                                            
      do 11231 j=1,ni                                                   
      if(ju(j).eq.0)goto 11231                                          
      xm(j)=0.0                                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                       
      if(isd .eq. 0)goto 11251                                          
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                         
      vc=xv(j)-xbq                                                      
      xs(j)=sqrt(vc)                                                    
      xv(j)=1.0+xbq/vc                                                  
      goto 11261                                                        
11251 continue                                                          
      xs(j)=1.0                                                         
11261 continue                                                          
      continue                                                          
11231 continue                                                          
      continue                                                          
      goto 11271                                                        
11221 continue                                                          
      do 11281 j=1,ni                                                   
      if(ju(j).eq.0)goto 11281                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2              
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                    
11281 continue                                                          
      continue                                                          
      if(isd .ne. 0)goto 11301                                          
      xs=1.0                                                            
      goto 11311                                                        
11301 continue                                                          
      xv=1.0                                                            
11311 continue                                                          
      continue                                                          
      ym=dot_product(w,y)                                               
      y=y-ym                                                            
      ys=sqrt(dot_product(w,y**2))                                      
      y=y/ys                                                            
11271 continue                                                          
      continue                                                          
      g=0.0                                                             
      do 11321 j=1,ni                                                   
      if(ju(j).eq.0)goto 11321                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)        
11321 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision g(ni),vp(ni),x(*),ulam(nlam),w(no)               
      double precision ao(nx,nlam),rsqo(nlam),almo(nlam)                
      double precision xm(ni),xs(ni),xv(ni),cl(2,ni)                    
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                       
      double precision, dimension (:), allocatable :: a,da              
      integer, dimension (:), allocatable :: mm                         
      double precision, dimension (:,:), allocatable :: c               
      allocate(c(1:ni,1:nx),stat=jerr)                                  
      if(jerr.ne.0) return;                                             
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx,itrace)     
      allocate(a(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(da(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      bta=beta                                                          
      omb=1.0-bta                                                       
      alm=0.0                                                           
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 11341                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
11341 continue                                                          
      rsq=0.0                                                           
      a=0.0                                                             
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      iz=0                                                              
      mnl=min(mnlam,nlam)                                               
      do 11351 m=1,nlam                                                 
      if(itrace.ne.0) call setpb(m-1)                                   
      if(flmin .lt. 1.0)goto 11371                                      
      alm=ulam(m)                                                       
      goto 11361                                                        
11371 if(m .le. 2)goto 11381                                            
      alm=alm*alf                                                       
      goto 11361                                                        
11381 if(m .ne. 1)goto 11391                                            
      alm=big                                                           
      goto 11401                                                        
11391 continue                                                          
      alm=0.0                                                           
      do 11411 j=1,ni                                                   
      if(ju(j).eq.0)goto 11411                                          
      if(vp(j).le.0.0)goto 11411                                        
      alm=max(alm,abs(g(j))/vp(j))                                      
11411 continue                                                          
      continue                                                          
      alm=alf*alm/max(bta,1.0d-3)                                       
11401 continue                                                          
11361 continue                                                          
      dem=alm*omb                                                       
      ab=alm*bta                                                        
      rsq0=rsq                                                          
      jz=1                                                              
      continue                                                          
11421 continue                                                          
      if(iz*jz.ne.0) go to 10360                                        
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 11431 k=1,ni                                                   
      if(ju(k).eq.0)goto 11431                                          
      ak=a(k)                                                           
      u=g(k)+ak*xv(k)                                                   
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 11431                                          
      if(mm(k) .ne. 0)goto 11451                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 11432                                           
      do 11461 j=1,ni                                                   
      if(ju(j).eq.0)goto 11461                                          
      if(mm(j) .eq. 0)goto 11481                                        
      c(j,nin)=c(k,mm(j))                                               
      goto 11461                                                        
11481 continue                                                          
      if(j .ne. k)goto 11501                                            
      c(j,nin)=xv(j)                                                    
      goto 11461                                                        
11501 continue                                                          
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))    
11461 continue                                                          
      continue                                                          
      mm(k)=nin                                                         
      ia(nin)=k                                                         
11451 continue                                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                  
      dlx=max(xv(k)*del**2,dlx)                                         
      do 11511 j=1,ni                                                   
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                           
11511 continue                                                          
      continue                                                          
11431 continue                                                          
11432 continue                                                          
      if(dlx.lt.thr)goto 11422                                          
      if(nin.gt.nx)goto 11422                                           
      if(nlp .le. maxit)goto 11531                                      
      jerr=-m                                                           
      return                                                            
11531 continue                                                          
10360 continue                                                          
      iz=1                                                              
      da(1:nin)=a(ia(1:nin))                                            
      continue                                                          
11541 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 11551 l=1,nin                                                  
      k=ia(l)                                                           
      ak=a(k)                                                           
      u=g(k)+ak*xv(k)                                                   
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 11551                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                  
      dlx=max(xv(k)*del**2,dlx)                                         
      do 11561 j=1,nin                                                  
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                              
11561 continue                                                          
      continue                                                          
11551 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 11542                                          
      if(nlp .le. maxit)goto 11581                                      
      jerr=-m                                                           
      return                                                            
11581 continue                                                          
      goto 11541                                                        
11542 continue                                                          
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                  
      do 11591 j=1,ni                                                   
      if(mm(j).ne.0)goto 11591                                          
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))        
11591 continue                                                          
      continue                                                          
      jz=0                                                              
      goto 11421                                                        
11422 continue                                                          
      if(nin .le. nx)goto 11611                                         
      jerr=-10000-m                                                     
      goto 11352                                                        
11611 continue                                                          
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                             
      kin(m)=nin                                                        
      rsqo(m)=rsq                                                       
      almo(m)=alm                                                       
      lmu=m                                                             
      if(m.lt.mnl)goto 11351                                            
      if(flmin.ge.1.0)goto 11351                                        
      me=0                                                              
      do 11621 j=1,nin                                                  
      if(ao(j,m).ne.0.0) me=me+1                                        
11621 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 11352                                            
      if(rsq-rsq0.lt.sml*rsq)goto 11352                                 
      if(rsq.gt.rsqmax)goto 11352                                       
11351 continue                                                          
11352 continue                                                          
      deallocate(a,mm,c,da)                                             
      return                                                            
      end                                                               
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flm
     *in,ulam,  thr,isd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),vp(ni),y(no),w(no),ulam(nlam),cl(2,ni)      
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)         
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                        
      double precision, dimension (:), allocatable :: xm,xs,xv,vlam     
      integer, dimension (:), allocatable :: ju                         
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vlam(1:nlam),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      call spchkvars(no,ni,x,ix,ju)                                     
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 11641                                   
      jerr=7777                                                         
      return                                                            
11641 continue                                                          
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,xv,jerr
     *)
      if(jerr.ne.0) return                                              
      cl=cl/ys                                                          
      if(isd .le. 0)goto 11661                                          
      do 11671 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
11671 continue                                                          
      continue                                                          
11661 continue                                                          
      if(flmin.ge.1.0) vlam=ulam/ys                                     
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                              
      do 11681 k=1,lmu                                                  
      alm(k)=ys*alm(k)                                                  
      nk=nin(k)                                                         
      do 11691 l=1,nk                                                   
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                      
11691 continue                                                          
      continue                                                          
      a0(k)=0.0                                                         
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))       
11681 continue                                                          
      continue                                                          
      deallocate(xm,xs,ju,xv,vlam)                                      
      return                                                            
      end                                                               
      subroutine spstandard1(no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,x
     *v,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)            
      integer ix(*),jx(*),ju(ni)                                        
      jerr = jerr*1                                                     
      w=w/sum(w)                                                        
      if(intr .ne. 0)goto 11711                                         
      ym=0.0                                                            
      ys=sqrt(dot_product(w,y**2))                                      
      y=y/ys                                                            
      do 11721 j=1,ni                                                   
      if(ju(j).eq.0)goto 11721                                          
      xm(j)=0.0                                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                       
      if(isd .eq. 0)goto 11741                                          
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                         
      vc=xv(j)-xbq                                                      
      xs(j)=sqrt(vc)                                                    
      xv(j)=1.0+xbq/vc                                                  
      goto 11751                                                        
11741 continue                                                          
      xs(j)=1.0                                                         
11751 continue                                                          
      continue                                                          
11721 continue                                                          
      continue                                                          
      return                                                            
11711 continue                                                          
      do 11761 j=1,ni                                                   
      if(ju(j).eq.0)goto 11761                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2              
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                    
11761 continue                                                          
      continue                                                          
      if(isd .ne. 0)goto 11781                                          
      xs=1.0                                                            
      goto 11791                                                        
11781 continue                                                          
      xv=1.0                                                            
11791 continue                                                          
      continue                                                          
      ym=dot_product(w,y)                                               
      y=y-ym                                                            
      ys=sqrt(dot_product(w,y**2))                                      
      y=y/ys                                                            
      return                                                            
      end                                                               
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision y(no),w(no),x(*),vp(ni),ulam(nlam),cl(2,ni)      
      double precision ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),x
     *v(ni)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                       
      double precision, dimension (:), allocatable :: a,g               
      integer, dimension (:), allocatable :: mm,iy                      
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx,itrace)     
      allocate(a(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(g(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(iy(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      bta=beta                                                          
      omb=1.0-bta                                                       
      alm=0.0                                                           
      iy=0                                                              
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 11811                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
11811 continue                                                          
      rsq=0.0                                                           
      a=0.0                                                             
      mm=0                                                              
      o=0.0                                                             
      nlp=0                                                             
      nin=nlp                                                           
      iz=0                                                              
      mnl=min(mnlam,nlam)                                               
      do 11821 j=1,ni                                                   
      if(ju(j).eq.0)goto 11821                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j)) 
11821 continue                                                          
      continue                                                          
      do 11831 m=1,nlam                                                 
      if(itrace.ne.0) call setpb(m-1)                                   
      alm0=alm                                                          
      if(flmin .lt. 1.0)goto 11851                                      
      alm=ulam(m)                                                       
      goto 11841                                                        
11851 if(m .le. 2)goto 11861                                            
      alm=alm*alf                                                       
      goto 11841                                                        
11861 if(m .ne. 1)goto 11871                                            
      alm=big                                                           
      goto 11881                                                        
11871 continue                                                          
      alm0=0.0                                                          
      do 11891 j=1,ni                                                   
      if(ju(j).eq.0)goto 11891                                          
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                        
11891 continue                                                          
      continue                                                          
      alm0=alm0/max(bta,1.0d-3)                                         
      alm=alf*alm0                                                      
11881 continue                                                          
11841 continue                                                          
      dem=alm*omb                                                       
      ab=alm*bta                                                        
      rsq0=rsq                                                          
      jz=1                                                              
      tlam=bta*(2.0*alm-alm0)                                           
      do 11901 k=1,ni                                                   
      if(iy(k).eq.1)goto 11901                                          
      if(ju(k).eq.0)goto 11901                                          
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                    
11901 continue                                                          
      continue                                                          
      continue                                                          
11911 continue                                                          
      if(iz*jz.ne.0) go to 10360                                        
10880 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 11921 k=1,ni                                                   
      if(iy(k).eq.0)goto 11921                                          
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)        
      ak=a(k)                                                           
      u=gk+ak*xv(k)                                                     
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 11921                                          
      if(mm(k) .ne. 0)goto 11941                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 11922                                           
      mm(k)=nin                                                         
      ia(nin)=k                                                         
11941 continue                                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*gk-del*xv(k))                                    
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                      
      o=o+del*xm(k)/xs(k)                                               
      dlx=max(xv(k)*del**2,dlx)                                         
11921 continue                                                          
11922 continue                                                          
      if(nin.gt.nx)goto 11912                                           
      if(dlx .ge. thr)goto 11961                                        
      ixx=0                                                             
      do 11971 j=1,ni                                                   
      if(iy(j).eq.1)goto 11971                                          
      if(ju(j).eq.0)goto 11971                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j)) 
      if(g(j) .le. ab*vp(j))goto 11991                                  
      iy(j)=1                                                           
      ixx=1                                                             
11991 continue                                                          
11971 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10880                                          
      goto 11912                                                        
11961 continue                                                          
      if(nlp .le. maxit)goto 12011                                      
      jerr=-m                                                           
      return                                                            
12011 continue                                                          
10360 continue                                                          
      iz=1                                                              
      continue                                                          
12021 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 12031 l=1,nin                                                  
      k=ia(l)                                                           
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)        
      ak=a(k)                                                           
      u=gk+ak*xv(k)                                                     
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 12031                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*gk-del*xv(k))                                    
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                      
      o=o+del*xm(k)/xs(k)                                               
      dlx=max(xv(k)*del**2,dlx)                                         
12031 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 12022                                          
      if(nlp .le. maxit)goto 12051                                      
      jerr=-m                                                           
      return                                                            
12051 continue                                                          
      goto 12021                                                        
12022 continue                                                          
      jz=0                                                              
      goto 11911                                                        
11912 continue                                                          
      if(nin .le. nx)goto 12071                                         
      jerr=-10000-m                                                     
      goto 11832                                                        
12071 continue                                                          
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                             
      kin(m)=nin                                                        
      rsqo(m)=rsq                                                       
      almo(m)=alm                                                       
      lmu=m                                                             
      if(m.lt.mnl)goto 11831                                            
      if(flmin.ge.1.0)goto 11831                                        
      me=0                                                              
      do 12081 j=1,nin                                                  
      if(ao(j,m).ne.0.0) me=me+1                                        
12081 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 11832                                            
      if(rsq-rsq0.lt.sml*rsq)goto 11832                                 
      if(rsq.gt.rsqmax)goto 11832                                       
11831 continue                                                          
11832 continue                                                          
      deallocate(a,mm,g,iy)                                             
      return                                                            
      end                                                               
      subroutine spchkvars(no,ni,x,ix,ju)                               
      implicit double precision(a-h,o-z)                                
      double precision x(*)                                             
      integer ix(*),ju(ni)                                              
      do 12091 j=1,ni                                                   
      ju(j)=0                                                           
      jb=ix(j)                                                          
      nj=ix(j+1)-jb                                                     
      if(nj.eq.0)goto 12091                                             
      je=ix(j+1)-1                                                      
      if(nj .ge. no)goto 12111                                          
      do 12121 i=jb,je                                                  
      if(x(i).eq.0.0)goto 12121                                         
      ju(j)=1                                                           
      goto 12122                                                        
12121 continue                                                          
12122 continue                                                          
      goto 12131                                                        
12111 continue                                                          
      t=x(jb)                                                           
      do 12141 i=jb+1,je                                                
      if(x(i).eq.t)goto 12141                                           
      ju(j)=1                                                           
      goto 12142                                                        
12141 continue                                                          
12142 continue                                                          
12131 continue                                                          
      continue                                                          
12091 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                      
      implicit double precision(a-h,o-z)                                
      double precision ca(*),x(*),f(n)                                  
      integer ia(*),ix(*),jx(*)                                         
      f=a0                                                              
      do 12151 j=1,nin                                                  
      k=ia(j)                                                           
      kb=ix(k)                                                          
      ke=ix(k+1)-1                                                      
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                          
12151 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      function row_prod(i,j,ia,ja,ra,w)                                 
      implicit double precision(a-h,o-z)                                
      integer ia(*),ja(*)                                               
      double precision ra(*),w(*)                                       
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(
     *i),ia(j+1)-ia(j),w)
      return                                                            
      end                                                               
      function dot(x,y,mx,my,nx,ny,w)                                   
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(*),w(*)                                   
      integer mx(*),my(*)                                               
      i=1                                                               
      j=i                                                               
      s=0.0                                                             
      continue                                                          
12161 continue                                                          
      continue                                                          
12171 if(mx(i).ge.my(j))goto 12172                                      
      i=i+1                                                             
      if(i.gt.nx) go to 12180                                           
      goto 12171                                                        
12172 continue                                                          
      if(mx(i).eq.my(j)) go to 12190                                    
      continue                                                          
12201 if(my(j).ge.mx(i))goto 12202                                      
      j=j+1                                                             
      if(j.gt.ny) go to 12180                                           
      goto 12201                                                        
12202 continue                                                          
      if(mx(i).eq.my(j)) go to 12190                                    
      goto 12161                                                        
12190 continue                                                          
      s=s+w(mx(i))*x(i)*y(j)                                            
      i=i+1                                                             
      if(i.gt.nx)goto 12162                                             
      j=j+1                                                             
      if(j.gt.ny)goto 12162                                             
      goto 12161                                                        
12162 continue                                                          
12180 continue                                                          
      dot=s                                                             
      return                                                            
      end                                                               
      subroutine lognet(parm,no,ni,nc,x,y,g,jd,vp,cl,ne,nx,nlam,flmin,ul
     *am,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nla
     *m)
      double precision ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl
     *(2,ni)
      integer jd(*),ia(nx),nin(nlam)                                    
      double precision, dimension (:), allocatable :: xm,xs,ww,vq,xv    
      integer, dimension (:), allocatable :: ju                         
      if(maxval(vp) .gt. 0.0)goto 12221                                 
      jerr=10000                                                        
      return                                                            
12221 continue                                                          
      allocate(ww(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      if(kopt .ne. 2)goto 12241                                         
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
12241 continue                                                          
      if(isd .le. 0)goto 12261                                          
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
12261 continue                                                          
      call chkvars(no,ni,x,ju)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 12281                                   
      jerr=7777                                                         
      return                                                            
12281 continue                                                          
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      do 12291 i=1,no                                                   
      ww(i)=sum(y(i,:))                                                 
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                              
12291 continue                                                          
      continue                                                          
      sw=sum(ww)                                                        
      ww=ww/sw                                                          
      if(nc .ne. 1)goto 12311                                           
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                     
      if(isd .le. 0)goto 12331                                          
      do 12341 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
12341 continue                                                          
      continue                                                          
12331 continue                                                          
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,cl,ne,nx,nlam,fl
     *min,ulam,  thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,n
     *lp,jerr)
      goto 12301                                                        
12311 if(kopt .ne. 2)goto 12351                                         
      call multlstandard1(no,ni,x,ww,ju,isd,intr,xm,xs,xv)              
      if(isd .le. 0)goto 12371                                          
      do 12381 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
12381 continue                                                          
      continue                                                          
12371 continue                                                          
      call multlognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,
     *ulam,thr,  intr,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 12391                                                        
12351 continue                                                          
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                     
      if(isd .le. 0)goto 12411                                          
      do 12421 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
12421 continue                                                          
      continue                                                          
12411 continue                                                          
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam
     *,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
12391 continue                                                          
12301 continue                                                          
      if(jerr.gt.0) return                                              
      dev0=2.0*sw*dev0                                                  
      do 12431 k=1,lmu                                                  
      nk=nin(k)                                                         
      do 12441 ic=1,nc                                                  
      if(isd .le. 0)goto 12461                                          
      do 12471 l=1,nk                                                   
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                   
12471 continue                                                          
      continue                                                          
12461 continue                                                          
      if(intr .ne. 0)goto 12491                                         
      a0(ic,k)=0.0                                                      
      goto 12501                                                        
12491 continue                                                          
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))         
12501 continue                                                          
      continue                                                          
12441 continue                                                          
      continue                                                          
12431 continue                                                          
      continue                                                          
      deallocate(ww,ju,vq,xm)                                           
      if(isd.gt.0) deallocate(xs)                                       
      if(kopt.eq.2) deallocate(xv)                                      
      return                                                            
      end                                                               
      subroutine lstandard1(no,ni,x,w,ju,isd,intr,xm,xs)                
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),w(no),xm(ni),xs(ni)                     
      integer ju(ni)                                                    
      if(intr .ne. 0)goto 12521                                         
      do 12531 j=1,ni                                                   
      if(ju(j).eq.0)goto 12531                                          
      xm(j)=0.0                                                         
      if(isd .eq. 0)goto 12551                                          
      vc=dot_product(w,x(:,j)**2)-dot_product(w,x(:,j))**2              
      xs(j)=sqrt(vc)                                                    
      x(:,j)=x(:,j)/xs(j)                                               
12551 continue                                                          
12531 continue                                                          
      continue                                                          
      return                                                            
12521 continue                                                          
      do 12561 j=1,ni                                                   
      if(ju(j).eq.0)goto 12561                                          
      xm(j)=dot_product(w,x(:,j))                                       
      x(:,j)=x(:,j)-xm(j)                                               
      if(isd .le. 0)goto 12581                                          
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                              
      x(:,j)=x(:,j)/xs(j)                                               
12581 continue                                                          
12561 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine multlstandard1(no,ni,x,w,ju,isd,intr,xm,xs,xv)         
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),w(no),xm(ni),xs(ni),xv(ni)              
      integer ju(ni)                                                    
      if(intr .ne. 0)goto 12601                                         
      do 12611 j=1,ni                                                   
      if(ju(j).eq.0)goto 12611                                          
      xm(j)=0.0                                                         
      xv(j)=dot_product(w,x(:,j)**2)                                    
      if(isd .eq. 0)goto 12631                                          
      xbq=dot_product(w,x(:,j))**2                                      
      vc=xv(j)-xbq                                                      
      xs(j)=sqrt(vc)                                                    
      x(:,j)=x(:,j)/xs(j)                                               
      xv(j)=1.0+xbq/vc                                                  
12631 continue                                                          
12611 continue                                                          
      continue                                                          
      return                                                            
12601 continue                                                          
      do 12641 j=1,ni                                                   
      if(ju(j).eq.0)goto 12641                                          
      xm(j)=dot_product(w,x(:,j))                                       
      x(:,j)=x(:,j)-xm(j)                                               
      xv(j)=dot_product(w,x(:,j)**2)                                    
      if(isd .le. 0)goto 12661                                          
      xs(j)=sqrt(xv(j))                                                 
      x(:,j)=x(:,j)/xs(j)                                               
      xv(j)=1.0                                                         
12661 continue                                                          
12641 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin,u
     *lam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2
     *,ni)
      double precision a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)          
      integer ju(ni),m(nx),kin(nlam)                                    
      double precision, dimension (:), allocatable :: b,bs,v,r,xv,q,ga  
      integer, dimension (:), allocatable :: mm,ixx                     
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx,itrace)     
      allocate(b(0:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ga(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(bs(0:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ixx(1:ni),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(r(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(v(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(q(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      fmax=log(1.0/pmin-1.0)                                            
      fmin=-fmax                                                        
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                   
      bta=parm                                                          
      omb=1.0-bta                                                       
      q0=dot_product(w,y)                                               
      if(q0 .gt. pmin)goto 12681                                        
      jerr=8001                                                         
      return                                                            
12681 continue                                                          
      if(q0 .lt. 1.0-pmin)goto 12701                                    
      jerr=9001                                                         
      return                                                            
12701 continue                                                          
      if(intr.eq.0.0) q0=0.5                                            
      ixx=0                                                             
      al=0.0                                                            
      bz=0.0                                                            
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                 
      if(nonzero(no,g) .ne. 0)goto 12721                                
      vi=q0*(1.0-q0)                                                    
      b(0)=bz                                                           
      v=vi*w                                                            
      r=w*(y-q0)                                                        
      q=q0                                                              
      xmz=vi                                                            
      dev1=-(bz*q0+log(1.0-q0))                                         
      goto 12731                                                        
12721 continue                                                          
      b(0)=0.0                                                          
      if(intr .eq. 0)goto 12751                                         
      b(0)=azero(no,y,g,w,jerr)                                         
      if(jerr.ne.0) return                                              
12751 continue                                                          
      q=1.0/(1.0+exp(-b(0)-g))                                          
      v=w*q*(1.0-q)                                                     
      r=w*(y-q)                                                         
      xmz=sum(v)                                                        
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                     
12731 continue                                                          
      continue                                                          
      if(kopt .le. 0)goto 12771                                         
      if(isd .le. 0 .or. intr .eq. 0)goto 12791                         
      xv=0.25                                                           
      goto 12801                                                        
12791 continue                                                          
      do 12811 j=1,ni                                                   
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                
12811 continue                                                          
      continue                                                          
12801 continue                                                          
      continue                                                          
12771 continue                                                          
      dev0=dev1                                                         
      do 12821 i=1,no                                                   
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                     
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))           
12821 continue                                                          
      continue                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 12841                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
12841 continue                                                          
      m=0                                                               
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      mnl=min(mnlam,nlam)                                               
      bs=0.0                                                            
      b(1:ni)=0.0                                                       
      shr=shri*dev0                                                     
      do 12851 j=1,ni                                                   
      if(ju(j).eq.0)goto 12851                                          
      ga(j)=abs(dot_product(r,x(:,j)))                                  
12851 continue                                                          
      continue                                                          
      do 12861 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 12881                                      
      al=ulam(ilm)                                                      
      goto 12871                                                        
12881 if(ilm .le. 2)goto 12891                                          
      al=al*alf                                                         
      goto 12871                                                        
12891 if(ilm .ne. 1)goto 12901                                          
      al=big                                                            
      goto 12911                                                        
12901 continue                                                          
      al0=0.0                                                           
      do 12921 j=1,ni                                                   
      if(ju(j).eq.0)goto 12921                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
12921 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
12911 continue                                                          
12871 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 12931 k=1,ni                                                   
      if(ixx(k).eq.1)goto 12931                                         
      if(ju(k).eq.0)goto 12931                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
12931 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
12941 continue                                                          
      bs(0)=b(0)                                                        
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                             
      if(kopt .ne. 0)goto 12961                                         
      do 12971 j=1,ni                                                   
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                    
12971 continue                                                          
      continue                                                          
12961 continue                                                          
      continue                                                          
12981 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 12991 k=1,ni                                                   
      if(ixx(k).eq.0)goto 12991                                         
      bk=b(k)                                                           
      gk=dot_product(r,x(:,k))                                          
      u=gk+xv(k)*b(k)                                                   
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 13011                                         
      b(k)=0.0                                                          
      goto 13021                                                        
13011 continue                                                          
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))       
13021 continue                                                          
      continue                                                          
      d=b(k)-bk                                                         
      if(abs(d).le.0.0)goto 12991                                       
      dlx=max(dlx,xv(k)*d**2)                                           
      r=r-d*v*x(:,k)                                                    
      if(mm(k) .ne. 0)goto 13041                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 12992                                           
      mm(k)=nin                                                         
      m(nin)=k                                                          
13041 continue                                                          
12991 continue                                                          
12992 continue                                                          
      if(nin.gt.nx)goto 12982                                           
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 13061                                          
      b(0)=b(0)+d                                                       
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
13061 continue                                                          
      if(dlx.lt.shr)goto 12982                                          
      if(nlp .le. maxit)goto 13081                                      
      jerr=-ilm                                                         
      return                                                            
13081 continue                                                          
      continue                                                          
13091 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 13101 l=1,nin                                                  
      k=m(l)                                                            
      bk=b(k)                                                           
      gk=dot_product(r,x(:,k))                                          
      u=gk+xv(k)*b(k)                                                   
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 13121                                         
      b(k)=0.0                                                          
      goto 13131                                                        
13121 continue                                                          
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))       
13131 continue                                                          
      continue                                                          
      d=b(k)-bk                                                         
      if(abs(d).le.0.0)goto 13101                                       
      dlx=max(dlx,xv(k)*d**2)                                           
      r=r-d*v*x(:,k)                                                    
13101 continue                                                          
      continue                                                          
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 13151                                          
      b(0)=b(0)+d                                                       
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
13151 continue                                                          
      if(dlx.lt.shr)goto 13092                                          
      if(nlp .le. maxit)goto 13171                                      
      jerr=-ilm                                                         
      return                                                            
13171 continue                                                          
      goto 13091                                                        
13092 continue                                                          
      goto 12981                                                        
12982 continue                                                          
      if(nin.gt.nx)goto 12942                                           
      do 13181 i=1,no                                                   
      fi=b(0)+g(i)                                                      
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))         
      if(fi .ge. fmin)goto 13201                                        
      q(i)=0.0                                                          
      goto 13191                                                        
13201 if(fi .le. fmax)goto 13211                                        
      q(i)=1.0                                                          
      goto 13221                                                        
13211 continue                                                          
      q(i)=1.0/(1.0+exp(-fi))                                           
13221 continue                                                          
13191 continue                                                          
13181 continue                                                          
      continue                                                          
      v=w*q*(1.0-q)                                                     
      xmz=sum(v)                                                        
      if(xmz.le.vmin)goto 12942                                         
      r=w*(y-q)                                                         
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 13241                        
      ix=0                                                              
      do 13251 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 13251                        
      ix=1                                                              
      goto 13252                                                        
13251 continue                                                          
13252 continue                                                          
      if(ix .ne. 0)goto 13271                                           
      do 13281 k=1,ni                                                   
      if(ixx(k).eq.1)goto 13281                                         
      if(ju(k).eq.0)goto 13281                                          
      ga(k)=abs(dot_product(r,x(:,k)))                                  
      if(ga(k) .le. al1*vp(k))goto 13301                                
      ixx(k)=1                                                          
      ix=1                                                              
13301 continue                                                          
13281 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10880                                           
      goto 12942                                                        
13271 continue                                                          
13241 continue                                                          
      goto 12941                                                        
12942 continue                                                          
      if(nin .le. nx)goto 13321                                         
      jerr=-10000-ilm                                                   
      goto 12862                                                        
13321 continue                                                          
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                             
      kin(ilm)=nin                                                      
      a0(ilm)=b(0)                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      devi=dev2(no,w,y,q,pmin)                                          
      dev(ilm)=(dev1-devi)/dev0                                         
      if(xmz.le.vmin)goto 12862                                         
      if(ilm.lt.mnl)goto 12861                                          
      if(flmin.ge.1.0)goto 12861                                        
      me=0                                                              
      do 13331 j=1,nin                                                  
      if(a(j,ilm).ne.0.0) me=me+1                                       
13331 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 12862                                            
      if(dev(ilm).gt.devmax)goto 12862                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12862                          
12861 continue                                                          
12862 continue                                                          
      g=log(q/(1.0-q))                                                  
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                               
      return                                                            
      end                                                               
      function dev2(n,w,y,p,pmin)                                       
      implicit double precision(a-h,o-z)                                
      double precision w(n),y(n),p(n)                                   
      pmax=1.0-pmin                                                     
      s=0.0                                                             
      do 13341 i=1,n                                                    
      pi=min(max(pmin,p(i)),pmax)                                       
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                    
13341 continue                                                          
      continue                                                          
      dev2=s                                                            
      return                                                            
      end                                                               
      function azero(n,y,g,q,jerr)                                      
      implicit double precision(a-h,o-z)                                
      parameter(eps=1.0d-7)                                             
      double precision y(n),g(n),q(n)                                   
      double precision, dimension (:), allocatable :: e,p,w             
      azero = 0.0                                                       
      allocate(e(1:n),stat=jerr)                                        
      if(jerr.ne.0) return                                              
      allocate(p(1:n),stat=jerr)                                        
      if(jerr.ne.0) return                                              
      allocate(w(1:n),stat=jerr)                                        
      if(jerr.ne.0) return                                              
      az=0.0                                                            
      e=exp(-g)                                                         
      qy=dot_product(q,y)                                               
      p=1.0/(1.0+e)                                                     
      continue                                                          
13351 continue                                                          
      w=q*p*(1.0-p)                                                     
      d=(qy-dot_product(q,p))/sum(w)                                    
      az=az+d                                                           
      if(abs(d).lt.eps)goto 13352                                       
      ea0=exp(-az)                                                      
      p=1.0/(1.0+ea0*e)                                                 
      goto 13351                                                        
13352 continue                                                          
      azero=az                                                          
      deallocate(e,p,w)                                                 
      return                                                            
      end                                                               
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin
     *,ulam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,j
     *err)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam
     *)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(
     *2,ni)
      integer ju(ni),m(nx),kin(nlam)                                    
      double precision, dimension (:,:), allocatable :: q               
      double precision, dimension (:), allocatable :: sxp,sxpl          
      double precision, dimension (:), allocatable :: di,v,r,ga         
      double precision, dimension (:,:), allocatable :: b,bs,xv         
      integer, dimension (:), allocatable :: mm,is,ixx                  
      allocate(b(0:ni,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni,1:nc),stat=jerr)                                 
      if(jerr.ne.0) return                                              
      allocate(bs(0:ni,1:nc),stat=jerr)                                 
      if(jerr.ne.0) return                                              
      allocate(q(1:no,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx,itrace)     
      exmn=-exmx                                                        
      allocate(r(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(v(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(is(1:max(nc,ni)),stat=jerr)                              
      if(jerr.ne.0) return                                              
      allocate(sxp(1:no),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(sxpl(1:no),stat=jerr)                                    
      if(jerr.ne.0) return                                              
      allocate(di(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ga(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ixx(1:ni),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      pmax=1.0-pmin                                                     
      emin=pmin/pmax                                                    
      emax=1.0/emin                                                     
      pfm=(1.0+pmin)*pmin                                               
      pfx=(1.0-pmin)*pmax                                               
      vmin=pfm*pmax                                                     
      bta=parm                                                          
      omb=1.0-bta                                                       
      dev1=0.0                                                          
      dev0=0.0                                                          
      do 13361 ic=1,nc                                                  
      q0=dot_product(w,y(:,ic))                                         
      if(q0 .gt. pmin)goto 13381                                        
      jerr =8000+ic                                                     
      return                                                            
13381 continue                                                          
      if(q0 .lt. 1.0-pmin)goto 13401                                    
      jerr =9000+ic                                                     
      return                                                            
13401 continue                                                          
      if(intr .ne. 0)goto 13421                                         
      q0=1.0/nc                                                         
      b(0,ic)=0.0                                                       
      goto 13431                                                        
13421 continue                                                          
      b(0,ic)=log(q0)                                                   
      dev1=dev1-q0*b(0,ic)                                              
13431 continue                                                          
      continue                                                          
      b(1:ni,ic)=0.0                                                    
13361 continue                                                          
      continue                                                          
      if(intr.eq.0) dev1=log(float(nc))                                 
      ixx=0                                                             
      al=0.0                                                            
      if(nonzero(no*nc,g) .ne. 0)goto 13451                             
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                      
      sxp=0.0                                                           
      do 13461 ic=1,nc                                                  
      q(:,ic)=exp(b(0,ic))                                              
      sxp=sxp+q(:,ic)                                                   
13461 continue                                                          
      continue                                                          
      goto 13471                                                        
13451 continue                                                          
      do 13481 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
13481 continue                                                          
      continue                                                          
      sxp=0.0                                                           
      if(intr .ne. 0)goto 13501                                         
      b(0,:)=0.0                                                        
      goto 13511                                                        
13501 continue                                                          
      call kazero(nc,no,y,g,w,b(0,:),jerr)                              
      if(jerr.ne.0) return                                              
13511 continue                                                          
      continue                                                          
      dev1=0.0                                                          
      do 13521 ic=1,nc                                                  
      q(:,ic)=b(0,ic)+g(:,ic)                                           
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                          
      q(:,ic)=exp(q(:,ic))                                              
      sxp=sxp+q(:,ic)                                                   
13521 continue                                                          
      continue                                                          
      sxpl=w*log(sxp)                                                   
      do 13531 ic=1,nc                                                  
      dev1=dev1+dot_product(y(:,ic),sxpl)                               
13531 continue                                                          
      continue                                                          
13471 continue                                                          
      continue                                                          
      do 13541 ic=1,nc                                                  
      do 13551 i=1,no                                                   
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))            
13551 continue                                                          
      continue                                                          
13541 continue                                                          
      continue                                                          
      dev0=dev0+dev1                                                    
      if(kopt .le. 0)goto 13571                                         
      if(isd .le. 0 .or. intr .eq. 0)goto 13591                         
      xv=0.25                                                           
      goto 13601                                                        
13591 continue                                                          
      do 13611 j=1,ni                                                   
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)              
13611 continue                                                          
      continue                                                          
13601 continue                                                          
      continue                                                          
13571 continue                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 13631                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
13631 continue                                                          
      m=0                                                               
      mm=0                                                              
      nin=0                                                             
      nlp=0                                                             
      mnl=min(mnlam,nlam)                                               
      bs=0.0                                                            
      shr=shri*dev0                                                     
      ga=0.0                                                            
      do 13641 ic=1,nc                                                  
      r=w*(y(:,ic)-q(:,ic)/sxp)                                         
      do 13651 j=1,ni                                                   
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))        
13651 continue                                                          
      continue                                                          
13641 continue                                                          
      continue                                                          
      do 13661 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 13681                                      
      al=ulam(ilm)                                                      
      goto 13671                                                        
13681 if(ilm .le. 2)goto 13691                                          
      al=al*alf                                                         
      goto 13671                                                        
13691 if(ilm .ne. 1)goto 13701                                          
      al=big                                                            
      goto 13711                                                        
13701 continue                                                          
      al0=0.0                                                           
      do 13721 j=1,ni                                                   
      if(ju(j).eq.0)goto 13721                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
13721 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
13711 continue                                                          
13671 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 13731 k=1,ni                                                   
      if(ixx(k).eq.1)goto 13731                                         
      if(ju(k).eq.0)goto 13731                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
13731 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
13741 continue                                                          
      ix=0                                                              
      jx=ix                                                             
      ig=0                                                              
      do 13751 ic=1,nc                                                  
      bs(0,ic)=b(0,ic)                                                  
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                       
      xmz=0.0                                                           
      do 13761 i=1,no                                                   
      pic=q(i,ic)/sxp(i)                                                
      if(pic .ge. pfm)goto 13781                                        
      pic=0.0                                                           
      v(i)=0.0                                                          
      goto 13771                                                        
13781 if(pic .le. pfx)goto 13791                                        
      pic=1.0                                                           
      v(i)=0.0                                                          
      goto 13801                                                        
13791 continue                                                          
      v(i)=w(i)*pic*(1.0-pic)                                           
      xmz=xmz+v(i)                                                      
13801 continue                                                          
13771 continue                                                          
      r(i)=w(i)*(y(i,ic)-pic)                                           
13761 continue                                                          
      continue                                                          
      if(xmz.le.vmin)goto 13751                                         
      ig=1                                                              
      if(kopt .ne. 0)goto 13821                                         
      do 13831 j=1,ni                                                   
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                 
13831 continue                                                          
      continue                                                          
13821 continue                                                          
      continue                                                          
13841 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 13851 k=1,ni                                                   
      if(ixx(k).eq.0)goto 13851                                         
      bk=b(k,ic)                                                        
      gk=dot_product(r,x(:,k))                                          
      u=gk+xv(k,ic)*b(k,ic)                                             
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 13871                                         
      b(k,ic)=0.0                                                       
      goto 13881                                                        
13871 continue                                                          
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))
     *)
13881 continue                                                          
      continue                                                          
      d=b(k,ic)-bk                                                      
      if(abs(d).le.0.0)goto 13851                                       
      dlx=max(dlx,xv(k,ic)*d**2)                                        
      r=r-d*v*x(:,k)                                                    
      if(mm(k) .ne. 0)goto 13901                                        
      nin=nin+1                                                         
      if(nin .le. nx)goto 13921                                         
      jx=1                                                              
      goto 13852                                                        
13921 continue                                                          
      mm(k)=nin                                                         
      m(nin)=k                                                          
13901 continue                                                          
13851 continue                                                          
13852 continue                                                          
      if(jx.gt.0)goto 13842                                             
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 13941                                          
      b(0,ic)=b(0,ic)+d                                                 
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
13941 continue                                                          
      if(dlx.lt.shr)goto 13842                                          
      if(nlp .le. maxit)goto 13961                                      
      jerr=-ilm                                                         
      return                                                            
13961 continue                                                          
      continue                                                          
13971 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 13981 l=1,nin                                                  
      k=m(l)                                                            
      bk=b(k,ic)                                                        
      gk=dot_product(r,x(:,k))                                          
      u=gk+xv(k,ic)*b(k,ic)                                             
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 14001                                         
      b(k,ic)=0.0                                                       
      goto 14011                                                        
14001 continue                                                          
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))
     *)
14011 continue                                                          
      continue                                                          
      d=b(k,ic)-bk                                                      
      if(abs(d).le.0.0)goto 13981                                       
      dlx=max(dlx,xv(k,ic)*d**2)                                        
      r=r-d*v*x(:,k)                                                    
13981 continue                                                          
      continue                                                          
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 14031                                          
      b(0,ic)=b(0,ic)+d                                                 
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
14031 continue                                                          
      if(dlx.lt.shr)goto 13972                                          
      if(nlp .le. maxit)goto 14051                                      
      jerr=-ilm                                                         
      return                                                            
14051 continue                                                          
      goto 13971                                                        
13972 continue                                                          
      goto 13841                                                        
13842 continue                                                          
      if(jx.gt.0)goto 13752                                             
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                         
      if(ix .ne. 0)goto 14071                                           
      do 14081 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 14101             
      ix=1                                                              
      goto 14082                                                        
14101 continue                                                          
14081 continue                                                          
14082 continue                                                          
14071 continue                                                          
      do 14111 i=1,no                                                   
      fi=b(0,ic)+g(i,ic)                                                
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))      
      fi=min(max(exmn,fi),exmx)                                         
      sxp(i)=sxp(i)-q(i,ic)                                             
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                 
      sxp(i)=sxp(i)+q(i,ic)                                             
14111 continue                                                          
      continue                                                          
13751 continue                                                          
13752 continue                                                          
      s=-sum(b(0,:))/nc                                                 
      b(0,:)=b(0,:)+s                                                   
      di=s                                                              
      do 14121 j=1,nin                                                  
      l=m(j)                                                            
      if(vp(l) .gt. 0.0)goto 14141                                      
      s=sum(b(l,:))/nc                                                  
      goto 14151                                                        
14141 continue                                                          
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                  
14151 continue                                                          
      continue                                                          
      b(l,:)=b(l,:)-s                                                   
      di=di-s*x(:,l)                                                    
14121 continue                                                          
      continue                                                          
      di=exp(di)                                                        
      sxp=sxp*di                                                        
      do 14161 ic=1,nc                                                  
      q(:,ic)=q(:,ic)*di                                                
14161 continue                                                          
      continue                                                          
      if(jx.gt.0)goto 13742                                             
      if(ig.eq.0)goto 13742                                             
      if(ix .ne. 0)goto 14181                                           
      do 14191 k=1,ni                                                   
      if(ixx(k).eq.1)goto 14191                                         
      if(ju(k).eq.0)goto 14191                                          
      ga(k)=0.0                                                         
14191 continue                                                          
      continue                                                          
      do 14201 ic=1,nc                                                  
      r=w*(y(:,ic)-q(:,ic)/sxp)                                         
      do 14211 k=1,ni                                                   
      if(ixx(k).eq.1)goto 14211                                         
      if(ju(k).eq.0)goto 14211                                          
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                       
14211 continue                                                          
      continue                                                          
14201 continue                                                          
      continue                                                          
      do 14221 k=1,ni                                                   
      if(ixx(k).eq.1)goto 14221                                         
      if(ju(k).eq.0)goto 14221                                          
      if(ga(k) .le. al1*vp(k))goto 14241                                
      ixx(k)=1                                                          
      ix=1                                                              
14241 continue                                                          
14221 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10880                                           
      goto 13742                                                        
14181 continue                                                          
      goto 13741                                                        
13742 continue                                                          
      if(jx .le. 0)goto 14261                                           
      jerr=-10000-ilm                                                   
      goto 13662                                                        
14261 continue                                                          
      devi=0.0                                                          
      do 14271 ic=1,nc                                                  
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                       
      a0(ic,ilm)=b(0,ic)                                                
      do 14281 i=1,no                                                   
      if(y(i,ic).le.0.0)goto 14281                                      
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                        
14281 continue                                                          
      continue                                                          
14271 continue                                                          
      continue                                                          
      kin(ilm)=nin                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dev1-devi)/dev0                                         
      if(ig.eq.0)goto 13662                                             
      if(ilm.lt.mnl)goto 13661                                          
      if(flmin.ge.1.0)goto 13661                                        
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 13662          
      if(dev(ilm).gt.devmax)goto 13662                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13662                          
13661 continue                                                          
13662 continue                                                          
      g=log(q)                                                          
      do 14291 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
14291 continue                                                          
      continue                                                          
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                        
      return                                                            
      end                                                               
      subroutine kazero(kk,n,y,g,q,az,jerr)                             
      implicit double precision(a-h,o-z)                                
      parameter(eps=1.0d-7)                                             
      double precision y(n,kk),g(n,kk),q(n),az(kk)                      
      double precision, dimension (:), allocatable :: s                 
      double precision, dimension (:,:), allocatable :: e               
      allocate(e(1:n,1:kk),stat=jerr)                                   
      if(jerr.ne.0) return                                              
      allocate(s(1:n),stat=jerr)                                        
      if(jerr.ne.0) return                                              
      az=0.0                                                            
      e=exp(g)                                                          
      do 14301 i=1,n                                                    
      s(i)=sum(e(i,:))                                                  
14301 continue                                                          
      continue                                                          
      continue                                                          
14311 continue                                                          
      dm=0.0                                                            
      do 14321 k=1,kk                                                   
      t=0.0                                                             
      u=t                                                               
      do 14331 i=1,n                                                    
      pik=e(i,k)/s(i)                                                   
      t=t+q(i)*(y(i,k)-pik)                                             
      u=u+q(i)*pik*(1.0-pik)                                            
14331 continue                                                          
      continue                                                          
      d=t/u                                                             
      az(k)=az(k)+d                                                     
      ed=exp(d)                                                         
      dm=max(dm,abs(d))                                                 
      do 14341 i=1,n                                                    
      z=e(i,k)                                                          
      e(i,k)=z*ed                                                       
      s(i)=s(i)-z+e(i,k)                                                
14341 continue                                                          
      continue                                                          
14321 continue                                                          
      continue                                                          
      if(dm.lt.eps)goto 14312                                           
      goto 14311                                                        
14312 continue                                                          
      az=az-sum(az)/kk                                                  
      deallocate(e,s)                                                   
      return                                                            
      end                                                               
      function elc(parm,n,cl,a,m)                                       
      implicit double precision(a-h,o-z)                                
      double precision a(n),cl(2)                                       
      integer m(n)                                                      
      fn=n                                                              
      am=sum(a)/fn                                                      
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 14361                    
      elc=am                                                            
      go to 14370                                                       
14361 continue                                                          
      do 14381 i=1,n                                                    
      m(i)=i                                                            
14381 continue                                                          
      continue                                                          
      call psort7(a,m,1,n)                                              
      if(a(m(1)) .ne. a(m(n)))goto 14401                                
      elc=a(1)                                                          
      go to 14370                                                       
14401 continue                                                          
      if(mod(n,2) .ne. 1)goto 14421                                     
      ad=a(m(n/2+1))                                                    
      goto 14431                                                        
14421 continue                                                          
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                    
14431 continue                                                          
      continue                                                          
      if(parm .ne. 1.0)goto 14451                                       
      elc=ad                                                            
      go to 14370                                                       
14451 continue                                                          
      b1=min(am,ad)                                                     
      b2=max(am,ad)                                                     
      k2=1                                                              
      continue                                                          
14461 if(a(m(k2)).gt.b1)goto 14462                                      
      k2=k2+1                                                           
      goto 14461                                                        
14462 continue                                                          
      k1=k2-1                                                           
      continue                                                          
14471 if(a(m(k2)).ge.b2)goto 14472                                      
      k2=k2+1                                                           
      goto 14471                                                        
14472 continue                                                          
      r=parm/((1.0-parm)*fn)                                            
      is=0                                                              
      sm=n-2*(k1-1)                                                     
      do 14481 k=k1,k2-1                                                
      sm=sm-2.0                                                         
      s=r*sm+am                                                         
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 14501                
      is=k                                                              
      goto 14482                                                        
14501 continue                                                          
14481 continue                                                          
14482 continue                                                          
      if(is .eq. 0)goto 14521                                           
      elc=s                                                             
      go to 14370                                                       
14521 continue                                                          
      r2=2.0*r                                                          
      s1=a(m(k1))                                                       
      am2=2.0*am                                                        
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                 
      elc=s1                                                            
      do 14531 k=k1+1,k2                                                
      s=a(m(k))                                                         
      if(s.eq.s1)goto 14531                                             
      c=r2*sum(abs(a-s))+s*(s-am2)                                      
      if(c .ge. cri)goto 14551                                          
      cri=c                                                             
      elc=s                                                             
14551 continue                                                          
      s1=s                                                              
14531 continue                                                          
      continue                                                          
14370 continue                                                          
      elc=max(maxval(a-cl(2)),min(minval(a-cl(1)),elc))                 
      return                                                            
      end                                                               
      function nintot(ni,nx,nc,a,m,nin,is)                              
      implicit double precision(a-h,o-z)                                
      double precision a(nx,nc)                                         
      integer m(nx),is(ni)                                              
      is=0                                                              
      nintot=0                                                          
      do 14561 ic=1,nc                                                  
      do 14571 j=1,nin                                                  
      k=m(j)                                                            
      if(is(k).ne.0)goto 14571                                          
      if(a(j,ic).eq.0.0)goto 14571                                      
      is(k)=k                                                           
      nintot=nintot+1                                                   
14571 continue                                                          
      continue                                                          
14561 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                          
      implicit double precision(a-h,o-z)                                
      double precision ca(nx,nc),a(ni,nc)                               
      integer ia(nx)                                                    
      a=0.0                                                             
      do 14581 ic=1,nc                                                  
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                         
14581 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                   
      implicit double precision(a-h,o-z)                                
      double precision a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)              
      integer ia(nx)                                                    
      do 14591 i=1,nt                                                   
      do 14601 ic=1,nc                                                  
      ans(ic,i)=a0(ic)                                                  
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1
     *:nin)))
14601 continue                                                          
      continue                                                          
14591 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine splognet(parm,no,ni,nc,x,ix,jx,y,g,jd,vp,cl,ne,nx,nlam,
     *flmin,  ulam,thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm
     *,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)  
      double precision ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl
     *(2,ni)
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                        
      double precision, dimension (:), allocatable :: xm,xs,ww,vq,xv    
      integer, dimension (:), allocatable :: ju                         
      if(maxval(vp) .gt. 0.0)goto 14621                                 
      jerr=10000                                                        
      return                                                            
14621 continue                                                          
      allocate(ww(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      if(kopt .ne. 2)goto 14641                                         
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
14641 continue                                                          
      call spchkvars(no,ni,x,ix,ju)                                     
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 14661                                   
      jerr=7777                                                         
      return                                                            
14661 continue                                                          
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      do 14671 i=1,no                                                   
      ww(i)=sum(y(i,:))                                                 
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                              
14671 continue                                                          
      continue                                                          
      sw=sum(ww)                                                        
      ww=ww/sw                                                          
      if(nc .ne. 1)goto 14691                                           
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)             
      if(isd .le. 0)goto 14711                                          
      do 14721 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
14721 continue                                                          
      continue                                                          
14711 continue                                                          
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,cl,ne,n
     *x,nlam,  flmin,ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin
     *,dev0,dev,  alm,nlp,jerr)
      goto 14681                                                        
14691 if(kopt .ne. 2)goto 14731                                         
      call multsplstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs,xv)      
      if(isd .le. 0)goto 14751                                          
      do 14761 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
14761 continue                                                          
      continue                                                          
14751 continue                                                          
      call multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nl
     *am,flmin,  ulam,thr,intr,maxit,xv,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,
     *alm,nlp,jerr)
      goto 14771                                                        
14731 continue                                                          
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)             
      if(isd .le. 0)goto 14791                                          
      do 14801 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
14801 continue                                                          
      continue                                                          
14791 continue                                                          
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,f
     *lmin,  ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,  ia,nin,dev0,
     *dev,alm,nlp,jerr)
14771 continue                                                          
14681 continue                                                          
      if(jerr.gt.0) return                                              
      dev0=2.0*sw*dev0                                                  
      do 14811 k=1,lmu                                                  
      nk=nin(k)                                                         
      do 14821 ic=1,nc                                                  
      if(isd .le. 0)goto 14841                                          
      do 14851 l=1,nk                                                   
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                   
14851 continue                                                          
      continue                                                          
14841 continue                                                          
      if(intr .ne. 0)goto 14871                                         
      a0(ic,k)=0.0                                                      
      goto 14881                                                        
14871 continue                                                          
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))         
14881 continue                                                          
      continue                                                          
14821 continue                                                          
      continue                                                          
14811 continue                                                          
      continue                                                          
      deallocate(ww,ju,vq,xm,xs)                                        
      if(kopt.eq.2) deallocate(xv)                                      
      return                                                            
      end                                                               
      subroutine multsplstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs,xv) 
      implicit double precision(a-h,o-z)                                
      double precision x(*),w(no),xm(ni),xs(ni),xv(ni)                  
      integer ix(*),jx(*),ju(ni)                                        
      if(intr .ne. 0)goto 14901                                         
      do 14911 j=1,ni                                                   
      if(ju(j).eq.0)goto 14911                                          
      xm(j)=0.0                                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                       
      if(isd .eq. 0)goto 14931                                          
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                         
      vc=xv(j)-xbq                                                      
      xs(j)=sqrt(vc)                                                    
      xv(j)=1.0+xbq/vc                                                  
      goto 14941                                                        
14931 continue                                                          
      xs(j)=1.0                                                         
14941 continue                                                          
      continue                                                          
14911 continue                                                          
      continue                                                          
      return                                                            
14901 continue                                                          
      do 14951 j=1,ni                                                   
      if(ju(j).eq.0)goto 14951                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2              
      if(isd .le. 0)goto 14971                                          
      xs(j)=sqrt(xv(j))                                                 
      xv(j)=1.0                                                         
14971 continue                                                          
14951 continue                                                          
      continue                                                          
      if(isd.eq.0) xs=1.0                                               
      return                                                            
      end                                                               
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs)        
      implicit double precision(a-h,o-z)                                
      double precision x(*),w(no),xm(ni),xs(ni)                         
      integer ix(*),jx(*),ju(ni)                                        
      if(intr .ne. 0)goto 14991                                         
      do 15001 j=1,ni                                                   
      if(ju(j).eq.0)goto 15001                                          
      xm(j)=0.0                                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      if(isd .eq. 0)goto 15021                                          
      vc=dot_product(w(jx(jb:je)),x(jb:je)**2)  -dot_product(w(jx(jb:je)
     *),x(jb:je))**2
      xs(j)=sqrt(vc)                                                    
      goto 15031                                                        
15021 continue                                                          
      xs(j)=1.0                                                         
15031 continue                                                          
      continue                                                          
15001 continue                                                          
      continue                                                          
      return                                                            
14991 continue                                                          
      do 15041 j=1,ni                                                   
      if(ju(j).eq.0)goto 15041                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j
     *)**2)
15041 continue                                                          
      continue                                                          
      if(isd.eq.0) xs=1.0                                               
      return                                                            
      end                                                               
      subroutine sprlognet2n(parm,no,ni,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,nla
     *m,  flmin,ulam,shri,isd,intr,maxit,kopt,xb,xs,  lmu,a0,a,m,kin,dev
     *0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)
      double precision a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)          
      double precision xb(ni),xs(ni)                                    
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                        
      double precision, dimension (:), allocatable :: xm,b,bs,v,r       
      double precision, dimension (:), allocatable :: sc,xv,q,ga        
      integer, dimension (:), allocatable :: mm,ixx                     
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx,itrace)     
      allocate(b(0:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(xm(0:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(bs(0:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ga(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ixx(1:ni),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(q(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(r(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(v(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(sc(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      fmax=log(1.0/pmin-1.0)                                            
      fmin=-fmax                                                        
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                   
      bta=parm                                                          
      omb=1.0-bta                                                       
      q0=dot_product(w,y)                                               
      if(q0 .gt. pmin)goto 15061                                        
      jerr=8001                                                         
      return                                                            
15061 continue                                                          
      if(q0 .lt. 1.0-pmin)goto 15081                                    
      jerr=9001                                                         
      return                                                            
15081 continue                                                          
      if(intr.eq.0) q0=0.5                                              
      bz=0.0                                                            
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                 
      if(nonzero(no,g) .ne. 0)goto 15101                                
      vi=q0*(1.0-q0)                                                    
      b(0)=bz                                                           
      v=vi*w                                                            
      r=w*(y-q0)                                                        
      q=q0                                                              
      xm(0)=vi                                                          
      dev1=-(bz*q0+log(1.0-q0))                                         
      goto 15111                                                        
15101 continue                                                          
      b(0)=0.0                                                          
      if(intr .eq. 0)goto 15131                                         
      b(0)=azero(no,y,g,w,jerr)                                         
      if(jerr.ne.0) return                                              
15131 continue                                                          
      q=1.0/(1.0+exp(-b(0)-g))                                          
      v=w*q*(1.0-q)                                                     
      r=w*(y-q)                                                         
      xm(0)=sum(v)                                                      
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                     
15111 continue                                                          
      continue                                                          
      if(kopt .le. 0)goto 15151                                         
      if(isd .le. 0 .or. intr .eq. 0)goto 15171                         
      xv=0.25                                                           
      goto 15181                                                        
15171 continue                                                          
      do 15191 j=1,ni                                                   
      if(ju(j).eq.0)goto 15191                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)       
15191 continue                                                          
      continue                                                          
15181 continue                                                          
      continue                                                          
15151 continue                                                          
      b(1:ni)=0.0                                                       
      dev0=dev1                                                         
      do 15201 i=1,no                                                   
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                     
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))           
15201 continue                                                          
      continue                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 15221                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
15221 continue                                                          
      m=0                                                               
      mm=0                                                              
      nin=0                                                             
      o=0.0                                                             
      svr=o                                                             
      mnl=min(mnlam,nlam)                                               
      bs=0.0                                                            
      nlp=0                                                             
      nin=nlp                                                           
      shr=shri*dev0                                                     
      al=0.0                                                            
      ixx=0                                                             
      do 15231 j=1,ni                                                   
      if(ju(j).eq.0)goto 15231                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      jn=ix(j+1)-ix(j)                                                  
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                              
      gj=dot_product(sc(1:jn),x(jb:je))                                 
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                   
15231 continue                                                          
      continue                                                          
      do 15241 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 15261                                      
      al=ulam(ilm)                                                      
      goto 15251                                                        
15261 if(ilm .le. 2)goto 15271                                          
      al=al*alf                                                         
      goto 15251                                                        
15271 if(ilm .ne. 1)goto 15281                                          
      al=big                                                            
      goto 15291                                                        
15281 continue                                                          
      al0=0.0                                                           
      do 15301 j=1,ni                                                   
      if(ju(j).eq.0)goto 15301                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
15301 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
15291 continue                                                          
15251 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 15311 k=1,ni                                                   
      if(ixx(k).eq.1)goto 15311                                         
      if(ju(k).eq.0)goto 15311                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
15311 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
15321 continue                                                          
      bs(0)=b(0)                                                        
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                             
      do 15331 j=1,ni                                                   
      if(ixx(j).eq.0)goto 15331                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      jn=ix(j+1)-ix(j)                                                  
      sc(1:jn)=v(jx(jb:je))                                             
      xm(j)=dot_product(sc(1:jn),x(jb:je))                              
      if(kopt .ne. 0)goto 15351                                         
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                           
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2             
15351 continue                                                          
15331 continue                                                          
      continue                                                          
      continue                                                          
15361 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 15371 k=1,ni                                                   
      if(ixx(k).eq.0)goto 15371                                         
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      jn=ix(k+1)-ix(k)                                                  
      bk=b(k)                                                           
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                              
      gk=dot_product(sc(1:jn),x(jb:je))                                 
      gk=(gk-svr*xb(k))/xs(k)                                           
      u=gk+xv(k)*b(k)                                                   
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 15391                                         
      b(k)=0.0                                                          
      goto 15401                                                        
15391 continue                                                          
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))       
15401 continue                                                          
      continue                                                          
      d=b(k)-bk                                                         
      if(abs(d).le.0.0)goto 15371                                       
      dlx=max(dlx,xv(k)*d**2)                                           
      if(mm(k) .ne. 0)goto 15421                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 15372                                           
      mm(k)=nin                                                         
      m(nin)=k                                                          
      sc(1:jn)=v(jx(jb:je))                                             
      xm(k)=dot_product(sc(1:jn),x(jb:je))                              
15421 continue                                                          
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)           
      o=o+d*(xb(k)/xs(k))                                               
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                               
15371 continue                                                          
15372 continue                                                          
      if(nin.gt.nx)goto 15362                                           
      d=0.0                                                             
      if(intr.ne.0) d=svr/xm(0)                                         
      if(d .eq. 0.0)goto 15441                                          
      b(0)=b(0)+d                                                       
      dlx=max(dlx,xm(0)*d**2)                                           
      r=r-d*v                                                           
      svr=svr-d*xm(0)                                                   
15441 continue                                                          
      if(dlx.lt.shr)goto 15362                                          
      if(nlp .le. maxit)goto 15461                                      
      jerr=-ilm                                                         
      return                                                            
15461 continue                                                          
      continue                                                          
15471 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 15481 l=1,nin                                                  
      k=m(l)                                                            
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      jn=ix(k+1)-ix(k)                                                  
      bk=b(k)                                                           
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                              
      gk=dot_product(sc(1:jn),x(jb:je))                                 
      gk=(gk-svr*xb(k))/xs(k)                                           
      u=gk+xv(k)*b(k)                                                   
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 15501                                         
      b(k)=0.0                                                          
      goto 15511                                                        
15501 continue                                                          
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))       
15511 continue                                                          
      continue                                                          
      d=b(k)-bk                                                         
      if(abs(d).le.0.0)goto 15481                                       
      dlx=max(dlx,xv(k)*d**2)                                           
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)           
      o=o+d*(xb(k)/xs(k))                                               
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                               
15481 continue                                                          
      continue                                                          
      d=0.0                                                             
      if(intr.ne.0) d=svr/xm(0)                                         
      if(d .eq. 0.0)goto 15531                                          
      b(0)=b(0)+d                                                       
      dlx=max(dlx,xm(0)*d**2)                                           
      r=r-d*v                                                           
      svr=svr-d*xm(0)                                                   
15531 continue                                                          
      if(dlx.lt.shr)goto 15472                                          
      if(nlp .le. maxit)goto 15551                                      
      jerr=-ilm                                                         
      return                                                            
15551 continue                                                          
      goto 15471                                                        
15472 continue                                                          
      goto 15361                                                        
15362 continue                                                          
      if(nin.gt.nx)goto 15322                                           
      sc=b(0)                                                           
      b0=0.0                                                            
      do 15561 j=1,nin                                                  
      l=m(j)                                                            
      jb=ix(l)                                                          
      je=ix(l+1)-1                                                      
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                   
      b0=b0-b(l)*xb(l)/xs(l)                                            
15561 continue                                                          
      continue                                                          
      sc=sc+b0                                                          
      do 15571 i=1,no                                                   
      fi=sc(i)+g(i)                                                     
      if(fi .ge. fmin)goto 15591                                        
      q(i)=0.0                                                          
      goto 15581                                                        
15591 if(fi .le. fmax)goto 15601                                        
      q(i)=1.0                                                          
      goto 15611                                                        
15601 continue                                                          
      q(i)=1.0/(1.0+exp(-fi))                                           
15611 continue                                                          
15581 continue                                                          
15571 continue                                                          
      continue                                                          
      v=w*q*(1.0-q)                                                     
      xm(0)=sum(v)                                                      
      if(xm(0).lt.vmin)goto 15322                                       
      r=w*(y-q)                                                         
      svr=sum(r)                                                        
      o=0.0                                                             
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 15631                      
      kx=0                                                              
      do 15641 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 15641                        
      kx=1                                                              
      goto 15642                                                        
15641 continue                                                          
15642 continue                                                          
      if(kx .ne. 0)goto 15661                                           
      do 15671 j=1,ni                                                   
      if(ixx(j).eq.1)goto 15671                                         
      if(ju(j).eq.0)goto 15671                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      jn=ix(j+1)-ix(j)                                                  
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                              
      gj=dot_product(sc(1:jn),x(jb:je))                                 
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                   
      if(ga(j) .le. al1*vp(j))goto 15691                                
      ixx(j)=1                                                          
      kx=1                                                              
15691 continue                                                          
15671 continue                                                          
      continue                                                          
      if(kx.eq.1) go to 10880                                           
      goto 15322                                                        
15661 continue                                                          
15631 continue                                                          
      goto 15321                                                        
15322 continue                                                          
      if(nin .le. nx)goto 15711                                         
      jerr=-10000-ilm                                                   
      goto 15242                                                        
15711 continue                                                          
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                             
      kin(ilm)=nin                                                      
      a0(ilm)=b(0)                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      devi=dev2(no,w,y,q,pmin)                                          
      dev(ilm)=(dev1-devi)/dev0                                         
      if(ilm.lt.mnl)goto 15241                                          
      if(flmin.ge.1.0)goto 15241                                        
      me=0                                                              
      do 15721 j=1,nin                                                  
      if(a(j,ilm).ne.0.0) me=me+1                                       
15721 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 15242                                            
      if(dev(ilm).gt.devmax)goto 15242                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 15242                          
      if(xm(0).lt.vmin)goto 15242                                       
15241 continue                                                          
15242 continue                                                          
      g=log(q/(1.0-q))                                                  
      deallocate(xm,b,bs,v,r,sc,xv,q,mm,ga,ixx)                         
      return                                                            
      end                                                               
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,n
     *lam,flmin,  ulam,shri,isd,intr,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev
     *0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb
     *(ni),xs(ni)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(
     *2,ni)
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                        
      double precision, dimension (:,:), allocatable :: q               
      double precision, dimension (:), allocatable :: sxp,sxpl          
      double precision, dimension (:), allocatable :: sc,xm,v,r,ga      
      double precision, dimension (:,:), allocatable :: b,bs,xv         
      integer, dimension (:), allocatable :: mm,is,iy                   
      allocate(b(0:ni,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni,1:nc),stat=jerr)                                 
      if(jerr.ne.0) return                                              
      allocate(bs(0:ni,1:nc),stat=jerr)                                 
      if(jerr.ne.0) return                                              
      allocate(q(1:no,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx,itrace)     
      exmn=-exmx                                                        
      allocate(xm(0:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(r(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(v(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ga(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(iy(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(is(1:max(nc,ni)),stat=jerr)                              
      if(jerr.ne.0) return                                              
      allocate(sxp(1:no),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(sxpl(1:no),stat=jerr)                                    
      if(jerr.ne.0) return                                              
      allocate(sc(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      pmax=1.0-pmin                                                     
      emin=pmin/pmax                                                    
      emax=1.0/emin                                                     
      pfm=(1.0+pmin)*pmin                                               
      pfx=(1.0-pmin)*pmax                                               
      vmin=pfm*pmax                                                     
      bta=parm                                                          
      omb=1.0-bta                                                       
      dev1=0.0                                                          
      dev0=0.0                                                          
      do 15731 ic=1,nc                                                  
      q0=dot_product(w,y(:,ic))                                         
      if(q0 .gt. pmin)goto 15751                                        
      jerr =8000+ic                                                     
      return                                                            
15751 continue                                                          
      if(q0 .lt. 1.0-pmin)goto 15771                                    
      jerr =9000+ic                                                     
      return                                                            
15771 continue                                                          
      if(intr.eq.0) q0=1.0/nc                                           
      b(1:ni,ic)=0.0                                                    
      b(0,ic)=0.0                                                       
      if(intr .eq. 0)goto 15791                                         
      b(0,ic)=log(q0)                                                   
      dev1=dev1-q0*b(0,ic)                                              
15791 continue                                                          
15731 continue                                                          
      continue                                                          
      if(intr.eq.0) dev1=log(float(nc))                                 
      iy=0                                                              
      al=0.0                                                            
      if(nonzero(no*nc,g) .ne. 0)goto 15811                             
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                      
      sxp=0.0                                                           
      do 15821 ic=1,nc                                                  
      q(:,ic)=exp(b(0,ic))                                              
      sxp=sxp+q(:,ic)                                                   
15821 continue                                                          
      continue                                                          
      goto 15831                                                        
15811 continue                                                          
      do 15841 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
15841 continue                                                          
      continue                                                          
      sxp=0.0                                                           
      if(intr .ne. 0)goto 15861                                         
      b(0,:)=0.0                                                        
      goto 15871                                                        
15861 continue                                                          
      call kazero(nc,no,y,g,w,b(0,:),jerr)                              
      if(jerr.ne.0) return                                              
15871 continue                                                          
      continue                                                          
      dev1=0.0                                                          
      do 15881 ic=1,nc                                                  
      q(:,ic)=b(0,ic)+g(:,ic)                                           
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                          
      q(:,ic)=exp(q(:,ic))                                              
      sxp=sxp+q(:,ic)                                                   
15881 continue                                                          
      continue                                                          
      sxpl=w*log(sxp)                                                   
      do 15891 ic=1,nc                                                  
      dev1=dev1+dot_product(y(:,ic),sxpl)                               
15891 continue                                                          
      continue                                                          
15831 continue                                                          
      continue                                                          
      do 15901 ic=1,nc                                                  
      do 15911 i=1,no                                                   
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))            
15911 continue                                                          
      continue                                                          
15901 continue                                                          
      continue                                                          
      dev0=dev0+dev1                                                    
      if(kopt .le. 0)goto 15931                                         
      if(isd .le. 0 .or. intr .eq. 0)goto 15951                         
      xv=0.25                                                           
      goto 15961                                                        
15951 continue                                                          
      do 15971 j=1,ni                                                   
      if(ju(j).eq.0)goto 15971                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)     
15971 continue                                                          
      continue                                                          
15961 continue                                                          
      continue                                                          
15931 continue                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 15991                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
15991 continue                                                          
      m=0                                                               
      mm=0                                                              
      nin=0                                                             
      nlp=0                                                             
      mnl=min(mnlam,nlam)                                               
      bs=0.0                                                            
      svr=0.0                                                           
      o=0.0                                                             
      shr=shri*dev0                                                     
      ga=0.0                                                            
      do 16001 ic=1,nc                                                  
      v=q(:,ic)/sxp                                                     
      r=w*(y(:,ic)-v)                                                   
      v=w*v*(1.0-v)                                                     
      do 16011 j=1,ni                                                   
      if(ju(j).eq.0)goto 16011                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      jn=ix(j+1)-ix(j)                                                  
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                              
      gj=dot_product(sc(1:jn),x(jb:je))                                 
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                          
16011 continue                                                          
      continue                                                          
16001 continue                                                          
      continue                                                          
      do 16021 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 16041                                      
      al=ulam(ilm)                                                      
      goto 16031                                                        
16041 if(ilm .le. 2)goto 16051                                          
      al=al*alf                                                         
      goto 16031                                                        
16051 if(ilm .ne. 1)goto 16061                                          
      al=big                                                            
      goto 16071                                                        
16061 continue                                                          
      al0=0.0                                                           
      do 16081 j=1,ni                                                   
      if(ju(j).eq.0)goto 16081                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
16081 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
16071 continue                                                          
16031 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 16091 k=1,ni                                                   
      if(iy(k).eq.1)goto 16091                                          
      if(ju(k).eq.0)goto 16091                                          
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                   
16091 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
16101 continue                                                          
      ixx=0                                                             
      jxx=ixx                                                           
      ig=0                                                              
      do 16111 ic=1,nc                                                  
      bs(0,ic)=b(0,ic)                                                  
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                       
      xm(0)=0.0                                                         
      svr=0.0                                                           
      o=0.0                                                             
      do 16121 i=1,no                                                   
      pic=q(i,ic)/sxp(i)                                                
      if(pic .ge. pfm)goto 16141                                        
      pic=0.0                                                           
      v(i)=0.0                                                          
      goto 16131                                                        
16141 if(pic .le. pfx)goto 16151                                        
      pic=1.0                                                           
      v(i)=0.0                                                          
      goto 16161                                                        
16151 continue                                                          
      v(i)=w(i)*pic*(1.0-pic)                                           
      xm(0)=xm(0)+v(i)                                                  
16161 continue                                                          
16131 continue                                                          
      r(i)=w(i)*(y(i,ic)-pic)                                           
      svr=svr+r(i)                                                      
16121 continue                                                          
      continue                                                          
      if(xm(0).le.vmin)goto 16111                                       
      ig=1                                                              
      do 16171 j=1,ni                                                   
      if(iy(j).eq.0)goto 16171                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                          
      if(kopt .ne. 0)goto 16191                                         
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                    
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2       
16191 continue                                                          
16171 continue                                                          
      continue                                                          
      continue                                                          
16201 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 16211 k=1,ni                                                   
      if(iy(k).eq.0)goto 16211                                          
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      jn=ix(k+1)-ix(k)                                                  
      bk=b(k,ic)                                                        
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                              
      gk=dot_product(sc(1:jn),x(jb:je))                                 
      gk=(gk-svr*xb(k))/xs(k)                                           
      u=gk+xv(k,ic)*b(k,ic)                                             
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 16231                                         
      b(k,ic)=0.0                                                       
      goto 16241                                                        
16231 continue                                                          
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))
     *)
16241 continue                                                          
      continue                                                          
      d=b(k,ic)-bk                                                      
      if(abs(d).le.0.0)goto 16211                                       
      dlx=max(dlx,xv(k,ic)*d**2)                                        
      if(mm(k) .ne. 0)goto 16261                                        
      nin=nin+1                                                         
      if(nin .le. nx)goto 16281                                         
      jxx=1                                                             
      goto 16212                                                        
16281 continue                                                          
      mm(k)=nin                                                         
      m(nin)=k                                                          
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                          
16261 continue                                                          
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)           
      o=o+d*(xb(k)/xs(k))                                               
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                               
16211 continue                                                          
16212 continue                                                          
      if(jxx.gt.0)goto 16202                                            
      d=0.0                                                             
      if(intr.ne.0) d=svr/xm(0)                                         
      if(d .eq. 0.0)goto 16301                                          
      b(0,ic)=b(0,ic)+d                                                 
      dlx=max(dlx,xm(0)*d**2)                                           
      r=r-d*v                                                           
      svr=svr-d*xm(0)                                                   
16301 continue                                                          
      if(dlx.lt.shr)goto 16202                                          
      if(nlp .le. maxit)goto 16321                                      
      jerr=-ilm                                                         
      return                                                            
16321 continue                                                          
      continue                                                          
16331 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 16341 l=1,nin                                                  
      k=m(l)                                                            
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      jn=ix(k+1)-ix(k)                                                  
      bk=b(k,ic)                                                        
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                              
      gk=dot_product(sc(1:jn),x(jb:je))                                 
      gk=(gk-svr*xb(k))/xs(k)                                           
      u=gk+xv(k,ic)*b(k,ic)                                             
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 16361                                         
      b(k,ic)=0.0                                                       
      goto 16371                                                        
16361 continue                                                          
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))
     *)
16371 continue                                                          
      continue                                                          
      d=b(k,ic)-bk                                                      
      if(abs(d).le.0.0)goto 16341                                       
      dlx=max(dlx,xv(k,ic)*d**2)                                        
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)           
      o=o+d*(xb(k)/xs(k))                                               
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                               
16341 continue                                                          
      continue                                                          
      d=0.0                                                             
      if(intr.ne.0) d=svr/xm(0)                                         
      if(d .eq. 0.0)goto 16391                                          
      b(0,ic)=b(0,ic)+d                                                 
      dlx=max(dlx,xm(0)*d**2)                                           
      r=r-d*v                                                           
      svr=svr-d*xm(0)                                                   
16391 continue                                                          
      if(dlx.lt.shr)goto 16332                                          
      if(nlp .le. maxit)goto 16411                                      
      jerr=-ilm                                                         
      return                                                            
16411 continue                                                          
      goto 16331                                                        
16332 continue                                                          
      goto 16201                                                        
16202 continue                                                          
      if(jxx.gt.0)goto 16112                                            
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                      
      if(ixx .ne. 0)goto 16431                                          
      do 16441 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 16461             
      ixx=1                                                             
      goto 16442                                                        
16461 continue                                                          
16441 continue                                                          
16442 continue                                                          
16431 continue                                                          
      sc=b(0,ic)+g(:,ic)                                                
      b0=0.0                                                            
      do 16471 j=1,nin                                                  
      l=m(j)                                                            
      jb=ix(l)                                                          
      je=ix(l+1)-1                                                      
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                
      b0=b0-b(l,ic)*xb(l)/xs(l)                                         
16471 continue                                                          
      continue                                                          
      sc=min(max(exmn,sc+b0),exmx)                                      
      sxp=sxp-q(:,ic)                                                   
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                       
      sxp=sxp+q(:,ic)                                                   
16111 continue                                                          
16112 continue                                                          
      s=-sum(b(0,:))/nc                                                 
      b(0,:)=b(0,:)+s                                                   
      sc=s                                                              
      b0=0.0                                                            
      do 16481 j=1,nin                                                  
      l=m(j)                                                            
      if(vp(l) .gt. 0.0)goto 16501                                      
      s=sum(b(l,:))/nc                                                  
      goto 16511                                                        
16501 continue                                                          
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                  
16511 continue                                                          
      continue                                                          
      b(l,:)=b(l,:)-s                                                   
      jb=ix(l)                                                          
      je=ix(l+1)-1                                                      
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                      
      b0=b0+s*xb(l)/xs(l)                                               
16481 continue                                                          
      continue                                                          
      sc=sc+b0                                                          
      sc=exp(sc)                                                        
      sxp=sxp*sc                                                        
      do 16521 ic=1,nc                                                  
      q(:,ic)=q(:,ic)*sc                                                
16521 continue                                                          
      continue                                                          
      if(jxx.gt.0)goto 16102                                            
      if(ig.eq.0)goto 16102                                             
      if(ixx .ne. 0)goto 16541                                          
      do 16551 j=1,ni                                                   
      if(iy(j).eq.1)goto 16551                                          
      if(ju(j).eq.0)goto 16551                                          
      ga(j)=0.0                                                         
16551 continue                                                          
      continue                                                          
      do 16561 ic=1,nc                                                  
      v=q(:,ic)/sxp                                                     
      r=w*(y(:,ic)-v)                                                   
      v=w*v*(1.0-v)                                                     
      do 16571 j=1,ni                                                   
      if(iy(j).eq.1)goto 16571                                          
      if(ju(j).eq.0)goto 16571                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      jn=ix(j+1)-ix(j)                                                  
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                              
      gj=dot_product(sc(1:jn),x(jb:je))                                 
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                          
16571 continue                                                          
      continue                                                          
16561 continue                                                          
      continue                                                          
      do 16581 k=1,ni                                                   
      if(iy(k).eq.1)goto 16581                                          
      if(ju(k).eq.0)goto 16581                                          
      if(ga(k) .le. al1*vp(k))goto 16601                                
      iy(k)=1                                                           
      ixx=1                                                             
16601 continue                                                          
16581 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10880                                          
      goto 16102                                                        
16541 continue                                                          
      goto 16101                                                        
16102 continue                                                          
      if(jxx .le. 0)goto 16621                                          
      jerr=-10000-ilm                                                   
      goto 16022                                                        
16621 continue                                                          
      devi=0.0                                                          
      do 16631 ic=1,nc                                                  
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                       
      a0(ic,ilm)=b(0,ic)                                                
      do 16641 i=1,no                                                   
      if(y(i,ic).le.0.0)goto 16641                                      
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                        
16641 continue                                                          
      continue                                                          
16631 continue                                                          
      continue                                                          
      kin(ilm)=nin                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dev1-devi)/dev0                                         
      if(ig.eq.0)goto 16022                                             
      if(ilm.lt.mnl)goto 16021                                          
      if(flmin.ge.1.0)goto 16021                                        
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 16022          
      if(dev(ilm).gt.devmax)goto 16022                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 16022                          
16021 continue                                                          
16022 continue                                                          
      g=log(q)                                                          
      do 16651 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
16651 continue                                                          
      continue                                                          
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy)                   
      return                                                            
      end                                                               
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)               
      implicit double precision(a-h,o-z)                                
      double precision a0(nc),ca(nx,nc),x(*),f(nc,n)                    
      integer ia(*),ix(*),jx(*)                                         
      do 16661 ic=1,nc                                                  
      f(ic,:)=a0(ic)                                                    
16661 continue                                                          
      continue                                                          
      do 16671 j=1,nin                                                  
      k=ia(j)                                                           
      kb=ix(k)                                                          
      ke=ix(k+1)-1                                                      
      do 16681 ic=1,nc                                                  
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                 
16681 continue                                                          
      continue                                                          
16671 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine coxnet(parm,no,ni,x,y,d,g,w,jd,vp,cl,ne,nx,nlam,flmin,u
     *lam,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam
     *)
      double precision ca(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)         
      integer jd(*),ia(nx),nin(nlam)                                    
      double precision, dimension (:), allocatable :: xs,ww,vq          
      integer, dimension (:), allocatable :: ju                         
      if(maxval(vp) .gt. 0.0)goto 16701                                 
      jerr=10000                                                        
      return                                                            
16701 continue                                                          
      allocate(ww(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      if(isd .le. 0)goto 16721                                          
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
16721 continue                                                          
      call chkvars(no,ni,x,ju)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 16741                                   
      jerr=7777                                                         
      return                                                            
16741 continue                                                          
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      ww=max(0d0,w)                                                     
      sw=sum(ww)                                                        
      if(sw .gt. 0.0)goto 16761                                         
      jerr=9999                                                         
      return                                                            
16761 continue                                                          
      ww=ww/sw                                                          
      call cstandard(no,ni,x,ww,ju,isd,xs)                              
      if(isd .le. 0)goto 16781                                          
      do 16791 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
16791 continue                                                          
      continue                                                          
16781 continue                                                          
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,
     *thr,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                              
      dev0=2.0*sw*dev0                                                  
      if(isd .le. 0)goto 16811                                          
      do 16821 k=1,lmu                                                  
      nk=nin(k)                                                         
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                
16821 continue                                                          
      continue                                                          
16811 continue                                                          
      deallocate(ww,ju,vq)                                              
      if(isd.gt.0) deallocate(xs)                                       
      return                                                            
      end                                                               
      subroutine cstandard(no,ni,x,w,ju,isd,xs)                         
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),w(no),xs(ni)                            
      integer ju(ni)                                                    
      do 16831 j=1,ni                                                   
      if(ju(j).eq.0)goto 16831                                          
      xm=dot_product(w,x(:,j))                                          
      x(:,j)=x(:,j)-xm                                                  
      if(isd .le. 0)goto 16851                                          
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                              
      x(:,j)=x(:,j)/xs(j)                                               
16851 continue                                                          
16831 continue                                                          
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
      if(jerr.ne.0) go to 12180                                         
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
      if(nonzero(no,g) .eq. 0)goto 16871                                
      f=g-dot_product(q,g)                                              
      e=q*exp(sign(min(abs(f),fmax),f))                                 
      goto 16881                                                        
16871 continue                                                          
      f=0.0                                                             
      e=q                                                               
16881 continue                                                          
      continue                                                          
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                              
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                      
      dev0=rr                                                           
      do 16891 i=1,no                                                   
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 16911                
      w(i)=0.0                                                          
      wr(i)=w(i)                                                        
16911 continue                                                          
16891 continue                                                          
      continue                                                          
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                      
      if(jerr.ne.0) go to 12180                                         
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 16931                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
16931 continue                                                          
      m=0                                                               
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      mnl=min(mnlam,nlam)                                               
      as=0.0                                                            
      cthr=cthri*dev0                                                   
      do 16941 j=1,ni                                                   
      if(ju(j).eq.0)goto 16941                                          
      ga(j)=abs(dot_product(wr,x(:,j)))                                 
16941 continue                                                          
      continue                                                          
      do 16951 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 16971                                      
      al=ulam(ilm)                                                      
      goto 16961                                                        
16971 if(ilm .le. 2)goto 16981                                          
      al=al*alf                                                         
      goto 16961                                                        
16981 if(ilm .ne. 1)goto 16991                                          
      al=big                                                            
      goto 17001                                                        
16991 continue                                                          
      al0=0.0                                                           
      do 17011 j=1,ni                                                   
      if(ju(j).eq.0)goto 17011                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
17011 continue                                                          
      continue                                                          
      al0=al0/max(parm,1.0d-3)                                          
      al=alf*al0                                                        
17001 continue                                                          
16961 continue                                                          
      sa=alpha*al                                                       
      omal=oma*al                                                       
      tlam=alpha*(2.0*al-al0)                                           
      do 17021 k=1,ni                                                   
      if(ixx(k).eq.1)goto 17021                                         
      if(ju(k).eq.0)goto 17021                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
17021 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
17031 continue                                                          
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                             
      call vars(no,ni,x,w,ixx,v)                                        
      continue                                                          
17041 continue                                                          
      nlp=nlp+1                                                         
      dli=0.0                                                           
      do 17051 j=1,ni                                                   
      if(ixx(j).eq.0)goto 17051                                         
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                
      if(abs(u) .gt. vp(j)*sa)goto 17071                                
      at=0.0                                                            
      goto 17081                                                        
17071 continue                                                          
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o
     *mal)))
17081 continue                                                          
      continue                                                          
      if(at .eq. a(j))goto 17101                                        
      del=at-a(j)                                                       
      a(j)=at                                                           
      dli=max(dli,v(j)*del**2)                                          
      wr=wr-del*w*x(:,j)                                                
      f=f+del*x(:,j)                                                    
      if(mm(j) .ne. 0)goto 17121                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 17052                                           
      mm(j)=nin                                                         
      m(nin)=j                                                          
17121 continue                                                          
17101 continue                                                          
17051 continue                                                          
17052 continue                                                          
      if(nin.gt.nx)goto 17042                                           
      if(dli.lt.cthr)goto 17042                                         
      if(nlp .le. maxit)goto 17141                                      
      jerr=-ilm                                                         
      return                                                            
17141 continue                                                          
      continue                                                          
17151 continue                                                          
      nlp=nlp+1                                                         
      dli=0.0                                                           
      do 17161 l=1,nin                                                  
      j=m(l)                                                            
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                
      if(abs(u) .gt. vp(j)*sa)goto 17181                                
      at=0.0                                                            
      goto 17191                                                        
17181 continue                                                          
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o
     *mal)))
17191 continue                                                          
      continue                                                          
      if(at .eq. a(j))goto 17211                                        
      del=at-a(j)                                                       
      a(j)=at                                                           
      dli=max(dli,v(j)*del**2)                                          
      wr=wr-del*w*x(:,j)                                                
      f=f+del*x(:,j)                                                    
17211 continue                                                          
17161 continue                                                          
      continue                                                          
      if(dli.lt.cthr)goto 17152                                         
      if(nlp .le. maxit)goto 17231                                      
      jerr=-ilm                                                         
      return                                                            
17231 continue                                                          
      goto 17151                                                        
17152 continue                                                          
      goto 17041                                                        
17042 continue                                                          
      if(nin.gt.nx)goto 17032                                           
      e=q*exp(sign(min(abs(f),fmax),f))                                 
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                      
      if(jerr .eq. 0)goto 17251                                         
      jerr=jerr-ilm                                                     
      go to 12180                                                       
17251 continue                                                          
      ix=0                                                              
      do 17261 j=1,nin                                                  
      k=m(j)                                                            
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 17261                        
      ix=1                                                              
      goto 17262                                                        
17261 continue                                                          
17262 continue                                                          
      if(ix .ne. 0)goto 17281                                           
      do 17291 k=1,ni                                                   
      if(ixx(k).eq.1)goto 17291                                         
      if(ju(k).eq.0)goto 17291                                          
      ga(k)=abs(dot_product(wr,x(:,k)))                                 
      if(ga(k) .le. sa*vp(k))goto 17311                                 
      ixx(k)=1                                                          
      ix=1                                                              
17311 continue                                                          
17291 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10880                                           
      goto 17032                                                        
17281 continue                                                          
      goto 17031                                                        
17032 continue                                                          
      if(nin .le. nx)goto 17331                                         
      jerr=-10000-ilm                                                   
      goto 16952                                                        
17331 continue                                                          
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                            
      kin(ilm)=nin                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                
      if(ilm.lt.mnl)goto 16951                                          
      if(flmin.ge.1.0)goto 16951                                        
      me=0                                                              
      do 17341 j=1,nin                                                  
      if(ao(j,ilm).ne.0.0) me=me+1                                      
17341 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 16952                                            
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 16952           
      if(dev(ilm).gt.devmax)goto 16952                                  
16951 continue                                                          
16952 continue                                                          
      g=f                                                               
12180 continue                                                          
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)           
      return                                                            
      end                                                               
      subroutine cxmodval(ca,ia,nin,n,x,f)                              
      implicit double precision(a-h,o-z)                                
      double precision ca(nin),x(n,*),f(n)                              
      integer ia(nin)                                                   
      f=0.0                                                             
      if(nin.le.0) return                                               
      do 17351 i=1,n                                                    
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                   
17351 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                      
      implicit double precision(a-h,o-z)                                
      double precision y(no),d(no),q(no)                                
      integer jp(no),kp(*)                                              
      do 17361 j=1,no                                                   
      jp(j)=j                                                           
17361 continue                                                          
      continue                                                          
      call psort7(y,jp,1,no)                                            
      nj=0                                                              
      do 17371 j=1,no                                                   
      if(q(jp(j)).le.0.0)goto 17371                                     
      nj=nj+1                                                           
      jp(nj)=jp(j)                                                      
17371 continue                                                          
      continue                                                          
      if(nj .ne. 0)goto 17391                                           
      jerr=20000                                                        
      return                                                            
17391 continue                                                          
      j=1                                                               
      continue                                                          
17401 if(d(jp(j)).gt.0.0)goto 17402                                     
      j=j+1                                                             
      if(j.gt.nj)goto 17402                                             
      goto 17401                                                        
17402 continue                                                          
      if(j .lt. nj-1)goto 17421                                         
      jerr=30000                                                        
      return                                                            
17421 continue                                                          
      t0=y(jp(j))                                                       
      j0=j-1                                                            
      if(j0 .le. 0)goto 17441                                           
      continue                                                          
17451 if(y(jp(j0)).lt.t0)goto 17452                                     
      j0=j0-1                                                           
      if(j0.eq.0)goto 17452                                             
      goto 17451                                                        
17452 continue                                                          
      if(j0 .le. 0)goto 17471                                           
      nj=nj-j0                                                          
      do 17481 j=1,nj                                                   
      jp(j)=jp(j+j0)                                                    
17481 continue                                                          
      continue                                                          
17471 continue                                                          
17441 continue                                                          
      jerr=0                                                            
      nk=0                                                              
      yk=t0                                                             
      j=2                                                               
      continue                                                          
17491 continue                                                          
      continue                                                          
17501 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 17502                  
      j=j+1                                                             
      if(j.gt.nj)goto 17502                                             
      goto 17501                                                        
17502 continue                                                          
      nk=nk+1                                                           
      kp(nk)=j-1                                                        
      if(j.gt.nj)goto 17492                                             
      if(j .ne. nj)goto 17521                                           
      nk=nk+1                                                           
      kp(nk)=nj                                                         
      goto 17492                                                        
17521 continue                                                          
      yk=y(jp(j))                                                       
      j=j+1                                                             
      goto 17491                                                        
17492 continue                                                          
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
      do 17531 j=1,kp(1)                                                
      i=jp(j)                                                           
      w(i)=e(i)*(b-e(i)*c)                                              
      if(w(i) .gt. 0.0)goto 17551                                       
      jerr=-30000                                                       
      return                                                            
17551 continue                                                          
      wr(i)=d(i)-e(i)*b                                                 
17531 continue                                                          
      continue                                                          
      do 17561 k=2,nk                                                   
      j1=kp(k-1)+1                                                      
      j2=kp(k)                                                          
      b=b+dk(k)/u(k)                                                    
      c=c+dk(k)/u(k)**2                                                 
      do 17571 j=j1,j2                                                  
      i=jp(j)                                                           
      w(i)=e(i)*(b-e(i)*c)                                              
      if(w(i) .gt. 0.0)goto 17591                                       
      jerr=-30000                                                       
      return                                                            
17591 continue                                                          
      wr(i)=d(i)-e(i)*b                                                 
17571 continue                                                          
      continue                                                          
17561 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine vars(no,ni,x,w,ixx,v)                                  
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),w(no),v(ni)                             
      integer ixx(ni)                                                   
      do 17601 j=1,ni                                                   
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                     
17601 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine died(no,nk,d,kp,jp,dk)                                 
      implicit double precision(a-h,o-z)                                
      double precision d(no),dk(nk)                                     
      integer kp(nk),jp(no)                                             
      dk(1)=sum(d(jp(1:kp(1))))                                         
      do 17611 k=2,nk                                                   
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                               
17611 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine usk(no,nk,kp,jp,e,u)                                   
      implicit double precision(a-h,o-z)                                
      double precision e(no),u(nk),h                                    
      integer kp(nk),jp(no)                                             
      h=0.0                                                             
      do 17621 k=nk,1,-1                                                
      j2=kp(k)                                                          
      j1=1                                                              
      if(k.gt.1) j1=kp(k-1)+1                                           
      do 17631 j=j2,j1,-1                                               
      h=h+e(jp(j))                                                      
17631 continue                                                          
      continue                                                          
      u(k)=h                                                            
17621 continue                                                          
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
      if(sw .gt. 0.0)goto 17651                                         
      jerr=9999                                                         
      go to 12180                                                       
17651 continue                                                          
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                            
      if(jerr.ne.0) go to 12180                                         
      fmax=log(huge(e(1))*0.1)                                          
      dq=d*q                                                            
      call died(no,nk,dq,kp,jp,dk)                                      
      gm=dot_product(q,g)/sw                                            
      do 17661 j=1,ni                                                   
      xm(j)=dot_product(q,x(:,j))/sw                                    
17661 continue                                                          
      continue                                                          
      do 17671 lam=1,nlam                                               
      do 17681 i=1,no                                                   
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                    
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                     
17681 continue                                                          
      continue                                                          
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                       
17671 continue                                                          
      continue                                                          
12180 continue                                                          
      deallocate(e,uu,dk,f,jp,kp,dq)                                    
      return                                                            
      end                                                               
      subroutine fishnet(parm,no,ni,x,y,g,w,jd,vp,cl,ne,nx,nlam,flmin,ul
     *am,thr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)     
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)
      integer jd(*),ia(nx),nin(nlam)                                    
      double precision, dimension (:), allocatable :: xm,xs,ww,vq       
      integer, dimension (:), allocatable :: ju                         
      if(maxval(vp) .gt. 0.0)goto 17701                                 
      jerr=10000                                                        
      return                                                            
17701 continue                                                          
      if(minval(y) .ge. 0.0)goto 17721                                  
      jerr=8888                                                         
      return                                                            
17721 continue                                                          
      allocate(ww(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      if(isd .le. 0)goto 17741                                          
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
17741 continue                                                          
      call chkvars(no,ni,x,ju)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 17761                                   
      jerr=7777                                                         
      go to 12180                                                       
17761 continue                                                          
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      ww=max(0d0,w)                                                     
      sw=sum(ww)                                                        
      if(sw .gt. 0.0)goto 17781                                         
      jerr=9999                                                         
      go to 12180                                                       
17781 continue                                                          
      ww=ww/sw                                                          
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                     
      if(isd .le. 0)goto 17801                                          
      do 17811 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
17811 continue                                                          
      continue                                                          
17801 continue                                                          
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,t
     *hr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 12180                                         
      dev0=2.0*sw*dev0                                                  
      do 17821 k=1,lmu                                                  
      nk=nin(k)                                                         
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                   
      if(intr .ne. 0)goto 17841                                         
      a0(k)=0.0                                                         
      goto 17851                                                        
17841 continue                                                          
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                  
17851 continue                                                          
      continue                                                          
17821 continue                                                          
      continue                                                          
12180 continue                                                          
      deallocate(ww,ju,vq,xm)                                           
      if(isd.gt.0) deallocate(xs)                                       
      return                                                            
      end                                                               
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,cl,ne,nx,nlam,flmin,u
     *lam,shri,  isd,intr,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)     
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)
      integer ju(ni),m(nx),kin(nlam)                                    
      double precision, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga
      integer, dimension (:), allocatable :: mm,ixx                     
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx,itrace)     
      sml=sml*10.0                                                      
      isd = isd*1                                                       
      allocate(a(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(as(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(t(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ga(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ixx(1:ni),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(wr(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(v(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(w(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(f(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      bta=parm                                                          
      omb=1.0-bta                                                       
      t=q*y                                                             
      yb=sum(t)                                                         
      fmax=log(huge(bta)*0.1)                                           
      if(nonzero(no,g) .ne. 0)goto 17871                                
      if(intr .eq. 0)goto 17891                                         
      w=q*yb                                                            
      az=log(yb)                                                        
      f=az                                                              
      dv0=yb*(az-1.0)                                                   
      goto 17901                                                        
17891 continue                                                          
      w=q                                                               
      az=0.0                                                            
      f=az                                                              
      dv0=-1.0                                                          
17901 continue                                                          
      continue                                                          
      goto 17911                                                        
17871 continue                                                          
      w=q*exp(sign(min(abs(g),fmax),g))                                 
      v0=sum(w)                                                         
      if(intr .eq. 0)goto 17931                                         
      eaz=yb/v0                                                         
      w=eaz*w                                                           
      az=log(eaz)                                                       
      f=az+g                                                            
      dv0=dot_product(t,g)-yb*(1.0-az)                                  
      goto 17941                                                        
17931 continue                                                          
      az=0.0                                                            
      f=g                                                               
      dv0=dot_product(t,g)-v0                                           
17941 continue                                                          
      continue                                                          
17911 continue                                                          
      continue                                                          
      a=0.0                                                             
      as=0.0                                                            
      wr=t-w                                                            
      v0=1.0                                                            
      if(intr.ne.0) v0=yb                                               
      dvr=-yb                                                           
      do 17951 i=1,no                                                   
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                            
17951 continue                                                          
      continue                                                          
      dvr=dvr-dv0                                                       
      dev0=dvr                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 17971                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
17971 continue                                                          
      m=0                                                               
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      mnl=min(mnlam,nlam)                                               
      shr=shri*dev0                                                     
      ixx=0                                                             
      al=0.0                                                            
      do 17981 j=1,ni                                                   
      if(ju(j).eq.0)goto 17981                                          
      ga(j)=abs(dot_product(wr,x(:,j)))                                 
17981 continue                                                          
      continue                                                          
      do 17991 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 18011                                      
      al=ulam(ilm)                                                      
      goto 18001                                                        
18011 if(ilm .le. 2)goto 18021                                          
      al=al*alf                                                         
      goto 18001                                                        
18021 if(ilm .ne. 1)goto 18031                                          
      al=big                                                            
      goto 18041                                                        
18031 continue                                                          
      al0=0.0                                                           
      do 18051 j=1,ni                                                   
      if(ju(j).eq.0)goto 18051                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
18051 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
18041 continue                                                          
18001 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 18061 k=1,ni                                                   
      if(ixx(k).eq.1)goto 18061                                         
      if(ju(k).eq.0)goto 18061                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
18061 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
18071 continue                                                          
      az0=az                                                            
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                             
      do 18081 j=1,ni                                                   
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                     
18081 continue                                                          
      continue                                                          
      continue                                                          
18091 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 18101 k=1,ni                                                   
      if(ixx(k).eq.0)goto 18101                                         
      ak=a(k)                                                           
      u=dot_product(wr,x(:,k))+v(k)*ak                                  
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 18121                                         
      a(k)=0.0                                                          
      goto 18131                                                        
18121 continue                                                          
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))        
18131 continue                                                          
      continue                                                          
      if(a(k).eq.ak)goto 18101                                          
      d=a(k)-ak                                                         
      dlx=max(dlx,v(k)*d**2)                                            
      wr=wr-d*w*x(:,k)                                                  
      f=f+d*x(:,k)                                                      
      if(mm(k) .ne. 0)goto 18151                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 18102                                           
      mm(k)=nin                                                         
      m(nin)=k                                                          
18151 continue                                                          
18101 continue                                                          
18102 continue                                                          
      if(nin.gt.nx)goto 18092                                           
      if(intr .eq. 0)goto 18171                                         
      d=sum(wr)/v0                                                      
      az=az+d                                                           
      dlx=max(dlx,v0*d**2)                                              
      wr=wr-d*w                                                         
      f=f+d                                                             
18171 continue                                                          
      if(dlx.lt.shr)goto 18092                                          
      if(nlp .le. maxit)goto 18191                                      
      jerr=-ilm                                                         
      return                                                            
18191 continue                                                          
      continue                                                          
18201 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 18211 l=1,nin                                                  
      k=m(l)                                                            
      ak=a(k)                                                           
      u=dot_product(wr,x(:,k))+v(k)*ak                                  
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 18231                                         
      a(k)=0.0                                                          
      goto 18241                                                        
18231 continue                                                          
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))        
18241 continue                                                          
      continue                                                          
      if(a(k).eq.ak)goto 18211                                          
      d=a(k)-ak                                                         
      dlx=max(dlx,v(k)*d**2)                                            
      wr=wr-d*w*x(:,k)                                                  
      f=f+d*x(:,k)                                                      
18211 continue                                                          
      continue                                                          
      if(intr .eq. 0)goto 18261                                         
      d=sum(wr)/v0                                                      
      az=az+d                                                           
      dlx=max(dlx,v0*d**2)                                              
      wr=wr-d*w                                                         
      f=f+d                                                             
18261 continue                                                          
      if(dlx.lt.shr)goto 18202                                          
      if(nlp .le. maxit)goto 18281                                      
      jerr=-ilm                                                         
      return                                                            
18281 continue                                                          
      goto 18201                                                        
18202 continue                                                          
      goto 18091                                                        
18092 continue                                                          
      if(nin.gt.nx)goto 18072                                           
      w=q*exp(sign(min(abs(f),fmax),f))                                 
      v0=sum(w)                                                         
      wr=t-w                                                            
      if(v0*(az-az0)**2 .ge. shr)goto 18301                             
      ix=0                                                              
      do 18311 j=1,nin                                                  
      k=m(j)                                                            
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 18311                         
      ix=1                                                              
      goto 18312                                                        
18311 continue                                                          
18312 continue                                                          
      if(ix .ne. 0)goto 18331                                           
      do 18341 k=1,ni                                                   
      if(ixx(k).eq.1)goto 18341                                         
      if(ju(k).eq.0)goto 18341                                          
      ga(k)=abs(dot_product(wr,x(:,k)))                                 
      if(ga(k) .le. al1*vp(k))goto 18361                                
      ixx(k)=1                                                          
      ix=1                                                              
18361 continue                                                          
18341 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10880                                           
      goto 18072                                                        
18331 continue                                                          
18301 continue                                                          
      goto 18071                                                        
18072 continue                                                          
      if(nin .le. nx)goto 18381                                         
      jerr=-10000-ilm                                                   
      goto 17992                                                        
18381 continue                                                          
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                            
      kin(ilm)=nin                                                      
      a0(ilm)=az                                                        
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                            
      if(ilm.lt.mnl)goto 17991                                          
      if(flmin.ge.1.0)goto 17991                                        
      me=0                                                              
      do 18391 j=1,nin                                                  
      if(ca(j,ilm).ne.0.0) me=me+1                                      
18391 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 17992                                            
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17992           
      if(dev(ilm).gt.devmax)goto 17992                                  
17991 continue                                                          
17992 continue                                                          
      g=f                                                               
      continue                                                          
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                             
      return                                                            
      end                                                               
      function nonzero(n,v)                                             
      implicit double precision(a-h,o-z)                                
      double precision v(n)                                             
      nonzero=0                                                         
      do 18401 i=1,n                                                    
      if(v(i) .eq. 0.0)goto 18421                                       
      nonzero=1                                                         
      return                                                            
18421 continue                                                          
18401 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                            
      implicit double precision(a-h,o-z)                                
      double precision a(nx,lmu),b(ni,lmu)                              
      integer ia(nx),nin(lmu)                                           
      do 18431 lam=1,lmu                                                
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                     
18431 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                        
      implicit double precision(a-h,o-z)                                
      double precision a(nx,nc,lmu),b(ni,nc,lmu)                        
      integer ia(nx),nin(lmu)                                           
      do 18441 lam=1,lmu                                                
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))          
18441 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)            
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),fl
     *og(nlam)
      double precision, dimension (:), allocatable :: w                 
      if(minval(y) .ge. 0.0)goto 18461                                  
      jerr=8888                                                         
      return                                                            
18461 continue                                                          
      allocate(w(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      w=max(0d0,q)                                                      
      sw=sum(w)                                                         
      if(sw .gt. 0.0)goto 18481                                         
      jerr=9999                                                         
      go to 12180                                                       
18481 continue                                                          
      yb=dot_product(w,y)/sw                                            
      fmax=log(huge(y(1))*0.1)                                          
      do 18491 lam=1,nlam                                               
      s=0.0                                                             
      do 18501 i=1,no                                                   
      if(w(i).le.0.0)goto 18501                                         
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                       
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                   
18501 continue                                                          
      continue                                                          
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                             
18491 continue                                                          
      continue                                                          
12180 continue                                                          
      deallocate(w)                                                     
      return                                                            
      end                                                               
      subroutine spfishnet(parm,no,ni,x,ix,jx,y,g,w,jd,vp,cl,ne,nx,nlam,
     *flmin,  ulam,thr,isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,
     *jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)         
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                        
      double precision, dimension (:), allocatable :: xm,xs,ww,vq       
      integer, dimension (:), allocatable :: ju                         
      if(maxval(vp) .gt. 0.0)goto 18521                                 
      jerr=10000                                                        
      return                                                            
18521 continue                                                          
      if(minval(y) .ge. 0.0)goto 18541                                  
      jerr=8888                                                         
      return                                                            
18541 continue                                                          
      allocate(ww(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      call spchkvars(no,ni,x,ix,ju)                                     
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 18561                                   
      jerr=7777                                                         
      go to 12180                                                       
18561 continue                                                          
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      ww=max(0d0,w)                                                     
      sw=sum(ww)                                                        
      if(sw .gt. 0.0)goto 18581                                         
      jerr=9999                                                         
      go to 12180                                                       
18581 continue                                                          
      ww=ww/sw                                                          
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)             
      if(isd .le. 0)goto 18601                                          
      do 18611 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
18611 continue                                                          
      continue                                                          
18601 continue                                                          
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,flmi
     *n,ulam,thr,  isd,intr,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nl
     *p,jerr)
      if(jerr.gt.0) go to 12180                                         
      dev0=2.0*sw*dev0                                                  
      do 18621 k=1,lmu                                                  
      nk=nin(k)                                                         
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                   
      if(intr .ne. 0)goto 18641                                         
      a0(k)=0.0                                                         
      goto 18651                                                        
18641 continue                                                          
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                  
18651 continue                                                          
      continue                                                          
18621 continue                                                          
      continue                                                          
12180 continue                                                          
      deallocate(ww,ju,vq,xm,xs)                                        
      return                                                            
      end                                                               
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,cl,ne,nx,nlam
     *,flmin,ulam,  shri,isd,intr,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,a
     *lm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),x
     *s(ni)
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                        
      double precision, dimension (:), allocatable :: qy,t,w,wr,v       
      double precision, dimension (:), allocatable :: a,as,xm,ga        
      integer, dimension (:), allocatable :: mm,ixx                     
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx,itrace)     
      sml=sml*10.0                                                      
      isd = isd*1                                                       
      allocate(a(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(as(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(t(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ga(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ixx(1:ni),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(wr(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(v(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(w(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(qy(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      bta=parm                                                          
      omb=1.0-bta                                                       
      fmax=log(huge(bta)*0.1)                                           
      qy=q*y                                                            
      yb=sum(qy)                                                        
      if(nonzero(no,g) .ne. 0)goto 18671                                
      t=0.0                                                             
      if(intr .eq. 0)goto 18691                                         
      w=q*yb                                                            
      az=log(yb)                                                        
      uu=az                                                             
      xm=yb*xb                                                          
      dv0=yb*(az-1.0)                                                   
      goto 18701                                                        
18691 continue                                                          
      w=q                                                               
      xm=0.0                                                            
      uu=0.0                                                            
      az=uu                                                             
      dv0=-1.0                                                          
18701 continue                                                          
      continue                                                          
      goto 18711                                                        
18671 continue                                                          
      w=q*exp(sign(min(abs(g),fmax),g))                                 
      ww=sum(w)                                                         
      t=g                                                               
      if(intr .eq. 0)goto 18731                                         
      eaz=yb/ww                                                         
      w=eaz*w                                                           
      az=log(eaz)                                                       
      uu=az                                                             
      dv0=dot_product(qy,g)-yb*(1.0-az)                                 
      goto 18741                                                        
18731 continue                                                          
      uu=0.0                                                            
      az=uu                                                             
      dv0=dot_product(qy,g)-ww                                          
18741 continue                                                          
      continue                                                          
      do 18751 j=1,ni                                                   
      if(ju(j).eq.0)goto 18751                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
18751 continue                                                          
      continue                                                          
18711 continue                                                          
      continue                                                          
      tt=yb*uu                                                          
      ww=1.0                                                            
      if(intr.ne.0) ww=yb                                               
      wr=qy-q*(yb*(1.0-uu))                                             
      a=0.0                                                             
      as=0.0                                                            
      dvr=-yb                                                           
      do 18761 i=1,no                                                   
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                          
18761 continue                                                          
      continue                                                          
      dvr=dvr-dv0                                                       
      dev0=dvr                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 18781                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
18781 continue                                                          
      m=0                                                               
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      mnl=min(mnlam,nlam)                                               
      shr=shri*dev0                                                     
      al=0.0                                                            
      ixx=0                                                             
      do 18791 j=1,ni                                                   
      if(ju(j).eq.0)goto 18791                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)
     *)-xb(j)*tt)/xs(j)
18791 continue                                                          
      continue                                                          
      do 18801 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 18821                                      
      al=ulam(ilm)                                                      
      goto 18811                                                        
18821 if(ilm .le. 2)goto 18831                                          
      al=al*alf                                                         
      goto 18811                                                        
18831 if(ilm .ne. 1)goto 18841                                          
      al=big                                                            
      goto 18851                                                        
18841 continue                                                          
      al0=0.0                                                           
      do 18861 j=1,ni                                                   
      if(ju(j).eq.0)goto 18861                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
18861 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
18851 continue                                                          
18811 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 18871 k=1,ni                                                   
      if(ixx(k).eq.1)goto 18871                                         
      if(ju(k).eq.0)goto 18871                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
18871 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
18881 continue                                                          
      az0=az                                                            
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                             
      do 18891 j=1,ni                                                   
      if(ixx(j).eq.0)goto 18891                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x
     *b(j)**2)/xs(j)**2
18891 continue                                                          
      continue                                                          
      continue                                                          
18901 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 18911 k=1,ni                                                   
      if(ixx(k).eq.0)goto 18911                                         
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      ak=a(k)                                                           
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 18931                                         
      a(k)=0.0                                                          
      goto 18941                                                        
18931 continue                                                          
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))        
18941 continue                                                          
      continue                                                          
      if(a(k).eq.ak)goto 18911                                          
      if(mm(k) .ne. 0)goto 18961                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 18912                                           
      mm(k)=nin                                                         
      m(nin)=k                                                          
18961 continue                                                          
      d=a(k)-ak                                                         
      dlx=max(dlx,v(k)*d**2)                                            
      dv=d/xs(k)                                                        
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)              
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                             
      uu=uu-dv*xb(k)                                                    
      tt=tt-dv*xm(k)                                                    
18911 continue                                                          
18912 continue                                                          
      if(nin.gt.nx)goto 18902                                           
      if(intr .eq. 0)goto 18981                                         
      d=tt/ww-uu                                                        
      az=az+d                                                           
      dlx=max(dlx,ww*d**2)                                              
      uu=uu+d                                                           
18981 continue                                                          
      if(dlx.lt.shr)goto 18902                                          
      if(nlp .le. maxit)goto 19001                                      
      jerr=-ilm                                                         
      return                                                            
19001 continue                                                          
      continue                                                          
19011 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 19021 l=1,nin                                                  
      k=m(l)                                                            
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      ak=a(k)                                                           
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 19041                                         
      a(k)=0.0                                                          
      goto 19051                                                        
19041 continue                                                          
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))        
19051 continue                                                          
      continue                                                          
      if(a(k).eq.ak)goto 19021                                          
      d=a(k)-ak                                                         
      dlx=max(dlx,v(k)*d**2)                                            
      dv=d/xs(k)                                                        
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)              
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                             
      uu=uu-dv*xb(k)                                                    
      tt=tt-dv*xm(k)                                                    
19021 continue                                                          
      continue                                                          
      if(intr .eq. 0)goto 19071                                         
      d=tt/ww-uu                                                        
      az=az+d                                                           
      dlx=max(dlx,ww*d**2)                                              
      uu=uu+d                                                           
19071 continue                                                          
      if(dlx.lt.shr)goto 19012                                          
      if(nlp .le. maxit)goto 19091                                      
      jerr=-ilm                                                         
      return                                                            
19091 continue                                                          
      goto 19011                                                        
19012 continue                                                          
      goto 18901                                                        
18902 continue                                                          
      if(nin.gt.nx)goto 18882                                           
      euu=exp(sign(min(abs(uu),fmax),uu))                               
      w=euu*q*exp(sign(min(abs(t),fmax),t))                             
      ww=sum(w)                                                         
      wr=qy-w*(1.0-uu)                                                  
      tt=sum(wr)                                                        
      if(ww*(az-az0)**2 .ge. shr)goto 19111                             
      kx=0                                                              
      do 19121 j=1,nin                                                  
      k=m(j)                                                            
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 19121                         
      kx=1                                                              
      goto 19122                                                        
19121 continue                                                          
19122 continue                                                          
      if(kx .ne. 0)goto 19141                                           
      do 19151 j=1,ni                                                   
      if(ixx(j).eq.1)goto 19151                                         
      if(ju(j).eq.0)goto 19151                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)
     *)-xb(j)*tt)/xs(j)
      if(ga(j) .le. al1*vp(j))goto 19171                                
      ixx(j)=1                                                          
      kx=1                                                              
19171 continue                                                          
19151 continue                                                          
      continue                                                          
      if(kx.eq.1) go to 10880                                           
      goto 18882                                                        
19141 continue                                                          
19111 continue                                                          
      goto 18881                                                        
18882 continue                                                          
      if(nin .le. nx)goto 19191                                         
      jerr=-10000-ilm                                                   
      goto 18802                                                        
19191 continue                                                          
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                            
      kin(ilm)=nin                                                      
      a0(ilm)=az                                                        
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                     
      if(ilm.lt.mnl)goto 18801                                          
      if(flmin.ge.1.0)goto 18801                                        
      me=0                                                              
      do 19201 j=1,nin                                                  
      if(ca(j,ilm).ne.0.0) me=me+1                                      
19201 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 18802                                            
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 18802           
      if(dev(ilm).gt.devmax)goto 18802                                  
18801 continue                                                          
18802 continue                                                          
      g=t+uu                                                            
      continue                                                          
      deallocate(t,w,wr,v,a,qy,xm,as,mm,ga,ixx)                         
      return                                                            
      end                                                               
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)    
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(n
     *lam)
      integer ix(*),jx(*)                                               
      double precision, dimension (:), allocatable :: w,f               
      if(minval(y) .ge. 0.0)goto 19221                                  
      jerr=8888                                                         
      return                                                            
19221 continue                                                          
      allocate(w(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(f(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      w=max(0d0,q)                                                      
      sw=sum(w)                                                         
      if(sw .gt. 0.0)goto 19241                                         
      jerr=9999                                                         
      go to 12180                                                       
19241 continue                                                          
      yb=dot_product(w,y)/sw                                            
      fmax=log(huge(y(1))*0.1)                                          
      do 19251 lam=1,nlam                                               
      f=a0(lam)                                                         
      do 19261 j=1,ni                                                   
      if(a(j,lam).eq.0.0)goto 19261                                     
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                       
19261 continue                                                          
      continue                                                          
      f=f+g                                                             
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                             
19251 continue                                                          
      continue                                                          
12180 continue                                                          
      deallocate(w,f)                                                   
      return                                                            
      end                                                               
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,
     *jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(
     *nlam)
      integer ix(*),jx(*),nin(nlam),ia(nx)                              
      double precision, dimension (:), allocatable :: w,f               
      if(minval(y) .ge. 0.0)goto 19281                                  
      jerr=8888                                                         
      return                                                            
19281 continue                                                          
      allocate(w(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(f(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      w=max(0d0,q)                                                      
      sw=sum(w)                                                         
      if(sw .gt. 0.0)goto 19301                                         
      jerr=9999                                                         
      go to 12180                                                       
19301 continue                                                          
      yb=dot_product(w,y)/sw                                            
      fmax=log(huge(y(1))*0.1)                                          
      do 19311 lam=1,nlam                                               
      f=a0(lam)                                                         
      do 19321 k=1,nin(lam)                                             
      j=ia(k)                                                           
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                      
19321 continue                                                          
      continue                                                          
      f=f+g                                                             
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                             
19311 continue                                                          
      continue                                                          
12180 continue                                                          
      deallocate(w,f)                                                   
      return                                                            
      end                                                               
      subroutine multelnet(parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,  flm
     *in,ulam,thr,isd,jsd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no,nr),w(no),vp(ni),ca(nx,nr,nlam)    
      double precision ulam(nlam),a0(nr,nlam),rsq(nlam),alm(nlam),cl(2,n
     *i)
      integer jd(*),ia(nx),nin(nlam)                                    
      double precision, dimension (:), allocatable :: vq;               
      if(maxval(vp) .gt. 0.0)goto 19341                                 
      jerr=10000                                                        
      return                                                            
19341 continue                                                          
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      call multelnetn(parm,no,ni,nr,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam
     *,thr,isd,  jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      deallocate(vq)                                                    
      return                                                            
      end                                                               
      subroutine multelnetn(parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flmi
     *n,ulam,thr,  isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision vp(ni),x(no,ni),y(no,nr),w(no),ulam(nlam),cl(2,ni
     *)
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)   
      integer jd(*),ia(nx),nin(nlam)                                    
      double precision, dimension (:), allocatable :: xm,xs,xv,ym,ys    
      integer, dimension (:), allocatable :: ju                         
      double precision, dimension (:,:,:), allocatable :: clt           
      allocate(clt(1:2,1:nr,1:ni),stat=jerr);                           
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ym(1:nr),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ys(1:nr),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      call chkvars(no,ni,x,ju)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 19361                                   
      jerr=7777                                                         
      return                                                            
19361 continue                                                          
      call multstandard1(no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym,ys,xv,y
     *s0,jerr)
      if(jerr.ne.0) return                                              
      do 19371 j=1,ni                                                   
      do 19381 k=1,nr                                                   
      do 19391 i=1,2                                                    
      clt(i,k,j)=cl(i,j)                                                
19391 continue                                                          
      continue                                                          
19381 continue                                                          
      continue                                                          
19371 continue                                                          
      continue                                                          
      if(isd .le. 0)goto 19411                                          
      do 19421 j=1,ni                                                   
      do 19431 k=1,nr                                                   
      do 19441 i=1,2                                                    
      clt(i,k,j)=clt(i,k,j)*xs(j)                                       
19441 continue                                                          
      continue                                                          
19431 continue                                                          
      continue                                                          
19421 continue                                                          
      continue                                                          
19411 continue                                                          
      if(jsd .le. 0)goto 19461                                          
      do 19471 j=1,ni                                                   
      do 19481 k=1,nr                                                   
      do 19491 i=1,2                                                    
      clt(i,k,j)=clt(i,k,j)/ys(k)                                       
19491 continue                                                          
      continue                                                          
19481 continue                                                          
      continue                                                          
19471 continue                                                          
      continue                                                          
19461 continue                                                          
      call multelnet2(parm,ni,nr,ju,vp,clt,y,no,ne,nx,x,nlam,flmin,ulam,
     *thr,maxit,xv,  ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                              
      do 19501 k=1,lmu                                                  
      nk=nin(k)                                                         
      do 19511 j=1,nr                                                   
      do 19521 l=1,nk                                                   
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                               
19521 continue                                                          
      continue                                                          
      if(intr .ne. 0)goto 19541                                         
      a0(j,k)=0.0                                                       
      goto 19551                                                        
19541 continue                                                          
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))              
19551 continue                                                          
      continue                                                          
19511 continue                                                          
      continue                                                          
19501 continue                                                          
      continue                                                          
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                 
      return                                                            
      end                                                               
      subroutine multstandard1(no,ni,nr,x,y,w,isd,jsd,intr,ju,  xm,xs,ym
     *,ys,xv,ys0,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(n
     *r),ys(nr)
      integer ju(ni)                                                    
      double precision, dimension (:), allocatable :: v                 
      allocate(v(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      w=w/sum(w)                                                        
      v=sqrt(w)                                                         
      if(intr .ne. 0)goto 19571                                         
      do 19581 j=1,ni                                                   
      if(ju(j).eq.0)goto 19581                                          
      xm(j)=0.0                                                         
      x(:,j)=v*x(:,j)                                                   
      z=dot_product(x(:,j),x(:,j))                                      
      if(isd .le. 0)goto 19601                                          
      xbq=dot_product(v,x(:,j))**2                                      
      vc=z-xbq                                                          
      xs(j)=sqrt(vc)                                                    
      x(:,j)=x(:,j)/xs(j)                                               
      xv(j)=1.0+xbq/vc                                                  
      goto 19611                                                        
19601 continue                                                          
      xs(j)=1.0                                                         
      xv(j)=z                                                           
19611 continue                                                          
      continue                                                          
19581 continue                                                          
      continue                                                          
      ys0=0.0                                                           
      do 19621 j=1,nr                                                   
      ym(j)=0.0                                                         
      y(:,j)=v*y(:,j)                                                   
      z=dot_product(y(:,j),y(:,j))                                      
      if(jsd .le. 0)goto 19641                                          
      u=z-dot_product(v,y(:,j))**2                                      
      ys0=ys0+z/u                                                       
      ys(j)=sqrt(u)                                                     
      y(:,j)=y(:,j)/ys(j)                                               
      goto 19651                                                        
19641 continue                                                          
      ys(j)=1.0                                                         
      ys0=ys0+z                                                         
19651 continue                                                          
      continue                                                          
19621 continue                                                          
      continue                                                          
      go to 10700                                                       
19571 continue                                                          
      do 19661 j=1,ni                                                   
      if(ju(j).eq.0)goto 19661                                          
      xm(j)=dot_product(w,x(:,j))                                       
      x(:,j)=v*(x(:,j)-xm(j))                                           
      xv(j)=dot_product(x(:,j),x(:,j))                                  
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                    
19661 continue                                                          
      continue                                                          
      if(isd .ne. 0)goto 19681                                          
      xs=1.0                                                            
      goto 19691                                                        
19681 continue                                                          
      do 19701 j=1,ni                                                   
      if(ju(j).eq.0)goto 19701                                          
      x(:,j)=x(:,j)/xs(j)                                               
19701 continue                                                          
      continue                                                          
      xv=1.0                                                            
19691 continue                                                          
      continue                                                          
      ys0=0.0                                                           
      do 19711 j=1,nr                                                   
      ym(j)=dot_product(w,y(:,j))                                       
      y(:,j)=v*(y(:,j)-ym(j))                                           
      z=dot_product(y(:,j),y(:,j))                                      
      if(jsd .le. 0)goto 19731                                          
      ys(j)=sqrt(z)                                                     
      y(:,j)=y(:,j)/ys(j)                                               
      goto 19741                                                        
19731 continue                                                          
      ys0=ys0+z                                                         
19741 continue                                                          
      continue                                                          
19711 continue                                                          
      continue                                                          
      if(jsd .ne. 0)goto 19761                                          
      ys=1.0                                                            
      goto 19771                                                        
19761 continue                                                          
      ys0=nr                                                            
19771 continue                                                          
      continue                                                          
10700 continue                                                          
      deallocate(v)                                                     
      return                                                            
      end                                                               
      subroutine multelnet2(beta,ni,nr,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,
     *ulam,thri,  maxit,xv,ys0,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision vp(ni),y(no,nr),x(no,ni),ulam(nlam),ao(nx,nr,nlam
     *)
      double precision rsqo(nlam),almo(nlam),xv(ni),cl(2,nr,ni)         
      integer ju(ni),ia(nx),kin(nlam)                                   
      double precision, dimension (:), allocatable :: g,gk,del,gj       
      integer, dimension (:), allocatable :: mm,ix,isc                  
      double precision, dimension (:,:), allocatable :: a               
      allocate(a(1:nr,1:ni),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx,itrace)     
      allocate(gj(1:nr),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(gk(1:nr),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(del(1:nr),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(g(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(ix(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(isc(1:nr),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      bta=beta                                                          
      omb=1.0-bta                                                       
      ix=0                                                              
      thr=thri*ys0/nr                                                   
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 19791                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
19791 continue                                                          
      rsq=ys0                                                           
      a=0.0                                                             
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      iz=0                                                              
      mnl=min(mnlam,nlam)                                               
      alm=0.0                                                           
      do 19801 j=1,ni                                                   
      if(ju(j).eq.0)goto 19801                                          
      g(j)=0.0                                                          
      do 19811 k=1,nr                                                   
      g(j)=g(j)+dot_product(y(:,k),x(:,j))**2                           
19811 continue                                                          
      continue                                                          
      g(j)=sqrt(g(j))                                                   
19801 continue                                                          
      continue                                                          
      do 19821 m=1,nlam                                                 
      if(itrace.ne.0) call setpb(m-1)                                   
      alm0=alm                                                          
      if(flmin .lt. 1.0)goto 19841                                      
      alm=ulam(m)                                                       
      goto 19831                                                        
19841 if(m .le. 2)goto 19851                                            
      alm=alm*alf                                                       
      goto 19831                                                        
19851 if(m .ne. 1)goto 19861                                            
      alm=big                                                           
      goto 19871                                                        
19861 continue                                                          
      alm0=0.0                                                          
      do 19881 j=1,ni                                                   
      if(ju(j).eq.0)goto 19881                                          
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                        
19881 continue                                                          
      continue                                                          
      alm0=alm0/max(bta,1.0d-3)                                         
      alm=alf*alm0                                                      
19871 continue                                                          
19831 continue                                                          
      dem=alm*omb                                                       
      ab=alm*bta                                                        
      rsq0=rsq                                                          
      jz=1                                                              
      tlam=bta*(2.0*alm-alm0)                                           
      do 19891 k=1,ni                                                   
      if(ix(k).eq.1)goto 19891                                          
      if(ju(k).eq.0)goto 19891                                          
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                    
19891 continue                                                          
      continue                                                          
      continue                                                          
19901 continue                                                          
      if(iz*jz.ne.0) go to 10360                                        
10880 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 19911 k=1,ni                                                   
      if(ix(k).eq.0)goto 19911                                          
      gkn=0.0                                                           
      do 19921 j=1,nr                                                   
      gj(j)=dot_product(y(:,j),x(:,k))                                  
      gk(j)=gj(j)+a(j,k)*xv(k)                                          
      gkn=gkn+gk(j)**2                                                  
19921 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-ab*vp(k)/gkn                                                
      del=a(:,k)                                                        
      if(u .gt. 0.0)goto 19941                                          
      a(:,k)=0.0                                                        
      goto 19951                                                        
19941 continue                                                          
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                   
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)
     *,isc,jerr)
      if(jerr.ne.0) return                                              
19951 continue                                                          
      continue                                                          
      del=a(:,k)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 19911                             
      do 19961 j=1,nr                                                   
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                           
      y(:,j)=y(:,j)-del(j)*x(:,k)                                       
      dlx=max(dlx,xv(k)*del(j)**2)                                      
19961 continue                                                          
      continue                                                          
      if(mm(k) .ne. 0)goto 19981                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 19912                                           
      mm(k)=nin                                                         
      ia(nin)=k                                                         
19981 continue                                                          
19911 continue                                                          
19912 continue                                                          
      if(nin.gt.nx)goto 19902                                           
      if(dlx .ge. thr)goto 20001                                        
      ixx=0                                                             
      do 20011 k=1,ni                                                   
      if(ix(k).eq.1)goto 20011                                          
      if(ju(k).eq.0)goto 20011                                          
      g(k)=0.0                                                          
      do 20021 j=1,nr                                                   
      g(k)=g(k)+dot_product(y(:,j),x(:,k))**2                           
20021 continue                                                          
      continue                                                          
      g(k)=sqrt(g(k))                                                   
      if(g(k) .le. ab*vp(k))goto 20041                                  
      ix(k)=1                                                           
      ixx=1                                                             
20041 continue                                                          
20011 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10880                                          
      goto 19902                                                        
20001 continue                                                          
      if(nlp .le. maxit)goto 20061                                      
      jerr=-m                                                           
      return                                                            
20061 continue                                                          
10360 continue                                                          
      iz=1                                                              
      continue                                                          
20071 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 20081 l=1,nin                                                  
      k=ia(l)                                                           
      gkn=0.0                                                           
      do 20091 j=1,nr                                                   
      gj(j)=dot_product(y(:,j),x(:,k))                                  
      gk(j)=gj(j)+a(j,k)*xv(k)                                          
      gkn=gkn+gk(j)**2                                                  
20091 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-ab*vp(k)/gkn                                                
      del=a(:,k)                                                        
      if(u .gt. 0.0)goto 20111                                          
      a(:,k)=0.0                                                        
      goto 20121                                                        
20111 continue                                                          
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                   
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)
     *,isc,jerr)
      if(jerr.ne.0) return                                              
20121 continue                                                          
      continue                                                          
      del=a(:,k)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 20081                             
      do 20131 j=1,nr                                                   
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                           
      y(:,j)=y(:,j)-del(j)*x(:,k)                                       
      dlx=max(dlx,xv(k)*del(j)**2)                                      
20131 continue                                                          
      continue                                                          
20081 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 20072                                          
      if(nlp .le. maxit)goto 20151                                      
      jerr=-m                                                           
      return                                                            
20151 continue                                                          
      goto 20071                                                        
20072 continue                                                          
      jz=0                                                              
      goto 19901                                                        
19902 continue                                                          
      if(nin .le. nx)goto 20171                                         
      jerr=-10000-m                                                     
      goto 19822                                                        
20171 continue                                                          
      if(nin .le. 0)goto 20191                                          
      do 20201 j=1,nr                                                   
      ao(1:nin,j,m)=a(j,ia(1:nin))                                      
20201 continue                                                          
      continue                                                          
20191 continue                                                          
      kin(m)=nin                                                        
      rsqo(m)=1.0-rsq/ys0                                               
      almo(m)=alm                                                       
      lmu=m                                                             
      if(m.lt.mnl)goto 19821                                            
      if(flmin.ge.1.0)goto 19821                                        
      me=0                                                              
      do 20211 j=1,nin                                                  
      if(ao(j,1,m).ne.0.0) me=me+1                                      
20211 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 19822                                            
      if(rsq0-rsq.lt.sml*rsq)goto 19822                                 
      if(rsqo(m).gt.rsqmax)goto 19822                                   
19821 continue                                                          
19822 continue                                                          
      deallocate(a,mm,g,ix,del,gj,gk)                                   
      return                                                            
      end                                                               
      subroutine chkbnds(nr,gk,gkn,xv,cl,al1,al2,a,isc,jerr)            
      implicit double precision(a-h,o-z)                                
      double precision gk(nr),cl(2,nr),a(nr)                            
      integer isc(nr)                                                   
      kerr=0                                                            
      al1p=1.0+al1/xv                                                   
      al2p=al2/xv                                                       
      isc=0                                                             
      gsq=gkn**2                                                        
      asq=dot_product(a,a)                                              
      usq=0.0                                                           
      u=0.0                                                             
      kn=-1                                                             
      continue                                                          
20221 continue                                                          
      vmx=0.0                                                           
      do 20231 k=1,nr                                                   
      v=max(a(k)-cl(2,k),cl(1,k)-a(k))                                  
      if(v .le. vmx)goto 20251                                          
      vmx=v                                                             
      kn=k                                                              
20251 continue                                                          
20231 continue                                                          
      continue                                                          
      if(vmx.le.0.0)goto 20222                                          
      if(isc(kn).ne.0)goto 20222                                        
      gsq=gsq-gk(kn)**2                                                 
      g=sqrt(gsq)/xv                                                    
      if(a(kn).lt.cl(1,kn)) u=cl(1,kn)                                  
      if(a(kn).gt.cl(2,kn)) u=cl(2,kn)                                  
      usq=usq+u**2                                                      
      if(usq .ne. 0.0)goto 20271                                        
      b=max(0d0,(g-al2p)/al1p)                                          
      goto 20281                                                        
20271 continue                                                          
      b0=sqrt(asq-a(kn)**2)                                             
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                  
      if(kerr.ne.0)goto 20222                                           
20281 continue                                                          
      continue                                                          
      asq=usq+b**2                                                      
      if(asq .gt. 0.0)goto 20301                                        
      a=0.0                                                             
      goto 20222                                                        
20301 continue                                                          
      a(kn)=u                                                           
      isc(kn)=1                                                         
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                  
      do 20311 j=1,nr                                                   
      if(isc(j).eq.0) a(j)=f*gk(j)                                      
20311 continue                                                          
      continue                                                          
      goto 20221                                                        
20222 continue                                                          
      if(kerr.ne.0) jerr=kerr                                           
      return                                                            
      end                                                               
      subroutine chkbnds1(nr,gk,gkn,xv,cl1,cl2,al1,al2,a,isc,jerr)      
      implicit double precision(a-h,o-z)                                
      double precision gk(nr),a(nr)                                     
      integer isc(nr)                                                   
      kerr=0                                                            
      al1p=1.0+al1/xv                                                   
      al2p=al2/xv                                                       
      isc=0                                                             
      gsq=gkn**2                                                        
      asq=dot_product(a,a)                                              
      usq=0.0                                                           
      u=0.0                                                             
      kn=-1                                                             
      continue                                                          
20321 continue                                                          
      vmx=0.0                                                           
      do 20331 k=1,nr                                                   
      v=max(a(k)-cl2,cl1-a(k))                                          
      if(v .le. vmx)goto 20351                                          
      vmx=v                                                             
      kn=k                                                              
20351 continue                                                          
20331 continue                                                          
      continue                                                          
      if(vmx.le.0.0)goto 20322                                          
      if(isc(kn).ne.0)goto 20322                                        
      gsq=gsq-gk(kn)**2                                                 
      g=sqrt(gsq)/xv                                                    
      if(a(kn).lt.cl1) u=cl1                                            
      if(a(kn).gt.cl2) u=cl2                                            
      usq=usq+u**2                                                      
      if(usq .ne. 0.0)goto 20371                                        
      b=max(0d0,(g-al2p)/al1p)                                          
      goto 20381                                                        
20371 continue                                                          
      b0=sqrt(asq-a(kn)**2)                                             
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                  
      if(kerr.ne.0)goto 20322                                           
20381 continue                                                          
      continue                                                          
      asq=usq+b**2                                                      
      if(asq .gt. 0.0)goto 20401                                        
      a=0.0                                                             
      goto 20322                                                        
20401 continue                                                          
      a(kn)=u                                                           
      isc(kn)=1                                                         
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                  
      do 20411 j=1,nr                                                   
      if(isc(j).eq.0) a(j)=f*gk(j)                                      
20411 continue                                                          
      continue                                                          
      goto 20321                                                        
20322 continue                                                          
      if(kerr.ne.0) jerr=kerr                                           
      return                                                            
      end                                                               
      function bnorm(b0,al1p,al2p,g,usq,jerr)                           
      implicit double precision(a-h,o-z)                                
      data thr,mxit /1.0d-10,100/                                       
      b=b0                                                              
      zsq=b**2+usq                                                      
      if(zsq .gt. 0.0)goto 20431                                        
      bnorm=0.0                                                         
      return                                                            
20431 continue                                                          
      z=sqrt(zsq)                                                       
      f=b*(al1p+al2p/z)-g                                               
      jerr=0                                                            
      do 20441 it=1,mxit                                                
      b=b-f/(al1p+al2p*usq/(z*zsq))                                     
      zsq=b**2+usq                                                      
      if(zsq .gt. 0.0)goto 20461                                        
      bnorm=0.0                                                         
      return                                                            
20461 continue                                                          
      z=sqrt(zsq)                                                       
      f=b*(al1p+al2p/z)-g                                               
      if(abs(f).le.thr)goto 20442                                       
      if(b .gt. 0.0)goto 20481                                          
      b=0.0                                                             
      goto 20442                                                        
20481 continue                                                          
20441 continue                                                          
20442 continue                                                          
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
      subroutine multsolns(ni,nx,nr,lmu,a,ia,nin,b)                     
      implicit double precision(a-h,o-z)                                
      double precision a(nx,nr,lmu),b(ni,nr,lmu)                        
      integer ia(nx),nin(lmu)                                           
      do 20491 lam=1,lmu                                                
      call multuncomp(ni,nr,nx,a(1,1,lam),ia,nin(lam),b(1,1,lam))       
20491 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine multuncomp(ni,nr,nx,ca,ia,nin,a)                       
      implicit double precision(a-h,o-z)                                
      double precision ca(nx,nr),a(ni,nr)                               
      integer ia(nx)                                                    
      a=0.0                                                             
      if(nin .le. 0)goto 20511                                          
      do 20521 j=1,nr                                                   
      a(ia(1:nin),j)=ca(1:nin,j)                                        
20521 continue                                                          
      continue                                                          
20511 continue                                                          
      return                                                            
      end                                                               
      subroutine multmodval(nx,nr,a0,ca,ia,nin,n,x,f)                   
      implicit double precision(a-h,o-z)                                
      double precision a0(nr),ca(nx,nr),x(n,*),f(nr,n)                  
      integer ia(nx)                                                    
      do 20531 i=1,n                                                    
      f(:,i)=a0                                                         
20531 continue                                                          
      continue                                                          
      if(nin.le.0) return                                               
      do 20541 i=1,n                                                    
      do 20551 j=1,nr                                                   
      f(j,i)=f(j,i)+dot_product(ca(1:nin,j),x(i,ia(1:nin)))             
20551 continue                                                          
      continue                                                          
20541 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine multspelnet(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,  
     *nlam,flmin,ulam,thr,isd,  jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,
     *nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no,nr),w(no),vp(ni),ulam(nlam),cl(2,ni)   
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)   
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                        
      double precision, dimension (:), allocatable :: vq;               
      if(maxval(vp) .gt. 0.0)goto 20571                                 
      jerr=10000                                                        
      return                                                            
20571 continue                                                          
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      call multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,fl
     *min,  ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jer
     *r)
      deallocate(vq)                                                    
      return                                                            
      end                                                               
      subroutine multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,n
     *lam,flmin,  ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,n
     *lp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),vp(ni),y(no,nr),w(no),ulam(nlam),cl(2,ni)   
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)   
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                        
      double precision, dimension (:), allocatable :: xm,xs,xv,ym,ys    
      integer, dimension (:), allocatable :: ju                         
      double precision, dimension (:,:,:), allocatable :: clt           
      allocate(clt(1:2,1:nr,1:ni),stat=jerr)                            
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ym(1:nr),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ys(1:nr),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      call spchkvars(no,ni,x,ix,ju)                                     
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 20591                                   
      jerr=7777                                                         
      return                                                            
20591 continue                                                          
      call multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,  xm,xs,
     *ym,ys,xv,ys0,jerr)
      if(jerr.ne.0) return                                              
      do 20601 j=1,ni                                                   
      do 20611 k=1,nr                                                   
      do 20621 i=1,2                                                    
      clt(i,k,j)=cl(i,j)                                                
20621 continue                                                          
      continue                                                          
20611 continue                                                          
      continue                                                          
20601 continue                                                          
      continue                                                          
      if(isd .le. 0)goto 20641                                          
      do 20651 j=1,ni                                                   
      do 20661 k=1,nr                                                   
      do 20671 i=1,2                                                    
      clt(i,k,j)=clt(i,k,j)*xs(j)                                       
20671 continue                                                          
      continue                                                          
20661 continue                                                          
      continue                                                          
20651 continue                                                          
      continue                                                          
20641 continue                                                          
      if(jsd .le. 0)goto 20691                                          
      do 20701 j=1,ni                                                   
      do 20711 k=1,nr                                                   
      do 20721 i=1,2                                                    
      clt(i,k,j)=clt(i,k,j)/ys(k)                                       
20721 continue                                                          
      continue                                                          
20711 continue                                                          
      continue                                                          
20701 continue                                                          
      continue                                                          
20691 continue                                                          
      call multspelnet2(parm,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,clt,nlam,f
     *lmin,  ulam,thr,maxit,xm,xs,xv,ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                              
      do 20731 k=1,lmu                                                  
      nk=nin(k)                                                         
      do 20741 j=1,nr                                                   
      do 20751 l=1,nk                                                   
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                               
20751 continue                                                          
      continue                                                          
      if(intr .ne. 0)goto 20771                                         
      a0(j,k)=0.0                                                       
      goto 20781                                                        
20771 continue                                                          
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))              
20781 continue                                                          
      continue                                                          
20741 continue                                                          
      continue                                                          
20731 continue                                                          
      continue                                                          
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                 
      return                                                            
      end                                                               
      subroutine multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,  
     *xm,xs,ym,ys,xv,ys0,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),y
     *s(nr)
      integer ix(*),jx(*),ju(ni)                                        
      jerr = jerr*1                                                     
      w=w/sum(w)                                                        
      if(intr .ne. 0)goto 20801                                         
      do 20811 j=1,ni                                                   
      if(ju(j).eq.0)goto 20811                                          
      xm(j)=0.0                                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      z=dot_product(w(jx(jb:je)),x(jb:je)**2)                           
      if(isd .le. 0)goto 20831                                          
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                         
      vc=z-xbq                                                          
      xs(j)=sqrt(vc)                                                    
      xv(j)=1.0+xbq/vc                                                  
      goto 20841                                                        
20831 continue                                                          
      xs(j)=1.0                                                         
      xv(j)=z                                                           
20841 continue                                                          
      continue                                                          
20811 continue                                                          
      continue                                                          
      ys0=0.0                                                           
      do 20851 j=1,nr                                                   
      ym(j)=0.0                                                         
      z=dot_product(w,y(:,j)**2)                                        
      if(jsd .le. 0)goto 20871                                          
      u=z-dot_product(w,y(:,j))**2                                      
      ys0=ys0+z/u                                                       
      ys(j)=sqrt(u)                                                     
      y(:,j)=y(:,j)/ys(j)                                               
      goto 20881                                                        
20871 continue                                                          
      ys(j)=1.0                                                         
      ys0=ys0+z                                                         
20881 continue                                                          
      continue                                                          
20851 continue                                                          
      continue                                                          
      return                                                            
20801 continue                                                          
      do 20891 j=1,ni                                                   
      if(ju(j).eq.0)goto 20891                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2              
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                    
20891 continue                                                          
      continue                                                          
      if(isd .ne. 0)goto 20911                                          
      xs=1.0                                                            
      goto 20921                                                        
20911 continue                                                          
      xv=1.0                                                            
20921 continue                                                          
      continue                                                          
      ys0=0.0                                                           
      do 20931 j=1,nr                                                   
      ym(j)=dot_product(w,y(:,j))                                       
      y(:,j)=y(:,j)-ym(j)                                               
      z=dot_product(w,y(:,j)**2)                                        
      if(jsd .le. 0)goto 20951                                          
      ys(j)=sqrt(z)                                                     
      y(:,j)=y(:,j)/ys(j)                                               
      goto 20961                                                        
20951 continue                                                          
      ys0=ys0+z                                                         
20961 continue                                                          
      continue                                                          
20931 continue                                                          
      continue                                                          
      if(jsd .ne. 0)goto 20981                                          
      ys=1.0                                                            
      goto 20991                                                        
20981 continue                                                          
      ys0=nr                                                            
20991 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine multspelnet2(beta,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,n
     *lam,flmin,  ulam,thri,maxit,xm,xs,xv,ys0,lmu,ao,ia,kin,rsqo,almo,n
     *lp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision y(no,nr),w(no),x(*),vp(ni),ulam(nlam),cl(2,nr,ni)
      double precision ao(nx,nr,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni
     *),xv(ni)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                       
      double precision, dimension (:), allocatable :: g,gj,gk,del,o     
      integer, dimension (:), allocatable :: mm,iy,isc                  
      double precision, dimension (:,:), allocatable :: a               
      allocate(a(1:nr,1:ni),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx,itrace)     
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(g(1:ni),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(gj(1:nr),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(gk(1:nr),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(del(1:nr),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(o(1:nr),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(iy(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(isc(1:nr),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      bta=beta                                                          
      omb=1.0-bta                                                       
      alm=0.0                                                           
      iy=0                                                              
      thr=thri*ys0/nr                                                   
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 21011                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
21011 continue                                                          
      rsq=ys0                                                           
      a=0.0                                                             
      mm=0                                                              
      o=0.0                                                             
      nlp=0                                                             
      nin=nlp                                                           
      iz=0                                                              
      mnl=min(mnlam,nlam)                                               
      do 21021 j=1,ni                                                   
      if(ju(j).eq.0)goto 21021                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      g(j)=0.0                                                          
      do 21031 k=1,nr                                                   
      g(j)=g(j)+(dot_product(y(jx(jb:je),k),w(jx(jb:je))*x(jb:je))/xs(j)
     *)**2
21031 continue                                                          
      continue                                                          
      g(j)=sqrt(g(j))                                                   
21021 continue                                                          
      continue                                                          
      do 21041 m=1,nlam                                                 
      if(itrace.ne.0) call setpb(m-1)                                   
      alm0=alm                                                          
      if(flmin .lt. 1.0)goto 21061                                      
      alm=ulam(m)                                                       
      goto 21051                                                        
21061 if(m .le. 2)goto 21071                                            
      alm=alm*alf                                                       
      goto 21051                                                        
21071 if(m .ne. 1)goto 21081                                            
      alm=big                                                           
      goto 21091                                                        
21081 continue                                                          
      alm0=0.0                                                          
      do 21101 j=1,ni                                                   
      if(ju(j).eq.0)goto 21101                                          
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                        
21101 continue                                                          
      continue                                                          
      alm0=alm0/max(bta,1.0d-3)                                         
      alm=alf*alm0                                                      
21091 continue                                                          
21051 continue                                                          
      dem=alm*omb                                                       
      ab=alm*bta                                                        
      rsq0=rsq                                                          
      jz=1                                                              
      tlam=bta*(2.0*alm-alm0)                                           
      do 21111 k=1,ni                                                   
      if(iy(k).eq.1)goto 21111                                          
      if(ju(k).eq.0)goto 21111                                          
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                    
21111 continue                                                          
      continue                                                          
      continue                                                          
21121 continue                                                          
      if(iz*jz.ne.0) go to 10360                                        
10880 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 21131 k=1,ni                                                   
      if(iy(k).eq.0)goto 21131                                          
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      gkn=0.0                                                           
      do 21141 j=1,nr                                                   
      gj(j)=dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(k)
      gk(j)=gj(j)+a(j,k)*xv(k)                                          
      gkn=gkn+gk(j)**2                                                  
21141 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-ab*vp(k)/gkn                                                
      del=a(:,k)                                                        
      if(u .gt. 0.0)goto 21161                                          
      a(:,k)=0.0                                                        
      goto 21171                                                        
21161 continue                                                          
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                   
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)
     *,isc,jerr)
      if(jerr.ne.0) return                                              
21171 continue                                                          
      continue                                                          
      del=a(:,k)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 21131                             
      if(mm(k) .ne. 0)goto 21191                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 21132                                           
      mm(k)=nin                                                         
      ia(nin)=k                                                         
21191 continue                                                          
      do 21201 j=1,nr                                                   
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                           
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)               
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                      
      dlx=max(xv(k)*del(j)**2,dlx)                                      
21201 continue                                                          
      continue                                                          
21131 continue                                                          
21132 continue                                                          
      if(nin.gt.nx)goto 21122                                           
      if(dlx .ge. thr)goto 21221                                        
      ixx=0                                                             
      do 21231 j=1,ni                                                   
      if(iy(j).eq.1)goto 21231                                          
      if(ju(j).eq.0)goto 21231                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      g(j)=0.0                                                          
      do 21241 k=1,nr                                                   
      g(j)=g(j)+  (dot_product(y(jx(jb:je),k)+o(k),w(jx(jb:je))*x(jb:je)
     *)/xs(j))**2
21241 continue                                                          
      continue                                                          
      g(j)=sqrt(g(j))                                                   
      if(g(j) .le. ab*vp(j))goto 21261                                  
      iy(j)=1                                                           
      ixx=1                                                             
21261 continue                                                          
21231 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10880                                          
      goto 21122                                                        
21221 continue                                                          
      if(nlp .le. maxit)goto 21281                                      
      jerr=-m                                                           
      return                                                            
21281 continue                                                          
10360 continue                                                          
      iz=1                                                              
      continue                                                          
21291 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 21301 l=1,nin                                                  
      k=ia(l)                                                           
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      gkn=0.0                                                           
      do 21311 j=1,nr                                                   
      gj(j)=  dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(
     *k)
      gk(j)=gj(j)+a(j,k)*xv(k)                                          
      gkn=gkn+gk(j)**2                                                  
21311 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-ab*vp(k)/gkn                                                
      del=a(:,k)                                                        
      if(u .gt. 0.0)goto 21331                                          
      a(:,k)=0.0                                                        
      goto 21341                                                        
21331 continue                                                          
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                   
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)
     *,isc,jerr)
      if(jerr.ne.0) return                                              
21341 continue                                                          
      continue                                                          
      del=a(:,k)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 21301                             
      do 21351 j=1,nr                                                   
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                           
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)               
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                      
      dlx=max(xv(k)*del(j)**2,dlx)                                      
21351 continue                                                          
      continue                                                          
21301 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 21292                                          
      if(nlp .le. maxit)goto 21371                                      
      jerr=-m                                                           
      return                                                            
21371 continue                                                          
      goto 21291                                                        
21292 continue                                                          
      jz=0                                                              
      goto 21121                                                        
21122 continue                                                          
      if(nin .le. nx)goto 21391                                         
      jerr=-10000-m                                                     
      goto 21042                                                        
21391 continue                                                          
      if(nin .le. 0)goto 21411                                          
      do 21421 j=1,nr                                                   
      ao(1:nin,j,m)=a(j,ia(1:nin))                                      
21421 continue                                                          
      continue                                                          
21411 continue                                                          
      kin(m)=nin                                                        
      rsqo(m)=1.0-rsq/ys0                                               
      almo(m)=alm                                                       
      lmu=m                                                             
      if(m.lt.mnl)goto 21041                                            
      if(flmin.ge.1.0)goto 21041                                        
      me=0                                                              
      do 21431 j=1,nin                                                  
      if(ao(j,1,m).ne.0.0) me=me+1                                      
21431 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 21042                                            
      if(rsq0-rsq.lt.sml*rsq)goto 21042                                 
      if(rsqo(m).gt.rsqmax)goto 21042                                   
21041 continue                                                          
21042 continue                                                          
      deallocate(a,mm,g,iy,gj,gk,del,o)                                 
      return                                                            
      end                                                               
      subroutine multlognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,f
     *lmin,ulam,  shri,intr,maxit,xv,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam
     *),cl(2,ni)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),xv(
     *ni)
      integer ju(ni),m(nx),kin(nlam)                                    
      double precision, dimension (:,:), allocatable :: q,r,b,bs        
      double precision, dimension (:), allocatable :: sxp,sxpl,ga,gk,del
      integer, dimension (:), allocatable :: mm,is,ixx,isc              
      allocate(b(0:ni,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      allocate(bs(0:ni,1:nc),stat=jerr)                                 
      if(jerr.ne.0) return                                              
      allocate(q(1:no,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      allocate(r(1:no,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx,itrace)     
      exmn=-exmx                                                        
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(is(1:max(nc,ni)),stat=jerr)                              
      if(jerr.ne.0) return                                              
      allocate(sxp(1:no),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(sxpl(1:no),stat=jerr)                                    
      if(jerr.ne.0) return                                              
      allocate(ga(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ixx(1:ni),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(gk(1:nc),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(del(1:nc),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(isc(1:nc),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      pmax=1.0-pmin                                                     
      emin=pmin/pmax                                                    
      emax=1.0/emin                                                     
      bta=parm                                                          
      omb=1.0-bta                                                       
      dev1=0.0                                                          
      dev0=0.0                                                          
      do 21441 ic=1,nc                                                  
      q0=dot_product(w,y(:,ic))                                         
      if(q0 .gt. pmin)goto 21461                                        
      jerr =8000+ic                                                     
      return                                                            
21461 continue                                                          
      if(q0 .lt. pmax)goto 21481                                        
      jerr =9000+ic                                                     
      return                                                            
21481 continue                                                          
      if(intr .ne. 0)goto 21501                                         
      q0=1.0/nc                                                         
      b(0,ic)=0.0                                                       
      goto 21511                                                        
21501 continue                                                          
      b(0,ic)=log(q0)                                                   
      dev1=dev1-q0*b(0,ic)                                              
21511 continue                                                          
      continue                                                          
      b(1:ni,ic)=0.0                                                    
21441 continue                                                          
      continue                                                          
      if(intr.eq.0) dev1=log(float(nc))                                 
      ixx=0                                                             
      al=0.0                                                            
      if(nonzero(no*nc,g) .ne. 0)goto 21531                             
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                      
      sxp=0.0                                                           
      do 21541 ic=1,nc                                                  
      q(:,ic)=exp(b(0,ic))                                              
      sxp=sxp+q(:,ic)                                                   
21541 continue                                                          
      continue                                                          
      goto 21551                                                        
21531 continue                                                          
      do 21561 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
21561 continue                                                          
      continue                                                          
      sxp=0.0                                                           
      if(intr .ne. 0)goto 21581                                         
      b(0,:)=0.0                                                        
      goto 21591                                                        
21581 continue                                                          
      call kazero(nc,no,y,g,w,b(0,:),jerr)                              
      if(jerr.ne.0) return                                              
21591 continue                                                          
      continue                                                          
      dev1=0.0                                                          
      do 21601 ic=1,nc                                                  
      q(:,ic)=b(0,ic)+g(:,ic)                                           
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                          
      q(:,ic)=exp(q(:,ic))                                              
      sxp=sxp+q(:,ic)                                                   
21601 continue                                                          
      continue                                                          
      sxpl=w*log(sxp)                                                   
      do 21611 ic=1,nc                                                  
      dev1=dev1+dot_product(y(:,ic),sxpl)                               
21611 continue                                                          
      continue                                                          
21551 continue                                                          
      continue                                                          
      do 21621 ic=1,nc                                                  
      do 21631 i=1,no                                                   
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))            
21631 continue                                                          
      continue                                                          
21621 continue                                                          
      continue                                                          
      dev0=dev0+dev1                                                    
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 21651                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
21651 continue                                                          
      m=0                                                               
      mm=0                                                              
      nin=0                                                             
      nlp=0                                                             
      mnl=min(mnlam,nlam)                                               
      bs=0.0                                                            
      shr=shri*dev0                                                     
      ga=0.0                                                            
      do 21661 ic=1,nc                                                  
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                   
      do 21671 j=1,ni                                                   
      if(ju(j).ne.0) ga(j)=ga(j)+dot_product(r(:,ic),x(:,j))**2         
21671 continue                                                          
      continue                                                          
21661 continue                                                          
      continue                                                          
      ga=sqrt(ga)                                                       
      do 21681 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 21701                                      
      al=ulam(ilm)                                                      
      goto 21691                                                        
21701 if(ilm .le. 2)goto 21711                                          
      al=al*alf                                                         
      goto 21691                                                        
21711 if(ilm .ne. 1)goto 21721                                          
      al=big                                                            
      goto 21731                                                        
21721 continue                                                          
      al0=0.0                                                           
      do 21741 j=1,ni                                                   
      if(ju(j).eq.0)goto 21741                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
21741 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
21731 continue                                                          
21691 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 21751 k=1,ni                                                   
      if(ixx(k).eq.1)goto 21751                                         
      if(ju(k).eq.0)goto 21751                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
21751 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
21761 continue                                                          
      ix=0                                                              
      jx=ix                                                             
      kx=jx                                                             
      t=0.0                                                             
      do 21771 ic=1,nc                                                  
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                    
21771 continue                                                          
      continue                                                          
      if(t .ge. eps)goto 21791                                          
      kx=1                                                              
      goto 21762                                                        
21791 continue                                                          
      t=2.0*t                                                           
      alt=al1/t                                                         
      al2t=al2/t                                                        
      do 21801 ic=1,nc                                                  
      bs(0,ic)=b(0,ic)                                                  
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                       
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                 
      d=0.0                                                             
      if(intr.ne.0) d=sum(r(:,ic))                                      
      if(d .eq. 0.0)goto 21821                                          
      b(0,ic)=b(0,ic)+d                                                 
      r(:,ic)=r(:,ic)-d*w                                               
      dlx=max(dlx,d**2)                                                 
21821 continue                                                          
21801 continue                                                          
      continue                                                          
      continue                                                          
21831 continue                                                          
      nlp=nlp+nc                                                        
      dlx=0.0                                                           
      do 21841 k=1,ni                                                   
      if(ixx(k).eq.0)goto 21841                                         
      gkn=0.0                                                           
      do 21851 ic=1,nc                                                  
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                  
      gkn=gkn+gk(ic)**2                                                 
21851 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-alt*vp(k)/gkn                                               
      del=b(k,:)                                                        
      if(u .gt. 0.0)goto 21871                                          
      b(k,:)=0.0                                                        
      goto 21881                                                        
21871 continue                                                          
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                  
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                              
21881 continue                                                          
      continue                                                          
      del=b(k,:)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 21841                             
      do 21891 ic=1,nc                                                  
      dlx=max(dlx,xv(k)*del(ic)**2)                                     
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                  
21891 continue                                                          
      continue                                                          
      if(mm(k) .ne. 0)goto 21911                                        
      nin=nin+1                                                         
      if(nin .le. nx)goto 21931                                         
      jx=1                                                              
      goto 21842                                                        
21931 continue                                                          
      mm(k)=nin                                                         
      m(nin)=k                                                          
21911 continue                                                          
21841 continue                                                          
21842 continue                                                          
      if(jx.gt.0)goto 21832                                             
      if(dlx.lt.shr)goto 21832                                          
      if(nlp .le. maxit)goto 21951                                      
      jerr=-ilm                                                         
      return                                                            
21951 continue                                                          
      continue                                                          
21961 continue                                                          
      nlp=nlp+nc                                                        
      dlx=0.0                                                           
      do 21971 l=1,nin                                                  
      k=m(l)                                                            
      gkn=0.0                                                           
      do 21981 ic=1,nc                                                  
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                  
      gkn=gkn+gk(ic)**2                                                 
21981 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-alt*vp(k)/gkn                                               
      del=b(k,:)                                                        
      if(u .gt. 0.0)goto 22001                                          
      b(k,:)=0.0                                                        
      goto 22011                                                        
22001 continue                                                          
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                  
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                              
22011 continue                                                          
      continue                                                          
      del=b(k,:)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 21971                             
      do 22021 ic=1,nc                                                  
      dlx=max(dlx,xv(k)*del(ic)**2)                                     
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                  
22021 continue                                                          
      continue                                                          
21971 continue                                                          
      continue                                                          
      if(dlx.lt.shr)goto 21962                                          
      if(nlp .le. maxit)goto 22041                                      
      jerr=-ilm                                                         
      return                                                            
22041 continue                                                          
      goto 21961                                                        
21962 continue                                                          
      goto 21831                                                        
21832 continue                                                          
      if(jx.gt.0)goto 21762                                             
      do 22051 ic=1,nc                                                  
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                             
      if(ix .ne. 0)goto 22071                                           
      do 22081 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 22101                
      ix=1                                                              
      goto 22082                                                        
22101 continue                                                          
22081 continue                                                          
22082 continue                                                          
22071 continue                                                          
      do 22111 i=1,no                                                   
      fi=b(0,ic)+g(i,ic)                                                
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))      
      fi=min(max(exmn,fi),exmx)                                         
      sxp(i)=sxp(i)-q(i,ic)                                             
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                 
      sxp(i)=sxp(i)+q(i,ic)                                             
22111 continue                                                          
      continue                                                          
22051 continue                                                          
      continue                                                          
      s=-sum(b(0,:))/nc                                                 
      b(0,:)=b(0,:)+s                                                   
      if(jx.gt.0)goto 21762                                             
      if(ix .ne. 0)goto 22131                                           
      do 22141 k=1,ni                                                   
      if(ixx(k).eq.1)goto 22141                                         
      if(ju(k).eq.0)goto 22141                                          
      ga(k)=0.0                                                         
22141 continue                                                          
      continue                                                          
      do 22151 ic=1,nc                                                  
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                   
      do 22161 k=1,ni                                                   
      if(ixx(k).eq.1)goto 22161                                         
      if(ju(k).eq.0)goto 22161                                          
      ga(k)=ga(k)+dot_product(r(:,ic),x(:,k))**2                        
22161 continue                                                          
      continue                                                          
22151 continue                                                          
      continue                                                          
      ga=sqrt(ga)                                                       
      do 22171 k=1,ni                                                   
      if(ixx(k).eq.1)goto 22171                                         
      if(ju(k).eq.0)goto 22171                                          
      if(ga(k) .le. al1*vp(k))goto 22191                                
      ixx(k)=1                                                          
      ix=1                                                              
22191 continue                                                          
22171 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10880                                           
      goto 21762                                                        
22131 continue                                                          
      goto 21761                                                        
21762 continue                                                          
      if(kx .le. 0)goto 22211                                           
      jerr=-20000-ilm                                                   
      goto 21682                                                        
22211 continue                                                          
      if(jx .le. 0)goto 22231                                           
      jerr=-10000-ilm                                                   
      goto 21682                                                        
22231 continue                                                          
      devi=0.0                                                          
      do 22241 ic=1,nc                                                  
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                       
      a0(ic,ilm)=b(0,ic)                                                
      do 22251 i=1,no                                                   
      if(y(i,ic).le.0.0)goto 22251                                      
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                        
22251 continue                                                          
      continue                                                          
22241 continue                                                          
      continue                                                          
      kin(ilm)=nin                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dev1-devi)/dev0                                         
      if(ilm.lt.mnl)goto 21681                                          
      if(flmin.ge.1.0)goto 21681                                        
      me=0                                                              
      do 22261 j=1,nin                                                  
      if(a(j,1,ilm).ne.0.0) me=me+1                                     
22261 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 21682                                            
      if(dev(ilm).gt.devmax)goto 21682                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 21682                          
21681 continue                                                          
21682 continue                                                          
      g=log(q)                                                          
      do 22271 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
22271 continue                                                          
      continue                                                          
      deallocate(sxp,b,bs,r,q,mm,is,ga,ixx,gk,del,sxpl)                 
      return                                                            
      end                                                               
      subroutine multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,
     *nx,nlam,  flmin,ulam,shri,intr,maxit,xv,xb,xs,lmu,a0,a,m,kin,dev0,
     *dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),y(no,nc),g(no,nc),w(no),vp(ni)              
      double precision ulam(nlam),xb(ni),xs(ni),xv(ni)                  
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(
     *2,ni)
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                        
      double precision, dimension (:,:), allocatable :: q,r,b,bs        
      double precision, dimension (:), allocatable :: sxp,sxpl,ga,gk    
      double precision, dimension (:), allocatable :: del,sc,svr        
      integer, dimension (:), allocatable :: mm,is,iy,isc               
      allocate(b(0:ni,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      allocate(bs(0:ni,1:nc),stat=jerr)                                 
      if(jerr.ne.0) return                                              
      allocate(q(1:no,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      allocate(r(1:no,1:nc),stat=jerr)                                  
      if(jerr.ne.0) return                                              
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx,itrace)     
      exmn=-exmx                                                        
      allocate(mm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ga(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(gk(1:nc),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(del(1:nc),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(iy(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(is(1:max(nc,ni)),stat=jerr)                              
      if(jerr.ne.0) return                                              
      allocate(sxp(1:no),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(sxpl(1:no),stat=jerr)                                    
      if(jerr.ne.0) return                                              
      allocate(svr(1:nc),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      allocate(sc(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(isc(1:nc),stat=jerr)                                     
      if(jerr.ne.0) return                                              
      pmax=1.0-pmin                                                     
      emin=pmin/pmax                                                    
      emax=1.0/emin                                                     
      bta=parm                                                          
      omb=1.0-bta                                                       
      dev1=0.0                                                          
      dev0=0.0                                                          
      do 22281 ic=1,nc                                                  
      q0=dot_product(w,y(:,ic))                                         
      if(q0 .gt. pmin)goto 22301                                        
      jerr =8000+ic                                                     
      return                                                            
22301 continue                                                          
      if(q0 .lt. pmax)goto 22321                                        
      jerr =9000+ic                                                     
      return                                                            
22321 continue                                                          
      b(1:ni,ic)=0.0                                                    
      if(intr .ne. 0)goto 22341                                         
      q0=1.0/nc                                                         
      b(0,ic)=0.0                                                       
      goto 22351                                                        
22341 continue                                                          
      b(0,ic)=log(q0)                                                   
      dev1=dev1-q0*b(0,ic)                                              
22351 continue                                                          
      continue                                                          
22281 continue                                                          
      continue                                                          
      if(intr.eq.0) dev1=log(float(nc))                                 
      iy=0                                                              
      al=0.0                                                            
      if(nonzero(no*nc,g) .ne. 0)goto 22371                             
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                      
      sxp=0.0                                                           
      do 22381 ic=1,nc                                                  
      q(:,ic)=exp(b(0,ic))                                              
      sxp=sxp+q(:,ic)                                                   
22381 continue                                                          
      continue                                                          
      goto 22391                                                        
22371 continue                                                          
      do 22401 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
22401 continue                                                          
      continue                                                          
      sxp=0.0                                                           
      if(intr .ne. 0)goto 22421                                         
      b(0,:)=0.0                                                        
      goto 22431                                                        
22421 continue                                                          
      call kazero(nc,no,y,g,w,b(0,:),jerr)                              
      if(jerr.ne.0) return                                              
22431 continue                                                          
      continue                                                          
      dev1=0.0                                                          
      do 22441 ic=1,nc                                                  
      q(:,ic)=b(0,ic)+g(:,ic)                                           
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                          
      q(:,ic)=exp(q(:,ic))                                              
      sxp=sxp+q(:,ic)                                                   
22441 continue                                                          
      continue                                                          
      sxpl=w*log(sxp)                                                   
      do 22451 ic=1,nc                                                  
      dev1=dev1+dot_product(y(:,ic),sxpl)                               
22451 continue                                                          
      continue                                                          
22391 continue                                                          
      continue                                                          
      do 22461 ic=1,nc                                                  
      do 22471 i=1,no                                                   
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))            
22471 continue                                                          
      continue                                                          
22461 continue                                                          
      continue                                                          
      dev0=dev0+dev1                                                    
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 22491                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
22491 continue                                                          
      m=0                                                               
      mm=0                                                              
      nin=0                                                             
      nlp=0                                                             
      mnl=min(mnlam,nlam)                                               
      bs=0.0                                                            
      shr=shri*dev0                                                     
      ga=0.0                                                            
      do 22501 ic=1,nc                                                  
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                   
      svr(ic)=sum(r(:,ic))                                              
      do 22511 j=1,ni                                                   
      if(ju(j).eq.0)goto 22511                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                          
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                         
22511 continue                                                          
      continue                                                          
22501 continue                                                          
      continue                                                          
      ga=sqrt(ga)                                                       
      do 22521 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 22541                                      
      al=ulam(ilm)                                                      
      goto 22531                                                        
22541 if(ilm .le. 2)goto 22551                                          
      al=al*alf                                                         
      goto 22531                                                        
22551 if(ilm .ne. 1)goto 22561                                          
      al=big                                                            
      goto 22571                                                        
22561 continue                                                          
      al0=0.0                                                           
      do 22581 j=1,ni                                                   
      if(ju(j).eq.0)goto 22581                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
22581 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
22571 continue                                                          
22531 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 22591 k=1,ni                                                   
      if(iy(k).eq.1)goto 22591                                          
      if(ju(k).eq.0)goto 22591                                          
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                   
22591 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
22601 continue                                                          
      ixx=0                                                             
      jxx=ixx                                                           
      kxx=jxx                                                           
      t=0.0                                                             
      do 22611 ic=1,nc                                                  
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                    
22611 continue                                                          
      continue                                                          
      if(t .ge. eps)goto 22631                                          
      kxx=1                                                             
      goto 22602                                                        
22631 continue                                                          
      t=2.0*t                                                           
      alt=al1/t                                                         
      al2t=al2/t                                                        
      do 22641 ic=1,nc                                                  
      bs(0,ic)=b(0,ic)                                                  
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                       
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                 
      svr(ic)=sum(r(:,ic))                                              
      if(intr .eq. 0)goto 22661                                         
      b(0,ic)=b(0,ic)+svr(ic)                                           
      r(:,ic)=r(:,ic)-svr(ic)*w                                         
      dlx=max(dlx,svr(ic)**2)                                           
22661 continue                                                          
22641 continue                                                          
      continue                                                          
      continue                                                          
22671 continue                                                          
      nlp=nlp+nc                                                        
      dlx=0.0                                                           
      do 22681 k=1,ni                                                   
      if(iy(k).eq.0)goto 22681                                          
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      del=b(k,:)                                                        
      gkn=0.0                                                           
      do 22691 ic=1,nc                                                  
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))-svr(ic)*xb(k))/xs(k)     
      gk(ic)=u+del(ic)*xv(k)                                            
      gkn=gkn+gk(ic)**2                                                 
22691 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-alt*vp(k)/gkn                                               
      if(u .gt. 0.0)goto 22711                                          
      b(k,:)=0.0                                                        
      goto 22721                                                        
22711 continue                                                          
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                  
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                              
22721 continue                                                          
      continue                                                          
      del=b(k,:)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 22681                             
      do 22731 ic=1,nc                                                  
      dlx=max(dlx,xv(k)*del(ic)**2)                                     
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x
     *b(k))/xs(k)
22731 continue                                                          
      continue                                                          
      if(mm(k) .ne. 0)goto 22751                                        
      nin=nin+1                                                         
      if(nin .le. nx)goto 22771                                         
      jxx=1                                                             
      goto 22682                                                        
22771 continue                                                          
      mm(k)=nin                                                         
      m(nin)=k                                                          
22751 continue                                                          
22681 continue                                                          
22682 continue                                                          
      if(jxx.gt.0)goto 22672                                            
      if(dlx.lt.shr)goto 22672                                          
      if(nlp .le. maxit)goto 22791                                      
      jerr=-ilm                                                         
      return                                                            
22791 continue                                                          
      continue                                                          
22801 continue                                                          
      nlp=nlp+nc                                                        
      dlx=0.0                                                           
      do 22811 l=1,nin                                                  
      k=m(l)                                                            
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      del=b(k,:)                                                        
      gkn=0.0                                                           
      do 22821 ic=1,nc                                                  
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))  -svr(ic)*xb(k))/xs(k)   
      gk(ic)=u+del(ic)*xv(k)                                            
      gkn=gkn+gk(ic)**2                                                 
22821 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-alt*vp(k)/gkn                                               
      if(u .gt. 0.0)goto 22841                                          
      b(k,:)=0.0                                                        
      goto 22851                                                        
22841 continue                                                          
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                  
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                              
22851 continue                                                          
      continue                                                          
      del=b(k,:)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 22811                             
      do 22861 ic=1,nc                                                  
      dlx=max(dlx,xv(k)*del(ic)**2)                                     
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x
     *b(k))/xs(k)
22861 continue                                                          
      continue                                                          
22811 continue                                                          
      continue                                                          
      if(dlx.lt.shr)goto 22802                                          
      if(nlp .le. maxit)goto 22881                                      
      jerr=-ilm                                                         
      return                                                            
22881 continue                                                          
      goto 22801                                                        
22802 continue                                                          
      goto 22671                                                        
22672 continue                                                          
      if(jxx.gt.0)goto 22602                                            
      do 22891 ic=1,nc                                                  
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                            
      if(ixx .ne. 0)goto 22911                                          
      do 22921 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 22941                
      ixx=1                                                             
      goto 22922                                                        
22941 continue                                                          
22921 continue                                                          
22922 continue                                                          
22911 continue                                                          
      sc=b(0,ic)+g(:,ic)                                                
      b0=0.0                                                            
      do 22951 j=1,nin                                                  
      l=m(j)                                                            
      jb=ix(l)                                                          
      je=ix(l+1)-1                                                      
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                
      b0=b0-b(l,ic)*xb(l)/xs(l)                                         
22951 continue                                                          
      continue                                                          
      sc=min(max(exmn,sc+b0),exmx)                                      
      sxp=sxp-q(:,ic)                                                   
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                       
      sxp=sxp+q(:,ic)                                                   
22891 continue                                                          
      continue                                                          
      s=sum(b(0,:))/nc                                                  
      b(0,:)=b(0,:)-s                                                   
      if(jxx.gt.0)goto 22602                                            
      if(ixx .ne. 0)goto 22971                                          
      do 22981 j=1,ni                                                   
      if(iy(j).eq.1)goto 22981                                          
      if(ju(j).eq.0)goto 22981                                          
      ga(j)=0.0                                                         
22981 continue                                                          
      continue                                                          
      do 22991 ic=1,nc                                                  
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                   
      do 23001 j=1,ni                                                   
      if(iy(j).eq.1)goto 23001                                          
      if(ju(j).eq.0)goto 23001                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                          
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                         
23001 continue                                                          
      continue                                                          
22991 continue                                                          
      continue                                                          
      ga=sqrt(ga)                                                       
      do 23011 k=1,ni                                                   
      if(iy(k).eq.1)goto 23011                                          
      if(ju(k).eq.0)goto 23011                                          
      if(ga(k) .le. al1*vp(k))goto 23031                                
      iy(k)=1                                                           
      ixx=1                                                             
23031 continue                                                          
23011 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10880                                          
      goto 22602                                                        
22971 continue                                                          
      goto 22601                                                        
22602 continue                                                          
      if(kxx .le. 0)goto 23051                                          
      jerr=-20000-ilm                                                   
      goto 22522                                                        
23051 continue                                                          
      if(jxx .le. 0)goto 23071                                          
      jerr=-10000-ilm                                                   
      goto 22522                                                        
23071 continue                                                          
      devi=0.0                                                          
      do 23081 ic=1,nc                                                  
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                       
      a0(ic,ilm)=b(0,ic)                                                
      do 23091 i=1,no                                                   
      if(y(i,ic).le.0.0)goto 23091                                      
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                        
23091 continue                                                          
      continue                                                          
23081 continue                                                          
      continue                                                          
      kin(ilm)=nin                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dev1-devi)/dev0                                         
      if(ilm.lt.mnl)goto 22521                                          
      if(flmin.ge.1.0)goto 22521                                        
      me=0                                                              
      do 23101 j=1,nin                                                  
      if(a(j,1,ilm).ne.0.0) me=me+1                                     
23101 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 22522                                            
      if(dev(ilm).gt.devmax)goto 22522                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 22522                          
22521 continue                                                          
22522 continue                                                          
      g=log(q)                                                          
      do 23111 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
23111 continue                                                          
      continue                                                          
      deallocate(sxp,b,bs,r,q,mm,is,sc,ga,iy,gk,del,sxpl)               
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
