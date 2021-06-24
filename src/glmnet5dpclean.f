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
10880 if(nlp .le. maxit)goto 10901                                      
      jerr=-m                                                           
      return                                                            
10901 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10911 k=1,ni                                                   
      if(ix(k).eq.0)goto 10911                                          
      gk=dot_product(y,x(:,k))                                          
      ak=a(k)                                                           
      u=gk+ak*xv(k)                                                     
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 10911                                          
      if(mm(k) .ne. 0)goto 10931                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 10912                                           
      mm(k)=nin                                                         
      ia(nin)=k                                                         
10931 continue                                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*gk-del*xv(k))                                    
      y=y-del*x(:,k)                                                    
      dlx=max(xv(k)*del**2,dlx)                                         
10911 continue                                                          
10912 continue                                                          
      if(nin.gt.nx)goto 10872                                           
      if(dlx .ge. thr)goto 10951                                        
      ixx=0                                                             
      do 10961 k=1,ni                                                   
      if(ix(k).eq.1)goto 10961                                          
      if(ju(k).eq.0)goto 10961                                          
      g(k)=abs(dot_product(y,x(:,k)))                                   
      if(g(k) .le. ab*vp(k))goto 10981                                  
      ix(k)=1                                                           
      ixx=1                                                             
10981 continue                                                          
10961 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10880                                          
      goto 10872                                                        
10951 continue                                                          
      if(nlp .le. maxit)goto 11001                                      
      jerr=-m                                                           
      return                                                            
11001 continue                                                          
10360 continue                                                          
      iz=1                                                              
      continue                                                          
11011 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 11021 l=1,nin                                                  
      k=ia(l)                                                           
      gk=dot_product(y,x(:,k))                                          
      ak=a(k)                                                           
      u=gk+ak*xv(k)                                                     
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 11021                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*gk-del*xv(k))                                    
      y=y-del*x(:,k)                                                    
      dlx=max(xv(k)*del**2,dlx)                                         
11021 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 11012                                          
      if(nlp .le. maxit)goto 11041                                      
      jerr=-m                                                           
      return                                                            
11041 continue                                                          
      goto 11011                                                        
11012 continue                                                          
      jz=0                                                              
      goto 10871                                                        
10872 continue                                                          
      if(nin .le. nx)goto 11061                                         
      jerr=-10000-m                                                     
      goto 10792                                                        
11061 continue                                                          
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                             
      kin(m)=nin                                                        
      rsqo(m)=rsq                                                       
      almo(m)=alm                                                       
      lmu=m                                                             
      if(m.lt.mnl)goto 10791                                            
      if(flmin.ge.1.0)goto 10791                                        
      me=0                                                              
      do 11071 j=1,nin                                                  
      if(ao(j,m).ne.0.0) me=me+1                                        
11071 continue                                                          
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
      do 11081 j=1,ni                                                   
      ju(j)=0                                                           
      t=x(1,j)                                                          
      do 11091 i=2,no                                                   
      if(x(i,j).eq.t)goto 11091                                         
      ju(j)=1                                                           
      goto 11092                                                        
11091 continue                                                          
11092 continue                                                          
11081 continue                                                          
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
      do 11101 i=1,n                                                    
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                   
11101 continue                                                          
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
      if(maxval(vp) .gt. 0.0)goto 11121                                 
      jerr=10000                                                        
      return                                                            
11121 continue                                                          
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      if(ka .ne. 1)goto 11141                                           
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,u
     *lam,thr,isd,  intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 11151                                                        
11141 continue                                                          
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,ul
     *am,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
11151 continue                                                          
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
      if(maxval(ju) .gt. 0)goto 11171                                   
      jerr=7777                                                         
      return                                                            
11171 continue                                                          
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys,xv,jer
     *r)
      if(jerr.ne.0) return                                              
      cl=cl/ys                                                          
      if(isd .le. 0)goto 11191                                          
      do 11201 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
11201 continue                                                          
      continue                                                          
11191 continue                                                          
      if(flmin.ge.1.0) vlam=ulam/ys                                     
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                              
      do 11211 k=1,lmu                                                  
      alm(k)=ys*alm(k)                                                  
      nk=nin(k)                                                         
      do 11221 l=1,nk                                                   
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                      
11221 continue                                                          
      continue                                                          
      a0(k)=0.0                                                         
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))       
11211 continue                                                          
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
      if(intr .ne. 0)goto 11241                                         
      ym=0.0                                                            
      ys=sqrt(dot_product(w,y**2))                                      
      y=y/ys                                                            
      do 11251 j=1,ni                                                   
      if(ju(j).eq.0)goto 11251                                          
      xm(j)=0.0                                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                       
      if(isd .eq. 0)goto 11271                                          
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                         
      vc=xv(j)-xbq                                                      
      xs(j)=sqrt(vc)                                                    
      xv(j)=1.0+xbq/vc                                                  
      goto 11281                                                        
11271 continue                                                          
      xs(j)=1.0                                                         
11281 continue                                                          
      continue                                                          
11251 continue                                                          
      continue                                                          
      goto 11291                                                        
11241 continue                                                          
      do 11301 j=1,ni                                                   
      if(ju(j).eq.0)goto 11301                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2              
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                    
11301 continue                                                          
      continue                                                          
      if(isd .ne. 0)goto 11321                                          
      xs=1.0                                                            
      goto 11331                                                        
11321 continue                                                          
      xv=1.0                                                            
11331 continue                                                          
      continue                                                          
      ym=dot_product(w,y)                                               
      y=y-ym                                                            
      ys=sqrt(dot_product(w,y**2))                                      
      y=y/ys                                                            
11291 continue                                                          
      continue                                                          
      g=0.0                                                             
      do 11341 j=1,ni                                                   
      if(ju(j).eq.0)goto 11341                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)        
11341 continue                                                          
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
      if(flmin .ge. 1.0)goto 11361                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
11361 continue                                                          
      rsq=0.0                                                           
      a=0.0                                                             
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      iz=0                                                              
      mnl=min(mnlam,nlam)                                               
      do 11371 m=1,nlam                                                 
      if(itrace.ne.0) call setpb(m-1)                                   
      if(flmin .lt. 1.0)goto 11391                                      
      alm=ulam(m)                                                       
      goto 11381                                                        
11391 if(m .le. 2)goto 11401                                            
      alm=alm*alf                                                       
      goto 11381                                                        
11401 if(m .ne. 1)goto 11411                                            
      alm=big                                                           
      goto 11421                                                        
11411 continue                                                          
      alm=0.0                                                           
      do 11431 j=1,ni                                                   
      if(ju(j).eq.0)goto 11431                                          
      if(vp(j).le.0.0)goto 11431                                        
      alm=max(alm,abs(g(j))/vp(j))                                      
11431 continue                                                          
      continue                                                          
      alm=alf*alm/max(bta,1.0d-3)                                       
11421 continue                                                          
11381 continue                                                          
      dem=alm*omb                                                       
      ab=alm*bta                                                        
      rsq0=rsq                                                          
      jz=1                                                              
      continue                                                          
11441 continue                                                          
      if(iz*jz.ne.0) go to 10360                                        
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 11451 k=1,ni                                                   
      if(ju(k).eq.0)goto 11451                                          
      ak=a(k)                                                           
      u=g(k)+ak*xv(k)                                                   
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 11451                                          
      if(mm(k) .ne. 0)goto 11471                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 11452                                           
      do 11481 j=1,ni                                                   
      if(ju(j).eq.0)goto 11481                                          
      if(mm(j) .eq. 0)goto 11501                                        
      c(j,nin)=c(k,mm(j))                                               
      goto 11481                                                        
11501 continue                                                          
      if(j .ne. k)goto 11521                                            
      c(j,nin)=xv(j)                                                    
      goto 11481                                                        
11521 continue                                                          
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))    
11481 continue                                                          
      continue                                                          
      mm(k)=nin                                                         
      ia(nin)=k                                                         
11471 continue                                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                  
      dlx=max(xv(k)*del**2,dlx)                                         
      do 11531 j=1,ni                                                   
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                           
11531 continue                                                          
      continue                                                          
11451 continue                                                          
11452 continue                                                          
      if(dlx.lt.thr)goto 11442                                          
      if(nin.gt.nx)goto 11442                                           
      if(nlp .le. maxit)goto 11551                                      
      jerr=-m                                                           
      return                                                            
11551 continue                                                          
10360 continue                                                          
      iz=1                                                              
      da(1:nin)=a(ia(1:nin))                                            
      continue                                                          
11561 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 11571 l=1,nin                                                  
      k=ia(l)                                                           
      ak=a(k)                                                           
      u=g(k)+ak*xv(k)                                                   
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 11571                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                  
      dlx=max(xv(k)*del**2,dlx)                                         
      do 11581 j=1,nin                                                  
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                              
11581 continue                                                          
      continue                                                          
11571 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 11562                                          
      if(nlp .le. maxit)goto 11601                                      
      jerr=-m                                                           
      return                                                            
11601 continue                                                          
      goto 11561                                                        
11562 continue                                                          
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                  
      do 11611 j=1,ni                                                   
      if(mm(j).ne.0)goto 11611                                          
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))        
11611 continue                                                          
      continue                                                          
      jz=0                                                              
      goto 11441                                                        
11442 continue                                                          
      if(nin .le. nx)goto 11631                                         
      jerr=-10000-m                                                     
      goto 11372                                                        
11631 continue                                                          
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                             
      kin(m)=nin                                                        
      rsqo(m)=rsq                                                       
      almo(m)=alm                                                       
      lmu=m                                                             
      if(m.lt.mnl)goto 11371                                            
      if(flmin.ge.1.0)goto 11371                                        
      me=0                                                              
      do 11641 j=1,nin                                                  
      if(ao(j,m).ne.0.0) me=me+1                                        
11641 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 11372                                            
      if(rsq-rsq0.lt.sml*rsq)goto 11372                                 
      if(rsq.gt.rsqmax)goto 11372                                       
11371 continue                                                          
11372 continue                                                          
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
      if(maxval(ju) .gt. 0)goto 11661                                   
      jerr=7777                                                         
      return                                                            
11661 continue                                                          
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,xv,jerr
     *)
      if(jerr.ne.0) return                                              
      cl=cl/ys                                                          
      if(isd .le. 0)goto 11681                                          
      do 11691 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
11691 continue                                                          
      continue                                                          
11681 continue                                                          
      if(flmin.ge.1.0) vlam=ulam/ys                                     
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                              
      do 11701 k=1,lmu                                                  
      alm(k)=ys*alm(k)                                                  
      nk=nin(k)                                                         
      do 11711 l=1,nk                                                   
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                      
11711 continue                                                          
      continue                                                          
      a0(k)=0.0                                                         
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))       
11701 continue                                                          
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
      if(intr .ne. 0)goto 11731                                         
      ym=0.0                                                            
      ys=sqrt(dot_product(w,y**2))                                      
      y=y/ys                                                            
      do 11741 j=1,ni                                                   
      if(ju(j).eq.0)goto 11741                                          
      xm(j)=0.0                                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                       
      if(isd .eq. 0)goto 11761                                          
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                         
      vc=xv(j)-xbq                                                      
      xs(j)=sqrt(vc)                                                    
      xv(j)=1.0+xbq/vc                                                  
      goto 11771                                                        
11761 continue                                                          
      xs(j)=1.0                                                         
11771 continue                                                          
      continue                                                          
11741 continue                                                          
      continue                                                          
      return                                                            
11731 continue                                                          
      do 11781 j=1,ni                                                   
      if(ju(j).eq.0)goto 11781                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2              
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                    
11781 continue                                                          
      continue                                                          
      if(isd .ne. 0)goto 11801                                          
      xs=1.0                                                            
      goto 11811                                                        
11801 continue                                                          
      xv=1.0                                                            
11811 continue                                                          
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
      if(flmin .ge. 1.0)goto 11831                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
11831 continue                                                          
      rsq=0.0                                                           
      a=0.0                                                             
      mm=0                                                              
      o=0.0                                                             
      nlp=0                                                             
      nin=nlp                                                           
      iz=0                                                              
      mnl=min(mnlam,nlam)                                               
      do 11841 j=1,ni                                                   
      if(ju(j).eq.0)goto 11841                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j)) 
11841 continue                                                          
      continue                                                          
      do 11851 m=1,nlam                                                 
      if(itrace.ne.0) call setpb(m-1)                                   
      alm0=alm                                                          
      if(flmin .lt. 1.0)goto 11871                                      
      alm=ulam(m)                                                       
      goto 11861                                                        
11871 if(m .le. 2)goto 11881                                            
      alm=alm*alf                                                       
      goto 11861                                                        
11881 if(m .ne. 1)goto 11891                                            
      alm=big                                                           
      goto 11901                                                        
11891 continue                                                          
      alm0=0.0                                                          
      do 11911 j=1,ni                                                   
      if(ju(j).eq.0)goto 11911                                          
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                        
11911 continue                                                          
      continue                                                          
      alm0=alm0/max(bta,1.0d-3)                                         
      alm=alf*alm0                                                      
11901 continue                                                          
11861 continue                                                          
      dem=alm*omb                                                       
      ab=alm*bta                                                        
      rsq0=rsq                                                          
      jz=1                                                              
      tlam=bta*(2.0*alm-alm0)                                           
      do 11921 k=1,ni                                                   
      if(iy(k).eq.1)goto 11921                                          
      if(ju(k).eq.0)goto 11921                                          
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                    
11921 continue                                                          
      continue                                                          
      continue                                                          
11931 continue                                                          
      if(iz*jz.ne.0) go to 10360                                        
10880 if(nlp .le. maxit)goto 11951                                      
      jerr=-m                                                           
      return                                                            
11951 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 11961 k=1,ni                                                   
      if(iy(k).eq.0)goto 11961                                          
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)        
      ak=a(k)                                                           
      u=gk+ak*xv(k)                                                     
      v=abs(u)-vp(k)*ab                                                 
      a(k)=0.0                                                          
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d
     *em)))
      if(a(k).eq.ak)goto 11961                                          
      if(mm(k) .ne. 0)goto 11981                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 11962                                           
      mm(k)=nin                                                         
      ia(nin)=k                                                         
11981 continue                                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*gk-del*xv(k))                                    
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                      
      o=o+del*xm(k)/xs(k)                                               
      dlx=max(xv(k)*del**2,dlx)                                         
11961 continue                                                          
11962 continue                                                          
      if(nin.gt.nx)goto 11932                                           
      if(dlx .ge. thr)goto 12001                                        
      ixx=0                                                             
      do 12011 j=1,ni                                                   
      if(iy(j).eq.1)goto 12011                                          
      if(ju(j).eq.0)goto 12011                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j)) 
      if(g(j) .le. ab*vp(j))goto 12031                                  
      iy(j)=1                                                           
      ixx=1                                                             
12031 continue                                                          
12011 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10880                                          
      goto 11932                                                        
12001 continue                                                          
      if(nlp .le. maxit)goto 12051                                      
      jerr=-m                                                           
      return                                                            
12051 continue                                                          
10360 continue                                                          
      iz=1                                                              
      continue                                                          
12061 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 12071 l=1,nin                                                  
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
      if(a(k).eq.ak)goto 12071                                          
      del=a(k)-ak                                                       
      rsq=rsq+del*(2.0*gk-del*xv(k))                                    
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                      
      o=o+del*xm(k)/xs(k)                                               
      dlx=max(xv(k)*del**2,dlx)                                         
12071 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 12062                                          
      if(nlp .le. maxit)goto 12091                                      
      jerr=-m                                                           
      return                                                            
12091 continue                                                          
      goto 12061                                                        
12062 continue                                                          
      jz=0                                                              
      goto 11931                                                        
11932 continue                                                          
      if(nin .le. nx)goto 12111                                         
      jerr=-10000-m                                                     
      goto 11852                                                        
12111 continue                                                          
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                             
      kin(m)=nin                                                        
      rsqo(m)=rsq                                                       
      almo(m)=alm                                                       
      lmu=m                                                             
      if(m.lt.mnl)goto 11851                                            
      if(flmin.ge.1.0)goto 11851                                        
      me=0                                                              
      do 12121 j=1,nin                                                  
      if(ao(j,m).ne.0.0) me=me+1                                        
12121 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 11852                                            
      if(rsq-rsq0.lt.sml*rsq)goto 11852                                 
      if(rsq.gt.rsqmax)goto 11852                                       
11851 continue                                                          
11852 continue                                                          
      deallocate(a,mm,g,iy)                                             
      return                                                            
      end                                                               
      subroutine spchkvars(no,ni,x,ix,ju)                               
      implicit double precision(a-h,o-z)                                
      double precision x(*)                                             
      integer ix(*),ju(ni)                                              
      do 12131 j=1,ni                                                   
      ju(j)=0                                                           
      jb=ix(j)                                                          
      nj=ix(j+1)-jb                                                     
      if(nj.eq.0)goto 12131                                             
      je=ix(j+1)-1                                                      
      if(nj .ge. no)goto 12151                                          
      do 12161 i=jb,je                                                  
      if(x(i).eq.0.0)goto 12161                                         
      ju(j)=1                                                           
      goto 12162                                                        
12161 continue                                                          
12162 continue                                                          
      goto 12171                                                        
12151 continue                                                          
      t=x(jb)                                                           
      do 12181 i=jb+1,je                                                
      if(x(i).eq.t)goto 12181                                           
      ju(j)=1                                                           
      goto 12182                                                        
12181 continue                                                          
12182 continue                                                          
12171 continue                                                          
      continue                                                          
12131 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                      
      implicit double precision(a-h,o-z)                                
      double precision ca(*),x(*),f(n)                                  
      integer ia(*),ix(*),jx(*)                                         
      f=a0                                                              
      do 12191 j=1,nin                                                  
      k=ia(j)                                                           
      kb=ix(k)                                                          
      ke=ix(k+1)-1                                                      
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                          
12191 continue                                                          
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
12201 continue                                                          
      continue                                                          
12211 if(mx(i).ge.my(j))goto 12212                                      
      i=i+1                                                             
      if(i.gt.nx) go to 12220                                           
      goto 12211                                                        
12212 continue                                                          
      if(mx(i).eq.my(j)) go to 12230                                    
      continue                                                          
12241 if(my(j).ge.mx(i))goto 12242                                      
      j=j+1                                                             
      if(j.gt.ny) go to 12220                                           
      goto 12241                                                        
12242 continue                                                          
      if(mx(i).eq.my(j)) go to 12230                                    
      goto 12201                                                        
12230 continue                                                          
      s=s+w(mx(i))*x(i)*y(j)                                            
      i=i+1                                                             
      if(i.gt.nx)goto 12202                                             
      j=j+1                                                             
      if(j.gt.ny)goto 12202                                             
      goto 12201                                                        
12202 continue                                                          
12220 continue                                                          
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
      if(maxval(vp) .gt. 0.0)goto 12261                                 
      jerr=10000                                                        
      return                                                            
12261 continue                                                          
      allocate(ww(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      if(kopt .ne. 2)goto 12281                                         
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
12281 continue                                                          
      if(isd .le. 0)goto 12301                                          
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
12301 continue                                                          
      call chkvars(no,ni,x,ju)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 12321                                   
      jerr=7777                                                         
      return                                                            
12321 continue                                                          
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      do 12331 i=1,no                                                   
      ww(i)=sum(y(i,:))                                                 
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                              
12331 continue                                                          
      continue                                                          
      sw=sum(ww)                                                        
      ww=ww/sw                                                          
      if(nc .ne. 1)goto 12351                                           
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                     
      if(isd .le. 0)goto 12371                                          
      do 12381 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
12381 continue                                                          
      continue                                                          
12371 continue                                                          
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,cl,ne,nx,nlam,fl
     *min,ulam,  thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,n
     *lp,jerr)
      goto 12341                                                        
12351 if(kopt .ne. 2)goto 12391                                         
      call multlstandard1(no,ni,x,ww,ju,isd,intr,xm,xs,xv)              
      if(isd .le. 0)goto 12411                                          
      do 12421 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
12421 continue                                                          
      continue                                                          
12411 continue                                                          
      call multlognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,
     *ulam,thr,  intr,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 12431                                                        
12391 continue                                                          
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                     
      if(isd .le. 0)goto 12451                                          
      do 12461 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
12461 continue                                                          
      continue                                                          
12451 continue                                                          
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam
     *,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
12431 continue                                                          
12341 continue                                                          
      if(jerr.gt.0) return                                              
      dev0=2.0*sw*dev0                                                  
      do 12471 k=1,lmu                                                  
      nk=nin(k)                                                         
      do 12481 ic=1,nc                                                  
      if(isd .le. 0)goto 12501                                          
      do 12511 l=1,nk                                                   
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                   
12511 continue                                                          
      continue                                                          
12501 continue                                                          
      if(intr .ne. 0)goto 12531                                         
      a0(ic,k)=0.0                                                      
      goto 12541                                                        
12531 continue                                                          
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))         
12541 continue                                                          
      continue                                                          
12481 continue                                                          
      continue                                                          
12471 continue                                                          
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
      if(intr .ne. 0)goto 12561                                         
      do 12571 j=1,ni                                                   
      if(ju(j).eq.0)goto 12571                                          
      xm(j)=0.0                                                         
      if(isd .eq. 0)goto 12591                                          
      vc=dot_product(w,x(:,j)**2)-dot_product(w,x(:,j))**2              
      xs(j)=sqrt(vc)                                                    
      x(:,j)=x(:,j)/xs(j)                                               
12591 continue                                                          
12571 continue                                                          
      continue                                                          
      return                                                            
12561 continue                                                          
      do 12601 j=1,ni                                                   
      if(ju(j).eq.0)goto 12601                                          
      xm(j)=dot_product(w,x(:,j))                                       
      x(:,j)=x(:,j)-xm(j)                                               
      if(isd .le. 0)goto 12621                                          
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                              
      x(:,j)=x(:,j)/xs(j)                                               
12621 continue                                                          
12601 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine multlstandard1(no,ni,x,w,ju,isd,intr,xm,xs,xv)         
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),w(no),xm(ni),xs(ni),xv(ni)              
      integer ju(ni)                                                    
      if(intr .ne. 0)goto 12641                                         
      do 12651 j=1,ni                                                   
      if(ju(j).eq.0)goto 12651                                          
      xm(j)=0.0                                                         
      xv(j)=dot_product(w,x(:,j)**2)                                    
      if(isd .eq. 0)goto 12671                                          
      xbq=dot_product(w,x(:,j))**2                                      
      vc=xv(j)-xbq                                                      
      xs(j)=sqrt(vc)                                                    
      x(:,j)=x(:,j)/xs(j)                                               
      xv(j)=1.0+xbq/vc                                                  
12671 continue                                                          
12651 continue                                                          
      continue                                                          
      return                                                            
12641 continue                                                          
      do 12681 j=1,ni                                                   
      if(ju(j).eq.0)goto 12681                                          
      xm(j)=dot_product(w,x(:,j))                                       
      x(:,j)=x(:,j)-xm(j)                                               
      xv(j)=dot_product(w,x(:,j)**2)                                    
      if(isd .le. 0)goto 12701                                          
      xs(j)=sqrt(xv(j))                                                 
      x(:,j)=x(:,j)/xs(j)                                               
      xv(j)=1.0                                                         
12701 continue                                                          
12681 continue                                                          
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
      if(q0 .gt. pmin)goto 12721                                        
      jerr=8001                                                         
      return                                                            
12721 continue                                                          
      if(q0 .lt. 1.0-pmin)goto 12741                                    
      jerr=9001                                                         
      return                                                            
12741 continue                                                          
      if(intr.eq.0.0) q0=0.5                                            
      ixx=0                                                             
      al=0.0                                                            
      bz=0.0                                                            
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                 
      if(nonzero(no,g) .ne. 0)goto 12761                                
      vi=q0*(1.0-q0)                                                    
      b(0)=bz                                                           
      v=vi*w                                                            
      r=w*(y-q0)                                                        
      q=q0                                                              
      xmz=vi                                                            
      dev1=-(bz*q0+log(1.0-q0))                                         
      goto 12771                                                        
12761 continue                                                          
      b(0)=0.0                                                          
      if(intr .eq. 0)goto 12791                                         
      b(0)=azero(no,y,g,w,jerr)                                         
      if(jerr.ne.0) return                                              
12791 continue                                                          
      q=1.0/(1.0+exp(-b(0)-g))                                          
      v=w*q*(1.0-q)                                                     
      r=w*(y-q)                                                         
      xmz=sum(v)                                                        
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                     
12771 continue                                                          
      continue                                                          
      if(kopt .le. 0)goto 12811                                         
      if(isd .le. 0 .or. intr .eq. 0)goto 12831                         
      xv=0.25                                                           
      goto 12841                                                        
12831 continue                                                          
      do 12851 j=1,ni                                                   
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                
12851 continue                                                          
      continue                                                          
12841 continue                                                          
      continue                                                          
12811 continue                                                          
      dev0=dev1                                                         
      do 12861 i=1,no                                                   
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                     
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))           
12861 continue                                                          
      continue                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 12881                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
12881 continue                                                          
      m=0                                                               
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      mnl=min(mnlam,nlam)                                               
      bs=0.0                                                            
      b(1:ni)=0.0                                                       
      shr=shri*dev0                                                     
      do 12891 j=1,ni                                                   
      if(ju(j).eq.0)goto 12891                                          
      ga(j)=abs(dot_product(r,x(:,j)))                                  
12891 continue                                                          
      continue                                                          
      do 12901 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 12921                                      
      al=ulam(ilm)                                                      
      goto 12911                                                        
12921 if(ilm .le. 2)goto 12931                                          
      al=al*alf                                                         
      goto 12911                                                        
12931 if(ilm .ne. 1)goto 12941                                          
      al=big                                                            
      goto 12951                                                        
12941 continue                                                          
      al0=0.0                                                           
      do 12961 j=1,ni                                                   
      if(ju(j).eq.0)goto 12961                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
12961 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
12951 continue                                                          
12911 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 12971 k=1,ni                                                   
      if(ixx(k).eq.1)goto 12971                                         
      if(ju(k).eq.0)goto 12971                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
12971 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
12981 continue                                                          
      if(nlp .le. maxit)goto 13001                                      
      jerr=-ilm                                                         
      return                                                            
13001 continue                                                          
      bs(0)=b(0)                                                        
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                             
      if(kopt .ne. 0)goto 13021                                         
      do 13031 j=1,ni                                                   
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                    
13031 continue                                                          
      continue                                                          
13021 continue                                                          
      continue                                                          
13041 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 13051 k=1,ni                                                   
      if(ixx(k).eq.0)goto 13051                                         
      bk=b(k)                                                           
      gk=dot_product(r,x(:,k))                                          
      u=gk+xv(k)*b(k)                                                   
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 13071                                         
      b(k)=0.0                                                          
      goto 13081                                                        
13071 continue                                                          
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))       
13081 continue                                                          
      continue                                                          
      d=b(k)-bk                                                         
      if(abs(d).le.0.0)goto 13051                                       
      dlx=max(dlx,xv(k)*d**2)                                           
      r=r-d*v*x(:,k)                                                    
      if(mm(k) .ne. 0)goto 13101                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 13052                                           
      mm(k)=nin                                                         
      m(nin)=k                                                          
13101 continue                                                          
13051 continue                                                          
13052 continue                                                          
      if(nin.gt.nx)goto 13042                                           
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 13121                                          
      b(0)=b(0)+d                                                       
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
13121 continue                                                          
      if(dlx.lt.shr)goto 13042                                          
      if(nlp .le. maxit)goto 13141                                      
      jerr=-ilm                                                         
      return                                                            
13141 continue                                                          
      continue                                                          
13151 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 13161 l=1,nin                                                  
      k=m(l)                                                            
      bk=b(k)                                                           
      gk=dot_product(r,x(:,k))                                          
      u=gk+xv(k)*b(k)                                                   
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 13181                                         
      b(k)=0.0                                                          
      goto 13191                                                        
13181 continue                                                          
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))       
13191 continue                                                          
      continue                                                          
      d=b(k)-bk                                                         
      if(abs(d).le.0.0)goto 13161                                       
      dlx=max(dlx,xv(k)*d**2)                                           
      r=r-d*v*x(:,k)                                                    
13161 continue                                                          
      continue                                                          
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 13211                                          
      b(0)=b(0)+d                                                       
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
13211 continue                                                          
      if(dlx.lt.shr)goto 13152                                          
      if(nlp .le. maxit)goto 13231                                      
      jerr=-ilm                                                         
      return                                                            
13231 continue                                                          
      goto 13151                                                        
13152 continue                                                          
      goto 13041                                                        
13042 continue                                                          
      if(nin.gt.nx)goto 12982                                           
      do 13241 i=1,no                                                   
      fi=b(0)+g(i)                                                      
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))         
      if(fi .ge. fmin)goto 13261                                        
      q(i)=0.0                                                          
      goto 13251                                                        
13261 if(fi .le. fmax)goto 13271                                        
      q(i)=1.0                                                          
      goto 13281                                                        
13271 continue                                                          
      q(i)=1.0/(1.0+exp(-fi))                                           
13281 continue                                                          
13251 continue                                                          
13241 continue                                                          
      continue                                                          
      v=w*q*(1.0-q)                                                     
      xmz=sum(v)                                                        
      if(xmz.le.vmin)goto 12982                                         
      r=w*(y-q)                                                         
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 13301                        
      ix=0                                                              
      do 13311 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 13311                        
      ix=1                                                              
      goto 13312                                                        
13311 continue                                                          
13312 continue                                                          
      if(ix .ne. 0)goto 13331                                           
      do 13341 k=1,ni                                                   
      if(ixx(k).eq.1)goto 13341                                         
      if(ju(k).eq.0)goto 13341                                          
      ga(k)=abs(dot_product(r,x(:,k)))                                  
      if(ga(k) .le. al1*vp(k))goto 13361                                
      ixx(k)=1                                                          
      ix=1                                                              
13361 continue                                                          
13341 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10880                                           
      goto 12982                                                        
13331 continue                                                          
13301 continue                                                          
      goto 12981                                                        
12982 continue                                                          
      if(nin .le. nx)goto 13381                                         
      jerr=-10000-ilm                                                   
      goto 12902                                                        
13381 continue                                                          
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                             
      kin(ilm)=nin                                                      
      a0(ilm)=b(0)                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      devi=dev2(no,w,y,q,pmin)                                          
      dev(ilm)=(dev1-devi)/dev0                                         
      if(xmz.le.vmin)goto 12902                                         
      if(ilm.lt.mnl)goto 12901                                          
      if(flmin.ge.1.0)goto 12901                                        
      me=0                                                              
      do 13391 j=1,nin                                                  
      if(a(j,ilm).ne.0.0) me=me+1                                       
13391 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 12902                                            
      if(dev(ilm).gt.devmax)goto 12902                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12902                          
12901 continue                                                          
12902 continue                                                          
      g=log(q/(1.0-q))                                                  
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                               
      return                                                            
      end                                                               
      function dev2(n,w,y,p,pmin)                                       
      implicit double precision(a-h,o-z)                                
      double precision w(n),y(n),p(n)                                   
      pmax=1.0-pmin                                                     
      s=0.0                                                             
      do 13401 i=1,n                                                    
      pi=min(max(pmin,p(i)),pmax)                                       
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                    
13401 continue                                                          
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
13411 continue                                                          
      w=q*p*(1.0-p)                                                     
      d=(qy-dot_product(q,p))/sum(w)                                    
      az=az+d                                                           
      if(abs(d).lt.eps)goto 13412                                       
      ea0=exp(-az)                                                      
      p=1.0/(1.0+ea0*e)                                                 
      goto 13411                                                        
13412 continue                                                          
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
      do 13421 ic=1,nc                                                  
      q0=dot_product(w,y(:,ic))                                         
      if(q0 .gt. pmin)goto 13441                                        
      jerr =8000+ic                                                     
      return                                                            
13441 continue                                                          
      if(q0 .lt. 1.0-pmin)goto 13461                                    
      jerr =9000+ic                                                     
      return                                                            
13461 continue                                                          
      if(intr .ne. 0)goto 13481                                         
      q0=1.0/nc                                                         
      b(0,ic)=0.0                                                       
      goto 13491                                                        
13481 continue                                                          
      b(0,ic)=log(q0)                                                   
      dev1=dev1-q0*b(0,ic)                                              
13491 continue                                                          
      continue                                                          
      b(1:ni,ic)=0.0                                                    
13421 continue                                                          
      continue                                                          
      if(intr.eq.0) dev1=log(float(nc))                                 
      ixx=0                                                             
      al=0.0                                                            
      if(nonzero(no*nc,g) .ne. 0)goto 13511                             
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                      
      sxp=0.0                                                           
      do 13521 ic=1,nc                                                  
      q(:,ic)=exp(b(0,ic))                                              
      sxp=sxp+q(:,ic)                                                   
13521 continue                                                          
      continue                                                          
      goto 13531                                                        
13511 continue                                                          
      do 13541 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
13541 continue                                                          
      continue                                                          
      sxp=0.0                                                           
      if(intr .ne. 0)goto 13561                                         
      b(0,:)=0.0                                                        
      goto 13571                                                        
13561 continue                                                          
      call kazero(nc,no,y,g,w,b(0,:),jerr)                              
      if(jerr.ne.0) return                                              
13571 continue                                                          
      continue                                                          
      dev1=0.0                                                          
      do 13581 ic=1,nc                                                  
      q(:,ic)=b(0,ic)+g(:,ic)                                           
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                          
      q(:,ic)=exp(q(:,ic))                                              
      sxp=sxp+q(:,ic)                                                   
13581 continue                                                          
      continue                                                          
      sxpl=w*log(sxp)                                                   
      do 13591 ic=1,nc                                                  
      dev1=dev1+dot_product(y(:,ic),sxpl)                               
13591 continue                                                          
      continue                                                          
13531 continue                                                          
      continue                                                          
      do 13601 ic=1,nc                                                  
      do 13611 i=1,no                                                   
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))            
13611 continue                                                          
      continue                                                          
13601 continue                                                          
      continue                                                          
      dev0=dev0+dev1                                                    
      if(kopt .le. 0)goto 13631                                         
      if(isd .le. 0 .or. intr .eq. 0)goto 13651                         
      xv=0.25                                                           
      goto 13661                                                        
13651 continue                                                          
      do 13671 j=1,ni                                                   
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)              
13671 continue                                                          
      continue                                                          
13661 continue                                                          
      continue                                                          
13631 continue                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 13691                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
13691 continue                                                          
      m=0                                                               
      mm=0                                                              
      nin=0                                                             
      nlp=0                                                             
      mnl=min(mnlam,nlam)                                               
      bs=0.0                                                            
      shr=shri*dev0                                                     
      ga=0.0                                                            
      do 13701 ic=1,nc                                                  
      r=w*(y(:,ic)-q(:,ic)/sxp)                                         
      do 13711 j=1,ni                                                   
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))        
13711 continue                                                          
      continue                                                          
13701 continue                                                          
      continue                                                          
      do 13721 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 13741                                      
      al=ulam(ilm)                                                      
      goto 13731                                                        
13741 if(ilm .le. 2)goto 13751                                          
      al=al*alf                                                         
      goto 13731                                                        
13751 if(ilm .ne. 1)goto 13761                                          
      al=big                                                            
      goto 13771                                                        
13761 continue                                                          
      al0=0.0                                                           
      do 13781 j=1,ni                                                   
      if(ju(j).eq.0)goto 13781                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
13781 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
13771 continue                                                          
13731 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 13791 k=1,ni                                                   
      if(ixx(k).eq.1)goto 13791                                         
      if(ju(k).eq.0)goto 13791                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
13791 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
13801 continue                                                          
      ix=0                                                              
      jx=ix                                                             
      ig=0                                                              
      if(nlp .le. maxit)goto 13821                                      
      jerr=-ilm                                                         
      return                                                            
13821 continue                                                          
      do 13831 ic=1,nc                                                  
      bs(0,ic)=b(0,ic)                                                  
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                       
      xmz=0.0                                                           
      do 13841 i=1,no                                                   
      pic=q(i,ic)/sxp(i)                                                
      if(pic .ge. pfm)goto 13861                                        
      pic=0.0                                                           
      v(i)=0.0                                                          
      goto 13851                                                        
13861 if(pic .le. pfx)goto 13871                                        
      pic=1.0                                                           
      v(i)=0.0                                                          
      goto 13881                                                        
13871 continue                                                          
      v(i)=w(i)*pic*(1.0-pic)                                           
      xmz=xmz+v(i)                                                      
13881 continue                                                          
13851 continue                                                          
      r(i)=w(i)*(y(i,ic)-pic)                                           
13841 continue                                                          
      continue                                                          
      if(xmz.le.vmin)goto 13831                                         
      ig=1                                                              
      if(kopt .ne. 0)goto 13901                                         
      do 13911 j=1,ni                                                   
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                 
13911 continue                                                          
      continue                                                          
13901 continue                                                          
      continue                                                          
13921 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 13931 k=1,ni                                                   
      if(ixx(k).eq.0)goto 13931                                         
      bk=b(k,ic)                                                        
      gk=dot_product(r,x(:,k))                                          
      u=gk+xv(k,ic)*b(k,ic)                                             
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 13951                                         
      b(k,ic)=0.0                                                       
      goto 13961                                                        
13951 continue                                                          
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))
     *)
13961 continue                                                          
      continue                                                          
      d=b(k,ic)-bk                                                      
      if(abs(d).le.0.0)goto 13931                                       
      dlx=max(dlx,xv(k,ic)*d**2)                                        
      r=r-d*v*x(:,k)                                                    
      if(mm(k) .ne. 0)goto 13981                                        
      nin=nin+1                                                         
      if(nin .le. nx)goto 14001                                         
      jx=1                                                              
      goto 13932                                                        
14001 continue                                                          
      mm(k)=nin                                                         
      m(nin)=k                                                          
13981 continue                                                          
13931 continue                                                          
13932 continue                                                          
      if(jx.gt.0)goto 13922                                             
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 14021                                          
      b(0,ic)=b(0,ic)+d                                                 
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
14021 continue                                                          
      if(dlx.lt.shr)goto 13922                                          
      if(nlp .le. maxit)goto 14041                                      
      jerr=-ilm                                                         
      return                                                            
14041 continue                                                          
      continue                                                          
14051 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 14061 l=1,nin                                                  
      k=m(l)                                                            
      bk=b(k,ic)                                                        
      gk=dot_product(r,x(:,k))                                          
      u=gk+xv(k,ic)*b(k,ic)                                             
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 14081                                         
      b(k,ic)=0.0                                                       
      goto 14091                                                        
14081 continue                                                          
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))
     *)
14091 continue                                                          
      continue                                                          
      d=b(k,ic)-bk                                                      
      if(abs(d).le.0.0)goto 14061                                       
      dlx=max(dlx,xv(k,ic)*d**2)                                        
      r=r-d*v*x(:,k)                                                    
14061 continue                                                          
      continue                                                          
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 14111                                          
      b(0,ic)=b(0,ic)+d                                                 
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
14111 continue                                                          
      if(dlx.lt.shr)goto 14052                                          
      if(nlp .le. maxit)goto 14131                                      
      jerr=-ilm                                                         
      return                                                            
14131 continue                                                          
      goto 14051                                                        
14052 continue                                                          
      goto 13921                                                        
13922 continue                                                          
      if(jx.gt.0)goto 13832                                             
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                         
      if(ix .ne. 0)goto 14151                                           
      do 14161 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 14181             
      ix=1                                                              
      goto 14162                                                        
14181 continue                                                          
14161 continue                                                          
14162 continue                                                          
14151 continue                                                          
      do 14191 i=1,no                                                   
      fi=b(0,ic)+g(i,ic)                                                
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))      
      fi=min(max(exmn,fi),exmx)                                         
      sxp(i)=sxp(i)-q(i,ic)                                             
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                 
      sxp(i)=sxp(i)+q(i,ic)                                             
14191 continue                                                          
      continue                                                          
13831 continue                                                          
13832 continue                                                          
      s=-sum(b(0,:))/nc                                                 
      b(0,:)=b(0,:)+s                                                   
      di=s                                                              
      do 14201 j=1,nin                                                  
      l=m(j)                                                            
      if(vp(l) .gt. 0.0)goto 14221                                      
      s=sum(b(l,:))/nc                                                  
      goto 14231                                                        
14221 continue                                                          
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                  
14231 continue                                                          
      continue                                                          
      b(l,:)=b(l,:)-s                                                   
      di=di-s*x(:,l)                                                    
14201 continue                                                          
      continue                                                          
      di=exp(di)                                                        
      sxp=sxp*di                                                        
      do 14241 ic=1,nc                                                  
      q(:,ic)=q(:,ic)*di                                                
14241 continue                                                          
      continue                                                          
      if(jx.gt.0)goto 13802                                             
      if(ig.eq.0)goto 13802                                             
      if(ix .ne. 0)goto 14261                                           
      do 14271 k=1,ni                                                   
      if(ixx(k).eq.1)goto 14271                                         
      if(ju(k).eq.0)goto 14271                                          
      ga(k)=0.0                                                         
14271 continue                                                          
      continue                                                          
      do 14281 ic=1,nc                                                  
      r=w*(y(:,ic)-q(:,ic)/sxp)                                         
      do 14291 k=1,ni                                                   
      if(ixx(k).eq.1)goto 14291                                         
      if(ju(k).eq.0)goto 14291                                          
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                       
14291 continue                                                          
      continue                                                          
14281 continue                                                          
      continue                                                          
      do 14301 k=1,ni                                                   
      if(ixx(k).eq.1)goto 14301                                         
      if(ju(k).eq.0)goto 14301                                          
      if(ga(k) .le. al1*vp(k))goto 14321                                
      ixx(k)=1                                                          
      ix=1                                                              
14321 continue                                                          
14301 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10880                                           
      goto 13802                                                        
14261 continue                                                          
      goto 13801                                                        
13802 continue                                                          
      if(jx .le. 0)goto 14341                                           
      jerr=-10000-ilm                                                   
      goto 13722                                                        
14341 continue                                                          
      devi=0.0                                                          
      do 14351 ic=1,nc                                                  
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                       
      a0(ic,ilm)=b(0,ic)                                                
      do 14361 i=1,no                                                   
      if(y(i,ic).le.0.0)goto 14361                                      
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                        
14361 continue                                                          
      continue                                                          
14351 continue                                                          
      continue                                                          
      kin(ilm)=nin                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dev1-devi)/dev0                                         
      if(ig.eq.0)goto 13722                                             
      if(ilm.lt.mnl)goto 13721                                          
      if(flmin.ge.1.0)goto 13721                                        
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 13722          
      if(dev(ilm).gt.devmax)goto 13722                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13722                          
13721 continue                                                          
13722 continue                                                          
      g=log(q)                                                          
      do 14371 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
14371 continue                                                          
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
      do 14381 i=1,n                                                    
      s(i)=sum(e(i,:))                                                  
14381 continue                                                          
      continue                                                          
      continue                                                          
14391 continue                                                          
      dm=0.0                                                            
      do 14401 k=1,kk                                                   
      t=0.0                                                             
      u=t                                                               
      do 14411 i=1,n                                                    
      pik=e(i,k)/s(i)                                                   
      t=t+q(i)*(y(i,k)-pik)                                             
      u=u+q(i)*pik*(1.0-pik)                                            
14411 continue                                                          
      continue                                                          
      d=t/u                                                             
      az(k)=az(k)+d                                                     
      ed=exp(d)                                                         
      dm=max(dm,abs(d))                                                 
      do 14421 i=1,n                                                    
      z=e(i,k)                                                          
      e(i,k)=z*ed                                                       
      s(i)=s(i)-z+e(i,k)                                                
14421 continue                                                          
      continue                                                          
14401 continue                                                          
      continue                                                          
      if(dm.lt.eps)goto 14392                                           
      goto 14391                                                        
14392 continue                                                          
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
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 14441                    
      elc=am                                                            
      go to 14450                                                       
14441 continue                                                          
      do 14461 i=1,n                                                    
      m(i)=i                                                            
14461 continue                                                          
      continue                                                          
      call psort7(a,m,1,n)                                              
      if(a(m(1)) .ne. a(m(n)))goto 14481                                
      elc=a(1)                                                          
      go to 14450                                                       
14481 continue                                                          
      if(mod(n,2) .ne. 1)goto 14501                                     
      ad=a(m(n/2+1))                                                    
      goto 14511                                                        
14501 continue                                                          
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                    
14511 continue                                                          
      continue                                                          
      if(parm .ne. 1.0)goto 14531                                       
      elc=ad                                                            
      go to 14450                                                       
14531 continue                                                          
      b1=min(am,ad)                                                     
      b2=max(am,ad)                                                     
      k2=1                                                              
      continue                                                          
14541 if(a(m(k2)).gt.b1)goto 14542                                      
      k2=k2+1                                                           
      goto 14541                                                        
14542 continue                                                          
      k1=k2-1                                                           
      continue                                                          
14551 if(a(m(k2)).ge.b2)goto 14552                                      
      k2=k2+1                                                           
      goto 14551                                                        
14552 continue                                                          
      r=parm/((1.0-parm)*fn)                                            
      is=0                                                              
      sm=n-2*(k1-1)                                                     
      do 14561 k=k1,k2-1                                                
      sm=sm-2.0                                                         
      s=r*sm+am                                                         
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 14581                
      is=k                                                              
      goto 14562                                                        
14581 continue                                                          
14561 continue                                                          
14562 continue                                                          
      if(is .eq. 0)goto 14601                                           
      elc=s                                                             
      go to 14450                                                       
14601 continue                                                          
      r2=2.0*r                                                          
      s1=a(m(k1))                                                       
      am2=2.0*am                                                        
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                 
      elc=s1                                                            
      do 14611 k=k1+1,k2                                                
      s=a(m(k))                                                         
      if(s.eq.s1)goto 14611                                             
      c=r2*sum(abs(a-s))+s*(s-am2)                                      
      if(c .ge. cri)goto 14631                                          
      cri=c                                                             
      elc=s                                                             
14631 continue                                                          
      s1=s                                                              
14611 continue                                                          
      continue                                                          
14450 continue                                                          
      elc=max(maxval(a-cl(2)),min(minval(a-cl(1)),elc))                 
      return                                                            
      end                                                               
      function nintot(ni,nx,nc,a,m,nin,is)                              
      implicit double precision(a-h,o-z)                                
      double precision a(nx,nc)                                         
      integer m(nx),is(ni)                                              
      is=0                                                              
      nintot=0                                                          
      do 14641 ic=1,nc                                                  
      do 14651 j=1,nin                                                  
      k=m(j)                                                            
      if(is(k).ne.0)goto 14651                                          
      if(a(j,ic).eq.0.0)goto 14651                                      
      is(k)=k                                                           
      nintot=nintot+1                                                   
14651 continue                                                          
      continue                                                          
14641 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                          
      implicit double precision(a-h,o-z)                                
      double precision ca(nx,nc),a(ni,nc)                               
      integer ia(nx)                                                    
      a=0.0                                                             
      do 14661 ic=1,nc                                                  
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                         
14661 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                   
      implicit double precision(a-h,o-z)                                
      double precision a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)              
      integer ia(nx)                                                    
      do 14671 i=1,nt                                                   
      do 14681 ic=1,nc                                                  
      ans(ic,i)=a0(ic)                                                  
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1
     *:nin)))
14681 continue                                                          
      continue                                                          
14671 continue                                                          
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
      if(maxval(vp) .gt. 0.0)goto 14701                                 
      jerr=10000                                                        
      return                                                            
14701 continue                                                          
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
      if(kopt .ne. 2)goto 14721                                         
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
14721 continue                                                          
      call spchkvars(no,ni,x,ix,ju)                                     
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 14741                                   
      jerr=7777                                                         
      return                                                            
14741 continue                                                          
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      do 14751 i=1,no                                                   
      ww(i)=sum(y(i,:))                                                 
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                              
14751 continue                                                          
      continue                                                          
      sw=sum(ww)                                                        
      ww=ww/sw                                                          
      if(nc .ne. 1)goto 14771                                           
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)             
      if(isd .le. 0)goto 14791                                          
      do 14801 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
14801 continue                                                          
      continue                                                          
14791 continue                                                          
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,cl,ne,n
     *x,nlam,  flmin,ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin
     *,dev0,dev,  alm,nlp,jerr)
      goto 14761                                                        
14771 if(kopt .ne. 2)goto 14811                                         
      call multsplstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs,xv)      
      if(isd .le. 0)goto 14831                                          
      do 14841 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
14841 continue                                                          
      continue                                                          
14831 continue                                                          
      call multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nl
     *am,flmin,  ulam,thr,intr,maxit,xv,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,
     *alm,nlp,jerr)
      goto 14851                                                        
14811 continue                                                          
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)             
      if(isd .le. 0)goto 14871                                          
      do 14881 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
14881 continue                                                          
      continue                                                          
14871 continue                                                          
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,f
     *lmin,  ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,  ia,nin,dev0,
     *dev,alm,nlp,jerr)
14851 continue                                                          
14761 continue                                                          
      if(jerr.gt.0) return                                              
      dev0=2.0*sw*dev0                                                  
      do 14891 k=1,lmu                                                  
      nk=nin(k)                                                         
      do 14901 ic=1,nc                                                  
      if(isd .le. 0)goto 14921                                          
      do 14931 l=1,nk                                                   
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                   
14931 continue                                                          
      continue                                                          
14921 continue                                                          
      if(intr .ne. 0)goto 14951                                         
      a0(ic,k)=0.0                                                      
      goto 14961                                                        
14951 continue                                                          
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))         
14961 continue                                                          
      continue                                                          
14901 continue                                                          
      continue                                                          
14891 continue                                                          
      continue                                                          
      deallocate(ww,ju,vq,xm,xs)                                        
      if(kopt.eq.2) deallocate(xv)                                      
      return                                                            
      end                                                               
      subroutine multsplstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs,xv) 
      implicit double precision(a-h,o-z)                                
      double precision x(*),w(no),xm(ni),xs(ni),xv(ni)                  
      integer ix(*),jx(*),ju(ni)                                        
      if(intr .ne. 0)goto 14981                                         
      do 14991 j=1,ni                                                   
      if(ju(j).eq.0)goto 14991                                          
      xm(j)=0.0                                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                       
      if(isd .eq. 0)goto 15011                                          
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                         
      vc=xv(j)-xbq                                                      
      xs(j)=sqrt(vc)                                                    
      xv(j)=1.0+xbq/vc                                                  
      goto 15021                                                        
15011 continue                                                          
      xs(j)=1.0                                                         
15021 continue                                                          
      continue                                                          
14991 continue                                                          
      continue                                                          
      return                                                            
14981 continue                                                          
      do 15031 j=1,ni                                                   
      if(ju(j).eq.0)goto 15031                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2              
      if(isd .le. 0)goto 15051                                          
      xs(j)=sqrt(xv(j))                                                 
      xv(j)=1.0                                                         
15051 continue                                                          
15031 continue                                                          
      continue                                                          
      if(isd.eq.0) xs=1.0                                               
      return                                                            
      end                                                               
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs)        
      implicit double precision(a-h,o-z)                                
      double precision x(*),w(no),xm(ni),xs(ni)                         
      integer ix(*),jx(*),ju(ni)                                        
      if(intr .ne. 0)goto 15071                                         
      do 15081 j=1,ni                                                   
      if(ju(j).eq.0)goto 15081                                          
      xm(j)=0.0                                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      if(isd .eq. 0)goto 15101                                          
      vc=dot_product(w(jx(jb:je)),x(jb:je)**2)  -dot_product(w(jx(jb:je)
     *),x(jb:je))**2
      xs(j)=sqrt(vc)                                                    
      goto 15111                                                        
15101 continue                                                          
      xs(j)=1.0                                                         
15111 continue                                                          
      continue                                                          
15081 continue                                                          
      continue                                                          
      return                                                            
15071 continue                                                          
      do 15121 j=1,ni                                                   
      if(ju(j).eq.0)goto 15121                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j
     *)**2)
15121 continue                                                          
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
      if(q0 .gt. pmin)goto 15141                                        
      jerr=8001                                                         
      return                                                            
15141 continue                                                          
      if(q0 .lt. 1.0-pmin)goto 15161                                    
      jerr=9001                                                         
      return                                                            
15161 continue                                                          
      if(intr.eq.0) q0=0.5                                              
      bz=0.0                                                            
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                 
      if(nonzero(no,g) .ne. 0)goto 15181                                
      vi=q0*(1.0-q0)                                                    
      b(0)=bz                                                           
      v=vi*w                                                            
      r=w*(y-q0)                                                        
      q=q0                                                              
      xm(0)=vi                                                          
      dev1=-(bz*q0+log(1.0-q0))                                         
      goto 15191                                                        
15181 continue                                                          
      b(0)=0.0                                                          
      if(intr .eq. 0)goto 15211                                         
      b(0)=azero(no,y,g,w,jerr)                                         
      if(jerr.ne.0) return                                              
15211 continue                                                          
      q=1.0/(1.0+exp(-b(0)-g))                                          
      v=w*q*(1.0-q)                                                     
      r=w*(y-q)                                                         
      xm(0)=sum(v)                                                      
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                     
15191 continue                                                          
      continue                                                          
      if(kopt .le. 0)goto 15231                                         
      if(isd .le. 0 .or. intr .eq. 0)goto 15251                         
      xv=0.25                                                           
      goto 15261                                                        
15251 continue                                                          
      do 15271 j=1,ni                                                   
      if(ju(j).eq.0)goto 15271                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)       
15271 continue                                                          
      continue                                                          
15261 continue                                                          
      continue                                                          
15231 continue                                                          
      b(1:ni)=0.0                                                       
      dev0=dev1                                                         
      do 15281 i=1,no                                                   
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                     
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))           
15281 continue                                                          
      continue                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 15301                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
15301 continue                                                          
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
      do 15311 j=1,ni                                                   
      if(ju(j).eq.0)goto 15311                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      jn=ix(j+1)-ix(j)                                                  
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                              
      gj=dot_product(sc(1:jn),x(jb:je))                                 
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                   
15311 continue                                                          
      continue                                                          
      do 15321 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 15341                                      
      al=ulam(ilm)                                                      
      goto 15331                                                        
15341 if(ilm .le. 2)goto 15351                                          
      al=al*alf                                                         
      goto 15331                                                        
15351 if(ilm .ne. 1)goto 15361                                          
      al=big                                                            
      goto 15371                                                        
15361 continue                                                          
      al0=0.0                                                           
      do 15381 j=1,ni                                                   
      if(ju(j).eq.0)goto 15381                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
15381 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
15371 continue                                                          
15331 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 15391 k=1,ni                                                   
      if(ixx(k).eq.1)goto 15391                                         
      if(ju(k).eq.0)goto 15391                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
15391 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
15401 continue                                                          
      if(nlp .le. maxit)goto 15421                                      
      jerr=-ilm                                                         
      return                                                            
15421 continue                                                          
      bs(0)=b(0)                                                        
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                             
      do 15431 j=1,ni                                                   
      if(ixx(j).eq.0)goto 15431                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      jn=ix(j+1)-ix(j)                                                  
      sc(1:jn)=v(jx(jb:je))                                             
      xm(j)=dot_product(sc(1:jn),x(jb:je))                              
      if(kopt .ne. 0)goto 15451                                         
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                           
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2             
15451 continue                                                          
15431 continue                                                          
      continue                                                          
      continue                                                          
15461 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 15471 k=1,ni                                                   
      if(ixx(k).eq.0)goto 15471                                         
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      jn=ix(k+1)-ix(k)                                                  
      bk=b(k)                                                           
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                              
      gk=dot_product(sc(1:jn),x(jb:je))                                 
      gk=(gk-svr*xb(k))/xs(k)                                           
      u=gk+xv(k)*b(k)                                                   
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 15491                                         
      b(k)=0.0                                                          
      goto 15501                                                        
15491 continue                                                          
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))       
15501 continue                                                          
      continue                                                          
      d=b(k)-bk                                                         
      if(abs(d).le.0.0)goto 15471                                       
      dlx=max(dlx,xv(k)*d**2)                                           
      if(mm(k) .ne. 0)goto 15521                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 15472                                           
      mm(k)=nin                                                         
      m(nin)=k                                                          
      sc(1:jn)=v(jx(jb:je))                                             
      xm(k)=dot_product(sc(1:jn),x(jb:je))                              
15521 continue                                                          
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)           
      o=o+d*(xb(k)/xs(k))                                               
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                               
15471 continue                                                          
15472 continue                                                          
      if(nin.gt.nx)goto 15462                                           
      d=0.0                                                             
      if(intr.ne.0) d=svr/xm(0)                                         
      if(d .eq. 0.0)goto 15541                                          
      b(0)=b(0)+d                                                       
      dlx=max(dlx,xm(0)*d**2)                                           
      r=r-d*v                                                           
      svr=svr-d*xm(0)                                                   
15541 continue                                                          
      if(dlx.lt.shr)goto 15462                                          
      if(nlp .le. maxit)goto 15561                                      
      jerr=-ilm                                                         
      return                                                            
15561 continue                                                          
      continue                                                          
15571 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 15581 l=1,nin                                                  
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
      if(au .gt. 0.0)goto 15601                                         
      b(k)=0.0                                                          
      goto 15611                                                        
15601 continue                                                          
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))       
15611 continue                                                          
      continue                                                          
      d=b(k)-bk                                                         
      if(abs(d).le.0.0)goto 15581                                       
      dlx=max(dlx,xv(k)*d**2)                                           
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)           
      o=o+d*(xb(k)/xs(k))                                               
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                               
15581 continue                                                          
      continue                                                          
      d=0.0                                                             
      if(intr.ne.0) d=svr/xm(0)                                         
      if(d .eq. 0.0)goto 15631                                          
      b(0)=b(0)+d                                                       
      dlx=max(dlx,xm(0)*d**2)                                           
      r=r-d*v                                                           
      svr=svr-d*xm(0)                                                   
15631 continue                                                          
      if(dlx.lt.shr)goto 15572                                          
      if(nlp .le. maxit)goto 15651                                      
      jerr=-ilm                                                         
      return                                                            
15651 continue                                                          
      goto 15571                                                        
15572 continue                                                          
      goto 15461                                                        
15462 continue                                                          
      if(nin.gt.nx)goto 15402                                           
      sc=b(0)                                                           
      b0=0.0                                                            
      do 15661 j=1,nin                                                  
      l=m(j)                                                            
      jb=ix(l)                                                          
      je=ix(l+1)-1                                                      
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                   
      b0=b0-b(l)*xb(l)/xs(l)                                            
15661 continue                                                          
      continue                                                          
      sc=sc+b0                                                          
      do 15671 i=1,no                                                   
      fi=sc(i)+g(i)                                                     
      if(fi .ge. fmin)goto 15691                                        
      q(i)=0.0                                                          
      goto 15681                                                        
15691 if(fi .le. fmax)goto 15701                                        
      q(i)=1.0                                                          
      goto 15711                                                        
15701 continue                                                          
      q(i)=1.0/(1.0+exp(-fi))                                           
15711 continue                                                          
15681 continue                                                          
15671 continue                                                          
      continue                                                          
      v=w*q*(1.0-q)                                                     
      xm(0)=sum(v)                                                      
      if(xm(0).lt.vmin)goto 15402                                       
      r=w*(y-q)                                                         
      svr=sum(r)                                                        
      o=0.0                                                             
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 15731                      
      kx=0                                                              
      do 15741 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 15741                        
      kx=1                                                              
      goto 15742                                                        
15741 continue                                                          
15742 continue                                                          
      if(kx .ne. 0)goto 15761                                           
      do 15771 j=1,ni                                                   
      if(ixx(j).eq.1)goto 15771                                         
      if(ju(j).eq.0)goto 15771                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      jn=ix(j+1)-ix(j)                                                  
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                              
      gj=dot_product(sc(1:jn),x(jb:je))                                 
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                   
      if(ga(j) .le. al1*vp(j))goto 15791                                
      ixx(j)=1                                                          
      kx=1                                                              
15791 continue                                                          
15771 continue                                                          
      continue                                                          
      if(kx.eq.1) go to 10880                                           
      goto 15402                                                        
15761 continue                                                          
15731 continue                                                          
      goto 15401                                                        
15402 continue                                                          
      if(nin .le. nx)goto 15811                                         
      jerr=-10000-ilm                                                   
      goto 15322                                                        
15811 continue                                                          
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                             
      kin(ilm)=nin                                                      
      a0(ilm)=b(0)                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      devi=dev2(no,w,y,q,pmin)                                          
      dev(ilm)=(dev1-devi)/dev0                                         
      if(ilm.lt.mnl)goto 15321                                          
      if(flmin.ge.1.0)goto 15321                                        
      me=0                                                              
      do 15821 j=1,nin                                                  
      if(a(j,ilm).ne.0.0) me=me+1                                       
15821 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 15322                                            
      if(dev(ilm).gt.devmax)goto 15322                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 15322                          
      if(xm(0).lt.vmin)goto 15322                                       
15321 continue                                                          
15322 continue                                                          
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
      do 15831 ic=1,nc                                                  
      q0=dot_product(w,y(:,ic))                                         
      if(q0 .gt. pmin)goto 15851                                        
      jerr =8000+ic                                                     
      return                                                            
15851 continue                                                          
      if(q0 .lt. 1.0-pmin)goto 15871                                    
      jerr =9000+ic                                                     
      return                                                            
15871 continue                                                          
      if(intr.eq.0) q0=1.0/nc                                           
      b(1:ni,ic)=0.0                                                    
      b(0,ic)=0.0                                                       
      if(intr .eq. 0)goto 15891                                         
      b(0,ic)=log(q0)                                                   
      dev1=dev1-q0*b(0,ic)                                              
15891 continue                                                          
15831 continue                                                          
      continue                                                          
      if(intr.eq.0) dev1=log(float(nc))                                 
      iy=0                                                              
      al=0.0                                                            
      if(nonzero(no*nc,g) .ne. 0)goto 15911                             
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                      
      sxp=0.0                                                           
      do 15921 ic=1,nc                                                  
      q(:,ic)=exp(b(0,ic))                                              
      sxp=sxp+q(:,ic)                                                   
15921 continue                                                          
      continue                                                          
      goto 15931                                                        
15911 continue                                                          
      do 15941 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
15941 continue                                                          
      continue                                                          
      sxp=0.0                                                           
      if(intr .ne. 0)goto 15961                                         
      b(0,:)=0.0                                                        
      goto 15971                                                        
15961 continue                                                          
      call kazero(nc,no,y,g,w,b(0,:),jerr)                              
      if(jerr.ne.0) return                                              
15971 continue                                                          
      continue                                                          
      dev1=0.0                                                          
      do 15981 ic=1,nc                                                  
      q(:,ic)=b(0,ic)+g(:,ic)                                           
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                          
      q(:,ic)=exp(q(:,ic))                                              
      sxp=sxp+q(:,ic)                                                   
15981 continue                                                          
      continue                                                          
      sxpl=w*log(sxp)                                                   
      do 15991 ic=1,nc                                                  
      dev1=dev1+dot_product(y(:,ic),sxpl)                               
15991 continue                                                          
      continue                                                          
15931 continue                                                          
      continue                                                          
      do 16001 ic=1,nc                                                  
      do 16011 i=1,no                                                   
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))            
16011 continue                                                          
      continue                                                          
16001 continue                                                          
      continue                                                          
      dev0=dev0+dev1                                                    
      if(kopt .le. 0)goto 16031                                         
      if(isd .le. 0 .or. intr .eq. 0)goto 16051                         
      xv=0.25                                                           
      goto 16061                                                        
16051 continue                                                          
      do 16071 j=1,ni                                                   
      if(ju(j).eq.0)goto 16071                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)     
16071 continue                                                          
      continue                                                          
16061 continue                                                          
      continue                                                          
16031 continue                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 16091                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
16091 continue                                                          
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
      do 16101 ic=1,nc                                                  
      v=q(:,ic)/sxp                                                     
      r=w*(y(:,ic)-v)                                                   
      v=w*v*(1.0-v)                                                     
      do 16111 j=1,ni                                                   
      if(ju(j).eq.0)goto 16111                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      jn=ix(j+1)-ix(j)                                                  
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                              
      gj=dot_product(sc(1:jn),x(jb:je))                                 
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                          
16111 continue                                                          
      continue                                                          
16101 continue                                                          
      continue                                                          
      do 16121 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 16141                                      
      al=ulam(ilm)                                                      
      goto 16131                                                        
16141 if(ilm .le. 2)goto 16151                                          
      al=al*alf                                                         
      goto 16131                                                        
16151 if(ilm .ne. 1)goto 16161                                          
      al=big                                                            
      goto 16171                                                        
16161 continue                                                          
      al0=0.0                                                           
      do 16181 j=1,ni                                                   
      if(ju(j).eq.0)goto 16181                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
16181 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
16171 continue                                                          
16131 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 16191 k=1,ni                                                   
      if(iy(k).eq.1)goto 16191                                          
      if(ju(k).eq.0)goto 16191                                          
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                   
16191 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
16201 continue                                                          
      ixx=0                                                             
      jxx=ixx                                                           
      ig=0                                                              
      if(nlp .le. maxit)goto 16221                                      
      jerr=-ilm                                                         
      return                                                            
16221 continue                                                          
      do 16231 ic=1,nc                                                  
      bs(0,ic)=b(0,ic)                                                  
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                       
      xm(0)=0.0                                                         
      svr=0.0                                                           
      o=0.0                                                             
      do 16241 i=1,no                                                   
      pic=q(i,ic)/sxp(i)                                                
      if(pic .ge. pfm)goto 16261                                        
      pic=0.0                                                           
      v(i)=0.0                                                          
      goto 16251                                                        
16261 if(pic .le. pfx)goto 16271                                        
      pic=1.0                                                           
      v(i)=0.0                                                          
      goto 16281                                                        
16271 continue                                                          
      v(i)=w(i)*pic*(1.0-pic)                                           
      xm(0)=xm(0)+v(i)                                                  
16281 continue                                                          
16251 continue                                                          
      r(i)=w(i)*(y(i,ic)-pic)                                           
      svr=svr+r(i)                                                      
16241 continue                                                          
      continue                                                          
      if(xm(0).le.vmin)goto 16231                                       
      ig=1                                                              
      do 16291 j=1,ni                                                   
      if(iy(j).eq.0)goto 16291                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                          
      if(kopt .ne. 0)goto 16311                                         
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                    
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2       
16311 continue                                                          
16291 continue                                                          
      continue                                                          
      continue                                                          
16321 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 16331 k=1,ni                                                   
      if(iy(k).eq.0)goto 16331                                          
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      jn=ix(k+1)-ix(k)                                                  
      bk=b(k,ic)                                                        
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                              
      gk=dot_product(sc(1:jn),x(jb:je))                                 
      gk=(gk-svr*xb(k))/xs(k)                                           
      u=gk+xv(k,ic)*b(k,ic)                                             
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 16351                                         
      b(k,ic)=0.0                                                       
      goto 16361                                                        
16351 continue                                                          
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))
     *)
16361 continue                                                          
      continue                                                          
      d=b(k,ic)-bk                                                      
      if(abs(d).le.0.0)goto 16331                                       
      dlx=max(dlx,xv(k,ic)*d**2)                                        
      if(mm(k) .ne. 0)goto 16381                                        
      nin=nin+1                                                         
      if(nin .le. nx)goto 16401                                         
      jxx=1                                                             
      goto 16332                                                        
16401 continue                                                          
      mm(k)=nin                                                         
      m(nin)=k                                                          
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                          
16381 continue                                                          
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)           
      o=o+d*(xb(k)/xs(k))                                               
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                               
16331 continue                                                          
16332 continue                                                          
      if(jxx.gt.0)goto 16322                                            
      d=0.0                                                             
      if(intr.ne.0) d=svr/xm(0)                                         
      if(d .eq. 0.0)goto 16421                                          
      b(0,ic)=b(0,ic)+d                                                 
      dlx=max(dlx,xm(0)*d**2)                                           
      r=r-d*v                                                           
      svr=svr-d*xm(0)                                                   
16421 continue                                                          
      if(dlx.lt.shr)goto 16322                                          
      if(nlp .le. maxit)goto 16441                                      
      jerr=-ilm                                                         
      return                                                            
16441 continue                                                          
      continue                                                          
16451 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 16461 l=1,nin                                                  
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
      if(au .gt. 0.0)goto 16481                                         
      b(k,ic)=0.0                                                       
      goto 16491                                                        
16481 continue                                                          
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))
     *)
16491 continue                                                          
      continue                                                          
      d=b(k,ic)-bk                                                      
      if(abs(d).le.0.0)goto 16461                                       
      dlx=max(dlx,xv(k,ic)*d**2)                                        
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)           
      o=o+d*(xb(k)/xs(k))                                               
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                               
16461 continue                                                          
      continue                                                          
      d=0.0                                                             
      if(intr.ne.0) d=svr/xm(0)                                         
      if(d .eq. 0.0)goto 16511                                          
      b(0,ic)=b(0,ic)+d                                                 
      dlx=max(dlx,xm(0)*d**2)                                           
      r=r-d*v                                                           
      svr=svr-d*xm(0)                                                   
16511 continue                                                          
      if(dlx.lt.shr)goto 16452                                          
      if(nlp .le. maxit)goto 16531                                      
      jerr=-ilm                                                         
      return                                                            
16531 continue                                                          
      goto 16451                                                        
16452 continue                                                          
      goto 16321                                                        
16322 continue                                                          
      if(jxx.gt.0)goto 16232                                            
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                      
      if(ixx .ne. 0)goto 16551                                          
      do 16561 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 16581             
      ixx=1                                                             
      goto 16562                                                        
16581 continue                                                          
16561 continue                                                          
16562 continue                                                          
16551 continue                                                          
      sc=b(0,ic)+g(:,ic)                                                
      b0=0.0                                                            
      do 16591 j=1,nin                                                  
      l=m(j)                                                            
      jb=ix(l)                                                          
      je=ix(l+1)-1                                                      
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                
      b0=b0-b(l,ic)*xb(l)/xs(l)                                         
16591 continue                                                          
      continue                                                          
      sc=min(max(exmn,sc+b0),exmx)                                      
      sxp=sxp-q(:,ic)                                                   
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                       
      sxp=sxp+q(:,ic)                                                   
16231 continue                                                          
16232 continue                                                          
      s=-sum(b(0,:))/nc                                                 
      b(0,:)=b(0,:)+s                                                   
      sc=s                                                              
      b0=0.0                                                            
      do 16601 j=1,nin                                                  
      l=m(j)                                                            
      if(vp(l) .gt. 0.0)goto 16621                                      
      s=sum(b(l,:))/nc                                                  
      goto 16631                                                        
16621 continue                                                          
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                  
16631 continue                                                          
      continue                                                          
      b(l,:)=b(l,:)-s                                                   
      jb=ix(l)                                                          
      je=ix(l+1)-1                                                      
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                      
      b0=b0+s*xb(l)/xs(l)                                               
16601 continue                                                          
      continue                                                          
      sc=sc+b0                                                          
      sc=exp(sc)                                                        
      sxp=sxp*sc                                                        
      do 16641 ic=1,nc                                                  
      q(:,ic)=q(:,ic)*sc                                                
16641 continue                                                          
      continue                                                          
      if(jxx.gt.0)goto 16202                                            
      if(ig.eq.0)goto 16202                                             
      if(ixx .ne. 0)goto 16661                                          
      do 16671 j=1,ni                                                   
      if(iy(j).eq.1)goto 16671                                          
      if(ju(j).eq.0)goto 16671                                          
      ga(j)=0.0                                                         
16671 continue                                                          
      continue                                                          
      do 16681 ic=1,nc                                                  
      v=q(:,ic)/sxp                                                     
      r=w*(y(:,ic)-v)                                                   
      v=w*v*(1.0-v)                                                     
      do 16691 j=1,ni                                                   
      if(iy(j).eq.1)goto 16691                                          
      if(ju(j).eq.0)goto 16691                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      jn=ix(j+1)-ix(j)                                                  
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                              
      gj=dot_product(sc(1:jn),x(jb:je))                                 
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                          
16691 continue                                                          
      continue                                                          
16681 continue                                                          
      continue                                                          
      do 16701 k=1,ni                                                   
      if(iy(k).eq.1)goto 16701                                          
      if(ju(k).eq.0)goto 16701                                          
      if(ga(k) .le. al1*vp(k))goto 16721                                
      iy(k)=1                                                           
      ixx=1                                                             
16721 continue                                                          
16701 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10880                                          
      goto 16202                                                        
16661 continue                                                          
      goto 16201                                                        
16202 continue                                                          
      if(jxx .le. 0)goto 16741                                          
      jerr=-10000-ilm                                                   
      goto 16122                                                        
16741 continue                                                          
      devi=0.0                                                          
      do 16751 ic=1,nc                                                  
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                       
      a0(ic,ilm)=b(0,ic)                                                
      do 16761 i=1,no                                                   
      if(y(i,ic).le.0.0)goto 16761                                      
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                        
16761 continue                                                          
      continue                                                          
16751 continue                                                          
      continue                                                          
      kin(ilm)=nin                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dev1-devi)/dev0                                         
      if(ig.eq.0)goto 16122                                             
      if(ilm.lt.mnl)goto 16121                                          
      if(flmin.ge.1.0)goto 16121                                        
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 16122          
      if(dev(ilm).gt.devmax)goto 16122                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 16122                          
16121 continue                                                          
16122 continue                                                          
      g=log(q)                                                          
      do 16771 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
16771 continue                                                          
      continue                                                          
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy)                   
      return                                                            
      end                                                               
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)               
      implicit double precision(a-h,o-z)                                
      double precision a0(nc),ca(nx,nc),x(*),f(nc,n)                    
      integer ia(*),ix(*),jx(*)                                         
      do 16781 ic=1,nc                                                  
      f(ic,:)=a0(ic)                                                    
16781 continue                                                          
      continue                                                          
      do 16791 j=1,nin                                                  
      k=ia(j)                                                           
      kb=ix(k)                                                          
      ke=ix(k+1)-1                                                      
      do 16801 ic=1,nc                                                  
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                 
16801 continue                                                          
      continue                                                          
16791 continue                                                          
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
      if(maxval(vp) .gt. 0.0)goto 16821                                 
      jerr=10000                                                        
      return                                                            
16821 continue                                                          
      allocate(ww(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      if(isd .le. 0)goto 16841                                          
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
16841 continue                                                          
      call chkvars(no,ni,x,ju)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 16861                                   
      jerr=7777                                                         
      return                                                            
16861 continue                                                          
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      ww=max(0d0,w)                                                     
      sw=sum(ww)                                                        
      if(sw .gt. 0.0)goto 16881                                         
      jerr=9999                                                         
      return                                                            
16881 continue                                                          
      ww=ww/sw                                                          
      call cstandard(no,ni,x,ww,ju,isd,xs)                              
      if(isd .le. 0)goto 16901                                          
      do 16911 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
16911 continue                                                          
      continue                                                          
16901 continue                                                          
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,
     *thr,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                              
      dev0=2.0*sw*dev0                                                  
      if(isd .le. 0)goto 16931                                          
      do 16941 k=1,lmu                                                  
      nk=nin(k)                                                         
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                
16941 continue                                                          
      continue                                                          
16931 continue                                                          
      deallocate(ww,ju,vq)                                              
      if(isd.gt.0) deallocate(xs)                                       
      return                                                            
      end                                                               
      subroutine cstandard(no,ni,x,w,ju,isd,xs)                         
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),w(no),xs(ni)                            
      integer ju(ni)                                                    
      do 16951 j=1,ni                                                   
      if(ju(j).eq.0)goto 16951                                          
      xm=dot_product(w,x(:,j))                                          
      x(:,j)=x(:,j)-xm                                                  
      if(isd .le. 0)goto 16971                                          
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                              
      x(:,j)=x(:,j)/xs(j)                                               
16971 continue                                                          
16951 continue                                                          
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
      if(jerr.ne.0) go to 12220                                         
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
      if(nonzero(no,g) .eq. 0)goto 16991                                
      f=g-dot_product(q,g)                                              
      e=q*exp(sign(min(abs(f),fmax),f))                                 
      goto 17001                                                        
16991 continue                                                          
      f=0.0                                                             
      e=q                                                               
17001 continue                                                          
      continue                                                          
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                              
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                      
      dev0=rr                                                           
      do 17011 i=1,no                                                   
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 17031                
      w(i)=0.0                                                          
      wr(i)=w(i)                                                        
17031 continue                                                          
17011 continue                                                          
      continue                                                          
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                      
      if(jerr.ne.0) go to 12220                                         
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 17051                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
17051 continue                                                          
      m=0                                                               
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      mnl=min(mnlam,nlam)                                               
      as=0.0                                                            
      cthr=cthri*dev0                                                   
      do 17061 j=1,ni                                                   
      if(ju(j).eq.0)goto 17061                                          
      ga(j)=abs(dot_product(wr,x(:,j)))                                 
17061 continue                                                          
      continue                                                          
      do 17071 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 17091                                      
      al=ulam(ilm)                                                      
      goto 17081                                                        
17091 if(ilm .le. 2)goto 17101                                          
      al=al*alf                                                         
      goto 17081                                                        
17101 if(ilm .ne. 1)goto 17111                                          
      al=big                                                            
      goto 17121                                                        
17111 continue                                                          
      al0=0.0                                                           
      do 17131 j=1,ni                                                   
      if(ju(j).eq.0)goto 17131                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
17131 continue                                                          
      continue                                                          
      al0=al0/max(parm,1.0d-3)                                          
      al=alf*al0                                                        
17121 continue                                                          
17081 continue                                                          
      sa=alpha*al                                                       
      omal=oma*al                                                       
      tlam=alpha*(2.0*al-al0)                                           
      do 17141 k=1,ni                                                   
      if(ixx(k).eq.1)goto 17141                                         
      if(ju(k).eq.0)goto 17141                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
17141 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
17151 continue                                                          
      if(nlp .le. maxit)goto 17171                                      
      jerr=-ilm                                                         
      return                                                            
17171 continue                                                          
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                             
      call vars(no,ni,x,w,ixx,v)                                        
      continue                                                          
17181 continue                                                          
      nlp=nlp+1                                                         
      dli=0.0                                                           
      do 17191 j=1,ni                                                   
      if(ixx(j).eq.0)goto 17191                                         
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                
      if(abs(u) .gt. vp(j)*sa)goto 17211                                
      at=0.0                                                            
      goto 17221                                                        
17211 continue                                                          
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o
     *mal)))
17221 continue                                                          
      continue                                                          
      if(at .eq. a(j))goto 17241                                        
      del=at-a(j)                                                       
      a(j)=at                                                           
      dli=max(dli,v(j)*del**2)                                          
      wr=wr-del*w*x(:,j)                                                
      f=f+del*x(:,j)                                                    
      if(mm(j) .ne. 0)goto 17261                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 17192                                           
      mm(j)=nin                                                         
      m(nin)=j                                                          
17261 continue                                                          
17241 continue                                                          
17191 continue                                                          
17192 continue                                                          
      if(nin.gt.nx)goto 17182                                           
      if(dli.lt.cthr)goto 17182                                         
      if(nlp .le. maxit)goto 17281                                      
      jerr=-ilm                                                         
      return                                                            
17281 continue                                                          
      continue                                                          
17291 continue                                                          
      nlp=nlp+1                                                         
      dli=0.0                                                           
      do 17301 l=1,nin                                                  
      j=m(l)                                                            
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                
      if(abs(u) .gt. vp(j)*sa)goto 17321                                
      at=0.0                                                            
      goto 17331                                                        
17321 continue                                                          
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o
     *mal)))
17331 continue                                                          
      continue                                                          
      if(at .eq. a(j))goto 17351                                        
      del=at-a(j)                                                       
      a(j)=at                                                           
      dli=max(dli,v(j)*del**2)                                          
      wr=wr-del*w*x(:,j)                                                
      f=f+del*x(:,j)                                                    
17351 continue                                                          
17301 continue                                                          
      continue                                                          
      if(dli.lt.cthr)goto 17292                                         
      if(nlp .le. maxit)goto 17371                                      
      jerr=-ilm                                                         
      return                                                            
17371 continue                                                          
      goto 17291                                                        
17292 continue                                                          
      goto 17181                                                        
17182 continue                                                          
      if(nin.gt.nx)goto 17152                                           
      e=q*exp(sign(min(abs(f),fmax),f))                                 
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                      
      if(jerr .eq. 0)goto 17391                                         
      jerr=jerr-ilm                                                     
      go to 12220                                                       
17391 continue                                                          
      ix=0                                                              
      do 17401 j=1,nin                                                  
      k=m(j)                                                            
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 17401                        
      ix=1                                                              
      goto 17402                                                        
17401 continue                                                          
17402 continue                                                          
      if(ix .ne. 0)goto 17421                                           
      do 17431 k=1,ni                                                   
      if(ixx(k).eq.1)goto 17431                                         
      if(ju(k).eq.0)goto 17431                                          
      ga(k)=abs(dot_product(wr,x(:,k)))                                 
      if(ga(k) .le. sa*vp(k))goto 17451                                 
      ixx(k)=1                                                          
      ix=1                                                              
17451 continue                                                          
17431 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10880                                           
      goto 17152                                                        
17421 continue                                                          
      goto 17151                                                        
17152 continue                                                          
      if(nin .le. nx)goto 17471                                         
      jerr=-10000-ilm                                                   
      goto 17072                                                        
17471 continue                                                          
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                            
      kin(ilm)=nin                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                
      if(ilm.lt.mnl)goto 17071                                          
      if(flmin.ge.1.0)goto 17071                                        
      me=0                                                              
      do 17481 j=1,nin                                                  
      if(ao(j,ilm).ne.0.0) me=me+1                                      
17481 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 17072                                            
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17072           
      if(dev(ilm).gt.devmax)goto 17072                                  
17071 continue                                                          
17072 continue                                                          
      g=f                                                               
12220 continue                                                          
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)           
      return                                                            
      end                                                               
      subroutine cxmodval(ca,ia,nin,n,x,f)                              
      implicit double precision(a-h,o-z)                                
      double precision ca(nin),x(n,*),f(n)                              
      integer ia(nin)                                                   
      f=0.0                                                             
      if(nin.le.0) return                                               
      do 17491 i=1,n                                                    
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                   
17491 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                      
      implicit double precision(a-h,o-z)                                
      double precision y(no),d(no),q(no)                                
      integer jp(no),kp(*)                                              
      do 17501 j=1,no                                                   
      jp(j)=j                                                           
17501 continue                                                          
      continue                                                          
      call psort7(y,jp,1,no)                                            
      nj=0                                                              
      do 17511 j=1,no                                                   
      if(q(jp(j)).le.0.0)goto 17511                                     
      nj=nj+1                                                           
      jp(nj)=jp(j)                                                      
17511 continue                                                          
      continue                                                          
      if(nj .ne. 0)goto 17531                                           
      jerr=20000                                                        
      return                                                            
17531 continue                                                          
      j=1                                                               
      continue                                                          
17541 if(d(jp(j)).gt.0.0)goto 17542                                     
      j=j+1                                                             
      if(j.gt.nj)goto 17542                                             
      goto 17541                                                        
17542 continue                                                          
      if(j .lt. nj-1)goto 17561                                         
      jerr=30000                                                        
      return                                                            
17561 continue                                                          
      t0=y(jp(j))                                                       
      j0=j-1                                                            
      if(j0 .le. 0)goto 17581                                           
      continue                                                          
17591 if(y(jp(j0)).lt.t0)goto 17592                                     
      j0=j0-1                                                           
      if(j0.eq.0)goto 17592                                             
      goto 17591                                                        
17592 continue                                                          
      if(j0 .le. 0)goto 17611                                           
      nj=nj-j0                                                          
      do 17621 j=1,nj                                                   
      jp(j)=jp(j+j0)                                                    
17621 continue                                                          
      continue                                                          
17611 continue                                                          
17581 continue                                                          
      jerr=0                                                            
      nk=0                                                              
      yk=t0                                                             
      j=2                                                               
      continue                                                          
17631 continue                                                          
      continue                                                          
17641 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 17642                  
      j=j+1                                                             
      if(j.gt.nj)goto 17642                                             
      goto 17641                                                        
17642 continue                                                          
      nk=nk+1                                                           
      kp(nk)=j-1                                                        
      if(j.gt.nj)goto 17632                                             
      if(j .ne. nj)goto 17661                                           
      nk=nk+1                                                           
      kp(nk)=nj                                                         
      goto 17632                                                        
17661 continue                                                          
      yk=y(jp(j))                                                       
      j=j+1                                                             
      goto 17631                                                        
17632 continue                                                          
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
      do 17671 j=1,kp(1)                                                
      i=jp(j)                                                           
      w(i)=e(i)*(b-e(i)*c)                                              
      if(w(i) .gt. 0.0)goto 17691                                       
      jerr=-30000                                                       
      return                                                            
17691 continue                                                          
      wr(i)=d(i)-e(i)*b                                                 
17671 continue                                                          
      continue                                                          
      do 17701 k=2,nk                                                   
      j1=kp(k-1)+1                                                      
      j2=kp(k)                                                          
      b=b+dk(k)/u(k)                                                    
      c=c+dk(k)/u(k)**2                                                 
      do 17711 j=j1,j2                                                  
      i=jp(j)                                                           
      w(i)=e(i)*(b-e(i)*c)                                              
      if(w(i) .gt. 0.0)goto 17731                                       
      jerr=-30000                                                       
      return                                                            
17731 continue                                                          
      wr(i)=d(i)-e(i)*b                                                 
17711 continue                                                          
      continue                                                          
17701 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine vars(no,ni,x,w,ixx,v)                                  
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),w(no),v(ni)                             
      integer ixx(ni)                                                   
      do 17741 j=1,ni                                                   
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                     
17741 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine died(no,nk,d,kp,jp,dk)                                 
      implicit double precision(a-h,o-z)                                
      double precision d(no),dk(nk)                                     
      integer kp(nk),jp(no)                                             
      dk(1)=sum(d(jp(1:kp(1))))                                         
      do 17751 k=2,nk                                                   
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                               
17751 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine usk(no,nk,kp,jp,e,u)                                   
      implicit double precision(a-h,o-z)                                
      double precision e(no),u(nk),h                                    
      integer kp(nk),jp(no)                                             
      h=0.0                                                             
      do 17761 k=nk,1,-1                                                
      j2=kp(k)                                                          
      j1=1                                                              
      if(k.gt.1) j1=kp(k-1)+1                                           
      do 17771 j=j2,j1,-1                                               
      h=h+e(jp(j))                                                      
17771 continue                                                          
      continue                                                          
      u(k)=h                                                            
17761 continue                                                          
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
      if(sw .gt. 0.0)goto 17791                                         
      jerr=9999                                                         
      go to 12220                                                       
17791 continue                                                          
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                            
      if(jerr.ne.0) go to 12220                                         
      fmax=log(huge(e(1))*0.1)                                          
      dq=d*q                                                            
      call died(no,nk,dq,kp,jp,dk)                                      
      gm=dot_product(q,g)/sw                                            
      do 17801 j=1,ni                                                   
      xm(j)=dot_product(q,x(:,j))/sw                                    
17801 continue                                                          
      continue                                                          
      do 17811 lam=1,nlam                                               
      do 17821 i=1,no                                                   
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                    
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                     
17821 continue                                                          
      continue                                                          
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                       
17811 continue                                                          
      continue                                                          
12220 continue                                                          
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
      if(maxval(vp) .gt. 0.0)goto 17841                                 
      jerr=10000                                                        
      return                                                            
17841 continue                                                          
      if(minval(y) .ge. 0.0)goto 17861                                  
      jerr=8888                                                         
      return                                                            
17861 continue                                                          
      allocate(ww(1:no),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(ju(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(vq(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      allocate(xm(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      if(isd .le. 0)goto 17881                                          
      allocate(xs(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
17881 continue                                                          
      call chkvars(no,ni,x,ju)                                          
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                              
      if(maxval(ju) .gt. 0)goto 17901                                   
      jerr=7777                                                         
      go to 12220                                                       
17901 continue                                                          
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      ww=max(0d0,w)                                                     
      sw=sum(ww)                                                        
      if(sw .gt. 0.0)goto 17921                                         
      jerr=9999                                                         
      go to 12220                                                       
17921 continue                                                          
      ww=ww/sw                                                          
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                     
      if(isd .le. 0)goto 17941                                          
      do 17951 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
17951 continue                                                          
      continue                                                          
17941 continue                                                          
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,t
     *hr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 12220                                         
      dev0=2.0*sw*dev0                                                  
      do 17961 k=1,lmu                                                  
      nk=nin(k)                                                         
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                   
      if(intr .ne. 0)goto 17981                                         
      a0(k)=0.0                                                         
      goto 17991                                                        
17981 continue                                                          
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                  
17991 continue                                                          
      continue                                                          
17961 continue                                                          
      continue                                                          
12220 continue                                                          
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
      if(nonzero(no,g) .ne. 0)goto 18011                                
      if(intr .eq. 0)goto 18031                                         
      w=q*yb                                                            
      az=log(yb)                                                        
      f=az                                                              
      dv0=yb*(az-1.0)                                                   
      goto 18041                                                        
18031 continue                                                          
      w=q                                                               
      az=0.0                                                            
      f=az                                                              
      dv0=-1.0                                                          
18041 continue                                                          
      continue                                                          
      goto 18051                                                        
18011 continue                                                          
      w=q*exp(sign(min(abs(g),fmax),g))                                 
      v0=sum(w)                                                         
      if(intr .eq. 0)goto 18071                                         
      eaz=yb/v0                                                         
      w=eaz*w                                                           
      az=log(eaz)                                                       
      f=az+g                                                            
      dv0=dot_product(t,g)-yb*(1.0-az)                                  
      goto 18081                                                        
18071 continue                                                          
      az=0.0                                                            
      f=g                                                               
      dv0=dot_product(t,g)-v0                                           
18081 continue                                                          
      continue                                                          
18051 continue                                                          
      continue                                                          
      a=0.0                                                             
      as=0.0                                                            
      wr=t-w                                                            
      v0=1.0                                                            
      if(intr.ne.0) v0=yb                                               
      dvr=-yb                                                           
      do 18091 i=1,no                                                   
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                            
18091 continue                                                          
      continue                                                          
      dvr=dvr-dv0                                                       
      dev0=dvr                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 18111                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
18111 continue                                                          
      m=0                                                               
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      mnl=min(mnlam,nlam)                                               
      shr=shri*dev0                                                     
      ixx=0                                                             
      al=0.0                                                            
      do 18121 j=1,ni                                                   
      if(ju(j).eq.0)goto 18121                                          
      ga(j)=abs(dot_product(wr,x(:,j)))                                 
18121 continue                                                          
      continue                                                          
      do 18131 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 18151                                      
      al=ulam(ilm)                                                      
      goto 18141                                                        
18151 if(ilm .le. 2)goto 18161                                          
      al=al*alf                                                         
      goto 18141                                                        
18161 if(ilm .ne. 1)goto 18171                                          
      al=big                                                            
      goto 18181                                                        
18171 continue                                                          
      al0=0.0                                                           
      do 18191 j=1,ni                                                   
      if(ju(j).eq.0)goto 18191                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
18191 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
18181 continue                                                          
18141 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 18201 k=1,ni                                                   
      if(ixx(k).eq.1)goto 18201                                         
      if(ju(k).eq.0)goto 18201                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
18201 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
18211 continue                                                          
      if(nlp .le. maxit)goto 18231                                      
      jerr=-ilm                                                         
      return                                                            
18231 continue                                                          
      az0=az                                                            
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                             
      do 18241 j=1,ni                                                   
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                     
18241 continue                                                          
      continue                                                          
      continue                                                          
18251 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 18261 k=1,ni                                                   
      if(ixx(k).eq.0)goto 18261                                         
      ak=a(k)                                                           
      u=dot_product(wr,x(:,k))+v(k)*ak                                  
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 18281                                         
      a(k)=0.0                                                          
      goto 18291                                                        
18281 continue                                                          
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))        
18291 continue                                                          
      continue                                                          
      if(a(k).eq.ak)goto 18261                                          
      d=a(k)-ak                                                         
      dlx=max(dlx,v(k)*d**2)                                            
      wr=wr-d*w*x(:,k)                                                  
      f=f+d*x(:,k)                                                      
      if(mm(k) .ne. 0)goto 18311                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 18262                                           
      mm(k)=nin                                                         
      m(nin)=k                                                          
18311 continue                                                          
18261 continue                                                          
18262 continue                                                          
      if(nin.gt.nx)goto 18252                                           
      if(intr .eq. 0)goto 18331                                         
      d=sum(wr)/v0                                                      
      az=az+d                                                           
      dlx=max(dlx,v0*d**2)                                              
      wr=wr-d*w                                                         
      f=f+d                                                             
18331 continue                                                          
      if(dlx.lt.shr)goto 18252                                          
      if(nlp .le. maxit)goto 18351                                      
      jerr=-ilm                                                         
      return                                                            
18351 continue                                                          
      continue                                                          
18361 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 18371 l=1,nin                                                  
      k=m(l)                                                            
      ak=a(k)                                                           
      u=dot_product(wr,x(:,k))+v(k)*ak                                  
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 18391                                         
      a(k)=0.0                                                          
      goto 18401                                                        
18391 continue                                                          
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))        
18401 continue                                                          
      continue                                                          
      if(a(k).eq.ak)goto 18371                                          
      d=a(k)-ak                                                         
      dlx=max(dlx,v(k)*d**2)                                            
      wr=wr-d*w*x(:,k)                                                  
      f=f+d*x(:,k)                                                      
18371 continue                                                          
      continue                                                          
      if(intr .eq. 0)goto 18421                                         
      d=sum(wr)/v0                                                      
      az=az+d                                                           
      dlx=max(dlx,v0*d**2)                                              
      wr=wr-d*w                                                         
      f=f+d                                                             
18421 continue                                                          
      if(dlx.lt.shr)goto 18362                                          
      if(nlp .le. maxit)goto 18441                                      
      jerr=-ilm                                                         
      return                                                            
18441 continue                                                          
      goto 18361                                                        
18362 continue                                                          
      goto 18251                                                        
18252 continue                                                          
      if(nin.gt.nx)goto 18212                                           
      w=q*exp(sign(min(abs(f),fmax),f))                                 
      v0=sum(w)                                                         
      wr=t-w                                                            
      if(v0*(az-az0)**2 .ge. shr)goto 18461                             
      ix=0                                                              
      do 18471 j=1,nin                                                  
      k=m(j)                                                            
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 18471                         
      ix=1                                                              
      goto 18472                                                        
18471 continue                                                          
18472 continue                                                          
      if(ix .ne. 0)goto 18491                                           
      do 18501 k=1,ni                                                   
      if(ixx(k).eq.1)goto 18501                                         
      if(ju(k).eq.0)goto 18501                                          
      ga(k)=abs(dot_product(wr,x(:,k)))                                 
      if(ga(k) .le. al1*vp(k))goto 18521                                
      ixx(k)=1                                                          
      ix=1                                                              
18521 continue                                                          
18501 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10880                                           
      goto 18212                                                        
18491 continue                                                          
18461 continue                                                          
      goto 18211                                                        
18212 continue                                                          
      if(nin .le. nx)goto 18541                                         
      jerr=-10000-ilm                                                   
      goto 18132                                                        
18541 continue                                                          
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                            
      kin(ilm)=nin                                                      
      a0(ilm)=az                                                        
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                            
      if(ilm.lt.mnl)goto 18131                                          
      if(flmin.ge.1.0)goto 18131                                        
      me=0                                                              
      do 18551 j=1,nin                                                  
      if(ca(j,ilm).ne.0.0) me=me+1                                      
18551 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 18132                                            
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 18132           
      if(dev(ilm).gt.devmax)goto 18132                                  
18131 continue                                                          
18132 continue                                                          
      g=f                                                               
      continue                                                          
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                             
      return                                                            
      end                                                               
      function nonzero(n,v)                                             
      implicit double precision(a-h,o-z)                                
      double precision v(n)                                             
      nonzero=0                                                         
      do 18561 i=1,n                                                    
      if(v(i) .eq. 0.0)goto 18581                                       
      nonzero=1                                                         
      return                                                            
18581 continue                                                          
18561 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                            
      implicit double precision(a-h,o-z)                                
      double precision a(nx,lmu),b(ni,lmu)                              
      integer ia(nx),nin(lmu)                                           
      do 18591 lam=1,lmu                                                
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                     
18591 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                        
      implicit double precision(a-h,o-z)                                
      double precision a(nx,nc,lmu),b(ni,nc,lmu)                        
      integer ia(nx),nin(lmu)                                           
      do 18601 lam=1,lmu                                                
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))          
18601 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)            
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),fl
     *og(nlam)
      double precision, dimension (:), allocatable :: w                 
      if(minval(y) .ge. 0.0)goto 18621                                  
      jerr=8888                                                         
      return                                                            
18621 continue                                                          
      allocate(w(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      w=max(0d0,q)                                                      
      sw=sum(w)                                                         
      if(sw .gt. 0.0)goto 18641                                         
      jerr=9999                                                         
      go to 12220                                                       
18641 continue                                                          
      yb=dot_product(w,y)/sw                                            
      fmax=log(huge(y(1))*0.1)                                          
      do 18651 lam=1,nlam                                               
      s=0.0                                                             
      do 18661 i=1,no                                                   
      if(w(i).le.0.0)goto 18661                                         
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                       
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                   
18661 continue                                                          
      continue                                                          
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                             
18651 continue                                                          
      continue                                                          
12220 continue                                                          
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
      if(maxval(vp) .gt. 0.0)goto 18681                                 
      jerr=10000                                                        
      return                                                            
18681 continue                                                          
      if(minval(y) .ge. 0.0)goto 18701                                  
      jerr=8888                                                         
      return                                                            
18701 continue                                                          
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
      if(maxval(ju) .gt. 0)goto 18721                                   
      jerr=7777                                                         
      go to 12220                                                       
18721 continue                                                          
      vq=max(0d0,vp)                                                    
      vq=vq*ni/sum(vq)                                                  
      ww=max(0d0,w)                                                     
      sw=sum(ww)                                                        
      if(sw .gt. 0.0)goto 18741                                         
      jerr=9999                                                         
      go to 12220                                                       
18741 continue                                                          
      ww=ww/sw                                                          
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)             
      if(isd .le. 0)goto 18761                                          
      do 18771 j=1,ni                                                   
      cl(:,j)=cl(:,j)*xs(j)                                             
18771 continue                                                          
      continue                                                          
18761 continue                                                          
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,flmi
     *n,ulam,thr,  isd,intr,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nl
     *p,jerr)
      if(jerr.gt.0) go to 12220                                         
      dev0=2.0*sw*dev0                                                  
      do 18781 k=1,lmu                                                  
      nk=nin(k)                                                         
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                   
      if(intr .ne. 0)goto 18801                                         
      a0(k)=0.0                                                         
      goto 18811                                                        
18801 continue                                                          
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                  
18811 continue                                                          
      continue                                                          
18781 continue                                                          
      continue                                                          
12220 continue                                                          
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
      if(nonzero(no,g) .ne. 0)goto 18831                                
      t=0.0                                                             
      if(intr .eq. 0)goto 18851                                         
      w=q*yb                                                            
      az=log(yb)                                                        
      uu=az                                                             
      xm=yb*xb                                                          
      dv0=yb*(az-1.0)                                                   
      goto 18861                                                        
18851 continue                                                          
      w=q                                                               
      xm=0.0                                                            
      uu=0.0                                                            
      az=uu                                                             
      dv0=-1.0                                                          
18861 continue                                                          
      continue                                                          
      goto 18871                                                        
18831 continue                                                          
      w=q*exp(sign(min(abs(g),fmax),g))                                 
      ww=sum(w)                                                         
      t=g                                                               
      if(intr .eq. 0)goto 18891                                         
      eaz=yb/ww                                                         
      w=eaz*w                                                           
      az=log(eaz)                                                       
      uu=az                                                             
      dv0=dot_product(qy,g)-yb*(1.0-az)                                 
      goto 18901                                                        
18891 continue                                                          
      uu=0.0                                                            
      az=uu                                                             
      dv0=dot_product(qy,g)-ww                                          
18901 continue                                                          
      continue                                                          
      do 18911 j=1,ni                                                   
      if(ju(j).eq.0)goto 18911                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
18911 continue                                                          
      continue                                                          
18871 continue                                                          
      continue                                                          
      tt=yb*uu                                                          
      ww=1.0                                                            
      if(intr.ne.0) ww=yb                                               
      wr=qy-q*(yb*(1.0-uu))                                             
      a=0.0                                                             
      as=0.0                                                            
      dvr=-yb                                                           
      do 18921 i=1,no                                                   
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                          
18921 continue                                                          
      continue                                                          
      dvr=dvr-dv0                                                       
      dev0=dvr                                                          
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 18941                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
18941 continue                                                          
      m=0                                                               
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      mnl=min(mnlam,nlam)                                               
      shr=shri*dev0                                                     
      al=0.0                                                            
      ixx=0                                                             
      do 18951 j=1,ni                                                   
      if(ju(j).eq.0)goto 18951                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)
     *)-xb(j)*tt)/xs(j)
18951 continue                                                          
      continue                                                          
      do 18961 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 18981                                      
      al=ulam(ilm)                                                      
      goto 18971                                                        
18981 if(ilm .le. 2)goto 18991                                          
      al=al*alf                                                         
      goto 18971                                                        
18991 if(ilm .ne. 1)goto 19001                                          
      al=big                                                            
      goto 19011                                                        
19001 continue                                                          
      al0=0.0                                                           
      do 19021 j=1,ni                                                   
      if(ju(j).eq.0)goto 19021                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
19021 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
19011 continue                                                          
18971 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 19031 k=1,ni                                                   
      if(ixx(k).eq.1)goto 19031                                         
      if(ju(k).eq.0)goto 19031                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
19031 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
19041 continue                                                          
      if(nlp .le. maxit)goto 19061                                      
      jerr=-ilm                                                         
      return                                                            
19061 continue                                                          
      az0=az                                                            
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                             
      do 19071 j=1,ni                                                   
      if(ixx(j).eq.0)goto 19071                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x
     *b(j)**2)/xs(j)**2
19071 continue                                                          
      continue                                                          
      continue                                                          
19081 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 19091 k=1,ni                                                   
      if(ixx(k).eq.0)goto 19091                                         
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      ak=a(k)                                                           
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 19111                                         
      a(k)=0.0                                                          
      goto 19121                                                        
19111 continue                                                          
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))        
19121 continue                                                          
      continue                                                          
      if(a(k).eq.ak)goto 19091                                          
      if(mm(k) .ne. 0)goto 19141                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 19092                                           
      mm(k)=nin                                                         
      m(nin)=k                                                          
19141 continue                                                          
      d=a(k)-ak                                                         
      dlx=max(dlx,v(k)*d**2)                                            
      dv=d/xs(k)                                                        
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)              
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                             
      uu=uu-dv*xb(k)                                                    
      tt=tt-dv*xm(k)                                                    
19091 continue                                                          
19092 continue                                                          
      if(nin.gt.nx)goto 19082                                           
      if(intr .eq. 0)goto 19161                                         
      d=tt/ww-uu                                                        
      az=az+d                                                           
      dlx=max(dlx,ww*d**2)                                              
      uu=uu+d                                                           
19161 continue                                                          
      if(dlx.lt.shr)goto 19082                                          
      if(nlp .le. maxit)goto 19181                                      
      jerr=-ilm                                                         
      return                                                            
19181 continue                                                          
      continue                                                          
19191 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 19201 l=1,nin                                                  
      k=m(l)                                                            
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      ak=a(k)                                                           
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                               
      if(au .gt. 0.0)goto 19221                                         
      a(k)=0.0                                                          
      goto 19231                                                        
19221 continue                                                          
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))        
19231 continue                                                          
      continue                                                          
      if(a(k).eq.ak)goto 19201                                          
      d=a(k)-ak                                                         
      dlx=max(dlx,v(k)*d**2)                                            
      dv=d/xs(k)                                                        
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)              
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                             
      uu=uu-dv*xb(k)                                                    
      tt=tt-dv*xm(k)                                                    
19201 continue                                                          
      continue                                                          
      if(intr .eq. 0)goto 19251                                         
      d=tt/ww-uu                                                        
      az=az+d                                                           
      dlx=max(dlx,ww*d**2)                                              
      uu=uu+d                                                           
19251 continue                                                          
      if(dlx.lt.shr)goto 19192                                          
      if(nlp .le. maxit)goto 19271                                      
      jerr=-ilm                                                         
      return                                                            
19271 continue                                                          
      goto 19191                                                        
19192 continue                                                          
      goto 19081                                                        
19082 continue                                                          
      if(nin.gt.nx)goto 19042                                           
      euu=exp(sign(min(abs(uu),fmax),uu))                               
      w=euu*q*exp(sign(min(abs(t),fmax),t))                             
      ww=sum(w)                                                         
      wr=qy-w*(1.0-uu)                                                  
      tt=sum(wr)                                                        
      if(ww*(az-az0)**2 .ge. shr)goto 19291                             
      kx=0                                                              
      do 19301 j=1,nin                                                  
      k=m(j)                                                            
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 19301                         
      kx=1                                                              
      goto 19302                                                        
19301 continue                                                          
19302 continue                                                          
      if(kx .ne. 0)goto 19321                                           
      do 19331 j=1,ni                                                   
      if(ixx(j).eq.1)goto 19331                                         
      if(ju(j).eq.0)goto 19331                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)
     *)-xb(j)*tt)/xs(j)
      if(ga(j) .le. al1*vp(j))goto 19351                                
      ixx(j)=1                                                          
      kx=1                                                              
19351 continue                                                          
19331 continue                                                          
      continue                                                          
      if(kx.eq.1) go to 10880                                           
      goto 19042                                                        
19321 continue                                                          
19291 continue                                                          
      goto 19041                                                        
19042 continue                                                          
      if(nin .le. nx)goto 19371                                         
      jerr=-10000-ilm                                                   
      goto 18962                                                        
19371 continue                                                          
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                            
      kin(ilm)=nin                                                      
      a0(ilm)=az                                                        
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                     
      if(ilm.lt.mnl)goto 18961                                          
      if(flmin.ge.1.0)goto 18961                                        
      me=0                                                              
      do 19381 j=1,nin                                                  
      if(ca(j,ilm).ne.0.0) me=me+1                                      
19381 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 18962                                            
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 18962           
      if(dev(ilm).gt.devmax)goto 18962                                  
18961 continue                                                          
18962 continue                                                          
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
      if(minval(y) .ge. 0.0)goto 19401                                  
      jerr=8888                                                         
      return                                                            
19401 continue                                                          
      allocate(w(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(f(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      w=max(0d0,q)                                                      
      sw=sum(w)                                                         
      if(sw .gt. 0.0)goto 19421                                         
      jerr=9999                                                         
      go to 12220                                                       
19421 continue                                                          
      yb=dot_product(w,y)/sw                                            
      fmax=log(huge(y(1))*0.1)                                          
      do 19431 lam=1,nlam                                               
      f=a0(lam)                                                         
      do 19441 j=1,ni                                                   
      if(a(j,lam).eq.0.0)goto 19441                                     
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                       
19441 continue                                                          
      continue                                                          
      f=f+g                                                             
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                             
19431 continue                                                          
      continue                                                          
12220 continue                                                          
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
      if(minval(y) .ge. 0.0)goto 19461                                  
      jerr=8888                                                         
      return                                                            
19461 continue                                                          
      allocate(w(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      allocate(f(1:no),stat=jerr)                                       
      if(jerr.ne.0) return                                              
      w=max(0d0,q)                                                      
      sw=sum(w)                                                         
      if(sw .gt. 0.0)goto 19481                                         
      jerr=9999                                                         
      go to 12220                                                       
19481 continue                                                          
      yb=dot_product(w,y)/sw                                            
      fmax=log(huge(y(1))*0.1)                                          
      do 19491 lam=1,nlam                                               
      f=a0(lam)                                                         
      do 19501 k=1,nin(lam)                                             
      j=ia(k)                                                           
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                      
19501 continue                                                          
      continue                                                          
      f=f+g                                                             
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                             
19491 continue                                                          
      continue                                                          
12220 continue                                                          
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
      if(maxval(vp) .gt. 0.0)goto 19521                                 
      jerr=10000                                                        
      return                                                            
19521 continue                                                          
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
      if(maxval(ju) .gt. 0)goto 19541                                   
      jerr=7777                                                         
      return                                                            
19541 continue                                                          
      call multstandard1(no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym,ys,xv,y
     *s0,jerr)
      if(jerr.ne.0) return                                              
      do 19551 j=1,ni                                                   
      do 19561 k=1,nr                                                   
      do 19571 i=1,2                                                    
      clt(i,k,j)=cl(i,j)                                                
19571 continue                                                          
      continue                                                          
19561 continue                                                          
      continue                                                          
19551 continue                                                          
      continue                                                          
      if(isd .le. 0)goto 19591                                          
      do 19601 j=1,ni                                                   
      do 19611 k=1,nr                                                   
      do 19621 i=1,2                                                    
      clt(i,k,j)=clt(i,k,j)*xs(j)                                       
19621 continue                                                          
      continue                                                          
19611 continue                                                          
      continue                                                          
19601 continue                                                          
      continue                                                          
19591 continue                                                          
      if(jsd .le. 0)goto 19641                                          
      do 19651 j=1,ni                                                   
      do 19661 k=1,nr                                                   
      do 19671 i=1,2                                                    
      clt(i,k,j)=clt(i,k,j)/ys(k)                                       
19671 continue                                                          
      continue                                                          
19661 continue                                                          
      continue                                                          
19651 continue                                                          
      continue                                                          
19641 continue                                                          
      call multelnet2(parm,ni,nr,ju,vp,clt,y,no,ne,nx,x,nlam,flmin,ulam,
     *thr,maxit,xv,  ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                              
      do 19681 k=1,lmu                                                  
      nk=nin(k)                                                         
      do 19691 j=1,nr                                                   
      do 19701 l=1,nk                                                   
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                               
19701 continue                                                          
      continue                                                          
      if(intr .ne. 0)goto 19721                                         
      a0(j,k)=0.0                                                       
      goto 19731                                                        
19721 continue                                                          
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))              
19731 continue                                                          
      continue                                                          
19691 continue                                                          
      continue                                                          
19681 continue                                                          
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
      if(intr .ne. 0)goto 19751                                         
      do 19761 j=1,ni                                                   
      if(ju(j).eq.0)goto 19761                                          
      xm(j)=0.0                                                         
      x(:,j)=v*x(:,j)                                                   
      z=dot_product(x(:,j),x(:,j))                                      
      if(isd .le. 0)goto 19781                                          
      xbq=dot_product(v,x(:,j))**2                                      
      vc=z-xbq                                                          
      xs(j)=sqrt(vc)                                                    
      x(:,j)=x(:,j)/xs(j)                                               
      xv(j)=1.0+xbq/vc                                                  
      goto 19791                                                        
19781 continue                                                          
      xs(j)=1.0                                                         
      xv(j)=z                                                           
19791 continue                                                          
      continue                                                          
19761 continue                                                          
      continue                                                          
      ys0=0.0                                                           
      do 19801 j=1,nr                                                   
      ym(j)=0.0                                                         
      y(:,j)=v*y(:,j)                                                   
      z=dot_product(y(:,j),y(:,j))                                      
      if(jsd .le. 0)goto 19821                                          
      u=z-dot_product(v,y(:,j))**2                                      
      ys0=ys0+z/u                                                       
      ys(j)=sqrt(u)                                                     
      y(:,j)=y(:,j)/ys(j)                                               
      goto 19831                                                        
19821 continue                                                          
      ys(j)=1.0                                                         
      ys0=ys0+z                                                         
19831 continue                                                          
      continue                                                          
19801 continue                                                          
      continue                                                          
      go to 10700                                                       
19751 continue                                                          
      do 19841 j=1,ni                                                   
      if(ju(j).eq.0)goto 19841                                          
      xm(j)=dot_product(w,x(:,j))                                       
      x(:,j)=v*(x(:,j)-xm(j))                                           
      xv(j)=dot_product(x(:,j),x(:,j))                                  
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                    
19841 continue                                                          
      continue                                                          
      if(isd .ne. 0)goto 19861                                          
      xs=1.0                                                            
      goto 19871                                                        
19861 continue                                                          
      do 19881 j=1,ni                                                   
      if(ju(j).eq.0)goto 19881                                          
      x(:,j)=x(:,j)/xs(j)                                               
19881 continue                                                          
      continue                                                          
      xv=1.0                                                            
19871 continue                                                          
      continue                                                          
      ys0=0.0                                                           
      do 19891 j=1,nr                                                   
      ym(j)=dot_product(w,y(:,j))                                       
      y(:,j)=v*(y(:,j)-ym(j))                                           
      z=dot_product(y(:,j),y(:,j))                                      
      if(jsd .le. 0)goto 19911                                          
      ys(j)=sqrt(z)                                                     
      y(:,j)=y(:,j)/ys(j)                                               
      goto 19921                                                        
19911 continue                                                          
      ys0=ys0+z                                                         
19921 continue                                                          
      continue                                                          
19891 continue                                                          
      continue                                                          
      if(jsd .ne. 0)goto 19941                                          
      ys=1.0                                                            
      goto 19951                                                        
19941 continue                                                          
      ys0=nr                                                            
19951 continue                                                          
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
      if(flmin .ge. 1.0)goto 19971                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
19971 continue                                                          
      rsq=ys0                                                           
      a=0.0                                                             
      mm=0                                                              
      nlp=0                                                             
      nin=nlp                                                           
      iz=0                                                              
      mnl=min(mnlam,nlam)                                               
      alm=0.0                                                           
      do 19981 j=1,ni                                                   
      if(ju(j).eq.0)goto 19981                                          
      g(j)=0.0                                                          
      do 19991 k=1,nr                                                   
      g(j)=g(j)+dot_product(y(:,k),x(:,j))**2                           
19991 continue                                                          
      continue                                                          
      g(j)=sqrt(g(j))                                                   
19981 continue                                                          
      continue                                                          
      do 20001 m=1,nlam                                                 
      if(itrace.ne.0) call setpb(m-1)                                   
      alm0=alm                                                          
      if(flmin .lt. 1.0)goto 20021                                      
      alm=ulam(m)                                                       
      goto 20011                                                        
20021 if(m .le. 2)goto 20031                                            
      alm=alm*alf                                                       
      goto 20011                                                        
20031 if(m .ne. 1)goto 20041                                            
      alm=big                                                           
      goto 20051                                                        
20041 continue                                                          
      alm0=0.0                                                          
      do 20061 j=1,ni                                                   
      if(ju(j).eq.0)goto 20061                                          
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                        
20061 continue                                                          
      continue                                                          
      alm0=alm0/max(bta,1.0d-3)                                         
      alm=alf*alm0                                                      
20051 continue                                                          
20011 continue                                                          
      dem=alm*omb                                                       
      ab=alm*bta                                                        
      rsq0=rsq                                                          
      jz=1                                                              
      tlam=bta*(2.0*alm-alm0)                                           
      do 20071 k=1,ni                                                   
      if(ix(k).eq.1)goto 20071                                          
      if(ju(k).eq.0)goto 20071                                          
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                    
20071 continue                                                          
      continue                                                          
      continue                                                          
20081 continue                                                          
      if(iz*jz.ne.0) go to 10360                                        
10880 if(nlp .le. maxit)goto 20101                                      
      jerr=-m                                                           
      return                                                            
20101 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 20111 k=1,ni                                                   
      if(ix(k).eq.0)goto 20111                                          
      gkn=0.0                                                           
      do 20121 j=1,nr                                                   
      gj(j)=dot_product(y(:,j),x(:,k))                                  
      gk(j)=gj(j)+a(j,k)*xv(k)                                          
      gkn=gkn+gk(j)**2                                                  
20121 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-ab*vp(k)/gkn                                                
      del=a(:,k)                                                        
      if(u .gt. 0.0)goto 20141                                          
      a(:,k)=0.0                                                        
      goto 20151                                                        
20141 continue                                                          
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                   
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)
     *,isc,jerr)
      if(jerr.ne.0) return                                              
20151 continue                                                          
      continue                                                          
      del=a(:,k)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 20111                             
      do 20161 j=1,nr                                                   
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                           
      y(:,j)=y(:,j)-del(j)*x(:,k)                                       
      dlx=max(dlx,xv(k)*del(j)**2)                                      
20161 continue                                                          
      continue                                                          
      if(mm(k) .ne. 0)goto 20181                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 20112                                           
      mm(k)=nin                                                         
      ia(nin)=k                                                         
20181 continue                                                          
20111 continue                                                          
20112 continue                                                          
      if(nin.gt.nx)goto 20082                                           
      if(dlx .ge. thr)goto 20201                                        
      ixx=0                                                             
      do 20211 k=1,ni                                                   
      if(ix(k).eq.1)goto 20211                                          
      if(ju(k).eq.0)goto 20211                                          
      g(k)=0.0                                                          
      do 20221 j=1,nr                                                   
      g(k)=g(k)+dot_product(y(:,j),x(:,k))**2                           
20221 continue                                                          
      continue                                                          
      g(k)=sqrt(g(k))                                                   
      if(g(k) .le. ab*vp(k))goto 20241                                  
      ix(k)=1                                                           
      ixx=1                                                             
20241 continue                                                          
20211 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10880                                          
      goto 20082                                                        
20201 continue                                                          
      if(nlp .le. maxit)goto 20261                                      
      jerr=-m                                                           
      return                                                            
20261 continue                                                          
10360 continue                                                          
      iz=1                                                              
      continue                                                          
20271 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 20281 l=1,nin                                                  
      k=ia(l)                                                           
      gkn=0.0                                                           
      do 20291 j=1,nr                                                   
      gj(j)=dot_product(y(:,j),x(:,k))                                  
      gk(j)=gj(j)+a(j,k)*xv(k)                                          
      gkn=gkn+gk(j)**2                                                  
20291 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-ab*vp(k)/gkn                                                
      del=a(:,k)                                                        
      if(u .gt. 0.0)goto 20311                                          
      a(:,k)=0.0                                                        
      goto 20321                                                        
20311 continue                                                          
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                   
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)
     *,isc,jerr)
      if(jerr.ne.0) return                                              
20321 continue                                                          
      continue                                                          
      del=a(:,k)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 20281                             
      do 20331 j=1,nr                                                   
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                           
      y(:,j)=y(:,j)-del(j)*x(:,k)                                       
      dlx=max(dlx,xv(k)*del(j)**2)                                      
20331 continue                                                          
      continue                                                          
20281 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 20272                                          
      if(nlp .le. maxit)goto 20351                                      
      jerr=-m                                                           
      return                                                            
20351 continue                                                          
      goto 20271                                                        
20272 continue                                                          
      jz=0                                                              
      goto 20081                                                        
20082 continue                                                          
      if(nin .le. nx)goto 20371                                         
      jerr=-10000-m                                                     
      goto 20002                                                        
20371 continue                                                          
      if(nin .le. 0)goto 20391                                          
      do 20401 j=1,nr                                                   
      ao(1:nin,j,m)=a(j,ia(1:nin))                                      
20401 continue                                                          
      continue                                                          
20391 continue                                                          
      kin(m)=nin                                                        
      rsqo(m)=1.0-rsq/ys0                                               
      almo(m)=alm                                                       
      lmu=m                                                             
      if(m.lt.mnl)goto 20001                                            
      if(flmin.ge.1.0)goto 20001                                        
      me=0                                                              
      do 20411 j=1,nin                                                  
      if(ao(j,1,m).ne.0.0) me=me+1                                      
20411 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 20002                                            
      if(rsq0-rsq.lt.sml*rsq)goto 20002                                 
      if(rsqo(m).gt.rsqmax)goto 20002                                   
20001 continue                                                          
20002 continue                                                          
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
20421 continue                                                          
      vmx=0.0                                                           
      do 20431 k=1,nr                                                   
      v=max(a(k)-cl(2,k),cl(1,k)-a(k))                                  
      if(v .le. vmx)goto 20451                                          
      vmx=v                                                             
      kn=k                                                              
20451 continue                                                          
20431 continue                                                          
      continue                                                          
      if(vmx.le.0.0)goto 20422                                          
      if(isc(kn).ne.0)goto 20422                                        
      gsq=gsq-gk(kn)**2                                                 
      g=sqrt(gsq)/xv                                                    
      if(a(kn).lt.cl(1,kn)) u=cl(1,kn)                                  
      if(a(kn).gt.cl(2,kn)) u=cl(2,kn)                                  
      usq=usq+u**2                                                      
      if(usq .ne. 0.0)goto 20471                                        
      b=max(0d0,(g-al2p)/al1p)                                          
      goto 20481                                                        
20471 continue                                                          
      b0=sqrt(asq-a(kn)**2)                                             
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                  
      if(kerr.ne.0)goto 20422                                           
20481 continue                                                          
      continue                                                          
      asq=usq+b**2                                                      
      if(asq .gt. 0.0)goto 20501                                        
      a=0.0                                                             
      goto 20422                                                        
20501 continue                                                          
      a(kn)=u                                                           
      isc(kn)=1                                                         
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                  
      do 20511 j=1,nr                                                   
      if(isc(j).eq.0) a(j)=f*gk(j)                                      
20511 continue                                                          
      continue                                                          
      goto 20421                                                        
20422 continue                                                          
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
20521 continue                                                          
      vmx=0.0                                                           
      do 20531 k=1,nr                                                   
      v=max(a(k)-cl2,cl1-a(k))                                          
      if(v .le. vmx)goto 20551                                          
      vmx=v                                                             
      kn=k                                                              
20551 continue                                                          
20531 continue                                                          
      continue                                                          
      if(vmx.le.0.0)goto 20522                                          
      if(isc(kn).ne.0)goto 20522                                        
      gsq=gsq-gk(kn)**2                                                 
      g=sqrt(gsq)/xv                                                    
      if(a(kn).lt.cl1) u=cl1                                            
      if(a(kn).gt.cl2) u=cl2                                            
      usq=usq+u**2                                                      
      if(usq .ne. 0.0)goto 20571                                        
      b=max(0d0,(g-al2p)/al1p)                                          
      goto 20581                                                        
20571 continue                                                          
      b0=sqrt(asq-a(kn)**2)                                             
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                  
      if(kerr.ne.0)goto 20522                                           
20581 continue                                                          
      continue                                                          
      asq=usq+b**2                                                      
      if(asq .gt. 0.0)goto 20601                                        
      a=0.0                                                             
      goto 20522                                                        
20601 continue                                                          
      a(kn)=u                                                           
      isc(kn)=1                                                         
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                  
      do 20611 j=1,nr                                                   
      if(isc(j).eq.0) a(j)=f*gk(j)                                      
20611 continue                                                          
      continue                                                          
      goto 20521                                                        
20522 continue                                                          
      if(kerr.ne.0) jerr=kerr                                           
      return                                                            
      end                                                               
      function bnorm(b0,al1p,al2p,g,usq,jerr)                           
      implicit double precision(a-h,o-z)                                
      data thr,mxit /1.0d-10,100/                                       
      b=b0                                                              
      zsq=b**2+usq                                                      
      if(zsq .gt. 0.0)goto 20631                                        
      bnorm=0.0                                                         
      return                                                            
20631 continue                                                          
      z=sqrt(zsq)                                                       
      f=b*(al1p+al2p/z)-g                                               
      jerr=0                                                            
      do 20641 it=1,mxit                                                
      b=b-f/(al1p+al2p*usq/(z*zsq))                                     
      zsq=b**2+usq                                                      
      if(zsq .gt. 0.0)goto 20661                                        
      bnorm=0.0                                                         
      return                                                            
20661 continue                                                          
      z=sqrt(zsq)                                                       
      f=b*(al1p+al2p/z)-g                                               
      if(abs(f).le.thr)goto 20642                                       
      if(b .gt. 0.0)goto 20681                                          
      b=0.0                                                             
      goto 20642                                                        
20681 continue                                                          
20641 continue                                                          
20642 continue                                                          
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
      do 20691 lam=1,lmu                                                
      call multuncomp(ni,nr,nx,a(1,1,lam),ia,nin(lam),b(1,1,lam))       
20691 continue                                                          
      continue                                                          
      return                                                            
      end                                                               
      subroutine multuncomp(ni,nr,nx,ca,ia,nin,a)                       
      implicit double precision(a-h,o-z)                                
      double precision ca(nx,nr),a(ni,nr)                               
      integer ia(nx)                                                    
      a=0.0                                                             
      if(nin .le. 0)goto 20711                                          
      do 20721 j=1,nr                                                   
      a(ia(1:nin),j)=ca(1:nin,j)                                        
20721 continue                                                          
      continue                                                          
20711 continue                                                          
      return                                                            
      end                                                               
      subroutine multmodval(nx,nr,a0,ca,ia,nin,n,x,f)                   
      implicit double precision(a-h,o-z)                                
      double precision a0(nr),ca(nx,nr),x(n,*),f(nr,n)                  
      integer ia(nx)                                                    
      do 20731 i=1,n                                                    
      f(:,i)=a0                                                         
20731 continue                                                          
      continue                                                          
      if(nin.le.0) return                                               
      do 20741 i=1,n                                                    
      do 20751 j=1,nr                                                   
      f(j,i)=f(j,i)+dot_product(ca(1:nin,j),x(i,ia(1:nin)))             
20751 continue                                                          
      continue                                                          
20741 continue                                                          
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
      if(maxval(vp) .gt. 0.0)goto 20771                                 
      jerr=10000                                                        
      return                                                            
20771 continue                                                          
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
      if(maxval(ju) .gt. 0)goto 20791                                   
      jerr=7777                                                         
      return                                                            
20791 continue                                                          
      call multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,  xm,xs,
     *ym,ys,xv,ys0,jerr)
      if(jerr.ne.0) return                                              
      do 20801 j=1,ni                                                   
      do 20811 k=1,nr                                                   
      do 20821 i=1,2                                                    
      clt(i,k,j)=cl(i,j)                                                
20821 continue                                                          
      continue                                                          
20811 continue                                                          
      continue                                                          
20801 continue                                                          
      continue                                                          
      if(isd .le. 0)goto 20841                                          
      do 20851 j=1,ni                                                   
      do 20861 k=1,nr                                                   
      do 20871 i=1,2                                                    
      clt(i,k,j)=clt(i,k,j)*xs(j)                                       
20871 continue                                                          
      continue                                                          
20861 continue                                                          
      continue                                                          
20851 continue                                                          
      continue                                                          
20841 continue                                                          
      if(jsd .le. 0)goto 20891                                          
      do 20901 j=1,ni                                                   
      do 20911 k=1,nr                                                   
      do 20921 i=1,2                                                    
      clt(i,k,j)=clt(i,k,j)/ys(k)                                       
20921 continue                                                          
      continue                                                          
20911 continue                                                          
      continue                                                          
20901 continue                                                          
      continue                                                          
20891 continue                                                          
      call multspelnet2(parm,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,clt,nlam,f
     *lmin,  ulam,thr,maxit,xm,xs,xv,ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                              
      do 20931 k=1,lmu                                                  
      nk=nin(k)                                                         
      do 20941 j=1,nr                                                   
      do 20951 l=1,nk                                                   
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                               
20951 continue                                                          
      continue                                                          
      if(intr .ne. 0)goto 20971                                         
      a0(j,k)=0.0                                                       
      goto 20981                                                        
20971 continue                                                          
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))              
20981 continue                                                          
      continue                                                          
20941 continue                                                          
      continue                                                          
20931 continue                                                          
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
      if(intr .ne. 0)goto 21001                                         
      do 21011 j=1,ni                                                   
      if(ju(j).eq.0)goto 21011                                          
      xm(j)=0.0                                                         
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      z=dot_product(w(jx(jb:je)),x(jb:je)**2)                           
      if(isd .le. 0)goto 21031                                          
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                         
      vc=z-xbq                                                          
      xs(j)=sqrt(vc)                                                    
      xv(j)=1.0+xbq/vc                                                  
      goto 21041                                                        
21031 continue                                                          
      xs(j)=1.0                                                         
      xv(j)=z                                                           
21041 continue                                                          
      continue                                                          
21011 continue                                                          
      continue                                                          
      ys0=0.0                                                           
      do 21051 j=1,nr                                                   
      ym(j)=0.0                                                         
      z=dot_product(w,y(:,j)**2)                                        
      if(jsd .le. 0)goto 21071                                          
      u=z-dot_product(w,y(:,j))**2                                      
      ys0=ys0+z/u                                                       
      ys(j)=sqrt(u)                                                     
      y(:,j)=y(:,j)/ys(j)                                               
      goto 21081                                                        
21071 continue                                                          
      ys(j)=1.0                                                         
      ys0=ys0+z                                                         
21081 continue                                                          
      continue                                                          
21051 continue                                                          
      continue                                                          
      return                                                            
21001 continue                                                          
      do 21091 j=1,ni                                                   
      if(ju(j).eq.0)goto 21091                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                          
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2              
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                    
21091 continue                                                          
      continue                                                          
      if(isd .ne. 0)goto 21111                                          
      xs=1.0                                                            
      goto 21121                                                        
21111 continue                                                          
      xv=1.0                                                            
21121 continue                                                          
      continue                                                          
      ys0=0.0                                                           
      do 21131 j=1,nr                                                   
      ym(j)=dot_product(w,y(:,j))                                       
      y(:,j)=y(:,j)-ym(j)                                               
      z=dot_product(w,y(:,j)**2)                                        
      if(jsd .le. 0)goto 21151                                          
      ys(j)=sqrt(z)                                                     
      y(:,j)=y(:,j)/ys(j)                                               
      goto 21161                                                        
21151 continue                                                          
      ys0=ys0+z                                                         
21161 continue                                                          
      continue                                                          
21131 continue                                                          
      continue                                                          
      if(jsd .ne. 0)goto 21181                                          
      ys=1.0                                                            
      goto 21191                                                        
21181 continue                                                          
      ys0=nr                                                            
21191 continue                                                          
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
      if(flmin .ge. 1.0)goto 21211                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
21211 continue                                                          
      rsq=ys0                                                           
      a=0.0                                                             
      mm=0                                                              
      o=0.0                                                             
      nlp=0                                                             
      nin=nlp                                                           
      iz=0                                                              
      mnl=min(mnlam,nlam)                                               
      do 21221 j=1,ni                                                   
      if(ju(j).eq.0)goto 21221                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      g(j)=0.0                                                          
      do 21231 k=1,nr                                                   
      g(j)=g(j)+(dot_product(y(jx(jb:je),k),w(jx(jb:je))*x(jb:je))/xs(j)
     *)**2
21231 continue                                                          
      continue                                                          
      g(j)=sqrt(g(j))                                                   
21221 continue                                                          
      continue                                                          
      do 21241 m=1,nlam                                                 
      if(itrace.ne.0) call setpb(m-1)                                   
      alm0=alm                                                          
      if(flmin .lt. 1.0)goto 21261                                      
      alm=ulam(m)                                                       
      goto 21251                                                        
21261 if(m .le. 2)goto 21271                                            
      alm=alm*alf                                                       
      goto 21251                                                        
21271 if(m .ne. 1)goto 21281                                            
      alm=big                                                           
      goto 21291                                                        
21281 continue                                                          
      alm0=0.0                                                          
      do 21301 j=1,ni                                                   
      if(ju(j).eq.0)goto 21301                                          
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                        
21301 continue                                                          
      continue                                                          
      alm0=alm0/max(bta,1.0d-3)                                         
      alm=alf*alm0                                                      
21291 continue                                                          
21251 continue                                                          
      dem=alm*omb                                                       
      ab=alm*bta                                                        
      rsq0=rsq                                                          
      jz=1                                                              
      tlam=bta*(2.0*alm-alm0)                                           
      do 21311 k=1,ni                                                   
      if(iy(k).eq.1)goto 21311                                          
      if(ju(k).eq.0)goto 21311                                          
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                    
21311 continue                                                          
      continue                                                          
      continue                                                          
21321 continue                                                          
      if(iz*jz.ne.0) go to 10360                                        
10880 if(nlp .le. maxit)goto 21341                                      
      jerr=-m                                                           
      return                                                            
21341 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 21351 k=1,ni                                                   
      if(iy(k).eq.0)goto 21351                                          
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      gkn=0.0                                                           
      do 21361 j=1,nr                                                   
      gj(j)=dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(k)
      gk(j)=gj(j)+a(j,k)*xv(k)                                          
      gkn=gkn+gk(j)**2                                                  
21361 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-ab*vp(k)/gkn                                                
      del=a(:,k)                                                        
      if(u .gt. 0.0)goto 21381                                          
      a(:,k)=0.0                                                        
      goto 21391                                                        
21381 continue                                                          
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                   
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)
     *,isc,jerr)
      if(jerr.ne.0) return                                              
21391 continue                                                          
      continue                                                          
      del=a(:,k)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 21351                             
      if(mm(k) .ne. 0)goto 21411                                        
      nin=nin+1                                                         
      if(nin.gt.nx)goto 21352                                           
      mm(k)=nin                                                         
      ia(nin)=k                                                         
21411 continue                                                          
      do 21421 j=1,nr                                                   
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                           
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)               
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                      
      dlx=max(xv(k)*del(j)**2,dlx)                                      
21421 continue                                                          
      continue                                                          
21351 continue                                                          
21352 continue                                                          
      if(nin.gt.nx)goto 21322                                           
      if(dlx .ge. thr)goto 21441                                        
      ixx=0                                                             
      do 21451 j=1,ni                                                   
      if(iy(j).eq.1)goto 21451                                          
      if(ju(j).eq.0)goto 21451                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      g(j)=0.0                                                          
      do 21461 k=1,nr                                                   
      g(j)=g(j)+  (dot_product(y(jx(jb:je),k)+o(k),w(jx(jb:je))*x(jb:je)
     *)/xs(j))**2
21461 continue                                                          
      continue                                                          
      g(j)=sqrt(g(j))                                                   
      if(g(j) .le. ab*vp(j))goto 21481                                  
      iy(j)=1                                                           
      ixx=1                                                             
21481 continue                                                          
21451 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10880                                          
      goto 21322                                                        
21441 continue                                                          
      if(nlp .le. maxit)goto 21501                                      
      jerr=-m                                                           
      return                                                            
21501 continue                                                          
10360 continue                                                          
      iz=1                                                              
      continue                                                          
21511 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 21521 l=1,nin                                                  
      k=ia(l)                                                           
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      gkn=0.0                                                           
      do 21531 j=1,nr                                                   
      gj(j)=  dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(
     *k)
      gk(j)=gj(j)+a(j,k)*xv(k)                                          
      gkn=gkn+gk(j)**2                                                  
21531 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-ab*vp(k)/gkn                                                
      del=a(:,k)                                                        
      if(u .gt. 0.0)goto 21551                                          
      a(:,k)=0.0                                                        
      goto 21561                                                        
21551 continue                                                          
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                   
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)
     *,isc,jerr)
      if(jerr.ne.0) return                                              
21561 continue                                                          
      continue                                                          
      del=a(:,k)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 21521                             
      do 21571 j=1,nr                                                   
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                           
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)               
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                      
      dlx=max(xv(k)*del(j)**2,dlx)                                      
21571 continue                                                          
      continue                                                          
21521 continue                                                          
      continue                                                          
      if(dlx.lt.thr)goto 21512                                          
      if(nlp .le. maxit)goto 21591                                      
      jerr=-m                                                           
      return                                                            
21591 continue                                                          
      goto 21511                                                        
21512 continue                                                          
      jz=0                                                              
      goto 21321                                                        
21322 continue                                                          
      if(nin .le. nx)goto 21611                                         
      jerr=-10000-m                                                     
      goto 21242                                                        
21611 continue                                                          
      if(nin .le. 0)goto 21631                                          
      do 21641 j=1,nr                                                   
      ao(1:nin,j,m)=a(j,ia(1:nin))                                      
21641 continue                                                          
      continue                                                          
21631 continue                                                          
      kin(m)=nin                                                        
      rsqo(m)=1.0-rsq/ys0                                               
      almo(m)=alm                                                       
      lmu=m                                                             
      if(m.lt.mnl)goto 21241                                            
      if(flmin.ge.1.0)goto 21241                                        
      me=0                                                              
      do 21651 j=1,nin                                                  
      if(ao(j,1,m).ne.0.0) me=me+1                                      
21651 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 21242                                            
      if(rsq0-rsq.lt.sml*rsq)goto 21242                                 
      if(rsqo(m).gt.rsqmax)goto 21242                                   
21241 continue                                                          
21242 continue                                                          
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
      do 21661 ic=1,nc                                                  
      q0=dot_product(w,y(:,ic))                                         
      if(q0 .gt. pmin)goto 21681                                        
      jerr =8000+ic                                                     
      return                                                            
21681 continue                                                          
      if(q0 .lt. pmax)goto 21701                                        
      jerr =9000+ic                                                     
      return                                                            
21701 continue                                                          
      if(intr .ne. 0)goto 21721                                         
      q0=1.0/nc                                                         
      b(0,ic)=0.0                                                       
      goto 21731                                                        
21721 continue                                                          
      b(0,ic)=log(q0)                                                   
      dev1=dev1-q0*b(0,ic)                                              
21731 continue                                                          
      continue                                                          
      b(1:ni,ic)=0.0                                                    
21661 continue                                                          
      continue                                                          
      if(intr.eq.0) dev1=log(float(nc))                                 
      ixx=0                                                             
      al=0.0                                                            
      if(nonzero(no*nc,g) .ne. 0)goto 21751                             
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                      
      sxp=0.0                                                           
      do 21761 ic=1,nc                                                  
      q(:,ic)=exp(b(0,ic))                                              
      sxp=sxp+q(:,ic)                                                   
21761 continue                                                          
      continue                                                          
      goto 21771                                                        
21751 continue                                                          
      do 21781 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
21781 continue                                                          
      continue                                                          
      sxp=0.0                                                           
      if(intr .ne. 0)goto 21801                                         
      b(0,:)=0.0                                                        
      goto 21811                                                        
21801 continue                                                          
      call kazero(nc,no,y,g,w,b(0,:),jerr)                              
      if(jerr.ne.0) return                                              
21811 continue                                                          
      continue                                                          
      dev1=0.0                                                          
      do 21821 ic=1,nc                                                  
      q(:,ic)=b(0,ic)+g(:,ic)                                           
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                          
      q(:,ic)=exp(q(:,ic))                                              
      sxp=sxp+q(:,ic)                                                   
21821 continue                                                          
      continue                                                          
      sxpl=w*log(sxp)                                                   
      do 21831 ic=1,nc                                                  
      dev1=dev1+dot_product(y(:,ic),sxpl)                               
21831 continue                                                          
      continue                                                          
21771 continue                                                          
      continue                                                          
      do 21841 ic=1,nc                                                  
      do 21851 i=1,no                                                   
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))            
21851 continue                                                          
      continue                                                          
21841 continue                                                          
      continue                                                          
      dev0=dev0+dev1                                                    
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 21871                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
21871 continue                                                          
      m=0                                                               
      mm=0                                                              
      nin=0                                                             
      nlp=0                                                             
      mnl=min(mnlam,nlam)                                               
      bs=0.0                                                            
      shr=shri*dev0                                                     
      ga=0.0                                                            
      do 21881 ic=1,nc                                                  
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                   
      do 21891 j=1,ni                                                   
      if(ju(j).ne.0) ga(j)=ga(j)+dot_product(r(:,ic),x(:,j))**2         
21891 continue                                                          
      continue                                                          
21881 continue                                                          
      continue                                                          
      ga=sqrt(ga)                                                       
      do 21901 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 21921                                      
      al=ulam(ilm)                                                      
      goto 21911                                                        
21921 if(ilm .le. 2)goto 21931                                          
      al=al*alf                                                         
      goto 21911                                                        
21931 if(ilm .ne. 1)goto 21941                                          
      al=big                                                            
      goto 21951                                                        
21941 continue                                                          
      al0=0.0                                                           
      do 21961 j=1,ni                                                   
      if(ju(j).eq.0)goto 21961                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
21961 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
21951 continue                                                          
21911 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 21971 k=1,ni                                                   
      if(ixx(k).eq.1)goto 21971                                         
      if(ju(k).eq.0)goto 21971                                          
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                  
21971 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
21981 continue                                                          
      ix=0                                                              
      jx=ix                                                             
      kx=jx                                                             
      t=0.0                                                             
      if(nlp .le. maxit)goto 22001                                      
      jerr=-ilm                                                         
      return                                                            
22001 continue                                                          
      do 22011 ic=1,nc                                                  
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                    
22011 continue                                                          
      continue                                                          
      if(t .ge. eps)goto 22031                                          
      kx=1                                                              
      goto 21982                                                        
22031 continue                                                          
      t=2.0*t                                                           
      alt=al1/t                                                         
      al2t=al2/t                                                        
      do 22041 ic=1,nc                                                  
      bs(0,ic)=b(0,ic)                                                  
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                       
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                 
      d=0.0                                                             
      if(intr.ne.0) d=sum(r(:,ic))                                      
      if(d .eq. 0.0)goto 22061                                          
      b(0,ic)=b(0,ic)+d                                                 
      r(:,ic)=r(:,ic)-d*w                                               
      dlx=max(dlx,d**2)                                                 
22061 continue                                                          
22041 continue                                                          
      continue                                                          
      continue                                                          
22071 continue                                                          
      nlp=nlp+nc                                                        
      dlx=0.0                                                           
      do 22081 k=1,ni                                                   
      if(ixx(k).eq.0)goto 22081                                         
      gkn=0.0                                                           
      do 22091 ic=1,nc                                                  
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                  
      gkn=gkn+gk(ic)**2                                                 
22091 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-alt*vp(k)/gkn                                               
      del=b(k,:)                                                        
      if(u .gt. 0.0)goto 22111                                          
      b(k,:)=0.0                                                        
      goto 22121                                                        
22111 continue                                                          
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                  
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                              
22121 continue                                                          
      continue                                                          
      del=b(k,:)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 22081                             
      do 22131 ic=1,nc                                                  
      dlx=max(dlx,xv(k)*del(ic)**2)                                     
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                  
22131 continue                                                          
      continue                                                          
      if(mm(k) .ne. 0)goto 22151                                        
      nin=nin+1                                                         
      if(nin .le. nx)goto 22171                                         
      jx=1                                                              
      goto 22082                                                        
22171 continue                                                          
      mm(k)=nin                                                         
      m(nin)=k                                                          
22151 continue                                                          
22081 continue                                                          
22082 continue                                                          
      if(jx.gt.0)goto 22072                                             
      if(dlx.lt.shr)goto 22072                                          
      if(nlp .le. maxit)goto 22191                                      
      jerr=-ilm                                                         
      return                                                            
22191 continue                                                          
      continue                                                          
22201 continue                                                          
      nlp=nlp+nc                                                        
      dlx=0.0                                                           
      do 22211 l=1,nin                                                  
      k=m(l)                                                            
      gkn=0.0                                                           
      do 22221 ic=1,nc                                                  
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                  
      gkn=gkn+gk(ic)**2                                                 
22221 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-alt*vp(k)/gkn                                               
      del=b(k,:)                                                        
      if(u .gt. 0.0)goto 22241                                          
      b(k,:)=0.0                                                        
      goto 22251                                                        
22241 continue                                                          
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                  
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                              
22251 continue                                                          
      continue                                                          
      del=b(k,:)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 22211                             
      do 22261 ic=1,nc                                                  
      dlx=max(dlx,xv(k)*del(ic)**2)                                     
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                  
22261 continue                                                          
      continue                                                          
22211 continue                                                          
      continue                                                          
      if(dlx.lt.shr)goto 22202                                          
      if(nlp .le. maxit)goto 22281                                      
      jerr=-ilm                                                         
      return                                                            
22281 continue                                                          
      goto 22201                                                        
22202 continue                                                          
      goto 22071                                                        
22072 continue                                                          
      if(jx.gt.0)goto 21982                                             
      do 22291 ic=1,nc                                                  
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                             
      if(ix .ne. 0)goto 22311                                           
      do 22321 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 22341                
      ix=1                                                              
      goto 22322                                                        
22341 continue                                                          
22321 continue                                                          
22322 continue                                                          
22311 continue                                                          
      do 22351 i=1,no                                                   
      fi=b(0,ic)+g(i,ic)                                                
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))      
      fi=min(max(exmn,fi),exmx)                                         
      sxp(i)=sxp(i)-q(i,ic)                                             
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                 
      sxp(i)=sxp(i)+q(i,ic)                                             
22351 continue                                                          
      continue                                                          
22291 continue                                                          
      continue                                                          
      s=-sum(b(0,:))/nc                                                 
      b(0,:)=b(0,:)+s                                                   
      if(jx.gt.0)goto 21982                                             
      if(ix .ne. 0)goto 22371                                           
      do 22381 k=1,ni                                                   
      if(ixx(k).eq.1)goto 22381                                         
      if(ju(k).eq.0)goto 22381                                          
      ga(k)=0.0                                                         
22381 continue                                                          
      continue                                                          
      do 22391 ic=1,nc                                                  
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                   
      do 22401 k=1,ni                                                   
      if(ixx(k).eq.1)goto 22401                                         
      if(ju(k).eq.0)goto 22401                                          
      ga(k)=ga(k)+dot_product(r(:,ic),x(:,k))**2                        
22401 continue                                                          
      continue                                                          
22391 continue                                                          
      continue                                                          
      ga=sqrt(ga)                                                       
      do 22411 k=1,ni                                                   
      if(ixx(k).eq.1)goto 22411                                         
      if(ju(k).eq.0)goto 22411                                          
      if(ga(k) .le. al1*vp(k))goto 22431                                
      ixx(k)=1                                                          
      ix=1                                                              
22431 continue                                                          
22411 continue                                                          
      continue                                                          
      if(ix.eq.1) go to 10880                                           
      goto 21982                                                        
22371 continue                                                          
      goto 21981                                                        
21982 continue                                                          
      if(kx .le. 0)goto 22451                                           
      jerr=-20000-ilm                                                   
      goto 21902                                                        
22451 continue                                                          
      if(jx .le. 0)goto 22471                                           
      jerr=-10000-ilm                                                   
      goto 21902                                                        
22471 continue                                                          
      devi=0.0                                                          
      do 22481 ic=1,nc                                                  
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                       
      a0(ic,ilm)=b(0,ic)                                                
      do 22491 i=1,no                                                   
      if(y(i,ic).le.0.0)goto 22491                                      
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                        
22491 continue                                                          
      continue                                                          
22481 continue                                                          
      continue                                                          
      kin(ilm)=nin                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dev1-devi)/dev0                                         
      if(ilm.lt.mnl)goto 21901                                          
      if(flmin.ge.1.0)goto 21901                                        
      me=0                                                              
      do 22501 j=1,nin                                                  
      if(a(j,1,ilm).ne.0.0) me=me+1                                     
22501 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 21902                                            
      if(dev(ilm).gt.devmax)goto 21902                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 21902                          
21901 continue                                                          
21902 continue                                                          
      g=log(q)                                                          
      do 22511 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
22511 continue                                                          
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
      do 22521 ic=1,nc                                                  
      q0=dot_product(w,y(:,ic))                                         
      if(q0 .gt. pmin)goto 22541                                        
      jerr =8000+ic                                                     
      return                                                            
22541 continue                                                          
      if(q0 .lt. pmax)goto 22561                                        
      jerr =9000+ic                                                     
      return                                                            
22561 continue                                                          
      b(1:ni,ic)=0.0                                                    
      if(intr .ne. 0)goto 22581                                         
      q0=1.0/nc                                                         
      b(0,ic)=0.0                                                       
      goto 22591                                                        
22581 continue                                                          
      b(0,ic)=log(q0)                                                   
      dev1=dev1-q0*b(0,ic)                                              
22591 continue                                                          
      continue                                                          
22521 continue                                                          
      continue                                                          
      if(intr.eq.0) dev1=log(float(nc))                                 
      iy=0                                                              
      al=0.0                                                            
      if(nonzero(no*nc,g) .ne. 0)goto 22611                             
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                      
      sxp=0.0                                                           
      do 22621 ic=1,nc                                                  
      q(:,ic)=exp(b(0,ic))                                              
      sxp=sxp+q(:,ic)                                                   
22621 continue                                                          
      continue                                                          
      goto 22631                                                        
22611 continue                                                          
      do 22641 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
22641 continue                                                          
      continue                                                          
      sxp=0.0                                                           
      if(intr .ne. 0)goto 22661                                         
      b(0,:)=0.0                                                        
      goto 22671                                                        
22661 continue                                                          
      call kazero(nc,no,y,g,w,b(0,:),jerr)                              
      if(jerr.ne.0) return                                              
22671 continue                                                          
      continue                                                          
      dev1=0.0                                                          
      do 22681 ic=1,nc                                                  
      q(:,ic)=b(0,ic)+g(:,ic)                                           
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                          
      q(:,ic)=exp(q(:,ic))                                              
      sxp=sxp+q(:,ic)                                                   
22681 continue                                                          
      continue                                                          
      sxpl=w*log(sxp)                                                   
      do 22691 ic=1,nc                                                  
      dev1=dev1+dot_product(y(:,ic),sxpl)                               
22691 continue                                                          
      continue                                                          
22631 continue                                                          
      continue                                                          
      do 22701 ic=1,nc                                                  
      do 22711 i=1,no                                                   
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))            
22711 continue                                                          
      continue                                                          
22701 continue                                                          
      continue                                                          
      dev0=dev0+dev1                                                    
      alf=1.0                                                           
      if(flmin .ge. 1.0)goto 22731                                      
      eqs=max(eps,flmin)                                                
      alf=eqs**(1.0/(nlam-1))                                           
22731 continue                                                          
      m=0                                                               
      mm=0                                                              
      nin=0                                                             
      nlp=0                                                             
      mnl=min(mnlam,nlam)                                               
      bs=0.0                                                            
      shr=shri*dev0                                                     
      ga=0.0                                                            
      do 22741 ic=1,nc                                                  
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                   
      svr(ic)=sum(r(:,ic))                                              
      do 22751 j=1,ni                                                   
      if(ju(j).eq.0)goto 22751                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                          
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                         
22751 continue                                                          
      continue                                                          
22741 continue                                                          
      continue                                                          
      ga=sqrt(ga)                                                       
      do 22761 ilm=1,nlam                                               
      if(itrace.ne.0) call setpb(ilm-1)                                 
      al0=al                                                            
      if(flmin .lt. 1.0)goto 22781                                      
      al=ulam(ilm)                                                      
      goto 22771                                                        
22781 if(ilm .le. 2)goto 22791                                          
      al=al*alf                                                         
      goto 22771                                                        
22791 if(ilm .ne. 1)goto 22801                                          
      al=big                                                            
      goto 22811                                                        
22801 continue                                                          
      al0=0.0                                                           
      do 22821 j=1,ni                                                   
      if(ju(j).eq.0)goto 22821                                          
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                         
22821 continue                                                          
      continue                                                          
      al0=al0/max(bta,1.0d-3)                                           
      al=alf*al0                                                        
22811 continue                                                          
22771 continue                                                          
      al2=al*omb                                                        
      al1=al*bta                                                        
      tlam=bta*(2.0*al-al0)                                             
      do 22831 k=1,ni                                                   
      if(iy(k).eq.1)goto 22831                                          
      if(ju(k).eq.0)goto 22831                                          
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                   
22831 continue                                                          
      continue                                                          
10880 continue                                                          
      continue                                                          
22841 continue                                                          
      ixx=0                                                             
      jxx=ixx                                                           
      kxx=jxx                                                           
      t=0.0                                                             
      if(nlp .le. maxit)goto 22861                                      
      jerr=-ilm                                                         
      return                                                            
22861 continue                                                          
      do 22871 ic=1,nc                                                  
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                    
22871 continue                                                          
      continue                                                          
      if(t .ge. eps)goto 22891                                          
      kxx=1                                                             
      goto 22842                                                        
22891 continue                                                          
      t=2.0*t                                                           
      alt=al1/t                                                         
      al2t=al2/t                                                        
      do 22901 ic=1,nc                                                  
      bs(0,ic)=b(0,ic)                                                  
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                       
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                 
      svr(ic)=sum(r(:,ic))                                              
      if(intr .eq. 0)goto 22921                                         
      b(0,ic)=b(0,ic)+svr(ic)                                           
      r(:,ic)=r(:,ic)-svr(ic)*w                                         
      dlx=max(dlx,svr(ic)**2)                                           
22921 continue                                                          
22901 continue                                                          
      continue                                                          
      continue                                                          
22931 continue                                                          
      nlp=nlp+nc                                                        
      dlx=0.0                                                           
      do 22941 k=1,ni                                                   
      if(iy(k).eq.0)goto 22941                                          
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      del=b(k,:)                                                        
      gkn=0.0                                                           
      do 22951 ic=1,nc                                                  
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))-svr(ic)*xb(k))/xs(k)     
      gk(ic)=u+del(ic)*xv(k)                                            
      gkn=gkn+gk(ic)**2                                                 
22951 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-alt*vp(k)/gkn                                               
      if(u .gt. 0.0)goto 22971                                          
      b(k,:)=0.0                                                        
      goto 22981                                                        
22971 continue                                                          
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                  
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                              
22981 continue                                                          
      continue                                                          
      del=b(k,:)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 22941                             
      do 22991 ic=1,nc                                                  
      dlx=max(dlx,xv(k)*del(ic)**2)                                     
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x
     *b(k))/xs(k)
22991 continue                                                          
      continue                                                          
      if(mm(k) .ne. 0)goto 23011                                        
      nin=nin+1                                                         
      if(nin .le. nx)goto 23031                                         
      jxx=1                                                             
      goto 22942                                                        
23031 continue                                                          
      mm(k)=nin                                                         
      m(nin)=k                                                          
23011 continue                                                          
22941 continue                                                          
22942 continue                                                          
      if(jxx.gt.0)goto 22932                                            
      if(dlx.lt.shr)goto 22932                                          
      if(nlp .le. maxit)goto 23051                                      
      jerr=-ilm                                                         
      return                                                            
23051 continue                                                          
      continue                                                          
23061 continue                                                          
      nlp=nlp+nc                                                        
      dlx=0.0                                                           
      do 23071 l=1,nin                                                  
      k=m(l)                                                            
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      del=b(k,:)                                                        
      gkn=0.0                                                           
      do 23081 ic=1,nc                                                  
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))  -svr(ic)*xb(k))/xs(k)   
      gk(ic)=u+del(ic)*xv(k)                                            
      gkn=gkn+gk(ic)**2                                                 
23081 continue                                                          
      continue                                                          
      gkn=sqrt(gkn)                                                     
      u=1.0-alt*vp(k)/gkn                                               
      if(u .gt. 0.0)goto 23101                                          
      b(k,:)=0.0                                                        
      goto 23111                                                        
23101 continue                                                          
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                  
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                              
23111 continue                                                          
      continue                                                          
      del=b(k,:)-del                                                    
      if(maxval(abs(del)).le.0.0)goto 23071                             
      do 23121 ic=1,nc                                                  
      dlx=max(dlx,xv(k)*del(ic)**2)                                     
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x
     *b(k))/xs(k)
23121 continue                                                          
      continue                                                          
23071 continue                                                          
      continue                                                          
      if(dlx.lt.shr)goto 23062                                          
      if(nlp .le. maxit)goto 23141                                      
      jerr=-ilm                                                         
      return                                                            
23141 continue                                                          
      goto 23061                                                        
23062 continue                                                          
      goto 22931                                                        
22932 continue                                                          
      if(jxx.gt.0)goto 22842                                            
      do 23151 ic=1,nc                                                  
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                            
      if(ixx .ne. 0)goto 23171                                          
      do 23181 j=1,nin                                                  
      k=m(j)                                                            
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 23201                
      ixx=1                                                             
      goto 23182                                                        
23201 continue                                                          
23181 continue                                                          
23182 continue                                                          
23171 continue                                                          
      sc=b(0,ic)+g(:,ic)                                                
      b0=0.0                                                            
      do 23211 j=1,nin                                                  
      l=m(j)                                                            
      jb=ix(l)                                                          
      je=ix(l+1)-1                                                      
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                
      b0=b0-b(l,ic)*xb(l)/xs(l)                                         
23211 continue                                                          
      continue                                                          
      sc=min(max(exmn,sc+b0),exmx)                                      
      sxp=sxp-q(:,ic)                                                   
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                       
      sxp=sxp+q(:,ic)                                                   
23151 continue                                                          
      continue                                                          
      s=sum(b(0,:))/nc                                                  
      b(0,:)=b(0,:)-s                                                   
      if(jxx.gt.0)goto 22842                                            
      if(ixx .ne. 0)goto 23231                                          
      do 23241 j=1,ni                                                   
      if(iy(j).eq.1)goto 23241                                          
      if(ju(j).eq.0)goto 23241                                          
      ga(j)=0.0                                                         
23241 continue                                                          
      continue                                                          
      do 23251 ic=1,nc                                                  
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                   
      do 23261 j=1,ni                                                   
      if(iy(j).eq.1)goto 23261                                          
      if(ju(j).eq.0)goto 23261                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                          
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                         
23261 continue                                                          
      continue                                                          
23251 continue                                                          
      continue                                                          
      ga=sqrt(ga)                                                       
      do 23271 k=1,ni                                                   
      if(iy(k).eq.1)goto 23271                                          
      if(ju(k).eq.0)goto 23271                                          
      if(ga(k) .le. al1*vp(k))goto 23291                                
      iy(k)=1                                                           
      ixx=1                                                             
23291 continue                                                          
23271 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10880                                          
      goto 22842                                                        
23231 continue                                                          
      goto 22841                                                        
22842 continue                                                          
      if(kxx .le. 0)goto 23311                                          
      jerr=-20000-ilm                                                   
      goto 22762                                                        
23311 continue                                                          
      if(jxx .le. 0)goto 23331                                          
      jerr=-10000-ilm                                                   
      goto 22762                                                        
23331 continue                                                          
      devi=0.0                                                          
      do 23341 ic=1,nc                                                  
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                       
      a0(ic,ilm)=b(0,ic)                                                
      do 23351 i=1,no                                                   
      if(y(i,ic).le.0.0)goto 23351                                      
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                        
23351 continue                                                          
      continue                                                          
23341 continue                                                          
      continue                                                          
      kin(ilm)=nin                                                      
      alm(ilm)=al                                                       
      lmu=ilm                                                           
      dev(ilm)=(dev1-devi)/dev0                                         
      if(ilm.lt.mnl)goto 22761                                          
      if(flmin.ge.1.0)goto 22761                                        
      me=0                                                              
      do 23361 j=1,nin                                                  
      if(a(j,1,ilm).ne.0.0) me=me+1                                     
23361 continue                                                          
      continue                                                          
      if(me.gt.ne)goto 22762                                            
      if(dev(ilm).gt.devmax)goto 22762                                  
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 22762                          
22761 continue                                                          
22762 continue                                                          
      g=log(q)                                                          
      do 23371 i=1,no                                                   
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                      
23371 continue                                                          
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
