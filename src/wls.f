c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))             
      subroutine wls(alm0,almc,alpha,m,no,ni,x,r,v,intr,ju,vp,cl,nx,thr,
     *maxit,  a,aint,g,ia,iy,iz,mm,nino,rsqc,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(no,ni),r(no),a(ni),vp(ni),cl(2,ni)             
      double precision v(no),g(ni)                                      
      integer iy(ni),ia(nx),ju(ni),mm(ni)                               
      double precision, dimension (:), allocatable :: xv                
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      do 10011 j=1,ni                                                   
      if(ju(j).eq.0)goto 10011                                          
      g(j)=abs(dot_product(r,x(:,j)))                                   
10011 continue                                                          
      continue                                                          
      do 10021 j=1,ni                                                   
      if(iy(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                     
10021 continue                                                          
      continue                                                          
      xmz = sum(v)                                                      
      ab=almc*alpha                                                     
      dem=almc*(1.0-alpha)                                              
      tlam=alpha*(2.0*almc-alm0)                                        
      do 10031 k=1,ni                                                   
      if(iy(k).eq.1)goto 10031                                          
      if(ju(k).eq.0)goto 10031                                          
      if(g(k) .le. tlam*vp(k))goto 10051                                
      iy(k)=1                                                           
      xv(k)=dot_product(v,x(:,k)**2)                                    
10051 continue                                                          
10031 continue                                                          
      continue                                                          
      jz = 1                                                            
      continue                                                          
10061 continue                                                          
      if(iz*jz.ne.0) go to 10070                                        
10080 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10091 k=1,ni                                                   
      if(iy(k).eq.0)goto 10091                                          
      gk=dot_product(r,x(:,k))                                          
      ak=a(k)                                                           
      u=gk+ak*xv(k)                                                     
      au=abs(u)-vp(k)*ab                                                
      if(au .gt. 0.0)goto 10111                                         
      a(k)=0.0                                                          
      goto 10121                                                        
10111 continue                                                          
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*dem)))       
10121 continue                                                          
      continue                                                          
      if(a(k).eq.ak)goto 10091                                          
      if(mm(k) .ne. 0)goto 10141                                        
      nino=nino+1                                                       
      if(nino.gt.nx)goto 10092                                          
      mm(k)=nino                                                        
      ia(nino)=k                                                        
10141 continue                                                          
      d=a(k)-ak                                                         
      rsqc=rsqc+d*(2.0*gk-d*xv(k))                                      
      r=r-d*v*x(:,k)                                                    
      dlx=max(xv(k)*d**2,dlx)                                           
10091 continue                                                          
10092 continue                                                          
      if(nino.gt.nx)goto 10062                                          
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 10161                                          
      aint=aint+d                                                       
      rsqc=rsqc+d*(2.0*sum(r)-d*xmz)                                    
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
10161 continue                                                          
      if(dlx .ge. thr)goto 10181                                        
      ixx=0                                                             
      do 10191 k=1,ni                                                   
      if(iy(k).eq.1)goto 10191                                          
      if(ju(k).eq.0)goto 10191                                          
      g(k)=abs(dot_product(r,x(:,k)))                                   
      if(g(k) .le. ab*vp(k))goto 10211                                  
      iy(k)=1                                                           
      xv(k)=dot_product(v,x(:,k)**2)                                    
      ixx=1                                                             
10211 continue                                                          
10191 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10080                                          
      goto 10062                                                        
10181 continue                                                          
      if(nlp .le. maxit)goto 10231                                      
      jerr=-m                                                           
      return                                                            
10231 continue                                                          
10070 continue                                                          
      iz = 1                                                            
      continue                                                          
10241 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10251 l=1,nino                                                 
      k=ia(l)                                                           
      gk=dot_product(r,x(:,k))                                          
      ak=a(k)                                                           
      u=gk+ak*xv(k)                                                     
      au=abs(u)-vp(k)*ab                                                
      if(au .gt. 0.0)goto 10271                                         
      a(k)=0.0                                                          
      goto 10281                                                        
10271 continue                                                          
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*dem)))       
10281 continue                                                          
      continue                                                          
      if(a(k).eq.ak)goto 10251                                          
      d=a(k)-ak                                                         
      rsqc=rsqc+d*(2.0*gk-d*xv(k))                                      
      r=r-d*v*x(:,k)                                                    
      dlx=max(xv(k)*d**2,dlx)                                           
10251 continue                                                          
      continue                                                          
      d=0.0                                                             
      if(intr.ne.0) d=sum(r)/xmz                                        
      if(d .eq. 0.0)goto 10301                                          
      aint=aint+d                                                       
      rsqc=rsqc+d*(2.0*sum(r)-d*xmz)                                    
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
10301 continue                                                          
      if(dlx.lt.thr)goto 10242                                          
      if(nlp .le. maxit)goto 10321                                      
      jerr=-m                                                           
      return                                                            
10321 continue                                                          
      goto 10241                                                        
10242 continue                                                          
      jz=0                                                              
      goto 10061                                                        
10062 continue                                                          
      deallocate(xv)                                                    
      return                                                            
      end                                                               
      subroutine spwls(alm0,almc,alpha,m,no,ni,x,ix,jx,xm,xs,r,v,intr,ju
     *,  vp,cl,nx,thr,maxit,a,aint,g,ia,iy,iz,mm,nino,rsqc,nlp,jerr)
      implicit double precision(a-h,o-z)                                
      double precision x(*),xm(ni),xs(ni),r(no),a(ni),vp(ni),cl(2,ni)   
      double precision v(no),g(ni)                                      
      integer ix(*),jx(*),iy(ni),ia(nx),ju(ni),mm(ni)                   
      double precision, dimension (:), allocatable :: xv                
      allocate(xv(1:ni),stat=jerr)                                      
      if(jerr.ne.0) return                                              
      xmz = sum(v)                                                      
      rsum = sum(r)                                                     
      do 10331 j=1,ni                                                   
      if(ju(j).eq.0)goto 10331                                          
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      g(j)=abs(dot_product(r(jx(jb:je)),x(jb:je))-rsum*xm(j))/xs(j)     
10331 continue                                                          
      continue                                                          
      do 10341 j=1,ni                                                   
      if(iy(j) .le. 0)goto 10361                                        
      jb=ix(j)                                                          
      je=ix(j+1)-1                                                      
      xv(j)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       
      xv(j)=xv(j)-2*xm(j)*dot_product(v(jx(jb:je)),x(jb:je))            
      xv(j)=(xv(j)+xmz*xm(j)**2)/xs(j)**2                               
10361 continue                                                          
10341 continue                                                          
      continue                                                          
      ab=almc*alpha                                                     
      dem=almc*(1.0-alpha)                                              
      tlam=alpha*(2.0*almc-alm0)                                        
      do 10371 k=1,ni                                                   
      if(iy(k).eq.1)goto 10371                                          
      if(ju(k).eq.0)goto 10371                                          
      if(g(k) .le. tlam*vp(k))goto 10391                                
      iy(k)=1                                                           
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      xv(k)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       
      xv(k)=xv(k)-2*xm(k)*dot_product(v(jx(jb:je)),x(jb:je))            
      xv(k)=(xv(k)+xmz*xm(k)**2)/xs(k)**2                               
10391 continue                                                          
10371 continue                                                          
      continue                                                          
      jz = 1                                                            
      continue                                                          
10401 continue                                                          
      if(iz*jz.ne.0) go to 10070                                        
10080 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10411 k=1,ni                                                   
      if(iy(k).eq.0)goto 10411                                          
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      gk=(dot_product(r(jx(jb:je)),x(jb:je))-rsum*xm(k))/xs(k)          
      ak=a(k)                                                           
      u=gk+ak*xv(k)                                                     
      au=abs(u)-vp(k)*ab                                                
      if(au .gt. 0.0)goto 10431                                         
      a(k)=0.0                                                          
      goto 10441                                                        
10431 continue                                                          
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*dem)))       
10441 continue                                                          
      continue                                                          
      if(a(k).eq.ak)goto 10411                                          
      if(mm(k) .ne. 0)goto 10461                                        
      nino=nino+1                                                       
      if(nino.gt.nx)goto 10412                                          
      mm(k)=nino                                                        
      ia(nino)=k                                                        
10461 continue                                                          
      d=a(k)-ak                                                         
      rsqc=rsqc+d*(2.0*gk-d*xv(k))                                      
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)           
      r=r+d*v*xm(k)/xs(k)                                               
      rsum=sum(r)                                                       
      dlx=max(xv(k)*d**2,dlx)                                           
10411 continue                                                          
10412 continue                                                          
      if(nino.gt.nx)goto 10402                                          
      d=0.0                                                             
      if(intr.ne.0) d=rsum/xmz                                          
      if(d .eq. 0.0)goto 10481                                          
      aint=aint+d                                                       
      rsqc=rsqc+d*(2.0*rsum-d*xmz)                                      
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
      rsum=sum(r)                                                       
10481 continue                                                          
      if(dlx .ge. thr)goto 10501                                        
      ixx=0                                                             
      do 10511 k=1,ni                                                   
      if(iy(k).eq.1)goto 10511                                          
      if(ju(k).eq.0)goto 10511                                          
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      g(k)=dot_product(r(jx(jb:je)),x(jb:je))                           
      g(k)=abs(g(k)-rsum*xm(k))/xs(k)                                   
      if(g(k) .le. ab*vp(k))goto 10531                                  
      iy(k)=1                                                           
      xv(k)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       
      vx=dot_product(v(jx(jb:je)),x(jb:je))                             
      xv(k)=xv(k)-2*xm(k)*vx                                            
      xv(k)=(xv(k)+xmz*xm(k)**2)/xs(k)**2                               
      ixx=1                                                             
10531 continue                                                          
10511 continue                                                          
      continue                                                          
      if(ixx.eq.1) go to 10080                                          
      goto 10402                                                        
10501 continue                                                          
      if(nlp .le. maxit)goto 10551                                      
      jerr=-m                                                           
      return                                                            
10551 continue                                                          
10070 continue                                                          
      iz = 1                                                            
      continue                                                          
10561 continue                                                          
      nlp=nlp+1                                                         
      dlx=0.0                                                           
      do 10571 l=1,nino                                                 
      k=ia(l)                                                           
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      gk=(dot_product(r(jx(jb:je)),x(jb:je))-rsum*xm(k))/xs(k)          
      ak=a(k)                                                           
      u=gk+ak*xv(k)                                                     
      au=abs(u)-vp(k)*ab                                                
      if(au .gt. 0.0)goto 10591                                         
      a(k)=0.0                                                          
      goto 10601                                                        
10591 continue                                                          
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*dem)))       
10601 continue                                                          
      continue                                                          
      if(a(k).eq.ak)goto 10571                                          
      d=a(k)-ak                                                         
      rsqc=rsqc+d*(2.0*gk-d*xv(k))                                      
      jb=ix(k)                                                          
      je=ix(k+1)-1                                                      
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)           
      r=r+d*v*xm(k)/xs(k)                                               
      rsum=sum(r)                                                       
      dlx=max(xv(k)*d**2,dlx)                                           
10571 continue                                                          
      continue                                                          
      d=0.0                                                             
      if(intr.ne.0) d=rsum/xmz                                          
      if(d .eq. 0.0)goto 10621                                          
      aint=aint+d                                                       
      rsqc=rsqc+d*(2.0*rsum-d*xmz)                                      
      dlx=max(dlx,xmz*d**2)                                             
      r=r-d*v                                                           
      rsum=rsum-d*xmz                                                   
10621 continue                                                          
      if(dlx.lt.thr)goto 10562                                          
      if(nlp .le. maxit)goto 10641                                      
      jerr=-m                                                           
      return                                                            
10641 continue                                                          
      goto 10561                                                        
10562 continue                                                          
      jz=0                                                              
      goto 10401                                                        
10402 continue                                                          
      deallocate(xv)                                                    
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
