c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))              
      subroutine get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)          772
      implicit double precision(a-h,o-z)                                    773
      data sml0,eps0,big0,mnlam0,rsqmax0,pmin0,exmx0  /1.0d-5,1.0d-6,9.9    775 
     *d35,5,0.999,1.0d-9,250.0/
      sml=sml0                                                              775
      eps=eps0                                                              775
      big=big0                                                              775
      mnlam=mnlam0                                                          775
      rsqmax=rsqmax0                                                        776
      pmin=pmin0                                                            776
      exmx=exmx0                                                            777
      return                                                                778
      entry chg_fract_dev(arg)                                              778
      sml0=arg                                                              778
      return                                                                779
      entry chg_dev_max(arg)                                                779
      rsqmax0=arg                                                           779
      return                                                                780
      entry chg_min_flmin(arg)                                              780
      eps0=arg                                                              780
      return                                                                781
      entry chg_big(arg)                                                    781
      big0=arg                                                              781
      return                                                                782
      entry chg_min_lambdas(irg)                                            782
      mnlam0=irg                                                            782
      return                                                                783
      entry chg_min_null_prob(arg)                                          783
      pmin0=arg                                                             783
      return                                                                784
      entry chg_max_exp(arg)                                                784
      exmx0=arg                                                             784
      return                                                                785
      end                                                                   786
      subroutine elnet  (ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,u    789 
     *lam,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                    790
      double precision x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam),cl(2,ni)     791
      double precision ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)              792
      integer jd(*),ia(nx),nin(nlam)                                        793
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 10021                                     796
      jerr=10000                                                            796
      return                                                                796
10021 continue                                                              797
      allocate(vq(1:ni),stat=jerr)                                          797
      if(jerr.ne.0) return                                                  798
      vq=max(0d0,vp)                                                        798
      vq=vq*ni/sum(vq)                                                      799
      if(ka .ne. 1)goto 10041                                               800
      call elnetu  (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,    803 
     *isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                            804
10041 continue                                                              805
      call elnetn (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,i    808 
     *sd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                              809
10031 continue                                                              809
      deallocate(vq)                                                        810
      return                                                                811
      end                                                                   812
      subroutine elnetu  (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ula    815 
     *m,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                    816
      double precision x(no,ni),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)      817
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)             818
      integer jd(*),ia(nx),nin(nlam)                                        819
      double precision, dimension (:), allocatable :: xm,xs,g,xv,vlam           
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           824
      if(jerr.ne.0) return                                                  825
      allocate(xm(1:ni),stat=jerr)                                          826
      if(jerr.ne.0) return                                                  827
      allocate(xs(1:ni),stat=jerr)                                          828
      if(jerr.ne.0) return                                                  829
      allocate(ju(1:ni),stat=jerr)                                          830
      if(jerr.ne.0) return                                                  831
      allocate(xv(1:ni),stat=jerr)                                          832
      if(jerr.ne.0) return                                                  833
      allocate(vlam(1:nlam),stat=jerr)                                      834
      if(jerr.ne.0) return                                                  835
      call chkvars(no,ni,x,ju)                                              836
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  837
      if(maxval(ju) .gt. 0)goto 10071                                       837
      jerr=7777                                                             837
      return                                                                837
10071 continue                                                              838
      call standard(no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr)          839
      if(jerr.ne.0) return                                                  840
      cl=cl/ys                                                              840
      if(isd .le. 0)goto 10091                                              840
10100 do 10101 j=1,ni                                                       840
      cl(:,j)=cl(:,j)*xs(j)                                                 840
10101 continue                                                              840
10102 continue                                                              840
10091 continue                                                              841
      if(flmin.ge.1.0) vlam=ulam/ys                                         842
      call elnet1(parm,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxi    844 
     *t,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  845
10110 do 10111 k=1,lmu                                                      845
      alm(k)=ys*alm(k)                                                      845
      nk=nin(k)                                                             846
10120 do 10121 l=1,nk                                                       846
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          846
10121 continue                                                              846
10122 continue                                                              846
      a0(k)=0.0                                                             847
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))           848
10111 continue                                                              849
10112 continue                                                              849
      deallocate(xm,xs,g,ju,xv,vlam)                                        850
      return                                                                851
      end                                                                   852
      subroutine standard (no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr    853 
     *)
      implicit double precision(a-h,o-z)                                    854
      double precision x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)      855
      integer ju(ni)                                                        856
      double precision, dimension (:), allocatable :: v                         
      allocate(v(1:no),stat=jerr)                                           859
      if(jerr.ne.0) return                                                  860
      w=w/sum(w)                                                            860
      v=sqrt(w)                                                             861
      if(intr .ne. 0)goto 10141                                             861
      ym=0.0                                                                861
      y=v*y                                                                 862
      ys=sqrt(dot_product(y,y)-dot_product(v,y)**2)                         862
      y=y/ys                                                                863
10150 do 10151 j=1,ni                                                       863
      if(ju(j).eq.0)goto 10151                                              863
      xm(j)=0.0                                                             863
      x(:,j)=v*x(:,j)                                                       864
      xv(j)=dot_product(x(:,j),x(:,j))                                      865
      if(isd .eq. 0)goto 10171                                              865
      xbq=dot_product(v,x(:,j))**2                                          865
      vc=xv(j)-xbq                                                          866
      xs(j)=sqrt(vc)                                                        866
      x(:,j)=x(:,j)/xs(j)                                                   866
      xv(j)=1.0+xbq/vc                                                      867
      goto 10181                                                            868
10171 continue                                                              868
      xs(j)=1.0                                                             868
10181 continue                                                              869
10161 continue                                                              869
10151 continue                                                              870
10152 continue                                                              870
      goto 10191                                                            871
10141 continue                                                              872
10200 do 10201 j=1,ni                                                       872
      if(ju(j).eq.0)goto 10201                                              873
      xm(j)=dot_product(w,x(:,j))                                           873
      x(:,j)=v*(x(:,j)-xm(j))                                               874
      xv(j)=dot_product(x(:,j),x(:,j))                                      874
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        875
10201 continue                                                              876
10202 continue                                                              876
      if(isd .ne. 0)goto 10221                                              876
      xs=1.0                                                                876
      goto 10231                                                            877
10221 continue                                                              878
10240 do 10241 j=1,ni                                                       878
      if(ju(j).eq.0)goto 10241                                              878
      x(:,j)=x(:,j)/xs(j)                                                   878
10241 continue                                                              879
10242 continue                                                              879
      xv=1.0                                                                880
10231 continue                                                              881
10211 continue                                                              881
      ym=dot_product(w,y)                                                   881
      y=v*(y-ym)                                                            881
      ys=sqrt(dot_product(y,y))                                             881
      y=y/ys                                                                882
10191 continue                                                              883
10131 continue                                                              883
      g=0.0                                                                 883
10250 do 10251 j=1,ni                                                       883
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                             883
10251 continue                                                              884
10252 continue                                                              884
      deallocate(v)                                                         885
      return                                                                886
      end                                                                   887
      subroutine elnet1 (beta,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,ulam,t    889 
     *hr,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                    890
      double precision vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam)         891
      double precision rsqo(nlam),almo(nlam),xv(ni)                         892
      double precision cl(2,ni)                                             893
      integer ju(ni),ia(nx),kin(nlam)                                       894
      double precision, dimension (:), allocatable :: a,da                      
      integer, dimension (:), allocatable :: mm                                 
      double precision, dimension (:,:), allocatable :: c                       
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      if(jerr.ne.0) return;                                                     
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)                902
      allocate(a(1:ni),stat=jerr)                                           903
      if(jerr.ne.0) return                                                  904
      allocate(mm(1:ni),stat=jerr)                                          905
      if(jerr.ne.0) return                                                  906
      allocate(da(1:ni),stat=jerr)                                          907
      if(jerr.ne.0) return                                                  908
      bta=beta                                                              908
      omb=1.0-bta                                                           909
      if(flmin .ge. 1.0)goto 10271                                          909
      eqs=max(eps,flmin)                                                    909
      alf=eqs**(1.0/(nlam-1))                                               909
10271 continue                                                              910
      rsq=0.0                                                               910
      a=0.0                                                                 910
      mm=0                                                                  910
      nlp=0                                                                 910
      nin=nlp                                                               910
      iz=0                                                                  910
      mnl=min(mnlam,nlam)                                                   912
      alm=0.0                                                               914
10280 do 10281 m=1,nlam                                                     915
      if(flmin .lt. 1.0)goto 10301                                          915
      alm=ulam(m)                                                           915
      goto 10291                                                            916
10301 if(m .le. 2)goto 10311                                                916
      alm=alm*alf                                                           916
      goto 10291                                                            917
10311 if(m .ne. 1)goto 10321                                                917
      alm=big                                                               917
      goto 10331                                                            918
10321 continue                                                              918
      alm=0.0                                                               919
10340 do 10341 j=1,ni                                                       919
      if(ju(j).eq.0)goto 10341                                              919
      if(vp(j).le.0.0)goto 10341                                            920
      alm=max(alm,abs(g(j))/vp(j))                                          921
10341 continue                                                              922
10342 continue                                                              922
      alm=alf*alm/max(bta,1.0d-3)                                           923
10331 continue                                                              924
10291 continue                                                              924
      dem=alm*omb                                                           924
      ab=alm*bta                                                            924
      rsq0=rsq                                                              924
      jz=1                                                                  925
10350 continue                                                              925
10351 continue                                                              925
      if(iz*jz.ne.0) go to 10360                                            925
      nlp=nlp+1                                                             925
      dlx=0.0                                                               926
10370 do 10371 k=1,ni                                                       926
      if(ju(k).eq.0)goto 10371                                              927
      ak=a(k)                                                               927
      u=g(k)+ak*xv(k)                                                       927
      v=abs(u)-vp(k)*ab                                                     927
      a(k)=0.0                                                              929
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d    930 
     *em)))
      if(a(k).eq.ak)goto 10371                                              931
      if(mm(k) .ne. 0)goto 10391                                            931
      nin=nin+1                                                             931
      if(nin.gt.nx)goto 10372                                               932
10400 do 10401 j=1,ni                                                       932
      if(ju(j).eq.0)goto 10401                                              933
      if(mm(j) .eq. 0)goto 10421                                            933
      c(j,nin)=c(k,mm(j))                                                   933
      goto 10401                                                            933
10421 continue                                                              934
      if(j .ne. k)goto 10441                                                934
      c(j,nin)=xv(j)                                                        934
      goto 10401                                                            934
10441 continue                                                              935
      c(j,nin)=dot_product(x(:,j),x(:,k))                                   936
10401 continue                                                              937
10402 continue                                                              937
      mm(k)=nin                                                             937
      ia(nin)=k                                                             938
10391 continue                                                              939
      del=a(k)-ak                                                           939
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      940
      dlx=max(xv(k)*del**2,dlx)                                             941
10450 do 10451 j=1,ni                                                       941
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                               941
10451 continue                                                              942
10452 continue                                                              942
10371 continue                                                              943
10372 continue                                                              943
      if(dlx.lt.thr)goto 10352                                              943
      if(nin.gt.nx)goto 10352                                               944
      if(nlp .le. maxit)goto 10471                                          944
      jerr=-m                                                               944
      return                                                                944
10471 continue                                                              945
10360 continue                                                              945
      iz=1                                                                  945
      da(1:nin)=a(ia(1:nin))                                                946
10480 continue                                                              946
10481 continue                                                              946
      nlp=nlp+1                                                             946
      dlx=0.0                                                               947
10490 do 10491 l=1,nin                                                      947
      k=ia(l)                                                               947
      ak=a(k)                                                               947
      u=g(k)+ak*xv(k)                                                       947
      v=abs(u)-vp(k)*ab                                                     948
      a(k)=0.0                                                              950
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d    951 
     *em)))
      if(a(k).eq.ak)goto 10491                                              952
      del=a(k)-ak                                                           952
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      953
      dlx=max(xv(k)*del**2,dlx)                                             954
10500 do 10501 j=1,nin                                                      954
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                  954
10501 continue                                                              955
10502 continue                                                              955
10491 continue                                                              956
10492 continue                                                              956
      if(dlx.lt.thr)goto 10482                                              956
      if(nlp .le. maxit)goto 10521                                          956
      jerr=-m                                                               956
      return                                                                956
10521 continue                                                              957
      goto 10481                                                            958
10482 continue                                                              958
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                      959
10530 do 10531 j=1,ni                                                       959
      if(mm(j).ne.0)goto 10531                                              960
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))            961
10531 continue                                                              962
10532 continue                                                              962
      jz=0                                                                  963
      goto 10351                                                            964
10352 continue                                                              964
      if(nin .le. nx)goto 10551                                             964
      jerr=-10000-m                                                         964
      goto 10282                                                            964
10551 continue                                                              965
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 965
      kin(m)=nin                                                            966
      rsqo(m)=rsq                                                           966
      almo(m)=alm                                                           966
      lmu=m                                                                 967
      if(m.lt.mnl)goto 10281                                                967
      if(flmin.ge.1.0)goto 10281                                            968
      me=0                                                                  968
10560 do 10561 j=1,nin                                                      968
      if(ao(j,m).ne.0.0) me=me+1                                            968
10561 continue                                                              968
10562 continue                                                              968
      if(me.gt.ne)goto 10282                                                969
      if(rsq-rsq0.lt.sml*rsq)goto 10282                                     969
      if(rsq.gt.rsqmax)goto 10282                                           970
10281 continue                                                              971
10282 continue                                                              971
      deallocate(a,mm,c,da)                                                 972
      return                                                                973
      end                                                                   974
      subroutine elnetn (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam    976 
     *,thr,isd,  intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                    977
      double precision vp(ni),x(no,ni),y(no),w(no),ulam(nlam),cl(2,ni)      978
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)             979
      integer jd(*),ia(nx),nin(nlam)                                        980
      double precision, dimension (:), allocatable :: xm,xs,xv,vlam             
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          985
      if(jerr.ne.0) return                                                  986
      allocate(xs(1:ni),stat=jerr)                                          987
      if(jerr.ne.0) return                                                  988
      allocate(ju(1:ni),stat=jerr)                                          989
      if(jerr.ne.0) return                                                  990
      allocate(xv(1:ni),stat=jerr)                                          991
      if(jerr.ne.0) return                                                  992
      allocate(vlam(1:nlam),stat=jerr)                                      993
      if(jerr.ne.0) return                                                  994
      call chkvars(no,ni,x,ju)                                              995
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  996
      if(maxval(ju) .gt. 0)goto 10581                                       996
      jerr=7777                                                             996
      return                                                                996
10581 continue                                                              997
      call standard1(no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)           998
      if(jerr.ne.0) return                                                  999
      cl=cl/ys                                                              999
      if(isd .le. 0)goto 10601                                              999
10610 do 10611 j=1,ni                                                       999
      cl(:,j)=cl(:,j)*xs(j)                                                 999
10611 continue                                                              999
10612 continue                                                              999
10601 continue                                                             1000
      if(flmin.ge.1.0) vlam=ulam/ys                                        1001
      call elnet2(parm,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,vlam,thr,maxi   1003 
     *t,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1004
10620 do 10621 k=1,lmu                                                     1004
      alm(k)=ys*alm(k)                                                     1004
      nk=nin(k)                                                            1005
10630 do 10631 l=1,nk                                                      1005
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1005
10631 continue                                                             1005
10632 continue                                                             1005
      a0(k)=0.0                                                            1006
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))          1007
10621 continue                                                             1008
10622 continue                                                             1008
      deallocate(xm,xs,ju,xv,vlam)                                         1009
      return                                                               1010
      end                                                                  1011
      subroutine standard1 (no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)   1012
      implicit double precision(a-h,o-z)                                   1013
      double precision x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)           1013
      integer ju(ni)                                                       1014
      double precision, dimension (:), allocatable :: v                         
      allocate(v(1:no),stat=jerr)                                          1017
      if(jerr.ne.0) return                                                 1018
      w=w/sum(w)                                                           1018
      v=sqrt(w)                                                            1019
      if(intr .ne. 0)goto 10651                                            1019
      ym=0.0                                                               1019
      y=v*y                                                                1020
      ys=sqrt(dot_product(y,y)-dot_product(v,y)**2)                        1020
      y=y/ys                                                               1021
10660 do 10661 j=1,ni                                                      1021
      if(ju(j).eq.0)goto 10661                                             1021
      xm(j)=0.0                                                            1021
      x(:,j)=v*x(:,j)                                                      1022
      xv(j)=dot_product(x(:,j),x(:,j))                                     1023
      if(isd .eq. 0)goto 10681                                             1023
      xbq=dot_product(v,x(:,j))**2                                         1023
      vc=xv(j)-xbq                                                         1024
      xs(j)=sqrt(vc)                                                       1024
      x(:,j)=x(:,j)/xs(j)                                                  1024
      xv(j)=1.0+xbq/vc                                                     1025
      goto 10691                                                           1026
10681 continue                                                             1026
      xs(j)=1.0                                                            1026
10691 continue                                                             1027
10671 continue                                                             1027
10661 continue                                                             1028
10662 continue                                                             1028
      go to 10700                                                          1029
10651 continue                                                             1030
10710 do 10711 j=1,ni                                                      1030
      if(ju(j).eq.0)goto 10711                                             1031
      xm(j)=dot_product(w,x(:,j))                                          1031
      x(:,j)=v*(x(:,j)-xm(j))                                              1032
      xv(j)=dot_product(x(:,j),x(:,j))                                     1032
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1033
10711 continue                                                             1034
10712 continue                                                             1034
      if(isd .ne. 0)goto 10731                                             1034
      xs=1.0                                                               1034
      goto 10741                                                           1035
10731 continue                                                             1035
10750 do 10751 j=1,ni                                                      1035
      if(ju(j).eq.0)goto 10751                                             1035
      x(:,j)=x(:,j)/xs(j)                                                  1035
10751 continue                                                             1036
10752 continue                                                             1036
      xv=1.0                                                               1037
10741 continue                                                             1038
10721 continue                                                             1038
      ym=dot_product(w,y)                                                  1038
      y=v*(y-ym)                                                           1038
      ys=sqrt(dot_product(y,y))                                            1038
      y=y/ys                                                               1039
10700 continue                                                             1039
      deallocate(v)                                                        1040
      return                                                               1041
      end                                                                  1042
      subroutine elnet2(beta,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,ulam,th   1044 
     *r,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1045
      double precision vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam)        1046
      double precision rsqo(nlam),almo(nlam),xv(ni)                        1047
      double precision cl(2,ni)                                            1048
      integer ju(ni),ia(nx),kin(nlam)                                      1049
      double precision, dimension (:), allocatable :: a,g                       
      integer, dimension (:), allocatable :: mm,ix                              
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1054
      allocate(a(1:ni),stat=jerr)                                          1055
      if(jerr.ne.0) return                                                 1056
      allocate(mm(1:ni),stat=jerr)                                         1057
      if(jerr.ne.0) return                                                 1058
      allocate(g(1:ni),stat=jerr)                                          1059
      if(jerr.ne.0) return                                                 1060
      allocate(ix(1:ni),stat=jerr)                                         1061
      if(jerr.ne.0) return                                                 1062
      bta=beta                                                             1062
      omb=1.0-bta                                                          1062
      ix=0                                                                 1064
      alf=1.0                                                              1066
      if(flmin .ge. 1.0)goto 10771                                         1066
      eqs=max(eps,flmin)                                                   1066
      alf=eqs**(1.0/(nlam-1))                                              1066
10771 continue                                                             1067
      rsq=0.0                                                              1067
      a=0.0                                                                1067
      mm=0                                                                 1067
      nlp=0                                                                1067
      nin=nlp                                                              1067
      iz=0                                                                 1067
      mnl=min(mnlam,nlam)                                                  1067
      alm=0.0                                                              1068
10780 do 10781 j=1,ni                                                      1068
      if(ju(j).eq.0)goto 10781                                             1068
      g(j)=abs(dot_product(y,x(:,j)))                                      1068
10781 continue                                                             1069
10782 continue                                                             1069
10790 do 10791 m=1,nlam                                                    1069
      alm0=alm                                                             1070
      if(flmin .lt. 1.0)goto 10811                                         1070
      alm=ulam(m)                                                          1070
      goto 10801                                                           1071
10811 if(m .le. 2)goto 10821                                               1071
      alm=alm*alf                                                          1071
      goto 10801                                                           1072
10821 if(m .ne. 1)goto 10831                                               1072
      alm=big                                                              1072
      goto 10841                                                           1073
10831 continue                                                             1073
      alm0=0.0                                                             1074
10850 do 10851 j=1,ni                                                      1074
      if(ju(j).eq.0)goto 10851                                             1074
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           1074
10851 continue                                                             1075
10852 continue                                                             1075
      alm0=alm0/max(bta,1.0d-3)                                            1075
      alm=alf*alm0                                                         1076
10841 continue                                                             1077
10801 continue                                                             1077
      dem=alm*omb                                                          1077
      ab=alm*bta                                                           1077
      rsq0=rsq                                                             1077
      jz=1                                                                 1078
      tlam=bta*(2.0*alm-alm0)                                              1079
10860 do 10861 k=1,ni                                                      1079
      if(ix(k).eq.1)goto 10861                                             1079
      if(ju(k).eq.0)goto 10861                                             1080
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       1081
10861 continue                                                             1082
10862 continue                                                             1082
10870 continue                                                             1082
10871 continue                                                             1082
      if(iz*jz.ne.0) go to 10360                                           1083
10880 continue                                                             1083
      nlp=nlp+1                                                            1083
      dlx=0.0                                                              1084
10890 do 10891 k=1,ni                                                      1084
      if(ix(k).eq.0)goto 10891                                             1084
      gk=dot_product(y,x(:,k))                                             1085
      ak=a(k)                                                              1085
      u=gk+ak*xv(k)                                                        1085
      v=abs(u)-vp(k)*ab                                                    1085
      a(k)=0.0                                                             1087
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1088 
     *em)))
      if(a(k).eq.ak)goto 10891                                             1089
      if(mm(k) .ne. 0)goto 10911                                           1089
      nin=nin+1                                                            1089
      if(nin.gt.nx)goto 10892                                              1090
      mm(k)=nin                                                            1090
      ia(nin)=k                                                            1091
10911 continue                                                             1092
      del=a(k)-ak                                                          1092
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1093
      y=y-del*x(:,k)                                                       1093
      dlx=max(xv(k)*del**2,dlx)                                            1094
10891 continue                                                             1095
10892 continue                                                             1095
      if(nin.gt.nx)goto 10872                                              1096
      if(dlx .ge. thr)goto 10931                                           1096
      ixx=0                                                                1097
10940 do 10941 k=1,ni                                                      1097
      if(ix(k).eq.1)goto 10941                                             1097
      if(ju(k).eq.0)goto 10941                                             1098
      g(k)=abs(dot_product(y,x(:,k)))                                      1099
      if(g(k) .le. ab*vp(k))goto 10961                                     1099
      ix(k)=1                                                              1099
      ixx=1                                                                1099
10961 continue                                                             1100
10941 continue                                                             1101
10942 continue                                                             1101
      if(ixx.eq.1) go to 10880                                             1102
      goto 10872                                                           1103
10931 continue                                                             1104
      if(nlp .le. maxit)goto 10981                                         1104
      jerr=-m                                                              1104
      return                                                               1104
10981 continue                                                             1105
10360 continue                                                             1105
      iz=1                                                                 1106
10990 continue                                                             1106
10991 continue                                                             1106
      nlp=nlp+1                                                            1106
      dlx=0.0                                                              1107
11000 do 11001 l=1,nin                                                     1107
      k=ia(l)                                                              1107
      gk=dot_product(y,x(:,k))                                             1108
      ak=a(k)                                                              1108
      u=gk+ak*xv(k)                                                        1108
      v=abs(u)-vp(k)*ab                                                    1108
      a(k)=0.0                                                             1110
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1111 
     *em)))
      if(a(k).eq.ak)goto 11001                                             1112
      del=a(k)-ak                                                          1112
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1113
      y=y-del*x(:,k)                                                       1113
      dlx=max(xv(k)*del**2,dlx)                                            1114
11001 continue                                                             1115
11002 continue                                                             1115
      if(dlx.lt.thr)goto 10992                                             1115
      if(nlp .le. maxit)goto 11021                                         1115
      jerr=-m                                                              1115
      return                                                               1115
11021 continue                                                             1116
      goto 10991                                                           1117
10992 continue                                                             1117
      jz=0                                                                 1118
      goto 10871                                                           1119
10872 continue                                                             1119
      if(nin .le. nx)goto 11041                                            1119
      jerr=-10000-m                                                        1119
      goto 10792                                                           1119
11041 continue                                                             1120
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1120
      kin(m)=nin                                                           1121
      rsqo(m)=rsq                                                          1121
      almo(m)=alm                                                          1121
      lmu=m                                                                1122
      if(m.lt.mnl)goto 10791                                               1122
      if(flmin.ge.1.0)goto 10791                                           1123
      me=0                                                                 1123
11050 do 11051 j=1,nin                                                     1123
      if(ao(j,m).ne.0.0) me=me+1                                           1123
11051 continue                                                             1123
11052 continue                                                             1123
      if(me.gt.ne)goto 10792                                               1124
      if(rsq-rsq0.lt.sml*rsq)goto 10792                                    1124
      if(rsq.gt.rsqmax)goto 10792                                          1125
10791 continue                                                             1126
10792 continue                                                             1126
      deallocate(a,mm,g,ix)                                                1127
      return                                                               1128
      end                                                                  1129
      subroutine chkvars(no,ni,x,ju)                                       1130
      implicit double precision(a-h,o-z)                                   1131
      double precision x(no,ni)                                            1131
      integer ju(ni)                                                       1132
11060 do 11061 j=1,ni                                                      1132
      ju(j)=0                                                              1132
      t=x(1,j)                                                             1133
11070 do 11071 i=2,no                                                      1133
      if(x(i,j).eq.t)goto 11071                                            1133
      ju(j)=1                                                              1133
      goto 11072                                                           1133
11071 continue                                                             1134
11072 continue                                                             1134
11061 continue                                                             1135
11062 continue                                                             1135
      return                                                               1136
      end                                                                  1137
      subroutine uncomp(ni,ca,ia,nin,a)                                    1138
      implicit double precision(a-h,o-z)                                   1139
      double precision ca(*),a(ni)                                         1139
      integer ia(*)                                                        1140
      a=0.0                                                                1140
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                  1141
      return                                                               1142
      end                                                                  1143
      subroutine modval(a0,ca,ia,nin,n,x,f)                                1144
      implicit double precision(a-h,o-z)                                   1145
      double precision ca(nin),x(n,*),f(n)                                 1145
      integer ia(nin)                                                      1146
      f=a0                                                                 1146
      if(nin.le.0) return                                                  1147
11080 do 11081 i=1,n                                                       1147
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      1147
11081 continue                                                             1148
11082 continue                                                             1148
      return                                                               1149
      end                                                                  1150
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam   1153 
     *,flmin,ulam,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      implicit double precision(a-h,o-z)                                   1154
      double precision x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)         1155
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)            1156
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1157
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 11101                                    1160
      jerr=10000                                                           1160
      return                                                               1160
11101 continue                                                             1161
      allocate(vq(1:ni),stat=jerr)                                         1161
      if(jerr.ne.0) return                                                 1162
      vq=max(0d0,vp)                                                       1162
      vq=vq*ni/sum(vq)                                                     1163
      if(ka .ne. 1)goto 11121                                              1164
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,u   1167 
     *lam,thr,isd,  intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 11131                                                           1168
11121 continue                                                             1169
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,ul   1172 
     *am,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
11131 continue                                                             1173
11111 continue                                                             1173
      deallocate(vq)                                                       1174
      return                                                               1175
      end                                                                  1176
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,f   1179 
     *lmin,ulam,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1180
      double precision x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)         1181
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)            1182
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1183
      double precision, dimension (:), allocatable :: xm,xs,g,xv,vlam           
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                          1188
      if(jerr.ne.0) return                                                 1189
      allocate(xm(1:ni),stat=jerr)                                         1190
      if(jerr.ne.0) return                                                 1191
      allocate(xs(1:ni),stat=jerr)                                         1192
      if(jerr.ne.0) return                                                 1193
      allocate(ju(1:ni),stat=jerr)                                         1194
      if(jerr.ne.0) return                                                 1195
      allocate(xv(1:ni),stat=jerr)                                         1196
      if(jerr.ne.0) return                                                 1197
      allocate(vlam(1:nlam),stat=jerr)                                     1198
      if(jerr.ne.0) return                                                 1199
      call spchkvars(no,ni,x,ix,ju)                                        1200
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1201
      if(maxval(ju) .gt. 0)goto 11151                                      1201
      jerr=7777                                                            1201
      return                                                               1201
11151 continue                                                             1202
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys,xv,jer   1203 
     *r)
      if(jerr.ne.0) return                                                 1204
      cl=cl/ys                                                             1204
      if(isd .le. 0)goto 11171                                             1204
11180 do 11181 j=1,ni                                                      1204
      cl(:,j)=cl(:,j)*xs(j)                                                1204
11181 continue                                                             1204
11182 continue                                                             1204
11171 continue                                                             1205
      if(flmin.ge.1.0) vlam=ulam/ys                                        1206
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla   1208 
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1209
11190 do 11191 k=1,lmu                                                     1209
      alm(k)=ys*alm(k)                                                     1209
      nk=nin(k)                                                            1210
11200 do 11201 l=1,nk                                                      1210
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1210
11201 continue                                                             1210
11202 continue                                                             1210
      a0(k)=0.0                                                            1211
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))          1212
11191 continue                                                             1213
11192 continue                                                             1213
      deallocate(xm,xs,g,ju,xv,vlam)                                       1214
      return                                                               1215
      end                                                                  1216
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys   1217 
     *,xv,jerr)
      implicit double precision(a-h,o-z)                                   1218
      double precision x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)         1219
      integer ix(*),jx(*),ju(ni)                                           1220
      w=w/sum(w)                                                           1221
      if(intr .ne. 0)goto 11221                                            1221
      ym=0.0                                                               1222
      ys=sqrt(dot_product(w,y**2)-dot_product(w,y)**2)                     1222
      y=y/ys                                                               1223
11230 do 11231 j=1,ni                                                      1223
      if(ju(j).eq.0)goto 11231                                             1223
      xm(j)=0.0                                                            1223
      jb=ix(j)                                                             1223
      je=ix(j+1)-1                                                         1224
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          1225
      if(isd .eq. 0)goto 11251                                             1225
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            1225
      vc=xv(j)-xbq                                                         1226
      xs(j)=sqrt(vc)                                                       1226
      xv(j)=1.0+xbq/vc                                                     1227
      goto 11261                                                           1228
11251 continue                                                             1228
      xs(j)=1.0                                                            1228
11261 continue                                                             1229
11241 continue                                                             1229
11231 continue                                                             1230
11232 continue                                                             1230
      goto 11271                                                           1231
11221 continue                                                             1232
11280 do 11281 j=1,ni                                                      1232
      if(ju(j).eq.0)goto 11281                                             1233
      jb=ix(j)                                                             1233
      je=ix(j+1)-1                                                         1233
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1234
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1235
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1236
11281 continue                                                             1237
11282 continue                                                             1237
      if(isd .ne. 0)goto 11301                                             1237
      xs=1.0                                                               1237
      goto 11311                                                           1237
11301 continue                                                             1237
      xv=1.0                                                               1237
11311 continue                                                             1238
11291 continue                                                             1238
      ym=dot_product(w,y)                                                  1238
      y=y-ym                                                               1238
      ys=sqrt(dot_product(w,y**2))                                         1238
      y=y/ys                                                               1239
11271 continue                                                             1240
11211 continue                                                             1240
      g=0.0                                                                1241
11320 do 11321 j=1,ni                                                      1241
      if(ju(j).eq.0)goto 11321                                             1241
      jb=ix(j)                                                             1241
      je=ix(j+1)-1                                                         1242
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)           1243
11321 continue                                                             1244
11322 continue                                                             1244
      return                                                               1245
      end                                                                  1246
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm   1248 
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1249
      double precision g(ni),vp(ni),x(*),ulam(nlam),w(no)                  1250
      double precision ao(nx,nlam),rsqo(nlam),almo(nlam)                   1251
      double precision xm(ni),xs(ni),xv(ni),cl(2,ni)                       1252
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1253
      double precision, dimension (:), allocatable :: a,da                      
      integer, dimension (:), allocatable :: mm                                 
      double precision, dimension (:,:), allocatable :: c                       
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      if(jerr.ne.0) return;                                                     
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1261
      allocate(a(1:ni),stat=jerr)                                          1262
      if(jerr.ne.0) return                                                 1263
      allocate(mm(1:ni),stat=jerr)                                         1264
      if(jerr.ne.0) return                                                 1265
      allocate(da(1:ni),stat=jerr)                                         1266
      if(jerr.ne.0) return                                                 1267
      bta=beta                                                             1267
      omb=1.0-bta                                                          1269
      alm=0.0                                                              1269
      alf=1.0                                                              1271
      if(flmin .ge. 1.0)goto 11341                                         1271
      eqs=max(eps,flmin)                                                   1271
      alf=eqs**(1.0/(nlam-1))                                              1271
11341 continue                                                             1272
      rsq=0.0                                                              1272
      a=0.0                                                                1272
      mm=0                                                                 1272
      nlp=0                                                                1272
      nin=nlp                                                              1272
      iz=0                                                                 1272
      mnl=min(mnlam,nlam)                                                  1273
11350 do 11351 m=1,nlam                                                    1274
      if(flmin .lt. 1.0)goto 11371                                         1274
      alm=ulam(m)                                                          1274
      goto 11361                                                           1275
11371 if(m .le. 2)goto 11381                                               1275
      alm=alm*alf                                                          1275
      goto 11361                                                           1276
11381 if(m .ne. 1)goto 11391                                               1276
      alm=big                                                              1276
      goto 11401                                                           1277
11391 continue                                                             1277
      alm=0.0                                                              1278
11410 do 11411 j=1,ni                                                      1278
      if(ju(j).eq.0)goto 11411                                             1278
      if(vp(j).le.0.0)goto 11411                                           1279
      alm=max(alm,abs(g(j))/vp(j))                                         1280
11411 continue                                                             1281
11412 continue                                                             1281
      alm=alf*alm/max(bta,1.0d-3)                                          1282
11401 continue                                                             1283
11361 continue                                                             1283
      dem=alm*omb                                                          1283
      ab=alm*bta                                                           1283
      rsq0=rsq                                                             1283
      jz=1                                                                 1284
11420 continue                                                             1284
11421 continue                                                             1284
      if(iz*jz.ne.0) go to 10360                                           1284
      nlp=nlp+1                                                            1284
      dlx=0.0                                                              1285
11430 do 11431 k=1,ni                                                      1285
      if(ju(k).eq.0)goto 11431                                             1286
      ak=a(k)                                                              1286
      u=g(k)+ak*xv(k)                                                      1286
      v=abs(u)-vp(k)*ab                                                    1286
      a(k)=0.0                                                             1288
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1289 
     *em)))
      if(a(k).eq.ak)goto 11431                                             1290
      if(mm(k) .ne. 0)goto 11451                                           1290
      nin=nin+1                                                            1290
      if(nin.gt.nx)goto 11432                                              1291
11460 do 11461 j=1,ni                                                      1291
      if(ju(j).eq.0)goto 11461                                             1292
      if(mm(j) .eq. 0)goto 11481                                           1292
      c(j,nin)=c(k,mm(j))                                                  1292
      goto 11461                                                           1292
11481 continue                                                             1293
      if(j .ne. k)goto 11501                                               1293
      c(j,nin)=xv(j)                                                       1293
      goto 11461                                                           1293
11501 continue                                                             1294
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))       1296
11461 continue                                                             1297
11462 continue                                                             1297
      mm(k)=nin                                                            1297
      ia(nin)=k                                                            1298
11451 continue                                                             1299
      del=a(k)-ak                                                          1299
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1300
      dlx=max(xv(k)*del**2,dlx)                                            1301
11510 do 11511 j=1,ni                                                      1301
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                              1301
11511 continue                                                             1302
11512 continue                                                             1302
11431 continue                                                             1303
11432 continue                                                             1303
      if(dlx.lt.thr)goto 11422                                             1303
      if(nin.gt.nx)goto 11422                                              1304
      if(nlp .le. maxit)goto 11531                                         1304
      jerr=-m                                                              1304
      return                                                               1304
11531 continue                                                             1305
10360 continue                                                             1305
      iz=1                                                                 1305
      da(1:nin)=a(ia(1:nin))                                               1306
11540 continue                                                             1306
11541 continue                                                             1306
      nlp=nlp+1                                                            1306
      dlx=0.0                                                              1307
11550 do 11551 l=1,nin                                                     1307
      k=ia(l)                                                              1308
      ak=a(k)                                                              1308
      u=g(k)+ak*xv(k)                                                      1308
      v=abs(u)-vp(k)*ab                                                    1308
      a(k)=0.0                                                             1310
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1311 
     *em)))
      if(a(k).eq.ak)goto 11551                                             1312
      del=a(k)-ak                                                          1312
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     1313
      dlx=max(xv(k)*del**2,dlx)                                            1314
11560 do 11561 j=1,nin                                                     1314
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                 1314
11561 continue                                                             1315
11562 continue                                                             1315
11551 continue                                                             1316
11552 continue                                                             1316
      if(dlx.lt.thr)goto 11542                                             1316
      if(nlp .le. maxit)goto 11581                                         1316
      jerr=-m                                                              1316
      return                                                               1316
11581 continue                                                             1317
      goto 11541                                                           1318
11542 continue                                                             1318
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                     1319
11590 do 11591 j=1,ni                                                      1319
      if(mm(j).ne.0)goto 11591                                             1320
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))           1321
11591 continue                                                             1322
11592 continue                                                             1322
      jz=0                                                                 1323
      goto 11421                                                           1324
11422 continue                                                             1324
      if(nin .le. nx)goto 11611                                            1324
      jerr=-10000-m                                                        1324
      goto 11352                                                           1324
11611 continue                                                             1325
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1325
      kin(m)=nin                                                           1326
      rsqo(m)=rsq                                                          1326
      almo(m)=alm                                                          1326
      lmu=m                                                                1327
      if(m.lt.mnl)goto 11351                                               1327
      if(flmin.ge.1.0)goto 11351                                           1328
      me=0                                                                 1328
11620 do 11621 j=1,nin                                                     1328
      if(ao(j,m).ne.0.0) me=me+1                                           1328
11621 continue                                                             1328
11622 continue                                                             1328
      if(me.gt.ne)goto 11352                                               1329
      if(rsq-rsq0.lt.sml*rsq)goto 11352                                    1329
      if(rsq.gt.rsqmax)goto 11352                                          1330
11351 continue                                                             1331
11352 continue                                                             1331
      deallocate(a,mm,c,da)                                                1332
      return                                                               1333
      end                                                                  1334
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flm   1336 
     *in,ulam,  thr,isd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1337
      double precision x(*),vp(ni),y(no),w(no),ulam(nlam),cl(2,ni)         1338
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)            1339
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           1340
      double precision, dimension (:), allocatable :: xm,xs,xv,vlam             
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         1345
      if(jerr.ne.0) return                                                 1346
      allocate(xs(1:ni),stat=jerr)                                         1347
      if(jerr.ne.0) return                                                 1348
      allocate(ju(1:ni),stat=jerr)                                         1349
      if(jerr.ne.0) return                                                 1350
      allocate(xv(1:ni),stat=jerr)                                         1351
      if(jerr.ne.0) return                                                 1352
      allocate(vlam(1:nlam),stat=jerr)                                     1353
      if(jerr.ne.0) return                                                 1354
      call spchkvars(no,ni,x,ix,ju)                                        1355
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1356
      if(maxval(ju) .gt. 0)goto 11641                                      1356
      jerr=7777                                                            1356
      return                                                               1356
11641 continue                                                             1357
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,xv,jerr   1358 
     *)
      if(jerr.ne.0) return                                                 1359
      cl=cl/ys                                                             1359
      if(isd .le. 0)goto 11661                                             1359
11670 do 11671 j=1,ni                                                      1359
      cl(:,j)=cl(:,j)*xs(j)                                                1359
11671 continue                                                             1359
11672 continue                                                             1359
11661 continue                                                             1360
      if(flmin.ge.1.0) vlam=ulam/ys                                        1361
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla   1363 
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 1364
11680 do 11681 k=1,lmu                                                     1364
      alm(k)=ys*alm(k)                                                     1364
      nk=nin(k)                                                            1365
11690 do 11691 l=1,nk                                                      1365
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         1365
11691 continue                                                             1365
11692 continue                                                             1365
      a0(k)=0.0                                                            1366
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))          1367
11681 continue                                                             1368
11682 continue                                                             1368
      deallocate(xm,xs,ju,xv,vlam)                                         1369
      return                                                               1370
      end                                                                  1371
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,   1372 
     *xv,jerr)
      implicit double precision(a-h,o-z)                                   1373
      double precision x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)               1374
      integer ix(*),jx(*),ju(ni)                                           1375
      w=w/sum(w)                                                           1376
      if(intr .ne. 0)goto 11711                                            1376
      ym=0.0                                                               1377
      ys=sqrt(dot_product(w,y**2)-dot_product(w,y)**2)                     1377
      y=y/ys                                                               1378
11720 do 11721 j=1,ni                                                      1378
      if(ju(j).eq.0)goto 11721                                             1378
      xm(j)=0.0                                                            1378
      jb=ix(j)                                                             1378
      je=ix(j+1)-1                                                         1379
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          1380
      if(isd .eq. 0)goto 11741                                             1380
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            1380
      vc=xv(j)-xbq                                                         1381
      xs(j)=sqrt(vc)                                                       1381
      xv(j)=1.0+xbq/vc                                                     1382
      goto 11751                                                           1383
11741 continue                                                             1383
      xs(j)=1.0                                                            1383
11751 continue                                                             1384
11731 continue                                                             1384
11721 continue                                                             1385
11722 continue                                                             1385
      return                                                               1386
11711 continue                                                             1387
11760 do 11761 j=1,ni                                                      1387
      if(ju(j).eq.0)goto 11761                                             1388
      jb=ix(j)                                                             1388
      je=ix(j+1)-1                                                         1388
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             1389
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 1390
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       1391
11761 continue                                                             1392
11762 continue                                                             1392
      if(isd .ne. 0)goto 11781                                             1392
      xs=1.0                                                               1392
      goto 11791                                                           1392
11781 continue                                                             1392
      xv=1.0                                                               1392
11791 continue                                                             1393
11771 continue                                                             1393
      ym=dot_product(w,y)                                                  1393
      y=y-ym                                                               1393
      ys=sqrt(dot_product(w,y**2))                                         1393
      y=y/ys                                                               1394
      return                                                               1395
      end                                                                  1396
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm   1398 
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   1399
      double precision y(no),w(no),x(*),vp(ni),ulam(nlam),cl(2,ni)         1400
      double precision ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),x   1401 
     *v(ni)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          1402
      double precision, dimension (:), allocatable :: a,g                       
      integer, dimension (:), allocatable :: mm,iy                              
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               1407
      allocate(a(1:ni),stat=jerr)                                          1408
      if(jerr.ne.0) return                                                 1409
      allocate(mm(1:ni),stat=jerr)                                         1410
      if(jerr.ne.0) return                                                 1411
      allocate(g(1:ni),stat=jerr)                                          1412
      if(jerr.ne.0) return                                                 1413
      allocate(iy(1:ni),stat=jerr)                                         1414
      if(jerr.ne.0) return                                                 1415
      bta=beta                                                             1415
      omb=1.0-bta                                                          1415
      alm=0.0                                                              1415
      iy=0                                                                 1417
      alf=1.0                                                              1419
      if(flmin .ge. 1.0)goto 11811                                         1419
      eqs=max(eps,flmin)                                                   1419
      alf=eqs**(1.0/(nlam-1))                                              1419
11811 continue                                                             1420
      rsq=0.0                                                              1420
      a=0.0                                                                1420
      mm=0                                                                 1420
      o=0.0                                                                1420
      nlp=0                                                                1420
      nin=nlp                                                              1420
      iz=0                                                                 1420
      mnl=min(mnlam,nlam)                                                  1421
11820 do 11821 j=1,ni                                                      1421
      if(ju(j).eq.0)goto 11821                                             1422
      jb=ix(j)                                                             1422
      je=ix(j+1)-1                                                         1423
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1424
11821 continue                                                             1425
11822 continue                                                             1425
11830 do 11831 m=1,nlam                                                    1425
      alm0=alm                                                             1426
      if(flmin .lt. 1.0)goto 11851                                         1426
      alm=ulam(m)                                                          1426
      goto 11841                                                           1427
11851 if(m .le. 2)goto 11861                                               1427
      alm=alm*alf                                                          1427
      goto 11841                                                           1428
11861 if(m .ne. 1)goto 11871                                               1428
      alm=big                                                              1428
      goto 11881                                                           1429
11871 continue                                                             1429
      alm0=0.0                                                             1430
11890 do 11891 j=1,ni                                                      1430
      if(ju(j).eq.0)goto 11891                                             1430
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           1430
11891 continue                                                             1431
11892 continue                                                             1431
      alm0=alm0/max(bta,1.0d-3)                                            1431
      alm=alf*alm0                                                         1432
11881 continue                                                             1433
11841 continue                                                             1433
      dem=alm*omb                                                          1433
      ab=alm*bta                                                           1433
      rsq0=rsq                                                             1433
      jz=1                                                                 1434
      tlam=bta*(2.0*alm-alm0)                                              1435
11900 do 11901 k=1,ni                                                      1435
      if(iy(k).eq.1)goto 11901                                             1435
      if(ju(k).eq.0)goto 11901                                             1436
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       1437
11901 continue                                                             1438
11902 continue                                                             1438
11910 continue                                                             1438
11911 continue                                                             1438
      if(iz*jz.ne.0) go to 10360                                           1439
10880 continue                                                             1439
      nlp=nlp+1                                                            1439
      dlx=0.0                                                              1440
11920 do 11921 k=1,ni                                                      1440
      if(iy(k).eq.0)goto 11921                                             1440
      jb=ix(k)                                                             1440
      je=ix(k+1)-1                                                         1441
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1442
      ak=a(k)                                                              1442
      u=gk+ak*xv(k)                                                        1442
      v=abs(u)-vp(k)*ab                                                    1442
      a(k)=0.0                                                             1444
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1445 
     *em)))
      if(a(k).eq.ak)goto 11921                                             1446
      if(mm(k) .ne. 0)goto 11941                                           1446
      nin=nin+1                                                            1446
      if(nin.gt.nx)goto 11922                                              1447
      mm(k)=nin                                                            1447
      ia(nin)=k                                                            1448
11941 continue                                                             1449
      del=a(k)-ak                                                          1449
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1450
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1451
      o=o+del*xm(k)/xs(k)                                                  1451
      dlx=max(xv(k)*del**2,dlx)                                            1452
11921 continue                                                             1453
11922 continue                                                             1453
      if(nin.gt.nx)goto 11912                                              1454
      if(dlx .ge. thr)goto 11961                                           1454
      ixx=0                                                                1455
11970 do 11971 j=1,ni                                                      1455
      if(iy(j).eq.1)goto 11971                                             1455
      if(ju(j).eq.0)goto 11971                                             1456
      jb=ix(j)                                                             1456
      je=ix(j+1)-1                                                         1457
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    1458
      if(g(j) .le. ab*vp(j))goto 11991                                     1458
      iy(j)=1                                                              1458
      ixx=1                                                                1458
11991 continue                                                             1459
11971 continue                                                             1460
11972 continue                                                             1460
      if(ixx.eq.1) go to 10880                                             1461
      goto 11912                                                           1462
11961 continue                                                             1463
      if(nlp .le. maxit)goto 12011                                         1463
      jerr=-m                                                              1463
      return                                                               1463
12011 continue                                                             1464
10360 continue                                                             1464
      iz=1                                                                 1465
12020 continue                                                             1465
12021 continue                                                             1465
      nlp=nlp+1                                                            1465
      dlx=0.0                                                              1466
12030 do 12031 l=1,nin                                                     1466
      k=ia(l)                                                              1466
      jb=ix(k)                                                             1466
      je=ix(k+1)-1                                                         1467
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           1468
      ak=a(k)                                                              1468
      u=gk+ak*xv(k)                                                        1468
      v=abs(u)-vp(k)*ab                                                    1468
      a(k)=0.0                                                             1470
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   1471 
     *em)))
      if(a(k).eq.ak)goto 12031                                             1472
      del=a(k)-ak                                                          1472
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       1473
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         1474
      o=o+del*xm(k)/xs(k)                                                  1474
      dlx=max(xv(k)*del**2,dlx)                                            1475
12031 continue                                                             1476
12032 continue                                                             1476
      if(dlx.lt.thr)goto 12022                                             1476
      if(nlp .le. maxit)goto 12051                                         1476
      jerr=-m                                                              1476
      return                                                               1476
12051 continue                                                             1477
      goto 12021                                                           1478
12022 continue                                                             1478
      jz=0                                                                 1479
      goto 11911                                                           1480
11912 continue                                                             1480
      if(nin .le. nx)goto 12071                                            1480
      jerr=-10000-m                                                        1480
      goto 11832                                                           1480
12071 continue                                                             1481
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                1481
      kin(m)=nin                                                           1482
      rsqo(m)=rsq                                                          1482
      almo(m)=alm                                                          1482
      lmu=m                                                                1483
      if(m.lt.mnl)goto 11831                                               1483
      if(flmin.ge.1.0)goto 11831                                           1484
      me=0                                                                 1484
12080 do 12081 j=1,nin                                                     1484
      if(ao(j,m).ne.0.0) me=me+1                                           1484
12081 continue                                                             1484
12082 continue                                                             1484
      if(me.gt.ne)goto 11832                                               1485
      if(rsq-rsq0.lt.sml*rsq)goto 11832                                    1485
      if(rsq.gt.rsqmax)goto 11832                                          1486
11831 continue                                                             1487
11832 continue                                                             1487
      deallocate(a,mm,g,iy)                                                1488
      return                                                               1489
      end                                                                  1490
      subroutine spchkvars(no,ni,x,ix,ju)                                  1491
      implicit double precision(a-h,o-z)                                   1492
      double precision x(*)                                                1492
      integer ix(*),ju(ni)                                                 1493
12090 do 12091 j=1,ni                                                      1493
      ju(j)=0                                                              1493
      jb=ix(j)                                                             1493
      nj=ix(j+1)-jb                                                        1493
      if(nj.eq.0)goto 12091                                                1494
      je=ix(j+1)-1                                                         1495
      if(nj .ge. no)goto 12111                                             1495
12120 do 12121 i=jb,je                                                     1495
      if(x(i).eq.0.0)goto 12121                                            1495
      ju(j)=1                                                              1495
      goto 12122                                                           1495
12121 continue                                                             1495
12122 continue                                                             1495
      goto 12131                                                           1496
12111 continue                                                             1496
      t=x(jb)                                                              1496
12140 do 12141 i=jb+1,je                                                   1496
      if(x(i).eq.t)goto 12141                                              1496
      ju(j)=1                                                              1496
      goto 12142                                                           1496
12141 continue                                                             1496
12142 continue                                                             1496
12131 continue                                                             1497
12101 continue                                                             1497
12091 continue                                                             1498
12092 continue                                                             1498
      return                                                               1499
      end                                                                  1500
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                         1501
      implicit double precision(a-h,o-z)                                   1502
      double precision ca(*),x(*),f(n)                                     1502
      integer ia(*),ix(*),jx(*)                                            1503
      f=a0                                                                 1504
12150 do 12151 j=1,nin                                                     1504
      k=ia(j)                                                              1504
      kb=ix(k)                                                             1504
      ke=ix(k+1)-1                                                         1505
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                             1506
12151 continue                                                             1507
12152 continue                                                             1507
      return                                                               1508
      end                                                                  1509
      function row_prod(i,j,ia,ja,ra,w)                                    1510
      implicit double precision(a-h,o-z)                                   1511
      integer ia(*),ja(*)                                                  1511
      double precision ra(*),w(*)                                          1512
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(   1514 
     *i),ia(j+1)-ia(j),w)
      return                                                               1515
      end                                                                  1516
      function dot(x,y,mx,my,nx,ny,w)                                      1517
      implicit double precision(a-h,o-z)                                   1518
      double precision x(*),y(*),w(*)                                      1518
      integer mx(*),my(*)                                                  1519
      i=1                                                                  1519
      j=i                                                                  1519
      s=0.0                                                                1520
12160 continue                                                             1520
12161 continue                                                             1520
12170 continue                                                             1521
12171 if(mx(i).ge.my(j))goto 12172                                         1521
      i=i+1                                                                1521
      if(i.gt.nx) go to 12180                                              1521
      goto 12171                                                           1522
12172 continue                                                             1522
      if(mx(i).eq.my(j)) go to 12190                                       1523
12200 continue                                                             1523
12201 if(my(j).ge.mx(i))goto 12202                                         1523
      j=j+1                                                                1523
      if(j.gt.ny) go to 12180                                              1523
      goto 12201                                                           1524
12202 continue                                                             1524
      if(mx(i).eq.my(j)) go to 12190                                       1524
      goto 12161                                                           1525
12190 continue                                                             1525
      s=s+w(mx(i))*x(i)*y(j)                                               1526
      i=i+1                                                                1526
      if(i.gt.nx)goto 12162                                                1526
      j=j+1                                                                1526
      if(j.gt.ny)goto 12162                                                1527
      goto 12161                                                           1528
12162 continue                                                             1528
12180 continue                                                             1528
      dot=s                                                                1529
      return                                                               1530
      end                                                                  1531
      subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,cl,ne,nx,nlam,flmin,u   1533 
     *lam,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      implicit double precision(a-h,o-z)                                   1534
      double precision x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nla   1535 
     *m)
      double precision ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl   1536 
     *(2,ni)
      integer jd(*),ia(nx),nin(nlam)                                       1537
      double precision, dimension (:), allocatable :: xm,xs,ww,vq,xv            
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 12221                                    1541
      jerr=10000                                                           1541
      return                                                               1541
12221 continue                                                             1542
      allocate(ww(1:no),stat=jerr)                                         1543
      if(jerr.ne.0) return                                                 1544
      allocate(ju(1:ni),stat=jerr)                                         1545
      if(jerr.ne.0) return                                                 1546
      allocate(vq(1:ni),stat=jerr)                                         1547
      if(jerr.ne.0) return                                                 1548
      allocate(xm(1:ni),stat=jerr)                                         1549
      if(jerr.ne.0) return                                                 1550
      if(kopt .ne. 2)goto 12241                                            1550
      allocate(xv(1:ni),stat=jerr)                                         1550
      if(jerr.ne.0) return                                                 1550
12241 continue                                                             1551
      if(isd .le. 0)goto 12261                                             1551
      allocate(xs(1:ni),stat=jerr)                                         1551
      if(jerr.ne.0) return                                                 1551
12261 continue                                                             1553
      call chkvars(no,ni,x,ju)                                             1554
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1555
      if(maxval(ju) .gt. 0)goto 12281                                      1555
      jerr=7777                                                            1555
      return                                                               1555
12281 continue                                                             1556
      vq=max(0d0,vp)                                                       1556
      vq=vq*ni/sum(vq)                                                     1557
12290 do 12291 i=1,no                                                      1557
      ww(i)=sum(y(i,:))                                                    1557
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                                 1557
12291 continue                                                             1558
12292 continue                                                             1558
      sw=sum(ww)                                                           1558
      ww=ww/sw                                                             1559
      if(nc .ne. 1)goto 12311                                              1559
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        1560
      if(isd .le. 0)goto 12331                                             1560
12340 do 12341 j=1,ni                                                      1560
      cl(:,j)=cl(:,j)*xs(j)                                                1560
12341 continue                                                             1560
12342 continue                                                             1560
12331 continue                                                             1561
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,cl,ne,nx,nlam,fl   1563 
     *min,ulam,  thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,n
     *lp,jerr)
      goto 12301                                                           1564
12311 if(kopt .ne. 2)goto 12351                                            1564
      call multlstandard1(no,ni,x,ww,ju,isd,intr,xm,xs,xv)                 1565
      if(isd .le. 0)goto 12371                                             1565
12380 do 12381 j=1,ni                                                      1565
      cl(:,j)=cl(:,j)*xs(j)                                                1565
12381 continue                                                             1565
12382 continue                                                             1565
12371 continue                                                             1566
      call multlognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,   1568 
     *ulam,thr,  intr,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 12391                                                           1569
12351 continue                                                             1569
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        1570
      if(isd .le. 0)goto 12411                                             1570
12420 do 12421 j=1,ni                                                      1570
      cl(:,j)=cl(:,j)*xs(j)                                                1570
12421 continue                                                             1570
12422 continue                                                             1570
12411 continue                                                             1571
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam   1573 
     *,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
12391 continue                                                             1574
12301 continue                                                             1574
      if(jerr.gt.0) return                                                 1574
      dev0=2.0*sw*dev0                                                     1575
12430 do 12431 k=1,lmu                                                     1575
      nk=nin(k)                                                            1576
12440 do 12441 ic=1,nc                                                     1576
      if(isd .le. 0)goto 12461                                             1576
12470 do 12471 l=1,nk                                                      1576
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1576
12471 continue                                                             1576
12472 continue                                                             1576
12461 continue                                                             1577
      if(intr .ne. 0)goto 12491                                            1577
      a0(ic,k)=0.0                                                         1577
      goto 12501                                                           1578
12491 continue                                                             1578
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1578
12501 continue                                                             1579
12481 continue                                                             1579
12441 continue                                                             1580
12442 continue                                                             1580
12431 continue                                                             1581
12432 continue                                                             1581
      deallocate(ww,ju,vq,xm)                                              1581
      if(isd.gt.0) deallocate(xs)                                          1582
      if(kopt.eq.2) deallocate(xv)                                         1583
      return                                                               1584
      end                                                                  1585
      subroutine lstandard1 (no,ni,x,w,ju,isd,intr,xm,xs)                  1586
      implicit double precision(a-h,o-z)                                   1587
      double precision x(no,ni),w(no),xm(ni),xs(ni)                        1587
      integer ju(ni)                                                       1588
      if(intr .ne. 0)goto 12521                                            1589
12530 do 12531 j=1,ni                                                      1589
      if(ju(j).eq.0)goto 12531                                             1589
      xm(j)=0.0                                                            1590
      if(isd .eq. 0)goto 12551                                             1590
      vc=dot_product(w,x(:,j)**2)-dot_product(w,x(:,j))**2                 1591
      xs(j)=sqrt(vc)                                                       1591
      x(:,j)=x(:,j)/xs(j)                                                  1592
12551 continue                                                             1593
12531 continue                                                             1594
12532 continue                                                             1594
      return                                                               1595
12521 continue                                                             1596
12560 do 12561 j=1,ni                                                      1596
      if(ju(j).eq.0)goto 12561                                             1597
      xm(j)=dot_product(w,x(:,j))                                          1597
      x(:,j)=x(:,j)-xm(j)                                                  1598
      if(isd .le. 0)goto 12581                                             1598
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 1598
      x(:,j)=x(:,j)/xs(j)                                                  1598
12581 continue                                                             1599
12561 continue                                                             1600
12562 continue                                                             1600
      return                                                               1601
      end                                                                  1602
      subroutine multlstandard1 (no,ni,x,w,ju,isd,intr,xm,xs,xv)           1603
      implicit double precision(a-h,o-z)                                   1604
      double precision x(no,ni),w(no),xm(ni),xs(ni),xv(ni)                 1604
      integer ju(ni)                                                       1605
      if(intr .ne. 0)goto 12601                                            1606
12610 do 12611 j=1,ni                                                      1606
      if(ju(j).eq.0)goto 12611                                             1606
      xm(j)=0.0                                                            1607
      xv(j)=dot_product(w,x(:,j)**2)                                       1608
      if(isd .eq. 0)goto 12631                                             1608
      xbq=dot_product(w,x(:,j))**2                                         1608
      vc=xv(j)-xbq                                                         1609
      xs(j)=sqrt(vc)                                                       1609
      x(:,j)=x(:,j)/xs(j)                                                  1609
      xv(j)=1.0+xbq/vc                                                     1610
12631 continue                                                             1611
12611 continue                                                             1612
12612 continue                                                             1612
      return                                                               1613
12601 continue                                                             1614
12640 do 12641 j=1,ni                                                      1614
      if(ju(j).eq.0)goto 12641                                             1615
      xm(j)=dot_product(w,x(:,j))                                          1615
      x(:,j)=x(:,j)-xm(j)                                                  1616
      xv(j)=dot_product(w,x(:,j)**2)                                       1617
      if(isd .le. 0)goto 12661                                             1617
      xs(j)=sqrt(xv(j))                                                    1617
      x(:,j)=x(:,j)/xs(j)                                                  1617
      xv(j)=1.0                                                            1617
12661 continue                                                             1618
12641 continue                                                             1619
12642 continue                                                             1619
      return                                                               1620
      end                                                                  1621
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin,u   1623 
     *lam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                   1624
      double precision x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2   1625 
     *,ni)
      double precision a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)             1626
      integer ju(ni),m(nx),kin(nlam)                                       1627
      double precision, dimension (:), allocatable :: b,bs,v,r,xv,q,ga          
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               1632
      allocate(b(0:ni),stat=jerr)                                          1633
      if(jerr.ne.0) return                                                 1634
      allocate(xv(1:ni),stat=jerr)                                         1635
      if(jerr.ne.0) return                                                 1636
      allocate(ga(1:ni),stat=jerr)                                         1637
      if(jerr.ne.0) return                                                 1638
      allocate(bs(0:ni),stat=jerr)                                         1639
      if(jerr.ne.0) return                                                 1640
      allocate(mm(1:ni),stat=jerr)                                         1641
      if(jerr.ne.0) return                                                 1642
      allocate(ixx(1:ni),stat=jerr)                                        1643
      if(jerr.ne.0) return                                                 1644
      allocate(r(1:no),stat=jerr)                                          1645
      if(jerr.ne.0) return                                                 1646
      allocate(v(1:no),stat=jerr)                                          1647
      if(jerr.ne.0) return                                                 1648
      allocate(q(1:no),stat=jerr)                                          1649
      if(jerr.ne.0) return                                                 1650
      fmax=log(1.0/pmin-1.0)                                               1650
      fmin=-fmax                                                           1650
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1651
      bta=parm                                                             1651
      omb=1.0-bta                                                          1652
      q0=dot_product(w,y)                                                  1652
      if(q0 .gt. pmin)goto 12681                                           1652
      jerr=8001                                                            1652
      return                                                               1652
12681 continue                                                             1653
      if(q0 .lt. 1.0-pmin)goto 12701                                       1653
      jerr=9001                                                            1653
      return                                                               1653
12701 continue                                                             1654
      if(intr.eq.0.0) q0=0.5                                               1655
      ixx=0                                                                1655
      al=0.0                                                               1655
      bz=0.0                                                               1655
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                    1656
      if(nonzero(no,g) .ne. 0)goto 12721                                   1656
      vi=q0*(1.0-q0)                                                       1656
      b(0)=bz                                                              1656
      v=vi*w                                                               1657
      r=w*(y-q0)                                                           1657
      q=q0                                                                 1657
      xmz=vi                                                               1657
      dev1=-(bz*q0+log(1.0-q0))                                            1658
      goto 12731                                                           1659
12721 continue                                                             1659
      b(0)=0.0                                                             1660
      if(intr .eq. 0)goto 12751                                            1660
      b(0)=azero(no,y,g,w,jerr)                                            1660
      if(jerr.ne.0) return                                                 1660
12751 continue                                                             1661
      q=1.0/(1.0+exp(-b(0)-g))                                             1661
      v=w*q*(1.0-q)                                                        1661
      r=w*(y-q)                                                            1661
      xmz=sum(v)                                                           1662
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1663
12731 continue                                                             1664
12711 continue                                                             1664
      if(kopt .le. 0)goto 12771                                            1665
      if(isd .le. 0 .or. intr .eq. 0)goto 12791                            1665
      xv=0.25                                                              1665
      goto 12801                                                           1666
12791 continue                                                             1666
12810 do 12811 j=1,ni                                                      1666
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1666
12811 continue                                                             1666
12812 continue                                                             1666
12801 continue                                                             1667
12781 continue                                                             1667
12771 continue                                                             1668
      dev0=dev1                                                            1669
12820 do 12821 i=1,no                                                      1669
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1670
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1671
12821 continue                                                             1673
12822 continue                                                             1673
      alf=1.0                                                              1675
      if(flmin .ge. 1.0)goto 12841                                         1675
      eqs=max(eps,flmin)                                                   1675
      alf=eqs**(1.0/(nlam-1))                                              1675
12841 continue                                                             1676
      m=0                                                                  1676
      mm=0                                                                 1676
      nlp=0                                                                1676
      nin=nlp                                                              1676
      mnl=min(mnlam,nlam)                                                  1676
      bs=0.0                                                               1676
      b(1:ni)=0.0                                                          1677
      shr=shri*dev0                                                        1678
12850 do 12851 j=1,ni                                                      1678
      if(ju(j).eq.0)goto 12851                                             1678
      ga(j)=abs(dot_product(r,x(:,j)))                                     1678
12851 continue                                                             1679
12852 continue                                                             1679
12860 do 12861 ilm=1,nlam                                                  1679
      al0=al                                                               1680
      if(flmin .lt. 1.0)goto 12881                                         1680
      al=ulam(ilm)                                                         1680
      goto 12871                                                           1681
12881 if(ilm .le. 2)goto 12891                                             1681
      al=al*alf                                                            1681
      goto 12871                                                           1682
12891 if(ilm .ne. 1)goto 12901                                             1682
      al=big                                                               1682
      goto 12911                                                           1683
12901 continue                                                             1683
      al0=0.0                                                              1684
12920 do 12921 j=1,ni                                                      1684
      if(ju(j).eq.0)goto 12921                                             1684
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1684
12921 continue                                                             1685
12922 continue                                                             1685
      al0=al0/max(bta,1.0d-3)                                              1685
      al=alf*al0                                                           1686
12911 continue                                                             1687
12871 continue                                                             1687
      al2=al*omb                                                           1687
      al1=al*bta                                                           1687
      tlam=bta*(2.0*al-al0)                                                1688
12930 do 12931 k=1,ni                                                      1688
      if(ixx(k).eq.1)goto 12931                                            1688
      if(ju(k).eq.0)goto 12931                                             1689
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1690
12931 continue                                                             1691
12932 continue                                                             1691
10880 continue                                                             1692
12940 continue                                                             1692
12941 continue                                                             1692
      bs(0)=b(0)                                                           1692
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1693
      if(kopt .ne. 0)goto 12961                                            1694
12970 do 12971 j=1,ni                                                      1694
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                       1694
12971 continue                                                             1695
12972 continue                                                             1695
12961 continue                                                             1696
12980 continue                                                             1696
12981 continue                                                             1696
      nlp=nlp+1                                                            1696
      dlx=0.0                                                              1697
12990 do 12991 k=1,ni                                                      1697
      if(ixx(k).eq.0)goto 12991                                            1698
      bk=b(k)                                                              1698
      gk=dot_product(r,x(:,k))                                             1699
      u=gk+xv(k)*b(k)                                                      1699
      au=abs(u)-vp(k)*al1                                                  1700
      if(au .gt. 0.0)goto 13011                                            1700
      b(k)=0.0                                                             1700
      goto 13021                                                           1701
13011 continue                                                             1702
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          1703
13021 continue                                                             1704
13001 continue                                                             1704
      d=b(k)-bk                                                            1704
      if(abs(d).le.0.0)goto 12991                                          1704
      dlx=max(dlx,xv(k)*d**2)                                              1705
      r=r-d*v*x(:,k)                                                       1706
      if(mm(k) .ne. 0)goto 13041                                           1706
      nin=nin+1                                                            1706
      if(nin.gt.nx)goto 12992                                              1707
      mm(k)=nin                                                            1707
      m(nin)=k                                                             1708
13041 continue                                                             1709
12991 continue                                                             1710
12992 continue                                                             1710
      if(nin.gt.nx)goto 12982                                              1711
      d=0.0                                                                1711
      if(intr.ne.0) d=sum(r)/xmz                                           1712
      if(d .eq. 0.0)goto 13061                                             1712
      b(0)=b(0)+d                                                          1712
      dlx=max(dlx,xmz*d**2)                                                1712
      r=r-d*v                                                              1712
13061 continue                                                             1713
      if(dlx.lt.shr)goto 12982                                             1713
      if(nlp .le. maxit)goto 13081                                         1713
      jerr=-ilm                                                            1713
      return                                                               1713
13081 continue                                                             1714
13090 continue                                                             1714
13091 continue                                                             1714
      nlp=nlp+1                                                            1714
      dlx=0.0                                                              1715
13100 do 13101 l=1,nin                                                     1715
      k=m(l)                                                               1715
      bk=b(k)                                                              1716
      gk=dot_product(r,x(:,k))                                             1717
      u=gk+xv(k)*b(k)                                                      1717
      au=abs(u)-vp(k)*al1                                                  1718
      if(au .gt. 0.0)goto 13121                                            1718
      b(k)=0.0                                                             1718
      goto 13131                                                           1719
13121 continue                                                             1720
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          1721
13131 continue                                                             1722
13111 continue                                                             1722
      d=b(k)-bk                                                            1722
      if(abs(d).le.0.0)goto 13101                                          1722
      dlx=max(dlx,xv(k)*d**2)                                              1723
      r=r-d*v*x(:,k)                                                       1724
13101 continue                                                             1725
13102 continue                                                             1725
      d=0.0                                                                1725
      if(intr.ne.0) d=sum(r)/xmz                                           1726
      if(d .eq. 0.0)goto 13151                                             1726
      b(0)=b(0)+d                                                          1726
      dlx=max(dlx,xmz*d**2)                                                1726
      r=r-d*v                                                              1726
13151 continue                                                             1727
      if(dlx.lt.shr)goto 13092                                             1727
      if(nlp .le. maxit)goto 13171                                         1727
      jerr=-ilm                                                            1727
      return                                                               1727
13171 continue                                                             1728
      goto 13091                                                           1729
13092 continue                                                             1729
      goto 12981                                                           1730
12982 continue                                                             1730
      if(nin.gt.nx)goto 12942                                              1731
13180 do 13181 i=1,no                                                      1731
      fi=b(0)+g(i)                                                         1732
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1733
      if(fi .ge. fmin)goto 13201                                           1733
      q(i)=0.0                                                             1733
      goto 13191                                                           1733
13201 if(fi .le. fmax)goto 13211                                           1733
      q(i)=1.0                                                             1733
      goto 13221                                                           1734
13211 continue                                                             1734
      q(i)=1.0/(1.0+exp(-fi))                                              1734
13221 continue                                                             1735
13191 continue                                                             1735
13181 continue                                                             1736
13182 continue                                                             1736
      v=w*q*(1.0-q)                                                        1736
      xmz=sum(v)                                                           1736
      if(xmz.le.vmin)goto 12942                                            1736
      r=w*(y-q)                                                            1737
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 13241                           1737
      ix=0                                                                 1738
13250 do 13251 j=1,nin                                                     1738
      k=m(j)                                                               1739
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 13251                           1739
      ix=1                                                                 1739
      goto 13252                                                           1740
13251 continue                                                             1741
13252 continue                                                             1741
      if(ix .ne. 0)goto 13271                                              1742
13280 do 13281 k=1,ni                                                      1742
      if(ixx(k).eq.1)goto 13281                                            1742
      if(ju(k).eq.0)goto 13281                                             1743
      ga(k)=abs(dot_product(r,x(:,k)))                                     1744
      if(ga(k) .le. al1*vp(k))goto 13301                                   1744
      ixx(k)=1                                                             1744
      ix=1                                                                 1744
13301 continue                                                             1745
13281 continue                                                             1746
13282 continue                                                             1746
      if(ix.eq.1) go to 10880                                              1747
      goto 12942                                                           1748
13271 continue                                                             1749
13241 continue                                                             1750
      goto 12941                                                           1751
12942 continue                                                             1751
      if(nin .le. nx)goto 13321                                            1751
      jerr=-10000-ilm                                                      1751
      goto 12862                                                           1751
13321 continue                                                             1752
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1752
      kin(ilm)=nin                                                         1753
      a0(ilm)=b(0)                                                         1753
      alm(ilm)=al                                                          1753
      lmu=ilm                                                              1754
      devi=dev2(no,w,y,q,pmin)                                             1755
      dev(ilm)=(dev1-devi)/dev0                                            1755
      if(xmz.le.vmin)goto 12862                                            1756
      if(ilm.lt.mnl)goto 12861                                             1756
      if(flmin.ge.1.0)goto 12861                                           1757
      me=0                                                                 1757
13330 do 13331 j=1,nin                                                     1757
      if(a(j,ilm).ne.0.0) me=me+1                                          1757
13331 continue                                                             1757
13332 continue                                                             1757
      if(me.gt.ne)goto 12862                                               1758
      if(dev(ilm).gt.devmax)goto 12862                                     1758
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12862                             1759
12861 continue                                                             1760
12862 continue                                                             1760
      g=log(q/(1.0-q))                                                     1761
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                                  1762
      return                                                               1763
      end                                                                  1764
      function dev2(n,w,y,p,pmin)                                          1765
      implicit double precision(a-h,o-z)                                   1766
      double precision w(n),y(n),p(n)                                      1767
      pmax=1.0-pmin                                                        1767
      s=0.0                                                                1768
13340 do 13341 i=1,n                                                       1768
      pi=min(max(pmin,p(i)),pmax)                                          1769
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1770
13341 continue                                                             1771
13342 continue                                                             1771
      dev2=s                                                               1772
      return                                                               1773
      end                                                                  1774
      function azero(n,y,g,q,jerr)                                         1775
      implicit double precision(a-h,o-z)                                   1776
      parameter(eps=1.0d-7)                                                1777
      double precision y(n),g(n),q(n)                                      1778
      double precision, dimension (:), allocatable :: e,p,w                     
      azero = 0.0                                                          1782
      allocate(e(1:n),stat=jerr)                                           1783
      if(jerr.ne.0) return                                                 1784
      allocate(p(1:n),stat=jerr)                                           1785
      if(jerr.ne.0) return                                                 1786
      allocate(w(1:n),stat=jerr)                                           1787
      if(jerr.ne.0) return                                                 1788
      az=0.0                                                               1788
      e=exp(-g)                                                            1788
      qy=dot_product(q,y)                                                  1788
      p=1.0/(1.0+e)                                                        1789
13350 continue                                                             1789
13351 continue                                                             1789
      w=q*p*(1.0-p)                                                        1790
      d=(qy-dot_product(q,p))/sum(w)                                       1790
      az=az+d                                                              1790
      if(abs(d).lt.eps)goto 13352                                          1791
      ea0=exp(-az)                                                         1791
      p=1.0/(1.0+ea0*e)                                                    1792
      goto 13351                                                           1793
13352 continue                                                             1793
      azero=az                                                             1794
      deallocate(e,p,w)                                                    1795
      return                                                               1796
      end                                                                  1797
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin   1799 
     *,ulam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,j
     *err)
      implicit double precision(a-h,o-z)                                   1800
      double precision x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam   1801 
     *)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(   1802 
     *2,ni)
      integer ju(ni),m(nx),kin(nlam)                                       1803
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
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               1818
      exmn=-exmx                                                           1819
      allocate(r(1:no),stat=jerr)                                          1820
      if(jerr.ne.0) return                                                 1821
      allocate(v(1:no),stat=jerr)                                          1822
      if(jerr.ne.0) return                                                 1823
      allocate(mm(1:ni),stat=jerr)                                         1824
      if(jerr.ne.0) return                                                 1825
      allocate(is(1:max(nc,ni)),stat=jerr)                                 1826
      if(jerr.ne.0) return                                                 1827
      allocate(sxp(1:no),stat=jerr)                                        1828
      if(jerr.ne.0) return                                                 1829
      allocate(sxpl(1:no),stat=jerr)                                       1830
      if(jerr.ne.0) return                                                 1831
      allocate(di(1:no),stat=jerr)                                         1832
      if(jerr.ne.0) return                                                 1833
      allocate(ga(1:ni),stat=jerr)                                         1834
      if(jerr.ne.0) return                                                 1835
      allocate(ixx(1:ni),stat=jerr)                                        1836
      if(jerr.ne.0) return                                                 1837
      pmax=1.0-pmin                                                        1837
      emin=pmin/pmax                                                       1837
      emax=1.0/emin                                                        1838
      pfm=(1.0+pmin)*pmin                                                  1838
      pfx=(1.0-pmin)*pmax                                                  1838
      vmin=pfm*pmax                                                        1839
      bta=parm                                                             1839
      omb=1.0-bta                                                          1839
      dev1=0.0                                                             1839
      dev0=0.0                                                             1840
13360 do 13361 ic=1,nc                                                     1840
      q0=dot_product(w,y(:,ic))                                            1841
      if(q0 .gt. pmin)goto 13381                                           1841
      jerr =8000+ic                                                        1841
      return                                                               1841
13381 continue                                                             1842
      if(q0 .lt. 1.0-pmin)goto 13401                                       1842
      jerr =9000+ic                                                        1842
      return                                                               1842
13401 continue                                                             1843
      if(intr .ne. 0)goto 13421                                            1843
      q0=1.0/nc                                                            1843
      b(0,ic)=0.0                                                          1843
      goto 13431                                                           1844
13421 continue                                                             1844
      b(0,ic)=log(q0)                                                      1844
      dev1=dev1-q0*b(0,ic)                                                 1844
13431 continue                                                             1845
13411 continue                                                             1845
      b(1:ni,ic)=0.0                                                       1846
13361 continue                                                             1847
13362 continue                                                             1847
      if(intr.eq.0) dev1=log(float(nc))                                    1847
      ixx=0                                                                1847
      al=0.0                                                               1848
      if(nonzero(no*nc,g) .ne. 0)goto 13451                                1849
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1849
      sxp=0.0                                                              1850
13460 do 13461 ic=1,nc                                                     1850
      q(:,ic)=exp(b(0,ic))                                                 1850
      sxp=sxp+q(:,ic)                                                      1850
13461 continue                                                             1851
13462 continue                                                             1851
      goto 13471                                                           1852
13451 continue                                                             1852
13480 do 13481 i=1,no                                                      1852
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1852
13481 continue                                                             1852
13482 continue                                                             1852
      sxp=0.0                                                              1853
      if(intr .ne. 0)goto 13501                                            1853
      b(0,:)=0.0                                                           1853
      goto 13511                                                           1854
13501 continue                                                             1854
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1854
      if(jerr.ne.0) return                                                 1854
13511 continue                                                             1855
13491 continue                                                             1855
      dev1=0.0                                                             1856
13520 do 13521 ic=1,nc                                                     1856
      q(:,ic)=b(0,ic)+g(:,ic)                                              1857
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             1858
      q(:,ic)=exp(q(:,ic))                                                 1858
      sxp=sxp+q(:,ic)                                                      1859
13521 continue                                                             1860
13522 continue                                                             1860
      sxpl=w*log(sxp)                                                      1860
13530 do 13531 ic=1,nc                                                     1860
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  1860
13531 continue                                                             1861
13532 continue                                                             1861
13471 continue                                                             1862
13441 continue                                                             1862
13540 do 13541 ic=1,nc                                                     1862
13550 do 13551 i=1,no                                                      1862
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1862
13551 continue                                                             1862
13552 continue                                                             1862
13541 continue                                                             1863
13542 continue                                                             1863
      dev0=dev0+dev1                                                       1864
      if(kopt .le. 0)goto 13571                                            1865
      if(isd .le. 0 .or. intr .eq. 0)goto 13591                            1865
      xv=0.25                                                              1865
      goto 13601                                                           1866
13591 continue                                                             1866
13610 do 13611 j=1,ni                                                      1866
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1866
13611 continue                                                             1866
13612 continue                                                             1866
13601 continue                                                             1867
13581 continue                                                             1867
13571 continue                                                             1869
      alf=1.0                                                              1871
      if(flmin .ge. 1.0)goto 13631                                         1871
      eqs=max(eps,flmin)                                                   1871
      alf=eqs**(1.0/(nlam-1))                                              1871
13631 continue                                                             1872
      m=0                                                                  1872
      mm=0                                                                 1872
      nin=0                                                                1872
      nlp=0                                                                1872
      mnl=min(mnlam,nlam)                                                  1872
      bs=0.0                                                               1872
      shr=shri*dev0                                                        1873
      ga=0.0                                                               1874
13640 do 13641 ic=1,nc                                                     1874
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1875
13650 do 13651 j=1,ni                                                      1875
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))           1875
13651 continue                                                             1876
13652 continue                                                             1876
13641 continue                                                             1877
13642 continue                                                             1877
13660 do 13661 ilm=1,nlam                                                  1877
      al0=al                                                               1878
      if(flmin .lt. 1.0)goto 13681                                         1878
      al=ulam(ilm)                                                         1878
      goto 13671                                                           1879
13681 if(ilm .le. 2)goto 13691                                             1879
      al=al*alf                                                            1879
      goto 13671                                                           1880
13691 if(ilm .ne. 1)goto 13701                                             1880
      al=big                                                               1880
      goto 13711                                                           1881
13701 continue                                                             1881
      al0=0.0                                                              1882
13720 do 13721 j=1,ni                                                      1882
      if(ju(j).eq.0)goto 13721                                             1882
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1882
13721 continue                                                             1883
13722 continue                                                             1883
      al0=al0/max(bta,1.0d-3)                                              1883
      al=alf*al0                                                           1884
13711 continue                                                             1885
13671 continue                                                             1885
      al2=al*omb                                                           1885
      al1=al*bta                                                           1885
      tlam=bta*(2.0*al-al0)                                                1886
13730 do 13731 k=1,ni                                                      1886
      if(ixx(k).eq.1)goto 13731                                            1886
      if(ju(k).eq.0)goto 13731                                             1887
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1888
13731 continue                                                             1889
13732 continue                                                             1889
10880 continue                                                             1890
13740 continue                                                             1890
13741 continue                                                             1890
      ix=0                                                                 1890
      jx=ix                                                                1890
      ig=0                                                                 1891
13750 do 13751 ic=1,nc                                                     1891
      bs(0,ic)=b(0,ic)                                                     1892
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1893
      xmz=0.0                                                              1894
13760 do 13761 i=1,no                                                      1894
      pic=q(i,ic)/sxp(i)                                                   1895
      if(pic .ge. pfm)goto 13781                                           1895
      pic=0.0                                                              1895
      v(i)=0.0                                                             1895
      goto 13771                                                           1896
13781 if(pic .le. pfx)goto 13791                                           1896
      pic=1.0                                                              1896
      v(i)=0.0                                                             1896
      goto 13801                                                           1897
13791 continue                                                             1897
      v(i)=w(i)*pic*(1.0-pic)                                              1897
      xmz=xmz+v(i)                                                         1897
13801 continue                                                             1898
13771 continue                                                             1898
      r(i)=w(i)*(y(i,ic)-pic)                                              1899
13761 continue                                                             1900
13762 continue                                                             1900
      if(xmz.le.vmin)goto 13751                                            1900
      ig=1                                                                 1901
      if(kopt .ne. 0)goto 13821                                            1902
13830 do 13831 j=1,ni                                                      1902
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                    1902
13831 continue                                                             1903
13832 continue                                                             1903
13821 continue                                                             1904
13840 continue                                                             1904
13841 continue                                                             1904
      nlp=nlp+1                                                            1904
      dlx=0.0                                                              1905
13850 do 13851 k=1,ni                                                      1905
      if(ixx(k).eq.0)goto 13851                                            1906
      bk=b(k,ic)                                                           1906
      gk=dot_product(r,x(:,k))                                             1907
      u=gk+xv(k,ic)*b(k,ic)                                                1907
      au=abs(u)-vp(k)*al1                                                  1908
      if(au .gt. 0.0)goto 13871                                            1908
      b(k,ic)=0.0                                                          1908
      goto 13881                                                           1909
13871 continue                                                             1910
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   1912 
     *)
13881 continue                                                             1913
13861 continue                                                             1913
      d=b(k,ic)-bk                                                         1913
      if(abs(d).le.0.0)goto 13851                                          1914
      dlx=max(dlx,xv(k,ic)*d**2)                                           1914
      r=r-d*v*x(:,k)                                                       1915
      if(mm(k) .ne. 0)goto 13901                                           1915
      nin=nin+1                                                            1916
      if(nin .le. nx)goto 13921                                            1916
      jx=1                                                                 1916
      goto 13852                                                           1916
13921 continue                                                             1917
      mm(k)=nin                                                            1917
      m(nin)=k                                                             1918
13901 continue                                                             1919
13851 continue                                                             1920
13852 continue                                                             1920
      if(jx.gt.0)goto 13842                                                1921
      d=0.0                                                                1921
      if(intr.ne.0) d=sum(r)/xmz                                           1922
      if(d .eq. 0.0)goto 13941                                             1922
      b(0,ic)=b(0,ic)+d                                                    1922
      dlx=max(dlx,xmz*d**2)                                                1922
      r=r-d*v                                                              1922
13941 continue                                                             1923
      if(dlx.lt.shr)goto 13842                                             1924
      if(nlp .le. maxit)goto 13961                                         1924
      jerr=-ilm                                                            1924
      return                                                               1924
13961 continue                                                             1925
13970 continue                                                             1925
13971 continue                                                             1925
      nlp=nlp+1                                                            1925
      dlx=0.0                                                              1926
13980 do 13981 l=1,nin                                                     1926
      k=m(l)                                                               1926
      bk=b(k,ic)                                                           1927
      gk=dot_product(r,x(:,k))                                             1928
      u=gk+xv(k,ic)*b(k,ic)                                                1928
      au=abs(u)-vp(k)*al1                                                  1929
      if(au .gt. 0.0)goto 14001                                            1929
      b(k,ic)=0.0                                                          1929
      goto 14011                                                           1930
14001 continue                                                             1931
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   1933 
     *)
14011 continue                                                             1934
13991 continue                                                             1934
      d=b(k,ic)-bk                                                         1934
      if(abs(d).le.0.0)goto 13981                                          1935
      dlx=max(dlx,xv(k,ic)*d**2)                                           1935
      r=r-d*v*x(:,k)                                                       1936
13981 continue                                                             1937
13982 continue                                                             1937
      d=0.0                                                                1937
      if(intr.ne.0) d=sum(r)/xmz                                           1938
      if(d .eq. 0.0)goto 14031                                             1938
      b(0,ic)=b(0,ic)+d                                                    1939
      dlx=max(dlx,xmz*d**2)                                                1939
      r=r-d*v                                                              1940
14031 continue                                                             1941
      if(dlx.lt.shr)goto 13972                                             1941
      if(nlp .le. maxit)goto 14051                                         1941
      jerr=-ilm                                                            1941
      return                                                               1941
14051 continue                                                             1942
      goto 13971                                                           1943
13972 continue                                                             1943
      goto 13841                                                           1944
13842 continue                                                             1944
      if(jx.gt.0)goto 13752                                                1945
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                            1946
      if(ix .ne. 0)goto 14071                                              1947
14080 do 14081 j=1,nin                                                     1947
      k=m(j)                                                               1948
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 14101                1948
      ix=1                                                                 1948
      goto 14082                                                           1948
14101 continue                                                             1949
14081 continue                                                             1950
14082 continue                                                             1950
14071 continue                                                             1951
14110 do 14111 i=1,no                                                      1951
      fi=b(0,ic)+g(i,ic)                                                   1953
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1954
      fi=min(max(exmn,fi),exmx)                                            1954
      sxp(i)=sxp(i)-q(i,ic)                                                1955
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    1956
      sxp(i)=sxp(i)+q(i,ic)                                                1957
14111 continue                                                             1958
14112 continue                                                             1958
13751 continue                                                             1959
13752 continue                                                             1959
      s=-sum(b(0,:))/nc                                                    1959
      b(0,:)=b(0,:)+s                                                      1959
      di=s                                                                 1960
14120 do 14121 j=1,nin                                                     1960
      l=m(j)                                                               1961
      if(vp(l) .gt. 0.0)goto 14141                                         1961
      s=sum(b(l,:))/nc                                                     1961
      goto 14151                                                           1962
14141 continue                                                             1962
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                     1962
14151 continue                                                             1963
14131 continue                                                             1963
      b(l,:)=b(l,:)-s                                                      1963
      di=di-s*x(:,l)                                                       1964
14121 continue                                                             1965
14122 continue                                                             1965
      di=exp(di)                                                           1965
      sxp=sxp*di                                                           1965
14160 do 14161 ic=1,nc                                                     1965
      q(:,ic)=q(:,ic)*di                                                   1965
14161 continue                                                             1966
14162 continue                                                             1966
      if(jx.gt.0)goto 13742                                                1966
      if(ig.eq.0)goto 13742                                                1967
      if(ix .ne. 0)goto 14181                                              1968
14190 do 14191 k=1,ni                                                      1968
      if(ixx(k).eq.1)goto 14191                                            1968
      if(ju(k).eq.0)goto 14191                                             1968
      ga(k)=0.0                                                            1968
14191 continue                                                             1969
14192 continue                                                             1969
14200 do 14201 ic=1,nc                                                     1969
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1970
14210 do 14211 k=1,ni                                                      1970
      if(ixx(k).eq.1)goto 14211                                            1970
      if(ju(k).eq.0)goto 14211                                             1971
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                          1972
14211 continue                                                             1973
14212 continue                                                             1973
14201 continue                                                             1974
14202 continue                                                             1974
14220 do 14221 k=1,ni                                                      1974
      if(ixx(k).eq.1)goto 14221                                            1974
      if(ju(k).eq.0)goto 14221                                             1975
      if(ga(k) .le. al1*vp(k))goto 14241                                   1975
      ixx(k)=1                                                             1975
      ix=1                                                                 1975
14241 continue                                                             1976
14221 continue                                                             1977
14222 continue                                                             1977
      if(ix.eq.1) go to 10880                                              1978
      goto 13742                                                           1979
14181 continue                                                             1980
      goto 13741                                                           1981
13742 continue                                                             1981
      if(jx .le. 0)goto 14261                                              1981
      jerr=-10000-ilm                                                      1981
      goto 13662                                                           1981
14261 continue                                                             1981
      devi=0.0                                                             1982
14270 do 14271 ic=1,nc                                                     1983
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1983
      a0(ic,ilm)=b(0,ic)                                                   1984
14280 do 14281 i=1,no                                                      1984
      if(y(i,ic).le.0.0)goto 14281                                         1985
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1986
14281 continue                                                             1987
14282 continue                                                             1987
14271 continue                                                             1988
14272 continue                                                             1988
      kin(ilm)=nin                                                         1988
      alm(ilm)=al                                                          1988
      lmu=ilm                                                              1989
      dev(ilm)=(dev1-devi)/dev0                                            1989
      if(ig.eq.0)goto 13662                                                1990
      if(ilm.lt.mnl)goto 13661                                             1990
      if(flmin.ge.1.0)goto 13661                                           1991
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 13662             1992
      if(dev(ilm).gt.devmax)goto 13662                                     1992
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13662                             1993
13661 continue                                                             1994
13662 continue                                                             1994
      g=log(q)                                                             1994
14290 do 14291 i=1,no                                                      1994
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1994
14291 continue                                                             1995
14292 continue                                                             1995
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                           1996
      return                                                               1997
      end                                                                  1998
      subroutine kazero(kk,n,y,g,q,az,jerr)                                1999
      implicit double precision(a-h,o-z)                                   2000
      parameter(eps=1.0d-7)                                                2001
      double precision y(n,kk),g(n,kk),q(n),az(kk)                         2002
      double precision, dimension (:), allocatable :: s                         
      double precision, dimension (:,:), allocatable :: e                       
      allocate(e(1:n,1:kk),stat=jerr)                                           
      if(jerr.ne.0) return                                                      
      allocate(s(1:n),stat=jerr)                                           2009
      if(jerr.ne.0) return                                                 2010
      az=0.0                                                               2010
      e=exp(g)                                                             2010
14300 do 14301 i=1,n                                                       2010
      s(i)=sum(e(i,:))                                                     2010
14301 continue                                                             2011
14302 continue                                                             2011
14310 continue                                                             2011
14311 continue                                                             2011
      dm=0.0                                                               2012
14320 do 14321 k=1,kk                                                      2012
      t=0.0                                                                2012
      u=t                                                                  2013
14330 do 14331 i=1,n                                                       2013
      pik=e(i,k)/s(i)                                                      2014
      t=t+q(i)*(y(i,k)-pik)                                                2014
      u=u+q(i)*pik*(1.0-pik)                                               2015
14331 continue                                                             2016
14332 continue                                                             2016
      d=t/u                                                                2016
      az(k)=az(k)+d                                                        2016
      ed=exp(d)                                                            2016
      dm=max(dm,abs(d))                                                    2017
14340 do 14341 i=1,n                                                       2017
      z=e(i,k)                                                             2017
      e(i,k)=z*ed                                                          2017
      s(i)=s(i)-z+e(i,k)                                                   2017
14341 continue                                                             2018
14342 continue                                                             2018
14321 continue                                                             2019
14322 continue                                                             2019
      if(dm.lt.eps)goto 14312                                              2019
      goto 14311                                                           2020
14312 continue                                                             2020
      az=az-sum(az)/kk                                                     2021
      deallocate(e,s)                                                      2022
      return                                                               2023
      end                                                                  2024
      function elc(parm,n,cl,a,m)                                          2025
      implicit double precision(a-h,o-z)                                   2026
      double precision a(n),cl(2)                                          2026
      integer m(n)                                                         2027
      fn=n                                                                 2027
      am=sum(a)/fn                                                         2028
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 14361                       2028
      elc=am                                                               2028
      go to 14370                                                          2028
14361 continue                                                             2029
14380 do 14381 i=1,n                                                       2029
      m(i)=i                                                               2029
14381 continue                                                             2029
14382 continue                                                             2029
      call psort7(a,m,1,n)                                                 2030
      if(a(m(1)) .ne. a(m(n)))goto 14401                                   2030
      elc=a(1)                                                             2030
      go to 14370                                                          2030
14401 continue                                                             2031
      if(mod(n,2) .ne. 1)goto 14421                                        2031
      ad=a(m(n/2+1))                                                       2031
      goto 14431                                                           2032
14421 continue                                                             2032
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       2032
14431 continue                                                             2033
14411 continue                                                             2033
      if(parm .ne. 1.0)goto 14451                                          2033
      elc=ad                                                               2033
      go to 14370                                                          2033
14451 continue                                                             2034
      b1=min(am,ad)                                                        2034
      b2=max(am,ad)                                                        2034
      k2=1                                                                 2035
14460 continue                                                             2035
14461 if(a(m(k2)).gt.b1)goto 14462                                         2035
      k2=k2+1                                                              2035
      goto 14461                                                           2035
14462 continue                                                             2035
      k1=k2-1                                                              2036
14470 continue                                                             2036
14471 if(a(m(k2)).ge.b2)goto 14472                                         2036
      k2=k2+1                                                              2036
      goto 14471                                                           2037
14472 continue                                                             2037
      r=parm/((1.0-parm)*fn)                                               2037
      is=0                                                                 2037
      sm=n-2*(k1-1)                                                        2038
14480 do 14481 k=k1,k2-1                                                   2038
      sm=sm-2.0                                                            2038
      s=r*sm+am                                                            2039
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 14501                   2039
      is=k                                                                 2039
      goto 14482                                                           2039
14501 continue                                                             2040
14481 continue                                                             2041
14482 continue                                                             2041
      if(is .eq. 0)goto 14521                                              2041
      elc=s                                                                2041
      go to 14370                                                          2041
14521 continue                                                             2041
      r2=2.0*r                                                             2041
      s1=a(m(k1))                                                          2041
      am2=2.0*am                                                           2042
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    2042
      elc=s1                                                               2043
14530 do 14531 k=k1+1,k2                                                   2043
      s=a(m(k))                                                            2043
      if(s.eq.s1)goto 14531                                                2044
      c=r2*sum(abs(a-s))+s*(s-am2)                                         2045
      if(c .ge. cri)goto 14551                                             2045
      cri=c                                                                2045
      elc=s                                                                2045
14551 continue                                                             2045
      s1=s                                                                 2046
14531 continue                                                             2047
14532 continue                                                             2047
14370 continue                                                             2047
      elc=max(maxval(a-cl(2)),min(minval(a-cl(1)),elc))                    2048
      return                                                               2049
      end                                                                  2050
      function nintot(ni,nx,nc,a,m,nin,is)                                 2051
      implicit double precision(a-h,o-z)                                   2052
      double precision a(nx,nc)                                            2052
      integer m(nx),is(ni)                                                 2053
      is=0                                                                 2053
      nintot=0                                                             2054
14560 do 14561 ic=1,nc                                                     2054
14570 do 14571 j=1,nin                                                     2054
      k=m(j)                                                               2054
      if(is(k).ne.0)goto 14571                                             2055
      if(a(j,ic).eq.0.0)goto 14571                                         2055
      is(k)=k                                                              2055
      nintot=nintot+1                                                      2056
14571 continue                                                             2056
14572 continue                                                             2056
14561 continue                                                             2057
14562 continue                                                             2057
      return                                                               2058
      end                                                                  2059
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             2060
      implicit double precision(a-h,o-z)                                   2061
      double precision ca(nx,nc),a(ni,nc)                                  2061
      integer ia(nx)                                                       2062
      a=0.0                                                                2063
14580 do 14581 ic=1,nc                                                     2063
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            2063
14581 continue                                                             2064
14582 continue                                                             2064
      return                                                               2065
      end                                                                  2066
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      2067
      implicit double precision(a-h,o-z)                                   2068
      double precision a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                 2068
      integer ia(nx)                                                       2069
14590 do 14591 i=1,nt                                                      2069
14600 do 14601 ic=1,nc                                                     2069
      ans(ic,i)=a0(ic)                                                     2071
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   2072 
     *:nin)))
14601 continue                                                             2072
14602 continue                                                             2072
14591 continue                                                             2073
14592 continue                                                             2073
      return                                                               2074
      end                                                                  2075
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,cl,ne,nx,nlam   2077 
     *,flmin,  ulam,thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,al
     *m,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2078
      double precision x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)     2079
      double precision ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl   2080 
     *(2,ni)
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           2081
      double precision, dimension (:), allocatable :: xm,xs,ww,vq,xv            
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 14621                                    2085
      jerr=10000                                                           2085
      return                                                               2085
14621 continue                                                             2086
      allocate(ww(1:no),stat=jerr)                                         2087
      if(jerr.ne.0) return                                                 2088
      allocate(ju(1:ni),stat=jerr)                                         2089
      if(jerr.ne.0) return                                                 2090
      allocate(vq(1:ni),stat=jerr)                                         2091
      if(jerr.ne.0) return                                                 2092
      allocate(xm(1:ni),stat=jerr)                                         2093
      if(jerr.ne.0) return                                                 2094
      allocate(xs(1:ni),stat=jerr)                                         2095
      if(jerr.ne.0) return                                                 2096
      if(kopt .ne. 2)goto 14641                                            2096
      allocate(xv(1:ni),stat=jerr)                                         2096
      if(jerr.ne.0) return                                                 2096
14641 continue                                                             2098
      call spchkvars(no,ni,x,ix,ju)                                        2099
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2100
      if(maxval(ju) .gt. 0)goto 14661                                      2100
      jerr=7777                                                            2100
      return                                                               2100
14661 continue                                                             2101
      vq=max(0d0,vp)                                                       2101
      vq=vq*ni/sum(vq)                                                     2102
14670 do 14671 i=1,no                                                      2102
      ww(i)=sum(y(i,:))                                                    2102
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                                 2102
14671 continue                                                             2103
14672 continue                                                             2103
      sw=sum(ww)                                                           2103
      ww=ww/sw                                                             2104
      if(nc .ne. 1)goto 14691                                              2104
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                2105
      if(isd .le. 0)goto 14711                                             2105
14720 do 14721 j=1,ni                                                      2105
      cl(:,j)=cl(:,j)*xs(j)                                                2105
14721 continue                                                             2105
14722 continue                                                             2105
14711 continue                                                             2106
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,cl,ne,n   2109 
     *x,nlam,  flmin,ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin
     *,dev0,dev,  alm,nlp,jerr)
      goto 14681                                                           2110
14691 if(kopt .ne. 2)goto 14731                                            2111
      call multsplstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs,xv)         2112
      if(isd .le. 0)goto 14751                                             2112
14760 do 14761 j=1,ni                                                      2112
      cl(:,j)=cl(:,j)*xs(j)                                                2112
14761 continue                                                             2112
14762 continue                                                             2112
14751 continue                                                             2113
      call multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nl   2115 
     *am,flmin,  ulam,thr,intr,maxit,xv,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,
     *alm,nlp,jerr)
      goto 14771                                                           2116
14731 continue                                                             2116
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                2117
      if(isd .le. 0)goto 14791                                             2117
14800 do 14801 j=1,ni                                                      2117
      cl(:,j)=cl(:,j)*xs(j)                                                2117
14801 continue                                                             2117
14802 continue                                                             2117
14791 continue                                                             2118
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,f   2121 
     *lmin,  ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,  ia,nin,dev0,
     *dev,alm,nlp,jerr)
14771 continue                                                             2122
14681 continue                                                             2122
      if(jerr.gt.0) return                                                 2122
      dev0=2.0*sw*dev0                                                     2123
14810 do 14811 k=1,lmu                                                     2123
      nk=nin(k)                                                            2124
14820 do 14821 ic=1,nc                                                     2124
      if(isd .le. 0)goto 14841                                             2124
14850 do 14851 l=1,nk                                                      2124
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      2124
14851 continue                                                             2124
14852 continue                                                             2124
14841 continue                                                             2125
      if(intr .ne. 0)goto 14871                                            2125
      a0(ic,k)=0.0                                                         2125
      goto 14881                                                           2126
14871 continue                                                             2126
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            2126
14881 continue                                                             2127
14861 continue                                                             2127
14821 continue                                                             2128
14822 continue                                                             2128
14811 continue                                                             2129
14812 continue                                                             2129
      deallocate(ww,ju,vq,xm,xs)                                           2129
      if(kopt.eq.2) deallocate(xv)                                         2130
      return                                                               2131
      end                                                                  2132
      subroutine multsplstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs,xv)    2133
      implicit double precision(a-h,o-z)                                   2134
      double precision x(*),w(no),xm(ni),xs(ni),xv(ni)                     2134
      integer ix(*),jx(*),ju(ni)                                           2135
      if(intr .ne. 0)goto 14901                                            2136
14910 do 14911 j=1,ni                                                      2136
      if(ju(j).eq.0)goto 14911                                             2136
      xm(j)=0.0                                                            2136
      jb=ix(j)                                                             2136
      je=ix(j+1)-1                                                         2137
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          2138
      if(isd .eq. 0)goto 14931                                             2138
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            2138
      vc=xv(j)-xbq                                                         2139
      xs(j)=sqrt(vc)                                                       2139
      xv(j)=1.0+xbq/vc                                                     2140
      goto 14941                                                           2141
14931 continue                                                             2141
      xs(j)=1.0                                                            2141
14941 continue                                                             2142
14921 continue                                                             2142
14911 continue                                                             2143
14912 continue                                                             2143
      return                                                               2144
14901 continue                                                             2145
14950 do 14951 j=1,ni                                                      2145
      if(ju(j).eq.0)goto 14951                                             2145
      jb=ix(j)                                                             2145
      je=ix(j+1)-1                                                         2146
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2147
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 2148
      if(isd .le. 0)goto 14971                                             2148
      xs(j)=sqrt(xv(j))                                                    2148
      xv(j)=1.0                                                            2148
14971 continue                                                             2149
14951 continue                                                             2150
14952 continue                                                             2150
      if(isd.eq.0) xs=1.0                                                  2151
      return                                                               2152
      end                                                                  2153
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs)           2154
      implicit double precision(a-h,o-z)                                   2155
      double precision x(*),w(no),xm(ni),xs(ni)                            2155
      integer ix(*),jx(*),ju(ni)                                           2156
      if(intr .ne. 0)goto 14991                                            2157
15000 do 15001 j=1,ni                                                      2157
      if(ju(j).eq.0)goto 15001                                             2157
      xm(j)=0.0                                                            2157
      jb=ix(j)                                                             2157
      je=ix(j+1)-1                                                         2158
      if(isd .eq. 0)goto 15021                                             2159
      vc=dot_product(w(jx(jb:je)),x(jb:je)**2)  -dot_product(w(jx(jb:je)   2161 
     *),x(jb:je))**2
      xs(j)=sqrt(vc)                                                       2162
      goto 15031                                                           2163
15021 continue                                                             2163
      xs(j)=1.0                                                            2163
15031 continue                                                             2164
15011 continue                                                             2164
15001 continue                                                             2165
15002 continue                                                             2165
      return                                                               2166
14991 continue                                                             2167
15040 do 15041 j=1,ni                                                      2167
      if(ju(j).eq.0)goto 15041                                             2167
      jb=ix(j)                                                             2167
      je=ix(j+1)-1                                                         2168
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             2169
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   2170 
     *)**2)
15041 continue                                                             2171
15042 continue                                                             2171
      if(isd.eq.0) xs=1.0                                                  2172
      return                                                               2173
      end                                                                  2174
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,nl   2177 
     *am,  flmin,ulam,shri,isd,intr,maxit,kopt,xb,xs,  lmu,a0,a,m,kin,de
     *v0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2178
      double precision x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)   2179
      double precision a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)             2180
      double precision xb(ni),xs(ni)                                       2180
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2181
      double precision, dimension (:), allocatable :: xm,b,bs,v,r               
      double precision, dimension (:), allocatable :: sc,xv,q,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2187
      allocate(b(0:ni),stat=jerr)                                          2188
      if(jerr.ne.0) return                                                 2189
      allocate(xm(0:ni),stat=jerr)                                         2190
      if(jerr.ne.0) return                                                 2191
      allocate(xv(1:ni),stat=jerr)                                         2192
      if(jerr.ne.0) return                                                 2193
      allocate(bs(0:ni),stat=jerr)                                         2194
      if(jerr.ne.0) return                                                 2195
      allocate(ga(1:ni),stat=jerr)                                         2196
      if(jerr.ne.0) return                                                 2197
      allocate(mm(1:ni),stat=jerr)                                         2198
      if(jerr.ne.0) return                                                 2199
      allocate(ixx(1:ni),stat=jerr)                                        2200
      if(jerr.ne.0) return                                                 2201
      allocate(q(1:no),stat=jerr)                                          2202
      if(jerr.ne.0) return                                                 2203
      allocate(r(1:no),stat=jerr)                                          2204
      if(jerr.ne.0) return                                                 2205
      allocate(v(1:no),stat=jerr)                                          2206
      if(jerr.ne.0) return                                                 2207
      allocate(sc(1:no),stat=jerr)                                         2208
      if(jerr.ne.0) return                                                 2209
      fmax=log(1.0/pmin-1.0)                                               2209
      fmin=-fmax                                                           2209
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      2210
      bta=parm                                                             2210
      omb=1.0-bta                                                          2211
      q0=dot_product(w,y)                                                  2211
      if(q0 .gt. pmin)goto 15061                                           2211
      jerr=8001                                                            2211
      return                                                               2211
15061 continue                                                             2212
      if(q0 .lt. 1.0-pmin)goto 15081                                       2212
      jerr=9001                                                            2212
      return                                                               2212
15081 continue                                                             2213
      if(intr.eq.0) q0=0.5                                                 2213
      bz=0.0                                                               2213
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                    2214
      if(nonzero(no,g) .ne. 0)goto 15101                                   2214
      vi=q0*(1.0-q0)                                                       2214
      b(0)=bz                                                              2214
      v=vi*w                                                               2215
      r=w*(y-q0)                                                           2215
      q=q0                                                                 2215
      xm(0)=vi                                                             2215
      dev1=-(bz*q0+log(1.0-q0))                                            2216
      goto 15111                                                           2217
15101 continue                                                             2217
      b(0)=0.0                                                             2218
      if(intr .eq. 0)goto 15131                                            2218
      b(0)=azero(no,y,g,w,jerr)                                            2218
      if(jerr.ne.0) return                                                 2218
15131 continue                                                             2219
      q=1.0/(1.0+exp(-b(0)-g))                                             2219
      v=w*q*(1.0-q)                                                        2219
      r=w*(y-q)                                                            2219
      xm(0)=sum(v)                                                         2220
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        2221
15111 continue                                                             2222
15091 continue                                                             2222
      if(kopt .le. 0)goto 15151                                            2223
      if(isd .le. 0 .or. intr .eq. 0)goto 15171                            2223
      xv=0.25                                                              2223
      goto 15181                                                           2224
15171 continue                                                             2225
15190 do 15191 j=1,ni                                                      2225
      if(ju(j).eq.0)goto 15191                                             2225
      jb=ix(j)                                                             2225
      je=ix(j+1)-1                                                         2226
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          2227
15191 continue                                                             2228
15192 continue                                                             2228
15181 continue                                                             2229
15161 continue                                                             2229
15151 continue                                                             2230
      b(1:ni)=0.0                                                          2230
      dev0=dev1                                                            2231
15200 do 15201 i=1,no                                                      2231
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        2232
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              2233
15201 continue                                                             2235
15202 continue                                                             2235
      alf=1.0                                                              2237
      if(flmin .ge. 1.0)goto 15221                                         2237
      eqs=max(eps,flmin)                                                   2237
      alf=eqs**(1.0/(nlam-1))                                              2237
15221 continue                                                             2238
      m=0                                                                  2238
      mm=0                                                                 2238
      nin=0                                                                2238
      o=0.0                                                                2238
      svr=o                                                                2238
      mnl=min(mnlam,nlam)                                                  2238
      bs=0.0                                                               2238
      nlp=0                                                                2238
      nin=nlp                                                              2239
      shr=shri*dev0                                                        2239
      al=0.0                                                               2239
      ixx=0                                                                2240
15230 do 15231 j=1,ni                                                      2240
      if(ju(j).eq.0)goto 15231                                             2241
      jb=ix(j)                                                             2241
      je=ix(j+1)-1                                                         2241
      jn=ix(j+1)-ix(j)                                                     2242
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2243
      gj=dot_product(sc(1:jn),x(jb:je))                                    2244
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      2245
15231 continue                                                             2246
15232 continue                                                             2246
15240 do 15241 ilm=1,nlam                                                  2246
      al0=al                                                               2247
      if(flmin .lt. 1.0)goto 15261                                         2247
      al=ulam(ilm)                                                         2247
      goto 15251                                                           2248
15261 if(ilm .le. 2)goto 15271                                             2248
      al=al*alf                                                            2248
      goto 15251                                                           2249
15271 if(ilm .ne. 1)goto 15281                                             2249
      al=big                                                               2249
      goto 15291                                                           2250
15281 continue                                                             2250
      al0=0.0                                                              2251
15300 do 15301 j=1,ni                                                      2251
      if(ju(j).eq.0)goto 15301                                             2251
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2251
15301 continue                                                             2252
15302 continue                                                             2252
      al0=al0/max(bta,1.0d-3)                                              2252
      al=alf*al0                                                           2253
15291 continue                                                             2254
15251 continue                                                             2254
      al2=al*omb                                                           2254
      al1=al*bta                                                           2254
      tlam=bta*(2.0*al-al0)                                                2255
15310 do 15311 k=1,ni                                                      2255
      if(ixx(k).eq.1)goto 15311                                            2255
      if(ju(k).eq.0)goto 15311                                             2256
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2257
15311 continue                                                             2258
15312 continue                                                             2258
10880 continue                                                             2259
15320 continue                                                             2259
15321 continue                                                             2259
      bs(0)=b(0)                                                           2259
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                2260
15330 do 15331 j=1,ni                                                      2260
      if(ixx(j).eq.0)goto 15331                                            2261
      jb=ix(j)                                                             2261
      je=ix(j+1)-1                                                         2261
      jn=ix(j+1)-ix(j)                                                     2262
      sc(1:jn)=v(jx(jb:je))                                                2263
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 2264
      if(kopt .ne. 0)goto 15351                                            2265
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              2266
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                2267
15351 continue                                                             2268
15331 continue                                                             2269
15332 continue                                                             2269
15360 continue                                                             2269
15361 continue                                                             2269
      nlp=nlp+1                                                            2269
      dlx=0.0                                                              2270
15370 do 15371 k=1,ni                                                      2270
      if(ixx(k).eq.0)goto 15371                                            2271
      jb=ix(k)                                                             2271
      je=ix(k+1)-1                                                         2271
      jn=ix(k+1)-ix(k)                                                     2271
      bk=b(k)                                                              2272
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2273
      gk=dot_product(sc(1:jn),x(jb:je))                                    2274
      gk=(gk-svr*xb(k))/xs(k)                                              2275
      u=gk+xv(k)*b(k)                                                      2275
      au=abs(u)-vp(k)*al1                                                  2276
      if(au .gt. 0.0)goto 15391                                            2276
      b(k)=0.0                                                             2276
      goto 15401                                                           2277
15391 continue                                                             2278
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          2279
15401 continue                                                             2280
15381 continue                                                             2280
      d=b(k)-bk                                                            2280
      if(abs(d).le.0.0)goto 15371                                          2280
      dlx=max(dlx,xv(k)*d**2)                                              2281
      if(mm(k) .ne. 0)goto 15421                                           2281
      nin=nin+1                                                            2281
      if(nin.gt.nx)goto 15372                                              2282
      mm(k)=nin                                                            2282
      m(nin)=k                                                             2282
      sc(1:jn)=v(jx(jb:je))                                                2283
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 2284
15421 continue                                                             2285
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2286
      o=o+d*(xb(k)/xs(k))                                                  2287
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2288
15371 continue                                                             2289
15372 continue                                                             2289
      if(nin.gt.nx)goto 15362                                              2290
      d=0.0                                                                2290
      if(intr.ne.0) d=svr/xm(0)                                            2291
      if(d .eq. 0.0)goto 15441                                             2291
      b(0)=b(0)+d                                                          2291
      dlx=max(dlx,xm(0)*d**2)                                              2291
      r=r-d*v                                                              2292
      svr=svr-d*xm(0)                                                      2293
15441 continue                                                             2294
      if(dlx.lt.shr)goto 15362                                             2295
      if(nlp .le. maxit)goto 15461                                         2295
      jerr=-ilm                                                            2295
      return                                                               2295
15461 continue                                                             2296
15470 continue                                                             2296
15471 continue                                                             2296
      nlp=nlp+1                                                            2296
      dlx=0.0                                                              2297
15480 do 15481 l=1,nin                                                     2297
      k=m(l)                                                               2297
      jb=ix(k)                                                             2297
      je=ix(k+1)-1                                                         2298
      jn=ix(k+1)-ix(k)                                                     2298
      bk=b(k)                                                              2299
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2300
      gk=dot_product(sc(1:jn),x(jb:je))                                    2301
      gk=(gk-svr*xb(k))/xs(k)                                              2302
      u=gk+xv(k)*b(k)                                                      2302
      au=abs(u)-vp(k)*al1                                                  2303
      if(au .gt. 0.0)goto 15501                                            2303
      b(k)=0.0                                                             2303
      goto 15511                                                           2304
15501 continue                                                             2305
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          2306
15511 continue                                                             2307
15491 continue                                                             2307
      d=b(k)-bk                                                            2307
      if(abs(d).le.0.0)goto 15481                                          2307
      dlx=max(dlx,xv(k)*d**2)                                              2308
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2309
      o=o+d*(xb(k)/xs(k))                                                  2310
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2311
15481 continue                                                             2312
15482 continue                                                             2312
      d=0.0                                                                2312
      if(intr.ne.0) d=svr/xm(0)                                            2313
      if(d .eq. 0.0)goto 15531                                             2313
      b(0)=b(0)+d                                                          2313
      dlx=max(dlx,xm(0)*d**2)                                              2313
      r=r-d*v                                                              2314
      svr=svr-d*xm(0)                                                      2315
15531 continue                                                             2316
      if(dlx.lt.shr)goto 15472                                             2317
      if(nlp .le. maxit)goto 15551                                         2317
      jerr=-ilm                                                            2317
      return                                                               2317
15551 continue                                                             2318
      goto 15471                                                           2319
15472 continue                                                             2319
      goto 15361                                                           2320
15362 continue                                                             2320
      if(nin.gt.nx)goto 15322                                              2321
      sc=b(0)                                                              2321
      b0=0.0                                                               2322
15560 do 15561 j=1,nin                                                     2322
      l=m(j)                                                               2322
      jb=ix(l)                                                             2322
      je=ix(l+1)-1                                                         2323
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      2324
      b0=b0-b(l)*xb(l)/xs(l)                                               2325
15561 continue                                                             2326
15562 continue                                                             2326
      sc=sc+b0                                                             2327
15570 do 15571 i=1,no                                                      2327
      fi=sc(i)+g(i)                                                        2328
      if(fi .ge. fmin)goto 15591                                           2328
      q(i)=0.0                                                             2328
      goto 15581                                                           2328
15591 if(fi .le. fmax)goto 15601                                           2328
      q(i)=1.0                                                             2328
      goto 15611                                                           2329
15601 continue                                                             2329
      q(i)=1.0/(1.0+exp(-fi))                                              2329
15611 continue                                                             2330
15581 continue                                                             2330
15571 continue                                                             2331
15572 continue                                                             2331
      v=w*q*(1.0-q)                                                        2331
      xm(0)=sum(v)                                                         2331
      if(xm(0).lt.vmin)goto 15322                                          2332
      r=w*(y-q)                                                            2332
      svr=sum(r)                                                           2332
      o=0.0                                                                2333
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 15631                         2333
      kx=0                                                                 2334
15640 do 15641 j=1,nin                                                     2334
      k=m(j)                                                               2335
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 15641                           2335
      kx=1                                                                 2335
      goto 15642                                                           2336
15641 continue                                                             2337
15642 continue                                                             2337
      if(kx .ne. 0)goto 15661                                              2338
15670 do 15671 j=1,ni                                                      2338
      if(ixx(j).eq.1)goto 15671                                            2338
      if(ju(j).eq.0)goto 15671                                             2339
      jb=ix(j)                                                             2339
      je=ix(j+1)-1                                                         2339
      jn=ix(j+1)-ix(j)                                                     2340
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 2341
      gj=dot_product(sc(1:jn),x(jb:je))                                    2342
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      2343
      if(ga(j) .le. al1*vp(j))goto 15691                                   2343
      ixx(j)=1                                                             2343
      kx=1                                                                 2343
15691 continue                                                             2344
15671 continue                                                             2345
15672 continue                                                             2345
      if(kx.eq.1) go to 10880                                              2346
      goto 15322                                                           2347
15661 continue                                                             2348
15631 continue                                                             2349
      goto 15321                                                           2350
15322 continue                                                             2350
      if(nin .le. nx)goto 15711                                            2350
      jerr=-10000-ilm                                                      2350
      goto 15242                                                           2350
15711 continue                                                             2351
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                2351
      kin(ilm)=nin                                                         2352
      a0(ilm)=b(0)                                                         2352
      alm(ilm)=al                                                          2352
      lmu=ilm                                                              2353
      devi=dev2(no,w,y,q,pmin)                                             2354
      dev(ilm)=(dev1-devi)/dev0                                            2355
      if(ilm.lt.mnl)goto 15241                                             2355
      if(flmin.ge.1.0)goto 15241                                           2356
      me=0                                                                 2356
15720 do 15721 j=1,nin                                                     2356
      if(a(j,ilm).ne.0.0) me=me+1                                          2356
15721 continue                                                             2356
15722 continue                                                             2356
      if(me.gt.ne)goto 15242                                               2357
      if(dev(ilm).gt.devmax)goto 15242                                     2357
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 15242                             2358
      if(xm(0).lt.vmin)goto 15242                                          2359
15241 continue                                                             2360
15242 continue                                                             2360
      g=log(q/(1.0-q))                                                     2361
      deallocate(xm,b,bs,v,r,sc,xv,q,mm,ga,ixx)                            2362
      return                                                               2363
      end                                                                  2364
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,n   2366 
     *lam,flmin,  ulam,shri,isd,intr,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev
     *0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2367
      double precision x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb   2368 
     *(ni),xs(ni)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(   2369 
     *2,ni)
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           2370
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
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2385
      exmn=-exmx                                                           2386
      allocate(xm(0:ni),stat=jerr)                                         2387
      if(jerr.ne.0) return                                                 2388
      allocate(r(1:no),stat=jerr)                                          2389
      if(jerr.ne.0) return                                                 2390
      allocate(v(1:no),stat=jerr)                                          2391
      if(jerr.ne.0) return                                                 2392
      allocate(mm(1:ni),stat=jerr)                                         2393
      if(jerr.ne.0) return                                                 2394
      allocate(ga(1:ni),stat=jerr)                                         2395
      if(jerr.ne.0) return                                                 2396
      allocate(iy(1:ni),stat=jerr)                                         2397
      if(jerr.ne.0) return                                                 2398
      allocate(is(1:max(nc,ni)),stat=jerr)                                 2399
      if(jerr.ne.0) return                                                 2400
      allocate(sxp(1:no),stat=jerr)                                        2401
      if(jerr.ne.0) return                                                 2402
      allocate(sxpl(1:no),stat=jerr)                                       2403
      if(jerr.ne.0) return                                                 2404
      allocate(sc(1:no),stat=jerr)                                         2405
      if(jerr.ne.0) return                                                 2406
      pmax=1.0-pmin                                                        2406
      emin=pmin/pmax                                                       2406
      emax=1.0/emin                                                        2407
      pfm=(1.0+pmin)*pmin                                                  2407
      pfx=(1.0-pmin)*pmax                                                  2407
      vmin=pfm*pmax                                                        2408
      bta=parm                                                             2408
      omb=1.0-bta                                                          2408
      dev1=0.0                                                             2408
      dev0=0.0                                                             2409
15730 do 15731 ic=1,nc                                                     2409
      q0=dot_product(w,y(:,ic))                                            2410
      if(q0 .gt. pmin)goto 15751                                           2410
      jerr =8000+ic                                                        2410
      return                                                               2410
15751 continue                                                             2411
      if(q0 .lt. 1.0-pmin)goto 15771                                       2411
      jerr =9000+ic                                                        2411
      return                                                               2411
15771 continue                                                             2412
      if(intr.eq.0) q0=1.0/nc                                              2413
      b(1:ni,ic)=0.0                                                       2413
      b(0,ic)=0.0                                                          2414
      if(intr .eq. 0)goto 15791                                            2414
      b(0,ic)=log(q0)                                                      2414
      dev1=dev1-q0*b(0,ic)                                                 2414
15791 continue                                                             2415
15731 continue                                                             2416
15732 continue                                                             2416
      if(intr.eq.0) dev1=log(float(nc))                                    2416
      iy=0                                                                 2416
      al=0.0                                                               2417
      if(nonzero(no*nc,g) .ne. 0)goto 15811                                2418
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         2418
      sxp=0.0                                                              2419
15820 do 15821 ic=1,nc                                                     2419
      q(:,ic)=exp(b(0,ic))                                                 2419
      sxp=sxp+q(:,ic)                                                      2419
15821 continue                                                             2420
15822 continue                                                             2420
      goto 15831                                                           2421
15811 continue                                                             2421
15840 do 15841 i=1,no                                                      2421
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2421
15841 continue                                                             2421
15842 continue                                                             2421
      sxp=0.0                                                              2422
      if(intr .ne. 0)goto 15861                                            2422
      b(0,:)=0.0                                                           2422
      goto 15871                                                           2423
15861 continue                                                             2423
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 2423
      if(jerr.ne.0) return                                                 2423
15871 continue                                                             2424
15851 continue                                                             2424
      dev1=0.0                                                             2425
15880 do 15881 ic=1,nc                                                     2425
      q(:,ic)=b(0,ic)+g(:,ic)                                              2426
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             2427
      q(:,ic)=exp(q(:,ic))                                                 2427
      sxp=sxp+q(:,ic)                                                      2428
15881 continue                                                             2429
15882 continue                                                             2429
      sxpl=w*log(sxp)                                                      2429
15890 do 15891 ic=1,nc                                                     2429
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  2429
15891 continue                                                             2430
15892 continue                                                             2430
15831 continue                                                             2431
15801 continue                                                             2431
15900 do 15901 ic=1,nc                                                     2431
15910 do 15911 i=1,no                                                      2431
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               2431
15911 continue                                                             2431
15912 continue                                                             2431
15901 continue                                                             2432
15902 continue                                                             2432
      dev0=dev0+dev1                                                       2433
      if(kopt .le. 0)goto 15931                                            2434
      if(isd .le. 0 .or. intr .eq. 0)goto 15951                            2434
      xv=0.25                                                              2434
      goto 15961                                                           2435
15951 continue                                                             2436
15970 do 15971 j=1,ni                                                      2436
      if(ju(j).eq.0)goto 15971                                             2436
      jb=ix(j)                                                             2436
      je=ix(j+1)-1                                                         2437
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        2438
15971 continue                                                             2439
15972 continue                                                             2439
15961 continue                                                             2440
15941 continue                                                             2440
15931 continue                                                             2442
      alf=1.0                                                              2444
      if(flmin .ge. 1.0)goto 15991                                         2444
      eqs=max(eps,flmin)                                                   2444
      alf=eqs**(1.0/(nlam-1))                                              2444
15991 continue                                                             2445
      m=0                                                                  2445
      mm=0                                                                 2445
      nin=0                                                                2445
      nlp=0                                                                2445
      mnl=min(mnlam,nlam)                                                  2445
      bs=0.0                                                               2445
      svr=0.0                                                              2445
      o=0.0                                                                2446
      shr=shri*dev0                                                        2446
      ga=0.0                                                               2447
16000 do 16001 ic=1,nc                                                     2447
      v=q(:,ic)/sxp                                                        2447
      r=w*(y(:,ic)-v)                                                      2447
      v=w*v*(1.0-v)                                                        2448
16010 do 16011 j=1,ni                                                      2448
      if(ju(j).eq.0)goto 16011                                             2449
      jb=ix(j)                                                             2449
      je=ix(j+1)-1                                                         2449
      jn=ix(j+1)-ix(j)                                                     2450
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2451
      gj=dot_product(sc(1:jn),x(jb:je))                                    2452
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2453
16011 continue                                                             2454
16012 continue                                                             2454
16001 continue                                                             2455
16002 continue                                                             2455
16020 do 16021 ilm=1,nlam                                                  2455
      al0=al                                                               2456
      if(flmin .lt. 1.0)goto 16041                                         2456
      al=ulam(ilm)                                                         2456
      goto 16031                                                           2457
16041 if(ilm .le. 2)goto 16051                                             2457
      al=al*alf                                                            2457
      goto 16031                                                           2458
16051 if(ilm .ne. 1)goto 16061                                             2458
      al=big                                                               2458
      goto 16071                                                           2459
16061 continue                                                             2459
      al0=0.0                                                              2460
16080 do 16081 j=1,ni                                                      2460
      if(ju(j).eq.0)goto 16081                                             2460
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2460
16081 continue                                                             2461
16082 continue                                                             2461
      al0=al0/max(bta,1.0d-3)                                              2461
      al=alf*al0                                                           2462
16071 continue                                                             2463
16031 continue                                                             2463
      al2=al*omb                                                           2463
      al1=al*bta                                                           2463
      tlam=bta*(2.0*al-al0)                                                2464
16090 do 16091 k=1,ni                                                      2464
      if(iy(k).eq.1)goto 16091                                             2464
      if(ju(k).eq.0)goto 16091                                             2465
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      2466
16091 continue                                                             2467
16092 continue                                                             2467
10880 continue                                                             2468
16100 continue                                                             2468
16101 continue                                                             2468
      ixx=0                                                                2468
      jxx=ixx                                                              2468
      ig=0                                                                 2469
16110 do 16111 ic=1,nc                                                     2469
      bs(0,ic)=b(0,ic)                                                     2470
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          2471
      xm(0)=0.0                                                            2471
      svr=0.0                                                              2471
      o=0.0                                                                2472
16120 do 16121 i=1,no                                                      2472
      pic=q(i,ic)/sxp(i)                                                   2473
      if(pic .ge. pfm)goto 16141                                           2473
      pic=0.0                                                              2473
      v(i)=0.0                                                             2473
      goto 16131                                                           2474
16141 if(pic .le. pfx)goto 16151                                           2474
      pic=1.0                                                              2474
      v(i)=0.0                                                             2474
      goto 16161                                                           2475
16151 continue                                                             2475
      v(i)=w(i)*pic*(1.0-pic)                                              2475
      xm(0)=xm(0)+v(i)                                                     2475
16161 continue                                                             2476
16131 continue                                                             2476
      r(i)=w(i)*(y(i,ic)-pic)                                              2476
      svr=svr+r(i)                                                         2477
16121 continue                                                             2478
16122 continue                                                             2478
      if(xm(0).le.vmin)goto 16111                                          2478
      ig=1                                                                 2479
16170 do 16171 j=1,ni                                                      2479
      if(iy(j).eq.0)goto 16171                                             2480
      jb=ix(j)                                                             2480
      je=ix(j+1)-1                                                         2481
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             2482
      if(kopt .ne. 0)goto 16191                                            2483
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       2484
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          2485
16191 continue                                                             2486
16171 continue                                                             2487
16172 continue                                                             2487
16200 continue                                                             2487
16201 continue                                                             2487
      nlp=nlp+1                                                            2487
      dlx=0.0                                                              2488
16210 do 16211 k=1,ni                                                      2488
      if(iy(k).eq.0)goto 16211                                             2489
      jb=ix(k)                                                             2489
      je=ix(k+1)-1                                                         2489
      jn=ix(k+1)-ix(k)                                                     2489
      bk=b(k,ic)                                                           2490
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2491
      gk=dot_product(sc(1:jn),x(jb:je))                                    2492
      gk=(gk-svr*xb(k))/xs(k)                                              2493
      u=gk+xv(k,ic)*b(k,ic)                                                2493
      au=abs(u)-vp(k)*al1                                                  2494
      if(au .gt. 0.0)goto 16231                                            2494
      b(k,ic)=0.0                                                          2494
      goto 16241                                                           2495
16231 continue                                                             2496
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   2498 
     *)
16241 continue                                                             2499
16221 continue                                                             2499
      d=b(k,ic)-bk                                                         2499
      if(abs(d).le.0.0)goto 16211                                          2500
      dlx=max(dlx,xv(k,ic)*d**2)                                           2501
      if(mm(k) .ne. 0)goto 16261                                           2501
      nin=nin+1                                                            2502
      if(nin .le. nx)goto 16281                                            2502
      jxx=1                                                                2502
      goto 16212                                                           2502
16281 continue                                                             2503
      mm(k)=nin                                                            2503
      m(nin)=k                                                             2504
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             2505
16261 continue                                                             2506
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2507
      o=o+d*(xb(k)/xs(k))                                                  2508
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2509
16211 continue                                                             2510
16212 continue                                                             2510
      if(jxx.gt.0)goto 16202                                               2511
      d=0.0                                                                2511
      if(intr.ne.0) d=svr/xm(0)                                            2512
      if(d .eq. 0.0)goto 16301                                             2512
      b(0,ic)=b(0,ic)+d                                                    2512
      dlx=max(dlx,xm(0)*d**2)                                              2513
      r=r-d*v                                                              2513
      svr=svr-d*xm(0)                                                      2514
16301 continue                                                             2515
      if(dlx.lt.shr)goto 16202                                             2515
      if(nlp .le. maxit)goto 16321                                         2515
      jerr=-ilm                                                            2515
      return                                                               2515
16321 continue                                                             2516
16330 continue                                                             2516
16331 continue                                                             2516
      nlp=nlp+1                                                            2516
      dlx=0.0                                                              2517
16340 do 16341 l=1,nin                                                     2517
      k=m(l)                                                               2517
      jb=ix(k)                                                             2517
      je=ix(k+1)-1                                                         2518
      jn=ix(k+1)-ix(k)                                                     2518
      bk=b(k,ic)                                                           2519
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2520
      gk=dot_product(sc(1:jn),x(jb:je))                                    2521
      gk=(gk-svr*xb(k))/xs(k)                                              2522
      u=gk+xv(k,ic)*b(k,ic)                                                2522
      au=abs(u)-vp(k)*al1                                                  2523
      if(au .gt. 0.0)goto 16361                                            2523
      b(k,ic)=0.0                                                          2523
      goto 16371                                                           2524
16361 continue                                                             2525
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   2527 
     *)
16371 continue                                                             2528
16351 continue                                                             2528
      d=b(k,ic)-bk                                                         2528
      if(abs(d).le.0.0)goto 16341                                          2529
      dlx=max(dlx,xv(k,ic)*d**2)                                           2530
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              2531
      o=o+d*(xb(k)/xs(k))                                                  2532
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  2533
16341 continue                                                             2534
16342 continue                                                             2534
      d=0.0                                                                2534
      if(intr.ne.0) d=svr/xm(0)                                            2535
      if(d .eq. 0.0)goto 16391                                             2535
      b(0,ic)=b(0,ic)+d                                                    2535
      dlx=max(dlx,xm(0)*d**2)                                              2536
      r=r-d*v                                                              2536
      svr=svr-d*xm(0)                                                      2537
16391 continue                                                             2538
      if(dlx.lt.shr)goto 16332                                             2538
      if(nlp .le. maxit)goto 16411                                         2538
      jerr=-ilm                                                            2538
      return                                                               2538
16411 continue                                                             2539
      goto 16331                                                           2540
16332 continue                                                             2540
      goto 16201                                                           2541
16202 continue                                                             2541
      if(jxx.gt.0)goto 16112                                               2542
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                         2543
      if(ixx .ne. 0)goto 16431                                             2544
16440 do 16441 j=1,nin                                                     2544
      k=m(j)                                                               2545
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 16461                2545
      ixx=1                                                                2545
      goto 16442                                                           2545
16461 continue                                                             2546
16441 continue                                                             2547
16442 continue                                                             2547
16431 continue                                                             2548
      sc=b(0,ic)+g(:,ic)                                                   2548
      b0=0.0                                                               2549
16470 do 16471 j=1,nin                                                     2549
      l=m(j)                                                               2549
      jb=ix(l)                                                             2549
      je=ix(l+1)-1                                                         2550
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   2551
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            2552
16471 continue                                                             2553
16472 continue                                                             2553
      sc=min(max(exmn,sc+b0),exmx)                                         2554
      sxp=sxp-q(:,ic)                                                      2555
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          2556
      sxp=sxp+q(:,ic)                                                      2557
16111 continue                                                             2558
16112 continue                                                             2558
      s=-sum(b(0,:))/nc                                                    2558
      b(0,:)=b(0,:)+s                                                      2558
      sc=s                                                                 2558
      b0=0.0                                                               2559
16480 do 16481 j=1,nin                                                     2559
      l=m(j)                                                               2560
      if(vp(l) .gt. 0.0)goto 16501                                         2560
      s=sum(b(l,:))/nc                                                     2560
      goto 16511                                                           2561
16501 continue                                                             2561
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                     2561
16511 continue                                                             2562
16491 continue                                                             2562
      b(l,:)=b(l,:)-s                                                      2563
      jb=ix(l)                                                             2563
      je=ix(l+1)-1                                                         2564
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         2565
      b0=b0+s*xb(l)/xs(l)                                                  2566
16481 continue                                                             2567
16482 continue                                                             2567
      sc=sc+b0                                                             2567
      sc=exp(sc)                                                           2567
      sxp=sxp*sc                                                           2567
16520 do 16521 ic=1,nc                                                     2567
      q(:,ic)=q(:,ic)*sc                                                   2567
16521 continue                                                             2568
16522 continue                                                             2568
      if(jxx.gt.0)goto 16102                                               2568
      if(ig.eq.0)goto 16102                                                2569
      if(ixx .ne. 0)goto 16541                                             2570
16550 do 16551 j=1,ni                                                      2570
      if(iy(j).eq.1)goto 16551                                             2570
      if(ju(j).eq.0)goto 16551                                             2570
      ga(j)=0.0                                                            2570
16551 continue                                                             2571
16552 continue                                                             2571
16560 do 16561 ic=1,nc                                                     2571
      v=q(:,ic)/sxp                                                        2571
      r=w*(y(:,ic)-v)                                                      2571
      v=w*v*(1.0-v)                                                        2572
16570 do 16571 j=1,ni                                                      2572
      if(iy(j).eq.1)goto 16571                                             2572
      if(ju(j).eq.0)goto 16571                                             2573
      jb=ix(j)                                                             2573
      je=ix(j+1)-1                                                         2573
      jn=ix(j+1)-ix(j)                                                     2574
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 2575
      gj=dot_product(sc(1:jn),x(jb:je))                                    2576
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             2577
16571 continue                                                             2578
16572 continue                                                             2578
16561 continue                                                             2579
16562 continue                                                             2579
16580 do 16581 k=1,ni                                                      2579
      if(iy(k).eq.1)goto 16581                                             2579
      if(ju(k).eq.0)goto 16581                                             2580
      if(ga(k) .le. al1*vp(k))goto 16601                                   2580
      iy(k)=1                                                              2580
      ixx=1                                                                2580
16601 continue                                                             2581
16581 continue                                                             2582
16582 continue                                                             2582
      if(ixx.eq.1) go to 10880                                             2583
      goto 16102                                                           2584
16541 continue                                                             2585
      goto 16101                                                           2586
16102 continue                                                             2586
      if(jxx .le. 0)goto 16621                                             2586
      jerr=-10000-ilm                                                      2586
      goto 16022                                                           2586
16621 continue                                                             2586
      devi=0.0                                                             2587
16630 do 16631 ic=1,nc                                                     2588
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          2588
      a0(ic,ilm)=b(0,ic)                                                   2589
16640 do 16641 i=1,no                                                      2589
      if(y(i,ic).le.0.0)goto 16641                                         2590
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           2591
16641 continue                                                             2592
16642 continue                                                             2592
16631 continue                                                             2593
16632 continue                                                             2593
      kin(ilm)=nin                                                         2593
      alm(ilm)=al                                                          2593
      lmu=ilm                                                              2594
      dev(ilm)=(dev1-devi)/dev0                                            2594
      if(ig.eq.0)goto 16022                                                2595
      if(ilm.lt.mnl)goto 16021                                             2595
      if(flmin.ge.1.0)goto 16021                                           2596
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 16022             2597
      if(dev(ilm).gt.devmax)goto 16022                                     2597
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 16022                             2598
16021 continue                                                             2599
16022 continue                                                             2599
      g=log(q)                                                             2599
16650 do 16651 i=1,no                                                      2599
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         2599
16651 continue                                                             2600
16652 continue                                                             2600
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy)                      2601
      return                                                               2602
      end                                                                  2603
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  2604
      implicit double precision(a-h,o-z)                                   2605
      double precision a0(nc),ca(nx,nc),x(*),f(nc,n)                       2605
      integer ia(*),ix(*),jx(*)                                            2606
16660 do 16661 ic=1,nc                                                     2606
      f(ic,:)=a0(ic)                                                       2606
16661 continue                                                             2607
16662 continue                                                             2607
16670 do 16671 j=1,nin                                                     2607
      k=ia(j)                                                              2607
      kb=ix(k)                                                             2607
      ke=ix(k+1)-1                                                         2608
16680 do 16681 ic=1,nc                                                     2608
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    2608
16681 continue                                                             2609
16682 continue                                                             2609
16671 continue                                                             2610
16672 continue                                                             2610
      return                                                               2611
      end                                                                  2612
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,cl,ne,nx,nlam,flmin,   2614 
     *ulam,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2615
      double precision x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam   2616 
     *)
      double precision ca(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)            2617
      integer jd(*),ia(nx),nin(nlam)                                       2618
      double precision, dimension (:), allocatable :: xs,ww,vq                  
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16701                                    2622
      jerr=10000                                                           2622
      return                                                               2622
16701 continue                                                             2623
      allocate(ww(1:no),stat=jerr)                                         2624
      if(jerr.ne.0) return                                                 2625
      allocate(ju(1:ni),stat=jerr)                                         2626
      if(jerr.ne.0) return                                                 2627
      allocate(vq(1:ni),stat=jerr)                                         2628
      if(jerr.ne.0) return                                                 2629
      if(isd .le. 0)goto 16721                                             2629
      allocate(xs(1:ni),stat=jerr)                                         2629
      if(jerr.ne.0) return                                                 2629
16721 continue                                                             2631
      call chkvars(no,ni,x,ju)                                             2632
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2633
      if(maxval(ju) .gt. 0)goto 16741                                      2633
      jerr=7777                                                            2633
      return                                                               2633
16741 continue                                                             2634
      vq=max(0d0,vp)                                                       2634
      vq=vq*ni/sum(vq)                                                     2635
      ww=max(0d0,w)                                                        2635
      sw=sum(ww)                                                           2636
      if(sw .gt. 0.0)goto 16761                                            2636
      jerr=9999                                                            2636
      return                                                               2636
16761 continue                                                             2636
      ww=ww/sw                                                             2637
      call cstandard(no,ni,x,ww,ju,isd,xs)                                 2638
      if(isd .le. 0)goto 16781                                             2638
16790 do 16791 j=1,ni                                                      2638
      cl(:,j)=cl(:,j)*xs(j)                                                2638
16791 continue                                                             2638
16792 continue                                                             2638
16781 continue                                                             2639
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,   2641 
     *thr,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 2641
      dev0=2.0*sw*dev0                                                     2642
      if(isd .le. 0)goto 16811                                             2642
16820 do 16821 k=1,lmu                                                     2642
      nk=nin(k)                                                            2642
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   2642
16821 continue                                                             2642
16822 continue                                                             2642
16811 continue                                                             2643
      deallocate(ww,ju,vq)                                                 2643
      if(isd.gt.0) deallocate(xs)                                          2644
      return                                                               2645
      end                                                                  2646
      subroutine cstandard (no,ni,x,w,ju,isd,xs)                           2647
      implicit double precision(a-h,o-z)                                   2648
      double precision x(no,ni),w(no),xs(ni)                               2648
      integer ju(ni)                                                       2649
16830 do 16831 j=1,ni                                                      2649
      if(ju(j).eq.0)goto 16831                                             2650
      xm=dot_product(w,x(:,j))                                             2650
      x(:,j)=x(:,j)-xm                                                     2651
      if(isd .le. 0)goto 16851                                             2651
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 2651
      x(:,j)=x(:,j)/xs(j)                                                  2651
16851 continue                                                             2652
16831 continue                                                             2653
16832 continue                                                             2653
      return                                                               2654
      end                                                                  2655
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,cl,ne,nx,nlam,flmin,   2657 
     *ulam,cthri,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2658
      double precision x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam   2659 
     *)
      double precision ao(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)            2660
      integer ju(ni),m(nx),kin(nlam)                                       2661
      double precision, dimension (:), allocatable :: w,dk,v,xs,wr              
      double precision, dimension (:), allocatable :: a,as,f,dq                 
      double precision, dimension (:), allocatable :: e,uu,ga                   
      integer, dimension (:), allocatable :: jp,kp,mm,ixx                       
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2668
      sml=sml*100.0                                                        2668
      devmax=devmax*0.99/0.999                                             2669
      allocate(e(1:no),stat=jerr)                                          2670
      if(jerr.ne.0)go to 12180                                             2671
      allocate(uu(1:no),stat=jerr)                                         2672
      if(jerr.ne.0)go to 12180                                             2673
      allocate(f(1:no),stat=jerr)                                          2674
      if(jerr.ne.0)go to 12180                                             2675
      allocate(w(1:no),stat=jerr)                                          2676
      if(jerr.ne.0)go to 12180                                             2677
      allocate(v(1:ni),stat=jerr)                                          2678
      if(jerr.ne.0)go to 12180                                             2679
      allocate(a(1:ni),stat=jerr)                                          2680
      if(jerr.ne.0)go to 12180                                             2681
      allocate(as(1:ni),stat=jerr)                                         2682
      if(jerr.ne.0)go to 12180                                             2683
      allocate(xs(1:ni),stat=jerr)                                         2684
      if(jerr.ne.0)go to 12180                                             2685
      allocate(ga(1:ni),stat=jerr)                                         2686
      if(jerr.ne.0)go to 12180                                             2687
      allocate(ixx(1:ni),stat=jerr)                                        2688
      if(jerr.ne.0)go to 12180                                             2689
      allocate(jp(1:no),stat=jerr)                                         2690
      if(jerr.ne.0)go to 12180                                             2691
      allocate(kp(1:no),stat=jerr)                                         2692
      if(jerr.ne.0)go to 12180                                             2693
      allocate(dk(1:no),stat=jerr)                                         2694
      if(jerr.ne.0)go to 12180                                             2695
      allocate(wr(1:no),stat=jerr)                                         2696
      if(jerr.ne.0)go to 12180                                             2697
      allocate(dq(1:no),stat=jerr)                                         2698
      if(jerr.ne.0)go to 12180                                             2699
      allocate(mm(1:ni),stat=jerr)                                         2700
      if(jerr.ne.0)go to 12180                                             2701
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2702
      if(jerr.ne.0) go to 12180                                            2702
      alpha=parm                                                           2703
      oma=1.0-alpha                                                        2703
      nlm=0                                                                2703
      ixx=0                                                                2703
      al=0.0                                                               2704
      dq=d*q                                                               2704
      call died(no,nk,dq,kp,jp,dk)                                         2705
      a=0.0                                                                2705
      f(1)=0.0                                                             2705
      fmax=log(huge(f(1))*0.1)                                             2706
      if(nonzero(no,g) .eq. 0)goto 16871                                   2706
      f=g-dot_product(q,g)                                                 2707
      e=q*exp(sign(min(abs(f),fmax),f))                                    2708
      goto 16881                                                           2709
16871 continue                                                             2709
      f=0.0                                                                2709
      e=q                                                                  2709
16881 continue                                                             2710
16861 continue                                                             2710
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 2711
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         2711
      dev0=rr                                                              2712
16890 do 16891 i=1,no                                                      2712
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 16911                   2712
      w(i)=0.0                                                             2712
      wr(i)=w(i)                                                           2712
16911 continue                                                             2712
16891 continue                                                             2713
16892 continue                                                             2713
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2714
      if(jerr.ne.0) go to 12180                                            2716
      alf=1.0                                                              2718
      if(flmin .ge. 1.0)goto 16931                                         2718
      eqs=max(eps,flmin)                                                   2718
      alf=eqs**(1.0/(nlam-1))                                              2718
16931 continue                                                             2719
      m=0                                                                  2719
      mm=0                                                                 2719
      nlp=0                                                                2719
      nin=nlp                                                              2719
      mnl=min(mnlam,nlam)                                                  2719
      as=0.0                                                               2719
      cthr=cthri*dev0                                                      2720
16940 do 16941 j=1,ni                                                      2720
      if(ju(j).eq.0)goto 16941                                             2720
      ga(j)=abs(dot_product(wr,x(:,j)))                                    2720
16941 continue                                                             2721
16942 continue                                                             2721
16950 do 16951 ilm=1,nlam                                                  2721
      al0=al                                                               2722
      if(flmin .lt. 1.0)goto 16971                                         2722
      al=ulam(ilm)                                                         2722
      goto 16961                                                           2723
16971 if(ilm .le. 2)goto 16981                                             2723
      al=al*alf                                                            2723
      goto 16961                                                           2724
16981 if(ilm .ne. 1)goto 16991                                             2724
      al=big                                                               2724
      goto 17001                                                           2725
16991 continue                                                             2725
      al0=0.0                                                              2726
17010 do 17011 j=1,ni                                                      2726
      if(ju(j).eq.0)goto 17011                                             2726
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            2726
17011 continue                                                             2727
17012 continue                                                             2727
      al0=al0/max(parm,1.0d-3)                                             2727
      al=alf*al0                                                           2728
17001 continue                                                             2729
16961 continue                                                             2729
      sa=alpha*al                                                          2729
      omal=oma*al                                                          2729
      tlam=alpha*(2.0*al-al0)                                              2730
17020 do 17021 k=1,ni                                                      2730
      if(ixx(k).eq.1)goto 17021                                            2730
      if(ju(k).eq.0)goto 17021                                             2731
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     2732
17021 continue                                                             2733
17022 continue                                                             2733
10880 continue                                                             2734
17030 continue                                                             2734
17031 continue                                                             2734
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                2735
      call vars(no,ni,x,w,ixx,v)                                           2736
17040 continue                                                             2736
17041 continue                                                             2736
      nlp=nlp+1                                                            2736
      dli=0.0                                                              2737
17050 do 17051 j=1,ni                                                      2737
      if(ixx(j).eq.0)goto 17051                                            2738
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2739
      if(abs(u) .gt. vp(j)*sa)goto 17071                                   2739
      at=0.0                                                               2739
      goto 17081                                                           2740
17071 continue                                                             2740
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o   2742 
     *mal)))
17081 continue                                                             2743
17061 continue                                                             2743
      if(at .eq. a(j))goto 17101                                           2743
      del=at-a(j)                                                          2743
      a(j)=at                                                              2743
      dli=max(dli,v(j)*del**2)                                             2744
      wr=wr-del*w*x(:,j)                                                   2744
      f=f+del*x(:,j)                                                       2745
      if(mm(j) .ne. 0)goto 17121                                           2745
      nin=nin+1                                                            2745
      if(nin.gt.nx)goto 17052                                              2746
      mm(j)=nin                                                            2746
      m(nin)=j                                                             2747
17121 continue                                                             2748
17101 continue                                                             2749
17051 continue                                                             2750
17052 continue                                                             2750
      if(nin.gt.nx)goto 17042                                              2750
      if(dli.lt.cthr)goto 17042                                            2751
      if(nlp .le. maxit)goto 17141                                         2751
      jerr=-ilm                                                            2751
      return                                                               2751
17141 continue                                                             2752
17150 continue                                                             2752
17151 continue                                                             2752
      nlp=nlp+1                                                            2752
      dli=0.0                                                              2753
17160 do 17161 l=1,nin                                                     2753
      j=m(l)                                                               2754
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   2755
      if(abs(u) .gt. vp(j)*sa)goto 17181                                   2755
      at=0.0                                                               2755
      goto 17191                                                           2756
17181 continue                                                             2756
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o   2758 
     *mal)))
17191 continue                                                             2759
17171 continue                                                             2759
      if(at .eq. a(j))goto 17211                                           2759
      del=at-a(j)                                                          2759
      a(j)=at                                                              2759
      dli=max(dli,v(j)*del**2)                                             2760
      wr=wr-del*w*x(:,j)                                                   2760
      f=f+del*x(:,j)                                                       2761
17211 continue                                                             2762
17161 continue                                                             2763
17162 continue                                                             2763
      if(dli.lt.cthr)goto 17152                                            2763
      if(nlp .le. maxit)goto 17231                                         2763
      jerr=-ilm                                                            2763
      return                                                               2763
17231 continue                                                             2764
      goto 17151                                                           2765
17152 continue                                                             2765
      goto 17041                                                           2766
17042 continue                                                             2766
      if(nin.gt.nx)goto 17032                                              2767
      e=q*exp(sign(min(abs(f),fmax),f))                                    2768
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         2769
      if(jerr .eq. 0)goto 17251                                            2769
      jerr=jerr-ilm                                                        2769
      go to 12180                                                          2769
17251 continue                                                             2770
      ix=0                                                                 2771
17260 do 17261 j=1,nin                                                     2771
      k=m(j)                                                               2772
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 17261                           2772
      ix=1                                                                 2772
      goto 17262                                                           2772
17261 continue                                                             2773
17262 continue                                                             2773
      if(ix .ne. 0)goto 17281                                              2774
17290 do 17291 k=1,ni                                                      2774
      if(ixx(k).eq.1)goto 17291                                            2774
      if(ju(k).eq.0)goto 17291                                             2775
      ga(k)=abs(dot_product(wr,x(:,k)))                                    2776
      if(ga(k) .le. sa*vp(k))goto 17311                                    2776
      ixx(k)=1                                                             2776
      ix=1                                                                 2776
17311 continue                                                             2777
17291 continue                                                             2778
17292 continue                                                             2778
      if(ix.eq.1) go to 10880                                              2779
      goto 17032                                                           2780
17281 continue                                                             2781
      goto 17031                                                           2782
17032 continue                                                             2782
      if(nin .le. nx)goto 17331                                            2782
      jerr=-10000-ilm                                                      2782
      goto 16952                                                           2782
17331 continue                                                             2783
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               2783
      kin(ilm)=nin                                                         2784
      alm(ilm)=al                                                          2784
      lmu=ilm                                                              2785
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   2786
      if(ilm.lt.mnl)goto 16951                                             2786
      if(flmin.ge.1.0)goto 16951                                           2787
      me=0                                                                 2787
17340 do 17341 j=1,nin                                                     2787
      if(ao(j,ilm).ne.0.0) me=me+1                                         2787
17341 continue                                                             2787
17342 continue                                                             2787
      if(me.gt.ne)goto 16952                                               2788
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 16952              2789
      if(dev(ilm).gt.devmax)goto 16952                                     2790
16951 continue                                                             2791
16952 continue                                                             2791
      g=f                                                                  2792
12180 continue                                                             2792
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)              2793
      return                                                               2794
      end                                                                  2795
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 2796
      implicit double precision(a-h,o-z)                                   2797
      double precision ca(nin),x(n,*),f(n)                                 2797
      integer ia(nin)                                                      2798
      f=0.0                                                                2798
      if(nin.le.0) return                                                  2799
17350 do 17351 i=1,n                                                       2799
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      2799
17351 continue                                                             2800
17352 continue                                                             2800
      return                                                               2801
      end                                                                  2802
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         2803
      implicit double precision(a-h,o-z)                                   2804
      double precision y(no),d(no),q(no)                                   2804
      integer jp(no),kp(*)                                                 2805
17360 do 17361 j=1,no                                                      2805
      jp(j)=j                                                              2805
17361 continue                                                             2805
17362 continue                                                             2805
      call psort7(y,jp,1,no)                                               2806
      nj=0                                                                 2806
17370 do 17371 j=1,no                                                      2806
      if(q(jp(j)).le.0.0)goto 17371                                        2806
      nj=nj+1                                                              2806
      jp(nj)=jp(j)                                                         2806
17371 continue                                                             2807
17372 continue                                                             2807
      if(nj .ne. 0)goto 17391                                              2807
      jerr=20000                                                           2807
      return                                                               2807
17391 continue                                                             2808
      j=1                                                                  2808
17400 continue                                                             2808
17401 if(d(jp(j)).gt.0.0)goto 17402                                        2808
      j=j+1                                                                2808
      if(j.gt.nj)goto 17402                                                2808
      goto 17401                                                           2809
17402 continue                                                             2809
      if(j .lt. nj-1)goto 17421                                            2809
      jerr=30000                                                           2809
      return                                                               2809
17421 continue                                                             2810
      t0=y(jp(j))                                                          2810
      j0=j-1                                                               2811
      if(j0 .le. 0)goto 17441                                              2812
17450 continue                                                             2812
17451 if(y(jp(j0)).lt.t0)goto 17452                                        2812
      j0=j0-1                                                              2812
      if(j0.eq.0)goto 17452                                                2812
      goto 17451                                                           2813
17452 continue                                                             2813
      if(j0 .le. 0)goto 17471                                              2813
      nj=nj-j0                                                             2813
17480 do 17481 j=1,nj                                                      2813
      jp(j)=jp(j+j0)                                                       2813
17481 continue                                                             2813
17482 continue                                                             2813
17471 continue                                                             2814
17441 continue                                                             2815
      jerr=0                                                               2815
      nk=0                                                                 2815
      yk=t0                                                                2815
      j=2                                                                  2816
17490 continue                                                             2816
17491 continue                                                             2816
17500 continue                                                             2817
17501 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 17502                     2817
      j=j+1                                                                2817
      if(j.gt.nj)goto 17502                                                2817
      goto 17501                                                           2818
17502 continue                                                             2818
      nk=nk+1                                                              2818
      kp(nk)=j-1                                                           2818
      if(j.gt.nj)goto 17492                                                2819
      if(j .ne. nj)goto 17521                                              2819
      nk=nk+1                                                              2819
      kp(nk)=nj                                                            2819
      goto 17492                                                           2819
17521 continue                                                             2820
      yk=y(jp(j))                                                          2820
      j=j+1                                                                2821
      goto 17491                                                           2822
17492 continue                                                             2822
      return                                                               2823
      end                                                                  2824
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     2825
      implicit double precision(a-h,o-z)                                   2826
      double precision d(no),dk(nk),wr(no),w(no)                           2827
      double precision e(no),u(no),b,c                                     2827
      integer kp(nk),jp(no)                                                2828
      call usk(no,nk,kp,jp,e,u)                                            2829
      b=dk(1)/u(1)                                                         2829
      c=dk(1)/u(1)**2                                                      2829
      jerr=0                                                               2830
17530 do 17531 j=1,kp(1)                                                   2830
      i=jp(j)                                                              2831
      w(i)=e(i)*(b-e(i)*c)                                                 2831
      if(w(i) .gt. 0.0)goto 17551                                          2831
      jerr=-30000                                                          2831
      return                                                               2831
17551 continue                                                             2832
      wr(i)=d(i)-e(i)*b                                                    2833
17531 continue                                                             2834
17532 continue                                                             2834
17560 do 17561 k=2,nk                                                      2834
      j1=kp(k-1)+1                                                         2834
      j2=kp(k)                                                             2835
      b=b+dk(k)/u(k)                                                       2835
      c=c+dk(k)/u(k)**2                                                    2836
17570 do 17571 j=j1,j2                                                     2836
      i=jp(j)                                                              2837
      w(i)=e(i)*(b-e(i)*c)                                                 2837
      if(w(i) .gt. 0.0)goto 17591                                          2837
      jerr=-30000                                                          2837
      return                                                               2837
17591 continue                                                             2838
      wr(i)=d(i)-e(i)*b                                                    2839
17571 continue                                                             2840
17572 continue                                                             2840
17561 continue                                                             2841
17562 continue                                                             2841
      return                                                               2842
      end                                                                  2843
      subroutine vars(no,ni,x,w,ixx,v)                                     2844
      implicit double precision(a-h,o-z)                                   2845
      double precision x(no,ni),w(no),v(ni)                                2845
      integer ixx(ni)                                                      2846
17600 do 17601 j=1,ni                                                      2846
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                        2846
17601 continue                                                             2847
17602 continue                                                             2847
      return                                                               2848
      end                                                                  2849
      subroutine died(no,nk,d,kp,jp,dk)                                    2850
      implicit double precision(a-h,o-z)                                   2851
      double precision d(no),dk(nk)                                        2851
      integer kp(nk),jp(no)                                                2852
      dk(1)=sum(d(jp(1:kp(1))))                                            2853
17610 do 17611 k=2,nk                                                      2853
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  2853
17611 continue                                                             2854
17612 continue                                                             2854
      return                                                               2855
      end                                                                  2856
      subroutine usk(no,nk,kp,jp,e,u)                                      2857
      implicit double precision(a-h,o-z)                                   2858
      double precision e(no),u(nk),h                                       2858
      integer kp(nk),jp(no)                                                2859
      h=0.0                                                                2860
17620 do 17621 k=nk,1,-1                                                   2860
      j2=kp(k)                                                             2861
      j1=1                                                                 2861
      if(k.gt.1) j1=kp(k-1)+1                                              2862
17630 do 17631 j=j2,j1,-1                                                  2862
      h=h+e(jp(j))                                                         2862
17631 continue                                                             2863
17632 continue                                                             2863
      u(k)=h                                                               2864
17621 continue                                                             2865
17622 continue                                                             2865
      return                                                               2866
      end                                                                  2867
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             2868
      implicit double precision(a-h,o-z)                                   2869
      double precision d(no),dk(nk),f(no)                                  2870
      integer kp(nk),jp(no)                                                2870
      double precision e(no),u(nk),s                                       2871
      call usk(no,nk,kp,jp,e,u)                                            2871
      u=log(u)                                                             2872
      risk=dot_product(d,f)-dot_product(dk,u)                              2873
      return                                                               2874
      end                                                                  2875
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 2876
      implicit double precision(a-h,o-z)                                   2877
      double precision x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(   2878 
     *nlam)
      double precision, dimension (:), allocatable :: dk,f,xm,dq,q              
      double precision, dimension (:), allocatable :: e,uu                      
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          2884
      if(jerr.ne.0) go to 12180                                            2885
      allocate(q(1:no),stat=jerr)                                          2886
      if(jerr.ne.0) go to 12180                                            2887
      allocate(uu(1:no),stat=jerr)                                         2888
      if(jerr.ne.0) go to 12180                                            2889
      allocate(f(1:no),stat=jerr)                                          2890
      if(jerr.ne.0) go to 12180                                            2891
      allocate(dk(1:no),stat=jerr)                                         2892
      if(jerr.ne.0) go to 12180                                            2893
      allocate(jp(1:no),stat=jerr)                                         2894
      if(jerr.ne.0) go to 12180                                            2895
      allocate(kp(1:no),stat=jerr)                                         2896
      if(jerr.ne.0) go to 12180                                            2897
      allocate(dq(1:no),stat=jerr)                                         2898
      if(jerr.ne.0) go to 12180                                            2899
      allocate(xm(1:ni),stat=jerr)                                         2900
      if(jerr.ne.0) go to 12180                                            2901
      q=max(0d0,w)                                                         2901
      sw=sum(q)                                                            2902
      if(sw .gt. 0.0)goto 17651                                            2902
      jerr=9999                                                            2902
      go to 12180                                                          2902
17651 continue                                                             2903
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               2904
      if(jerr.ne.0) go to 12180                                            2904
      fmax=log(huge(e(1))*0.1)                                             2905
      dq=d*q                                                               2905
      call died(no,nk,dq,kp,jp,dk)                                         2905
      gm=dot_product(q,g)/sw                                               2906
17660 do 17661 j=1,ni                                                      2906
      xm(j)=dot_product(q,x(:,j))/sw                                       2906
17661 continue                                                             2907
17662 continue                                                             2907
17670 do 17671 lam=1,nlam                                                  2908
17680 do 17681 i=1,no                                                      2908
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       2909
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        2910
17681 continue                                                             2911
17682 continue                                                             2911
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          2912
17671 continue                                                             2913
17672 continue                                                             2913
12180 continue                                                             2913
      deallocate(e,uu,dk,f,jp,kp,dq)                                       2914
      return                                                               2915
      end                                                                  2916
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,cl,ne,nx,nlam,flmin,u   2918 
     *lam,thr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2919
      double precision x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)        2920
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   2921
      integer jd(*),ia(nx),nin(nlam)                                       2922
      double precision, dimension (:), allocatable :: xm,xs,ww,vq               
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 17701                                    2926
      jerr=10000                                                           2926
      return                                                               2926
17701 continue                                                             2927
      if(minval(y) .ge. 0.0)goto 17721                                     2927
      jerr=8888                                                            2927
      return                                                               2927
17721 continue                                                             2928
      allocate(ww(1:no),stat=jerr)                                         2929
      if(jerr.ne.0) return                                                 2930
      allocate(ju(1:ni),stat=jerr)                                         2931
      if(jerr.ne.0) return                                                 2932
      allocate(vq(1:ni),stat=jerr)                                         2933
      if(jerr.ne.0) return                                                 2934
      allocate(xm(1:ni),stat=jerr)                                         2935
      if(jerr.ne.0) return                                                 2936
      if(isd .le. 0)goto 17741                                             2936
      allocate(xs(1:ni),stat=jerr)                                         2936
      if(jerr.ne.0) return                                                 2936
17741 continue                                                             2937
      call chkvars(no,ni,x,ju)                                             2938
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2939
      if(maxval(ju) .gt. 0)goto 17761                                      2939
      jerr=7777                                                            2939
      go to 12180                                                          2939
17761 continue                                                             2940
      vq=max(0d0,vp)                                                       2940
      vq=vq*ni/sum(vq)                                                     2941
      ww=max(0d0,w)                                                        2941
      sw=sum(ww)                                                           2941
      if(sw .gt. 0.0)goto 17781                                            2941
      jerr=9999                                                            2941
      go to 12180                                                          2941
17781 continue                                                             2942
      ww=ww/sw                                                             2943
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        2944
      if(isd .le. 0)goto 17801                                             2944
17810 do 17811 j=1,ni                                                      2944
      cl(:,j)=cl(:,j)*xs(j)                                                2944
17811 continue                                                             2944
17812 continue                                                             2944
17801 continue                                                             2945
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,t   2947 
     *hr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 12180                                            2947
      dev0=2.0*sw*dev0                                                     2948
17820 do 17821 k=1,lmu                                                     2948
      nk=nin(k)                                                            2949
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2950
      if(intr .ne. 0)goto 17841                                            2950
      a0(k)=0.0                                                            2950
      goto 17851                                                           2951
17841 continue                                                             2951
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2951
17851 continue                                                             2952
17831 continue                                                             2952
17821 continue                                                             2953
17822 continue                                                             2953
12180 continue                                                             2953
      deallocate(ww,ju,vq,xm)                                              2953
      if(isd.gt.0) deallocate(xs)                                          2954
      return                                                               2955
      end                                                                  2956
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,cl,ne,nx,nlam,flmin,u   2958 
     *lam,shri,  isd,intr,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2959
      double precision x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)        2960
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   2961
      integer ju(ni),m(nx),kin(nlam)                                       2962
      double precision, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga        
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2966
      sml=sml*10.0                                                         2967
      allocate(a(1:ni),stat=jerr)                                          2968
      if(jerr.ne.0) return                                                 2969
      allocate(as(1:ni),stat=jerr)                                         2970
      if(jerr.ne.0) return                                                 2971
      allocate(t(1:no),stat=jerr)                                          2972
      if(jerr.ne.0) return                                                 2973
      allocate(mm(1:ni),stat=jerr)                                         2974
      if(jerr.ne.0) return                                                 2975
      allocate(ga(1:ni),stat=jerr)                                         2976
      if(jerr.ne.0) return                                                 2977
      allocate(ixx(1:ni),stat=jerr)                                        2978
      if(jerr.ne.0) return                                                 2979
      allocate(wr(1:no),stat=jerr)                                         2980
      if(jerr.ne.0) return                                                 2981
      allocate(v(1:ni),stat=jerr)                                          2982
      if(jerr.ne.0) return                                                 2983
      allocate(w(1:no),stat=jerr)                                          2984
      if(jerr.ne.0) return                                                 2985
      allocate(f(1:no),stat=jerr)                                          2986
      if(jerr.ne.0) return                                                 2987
      bta=parm                                                             2987
      omb=1.0-bta                                                          2988
      t=q*y                                                                2988
      yb=sum(t)                                                            2988
      fmax=log(huge(bta)*0.1)                                              2989
      if(nonzero(no,g) .ne. 0)goto 17871                                   2990
      if(intr .eq. 0)goto 17891                                            2990
      w=q*yb                                                               2990
      az=log(yb)                                                           2990
      f=az                                                                 2990
      dv0=yb*(az-1.0)                                                      2990
      goto 17901                                                           2991
17891 continue                                                             2991
      w=q                                                                  2991
      az=0.0                                                               2991
      f=az                                                                 2991
      dv0=-1.0                                                             2991
17901 continue                                                             2992
17881 continue                                                             2992
      goto 17911                                                           2993
17871 continue                                                             2993
      w=q*exp(sign(min(abs(g),fmax),g))                                    2993
      v0=sum(w)                                                            2994
      if(intr .eq. 0)goto 17931                                            2994
      eaz=yb/v0                                                            2994
      w=eaz*w                                                              2994
      az=log(eaz)                                                          2994
      f=az+g                                                               2995
      dv0=dot_product(t,g)-yb*(1.0-az)                                     2996
      goto 17941                                                           2997
17931 continue                                                             2997
      az=0.0                                                               2997
      f=g                                                                  2997
      dv0=dot_product(t,g)-v0                                              2997
17941 continue                                                             2998
17921 continue                                                             2998
17911 continue                                                             2999
17861 continue                                                             2999
      a=0.0                                                                2999
      as=0.0                                                               2999
      wr=t-w                                                               2999
      v0=1.0                                                               2999
      if(intr.ne.0) v0=yb                                                  2999
      dvr=-yb                                                              3000
17950 do 17951 i=1,no                                                      3000
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               3000
17951 continue                                                             3000
17952 continue                                                             3000
      dvr=dvr-dv0                                                          3000
      dev0=dvr                                                             3002
      alf=1.0                                                              3004
      if(flmin .ge. 1.0)goto 17971                                         3004
      eqs=max(eps,flmin)                                                   3004
      alf=eqs**(1.0/(nlam-1))                                              3004
17971 continue                                                             3005
      m=0                                                                  3005
      mm=0                                                                 3005
      nlp=0                                                                3005
      nin=nlp                                                              3005
      mnl=min(mnlam,nlam)                                                  3005
      shr=shri*dev0                                                        3005
      ixx=0                                                                3005
      al=0.0                                                               3006
17980 do 17981 j=1,ni                                                      3006
      if(ju(j).eq.0)goto 17981                                             3006
      ga(j)=abs(dot_product(wr,x(:,j)))                                    3006
17981 continue                                                             3007
17982 continue                                                             3007
17990 do 17991 ilm=1,nlam                                                  3007
      al0=al                                                               3008
      if(flmin .lt. 1.0)goto 18011                                         3008
      al=ulam(ilm)                                                         3008
      goto 18001                                                           3009
18011 if(ilm .le. 2)goto 18021                                             3009
      al=al*alf                                                            3009
      goto 18001                                                           3010
18021 if(ilm .ne. 1)goto 18031                                             3010
      al=big                                                               3010
      goto 18041                                                           3011
18031 continue                                                             3011
      al0=0.0                                                              3012
18050 do 18051 j=1,ni                                                      3012
      if(ju(j).eq.0)goto 18051                                             3012
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3012
18051 continue                                                             3013
18052 continue                                                             3013
      al0=al0/max(bta,1.0d-3)                                              3013
      al=alf*al0                                                           3014
18041 continue                                                             3015
18001 continue                                                             3015
      al2=al*omb                                                           3015
      al1=al*bta                                                           3015
      tlam=bta*(2.0*al-al0)                                                3016
18060 do 18061 k=1,ni                                                      3016
      if(ixx(k).eq.1)goto 18061                                            3016
      if(ju(k).eq.0)goto 18061                                             3017
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     3018
18061 continue                                                             3019
18062 continue                                                             3019
10880 continue                                                             3020
18070 continue                                                             3020
18071 continue                                                             3020
      az0=az                                                               3021
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                3022
18080 do 18081 j=1,ni                                                      3022
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                        3022
18081 continue                                                             3023
18082 continue                                                             3023
18090 continue                                                             3023
18091 continue                                                             3023
      nlp=nlp+1                                                            3023
      dlx=0.0                                                              3024
18100 do 18101 k=1,ni                                                      3024
      if(ixx(k).eq.0)goto 18101                                            3024
      ak=a(k)                                                              3025
      u=dot_product(wr,x(:,k))+v(k)*ak                                     3025
      au=abs(u)-vp(k)*al1                                                  3026
      if(au .gt. 0.0)goto 18121                                            3026
      a(k)=0.0                                                             3026
      goto 18131                                                           3027
18121 continue                                                             3028
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3029
18131 continue                                                             3030
18111 continue                                                             3030
      if(a(k).eq.ak)goto 18101                                             3030
      d=a(k)-ak                                                            3030
      dlx=max(dlx,v(k)*d**2)                                               3031
      wr=wr-d*w*x(:,k)                                                     3031
      f=f+d*x(:,k)                                                         3032
      if(mm(k) .ne. 0)goto 18151                                           3032
      nin=nin+1                                                            3032
      if(nin.gt.nx)goto 18102                                              3033
      mm(k)=nin                                                            3033
      m(nin)=k                                                             3034
18151 continue                                                             3035
18101 continue                                                             3036
18102 continue                                                             3036
      if(nin.gt.nx)goto 18092                                              3037
      if(intr .eq. 0)goto 18171                                            3037
      d=sum(wr)/v0                                                         3038
      az=az+d                                                              3038
      dlx=max(dlx,v0*d**2)                                                 3038
      wr=wr-d*w                                                            3038
      f=f+d                                                                3039
18171 continue                                                             3040
      if(dlx.lt.shr)goto 18092                                             3040
      if(nlp .le. maxit)goto 18191                                         3040
      jerr=-ilm                                                            3040
      return                                                               3040
18191 continue                                                             3041
18200 continue                                                             3041
18201 continue                                                             3041
      nlp=nlp+1                                                            3041
      dlx=0.0                                                              3042
18210 do 18211 l=1,nin                                                     3042
      k=m(l)                                                               3042
      ak=a(k)                                                              3043
      u=dot_product(wr,x(:,k))+v(k)*ak                                     3043
      au=abs(u)-vp(k)*al1                                                  3044
      if(au .gt. 0.0)goto 18231                                            3044
      a(k)=0.0                                                             3044
      goto 18241                                                           3045
18231 continue                                                             3046
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3047
18241 continue                                                             3048
18221 continue                                                             3048
      if(a(k).eq.ak)goto 18211                                             3048
      d=a(k)-ak                                                            3048
      dlx=max(dlx,v(k)*d**2)                                               3049
      wr=wr-d*w*x(:,k)                                                     3049
      f=f+d*x(:,k)                                                         3051
18211 continue                                                             3051
18212 continue                                                             3051
      if(intr .eq. 0)goto 18261                                            3051
      d=sum(wr)/v0                                                         3051
      az=az+d                                                              3052
      dlx=max(dlx,v0*d**2)                                                 3052
      wr=wr-d*w                                                            3052
      f=f+d                                                                3053
18261 continue                                                             3054
      if(dlx.lt.shr)goto 18202                                             3054
      if(nlp .le. maxit)goto 18281                                         3054
      jerr=-ilm                                                            3054
      return                                                               3054
18281 continue                                                             3055
      goto 18201                                                           3056
18202 continue                                                             3056
      goto 18091                                                           3057
18092 continue                                                             3057
      if(nin.gt.nx)goto 18072                                              3058
      w=q*exp(sign(min(abs(f),fmax),f))                                    3058
      v0=sum(w)                                                            3058
      wr=t-w                                                               3059
      if(v0*(az-az0)**2 .ge. shr)goto 18301                                3059
      ix=0                                                                 3060
18310 do 18311 j=1,nin                                                     3060
      k=m(j)                                                               3061
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 18311                            3061
      ix=1                                                                 3061
      goto 18312                                                           3062
18311 continue                                                             3063
18312 continue                                                             3063
      if(ix .ne. 0)goto 18331                                              3064
18340 do 18341 k=1,ni                                                      3064
      if(ixx(k).eq.1)goto 18341                                            3064
      if(ju(k).eq.0)goto 18341                                             3065
      ga(k)=abs(dot_product(wr,x(:,k)))                                    3066
      if(ga(k) .le. al1*vp(k))goto 18361                                   3066
      ixx(k)=1                                                             3066
      ix=1                                                                 3066
18361 continue                                                             3067
18341 continue                                                             3068
18342 continue                                                             3068
      if(ix.eq.1) go to 10880                                              3069
      goto 18072                                                           3070
18331 continue                                                             3071
18301 continue                                                             3072
      goto 18071                                                           3073
18072 continue                                                             3073
      if(nin .le. nx)goto 18381                                            3073
      jerr=-10000-ilm                                                      3073
      goto 17992                                                           3073
18381 continue                                                             3074
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               3074
      kin(ilm)=nin                                                         3075
      a0(ilm)=az                                                           3075
      alm(ilm)=al                                                          3075
      lmu=ilm                                                              3076
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               3077
      if(ilm.lt.mnl)goto 17991                                             3077
      if(flmin.ge.1.0)goto 17991                                           3078
      me=0                                                                 3078
18390 do 18391 j=1,nin                                                     3078
      if(ca(j,ilm).ne.0.0) me=me+1                                         3078
18391 continue                                                             3078
18392 continue                                                             3078
      if(me.gt.ne)goto 17992                                               3079
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17992              3080
      if(dev(ilm).gt.devmax)goto 17992                                     3081
17991 continue                                                             3082
17992 continue                                                             3082
      g=f                                                                  3083
12180 continue                                                             3083
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                                3084
      return                                                               3085
      end                                                                  3086
      function nonzero(n,v)                                                3087
      implicit double precision(a-h,o-z)                                   3088
      double precision v(n)                                                3089
      nonzero=0                                                            3089
18400 do 18401 i=1,n                                                       3089
      if(v(i) .eq. 0.0)goto 18421                                          3089
      nonzero=1                                                            3089
      return                                                               3089
18421 continue                                                             3089
18401 continue                                                             3090
18402 continue                                                             3090
      return                                                               3091
      end                                                                  3092
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               3093
      implicit double precision(a-h,o-z)                                   3094
      double precision a(nx,lmu),b(ni,lmu)                                 3094
      integer ia(nx),nin(lmu)                                              3095
18430 do 18431 lam=1,lmu                                                   3095
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        3095
18431 continue                                                             3096
18432 continue                                                             3096
      return                                                               3097
      end                                                                  3098
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           3099
      implicit double precision(a-h,o-z)                                   3100
      double precision a(nx,nc,lmu),b(ni,nc,lmu)                           3100
      integer ia(nx),nin(lmu)                                              3101
18440 do 18441 lam=1,lmu                                                   3101
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             3101
18441 continue                                                             3102
18442 continue                                                             3102
      return                                                               3103
      end                                                                  3104
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               3105
      implicit double precision(a-h,o-z)                                   3106
      double precision x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),fl   3107 
     *og(nlam)
      double precision, dimension (:), allocatable :: w                         
      if(minval(y) .ge. 0.0)goto 18461                                     3110
      jerr=8888                                                            3110
      return                                                               3110
18461 continue                                                             3111
      allocate(w(1:no),stat=jerr)                                          3111
      if(jerr.ne.0) return                                                 3112
      w=max(0d0,q)                                                         3112
      sw=sum(w)                                                            3112
      if(sw .gt. 0.0)goto 18481                                            3112
      jerr=9999                                                            3112
      go to 12180                                                          3112
18481 continue                                                             3113
      yb=dot_product(w,y)/sw                                               3113
      fmax=log(huge(y(1))*0.1)                                             3114
18490 do 18491 lam=1,nlam                                                  3114
      s=0.0                                                                3115
18500 do 18501 i=1,no                                                      3115
      if(w(i).le.0.0)goto 18501                                            3116
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          3117
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      3118
18501 continue                                                             3119
18502 continue                                                             3119
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                3120
18491 continue                                                             3121
18492 continue                                                             3121
12180 continue                                                             3121
      deallocate(w)                                                        3122
      return                                                               3123
      end                                                                  3124
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,cl,ne,nx,nlam   3126 
     *,flmin,  ulam,thr,isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp
     *,jerr)
      implicit double precision(a-h,o-z)                                   3127
      double precision x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)   3128
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)            3129
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3130
      double precision, dimension (:), allocatable :: xm,xs,ww,vq               
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 18521                                    3134
      jerr=10000                                                           3134
      return                                                               3134
18521 continue                                                             3135
      if(minval(y) .ge. 0.0)goto 18541                                     3135
      jerr=8888                                                            3135
      return                                                               3135
18541 continue                                                             3136
      allocate(ww(1:no),stat=jerr)                                         3137
      if(jerr.ne.0) return                                                 3138
      allocate(ju(1:ni),stat=jerr)                                         3139
      if(jerr.ne.0) return                                                 3140
      allocate(vq(1:ni),stat=jerr)                                         3141
      if(jerr.ne.0) return                                                 3142
      allocate(xm(1:ni),stat=jerr)                                         3143
      if(jerr.ne.0) return                                                 3144
      allocate(xs(1:ni),stat=jerr)                                         3145
      if(jerr.ne.0) return                                                 3146
      call spchkvars(no,ni,x,ix,ju)                                        3147
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3148
      if(maxval(ju) .gt. 0)goto 18561                                      3148
      jerr=7777                                                            3148
      go to 12180                                                          3148
18561 continue                                                             3149
      vq=max(0d0,vp)                                                       3149
      vq=vq*ni/sum(vq)                                                     3150
      ww=max(0d0,w)                                                        3150
      sw=sum(ww)                                                           3150
      if(sw .gt. 0.0)goto 18581                                            3150
      jerr=9999                                                            3150
      go to 12180                                                          3150
18581 continue                                                             3151
      ww=ww/sw                                                             3152
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                3153
      if(isd .le. 0)goto 18601                                             3153
18610 do 18611 j=1,ni                                                      3153
      cl(:,j)=cl(:,j)*xs(j)                                                3153
18611 continue                                                             3153
18612 continue                                                             3153
18601 continue                                                             3154
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,flmi   3156 
     *n,ulam,thr,  isd,intr,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nl
     *p,jerr)
      if(jerr.gt.0) go to 12180                                            3156
      dev0=2.0*sw*dev0                                                     3157
18620 do 18621 k=1,lmu                                                     3157
      nk=nin(k)                                                            3158
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      3159
      if(intr .ne. 0)goto 18641                                            3159
      a0(k)=0.0                                                            3159
      goto 18651                                                           3160
18641 continue                                                             3160
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     3160
18651 continue                                                             3161
18631 continue                                                             3161
18621 continue                                                             3162
18622 continue                                                             3162
12180 continue                                                             3162
      deallocate(ww,ju,vq,xm,xs)                                           3163
      return                                                               3164
      end                                                                  3165
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,cl,ne,nx,nlam   3167 
     *,flmin,ulam,  shri,isd,intr,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,a
     *lm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   3168
      double precision x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),x   3169 
     *s(ni)
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   3170
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           3171
      double precision, dimension (:), allocatable :: qy,t,w,wr,v               
      double precision, dimension (:), allocatable :: a,as,xm,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               3176
      sml=sml*10.0                                                         3177
      allocate(a(1:ni),stat=jerr)                                          3178
      if(jerr.ne.0) return                                                 3179
      allocate(as(1:ni),stat=jerr)                                         3180
      if(jerr.ne.0) return                                                 3181
      allocate(t(1:no),stat=jerr)                                          3182
      if(jerr.ne.0) return                                                 3183
      allocate(mm(1:ni),stat=jerr)                                         3184
      if(jerr.ne.0) return                                                 3185
      allocate(ga(1:ni),stat=jerr)                                         3186
      if(jerr.ne.0) return                                                 3187
      allocate(ixx(1:ni),stat=jerr)                                        3188
      if(jerr.ne.0) return                                                 3189
      allocate(wr(1:no),stat=jerr)                                         3190
      if(jerr.ne.0) return                                                 3191
      allocate(v(1:ni),stat=jerr)                                          3192
      if(jerr.ne.0) return                                                 3193
      allocate(xm(1:ni),stat=jerr)                                         3194
      if(jerr.ne.0) return                                                 3195
      allocate(w(1:no),stat=jerr)                                          3196
      if(jerr.ne.0) return                                                 3197
      allocate(qy(1:no),stat=jerr)                                         3198
      if(jerr.ne.0) return                                                 3199
      bta=parm                                                             3199
      omb=1.0-bta                                                          3199
      fmax=log(huge(bta)*0.1)                                              3200
      qy=q*y                                                               3200
      yb=sum(qy)                                                           3201
      if(nonzero(no,g) .ne. 0)goto 18671                                   3201
      t=0.0                                                                3202
      if(intr .eq. 0)goto 18691                                            3202
      w=q*yb                                                               3202
      az=log(yb)                                                           3202
      uu=az                                                                3203
      xm=yb*xb                                                             3203
      dv0=yb*(az-1.0)                                                      3204
      goto 18701                                                           3205
18691 continue                                                             3205
      w=q                                                                  3205
      xm=0.0                                                               3205
      uu=0.0                                                               3205
      az=uu                                                                3205
      dv0=-1.0                                                             3205
18701 continue                                                             3206
18681 continue                                                             3206
      goto 18711                                                           3207
18671 continue                                                             3207
      w=q*exp(sign(min(abs(g),fmax),g))                                    3207
      ww=sum(w)                                                            3207
      t=g                                                                  3208
      if(intr .eq. 0)goto 18731                                            3208
      eaz=yb/ww                                                            3209
      w=eaz*w                                                              3209
      az=log(eaz)                                                          3209
      uu=az                                                                3209
      dv0=dot_product(qy,g)-yb*(1.0-az)                                    3210
      goto 18741                                                           3211
18731 continue                                                             3211
      uu=0.0                                                               3211
      az=uu                                                                3211
      dv0=dot_product(qy,g)-ww                                             3211
18741 continue                                                             3212
18721 continue                                                             3212
18750 do 18751 j=1,ni                                                      3212
      if(ju(j).eq.0)goto 18751                                             3212
      jb=ix(j)                                                             3212
      je=ix(j+1)-1                                                         3213
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3214
18751 continue                                                             3215
18752 continue                                                             3215
18711 continue                                                             3216
18661 continue                                                             3216
      tt=yb*uu                                                             3216
      ww=1.0                                                               3216
      if(intr.ne.0) ww=yb                                                  3216
      wr=qy-q*(yb*(1.0-uu))                                                3216
      a=0.0                                                                3216
      as=0.0                                                               3217
      dvr=-yb                                                              3218
18760 do 18761 i=1,no                                                      3218
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             3218
18761 continue                                                             3218
18762 continue                                                             3218
      dvr=dvr-dv0                                                          3218
      dev0=dvr                                                             3220
      alf=1.0                                                              3222
      if(flmin .ge. 1.0)goto 18781                                         3222
      eqs=max(eps,flmin)                                                   3222
      alf=eqs**(1.0/(nlam-1))                                              3222
18781 continue                                                             3223
      m=0                                                                  3223
      mm=0                                                                 3223
      nlp=0                                                                3223
      nin=nlp                                                              3223
      mnl=min(mnlam,nlam)                                                  3223
      shr=shri*dev0                                                        3223
      al=0.0                                                               3223
      ixx=0                                                                3224
18790 do 18791 j=1,ni                                                      3224
      if(ju(j).eq.0)goto 18791                                             3225
      jb=ix(j)                                                             3225
      je=ix(j+1)-1                                                         3226
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   3228 
     *)-xb(j)*tt)/xs(j)
18791 continue                                                             3229
18792 continue                                                             3229
18800 do 18801 ilm=1,nlam                                                  3229
      al0=al                                                               3230
      if(flmin .lt. 1.0)goto 18821                                         3230
      al=ulam(ilm)                                                         3230
      goto 18811                                                           3231
18821 if(ilm .le. 2)goto 18831                                             3231
      al=al*alf                                                            3231
      goto 18811                                                           3232
18831 if(ilm .ne. 1)goto 18841                                             3232
      al=big                                                               3232
      goto 18851                                                           3233
18841 continue                                                             3233
      al0=0.0                                                              3234
18860 do 18861 j=1,ni                                                      3234
      if(ju(j).eq.0)goto 18861                                             3234
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3234
18861 continue                                                             3235
18862 continue                                                             3235
      al0=al0/max(bta,1.0d-3)                                              3235
      al=alf*al0                                                           3236
18851 continue                                                             3237
18811 continue                                                             3237
      al2=al*omb                                                           3237
      al1=al*bta                                                           3237
      tlam=bta*(2.0*al-al0)                                                3238
18870 do 18871 k=1,ni                                                      3238
      if(ixx(k).eq.1)goto 18871                                            3238
      if(ju(k).eq.0)goto 18871                                             3239
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     3240
18871 continue                                                             3241
18872 continue                                                             3241
10880 continue                                                             3242
18880 continue                                                             3242
18881 continue                                                             3242
      az0=az                                                               3243
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                3244
18890 do 18891 j=1,ni                                                      3244
      if(ixx(j).eq.0)goto 18891                                            3244
      jb=ix(j)                                                             3244
      je=ix(j+1)-1                                                         3245
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3246
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   3248 
     *b(j)**2)/xs(j)**2
18891 continue                                                             3249
18892 continue                                                             3249
18900 continue                                                             3249
18901 continue                                                             3249
      nlp=nlp+1                                                            3250
      dlx=0.0                                                              3251
18910 do 18911 k=1,ni                                                      3251
      if(ixx(k).eq.0)goto 18911                                            3251
      jb=ix(k)                                                             3251
      je=ix(k+1)-1                                                         3251
      ak=a(k)                                                              3252
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   3254 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  3255
      if(au .gt. 0.0)goto 18931                                            3255
      a(k)=0.0                                                             3255
      goto 18941                                                           3256
18931 continue                                                             3257
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3258
18941 continue                                                             3259
18921 continue                                                             3259
      if(a(k).eq.ak)goto 18911                                             3260
      if(mm(k) .ne. 0)goto 18961                                           3260
      nin=nin+1                                                            3260
      if(nin.gt.nx)goto 18912                                              3261
      mm(k)=nin                                                            3261
      m(nin)=k                                                             3262
18961 continue                                                             3263
      d=a(k)-ak                                                            3263
      dlx=max(dlx,v(k)*d**2)                                               3263
      dv=d/xs(k)                                                           3264
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 3265
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                3266
      uu=uu-dv*xb(k)                                                       3266
      tt=tt-dv*xm(k)                                                       3267
18911 continue                                                             3268
18912 continue                                                             3268
      if(nin.gt.nx)goto 18902                                              3269
      if(intr .eq. 0)goto 18981                                            3269
      d=tt/ww-uu                                                           3270
      az=az+d                                                              3270
      dlx=max(dlx,ww*d**2)                                                 3270
      uu=uu+d                                                              3271
18981 continue                                                             3272
      if(dlx.lt.shr)goto 18902                                             3272
      if(nlp .le. maxit)goto 19001                                         3272
      jerr=-ilm                                                            3272
      return                                                               3272
19001 continue                                                             3273
19010 continue                                                             3273
19011 continue                                                             3273
      nlp=nlp+1                                                            3273
      dlx=0.0                                                              3274
19020 do 19021 l=1,nin                                                     3274
      k=m(l)                                                               3275
      jb=ix(k)                                                             3275
      je=ix(k+1)-1                                                         3275
      ak=a(k)                                                              3276
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   3278 
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  3279
      if(au .gt. 0.0)goto 19041                                            3279
      a(k)=0.0                                                             3279
      goto 19051                                                           3280
19041 continue                                                             3281
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3282
19051 continue                                                             3283
19031 continue                                                             3283
      if(a(k).eq.ak)goto 19021                                             3283
      d=a(k)-ak                                                            3283
      dlx=max(dlx,v(k)*d**2)                                               3284
      dv=d/xs(k)                                                           3284
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 3285
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                3286
      uu=uu-dv*xb(k)                                                       3286
      tt=tt-dv*xm(k)                                                       3287
19021 continue                                                             3288
19022 continue                                                             3288
      if(intr .eq. 0)goto 19071                                            3288
      d=tt/ww-uu                                                           3288
      az=az+d                                                              3289
      dlx=max(dlx,ww*d**2)                                                 3289
      uu=uu+d                                                              3290
19071 continue                                                             3291
      if(dlx.lt.shr)goto 19012                                             3291
      if(nlp .le. maxit)goto 19091                                         3291
      jerr=-ilm                                                            3291
      return                                                               3291
19091 continue                                                             3292
      goto 19011                                                           3293
19012 continue                                                             3293
      goto 18901                                                           3294
18902 continue                                                             3294
      if(nin.gt.nx)goto 18882                                              3295
      euu=exp(sign(min(abs(uu),fmax),uu))                                  3296
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                3296
      ww=sum(w)                                                            3297
      wr=qy-w*(1.0-uu)                                                     3297
      tt=sum(wr)                                                           3298
      if(ww*(az-az0)**2 .ge. shr)goto 19111                                3298
      kx=0                                                                 3299
19120 do 19121 j=1,nin                                                     3299
      k=m(j)                                                               3300
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 19121                            3300
      kx=1                                                                 3300
      goto 19122                                                           3301
19121 continue                                                             3302
19122 continue                                                             3302
      if(kx .ne. 0)goto 19141                                              3303
19150 do 19151 j=1,ni                                                      3303
      if(ixx(j).eq.1)goto 19151                                            3303
      if(ju(j).eq.0)goto 19151                                             3304
      jb=ix(j)                                                             3304
      je=ix(j+1)-1                                                         3305
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3306
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   3308 
     *)-xb(j)*tt)/xs(j)
      if(ga(j) .le. al1*vp(j))goto 19171                                   3308
      ixx(j)=1                                                             3308
      kx=1                                                                 3308
19171 continue                                                             3309
19151 continue                                                             3310
19152 continue                                                             3310
      if(kx.eq.1) go to 10880                                              3311
      goto 18882                                                           3312
19141 continue                                                             3313
19111 continue                                                             3314
      goto 18881                                                           3315
18882 continue                                                             3315
      if(nin .le. nx)goto 19191                                            3315
      jerr=-10000-ilm                                                      3315
      goto 18802                                                           3315
19191 continue                                                             3316
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               3316
      kin(ilm)=nin                                                         3317
      a0(ilm)=az                                                           3317
      alm(ilm)=al                                                          3317
      lmu=ilm                                                              3318
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        3319
      if(ilm.lt.mnl)goto 18801                                             3319
      if(flmin.ge.1.0)goto 18801                                           3320
      me=0                                                                 3320
19200 do 19201 j=1,nin                                                     3320
      if(ca(j,ilm).ne.0.0) me=me+1                                         3320
19201 continue                                                             3320
19202 continue                                                             3320
      if(me.gt.ne)goto 18802                                               3321
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 18802              3322
      if(dev(ilm).gt.devmax)goto 18802                                     3323
18801 continue                                                             3324
18802 continue                                                             3324
      g=t+uu                                                               3325
12180 continue                                                             3325
      deallocate(t,w,wr,v,a,qy,xm,as,mm,ga,ixx)                            3326
      return                                                               3327
      end                                                                  3328
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       3329
      implicit double precision(a-h,o-z)                                   3330
      double precision x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(n   3331 
     *lam)
      integer ix(*),jx(*)                                                  3332
      double precision, dimension (:), allocatable :: w,f                       
      if(minval(y) .ge. 0.0)goto 19221                                     3335
      jerr=8888                                                            3335
      return                                                               3335
19221 continue                                                             3336
      allocate(w(1:no),stat=jerr)                                          3337
      if(jerr.ne.0) return                                                 3338
      allocate(f(1:no),stat=jerr)                                          3339
      if(jerr.ne.0) return                                                 3340
      w=max(0d0,q)                                                         3340
      sw=sum(w)                                                            3340
      if(sw .gt. 0.0)goto 19241                                            3340
      jerr=9999                                                            3340
      go to 12180                                                          3340
19241 continue                                                             3341
      yb=dot_product(w,y)/sw                                               3341
      fmax=log(huge(y(1))*0.1)                                             3342
19250 do 19251 lam=1,nlam                                                  3342
      f=a0(lam)                                                            3343
19260 do 19261 j=1,ni                                                      3343
      if(a(j,lam).eq.0.0)goto 19261                                        3343
      jb=ix(j)                                                             3343
      je=ix(j+1)-1                                                         3344
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          3345
19261 continue                                                             3346
19262 continue                                                             3346
      f=f+g                                                                3347
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   3348
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                3349
19251 continue                                                             3350
19252 continue                                                             3350
12180 continue                                                             3350
      deallocate(w,f)                                                      3351
      return                                                               3352
      end                                                                  3353
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   3354 
     *jerr)
      implicit double precision(a-h,o-z)                                   3355
      double precision x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(   3356 
     *nlam)
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 3357
      double precision, dimension (:), allocatable :: w,f                       
      if(minval(y) .ge. 0.0)goto 19281                                     3360
      jerr=8888                                                            3360
      return                                                               3360
19281 continue                                                             3361
      allocate(w(1:no),stat=jerr)                                          3362
      if(jerr.ne.0) return                                                 3363
      allocate(f(1:no),stat=jerr)                                          3364
      if(jerr.ne.0) return                                                 3365
      w=max(0d0,q)                                                         3365
      sw=sum(w)                                                            3365
      if(sw .gt. 0.0)goto 19301                                            3365
      jerr=9999                                                            3365
      go to 12180                                                          3365
19301 continue                                                             3366
      yb=dot_product(w,y)/sw                                               3366
      fmax=log(huge(y(1))*0.1)                                             3367
19310 do 19311 lam=1,nlam                                                  3367
      f=a0(lam)                                                            3368
19320 do 19321 k=1,nin(lam)                                                3368
      j=ia(k)                                                              3368
      jb=ix(j)                                                             3368
      je=ix(j+1)-1                                                         3369
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         3370
19321 continue                                                             3371
19322 continue                                                             3371
      f=f+g                                                                3372
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   3373
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                3374
19311 continue                                                             3375
19312 continue                                                             3375
12180 continue                                                             3375
      deallocate(w,f)                                                      3376
      return                                                               3377
      end                                                                  3378
      subroutine multelnet  (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flm   3381 
     *in,ulam,thr,isd,jsd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      implicit double precision(a-h,o-z)                                   3382
      double precision x(no,ni),y(no,nr),w(no),vp(ni),ca(nx,nr,nlam)       3383
      double precision ulam(nlam),a0(nr,nlam),rsq(nlam),alm(nlam),cl(2,n   3384 
     *i)
      integer jd(*),ia(nx),nin(nlam)                                       3385
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 19341                                    3388
      jerr=10000                                                           3388
      return                                                               3388
19341 continue                                                             3389
      allocate(vq(1:ni),stat=jerr)                                         3389
      if(jerr.ne.0) return                                                 3390
      vq=max(0d0,vp)                                                       3390
      vq=vq*ni/sum(vq)                                                     3391
      call multelnetn(parm,no,ni,nr,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam   3393 
     *,thr,isd,  jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      deallocate(vq)                                                       3394
      return                                                               3395
      end                                                                  3396
      subroutine multelnetn (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flm   3398 
     *in,ulam,thr,  isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      implicit double precision(a-h,o-z)                                   3399
      double precision vp(ni),x(no,ni),y(no,nr),w(no),ulam(nlam),cl(2,ni   3400 
     *)
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)      3401
      integer jd(*),ia(nx),nin(nlam)                                       3402
      double precision, dimension (:), allocatable :: xm,xs,xv,ym,ys            
      integer, dimension (:), allocatable :: ju                                 
      double precision, dimension (:,:,:), allocatable :: clt                   
      allocate(clt(1:2,1:nr,1:ni),stat=jerr);                                   
      if(jerr.ne.0) return                                                      
      allocate(xm(1:ni),stat=jerr)                                         3410
      if(jerr.ne.0) return                                                 3411
      allocate(xs(1:ni),stat=jerr)                                         3412
      if(jerr.ne.0) return                                                 3413
      allocate(ym(1:nr),stat=jerr)                                         3414
      if(jerr.ne.0) return                                                 3415
      allocate(ys(1:nr),stat=jerr)                                         3416
      if(jerr.ne.0) return                                                 3417
      allocate(ju(1:ni),stat=jerr)                                         3418
      if(jerr.ne.0) return                                                 3419
      allocate(xv(1:ni),stat=jerr)                                         3420
      if(jerr.ne.0) return                                                 3421
      call chkvars(no,ni,x,ju)                                             3422
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3423
      if(maxval(ju) .gt. 0)goto 19361                                      3423
      jerr=7777                                                            3423
      return                                                               3423
19361 continue                                                             3424
      call multstandard1(no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym,ys,xv,y   3425 
     *s0,jerr)
      if(jerr.ne.0) return                                                 3426
19370 do 19371 j=1,ni                                                      3426
19380 do 19381 k=1,nr                                                      3426
19390 do 19391 i=1,2                                                       3426
      clt(i,k,j)=cl(i,j)                                                   3426
19391 continue                                                             3426
19392 continue                                                             3426
19381 continue                                                             3426
19382 continue                                                             3426
19371 continue                                                             3427
19372 continue                                                             3427
      if(isd .le. 0)goto 19411                                             3427
19420 do 19421 j=1,ni                                                      3427
19430 do 19431 k=1,nr                                                      3427
19440 do 19441 i=1,2                                                       3427
      clt(i,k,j)=clt(i,k,j)*xs(j)                                          3427
19441 continue                                                             3427
19442 continue                                                             3427
19431 continue                                                             3427
19432 continue                                                             3427
19421 continue                                                             3427
19422 continue                                                             3427
19411 continue                                                             3428
      if(jsd .le. 0)goto 19461                                             3428
19470 do 19471 j=1,ni                                                      3428
19480 do 19481 k=1,nr                                                      3428
19490 do 19491 i=1,2                                                       3428
      clt(i,k,j)=clt(i,k,j)/ys(k)                                          3428
19491 continue                                                             3428
19492 continue                                                             3428
19481 continue                                                             3428
19482 continue                                                             3428
19471 continue                                                             3428
19472 continue                                                             3428
19461 continue                                                             3429
      call multelnet2(parm,ni,nr,ju,vp,clt,y,no,ne,nx,x,nlam,flmin,ulam,   3431 
     *thr,maxit,xv,  ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 3432
19500 do 19501 k=1,lmu                                                     3432
      nk=nin(k)                                                            3433
19510 do 19511 j=1,nr                                                      3434
19520 do 19521 l=1,nk                                                      3434
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  3434
19521 continue                                                             3435
19522 continue                                                             3435
      if(intr .ne. 0)goto 19541                                            3435
      a0(j,k)=0.0                                                          3435
      goto 19551                                                           3436
19541 continue                                                             3436
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 3436
19551 continue                                                             3437
19531 continue                                                             3437
19511 continue                                                             3438
19512 continue                                                             3438
19501 continue                                                             3439
19502 continue                                                             3439
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                    3440
      return                                                               3441
      end                                                                  3442
      subroutine multstandard1  (no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym   3444 
     *,ys,xv,ys0,jerr)
      implicit double precision(a-h,o-z)                                   3445
      double precision x(no,ni),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(n   3446 
     *r),ys(nr)
      integer ju(ni)                                                       3447
      double precision, dimension (:), allocatable :: v                         
      allocate(v(1:no),stat=jerr)                                          3450
      if(jerr.ne.0) return                                                 3451
      w=w/sum(w)                                                           3451
      v=sqrt(w)                                                            3452
      if(intr .ne. 0)goto 19571                                            3453
19580 do 19581 j=1,ni                                                      3453
      if(ju(j).eq.0)goto 19581                                             3453
      xm(j)=0.0                                                            3453
      x(:,j)=v*x(:,j)                                                      3454
      z=dot_product(x(:,j),x(:,j))                                         3455
      if(isd .le. 0)goto 19601                                             3455
      xbq=dot_product(v,x(:,j))**2                                         3455
      vc=z-xbq                                                             3456
      xs(j)=sqrt(vc)                                                       3456
      x(:,j)=x(:,j)/xs(j)                                                  3456
      xv(j)=1.0+xbq/vc                                                     3457
      goto 19611                                                           3458
19601 continue                                                             3458
      xs(j)=1.0                                                            3458
      xv(j)=z                                                              3458
19611 continue                                                             3459
19591 continue                                                             3459
19581 continue                                                             3460
19582 continue                                                             3460
      ys0=0.0                                                              3461
19620 do 19621 j=1,nr                                                      3461
      ym(j)=0.0                                                            3461
      y(:,j)=v*y(:,j)                                                      3462
      z=dot_product(y(:,j),y(:,j))                                         3463
      if(jsd .le. 0)goto 19641                                             3463
      u=z-dot_product(v,y(:,j))**2                                         3463
      ys0=ys0+z/u                                                          3464
      ys(j)=sqrt(u)                                                        3464
      y(:,j)=y(:,j)/ys(j)                                                  3465
      goto 19651                                                           3466
19641 continue                                                             3466
      ys(j)=1.0                                                            3466
      ys0=ys0+z                                                            3466
19651 continue                                                             3467
19631 continue                                                             3467
19621 continue                                                             3468
19622 continue                                                             3468
      go to 10700                                                          3469
19571 continue                                                             3470
19660 do 19661 j=1,ni                                                      3470
      if(ju(j).eq.0)goto 19661                                             3471
      xm(j)=dot_product(w,x(:,j))                                          3471
      x(:,j)=v*(x(:,j)-xm(j))                                              3472
      xv(j)=dot_product(x(:,j),x(:,j))                                     3472
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       3473
19661 continue                                                             3474
19662 continue                                                             3474
      if(isd .ne. 0)goto 19681                                             3474
      xs=1.0                                                               3474
      goto 19691                                                           3475
19681 continue                                                             3475
19700 do 19701 j=1,ni                                                      3475
      if(ju(j).eq.0)goto 19701                                             3475
      x(:,j)=x(:,j)/xs(j)                                                  3475
19701 continue                                                             3476
19702 continue                                                             3476
      xv=1.0                                                               3477
19691 continue                                                             3478
19671 continue                                                             3478
      ys0=0.0                                                              3479
19710 do 19711 j=1,nr                                                      3480
      ym(j)=dot_product(w,y(:,j))                                          3480
      y(:,j)=v*(y(:,j)-ym(j))                                              3481
      z=dot_product(y(:,j),y(:,j))                                         3482
      if(jsd .le. 0)goto 19731                                             3482
      ys(j)=sqrt(z)                                                        3482
      y(:,j)=y(:,j)/ys(j)                                                  3482
      goto 19741                                                           3483
19731 continue                                                             3483
      ys0=ys0+z                                                            3483
19741 continue                                                             3484
19721 continue                                                             3484
19711 continue                                                             3485
19712 continue                                                             3485
      if(jsd .ne. 0)goto 19761                                             3485
      ys=1.0                                                               3485
      goto 19771                                                           3485
19761 continue                                                             3485
      ys0=nr                                                               3485
19771 continue                                                             3486
19751 continue                                                             3486
10700 continue                                                             3486
      deallocate(v)                                                        3487
      return                                                               3488
      end                                                                  3489
      subroutine multelnet2(beta,ni,nr,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,   3491 
     *ulam,thri,  maxit,xv,ys0,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   3492
      double precision vp(ni),y(no,nr),x(no,ni),ulam(nlam),ao(nx,nr,nlam   3493 
     *)
      double precision rsqo(nlam),almo(nlam),xv(ni),cl(2,nr,ni)            3494
      integer ju(ni),ia(nx),kin(nlam)                                      3495
      double precision, dimension (:), allocatable :: g,gk,del,gj               
      integer, dimension (:), allocatable :: mm,ix,isc                          
      double precision, dimension (:,:), allocatable :: a                       
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               3503
      allocate(gj(1:nr),stat=jerr)                                         3504
      if(jerr.ne.0) return                                                 3505
      allocate(gk(1:nr),stat=jerr)                                         3506
      if(jerr.ne.0) return                                                 3507
      allocate(del(1:nr),stat=jerr)                                        3508
      if(jerr.ne.0) return                                                 3509
      allocate(mm(1:ni),stat=jerr)                                         3510
      if(jerr.ne.0) return                                                 3511
      allocate(g(1:ni),stat=jerr)                                          3512
      if(jerr.ne.0) return                                                 3513
      allocate(ix(1:ni),stat=jerr)                                         3514
      if(jerr.ne.0) return                                                 3515
      allocate(isc(1:nr),stat=jerr)                                        3516
      if(jerr.ne.0) return                                                 3517
      bta=beta                                                             3517
      omb=1.0-bta                                                          3517
      ix=0                                                                 3517
      thr=thri*ys0/nr                                                      3519
      alf=1.0                                                              3521
      if(flmin .ge. 1.0)goto 19791                                         3521
      eqs=max(eps,flmin)                                                   3521
      alf=eqs**(1.0/(nlam-1))                                              3521
19791 continue                                                             3522
      rsq=ys0                                                              3522
      a=0.0                                                                3522
      mm=0                                                                 3522
      nlp=0                                                                3522
      nin=nlp                                                              3522
      iz=0                                                                 3522
      mnl=min(mnlam,nlam)                                                  3522
      alm=0.0                                                              3523
19800 do 19801 j=1,ni                                                      3523
      if(ju(j).eq.0)goto 19801                                             3523
      g(j)=0.0                                                             3524
19810 do 19811 k=1,nr                                                      3524
      g(j)=g(j)+dot_product(y(:,k),x(:,j))**2                              3524
19811 continue                                                             3525
19812 continue                                                             3525
      g(j)=sqrt(g(j))                                                      3526
19801 continue                                                             3527
19802 continue                                                             3527
19820 do 19821 m=1,nlam                                                    3527
      alm0=alm                                                             3528
      if(flmin .lt. 1.0)goto 19841                                         3528
      alm=ulam(m)                                                          3528
      goto 19831                                                           3529
19841 if(m .le. 2)goto 19851                                               3529
      alm=alm*alf                                                          3529
      goto 19831                                                           3530
19851 if(m .ne. 1)goto 19861                                               3530
      alm=big                                                              3530
      goto 19871                                                           3531
19861 continue                                                             3531
      alm0=0.0                                                             3532
19880 do 19881 j=1,ni                                                      3532
      if(ju(j).eq.0)goto 19881                                             3533
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           3534
19881 continue                                                             3535
19882 continue                                                             3535
      alm0=alm0/max(bta,1.0d-3)                                            3535
      alm=alf*alm0                                                         3536
19871 continue                                                             3537
19831 continue                                                             3537
      dem=alm*omb                                                          3537
      ab=alm*bta                                                           3537
      rsq0=rsq                                                             3537
      jz=1                                                                 3538
      tlam=bta*(2.0*alm-alm0)                                              3539
19890 do 19891 k=1,ni                                                      3539
      if(ix(k).eq.1)goto 19891                                             3539
      if(ju(k).eq.0)goto 19891                                             3540
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       3541
19891 continue                                                             3542
19892 continue                                                             3542
19900 continue                                                             3542
19901 continue                                                             3542
      if(iz*jz.ne.0) go to 10360                                           3543
10880 continue                                                             3543
      nlp=nlp+1                                                            3543
      dlx=0.0                                                              3544
19910 do 19911 k=1,ni                                                      3544
      if(ix(k).eq.0)goto 19911                                             3544
      gkn=0.0                                                              3545
19920 do 19921 j=1,nr                                                      3545
      gj(j)=dot_product(y(:,j),x(:,k))                                     3546
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3546
      gkn=gkn+gk(j)**2                                                     3548
19921 continue                                                             3548
19922 continue                                                             3548
      gkn=sqrt(gkn)                                                        3548
      u=1.0-ab*vp(k)/gkn                                                   3548
      del=a(:,k)                                                           3549
      if(u .gt. 0.0)goto 19941                                             3549
      a(:,k)=0.0                                                           3549
      goto 19951                                                           3550
19941 continue                                                             3550
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3551
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3553 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3554
19951 continue                                                             3555
19931 continue                                                             3555
      del=a(:,k)-del                                                       3555
      if(maxval(abs(del)).le.0.0)goto 19911                                3556
19960 do 19961 j=1,nr                                                      3556
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3557
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3557
      dlx=max(dlx,xv(k)*del(j)**2)                                         3558
19961 continue                                                             3559
19962 continue                                                             3559
      if(mm(k) .ne. 0)goto 19981                                           3559
      nin=nin+1                                                            3559
      if(nin.gt.nx)goto 19912                                              3560
      mm(k)=nin                                                            3560
      ia(nin)=k                                                            3561
19981 continue                                                             3562
19911 continue                                                             3563
19912 continue                                                             3563
      if(nin.gt.nx)goto 19902                                              3564
      if(dlx .ge. thr)goto 20001                                           3564
      ixx=0                                                                3565
20010 do 20011 k=1,ni                                                      3565
      if(ix(k).eq.1)goto 20011                                             3565
      if(ju(k).eq.0)goto 20011                                             3565
      g(k)=0.0                                                             3566
20020 do 20021 j=1,nr                                                      3566
      g(k)=g(k)+dot_product(y(:,j),x(:,k))**2                              3566
20021 continue                                                             3567
20022 continue                                                             3567
      g(k)=sqrt(g(k))                                                      3568
      if(g(k) .le. ab*vp(k))goto 20041                                     3568
      ix(k)=1                                                              3568
      ixx=1                                                                3568
20041 continue                                                             3569
20011 continue                                                             3570
20012 continue                                                             3570
      if(ixx.eq.1) go to 10880                                             3571
      goto 19902                                                           3572
20001 continue                                                             3573
      if(nlp .le. maxit)goto 20061                                         3573
      jerr=-m                                                              3573
      return                                                               3573
20061 continue                                                             3574
10360 continue                                                             3574
      iz=1                                                                 3575
20070 continue                                                             3575
20071 continue                                                             3575
      nlp=nlp+1                                                            3575
      dlx=0.0                                                              3576
20080 do 20081 l=1,nin                                                     3576
      k=ia(l)                                                              3576
      gkn=0.0                                                              3577
20090 do 20091 j=1,nr                                                      3577
      gj(j)=dot_product(y(:,j),x(:,k))                                     3578
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3578
      gkn=gkn+gk(j)**2                                                     3580
20091 continue                                                             3580
20092 continue                                                             3580
      gkn=sqrt(gkn)                                                        3580
      u=1.0-ab*vp(k)/gkn                                                   3580
      del=a(:,k)                                                           3581
      if(u .gt. 0.0)goto 20111                                             3581
      a(:,k)=0.0                                                           3581
      goto 20121                                                           3582
20111 continue                                                             3582
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3583
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3585 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3586
20121 continue                                                             3587
20101 continue                                                             3587
      del=a(:,k)-del                                                       3587
      if(maxval(abs(del)).le.0.0)goto 20081                                3588
20130 do 20131 j=1,nr                                                      3588
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3589
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          3589
      dlx=max(dlx,xv(k)*del(j)**2)                                         3590
20131 continue                                                             3591
20132 continue                                                             3591
20081 continue                                                             3592
20082 continue                                                             3592
      if(dlx.lt.thr)goto 20072                                             3592
      if(nlp .le. maxit)goto 20151                                         3592
      jerr=-m                                                              3592
      return                                                               3592
20151 continue                                                             3593
      goto 20071                                                           3594
20072 continue                                                             3594
      jz=0                                                                 3595
      goto 19901                                                           3596
19902 continue                                                             3596
      if(nin .le. nx)goto 20171                                            3596
      jerr=-10000-m                                                        3596
      goto 19822                                                           3596
20171 continue                                                             3597
      if(nin .le. 0)goto 20191                                             3597
20200 do 20201 j=1,nr                                                      3597
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3597
20201 continue                                                             3597
20202 continue                                                             3597
20191 continue                                                             3598
      kin(m)=nin                                                           3599
      rsqo(m)=1.0-rsq/ys0                                                  3599
      almo(m)=alm                                                          3599
      lmu=m                                                                3600
      if(m.lt.mnl)goto 19821                                               3600
      if(flmin.ge.1.0)goto 19821                                           3601
      me=0                                                                 3601
20210 do 20211 j=1,nin                                                     3601
      if(ao(j,1,m).ne.0.0) me=me+1                                         3601
20211 continue                                                             3601
20212 continue                                                             3601
      if(me.gt.ne)goto 19822                                               3602
      if(rsq0-rsq.lt.sml*rsq)goto 19822                                    3602
      if(rsqo(m).gt.rsqmax)goto 19822                                      3603
19821 continue                                                             3604
19822 continue                                                             3604
      deallocate(a,mm,g,ix,del,gj,gk)                                      3605
      return                                                               3606
      end                                                                  3607
      subroutine chkbnds(nr,gk,gkn,xv,cl,al1,al2,a,isc,jerr)               3608
      implicit double precision(a-h,o-z)                                   3609
      double precision gk(nr),cl(2,nr),a(nr)                               3609
      integer isc(nr)                                                      3610
      kerr=0                                                               3610
      al1p=1.0+al1/xv                                                      3610
      al2p=al2/xv                                                          3610
      isc=0                                                                3611
      gsq=gkn**2                                                           3611
      asq=dot_product(a,a)                                                 3611
      usq=0.0                                                              3613
      u=0.0                                                                3613
      kn=-1                                                                3615
20220 continue                                                             3615
20221 continue                                                             3615
      vmx=0.0                                                              3616
20230 do 20231 k=1,nr                                                      3616
      v=max(a(k)-cl(2,k),cl(1,k)-a(k))                                     3617
      if(v .le. vmx)goto 20251                                             3617
      vmx=v                                                                3617
      kn=k                                                                 3617
20251 continue                                                             3618
20231 continue                                                             3619
20232 continue                                                             3619
      if(vmx.le.0.0)goto 20222                                             3619
      if(isc(kn).ne.0)goto 20222                                           3620
      gsq=gsq-gk(kn)**2                                                    3620
      g=sqrt(gsq)/xv                                                       3621
      if(a(kn).lt.cl(1,kn)) u=cl(1,kn)                                     3621
      if(a(kn).gt.cl(2,kn)) u=cl(2,kn)                                     3622
      usq=usq+u**2                                                         3623
      if(usq .ne. 0.0)goto 20271                                           3623
      b=max(0d0,(g-al2p)/al1p)                                             3623
      goto 20281                                                           3624
20271 continue                                                             3624
      b0=sqrt(asq-a(kn)**2)                                                3625
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                     3625
      if(kerr.ne.0)goto 20222                                              3626
20281 continue                                                             3627
20261 continue                                                             3627
      asq=usq+b**2                                                         3627
      if(asq .gt. 0.0)goto 20301                                           3627
      a=0.0                                                                3627
      goto 20222                                                           3627
20301 continue                                                             3628
      a(kn)=u                                                              3628
      isc(kn)=1                                                            3628
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                     3629
20310 do 20311 j=1,nr                                                      3629
      if(isc(j).eq.0) a(j)=f*gk(j)                                         3629
20311 continue                                                             3630
20312 continue                                                             3630
      goto 20221                                                           3631
20222 continue                                                             3631
      if(kerr.ne.0) jerr=kerr                                              3632
      return                                                               3633
      end                                                                  3634
      subroutine chkbnds1(nr,gk,gkn,xv,cl1,cl2,al1,al2,a,isc,jerr)         3635
      implicit double precision(a-h,o-z)                                   3636
      double precision gk(nr),a(nr)                                        3636
      integer isc(nr)                                                      3637
      kerr=0                                                               3637
      al1p=1.0+al1/xv                                                      3637
      al2p=al2/xv                                                          3637
      isc=0                                                                3638
      gsq=gkn**2                                                           3638
      asq=dot_product(a,a)                                                 3638
      usq=0.0                                                              3640
      u=0.0                                                                3640
      kn=-1                                                                3642
20320 continue                                                             3642
20321 continue                                                             3642
      vmx=0.0                                                              3643
20330 do 20331 k=1,nr                                                      3643
      v=max(a(k)-cl2,cl1-a(k))                                             3644
      if(v .le. vmx)goto 20351                                             3644
      vmx=v                                                                3644
      kn=k                                                                 3644
20351 continue                                                             3645
20331 continue                                                             3646
20332 continue                                                             3646
      if(vmx.le.0.0)goto 20322                                             3646
      if(isc(kn).ne.0)goto 20322                                           3647
      gsq=gsq-gk(kn)**2                                                    3647
      g=sqrt(gsq)/xv                                                       3648
      if(a(kn).lt.cl1) u=cl1                                               3648
      if(a(kn).gt.cl2) u=cl2                                               3649
      usq=usq+u**2                                                         3650
      if(usq .ne. 0.0)goto 20371                                           3650
      b=max(0d0,(g-al2p)/al1p)                                             3650
      goto 20381                                                           3651
20371 continue                                                             3651
      b0=sqrt(asq-a(kn)**2)                                                3652
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                     3652
      if(kerr.ne.0)goto 20322                                              3653
20381 continue                                                             3654
20361 continue                                                             3654
      asq=usq+b**2                                                         3654
      if(asq .gt. 0.0)goto 20401                                           3654
      a=0.0                                                                3654
      goto 20322                                                           3654
20401 continue                                                             3655
      a(kn)=u                                                              3655
      isc(kn)=1                                                            3655
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                     3656
20410 do 20411 j=1,nr                                                      3656
      if(isc(j).eq.0) a(j)=f*gk(j)                                         3656
20411 continue                                                             3657
20412 continue                                                             3657
      goto 20321                                                           3658
20322 continue                                                             3658
      if(kerr.ne.0) jerr=kerr                                              3659
      return                                                               3660
      end                                                                  3661
      function bnorm(b0,al1p,al2p,g,usq,jerr)                              3662
      implicit double precision(a-h,o-z)                                   3663
      data thr,mxit /1.0d-10,100/                                          3664
      b=b0                                                                 3664
      zsq=b**2+usq                                                         3664
      if(zsq .gt. 0.0)goto 20431                                           3664
      bnorm=0.0                                                            3664
      return                                                               3664
20431 continue                                                             3665
      z=sqrt(zsq)                                                          3665
      f=b*(al1p+al2p/z)-g                                                  3665
      jerr=0                                                               3666
20440 do 20441 it=1,mxit                                                   3666
      b=b-f/(al1p+al2p*usq/(z*zsq))                                        3667
      zsq=b**2+usq                                                         3667
      if(zsq .gt. 0.0)goto 20461                                           3667
      bnorm=0.0                                                            3667
      return                                                               3667
20461 continue                                                             3668
      z=sqrt(zsq)                                                          3668
      f=b*(al1p+al2p/z)-g                                                  3669
      if(abs(f).le.thr)goto 20442                                          3669
      if(b .gt. 0.0)goto 20481                                             3669
      b=0.0                                                                3669
      goto 20442                                                           3669
20481 continue                                                             3670
20441 continue                                                             3671
20442 continue                                                             3671
      bnorm=b                                                              3671
      if(it.ge.mxit) jerr=90000                                            3672
      return                                                               3674
      entry chg_bnorm(arg,irg)                                             3674
      bnorm = 0.0                                                          3674
      thr=arg                                                              3674
      mxit=irg                                                             3674
      return                                                               3675
      entry get_bnorm(arg,irg)                                             3675
      bnorm = 0.0                                                          3675
      arg=thr                                                              3675
      irg=mxit                                                             3675
      return                                                               3677
      end                                                                  3678
      subroutine multsolns(ni,nx,nr,lmu,a,ia,nin,b)                        3679
      implicit double precision(a-h,o-z)                                   3680
      double precision a(nx,nr,lmu),b(ni,nr,lmu)                           3680
      integer ia(nx),nin(lmu)                                              3681
20490 do 20491 lam=1,lmu                                                   3681
      call multuncomp(ni,nr,nx,a(1,1,lam),ia,nin(lam),b(1,1,lam))          3681
20491 continue                                                             3682
20492 continue                                                             3682
      return                                                               3683
      end                                                                  3684
      subroutine multuncomp(ni,nr,nx,ca,ia,nin,a)                          3685
      implicit double precision(a-h,o-z)                                   3686
      double precision ca(nx,nr),a(ni,nr)                                  3686
      integer ia(nx)                                                       3687
      a=0.0                                                                3688
      if(nin .le. 0)goto 20511                                             3688
20520 do 20521 j=1,nr                                                      3688
      a(ia(1:nin),j)=ca(1:nin,j)                                           3688
20521 continue                                                             3688
20522 continue                                                             3688
20511 continue                                                             3689
      return                                                               3690
      end                                                                  3691
      subroutine multmodval(nx,nr,a0,ca,ia,nin,n,x,f)                      3692
      implicit double precision(a-h,o-z)                                   3693
      double precision a0(nr),ca(nx,nr),x(n,*),f(nr,n)                     3693
      integer ia(nx)                                                       3694
20530 do 20531 i=1,n                                                       3694
      f(:,i)=a0                                                            3694
20531 continue                                                             3694
20532 continue                                                             3694
      if(nin.le.0) return                                                  3695
20540 do 20541 i=1,n                                                       3695
20550 do 20551 j=1,nr                                                      3695
      f(j,i)=f(j,i)+dot_product(ca(1:nin,j),x(i,ia(1:nin)))                3695
20551 continue                                                             3695
20552 continue                                                             3695
20541 continue                                                             3696
20542 continue                                                             3696
      return                                                               3697
      end                                                                  3698
      subroutine multspelnet  (parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,   3701 
     *nlam,flmin,ulam,thr,isd,  jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,
     *nlp,jerr)
      implicit double precision(a-h,o-z)                                   3702
      double precision x(*),y(no,nr),w(no),vp(ni),ulam(nlam),cl(2,ni)      3703
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)      3704
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3705
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 20571                                    3708
      jerr=10000                                                           3708
      return                                                               3708
20571 continue                                                             3709
      allocate(vq(1:ni),stat=jerr)                                         3709
      if(jerr.ne.0) return                                                 3710
      vq=max(0d0,vp)                                                       3710
      vq=vq*ni/sum(vq)                                                     3711
      call multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,fl   3713 
     *min,  ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jer
     *r)
      deallocate(vq)                                                       3714
      return                                                               3715
      end                                                                  3716
      subroutine multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,n   3718 
     *lam,flmin,  ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,n
     *lp,jerr)
      implicit double precision(a-h,o-z)                                   3719
      double precision x(*),vp(ni),y(no,nr),w(no),ulam(nlam),cl(2,ni)      3720
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)      3721
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           3722
      double precision, dimension (:), allocatable :: xm,xs,xv,ym,ys            
      integer, dimension (:), allocatable :: ju                                 
      double precision, dimension (:,:,:), allocatable :: clt                   
      allocate(clt(1:2,1:nr,1:ni),stat=jerr)                                    
      if(jerr.ne.0) return                                                      
      allocate(xm(1:ni),stat=jerr)                                         3730
      if(jerr.ne.0) return                                                 3731
      allocate(xs(1:ni),stat=jerr)                                         3732
      if(jerr.ne.0) return                                                 3733
      allocate(ym(1:nr),stat=jerr)                                         3734
      if(jerr.ne.0) return                                                 3735
      allocate(ys(1:nr),stat=jerr)                                         3736
      if(jerr.ne.0) return                                                 3737
      allocate(ju(1:ni),stat=jerr)                                         3738
      if(jerr.ne.0) return                                                 3739
      allocate(xv(1:ni),stat=jerr)                                         3740
      if(jerr.ne.0) return                                                 3741
      call spchkvars(no,ni,x,ix,ju)                                        3742
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 3743
      if(maxval(ju) .gt. 0)goto 20591                                      3743
      jerr=7777                                                            3743
      return                                                               3743
20591 continue                                                             3744
      call multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,  xm,xs,   3746 
     *ym,ys,xv,ys0,jerr)
      if(jerr.ne.0) return                                                 3747
20600 do 20601 j=1,ni                                                      3747
20610 do 20611 k=1,nr                                                      3747
20620 do 20621 i=1,2                                                       3747
      clt(i,k,j)=cl(i,j)                                                   3747
20621 continue                                                             3747
20622 continue                                                             3747
20611 continue                                                             3747
20612 continue                                                             3747
20601 continue                                                             3748
20602 continue                                                             3748
      if(isd .le. 0)goto 20641                                             3748
20650 do 20651 j=1,ni                                                      3748
20660 do 20661 k=1,nr                                                      3748
20670 do 20671 i=1,2                                                       3748
      clt(i,k,j)=clt(i,k,j)*xs(j)                                          3748
20671 continue                                                             3748
20672 continue                                                             3748
20661 continue                                                             3748
20662 continue                                                             3748
20651 continue                                                             3748
20652 continue                                                             3748
20641 continue                                                             3749
      if(jsd .le. 0)goto 20691                                             3749
20700 do 20701 j=1,ni                                                      3749
20710 do 20711 k=1,nr                                                      3749
20720 do 20721 i=1,2                                                       3749
      clt(i,k,j)=clt(i,k,j)/ys(k)                                          3749
20721 continue                                                             3749
20722 continue                                                             3749
20711 continue                                                             3749
20712 continue                                                             3749
20701 continue                                                             3749
20702 continue                                                             3749
20691 continue                                                             3750
      call multspelnet2(parm,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,clt,nlam,f   3752 
     *lmin,  ulam,thr,maxit,xm,xs,xv,ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 3753
20730 do 20731 k=1,lmu                                                     3753
      nk=nin(k)                                                            3754
20740 do 20741 j=1,nr                                                      3755
20750 do 20751 l=1,nk                                                      3755
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  3755
20751 continue                                                             3756
20752 continue                                                             3756
      if(intr .ne. 0)goto 20771                                            3756
      a0(j,k)=0.0                                                          3756
      goto 20781                                                           3757
20771 continue                                                             3757
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 3757
20781 continue                                                             3758
20761 continue                                                             3758
20741 continue                                                             3759
20742 continue                                                             3759
20731 continue                                                             3760
20732 continue                                                             3760
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                    3761
      return                                                               3762
      end                                                                  3763
      subroutine multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,     3765 
     *xm,xs,ym,ys,xv,ys0,jerr)
      implicit double precision(a-h,o-z)                                   3766
      double precision x(*),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),y   3767 
     *s(nr)
      integer ix(*),jx(*),ju(ni)                                           3768
      w=w/sum(w)                                                           3769
      if(intr .ne. 0)goto 20801                                            3770
20810 do 20811 j=1,ni                                                      3770
      if(ju(j).eq.0)goto 20811                                             3770
      xm(j)=0.0                                                            3770
      jb=ix(j)                                                             3770
      je=ix(j+1)-1                                                         3771
      z=dot_product(w(jx(jb:je)),x(jb:je)**2)                              3772
      if(isd .le. 0)goto 20831                                             3772
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            3772
      vc=z-xbq                                                             3773
      xs(j)=sqrt(vc)                                                       3773
      xv(j)=1.0+xbq/vc                                                     3774
      goto 20841                                                           3775
20831 continue                                                             3775
      xs(j)=1.0                                                            3775
      xv(j)=z                                                              3775
20841 continue                                                             3776
20821 continue                                                             3776
20811 continue                                                             3777
20812 continue                                                             3777
      ys0=0.0                                                              3778
20850 do 20851 j=1,nr                                                      3778
      ym(j)=0.0                                                            3778
      z=dot_product(w,y(:,j)**2)                                           3779
      if(jsd .le. 0)goto 20871                                             3779
      u=z-dot_product(w,y(:,j))**2                                         3779
      ys0=ys0+z/u                                                          3780
      ys(j)=sqrt(u)                                                        3780
      y(:,j)=y(:,j)/ys(j)                                                  3781
      goto 20881                                                           3782
20871 continue                                                             3782
      ys(j)=1.0                                                            3782
      ys0=ys0+z                                                            3782
20881 continue                                                             3783
20861 continue                                                             3783
20851 continue                                                             3784
20852 continue                                                             3784
      return                                                               3785
20801 continue                                                             3786
20890 do 20891 j=1,ni                                                      3786
      if(ju(j).eq.0)goto 20891                                             3787
      jb=ix(j)                                                             3787
      je=ix(j+1)-1                                                         3787
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             3788
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 3789
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       3790
20891 continue                                                             3791
20892 continue                                                             3791
      if(isd .ne. 0)goto 20911                                             3791
      xs=1.0                                                               3791
      goto 20921                                                           3791
20911 continue                                                             3791
      xv=1.0                                                               3791
20921 continue                                                             3792
20901 continue                                                             3792
      ys0=0.0                                                              3793
20930 do 20931 j=1,nr                                                      3794
      ym(j)=dot_product(w,y(:,j))                                          3794
      y(:,j)=y(:,j)-ym(j)                                                  3795
      z=dot_product(w,y(:,j)**2)                                           3796
      if(jsd .le. 0)goto 20951                                             3796
      ys(j)=sqrt(z)                                                        3796
      y(:,j)=y(:,j)/ys(j)                                                  3796
      goto 20961                                                           3797
20951 continue                                                             3797
      ys0=ys0+z                                                            3797
20961 continue                                                             3798
20941 continue                                                             3798
20931 continue                                                             3799
20932 continue                                                             3799
      if(jsd .ne. 0)goto 20981                                             3799
      ys=1.0                                                               3799
      goto 20991                                                           3799
20981 continue                                                             3799
      ys0=nr                                                               3799
20991 continue                                                             3800
20971 continue                                                             3800
      return                                                               3801
      end                                                                  3802
      subroutine multspelnet2(beta,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,n   3804 
     *lam,flmin,  ulam,thri,maxit,xm,xs,xv,ys0,lmu,ao,ia,kin,rsqo,almo,n
     *lp,jerr)
      implicit double precision(a-h,o-z)                                   3805
      double precision y(no,nr),w(no),x(*),vp(ni),ulam(nlam),cl(2,nr,ni)   3806
      double precision ao(nx,nr,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni   3807 
     *),xv(ni)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          3808
      double precision, dimension (:), allocatable :: g,gj,gk,del,o             
      integer, dimension (:), allocatable :: mm,iy,isc                          
      double precision, dimension (:,:), allocatable :: a                       
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               3816
      allocate(mm(1:ni),stat=jerr)                                         3817
      if(jerr.ne.0) return                                                 3818
      allocate(g(1:ni),stat=jerr)                                          3819
      if(jerr.ne.0) return                                                 3820
      allocate(gj(1:nr),stat=jerr)                                         3821
      if(jerr.ne.0) return                                                 3822
      allocate(gk(1:nr),stat=jerr)                                         3823
      if(jerr.ne.0) return                                                 3824
      allocate(del(1:nr),stat=jerr)                                        3825
      if(jerr.ne.0) return                                                 3826
      allocate(o(1:nr),stat=jerr)                                          3827
      if(jerr.ne.0) return                                                 3828
      allocate(iy(1:ni),stat=jerr)                                         3829
      if(jerr.ne.0) return                                                 3830
      allocate(isc(1:nr),stat=jerr)                                        3831
      if(jerr.ne.0) return                                                 3832
      bta=beta                                                             3832
      omb=1.0-bta                                                          3832
      alm=0.0                                                              3832
      iy=0                                                                 3832
      thr=thri*ys0/nr                                                      3834
      alf=1.0                                                              3836
      if(flmin .ge. 1.0)goto 21011                                         3836
      eqs=max(eps,flmin)                                                   3836
      alf=eqs**(1.0/(nlam-1))                                              3836
21011 continue                                                             3837
      rsq=ys0                                                              3837
      a=0.0                                                                3837
      mm=0                                                                 3837
      o=0.0                                                                3837
      nlp=0                                                                3837
      nin=nlp                                                              3837
      iz=0                                                                 3837
      mnl=min(mnlam,nlam)                                                  3838
21020 do 21021 j=1,ni                                                      3838
      if(ju(j).eq.0)goto 21021                                             3838
      jb=ix(j)                                                             3838
      je=ix(j+1)-1                                                         3838
      g(j)=0.0                                                             3839
21030 do 21031 k=1,nr                                                      3840
      g(j)=g(j)+(dot_product(y(jx(jb:je),k),w(jx(jb:je))*x(jb:je))/xs(j)   3841 
     *)**2
21031 continue                                                             3842
21032 continue                                                             3842
      g(j)=sqrt(g(j))                                                      3843
21021 continue                                                             3844
21022 continue                                                             3844
21040 do 21041 m=1,nlam                                                    3844
      alm0=alm                                                             3845
      if(flmin .lt. 1.0)goto 21061                                         3845
      alm=ulam(m)                                                          3845
      goto 21051                                                           3846
21061 if(m .le. 2)goto 21071                                               3846
      alm=alm*alf                                                          3846
      goto 21051                                                           3847
21071 if(m .ne. 1)goto 21081                                               3847
      alm=big                                                              3847
      goto 21091                                                           3848
21081 continue                                                             3848
      alm0=0.0                                                             3849
21100 do 21101 j=1,ni                                                      3849
      if(ju(j).eq.0)goto 21101                                             3850
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           3851
21101 continue                                                             3852
21102 continue                                                             3852
      alm0=alm0/max(bta,1.0d-3)                                            3852
      alm=alf*alm0                                                         3853
21091 continue                                                             3854
21051 continue                                                             3854
      dem=alm*omb                                                          3854
      ab=alm*bta                                                           3854
      rsq0=rsq                                                             3854
      jz=1                                                                 3855
      tlam=bta*(2.0*alm-alm0)                                              3856
21110 do 21111 k=1,ni                                                      3856
      if(iy(k).eq.1)goto 21111                                             3856
      if(ju(k).eq.0)goto 21111                                             3857
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       3858
21111 continue                                                             3859
21112 continue                                                             3859
21120 continue                                                             3859
21121 continue                                                             3859
      if(iz*jz.ne.0) go to 10360                                           3860
10880 continue                                                             3860
      nlp=nlp+1                                                            3860
      dlx=0.0                                                              3861
21130 do 21131 k=1,ni                                                      3861
      if(iy(k).eq.0)goto 21131                                             3861
      jb=ix(k)                                                             3861
      je=ix(k+1)-1                                                         3861
      gkn=0.0                                                              3862
21140 do 21141 j=1,nr                                                      3863
      gj(j)=dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(k)   3864
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3864
      gkn=gkn+gk(j)**2                                                     3865
21141 continue                                                             3866
21142 continue                                                             3866
      gkn=sqrt(gkn)                                                        3866
      u=1.0-ab*vp(k)/gkn                                                   3866
      del=a(:,k)                                                           3867
      if(u .gt. 0.0)goto 21161                                             3867
      a(:,k)=0.0                                                           3867
      goto 21171                                                           3868
21161 continue                                                             3868
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3869
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3871 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3872
21171 continue                                                             3873
21151 continue                                                             3873
      del=a(:,k)-del                                                       3873
      if(maxval(abs(del)).le.0.0)goto 21131                                3874
      if(mm(k) .ne. 0)goto 21191                                           3874
      nin=nin+1                                                            3874
      if(nin.gt.nx)goto 21132                                              3875
      mm(k)=nin                                                            3875
      ia(nin)=k                                                            3876
21191 continue                                                             3877
21200 do 21201 j=1,nr                                                      3877
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3878
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3879
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3879
      dlx=max(xv(k)*del(j)**2,dlx)                                         3880
21201 continue                                                             3881
21202 continue                                                             3881
21131 continue                                                             3882
21132 continue                                                             3882
      if(nin.gt.nx)goto 21122                                              3883
      if(dlx .ge. thr)goto 21221                                           3883
      ixx=0                                                                3884
21230 do 21231 j=1,ni                                                      3884
      if(iy(j).eq.1)goto 21231                                             3884
      if(ju(j).eq.0)goto 21231                                             3885
      jb=ix(j)                                                             3885
      je=ix(j+1)-1                                                         3885
      g(j)=0.0                                                             3886
21240 do 21241 k=1,nr                                                      3886
      g(j)=g(j)+  (dot_product(y(jx(jb:je),k)+o(k),w(jx(jb:je))*x(jb:je)   3888 
     *)/xs(j))**2
21241 continue                                                             3889
21242 continue                                                             3889
      g(j)=sqrt(g(j))                                                      3890
      if(g(j) .le. ab*vp(j))goto 21261                                     3890
      iy(j)=1                                                              3890
      ixx=1                                                                3890
21261 continue                                                             3891
21231 continue                                                             3892
21232 continue                                                             3892
      if(ixx.eq.1) go to 10880                                             3893
      goto 21122                                                           3894
21221 continue                                                             3895
      if(nlp .le. maxit)goto 21281                                         3895
      jerr=-m                                                              3895
      return                                                               3895
21281 continue                                                             3896
10360 continue                                                             3896
      iz=1                                                                 3897
21290 continue                                                             3897
21291 continue                                                             3897
      nlp=nlp+1                                                            3897
      dlx=0.0                                                              3898
21300 do 21301 l=1,nin                                                     3898
      k=ia(l)                                                              3898
      jb=ix(k)                                                             3898
      je=ix(k+1)-1                                                         3898
      gkn=0.0                                                              3899
21310 do 21311 j=1,nr                                                      3899
      gj(j)=  dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(   3901 
     *k)
      gk(j)=gj(j)+a(j,k)*xv(k)                                             3901
      gkn=gkn+gk(j)**2                                                     3902
21311 continue                                                             3903
21312 continue                                                             3903
      gkn=sqrt(gkn)                                                        3903
      u=1.0-ab*vp(k)/gkn                                                   3903
      del=a(:,k)                                                           3904
      if(u .gt. 0.0)goto 21331                                             3904
      a(:,k)=0.0                                                           3904
      goto 21341                                                           3905
21331 continue                                                             3905
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      3906
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   3908 
     *,isc,jerr)
      if(jerr.ne.0) return                                                 3909
21341 continue                                                             3910
21321 continue                                                             3910
      del=a(:,k)-del                                                       3910
      if(maxval(abs(del)).le.0.0)goto 21301                                3911
21350 do 21351 j=1,nr                                                      3911
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              3912
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  3913
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         3913
      dlx=max(xv(k)*del(j)**2,dlx)                                         3914
21351 continue                                                             3915
21352 continue                                                             3915
21301 continue                                                             3916
21302 continue                                                             3916
      if(dlx.lt.thr)goto 21292                                             3916
      if(nlp .le. maxit)goto 21371                                         3916
      jerr=-m                                                              3916
      return                                                               3916
21371 continue                                                             3917
      goto 21291                                                           3918
21292 continue                                                             3918
      jz=0                                                                 3919
      goto 21121                                                           3920
21122 continue                                                             3920
      if(nin .le. nx)goto 21391                                            3920
      jerr=-10000-m                                                        3920
      goto 21042                                                           3920
21391 continue                                                             3921
      if(nin .le. 0)goto 21411                                             3921
21420 do 21421 j=1,nr                                                      3921
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         3921
21421 continue                                                             3921
21422 continue                                                             3921
21411 continue                                                             3922
      kin(m)=nin                                                           3923
      rsqo(m)=1.0-rsq/ys0                                                  3923
      almo(m)=alm                                                          3923
      lmu=m                                                                3924
      if(m.lt.mnl)goto 21041                                               3924
      if(flmin.ge.1.0)goto 21041                                           3925
      me=0                                                                 3925
21430 do 21431 j=1,nin                                                     3925
      if(ao(j,1,m).ne.0.0) me=me+1                                         3925
21431 continue                                                             3925
21432 continue                                                             3925
      if(me.gt.ne)goto 21042                                               3926
      if(rsq0-rsq.lt.sml*rsq)goto 21042                                    3926
      if(rsqo(m).gt.rsqmax)goto 21042                                      3927
21041 continue                                                             3928
21042 continue                                                             3928
      deallocate(a,mm,g,iy,gj,gk,del,o)                                    3929
      return                                                               3930
      end                                                                  3931
      subroutine multlognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,f   3933 
     *lmin,ulam,  shri,intr,maxit,xv,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                   3934
      double precision x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam   3935 
     *),cl(2,ni)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),xv(   3936 
     *ni)
      integer ju(ni),m(nx),kin(nlam)                                       3937
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
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               3950
      exmn=-exmx                                                           3951
      allocate(mm(1:ni),stat=jerr)                                         3952
      if(jerr.ne.0) return                                                 3953
      allocate(is(1:max(nc,ni)),stat=jerr)                                 3954
      if(jerr.ne.0) return                                                 3955
      allocate(sxp(1:no),stat=jerr)                                        3956
      if(jerr.ne.0) return                                                 3957
      allocate(sxpl(1:no),stat=jerr)                                       3958
      if(jerr.ne.0) return                                                 3959
      allocate(ga(1:ni),stat=jerr)                                         3960
      if(jerr.ne.0) return                                                 3961
      allocate(ixx(1:ni),stat=jerr)                                        3962
      if(jerr.ne.0) return                                                 3963
      allocate(gk(1:nc),stat=jerr)                                         3964
      if(jerr.ne.0) return                                                 3965
      allocate(del(1:nc),stat=jerr)                                        3966
      if(jerr.ne.0) return                                                 3967
      allocate(isc(1:nc),stat=jerr)                                        3968
      if(jerr.ne.0) return                                                 3969
      pmax=1.0-pmin                                                        3969
      emin=pmin/pmax                                                       3969
      emax=1.0/emin                                                        3970
      bta=parm                                                             3970
      omb=1.0-bta                                                          3970
      dev1=0.0                                                             3970
      dev0=0.0                                                             3971
21440 do 21441 ic=1,nc                                                     3971
      q0=dot_product(w,y(:,ic))                                            3972
      if(q0 .gt. pmin)goto 21461                                           3972
      jerr =8000+ic                                                        3972
      return                                                               3972
21461 continue                                                             3973
      if(q0 .lt. pmax)goto 21481                                           3973
      jerr =9000+ic                                                        3973
      return                                                               3973
21481 continue                                                             3974
      if(intr .ne. 0)goto 21501                                            3974
      q0=1.0/nc                                                            3974
      b(0,ic)=0.0                                                          3974
      goto 21511                                                           3975
21501 continue                                                             3975
      b(0,ic)=log(q0)                                                      3975
      dev1=dev1-q0*b(0,ic)                                                 3975
21511 continue                                                             3976
21491 continue                                                             3976
      b(1:ni,ic)=0.0                                                       3977
21441 continue                                                             3978
21442 continue                                                             3978
      if(intr.eq.0) dev1=log(float(nc))                                    3978
      ixx=0                                                                3978
      al=0.0                                                               3979
      if(nonzero(no*nc,g) .ne. 0)goto 21531                                3980
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         3980
      sxp=0.0                                                              3981
21540 do 21541 ic=1,nc                                                     3981
      q(:,ic)=exp(b(0,ic))                                                 3981
      sxp=sxp+q(:,ic)                                                      3981
21541 continue                                                             3982
21542 continue                                                             3982
      goto 21551                                                           3983
21531 continue                                                             3983
21560 do 21561 i=1,no                                                      3983
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         3983
21561 continue                                                             3983
21562 continue                                                             3983
      sxp=0.0                                                              3984
      if(intr .ne. 0)goto 21581                                            3984
      b(0,:)=0.0                                                           3984
      goto 21591                                                           3985
21581 continue                                                             3985
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 3985
      if(jerr.ne.0) return                                                 3985
21591 continue                                                             3986
21571 continue                                                             3986
      dev1=0.0                                                             3987
21600 do 21601 ic=1,nc                                                     3987
      q(:,ic)=b(0,ic)+g(:,ic)                                              3988
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             3989
      q(:,ic)=exp(q(:,ic))                                                 3989
      sxp=sxp+q(:,ic)                                                      3990
21601 continue                                                             3991
21602 continue                                                             3991
      sxpl=w*log(sxp)                                                      3991
21610 do 21611 ic=1,nc                                                     3991
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  3991
21611 continue                                                             3992
21612 continue                                                             3992
21551 continue                                                             3993
21521 continue                                                             3993
21620 do 21621 ic=1,nc                                                     3993
21630 do 21631 i=1,no                                                      3993
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               3993
21631 continue                                                             3993
21632 continue                                                             3993
21621 continue                                                             3994
21622 continue                                                             3994
      dev0=dev0+dev1                                                       3996
      alf=1.0                                                              3998
      if(flmin .ge. 1.0)goto 21651                                         3998
      eqs=max(eps,flmin)                                                   3998
      alf=eqs**(1.0/(nlam-1))                                              3998
21651 continue                                                             3999
      m=0                                                                  3999
      mm=0                                                                 3999
      nin=0                                                                3999
      nlp=0                                                                3999
      mnl=min(mnlam,nlam)                                                  3999
      bs=0.0                                                               3999
      shr=shri*dev0                                                        4000
      ga=0.0                                                               4001
21660 do 21661 ic=1,nc                                                     4001
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      4002
21670 do 21671 j=1,ni                                                      4002
      if(ju(j).ne.0) ga(j)=ga(j)+dot_product(r(:,ic),x(:,j))**2            4002
21671 continue                                                             4003
21672 continue                                                             4003
21661 continue                                                             4004
21662 continue                                                             4004
      ga=sqrt(ga)                                                          4005
21680 do 21681 ilm=1,nlam                                                  4005
      al0=al                                                               4006
      if(flmin .lt. 1.0)goto 21701                                         4006
      al=ulam(ilm)                                                         4006
      goto 21691                                                           4007
21701 if(ilm .le. 2)goto 21711                                             4007
      al=al*alf                                                            4007
      goto 21691                                                           4008
21711 if(ilm .ne. 1)goto 21721                                             4008
      al=big                                                               4008
      goto 21731                                                           4009
21721 continue                                                             4009
      al0=0.0                                                              4010
21740 do 21741 j=1,ni                                                      4010
      if(ju(j).eq.0)goto 21741                                             4010
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            4010
21741 continue                                                             4011
21742 continue                                                             4011
      al0=al0/max(bta,1.0d-3)                                              4011
      al=alf*al0                                                           4012
21731 continue                                                             4013
21691 continue                                                             4013
      al2=al*omb                                                           4013
      al1=al*bta                                                           4013
      tlam=bta*(2.0*al-al0)                                                4014
21750 do 21751 k=1,ni                                                      4014
      if(ixx(k).eq.1)goto 21751                                            4014
      if(ju(k).eq.0)goto 21751                                             4015
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     4016
21751 continue                                                             4017
21752 continue                                                             4017
10880 continue                                                             4018
21760 continue                                                             4018
21761 continue                                                             4018
      ix=0                                                                 4018
      jx=ix                                                                4018
      kx=jx                                                                4018
      t=0.0                                                                4019
21770 do 21771 ic=1,nc                                                     4019
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       4019
21771 continue                                                             4020
21772 continue                                                             4020
      if(t .ge. eps)goto 21791                                             4020
      kx=1                                                                 4020
      goto 21762                                                           4020
21791 continue                                                             4020
      t=2.0*t                                                              4020
      alt=al1/t                                                            4020
      al2t=al2/t                                                           4021
21800 do 21801 ic=1,nc                                                     4022
      bs(0,ic)=b(0,ic)                                                     4022
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          4023
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    4024
      d=0.0                                                                4024
      if(intr.ne.0) d=sum(r(:,ic))                                         4025
      if(d .eq. 0.0)goto 21821                                             4026
      b(0,ic)=b(0,ic)+d                                                    4026
      r(:,ic)=r(:,ic)-d*w                                                  4026
      dlx=max(dlx,d**2)                                                    4027
21821 continue                                                             4028
21801 continue                                                             4029
21802 continue                                                             4029
21830 continue                                                             4029
21831 continue                                                             4029
      nlp=nlp+nc                                                           4029
      dlx=0.0                                                              4030
21840 do 21841 k=1,ni                                                      4030
      if(ixx(k).eq.0)goto 21841                                            4030
      gkn=0.0                                                              4031
21850 do 21851 ic=1,nc                                                     4031
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     4032
      gkn=gkn+gk(ic)**2                                                    4033
21851 continue                                                             4034
21852 continue                                                             4034
      gkn=sqrt(gkn)                                                        4034
      u=1.0-alt*vp(k)/gkn                                                  4034
      del=b(k,:)                                                           4035
      if(u .gt. 0.0)goto 21871                                             4035
      b(k,:)=0.0                                                           4035
      goto 21881                                                           4036
21871 continue                                                             4036
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     4037
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(   4039 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 4040
21881 continue                                                             4041
21861 continue                                                             4041
      del=b(k,:)-del                                                       4041
      if(maxval(abs(del)).le.0.0)goto 21841                                4042
21890 do 21891 ic=1,nc                                                     4042
      dlx=max(dlx,xv(k)*del(ic)**2)                                        4043
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     4044
21891 continue                                                             4045
21892 continue                                                             4045
      if(mm(k) .ne. 0)goto 21911                                           4045
      nin=nin+1                                                            4046
      if(nin .le. nx)goto 21931                                            4046
      jx=1                                                                 4046
      goto 21842                                                           4046
21931 continue                                                             4047
      mm(k)=nin                                                            4047
      m(nin)=k                                                             4048
21911 continue                                                             4049
21841 continue                                                             4050
21842 continue                                                             4050
      if(jx.gt.0)goto 21832                                                4050
      if(dlx.lt.shr)goto 21832                                             4051
      if(nlp .le. maxit)goto 21951                                         4051
      jerr=-ilm                                                            4051
      return                                                               4051
21951 continue                                                             4052
21960 continue                                                             4052
21961 continue                                                             4052
      nlp=nlp+nc                                                           4052
      dlx=0.0                                                              4053
21970 do 21971 l=1,nin                                                     4053
      k=m(l)                                                               4053
      gkn=0.0                                                              4054
21980 do 21981 ic=1,nc                                                     4054
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     4055
      gkn=gkn+gk(ic)**2                                                    4056
21981 continue                                                             4057
21982 continue                                                             4057
      gkn=sqrt(gkn)                                                        4057
      u=1.0-alt*vp(k)/gkn                                                  4057
      del=b(k,:)                                                           4058
      if(u .gt. 0.0)goto 22001                                             4058
      b(k,:)=0.0                                                           4058
      goto 22011                                                           4059
22001 continue                                                             4059
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     4060
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(   4062 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 4063
22011 continue                                                             4064
21991 continue                                                             4064
      del=b(k,:)-del                                                       4064
      if(maxval(abs(del)).le.0.0)goto 21971                                4065
22020 do 22021 ic=1,nc                                                     4065
      dlx=max(dlx,xv(k)*del(ic)**2)                                        4066
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     4067
22021 continue                                                             4068
22022 continue                                                             4068
21971 continue                                                             4069
21972 continue                                                             4069
      if(dlx.lt.shr)goto 21962                                             4069
      if(nlp .le. maxit)goto 22041                                         4069
      jerr=-ilm                                                            4069
      return                                                               4069
22041 continue                                                             4071
      goto 21961                                                           4072
21962 continue                                                             4072
      goto 21831                                                           4073
21832 continue                                                             4073
      if(jx.gt.0)goto 21762                                                4074
22050 do 22051 ic=1,nc                                                     4075
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                                4076
      if(ix .ne. 0)goto 22071                                              4077
22080 do 22081 j=1,nin                                                     4077
      k=m(j)                                                               4078
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 22101                   4078
      ix=1                                                                 4078
      goto 22082                                                           4078
22101 continue                                                             4080
22081 continue                                                             4081
22082 continue                                                             4081
22071 continue                                                             4082
22110 do 22111 i=1,no                                                      4082
      fi=b(0,ic)+g(i,ic)                                                   4084
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         4085
      fi=min(max(exmn,fi),exmx)                                            4085
      sxp(i)=sxp(i)-q(i,ic)                                                4086
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    4087
      sxp(i)=sxp(i)+q(i,ic)                                                4088
22111 continue                                                             4089
22112 continue                                                             4089
22051 continue                                                             4090
22052 continue                                                             4090
      s=-sum(b(0,:))/nc                                                    4090
      b(0,:)=b(0,:)+s                                                      4091
      if(jx.gt.0)goto 21762                                                4092
      if(ix .ne. 0)goto 22131                                              4093
22140 do 22141 k=1,ni                                                      4093
      if(ixx(k).eq.1)goto 22141                                            4093
      if(ju(k).eq.0)goto 22141                                             4093
      ga(k)=0.0                                                            4093
22141 continue                                                             4094
22142 continue                                                             4094
22150 do 22151 ic=1,nc                                                     4094
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      4095
22160 do 22161 k=1,ni                                                      4095
      if(ixx(k).eq.1)goto 22161                                            4095
      if(ju(k).eq.0)goto 22161                                             4096
      ga(k)=ga(k)+dot_product(r(:,ic),x(:,k))**2                           4097
22161 continue                                                             4098
22162 continue                                                             4098
22151 continue                                                             4099
22152 continue                                                             4099
      ga=sqrt(ga)                                                          4100
22170 do 22171 k=1,ni                                                      4100
      if(ixx(k).eq.1)goto 22171                                            4100
      if(ju(k).eq.0)goto 22171                                             4101
      if(ga(k) .le. al1*vp(k))goto 22191                                   4101
      ixx(k)=1                                                             4101
      ix=1                                                                 4101
22191 continue                                                             4102
22171 continue                                                             4103
22172 continue                                                             4103
      if(ix.eq.1) go to 10880                                              4104
      goto 21762                                                           4105
22131 continue                                                             4106
      goto 21761                                                           4107
21762 continue                                                             4107
      if(kx .le. 0)goto 22211                                              4107
      jerr=-20000-ilm                                                      4107
      goto 21682                                                           4107
22211 continue                                                             4108
      if(jx .le. 0)goto 22231                                              4108
      jerr=-10000-ilm                                                      4108
      goto 21682                                                           4108
22231 continue                                                             4108
      devi=0.0                                                             4109
22240 do 22241 ic=1,nc                                                     4110
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          4110
      a0(ic,ilm)=b(0,ic)                                                   4111
22250 do 22251 i=1,no                                                      4111
      if(y(i,ic).le.0.0)goto 22251                                         4112
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           4113
22251 continue                                                             4114
22252 continue                                                             4114
22241 continue                                                             4115
22242 continue                                                             4115
      kin(ilm)=nin                                                         4115
      alm(ilm)=al                                                          4115
      lmu=ilm                                                              4116
      dev(ilm)=(dev1-devi)/dev0                                            4117
      if(ilm.lt.mnl)goto 21681                                             4117
      if(flmin.ge.1.0)goto 21681                                           4118
      me=0                                                                 4118
22260 do 22261 j=1,nin                                                     4118
      if(a(j,1,ilm).ne.0.0) me=me+1                                        4118
22261 continue                                                             4118
22262 continue                                                             4118
      if(me.gt.ne)goto 21682                                               4119
      if(dev(ilm).gt.devmax)goto 21682                                     4119
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 21682                             4120
21681 continue                                                             4121
21682 continue                                                             4121
      g=log(q)                                                             4121
22270 do 22271 i=1,no                                                      4121
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         4121
22271 continue                                                             4122
22272 continue                                                             4122
      deallocate(sxp,b,bs,r,q,mm,is,ga,ixx,gk,del,sxpl)                    4123
      return                                                               4124
      end                                                                  4125
      subroutine multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,   4127 
     *nx,nlam,  flmin,ulam,shri,intr,maxit,xv,xb,xs,lmu,a0,a,m,kin,dev0,
     *dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   4128
      double precision x(*),y(no,nc),g(no,nc),w(no),vp(ni)                 4129
      double precision ulam(nlam),xb(ni),xs(ni),xv(ni)                     4130
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(   4131 
     *2,ni)
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           4132
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
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               4146
      exmn=-exmx                                                           4147
      allocate(mm(1:ni),stat=jerr)                                         4148
      if(jerr.ne.0) return                                                 4149
      allocate(ga(1:ni),stat=jerr)                                         4150
      if(jerr.ne.0) return                                                 4151
      allocate(gk(1:nc),stat=jerr)                                         4152
      if(jerr.ne.0) return                                                 4153
      allocate(del(1:nc),stat=jerr)                                        4154
      if(jerr.ne.0) return                                                 4155
      allocate(iy(1:ni),stat=jerr)                                         4156
      if(jerr.ne.0) return                                                 4157
      allocate(is(1:max(nc,ni)),stat=jerr)                                 4158
      if(jerr.ne.0) return                                                 4159
      allocate(sxp(1:no),stat=jerr)                                        4160
      if(jerr.ne.0) return                                                 4161
      allocate(sxpl(1:no),stat=jerr)                                       4162
      if(jerr.ne.0) return                                                 4163
      allocate(svr(1:nc),stat=jerr)                                        4164
      if(jerr.ne.0) return                                                 4165
      allocate(sc(1:no),stat=jerr)                                         4166
      if(jerr.ne.0) return                                                 4167
      allocate(isc(1:nc),stat=jerr)                                        4168
      if(jerr.ne.0) return                                                 4169
      pmax=1.0-pmin                                                        4169
      emin=pmin/pmax                                                       4169
      emax=1.0/emin                                                        4170
      bta=parm                                                             4170
      omb=1.0-bta                                                          4170
      dev1=0.0                                                             4170
      dev0=0.0                                                             4171
22280 do 22281 ic=1,nc                                                     4171
      q0=dot_product(w,y(:,ic))                                            4172
      if(q0 .gt. pmin)goto 22301                                           4172
      jerr =8000+ic                                                        4172
      return                                                               4172
22301 continue                                                             4173
      if(q0 .lt. pmax)goto 22321                                           4173
      jerr =9000+ic                                                        4173
      return                                                               4173
22321 continue                                                             4174
      b(1:ni,ic)=0.0                                                       4175
      if(intr .ne. 0)goto 22341                                            4175
      q0=1.0/nc                                                            4175
      b(0,ic)=0.0                                                          4175
      goto 22351                                                           4176
22341 continue                                                             4176
      b(0,ic)=log(q0)                                                      4176
      dev1=dev1-q0*b(0,ic)                                                 4176
22351 continue                                                             4177
22331 continue                                                             4177
22281 continue                                                             4178
22282 continue                                                             4178
      if(intr.eq.0) dev1=log(float(nc))                                    4178
      iy=0                                                                 4178
      al=0.0                                                               4179
      if(nonzero(no*nc,g) .ne. 0)goto 22371                                4180
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         4180
      sxp=0.0                                                              4181
22380 do 22381 ic=1,nc                                                     4181
      q(:,ic)=exp(b(0,ic))                                                 4181
      sxp=sxp+q(:,ic)                                                      4181
22381 continue                                                             4182
22382 continue                                                             4182
      goto 22391                                                           4183
22371 continue                                                             4183
22400 do 22401 i=1,no                                                      4183
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         4183
22401 continue                                                             4183
22402 continue                                                             4183
      sxp=0.0                                                              4184
      if(intr .ne. 0)goto 22421                                            4184
      b(0,:)=0.0                                                           4184
      goto 22431                                                           4185
22421 continue                                                             4185
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 4185
      if(jerr.ne.0) return                                                 4185
22431 continue                                                             4186
22411 continue                                                             4186
      dev1=0.0                                                             4187
22440 do 22441 ic=1,nc                                                     4187
      q(:,ic)=b(0,ic)+g(:,ic)                                              4188
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             4189
      q(:,ic)=exp(q(:,ic))                                                 4189
      sxp=sxp+q(:,ic)                                                      4190
22441 continue                                                             4191
22442 continue                                                             4191
      sxpl=w*log(sxp)                                                      4191
22450 do 22451 ic=1,nc                                                     4191
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  4191
22451 continue                                                             4192
22452 continue                                                             4192
22391 continue                                                             4193
22361 continue                                                             4193
22460 do 22461 ic=1,nc                                                     4193
22470 do 22471 i=1,no                                                      4193
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               4193
22471 continue                                                             4193
22472 continue                                                             4193
22461 continue                                                             4194
22462 continue                                                             4194
      dev0=dev0+dev1                                                       4196
      alf=1.0                                                              4198
      if(flmin .ge. 1.0)goto 22491                                         4198
      eqs=max(eps,flmin)                                                   4198
      alf=eqs**(1.0/(nlam-1))                                              4198
22491 continue                                                             4199
      m=0                                                                  4199
      mm=0                                                                 4199
      nin=0                                                                4199
      nlp=0                                                                4199
      mnl=min(mnlam,nlam)                                                  4199
      bs=0.0                                                               4200
      shr=shri*dev0                                                        4200
      ga=0.0                                                               4201
22500 do 22501 ic=1,nc                                                     4201
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      4201
      svr(ic)=sum(r(:,ic))                                                 4202
22510 do 22511 j=1,ni                                                      4202
      if(ju(j).eq.0)goto 22511                                             4203
      jb=ix(j)                                                             4203
      je=ix(j+1)-1                                                         4204
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             4205
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            4206
22511 continue                                                             4207
22512 continue                                                             4207
22501 continue                                                             4208
22502 continue                                                             4208
      ga=sqrt(ga)                                                          4209
22520 do 22521 ilm=1,nlam                                                  4209
      al0=al                                                               4210
      if(flmin .lt. 1.0)goto 22541                                         4210
      al=ulam(ilm)                                                         4210
      goto 22531                                                           4211
22541 if(ilm .le. 2)goto 22551                                             4211
      al=al*alf                                                            4211
      goto 22531                                                           4212
22551 if(ilm .ne. 1)goto 22561                                             4212
      al=big                                                               4212
      goto 22571                                                           4213
22561 continue                                                             4213
      al0=0.0                                                              4214
22580 do 22581 j=1,ni                                                      4214
      if(ju(j).eq.0)goto 22581                                             4214
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            4214
22581 continue                                                             4215
22582 continue                                                             4215
      al0=al0/max(bta,1.0d-3)                                              4215
      al=alf*al0                                                           4216
22571 continue                                                             4217
22531 continue                                                             4217
      al2=al*omb                                                           4217
      al1=al*bta                                                           4217
      tlam=bta*(2.0*al-al0)                                                4218
22590 do 22591 k=1,ni                                                      4218
      if(iy(k).eq.1)goto 22591                                             4218
      if(ju(k).eq.0)goto 22591                                             4219
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      4220
22591 continue                                                             4221
22592 continue                                                             4221
10880 continue                                                             4222
22600 continue                                                             4222
22601 continue                                                             4222
      ixx=0                                                                4222
      jxx=ixx                                                              4222
      kxx=jxx                                                              4222
      t=0.0                                                                4223
22610 do 22611 ic=1,nc                                                     4223
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       4223
22611 continue                                                             4224
22612 continue                                                             4224
      if(t .ge. eps)goto 22631                                             4224
      kxx=1                                                                4224
      goto 22602                                                           4224
22631 continue                                                             4224
      t=2.0*t                                                              4224
      alt=al1/t                                                            4224
      al2t=al2/t                                                           4225
22640 do 22641 ic=1,nc                                                     4225
      bs(0,ic)=b(0,ic)                                                     4225
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          4226
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    4226
      svr(ic)=sum(r(:,ic))                                                 4227
      if(intr .eq. 0)goto 22661                                            4227
      b(0,ic)=b(0,ic)+svr(ic)                                              4227
      r(:,ic)=r(:,ic)-svr(ic)*w                                            4228
      dlx=max(dlx,svr(ic)**2)                                              4229
22661 continue                                                             4230
22641 continue                                                             4231
22642 continue                                                             4231
22670 continue                                                             4231
22671 continue                                                             4231
      nlp=nlp+nc                                                           4231
      dlx=0.0                                                              4232
22680 do 22681 k=1,ni                                                      4232
      if(iy(k).eq.0)goto 22681                                             4233
      jb=ix(k)                                                             4233
      je=ix(k+1)-1                                                         4233
      del=b(k,:)                                                           4233
      gkn=0.0                                                              4234
22690 do 22691 ic=1,nc                                                     4235
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))-svr(ic)*xb(k))/xs(k)        4236
      gk(ic)=u+del(ic)*xv(k)                                               4236
      gkn=gkn+gk(ic)**2                                                    4237
22691 continue                                                             4238
22692 continue                                                             4238
      gkn=sqrt(gkn)                                                        4238
      u=1.0-alt*vp(k)/gkn                                                  4239
      if(u .gt. 0.0)goto 22711                                             4239
      b(k,:)=0.0                                                           4239
      goto 22721                                                           4240
22711 continue                                                             4241
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     4242
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(   4244 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 4245
22721 continue                                                             4246
22701 continue                                                             4246
      del=b(k,:)-del                                                       4246
      if(maxval(abs(del)).le.0.0)goto 22681                                4247
22730 do 22731 ic=1,nc                                                     4247
      dlx=max(dlx,xv(k)*del(ic)**2)                                        4248
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   4250 
     *b(k))/xs(k)
22731 continue                                                             4251
22732 continue                                                             4251
      if(mm(k) .ne. 0)goto 22751                                           4251
      nin=nin+1                                                            4252
      if(nin .le. nx)goto 22771                                            4252
      jxx=1                                                                4252
      goto 22682                                                           4252
22771 continue                                                             4253
      mm(k)=nin                                                            4253
      m(nin)=k                                                             4254
22751 continue                                                             4255
22681 continue                                                             4256
22682 continue                                                             4256
      if(jxx.gt.0)goto 22672                                               4257
      if(dlx.lt.shr)goto 22672                                             4257
      if(nlp .le. maxit)goto 22791                                         4257
      jerr=-ilm                                                            4257
      return                                                               4257
22791 continue                                                             4258
22800 continue                                                             4258
22801 continue                                                             4258
      nlp=nlp+nc                                                           4258
      dlx=0.0                                                              4259
22810 do 22811 l=1,nin                                                     4259
      k=m(l)                                                               4259
      jb=ix(k)                                                             4259
      je=ix(k+1)-1                                                         4259
      del=b(k,:)                                                           4259
      gkn=0.0                                                              4260
22820 do 22821 ic=1,nc                                                     4261
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))  -svr(ic)*xb(k))/xs(k)      4263
      gk(ic)=u+del(ic)*xv(k)                                               4263
      gkn=gkn+gk(ic)**2                                                    4264
22821 continue                                                             4265
22822 continue                                                             4265
      gkn=sqrt(gkn)                                                        4265
      u=1.0-alt*vp(k)/gkn                                                  4266
      if(u .gt. 0.0)goto 22841                                             4266
      b(k,:)=0.0                                                           4266
      goto 22851                                                           4267
22841 continue                                                             4268
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     4269
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(   4271 
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 4272
22851 continue                                                             4273
22831 continue                                                             4273
      del=b(k,:)-del                                                       4273
      if(maxval(abs(del)).le.0.0)goto 22811                                4274
22860 do 22861 ic=1,nc                                                     4274
      dlx=max(dlx,xv(k)*del(ic)**2)                                        4275
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   4277 
     *b(k))/xs(k)
22861 continue                                                             4278
22862 continue                                                             4278
22811 continue                                                             4279
22812 continue                                                             4279
      if(dlx.lt.shr)goto 22802                                             4279
      if(nlp .le. maxit)goto 22881                                         4279
      jerr=-ilm                                                            4279
      return                                                               4279
22881 continue                                                             4281
      goto 22801                                                           4282
22802 continue                                                             4282
      goto 22671                                                           4283
22672 continue                                                             4283
      if(jxx.gt.0)goto 22602                                               4284
22890 do 22891 ic=1,nc                                                     4285
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                               4286
      if(ixx .ne. 0)goto 22911                                             4287
22920 do 22921 j=1,nin                                                     4287
      k=m(j)                                                               4288
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 22941                   4288
      ixx=1                                                                4288
      goto 22922                                                           4288
22941 continue                                                             4290
22921 continue                                                             4291
22922 continue                                                             4291
22911 continue                                                             4292
      sc=b(0,ic)+g(:,ic)                                                   4292
      b0=0.0                                                               4293
22950 do 22951 j=1,nin                                                     4293
      l=m(j)                                                               4293
      jb=ix(l)                                                             4293
      je=ix(l+1)-1                                                         4294
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   4295
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            4296
22951 continue                                                             4297
22952 continue                                                             4297
      sc=min(max(exmn,sc+b0),exmx)                                         4298
      sxp=sxp-q(:,ic)                                                      4299
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          4300
      sxp=sxp+q(:,ic)                                                      4301
22891 continue                                                             4302
22892 continue                                                             4302
      s=sum(b(0,:))/nc                                                     4302
      b(0,:)=b(0,:)-s                                                      4303
      if(jxx.gt.0)goto 22602                                               4304
      if(ixx .ne. 0)goto 22971                                             4305
22980 do 22981 j=1,ni                                                      4305
      if(iy(j).eq.1)goto 22981                                             4305
      if(ju(j).eq.0)goto 22981                                             4305
      ga(j)=0.0                                                            4305
22981 continue                                                             4306
22982 continue                                                             4306
22990 do 22991 ic=1,nc                                                     4306
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      4307
23000 do 23001 j=1,ni                                                      4307
      if(iy(j).eq.1)goto 23001                                             4307
      if(ju(j).eq.0)goto 23001                                             4308
      jb=ix(j)                                                             4308
      je=ix(j+1)-1                                                         4309
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             4310
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            4311
23001 continue                                                             4312
23002 continue                                                             4312
22991 continue                                                             4313
22992 continue                                                             4313
      ga=sqrt(ga)                                                          4314
23010 do 23011 k=1,ni                                                      4314
      if(iy(k).eq.1)goto 23011                                             4314
      if(ju(k).eq.0)goto 23011                                             4315
      if(ga(k) .le. al1*vp(k))goto 23031                                   4315
      iy(k)=1                                                              4315
      ixx=1                                                                4315
23031 continue                                                             4316
23011 continue                                                             4317
23012 continue                                                             4317
      if(ixx.eq.1) go to 10880                                             4318
      goto 22602                                                           4319
22971 continue                                                             4320
      goto 22601                                                           4321
22602 continue                                                             4321
      if(kxx .le. 0)goto 23051                                             4321
      jerr=-20000-ilm                                                      4321
      goto 22522                                                           4321
23051 continue                                                             4322
      if(jxx .le. 0)goto 23071                                             4322
      jerr=-10000-ilm                                                      4322
      goto 22522                                                           4322
23071 continue                                                             4322
      devi=0.0                                                             4323
23080 do 23081 ic=1,nc                                                     4324
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          4324
      a0(ic,ilm)=b(0,ic)                                                   4325
23090 do 23091 i=1,no                                                      4325
      if(y(i,ic).le.0.0)goto 23091                                         4326
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           4327
23091 continue                                                             4328
23092 continue                                                             4328
23081 continue                                                             4329
23082 continue                                                             4329
      kin(ilm)=nin                                                         4329
      alm(ilm)=al                                                          4329
      lmu=ilm                                                              4330
      dev(ilm)=(dev1-devi)/dev0                                            4331
      if(ilm.lt.mnl)goto 22521                                             4331
      if(flmin.ge.1.0)goto 22521                                           4332
      me=0                                                                 4332
23100 do 23101 j=1,nin                                                     4332
      if(a(j,1,ilm).ne.0.0) me=me+1                                        4332
23101 continue                                                             4332
23102 continue                                                             4332
      if(me.gt.ne)goto 22522                                               4333
      if(dev(ilm).gt.devmax)goto 22522                                     4333
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 22522                             4334
22521 continue                                                             4335
22522 continue                                                             4335
      g=log(q)                                                             4335
23110 do 23111 i=1,no                                                      4335
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         4335
23111 continue                                                             4336
23112 continue                                                             4336
      deallocate(sxp,b,bs,r,q,mm,is,sc,ga,iy,gk,del,sxpl)                  4337
      return                                                               4338
      end                                                                  4339
      subroutine psort7 (v,a,ii,jj)                                             
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
