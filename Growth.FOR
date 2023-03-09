c***********************************************************************
c                         G R O W T H - 23            11/1/2023        *
c                                                                      *
c        Compilation with Intel Visual Fortran,                        *
c            (Quickwin Application Project)                            *
c                                                                      *
c  Aim:      3-D gravity inversion for modelling gravity dg or CBA     *
c               as aggregated growing bodies for anomalous density     *
c                                                                      *
c  Reference: A free-geometry geodynamic modelling of surface gravity  *
c                changes using Growth-dg software                      *
c             A.G.Camacho, P.Vajda, J.Fernández                        *
c ---------------------------------------------------------------------*
c                                                                      *
c  INPUT                                                               *
c  -----                                                               *
c  file "Gra.txt": coordinates, altitudes, observed gravity values     *
c         xp,yp,zp : planar coordinates (metres) and altitude (metres) *
c         dg: observed gravity values, dg or CBA (in microgals)        *
c                                                                      *
c  OUTPUT                                                              *
c  ------                                                              *
c  file "mod.txt": Resulting 3-D model of the subsurface anomalous     *
c        structures as given by means of the 3-D grid of cells         *
c        Values for several characteristic parameters of the model     *
c                                                                      *
c  file "FIL.txt": regional, local, calculated and residual gravity    *
c        values resulting from inversion                               *
c                                                                      *
c ---------------------------------------------------------------------*
c  DIMENSIONS                                                          *
c  ----------                                                          *
c   ms:  max. number of gravity stations                               *
c   mc:  max. number of prismatic cells for the 3-D grid               *
c   mb:  number of points in file map.bln for map lines                *
c                                                                      *
c***********************************************************************

      use IFQWIN
      use IFPORT

c  Dimensions

      parameter (ms=5000,mc=99999,mb=2500,ml=35,mr=50,mpo=200)
      type (xycoord) xy
      dimension 
     & xc(mc),yc(mc),zc(mc),tx(mc),ty(mc),tz(mc),df(mc),far(mc),
     & qm(mc),isen(mc),vol(mc),lev(mc),lay(mc), iup(mpo,mc),ido(mpo,mc), 
     & rr(ms),g(ms),gc(ms),x(ms),y(ms),z(ms),re(ms),w(ms),
     & xn(ms),yn(ms),zn(ms),xk(ms),yk(ms), klev(2,ml) ,
     & cov(mr),swc(mr),dla(mr), xb(mb),yb(mb),ib(mb)
      integer col(21),  nada
      integer*2 u1,u2,v1,v2,up,vp,up2,vp2
      real*8 gxx,gyy,gxg,gyg,gxz,gyz,gzz,gzg,ryy,ryg,
     & px,py,pz,f1g,f11,srg,su,dg,f  
      character*50 obs,mod,fil
      character*120 texto
      character*9 hoy
      character*1 ts1,ts2,ts3,ts4,ts5,ts6,bl,as,le(15)
      data bl/' '/,as/'*'/,
     & le/'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o'/

      call seed(2016)
      nada=setbkcolorrgb(#f0f0f0)
      call clearscreen($gclearscreen)
      nada=initializefonts()
      nada=settextcolorrgb(#000000)

c    Default values for parameters

      ig0 =   0
      igxy=   0
      bal =  30
      nra =   4
      blun= 2.5
      mlev=   4
      tpc=    3
      ncer=   0
      kup=    0
      seli=  50  
      fdib= 106  
      pp  = 750
      base=-14000
      nlim=30

c   Reading file names

      call dial1(obs,mod,fil,nra,bal,tpc,kup,ncer,mlev)

c   Reading position of bechmarks and gravity value

      open(1,file=obs)
      ax= 9.d9
      bx=-ax
      ay= ax
      by=-ay
      xd=0.
      yd=0.
      zd=0.
      gm=0.
      i=0
      ns=0
    1 read(1,*,err=2,end=2) xp,yp,zp,dg 
      eg=1
      if(xp.eq.0.and.yp.eq.0.and.zp.eq.0) go to 2
      if(xp.lt.ax) ax=xp
      if(xp.gt.bx) bx=xp
      if(yp.lt.ay) ay=yp
      if(yp.gt.by) by=yp
      i=i+1
      if(i.ge.ms) call dim(ms)
      x(i)=xp
      y(i)=yp
      z(i)=zp  
      g(i)=dg
      w(i)=eg   
      xd=xd+xp
      yd=yd+yp
      zd=zd+zp
      gm=gm+dg
      go to 1
    2 ns=i 
      gm=gm/ns 
      ig0 =0
      igxy=0 
      if(ns.le.2) then
        write(*,'(//a,a)') ' ***** Error in data file ',obs
        stop
      endif
      close(1)

c    Normalized coordinates
      
      fe=(bx-ax+by-ay)/2./30000. 
      ixm= (xd/ns+ax+bx)/3.
      iym= (yd/ns+ay+by)/3.
      izm= zd/ns
      jmi=1
      if(fe.gt.0.1) then
        ixm= nint(ixm/100.)*100.
        iym= nint(iym/100.)*100.
        izm= nint(izm/10.)*10.
        jmi=0
      endif     

      tic=3000.*fe     
      call round(tic,r)   
      tic=r
    
      techo=-9.d9

      do 3 i=1,ns
      x(i)=x(i)-ixm         
      y(i)=y(i)-iym
      z(i)=z(i)-izm
      zn(i)=z(i)/fe
      if(zn(i).gt.techo) techo=zn(i)  
      w(i)=1.
      xk(i)=x(i)/1000.
      yk(i)=y(i)/1000.  
      xn(i)=x(i)/fe
    3 yn(i)=y(i)/fe

c    Coordinate limits

      ax=(ax-ixm)/fe
      bx=(bx-ixm)/fe
      ay=(ay-iym)/fe
      by=(by-iym)/fe

c    Read geometrical parameters for the 3D grid

      iex=3
      cell=pp*fe
   10 call dial2(obs,mod,fil,nra,ns,bal,tpc,kup,ncer,mlev,cell,
     - jin,jmi)
      pp=cell/fe*0.55                                       
     
c   Creating a 3D grid of cells (for normalized coordinates)

      avm=pp*pp*pp*pp*1.0d-10   
      nc=0
      kr=0
      write(*,'(//////10x,a\)') 'Please, wait ...'
      nlay=1
    7 ncr=nc
      kr=kr+1
      write(*,'(a1\)') ' '
      zs=techo
      ax=ax-4000
      ay=ay-4000
      bx=bx+4000
      by=by+4000
    8 dz=pp*0.3                
      call especapa(zs,avm,ms,ns,xn,yn,zn,dz) 
      if(dz.gt.1) go to 9
        zs=zs-pp*0.5
        if(zs.gt.base) go to 8
    9 if(kr.eq.2) dz=dz*1.7   
      if(kr.eq.3) dz=dz*2.0  
      write(*,'(a1\)') '.'
      if(kr.eq.1.and.nc.gt.0) nlay=nlay+1
      j=nlay
      if(kr.gt.1) j=0   
      call celdascapa(ax,bx,ay,by,zs,dz,avm,ms,ns,xn,yn,zn,
     & mc,nc,xc,yc,zc,tx,ty,tz,df,seli,ncr,lay,j)
      if(nc.ge.mc.or.nc.lt.0) then 
        write(*,'(////80x,a,i6)') '***** Warning: Num. of cells >',mc
        write(*,'(/84x,a,i6)') '>> Try a smaller length for cells'
        pp=680
        go to 10
      endif
      zs=zs-dz                      
      if(zs.gt.base) go to 8  
      if(kr.lt.iex) go to 7   
      c=0
      d=0
      mlay=0
      do 4 i=1,nc
      xc(i)=xc(i)*fe
      yc(i)=yc(i)*fe
      zc(i)=zc(i)*fe
      tx(i)=tx(i)*fe
      ty(i)=ty(i)*fe
      tz(i)=tz(i)*fe
      vol(i)=tx(i)*ty(i)*tz(i)
      far(i)=0

      if(lay(i).eq.0) far(i)=-999
      c=c+tx(i)+ty(i)+tz(i)
      d=d+abs(tx(i)-ty(i))+abs(tz(i)-ty(i))+abs(tz(i)-tx(i))
      if(lay(i).gt.mlay) mlay=lay(i)
    4 continue
      c=c/3./nc
      d=d/3./nc
      siz=nint(c)
      techo=techo*fe 

c   Inner parameters

      if(jin.eq.1) call dial4(blun,nlim)

      if(ns.gt.nlim) then
        ig0 =  1
        igxy=  1
        bal = 50
        nra = 80
        tpc = 60
        ncer= 10 
        kup =  5
        mlev= 30
      endif 

c  Reading map.bln
      
      open(1,file='map.bln') 
      k=1
      l=2
      lx=192
    6 read(1,*,end=11) j
      do 5 i=k,(j+k-1)
      read(1,*) xp,yp
      ix=(xp-ixm)/fe
      iy=(yp-iym)/fe
      xb(i)=ix
      yb(i)=iy
      ix=ix/fdib
      iy=iy/fdib
      if(ix.gt.lx.or.ix.lt.-lx.or.iy.gt.lx.or.iy.lt.-lx) l=1
      if(i.ge.mb) call dim(mb)
      ib(i)=l
      l=0
    5 continue
      k=k+j
      l=2
      go to 6
   11 close(1)
      nb=k-1

c    Reading options and parameters for inversion

      pp=cell
      call dial3(obs,mod,fil,nra,nc,ns,bal,mlev,pp,
     & ig0,igxy,tpc,kup,ncer,jmi,jin)

c    Some initial values 

      fbal=bal*bal*fe*fe*1000   
      gmi=9.d9
      gma=-9.d9
      disg=0
      do 12 i=1,ns         
      if(g(i).lt.gmi) gmi=g(i)
      if(g(i).gt.gma) gma=g(i)
      disg=disg+(g(i)-gm)**2
   12 continue
      disg=sqrt(disg/ns)
      kdg=10*nint(disg/20./10.)
      if(disg.gt.1000) kdg=20.*nint(disg/20./20.)
      if(ns.le.30) kdg=disg/10.
      if(kdg.lt.1) kdg=1
      rag=gma-gmi
      deg=gma*gmi
      if(ncer.le.0) mlev=4

c    Contiguous cells
 
      if(ncer.gt.0) then
            write(*,'(//10x,a\)') 'Please, wait again '
            call updown(mc,nc,xc,yc,tz,lay,mpo,iup,ido)
            if(ncer.gt.mpo) ncer=mpo
      endif

c  Calculating design matrix and sensitivity

      avm=siz/50.  
c      c=1.*mlay/mlev
      do 14 j=1,nc
      sv=0.
      pen=vol(j)*6.672e-3
      do 13 i=1,ns
      xd=x(i)-xc(j)
      yd=y(i)-yc(j)
      zd=z(i)-zc(j)                            
      d=xd*xd+yd*yd+zd*zd + 1.e2*fe*fe 
   13 sv= sv+zd*zd/d/d/d
      sv=sv/ns
      qm(j)=sv*pen*pen/avm/avm    
c      if(lay(j).gt.0) lay(j)=lay(j)/c+1.
   14 isen(j)=1.d-7/sqrt(sv)/fe/fe     


c--------Initial graphics and headings ----------------------------------
      jx2=1025                                                          !
      jx3=1433                                                          !
      jx4=1761                                                          !
      jy4=555                                                           !
      jy5=760                                                           !
      jy6=140                                                           !
      jy7=370                                                           !
      call clearscreen($gclearscreen)                                   !
      nada=setcolorrgb(#000000)                                         !
      nada=initializefonts()                                            !
      write(texto,'(a)') 'GROWTH-23'                                    !
      nada=setfont('t''Times New Roman''h25w12')                        !
      call moveto(130,15,xy)                                            !
      call outgtext(texto)                                              !
      kdib=0                                                            !
      write(texto,'(a)')                                                !
     & 'Inversion of gravity values as 3D mass sources '                !
      nada=setfont('t''Times New Roman''h19w8')                         !
      call moveto(26,42,xy)                                             !
      call outgtext(texto)                                              !
      call date(hoy)                                                    !
      write(texto,'(a)') hoy                                            !
      call ventana(160,252,240,265)                                     !
      call moveto(170,250,xy)                                           !
      nada=setcolorrgb(#000000)                                         !
      nada=setfont('t''Courier''h17w7')                                 !
      call outgtext(texto)                                              !

      call ventana(18,124,380,142)                                      !
      nada=setcolorrgb(#a0ffff)                                         !
      nada=rectangle($gfillinterior,302,125,342,141)                    !
      nada=setcolorrgb(#000000)                                         !
      write(texto,*)                                                    !
     & 'Num    Num   Side  Blun  Bala  Dep  Fla  Size  Rand'            !
      call moveto(18,90,xy)                                             !
      call outgtext(texto)                                              !
      write(texto,*)                                                    !
     & 'data   cells  (m)   val   val  val  val   %     val'            !
      call moveto(18,103,xy)                                            !
      call outgtext(texto)                                              !
      if(jmi.eq.0) write(texto,'(i5,i7,i6,f6.1,3i5,f7.1,i5)')           !
     & ns,nc,nint(cell),blun,nint(bal),kup,ncer,tpc,nra                 !
      if(jmi.eq.1) write(texto,'(i5,i7,f6.1,f6.1,3i5,f7.1,i5)')         !
     & ns,nc,cell,blun,nint(bal),kup,ncer,tpc,nra                       !
      call moveto(20,126,xy)                                            !
      call outgtext(texto)                                              !

      write(texto,'(a,3x,a)') 'Number of',                              !
     -' Dens   Resid  Filled  Off.G    Reg.G '                          !
      call moveto(25,165,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a,a)') 'Steps  Cells',                              !
     -'  kg/m3   uGal     %     uG     uG/km '                          !
      call moveto(20,180,xy)                                            !
      call outgtext(texto)                                              !
      call ventana(18,202,380,220)                                      !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)') 'Initial steps . . . . . .'                    !
      call moveto(40,204,xy)                                            !
      call outgtext(texto)                                              !
      
      ix=40                                                             !
      iy=300                                                            !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a,i4,a)') 'Aver.reliev ',izm,' m'                   !
      call moveto(ix+53,iy-2,xy)                                        !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Sea level'                                    !
      call moveto(ix+53,iy+22,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Neg.dens.contr.'                              !
      call moveto(ix+53,iy+46,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Pos.dens.contr.'                              !
      call moveto(ix+53,iy+69,xy)                                       !
      call outgtext(texto)                                              !
      call ventana(ix-20,iy-15,ix+42,iy+100)                            !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=ix-20                                                          !
      u2=ix+42                                                          !
      v1=iy+5                                                           !
      v2=iy+102                                                         !
      nada=rectangle($gfillinterior,u1,v1 ,u2,v2)                       !
      nada=setcolorrgb(#f57070)                                         !
      u1=ix-11                                                          !
      u2=ix-3                                                           !
      v1=iy+50                                                          !
      v2=iy+57                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#f06060)                                         !
      u1=ix+1                                                           !
      v1=iy+50                                                          !
      u2=ix+10                                                          !
      v2=iy+57                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#c05050)                                         !
      u1=ix+14                                                          !
      v1=iy+50                                                          !
      u2=ix+22                                                          !
      v2=iy+57                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#804040)                                         !
      u1=ix+26                                                          !
      v1=iy+50                                                          !
      u2=ix+34                                                          !
      v2=iy+57                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#7070f7)                                         !
      u1=ix-11                                                          !
      v1=iy+74                                                          !
      u2=ix-3                                                           !
      v2=iy+81                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#6060f0)                                         !
      u1=ix+1                                                           !
      v1=iy+74                                                          !
      u2=ix+10                                                          !
      v2=iy+81                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#5050c0)                                         !
      u1=ix+14                                                          !
      v1=iy+74                                                          !
      u2=ix+22                                                          !
      v2=iy+81                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#404080)                                         !
      u1=ix+26                                                          !
      v1=iy+74                                                          !
      u2=ix+34                                                          !
      v2=iy+81                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#f0f000)                                         !
      call moveto(ix-20,iy+30,xy)                                       !
      u1=ix+42                                                          !
      v1=iy+30                                                          !
      nada=lineto(u1,v1)                                                !

      jx1=612                                                           !
      nada=setcolorrgb(#000000)                                         !
      call setgtextrotation(900)                                        !
      write(texto,'(a,i6,a)') 'Tic=',nint(tic),' m'                     !
      call moveto(jx1-215,610,xy)                                       !
      call outgtext(texto)                                              !
      call moveto(jx1+196,220,xy)                                       !
      call outgtext(texto)                                              !
      call moveto(jx2+195,230,xy)                                       !
      call outgtext(texto)                                              !
      call moveto(jx3+195,230,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'A.G. Camacho 2023'                            !
      call moveto(jx1+195,410,xy)                                       !
      call outgtext(texto)                                              !
      call setgtextrotation(000)                                        !

      jz=445                                                            ! 
      jx=jx1                                                            !
      jy0=220                                                           !
      lx=190                                                            !
      ly=190                                                            !
      lz=172                                                            !
      call ventana(jx-lx,jz,jx+lx,jz+lz)                                !
      call ventana(jx-lx,jy0-ly,jx+lx,jy0+ly)                           !
c-------------sensitivity scale-----------------------------------------!
      jx0=205                                                           !
      call ventana(jx0-lx,jz,jx0+lx,jz+lz)                              !
      col( 1)=#76a076                                                   !
      col( 2)=#80a680                                                   !
      col( 3)=#86b086                                                   !
      col( 4)=#90b690                                                   !
      col( 5)=#96c096                                                   !
      col( 6)=#a0c6a0                                                   !
      col( 7)=#a6d0a6                                                   !
      col( 8)=#b0d6b0                                                   !
      col( 9)=#b6e0b6                                                   !
      col(10)=#c0e6c0                                                   !
      col(11)=#c6f0c6                                                   !
      col(12)=#d0f6d0                                                   !
      jx=340                                                            !
      jy=270                                                            !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)') 'Sens.'                                        !
      call moveto(jx,jy-25,xy)                                          !
      call outgtext(texto)                                              !
      c=1.e-8/fe/fe/fe/2.                                               !
      call setgtextrotation(900)                                        !
      write(texto,'(a,e8.1,a)') 'x',c,' kg/uGal'                        !
      call moveto(jx+27,jy+130,xy)                                      !
      call outgtext(texto)                                              !
      call setgtextrotation(000)                                        !
      call ventana(jx+15,jy,jx+25,jy+142)                               !
      do 21 i=1,12                                                      !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(i5)') i                                             !
      call moveto(jx-25,jy+(i-1)*12-3,xy)                               !
      call outgtext(texto)                                              !
      nada=setcolorrgb(col(i))                                          !
      u1=jx+15                                                          !
      v1=jy+(i-1)*12                                                    !
      u2=jx+25                                                          !
      v2=jy+i*12                                                        !
   21 nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      
      jx=jx1                                                            !
      jy=jy0                                                            !
      nada=setcolorrgb(#c0d0c0)                                         !
      u1=jx-lx                                                          !
      v1=jy-ly                                                          !
      u2=jx+lx                                                          !
      v2=jy+ly                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !

      zz=techo-8000*fe                                                  !
      do 22 i=1,39                                                      !
      ix=-lx+lx/20.*i                                                   !
      xx=ix*fe*fdib                                                     !
      do 22 l=1,39                                                      !
      iy=-ly+ly/20.*l                                                   !
      yy=iy*fe*fdib                                                     !
      call ss2(ms,ns,x,y,z,fe,xx,yy,zz,sen,k)                           !
      nada=setcolorrgb(col(k))                                          !
      k=10                                                              !
      u1=jx+ix-k                                                        !
      v1=jy-iy-k                                                        !
      u2=jx+ix+k                                                        !
      v2=jy-iy+k                                                        !
   22 nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !

      nada=setcolorrgb(#c0d0c0)                                         !
      u1=jx-lx                                                          !
      v1=jz+30                                                          !
      u2=jx+lx                                                          !
      v2=jz+lz                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !

      yy=0                                                              !
      do 23 i=1,39                                                      !
      ix=-lx+lx/20.*i                                                   !
      xx=ix*fe*fdib                                                     !
      do 23 l=1,14                                                      !
      iz= ly/20.*l                                                      !
      zz=techo-(iz+30)*fe*fdib                                          !
      call ss2(ms,ns,x,y,z,fe,xx,yy,zz,sen,k)                           !
      nada=setcolorrgb(col(k))                                          !
      k=10                                                              !
      u1=jx+ix-k                                                        !
      v1=jz+30+iz-k                                                     !
      u2=jx+ix+k                                                        !
      v2=jz+30+iz+k                                                     !
   23 nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#c0d0c0)                                         !
      jy=jx0                                                            !
      u1=jy-ly                                                          !
      v1=jz+30                                                          !
      u2=jy+ly                                                          !
      v2=jz+lz                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      xx=0                                                              !
      do 24 i=1,39                                                      !
      iy=-ly+ly/20.*i                                                   !
      yy=iy*fe*fdib                                                     !
      do 24 l=1,14                                                      !
      iz= ly/20.*l                                                      !
      zz=techo-(iz+30)*fe*fdib                                          !
      call ss2(ms,ns,x,y,z,fe,xx,yy,zz,sen,k)                           !
      nada=setcolorrgb(col(k))                                          !
      k=10                                                              !
      u1=jy+iy-k                                                        !
      v1=jz+30+iz-k                                                     !
      u2=jy+iy+k                                                        !
      v2=jz+30+iz+k                                                     !
   24 nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
c------------------Offset, misfit evoluton------------------------------!
      call ventana(20,868,590,962)                                      !
      call ventana(615,868,1215,962)                                    !
      nada=setcolorrgb(#202020)                                         !
      nada=rectangle($gfillinterior,80,868,80,962)                      !
      nada=rectangle($gfillinterior,80,920,590,920)                     !
      write(texto,'(a)') ' Grav.'                                       !
      call moveto(30,873,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'offset'                                       !
      call moveto(30,890,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a)')                                                !
     -'500   1000  1500  2000  2500  3000  3500  4000  4500  5000  5500'!
      call moveto(108,930,xy)                                           !
      call outgtext(texto)                                              !
      write(texto,'(a,a)') '5000 10000 ',                               !
     -'20000 25000 30000 35000 40000 45000 50000 55000 60000 65000'     !
      call moveto(703,930,xy)                                           !
      call outgtext(texto)                                              !
      do 27 i=1,12                                                      !
      c= 500*i/12.                                                      !
      nada=rectangle($gfillinterior,80+c,920,80+c,925)                  !
   27 nada=rectangle($gfillinterior,675+c,920,675+c,925)                !
      write(texto,'(a)') 'Growth steps'                                 !
      call moveto(500,892,xy)                                           !
      call outgtext(texto)                                              !
      call moveto(1120,880,xy)                                          !
      call outgtext(texto)                                              !
      write(texto,'(i4,a)') kdg,' uG'                                   !
      call moveto(23,927,xy)                                            !
      call outgtext(texto)                                              !
      nada=setcolorrgb(#008000)                                         !
      nada=rectangle($gfillinterior,79,920,81,950)                      !
      nada=setcolorrgb(#202020)                                         !
      nada=rectangle($gfillinterior,680,868,680,962)                    !
      nada=rectangle($gfillinterior,680,920,1217,920)                   !
      write(texto,'(a)') ' R.M.S.'                                      !
      call moveto(623,873,xy)                                           !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'residuals'                                    !
      call moveto(616,890,xy)                                           !
      call outgtext(texto)                                              !
      j=3*kdg                                                           !
      if(ns.le.30) j=6*kdg                                              !
      write(texto,'(i4,a)') j,' uG'                                     !
      call moveto(623,927,xy)                                           !
      call outgtext(texto)                                              !
      nada=setcolorrgb(#008000)                                         !
      nada=rectangle($gfillinterior,679,920,681,950)                    !
c--------------------------Obser,model,res------------------------------!
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a,50x,a,34x,a)')                                    !
     &   'Obser dg','Model dg','Regional & Local'                       !
      call moveto(jx2-20,10,xy)                                         !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Residual values'                              !
      call moveto(jx3-32,420,xy)                                        !
      call outgtext(texto)                                              !
      call moveto(jx2-55,jy5-85,xy)                                     !
      call outgtext(texto)                                              !
      jy=jy0                                                            !
      call ventana(jx2-lx,jy-ly,jx2+lx,jy+ly)                           !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx2-lx                                                         !
      v1=jy-ly                                                          !
      u2=jx2+lx                                                         !
      v2=jy+ly                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      call ventana(jx3-lx,jy-ly,jx3+lx,jy+ly)                           !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx3-lx                                                         !
      v1=jy-ly                                                          !
      u2=jx3+lx                                                         !
      v2=jy+ly                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      jy= jz+ly                                                         !
      call ventana(jx3-lx,jy-ly,jx3+lx,jy+ly)                           !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx3-lx                                                         !
      v1=jy-ly                                                          !
      u2=jx3+lx                                                         !
      v2=jy+ly                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
c--------------------------Regional and local---------------------------!
      lxr=110                                                           !
      call ventana(jx4-lxr,jy6-lxr,jx4+lxr,jy6+lxr)                     !
      call hor0(jx4,jy6,lxr,lxr,tic,fe,fdib*1.7,mb,nb,xb,yb,ib)         !
      call ventana(jx4-lxr,jy7-lxr,jx4+lxr,jy7+lxr)                     !
      call hor0(jx4,jy7,lxr,lxr,tic,fe,fdib*1.7,mb,nb,xb,yb,ib)         !
c------------------------- sea level -----------------------------------!     
      v1=jz+30+izm/fe/fdib                                              !
      nada=setcolorrgb(#f0f000)                                         !
      jx=jx1                                                            !
      call moveto(jx-lx,v1,xy)                                          !
      u1=jx+lx                                                          !
      nada=lineto(u1,v1)                                                !
      jx=jx0                                                            !
      call moveto(jx-lx,v1,xy)                                          !
      u1=jx+lx                                                          !
      nada=lineto(u1,v1)                                                !
c----------------------------Graphical scale----------------------------!
      ix=240                                                            !
      iy=355                                                            !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(i5,a)') nint(tic),' m'                              !
      call moveto(ix-3,iy-2,xy)                                         !
      call outgtext(texto)                                              !
      v1=iy+28                                                          !
      call moveto(ix+3,v1,xy)                                           !
      u1=ix+47                                                          !
      nada=lineto(u1,v1)                                                !
      u1=ix+10                                                          !
      call moveto(u1,v1-5,xy)                                           !
      nada=lineto(u1,v1)                                                !
      u1=ix+10+tic/fe/fdib                                              !
      call moveto(u1,v1-5,xy)                                           !
      nada=lineto(u1,v1)                                                !
      u1=ix-6                                                           !
      call moveto(u1,iy-17,xy)                                          !
      v1=iy+24                                                          !
      nada=lineto(u1,v1)                                                !
      v1=iy+20                                                          !
      call moveto(u1+5,v1,xy)                                           !
      nada=lineto(u1,v1)                                                !
      v1=iy+20-tic/fe/fdib                                              !
      call moveto(u1+5,v1,xy)                                           !
      nada=lineto(u1,v1)                                                !
c--------------------------axes tics------------------------------------!
      j=izm/tic+2                                                       !
      do 25 i=1,20                                                      !

      jx=jx1                                                            !
      jy=jy0                                                            !
      ix=jx-lx                                                          !
      k=(i-1)*tic/fe/fdib                                               !
      v1=jy-ly+k                                                        !
      call moveto(ix,v1,xy)                                             !
      u1=ix+4                                                           !
      if(k.lt.2*ly) nada=lineto(u1,v1)                                  !
      ix=jx+lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix-4                                                           !
      if(k.lt.2*ly) nada=lineto(u1,v1)                                  !
      ix=jx-lx+k                                                        !
      iy=jy-ly                                                          !
      call moveto(ix,iy,xy)                                             !
      u1=ix                                                             !
      v1=iy+4                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
      iy=jy+ly                                                          !
      call moveto(ix,iy,xy)                                             !
      v1=iy-4                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !

      jx=jx2                                                            !
      jy=jy0                                                            !
      ix=jx-lx                                                          !
      k=(i-1)*tic/fe/fdib                                               !
      v1=jy-ly+k                                                        !
      call moveto(ix,v1,xy)                                             !
      u1=ix+4                                                           !
      if(k.lt.2*ly) nada=lineto(u1,v1)                                  !
      ix=jx+lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix-4                                                           !
      if(k.lt.2*ly) nada=lineto(u1,v1)                                  !
      ix=jx-lx+k                                                        !
      iy=jy-ly                                                          !
      call moveto(ix,iy,xy)                                             !
      u1=ix                                                             !
      v1=iy+4                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
      iy=jy+ly                                                          !
      call moveto(ix,iy,xy)                                             !
      v1=iy-4                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
     
      jx=jx3                                                            !
      jy=jy0                                                            !
      ix=jx-lx                                                          !
      k=(i-1)*tic/fe/fdib                                               !
      v1=jy-ly+k                                                        !
      u1=ix                                                             !
      call moveto(u1,v1,xy)                                             !
      u1=ix+4                                                           !
      if(k.lt.2*ly) nada=lineto(u1,v1)                                  !
      ix=jx+lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix-4                                                           !
      if(k.lt.2*ly) nada=lineto(u1,v1)                                  !
      ix=jx-lx+k                                                        !
      iy=jy-ly                                                          !
      call moveto(ix,iy,xy)                                             !
      u1=ix                                                             !
      v1=iy+4                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
      iy=jy+ly                                                          !
      call moveto(ix,iy,xy)                                             !
      v1=iy-4                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
      
      k=((i-j)*tic+izm)/fe/fdib                                         !
      v1=jz+30+k                                                        !
      jx=jx1                                                            !
      ix=jx-lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix+4                                                           !
      if(v1.lt.(jz+lz)) nada=lineto(u1,v1)                              !
      ix=jx+lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix-4                                                           !
      if(v1.lt.(jz+lz)) nada=lineto(u1,v1)                              !
      k=(i-1)*tic/fe/fdib                                               !
      ix=jx-lx+k                                                        !
      call moveto(ix,jz,xy)                                             !
      u1=ix                                                             !
      v1=jz+5                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !

      k=((i-j)*tic+izm)/fe/fdib                                         !
      v1=jz+30+k                                                        !
      jx=jx0                                                            !
      ix=jx-lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix+4                                                           !
      if(v1.lt.(jz+lz)) nada=lineto(u1,v1)                              !
      ix=jx+lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix-4                                                           !
      if(v1.lt.(jz+lz)) nada=lineto(u1,v1)                              !
      k=(i-1)*tic/fe/fdib                                               !
      ix=jx-lx+k                                                        !
      call moveto(ix,jz,xy)                                             !
      u1=ix                                                             !
      v1=jz+5                                                           !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
      
   25 continue                                                          !
c-----------------------------------------------------------------------!      
      write(texto,'(a)') 'Plan view'                                    !
      call moveto(jx1-20,10,xy)                                         !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Elevation view W - E'                         !
      call moveto(jx1-60,420,xy)                                        !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Elevation view S - N'                         !
      call moveto(jx0-50,420,xy)                                        !
      call outgtext(texto)                                              !
c---------------color scales--------------------------------------------!
      col( 1)=#600000                                                   !
      col( 2)=#a00000                                                   !
      col( 3)=#d00000                                                   !
      col( 4)=#f06000                                                   !
      col( 5)=#f09000                                                   !
      col( 6)=#f0c000                                                   !
      col( 7)=#f0e000                                                   !
      col( 8)=#f0f060                                                   !
      col( 9)=#f0f0b0                                                   !
      col(10)=#f3f3d0                                                   !
      col(11)=#f6f6f6                                                   !
      col(12)=#a0f2f2                                                   !
      col(13)=#70e9f0                                                   !
      col(14)=#00e0f0                                                   !
      col(15)=#00c0f0                                                   !
      col(16)=#0090f0                                                   !
      col(17)=#0060f0                                                   !
      col(18)=#0000f0                                                   !
      col(19)=#0000d0                                                   !
      col(20)=#0000a0                                                   !
      col(21)=#000060                                                   !
      ix =jx3+320                                                       !
      ix2=ix                                                            !
      iy2=520                                                           ! 
      iy =iy2                                                           !
      pg=rag/20                                                         !
      dg=0                                                              !
      if(abs(gmi).gt.100) then                                          !
       dg=gmi+10*pg                                                     !
       write(texto,'(a)') '+'                                           !
       call moveto(ix-20,iy+268,xy)                                     !
       call outgtext(texto)                                             !
       write(texto,'(i7)') nint(dg)                                     !
       call moveto(ix-55,iy+284,xy)                                     !
       call outgtext(texto)                                             !
       write(texto,'(a)') 'uGal'                                        !
       call moveto(ix-32,iy+297,xy)                                     !
       call outgtext(texto)                                             !
      endif                                                             ! 
      do 26 i=1,21                                                      !
      k=gmi+(i-1)*pg-dg                                                 !
      nada=setcolorrgb(col(i))                                          !
      u1=ix+5                                                           !
      v1=iy +i*12                                                       !
      u2=ix+18                                                          !
      v2=iy +12+i*12                                                    !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      nada=setcolorrgb(#000000)                                         !
      nada=ellipse($gborder      ,u1,v1,u2,v2)                          !
      write(texto,'(i6)') k                                             !
      call moveto(ix-50,iy+i*12-2,xy)                                   !
      call outgtext(texto)                                              !
   26 continue                                                          !
      write(texto,'(a)') 'Obs.'                                         !
      call moveto(ix-35,iy-20,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Res.'                                         !
      call moveto(ix+40,iy2-20,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'uGal'                                         !
      call moveto(ix-35,iy-7,xy)                                        !
      call outgtext(texto)                                              !
      call moveto(ix+40,iy2-7,xy)                                       !
      call outgtext(texto)                                              !
c-----------------------Map.bln line------------------------------------!
      call cont(mb,nb,xb,yb,ib,jx2,jy0,fdib)                            !
      call cont(mb,nb,xb,yb,ib,jx1,jy0,fdib)                            !
c-----------------------Observed values---------------------------------!
      l=12                                                              !
      if(ns.gt.100) l=4                                                 !
      do 28 i=1,ns                                                      !
      ix=jx2 + x(i)/fe/fdib                                             !
      iy=jy0 - y(i)/fe/fdib                                             !
      k=(g(i)-gmi)/pg+1                                                 !
      nada=setcolorrgb(col(k))                                          !
      u1=ix-l                                                           !
      v1=iy-l                                                           !
      u2=ix+l                                                           !
      v2=iy+l                                                           !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      nada=setcolorrgb(#104010)                                         !
      nada=ellipse($gborder,u1,v1,u2,v2)                                !
   28 continue                                                          !
c---------------------------Observed-Modelled, residuals----------------!
      write(texto,'(a)') 'Observed     Modelled'                        !
      ly=110                                                            !
      call moveto(jx2-53,420,xy)                                        !
      call outgtext(texto)                                              !
      call ventana(jx2-lx,jy4-ly,jx2+lx,jy4+ly)                         !
      lyr=64                                                            !
      call ventana(jx2-lx,jy5-lyr,jx2+lx,jy5+lyr)                       !
      u1=jx2-lx+120                                                     !
      u2=u1+9                                                           !
      v1=jy4-125                                                        !
      v2=v1-5                                                           !
      nada=setcolorrgb(col(5))                                          !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      u1=u1+90                                                          !
      u2=u1+9                                                           !
      nada=setcolorrgb(col(16))                                         !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      u1=jx2-lx                                                         !
      u2=u1+2                                                           !
      v1=jy4                                                            !
      v2=v1-ly/3.                                                       !
      nada=setcolorrgb(#009000)                                         !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      v1=jy5                                                            !
      v2=v1-64/3.                                                       !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#000000)                                         !
      call setgtextrotation(900)                                        !
      write(texto,'(i5,a)') nint(rag/2/3.), ' uGal'                     !
      u1=jx2-lx-19                                                      !
      call moveto(u1,jy4+4,xy)                                          !
      call outgtext(texto)                                              !
      lyr=64                                                            !
      call ventana(jx2-lx,jy5-lyr,jx2+lx,jy5+lyr)                       !
      u1=jx2-lx                                                         !
      u2=u1+2                                                           !
      v1=jy5                                                            !
      v2=v1-64/3.                                                       !
      nada=setcolorrgb(#009000)                                         !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(i5,a)') nint(rag/2/3./3.), ' uGal'                  !
      u1=jx2-lx-19                                                      !
      call moveto(u1,jy5+10,xy)                                         !
      call outgtext(texto)                                              !

      call setgtextrotation(000)                                        !
                                                                        ! 
!-----------------------------------------------------------------------!
      
c    More initialization
      
      call step(ms,ns,x,y,rr,dr)
      do 15 i=1,ns
      w(i)=1.
   15 gc(i)=0.
      do 18 j=1,nc
      df(j)=0.
   18 lev(j)=0
      jend=0
      nlev=0
      nu=0
      ncf=0
      sm=0.
      fit0=9.d29
      ncc=nc/nra
      dgs=0
      pxs=0
      pys=0
      ddg=9.d5
      dg0=9.d5
      flev=1.1-ncer/10.

c-----------------------------------------------------------------------!
c                                                                       !     
c       Step by step growth of the anom anomalous structure             !
c                                                                       !
c-----------------------------------------------------------------------!
  
   30 continue
      c=1.00005 
      if(c.lt.3000) c=1.0001
      if(c.lt.2000) c=1.001
      if(c.lt.1000) c=1.01
      if(c.lt.500) c=1.1
      if(c.lt.100) c=2. 
      fit=c*fit0

c    Smoothing lambda value 

      cr=1.-0.5*exp(-nu*100./nc)    
      fl=fbal *cr 

c   Normal equations

      s11=0.
      s1x=0.
      s1y=0.
      s1g=0.
      sxx=0.
      sxy=0.
      syy=0.
      sxg=0.
      syg=0.      
      do 31 i=1,ns
      p2=w(i)
      xd=xk(i)
      yd=yk(i)
      dg=g(i)
      s11=s11+p2
      s1g=s1g+dg*p2
        s1x=s1x+xd*p2
        s1y=s1y+yd*p2
        sxx=sxx+xd*xd*p2
        sxy=sxy+xd*yd*p2
        sxg=sxg+xd*dg*p2
        syy=syy+yd*yd*p2
        syg=syg+yd*dg*p2
   31 continue

c  Exploring model space

      js=0
      do 35 j1=1,ncc
      j=j1
      if(ncc.lt.nc) then
        call RANDOM_NUMBER(ranval)
        j=(nc-1)*ranval+1
      endif

c  Number of revisiting for the selected cell

      klm=lev(j)                      
      kla=lay(j)
      if(klm.eq.mlev) go to 35        
        
      q0=qm(j) * ( 1.-zc(j)/fe*kup*1.e-4)

      tt=6.672e-3*vol(j) 
      do 34 k=1,2          
      if(klm.gt.1.and.klev(k,klm).lt.klev(k,1)*0.1) go to 35
      if(ncer.le.0) then
        if(klm.ge.1.and.klev(k,klm+1).gt.flev*klev(k,klm)) go to 35
      endif

      db=2*k-3                               
      fq=1
      if(ncer.gt.0.and.kla.gt.0) then
        ic=0 
        if(k.eq.1) then
          do 16 i=1,ncer     
          ii=iup(i,j)
          if(lay(ii).lt.kla) ic=ic+1
   16     continue
        endif
        if(k.eq.2) then
          do 17 i=1,ncer    
          ii=ido(i,j)
          if(lay(ii).gt.kla) ic=ic+1
   17   continue
        endif
        if(ic.lt.1) go to 34
        fq=(ncer*1./ic)**3            ! <<<<<<< 3  volumetrico
      endif

      qq=q0*fq
      srr=0.
      srg=0.
      s1r=0.
         px=pxs
         py=pys
         dg=dgs
         srx=0.
         sry=0.
         srz=0
      c=1.e4 *fe*fe
      do 32 i=1,ns
      xd=x(i)-xc(j)
      yd=y(i)-yc(j)
      zd=z(i)-zc(j)                 
      d=xd*xd+yd*yd+zd*zd + c 
      av=sqrt(zd*zd/d/d/d)*tt     
      rr(i)=gc(i)+av*db       
      zb=rr(i)*w(i)
      s1r=s1r+zb
      srg=srg+zb*g(i)
      if(igxy.ne.1) go to 32
        srx=srx+zb*xk(i)
        sry=sry+zb*yk(i)
   32 srr=srr+zb*rr(i)

      sv=sm+db*db*qq 

      su=srr+fl*sv  

         f11=s11*su-s1r*s1r
         f1x=s1x*su-s1r*srx
         f1y=s1y*su-s1r*sry
         f1g=s1g*su-s1r*srg
         if(igxy.eq.1) then
           fxx=sxx*su-srx*srx
           fxy=sxy*su-srx*sry
           fxg=sxg*su-srx*srg
           fyy=syy*su-sry*sry
           fyg=syg*su-sry*srg
           gxx=f11*fxx-f1x*f1x
           gxy=f11*fxy-f1x*f1y
           gxg=f11*fxg-f1x*f1g
           gyy=f11*fyy-f1y*f1y
           gyg=f11*fyg-f1y*f1g
           ryy=gyy*gxx-gxy*gxy
           ryg=gyg*gxx-gxy*gxg
           py=ryg/ryy
           px=(gxg-py*gxy)/gxx
         endif 

c   Solutions: offset and scale factor

         if(ig0.eq.1) dg=(f1g-px*f1x-py*f1y)/f11
         f=(srg-dg*s1r-px*srx-py*sry)/su
         if(f.lt.0.00001) go to 34 

c   Evaluating the solution: misfit value

      sv=f*f*sv*fl
      em=0
      do 33 i=1,ns
      c=g(i)-dg-px*xk(i)-py*yk(i)-rr(i)*f     
   33 em=em+c*c*w(i)                            
      c=sv+em  

c   Keeping best values

      if(c.gt.fit) go to 34
        fit=c
        fs=f
        js=j        
        ks=k
        dbs=db 
        qs=qq  
        ems=em
        svs=sv 
        dgs=dg
        pxs=px
        pys=py
   34 continue
   35 continue      
               
c   Growth does stop
      
      iend=1 
      if(js.eq.0) go to 50
      iend=2

c   Control of offset

c       c=abs(dgs-dg0)*0.02+ddg*0.98      
c       if(c.lt.ddg) ddg=c
c       if(ddg.lt.0.1) ig0=0
       if(nu.gt.6000) ig0=0

c   Levels of density. Average density

      lev(js)=lev(js)+1
      iz=lay(js)
      if(iz.gt.0.and.lev(js).gt.nlev) nlev=lev(js)
      dme=0
      k=0
      klm=0
      do 36 i=1,nlev
      if(lev(js).eq.i.and.iz.gt.0) then
        if(i.gt.ml) write(*,*) 'superado ml'
        klev(ks,i)=klev(ks,i)+1    
        if(i.gt.1) klev(ks,i-1)=klev(ks,i-1)-1
      endif
      if(klev(ks,i).eq.0) go to 36
        l=klev(1,i)+klev(2,i)
        if(l.gt.klm) klm=l
        k=k+klev(ks,i)
        dme=dme+klev(ks,i)*i
   36 continue
      dme=dme/k*fs                    
      c=fs*mlev/2.
      if(dme.lt.c) dme=c
     
c   Growth stop

      iend=2
      if(df(js).eq.0) ncf=ncf+1
      fin=ncf*100./nc
      if(nu.gt.6.and.tpc.le.fin) go to 50
      iend=0

c   Selected values    

      nu=nu+1
      if(fit.lt.fit0) fit0=fit
      df(js)=df(js)+dbs                           
      sm=sm+dbs*dbs*qs     
      do 37 i=1,ns
      xd=x(i)-xc(js)
      yd=y(i)-yc(js)
      zd=z(i)-zc(js)                       
      d=xd*xd+yd*yd+zd*zd + 1.e4*fe*fe  
      av=zd/sqrt(d*d*d)*6.672e-3*vol(js)  
   37 gc(i)=gc(i)+av*dbs                           
      
      far(js)=fs
      if(iz.le.0) far(js)=-fs
      if(far(js).gt.99999.) far(js)=99999.

      if(ncer.gt.0.and.iz.gt.0) lay(js)=lay(js)+ 2*ks-3 

      fsa=fs

c     Modelled and residual values 
      
      em=0.
      c=0
      rma=0
      do 38 i=1,ns
      tt=g(i)-dgs-xk(i)*pxs-yk(i)*pys
      re(i)= tt-gc(i)*fs
      rr(i)=abs(re(i))
      if(rr(i).gt.rma) rma=rr(i)   
      c=c  +w(i)
   38 em=em+w(i)*rr(i)*rr(i) 
      if(em.lt.0) em=0
      em=sqrt(em/c)
      rma=rma/10

c    New weigthing according previous residual values 

      bls=blun * sqrt(tpc/fin)   
      call rewe(ms,ns,re,bls,w)

c-----------------------Graphics----------------------------------------!
                                                                        !
      write(texto,'(i7,i6,f8.1,i6,f7.1,i7,i6,i5)')                      !
     & nu,ncf,dme*(2*ks-3),nint(em),fin,nint(dgs),                      !
     & nint(pxs),nint(pys)                                              !
      ly=204                                                            !
      nada=setcolorrgb(#ffffff)                                         !
      v1=ly                                                             !
      v2=ly+17                                                          !
      nada=rectangle($gfillinterior,18,v1,380,v2)                       !
      nada=setcolorrgb(#a0ffff)                                         !
      nada=rectangle($gfillinterior,210,v1,252,v2)                      !
      nada=setcolorrgb(#000000)                                         !
      call moveto(12,ly,xy)                                             !
      call outgtext(texto)                                              !
c------------------filled cell------------------------------------------!
      lx=190                                                            !
      ly=190                                                            !
      ix=xc(js)/fe/fdib                                                 !
      j=lx-10                                                           !
      if(ix.lt.-j.or.ix.gt.j) go to 30                                  !
      iy=yc(js)/fe/fdib                                                 !
      j=ly-10                                                           !
      if(iy.lt.-j.or.iy.gt.j) go to 30                                  !
      jx=jx1                                                            !
      jy=jy0                                                            !
      ix=ix+jx                                                          !
      iy=jy-iy                                                          !
      c=mlev                                                            !
      if(mlev.gt.12) c=12                                               !
      c=lev(js)/c                                                       !
      if(ks.eq.1) then                                                  !
        if(c.le.0.00) nada=setcolorrgb(#ffb070)                         !   
        if(c.ge.0.00) nada=setcolorrgb(#f57070)                         !
        if(c.ge.0.25) nada=setcolorrgb(#f06060)                         !
        if(c.ge.0.50) nada=setcolorrgb(#c05050)                         !
        if(c.ge.0.75) nada=setcolorrgb(#804040)                         !
      endif                                                             !
      if(ks.eq.2) then                                                  !
        if(c.le.0.00) nada=setcolorrgb(#70b0ff)                         !   
        if(c.ge.0.00) nada=setcolorrgb(#7070f7)                         !
        if(c.ge.0.25) nada=setcolorrgb(#6060f0)                         !
        if(c.ge.0.50) nada=setcolorrgb(#5050c0)                         !
        if(c.ge.0.75) nada=setcolorrgb(#404080)                         !
      endif                                                             !
      i=nint(tx(js)/2./fe/fdib)                                         !
      u1=ix-i                                                           !
      if(u1.lt.122) u1=122                                              !
      u2=ix+i                                                           !
      if(u2.gt.1028) u2=528                                             !
      j=nint(ty(js)/2./fe/fdib)                                         !
      v1=iy-j                                                           !
      if(v1.lt.22) v1=22                                                !
      v2=iy+j                                                           !
      if(v2.gt.418) v2=418                                              !
      if(lay(js).ne.0) nada=rectangle($gfillinterior,u1,v1,u2,v2)       !
      if(lay(js).eq.0) nada=rectangle($gborder,u1,v1,u2,v2)             !
      iz=-zc(js)/fe/fdib+jz+30                                          !
      k=nint(tz(js)/2./fe/fdib)                                         !
      v1=iz-k                                                           !
      v2=iz+k                                                           !
      if(lay(js).ne.0) nada=rectangle($gfillinterior,u1,v1,u2,v2)       !
      if(lay(js).eq.0) nada=rectangle($gborder,u1,v1,u2,v2)             !
      jy=jx0                                                            !
      iy=jy+yc(js)/fe/fdib                                              !
      u1=iy-i                                                           !
      u2=iy+i                                                           !
      if(lay(js).ne.0) nada=rectangle($gfillinterior,u1,v1,u2,v2)       !
      if(lay(js).eq.0) nada=rectangle($gborder,u1,v1,u2,v2)             !
c-------------------------Progres---------------------------------------!
      ix=jx3+210                                                        !
      iy=544                                                            !
      lx=20                                                             !
      ly=280                                                            !
      call ventana(ix+10,iy-30,ix+lx,iy+ly)                             !
      nada=setcolorrgb(#00c0f0)                                         !
      c=fin/tpc                                                         !
      j=c*(ly+32)                                                       !
      u1=ix+10                                                          !
      v1=iy+ly-j                                                        !
      u2=ix+lx                                                          !
      v2=iy+ly                                                          !
      if(v1.gt.(iy+ly)) v1=iy+ly                                        !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
c---------------------complete/partial graphic--------------------------!
      c=nu/30.                                                          !
      if(nu.gt.200) c=nu/60.                                            !
      if(nu.gt.400) c=nu/120.                                           !
      d=c-int(c)                                                        !
      if(d.ne.0) go to 30                                               !
c--------------------------scale residuals------------------------------!
      nada=setcolorrgb(#f0f0f0)                                         !
      u1=ix2+30                                                         !
      v1=iy2+12                                                         !
      u2=u1+45                                                          !
      v2=iy2+266                                                        !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#000000)                                         !
      do 41 i=1,21                                                      !
      k=(i-11)*rma                                                      !
      write(texto,'(i6)') k                                             !
      call moveto(ix2+30,iy2+i*12,xy)                                   !
      call outgtext(texto)                                              !
   41 continue                                                          !
c----------------------------offset-------------------------------------!   
      if(ig0.eq.1) then                                                 !
       u1=nu/12.+80                                                     !
       v1=920+(dgs-dg0)*30./kdg                                         !
       call moveto(up,vp,xy)                                            !
       nada=setcolorrgb(#0000b0)                                        !
       if(abs(vp-920).lt.300) nada=lineto(u1,v1)                        !
       up=u1                                                            !
       vp=v1                                                            !
      endif                                                             !
c----------------------------misfit-------------------------------------!
      if(nu.le.65000) then                                              !
       u1=nu/120.+680                                                   !
       j=kdg*3                                                          !
       if(ns.le.30) j=kdg*6                                             !
       v1=920- em*30./j                                                 !
       call moveto(up2,vp2,xy)                                          !
       nada=setcolorrgb(#b00000)                                        !
       if(abs(vp2-920).lt.300) nada=lineto(u1,v1)                       !
       up2=u1                                                           !
       vp2=v1                                                           !
      endif                                                             !
c----------------------- Residual and Model-----------------------------!
      lx=190                                                            !
      ly=190                                                            !
      jy=jz+ly                                                          !
      call hor0(jx3,jy0,lx,lx,tic,fe,fdib,mb,nb,xb,yb,ib)               !
      call hor0(jx3, jy,lx,lx,tic,fe,fdib,mb,nb,xb,yb,ib)               !
      do 42 i=1,ns                                                      !
      ix=jx3  + x(i)/fe/fdib                                            !
      iy= jy  - y(i)/fe/fdib                                            !
      k=11+re(i)/rma                                                    !
      if(k.gt.21) k=21                                                  !
      if(k.lt.1) k=1                                                    !
      nada=setcolorrgb(col(k))                                          !
      l=w(i)*14                                                         !
      if(l.lt.2) l=2                                                    !
      if(ns.gt.100) l=4                                                 !
      u1=ix-l                                                           !
      v1=iy-l                                                           !
      u2=ix+l                                                           !
      v2=iy+l                                                           !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      nada=setcolorrgb(#104010)                                         !
      nada=ellipse($gborder,u1,v1,u2,v2)                                !
      ix=jx3  + x(i)/fe/fdib                                            !
      iy=jy0  - y(i)/fe/fdib                                            !
      d=gc(i)*fs+xk(i)*pxs+yk(i)*pys+dgs                                !
      k=(d-gmi)/pg+1                                                    !
      if(k.gt.21) k=21                                                  !
      if(k.lt.1) k=1                                                    !
      nada=setcolorrgb(col(k))                                          !
      u1=ix-l                                                           !
      v1=iy-l                                                           !
      u2=ix+l                                                           !
      v2=iy+l                                                           !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      nada=setcolorrgb(#104010)                                         !
      nada=ellipse($gborder,u1,v1,u2,v2)                                !
c--------------------------Regional and local---------------------------!
      l=w(i)*12                                                         !
      if(l.lt.2) l=2                                                    !
      if(ns.gt.100) l=3                                                 !
      ix=jx4  + x(i)/fe/fdib/1.7                                        !
      iy=jy6  - y(i)/fe/fdib/1.7                                        !
      d=xk(i)*pxs+yk(i)*pys+dgs                                         !
      k=(d-gmi)/pg+1                                                    ! 
      if(k.gt.21) k=21                                                  !
      if(k.lt.1) k=1                                                    !
      nada=setcolorrgb(col(k))                                          !
      u1=ix-l                                                           !
      v1=iy-l                                                           !
      u2=ix+l                                                           !
      v2=iy+l                                                           !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      nada=setcolorrgb(#104010)                                         !
      nada=ellipse($gborder,u1,v1,u2,v2)                                !
      ix=jx4  + x(i)/fe/fdib/1.7                                        !
      iy=jy7  - y(i)/fe/fdib/1.7                                        !
      d=gc(i)*fs                                                        !
      k=11+d/pg                                                         !
      if(k.gt.21) k=21                                                  !
      if(k.lt.1) k=1                                                    !
      nada=setcolorrgb(col(k))                                          !
      u1=ix-l                                                           !
      v1=iy-l                                                           !
      u2=ix+l                                                           !
      v2=iy+l                                                           !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
      nada=setcolorrgb(#104010)                                         !
      nada=ellipse($gborder,u1,v1,u2,v2)                                !
   42 continue                                                          !
c---------------------Observed-modelled and residuals--      -----------!
      lx=190                                                            !
      ly=110                                                            !
      u1=jx2-lx+10                                                      !
      u2=jx2+lx                                                         !
      v1= jy4-ly                                                        !
      v2= jy4+ly                                                        !
      nada=setcolorrgb(#ffffff)                                         !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      lyr=64                                                            !
      v1= jy5-lyr                                                       !
      v2= jy5+lyr                                                       !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#808080)                                         !
      c=1.96*ly/rag                                                     !
      if(deg.lt.0) then                                                 !
         iy= jy4+ly-(dgs-gmi)*c                                         !
         u1=jx2-lx+15                                                   !
         v1=iy                                                          !
         call moveto(jx2+lx,iy,xy)                                      !
         nada=lineto(u1,v1)                                             !
      endif                                                             !
      do 44 i=1,ns                                                      !
      tt=gc(i)*fs+xk(i)*pxs+yk(i)*pys+dgs                               !
      ix=jx2 +x(i)/fe/fdib                                              !
      qq=ly-(g(i)-gmi)*c                                                !
      if(abs(qq).le.ly) then                                            !
       iy= jy4+qq                                                       !
       u1=ix-4                                                          !
       v1=iy-2                                                          !
       u2=ix+4                                                          !
       v2=iy+2                                                          !
       nada=setcolorrgb(col(5))                                         !
       nada=ellipse($gfillinterior,u1,v1,u2,v2)                         !
      endif                                                             !
      qq=ly-(tt-gmi)*c                                                  !
      if(abs(qq).le.ly) then                                            !
       iy= jy4+ly-(tt-gmi)*c                                            !
       u1=ix-4                                                          !
       v1=iy-2                                                          !
       u2=ix+4                                                          !
       v2=iy+2                                                          !
       nada=setcolorrgb(col(16))                                        !
       nada=ellipse($gfillinterior,u1,v1,u2,v2)                         !
      endif                                                             !

      qq=(g(i)-tt)*c*3                                                  !
      if(abs(qq).le.lyr) then                                           !
       u1=ix-1                                                          !
       v1=jy5                                                           !
       u2=ix+1                                                          !
       v2=jy5-qq                                                        !
       nada=setcolorrgb(#00b000)                                        !
       nada=rectangle($gfillinterior,u1,v1,u2,v2)                       !
      endif                                                             !
   44 continue                                                          !
      nada=setcolorrgb(#000000)                                         !
      nada=rectangle($gfillinterior,jx2-lx,jy5,jx2+lx,jy5)              !

c-------------------------Residual histog-------------------------------!
      iy=805                                                            !
      ix=258                                                            !
      lx=120                                                            !
      ly=138                                                            !
      call ventana(ix-20,iy-ly-8,ix+lx+17,iy+22)                        !
      call histog(ms,ns,w,re,kdib,ix,lx,iy,ly)                          !
c------------------------Density levels---------------------------------!
      ix=54                                                             !
      iy=805                                                            !
      lx=150/nlev                                                       !
      if(mlev.eq.1) lx=90                                               !    
      call ventana(ix-34,iy-146,ix+163,iy+23)                           !
      nada=setcolorrgb(#101040)                                         !
      do 45 i=1,nlev                                                    !
      v2=iy                                                             !
      if(klev(2,i).le.1) go to 43                                       !
      l=klev(2,i)*110./klm                                              !
      if(l.gt.110) l=110                                                !
      nada=setcolorrgb(#8080f0)                                         !
      u1=ix+2+(i-1)*lx                                                  !
      v1=iy                                                             !
      u2=ix+i*lx                                                        !
      v2=iy-l                                                           !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
   43 if(klev(1,i).le.1) go to 45                                       !
      l=klev(1,i)*110./klm                                              !
      if(l.gt.110) l=110                                                !
      nada=setcolorrgb(#f08080)                                         !
      u1=ix+2+(i-1)*lx                                                  !
      v1=v2                                                             !
      u2=ix+i*lx                                                        !
      v2=v2-l                                                           !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      write(texto,'(i4)') nint(i*fs)                                    !
      nada=setcolorrgb(#000000)                                         !
      k=ix-16+(i-0.5)*lx                                                !
      call moveto(k,iy+2,xy)                                            !
      if(nlev.le.6) call outgtext(texto)                                !
   45 continue                                                          !
      nada=setcolorrgb(#000000)                                         !
      call moveto(ix,iy,xy)                                             !
      u1=ix+165                                                         !
      v1=iy                                                             !
      nada=lineto(u1,v1)                                                !
      call moveto(ix,iy+5,xy)                                           !
      u1=ix                                                             !
      v1=iy-120                                                         !
      nada=lineto(u1,v1)                                                !
      write(texto,'(i5)') klm                                           !
      call moveto(ix-37,iy-115,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(i5)') nint(klm*2./3.)                               !
      call moveto(ix-37,iy-80,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(i4)') nint(klm/3.)                                  !
      call moveto(ix-30,iy-45,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(i1)') 0                                             !
      call moveto(ix-10,iy-8,xy)                                        !
      call outgtext(texto)                                              !
      write(texto,*) 'Num.Cells'                                        !
      call moveto(ix-25,iy-140,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Dens.con.'                                    !
      call moveto(ix+103,iy-113,xy)                                     !
      call outgtext(texto)                                              !
      write(texto,'(a,i3)') 'Ave=',nint(nlev*fs)                        !
      call moveto(ix+114,iy-93,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a,i3)') 'Max=',nint(dme)                            !
      call moveto(ix+114,iy-73,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a,f9.0)') 'kg/m3'                                   !
      call moveto(ix+127,iy-53,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a,i2)') 'nl=',nlev                                  !
      call moveto(ix+123,iy-138,xy)                                     !
      call outgtext(texto)                                              !
      if(kdib.eq.0) then                                                !
       write(texto,'(a)') 'Rep.cells on nl dens.levels'                 !
       call moveto(ix-28,iy-170,xy)                                     !
       call outgtext(texto)                                             !
      endif                                                             !
c----------------------------Autocorrel---------------------------------!
      call cov2(ms,ns,x,y,re,w,dr,mr,cov,swc)                           !
      ix=461                                                            !
      iy=770                                                            !
      lx=340                                                            !
      call ventana(ix-45,iy-110,ix+lx,iy+55)                            !
      nada=setcolorrgb(#000000)                                         !
      call moveto(ix,iy,xy)                                             !
      u1=ix+lx-2                                                        !
      v1=iy                                                             !
      nada=lineto(u1,v1)                                                !
      call moveto(ix,iy+40,xy)                                          !
      u1=ix                                                             !
      v1=iy-80                                                          !
      nada=lineto(u1,v1)                                                !
      if(kdib.eq.0) then                                                !
       write(texto,'(a,i5,a)')                                          ! 
     - 'Autocorrel. of Residues  (Step=',nint(dr),' m)'                 !
       call moveto(ix+30,iy-135,xy)                                     !
       call outgtext(texto)                                             !
      endif                                                             !
      write(texto,'(a,f5.2)') 'Autocorrel (1) =',cov(1)                 !
      call moveto(ix+100,iy-105,xy)                                     !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Horiz.'                                       !
      call moveto(ix+lx-42,iy-45,xy)                                    !
      call outgtext(texto)                                              !
      write(texto,'(a)') 'Dist.'                                        !
      call moveto(ix+lx-35,iy-25,xy)                                    !
      call outgtext(texto)                                              !
      do 46 i=1,10                                                      !
      u1=ix+i*tic/fe/fdib                                               !
      call moveto(u1,iy,xy)                                             !
      v1=iy+6                                                           !
      nada=lineto(u1,v1)                                                !
      if(i.gt.7) go to 46                                               !
      v1=iy-80*(i-3)*0.25                                               !
      u1=ix-5                                                           !
      call moveto(u1,v1,xy)                                             !
      u1=ix                                                             !
      nada=lineto(u1,v1)                                                !
      write(texto,'(f5.2)') -0.75+i*0.25                                !
      call moveto(ix-42,v1-10,xy)                                       !
      call outgtext(texto)                                              !
   46 continue                                                          !
      u1=ix+0.5*dr/fe/fdib                                              !
      v1=iy-cov(1)*80                                                   !
      nada=setcolorrgb(#009000)                                         !
      call moveto(u1,v1,xy)                                             !
      do 47 i=1,mr                                                      !
      u1=ix+(i-0.5)*dr/fe/fdib                                          !
      v1=iy-cov(i)*80                                                   !
      if(abs(cov(i)).eq.0) go to 47                                     !
      k=8*sqrt(swc(i))                                                  !
      if((v1-k).lt.(iy-110).or.(v1+k).gt.(iy+55)) go to 47              !
      if(u1.gt.(ix+lx)) go to 47                                        !
      if(k.eq.0) k=1                                                    !
      nada=lineto(u1,v1)                                                !
      u1=u1-k                                                           !
      v1=v1-k                                                           !
      u2=u1+2*k                                                         !
      v2=v1+2*k                                                         !
      nada=ellipse($gfillinterior,u1,v1,u2,v2)                          !
   47 continue                                                          !
      kdib=1                                                            !
!-----------------------------------------------------------------------

      dg0=dgs
      if(jend.eq.0) go to 30
      go to 60     

c-----------------------------------------------------------------------!
c           End of model growth                                         !
c-----------------------------------------------------------------------!

   50 continue

c-----------------------------------------------------------------------!
c           End of process                                              !
c-----------------------------------------------------------------------!

   60 fs=fsa
      fit0=sqrt(fit0) 

c Final residual and absolute residuals

      c=0
      do 51 i=1,ns
      re(i)=g(i)-gc(i)*fs
      c=c+re(i)
   51 gc(i)=abs(re(i))
      c=c/ns
      syg=c/0.6745


c  Writing modelled and residual values

      open(2,file=obs)
      open(3,file=fil)
      write(3,'(a,a)') '    X       Y       Obs       Off       ',
     & 'dg        Mod      Mod-dg    S-Res   Test'
      nms=0
      j1=9999
      j2=-j1
      disr=0
      disl=0
      disv=0
   56 read(2,*,end=58) xp,yp,zp,dg,eg
      if(xp.eq.0.and.yp.eq.0.and.zp.eq.0) go to 58
      xp=xp-ixm
      yp=yp-iym
      zp=zp-izm
      greg=dgs+xp/1000.*pxs+yp/1000.*pys
      gloc=dg-greg
      nms=nms+1
      gcal=0
      do 57 j=1,nc
      xd=xp-xc(j)
      yd=yp-yc(j)
      zd=zp-zc(j)                        
      d=xd*xd+yd*yd+zd*zd + 1.e4*fe*fe 
      av=zd/sqrt(d*d*d)*6.672e-3*vol(j)  
   57 gcal=gcal+av*df(j)*fs
      gres=gloc-gcal
      disl=disl+gloc*gloc
      disr=disr+greg*greg
      disv=disv+gres*gres
      gxg=abs(gres)/syg
      ts1=bl
      ts2=bl
      ts3=bl
      ts4=bl
      ts5=bl 
      ts6=bl
      if(gxg.gt.1.0*blun) ts1=as
      if(gxg.gt.1.5*blun) ts2=as
      if(gxg.gt.2.0*blun) ts3=as
      if(gxg.gt.2.5*blun) ts4=as      
      if(gxg.gt.3.0*blun) ts5=as      
      if(gxg.gt.4.0*blun) ts6=as            
      write(3,203) xp+ixm,yp+iym,dg,greg,gloc,gcal,gres,gxg,
     -  ts1,ts2,ts3,ts4,ts5,ts6
  203 format(f8.0,f9.0,2f9.1,3f10.1,f9.2,2x,6a1)
      if(gres.gt.j2) j2=gres
      if(gres.lt.j1) j1=gres
      go to 56
   58 close(3)
      close(2)
      disv=sqrt(disv/nms)
      disl=sqrt(disl/nms)
      disr=sqrt(disr/nms)


c   Writing resulting model in Mod.txt

      open(1,file=mod)
      nn=0
      zmin=0
      do 52 j=1,nc
      if(df(j).eq.0) go to 52
      nn=nn+1
      if(zc(j).lt.zmin) zmin=zc(j)
   52 continue 
      zmi=zmin -1
      zmin=zmin-500*fe
      if(ns.gt.30) nn=nc
      write(1,'(i6,5x,a,f6.0,3x,a,f6.1,3x,a,f4.1,
     & 3x,a,i4,3x,a,i2)') 
     & nn,'Cell side=',cell,'Balance=',bal,'%-model=',tpc,
     & 'Dep.val=',kup,'Flatt=',ncer
      write(1,*) '  X(m)   Y(m)   Z(m asl)   sx(m)  sy(m)   sz(m)',
     & '  den(kg/m3)   sen    ord'
      srx=0.
      sry=0.
      px=0.
      py=0.
      srg=0.
      tt=0
      zd=fs*(mlev-0.5)         !<<<<<
      if(mlev.eq.1) zd=fs
      do 54 j=1,nc
      df(j)=df(j)*fs             
      if(ns.lt.30.and.df(j).eq.0) go to 53
      dg=vol(j)*df(j)                                     ! mass
      c=abs(dg)                                           ! absol mass
      srg=srg+c                                           ! addition absol. mass
      if(df(j).lt.0) then
        srx=srx+c                                         ! addition neg. mass
        px=px+c*zc(j)                                     ! mean depth
      endif
      if(df(j).gt.0) then
        sry=sry+c                                         ! addition pos. mass
        py=py+c*zc(j)                                     ! mean depth
      endif
      if(ncer.gt.0.and.lay(j).gt.0) then
        if(zc(j).ge.zmi) df(j)=  lay(j)*fs   !*(1.-lay(j)/40.)
        if(df(j).gt.tt) tt=df(j)   
        if(zc(j).lt.zmi) df(j)=tt  
      endif 
      call ss2(ms,ns,x,y,z,fe,xc(j),yc(j),zc(j),sen,k)
      write(1,'(i7,2i8,2x,3i7,f9.1,2x,f9.1,f7.1 , 2i6)') 
     -nint(xc(j)+ixm),nint(yc(j)+iym),nint(zc(j)+izm),
     -nint(tx(j)),nint(ty(j)),nint(tz(j)),df(j),
     -sen,far(j)   ,lay(j),lev(j)
      if(lay(j).le.0) go to 53
        k=lay(j)
        swc(k)=zc(j)
        dla(k)=dla(k)+df(j)*tx(j)*ty(j)  
   53 continue
      zd=zc(j)
      if(zd.lt.zmin) zd=zmin
   54 continue 
      close(1)    
      zme=(px+py)/(srx+sry+1.)+izm
      if(srx.gt.0.) px=px/srx+izm                         ! weighted neg. depth
      if(sry.gt.0.) py=py/sry+izm                         ! weighted pos. depth   

c  Maximun density vertical variation for contiguos layers
    
      cmax=0
      do 55 i=1,nlay
      dla(i)=dla(i)/fe/fe*1.0e-6
      if(i.eq.1) go to 55
        c=dla(i-1)-dla(i)
        if(c.gt.cmax) cmax=c
   55 continue
     
c----------------------- Final drawing ----------------------------------
      iy=835                                                            !
      ix=jx3-205                                                        !
      call ventana(ix+130,iy+30,ix+330,iy+70)                           !
      call ventana(ix+130,iy+110,ix+365,iy+126)                         !
      nada=setcolorrgb(#000000)                                         !
      write(texto,*) 'Tot.      Neg.    Posit.'                         !
      call moveto(ix+150,iy+10,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a,3e10.2)') '    Anom.mass(kg)',srg,srx,sry         !
      call moveto(ix,iy+30,xy)                                          !
      call outgtext(texto)                                              !
      write(texto,'(a,3f10.0)') '    Aver.depth(m)',zme,px,py           !
      call moveto(ix,iy+50,xy)                                          !
      call outgtext(texto)                                              !
      write(texto,'(a)') '  Tot.   Off.  Loc.  Resid.  Re/Lo'           !
      call moveto(ix+135,iy+90,xy)                                      !
      call outgtext(texto)                                              !
      write(texto,'(a,f8.0,2f7.0,f6.0,f6.2)') '    Scatter.(uGal)',     !
     & disg,disr,disl,em, em/disv                                       !
      call moveto(ix,iy+110,xy)                                         !
      call outgtext(texto)                                              !
      ix=jx3+220                                                        !
      call ventana(ix,iy+60,ix+220,iy+80)                               !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)') 'Process completed'                            !
      call moveto(ix+45,iy+30,xy)                                       !
      call outgtext(texto)                                              !
      if(iend.eq.1) write(texto,'(a)') 'No further suitable cells'      !
      if(iend.eq.2) write(texto,'(a)') 'OK  End value (den,%) reached'  !
      call moveto(ix+5,iy+60,xy)                                        !
      call outgtext(texto)                                              !
      write(texto,'(a)') ' ===> Press (Enter)'                          !
      call moveto(ix+40,iy+100,xy)                                      !
      call outgtext(texto)                                              !

c-----------------------------------------------------------------------!
      iv=0                                                              !

   70 read(*,'(i1)') i
c------------------- Drawing several sections --------------------------!
      lx=190                                                            !
      ly=190                                                            !
      lz=162                                                            !
      seli=13.                                                          !
      zi=9.d9                                                           !
      zs=-zi                                                            !
      do 71 i=1,nc                                                      !
      if(lay(i).eq.0.or.lev(i).lt.1) go to 71                           !
      j=zc(i)/fe/fdib                                                   !
      if(j.lt.-lz) go to 71                                             !
      c=df(i)                                                           !
      if(c.lt.zi) zi=c                                                  !
      if(c.gt.zs) zs=c                                                  !
   71 continue                                                          !
      zi=zi-(zs-zi)/30.                                                 !
      pa=(zs-zi)/20.                                                    !
      call clearscreen($gclearscreen)                                   !
      nada=setcolorrgb(#000000)                                         !
      nada=initializefonts()                                            !
      write(texto,'(a)') 'GROWTH-23'                                    !
      nada=setfont('t''Times New Roman''h23w11')                        !
      call moveto(160,15,xy)                                            !
      call outgtext(texto)                                              !
      write(texto,'(a)')                                                !
     & 'Inversion of gravity changes as 3D mass sources'                !
      nada=setfont('t''Times New Roman''h19w9')                         !
      call moveto(40,42,xy)                                             !
      call outgtext(texto)                                              !
      nada=setfont('t''Times New Roman''h16w07')                        !
      write(texto,'(a)') hoy                                            !
      call moveto(200,70,xy)                                            !
      call outgtext(texto)                                              !
                                                                        !
      ix=360                                                            !
      iy=800                                                            !
      nada=setfont('t''Times New Roman''h14w06')                        !
      i=0                                                               !
      j=0                                                               !
      do 72 k=1,21                                                      !
      c=zi+(k-1)*pa                                                     !
      if(k.gt.11.and.j.eq.0) i=0                                        !
      if(k.gt.11.and.j.eq.0) j=1                                        !
      i=i+1                                                             !
      nada=setcolorrgb(col(k))                                          !
      u1=ix-3                                                           !
      if(k.gt.11) u1=ix+65                                              !
      v1=iy+i*14                                                        !
      u2=u1+18                                                          !
      v2=v1+13                                                          !
      nada=setcolorrgb(#d0f0d0)                                         !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(col(k))                                          !
      u1=u1+3                                                           !
      v1=v1+1                                                           !
      u2=u2-3                                                           !
      v2=v2-1                                                           !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(f6.0)') c                                           !
      call moveto(u1-37,v1,xy)                                          !
      call outgtext(texto)                                              !
   72 continue                                                          !
      nada=setfont('t''Times New Roman''h15w07')                        !
      write(texto,'(a)') 'kg/m3'                                        !
      call moveto(ix-95,iy+45,xy)                                       !
      call outgtext(texto)                                              !
      ix=145                                                            !
      iy=860                                                            !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(i6,a)') nint(tic),' m'                              !
      call moveto(ix,iy-2,xy)                                           !
      call outgtext(texto)                                              !
      v1=iy+28                                                          !
      call moveto(ix+3,v1,xy)                                           !
      u1=ix+47                                                          !
      nada=lineto(u1,v1)                                                !
      u1=ix+10                                                          !
      call moveto(u1,v1-5,xy)                                           !
      nada=lineto(u1,v1)                                                !
      u1=ix+10+tic/fe/fdib                                              !
      call moveto(u1,v1-5,xy)                                           !
      nada=lineto(u1,v1)                                                !
      u1=ix-6                                                           !
      call moveto(u1,iy-17,xy)                                          !
      v1=iy+24                                                          !
      nada=lineto(u1,v1)                                                !
      v1=iy+20                                                          !
      call moveto(u1+5,v1,xy)                                           !
      nada=lineto(u1,v1)                                                !
      v1=iy+20-tic/fe/fdib                                              !
      call moveto(u1+5,v1,xy)                                           !
      nada=lineto(u1,v1)                                                !
      nada=setfont('t''Times New Roman''h16w07')                        !
      jx=232                                                            !
      jy1=325                                                           !
      call ventana(jx-lx,jy1-ly,jx+lx,jy1+ly)                           !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx-lx                                                          !
      u2=jx+lx                                                          !
      v1=jy1-ly                                                         !
      v2=jy1+ly                                                         !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)') 'Plan view and profiles'                       !
      call moveto(jx-70,jy-ly-25,xy)                                    !
      call outgtext(texto)                                              !
      call hor0(jx,jy1,lx,ly,tic,fe,fdib,mb,nb,xb,yb,ib)                !
      do 75 ix=1,nc                                                     !
      i=nc-ix+1                                                         !
      if(df(i).eq.0.and.ns.lt.30) go to 75                              !
      if(lay(i).le.0) go to 75                                          !
      if(isen(i).gt.4.5*seli) go to 75                                  !
      if(isen(i).gt.       0) r=2.0                                     !
      if(isen(i).gt.1.5*seli) r=2.6                                     !
      if(isen(i).gt.2.0*seli) r=3.2                                     !
      if(isen(i).gt.2.5*seli) r=3.8                                     !
      if(isen(i).gt.3.0*seli) r=4.4                                     !
      if(isen(i).gt.3.5*seli) r=5.0                                     !
      if(isen(i).gt.4.0*seli) r=5.6                                     !
      j=(xc(i)+tx(i)/r)/fe/fdib                                         !
      if(j.gt.lx) go to 75                                              !
      u2=jx+j                                                           !
      j=(xc(i)-tx(i)/r)/fe/fdib                                         !
      if(j.lt.-lx) go to 75                                             !
      u1=jx+j                                                           !
      j=(yc(i)+ty(i)/r)/fe/fdib                                         !
      if(j.gt.ly) go to 75                                              !
      v2=jy1-j                                                          !
      j=(yc(i)-ty(i)/r)/fe/fdib                                         !
      if(j.lt.-ly) go to 75                                             !
      v1=jy1-j                                                          !
      k=(df(i)-zi)/pa +1                                                !
      if(k.gt.21) k=21                                                  !
      if(k.lt.1) k=1                                                    !
      nada=setcolorrgb(col(k))                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
   75 continue                                                          !
      call cont(mb,nb,xb,yb,ib,jx,jy1,fdib)                             !
      jx1=jx                                                            !
      nada=setfont('t''Times New Roman''h14w06')                        !
      nada=setcolorrgb(#000000)                                         !
      do 74 i=-10,10                                                    !
      k=i*tic/fe/fdib                                                   !
      if(abs(k).gt.lx) go to 74                                         !
      v1=jy1+k-6                                                        !
      write(texto,'(i7)') nint(iym-i*tic)                               !
      call moveto(jx+lx+6,v1,xy)                                        !
      nada=setfont('t''Times New Roman''h14w06')                        !
      call outgtext(texto)                                              !
      u1=jx+k-6                                                         !
      write(texto,'(i6)') nint((i-1)*tic+ixm)                           !
      call moveto(u1,jy1+ly+45,xy)                                      !
      call setgtextrotation(900)                                        !
      nada=setfont('t''Times New Roman''h12w05')                        !
      call outgtext(texto)                                              !
      call setgtextrotation(000)                                        !
   74 continue                                                          !
      nada=setfont('t''Times New Roman''h16w07')                        !
      jy2=620                                                           !
      call ventana(jx-lx,jy2,jx+lx,jy2+lz)                              !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx-lx                                                          !
      u2=jx+lx                                                          !
      v1=jy2+30                                                         !
      v2=jy2+lz                                                         !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)') 'WE elevation view and profiles'               !
      call moveto(jx1-100,jy2-25,xy)                                    !
      call outgtext(texto)                                              !
      do 76 i=1,nc                                                      !
      if(df(i).eq.0.and.ns.lt.30) go to 76                              !
      if(lay(i).le.0) go to 76                                          !
      j=zc(i)/fe/fdib                                                   !
      if(j.lt.-lz) go to 76                                             !
      c=ty(i)/2.                                                        !
      if(isen(i).gt.4.5*seli) go to 76                                  !
      if(isen(i).gt.       0) r=2.0                                     !
      if(isen(i).gt.1.5*seli) r=2.6                                     !
      if(isen(i).gt.2.0*seli) r=3.2                                     !
      if(isen(i).gt.2.5*seli) r=3.8                                     !
      if(isen(i).gt.3.0*seli) r=4.4                                     !
      if(isen(i).gt.3.5*seli) r=5.0                                     !
      if(isen(i).gt.4.0*seli) r=5.6                                     !
      j=(xc(i)+tx(i)/r)/fe/fdib                                         !
      if(j.gt.lx) go to 76                                              !
      u2=jx+j                                                           !
      j=(xc(i)-tx(i)/r)/fe/fdib                                         !
      if(j.lt.-lx) go to 76                                             !
      u1=jx+j                                                           !
      j=(zc(i)+tz(i)/r)/fe/fdib                                         !
      if(j.gt.lz) go to 76                                              !
      v2=jy2-j  +30                                                     !
      j=(zc(i)-tz(i)/r)/fe/fdib                                         !
      if(j.lt.-lz) go to 76                                             !
      v1=jy2-j  +30                                                     !
      k=(df(i)-zi)/pa +1                                                !
      if(k.gt.21) k=21                                                  !
      if(k.lt.1) k=1                                                    !
      nada=setcolorrgb(col(k))                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
   76 continue                                                          !
      call ver0(jx1,jy2,lx,ly,lz,izm,tic,fe,fdib)                       !
      nada=setfont('t''Times New Roman''h14w06')                        !
      nada=setcolorrgb(#000000)                                         !
      j=izm/tic+2                                                       !
      do 77 i=2,20                                                      !
      k=((i-j)*tic+izm)/fe/fdib                                         !
      v1=jy2+23+k                                                       !
      if(v1.gt.(jy2+lz)) go to 77                                       !
      write(texto,'(i6)') -nint((i-j)*tic)                              !
      call moveto(jx+lx+6,v1,xy)                                        !
      call outgtext(texto)                                              !
   77 continue                                                          !
      nada=setfont('t''Times New Roman''h16w07')                        !
      l=3000*fe                                                         !
      lk=0                                                              !
      do 79 i=1,6                                                       !
      k=-(i+1)*(i+1)*l/17.                                              !
      jy=250                                                            ! 
      if(i.ge.4) jy=700                                                 !
      if(i.eq.1.or.i.eq.4) jx=700                                       !
      if(i.eq.2.or.i.eq.5) jx=1110                                      !
      if(i.eq.3.or.i.eq.6) jx=1520                                      !
      lk=lk+1                                                           !
      call sect(jx,jy,mc,nc,mb,nb,xc,yc,zc,tx,ty,tz,df,isen,pa,zi,lay,  !
     &  xb,yb,ib,col,fe,iym,izm,fdib,tic,seli,jx1,jy2,1,k,le,lk)        !
      call cont(mb,nb,xb,yb,ib,jx,jy,fdib)                              !
   79 continue                                                          !

      write(texto,'(a)') ' ===> Press (Enter)'                          !
      call moveto(1200,950,xy)                                          !
      nada=setcolorrgb(#000000)                                         !
      call outgtext(texto)                                              !
      read(*,'(i1)') i                                                  !
      nada=setcolorrgb(#f0f0f0)                                         !
      nada=rectangle($gfillinterior,500,10,1860,970)                    !
      lk=0                                                              !
      do 78 i=1,15                                                      !
      if(i.ge. 1) jx=700                                                !
      if(i.ge. 6) jx=1110                                               !
      if(i.ge.11) jx=1520                                               !
      if(i.ge. 1) jy=50+(i- 1)*182                                      !  180
      if(i.ge. 6) jy=50+(i- 6)*182                                      !
      if(i.ge.11) jy=50+(i-11)*182                                      !
      k=(8-i)*l*0.8                                                     !
      lk=lk+1                                                           !
      call sect(jx,jy,mc,nc,mb,nb,xc,yc,zc,tx,ty,tz,df,isen,pa,zi,lay,  !
     &  xb,yb,ib,col,fe,iym,izm,fdib,tic,seli,jx1,jy1,2,k,le,lk)        !
   78 continue                                                          !
c-----------------------------------------------------------------------!
         go to 99

C  Some messages for error or warning

   61 nada=setcolorrgb(#000000)                                         !
      write(texto,'(a)')  '*** Error: Unexpected execution end !!!'     !
      call moveto(650,740,xy)                                           !
      call outgtext(texto)
      if(nu.le.1) write(texto,'(a)')
     &'Verify data suitability (units, magnitudes) and model hypothesis'
      call moveto(660,760,xy)
      if(nu.le.1) call outgtext(texto)
      go to 99

   62 write(texto,'(3a)') '*** Error reading file ',obs,' !!'
      call moveto(650,740,xy)
      call outgtext(texto)
      go to 99

C   End execution

   99 write(texto,'(a)') ' ===> Press (Enter)'                          !
      call moveto(1200,950,xy)                                          !
      nada=setcolorrgb(#000000)                                         !
      call outgtext(texto)                                              !
      read(*,'(i1)') i                                                  !

      stop
      end
      
c***********************************************************************
c     Autocorrelation of n<m data v(x,y) for a distance dr             c
c***********************************************************************
      subroutine cov2(m,n,x,y,v,w,dr,mr,cov,sw)
      dimension x(m),y(m),v(m),w(m),cov(mr),sw(mr)
        
      vm=0.
      suw=0.
      do 1 i=1,n
      suw=suw+w(i)
    1 vm=vm+v(i)*w(i)
      vm=vm/suw
      sd2=0
      do 2 i=1,n
      vi=v(i)-vm
    2 sd2=sd2+vi*vi*w(i)
      sd2=sd2/suw
      
      do 7 k=1,mr
      cov(k)=0.
    7 sw(k)=0. 
      dr2=dr*dr
     
      do 3 i=1,n-1
      do 4 j=i+1,n
      xx=x(i)-x(j)
      yy=y(i)-y(j)
      d=sqrt(xx*xx+yy*yy)
      do 5 k=1,mr
      if(d.ge.dr*(k-1).and.d.le.dr*k) then
        cov(k)=cov(k)+(v(j)-vm)*(v(i)-vm)*w(j)*w(i)
        sw(k)=sw(k)   +w(j)*w(i)
      endif
    5 continue  
    4 continue
    3 continue
      vm=0
      do 6 k=1,mr
      if(sw(k).gt.vm) vm=sw(k)    
    6 if(sw(k).ne.0) cov(k)=cov(k)/sw(k)/sd2
      do 8 k=1,mr
    8 sw(k)=sw(k)/vm
      
      return
      end
      
c***********************************************************************
c     Execution stops when number of data is higher than dimensions    *
c***********************************************************************
      subroutine dim(m)
      write(*,200) m
  200 format(/6x,'*** Error: data size > max. dimension',i6,' !!!'/
     -10x,'Reduce the data side or change dimension in the code'/)
      stop
      end
      
c***********************************************************************
c     Sensitivity calculus
c***********************************************************************
      subroutine at2(ms,ns,x,y,z,xx,yy,zz,dx,dy,dz,sen,av,zp)
      dimension x(ms),y(ms),z(ms)
      av=0.
      t=(dx+dy+dz)/3.  
      t2=t*t
      zp=0
      sp=0 
      do 1 i=1,ns
      zr=zz-z(i)
      xr=xx-x(i)
      yr=yy-y(i)
      d=xr*xr+yr*yr+zr*zr + 1.e4
      av=av+zr*zr/d/d/d   
      d=1./d/d
      sp=sp+d
    1 zp=zp+z(i)*d
      av=av/ns
      sen=3.d-7/sqrt(av)
      t2=dx*dy*dz
      av=av*t2*t2  
      zp=zp/sp
      return
      end

c***********************************************************************
c     New cells for a layer
c           zs: top of the layer    tz: width of the layer 
c***********************************************************************
      subroutine celdascapa(ax,bx,ay,by,zs,tz,avm,ms,ns,x,y,z,
     -mc,nc,xc,yc,zc,dx,dy,dz,at,seli,ncr,lay,nlay)
      dimension xc(mc),yc(mc),zc(mc),dx(mc),dy(mc),dz(mc),
     -x(ms),y(ms),z(ms),lay(mc)
      character*50 texto
      real*4 at(mc)
      zi=zs-tz                           
      zz=zs-tz/2.
      tx=tz                             
      ty1=tz/3                          
      ty2=tz*3                          
      tt=tz*tz 
      kk=nc
      
      do 7 i=1,999                       
      xx=ax+tx*(i-1)
      if(xx.gt.bx) go to 9
      ty=tz
      yr=ay
      do 6 j=1,999                        
      if(yr.ge.by) go to 7
      c=ty
      l=0
    3 yy=yr+c/2
      call at2(ms,ns,x,y,z,xx,yy,zz,tx,c,tz,sen,am,zp)
      if(sen.gt.seli*10.) go to 6             
      if((zp-100).lt.zs.or.am.le.0) go to 6   
      av=am/avm
      if(abs(av-1.).le.0.1) go to 4        
      c=c/av**0.41
      if(c.gt.ty2) go to 6  
      if(c.lt.ty1) go to 6           
      l=l+1
      if(l.le.10) go to 3
      go to 6
    4 if(c.gt.9999) go to 6
    5 ty=c                              

c   Eliminacion celdas superpuestas

      if(ncr.eq.0) go to 12
      do 11 k=1,ncr
      if(abs(xx-xc(k)).gt.(dx(k)+tx)/2.6) go to 11  
      if(abs(yy-yc(k)).gt.(dy(k)+ty)/2.6) go to 11  
      if(abs(zz-zc(k)).gt.(dz(k)+tz)/2.6) go to 11  
      go to 6
   11 continue

c  New cell  

   12 kk=kk+1
      xc(kk)=xx
      yc(kk)=yy
      zc(kk)=zz
      dx(kk)=tx
      dy(kk)=ty
      dz(kk)=tz
      at(kk)=av
      lay(kk)=nlay
      if(kk.lt.mc) go to 6 
         kk=-9999
         go to 9
    6 yr=yr+ty
    7 continue
      
    9 nc=kk
      return
      end
      
c***********************************************************************
c                 Width of a layer                                     *
c***********************************************************************
      subroutine especapa(zs,avm,ms,ns,x,y,z,tz)
      dimension x(ms),y(ms),z(ms)
      l=0
      c=tz                           
    1 zz=zs-c*0.5                    
      k=0
      ama=-9.d9
      ame=0
      do 2 i=1,ns
      call at2(ms,ns,x,y,z,x(i)*1.5,y(i)*1.5,zz,c,c,c,se,aa,zp)
      if((zp-zz).le.100) go to 2
      if(aa.gt.ama) ama=aa
      ame=ame+aa
      k=k+1
    2 continue
      if(k.le.1) go to 9
      am=ame/k/avm
      c=c/am**0.13
      l=l+1
      if(abs(am-1.).gt.0.1.and.l.lt.12) go to 1
      
      if(c.gt.tz) tz=c
    9 return
      end
      
c***********************************************************************
c            Step for autocorrelation                                  *
c***********************************************************************
      subroutine step(m,n,x,y,d,dr)
      dimension x(m),y(m),d(m)
      dr=0
      do 1 i=1,n
      d4=9.e9
      d3=9.e9
      d2=9.e9
      de=9.e9
      do 2 j=1,n
      xx=x(i)-x(j)
      yy=y(i)-y(j)
      dd= xx*xx+yy*yy
      if(dd.gt.d4.or.dd.le.1.) go to 2
      if(dd.lt.de) de=dd
      if(dd.lt.d2.and.dd.gt.de) d2=dd
      if(dd.lt.d3.and.dd.gt.d2) d3=dd
      if(dd.lt.d4.and.dd.gt.d3) d4=dd
    2 continue
      de=sqrt(de)
      d2=sqrt(d2)
      d3=sqrt(d3)
      d4=sqrt(d4)
      d(i)=de*0.4+d2*0.3+d3*0.2+d4*0.1  
    1 dr=dr+d(i)
      dr=dr/n
      return
      end
      
c********************************************************************
c         Square window  in drawing                                 *
c********************************************************************
      subroutine ventana(i1,j1,i2,j2)
      use IFQWIN
      integer*2 x1,x2,y1,y2
      nada=setcolorrgb(#909090)
      x1=i1-2
      y1=j1-2
      x2=i1-1
      y2=j2+1
      nada=rectangle($gfillinterior,x1,y1,x2,y2)
      x2=i2+1
      y2=j1-1
      nada=rectangle($gfillinterior,x1,y1,x2,y2)
      nada=setcolorrgb(#fefefe)
      x1=i1-1
      y1=j2+2 
      x2=i2+2
      y2=j2+2
      nada=rectangle($gfillinterior,x1,y1,x2,y2)
      y1=j1
      nada=rectangle($gfillinterior,x2,y1,x2,y2)
      nada=setcolorrgb(#ffffff)
      x1=i1
      x2=i2
      y2=j2
      nada=rectangle($gfillinterior,x1,y1,x2,y2)
      return
      end
      
c**********************************************************************
c       Rounding numbers                                              *
c**********************************************************************    
      subroutine round(c,r)
      j=1
      if(c.lt.0) j=-1
      r=abs(c)
      i=log10(r)                                               
      r=j*nint(r/10**i)*10**i  
      return 
      end
      
c***********************************************************************
c     Dialog interface 1 (data files)                                  *
c***********************************************************************
      subroutine dial1(obs,mod,fil,na,bal,tpc,kup,ncer,mlev)
      use dflogm
      include 'Growth.fd'
      logical retlog,nada
      character*50 texto,obs,mod,fil
      external fin
      type(dialog) dlg
      retlog=DlgInit( IDD_Growth, dlg)

      retlog = DlgSet(dlg,IDC_EDIT_e1 ,"Gra.txt" )
      retlog = DlgSet(dlg,IDC_EDIT_e3 ,"Mod.txt" )
      retlog = DlgSet(dlg,IDC_EDIT_e4 ,"Fil.txt" )
      write(texto,'(i3)') na            
      retlog = DlgSet(dlg,IDC_EDIT_e5,texto)
      write(texto,'(i3)') nint(bal)      
      retlog = DlgSet(dlg,IDC_EDIT_e6,texto)
      write(texto,'(f3.0)') tpc 
      retlog = DlgSet(dlg,IDC_EDIT_e7,texto)
      retlog = DlgSet(dlg,IDC_CHECK_c3,.false.,dlg_enable)
      write(texto,'(i2)') kup    
      retlog = DlgSet(dlg,IDC_EDIT_e11,texto) 
      write(texto,'(i2)') ncer    
      retlog = DlgSet(dlg,IDC_EDIT_e13,texto) 
      retlog = DlgSet(dlg,IDC_CHECK_c3,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_CHECK_c2,.false.)
      retlog = DlgSet(dlg,IDC_CHECK_c1,.false.)   
      retlog = DlgSet(dlg,IDC_CHECK_c3,.false.)   
  
      retlog=DlgSet(dlg,IDC_TEXT_t2 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t5 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t6 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t7 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t8 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t9,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t11,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t12,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t13,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t14,.false.,dlg_enable)

      retlog=DlgSet(dlg,IDC_EDIT_e2 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e5 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e6 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e7 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e11,.false.,dlg_enable) 
      retlog=DlgSet(dlg,IDC_EDIT_e13,.false.,dlg_enable) 

      retlog=DlgSet(dlg,IDC_BOX_b3  ,.false.,dlg_enable) 
      retlog=DlgSet(dlg,IDC_BOX_b4  ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_CHECK_c2,.false.,dlg_enable)      
      retlog=DlgSet(dlg,IDC_CHECK_c1,.false.,dlg_enable)
      retlog = DlgSetSub(dlg,IDCANCEL,fin)

      retint = DlgModal( dlg )

      nada=DlgGet(dlg,IDC_EDIT_e1, obs )
      nada=DlgGet(dlg,IDC_EDIT_e3, mod )
      nada=DlgGet(dlg,IDC_EDIT_e4, fil )

      return
      end      
           
c***********************************************************************
c     Dialog interface 2  (average cell size)                          *
c***********************************************************************
      subroutine dial2(obs,mod,fil,na,ns,bal,tpc,kup,ncer,mlev,pp,
     - jin,jmi)

      use dflogm
      include 'Growth.fd'
      logical retlog,nada
      character*50 texto,obs,mod,fil
      external fin
      type(dialog) dlg
      retlog=DlgInit( IDD_Growth, dlg)
                                                    
      retlog = DlgSet(dlg,IDC_EDIT_e1,obs )
      retlog = DlgSet(dlg,IDC_EDIT_e3,mod )
      retlog = DlgSet(dlg,IDC_EDIT_e4,fil )

      retlog = DlgSet(dlg,IDC_TEXT_t1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t3,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t4,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_BOX_b1 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e3,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e4,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e7,.false.,dlg_enable)

      retlog = DlgSet(dlg,IDC_CHECK_c1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_CHECK_c2,.false.,dlg_enable)

      write(texto,'(i6)') ns
      retlog = DlgSet(dlg,IDC_EDIT_e30,texto)
      write(texto,'(i3)') na            
      retlog = DlgSet(dlg,IDC_EDIT_e5,texto)
      write(texto,'(i3)') nint(bal)      
      retlog = DlgSet(dlg,IDC_EDIT_e6,texto)
      write(texto,'(f3.0)') tpc     
      retlog = DlgSet(dlg,IDC_EDIT_e7,texto) 
      write(texto,'(i2)') ncer  
      retlog = DlgSet(dlg,IDC_EDIT_e13,texto) 
      write(texto,'(i2)') kup     
      retlog = DlgSet(dlg,IDC_EDIT_e11,texto) 
      retlog = DlgSet(dlg,IDC_CHECK_c2,.false.)
      retlog = DlgSet(dlg,IDC_CHECK_c1,.false.) 
      retlog = DlgSet(dlg,IDC_CHECK_c3,.false.)  

      retlog=DlgSet(dlg,IDC_TEXT_t2 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t5 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t6 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t7 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t8 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t8 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t11,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t13,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_TEXT_t14,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e5 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e6 ,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e11,.false.,dlg_enable)      
      retlog=DlgSet(dlg,IDC_EDIT_e12,.false.,dlg_enable)
      retlog=DlgSet(dlg,IDC_EDIT_e13,.false.,dlg_enable)      
      retlog=DlgSet(dlg,IDC_BOX_b4  ,.false.,dlg_enable)
     
      write(texto,'(i6)') nint(pp)
      if(jmi.eq.1) write(texto,'(f6.1)') pp
      nada=DlgSet(dlg,IDC_EDIT_e2,texto)
      
      retlog = DlgSetSub(dlg,IDCANCEL,fin)

      retint = DlgModal( dlg )

      retlog = DlgGet( dlg, IDC_EDIT_e2,texto)
      read(texto,*) pp

      jin=0
      retlog = DlgGet( dlg, IDC_CHECK_c3, nada)
      if(nada) jin=1
      
      return
      end
      
c***********************************************************************
c     Dialog interface 3  (inversion parameters)                       *
c***********************************************************************
      subroutine dial3(obs,mod,fil,na,nc,ns,bal,mlev,
     & pp,ig0,igxy,tpc,kup,ncer,jmi,jin)
      use dflogm
      include 'Growth.fd'
      logical retlog,nada
      character*50 texto,obs,mod,fil
      external fin
      type(dialog) dlg
      retlog=DlgInit( IDD_Growth, dlg)
                                                    
      retlog = DlgSet(dlg,IDC_EDIT_e1,obs )
      retlog = DlgSet(dlg,IDC_EDIT_e3,mod )
      retlog = DlgSet(dlg,IDC_EDIT_e4,fil )

      retlog=DlgSet(dlg,IDC_CHECK_c3,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t3,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t4,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_TEXT_t9,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_BOX_b1 ,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e1,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e2,.false.,dlg_enable)      
      retlog = DlgSet(dlg,IDC_EDIT_e3,.false.,dlg_enable)
      retlog = DlgSet(dlg,IDC_EDIT_e4,.false.,dlg_enable)
      write(texto,'(i6)') ns
      retlog = DlgSet(dlg,IDC_EDIT_e30,texto)
      write(texto,'(i6)') nc
      retlog = DlgSet(dlg,IDC_EDIT_e31,texto)
      write(texto,'(i3)') na            
      retlog = DlgSet(dlg,IDC_EDIT_e5,texto)
      write(texto,'(i3)') nint(bal)      
      retlog = DlgSet(dlg,IDC_EDIT_e6,texto)
      write(texto,'(f3.0)') tpc        
      retlog = DlgSet(dlg,IDC_EDIT_e7,texto)
      write(texto,'(i2)') kup     
      retlog = DlgSet(dlg,IDC_EDIT_e11,texto)
      write(texto,'(i2)') ncer     
      retlog = DlgSet(dlg,IDC_EDIT_e13,texto)

      if(ig0.eq.0) retlog = DlgSet(dlg,IDC_CHECK_c1,.false.)
      if(ig0.eq.1) retlog = DlgSet(dlg,IDC_CHECK_c1,.true.)
                   retlog = DlgSet(dlg,IDC_CHECK_c3,.false.)
      if(jin.eq.1) retlog = DlgSet(dlg,IDC_CHECK_c3,.true.)
      if(igxy.eq.0) retlog = DlgSet(dlg,IDC_CHECK_c2,.false.)
      if(igxy.eq.1) retlog = DlgSet(dlg,IDC_CHECK_c2,.true.)
      write(texto,'(i6)') nint(pp)
      if(jmi.eq.1) write(texto,'(f6.1)') pp
      nada=DlgSet(dlg,IDC_EDIT_e2,texto)
      retlog = DlgSetSub(dlg,IDCANCEL,fin)
      
      retint = DlgModal( dlg )
                                                  
      retlog = DlgGet( dlg, IDC_EDIT_e5,texto)
      read (texto,*) na
      if(na.lt.1) na=1
      retlog = DlgGet( dlg, IDC_EDIT_e6, texto)
      read(texto,*) bal
      retlog = DlgGet( dlg, IDC_EDIT_e7, texto)
      read(texto,*) tpc
      retlog = DlgGet( dlg, IDC_EDIT_e11, texto)
      read(texto,*) kup 
      retlog = DlgGet( dlg, IDC_EDIT_e13, texto)
      read(texto,*) ncer 
      ig0=0
      retlog = DlgGet( dlg, IDC_CHECK_c1, nada)
      if(nada) ig0=1
      igxy=0
      retlog = DlgGet( dlg, IDC_CHECK_c2, nada)
      if(nada) igxy=1

      return
      end
      
c***********************************************************************
c     Dialog interface 4  (inner parameters)                           *
c***********************************************************************
      subroutine dial4(blun,nlim)
      use dflogm
      include 'Growth.fd'
      logical retlog,nada
      external fin
      character*50 texto
      type(dialog) dlg
      retlog=DlgInit( IDD_Param, dlg)                         
   
      write(texto,'(f5.2)') blun
      retlog = DlgSet(dlg,IDC_EDIT_e1,texto)
      write(texto,'(i4)') nlim
      retlog = DlgSet(dlg,IDC_EDIT_e2,texto)

      retlog = DlgSetSub(dlg,IDCANCEL,fin)
      
      retint = DlgModal( dlg )
                                                  
      retlog = DlgGet( dlg, IDC_EDIT_e1,texto)
      read (texto,*) blun
      retlog = DlgGet( dlg, IDC_EDIT_e2, texto)
      read(texto,*) nlim

      call DlgUninit(dlg )
      return
      end


c***********************************************************************
c     Simple sensitivity calculus  Ton/uGal
c***********************************************************************
      subroutine ss2(ms,ns,x,y,z,fe,xx,yy,zz,sen,k)
      dimension x(ms),y(ms),z(ms)
      av=0
      r=1.e4*fe*fe
      do 1 j=1,ns                          
      xd=xx-x(j)
      yd=yy-y(j)
      zd=zz-z(j)                            
      d=xd*xd+yd*yd+zd*zd + r
      c=zd*zd/d/d/d
    1 av=av+c
      av=av/ns
      av=6.672e-3 *sqrt(av)*1.e8   
      sen=1./av                    
      k=sen/fe/fe/50.

      if(k.lt.1) k=1
      if(k.gt.12) k=12

      return
      end

c***********************************************************************   
c     Histogram for residuals,  n columns
c***********************************************************************
      subroutine histog(ms,ns,w,re,kdib,ix,lx,iy,ly)
      use IFQWIN
      use IFPORT
      type (xycoord) xy
      integer*2 u1,u2,v1,v2
      character*140 texto
      dimension re(ms),w(ms),k(20)
      sw=0
      do 5 i=1,ns
    5 sw=sw+w(i)
      sw=sw/ns
      rmin=0
      rmax=0
      rmed=0
      do 1 i=1,ns
      c=re(i)  *w(i)/sw
      if(c.gt.rmax) rmax=c    
      if(c.lt.rmin) rmin=c   
      rmed=rmed+c
    1 continue
      rmed=rmed/ns   
      n=ns/4.7 
      if(n.gt.20) n=20
      p=(rmax-rmin)/n
      py=0
      xce=0
      do 2 j=1,n
      k(j)=0
      r1=rmin+(j-1)*p
      r2=r1+p
      do 3 i=1,ns
      c=re(i) *w(i)/sw
    3 if(c.gt.r1.and.c.le.r2) k(j)=k(j)+1   
      if(k(j).gt.py) py=k(j)
    2 continue                                                          !

      nada=setcolorrgb(#000000)                                         !
      write(texto,'(i9,8x,i9)') nint(rmin),nint(rmax)                   !
      call moveto(ix-44,iy+5,xy)                                        !
      call outgtext(texto)                                              !
      write(texto,'(i4)') 0                                             !
      call moveto(ix-27,iy-10,xy)                                       !
      call outgtext(texto)                                              !
      write(texto,'(i4)') nint(py/2)                                    !
      j=iy-ly/1.8                                                       !
      call moveto(ix-27,j,xy)                                           !
      call outgtext(texto)                                              !
      write(texto,'(i4)') nint(py)                                      !
      call moveto(ix-27,iy-ly,xy)                                       !
      call outgtext(texto)                                              !
      if(kdib.eq.0) then                                                !
       write(texto,*) 'Resid.(uGal)'                                    !
       call moveto(ix+10,iy-170,xy)                                     !
       call outgtext(texto)                                             !
      endif                                                             !
      py=ly/py*0.9                                                      ! 
      px=lx/n                                                           !
      do 4 j=1,n                                                        !
      l=py*k(j)                                                         !
      u1=ix+5+(j-1)*p/(rmax-rmin)*lx +1                                 !
      u2=ix+5+j*p   /(rmax-rmin)*lx  -2                                 !
      v1=iy                                                             !
      v2=iy-l                                                           !
      nada=setcolorrgb(#60d060)                                         !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
    4 continue                                                          !                  
      nada=setcolorrgb(#000000)                                         !
      call moveto(ix+5,iy,xy)                                           !
      u1=ix+5                                                           !
      v1=iy-ly                                                          !
      nada=lineto(u1,v1)                                                !
      i=(rmed-rmin)/(rmax-rmin)*lx                                      !
      call moveto(ix+5+i,iy,xy)                                         !
      u1=ix+5+i                                                         !
      nada=lineto(u1,v1)                                                !
      write(texto,'(i4)') nint(rmed)                                    !
      call moveto(u1-20,iy+5,xy)                                        !
      call outgtext(texto)                                              !
      call moveto(ix+5,iy,xy)                                           !
      u1=ix+5+lx                                                        !
      v1=iy                                                             !
      nada=lineto(u1,v1)                                                !
   99 return                                                            !
      end                                                               !

c***********************************************************************!
c     Draw a section of the 3D anomalous density model                  !
c***********************************************************************!
      subroutine sect(jx,jy,mc,nc,mb,nb,xc,yc,zc,tx,ty,tz,df,isen,pa,   !
     &  zi,lay,xb,yb,ib,col,fe,iym,izm,fdib,tic,seli,jx1,jy1,itipo,     !
     & ival,le,lk)                                                      !
      use IFQWIN                                                        !
      use IFPORT                                                        !
      type (xycoord) xy                                                 !
      dimension xb(mb),yb(mb),ib(mb),                                   !
     &xc(mc),yc(mc),zc(mc),tx(mc),ty(mc),tz(mc),df(mc),isen(mc),lay(mc) !
      integer col(21),nada                                              !  
      integer*2 u1,u2,v1,v2                                             !
      character*120 texto                                               !
      character*1 le(15)                                                !
      lx=190                                                            !
      ly=190                                                            !
      lz=147                                                            !
      if(itipo.eq.1) then                                               !
      call ventana(jx-lx,jy-ly,jx+lx,jy+ly)                             !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx-lx                                                          !
      u2=jx+lx                                                          !
      v1=jy-ly                                                          !
      v2=jy+ly                                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      call hor0(jx,jy,lx,ly,tic,fe,fdib,mb,nb,xb,yb,ib)                 !
      do 3 i=1,nc                                                       !
      c=tz(i)/2.                                                        !
      if(isen(i).gt.4.5*seli) go to 3                                   !
      if(isen(i).gt.     0) r=2.                                        !
      if(isen(i).gt.1.5*seli) r=2.9                                     !
      if(isen(i).gt.2.0*seli) r=3.8                                     !
      if(isen(i).gt.2.5*seli) r=4.7                                     !
      if(isen(i).gt.3.0*seli) r=5.6                                     !
      if(isen(i).gt.3.5*seli) r=6.5                                     !
      if(isen(i).gt.4.0*seli) r=7.4                                     !
      if((zc(i)+c).lt.ival) go to 3                                     !
      if((zc(i)-c).ge.ival) go to 3                                     !
      j=(xc(i)+tx(i)/r)/fe/fdib                                         !
      if(j.gt.lx) go to 3                                               !
      u2=jx+j                                                           !
      j=(xc(i)-tx(i)/r)/fe/fdib                                         !
      if(j.lt.-lx) go to 3                                              !
      u1=jx+j                                                           !
      j=(yc(i)+ty(i)/r)/fe/fdib                                         !
      if(j.gt.ly) go to 3                                               !
      v2=jy-j                                                           !
      j=(yc(i)-ty(i)/r)/fe/fdib                                         !
      if(j.lt.-ly) go to 3                                              !
      v1=jy-j                                                           !
      k= (df(i)-zi)/pa+1                                                !
      if(k.gt.21) k=21                                                  !
      if(k.lt.1) k=1                                                    !
      nada=setcolorrgb(col(k))                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
    3 continue                                                          !
      write(texto,'(a1,8x,a,i9)') le(lk),'Depth(m) bsl =',ival+izm      ! 
      call moveto(jx-90,jy-ly-20,xy)                                    !
      nada=setcolorrgb(#000000)                                         !
      call outgtext(texto)                                              !
      k=ival/fe/fdib                                                    !
      write(texto,'(a1)') le(lk)                                        ! 
      call moveto(jx1-lx-15,jy1+22-k,xy)                                ! 
      call outgtext(texto)                                              !
      nada=setcolorrgb(#c000c0)                                         !
      ix=jx-lx                                                          !
      call moveto(jx1-lx,jy1+30-k,xy)                                   !
      u1=jx1+lx                                                         !
      v1=jy1+30-k                                                       !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
      endif                                                             !

      if(itipo.eq.2) then                                               !
      call ventana(jx-lx,jy,jx+lx,jy+lz)                                !
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx-lx                                                          !
      u2=jx+lx                                                          !        
      v1=jy+30                                                          !
      v2=jy+lz                                                          !        
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
      do 6 i=1,nc                                                       !
      j=zc(i)/fe/fdib                                                   ! 
      if(j.lt.-lz+30) go to 6                                           !
      c=ty(i)/2.                                                        !
      if(isen(i).gt.4.5*seli) go to 6                                   !
      if(isen(i).gt.     0) r=2.                                        !                                          
      if(isen(i).gt.1.5*seli) r=2.9                                     !
      if(isen(i).gt.2.0*seli) r=3.8                                     !
      if(isen(i).gt.2.6*seli) r=4.7                                     !
      if(isen(i).gt.3.0*seli) r=5.6                                     !
      if(isen(i).gt.3.5*seli) r=6.5                                     !
      if(isen(i).gt.4.0*seli) r=7.4                                     !
      if((yc(i)+c).lt.ival) go to 6                                     !
      if((yc(i)-c).ge.ival) go to 6                                     !
      j=(xc(i)+tx(i)/r)/fe/fdib                                         !
      if(j.gt.lx) go to 6                                               !
      u2=jx+j                                                           !
      j=(xc(i)-tx(i)/r)/fe/fdib                                         !
      if(j.lt.-lx) go to 6                                              !
      u1=jx+j                                                           !
      j=(zc(i)+tz(i)/r)/fe/fdib                                         !
      if(j.gt.lz) go to 6                                               !
      v2=jy-j   +30                                                     !
      j=(zc(i)-tz(i)/r)/fe/fdib                                         !
      if(j.lt.-lz) go to 6                                              !
      v1=jy-j  +30                                                      !
      k=(df(i)-zi)/pa+1                                                 !
      if(k.gt.21) k=21                                                  !
      if(k.lt.1) k=1                                                    !
      nada=setcolorrgb(col(k))                                          !
      nada=rectangle($gfillinterior,u1,v1,u2,v2)                        !
    6 continue                                                          !
      call ver0(jx,jy,lx,ly,lz,izm,tic,fe,fdib)                         !
      write(texto,'(a1,8x,a,i9)') le(lk),'Y UTM (m) =',ival+iym         !
      call moveto(jx-90,jy-20,xy)                                       !
      nada=setcolorrgb(#000000)                                         !
      call outgtext(texto)                                              !
      k=ival/fe/fdib                                                    !
      write(texto,'(a1)') le(lk)                                        ! 
      call moveto(jx1-lx-15,jy1-k-10,xy)                                ! 
      call outgtext(texto)                                              !
      nada=setcolorrgb(#c000c0)                                         ! 
      ix=jx-lx                                                          !
      call moveto(jx1-lx,jy1-k,xy)                                      !
      u1=jx1+lx                                                         ! 
      v1=jy1-k                                                          !
      if(k.lt.2*lx) nada=lineto(u1,v1)                                  !
      endif                                                             !
c-----------------------------------------------------------------------!

   99 return
      end

c************************************************************************
c     Draw a background for a plan view of the model
c************************************************************************
      subroutine hor0(jx,jy,lx,ly,tic,fe,fdib,mb,nb,xb,yb,ib)           !
      use IFQWIN                                                        !
      use IFPORT                                                        !
      type (xycoord) xy                                                 !
      dimension xb(mb),yb(mb),ib(mb)                                    !
      integer*2 u1,u2,v1,v2                                             !  
      nada=setcolorrgb(#d0f0d0)                                         !
      u1=jx-lx                                                          !
      u2=jx+lx                                                          !
      v1=jy-ly                                                          !
      v2=jy+ly                                                          !
      nada=rectangle($gfillinterior,u1,v1 ,u2,v2)                       !
      nada=setcolorrgb(#000000)                                         !
      do 1 i=-10,10                                                     !
      ix=jx-lx                                                          !
      k=i*tic/fe/fdib                                                   !
      v1=jy+k                                                           !
      call moveto(ix,v1,xy)                                             !
      u1=ix+4                                                           !
      if(abs(k).lt.ly) nada=lineto(u1,v1)                               !
      ix=jx+lx                                                          !
      call moveto(ix,v1,xy)                                             !
      u1=ix-4                                                           !
      if(abs(k).lt.ly) nada=lineto(u1,v1)                               !
      ix=jx+k                                                           !
      iy=jy-ly                                                          !
      call moveto(ix,iy,xy)                                             !
      u1=ix                                                             !
      v1=iy+4                                                           !
      if(abs(k).lt.lx) nada=lineto(u1,v1)                               !
      iy=jy+ly                                                          !
      call moveto(ix,iy,xy)                                             !
      v1=iy-4                                                           !
      if(abs(k).lt.lx) nada=lineto(u1,v1)                               !
    1 continue                                                          !
      call cont(mb,nb,xb,yb,ib,jx,jy,fdib)                              !
      return                                                            !
      end                                                               !   

c***********************************************************************!
c     Draw a background for a elevation view of the model               !
c***********************************************************************!
      subroutine ver0(jx,jy,lx,ly,lz,izm,tic,fe,fdib)                   !
      use IFQWIN                                                        !
      use IFPORT
      type (xycoord) xy
      integer*2 u1,u2,v1,v2

      v1=jy+30+izm/fe/fdib                                         
      nada=setcolorrgb(#f0c000)                                        
      call moveto(jx-lx,v1,xy)  
      u1=jx+lx                                       
      nada=lineto(u1,v1) 

      nada=setcolorrgb(#000000)  
      j=izm/tic+2                                                      
      do 1 i=-10,10                                               
      k=((i-j)*tic+izm)/fe/fdib
      v1=jy+30+k            
      if(v1.lt.(jy+lz).and.v1.gt.jy+10) then
        ix=jx-lx            
        call moveto(ix,v1,xy) 
        u1=ix+4               
        nada=lineto(u1,v1)    
        ix=jx+lx              
        call moveto(ix,v1,xy) 
        u1=ix-4               
        nada=lineto(u1,v1)                              
      endif

      k=i*tic/fe/fdib
      ix=jx+k         
      call moveto(ix,jy,xy) 
      u1=ix                 
      v1=jy+5               
      if(abs(k).lt.lx) nada=lineto(u1,v1)   
    1 continue

      return
      end

c*********************************************************************
c       End
c*********************************************************************
      subroutine fin
      use dialogm
      include 'Growth.fd'
      stop
      end 

c*********************************************************************
c   Iterative reweighting of the data according to previous residuals
c*********************************************************************
      subroutine rewe(m,n,re,blun,w)
      dimension re(m),w(m)  
      rme=0
      do 1 i=1,n
      ra=abs(re(i))
      rme=rme+ra
    1 w(i)=ra
      rme=rme/n/0.6745 *blun
      c=4  !5
      sw=0
      do 2 i=1,n
      t=abs(re(i))/rme-1
      w(i)=1./(1+exp(c*t))
      if(w(i).lt.0.01) w(i)=0.01  
    2 sw=sw+w(i)
      sw=sw/n
      do 3 i=1,n
    3 w(i)=w(i)/sw
      return
      end

c************************************************************************
c     Draw contour line xb,yb in map.bln
c************************************************************************
      subroutine cont(mb,nb,xb,yb,ib,jx,jy,fdib)                        !
      use IFQWIN                                                        !
      use IFPORT                                                        !
      type (xycoord) xy                                                 !  G
      dimension xb(mb),yb(mb),ib(mb)                                    !  R
      integer*2 u1,v1                                                   !  A
      nada=setcolorrgb(#000000)                                         !  P
      do 2 i=2,nb                                                       !  H
      if(ib(i-1).eq.1) go to 2                                          !  I 
      if(ib(i).eq.1.or.ib(i).eq.2) go to 2                              !  C
      ix=xb(i-1)/fdib+jx                                                !  S
      iy=jy-yb(i-1)/fdib                                                !
      call moveto(ix,iy,xy)                                             !
      u1=xb(i)/fdib+jx                                                  !
      v1=jy-yb(i)/fdib                                                  !
      nada=lineto(u1,v1)                                                !
    2 continue                                                          !
      return                                                            !
      end                                                               !

c*********************************************************************
c   Up and down cells close to each cell
*********************************************************************
      subroutine updown(mc,nc,xc,yc,tz,lay,mpo,iup,ido)
      dimension xc(mc),yc(mc),tz(mc),lay(mc),iup(mpo,mc),ido(mpo,mc)
      dimension xd(200),ix(200),id(200)  

        num=0
        ndm=0
 
      do 1 i=1,nc
           j=i-int(i/1000)*1000
           if(j.eq.0) write(*,'(\a1)') '.'      
      l=lay(i)
      if(l.le.0) go to 1
      xx=xc(i)
      yy=yc(i)
            dm=4*tz(i)    
            dup=dm*dm
            ddo=dup*1.6

      do 2 k=1,mpo
    2 iup(k,i)=0
      nup=0
      li=l-1
      do 3 j=1,nc
      if(lay(j).ne.li) go to 3
        dx=xx-xc(j)   
        dy=yy-yc(j)
        d=dx*dx+dy*dy
        if(d.gt.dup) go to 3  
        nup=nup+1 
        if(nup.gt.mpo) go to 5
        xd(nup)=d
        id(nup)=j
    3 continue
    5 continue
      call orden(200,nup,xd,ix)
      do 4 j=1,nup 
    4 iup(j,i)=id(ix(j))

      do 6 k=1,mpo
    6 ido(k,i)=0
      ndo=0
      li=l+1
      do 7 j=1,nc
      if(lay(j).ne.li) go to 7
        dx=xx-xc(j)  
        dy=yy-yc(j)
        d=dx*dx+dy*dy
        if(d.gt.ddo) go to 7        
        ndo=ndo+1 
        if(ndo.gt.mpo) go to 8
        xd(ndo)=d
        id(ndo)=j
    7 continue
    8 continue
      call orden(200,ndo,xd,ix)
      do 9 j=1,ndo 
    9 ido(j,i)=id(ix(j))

         if(ndo.gt.ndm) ndm=ndo
         if(nup.gt.num) num=nup
    1 continue   

      return
      end
c*********************************************************************
c   Ordena numeros positivos
*********************************************************************
      subroutine orden(mx,nx,x,ix)
      dimension x(mx),ix(mx)

      do 2 i=1,nx
      xmin=9.d9
      do 3 j=1,nx
      if(x(j).lt.0) go to 3 
      if(x(j).ge.xmin) go to 3      
      xmin=x(j)
      jmin=j 
    3 continue 
      ix(i)=jmin
      x(jmin)=-x(jmin)
    2 continue

      return
      end
c*********************************End*************************************