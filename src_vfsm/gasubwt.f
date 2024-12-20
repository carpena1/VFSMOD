       SUBROUTINE GASUBWT(TIME,DT,L,R,RAIN,NEND,TRAI)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                            C
C     Authors: Rafael Muñoz-Carpena/Claire Lauvernet, 10/16/2010             C
C     Copyright 2010 University of Florida/CEMAGREF. All rights reserved.    C
C                                                                            C
C     units: m,s                                                             C
C                                                                            C
C     The subroutine solves the soil infiltration problem for unsteady rain  C
C     in the presence of a shallow water table using a modified Green-Ampt   C
C     infiltration model as proposed by Salvucci and Entekhabi (1995), Chu   C
C     (1996) and work by the authors of this oce. The method was extended    C
C     to include surface mass balance e as proposed by Skaggs (1982) and     C
C     Khaleel (1982) in Hydrologic modeling of small watersheds, ASAE mon.   C
C     no. 5, and Chu (1978, Water Resour. Res).                              C
C                                                                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                            C
C     THIS IS A MODIFIED VERSION OF THE ORIGINAL METHOD. IN THIS CASE        C
C     THE INFILTRATION IS ALLOWED TO BE ITS MAXIMUM POTENTIAL AFTER THE      C
C     FIRST PODING. THE IMPLICATION HERE IS THAT AFTER THE ORIGINAL          C
C     PONDING IS ACHIEVED, THE RUNOFF WATER MOVING AT THE SURFACE WILL       C
C     SUPPLY ENOUGH WATER TO SUSTAIN THE MAXIMUM POSSIBLE INFILTRATION       C
C     FOR THAT TIME STEP, IN OTHER WORDS, THE EFFECTIVE RAILFALL FEEDED      C
C     INTO THE MAIN FE-PROGRAM WILL BE R<0, INDICATING INFILTRATING PLANE    C
C                                                                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                            C
C     Parameters used:                                                       C
C     WTD: Water table depth (m                                              C
C     hb: Soil water bubling pressure (m)                                    C
C     ITHETATYPE: integer to select the soil water characteristic curve      C
C     type with values: 1=van Genuchten, 2=Brooks and Corey.                 C
C     Note that SAV and OI are now calculated internally, but some numbers   C
C     must be included to ensure backwards compatibility, even although      C
C     these numbers are ignored.                                             C
C     Equation   ITHETATYPE PAR(1) PAR(2)     PAR(3)    PAR(4)               C
C     van Genuchten 1      OR   VGALPHA  VGN     VGM                         C
C     Brooks&Corey  2      OR   BCALPHA  BCLAMBDA --                         C
C                                                                            C
C     IKUNSTYPE: an integer to select the unsaturated hydraulic conductivity C
C     curve type with values: 1=van Genuchten, 2=Brooks and Corey, 3=Gardner C
C     Equation   IKUNSTYPE PARK(1)   PARK(2)                                 C
C     van Genuchten 1      VGM     --                                        C
C     Brooks&Corey  2      BCETA   BCALPHA                                   C
C     Gardner       3      GDALPHA                                           C
C     Also, VKs (m/s) is read as part of the common soil properties in ISO   C
C                                                                            C
C     NPFORCE= New surface ponding forcing scheme, 1: when overland flow at  C
C     check node and zw>0, t<tw the infiltration is at capacity (CASE 3.2)   C
C     regardless of ponding or not from point excess calculation             C
C     WTSC=slope of water table (=surface slope) when t>tw (02/2012)         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      COMMON/GA2/LO,NPOND,NPFORCE,NPFORCE0
      COMMON/WTGA1/WTD,PARW(4),PARK(2),hb,RAINL,ZW,TW,RVH
      COMMON/WTGA2/ITHETATYPE,IKUNSTYPE,ITWBC
      COMMON/GRASSD/PART(3),SC,SS,VN1,VN2,VN,GSI,H,VLCM,POR,CI,ICO
      DIMENSION RAIN(200,2)

c-----Set some variable names for use in functions
      WTSC=SC
      VKS=AGA
      RAINL=RAIN(L,2)
      WTD=-dabs(WTD)
      aWTD=-WTD
      wdt=0.d0
      NPDIFF=NPFORCE-NPFORCE0
      NPFORCE0=NPFORCE

c----Select bottom BC: Three types of boundary conditions (ITWBC) are included:
c-------Type 1 (default): Dupuis-Forchheimer, f=Ksh.aWTD/VL.So
c-------Type 2: Vertical saturated flow (Vachaud), f=Ks
c-------Type 3: Simplified, f=Ksh.So
c-----The horizontal Ks calculated from input RVH=VKS/VKSH (default=1)
      VKSH=VKS/RVH
      IF(ITWBC.EQ.1) THEN
           FPTW=VKSH*WTSC*aWTD/(VLCM/100.d0)
      ELSEIF(ITWBC.EQ.2) THEN
           FPTW=VKS
       ELSE
           FPTW=VKSH*WTSC
      ENDIF

c-----Check if end of field runoff is reached and if so reset case as   -----
c-----no-ponding at beginning and end of the period (NPOND=0,CU<0)(Case 1b)--
      IF(NEND.EQ.1) NPFORCE=0

c------Check if time step belongs to same rainfall period as before(L=LO) ---
c------and if not reset the total rainfall value for the end of the period --
      IF(L.NE.LO) THEN
       PSOLD=PS
       PST = PS + RAIN(L,2)*(RAIN(L+1,1)-RAIN(L,1))
c---------rmc03/2011- missing case for NPOND=1 and new_rain<end_infiltration-
c---------of last period (fp,capacity since soil was ponded, i.e. NPOND=1)---
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   CASE 1. zw<0 (soil saturated by WT at "surface")                      C
C   CASE 2. zw>0,t>=tw (WT not on surface but time after column saturation)  C
C   CASE 3. zw>0,t<tw (WT not on surface and time before column saturation)  C
C         3.1. t<tp before ponding in the rainfall peeriod                C
C         3.2. t>=tp ater ponding in the rainfall period                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c---- Check shallow water table (zw=|L|-|hb|) position
      if(zw.le.0.d0) then
c--------CASE 1 (zw=<0, water table (or capillary fringe) at surface.
c----------Three types of boundary conditions (ITWBC) are considered for
c----------fw (FPTW) based on ITWBC selected by user (1 by default)
          IF(NPFORCE.eq.1) then
             FPI= FPTW
            ELSE
             FPI=MIN(FPTW,RAIN(L,2))
          ENDIF
          F=F+DT*FPI
          z=0.d0
          fcase=1.d0
        else
c--------CASES for zw>0, water table (capillary fringe) below surface for new
c--------rain (L.ne.LO) period calculate first singular times tp,to,tw
          if((L.ne.LO.or.NPDIFF.ne.0).and.(time.lt.tw.or.tw.eq.0.d0.
     &          or.LO.eq.0).and.rain(L,2).gt.0.d0)then
             fpiold=fpi
             if(rain(L,2).le.vks.and.NPFORCE.eq.0) then
c------------------get Fw(tw), dFw=Fw-F, tw=t+dFw/rain
                 call qgaus(3,-WTD-zw,-WTD,20,dFw)
                 Fw=Os*zw-dFw
                 dFzw=Fw-F
                 if (rain(L,2).gt.0.d0) tw=rain(L,1)+dFzw/rain(L,2)
c----------------- simple estimate for tp (large enough)
                tp=rain(L+1,1)*10.d0
               else
c------------------normal case, calculate singular times tp, to, tw
c------------------get potential FPI at the beginning of the new period to
c------------------compare with rain
                 zist=dabs(WTD)-dabs(z)
                 call qgaus(1,0.d0,zist,20,vkint)
                 fpr=vks+(1.d0/z)*vkint*vks
c------------------if new rainfall period after no-ponding calculate tp,
c------------------tpp, the cum. infiltration at tp (Ftp), and if there is
c------------------ponding at the end of period (CU>0)
                 if(RAIN(L,2).gt.FPR.OR.NPFORCE.eq.1) then
                    zp=z
                    tp=rain(L,1)
                    IF (NPDIFF.ne.0) tp=time
                    Ftp=F
                  else
                    tpold=tp
                    tppold=tpp
                    Ftpold=Ftp
                    call nrsol1(time,0.d0,dabs(WTD),zp,3)
                    call qgaus(3,-WTD-zp,-WTD,20,dFtp)
                    Ftp=Os*zp-dFtp
                    tp=(Ftp-PS+RO)/RAIN(L,2)+RAIN(L,1)
                 endif
c------------------find shifting time (tpp)
                 call qgaus(2,0.d0,zp,20,tpp)
c------------------find time to reach point in profile soil saturated by
c------------------shallow water table (zw=|L|-|hb|)
                 call qgaus(2,0.d0,zw,20,tw)
                 tw=tw+tp-tpp
             endif
             fpi=fpiold
          endif

c----start (z>zw) cases for all time steps based 1st) on tw; 2nd) on tp
          if(time.gt.tw.and.tw.ne.0.d0)then
c-RMC03/13-----CASE 2: t>=tw, z=zw>0 after wetting front reaches WT.
c--------------Three boundary conditions (ITWBC) are included:
c--------------Type BC1 (default): Dupuis-Forchheimer, f=Ksh.aWTD/VL.So
c--------------Type BC2: Vertical saturated flow (Salvucci & Entekabi),f=Ks
c--------------Type BC3: Simplified, f=Ksh.So
c--------------The horizontal Ks from input RVH=VKS/VKSH (default=1)
             IF(NPFORCE.eq.1) then
                  FPI= FPTW
               ELSE
                  FPI=MIN(FPTW,RAIN(L,2))
             ENDIF
             F=F+DT*FPI
             z=zw
             fcase=2.d0
            elseif(time.lt.tp.and.NPFORCE.eq.0)then
c--------------CASE 3.1: t<tp<tw,z<zw>0, before ponding (w/wo initial ponding)
             FPI=RAIN(L,2)
             F=F+DT*FPI
             if(rainL.gt.0.d0.and.(time+wdt*dt).lt.tw)then
                 call nrsol1(time,z,zw,zF,4)
                 z=zF
             endif
             fcase=3.1d0
            else
c--------------CASE 3.2: tp<t<tw,z<zw>0,after ponding OR NPFORCE=1 (GA infilt.
c--------------capacity)
             zist=dabs(WTD)-dabs(z)
             if(rain(L,2).gt.FPI.or.NPFORCE.eq.1)then
c------------------get Fw(tw), dFw=Fw-F, tw=t+dFw/rain
                 if(NPDIFF.ne.0) then
                    tp = time
c------------------find shifting time (tpp)
                    call qgaus(2,0.d0,z,20,tpp)
                    call qgaus(2,0.d0,zw,20,tw)
                    tw=tp-tpp+tw
                    fpi=0.d0
                  else
                    call nrsol1(time,z,zw,z,2)
                    call qgaus(1,0.d0,zist,20,vkint)
                    fpi=vks+(1.d0/z)*vkint*vks
                 endif
                 call qgaus(3,-WTD-z,-WTD,20,dF)
                 F=Os*z-dF
              else
                 FPI=RAIN(L,2)
                 F=F+DT*FPI
                 if(rainL.gt.0.d0.and.(time+wdt*dt).lt.tw)then
                   call nrsol1(time,z,zw,zF,4)
                   z=zF
                 endif
             endif
             fcase=3.2d0
          endif
      endif

c-----Left hand side of kinematic wave equation (re=i-f)
      R=RAIN(L,2)-FPI
c-----mass balance (Chu (1978),Skaggs and Khaleel, 1982)
      IF((TIME+DT).GT.RAIN(L+1,1)) THEN
         PS=PST
      ELSE
         PS= PSOLD + RAIN(L,2)*(TIME-RAIN(L,1))
      ENDIF
      EXCESS= PS -F - RO -STO
      STO =STO +EXCESS
      IF (STO.GT.SM) THEN
         ROI=STO-SM
         STO=SM
         RO =RO +ROI
      ENDIF

c--for output debug--gasubwt-----
cdb      write(*,100)time,tw,z,F,PS,RO,FPI,RAINL,tp,tpp,NPFORCE,fcase,
cdb     &            CU,zw,L,LO
cdb 100   format(2f12.1,6e11.3,2f12.4,i4,f5.2,2f10.4,2i4)
c--for debug--gasubwt-----

c-----Set rainfall period for next time step as the current one and total
c-----cummulative rain to the one just calculated in this time step
      LO=L
      TRAI=PS

      RETURN
      END


      subroutine qgaus(ieq,a,b,ngl,sst)
C
C  Numerical integration of a generic function using Gauss-Legendre quadrature
C   ieq: function to integrate (1:Kuns; 2: time integral; 3: water content)
C   a,b: limits of integration
C   ngl: orger of the gauss-Legendre quadrature
C   sst: value of integral
C
      implicit double precision (a-h,o-z)
      COMMON/WTGA2/ITHETATYPE,IKUNSTYPE,ITWBC
      COMMON/CINT/XI(20,20),W(20,20)
      external vKuns,func

c------Set selected soil hydraulic functions------
      ic=ithetatype
      ikc=ikunstype

      xm=0.5d0*(b+a)
      xr=0.5d0*(b-a)
      sst=0.d0
      do 11 j=1,ngl
        dx=xr*xi(j,ngl)
        if(ieq.eq.1) then
            sst=sst+w(j,ngl)*vKuns(xm+dx,ic,ikc)
         elseif(ieq.eq.2) then
            sst=sst+w(j,ngl)*func(xm+dx)
         else
            sst=sst+w(j,ngl)*swcc(xm+dx,ic)
        endif
11     continue
      sst=xr*sst

      return
      end


      subroutine nrsol1(time,x1,x2,xval,ieq)
C
C  Newton-Raphson method for time implicit equations
C   ieq: function to integrate (1-Kuns; 2-time integral equ.; 3-water content)
C   xval: value to optimize (implicit result)
C
C   Credit: modified from 1986-92 Numerical Recipes Software (C) #>0K!.

      implicit double precision (a-h,o-z)
      EXTERNAL gfGA,gfSvc,gtp,gFz

      iter = 0
      tol=1D-08
      error=100.d0
      maxiter=350
      if(ieq.eq.1) then
            call gfGA(time,x1,fl,dfval)
            call gfGA(time,x2,fh,dfval)
        elseif(ieq.eq.2) then
            call gfSvc(time,x1,fl,dfval)
            call gfSvc(time,x2,fh,dfval)
        elseif(ieq.eq.3) then
            call gtp(x1,fl,dfval)
            call gtp(x2,fh,dfval)
        else
            call gFz(x1,fl,dfval)
            call gFz(x2,fh,dfval)
      endif
      if((fl.gt.0.d0.and.fh.gt.0.d0).or.(fl.lt.0.d0.and.fh.
     &   lt.0.d0)) then
c         write(*,'(A39)')'WARNING:root must be bracketed in NRSOL'
         xval=x1
         return
      endif
      if(fl.eq.0.d0)then
          xval=x1
          return
        else if(fh.eq.0.d0)then
            xval=x2
            return
        else if(fl.lt.0.d0)then
            xl=x1
            xh=x2
        else
            xh=x1
            xl=x2
      endif
      xval=.5d0*(x1+x2)
      dxvalold=abs(x2-x1)
      dxval=dxvalold
      if(ieq.eq.1) then
            call gfGA(time,xval,fval,dfval)
        elseif(ieq.eq.2) then
            call gfSvc(time,xval,fval,dfval)
        elseif(ieq.eq.3) then
            call gtp(xval,fval,dfval)
        else
            call gFz(xval,fval,dfval)
      endif

      do 11 j=1,maxiter
          if(((xval-xh)*dfval-fval)*((xval-xl)*dfval-fval).ge.0.d0.
     &     or.abs(2.d0*fval).gt.abs(dxvalold*dfval) ) then
            dxvalold=dxval
            dxval=0.5d0*(xh-xl)
            xval=xl+dxval
            if(xl.eq.xval)return
            else
            dxvalold=dxval
            dxval=fval/dfval
            temp=xval
            xval=xval-dxval
            if(temp.eq.xval)return
          endif
          if(abs(dxval).lt.tol) return
          if(ieq.eq.1) then
            call gfGA(time,xval,fval,dfval)
            elseif(ieq.eq.2) then
            call gfSvc(time,xval,fval,dfval)
            elseif(ieq.eq.3) then
            call gtp(xval,fval,dfval)
            else
            call gFz(xval,fval,dfval)
          endif
          if(fval.lt.0.d0) then
            xl=xval
            else
            xh=xval
          endif
11     continue
      write(*,12)iter,error
      stop
12    format(2x,'ERROR: The solution did not converge: iter=',i4,
     2  '; error=',f8.4)
      end


      subroutine gFz(xval,fval,dfval)
C
C   Find the equivalent saturated front depth for a given value of F
      implicit double precision (a-h,o-z)
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      COMMON/WTGA1/WTD,PARW(4),PARK(2),hb,RAINL,ZW,TW,RVH
      COMMON/WTGA2/ITHETATYPE,IKUNSTYPE,ITWBC
      external swcc

      ic=ithetatype
      call qgaus(3,-WTD,-WTD-xval,20,F1)
      fval= F-OS*xval-F1
      dfval=-OS + swcc(dabs(WTD)-xval,ic)

      return
      end


      subroutine gtp(xval,fval,dfval)
C
C   Time to ponding calculation
      implicit double precision (a-h,o-z)
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      COMMON/WTGA1/WTD,PARW(4),PARK(2),hb,RAINL,ZW,TW,RVH
      COMMON/WTGA2/ITHETATYPE,IKUNSTYPE,ITWBC
      external vKuns

      ic=ithetatype
      ikc=ikunstype
      vks=AGA

      call qgaus(1,0.d0,dabs(-WTD-xval),20,vkint)
      fval= xval - 1.d0/(rainL-vKs)*vkint*vKs
      dfval=1.d0 + vKs/(rainL-vKs)*vKuns(-WTD-xval,ic,ikc)

      return
      end


      subroutine gfGA(time,bff,fval,dfval)
C
C   Regular Green-Ampt (1911) infiltration equation, no water table
      implicit double precision (a-h,o-z)
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS

      vks=AGA
      fval=bff-BGA/vks* log (1.d0 + (bff/(BGA/vks)))
     &      - vks * (time-tp+tpp)
      dfval= 1.d0 - ((BGA/vks) /(BGA/vks + bff))

      return
      end


      subroutine gfSvc(time,xval,fval,dfval)
C
C   Integral step in Salvucci and Enthekabi (1995) function for soil
C      infiltration with shallow water table
      implicit double precision (a-h,o-z)
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      external func

      a=0.d0
      b=xval

      call qgaus(2,a,b,20,fval)
      fval=time-(fval+tp-tpp)
      dfval=-func(xval)
      return
      end


      function vKuns(x,ic,ikc)
c   Unsaturated hydraulic conductivity. ic controls the SWCC user function
c     ikc=1, Mualem/van Genuchten
c     ikc=2, Brooks and Corey
c     ikc=3, Garner
      implicit double precision (a-h,o-z)
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      COMMON/WTGA1/WTD,PARW(4),PARK(2),hb,RAINL,ZW,TW,RVH

      WOR=PARW(1)
c---  Case 1:  Van Genuchten/Mualem
      IF(ikc.eq.1) THEN
           VGM=PARK(1)
           IF(dabs(x).gt.dabs(hb)) THEN
              Se=(swcc(x,ic)-WOR)/(OS-WOR)
            ELSE
              Se=1.d0
           ENDIF
           vKuns=(Se**(0.5d0))*(1.d0-(1.d0-Se**(1.d0/VGM))**VGM)**2.d0
      ENDIF
c--- Case 2: Brooks and Corey
      IF(ikc.eq.2) THEN
           BCETA=PARK(1)
           BCALPHA=PARK(2)
           IF(dabs(x).gt.dabs(hb)) THEN
              Se = (swcc(x,ic)-WOR)/(OS-WOR)
            ELSE
              Se=1.d0
           ENDIF
           vKuns=Se**BCETA
      ENDIF
c--- Case 3: Gardner
      IF(ikc.eq.3) THEN
         GNDALPHA=PARK(1)
         vKuns=dexp(-GNDALPHA*x)
      ENDIF

      RETURN
      END


      function swcc(h,ic)
c   Soil water characteristic curves. ic controls the function selected
c     ic=1, vn Genuchten
c     ic=2, Brooks and Corey
      implicit double precision (a-h,o-z)
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      COMMON/WTGA1/WTD,PARW(4),PARK(2),hb,RAINL,ZW,TW,RVH

c--- Case 1:  Van Genuchten/Mualem
      IF(ic.eq.1) THEN
           WOR=PARW(1)
           VGALPHA=PARW(2)
           VGN=PARW(3)
           VGM=PARW(4)
           IF(dabs(h).gt.dabs(hb)) THEN
                  Se=(1.d0+(dabs(VGALPHA*h))**VGN)**(-VGM)
            ELSE
                  Se=1.d0
           ENDIF
           swcc = WOR + Se*(OS-WOR)
c--- Case 2: Brooks and Corey
         ELSEIF(ic.eq.2) THEN
           WOR=PARW(1)
           BCALPHA=PARW(2)
           BCLAMBDA=PARW(3)
           IF(dabs(h).le.dabs(hb)) THEN
              swcc=OS
            ELSE
              swcc=WOR+(OS-WOR)*(BCALPHA*(dabs(h)))**(-BCLAMBDA)
           ENDIF
      ENDIF

      return
      end

      function func(zx)
c   Time integral equation (Chu, Salvucci function)
      implicit double precision (a-h,o-z)
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      COMMON/WTGA1/WTD,PARW(4),PARK(2),hb,RAINL,ZW,TW,RVH
      COMMON/WTGA2/ITHETATYPE,IKUNSTYPE,ITWBC

      ic=ithetatype
      vks=AGA

      zist=dabs(WTD)-dabs(zx)
      call qgaus(1,0.d0,zist,20,vkint)
      FPP=vks+(1.d0/zx)*vkint*vks
      OZ=swcc(zist,ic)
      func=(OS-OZ)/FPP

      return
      end
