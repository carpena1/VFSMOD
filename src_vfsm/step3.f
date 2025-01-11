      SUBROUTINE STEP3(GS2,TIME,NTRCAP,COARSE,QOUT,FWIDTH,ISCR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                          C
C     This program solves step 3 to find dfs, RSSs, and Se to go to step   C
C     4. After  Barfield et. al (1979)"Filtration of sediment by simulated C
C     vegetation II. Trans. ASAE, 22(5):540-548. and Hayes et. al (1984)   C
C     "Performance of Grass Filters under Laboratory and Field Conditions. C
C     Trans. ASAE, 27(5):1321-1331.                                        C
C ------------------------------------------------------------------------ C
C v1.05-12/2000 VFS wedge deposition on B(t) (US), based on App. 9C, Haan  C
C      et al. (1994), Prentice Hall.                                       C
C ------------------------------------------------------------------------ C
C     Known:                                                               C
C           dp: particle size, diameter (cm)                               C
C           Sc: filter main slope                                          C
C           gammaw, gammasb:  water and sediment weight density (g/cm3)    C
C           gacc= acceleration due to gravity (980 cm/s2)                  C
C           gsi: sediment load entering filter (g/s/cm)                    C
C           gsi_c: coarse sediment load fraction entering filter (g/s/cm)  C
C           gsi_f: fine sediment load fraction entering filter (g/s/cm)    C
C           gssi: sediment load entering B(t)(after upstream deposition)   C
C           gs2: sediment load into downstream section (D(t), gsd (g/s/cm) C
C           coarse: % of particles from incoming sediment with diameter    C
C                   > 0.0037 cm (coarse fraction that will be routed       C
C                   through wedge).                                        C
C     Unknown:                                                             C
C           Rss: hydraulic radius of the filter at B(t) (cm)               C
C           dfs: depth of flow at B(t) (cm)                                C
C           Vms: depth averaged velocity at B(t) (cm/s)                    C
C           Se: equilibrium slope at B(t)                                  C
C           f: fraccion trapped in the deposition wedge                    C
C           X1(t),Y(t),X2(t): sediment wedge geometry                      C
C           DEP: depth of deposited sediment at lower filter section       C
C           Totrap: Sediment trapping efficiency                           C
C                                                                          C
C     NOTE1: all units in CGS system (cm,g,s), including Manning's n       C
C     NOTE2: The QSED values are average values over the period, so the    C
C        integrals don't require to take the average value over interval   C
C                                                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/GRASSD/PART(3),SC,SS,VN1,VN2,VN,GSI,H,VLCM,POR,CI,ICO
      COMMON/GRASSD2/GSIMASS,GSOMASS,TTE,DEP,QSED(4),RS(3),DF(3),VM(3)
      COMMON/GRASSD3/SUSMASS,WEDGEMASS,NFUP
      COMMON/OLD/SEOLD,TOLD,XTOLD,YTOLD,CDEP,SE,VBTOLD,GSOLD

      FWID=FWIDTH*100.D0
      DTG=TIME-TOLD
      GAMMASB=(1.D0-POR)*PART(3)
      GAMMAW=1.D0
      GACC=981.D0
      GSI_C=COARSE*GSI
      GSI_F=(1.D0-COARSE)*GSI

c----(Ia)--If sediment transport capacity at the downstream side of sediment--
c-------wedge (gs2) is greater than the COARSE (diameter>0.0037cm) incoming---
c-------sediment load (GSI_C) all sediment goes through the wedge to the------
c-------lower part of the filter (ntrcap=1).----------------------------------
      IF(NTRCAP.EQ.1)THEN
            SE=SEOLD
            F=0.D0
            FI=0.D0
            XT=XTOLD
            YT=YTOLD
            X1=YT/SC
            GSSI=GSI
            GS2=GSI
            GS1=GSI
            GOTO 130
      ENDIF

c----(Ib)--If transport capacity is lower, then the FINE fraction goes through-
c-------the wedge and the COARSE (gsi_c) is filtered at the wedge ------------
      GS1=(GSI_C+GS2)/2.D0

C ---- Calculate df, Rs, Set on the center of zone B(t)
c ---- Newton-Raphson Method----------------------------
c------- start with US units to use algorithm in original units--
      DP=PART(1)*10.D0
      SSM=Ss/30.48D0
      GS1M=GS1*.067056D0
      QW1=QSED(1)/30.48D0**2.D0
      XK=6.462D07*PART(3)*(PART(3)-1.D0)**(-3.07D0)
      C1=VN*QW1*XK**.14D0/(1.5D0*GS1M**.14D0*DP**.2898D0)
      dfba=1.5D0
      MAXIT=500
      ERROR=1.D0
      NS=0
      DO 10 WHILE(NS.LT.MAXIT.AND.ERROR.GT.1D-8)
            dfbaO=dfba
            PHI= Ssm*dfba**7.D0 - C1**6.D0*(Ssm + 2.D0*dfba)
            DPHI=-2.D0*C1**6.D0 + 7.D0*Ssm*dfba**6.D0
            dfba=dfbaO-PHI/DPHI
            ERROR=DABS(dfbaO-dfba)
            NS=NS+1
10    CONTINUE
      RSm=dfba*Ssm/(2.D0*dfba+Ssm)
      Set=GS1M**.28D0*DP**.5798D0/(RSM*XK**.28D0)
c------- end of US units en conversion to CGS----
      dfba=dfba*30.48D0
      DF(1)=dfba
      RS(1)=dfba*Ss/(2.D0*dfba+Ss)
      VM(1)=QSED(1)/DF(1)
      Se=Set-Sc

C ---- Modify fraction of sediment trapped to account for upstream
c ---- deposition (After App. 9C in Haan et al, 1994)
      f=(GSI_C-GS2)/GSI_C
      alpha=(Set-Sc)/Sc
c ---- NOTE: there is a typo in Haan's reference: f multiplier to obtain fi
c ---- (fraction in lower end of wedge) must be 1/(1+alpha).
      fi=f*(1.D0/(1.D0+alpha))

C ---- Find wedge dimensions
      if(Ytold.lt.H) then
            FGSISE=fi*GSI_C*SE
            Yt=dsqrt(2.D0/gammasb*FGSISE*DTG+Ytold**2.D0)
        else
            fi=f
            Yt=H
      endif
      FGSISE=fi*GSI_C*SE
      if(Yt.ge.H) then
            Yt=H
            Xt=xtold+DTG*f*GSI_C/(H*gammasb)
        else
            Xt=dsqrt(2.D0/gammasb*fi*GSI_C/Se*DTG+xtold**2.D0)
      endif
      X1=Yt/Sc

C ---- Find upstream deposition (After App. 9C in Haan et al, 1994)
      if(Yt.lt.H) then
            gsu=gammasb*(Yt**2.D0-Ytold**2.D0)/(2.D0*DTG*Sc)
         else
            gsu=0.d0
      endif
      gssi=GSI_C-gsu
      fu=gsu/GSI_C

c ---- Calculate mass of sediment deposited and error in mass
c ---- from difference in sediment loads (GSI_C-GS2)
      SM2=(GSI_C-GS2)*DTG
      WEDGEMASS=WEDGEMASS+SM2
c----- v1.06-3/2002 recalculate sediment available at end of the wedge as coarse
c------amount transported through plus the fine fraction entering the wedge
      GS2=GS2+GSI_F

C ---- Check if top of vegetation is reached and also if strip has been
c ---- filled up and issue warnings -------------------------------
130   VLT=VLCM-Xt
      if(Yt.ge.H.and.Ytold.lt.H) WRITE(*,340)
      if(Xt.GE.VLCM) THEN
            WRITE(*,350)
            NFUP=1
            YT=H
            X1=YT/SC
            XT=VLCM
            VLT=0.D0
            FRAC=0.D0
            DEP=H
            CDEP=0.D0
            GOTO 140
      ENDIF

c---(II) Suspended sediment zone------------------------------------------
C=====> Tollner's (1977) Sediment Trapping Efficiency in lower section <==
c--------Note: Using Wilson's CDEP from last time step (DEP in inches)
      IF(TOLD.eq.0.D0)CDEP=1.D0
      IF (QSED(3).GT.0.D0)THEN
            RNUMBER=VM(3)*RS(3)/0.01D0
            PNF=(VLT*PART(2))/(VM(3)*DF(3))
            FRAC=CDEP*DEXP(-.00105D0*(RNUMBER**.82D0)*PNF**(-.91D0))
         ELSE
            FRAC=1.D0
      ENDIF

c-------Calculate subtotals for time step
140      GSO=GS2*(1.D0-FRAC)
c--05/08/03 If runoff outflow=0 force all sediment to deposit in lower part of the filter
      IF(QOUT.LE.0.D0)THEN
            FRAC=1.D0
            GSO=0.D0
      ENDIF
      GSOMASS=GSOMASS+0.5D0*(GSO+GSOLD)*DTG
      TOTRAP=(GSI-GSO)/GSI
      SUSMASS=(GS2-GSO)*DTG+SUSMASS

c-------On the assumption that the trapped sediment is uniformly distributed on
c-------bed of channel and the strip has not been filled up, calculate DEP,
c-------depth of sediment deposited for that dt and a multiplier for the next
c------- time step CDEP to reduce the actual GSO (Wilson et al. 1981).
      IF(NFUP.EQ.0) THEN
         GTRAP= GS2*FRAC
         GTRAPMASS=GTRAP*DTG
         VOLTRAP=GTRAPMASS/GAMMASB
         DEP=DEP+VOLTRAP/VLT
         DEPP=DEP/2.54D0
         CDEP=0.5D0*(DEXP(-3.D0*DEPP)+DEXP(15.D0*DEPP*(0.2D0-DEPP)))
      ENDIF

c-----Write output to .og1 file----
      IF(ISCR.EQ.0) THEN
         WRITE(13,80)TIME,YT,X1,XT,VLT,SE,GSI,GSSI,GS2,GSO,GSIMASS,
     &      WEDGEMASS,SUSMASS,GSOMASS,FI,FRAC,DEP,CDEP,TOTRAP
       ELSE
         WRITE(13,80)TIME,YT,X1,XT,VLT,SE,GSI*FWID,GSSI*FWID,GS2*FWID,
     &      GSO*FWID,GSIMASS*FWID,WEDGEMASS*FWID,SUSMASS*FWID,
     &      GSOMASS*FWID,FI,FRAC,DEP,CDEP,TOTRAP
      ENDIF

c-----Update values for next time step-------------
      SEOLD=SE
      GSOLD=GSO
      TOLD=TIME
      XTOLD=XT
      YTOLD=YT

c ---- Formats ---------
80    FORMAT(f7.0,4F10.3,F10.6,13E10.3)
340   format(/,66('*'),/,'*  WARNING: Top of vegetation reached',
     &     ' - trapezoidal wedge starts *',/,66('*'))
350   format(/,66('*'),/,'*  WARNING: Strip filled up!',37x,'*',
     &     /,66('*'))

      return
      end
