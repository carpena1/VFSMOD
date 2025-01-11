      SUBROUTINE WQPEST(JJ,VIN,VOUT,SMIN,SMOUT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C  Pesticide trapping, mass balance and degradation component                 C
C  Variables used:                                                            C
C  SAREA,VFSAREA (m2): Source (field) and filter areas                        C
C  VIN,VOUT(m3)= VFS runoff inflow and outflow volume hydrograph areas        C
C  SMIN,SMOUT(Kg)= VFS sediment inflow and outflow                            C
C  WAT_IN(m3)= Runoff inflow (VIN) + total rain on the filter for the event   C
C  VF(m3)= Total Infiltration in Filter                                       C
C  OS= soil saturated water content (porosity)                                C
C  TTE [-]= sediment reduction in filter                                      C
C  PDSED(%) = sediment reduction in filter (dE= TTE*100)                      C
C  DELTAP (%)= pesticide reduction in filter (dP)                             C
C  PDQ(%) = flow reduction in filter (dQ)                                     C
C  DGPIN(JJ)(mg/m2) = incoming pesticide in filter (per m2 of SOURCE area)    C
C  DGLD(JJ)(m)= ambda, dispersion length of chemical (m). This can be taken   C
C           as 0.05m from FOCUS-Pearl (Default)                               C
C  DGMfFd,DGMfFp,DGMfF(mg)= leached pesticide (dissolved, adsorbed, total)    C
C  DGMmld,DGMmlp,DGMml(mg)= mixing layer pest. (dissolved, adsorbed, total)   C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION   (A-H,O-Z)
      COMMON/PAR1/VL,FWIDTH,SWIDTH,SLENGTH
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      COMMON/GRASSD2/GSIMASS,GSOMASS,TTE,DEP,QSED(4),RS(3),DF(3),VM(3)
      COMMON/WQ1/VKD(10),CCP,CSAB(5),DGMRES0(10),DGMOL(10),DGFRAC(10,10)
      COMMON/IWQ2/NDGDAY,IDG,IWQ,IWQPRO,ICAT,IMOB
      COMMON/WQ3/DGKREF(10),FC,DGPIN(10),DGML,DGT(366),DGTHETA(366),DGLD(10),RF
      COMMON/CDE1/DGMfFd(10),DGMfFp(10),DGMfF(10),DGMmld(10),DGMmlp(10),DGMml(10)

      SAREA=SLENGTH*SWIDTH
      VFSAREA=VL*FWIDTH
      VF=F*VFSAREA
      TOTRAIN=PS*VFSAREA
      WAT_IN=Vin+TOTRAIN
c---- convert from m3 to L ----
      VINL=Vin*1000.D0
      TOTRAINL=TOTRAIN*1000.d0
      WAT_INL=WAT_IN*1000.d0
c---- calculate factors for pesticide trapping equations ----
      IF (SMIN.EQ.0.D0) THEN
          FPH=WAT_INL/(VKD(JJ)*0.001D0)
        ELSE
          FPH=WAT_INL/(VKD(JJ)*SMIN)
      ENDIF
      IF(WAT_IN.gt.0.d0.AND.VOUT.GT.0.d0)THEN
          PDQ=VF/WAT_IN*100.D0
        ELSE
          PDQ=100.D0
      ENDIF
      PDSED=TTE*100.D0

C-------Calculate pesticide trapping efficiency based on:--------------
C-------IWQPRO=1 - Sabbagh et al. (2009) semiempirical eq. ---------------
C-------IWQPRO=2 - Refitted Sabbagh et al. eq. w/user supplies coeff -----
C-------IWQPRO=3 - Based on Muñoz-Carpena et al. (2015) mechanistic eq. --
C-------IWQPRO=4 - Chen et al. (2017) empirical eq. ----------------------
      SELECT CASE (IWQPRO)
        CASE (1)
c-------Orifinal Sabbagh eq.
           DELTAP=CSAB(1)+CSAB(2)*PDQ+CSAB(3)*PDSED+CSAB(4)*
     &          DLOG(FPH+1.D0)+CSAB(5)*CCP
         CASE (2)
c--------Sabbagh refit eq.
           DELTAP=CSAB(1)+CSAB(2)*PDQ+CSAB(3)*PDSED+CSAB(4)*
     &          DLOG(FPH+1.D0)+CSAB(5)*CCP
         CASE (3)
c--------Mechanistic mass balance eq.
           DELTAP=100.D0*MIN(VINL+VKD(JJ)*SMIN,PDSED/100*SMIN*VKD(JJ)+
     &          PDQ/100*VINL)/(VINL+VKD(JJ)*SMIN)
         CASE DEFAULT
c--------Chen eq. (California)
           DELTAP=101.D0-(8.06D0-0.07D0*PDQ+0.02D0*PDSED+0.05D0*CCP-
     &          2.17D0*ICAT+0.02D0*PDQ*ICAT-0.0003D0*PDQ*PDSED)**2.D0
      END SELECT
      IF (DELTAP.GT.100.D0.OR.PDQ.GE.100.d0.OR.VOUT.LE.0.d0) DELTAP=100.D0
      IF (PDQ.EQ.0.d0.AND.PDSED.EQ.0.d0) DELTAP=0.D0
      IF (DELTAP.LT.0.d0) DELTAP=0.D0

      WRITE(18,*)
      WRITE(18,*)'Outputs for Water Quality'
      WRITE(18,3500)
      WRITE(18,4005)VIN
      WRITE(18,4105)SMIN
      IF(dabs(FPH).le.9999.d0) THEN
          WRITE(18,4200)FPH
        ELSE
          WRITE(18,4205)FPH
      ENDIF
      WRITE(18,4300)PDQ
      WRITE(18,4400)PDSED
      IF(Vin.gt.0.d0) THEN
          RRED=100.D0*(1-Vout/Vin)
        ELSE
          RRED=100.d0
      ENDIF
      IF(dabs(RRED).gt.200.d0) THEN
          WRITE(18,4455)RRED
        ELSE
          WRITE(18,4450)RRED
      ENDIF
      WRITE(18,*)
      WRITE(18,4500)DELTAP

C-------Calculate pesticide mass balance and degradation (IDG=1-4)-----
C-------based on Muñoz-Carpena et al. (2015) --------------------------
      DGROB=(1.d0-OS)*2.65d0
      OI=OS-DM
c---------Carry-over sorbed in mixing layer from previous event--------
      IF(IMOB.EQ.2) then
          DGMRESP0=0
        ELSE
          DGMRESP0=DGMRES0(JJ)*VKD(JJ)*DGROB/OI*SAREA
      ENDIF
c---------IN/OUT/SED/MIXING LAYER fractions----------------------------
      DGMI=DGPIN(JJ)*SAREA
      DGMO=DGMI*(1.d0-DELTAP/100.d0)
      IF(VIN.le.0.d0.and.SMIN.le.0.d0) THEN
          DGSI=0.d0
        ELSE
          DGSI=DGMI*VKD(JJ)/(VINL+SMIN*VKD(JJ))
      ENDIF
      DGMF=DGMI*DELTAP/100.d0
      DGMFSED=DGSI*SMIN*PDSED/100.d0
c-----For high sorption pesticides check if all pesticide is sediment-bonded
      IF(DGMFSED.GT.DGMF) DGMFSED=DGMF
      DGMFF=DGMF-DGMFSED
      DGMRES=DGMFSED+DGMml+DGMRESP0
      DGMRESFSED=DGMFSED
      DGMRESMML=DGMml
      DGTHETAN=OS

c-----Scheme for 4 degradation types (IDG=1: EU-FOCUS; 2: US-EPA K=Kref;
c-----3: k=k(T); 4:k=k(theta))
      IF(IDG.GT.0.AND.IDG.LE.4.AND.NDGDAY.GT.0) THEN
            DGTHETAN=DGTHETA(NDGDAY)
          DGBETA=-0.7d0
          DO 40 I=1,NDGDAY
              IF(IDG.EQ.1) THEN
c------------------For EU FOCUS reference temperature is 20 C, QT10=2.58,
c------------------EFSA Journal (2007):622:1-32
                   DGC1=65.4d0/(0.008314d0)
                   DGC2=1.d0/293.15d0
                   DGKTEMP=DEXP(DGC1*(DGC2-(1.d0/(DGT(I)+273.15d0))))
                   DGKTHETA=(DGTHETA(I)/FC)**DGBETA
                ELSEIF(IDG.EQ.2) THEN
                   DGKTEMP=1.D0
                   DGKTHETA=1.D0
                ELSEIF(IDG.EQ.3) THEN
c------------------For US EPA the reference temperature is 25 C, , QT10=2.0,
c----------------- EPA/PMRA PRZM Guidance v.1.0 (2012) (Appendix A, page 5)
                   DGC1=49.5d0/(8.314d0/1000.d0)
                   DGC2=1.d0/298.15d0
                   DGKTEMP=DEXP(DGC1*(DGC2-(1.d0/(DGT(I)+273.15d0))))
                   DGKTHETA=1.D0
                ELSEIF(IDG.EQ.4) THEN
                   DGKTEMP=1.D0
                   DGKTHETA=(DGTHETA(I)/FC)**DGBETA
               ENDIF
               DGKI=DGKREF(JJ)*DGKTEMP*DGKTHETA
               DGMRES=DGMRES*DEXP(-DGKI)
c              check for NaN
c old fortran  IF(ISNAN(DGMRES)) DGMRES=0.d0
               IF(DGMRES.NE.DGMRES) DGMRES=0.d0
c               DGMRESFSED=DGMRESFSED*DEXP(-DGKI)
c               DGMRESMML=DGMRESMML*DEXP(-DGKI)
40        CONTINUE
      ENDIF

c-----Solid/liquid phase partitioning between outgoing sediment (SMOUT)
c-----and water (VOUT)assuming linear equilibrium based on Kd (VKD(JJ)).
      IF(VOUT.GT.0.D0) THEN
            DGMOP=(DGMO*SMOUT*VKD(JJ))/(VOUT*1000.D0+SMOUT*VKD(JJ))
            DGMOD=DGMO-DGMOP
         ELSE
            DGMOP=0.D0
            DGMOD=0.D0
      ENDIF

c-----Calculate surface residue mass in the porewater based sorption
c-----equilibrium for the soil water content at the next event
      Vw=dgML*VFSAREA*DGTHETAN*10.d0
      DGMsoil=dgML*VFSAREA*DGROB*10.d0
      DGMRESD=DGMRES*Vw/(Vw+DGMsoil*VKD(JJ))
      DGMRESP=DGMRESD*VKD(JJ)*DGMsoil/Vw

c-----Remobilization (IMOB) options. Calculate
      SELECT CASE (IMOB)
c-------Full remobilization of surface residue, m'i= mml_res+mf,sed_res (IMOB=2)
        CASE (2)
          DGMRESI=DGMRES
c-------No remobilization of surface residue, m'i= 0 (IMOB=3)
        CASE (3)
          DGMRESI=0
c-------Partial remobilization of surface residue, m'i= mml_res (IMOB=1)
        CASE DEFAULT
          DGMRESI=DGMRESD
      END SELECT

C-----Write pesticide trapping outputs in OWQ file--------------
      ZEROPRT=1D-99
      WRITE(18,*)
      WRITE(18,3525)
      WRITE(18,3550)
      IF(DGMI.LT.ZEROPRT)DGMI=0.d0
      write(18,4600)DGMI
      IF(DGMO.LT.ZEROPRT)DGMO=0.d0
      write(18,4700)DGMO
      IF(DGMOP.LT.ZEROPRT)DGMOP=0.d0
      write(18,5800)DGMOP
      IF(DGMOD.LT.ZEROPRT)DGMOD=0.d0
      write(18,5900)DGMOD
      IF(DGMF.LT.ZEROPRT)DGMF=0.d0
      write(18,4800)DGMF
      IF(DGMFSED.LT.ZEROPRT)DGMFSED=0.d0
      write(18,4805)DGMFSED
      IF(DGMml.LT.ZEROPRT)DGMml=0.d0
      write(18,4900)DGMml
      IF(DGMRESP0.LT.ZEROPRT)DGMRESP0=0.d0
      write(18,4950)DGMRESP0
      write(18,5000)DGMml+DGMFSED+DGMRESP0
      IF(DGMRES.LT.ZEROPRT)DGMRES=0.d0
      write(18,5100)DGMRES,NDGDAY
      IF(DGMRESD.LT.ZEROPRT)DGMRESD=0.d0
      write(18,5105)DGMRESD,NDGDAY
      IF(DGMRESP.LT.ZEROPRT)DGMRESP=0.d0
      write(18,5110)DGMRESP,NDGDAY
      IF(DGMRESI.LT.ZEROPRT)DGMRESI=0.d0
      WRITE(18,6001)DGMRESI,IMOB
      WRITE(18,*)
      WRITE(18,*)'Normalized values by source area:'
      WRITE(18,*)
      WRITE(18,5150)SAREA
      write(18,5200)DGPIN(JJ)
      write(18,5300)DGMO/SAREA
      write(18,5850)DGMOP/SAREA
      write(18,5950)DGMOD/SAREA
      write(18,5400)DGMF/SAREA
      write(18,5405)DGMFSED/SAREA
      write(18,5500)DGMml/SAREA
      write(18,5550)DGMRESP0/SAREA
      write(18,5600)(DGMml+DGMFSED+DGMRESP0)/SAREA
      write(18,5700)DGMRES/SAREA,NDGDAY
      write(18,5705)DGMRESD/SAREA,NDGDAY
      write(18,5710)DGMRESP/SAREA,NDGDAY
      WRITE(18,6002)DGMRESI/SAREA,IMOB

3500  FORMAT(1X,26('-'))
3525  FORMAT(1X,'Pesticide mass balance, degradation & remobilization')
3550  FORMAT(1X,43('-'))
4005  FORMAT(E10.3,' m3 = Runoff inflow')
4105  FORMAT(E10.3,' Kg = Sediment inflow')
4200  FORMAT(F10.3,'    = Phase distribution, Fph')
4205  FORMAT(E10.3,'    = Phase distribution, Fph')
4300  FORMAT(F10.3,' %  = Infiltration (dQ)')
4400  FORMAT(F10.3,' %  = Sediment reduction (dE)')
4450  FORMAT(F10.3,' %  = Runoff inflow reduction')
4455  FORMAT(E10.3,' %  = Runoff inflow reduction')
4500  FORMAT(F10.3,' %  = Pesticide reduction (dP)')
4600  FORMAT(E15.6,' mg   = Pesticide input (mi)')
4700  FORMAT(E15.6,' mg   = Pesticide output (mo)')
4800  FORMAT(E15.6,' mg   = Pesticide trapped in VFS (mf)')
4805  FORMAT(E15.6,' mg   = Pesticide trapped with sediment (mfsed)')
4900  FORMAT(E15.6,' mg   = Pesticide trapped in mixing layer (mfml)')
4950  FORMAT(E15.6,' mg   = Pesticide in mixing layer from last event',
     &                 ' (mfml0)')
5000  FORMAT(E15.6,' mg   = Total surface residue ',
     &                 '(mres1=mfml+mfsed+mfml0)')
5100  FORMAT(E15.6,' mg   = Total surface residue after ',
     &                 'degradation (',I3,' days)')
5105  FORMAT(E15.6,' mg   = Dissolved surface residue after ',
     &                 'degradation (',I3,' days)')
5110  FORMAT(E15.6,' mg   = Sorbed surface residue after ',
     &                 'degradation (',I3,' days)')

5150  FORMAT(F15.2,' m^2  = Source Area (input)')
5200  FORMAT(E15.6,' mg/m2= Pesticide input (mi)')
5300  FORMAT(E15.6,' mg/m2= Pesticide output (mo)')
5400  FORMAT(E15.6,' mg/m2= Pesticide trapped in VFS (mf)')
5405  FORMAT(E15.6,' mg/m2= Pesticide trapped with sediment (mfsed)')
5500  FORMAT(E15.6,' mg/m2= Pesticide trapped in mixing layer (mfml)')
5550  FORMAT(E15.6,' mg/m2= Pesticide in mixing layer from last event',
     &                 ' (mfml0)')
5600  FORMAT(E15.6,' mg/m2= Total surface residue ',
     &                 '(mres1=mfml+mfsed+mfml0)')
5700  FORMAT(E15.6,' mg/m2= Total residue after ',
     &                 'degradation (',I3,' days)')
5705  FORMAT(E15.6,' mg/m2= Dissolved surface residue after ',
     &                 'degradation (',I3,' days)')
5710  FORMAT(E15.6,' mg/m2= Sorbed surface residue after ',
     &                 'degradation (',I3,' days)')
5800  FORMAT(E15.6,' mg   = Pesticide outflow in solid phase (mop)')
5850  FORMAT(E15.6,' mg/m2= Pesticide outflow in solid phase (mop)')
5900  FORMAT(E15.6,' mg   = Pesticide outflow in liquid phase (mod)')
5950  FORMAT(E15.6,' mg/m2= Pesticide outflow in liquid phase (mod)')
6001  FORMAT(E15.6,' mg   = Next event residue remobilization',
     & ' (mresn, IMOB=',I1,')')
6002  FORMAT(E15.6,' mg/m2= Next event residue remobilization',
     & ' (mresn, IMOB=',I1,')')

      RETURN
      END
