      SUBROUTINE WQPEST(VIN,VOUT,SMIN,SMOUT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                              C
C  Pesticide trapping, mass balance and degradation component                  C
C  Variables used:                                                             C 
C  SAREA,VFSAREA (m2): Source (field) and filter areas                         C
C  VIN,VOUT(m3)= VFS runoff inflow and outflow volume hydrograph areas         C
C  SMIN,SMOUT(Kg)= VFS sediment inflow and outflow                             C
C  WAT_IN(m3)= Runoff inflow (VIN) + total rain on the filter for the event    C
C  VF(m3)= Total Infiltration in Filter                                        C
C  OS= soil saturated water content (porosity)                                 C
C  TTE [-]= sediment reduction in filter                                       C
C  PDSED(%) = sediment reduction in filter (dE= TTE*100)                       C
C  DELTAP (%)= pesticide reduction in filter (dP)                              C
C  PDQ(%) = flow reduction in filter (dQ)                                      C
C  DGPIN(JJ)(mg/m2) = incoming pesticide in filter (per m2 of SOURCE area)     C
C  DGLD(JJ)(m)= ambda, dispersion length of chemical (m). This can be taken    C
C           as 0.05m from FOCUS-Pearl (Default)                                C
C  DGMfFd,DGMfFp,DGMfF(mg)= leached pesticide (dissolved, adsorbed, total)     C
C  DGMmld,DGMmlp,DGMml(mg)= mixing layer pest. (dissolved, adsorbed, total)    C
C  DGMRES(mg)= pesticide surface residues at the next storm (after degradation)C
C                                                                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION   (A-H,O-Z)
      COMMON/PAR1/VL,FWIDTH,SWIDTH,SLENGTH
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      COMMON/GRASSD2/GSIMASS,GSOMASS,TTE,DEP,QSED(4),RS(3),DF(3),VM(3)
      COMMON/WQ1/VKD(10),CCP,CSAB(5),DGMRES0(10),DGMOL(10),DGFRAC(10,10)
      COMMON/IWQ2/NDGDAY,IDG,IWQ,IWQPRO,ICAT,IMOB
      COMMON/WQ3/DGKREF(10),FC,DGPIN(10),DGML,DGT(366),DGTHETA(366),DGLD(10),RF
      COMMON/CDE1/DGMfFd(10),DGMfFp(10),DGMfF(10),DGMmld(10),DGMmlp(10),DGMml(10)
      
      DIMENSION DGMRES0P(10),DGMI(10),DGMO(10),DGMOD(10),DGMOP(10),DGMf(10),
     &      DGMfSED(10),DGMRES(10),DGMRESFSED(10),DGMRESMML(10),FPH(10),
     &      DGMRESP(10),DGMRESD(10), DGMRESI(10),DGMSurf(10)
      DIMENSION DGMRES0PAREA(10),DGMIAREA(10),DGMOAREA(10),DGMODAREA(10),
     &      DGMOPAREA(10),DGMfAREA(10), DGMfSEDAREA(10),DGMRESAREA(10),
     &      DGMRESFSEDAREA(10),DGMRESMMLAREA(10),DGMRESPAREA(10),
     &      DGMRESDAREA(10), DGMRESIAREA(10),DGMSurfAREA(10),DGMmlAREA(10)
      CHARACTER*200 LINE
      CHARACTER*3 SNDGDAY 
      CHARACTER*1 SIMOB 
      CHARACTER*58 DESCR1(13)
      DATA(DESCR1(I),I=1,13)/'Pesticide input (mi)',
     & 'Pesticide output (mo)',
     & 'Pesticide outflow in solid phase (mop)',
     & 'Pesticide outflow in liquid phase (mod)',
     & 'Pesticide trapped in VFS (mf)',
     & 'Pesticide trapped with sediment (mfsed)',
     & 'Pesticide trapped in mixing layer (mfml)',
     & 'Pesticide sorbed in mixing layer from last event (mfml0)',
     & 'Total surface residue (mres1=mfml+mfsed+mfml0)',
     & 'Total surface residue after degradation (',
     & 'Dissolved surface residue after degradation (',
     & 'Sorbed surface residue after degradation (',
     & 'Next event residue remobilization (mresn, IMOB='/

c---- Water in VFS
      SAREA=SLENGTH*SWIDTH
      VFSAREA=VL*FWIDTH
      VF=F*VFSAREA
      TOTRAIN=PS*VFSAREA
      WAT_IN=Vin+TOTRAIN
c---- convert from m3 to L ----
      VINL=Vin*1000.D0
      TOTRAINL=TOTRAIN*1000.d0
      WAT_INL=WAT_IN*1000.d0
c-----Auxiliary values for pesticide redistribution---
      DGROB=(1.d0-OS)*2.65d0
      OI=OS-DM
      DGTHETAN=DGTHETA(NDGDAY)
c---- Calculate flow and sediment reduction factors for pesticide trapping equations ----
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
      WRITE(18,3400)'Outputs for Water Quality'
      DO JJ=1,IWQ
            IF (SMIN.EQ.0.D0) THEN
                FPH(JJ)=WAT_INL/(VKD(JJ)*0.001D0)
              ELSE
                FPH(JJ)=WAT_INL/(VKD(JJ)*SMIN)
            ENDIF
            SELECT CASE (IWQPRO)
              CASE (1)
c-------------Original Sabbagh eq.
                DELTAP=CSAB(1)+CSAB(2)*PDQ+CSAB(3)*PDSED+CSAB(4)*
     &          DLOG(FPH(JJ)+1.D0)+CSAB(5)*CCP
              CASE (2)
c-------------Sabbagh refit eq.
                DELTAP=CSAB(1)+CSAB(2)*PDQ+CSAB(3)*PDSED+CSAB(4)*
     &          DLOG(FPH(JJ)+1.D0)+CSAB(5)*CCP
              CASE (3)
c-------------Mechanistic mass balance eq. - (DEFAULT)
                DELTAP=100.D0*MIN(VINL+VKD(JJ)*SMIN,PDSED/100*SMIN*VKD(JJ)+
     &          PDQ/100*VINL)/(VINL+VKD(JJ)*SMIN)
              CASE DEFAULT
c-------------Chen eq. (California)
                DELTAP=101.D0-(8.06D0-0.07D0*PDQ+0.02D0*PDSED+0.05D0*CCP-
     &          2.17D0*ICAT+0.02D0*PDQ*ICAT-0.0003D0*PDQ*PDSED)**2.D0
            END SELECT
            IF (DELTAP.GT.100.D0.OR.PDQ.GE.100.d0.OR.VOUT.LE.0.d0) DELTAP=100.D0
            IF (PDQ.EQ.0.d0.AND.PDSED.EQ.0.d0) DELTAP=0.D0
            IF (DELTAP.LT.0.d0) DELTAP=0.D0
            IF (Vin.EQ.0.d0) DELTAP=100.D0
C-------Print VFS efficiencies in .iwq file --------------
            IF(IWQ.GT.1) WRITE(18,3500)'COMPOUND',JJ,':'
            WRITE(18,4005)VIN
            WRITE(18,4105)SMIN
            IF(dabs(FPH(JJ)).le.9999.d0) THEN
              WRITE(18,4200)FPH(JJ)
             ELSE
              WRITE(18,4205)FPH(JJ)
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
c---------Carry-over sorbed in mixing layer from previous event--------
            IF(IMOB.EQ.2) then
              DGMRES0P(JJ)=0
             ELSE
              DGMRES0P(JJ)=DGMRES0(JJ)*VKD(JJ)*SAREA*DGROB/OI
            ENDIF
c---------IN/OUT/SED/MIXING LAYER fractions----------------------------
            DGMI(JJ)=DGPIN(JJ)*SAREA
            DGMO(JJ)=DGMI(JJ)*(1.d0-DELTAP/100.d0)
            IF(VIN.le.0.d0.and.SMIN.le.0.d0) THEN
              DGSI=0.d0
             ELSE
              DGSI=DGMI(JJ)*VKD(JJ)/(VINL+SMIN*VKD(JJ))
             ENDIF
            DGMf(JJ)=DGMI(JJ)*DELTAP/100.d0
            DGMfSED(JJ)=DGSI*SMIN*PDSED/100.d0
c---------Check approximation error on sediment bonded pesticide nad issue warning if >1%
            DGMfF1=DGMf(JJ)-DGMfSED(JJ)-DGMfF(JJ)
            IF(DABS(DGMfF1/DGMf(JJ)).gt.0.01) THEN
c              print*,'Warning: DGMfSED estimation error >1%)'
              DGMfSED(JJ)=DGMfSED(JJ)+DGMfF1
            ENDIF
c---------Pesticide residues on the VFS surface at the end of the runoff event
            DGMRES(JJ)=DGMfSED(JJ)+DGMml(JJ)+DGMRES0P(JJ)
            DGMRESFSED(JJ)=DGMfSED(JJ)
            DGMRESMML(JJ)=DGMml(JJ)
      END DO
c-----Pesticide degradation in days between runoff events for single compound 
c-----and multiple species reactions (parent-metabolites)     
      CALL multspcalc(DGMRES)   

c---- Final mass balance and remobilization
      Vw=dgML*VFSAREA*DGTHETAN*10.d0
      DGMsoil=dgML*VFSAREA*DGROB*10.d0
      DO JJ=1,IWQ
c-------Solid/liquid phase partitioning between outgoing sediment (SMOUT)
c-------and water (VOUT)assuming linear equilibrium based on Kd (VKD(JJ)).
        IF(VOUT.GT.0.D0) THEN
            DGMOP(JJ)=(DGMO(JJ)*SMOUT*VKD(JJ))/(VOUT*1000.D0+SMOUT*VKD(JJ))
            DGMOD(JJ)=DGMO(JJ)-DGMOP(JJ)
         ELSE
            DGMOP(JJ)=0.D0
            DGMOD(JJ)=0.D0
        ENDIF
c-----Calculate surface residue mass in the porewater based sorption
c-----equilibrium for the soil water content at the next event
        DGMRESD(JJ)=DGMRES(JJ)*Vw/(Vw+DGMsoil*VKD(JJ))
        DGMRESP(JJ)=DGMRESD(JJ)*VKD(JJ)*DGMsoil/Vw

c-----Remobilization (IMOB) options. Calculate
        SELECT CASE (IMOB)
c-------Full remobilization of surface residue, m'i= mml_res+mf,sed_res (IMOB=2)
          CASE (2)
            DGMRESI(JJ)=DGMRES(JJ)
c-------No remobilization of surface residue, m'i= 0 (IMOB=3)
          CASE (3)
            DGMRESI(JJ)=0
c-------Partial remobilization of surface residue, m'i= mml_res (IMOB=1)
          CASE DEFAULT
            DGMRESI(JJ)=DGMRESD(JJ)
        END SELECT
      END DO

c -----Write mass balance for all compounds in .owq file
      WRITE(18,*)
      WRITE(18,3525)
      ZEROPRT=1D-99
      WRITE (SNDGDAY, '(I3)') NDGDAY
      WRITE (SIMOB, '(I1)') IMOB
c------normalize values by area
      DO JJ=1,IWQ
            IF(DGMI(JJ).LT.ZEROPRT)DGMI(JJ)=0.d0
            DGMIAREA(JJ)=DGMI(JJ)/SAREA
            IF(DGMO(JJ).LT.ZEROPRT)DGMO(JJ)=0.d0
            DGMOAREA(JJ)=DGMO(JJ)/SAREA
            IF(DGMOP(JJ).LT.ZEROPRT)DGMOP(JJ)=0.d0
            DGMOPAREA(JJ)=DGMOP(JJ)/SAREA
            IF(DGMOD(JJ).LT.ZEROPRT)DGMOD(JJ)=0.d0
            DGMODAREA(JJ)=DGMOD(JJ)/SAREA
            IF(DGMf(JJ).LT.ZEROPRT)DGMf(JJ)=0.d0
            DGMfAREA(JJ)=DGMf(JJ)/SAREA
            IF(DGMfSED(JJ).LT.ZEROPRT)DGMfSED(JJ)=0.d0
            DGMfSEDAREA(JJ)=DGMfSED(JJ)/SAREA
            IF(DGMml(JJ).LT.ZEROPRT)DGMml(JJ)=0.d0
            DGMmlAREA(JJ)=DGMml(JJ)/SAREA
            IF(DGMRES0P(JJ).LT.ZEROPRT)DGMRES0P(JJ)=0.d0
            DGMRES0PAREA(JJ)=DGMRES0P(JJ)/SAREA
            DGMSurf(JJ)=DGMml(JJ)+DGMfSED(JJ)+DGMRES0P(JJ)
            DGMSurfAREA(JJ)=DGMSurf(JJ)/SAREA
            IF(DGMRES(JJ).LT.ZEROPRT)DGMRES(JJ)=0.d0
            DGMRESAREA(JJ)=DGMRES(JJ)/SAREA
            IF(DGMRESD(JJ).LT.ZEROPRT)DGMRESD(JJ)=0.d0
            DGMRESDAREA(JJ)=DGMRESD(JJ)/SAREA
            IF(DGMRESP(JJ).LT.ZEROPRT)DGMRESP(JJ)=0.d0
            DGMRESPAREA(JJ)=DGMRESP(JJ)/SAREA
            IF(DGMRESI(JJ).LT.ZEROPRT)DGMRESI(JJ)=0.d0
            DGMRESIAREA(JJ)=DGMRESI(JJ)/SAREA
          END DO
c-------Print mass balance for all compounds in .owq file
      IF(IWQ.GT.1) WRITE(18,3535)('Compound',JJ,JJ=1,IWQ)
      CALL COMP_STRING(0,DGMI,DESCR1(1))
      CALL COMP_STRING(0,DGMO,DESCR1(2))
      CALL COMP_STRING(0,DGMOP,DESCR1(3))
      CALL COMP_STRING(0,DGMOD,DESCR1(4))
      CALL COMP_STRING(0,DGMf,DESCR1(5))
      CALL COMP_STRING(0,DGMfSED,DESCR1(6))
      CALL COMP_STRING(0,DGMml,DESCR1(7))
      CALL COMP_STRING(0,DGMRES0P,DESCR1(8))
      CALL COMP_STRING(0,DGMSurf,DESCR1(9))
      CALL COMP_STRING(1,DGMRES,DESCR1(10))
      CALL COMP_STRING(1,DGMRESD,DESCR1(11))
      CALL COMP_STRING(1,DGMRESP,DESCR1(12))
      CALL COMP_STRING(2,DGMRESI,DESCR1(13))
      WRITE(18,*)
      WRITE(18,*)'Normalized values by source area:'
      WRITE(18,*)
      WRITE(18,5150)SAREA
      CALL COMP_STRING(3,DGMIAREA,DESCR1(1))
      CALL COMP_STRING(3,DGMOAREA,DESCR1(2))
      CALL COMP_STRING(3,DGMOPAREA,DESCR1(3))
      CALL COMP_STRING(3,DGMODAREA,DESCR1(4))
      CALL COMP_STRING(3,DGMfAREA,DESCR1(5))
      CALL COMP_STRING(3,DGMfSEDAREA,DESCR1(6))
      CALL COMP_STRING(3,DGMmlAREA,DESCR1(7))
      CALL COMP_STRING(3,DGMRES0PAREA,DESCR1(8))
      CALL COMP_STRING(3,DGMSurfAREA,DESCR1(9))
      CALL COMP_STRING(4,DGMRESAREA,DESCR1(10))
      CALL COMP_STRING(4,DGMRESDAREA,DESCR1(11))
      CALL COMP_STRING(4,DGMRESPAREA,DESCR1(12))
      CALL COMP_STRING(5,DGMRESIAREA,DESCR1(13))
      WRITE(18,*)

3400  FORMAT(/,43('-'),/,A26,/,43('-'))
3500  FORMAT(/,A8,I3,A1)
3525  FORMAT(100('-'),/,1x,'Pesticide mass balance, degradation & remobilization',/,100('-'))
3535   FORMAT(10(3x,A10,I2))
4005  FORMAT(E10.3,' m3 = Runoff inflow')
4105  FORMAT(E10.3,' Kg = Sediment inflow')
4200  FORMAT(F10.3,'    = Phase distribution, Fph')
4205  FORMAT(E10.3,'    = Phase distribution, Fph')
4300  FORMAT(F10.3,' %  = Infiltration (dQ)')
4400  FORMAT(F10.3,' %  = Sediment reduction (dE)')
4450  FORMAT(F10.3,' %  = Runoff inflow reduction')
4455  FORMAT(E10.3,' %  = Runoff inflow reduction')
4500  FORMAT(F10.3,' %  = Pesticide reduction (dP)')
5150  FORMAT(F15.2,' m^2  = Source Area (input)')

      RETURN
      END


      SUBROUTINE COMP_STRING(IFMT,DUM,DESCR)
c-----prepare lines for mass balance with multiple compounds in a single line
c-----DUM(JJ): array containing values for each compound (JJ=1,WQ)
c-----DESCR: Units and descriptor for value
c-----IFMT: 0 for outputs with no days or imob variable (mg), 1,2 for NDGDAY 
c-----      and IMOB (mg), 3 for outputs with no NDGDAYS or IMOB (mg/m2), 
c-----      4,5 for NDGDAY and IMOB (mg/m2), 

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/IWQ2/NDGDAY,IDG,IWQ,IWQPRO,ICAT,IMOB
      DIMENSION DUM(10)
      CHARACTER*200 LINE
      CHARACTER*3 SNDGDAY 
      CHARACTER*1 SIMOB 
      CHARACTER*58 DESCR     
      
      LINE = ' '
      DO JJ = 1, IWQ
            WRITE(LINE, '(A,E15.6,E15.6)') TRIM(LINE),DUM(JJ)
      END DO

      SELECT CASE (IFMT)
        CASE (1)
          WRITE (SNDGDAY, '(I3)') NDGDAY
          LINE=TRIM(LINE) // ' mg   = ' // TRIM(DESCR) // SNDGDAY // ' days)'
        CASE (2)
          WRITE (SIMOB, '(I1)') IMOB
          LINE=TRIM(LINE)  // ' mg   = ' // TRIM(DESCR) // SIMOB // ')'
        CASE (3)
          LINE=TRIM(LINE)  //' mg/m2= ' // TRIM(DESCR)
        CASE (4)
          WRITE (SNDGDAY, '(I3)') NDGDAY
          LINE=TRIM(LINE) // ' mg/m2= ' // TRIM(DESCR) // SNDGDAY // ' days)'
        CASE (5)
          WRITE (SIMOB, '(I1)') IMOB
          LINE=TRIM(LINE)  // ' mg/m2= ' // TRIM(DESCR) // SIMOB // ')'
        CASE DEFAULT
          LINE=TRIM(LINE)  //' mg   = ' // TRIM(DESCR)
      END SELECT

      LINE=TRIM(LINE)
      WRITE(18,'(A)'), LINE

      return
      END
      