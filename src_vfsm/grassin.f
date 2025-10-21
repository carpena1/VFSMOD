      SUBROUTINE GRASSIN(ICOARSE,COARSE,LISFIL,QSUM0,RMMH,BCROPEAK,
     &    ISCR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C     Get filter properties to be used in later from filename.igr             C
C     NOTE: all units in CGS system (cm,g,s), including modified Manning's n  C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/PAR1/VL,FWIDTH,SWIDTH,SLENGTH
      COMMON/GRASSD/PART(3),SC,SS,VN1,VN2,VN,GSI,H,VLCM,POR,CI,ICO
      DIMENSION PARTC(6,3)
      CHARACTER*19 PARTLBL(6),PARTLB
      CHARACTER*75 LISFIL(13)
      DATA (PARTC(I,1), I=1,6),(PARTC(I,2), I=1,6),(PARTC(I,3), I=1,6)
c-----particle diameters (cm)---------------
     &      /.0002D0, .001D0, .003D0, .03D0, .02D0, .0029D0,
c-----fall velocities (cm/s)-------------¬---
     & .0004D0, .0094D0, .0408D0 ,3.0625D0, 3.7431D0, .076D0,
c-----particle densities (g/cm3)------------
     &  2.6D0, 2.65D0, 1.80D0, 1.60D0, 2.65D0, 2.65D0/
      DATA (PARTLBL(I), I=1,6)/'    Clay','    Silt',' Sm.Aggregate',
     &  ' Lg.Aggregate','    Sand',' Silt (USDA)'/

c-----(.IGR file) Read vegetation inputs -----------------------------

      IF(ISCR.EQ.0) THEN
        WRITE(*,'(" ... Reading inputs from: ",A45)')LISFIL(7)
      ENDIF
      READ(12,*)SS, VN, H, VN2,ICO

c------(.ISD file) Read sediment inputs -----------------------------

      IF(ISCR.EQ.0) THEN
        WRITE(*,'(" ... Reading inputs from: ",A45)')LISFIL(11)
      ENDIF
      READ(16,*)NPART,COARSE,CI,POR

C--if the particle class is other than those defined above, calculate values--

      IF(NPART.LE.6)THEN
          DO 10 I=1,3
c-------NPART =1 to 6:sediment characteristics based on typical sediment classes
              PART(I)=PARTC(NPART,I)
10        CONTINUE
          PARTLB=PARTLBL(NPART)
        ELSE
          READ(16,*)dummy1,dummy2
          IF(NPART.EQ.7) THEN
c-------NPART =7:sediment characterstics provided by user (next line of ISD)
            DP=dummy1
            SG=dummy2
            PART(1)=DP
            PART(3)=SG
            PARTLB='    (user) coarse'
            IF(DP.LT.0.0037D0)PARTLB='    (user) fine'
           ELSE
c-------NPART =8:sediment characteristics calculated dynamically for each storm
c------ based on Reichengerger et al., (2022), where: ITILLAGE (0/1) from .isd;
c------ SILTFRAC (-) the soil silt fraction from .isd; SLENGTH (m) slope length
c------ of field from .iro; RMMH (mm/h) average rainfall for the event summed
c------ from .iro and divided over source area (SLENGTH*SWIDTH); Ei (sediment
c------ load into VFS, g/m2) is obtained from Ci in .isd and hydrograph in .iro
c------ over the source area; Qoin (mm/h) peak runoff into VFS from .iro.
            SILTFRAC=dummy1
            ITILLAGE=NINT(dummy2)
            IF((SILTFRAC.LT.0.d0).OR.(SILTFRAC.GT.1.d0)) STOP 'ERROR:
     & silt fraction in ISD file must be a fraction (0-1)'
            IF((ITILLAGE.NE.0).AND.(ITILLAGE.NE.1)) STOP 'ERROR:
     & itillage in ISD file must be 0: Conventional or 1: No Tillage'
            SAREA=SWIDTH*SLENGTH
            EIN=QSUM0*Ci*1.d6/SAREA
            QPIN=BCROPEAK*3.6d6/SAREA
c-------- dynamic d50 equation for each storm event ---------
           PART(1) = 5.38827d0+35.27513d0*ITILLAGE + 0.32195d0
     &       *DSQRT(Ein) + 0.06625d0*Rmmh + 0.11455d0*Qpin -29.57603d0*
     &      (1/SLENGTH) -1.51938d0*(1/sqrt(SILTFRAC))
            IF(PART(1).LT.1.1D0) PART(1)=1.1D0
            PART(1)=PART(1)/10000.d0
            PART(3)=2.65D0
            IF(PART(1).LE.0.0002D0) PART(3)=2.6D0
            DP=PART(1)
            SG=PART(3)
            PARTLB='   (storm) coarse'
            IF(DP.LT.0.0037D0)PARTLB='   (storm) fine'
          ENDIF
          
c----Particle fall velocity caculation:¬
c----a) Barfield et al. (1981) to calculate in GRASSF (dp in mm)-------
C            VC1=-0.34246727D0
C            VC2=0.98912185D0
C            VC3=1.1461280D0
C            CDP=DLOG(DP/10.D0)
C            PART(2)=10.D0**(VC1*CDP*CDP+VC2*CDP+VC3)
c----b) Fair and Geyer method (1954), based on Stokes (units in SI) ------
          AGRAV=9.8066352D0
          VISCOS=0.000001003352832D0
          MAXITER=100
          DIA=PART(1)/100.D0
          FV=AGRAV*(SG-1.D0)*DIA**2.D0/(18.D0*VISCOS)
          C1=DIA*VISCOS
          REYN=FV*C1
          IF(REYN.LE.0.1D0) GO TO 70
          C2=SQRT(4.D0*AGRAV*(SG-1.D0)*DIA/3.D0)
          ITER=0
          EPS=1.D0
          FVOLD=10000000.D0
          DO 69 WHILE(EPS.GT.0.00001D0.AND.ITER.LT.MAXITER)
              ITER=ITER+1
              CD=24.D0/REYN+3.D0/SQRT(REYN)+0.34D0
              FV=C2/SQRT(CD)
              REYN=FV*C1
              EPS=DABS(FV-FVOLD)
              FVOLD=FV
69        CONTINUE
70        PART(2)=FV*100
      ENDIF

c-------If particle is fine(<37 microns) don't run the wedge part (icoarse=0)-

      ICOARSE=1
      IF(COARSE.LE.0.D0)ICOARSE=0

c-------Output all input values for sediment transport-------------

      COARSEP= COARSE*100.D0
      WRITE(13,*)'Sediment parameters'
      WRITE(13,*)'-----------------------------------------'
      WRITE(13,200)'Sediment inflow concentration(Ci)=',CI,'g/cm3'
      WRITE(13,99)PARTLB
      WRITE(13,200)'     Particle size (diameter, dp)=',PART(1),'cm'
      WRITE(13,200)'      Particle fall velocity (Vf)=',PART(2),'cm/s'
      WRITE(13,200)'Particle weight density (gamma-s)=',PART(3),'g/cm3'
      WRITE(13,200)'% of particles with dp>0.0037 cm =',COARSEP,'%'
      WRITE(13,200)'Porosity of deposited sediment(P)=',POR*100.D0,'%'
      WRITE(13,*)
      WRITE(14,101)
      WRITE(14,102)
      WRITE(14,107)
      WRITE(13,*)'Filter parameters for sediment transport'
      WRITE(13,*)'-----------------------------------------'
      WRITE(13,200)'           Filter main slope (Sc)=',SC
      WRITE(13,200)'        Filter media spacing (Ss)=',SS,'cm'
      WRITE(13,200)'  Modified Manning coefficient (n)=',VN,'s.cm^-.33'
      WRITE(13,200)' Manning coeff. for bare soil (n2)=',VN2,'s.m^-.33'
      WRITE(13,200)'          Filter media height (H)=',H,'cm'
      WRITE(13,*)
      WRITE(13,1050)
      IF(ISCR.EQ.0) THEN
        WRITE(13,1060)
       ELSE
        WRITE(13,1061)
      ENDIF
      WRITE(13,1070)

99    FORMAT('                   Particle class =',a19)
101   FORMAT('   Time      qin      q1       Rs1       Vm1        ',
     & 'df1       q2       Rs2       Vm2       df2      q3        Rs3',
     & '        Vm3       df3      qout')
102   FORMAT('    (s)  (cm3/s/cm)(cm3/s/cm)  (cm)     (cm/s)     (cm)',
     & '    (cm3/s/cm)  (cm)     (cm/s)     (cm)   (cm3/s/cm)  (cm)',
     & '      (cm/s)      (cm) (cm3/s/cm)')
107   FORMAT(147('-'))
200   FORMAT(A35,F12.6,A11)
1050  FORMAT('  Time      Y(t)     X1(t)     X2(t)      L(t)',
     &'      Se       gsi       gsI       gs2       gso     Cum.gsi ',
     &  'Wedge_mass Lower_mass  Cum.gso     f        frac      DEP',
     &  '      CDEP       Tt')
1060  FORMAT('  (s)       (cm)     (cm)      (cm)       (cm)',
     &  '             (g/cm.s) (g/cm.s)   (g/cm.s)  (g/cm.s)   (g/cm)',
     &  '    (g/cm)    (g/cm)    (g/cm)                         (cm)')
1061  FORMAT('  (s)       (cm)     (cm)      (cm)       (cm)',
     &  '              (g/s)    (g/s)      (g/s)     (g/s)      (g) ',
     &  '      (g)       (g)       (g)                          (cm)')
1070  FORMAT(187('-'))


      RETURN
      END
