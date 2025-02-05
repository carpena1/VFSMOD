        SUBROUTINE EINSTEIN(GS2,NTRCAP,COARSE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                          C
C     This program solves Einstein's equation to find g_sd at the end of   C
C     C(t). After  Barfield et. al (1979)"Filtration of sediment by        C
C     simulated vegetation I.      Trans. ASAE, 22(5):540-548.             C
C                                                                          C
C     Known:                                                               C
C            dp: particle size, diameter (cm)                              C
C            Sc: filter main slope                                         C
C            Rs: hydraulic radius of the filter at D(t) (cm)               C
C            gammaw, gammasb:  water and sediment weight density (g/cm3)   C
C            gacc= acceleration due to gravity (980 cm/s2)                 C
C            coarse: % of particles from incoming sediment with diameter   C
C                    > 0.0037 cm (coarse fraction that will be routed      C
C                    through wedge).                                       C
C            gsi_c: coarse sediment load fraction entering filter (g/s/cm) C
C                                                                          C
C     Unknown:                                                             C
C            gs2: sediment load entering downstream section (g/s/cm)       C
C                                                                          C
C     NOTE: all units in CGS system (cm,g,s), including Manning's n        C
C                                                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/GRASSD/PART(3),SC,SS,VN1,VN2,VN,GSI,H,VLCM,POR,CI,ICO
      COMMON/GRASSD2/GSIMASS,GSOMASS,TTE,DEP,QSED(4),RS(3),DF(3),VM(3)

      IF (QSED(2).GT.0.d0) THEN
          GAMMAW=1.D0
          GACC=980.D0
          CC0=(PART(3)-GAMMAW)/GAMMAW
          CC1= PART(3)*DSQRT(CC0*GACC*PART(1)**3.D0)
          CC2= CC0*PART(1)/1.08D0/SC/RS(2)
          GS2=CC1 * CC2**(-1.D0/0.28D0)
c---------Check if the transport capacity is lower than concentration -------
c---------a) if lower, deposition at the wedge occurs (first part of step 3)
c---------b) if higher, there is enough energy to transport sediment through
C----------the wedge and no deposition occurs, all sediment is transported --
C----------to the suspended sediment zone (zones C(t) and D(t))(step3)

          GSI_C=COARSE*GSI
          IF(GS2.GT.GSI_C)THEN
               GS2=GSI
               NTRCAP=1
             ELSE
               NTRCAP=0
          ENDIF
        ELSE
          GS2=0.d0
          NTRCAP=0
      ENDIF

      RETURN
      END
