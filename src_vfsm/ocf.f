        SUBROUTINE OCF(NPLACE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C      This program solves the hydraulic properties to be used in later       C
C      steps of the problem by means of Manning's and Channel Flow            C
C      theory and a Newton-Raphson iterative method. It utilizes the method   C
C       proposed by Barfield et. al (1979)"Filtration of sediment by          C
C      simulated vegetation I.      Trans. ASAE, 22(5):540-548.               C
C                                                                             C
C      Known:                                                                 C
C            Ss: spacing of the filter media elments (cm)                     C
C            Sc: filter main slope                                            C
C            n: Manning's n= 0.0072 for cilindrical media (s/cm^1/3)          C
C            q_ocf: overland flow (cm2/s)                                     C
C      Unknown:                                                               C
C            df: depth of flow at D(t) (cm)                                   C
C            Vm: depth averaged velocity at D(t)(cm/s)                        C
C            Rs: hydraulic radius of the filter (cm)                          C
C                                                                             C
C      NOTE: all units in CGS system (cm,g,s), including Manning's n          C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/GRASSD/PART(3),SC,SS,VN1,VN2,VN,GSI,H,VLCM,POR,CI,ICO
      COMMON/GRASSD2/GSIMASS,GSOMASS,TTE,DEP,QSED(4),RS(3),DF(3),VM(3)

C------- flow depth and velocity set to zero for no flow at any given point --
      IF(QSED(NPLACE).LE.0.D0)THEN
            VM(NPLACE)=0.D0
            DF(NPLACE)=0.D0
            GOTO 130
      ENDIF
c------- otherwise, calculate Rs, Vm and df for the given point -----------
      ALPHA=DSQRT(SC)/VN
      QOCF=QSED(NPLACE)
      C0=(1.D0/ALPHA)**1.5D0
      C1= 2.D0*QOCF
      C2=SS
      C3=SS*QOCF

c------Newton-Raphson Method--------------------

      VMD=0.000001D0
      MAXIT=100
      ERROR=1.D0
      NS=0
      DO 120 WHILE(NS.LT.MAXIT.AND.ERROR.GT.1D-8)
            VMDO=VMD
            PHI=C0*(VMD**1.5D0)*(C1+C2*VMD)-C3
            DPHI=1.5D0*C0*DSQRT(VMD)*(C1+C2*VMD)+C0*(VMD**1.5D0)*C2
            VMD=VMDO-PHI/DPHI
            ERROR=DABS(VMDO-VMD)
            NS=NS+1
120   CONTINUE
      VM(NPLACE)=VMD
      DF(NPLACE)=QOCF/VMD
130   RS(NPLACE)=DF(NPLACE)*SS/(2*DF(NPLACE)+SS)

      RETURN
      END
