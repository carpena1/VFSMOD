      SUBROUTINE SHAPEF(XIS,PSI,DPSI,WF,PGPAR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C     SUBROUTINE SHAPEF CALCULATES THE VALUES OF THE WEIGHTING AND BASIS      C
C     FUNCTIONS PSI AND THEIR DERIVATIVES DPSI WITH RESPECT TO  THE MASTER    C
C     ELEMENT COORDINATES AT A SPECIFIED VALUE OF XIS.                        C
C                                                                             C
C     ANY TYPICAL ELEMENT   = [X1,Xk+1] CONSISTING OF k+1 NODES               C
C     X1, ..., Xk+1 IS ALWAYS NORMALIZED INTO THE MASTER                      C
C     ELEMENT   = [-1,1] BY THE TRANSFORMATION OVER A TYPICAL                 C
C     ELEMENT  [X1,Xk+1] THERE EXIST k+1 ELEMENT SHAPE FUNCTIONS PSI, EACH    C
C     IS A POLYNOMIAL OF DEGREE k,                                            C
C      XIS     :THE MASTER COORDINATE OF THE POINT AT WHICH VALUES OF PSI AND C
C      DPSI ARE DESIRED                                                       C
C      NPOL    :NUMBER OF NODES IN THE ELEMENT (= NUMBER OF SHAPE FUNCTIONS)  C
C      PSI(I)  :THE I-TH SHAPE FUNCTION AT XIS                                C
C      DPSI(I) :THE DERIVATIVE OF THE I-TH SHAPE FUNCTION AT XIS              C
C      WF(I)   :MODIFIED WEIGHTING FUNCTIONS                                  C
C      PGPAR(I):PETROV-GALERKIN PARAMETERS (I=1,4)                            C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/PAR/QK(1001),R,THETAW,DX,DT,NDT,NELEM,MAXISTER,NPOL,IOUT,NL
      DIMENSION PSI(4),DPSI(4),WF(4)
      DIMENSION PGPAR(4)

      IF (NPOL.LT.2.OR.NPOL.GT.4) GO TO 99
      GO TO (99,10,20,30) NPOL
C--------------------------------
C     LINEAR SHAPE FUNCTIONS
C--------------------------------
10    PSI(1) = 0.5D0*(1.D0-XIS)
      PSI(2) = 0.5D0*(1.D0+XIS)
      DPSI(1)= -0.5D0
      DPSI(2)= 0.5D0
      GO TO 80
C--------------------------------
C     QUADRATIC SHAPE FUNCTION
C--------------------------------
20    PSI(1) = XIS*(XIS-1.D0)*0.5D0
      PSI(2) = 1.D0-XIS**2.D0
      PSI(3) = XIS*(XIS+1.D0)*0.5D0
      DPSI(1)= XIS-0.5D0
      DPSI(2)= -2.D0*XIS
      DPSI(3)= XIS+0.5D0

C-----------Modified weighting functions-----------
      FCU= 5.D0/8.D0*XIS*(XIS+1.D0)*(XIS-1.D0)
      FQR= 21.D0/16.D0*(-XIS**4+XIS**2.D0)
      WF(1)= PSI(1) - PGPAR(1)*FCU - PGPAR(3)*FQR
      WF(2)= PSI(2) + 4.D0*PGPAR(2)*FCU + 4.D0*PGPAR(4)*FQR
      WF(3)= PSI(3) - PGPAR(1)*FCU - PGPAR(3)*FQR
      GO TO 80
C--------------------------------
C     CUBIC SHAPE FUNCTION
C--------------------------------
30    PSI(1) = 9.D0/16.D0*(1.D0/9.D0-XIS**2.D0)*(XIS-1.D0)
      PSI(2) = 27.D0/16.D0*(1.D0-XIS**2.D0)*(1.D0/3.D0-XIS)
      PSI(3) = 27.D0/16.D0*(1.D0-XIS**2.D0)*(1.D0/3.D0+XIS)
      PSI(4) = -9.D0/16.D0*(1.D0/9.D0-XIS**2.D0)*(1.D0+XIS)
      DPSI(1)= -9.D0/16.D0*(3.D0*XIS**2.D0-2.D0*XIS-1.D0/9.D0)
      DPSI(2)= 27.D0/16.D0*(3.D0*XIS**2.D0-2.D0/3.D0*XIS-1.D0)
      DPSI(3)= 27.D0/16.D0*(-3.D0*XIS**2.D0-2.D0/3.D0*XIS+1.D0)
      DPSI(4)= -9.D0/16.D0*(-3.D0*XIS**2.D0-2.D0*XIS+1.D0/9.D0)

80    RETURN

99    WRITE(6,*)'ERROR IN CALLING TO SHAPEF, NPOL= ', NPOL
      STOP
      END
