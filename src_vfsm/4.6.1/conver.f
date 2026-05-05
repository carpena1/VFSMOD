        SUBROUTINE CONVER(N,X,XM,MFLAG)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C               This subroutine checks for convergence as:                    C
C                       Max (deltaX)/Max (X^m+1 ) <= eps                      C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION   (A-H,O-Z)

      DIMENSION X(MAXEQN),XM(MAXEQN)

      EPS = .00000001D0
      BIG1=0.D0
      BIG2=0.D0
      DO 10 I=1,N
            C1= ABS(X(I)-XM(I))
            IF(X(I).LT.0.D0) THEN
                   C1=0.D0
                  X(I)=0.D0
            ENDIF
            C2= ABS(XM(I))
            BIG1=MAX(BIG1,C1)
            BIG2=MAX(BIG2,C2)
10    CONTINUE
      IF(BIG2.EQ.0.D0)THEN
            TOL=1.1D0*EPS
           ELSE
            TOL=BIG1/BIG2
      ENDIF

C-------If there is convergence return MFLAG=1, otherwise 0----

      MFLAG= 0
      IF(TOL.LE.EPS) MFLAG=1

      RETURN
      END
