        SUBROUTINE INI(A,B,X,XM,X0,Q0,QM,SSE,NODEX)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C      SET ALL MATRICES=0 TO START PROGRAM                                    C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION   (A-H,O-Z)

      DIMENSION A(MAXEQN,MAXBND),B(MAXEQN),NODEX(4)
      DIMENSION X(MAXEQN),X0(MAXEQN),XM(MAXEQN),QM(MAXEQN),Q0(MAXEQN)

      DO 10 I=1,MAXEQN
            B(I)= 0.D0
            X(I)= 0.D0
            X0(I)=0.D0
            XM(I)=0.D0
            QM(I)= 0.D0
            Q0(I)=0.D0
            DO 10 J=1,MAXBND
                A(I,J)= 0.D0
10    CONTINUE
            SSE=0.D0
      DO 20 I=1,4
            NODEX(I)=1
20    CONTINUE

      RETURN
      END
