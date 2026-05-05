        SUBROUTINE BCA(A,NBAND)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C     The boundary conditions for this program are (problem set 2, part 1)    C
C        C (x=0,t>0) = 1                                                      C
C        dC/dx (x=xL)= 0                                                      C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION   (A-H,O-Z)

      DIMENSION A(MAXEQN,MAXBND)

      NDIAG = NBAND/2 + 1

C-------Plug in first kind of BC (Dirichlet)------

      DO 10 I=1,NBAND
            A(1,I)=0.D0
10    CONTINUE
      A(1,NDIAG)=1.D0

      RETURN
      END
