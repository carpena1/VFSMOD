        SUBROUTINE UPDATE(N,X,X0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                              l    l+1                       C
C     Refresh values of the X vector, this is X  = X                          C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION   (A-H,O-Z)

      DIMENSION X(MAXEQN),X0(MAXEQN)
      DO 10 I=1,N
            X0(I)=X(I)
10    CONTINUE

      RETURN
      END
