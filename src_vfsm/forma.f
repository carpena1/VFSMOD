        SUBROUTINE FORMA(A,NBAND,PGPAR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C      This subroutine assembles the system matrix [A] as a banded matrix.    C
C      This procedure involves the calculation of element matrices EK and     C
C      their accumulation in the banded system matrix [A]. Finally we end     C
C      up plugging in the BC for the problem.                                 C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION   (A-H,O-Z)

      COMMON/PAR/QK(1001),R,THETAW,DX,DT,NDT,NELEM,MAXITER,NPOL,IOUT,NL
      DIMENSION A(MAXEQN,MAXBND),EK(4,4),PGPAR(4)

      DO 10 NEL=1,NELEM

C---------------Form the element matrices EK-----------------

          CALL ELEM(EK,PGPAR)

C---------------Assemble the matrix------------------

          CALL ASSM(A,EK,NBAND,NEL)
10    CONTINUE

C-------Apply boundary conditions over A----------

      CALL BCA(A,NBAND)

      RETURN
      END
