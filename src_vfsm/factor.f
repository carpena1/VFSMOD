        SUBROUTINE FACTOR (A,N,NBAND)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C    PERFORM THE LOWER AND UPPER DECOMPOSITON OVER SYSTEM MATRIX A  AND       C
C    AND STORE THE LOWER AND UPPER TRIANGULAR MATRICES ON THE OLD A MATRIX    C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION  (A-H, O-Z)

      DIMENSION A(MAXEQN,MAXBND)

      NMAX = NBAND/2
      NDIAG = NMAX + 1
      M = N -1
      DO 10 I = 1,M
         NA  = NMAX
         IF(N -I .LT. NA) NA  = N -I
         DO 10 J = 1,NA
            NEQN  = I+J
            A(NEQN ,NDIAG-J) = -A(NEQN ,NDIAG-J)/A(I,NDIAG)
            NLOW = NDIAG-J+1
            NHIGH = NLOW+NA-1
            DO 10 K = NLOW,NHIGH
               A(NEQN,K) =A(NEQN,K)+A(NEQN,NDIAG-J)*A(I,NDIAG+1-NLOW+K)
   10 CONTINUE

      RETURN
      END
