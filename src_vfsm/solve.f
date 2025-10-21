        SUBROUTINE SOLVE(A,B,X,N,NBAND)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C     SOLVE THE TRANSFORMED MATRIX A USING A BACKWARD AND FORWARD SUBTITUTION C
C      SUCH: [A] {x} = {b}                                                    C
C     Since [A]= [L].[U] then  [L][U]{x} = [L]([U]{x})= [L]{y}= {b}           C
C     solving  [L]{y}=b    (forward subtitution)                              C	
C	             [U]{x}=y    (backward substitution)                            C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION  (A-H, O-Z)

      DIMENSION A(MAXEQN,MAXBND), B(MAXEQN), X(MAXEQN)

      NDIAG = (NBAND/2)+1
      NMAX = NBAND-NDIAG
      M = N  -1
      DO 10 I = 1,M
         NA = NMAX
         IF(NA .GT. N  -I) NA = N  -I
         DO 10 J = 1,NA
            B(I+J) = B(I+J)+A(I+J,NDIAG-J)*B(I)
10    CONTINUE
      X(N) = B(N)/A(N,NDIAG)
      DO 20 J = M,1,-1
         NA = NMAX
         IF(NA .GT. N  -J) NA= N  -J
         DO 30 K = 1,NA
            B(J) = B(J)-X(J+K)*A(J,NDIAG+K)
30       CONTINUE
         X(J) = B(J)/A(J,NDIAG)
20    CONTINUE

      RETURN
      END
