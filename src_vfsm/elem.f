        SUBROUTINE ELEM(EK,PGPAR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C       SUBROUTINE ELEM EVALUATE THE COMPONENTS                               C
C       FOR THE ELEMENT STIFFNESS MATRIX k.                                   C
C                                                                             C
C       NPOL    - NUMBER OF NODAL POINTS IN THE ELEMENT                       C
C       THETAW  - TIME WEIGHTING FACTOR                                       C
C       EK(I,J) - ENTRY IN ELEMENT STIFFNESS MATRIX                           C
C       NL      - ORDER OF THE INTEGRATION RULE                               C
C       XI(L)   - c1, THE LOCATION OF THE L-TH GAUSS ABSCISSA                 C
C       W(L)    - w1, the l-th gauss quadrature  point                        C
C       DX1        - betwen nodes in element                                  C
C       PSI(I)      - Weighting functions                                     C
C       PSI(I),DPSI(I) -Basis functions and derivatives                       C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/PAR/QK(1001),R,THETAW,DX,DT,NDT,NELEM,MAXITER,NPOL,IOUT,NL
      COMMON/CINT/XI(20,20),W(20,20)
      DIMENSION PGPAR(4)
      DIMENSION EK(4,4),PSI(4),DPSI(4),WF(4)

C----------Inizialize element arrays-----------

      DO 10 I=1,NPOL
          DO 10 J=1,NPOL
            EK(I,J) = 0.D0
10    CONTINUE

C----------Begin integration point loop---------

      DX1=DX*(NPOL-1)
      DO 20 L=1,NL
          CALL SHAPEF(XI(L,NL),PSI,DPSI,WF,PGPAR)
            DO 20 J=1,NPOL
              DO 20 I=1,NPOL
                  EK(I,J)=EK(I,J)+(DX1*0.5D0*WF(I)*PSI(J))*W(L,NL)
20    CONTINUE

       RETURN
       END
