         SUBROUTINE MODIFY(QM,B,BCRO,PGPAR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C      In this subroutine the assembling of the right hand side part of the   C
C      equation (vector b) following the procedure discussed by Vieux et al   C
C      (1990)                                                                 C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION   (A-H,O-Z)

      COMMON/PAR/QK(1001),R,THETAW,DX,DT,NDT,NELEM,MAXITER,NPOL,IOUT,NL
      COMMON/CINT/XI(20,20),W(20,20)
      DIMENSION PGPAR(4)
      DIMENSION QM(MAXEQN),PSI(4),DPSI(4),WF(4),B(MAXEQN),BT(4)

C------Find dx1 for integration rule----------------------

      DX1=DX*(NPOL-1)

C-------Begin vector formation element by element----------

      VMN=5.D0/3.D0
      DO 60 NEL=1,NELEM
            K1=(NPOL-1)*NEL-NPOL+2
            K2=K1+1
            K3=K1+2
            BT(1)=0.D0
            BT(2)=0.D0
            BT(3)=0.D0

C-------Begin integration point----------------------------

            DO 30 L=1,NL
                CALL SHAPEF(XI(L,NL),PSI,DPSI,WF,PGPAR)
                DQM=DPSI(1)*QM(K1)+DPSI(2)*QM(K2)+DPSI(3)*QM(K3)
                BT(1)=BT(1)+WF(1)*DQM*W(L,NL)
                BT(2)=BT(2)+WF(2)*DQM*W(L,NL)
                BT(3)=BT(3)+WF(3)*DQM*W(L,NL)
30          CONTINUE

C-------Plug the element vector into the {b} vector--------


            B(K1)=B(K1)-THETAW*DT*BT(1)
            B(K2)=B(K2)-THETAW*DT*BT(2)
            B(K3)=B(K3)-THETAW*DT*BT(3)
60    CONTINUE

C------Plug in the boundary condition b(1)=0 --------------

      B(1) = BCRO

      RETURN
      END
