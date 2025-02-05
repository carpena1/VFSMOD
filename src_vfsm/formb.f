        SUBROUTINE FORMB(B0,X0,Q0,N,BCRO,PGPAR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C      In this subroutine the assembling of the right hand side part of the   C
C      equation (vector b).                                                   C
C      Where HOLD = is the depth of water at the the last time step           C
C            BO   = is a matrix of constant coefficients                      C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION   (A-H,O-Z)

      COMMON/PAR/QK(1001),R,THETAW,DX,DT,NDT,NELEM,MAXITER,NPOL,IOUT,NL
      COMMON/CINT/XI(20,20),W(20,20)
      DIMENSION PGPAR(4)
      DIMENSION PSI(4),DPSI(4),WF(4),BT0(4)
      DIMENSION X0(MAXEQN),Q0(MAXEQN),B0(MAXEQN)

C------Find dx1 for integration rule----------------------

      DX1=DX*(NPOL-1)

C-------Inizialize vector {b}------------------------------

      DO 10 I=1,N
            B0(I)=0.D0
10    CONTINUE

C-------Begin vector formation element by element----------

      VMN=5.D0/3.D0
      DO 60 NEL=1,NELEM

C-------Inizialize temporary vectors-----------------------

            K1 =(NPOL-1)*NEL-NPOL+2
            k2 =K1+1
            K3= K1+2
            BT0(1)=0.D0
            BT0(2)=0.D0
            BT0(3)=0.D0

C-------Begin integration point----------------------------

            DO 30 L=1,NL
                CALL SHAPEF(XI(L,NL),PSI,DPSI,WF,PGPAR)
                HOLD = PSI(1)*X0(K1)+PSI(2)*X0(K2)+PSI(3)*X0(K3)
                DQOLD = DPSI(1)*Q0(K1)+DPSI(2)*Q0(K2)+DPSI(3)*Q0(K3)
                RAIN= (PSI(1)+PSI(2)+PSI(3))*R
                BT0(1)=BT0(1)+(DX1*0.5D0*WF(1)*HOLD+WF(1)*RAIN*DX1
     &               *.5D0*DT -(1.D0-THETAW)*DT*WF(1)*DQOLD)*W(L,NL)
                BT0(2)=BT0(2)+(DX1*0.5D0*WF(2)*HOLD+WF(2)*RAIN*DX1
     &               *.5D0*DT -(1.D0-THETAW)*DT*WF(2)*DQOLD)*W(L,NL)
                BT0(3)=BT0(3)+(DX1*0.5D0*WF(3)*HOLD+WF(3)*RAIN*DX1
     &               *0.5D0*DT -(1.D0-THETAW)*DT*WF(3)*DQOLD)*W(L,NL)
30          CONTINUE

C-------Plug the element vector into the {b0} vector--------

            B0(K1)=B0(K1)+BT0(1)
            B0(K2)=B0(K2)+BT0(2)
            B0(K3)=B0(K3)+BT0(3)
60    CONTINUE

C------Plug in the boundary condition b(1)=0 --------------

      B0(1) = BCRO

      RETURN
      END
