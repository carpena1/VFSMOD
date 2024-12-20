        SUBROUTINE KWWRITE(N,L,M,QTEMP,X,BCROQ,FWIDTH)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C      Write hydrology outputs to the filename.ohy file                       C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      COMMON/WTGA1/WTD,PARW(4),PARK(2),hb,RAINL,ZW,TW,RVH
      COMMON/PAR/QK(1001),R,THETAW,DX,DT,NDT,NELEM,MAXITER,NPOL,IOUT,NL
      COMMON/KWW/TIMELAST,OFLOWLAST,OSUMFLOW,QFIELDLAST,FSUMFLOW,FP0

      DIMENSION X(MAXEQN)

c--10/20/08--TIME=DT*(L-NDT*.005D0)
      TIME=DT*L
      DMAX=R*TIME
      IF(IOUT.EQ.1) THEN

C-------Print values for nodes for each time step (IOUT=1)-------

            MID=N/2+1
            VMID=(N-1)*DX/2
            DO 10 I=1,N
               VNODE=(I-1)*DX
               IF(VNODE.EQ.VMID)WRITE(11,200)TIME,VMID,X(MID),M
10            CONTINUE
         ELSE

c------Print values of flowrate at the downstream end of the plane(IOUT><0)
            TIMEINCR=TIME-TIMELAST
            OFLOWRATE=FWIDTH*QTEMP
            OSUMFLOW=OSUMFLOW+0.5D0*(OFLOWRATE+OFLOWLAST)*TIMEINCR
            QFIELD=BCROQ
            FSUMFLOW=FSUMFLOW+0.5D0*(QFIELD+QFIELDLAST)*TIMEINCR
            TIMELAST=TIME
            OFLOWLAST=OFLOWRATE
            QFIELDLAST=QFIELD
            IF(WTD.EQ.999.d0) THEN
               dF=(FP0+FPI)*0.5d0*TIMEINCR
               IF(DM.LE.0.0001D0) then
c------------------For initial saturated soil conditions, f=Ks and z=Ks.time
                   Z=AGA*TIME/OS
                 ELSE
c------------------For initial unsaturated soil, G-A wetting front depth z=F/M
                   Z=Z+(FP0+FPI)*0.5d0*TIMEINCR/DM
               ENDIF
               FP0=FPI
            ENDIF
            WRITE(11,400)TIME,OFLOWRATE,OSUMFLOW,R,QFIELD,FSUMFLOW,
     &         FPI,Z,M
      ENDIF

200   FORMAT(E11.4,4E12.4,I6)
400   FORMAT(E11.4,7E12.4,I6)

      RETURN
      END
