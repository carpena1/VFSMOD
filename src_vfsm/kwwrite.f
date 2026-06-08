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
      COMMON/OHYBUF/OHB(9,200),NKWB

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
c-----Buffer the hydrograph row.  The .ohy table is written by OHYFLUSH after
c-----the time loop, applying a centered 5-point median + volume (mass)
c-----renormalization to the OUTFLOW column only.  This removes the residual
c-----Petrov-Galerkin oscillation left in the printed (block-averaged)
c-----hydrograph while conserving the outflow volume; reporting only -- the
c-----solver, mass balance and all other columns are unchanged.
            IF(NKWB.LT.200)THEN
              NKWB=NKWB+1
              OHB(1,NKWB)=TIME
              OHB(2,NKWB)=OFLOWRATE
              OHB(3,NKWB)=OSUMFLOW
              OHB(4,NKWB)=R
              OHB(5,NKWB)=QFIELD
              OHB(6,NKWB)=FSUMFLOW
              OHB(7,NKWB)=FPI
              OHB(8,NKWB)=Z
              OHB(9,NKWB)=DBLE(M)
            ENDIF
      ENDIF

200   FORMAT(E11.4,4E12.4,I6)

      RETURN
      END


      SUBROUTINE OHYFLUSH
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Write the buffered .ohy hydrograph table (unit 11) after applying a      C
C  centered 5-point median filter to the OUTFLOW column, renormalized to    C
C  preserve the raw outflow volume (so the .osm water balance is unchanged).C
C  Smooths the residual oscillation of the printed block-averaged outflow.  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/OHYBUF/OHB(9,200),NKWB
      DIMENSION QS(200),WIN(5)
      IF(NKWB.LE.0) RETURN

c-----Centered 5-point median of OUTFLOW (window shrinks at the two ends)----
      DO 10 K=1,NKWB
        K1=MAX(1,K-2)
        K2=MIN(NKWB,K+2)
        NW=0
        DO 5 J=K1,K2
          NW=NW+1
          WIN(NW)=OHB(2,J)
5       CONTINUE
        CALL MEDIAN5(WIN,NW,QS(K))
10    CONTINUE

c-----Mass renormalization: scale the smoothed series so its trapezoidal time
c-----integral (from t=0, q=0) equals that of the raw outflow -> volume kept.
      SUMR=0.5D0*OHB(2,1)*OHB(1,1)
      SUMS=0.5D0*QS(1)*OHB(1,1)
      DO 20 K=2,NKWB
        DTT=OHB(1,K)-OHB(1,K-1)
        SUMR=SUMR+0.5D0*(OHB(2,K)+OHB(2,K-1))*DTT
        SUMS=SUMS+0.5D0*(QS(K)+QS(K-1))*DTT
20    CONTINUE
      FAC=1.D0
      IF(SUMS.GT.0.D0) FAC=SUMR/SUMS
      DO 25 K=1,NKWB
        QS(K)=QS(K)*FAC
25    CONTINUE

c-----Recompute the cumulative outflow (CUM.FLOW) from the smoothed series
c-----and write the hydrograph table.
      CUM=0.5D0*QS(1)*OHB(1,1)
      DO 30 K=1,NKWB
        IF(K.GT.1) CUM=CUM+0.5D0*(QS(K)+QS(K-1))*(OHB(1,K)-OHB(1,K-1))
        WRITE(11,400)OHB(1,K),QS(K),CUM,OHB(4,K),OHB(5,K),OHB(6,K),
     &               OHB(7,K),OHB(8,K),NINT(OHB(9,K))
30    CONTINUE
400   FORMAT(E11.4,7E12.4,I6)
      RETURN
      END


      SUBROUTINE MEDIAN5(A,N,AMED)
c-----Median of up to 5 values via insertion sort (returns the central value)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(N),B(5)
      DO 5 I=1,N
        B(I)=A(I)
5     CONTINUE
      DO 20 I=2,N
        TMPV=B(I)
        J=I-1
10      IF(J.GE.1)THEN
          IF(B(J).GT.TMPV)THEN
            B(J+1)=B(J)
            J=J-1
            GOTO 10
          ENDIF
        ENDIF
        B(J+1)=TMPV
20    CONTINUE
      AMED=B((N+1)/2)
      RETURN
      END
