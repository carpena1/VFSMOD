        SUBROUTINE POINTS(N,XPOINTS,NODEX,VBT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                           C
C     Find qk, x1,x2,x3 to feed back to the overland flow. Definition of    C
C     variables:                                                            C
C                                                                           C
C     XPOINTS(1)= Mid-point of downface of sediment wedge (cm)              C 
C     XPOINTS(2)= Bottom point of downface of sediment wedge (cm)           C 
C     XPOINTS(3)= Mid-point of efective filter lnegth L(t) (cm)             C 
C     NODEX(1)= Node for mid-point of downface of sediment wedge            C 
C     NODEX(2)= Node for bottom point of downface of sediment wedge         C 
C     NODEX(3)= Node for mid-point of efective filter length L(t)           C 
C     NODEX(4)= Node for beginning of downface of sediment wedge            C 
C                                                                           C
C     Sections of the filter:                                               C
C                                                                           C
C     Section A= Top flat face of sediment wedge: slope=SC, n=VN2 (bare)    C
C     Section B= Downface of sediment wedge: slope=SET, n=VN2 (length VBT)  C
C     Section C & D= Effective filter (VLT): slope= unchanged, n= unchanged C
C                                                                           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/PAR/QK(1001),R,THETAW,DX,DT,NDT,NELEM,MAXITER,NPOL,IOUT,NL
      COMMON/GRASSD/PART(3),SC,SS,VN1,VN2,VN,GSI,H,VLCM,POR,CI,ICO
      COMMON/OLD/SEOLD,TOLD,XTOLD,YTOLD,CDEP,SE,VBTOLD,GSOLD
      DIMENSION XPOINTS(3),NODEX(4)

C-------Find points for each of the areas in filter---------------

      SET=SC+SE
      VLT=VLCM-XTOLD
      IF(YTOLD.EQ.0.D0.OR.SE.EQ.0.D0)THEN
            VBT=VBTOLD
         ELSE
            VBT=YTOLD/SE
      ENDIF
      IF(VBT.GT.XTOLD)VBT=XTOLD
      XPOINTS(1)=XTOLD-0.5D0*VBT
      XPOINTS(2)=XTOLD
      XPOINTS(3)=VLCM-VLT*0.5D0
      XPOINTBT=XTOLD-VBT
      NODEX(1)=IDNINT(XPOINTS(1)/DX/100.d0)+1
      NODEX(2)=IDNINT(XPOINTS(2)/DX/100.d0)+1
      NODEX(3)=IDNINT(XPOINTS(3)/DX/100.d0)+1
      NODEX(4)=IDNINT(XPOINTBT/DX/100.d0)+1
      
C-------If required reshape the surface topography of the filter.------

      IF(ICO.EQ.1)THEN
            ALPHA1=DSQRT(SC)/VN2
            ALPHA2=DSQRT(SET)/VN1
            DO 10 I=1,N
                IF(I.LE.NODEX(4))QK(I)=ALPHA1
                IF(I.GT.NODEX(4).AND.I.LE.NODEX(2))QK(I)=ALPHA2
10          CONTINUE
      ENDIF
      VBTOLD=VBT

      RETURN
      END
                  
