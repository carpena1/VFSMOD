        SUBROUTINE GASUB(TIME,DT,L,R,RAIN,NEND,TRAI)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C      The subroutine solves and infiltration problem for unsteady rainfall   C
C     case using the Green-Ampt infiltration model. The program is based on   C
C     Skaggs and Khaleel (1982) in Hydrologic modeling of small watersheds    C
C     ASAE monograph no. 5, and Chu (1978) Water Resour. Res.                 C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C      THIS IS A MODIFIED VERSION OF THE ORIGINAL METHOD. IN THIS CASE        C
C      THE INFILTRATION IS ALLOWED TO BE ITS MAXIMUM POTENTIAL AFTER THE      C
C      FIRST PODING. THE IMPLICATION HERE IS THAT AFTER THE ORIGINAL          C
C      PONDING IS ACHIEVED, THE RUNOFF WATER MOVING AT THE SURFACE WILL       C
C     SUPPLY ENOUGH WATER TO SUSTAIN THE MAXIMUM POSSIBLE INFILTRATION        C
C     FOR THAT TIME STEP, IN OTHER WORDS, THE EFFECTIVE RAILFALL FEEDED       C
C      INTO THE MAIN FE-PROGRAM (R) WILL BE IN MOST CASES A NEGATIVE VALUE    C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C      Ponding indicator: beginning of the period (NPOND=1 ponded;0 no-ponded)C
C       (Note:rmc,3/2011)  end of the period (CU>0 ponded;CU=<0 no-ponded)    C
C                                                                             C
C     units: m,s                                                              C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      COMMON/GA2/LO,NPOND,NPFORCE,NPFORCE0
      DIMENSION RAIN(200,2)

      MAXIT=100
      C2=BGA/AGA
      FPIO=FPI
c------Check if end of field runoff has been reached and if so reset case as---
c------no-ponding at beginning and end of the period (NPOND=0,CU<0)(Case 1b)---
c-----rmc03/2011-new surface ponding forcing scheme (NPFORCE=1) when overland flow
c-----sustains infiltration capacity regardless of point excess rainfall calculation

      IF(NPFORCE.EQ.1.and.NPOND.EQ.0) NPOND=1
      IF(NEND.EQ.1)THEN
            NPOND=0
            CU=-1.D0
      ENDIF

c------Check if time step belongs to the same rainfall period as before (L=LO)-
c------and if not reset the total rainfall value for the end of the period ----
      IF(L.NE.LO) THEN
            PSOLD=PS
            PST = PS + RAIN(L,2)*(RAIN(L+1,1)-RAIN(L,1))
c-----------rmc03/2011- missing case for NPOND=1 and new_rain<end_infiltration-
c-----------of last period (fp,capacity since soil was ponded, i.e. NPOND=1)---
            IF(NPOND.EQ.1.AND.FPI.GT.RAIN(L,2).AND.NPFORCE.EQ.0) NPOND=0
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Surface ponding at beginning? y(npond=1), n (npond=0):                  C
C      CASE 1. Without surface ponding at the beginning of the period         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(NPOND.EQ.0) THEN
c-----------If new rainfall period calculate the cum. infiltration at tp (FP)
c-----------to find out if at the end of period there is ponding (CU>0)
            IF(L.NE.LO) THEN
              FP= BGA/(RAIN(L,2)-AGA)
              IF(FP.LT.0.D0) FP=F
              CU= PST - RO - FP
              IF(AGA.GE.RAIN(L,2))CU=-1.D0
            ENDIF

C-----1-a) With ponding at the end of the period (Cu >0)-----------------------

            IF (CU.GT.0.D0) THEN
              fcase=1.1d0
c-------------Calculate ponding time (tp) and shift time (tpp) when we regain ponding
c-------------within the current the  time step (i.e. there was no ponding (npond=0)
c-------------at the beginning of time step but there is at the end (Cu>0).
              IF(L.NE.LO) THEN
                     TP=(FP-PS+RO)/RAIN(L,2)+RAIN(L,1)
                     IF (TP.LT.RAIN(L,1))THEN
                            TP=RAIN(L,1)
                            FP=F
                     ENDIF
                     TPP=(FP-C2*DLOG(1.D0+FP/(C2)))/AGA
              ENDIF
              PS= PSOLD + RAIN(L,2)*(TIME-RAIN(L,1))
c-------------Split time step in before tp (fpi=i) and after (fpi=f)
              IF(TIME.LT.TP)THEN
                     F= PS -RO
                     FPI=RAIN(L,2)
                ELSE
                     C1= AGA*(TIME-TP+TPP)
                     F=RNEWTON(C1,C2,F)
                     FPI=AGA+BGA/F
                     EXCESS= PS -F - RO -STO
                     STO =STO +EXCESS
                     IF (STO.GT.SM) THEN
                            ROI=STO-SM
                            STO=SM
                            RO =RO +ROI
                     ENDIF
              ENDIF
c-------------Left hand side of kinematic wave equation (re=i-f)
              R=RAIN(L,2)-FPI
c-------------If the end of time step falls in new rain period calculate f
c-------------for that time and reset initial condition to ponding
              IF((TIME+DT).GT.RAIN(L+1,1)) THEN
                     C1=AGA*(RAIN(L+1,1)-TP+TPP)
                     F=RNEWTON(C1,C2,F)
                     NPOND=1
              ENDIF

C-----1-b) No ponding at the end of the period (Cu <0, fpi=i)-----------
             ELSE
              fcase=1.2d0
              PS= PSOLD + RAIN(L,2)*(TIME-RAIN(L,1))
c-------------rmc/5/2011/ account for storage in mass balance at end of runoff
              IF(STO.GT.0.D0.AND.NEND.EQ.1)THEN
                     F=PS-RO+STO
                ELSE
                     F=PS-RO
              ENDIF
              FPI=RAIN(L,2)
c-------------Left hand side of kinematic wave equation (ie=i-f)
              R=RAIN(L,2)-FPI
              IF((TIME+DT).GT.RAIN(L+1,1))F = PST -RO
            ENDIF

C-----------Do mass balance at the end of this rainfall period----------
            IF((TIME+DT).GT.RAIN(L+1,1)) THEN
              PS=PST
              EXCESS= PS -F - RO -STO
              STO =STO +EXCESS
              IF (STO.GT.SM) THEN
                     ROI=STO-SM
                     STO=SM
                     RO =RO +ROI
              ENDIF
            ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      CASE 2. With surface ponding at the beginning of the period (NPOND=1)  C
C     Note: this assumes that once field inflow starts (at NCHCK) the soil    C
C     will be kept ponded by overland runoff until the inflow hydrograph ends C
C     (i.e. NEND=1). The case when ponding ceases at the end of the time step C
C     is not considered since the time steps are small and the error caused   C
C     by ignoring this is neglibible (i.e considering only ponding that step).C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        ELSE
            fcase=2.0d0
            PS= PSOLD + RAIN(L,2)*(TIME-RAIN(L,1))
            C1= AGA*(TIME-TP+TPP)
            F=RNEWTON(C1,C2,F)
            IF(ISNAN(F)) F=PS
            FPI=AGA+BGA/F
            ROI= PS -F - RO
            IF(ROI.LT.0.d0) ROI=0.d0
            RO =RO +ROI
            R=RAIN(L,2)-FPI
            IF((TIME+DT).GT.RAIN(L+1,1))THEN
              F=PST-RO
              PS=PST
              ROI= PS -F - RO -STO
              RO =RO +ROI
            ENDIF
      ENDIF

c-----Set rainfall period for next time step as the current one and total
c-----cummulative rain to the one just calculated in this time step
      LO=L
      TRAI=PS

c--for debug--gasub-----
c      write(19,100)time,tw,z,F,PS,RO,FPI,RAIN(L,2),tp,tpp,NPFORCE,
c     &             NPOND,NEND,fcase,CU
c100   format(2f10.1,6e11.3,2f12.4,3i2,f5.2,3f10.4)
c--for debug--gasub-----

      RETURN
      END



      FUNCTION RNEWTON(C1,C2,F)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

      VF=F+1.0D-6
      MAXIT=101000
      ERROR=1.D0
      NS=0
      DO 120 WHILE(NS.LT.MAXIT.AND.ERROR.GT.1D-8)
            VFO=VF
            PHI=C1-VF+C2*DLOG(1.D0+VF/C2)
            DPHI=C2/(C2+VF)-1.D0
            VF=VFO-PHI/DPHI
            ERROR=DABS(VFO-VF)
            NS=NS+1
120        CONTINUE
      RNEWTON=VF
      if(ns.ge.maxit) print*,'ns=',ns

      return
      END
