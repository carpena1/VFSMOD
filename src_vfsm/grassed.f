       SUBROUTINE GRASSED(TIME,N,QIN,QOUT,NODEX,ICOARSE,COARSE,
     &      FWIDTH,ISCR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C     This subroutine solves the sediment transport problem on a grass filter C
C     It utilizes the method proposed by:                                     C
C                                                                             C
C  1. Tollner et al. (1976). "Suspended sediment filtration capacity of       C
C     simulated vegetation". Trans. ASAE. 19(4):698-682.                      C
C  2. Tollner et al. (1977). "Sediment deposition patterns in simulated       C
C     grass filters". Trans. ASAE. 20(5):940-944.                             C
C  3. Barfield et. al (1979)"Filtration of sediment by simulated vegetation   C
C     I. Trans. ASAE, 22(3):540-548.                                          C
C  4. Hayes et. al (1979)"Filtration of sediment by simulated vegetation II"  C
C     Trans. ASAE, 22(5):1063-1067                                            C
C  5. Hayes et. al (1984)"Performance of grass filters under laboratory and   C
C     field Conditions".Trans. ASAE, 27(5):1321-1331, to account for          C
C     triangular upslope deposition and particle and size distribution        C
C  6. Wilson et al (1981)"A Hydrology and sedimentology model: Part I.        C
C     Modeling techniques. U.of Kentucky. Lexington. This is a major          C
C     rewrite of the prodedures involved.                                     C
C  7. Haan et al (1994)"Design Hydrology and Sedimentology for Small          C
C     Catchments". Prentice-Hall. Chapter 9C contains updated (and clearer)   C
C     procedures for sediment trapping and wedge formation                    C
C                                                                             C
C     This is just the sediment transport unit calling the subroutines:       C
C             OCF, EINSTEIN, STEP3                                            C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/GRASSD/PART(3),SC,SS,VN1,VN2,VN,GSI,H,VLCM,POR,CI,ICO
      COMMON/GRASSD2/GSIMASS,GSOMASS,TTE,DEP,QSED(4),RS(3),DF(3),VM(3)
      COMMON/GRASSD3/SUSMASS,WEDGEMASS,NFUP
      COMMON/OLD/SEOLD,TOLD,XTOLD,YTOLD,CDEP,SE,VBTOLD,GSOLD
      DIMENSION XPOINTS(3),NODEX(4)

C----- Select flow and sediment load at filter entry (incl. GSIMASS 05/06/03)--
      FWID=FWIDTH*100.D0
      J=1
      GSI=QIN*CI
      DTG=TIME-TOLD
      GSIMASS=GSIMASS+GSI*DTG

C----- (16/12/00) New way of handling a sediment filled strip. After this
c----- happens, all the sediment inflow is routed to the outflow end (no
c----- deposition) and without stopping the simulation. This way the
c----- simulation summary (*.osm) now shows realistic trpping efficiencies.
c----- Check if strip is filled-up (NFUP=1)
      IF(NFUP.EQ.1) THEN
            YT=H
            X1=YT/SC
            XT=VLCM
            VLT=0.D0
            SE=0.D0
            GSSI=GSI
            FI=0.D0
            GS2=GSI
            GSO=GSI
            FRAC=0.D0
            GSOMASS=GSOMASS+0.5D0*(GSO+GSOLD)*DTG
            GSOLD=GSO
            TOTRAP=0.D0
            DEP=H
            IF(ISCR.EQ.0) THEN
                  WRITE(13,180)TIME,YT,X1,XT,VLT,SE,GSI,GSSI,GS2,GSO,
     &            GSIMASS,WEDGEMASS,SUSMASS,GSOMASS,FI,FRAC,DEP,CDEP,
     &            TOTRAP
              ELSE
                  WRITE(13,180)TIME,YT,X1,XT,VLT,SE,GSI*FWID,GSSI*FWID,
     &            GS2*FWID,GSO*FWID,GSIMASS*FWID,WEDGEMASS*FWID,
     &            SUSMASS*FWID,GSOMASS*FWID,FI,FRAC,DEP,CDEP,TOTRAP
            ENDIF
            RETURN
      ENDIF

C-------STEP 1: Solves hydraulic properties at points (1), (2), (3) of ----
C-------        the filter to be used later on------------
      DO 10 NPLACE=1,3
          CALL OCF(NPLACE)
10    CONTINUE
      WRITE(14,200)TIME,QIN,(QSED(J),RS(J),VM(J),DF(J),J=1,3), QSED(4)
c-----Check and fix errors in incoming hydro/sedimentograph (values <0)
      IF((QIN.LE.0.D0.OR.QSED(2).LE.0.D0).AND.QOUT.GT.0.D0) THEN
          GS2=0.D0
          RETURN
        ELSEIF((QIN.LE.0.D0.OR.QSED(2).LE.0.D0).AND.QOUT.LE.0.D0) THEN
          YT=YTOLD
          X1=YT/SC
          XT=XTOLD
          VLT=VLCM-Xt
          GSSI=GSI
          FI=0.D0
          GS2=GSI
          GSO=0.D0
          FRAC=0.D0
          GSOMASS=GSOMASS+0.5D0*(GSO+GSOLD)*DTG
          GSOLD=GSO
          TOTRAP=0.D0
          IF(ISCR.EQ.0) THEN
                WRITE(13,180)TIME,YT,X1,XT,VLT,SE,GSI,GSSI,GS2,GSO,
     &          GSIMASS,WEDGEMASS,SUSMASS,GSOMASS,FI,FRAC,DEP,CDEP,
     &          TOTRAP
            ELSE
                WRITE(13,180)TIME,YT,X1,XT,VLT,SE,GSI*FWID,GSSI*FWID,
     &          GS2*FWID,GSO*FWID,GSIMASS*FWID,WEDGEMASS*FWID,
     &          SUSMASS*FWID,GSOMASS*FWID,FI,FRAC,DEP,CDEP,TOTRAP
          ENDIF
          RETURN
      ENDIF
C-------STEP 2: Solve Einstein's equation to find transport capacity (gs2)--
C-------        at the end of B(t) -----------------------------------------
      CALL EINSTEIN(GS2,NTRCAP,COARSE)
      IF(ICOARSE.EQ.0)NTRCAP=1

C-------STEP 3: Calculate shape of sediment wedge, sediment outflow, and ---
C-------        trapping efficiency for the filter  and finishes up---------
      CALL STEP3(GS2,TIME,NTRCAP,COARSE,QOUT,FWIDTH,ISCR)

C-------STEP 4: Position points (1), (2), (3) at system nodes so that flow --
C-------        rates can be read at those points at next time step ---------
      CALL POINTS(N,XPOINTS,NODEX,VBT)

180   FORMAT(f7.0,4F10.3,F10.6,13E10.3)
200   FORMAT(F7.0,14E10.3)
201   FORMAT(A32,F8.4,a4)

130   RETURN
      END
