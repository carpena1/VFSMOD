      SUBROUTINE multspcalc(DGMRES)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Subroutine for VFSMOD multispecies (parent-metabolites) degradation         C
C  11/2024, Rafael Muñoz-Carpena (c), UFL, carpena@ufl.edu                     C
C  Variables used:                                                             C 
C  DGMRES(mg)= pesticide surface residues                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/WQ1/VKD(10),CCP,CSAB(5),DGMRES0(10),DGMOL(10),DGFRAC(10,10)
      COMMON/IWQ2/NDGDAY,IDG,IWQ,IWQPRO,ICAT,IMOB
      COMMON/WQ3/DGKREF(10),FC,DGPIN(10),DGML,DGT(366),DGTHETA(366),DGLD(10),RF
      DIMENSION DGAI(10,10),DGF(10,10),DGDMASS(366,10),DGMRES(10)
      CHARACTER*75 LISFIL(13)
      CHARACTER*1 DASH

c------For EU FOCUS (IDG=1) reference temperature is 20 C, QT10=2.58,
c------EFSA Journal (2007):622:1-32
      DGBETA=-0.7d0
      DGC1=65.4d0/(0.008314d0)
      DGC2=1.d0/293.15d0
c------For US EPA (IDG=3) the reference temperature is 25 C, , QT10=2.0,
c------EPA/PMRA PRZM Guidance v.1.0 (2012) (Appendix A, page 5)
      DGC3=49.5d0/(8.314d0/1000.d0)
      DGC4=1.d0/298.15d0
c------Print section header in OWQ file (don;t print with IWQ=1 for version  compatibility)
      IF(IWQ.GT.1) THEN
        WRITE(18,160)'Product degradation between runoff events (mg)'
        WRITE(18,170)'Time(d)   ',('Product',II,II=1,IWQ)
      ENDIF
c-----A) Single species (IWQ=1) degradation (compatible ith previous version)
      IF(IWQ.EQ.1) THEN
        JJ=1
        IF(IDG.GT.0.AND.IDG.LE.4.AND.NDGDAY.GT.0) THEN
           DO 5 I=1,NDGDAY
            IF(IDG.EQ.1) THEN
c-------------- For EU FOCUS reference temperature is 20 C, QT10=2.58,
c-------------- EFSA Journal (2007):622:1-32
              DGKTEMP=DEXP(DGC1*(DGC2-(1.d0/(DGT(I)+273.15d0))))
              DGKTHETA=(DGTHETA(I)/FC)**DGBETA
             ELSEIF(IDG.EQ.2) THEN
              DGKTEMP=1.D0
              DGKTHETA=1.D0
             ELSEIF(IDG.EQ.3) THEN
c-------------- For US EPA the reference temperature is 25 C, , QT10=2.0,
c-------------- EPA/PMRA PRZM Guidance v.1.0 (2012) (Appendix A, page 5)
              DGKTEMP=DEXP(DGC3*(DGC4-(1.d0/(DGT(I)+273.15d0))))
              DGKTHETA=1.D0
             ELSEIF(IDG.EQ.4) THEN
              DGKTEMP=1.D0
              DGKTHETA=(DGTHETA(I)/FC)**DGBETA
            ENDIF
            DGKI=DGKREF(JJ)*DGKTEMP*DGKTHETA
            DGMRES(JJ)=DGMRES(JJ)*DEXP(-DGKI)
c           check for NaN
            IF(DGMRES(JJ).NE.DGMRES(JJ)) DGMRES(JJ)=0.d0
c            WRITE(18,180)I,DGMRES(JJ)       
5         CONTINUE
        ENDIF   
       ELSE
c-----(IWQ>1) Multispecies reaction (parent-metabilites) degration algorithm 
C-----Solve the daily mass degradation problem for multiple products with an explicit 
c-----finite difference backwards approximation. dt=0.01 days to minimize convergence errors.
        DGDT=0.01d0
        TIME=0.d0
C-----Assemble time invariant part of matrix A = [M]°[kref]°[f]---------
c------a) Since [f]=[f]^T (-), transpose adjacency matrix DGFRAC(II,JJ)
        DO 20 II=1,IWQ
          DO 10 JJ=1,IWQ 
            DGF(II,JJ)=DGFRAC(JJ,II)
10        CONTINUE
20      CONTINUE
c------b) Hadamard multiplication of elements to get time-invariant A matrix M, DGAI(II,JJ)
        DO 40 II=1,IWQ
          DO 30 JJ=1,IWQ
            DGAI(II,JJ)=DGF(II,JJ)*DGMOL(II)/DGMOL(JJ)*DGKREF(JJ)
30        CONTINUE
40      CONTINUE
c-----debug, print matrices for checking ------
!        DO 50 II=1,IWQ
!          write(*,'(10f10.6)') (DGAI(II,JJ),JJ=1,IWQ)
!50    CONTINUE
!        write(*,*)
!        DO 60 II=1,IWQ
!          write(*,'(10f10.6)') (DGF(II,JJ),JJ=1,IWQ)
!60    CONTINUE
!        write(*,*)
!        DO 70 II=1,IWQ
!           write(*,'(10f10.6)') (DGMOL(II)/DGMOL(JJ),JJ=1,IWQ)
!70    CONTINUE
!        write(*,*)
!        DO 80 II=1,IWQ
!        write(*,'(10f10.6)') (DGKREF(JJ),JJ=1,IWQ)
!80    CONTINUE
c----- Initial condition, residue in the VFS at the end of the runoff event
        DO 55 II=1,IWQ
          DGDMASS(1,II)=DGMRES(II)
55      CONTINUE
c-----Prepare temperature and moisture modifiers for first day based on IDG option
        NDAY=1
        IF(IDG.EQ.1) THEN
              DGKTEMP=DEXP(DGC1*(DGC2-(1.d0/(DGT(NDAY)+273.15d0))))
              DGKTHETA=(DGTHETA(NDAY)/FC)**DGBETA
          ELSEIF(IDG.EQ.2) THEN
              DGKTEMP=1.D0
              DGKTHETA=1.D0
          ELSEIF(IDG.EQ.3) THEN
              DGKTEMP=DEXP(DGC3*(DGC4-(1.d0/(DGT(NDAY)+273.15d0))))
              DGKTHETA=1.D0
          ELSEIF(IDG.EQ.4) THEN
              DGKTEMP=1.D0
              DGKTHETA=(DGTHETA(NDAY)/FC)**DGBETA
        ENDIF
c-----A) Multispecies case (IWQ>1): start time loop to calculate degradation along
C-----the NDGDAY days beteen runoff events. Get the daily temperature 
c-----and soil moisture modifiers, update [A] and solve daily---
        WRITE(18,210)0.d0,(DGDMASS(1,II),II=1,IWQ)
        DO 100 K=2,INT(NDGDAY/DGDT)+1
              TIME = (K-1) * DGDT
              DO 90 II=1,IWQ
                DGMULT=0.d0
                DO 85 JJ=1,IWQ
                  DGMULT=DGMULT+DGAI(II,JJ)*DGDMASS(K-1,JJ)
85              CONTINUE
c------- Modify by temperature (dgkTEMP) and moisture (dgkTHETA) as needed           
                DGMULT= DGMULT*DGKTEMP*DGKTHETA
                DGDMASS(K,II)=DGDMASS(K-1,II)+DGDT*DGMULT
90            CONTINUE
              IF (DMOD(TIME,1.d0).EQ.0.d0) THEN
                  WRITE(18,210)TIME,(DGDMASS(K,II),II=1,IWQ)
                  NEND=K
                  NDAY=NDAY+1
c---------------Change temperature and moisture modifiers based on IDG option
                  IF(IDG.EQ.1) THEN
                    DGKTEMP=DEXP(DGC1*(DGC2-(1.d0/(DGT(NDAY)+273.15d0))))
                    DGKTHETA=(DGTHETA(NDAY)/FC)**DGBETA
                   ELSEIF(IDG.EQ.2) THEN
                    DGKTEMP=1.D0
                    DGKTHETA=1.D0
                   ELSEIF(IDG.EQ.3) THEN
                    DGKTEMP=DEXP(DGC3*(DGC4-(1.d0/(DGT(NDAY)+273.15d0))))
                    DGKTHETA=1.D0
                   ELSEIF(IDG.EQ.4) THEN
                    DGKTEMP=1.D0
                    DGKTHETA=(DGTHETA(NDAY)/FC)**DGBETA
                  ENDIF
              ENDIF
100     CONTINUE
        DO II=1,IWQ
          DGMRES(II)=DGDMASS(NEND,II)
        END DO
      ENDIF

c-----Assign the final residues before the next storm to the DGMRES array

160   FORMAT(/,100('-'),/,4x,A46,/,100('-'))
170   FORMAT(4x,A7,10(A10,I2))
180   FORMAT(I11,E12.5)
200   FORMAT(4x,A7,10(A10,I2))
210   FORMAT(F12.0,10E12.5)

      RETURN
      END
