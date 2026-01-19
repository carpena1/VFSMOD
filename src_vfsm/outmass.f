      SUBROUTINE OUTMASS(TRAI,LISFIL,ISCR,NPRINT,VIN,VOUT,SMIN,SMOUT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C      This subroutine processes the output hydrograph and find components    C
C       of the water balance and hydrograph. The results go in "filename.osm" C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/PAR1/VL,FWIDTH,SWIDTH,SLENGTH
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      COMMON/WTGA1/WTD,PARW(4),PARK(2),hb,RAINL,ZW,TW,RVH
      COMMON/GRASSD/PART(3),SC,SS,VN1,VN2,VN,GSI,H,VLCM,POR,CI,ICO
      COMMON/GRASSD2/GSIMASS,GSOMASS,TTE,DEP,QSED(4),RS(3),DF(3),VM(3)
      COMMON/GRASSD3/SUSMASS,WEDGEMASS,NFUP
      COMMON/OLD/SEOLD,TOLD,XTOLD,YTOLD,CDEP,SE,VBTOLD,GSOLD
      COMMON/WQ1/VKD(10),CCP,CSAB(5),DGMRES0(10),DGMOL(10),DGFRAC(10,10)
      COMMON/IWQ2/NDGDAY,IDG,IWQ,IWQPRO,ICAT,IMOB
      COMMON/WQ3/DGKREF(10),FC,DGPIN(10),DGML,DGT(366),DGTHETA(366),DGLD(10),RF
      CHARACTER*75 DUMMY,DUMMY1
      CHARACTER*75 LISFIL(13)

C-------------- Summarize results in .osp from filename.ohy file -------------
      WRITE(10,*)'INPUTS'
      WRITE(10,*)'------'
c-----Close .OHY and .OSM and open again to output table (last NPRINT lines)--
      CLOSE(11)
      CLOSE(13)
      OPEN(11,FILE=LISFIL(6),STATUS='OLD')
      OPEN(13,FILE=LISFIL(8),STATUS='OLD')

c-----Echo hydrological inputs from .OHY and write into summary file .OSM-----
      VFSAREA=VL*FWIDTH
      IDX=0
      IDXR=0
      DO 5 WHILE (IDX.EQ.0)
            READ(11,'(A)')DUMMY
            IDX=INDEX(DUMMY,'OUTFLOW')
            IF(IDX.EQ.0) WRITE(10,*)DUMMY
            IDXR=INDEX(DUMMY,'Total rainfall')
            IF(IDXR.NE.0) THEN
                  BACKSPACE(11)
                  READ(11,100)DUMMY1,TRAIN
                  IF(TRAI.EQ.0) THEN
                        TRAI=TRAIN/1000.d0
                        PS=TRAI
                  ENDIF
             ENDIF
5     CONTINUE
      READ(11,'(A)')DUMMY
      READ(11,'(A)')DUMMY

      SUM0=0.D0
      SUM1=0.D0
      SUM2=0.D0
      SUM3=0.D0
      TIME0=0.D0
      BIG1=0.D0
      NZERO=1
      INI=0
      RAIN_E0=0.D0
      OUTF0=0.D0
      UPIN0=0.D0
      VF0=0.D0
      TINI=0.D0
      TEND=0.D0
      DO 10 I=1,NPRINT+1
            READ(11,*,END=16)TIME1,OUTF1,CUMFLOW,RAIN_E1,UPIN1,CUMIF,
     &             VF1
           TIMEINCR=TIME1-TIME0
            AREA0=TIMEINCR*(UPIN1+UPIN0)/2.D0
            SUM0=SUM0 + AREA0
            AREA1=TIMEINCR*(OUTF1+OUTF0)/2.D0
            SUM1=SUM1 + AREA1
            AREA2=TIMEINCR*(RAIN_E1+RAIN_E0)/2.D0
            SUM2=SUM2 + AREA2
            AREA3=TIMEINCR*(VF1+VF0)/2.D0
            SUM3=SUM3 + AREA3
            UPIN0=UPIN1
            OUTF0=OUTF1
            RAIN_E0=RAIN_E1
            VF0=VF1
            TIME0=TIME1
            C1= OUTF1
            BIG1=MAX(BIG1,C1)
            IF(C1.EQ.BIG1)THEN
                  TBIG=TIME1
                  QBIG=C1
            ENDIF
            IF(NZERO.EQ.1.AND.C1.GT.0.D0.AND.INI.EQ.0) THEN
                  TINI=TIME1
                  NZERO=0
                  INI=1
               ELSEIF(NZERO.EQ.0.AND.C1.EQ.0.D0) THEN
                  NZERO=1
                  TEND=TIME1
            ENDIF
            IF(C1.GT.0.D0)TEND=TIME1
10      CONTINUE
c-----Echo sediment inputs from .OG1 and write into summary file .OSM---------
12    READ(13,'(A)')DUMMY
      IDX=INDEX(DUMMY,'gsI')
      IF(IDX.NE.0) GOTO 16
      WRITE(10,*)DUMMY
      GOTO 12
16    CONTINUE

c-----Write summary outputs into summary file .OSM---------

      WRITE(10,*)'OUTPUTS'
      WRITE(10,*)'-------'
      WRITE(10,*)'Water balance'
      WRITE(10,*)'-------------'

C---------------Calculate totals for event in m3 ----------
      TOTRAIN=TRAI*VFSAREA
      WRITE(10,600)TOTRAIN
      Vout=SUM1
      Vie=SUM2*VFSAREA
      Vin=SUM0
      F=SUM3
      VF= F*VFSAREA

C---------------Water and sediment balance for the event----------

      WAT_IN=Vin+TOTRAIN
      VFE=WAT_IN - Vout
      if(AGA.LE.0.D0)THEN
          VF=0.d0
        ELSEIF(VF.LT.0.d0) then
          VF=VFE
      ENDIF
      WAT_OUT= Vout+VF
      WAT_BAL=WAT_IN-WAT_OUT
      IF (WAT_IN.GT.0.d0) THEN
              WAT_ERR=WAT_BAL/WAT_IN*100.d0
          ELSE
              WAT_ERR=0.0d0
      ENDIF
c-----Check water balance and adjust for leaching and pollutant calculations
      IF(DABS(WAT_ERR).GT.5.d0) THEN
          VF=VFE
          WAT_OUT= Vout+VF
          WAT_BAL=WAT_IN-WAT_OUT
          WAT_ERR=WAT_BAL/WAT_IN*100.d0
          F= VF/VFSAREA
          Z=F/DM
      ENDIF
c-----For shallow water table adjust rectangular wetting front of same depth for
c-----leaching and pollutant calculations -------------------------------------
      IF(WTD.NE.999.d0) DM=VF/(VFSAREA*Z)

c-(02/1999)--Write Filter performance summary output file .OSM -------------

      WRITE(10,300)Vin
      WRITE(10,150)Vout
      WRITE(10,400)VF
      WRITE(10,*)' '
      WRITE(10,*)'Hydrology'
      WRITE(10,*)'-----------'
      WRITE(10,700)TINI
      WRITE(10,800)TBIG,QBIG
      WRITE(10,900)TEND
      WRITE(10,*)' '
      WRITE(10,*)'Sediment'
      WRITE(10,*)'--------'
      TTE=1.D0
      IF(GSIMASS.GT.0.D0) TTE=(GSIMASS-GSOMASS)/GSIMASS
      FWID=FWIDTH*100.D0
      WRITE(10,1100)GSIMASS,GSIMASS*FWID
      WRITE(10,1200)GSOMASS,GSOMASS*FWID
      WRITE(10,1010)TTE*100

c-------Sediment wedge final shape ------------------------------

      WRITE(10,*)'                    ----------------'
      WRITE(10,1250)
      VLT=VLCM-XTOLD
      X1=YTOLD/SC
      WRITE(10,1300)YTOLD
      WRITE(10,1400)XTOLD
      WRITE(10,1500)VLT
      WRITE(10,1600)X1
      WRITE(10,1700)DEP

c--------Mass balance -------------------------------------------

      GAMMASB=(1.D0-POR)*PART(3)
c------------Bmass= buffer sediment in wedge + lower deposition
      IF(NFUP.EQ.0.AND.YTOLD.LT.H) THEN
c------------Triangular wedge
            BMASS=(DEP*VLT+(XTOLD+X1)*YTOLD*0.5D0)*GAMMASB
         ELSEIF(NFUP.EQ.0.AND.YTOLD.EQ.H) THEN
c------------Trapezoidal wedge
            X3=YTOLD/Se
            X2=XTOLD-X3
            TRBOT=X1+X2+X3
            TRTOP=X2
            BMASS=(DEP*VLT+(TRTOP+TRBOT)*YTOLD*0.5D0)*gammasb
         ELSE
c------------Filter filled up to top of vegetation
            BMASS=(H*VLCM+(H/SC)*H*0.5D0)*GAMMASB
      ENDIF
      IF (GSOMASS.GT.0.d0) THEN
            SED_ERR=(GSIMASS-BMASS-GSOMASS)/GSIMASS*100.D0
         ELSE
            SED_ERR=0.d0
      ENDIF
c      print*,'sed_err,wat_err=',(GSIMASS-(BMASS))/GSOMASS*100.D0,WAT_BAL/WAT_IN*100.d0
      WRITE(10,1800) SED_ERR
      IF (WAT_ERR.EQ.0.d0) THEN
            WRITE(10,1826) WAT_ERR
        ELSE
            WRITE(10,1825) WAT_ERR
      ENDIF
      IF(NFUP.EQ.1) THEN
            WRITE(10,1850)
        ELSEIF(YTOLD.EQ.H) THEN
            WRITE(10,1875)
      ENDIF

c-(02/1999)--Write Filter performance summary output file .OSP  -------------

      WRITE(15,*)' '
      WRITE(15,*)'      Summary of Buffer Performance Indicators:'
      WRITE(15,*)' '
      WRITE(15,1900)SWIDTH*SLENGTH
      WRITE(15,2000)SLENGTH
      WRITE(15,2100)SWIDTH
      WRITE(15,2200)VL
      WRITE(15,2225)FWIDTH
      WRITE(15,2250)VN1
      WRITE(15,*)' '
      WRITE(15,2300)VL/SLENGTH*100.d0
      WRITE(15,2400)TOTRAIN/(FWIDTH*VL)*1000.d0
      WRITE(15,2450)TOTRAIN
      WRITE(15,2500)Vin/(SWIDTH*SLENGTH)*1000.d0
      WRITE(15,2550)Vin
      WRITE(15,2600)Vout/(SWIDTH*SLENGTH+FWIDTH*VL)*1000.d0
      WRITE(15,2650)Vout
      WRITE(15,2675)VF
c--- rmc 04/20/03 --Fix for Vin(Q)=0 or Vout~0
      SMIN=GSIMASS*FWID/1000.d0
      SMOUT=GSOMASS*FWID/1000.d0
      if(Vin.le.0) then
          WRITE(15,2685)0.d0
          WRITE(15,*)' '
          WRITE(15,2700)SMIN
          WRITE(15,2800)0.d0
       else
          WRITE(15,2685)Vout/Vin
          WRITE(15,*)' '
          WRITE(15,2700)SMIN
          WRITE(15,2800)GSIMASS*FWID/(Vin*1000.d0)
      endif
      WRITE(15,2900)SMOUT
      if(Vout.le.0) then
            WRITE(15,3000)0.d0
            WRITE(15,3050)(GSIMASS-GSOMASS)*FWID/1000.d0
            WRITE(15,3075)0.d0
       else
            WRITE(15,3000)GSOMASS*FWID/(Vout*1000.d0)
            WRITE(15,3050)(GSIMASS-GSOMASS)*FWID/1000.d0
            WRITE(15,3075)SMOUT/SMIN
      endif
c--- rmc 04/20/03 --end of fix

      WRITE(15,*)' '
      WRITE(15,3100)VLT/100.d0
      WRITE(15,3200)XTOLD/100.d0
      IF(NFUP.EQ.1) THEN
            WRITE(15,1850)
        ELSEIF(YTOLD.EQ.H) THEN
            WRITE(15,1875)
      ENDIF

100   FORMAT(A31,F12.2)
150   FORMAT('Volume from outflow = ', E14.4,' m3')
200   FORMAT('Volume from i_e     = ', E14.4,' m3',F14.2,'%')
300   FORMAT('Volume from up-field= ', E14.4,' m3')
400   FORMAT('Volume infiltrated  = ', E14.4,' m3')
500   FORMAT(F8.2,3E12.4,I6)
600   FORMAT('Volume from rainfall= ', E14.4,' m3')
700   FORMAT('Time to beginning   = ', E14.4,' s')
800   FORMAT('Time and q at peak  = ', E14.4,' s ',E14.4,' m3/s')
900   FORMAT('Time to end runoff  = ', E14.4,' s')
1010  FORMAT('Trapping efficiency =      ', F5.1,' %')
1100  FORMAT('Sediment inflow     = ', E14.4,' g/cm',E14.4,' g')
1200  FORMAT('Sediment outflow    = ', E14.4,' g/cm',E14.4,' g')
1250  FORMAT('Sediment deposition :')
1300  FORMAT(8x,'- Sediment wedge depth,   Y(t) =',F8.2,' cm')
1400  FORMAT(8x,'- Sediment wedge length,  X2(t)=',F8.2,' cm')
1500  FORMAT(8x,'- Effective filter length,L(t) =',F8.2,' cm')
1600  FORMAT(8x,'- Sediment tail at field, X1(t)=',F8.2,' cm')
1700  FORMAT(8x,'- Sediment depth in low section=',F8.2,' cm')
1800  FORMAT(8x,'- Rough sediment balance error =',F8.2,' %')
1825  FORMAT(8x,'- Rough water balance error    =',F8.2,' %')
1826  FORMAT(8x,'- Rough water balance error    =',F8.2,' % (adjusted)')
1850  format(/,66('*'),/,'*  WARNING: Strip filled up!',37x,'*',
     &     /,66('*'))
1875  format(/,67('*'),/,'*  WARNING: Top of vegetation reached',
     &     ' - trapezoidal wedge started *',/,67('*'))
1900  FORMAT(1x,F10.1,' m^2 = Source Area (input)')
2000  FORMAT(1x,F10.2,' m   = Source Flow Length (input)')
2100  FORMAT(1x,F10.2,' m   = Source Area Width (input)')
2200  FORMAT(1x,F10.2,' m   = Filter Strip Length (input)')
2225  FORMAT(1x,F10.2,' m   = Filter Strip Width (input)')
2250  FORMAT(1x,F10.3,'     = Mean Filter Mannings Roughness (input)')
2300  FORMAT(1x,F10.2,' %   = ',
     1        'Ratio of Filter Length to Source Flow Length')
2400  FORMAT(1x,F10.3,' mm  = Total Rainfall')
2450  FORMAT(1x,F10.3,' m3  = Total Rainfall on Filter')
2500  FORMAT(1x,F10.3,' mm  = Total Runoff from Source (mm depth over',
     &    ' Source Area)')
2550  FORMAT(1x,F10.3,' m3  = Total Runoff from Source')
2600  FORMAT(1x,F10.3,' mm  = Total Runoff out from Filter (mm depth',
     &    ' over Source+Filter)')
2650  FORMAT(1x,F10.3,' m3  = Total Runoff out from Filter')
2675  FORMAT(1x,F10.3,' m3  = Total Infiltration in Filter')
2685  FORMAT(1x,F10.3,'     = Runoff Delivery Ratio')
2700  FORMAT(F12.3,' kg = Mass Sediment Input to Filter')
2800  FORMAT(F12.3,' g/L= Concentration Sediment in Runoff from',
     &    ' source Area')
2900  FORMAT(F12.3,' kg = Mass Sediment Output from Filter')
3000  FORMAT(F12.3,' g/L= Concentration Sediment in Runoff exiting',
     &   ' the Filter')
3050  FORMAT(F12.3,' kg = Mass Sediment retained in Filter')
3075  FORMAT(2x,F10.3,'    = Sediment Delivery Ratio')
3100  FORMAT(1x,F10.2,' m   = Effective Filter Length')
3200  FORMAT(1x,F10.2,' m   = Wedge Distance')
3500  FORMAT(26('-'))
3525  FORMAT(A40,I1,')')
3550  FORMAT(50('-'))

      RETURN
      END
