      SUBROUTINE CDE(JJ,VIN,VOUT,SMIN,SMOUT,TIME)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C      The subroutine solves Huang and van Genuchten (1995, Soil Sci. 159(4): C
C      217-223) CDE analytical solution with retardation (sorption),          C
C      initially unsaturated soil with piston flow folowing Green-Ampt, third C
C      (flux) BC on top and 0 flux on the bottom, and 0 initial concentration C
C      in the soil. Although solution can handle 2 pulses (Co for 0<t<to, 0   C
C      for t>to) only one is used here. Note that there is an error in Eq.22c C
C      of the original paper that is corrected here.                          C
C                                                                             C
C    Variables used:                                                          C
C      VIN,VOUT=runoff volume in and out of the VFS (m3), also Qi and Qo.     C
C      SMIN,SMOUT= Sediment in and out of the VFS (Kg), also Ei and Eo.       C
C      F= Infiltration during the event (m) = F0 for single pulse             C
C      PS= Ranfall during the event (m)                                       C
C      z= wetting front depth at the end of event (m)= Z1 for single pulse    C
C      SAREA,VFSAREA (m2): Source (field) and filter areas                    C
C      VIN,VOUT(m3)= VFS runoff inflow and outflow volume Ohydrograph areas)  C
C      SMIN,SMOUT(Kg)= VFS sediment inflow and outflow                        C
C      WAT_IN(m3)= Runoff inflow (VIN) + total rain on the filter             C
C      VF,VFL= Total Infiltration in Filter (m3 and L, respectively)          C
C      OS= soil saturated water content (porosity)                            C
C      DGPIN(mg/m2) = incoming pesticide in filter (per m2 of SOURCE area)    C
C      DGCIN (mg/L) = Chemical concentration of runoff                        C
C      DGCIN2 (mg/L) = Average chemical concentration of infiltrating water   C
C      DGMfFd,DGMfFp,DGMfF(mg)= leached pesticide (dissolved,adsorbed,total)  C
C      DGMmld,DGMmlp,DGMml(mg)= mixing layer pest. (dissolved,adsorbed,total) C
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/PAR1/VL,FWIDTH,SWIDTH,SLENGTH
      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      COMMON/GRASSD2/GSIMASS,GSOMASS,TTE,DEP,QSED(4),RS(3),DF(3),VM(3)
      COMMON/WTGA1/WTD,PARW(4),PARK(2),hb,RAINL,ZW,TW,RVH
      COMMON/WQ1/VKD(10),CCP,CSAB(5),DGMRES0(10),DGMOL(10),DGFRAC(10,10)
      COMMON/IWQ2/NDGDAY,IDG,IWQ,IWQPRO,ICAT,IMOB
      COMMON/WQ3/DGKREF(10),FC,DGPIN(10),DGML,DGT(366),DGTHETA(366),DGLD(10),RF
      COMMON/CINT/XI(20,20),W(20,20)
      COMMON/CDE1/DGMfFd(10),DGMfFp(10),DGMfF(10),DGMmld(10),DGMmlp(10),DGMml(10)

      MAXIT=100
c------Infiltration and water balance after event
      VFSAREA=VL*FWIDTH
      SAREA=SWIDTH*SLENGTH
      VFL=F*1000.d0*VFSAREA
      IF(DM.LE.0.0001d0.AND.WTD.EQ.999.D0) then
c---------For initial saturated soil conditions, f=Ks and z=Ks.time
          Z1=AGA*TIME/OS
        ELSE
c---------For initial unsaturated soil, G-A wetting front depth z=F/M
          Z1=F/DM
      ENDIF
      F0=F
      TOTRAIN=PS*VFSAREA
      TOTRAINL=TOTRAIN*1000.d0
      VinL=Vin*1000.d0
      WAT_INL=VinL+TOTRAINL

c------Average chemical concentration (DGCIN2) of infiltrating water. The incoming
c------mass is DGPIN plus remobilized residues from prev. event, DGPIN=DGPIN+DGM,FRES
c------(see inputs.f). The mass is partitioned in dissolved (mdi)+sorbed mass (mpi)
      Ei=SMIN
      IF(VinL.GT.0.d0) THEN
           DGPINd=DGPIN(JJ)*VinL/(VinL+Ei*VKD(JJ))
           DGPINp=DGPIN(JJ)*Ei*VKD(JJ)/(VinL+Ei*VKD(JJ))
           DGCIN=DGPINd/VinL
        ELSE
            DGPINd=0.d0
            DGPINp=0.d0
            DGPIN=0.d0
      ENDIF
      DGSIN=VKD(JJ)*DGCIN

C-------check distribution of incoming pesticide mass ------
      DGPIN_check=DGCIN*VinL+DGSIN*Ei
      IF(DGPIN_check-DGPIN(JJ).GT.0.000000000001d0) THEN
        print*,'WARNING: incoming mass distribution wrong:'
        print*,'         DGPIN_check,DGPIN=',DGPIN_check,DGPIN(JJ)
      ENDIF

c------A) Average infiltration resident concentration. The concentration is the
c------ runin dissolved mass over total surface water, wat_in = runoff+rain (L).
c------ The average residuent concentration is mass conservative and used here.
      DGCIN2=DGPINd*SAREA/WAT_INL
c------B) Average infiltration flux concentration. The concentration is the mass
c------trapped in filter (mf) less the mass in sediment (mfsed) over the water
c------infiltrated during the event (VF_L). This average flux concentrations
c------produces mass blance errors and it is not used here.
C-------First calculate pesticide trapping efficiency based on:--------------
C-------IWQPRO=1 - Sabbagh et al. (2009) semiempirical eq. ---------------
C-------IWQPRO=2 - Refitted Sabbagh et al. eq. w/user supplies coeff -----
C-------IWQPRO=3 - Based on Muñoz-Carpena et al. (2015) mechanistic eq. --
C-------IWQPRO=4 - Chen et al. (2017) empirical eq. ----------------------
c---- calculate factors for pesticide trapping equations ----
C      PDSED=TTE*100.D0
C      IF(WAT_INL.gt.0.d0)THEN
C          PDQ=VFL/WAT_INL*100.D0
C        ELSE
C          PDQ=100.D0
C      ENDIF
C      PDSED=TTE*100.D0
C      SELECT CASE (IWQPRO)
C        CASE (1)
C-------Original Sabbagh eq.
C           DELTAP=CSAB(1)+CSAB(2)*PDQ+CSAB(3)*PDSED+CSAB(4)*
C     &          DLOG(FPH+1.D0)+CSAB(5)*CCP
C         CASE (2)
C--------Sabbagh refit eq.
C           DELTAP=CSAB(1)+CSAB(2)*PDQ+CSAB(3)*PDSED+CSAB(4)*
C     &          DLOG(FPH+1.D0)+CSAB(5)*CCP
C         CASE (3)
C--------Mechanistic mass balance eq.
C           DELTAP=100.D0*MIN(VINL+VKD(JJ)*SMIN,PDSED/100*SMIN*VKD(JJ)+
C     &          PDQ/100*VINL)/(VINL+VKD(JJ)*SMIN)
C         CASE DEFAULT
C--------Chen eq. (California)
C           DELTAP=101.D0-(8.06D0-0.07D0*PDQ+0.02D0*PDSED+0.05D0*CCP-
C     &          2.17D0*ICAT+0.02D0*PDQ*ICAT-0.0003D0*PDQ*PDSED)**2.D0
C      END SELECT
C      IF (DELTAP.GT.100.D0.OR.PDQ.GE.100.d0) DELTAP=100.D0
C      IF (PDQ.EQ.0.d0.AND.PDSED.EQ.0.d0) DELTAP=0.D0
C      IF (DELTAP.LT.0.d0) DELTAP=0.D0
C      DGCIN2=(DGPIN*DELTAP/100.d0*SAREA-DGSIN*Ei*PDSED/100.d0)/VFL

C-------Prepare inputs for CDE transport solution ------
      C1=DGCIN2
      C2=0.d0
      DGROB=(1.d0-OS)*2.65d0
      RF=1.d0+VKD(JJ)*DGROB/OS
      TT=F/OS
      TT0=F0/OS
      Zf=Z1
      ZD=0.D0
      DGHALFJJ=DLOG(2.D0)/DGKREF(JJ)
 
C-------Output CDE transport inputs-------------------------------
      IF(JJ.EQ.1)WRITE(18,300)'Soil leaching and mixing layer calculations (CDE)'
      IF(IWQ.GT.1) WRITE(18,350)'COMPOUND',JJ,':'
      WRITE(18,*)'Huang & van Genuchten (1995) CDE analytical solution'
      WRITE(18,*)'Single pulse with C1=dissolved in mass/(runoff+rain)'
      WRITE(18,800)'Wetting front depth (z)=',Z1,'m'

      WRITE(18,800)'Incoming runoff volume (Vi)=',VinL,'L'
      WRITE(18,800)'Event rainfall (P)=',TOTRAINL,'L'
      WRITE(18,800)'Total infiltration (F)=',VFL,'L'
      WRITE(18,800)'Incoming sediment (Ei)=',Ei,'Kg'
      WRITE(18,800)'Incoming pesticide mass (mi)=',DGPIN(JJ),'mg/m2'
      WRITE(18,800)'Incoming dissolved mass (mdi)=',DGPINd,'mg/m2'
      WRITE(18,800)'Incoming sed-sorbed mass (mpi)=',DGPINp,'mg/m2'
      WRITE(18,800)'Incoming dissolved mass (mdi)=',DGPINd*SAREA,'mg'
      WRITE(18,800)'Incoming sed-sorbed mass (mpi)=',DGPINp*SAREA,'mg'
      WRITE(18,800)'Pulse concentration (C1)=',C1,'mg/L'
      WRITE(18,800)'Partition coefficient (Kd)=',VKD(JJ),'L/Kg'
      WRITE(18,800)'Dispersion length (l)=',DGLD(JJ),'m'
      WRITE(18,800)'Pesticide half-life (Ln2/Kref)=',DGHALFJJ,'days'
      WRITE(18,800)'Retardation factor (R)=',RF,'(-)'
      WRITE(18,800)'Transformed time (T)=',TT,'(-)'
      WRITE(18,*)
      WRITE(18,*)'     Z(m)      C(mg/L)      S(mg/mg)'
      WRITE(18,*)'-------------------------------------'

c-------Huang and Van Genuchten (1995) CDE analytical solution
      do while (ZD.le.Zf)
        C=CONC(JJ,C1,C2,TT,TT0,ZD)
        WRITE(18,500)ZD,C,VKD(JJ)*C
        ZD=ZD+0.01D0
      end do
      C=0.d0
      WRITE(18,500)ZD,C,VKD(JJ)*C

c-----Integrate mass in the soil profile for the event
      call qgausscde(JJ,C1,C2,TT,TT0,0.d0,Zf,20,CONCINTG)
      DGMfFd(JJ)=CONCINTG*OS*1000.d0*VFSAREA
      DGMfFp(JJ)=CONCINTG*VKD(JJ)*DGROB*1000.d0*VFSAREA
      DGMfF(JJ)=DGMfFd(JJ)+DGMfFp(JJ)
      IF(DGMfF(JJ).EQ.0.d0) THEN
          DGMfFdPCT= 0.d0
          DGMfFpPCT= 0.d0
          DGMfFPCT= 0.d0
        ELSE
          DGMfFdPCT=DGMfFd(JJ)/DGMfF(JJ)*100.d0
          DGMfFpPCT=DGMfFp(JJ)/DGMfF(JJ)*100.d0
          DGMfFPCT=DGMfF(JJ)/DGMfF(JJ)*100.d0
        END IF
      WRITE(18,*)
      WRITE(18,475)'Soil profile total mass (mfF mg)=',DGMfF(JJ),
     &  DGMfFPCT
      WRITE(18,475)'Soil profile dissolved mass (mfFd mg)=',DGMfFd(JJ),
     &  DGMfFdPCT
      WRITE(18,475)'Soil profile sorbed mass (mfFp mg)=',DGMfFp(JJ),
     &  DGMfFpPCT

c-----Integrate mass in mixing layer for the event
      call qgausscde(JJ,C1,C2,TT,TT0,0.d0,DGML/100.d0,20,CONCINTG1)
      IF(CONCINTG1.GT.CONCINTG) CONCINTG1=CONCINTG
      DGMmld(JJ)=CONCINTG1*OS*1000.d0*VFSAREA
      DGMmlp(JJ)=CONCINTG1*VKD(JJ)*DGROB*1000.d0*VFSAREA
      DGMml(JJ)= DGMmld(JJ)+DGMmlp(JJ)
      IF(DGMfF(JJ).EQ.0.d0) THEN
          DGMmldPCT= 0.d0
          DGMmlpPCT= 0.d0
          DGMmlPCT= 0.d0
        ELSE
          DGMmldPCT=DGMmld(JJ)/DGMfF(JJ)*100.d0
          DGMmlpPCT=DGMmlp(JJ)/DGMfF(JJ)*100.d0
          DGMmlPCT=DGMml(JJ)/DGMfF(JJ)*100.d0
      END IF
      WRITE(18,475)'Mixing layer total mass (mfml mg)=',DGMml(JJ),
     &  DGMmlPCT
      WRITE(18,475)'Mixing layer dissolved mass (mfmld mg)=',
     &  DGMmld(JJ),DGMmldPCT
      WRITE(18,475)'Mixing layer sorbed mass (mfmlp mg)=',DGMmlp(JJ),
     &  DGMmlpPCT

300   FORMAT(/,65('-'),/,A50,/,65('-'))
350   FORMAT(/,A8,I3,A1)
400   FORMAT(A41,3E12.4)
450   FORMAT(A20,F12.2)
475   FORMAT(A41,f12.4,F8.2,'%')
500   FORMAT(6F12.4)
600   FORMAT(A40,2F12.4)
800   FORMAT(A31,F12.4,A11)
c800   FORMAT(A31,F12.6,A20)



      return
      end

      FUNCTION CONC(JJ,C1,C2,TT,TT0,ZD)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Auxiliary function of van Genuchten and Alves (1980) A(x,t)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      if(TT.le.TT0) THEN
          CONC=C2+(C1-C2)*AUX(JJ,ZD,TT)
        ELSE
          CONC=C2+(C1-C2)*(AUX(JJ,ZD,TT)-AUX(JJ,ZD,TT-TT0))
      ENDIF

      return
      END


      FUNCTION AUX(JJ,ZD,TTF)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Auxiliary function of van Genuchten and Alves (1980) A(x,t)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/WQ3/DGKREF(10),FC,DGPIN(10),DGML,DGT(366),DGTHETA(366),DGLD(10),RF

      PI=3.14159265358979323846264338327950288419716939937510582d0
      D1=(RF*ZD-TTF)/(2.d0*dsqrt(RF*DGLD(JJ)*TTF))
      D2=(RF*ZD+TTF)/(2.d0*dsqrt(RF*DGLD(JJ)*TTF))
      A1= 0.5d0*erfcc(D1)
      A2= -0.5d0*(1.d0+ZD/DGLD(JJ)+TTF/(DGLD(JJ)*RF))*dexp(ZD/DGLD(JJ))*erfcc(D2)
      A3= dsqrt(TTF/(PI*DGLD(JJ)*RF))*dexp(-D1**2)
      AUX= A1 + A2 + A3

      return
      END

      FUNCTION erfcc(x)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Returns the complementary error function erfc(x) with fractional error    C
C   everywhere less than 1.2e−7, based on Chebyshev fitting.                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      x1=dabs(x)
      t=1.d0/(1.d0+0.5d0*x1)
      erfcc=t*dexp(-x1*x1-1.26551223d0+t*(1.00002368d0+t*(.37409196d0+t*
     & (.09678418d0+t*(-.18628806d0+t*(.27886807d0+t*(-1.13520398d0+t*
     & (1.48851587d0+t*(-.82215223d0+t*.17087277d0)))))))))
      if (x.lt.0.d0) erfcc=2.d0-erfcc
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #>0K!.

      FUNCTION erf(x)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   Returns the complementary error function erfc(x) with fractional error    C
C   everywhere less than 1.2e−7, based on Chebyshev fitting.                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      erf=1.d0-erfcc(x)

      return
      END

      subroutine qgausscde(JJ,C1,C2,TT,TT0,a,b,ngl,sst)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Numerical integration of a generic function using Gauss-Legendre quadrature
C   a,b: limits of integration
C   ngl: orger of the gauss-Legendre quadrature (up to 20)
C   sst: value of integral
C   Uses external function CONC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit double precision (a-h,o-z)
      COMMON/CINT/XI(20,20),W(20,20)
      external CONC

      xm=0.5d0*(b+a)
      xr=0.5d0*(b-a)
      sst=0.d0
      do 11 j=1,ngl
        dx=xr*xi(j,ngl)
        sst=sst+w(j,ngl)*CONC(JJ,C1,C2,TT,TT0,xm+dx)
11     continue
      sst=xr*sst

      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software #>0K!.
