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
      COMMON/WQ1/VKD(10),VKF(10),VKN(10),CCP,CSAB(5),DGMRES0(10),
     &           DGMOL(10),DGFRAC(10,10)
      COMMON/IWQ2/NDGDAY,IDG,IWQ,IWQPRO,ICAT,IMOB
      COMMON/WQ3/DGKREF(10),FC,DGPIN(10),DGML,DGT(366),DGTHETA(366),DGLD(10),RF,DGMPI(10)
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
c---------Phase partitioning of incoming pesticide mass:
c---------Linear (IKD=0,1): s=Kd*C; Freundlich (IKD=2): s=Kf*C^N
c---------For Freundlich use effective linearized Kd at inlet conc (v4.6.2)
      IF(VinL.GT.0.d0) THEN
c--------Effective linear Kd at inlet concentration for phase partitioning
c--------For Freundlich: Kd_eff = Kf*C^(N-1); for linear: Kd_eff = Kd
            IF(VKN(JJ).EQ.1.D0) THEN
               VKDEFF=VKD(JJ)
               DGPINd=DGPIN(JJ)*VinL/(VinL+Ei*VKDEFF)
               DGPINp=DGPIN(JJ)-DGPINd
               DGCIN=DGPINd/VinL
            ELSE
               IF(DGPIN(JJ).LE.0.D0) THEN
                  DGCIN=0.D0
               ELSEIF(Ei.LE.0.D0) THEN
                  DGCIN=DGPIN(JJ)/VinL
               ELSE
c-----------Solve exact Freundlich inlet partitioning: mi = C*Vi + Kf*C^N*Ei
                  CLO=0.D0
                  CHI=DMAX1(DGPIN(JJ)/VinL,1.D-10)
                  DO 71 IBR=1,60
                    FHI=CHI*VinL+VKF(JJ)*CHI**VKN(JJ)*Ei-DGPIN(JJ)
                    IF(FHI.GE.0.D0) GOTO 72
                    CHI=CHI*2.D0
 71               CONTINUE
 72               CONTINUE
                  DO 73 IBR=1,80
                    CMID=0.5D0*(CLO+CHI)
                    FMID=CMID*VinL+VKF(JJ)*CMID**VKN(JJ)*Ei-DGPIN(JJ)
                    IF(FMID.GT.0.D0) THEN
                      CHI=CMID
                    ELSE
                      CLO=CMID
                    ENDIF
                    IF((CHI-CLO)/DMAX1(CHI,1.D-10).LT.1.D-10) GOTO 74
 73               CONTINUE
 74               CONTINUE
                  DGCIN=0.5D0*(CLO+CHI)
               ENDIF
               IF(DGCIN.LE.0.D0) THEN
                  VKDEFF=0.D0
                  DGPINd=0.D0
                  DGPINp=0.D0
               ELSE
                  VKDEFF=VKF(JJ)*DGCIN**(VKN(JJ)-1.D0)
                  DGPINd=DGCIN*VinL
                  DGPINp=DGPIN(JJ)-DGPINd
               ENDIF
            ENDIF
         ELSE
             DGPINd=0.d0
             DGPINp=0.d0
             DGPIN=0.d0
            VKDEFF=VKD(JJ)
      ENDIF
c-----Incoming sediment-sorbed mass (mpi, mg/m2): used by WQPEST sediment trapping
      DGMPI(JJ)=DGPINp
c-----Sorbed concentration on incoming sediment: s = Kf*C^N (or Kd*C)
      IF(VKN(JJ).EQ.1.D0) THEN
         DGSIN=VKD(JJ)*DGCIN
      ELSE
         IF(DGCIN.GT.0.D0) THEN
            DGSIN=VKF(JJ)*DGCIN**VKN(JJ)
         ELSE
            DGSIN=0.D0
         ENDIF
      ENDIF

C-------check distribution of incoming pesticide mass ------
      DGPIN_check=DGCIN*VinL+DGSIN*Ei
      IF(DGPIN(JJ).GT.0.d0 .AND.
     &   (DGPIN_check-DGPIN(JJ))/DGPIN(JJ).GT.0.0001d0) THEN
        print*,'WARNING: incoming mass distribution wrong:'
        print*,'         DGPIN_check,DGPIN=',DGPIN_check,DGPIN(JJ)
      ENDIF

c------A) Average infiltration resident concentration. The concentration is the
c------ runin dissolved mass over total surface water, wat_in = runoff+rain (L).
c------ The average resident concentration is mass conservative (USED HERE!).
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
      TT=F/OS
      TT0=F0/OS
c---------Retardation factor R for CDE leaching solution (Eq. 3):
c---------Linear:     R = 1 + Kd*rho_b/theta_s (constant)
c---------Freundlich: R = 1 + N*Kf*C^(N-1)*rho_b/theta_s (concentration-
c---------            dependent, linearized at C1 via fixed-point iteration)
c---------            Munoz-Carpena et al. SETAC 2026 (v4.6.2)
c-----Retardation factor R: linear or Freundlich fixed-point iteration
c-----RF is in COMMON/WQ3 so CONC/AUX picks up each update automatically
      IF(VKN(JJ).EQ.1.D0) THEN
         RF=1.d0+VKD(JJ)*DGROB/OS
      ELSE
c--------Freundlich: R^(k+1) = 1 + N*Kf*[C1(R^k)]^(N-1)*rho_b/theta_s
         IF(C1.LE.0.D0) THEN
            RF=1.D0
         ELSE
c-----------Secant (chord) retardation for the loading shock front:
c-----------R = 1 + (rho/theta)*[q(C)/C] = 1 + (rho/theta)*Kf*C^(N-1).
c-----------The secant slope (no N factor) places the self-sharpening front by
c-----------mass balance and equals the equilibrium partition coeff, so transport
c-----------and sorbed-mass use one consistent coefficient (reduces to Kd at N=1).
            C1SURF=DMAX1(C1,1.D-10)
            RF=1.d0+VKF(JJ)*C1SURF**(VKN(JJ)-1.D0)*DGROB/OS
            DO 135 INRNR=1,50
               RFOLD=RF
               C1SURF=CONC(JJ,C1,C2,TT,TT0,0.D0)
               IF(C1SURF.LE.0.D0) C1SURF=1.D-10
               RF=1.D0+VKF(JJ)*C1SURF**(VKN(JJ)-1.D0)*DGROB/OS
               IF(DABS(RF-RFOLD)/DMAX1(RFOLD,1.D-10).LT.1.D-6) GOTO 136
135         CONTINUE
136         CONTINUE
         ENDIF
      ENDIF
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
      IF(VKN(JJ).EQ.1.D0) THEN
c--------Linear isotherm echo (v4.6.2)
         WRITE(18,800)'Linear sorption coeff. (Kd)=',VKD(JJ),'L/Kg'
      ELSE
c--------Non-linear Freundlich: print Kf and N on one line (v4.6.2)
         WRITE(18,811)'Non-linear Freundlich (Kf,Nf)=',VKF(JJ),
     &        '  L^N/Kg,',VKN(JJ),'  (-)'
      ENDIF
      WRITE(18,800)'Dispersion length (l)=',DGLD(JJ),'m'
      WRITE(18,800)'Pesticide half-life (Ln2/Kref)=',DGHALFJJ,'days'
      WRITE(18,800)'Retardation factor (R)=',RF,'(-)'
      WRITE(18,800)'Transformed time (T)=',TT,'(-)'
      WRITE(18,*)
      IF(VKN(JJ).EQ.1.D0) THEN
        WRITE(18,*)'     Z(m)   C_dis(mg/L)   S_lin(mg/mg)'
      ELSE
        WRITE(18,*)'     Z(m)   C_dis(mg/L)  S_Fre(mg/mg)'
      ENDIF
      WRITE(18,*)'-------------------------------------'

c---------Sorbed concentration output: s=Kd*C (linear) or s=Kf*C^N (Freundlich)
c-------Huang and Van Genuchten (1995) CDE analytical solution
      do while (ZD.le.Zf)
        C=CONC(JJ,C1,C2,TT,TT0,ZD)
        IF(VKN(JJ).EQ.1.D0) THEN
           WRITE(18,500)ZD,C,VKD(JJ)*C
        ELSE
           IF(C.GT.0.D0) THEN
              WRITE(18,500)ZD,C,VKF(JJ)*C**VKN(JJ)
           ELSE
              WRITE(18,500)ZD,C,0.D0
           ENDIF
        ENDIF
        ZD=ZD+0.01D0
      end do
      C=0.d0
      WRITE(18,500)ZD,C,0.D0

c---------Adsorbed mass in soil profile: integrate Kf*C^N*rho_b (Freundlich)
c---------or Kd*C*rho_b (linear) over depth (v4.6.2)
c-----Integrate mass in the soil profile for the event
      call qgausscde(JJ,C1,C2,TT,TT0,0.d0,Zf,20,CONCINTG)
      DGMfFd(JJ)=CONCINTG*OS*1000.d0*VFSAREA
c-----Sorbed profile mass consistent with the (secant) retardation used in
c-----transport: mfFp=(R-1)*mfFd. This conserves mass (mfF=R*mfFd=infiltrated
c-----dissolved mass) for any isotherm and equals CONCINTG*Kd*rho_b when N=1.
      DGMfFp(JJ)=(RF-1.D0)*DGMfFd(JJ)
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
c-----Mixing-layer sorbed mass, same (R-1)*dissolved consistency as the profile.
      DGMmlp(JJ)=(RF-1.D0)*DGMmld(JJ)
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
550   FORMAT(6E12.4)
600   FORMAT(A40,2F12.4)
800   FORMAT(A31,F12.4,A11)
811   FORMAT(A31,F12.4,A9,F8.6,A6)
c800   FORMAT(A31,F12.6,A20)



      return
      end

      FUNCTION CONC(JJ,C1,C2,TT,TT0,ZD)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Auxiliary function of van Genuchten and Alves (1980) A(x,t)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

c---------Heaviside unit function U(T-T0) of Eq. 20 (Huang & van Genuchten,
c---------1995): applies superposition for 2-pulse input. For single pulse
c---------with C2=0 and T>T0 this gives C=C1*(A(z,T)-A(z,T-T0)).
c---------Assumes zero initial pesticide concentration in the profile C(z,0)=0
c---------so G_i and E_i terms of Eq. 22b,c vanish (v4.6.2)
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
      COMMON/WQ3/DGKREF(10),FC,DGPIN(10),DGML,DGT(366),DGTHETA(366),DGLD(10),RF,DGMPI(10)

      PI=3.14159265358979323846264338327950288419716939937510582d0
      D1=(RF*ZD-TTF)/(2.d0*dsqrt(RF*DGLD(JJ)*TTF))
      D2=(RF*ZD+TTF)/(2.d0*dsqrt(RF*DGLD(JJ)*TTF))
      A1= 0.5d0*erfcc(D1)
      A2= -0.5d0*(1.d0+ZD/DGLD(JJ)+TTF/(DGLD(JJ)*RF))*dexp(ZD/DGLD(JJ))*erfcc(D2)
      A3= dsqrt(TTF/(PI*DGLD(JJ)*RF))*dexp(-D1**2)
      AUX= A1 + A2 + A3

      return
      END

      FUNCTION SCONC(JJ,C1,C2,TT,TT0,ZD)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Auxiliary sorbed concentration integrand for Freundlich mass integration    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/WQ1/VKD(10),VKF(10),VKN(10),CCP,CSAB(5),DGMRES0(10),
     &           DGMOL(10),DGFRAC(10,10)

      CVAL=CONC(JJ,C1,C2,TT,TT0,ZD)
      IF(CVAL.GT.0.D0) THEN
         SCONC=CVAL**VKN(JJ)
      ELSE
         SCONC=0.D0
      ENDIF

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

      subroutine qgausssorb(JJ,C1,C2,TT,TT0,a,b,ngl,sst)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Numerical integration of the Freundlich sorbed integrand C^N over depth     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit double precision (a-h,o-z)
      COMMON/CINT/XI(20,20),W(20,20)
      external SCONC

      xm=0.5d0*(b+a)
      xr=0.5d0*(b-a)
      sst=0.d0
      do 12 j=1,ngl
        dx=xr*xi(j,ngl)
        sst=sst+w(j,ngl)*SCONC(JJ,C1,C2,TT,TT0,xm+dx)
 12   continue
      sst=xr*sst

      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software #>0K!.
