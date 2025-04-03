      SUBROUTINE INPUTS(N,NBAND,NRAIN,RAIN,NBCROFF,BCROFF,TE,QMAX,
     & PGPAR,NCHK,LISFIL,ISCR,QSUM0,RMMH,BCROPEAK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Read data from input files (ikw,iso,igr,isd,irn,iro(,iwq)) in free  C
C  format.                                                             C
C  In file IKW, calculate the following parameters:                    C
C   1- N, NBAND,NELEM                                                  C
C   2- Maximum flow rate,depth at steady-state condition(QMAX,HMAX)    C
C   3- Celerity of the wave (C)                                        C
C   4- Courant time step (DTC)                                         C
C   5- Froude number (FR)                                              C
C   6- Kinematic flow number (FK)                                      C
C   7- Henderson's time to equilibrium (TE)                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (MAXEQN=1001,MAXBND=40)
      CHARACTER*30 PLABEL
      CHARACTER*75 LISFIL(13)
      CHARACTER*1 CWQ
      CHARACTER*10 CWTD,CIDG,CIMOB
      CHARACTER*200 line

      COMMON/GA1/PS,PSOLD,PST,F,RO,TP,TPP,FPI,STO,CU,AGA,BGA,DM,SM,Z,OS
      COMMON/WTGA1/WTD,PARW(4),PARK(2),hb,RAINL,ZW,TW,RVH
      COMMON/WTGA2/ITHETATYPE,IKUNSTYPE,ITWBC
      COMMON/PAR/QK(1001),R,THETAW,DX,DT,NDT,NELEM,MAXITER,NPOL,IOUT,NL
      COMMON/PAR1/VL,FWIDTH,SWIDTH,SLENGTH
      COMMON/GRASSD/PART(3),SC,SS,VN1,VN2,VN,GSI,H,VLCM,POR,CI,ICO
      COMMON/WQ1/VKD(10),CCP,CSAB(5),DGMRES0(10),DGMOL(10),DGFRAC(10,10)
      COMMON/IWQ2/NDGDAY,IDG,IWQ,IWQPRO,ICAT,IMOB
      COMMON/WQ3/DGKREF(10),FC,DGPIN(10),DGML,DGT(366),DGTHETA(366),DGLD(10),RF    
      
      DIMENSION BCROFF(200,2),RAIN(200,2)
      DIMENSION NODEP(MAXEQN),RNA(MAXEQN),SOA(MAXEQN),SX(MAXEQN)
      DIMENSION PGPAR(4),CSAB2(5),VKOC(10),DGHALF(10),test_array(100)

C-----(.IKW file) Read in main parameters of the program--------------

      READ(1,'(A30)')PLABEL
      READ(1,*)FWIDTH
      READ(1,*)VL,N,THETAW,CR,MAXITER,NPOL,IELOUT,KPG
      IOUT=0
      IF(ISCR.EQ.0) THEN
        WRITE(*,145)PLABEL
        WRITE(*,'(" ... Reading inputs from: ",A45)')LISFIL(1)
      ENDIF

c-----When N=-1 calculate N based on VL so that 11>=N>=101 for VL=1-100m
      IF(N.EQ.-1) N=2*nint(0.5d0*((90.d0/99.d0*(VL-1.d0)+11.d0)+1.d0))-1
 
C--[EVR-1998]-Check if N is compatible with type of shape funcion -----
      L=N-1
      M=NPOL-1
      IF (MOD (L,M).NE.0) N = N+M-MOD(L,M)
      SEPN = VL/(N-1)

c---------- Read surface properties of the filter -------------------

      READ(1,*)NPROP
      DO 5 IPROP=1,NPROP
        READ(1,*)SX(IPROP),RNA(IPROP),SOA(IPROP)
5     CONTINUE

c--[06/2008]--Read IWQ flag (0= no; >1= yes where IWQ sets the number of species)-
      IWQ=0
      read(1,'(i10)',IOSTAT=ios,ERR=8,END=8) IWQ
      IF(IOS.NE.0.OR.IWQ.LT.0)THEN
            IWQ=0
	END IF

C--[EVR-1998]-Assign nodes to the X-values where filter changes -----

8     XSEG=0.d0
      J = 1
      DO 80 I= 1, N
        DO WHILE (XSEG.GT.SX(J).AND.I.NE.N)
            J= J+1
            NODEP(J-1)=I
        END DO
        XSEG = XSEG + SEPN
80    CONTINUE
      NODEP(NPROP)=N

C -----Calculate alpha for Mannings equation------------------

      SMALLQK=1000.D0
      BIGQK=0.D0
      SOAVG=0.D0
      RNAVG=0.D0
      DO 15 I=1,N
        DO 10 IPROP=1,NPROP
            IF(I.LE.NODEP(IPROP))THEN
                RN=RNA(IPROP)
                SO=SOA(IPROP)
                GOTO 12
            ENDIF
10      CONTINUE
12      SOAVG=SOAVG+SO
        RNAVG=RNAVG+RN
        QK(I) = SO**0.5D0/RN
        BIGQK=DMAX1(BIGQK,QK(I))
        SMALLQK=DMIN1(SMALLQK,QK(I))
        IF (QK(I).eq.BIGQK)then
            nbig=i
           ELSEIF(QK(I).eq.SMALLQK)then
            nsmall=1
        ENDIF
15    CONTINUE

c-------Filter main slope and roughness for sediment calculations ------

      SC=SOAVG/N
      VN1=RNAVG/N

c-----(.IRN file) Read rainfall distribution [time[s]=RAIN(I,1);i[m/s]=RAIN(I,2)
      IF(ISCR.EQ.0) THEN
        WRITE(*,'(" ... Reading inputs from: ",A45)')LISFIL(2)
      ENDIF
      READ(2,*)NRAIN, RPEAK
      IF(NRAIN.GT.200) THEN
        WRITE(*,110)NRAIN
110     format(/,71('*'),/,'*  ERROR: Max NRAIN= 200, used ',I5,
     &     ' steps. Please correct and rerun  *',/,71('*'),/)
        STOP
      ENDIF
      TOTRAIN=0.d0
      STARTRAIN=0.d0
      ENDRAIN=0.d0
      DO 20 I=1,NRAIN
        READ(2,*,iostat=ios)(RAIN(I,J),J=1,2)
        IF((RAIN(I,2).gt.0.d0).and.(STARTRAIN.eq.0.d0)) STARTRAIN=
     &      RAIN(I,1)+1.0d-11
        IF (I.GT.1) THEN
            TOTRAIN=TOTRAIN+RAIN(I-1,2)*(RAIN(I,1)-RAIN(I-1,1))
            IF((RAIN(I,2).le.0.d0).and.(STARTRAIN.gt.0.d0).and.
     &        (ENDRAIN.eq.0.d0)) THEN
              ENDRAIN=RAIN(I,1)
             ELSEIF((RAIN(I,2).le.0.d0).and.(STARTRAIN.gt.0.d0).and.
     &        (RAIN(I-1,2).gt.0.d0)) THEN
              ENDRAIN=RAIN(I,1)
            ENDIF
        ENDIF
20    CONTINUE
      IF(IOS.NE.0) THEN
        WRITE(*,111)
111     format(/,71('*'),/,'*  NRAIN does not match rain steps.'
     &     ,' Please correct and rerun          *',/,71('*'),/)
        STOP
      ENDIF
      DR1=RAIN(NRAIN,1)
c----Average rain for d50 (NPART8) calculation
      IF(ENDRAIN.eq.0.d0) THEN
            RMMH=0.d0
        ELSE
            RMMH=TOTRAIN*1000.d0/(ENDRAIN-STARTRAIN)*3600.d0
      ENDIF

C-----(.ISO file) Read Green-Ampt infiltration parameters -------------------
      IF(ISCR.EQ.0) THEN
        WRITE(*,'(" ... Reading inputs from: ",A45)')LISFIL(4)
      ENDIF
      WTD=999.d0
      READ(7,*)VKS,Sav,OS,OI,SM,SCHK

C-----Get downslope node for flood checking (0=<schck=<1) or average along --
C-----filter Xavg (nchk= -1)
      IF(SCHK.EQ.-1.d0) THEN
            NCHK=-1
       ELSEIF(SCHK.LT.0.OR.SCHK.GT.1) THEN
         WRITE(*,109)
109      format(/,71('*'),/,'*  ERROR: SCHK is outside [0,1] range.
     &   Please correct and rerun',5(' '),'*',/,71('*'),/)
         STOP 
       ELSEIF(SCHK.EQ.0) THEN
         NCHK=1
       ELSE  
         NCHK=IDNINT(SCHK*N)
      ENDIF
c---------- check if there is shallow water table
      READ(7,*,END=25)CWTD
      BACKSPACE(7)
      READ(7,*,ERR=25)WTD

c---------- mod for shallow water table, rmc, 10/2013
!      WRITE(*,'(4x," Water table depth at (m):",F6.2)')WTD
!      WRITE(*,*) '... Reading shallow water table parameters'
! --------Soil water characteristic curve (ITHETATYPE= 1:vG; 2:BC)
      READ(7,*) ITHETATYPE
      BACKSPACE(7)
      IF(ITHETATYPE.eq.1) then
            READ(7,*) ITHETATYPE,PARW(1),PARW(2),PARW(3),PARW(4)
        elseif(ITHETATYPE.eq.2) then
            READ(7,*) ITHETATYPE,PARW(1),PARW(2),PARW(3)
        else
            READ(7,*) ITHETATYPE,PARW(1),PARW(2),PARW(3),PARW(4)
      ENDIF
! --------Hydraulic conductivity curve (IKUNSTYPE= 1:vG; 2:BC; 3:Gardner)
      READ(7,*) IKUNSTYPE
      BACKSPACE(7)
      IF(IKUNSTYPE.eq.2) then
            READ(7,*) IKUNSTYPE,PARK(1),PARK(2)
        else
            READ(7,*) IKUNSTYPE,PARK(1)
      ENDIF
! RMC----Set value of hb. By default hb=0 unless BC (ITHETATYPE=2)
!------- then hb=1/PAR2
      hb=0.d0
      IF(ITHETATYPE.EQ.2) hb=1/PARW(2)
      zw=dabs(WTD)-dabs(hb)
! CL-----Input checking for shallow water table
      IF(PARW(3).LT.1.d0.AND.ITHETATYPE.EQ.1)write(*,*)
     &     'WARNING:vnmvg<1 !!!'
      IF(PARW(2).EQ.0.d0) write(*,*)'warning: alpha = 0 !!!'
      IF(IKUNSTYPE.EQ.2.AND.PARK(2).EQ.0.d0) write(*,*)
     &    'WARNING: alpha BC = 0 !!!'
      IF (.NOT.ITHETATYPE.EQ.1.AND..NOT.ITHETATYPE.EQ.2) WRITE(*,*)
     &      'WARNING: missing code for soil water characteristic
     &      curve type !!'
      IF (.NOT. IKUNSTYPE.EQ.1 .AND. .NOT. IKUNSTYPE.EQ.2 .AND.
     &      .NOT.IKUNSTYPE.EQ.3) WRITE(*,*)'WARNING: missing code for
     & unsaturated hydraulic conductivity curve type!!'
      IF (WTD.LT.0.d0) WRITE(*,*)'WARNING: WTD must be positive!!
     & Please change in ISO file- Taken as positive'

c-RMC05/2020- Changed form previous versions. Three bottom boundary conditions
c(ITWBC) are included after wetting front reaches WT (t>=tw and z=zw>0),
c-----------Type BC1 (default): Dupuis-Forchheimer, f=Ksh.WTD/VL.So
c-----------Type BC2: Vertical saturated flow (Salvucci & Entekabi),f=Ks
c-----------Type BC3: Simplified, f=Ksh.So
c----NOTE: However, BC1 is hard-set as it is the one relevant to practical VFS
c----------applications (lateral draining to stream). The other 2 are kept for
c----------research and can be set in the code but not the user in the .ISO file
C      IF(ITWBC1.EQ.1) THEN
C         BACKSPACE(7)
C         READ(7,*)ITWBC,RVH
C       ELSEIF(ITWBC2.EQ.1) THEN
C         BACKSPACE(7)
C         READ(7,*)ITWBC,RVH
C       ELSEIF(ITWBC3.EQ.1) THEN
C         BACKSPACE(7)
C         READ(7,*)ITWBC,RVH
C      ENDIF
      ITWBC=1

c-RMC05/2020- Read if present RVH= Ksv/Ksh in ISO (last line). If not RVH=1.0
      RVH=1.d0
      READ(7,*,END=25)CWTD
      BACKSPACE(7)
      READ(7,*,ERR=25)RVH
      IF (RVH.LE.0.d0) RVH=1.d0

c----- mod for shallow water table, rmc, 10/2013
25    DM=OS-OI
      IF(DM.LE.0.d0)DM=0.0001d0
      IF(VKS.LT.0.d0) VKS=0.d0
      SavM=Sav*DM
      AGA= VKS
      BGA= VKS*SAVM

C-----(.IRO file) Read runoff inflow at upper side of strip (BC) in (m3/s) ---
      IF(ISCR.EQ.0) THEN
        WRITE(*,'(" ... Reading inputs from: ",A45)')LISFIL(3)
      ENDIF
      READ(3,*)SWIDTH,SLENGTH
      SAREA=SWIDTH*SLENGTH
      READ(3,*) NBCROFF,BCROPEAK
      IF(NBCROFF.GT.200) THEN
        WRITE(*,115)NBCROFF
115     format(/,71('*'),/,'*  ERROR: Max NBCROFF= 200, used ',I5,
     &     ' steps. Please correct and rerun      *',/,71('*'),/)
        STOP
      ENDIF
      QSUM0=0.D0
      TOTBCROFF=0.d0
      DO 30 I=1,NBCROFF
        READ(3,*,iostat=ios)(BCROFF(I,J),J=1,2)
        IF(IOS.NE.0) THEN
          WRITE(*,112)
112       format(/,71('*'),/,'*  NBCROFF does not match runoff steps'
     &     ,'. Please correct and rerun      *',/,71('*'),/)
          STOP
        ENDIF
        IF(I.GE.2) THEN
          TIMEINCR=BCROFF(I,1)-BCROFF(I-1,1)
          AREA0=TIMEINCR*(BCROFF(I,2)+BCROFF(I-1,2))/2.D0
          QSUM0=QSUM0 + AREA0
        ENDIF
30    CONTINUE
      DR2=BCROFF(NBCROFF,1)
      DR=DMAX1(DR1,DR2)
      IF(DR.GT.DR1) THEN
        NRAIN=NRAIN+1
        RAIN(NRAIN,1)=DR
        RAIN(NRAIN,2)=0.D0
      ENDIF
c-----Handle case with no water (rain and runoff) is provided
      IF(QSUM0.EQ.0.d0.and.TOTRAIN.EQ.0.d0) THEN
        WRITE(*,116)
116     format(/,64('*'),/,'*  NOTHING TO DO: NO RAIN OR',
     &     ' FIELD RUNOFF PROVIDED.            *',/,64('*'),/)
        STOP
      ENDIF

C-------Find the bandwidth for the matrix, #element, #nodes---
      NBAND=2*NPOL-1
      DX=VL/(N-1)
      NELEM=(N-1)/(NPOL-1)

C-------Calculate convergence and wave form parameters--------
C***** English Units
c        G=32.185D0
c       CMN=1.486D0
C***** Metric Units, PEAK & RPEAK (m/s), BCROPEAK (m3/s), QMAX(m2/s)----
      G=9.81D0
      CMN=1.D0
      VMN=5.D0/3.D0
      PEAK=RPEAK+BCROPEAK/(VL*FWIDTH)
      QMAX= VL*PEAK
      HMAX= (QMAX/BIGQK)**(1.D0/VMN)
      VMAX=QMAX/HMAX
      FR=VMAX/(G*HMAX)**0.5D0
      FK=(VL*SO*G)/VMAX**2.D0
      TE= HMAX/PEAK
c-- Old time step scheme based on Courant number
      C= VMN*BIGQK*HMAX**(VMN-1.D0)
      DTC= CR*DX/C
      
c-09/15 New scheme for dynamic time step (Jaber & Mohtar, 2002, Table 1,
c-----https://doi.org/10.1016/S0309-1708(02)00005-2.
c---- Time of concentration (h), approximation Ragan and Duru (1972) (RD)
      TIMEC1=0.0803D0*(VN1*VL)**0.6d0/(PEAK*3.6D5*SC**0.3d0)
c---- Time of concentration (h), method of characteristics (MOC)
      TIMEC2=PEAK**((1.d0-VMN)/VMN)*(VL/(DSQRT(SC)/VN1))**(1.d0/VMN)/
     &  3600.d0
c---- Dynamic time step scheme for varying node properies (Jaber & Mohtar, 2002)
c-----Use tc MOC (h)
      TIMEC=TIMEC2
c-----a1) Consistent scheme constant excess rainfall  (Table 1)
      DT=0.96d0*TIMEC/(NELEM**1.27d0)*3600.d0*CR
c-----a2) Consistent scheme variable excess rainfall  (Table 1)
c      DT=1.66d0*TIMEC/(NELEM**1.43d0)*3600.d0*CR
c-----b1) Lumped scheme constant excess rainfall (Table 2)
c      DT=1.27d0*TIMEC/(N**1.16d0)*3600.d0*CR
c-----b2) Lumped scheme variable excess rainfall (Table 2)
c      DT=0.79d0*(TIMEC/N)*3600.d0*CR
c-----choose conservative time step
      IF(DT.GT.DTC) DT=DTC
      CRR=C*DT/DX
c-- Find number of time steps in problem
      NDT=IDINT(DR/DT)
C-------Calculate the PG Parameters (n=50) (Muñoz-Carpena et al., 1993)-
c-09/15 Use new CRR calculated from dynamic time step
      IF(KPG.EQ.1)THEN
        PGPAR(1)=0.0215873D0 - 0.345217D0*CRR + 1.33259D0*CRR**2.D0 -
     &               1.62016D0*CRR**3.D0 + 0.670333D0*CRR**4.D0
        PGPAR(2)= 0.0592655D0 - 0.107237D0*CRR + 0.235216D0*CRR**2.D0 -
     &               0.426017D0*CRR**3.D0 + 0.222228D0*CRR**4.D0
        PGPAR(3)=0.0280422D0 + 0.175632D0*CRR - 0.592941D0*CRR**2.D0 -
     &               0.149698D0*CRR**3.D0 - 0.0704731D0*CRR**4.D0
        PGPAR(4)= -0.0456247D0 +0.00112745D0*CRR +0.420433D0*CRR**2.D0 -
     &               0.0935913D0*CRR**3.D0 - 0.0764558D0*CRR**4.D0
      ENDIF

C-------Set the order of the integration rule-------------------

      IF(KPG.EQ.0.OR.(PGPAR(4).EQ.0.D0.AND.PGPAR(3).EQ.0.D0))THEN
         NL=NPOL+1
        ELSE
         NL=5
      ENDIF

C-----(.IWQ file) [06/2008]-Read water quality parameters, ensure backwards compatibility ---

      IF(IWQ.GT.0) THEN
         IF (ISCR.EQ.0) THEN
           WRITE(*,'(" ... Reading inputs from: ",A45)')LISFIL(12)
         ENDIF
         OPEN(17,FILE=LISFIL(12),ERR=1500,STATUS='OLD')
         OPEN(18,FILE=LISFIL(13),STATUS='UNKNOWN')
         WRITE(18,220)LISFIL(13)
         write(15,225) 12,'iwq',lisfil(12)
         write(15,225) 13,'owq',lisfil(13)

c---------Read water quality problems (IWQPRO)
C----------IWQPRO=1 - Sabbagh et al. (2009) semiempirical eq.
C----------IWQPRO=2 - Refitted Sabbagh et al. eq. w/user supplies coeff.
C----------IWQPRO=3 - Based on Muñoz-Carpena et al. (2015) mechanistic eq.
C----------IWQPRO=4 - Chen et al. (2017) empirical eq.
         ICAT=0
         CSAB(1)= 24.79D0
         CSAB(2)= 0.54D0
         CSAB(3)= 0.52D0
         CSAB(4)= -2.42D0
         CSAB(5)= -.89D0
         READ(17,*)IWQPRO
         IF(IWQPRO.EQ.2) THEN
             BACKSPACE(17)
             READ(17,*,ERR=50)IWQPRO,(CSAB(I),I=1,5)
           ELSEIF(IWQPRO.LT.1.or.IWQPRO.GT.4) THEN
             IWQPRO=3
         ENDIF
50       IF(IWQPRO.GT.0.AND.IWQPRO.LE.4) THEN
           READ(17,*)IKD
           BACKSPACE(17)
           IF(IKD.EQ.1.and.IWQ.eq.1) THEN
              READ(17,*,iostat=ierr)IKD,VKOC(1),OCP
              if (ierr.ne.0) then
                  write(*,199)
                  print*,'ERROR: line 2 of .iwq file. With IKD=1, both Koc',
     &               '       and %OC must be provided. Please fix and rerun.'
                  write(*,199)
                 stop
               endif
            elseif(IKD.EQ.0.and.IWQ.eq.1) then
               READ(17,*,iostat=ierr)IKD,VKD(1)
              if (ierr.ne.0) then
                  write(*,199)
                  print*,'ERROR: line 2 of .iwq file. Missing IKD and/or Kd.',
     &               '       Please fix and rerun.'
                  write(*,199)
                  stop
               endif
            elseif(IKD.EQ.1.and.IWQ.GT.1) then
               READ(17,*,iostat=ierr)IKD,VKOC(1),OCP,(VKOC(JJ),JJ=2,IWQ)
               if (ierr.ne.0) then
                  write(*,199)
                  print*,'ERROR: line 2 of .iwq file. Missing Koc for IWQ',
     &               '     compounds set in .ikw. Please fix and rerun.'
                  write(*,199)
                  stop
               endif
               if (OCP.GT.100.d0) then
                  write(*,199)
                  print*,'ERROR: line 2 of .iwq file. %OC>100% not possible. ',
     &               '       Please fix and rerun.'
                  write(*,199)
                  stop
               endif
            else
               READ(17,*,iostat=ierr)IKD,VKD(1),(VKD(JJ),JJ=2,IWQ)
               if (ierr.ne.0) then
                  write(*,199)
                  print*,'ERROR: line 2 of .iwq file. Missing Kd for IWQ',
     &               '     compounds set in .ikw. Please fix and rerun.'
                  write(*,199)
                  stop
               endif
         ENDIF            
c---------Check when Chen et al. eq. was selected that KOC is provided
           IF(IWQPRO.EQ.4.AND.IKD.EQ.0) THEN
                  WRITE(*,120)
120    format(/,71('*'),/,'*  ERROR: Chen et al.(2016) eq. was',
     &     ' selected and KOC & %OC must be    *',/,'*  provided ',
     &     'not KD. Please correct & rerun.',22(' '),'*',/,
     &     71('*'),/)
                  STOP
           ENDIF
           IF(IKD.EQ.1) THEN
             DO 55 JJ=1,IWQ
                VKD(JJ)=VKOC(JJ)*OCP*.01D0
                IF(VKOC(JJ).GT.9000.D0) ICAT=1
55           CONTINUE
           ENDIF
           READ(17,*)CCP
c---------- check if pesticide degradation is requested (IDG=1 TO 4)
           IDG=0
           READ(17,*,END=32)CIDG
           BACKSPACE(17)
           READ(17,*,ERR=32)IDG
           IF(IDG.LT.0.OR.IDG.GT.4) THEN
             IDG=0
             write(*,199)
             write(*,*)'WARNING: IDG outside 0 to 4 range; set to 0',
     &                 ' (no degradation)'
             write(*,199)
           ENDIF
c----------When Pesticide degradation requested: IDG=1-4'         
           IF (IDG.GT.0) THEN
            READ(17,'(A)',iostat=io) line
            ncount = 0
            ipos = 1
            maxlen = len(trim(line))
c-------   Parse numbers from the line
            do while (ipos.LE.maxlen)
                  call getnumber(line, ipos, vnum, ios)
                  if (ios.ne.0) then
                    exit
                  else
                    ncount = ncount + 1
                    if (ncount.GT.100) then
                        print *, 'ERROR: .iwq file line 5. Exceeded max.',
     &                  'number compounds.'
                        stop
                    endif
                    test_array(ncount) = vnum
                  endif
            end do
            nexpect= 7+4*(IWQ-1)
c-------   debug error checking for compound properties ----
c           print *,'Numbers read and expected: ', ncount, nexpect
c           do i = 1, ncount
c                 print *, test_array(i)
c           end do
c-----------------------  
            BACKSPACE(17)
            IF(nexpect.EQ.ncount) then
                  READ(17,*,iostat=ierr)NDGDAY,DGHALF(1),FC,DGPIN(1),DGML,DGLD(1),
     &             DGMRES0(1),((dgHalf(JJ),dgPIN(JJ),dgLD(JJ),dgmres0(JJ)),JJ=2,IWQ)
                ELSEIF(IWQ.GT.1) THEN
                  write(*,199)
                  write(*,121)
                  write(*,199)
                  STOP
121   format('ERROR: .iwq file line 5. Number of inputs does not match number of',/
     &            ' chemicals (IWQ) selected in .ikw file. Please fix and rerun.')
               ELSEIF(ncount.EQ.6) THEN
                  READ(17,*,iostat=ierr)NDGDAY,DGHALF(1),FC,DGPIN(1),DGML,DGLD(1)
                  DGMRES0(1)=0.d0
                  write(*,199)
                  write(*,122)
                  write(*,199)
122   format('WARNING: .iwq file line 5. Missing DGMRES0 in IWQ Line 5. Set to 0 mg/m2')
               ELSEIF(ncount.EQ.5) THEN
                  READ(17,*,iostat=ierr)NDGDAY,DGHALF(1),FC,DGPIN(1),DGML
                  DGMRES0(1)=0.d0
                  DGLD(1)=0.05d0
                  write(*,199)
                  write(*,123)
                  write(*,199)
123   format('WARNING: .iwq file line 5. Missing DGLD and DGMRES0.',
     &                ' Set to 0.05 m, 0 mg/m2.')                        
               ELSE
                  write(*,199)
                  write(*,121)
                  write(*,199)
                  STOP
            ENDIF
c-rmc10/2021--Add surface redidues from last event to incoming mass---
            DO 127 JJ=1,IWQ
              DGPIN(JJ)=DGPIN(JJ)+DGMRES0(JJ)
              DGKREF(JJ)=DLOG(2.D0)/DGHALF(JJ)
127         continue
            IF(NDGDAY.LT.0.OR.NDGDAY.GT.365) THEN
              WRITE(*,128)NDGDAY
128   format(/,71('*'),/,'*  ERROR: NDGDAY= 0-365, used ',I5,
     &        ' days. Please correct and rerun    *',/,71('*'),/)
              STOP
            ENDIF
            READ(17,*,iostat=ierr)(DGT(I),I=1,NDGDAY)
            READ(17,*,iostat=ierr)(DGTHETA(I),I=1,NDGDAY)
            if (ierr.ne.0) then
              WRITE(*,129)
129    format(/,71('*'),/,'*  ERROR: Number of T and theta values in '
     &      ,'.IWQ is less than NDGDAY.   *',/,'*  Please correct and'
     &      ,' rerun',43(' '),'*',/,71('*'),/)
              STOP
            ENDIF
c-rmc01/2020--Check if remobilization scheme requested (IMOB=1 TO 3)
            READ(17,*,END=33)CIMOB
            BACKSPACE(17)
            READ(17,*,ERR=33)IMOB
c------------Residue remobilitation requested: IMOB
33          IF(IMOB.LT.1.OR.IMOB.GT.3) IMOB=1
c ----------Placeholder for other wq problems, i.e. TaRSE, 07/28/08 rmc
c
c-----------Multispecies (IWQ>2), read DGMOL(JJ), DGFRAC(JJ,JJ)
            IF(IWQ.GT.1) THEN
              READ(17,*,iostat=ierr)(DGMOL(JJ),JJ=1,IWQ)
               if (ierr.ne.0) then
                write(*,199)
                print*,'ERROR: line 9 of .iwq file. Missing DGMOL for',
     &         '      chemicals (required for IWQ>1). Please fix and rerun.'
                write(*,199)
                stop
              endif
              DO 130 II=1,IWQ
                  READ(17,*,iostat=ierr)(DGFRAC(II,JJ),JJ=1,IWQ)
                  if (ierr.ne.0) then
                        write(*,199)
                        print*,'ERROR: line 10 of .iwq file. Missing DGFRAC() for',
     &         '      chemicals (required for IWQ>1). Please fix and rerun.'
                        write(*,199)
                        stop
                  endif
130           CONTINUE                           
            ENDIF
          ELSE
             WRITE(*,131)
131    format(/,71('*'),/,'*     ERROR: WQ problem selected is not',
     &    ' available in this version     *',/,'*     (IWQPRO= 1-4 ',
     &    'only). IWQ set 0 (no water quality).',21(' '),'*',/,71('*'),/)
             IWQ=0
         ENDIF
        ELSE
           write(15,*)
           write(15,*)
      ENDIF
      ENDIF

C-------Output all the parameters-------------------------------

32    WRITE(11,*)'Storm parameters'
      WRITE(11,*)'----------------'
      WRITE(11,140)PLABEL
      WRITE(11,180)
      WRITE(11,150)
      WRITE(11,160)
      WRITE(11,180)
      DO 35 I=1,NRAIN-1
         WRITE(11,170)I,RAIN(I,1),RAIN(I+1,1),RAIN(I,2)
35    CONTINUE
      WRITE(11,180)
      WRITE(11,350)'Total rainfall (mm)=',TOTRAIN*1000.d0
      WRITE(11,400)'Peak rainfall intensity(m/s)=',RPEAK
      WRITE(11,400)'Peak inflow rate BC (m3/s)=',BCROPEAK
      WRITE(11,*)'        (The inflow hydrograph can be found in the',
     &              ' OUTPUTS)'
      WRITE(11,*)
      WRITE(11,*)'Filter parameters'
      WRITE(11,*)'-----------------'
              WRITE(11,200)'Length of the strip (m)=',VL
              WRITE(11,200)'Width of the strip (m)=',FWIDTH
              WRITE(11,200)'Surface characteristics='
      WRITE(11,525)
      WRITE(11,575)
      WRITE(11,525)
      DO 40 IPROP=1,NPROP
            IF(IPROP.EQ.1) THEN
                  EX1=0.d0
              ELSE
                  EX1=(NODEP(IPROP-1)-1)*DX
            ENDIF
            EX2=(NODEP(IPROP)-1)*DX
            WRITE(11,550)EX1,EX2,RNA(IPROP),SOA(IPROP)
40    CONTINUE
      WRITE(11,525)
      WRITE(11,*)
C-------Output nodal information if desired (ielout=1)----------------
      IF(IELOUT.EQ.1)THEN
        WRITE(11,*)' Elemental information follows (IELOUT=1):'
        WRITE(11,185)
        WRITE(11,186)
        WRITE(11,185)
        DO 45 NEL=1,NELEM
            K=(NPOL-1)*NEL-NPOL+1
            DO 44 I=1,NPOL
                K=K+1
                WRITE(11,600)NEL,K,I,QK(K),(k-1)*dx
44          CONTINUE
45        CONTINUE
        WRITE(11,185)
        WRITE(11,*)
      ENDIF

      WRITE(11,*)'Soil parameters'
      WRITE(11,*)'---------------------'
      WRITE(11,400)'Saturated hydraulic cond.(Ks)=',VKS
      IF(0.D0.LE.WTD.AND.WTD.LT.10.D0) THEN
         WRITE(11,201)'Water table depth at (m)=',WTD,'BC(t>tw),Rvh='
     &         ,ITWBC,RVH
         IF(ITHETATYPE.EQ.1) THEN
             WRITE(11,205)'van Genuchten'
             WRITE(11,200)'Or,alpha,n,m =',(PARW(k), k=1,4)
           ELSEIF(ITHETATYPE.EQ.2) THEN
             WRITE(11,205)'Brooks & Corey'
             WRITE(11,200)'Or,alpha,lambda =',(PARW(k), k=1,3)
         ENDIF
         IF(IKUNSTYPE.EQ.1) THEN
             WRITE(11,208)'van Genuchten'
             WRITE(11,200)'m=',PARK(1)
           ELSEIF(IKUNSTYPE.EQ.2) THEN
             WRITE(11,208)'Brooks & Corey'
             WRITE(11,200)'eta,alpha=',PARK(1),PARK(2)
           ELSEIF (IKUNSTYPE.EQ.3) THEN
             WRITE(11,208)'Gardner'
             WRITE(11,200)'alpha=',PARK(1)
         ENDIF
      ENDIF
      WRITE(11,200)'Sat. soil-water content(Os)=',Os
      IF(WTD.eq.999.d0)THEN
        WRITE(11,200)'Initial soil-water content(Oi)=',Oi
        WRITE(11,200)'Initial soil-water deficit (M)=',DM
        WRITE(11,200)'Avg. suction at wet front(Sav)=',Sav
        WRITE(11,400)'Green-Ampt parameters (A,B)=',AGA,BGA
      ENDIF
      IF(NCHK.LT.0) THEN
        WRITE(11,210)'Node number for flood checking=','Average'
       ELSE
        WRITE(11,209)'Node number for flood checking= ',NCHK
      ENDIF
      WRITE(11,*)
      WRITE(11,*)'Simulation parameters'
      WRITE(11,*)'---------------------'
      WRITE(11,400)'Length of simulation (s)=',DR
      WRITE(11,209)'Order of the basis functions=',NPOL-1
C      WRITE(11,*)' Output option (0=q@t;1= h@x)= ',IOUT
      WRITE(11,700)'Petrov-Galerkin parameters=',(PGPAR(I),I= 1,4)
      WRITE(11,200)'Time weighting parameter=',THETAW
      WRITE(11,200)'Space step, dx(m) =',DX
      WRITE(11,200)'Time step, dt (s) =',DT
      WRITE(11,*)'Number of nodes in system    =',N
      WRITE(11,*)'Number of elements in system =',NELEM
      WRITE(11,*)'Number of time steps        =',NDT
      WRITE(11,*)'Maximum number of iterations =',MAXITER
      WRITE(11,200)'Maximum flow rate and depth=',QMAX,HMAX
      WRITE(11,200)'Celerity of the wave=',C
      WRITE(11,200)'Courant time step=',DTC
      WRITE(11,200)'Froude number=',FR
      WRITE(11,350)'Kinematic wave number=',FK
      WRITE(11,200)'Courant number=',CRR
      IF(ICO.EQ.0)THEN
        WRITE(11,210)'Surface changes feedback=',' NO'
       ELSE
        WRITE(11,210)'Surface changes feedback=',' YES'
      ENDIF
      WRITE(11,*)

C-------Issue a warning if any of the criteria is not met-------

      IF(FK.LT.10.D0) THEN
        WRITE(*,*)'WARNING: Kinematic number smaller than 10'
       ELSE IF (FR.GT.1.5D0) THEN
        WRITE(*,*)'WARNING: Froude number greater than 2'
       ELSE IF (CRR.GT.1.D0) THEN
        WRITE(*,*)'WARNING: Courant number greater than 1'
      ENDIF
      WRITE(11,*)

C-------Print header for output values--------------------------

      IF(IOUT.EQ.0.D0) THEN
        WRITE(11,190)
        WRITE(11,192)
        WRITE(11,198)
        WRITE(11,500)0.d0,0.d0,0.d0,0.d0,BCROFF(1,2),0.d0,0.d0,
     &         0.d0,0
       ELSE
        WRITE(11,194)
        WRITE(11,196)
        WRITE(11,198)
      ENDIF

c-----------Output all input values for Water Quality (if IWQ>0)--------
      IF(IWQ.GT.0) THEN
            WRITE(18,803)'Parameters for Water Quality - No. Compounds=',IWQ
            SELECT CASE (IWQPRO)
              CASE (1)
                 WRITE(18,*)'Type of problem - ',
     &           'Pesticide trapping (Sabbagh et al., 2009)'
              CASE (2)
                 WRITE(18,799)'Type of problem - ',
     &           'Pesticide trapping (Sabbagh refit):',
     &           (CSAB(I),I=1,5)
              CASE (3)
                 WRITE(18,*)'Type of problem - ',
     &           'Pesticide trapping (mass balance (Reichenberger 2019)'
              CASE DEFAULT
                 WRITE(18,*)'Type of problem - ',
     &           'Pesticide trapping (Chen et al.,2017)'
            END SELECT
            WRITE(18,800)'Partition coefficient (Kd)=',VKD(1),'L/Kg (Compound 1) '
            WRITE(18,800)'% Clay in sediment (%CL)  =',CCP,'%                '
            WRITE(18,800)'Dispersion length (l)=',DGLD(1),'m (Compound 1)    '
            SELECT CASE (IMOB)
              CASE (2)
                WRITE(18,607)
              CASE (3)
                WRITE(18,608)
              CASE DEFAULT
                WRITE(18,606)
            END SELECT
            IF(IDG.GT.0.AND.IDG.LE.4) THEN
                SELECT CASE (IDG)
                  CASE (1)
                    WRITE(18,601)
                  CASE (2)
                    WRITE(18,602)
                  CASE (3)
                    WRITE(18,603)
                  CASE DEFAULT
                    WRITE(18,604)
                END SELECT
              ELSE
                WRITE(18,605)
            ENDIF
            WRITE(18,800)'Pesticide half-life (Ln2/Kref)=',
     &               DGHALF(1),'days (Compound 1) '
            WRITE(18,800)'  Soil field capacity (FC)=',FC,
     &               '(-)              '
            WRITE(18,800)'Incoming pesticide mass (mi)=',
     &               DGPIN(1)-DGMRES0(1),'mg/m2 (Compound 1)'
            WRITE(18,800)'Mixing layer thickness (dml)=',
     &               DGML,'cm               '
            WRITE(18,800)'Remob pest prev. event (mres0)=',
     &               DGMRES0(1),'mg/m2 (Compound 1)'
            WRITE(18,801)'No. of days between events=',NDGDAY
            WRITE(18,802)
            WRITE(18,910)(I,DGT(I),DGTHETA(I),I=1,NDGDAY)
            IF(IWQ.GT.1) THEN
                WRITE(18,200)'Pest. molar mass (dgMOL,g/mol)=',
     &               (DGMOL(JJ),JJ=1,IWQ)
                WRITE(18,200)'Pesticide half-life (Ln2/Kref)=',
     &               (DGHALF(JJ),JJ=1,IWQ)
                WRITE(18,*)'Molar transformation fractions of Compounds (dgFRACj)(-)='
                WRITE(18,'(31x,A7,10I8)')'Compound',(JJ,JJ=1,IWQ)
                DO 132 II=1,IWQ
                  WRITE(18,805)II,(DGFRAC(II,JJ),JJ=1,IWQ)
 132            CONTINUE
            ENDIF                           
      ENDIF

140   FORMAT(20x,'Storm data: ',A30)
145   FORMAT(1x,'Storm on: ',a30,14x,'...RUNNING...')
150   FORMAT(20x,'|Period|',4x,'Time interval',4x,'|  Rainfall  |')
160   FORMAT(20x,'|',6x,'|',9x,'(s)',9x,'|',3x,'(m/s)',4x,'|')
170   FORMAT(20x,'|',I5,' |',F8.1,'  to',F8.1,' |',E11.4,' |')
180   FORMAT(20x,'+',6('-'),'+',21('-'),'+',12('-'),'+')
185   FORMAT(20x,3('+',5('-')),'+',9('-'),'+',9('-'),'+')
186   FORMAT(20x,'| Elem| node|local| alpha   |   x(m)  |')
190   FORMAT(5x,'TIME',5x,'OUTFLOW',4x,'CUM.FLOW',5x,'ie =r-f',
     &     5x,'INFLOW',4x,'CUM.INFLOW',7x,'f',10x,'z',8x,'ITER')
192   FORMAT(5x,'(s)',6x,'(m3/s)',7x,'(m3)',8x,'(m/s)',6x,'(m3/s)',
     &     7x,'(m3)',8x,'(m/s)',7x,'(m)')
194   FORMAT('ITER',6x,'TIME',7x,'INFLOW',6x,'ie =r-f',
     &    5x,'DEPTH (X=L/2)')
196   FORMAT(12x,'(s)',10x,'(m)',9x,'(m/s)',9x,'(m)')
198   FORMAT(101('-'))
199   FORMAT(70('-'))
200   FORMAT(A31,10F12.6)
201   FORMAT(A31,F12.6,4x,A13,I4,F10.4)
205   FORMAT(4x,'SWCC curve in infiltration=',3x,A15)
208   FORMAT(2x,'KUNSAT curve in infiltration=',3x,A15)
209   FORMAT(A31,I12)
210   FORMAT(A31,A12)
220   FORMAT('File: ',A40,8x,'VFSMOD v4.6.0 11/2024')
225   format(3x,'File #=',i3,' code:',a3,'=',a)
350   FORMAT(A31,2F12.2)
400   FORMAT(A31,2E12.4)
500   FORMAT(E11.4,7E12.4,I6)
525   FORMAT(20x,'+',21('-'),'+',9('-'),'+',9('-'),'+')
550   FORMAT(20x,'|',F8.4,'  to',F8.4,' |',2(F8.4,' |'),E11.4,' |')
575   FORMAT(20x,'|',4x,'x(m) interval',4x,'|',4x,'n',4x,'|',
     &          4x,'So',3x,'|')
600   FORMAT(20x,'|',3(I5,'|'),2(F9.4,'|'))
601   FORMAT(8x,'Degradation type (IDG)=   1: EU-FOCUS k=Kref.k(T).',
     & 'k(theta)')
602   FORMAT(8x,'Degradation type (IDG)=   2: US-EPA k=Kref')
603   FORMAT(8x,'Degradation type (IDG)=   3: k=Kref.k(T)')
604   FORMAT(8x,'Degradation type (IDG)=   4: k=Kref.k(theta)')
605   FORMAT(8x,'Degradation type (IDG)=   0: no degradation')
606   FORMAT(1x,'Residue remobilitation (IMOB)=   1: partial (surface',
     & ' porewater)')
607   FORMAT(1x,'Residue remobilitation (IMOB)=   2: full (surface ',
     & 'porewater+sorbed)')
608   FORMAT(1x,'Residue remobilitation (IMOB)=   3: no remobilization')
700   FORMAT(A31,10F9.5)
799   FORMAT(A18,A34,5F7.3,')')
800   FORMAT(A31,F12.6,A20)
801   FORMAT(A31,I5)
802   FORMAT(40x,'day    T(C)  theta(-)')
803   FORMAT(51('-'),/,A45,I4,/,51('-'))
805   FORMAT(31x,I8,10F8.2)
910   FORMAT(365(35x,I6,2F10.5,/))

      RETURN

1500  WRITE(*,1600)'ERROR: Input file .iwq missing (check .ikw and project)'
1600  FORMAT(/,A55,/)
      STOP
      END

      subroutine getnumber(line, ipos, vnum, ios)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                 
C     Helper subroutine for error checking when reading. It parses the 
C     read line and identies the numbers to check if the number in the 
C     input file matches the number expected for the input
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	IMPLICIT DOUBLE PRECISION(a-h,o-z)
c     character*200 line
      CHARACTER*(*) line
      character*20 buffer

      ios = 0
      numstart = ipos
c     Skip non-numeric characters
      do while (numstart.LE.len(trim(line)))
         if (line(numstart:numstart) .GE. '0' .AND. line(numstart:numstart) .LE. '9') then
            if (numstart.EQ.1) then
                  exit   ! Skip first position on the line
              elseif(line(numstart-1:numstart-1) .EQ. ' ' .OR.
     &          line(numstart-1:numstart-1) .EQ. '(') then   ! Check if start of line or preceded by space or '('
               exit
            endif
         endif
         numstart = numstart + 1
      end do
      if (numstart > len(trim(line))) then
         ios = 1
         return
      endif
c     Find the end of the number
      numend = numstart
      numlen = len(trim(line))
      do while (numend <= numlen)
         if ((line(numend:numend) .GE. '0' .AND. line(numend:numend) .LE. '9') .OR.
     &       line(numend:numend) .EQ. '.' .OR. line(numend:numend) .EQ. 'E' .OR.
     &       line(numend:numend) .EQ. 'e' .OR.
     &       ((numend > numstart) .AND. (line(numend:numend) .EQ. '+' .OR. line(numend:numend) .EQ. '-'))) then
            numend = numend + 1
         elseif (line(numend:numend) .EQ. ' ' .OR. line(numend:numend) .EQ. ')') then  ! Check for termination condition
            exit
         else
            ios = 2
            return
         endif
      end do

      buffer = line(numstart:numend-1)
      read(buffer,*,iostat=ios) vnum
      if (ios /= 0) then
         ios = 3
         return
      endif
      ipos = numend

      return
      end