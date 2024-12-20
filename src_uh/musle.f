      subroutine musle(er,er1,erCoolm,ek,Y,pl,cfact,pfact,A,vol,tc,rain,
     1 D,soilty,dp,ieroty,sconc,om,aa1,b1,bigE,raimax30,qp)
c----------------------------------------------------
c    Calculation of soil loss and sediment concentration in runoff at 
c    the edge of the field (source) area using various modified USLE 
c    methods in metric units.
c    A0-A
c    A(Ha): Source area
c    vol(m3), volro (mm): runoff volume 
c    qp(m3/s),qpdepth(mm/h): runoff peak flow
c    ek: K erodibility factor in metric units [(kg.h)/(N.m^2)]
c    ekeng: K erodibility factor in English units [(0.1ton/ft-ton*)*(h/in)]
c----------------------------------------------------
      implicit double precision (a-h, o-z)
      character*20 soilty,musle_type
      character*20 types(21)
      dimension d50(21),sand(21),silt(21),Tf(21),Sf(21),Pf(21),sconc(6)

      data types/'Clay','Silty clay','Sandy clay','Silty clay loam',
     1  'Clay loam','Sandy clay loam','Silt','Silt loam','Loam',
     2  'Very fine sandy loam','Fine sandy loam','Sandy loam',
     3  'Coarse sandy loam','Loamy very fine sand','Loamy fine sand',
     4  'Loamy sand','Loamy coarse sand','Very fine sand',
     5  'Fine sand','Sand','Coarse sand'/
      data d50/23.d0,24.d0,66.d0,25.d0,
     1         18.d0,91.d0,19.d0,27.d0,35.d0,
     2         35.d0, 80.d0,98.d0,
     3         160.d0,90.d0,120.d0,
     4         135.d0,180.d0,140.d0,160.d0,170.d0,200.d0/
      data sand/20.d0,10.d0,50.d0,15.d0,35.d0,
     1          55.d0,5.d0,20.d0,45.d0,60.d0,
     2          60.d0,60.d0,60.d0,84.d0,84.d0,
     3          84.d0,84.d0,90.d0,90.d0,90.d0,
     4          90.d0/
      data silt/30.d0,45.d0,10.d0,50.d0,30.d0,
     1          20.d0,85.d0,60.d0,35.d0,25.d0,
     2          25.d0,25.d0,25.d0,8.d0,8.d0,
     3           8.d0,8.d0,5.d0,5.d0,5.d0,
     4           5.d0/
      data tf/0.01287d0,0.01870d0,0.01714d0,0.02606d0,0.0236d0,
     1        0.02778d0,0.05845d0,0.04259d0,0.03618d0,0.03877d0,
     2        0.03205d0,0.02549d0,0.01914d0,0.03726d0,0.02301d0,
     3        0.01624d0,0.00982d0,0.04401d0,0.02173d0,0.01481d0,
     4        0.00827d0/
      data sf/0.065d0,0.065d0,0.065d0,0.065d0,0.065d0,
     1        0.065d0,0.065d0,0.065d0,0.0325d0,-0.0325d0,
     2        0.d0,0.0325d0,0.0325d0,-0.0325d0,0.d0,
     3        0.0325d0,0.0325d0,-0.0325d0,0.d0,0.0325d0,
     4        0.0325d0/
      data pf/0.075d0,0.075d0,0.075d0,0.050d0,0.050d0,
     1        0.05d0,0.025d0,0.025d0,0.025d0,0.d0,
     2        0.d0,0.d0,0.d0,-0.025d0,-0.025d0,
     3        -0.025d0,-0.025d0,-0.05d0,-0.05d0,-0.05d0,
     4        -0.05d0/

c---------------------
c  compute R for musle (N/h)
c---------------------
c**convert bigE to SI metric - multiply by 1.702 / 100
c** units Rst=N/h, Runoff volume: A(Ha),vol(m3), volro (mm),qp(m3/s),qpdepth(mm/h)
      volro=vol/(A*10.d0)
      qpdepth=qp*360.d0/A
      Def=0.24d0*tc
      dtime=Def
      ndtime=D/dtime+1.d0
      bigEm=0.006700d0*bigE
      rst=1.702d0*(bigE/100.d0)*(raimax30/25.4d0/dtime)*dtime/0.5d0
c---------------------
c    er=> Foster et al. (1977b)
c---------------------
      rro=volro*(qpdepth)**(1.d0/3.d0)
      er=0.5d0* rst + 0.35d0*rro
c---------------------
c    er1=> Williams (1977)
c---------------------
      er1 = (9.d0*(vol*qp)**0.56d0)/A
      write(10,17)
17    format(/,/,1x,'RAINFALL ENERGY FACTOR R FOR EROSION CALCULATIONS')
c     1  ,'RAINFALL ENERGY FACTOR R FOR EROSION CALCULATIONS')
      write(10,22) bigE,bigEm,volro,qpdepth,rst,rro,er,A,vol,qp,er1
22    format(/,2x,'a) Foster et al. (1977)',
     1   /,4x,'E=  ',f10.3,' ft-tonf/acre =',f10.3,' MJ/ha',
     2   /,4x,'Volume Runoff=',f10.4,' mm; qpeak =',f10.4,' mm/h',
     3   /,4x,'Factors in Rm: Rstorm=',f10.4,'; Rrunoff=',f10.4,
     4   /,4x,'Rm (Foster)=',f10.4,' N/h',
     5   /,/,2x,'b) Williams (1975)',
     6   /,4x,'Watershed area=',f10.3,' ha',
     7   /,4x,'Volume runoff=',f10.4,' m3; qpeak=',f10.4,' m3/s',
     8   /,4x,'Rw (Williams)=',f10.4,' N/h')
c-------------------
c  erGLEAMS, from GLEAMS daily rain
c----------------------
      gei=7.87d0*(rain/25.4d0)**1.51d0
      geim=1.702d0*gei
      write(10,32) rain,gei,geim
32      format(/,2x,'c) GLEAMS/ daily CREAMS',
     1   /,4x,'Rain  =',f6.2,' mm',
     2   /,4x,'R_GLM =',f10.2,' From GLEAMS - Wischmeier',
     3   /,4x,'R_GLM =',f10.4,' N/h,  Converted to Metric')
c---------------------
c**   Cooley (1980) -> er for design storms, EI/100 = R ft tonsf/ac
c**                    1/0.67 * R for J/m^2
c---------------------
      erCooly=aa1*(rain/25.4d0)**(2.119d0*D**0.0086d0)/(D**b1)
      erCoolm=erCooly*1.702d0
      write(10,24) aa1,b1,rain,D,erCooly, erCoolm
24    format(/,2x,'d) Cooley (1980) - Design Storm ',
     1   /,4x,'a1=', f10.4, ' b1=',f10.4,
     2   /,4x,'Rain=', f10.3, ' mm   D=',f10.3,' hr',
     3   /,4x,'Rst =',f10.3,' ft-tonf/acre =',f10.3,' N/ha')

c--------------------
c  Soil type selection
c--------------------
      do i=1, 21
         if (soilty.eq.types(I)) then
            isoil=i
c            print *,'picked=>',isoil,soilty,types(i),d50(isoil)
         end if

      end do

c-------------
c K factor for MUSLE equations
c-------------
c -- Print a table of reference K values for different soil textures
      write(10,98)
98    format(/,/,1x,'ERODIBILITY K AND PARTICLE SIZE SELECTION',/,/,
     1 4x,'Table for computing Ksoil (from GLEAMS and KINEROS)-OM=2%',
     2 /,4x,' i',4x,'Soil Type',9x,'%Sand',2x,'%Silt',3x,'Tex.F.',
     3 4x,'Str.F.',4x,'Per.F.',4x,'D50(um)',1x,'K(kg-h/N-m^2)')
      om2=2.d0
      do i=1,21
c --  English units, ek0eng [(0.1ton/ft-ton*)*(h/in)] 
        ek0eng=tf(i)*(12.d0-om2)+sf(i)+pf(i)
c --  metric units, ek0 [kg/N*h/m^2] in 2 steps:
c --  a) (Eng->metric)  eK0=1.313* ek0eng  (USDA Ag. Handbook 537 Supp, 1981)
c --  b) (Tm/ha->Kg/m2) 0.1*ek0 
        ek0=1.313d0*0.1d0*ek0eng
        write(10,99)i,types(i),sand(i),silt(i),tf(i),sf(i),pf(i),
     1   d50(i),ek0
      end do
99    format(4x,i2,2x,a20,2x,f4.0,2x,f4.0,2x,f8.5,2x,f8.4,2x,f7.3,2x,
     1  f6.1,2x,f7.3)

c     calculate K internally if user sets K= -1 in .INP input file
      if (dp.le.0.d0) then
         dp=d50(isoil)
      endif
      if ((ek.lt.0.d0)) then
c --  English units, ekeng [(0.1ton/ft-ton*)*(h/in)] 
         ekeng=tf(isoil)*(12.d0-om)+sf(isoil)+pf(isoil)
c --  metric units, ek [kg/N*h/m^2] in 2 steps:
c --  a) (Eng->metric)  eK=1.313* ekeng  (USDA Ag. Handbook 537 Supp, 1981)
c --  b) (Tm/ha->Kg/m2) 0.1*ek 
         ek=1.313d0*0.1d0*ekeng
      end if

c-------------
c S factor (-)
c-------------
      theta=datan(Y)
      s=dsin(theta)
c ** Usle
c      bigS=65.4d0*s**2+4.56d0*s+0.065d0
c from haan p261
       if (s.lt.0.09) then
          bigS=10.8*s +0.03
       else
          bigS=16.8*s -0.5
       endif
       if (pl.lt.0.7) then
          bigS=3.0*s**0.8 +0.56
       endif

c-------------
c L factor after McCool, p262 Haan (-)
c-------------
c      if (x.lt.3.d0) then
c         x=0.3d0
c       elseif (x.eq.4.d0) then
c         x=0.4d0
c       else
c         x=0.5d0
c      endif
c * use distance not length along slope
      slopeL = pl*cos(theta)
      beta=11.16d0*s/(3.d0*s**0.8d0+0.56)
      x=beta/(1.d0+beta)
      bigL=(slopeL/22.d0)**x

c-------------
c Final soil loss calculations Ai(kg/m2)
c-------------
      A0=er*ek*bigL*bigS*cfact*pfact
      A1=er1*ek*bigL*bigS*cfact*pfact
      A2=geim*ek*bigL*bigS*cfact*pfact
      A3=erCoolm*ek*bigL*bigS*cfact*pfact
c----For PRZM erosion calculations English units are needed, convert K if provided by user
      ekeng=ek/(1.313d0*0.1d0)
c -- PRZM MUSLE (Williams,1975), Xew. Note: PRZM manual units error, should be Tm/ha-event
      Xew=1.586d0*(volro*qpdepth)**0.56d0*A**0.12d0*ekeng*bigL*bigS*
     & cfact*pfact
c --   convert Tm/ha to kg/m2 with 1 day=1 event
      A4=0.1d0*Xew
c -- PRZM MUSS (Kinnel,2004) Xem. Note: PRZM manual units error, should be Tm/ha-event
      Xem=0.79d0*(volro*qpdepth)**0.65d0*A**0.009d0*ekeng*bigL*bigS*
     & cfact*pfact
c --   convert Tm/ha to kg/m2 with 1 day=1 event
      A5=0.1d0*Xem

c-------------
c Runoff sediment concentration in (g/L) for VFSMOD .ISD file
c-------------
      if(volro.le.0)then
          sconc(1)=0.d0
          sconc(2)=sconc(1)
          sconc(3)=sconc(1)
          sconc(4)=sconc(1)
          sconc(5)=sconc(1)
          sconc(6)=sconc(1)
        else
          sconc(1)=A0*A*10000.d0/vol
          sconc(2)=A1*A*10000.d0/vol
          sconc(3)=A2*A*10000.d0/vol
          sconc(4)=A3*A*10000.d0/vol
          sconc(5)=A4*A*10000.d0/vol
          sconc(6)=A5*A*10000.d0/vol
      endif

      if (ieroty.eq.0) then
            musle_type='Foster Rw (MUSLE)'
        elseif (ieroty.eq.1) then
            musle_type='Williams Rw (MUSLE)'
        elseif (ieroty.eq.2) then
            musle_type='GLEAMS Rw (MUSLE)'
        elseif (ieroty.eq.3) then
            musle_type='Cooley Rw (MUSLE)'
        elseif (ieroty.eq.4) then
            musle_type='PRZM (MUSLE)'
        elseif (ieroty.eq.5) then
            musle_type='PRZM (MUSS)'
        else
            musle_type='Williams Rw (MUSLE)'
      endif

c-----------------------
c  print summary
c-----------------------
      write(10,19) soilty,ek,ekeng,om,dp
      write(10,21) slopeL,beta,bigL,bigS,cfact,pfact
      write(10,23)
c*** convert kg/m^2 to Eng Units ton/ac to show also on table
      ccc=0.00110231131d0/0.000247105381d0
      write(10,26) 'Rm (Foster)  ',A0,A0*ccc,sconc(1),er,ek,bigL,bigS,
     1                    cfact,pfact
      write(10,25) 'Rw (Williams)',A1,A1*ccc,sconc(2),er1
      write(10,25) 'Rm (GLEAMS)  ',A2,A2*ccc,sconc(3),geim
      write(10,25) 'Rst (Cooley) ',A3,A3*ccc,sconc(4),ercoolm
      write(10,27)
      write(10,28) 'PRZM (MUSLE) ',A4,A4*ccc,sconc(5),ekeng
      write(10,28) 'PRZM (MUSS)  ',A5,A5*ccc,sconc(6),ekeng
      write(10,30) musle_type
19    format(/,2x,'For the selected soil type: ',a20,
     1  /,5x,'K=',f10.3,' kg-h/N-m^2',2x,' Eng. K=',f10.3,
     2       2x,'% Org Matter=',f6.1,
     3  /,4x,' dp (d50)=',f12.2,' um')
21    format(/,/,1x,'MISCELLANEOUS CALCS:',/,
     5  /,4x,' SlopeL  =',f10.3,' m',2x,'Beta    =',f10.3,
     6  /,4x,' L-Factor=',f10.3,4x,'S-Factor=',f10.3,
     7  /,4x,' C-fact  =',f10.2,4x,'P-fact  =',f10.2)
23    format(/,/,1x,'FINAL CALCS:',
     1   /,/,2x,'Method',16x,'Soil Loss A',10x,'Sediment ',2x,'R-Factor'
     2   ,2x,'K-Factor',2x,'L-Factor',1x,'S-Factor',1x,'C-Factor',1x,
     3   'P-Factor',/,2x,18x,'kg/m^2',4x,'EngUnits t/ac',2x,'Conc g/l',
     4    3x,'N/ha',6x,'kg-N/N-m^2',/,110('-'))
25    format(
     1  1x,A13,2x,f10.2,4x,f10.2,4x,f10.2,2x,f10.2)
26    format(
     1  1x,A13,2x,f10.2,4x,f10.2,4x,f10.2,2x,f10.2,2x,f7.3,
     2    4(1x,f8.3))
27    format(110('-'),/,2x,18x,'kg/m^2',4x,'EngUnits t/ac',2x,
     1    'Conc g/l',3x,'N/ha',6x,'0.01ton.ac.h/(ac.ft-ton.in)')
28    format(
     1  1x,A13,2x,f10.2,4x,f10.2,4x,f10.2,11x,f10.3)
30    format(/,70('='),/,1x,'Method selected for VFSMOD .isd file:',
     1  2x,A20,/,70('='))

      return
      end
