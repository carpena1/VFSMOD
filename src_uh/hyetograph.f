      Subroutine hyetgh(jstype,P,D,volro,qdepth,vol,qp,A,xIa,rtpeak,er,
     1          er1,erCoolm,ti,nref,tc,a1,b1,bigE,raimax30,nhyet)
c !----------------------------------------------------
c ! Where P'(t)= is the cumulative hyetograph for the given duration
c !       Pd is the total rainfall for the given period (mm)
c !       D is the storm duration in hours
c !       t is time in hours
c !
c ! storm type II or III - from Haan
c ! hyetograph for 24 hour storms
c ! P(t)        T (  24.04  ) ^0.75
c ! ---- = 0.5+---(---------)
c ! P24        24 (2|T|+0.04)
c !
c !  where T=t -12  (with t in hours)
c !        P24 is the 24h storm
c !
c ! storm type I  - fitted from SCS tabular data - rmc 03/04/99
c !
c ! hyetograph for 24 hour storms
c !                   (     -0.1617    ) ^0.5853
c !       | 0.4511+ T (----------------)       ;for [-3.0163|T|+0.013]<0
c ! P(t)  |           (-3.0163|T|+0.013)
c ! ---- =|
c ! P24   |
c !       |   0.5129                           ;for [-3.0163|T|+0.013]>0
c !
c !  where T=t -9.995  (with t in hours)
c !
c ! storm type IA  - fitted from SCS tabular data - rmc 03/04/99
c !
c ! hyetograph for 24 hour storms
c ! P(t)             (     0.0843     ) ^0.4228
c ! ---- = 0.3919+ T (----------------)
c ! P24              (120.39|T|+0.3567)
c !
c !  where T=t -7.96  (with t in hours)
c !
c ! For any storm of any duration (from Haan et. al.(1994), eq. 3.7)
c !
c ! P'(t)    P(tmid+t-D/2) - P(tmid-D/2)
c ! ----- = --------------------------
c ! Pd       P(tmid+D/2) - P(tmid-D/2)
c !
c ! where where tmid=12. Alternatively (Munoz-Carpena and Parsons,2004),
c ! tmid=12.00 for storm type II & III, tmid=9.995 for storm type I, and
c ! tmid=7.960 for storm type IA
c !
c !----------------------------------------------------------------
c  Storm type I and IA - fitted equations from tabular data on Haan's
c !----------------------------------------------------------------

      implicit double precision (a-h, o-z)
      character*5 stype(6)
      common/rain/rfix,rti(10000),rfi(10000),rcum(10000,2),ref(10000,2),
     1 ncum
      dimension rainh(10000),rainh30(10000)
      data stype/'I','IA','II','III','user','user1'/

c----rmc 24/3/99 - hyetograph using unit hydrograph time step, Def
      pd = P
      Def=0.24d0*tc
      ndtime=nint(D/Def+0.5D0)
      dtime=D/ndtime
      ir30win=nint(0.5d0/dtime+0.5D0)
c-----
      pcumtot=0.d0
      refcum=0.d0
      iflag=0
      ti=0.d0
      bigE=0.d0
      raimax=0.d0
      raimax30=0.d0
      nref=0

c----rmc v3 09/2011 -- For each storm type (except for new user hyetographs)
c--- calculate scaling factors for D<24 h, based on eq. 3-7, Haan et.al. (1994)
c**> set Cooly (1980) a1,b1 coef. For user and user1 types it uses Type II coef.

      if (jstype.eq.1) then
         tmid=9.995d0
         a1 = 15.03d0
         b1 = 0.578d0
       else if(jstype.eq.2) then
         tmid=7.96d0
         a1 = 12.98d0
         b1 = 0.7488d0
       else if(jstype.eq.3) then
         tmid=11.8d0
         a1 = 17.9d0
         b1 = 0.4134d0
       else if(jstype.eq.4) then
         tmid=12.d0
         a1 = 21.51d0
         b1 = 0.2811d0
       else if(jstype.eq.5) then
         tmid=rcum(1,1)
         a1 = 17.9d0
         b1 = 0.4134d0
       else
         tmid=rcum(1,1)
         a1 = 17.9d0
         b1 = 0.4134d0
      endif
c---  scaling factors for other storm durations and volumen
      tminus=tmid-D*.5d0
      if(tminus.lt.0.d0) tminus=0.d0
      tplus=tmid+D*.5d0
      if(tplus.gt.24.d0) tplus=24.d0
      pm=SCStorm(jstype,tminus)
      pp=SCStorm(jstype,tplus)
c---------------------
c---  Main loop for rain time step calculations
c---  NOTE: the first position of the rain array is used to store {tmid, 0.5}
c---        so the actual hyetograph goes from i=2,ndtime+1
c---------------------
      do 3 i=2,ndtime+1
        smalle=0.d0
        rti(i)=(i-1)*dtime
        tsmall=tmid+rti(i)-D*.5d0
        if(tsmall.gt.24.d0) then
           tsmall=24.d0
         elseif(tsmall.lt.0.d0) then
           tsmall=0.d0
        endif
        ptp=SCStorm(jstype,tsmall)
c---------------------
c  Calculate cumulative hyetograph for any duration and volume
c---------------------
        cumtotal=pd*(ptp-pm)/(pp-pm)
c        STOP

        if(cumtotal.gt.xIa.and.iflag.eq.0) then
           iflag=1
           ti=(rti(i)-rti(i-1))/(cumtotal-pcumtot)*
     1            (xIa-pcumtot)+rti(i-1)
        endif
c----------------------------------------------------------------
c  Calculate instantaneous hyetograph and rainfall energy term for USLE
c----------------------------------------------------------------
        rainh(i)=cumtotal-pcumtot
        if (rainh(i).gt.0.d0) then
c ---> english units ft-tons/acre-inch
           if ((rainh(i)/25.4d0/dtime).gt.3.d0) then
             smalle=1074.d0
           else
             smalle=(rainh(i)/25.4d0)*
     1        (916.d0+331.d0*dlog10(rainh(i)/25.4d0/dtime))
           endif
c           smalle=
c     1      1099.d0 * (1.d0-0.72*exp(-1.27*(rainh(i)/25.4d0/dtime)))
c           print*,rti(i),smalle
           bigE=bigE+smalle
c ---> metric units
c           bigE=bigE+11.9d0+8.73d0*dlog10(rainh(i)/dtime)
        endif
        if (rainh(i).gt.raimax) then
           raimax=rainh(i)
           rtpeak=rti(i)
        endif
c---rmc - I30 calculation after Chow et al, 1987
        if(i.gt.ir30win) then
           rainh30(i)= rainh(i)
           do 4 i30max=1,ir30win
            rainh30(i)= rainh30(i)+rainh(i-i30max)
4          continue
        endif
        if (rainh30(i).gt.raimax30) then
           raimax30=rainh30(i)
           rtpeak30=rti(i)
        endif
        pcumtot=cumtotal
        rfi1=rainh(i)/dtime
        rfi(i)=rfi1/3600.d0/1000.d0
c---rmc-08/24/11-- excess rainfall hyetograph for tabular hydrograph
        if(cumtotal.ge.xIa) then
           nref=nref+1
           ref(nref,1)=rti(i)
           ref(nref,2)=(cumtotal-xIa)**2.d0/(cumtotal+4.d0*xIa)-refcum
           refcum=(cumtotal-xIa)**2.d0/(cumtotal+4.d0*xIa)
        endif
3     continue
c---rmc-08/24/11-- number of hyetograph steps
      nhyet=i-1
c----rmc 03/11/99---
      rfix=raimax/(dtime*3600.d0)/1000.d0
      rfix30=raimax30*2.d0
      rI30=rfix30/25.4d0

c----rmc 08/24/11-- done hyet - computing musle param
c---------------------
c  compute R for musle
c---------------------
c**convert bigE to SI metric - multiply by 1.702 / 100
c** units Rst=N/h
      bigEm=0.006700d0*bigE
      rst=1.702d0*(bigE/100.d0)*(raimax30/25.4d0/dtime)*dtime/0.5d0
      rro=volro*(qdepth)**(1.d0/3.d0)
c-- a)  er=> Foster et al. 1977b, units N/h
      er=0.5d0* rst + 0.35d0*rro
c--rmc 03/28/99- er1 = 9.05d0*(vol*qp)**0.56d0
c-- b)  er1=> Williams, units Mg h/ha N
      er1 = 9.05d0*(vol*qp)**0.56d0/A
c-- c) erCooly=> Cooley (1980), er for design storms
c--              EI/100 = R ft     tonsf/ac
c**              EI/100   R*1/0.67 J/m^2
        erCooly=a1*(P/25.4d0)**(2.119d0*rti(ndtime)**0.0086d0)
     1           /(rti(ndtime)**b1)
        erCoolm=erCooly*1.702d0
c---------------------
c results
c---------------------
      write(10,5)Def*60.d0,nhyet
5     format(/,1x,'SCS ',f4.1,'-MIN HYETOGRAPH (25 of',i5,
     1  1x,'steps printed)',/,/,
     2  2x,'No.',3x,'Time(hr)',3x,'Rainfall(mm)',1x,'Rain30(mm)',
     3  1x,'Eff.rain(mm)')

c---print hyetograph results 25 time steps only
      maxstep=24
      nwrite=nint(real(ndtime)/real(maxstep)+0.5D0)
      crainh=0.d0
      cref=0.d0
      iref=nhyet-nref
      do 8 i=2,ndtime+1
            crainh=crainh+rainh(i)
            if(i.gt.iref) cref=cref+ref((i-iref),2)
            do 7 k=1,maxstep
               if(i.eq.k*nwrite) then
                  write(10,10)k,rti(i),crainh,rainh30(i),cref
               endif
7           enddo
8     enddo
      write(10,10)k,rti(ndtime+1),crainh,rainh30(ndtime+1),cref
10    format(2x,i3,2x,f8.2,3x,3f10.2)

      write(10,15) cumtotal,P,refcum,raimax30,rfix30
15    format(/,2x,'Computed Total Rain =',f10.1,' mm'/,
     1   2x,'  Actual Total Rain =',f10.1,' mm',/,
     2   2x,'  Total Rain Excess =',f10.1,' mm',/,
     3   2x,'  raimax30          =',f10.1,' mm',/,
     4   2x,'  I30               =',f10.1,' mm/h')

      return
       end

       function SCStorm(Jstype,ptime)
c !----------------------------------------------------
c ! SCS design storm type equation using generalized coefficients
c ! (Munoz-Carpena and Parsons,2004) for 24 hour storms,
c !     P(t)       t-b (    d    ) ^g
c !     ---- = a + ----(---------)
c !     P24         c  (e|t-b|+f )
c !----------------------------------------------------
       implicit double precision (a-h,o-z)

      common/rain/rfix,rti(10000),rfi(10000),rcum(10000,2),ref(10000,2),
     1 ncum
      dimension cff(4,7)
      data cff/0.4511d0,0.3919d0,0.495d0,0.5d0,9.995d0,7.96d0,11.8d0,
     1   12.d0,1d0,1.d0,0.56d0,24.d0,-0.1617d0,0.843d0,10.6d0,24.04d0,
     2   -3.0163d0,120.39d0,130.d0,2.d0,0.013d0,0.3567d0,0.525d0,
     3   0.04d0,0.5853d0,0.4228d0,0.75d0,0.75d0/

      SCStorm=0.d0
      if(Jstype.le.4) then
          cffa=cff(Jstype,1)
          cffb=cff(Jstype,2)
          cffc=cff(Jstype,3)
          cffd=cff(Jstype,4)
          cffe=cff(Jstype,5)
          cfff=cff(Jstype,6)
          cffg=cff(Jstype,7)
          bigT= ptime-cffb
          denom=cffe*dabs(bigT)+cfff
          if(Jstype.eq.1.and.denom.ge.0.d0) then
            SCStorm=0.5129d0
           else
            SCStorm=cffa+(bigT/cffc)*(cffd/denom)**cffg
          endif
        else
          do 15 i=2,ncum-1
            t1=rcum(i,1)
            rcum1=rcum(i,2)
            t2=rcum(i+1,1)
            rcum2=rcum(i+1,2)
            if(ptime.ge.t1.and.ptime.le.t2) then
             SCStorm=(ptime-t1)/(t2-t1)*(rcum2-
     &         rcum1)+rcum1
            endif
15        continue
      endif
      
c100   format(7f9.4)

      return
      end
