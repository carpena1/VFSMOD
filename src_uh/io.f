      subroutine getinp(P,CN,A,jstype,D,pL,Y,ek,cfact,pfact,soilty,
     1                  ieroty,dp,om,dincr)
C---------------------------------------------------------------
C       Input parameters
C---------------------------------------------------------------
      implicit double precision (a-h, o-z)
      common/rain/rfix,rti(10000),rfi(10000),rcum(10000,2),ref(10000,2),
     1 ncum
      character*20 soilty

c----------------------------------
c read rainfall amount, watershed characteristics
c   P = rainfall amount (mm)
c   CN = SCS curve number
c   A = watershed area (ha)
c   jstype = storm type (1=I;2=II;3=III;4=Ia;5=user cum curve (P/P24);6=user1 real cum. curve (P))
c   D = storm duration (h)
c   pL = maximum flow path length (m)
c   Y = slope (%)
c   dincr =  (missing set=5) time step of hyet-/hydrographs in IRN, IRO (min)
c----------------------------------
      read(1,*,iostat=ierr)P,CN,A,jstype,D,pL,Y,dincr
      if (ierr.ne.0) then
            dincr=5.d0
            rewind(1)
            READ(1,*,iostat=ierr)P,CN,A,jstype,D,pL,Y
            if (ierr.ne.0) then
                 print*,'ierr=',ierr
                 print*,'ERROR: Check first line of INP file'
                 print*,''
                 STOP
            end if
      end if
      if(D.gt.24.d0) then
            print*,'ERROR: duration of the storm must be D =< 24 h.'
            print*,'       Please revise .inp file to fix this and',
     1                     ' rerun.'
            STOP
      endif
      read(1,'(A)')
c---------------------------------
c read inputs for soil erosion calculations
c   soilty = soil type (Character)
c   ek = soil erodibility
c   cfact= C factor
c   pfact= P factor
c------------------------------------
      read(1,'(A)')soilty
      read(1,*) ek, cfact, pfact, dp
c-      convert dp to um
       dp=dp*10000.d0
c----------------------------------
c-  ieroty = select method to estimate storm erosion
c-           Selections:
c-         0 = Foster's method for R-factor
c-         1 or not present = Using Williams R-factor (recommended)
c-         2 = Using R-factor from GLEAMS with daily rainfall
c-         3 = Using R-factor from Cooley with daily rainfall
c-         4 = Using R-factor from PRZM MUSLE
c-         5 = Using R-factor from PRZM MUSS
      read(1,*,END=22) ieroty
      if ((ieroty.lt.6).and.(ieroty.ge.0)) go to 24
  22  ieroty=1
  24  continue
c------------------------------------------------
c read om = % organic matter (always read even if not used when k=-1)
      read(1,*,end=32)om
  32  continue
c------------------------------------------------
c-- For user defined case, read "tmid" first and then 24-h P/P24 curve (jstype=5)
c-- or 24 P (mm) (jstype=6). "tmid" is stored in first position of the rcum(i,j)
c-- array. For jstype=6  "user1" option, tmid = 1/2 D (i.e. middle of hyetograph).
      if(jstype.gt.4) then
         read(1,*,end=40)rcum(1,1)
         if(jstype.eq.6) rcum(1,1)=D*0.5D0
         if (rcum(1,1).eq.0.or.rcum(1,1).gt.24) then
           if(jstype.eq.5) then
               print*,'ERROR: the first value must be tmid(h), followed'
     &        ,'in the next lines by the cumulative rainfall that must'
     &        ,'begin with (0,0) and end with (24,1)'
               print*,''
           endif
           STOP
         endif
         rcum(1,2)=0.5d0
         do 35 i=2,10000
           read(1,*,end=40)(rcum(i,j),j=1,2)
c rmc-- new 'user1' (jstype=6) option for reading real cum. hyetographs
           if(jstype.eq.6)rcum(i,2)=rcum(i,2)/P
35       continue
40       ncum=i-1
         if(jstype.eq.5.and.((rcum(ncum,1).ne.24).or.
     &       (rcum(ncum,2).ne.1))) then
            print*,'ERROR: the cumulative rainfall must begin with ',
     &      '(0,0) and end with (24,1)'
            print*,''
            STOP
           else if(jstype.eq.6.and.((rcum(ncum,1).ne.D).or.
     &       (rcum(ncum,2).ne.1))) then
            print*,'ERROR: the cumulative rainfall must begin with ',
     &      '(0,0) and end with (D,P)'
            print*,''
            STOP
         endif
      endif

      return
      end

      subroutine results(P,CN,Q,A,tc,xIa,jstype,D,pL,Y,qp,tp,qdepth,
     1                   ieroty)
C---------------------------------------------------------------
C      Output summary of hydrology results nicely
C---------------------------------------------------------------
      implicit double precision (a-h, o-z)
      character*5 stype(6)
      data stype/'I','IA','II','III','user','user1'/

      write(2,*)' '
      write(2,*)' HYDROGRAPH CALCULATION FOR WATERSHED-SCS METHOD'
      write(2,*)' '
      write(2,*)'INPUTS'
      write(2,*)'------'
      write(2,100)P
      write(2,200)stype(jstype)
      write(2,250)D
      write(2,300)CN
      write(2,400)A
      write(2,500)pL
      write(2,600)Y*100.d0
      write(2,610)ieroty
      write(2,*)' '
      write(2,*)'OUTPUTS'
      write(2,*)'-------'
      write(2,700)Q, Q*A*10.d0
      write(2,800)xIa
      write(2,900)tc,tc*60.d0
      write(2,*)' '

100   format('Storm Rainfall=',f8.2,' mm')
200   format('SCS storm type= ',a5)
250   format('Storm duration=',f6.1,' h')
300   format('SCS Curve number=',f6.1)
400   format('Watershed area=',f8.2,' ha')
500   format('Maximum flow path length=',f8.2,' m')
600   format('Average slope of flow path=',f8.2,' %')
610   format('MUSLE type=',i3,' where R-factor is (see manual)',/,
     1    2x,'0=Foster, 1=Williams,   2=GLEAMS',/,
     2    2x,'3=Cooley, 4=PRZM-MUSLE, 5=PRZM-MUSS')
700   format('Runoff volume=',f8.2,' mm=',f8.2,' m3')
800   format('Initial Abstraction=',f8.2,' mm')
900   format('Concentration time=',f8.2,' h= ',f8.2,' min')

      return
      end

      subroutine vfsout(dp,ieroty,sconc,A,pL,qp,tp,
     1                  tc,D,ti,nhyet,nhyd,dincr)
c--------------------------------
c      Output for VFSMOD input files
c--------------------------------
      implicit double precision (a-h,o-z)
      common/hydgph/u(10000,2),qh(10000,2)
      common/rain/rfix,rti(10000),rfi(10000),rcum(10000,2),ref(10000,2),
     1 ncum
      dimension sconc(6)
      logical ende
      character*72 dumline

c--------------------------------
c Output of VFSMOD input file: *.isd
c--------------------------------
      npart=7
      coarse=0.5d0
      cci= sconc(ieroty+1)/1000.d0

c--rmc  05/08/03 when soil loss for problem is too large with the erosion
c----- method selected, print a warning to suggest alternative methods known
c----- to produce lower values (Williams, MUSS) that consider runoff and
c----  typically avoids this problem.
      if(cci.ge.0.25d0.and.(ieroty.ne.1.and.ieroty.ne.5)) then
c         cci=sconc(2)/1000.d0
         write(*,160)ieroty
         write(10,*)
         write(10,160)ieroty
      endif
      por=0.434d0
      write (15,101) Npart,coarse,cci,por
      dpp=dp/10000.d0
      sg=2.65d0
      write (15,102)dpp,sg

c--------------------------------
c Output of VFSMOD runoff hydrograph: *.iro
c--------------------------------
c  Output hydrograph based on 5-, 15, 30-min or 60-min steps (dincr) selected by
c  user on INP file or if missing set = 5 min.
      dincr=dincr/60.d0
c  Use scratch file to store hyetograph and replace final peak rainfall later
      open(3,status='scratch')
      swidth=A*10000.d0/pL
      slength=pL
      write (3,103) swidth,slength
      nbcroff=nhyd
      bcropeak=qp
      dtstep=qh(nhyd,1)-qh(nhyd-1,1)
      nwrite1=int(dincr/dtstep)
      if (nwrite1.eq.0) nwrite1=1
      nstep1=int(nhyd/nwrite1)
      bcropeak=0.d0
      write(2,250)ti
      if(nhyd.le.nstep1) then
            nsteps=nbcroff+1
            write (3,104)nsteps,bcropeak
            do 20 ii=1, nbcroff-1
               tt=qh(ii,1)*3600.d0
               if (ii.eq.1) then
                  write(3,105)tt,qh(ii,2),dincr*60.d0
                 else
                  write(3,106)tt,qh(ii,2)
                endif
                bcropeak=max(bcropeak,qh(ii,2))
20          continue
         else
            nsteps=nstep1+3
            write(3,104)nsteps,bcropeak
            write(3,105)qh(1,1)*3600.d0,qh(1,2),dincr*60.d0
            nsum=0
            tsum=0.d0
            tmean=0.d0
            qsum=0.d0
            qmean=0.d0
            do 29 ii=2,nhyd
               tsum=tsum+qh(ii,1)*3600.d0
               qsum=qsum+qh(ii,2)
               nsum=nsum+1
               do 25 k=1,nstep1
                  if(ii.eq.k*nwrite1) then
c       print point values every nwrite1 steps
c                       write(3,106)tt,qh(ii,2)
c       print average values every nwrite1 steps
                        tmean=tsum/nsum
                        qmean=qsum/nsum
                        write(3,106)tmean,qmean
                        tsum=0.d0
                        qsum=0.d0
                        nsum=0
c       end of print selection
                  endif
                  bcropeak=max(bcropeak,qmean)
25             continue
29          continue
      endif
c---replace bcropeak from the final hydrograph at selected dincr steps
      rewind(3)
      read(3,*)dummy1,dummy2
      write (12,103)dummy1,dummy2
      read(3,*)dummy1,dummy2
      write(12,104)nsteps,bcropeak
      read(3,*)dummy1,dummy2
      write(12,105)dummy1,dummy2,dincr*60.d0
      ende = .FALSE.
      do while (.NOT.ende)
         read (3,*,end=30)dummy1, dummy2
         write(12,106)dummy1, dummy2
      enddo
30    continue
      close(3)
c---write 0 entry after last step
      tend1=qh(nhyd,1)*3600.d0
      write(12,106)tend1,qh(nhyd,2)
      write(12,107)tend1+300.d0,0.d0

c----------------------------------------
c Output VFSMOD rainfall hyetograph: *.irn
c----------------------------------------
c  Output hydrograph based on 5-, 15, 30-min or 60-min steps (dincr) selected as above
c  Use scratch file to store hyetograph and replace final peak rainfall later
      open(3,status='scratch')
      ndev=0
      dtstep=rti(nhyet)-rti(nhyet-1)
      nstep2=nint(rti(nhyet)/dincr)
      nwrite2=nint(real(nhyet)/real(nstep2))
      if(nstep2*nwrite2.gt.nhyet) then
        ndev=nint(real(nstep2*nwrite2-nhyet)/real(nwrite2))
        nstep2=nstep2-ndev
      endif
      if (nwrite2.eq.0) nwrite2=1
      rfix=0.d0
      if(nhyet.le.nstep2)then
            nsteps=nhyet+2
            write(3,201)nsteps,rfix
            write(3,203) rti(1)*3600.d0,rfi(1),dincr*60.d0
            do 31 ii=2,nhyet-1
               write(3,204)rti(ii)*3600.d0,rfi(ii)
               rfix=max(rfix,rfi(ii))
31          continue
         else
            nsteps=nstep2+3
            write(3,201)nsteps,rfix
            write(3,203) rti(1)*3600.d0,rfi(1),dincr*60.d0
            nsum=0
            tsum=0.d0
            tmean=0.d0
            rsum=0.d0
            rmean=0.d0
            do 33 ii=2,nhyet
               tsum=tsum+rti(ii)*3600.d0
               rsum=rsum+rfi(ii)
               nsum=nsum+1
               do 32 k=1,nstep2
c       print point values every nwrite2 steps
                  if(ii.eq.k*nwrite2) then
c       calculate moving average values every nwrite2 steps
                        tmean=tsum/nsum
                        rmean=rsum/nsum
                        write(3,204)tmean,rmean
                        tsum=0.d0
                        rsum=0.d0
                        nsum=0
c       end of print selection
                  endif
                  rfix=max(rfix,rmean)
32             continue
33          continue
c     write 0 entry after last step
            write(3,204)rti(nhyet)*3600.d0,0.d0
      endif
c     extend the duration of the run tend 5 more minutes (300 s)
      write(3,204)rti(nhyet)*3600.d0+300.d0,0.d0

c---replace peak rain from the final hyetograph at selected dincr steps
      rewind(3)
      read(3,*)dummy1,dummy2
      write(14,201)nsteps,rfix
      read(3,*)dummy1,dummy2
      write(14,203)dummy1,dummy2,dincr*60.d0
      ende = .FALSE.
      do while (.NOT.ende)
         read (3,*,end=40)dummy1, dummy2
         write(14,204)dummy1, dummy2
      end do
40    continue
      close(3)
c----------
c Output message at end of program -----------------
c----------
      WRITE(*,*)
      WRITE(*,*)'...FINISHED...','UH v3.0.8 09/2023'
      WRITE(*,*)

c-------------------
c  format statements
c-------------------
101   format (2x,i4,2x,f8.1,2x,f11.5,2x,f7.4,8x,
     1   'Npart, Coarse, Ci(g/cm3), Por')
102   format (2x,f10.7,2x,f7.1,21x,'Dp(cm), SG(g/cm3)')
103   format (2x,f7.1,2x,f7.1,21x,'Swidth(m), Slength(m)')
104   format (2x,i4,2x,e12.5,19x,'nbcroff, bcropeak (m3/s)')
105   format (2x,e12.5,2x,e12.5,10x,' time(s), ro(m3/s)-',f4.0,
     1   'min steps')
106   format (2x,e12.5,2x,e12.5)
107   format (2x,e12.5,2x,e12.5,/,30('-'))
160   format(72('-'),/,'WARNING: This case produces large sediment',
     1    ' concentration with',/,'the erosion method IEROTY#',i4,
     2    '  selected. Try Williams (1) or',/,'MUSS (5) and rerun',
     3    ' -- Please',' see manual',/,72('-'))
201   format(i4,2x,e12.5,20x,' NRAIN, RPEAK(m/s)')
203   format(2x,e12.5,3x,e12.5,10x,'time(s), rainfall(m/s)-',f4.0,
     1   'min steps')
204   format(2x,e12.5,3x,e12.5)
205   format(2x,e12.5,3x,e12.5,/,30('-'))
250   format('Time to ponding=',f8.3,' h')
260   format('Duration of rainfall excess=',f8.3,' h')
270   format('Time to peak after shifting=',f8.3,' h')
280   format('Time correction to match hyetograph=',f8.3,' h')

      return
      end
