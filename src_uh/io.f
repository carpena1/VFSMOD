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
c   Output for VFSMOD input files
c   dt = dincr*60 time step of hyet-/hydrographs in IRN, IRO (s)
c   qh(i,j)= hydrograph (i=time step, j=1 time (h), j=2 flow (m3/s))
c   rti(i),rfi(i)= hyetograph time (s) and instensity (m/s)(i=time step)
c   nhyd= total number of hydrograph steps
c   nhyet= total rnumber of hyetograph steps
c   predm(i,j)=interpolated hydro/hyetograph (i=time step, j=1 time (s), j=2 flow/rain)
c--------------------------------
      implicit double precision (a-h,o-z)
      common/hydgph/u(10000,2),qh(10000,2)
      common/rain/rfix,rti(10000),rfi(10000),rcum(10000,2),ref(10000,2),
     1 ncum
      dimension sconc(6),predm(10000,2)
      logical ende
      character*72 dumline
      character*1 adum

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
c  user on .INP input file (if missing default dincr= 5 min).

c  Interpolate original hydrograph from convolution to user desired time step 
      dt=dincr*60.d0
c--- debug: show the full synthetic hydrograph before interpolation ---
c      print*,dt,nhyd,qh(nhyd,1)*3600.d0
c      do 2 k=1,nhyd
c          write(*,106)qh(k,1)*3600.d0,qh(k,2)
c2     continue
c      write(*,206)   
c---end debug----
      time0=0.d0
      i=0
      bcropeak=0.d0
      do while (time0.le.qh(nhyd,1)*3600.d0)
            i=i+1
            if (i.eq.1) then
                  time0=qh(i,1)*3600.d0
                  predm(i,1)=time0
                  predm(i,2)=qh(i,2)
              else
                  do 10 j=1,nhyd-1
                     if(time0.gt.qh(j,1)*3600.d0) then
                        predm(i,2)=(time0/3600.d0-qh(j,1))*(qh(j+1,2)
     &		       -qh(j,2))/(qh(j+1,1)-qh(j,1)) + qh(j,2)
                        predm(i,1)=time0
                        goto 10
                     endif
10		      continue
		endif
            bcropeak=max(bcropeak,predm(i,2))
            time0=time0+dt
  	end do
c --to conserve the full graph, add last step of the original hydrograph 
      predm(i+1,1)=qh(nhyd,1)*3600.d0
      predm(i+1,2)=qh(nhyd,2)
      nsteps=i+1
 
c --Write *.iro file with the interpolated hydrograph
      swidth=A*10000.d0/pL
      slength=pL
      write (12,103) swidth,slength
      write (12,104)nsteps,bcropeak
      do 13 ii=1, nsteps
            if (ii.eq.1) then
               write(12,105)(predm(ii,j),j=1,2),dincr
              else
               write(12,106)(predm(ii,j),j=1,2)
             endif
 13   continue
      write(12,206)

c----------------------------------------
c Output VFSMOD rainfall hyetograph: *.irn
c----------------------------------------
c  Output hyetograph based on 5-, 15, 30-min or 60-min steps (dincr) selected by
c  user on .INP input file (if missing default dincr= 5 min).

c  Interpolate original hyetograph to the user desired time step 
c--- debug: show the full synthetic hyetograph before interpolation---
c      print*,dt,nhyet,rti(nhyet)*3600.d0
c      do 15 k=1,nhyet
c          write(*,106)rti(k)*3600.d0,rfi(k)
c15     continue
c      write(*,206)   
c---end debug---
      time0=0.d0
      i=0
      rpeak=0.d0
      do while (time0.le.rti(nhyet)*3600.d0)
            i=i+1
            if (i.eq.1) then
                  time0=rti(i)*3600.d0
                  predm(i,1)=time0
                  predm(i,2)=rfi(i)
              else
                  do 20 j=1,nhyet-1
                     if(time0.gt.rti(j)*3600.d0) then
                        predm(i,2)=(time0/3600.d0-rti(j))*(rfi(j+1)
     &		       -rfi(j))/(rti(j+1)-rti(j)) + rfi(j)
                        predm(i,1)=time0
                        goto 20
                     endif
20		      continue
		endif
            rpeak=max(rpeak,predm(i,2))
            time0=time0+dt
      end do
c --to conserve the full graph, add last step of the original hydrograph with rain=0 
      predm(i+1,1)=time0
      predm(i+1,2)=0.d0
c -- extend the duration of the run "tend" 5 more minutes (300 s) to allow for hydrograph to finish
      predm(i+2,1)=predm(i+1,1)+300.d0
      predm(i+2,2)=0.d0
      nsteps=i+2
 
c --Write *.iro file with the interpolated hydrograph
      write(14,201)nsteps,rpeak
      do 25 ii=1, nsteps
            if (ii.eq.1) then
               write(14,203)(predm(ii,j),j=1,2),dincr
              else
               write(14,106)(predm(ii,j),j=1,2)
             endif
25    continue
      write(14,206)

c----------
c Output message at end of program -----------------
c----------
      WRITE(*,*)
      WRITE(*,*)'...FINISHED...','UH v3.0.9 09/2024'
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
107   format (3e12.5)
160   format(72('-'),/,'WARNING: This case produces large sediment',
     1    ' concentration with',/,'the erosion method IEROTY#',i4,
     2    '  selected. Try Williams (1) or',/,'MUSS (5) and rerun',
     3    ' -- Please',' see manual',/,72('-'))
201   format(i4,2x,e12.5,20x,' NRAIN, RPEAK(m/s)')
203   format(2x,e12.5,2x,e12.5,10x,'time(s), rainfall(m/s)-',f4.0,
     1   ' min steps')
206   format(30('-'))

      return
      end