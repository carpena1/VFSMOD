      program uh
C---------------------------------------------------------------
c      version 3.0.9, Last Modified: See Modifications below
C      WRITTEN FOR: ASAE'99 Toronto paper, March 8, 2002
C      Written by: R. Munoz-Carpena (rmc)   &   J. E. Parsons, BAE (jep)
C                  University of Florida        BAE, NC State University
C                  Gainesville, FL 32611        Raleigh, NC 27695-7625(USA)
C                  e-mail: carpena@ufl.edu
C---------------------------------------------------------------
C Program to create input files for VFSmod, based on NRCS-TR55 and Haan
C et al, 1996, with additional work done on coefficients for unit peak
C flow calculation.
c
c    Date      Modification                                Initials
c   -------    ------------------------------                ------
c   2/17/99    Check for 0.1<Ia/P<0.5                            rmc
c   2/18/99    Added hyetograph output for 6 h storm             jep
c   2/18/99    Modify File Inputs for Erosion                    jep
c   2/20/99    Roughed in MUSLE                                  jep
c   3/01/99    Checked erosion parameters and units              rmc
c   3/02/99    Additional work on Musle - units close            jep
c   3/03/99    Added hyetographs for storm types I & IA          rmc
c   3/05/99    output irs file for VFSMOD                        jep
c   3/06/99    Input/Output files as in VFSMOD                   jep
c   3/10/99    Checked Input/Output files as in VFSMOD           rmc
c   3/10/99    Cleanup - created hydrograph.f for                jep
c              hydrograph subroutines, created io.f for          jep
c              input and output related processing               jep
c   3/28/99    Erosion part: fixes in I30 calculation            rmc
c              after Chow and checked for consistency in         rmc
c              units, clean up; Hydro: added delay time          rmc
c   8/27/99    Added option to select different methods          rmc
c              for applying MUSLE, default is Foster,            rmc
c              2=Williams, 3=GLEAMS                              rmc
c   10/01/99   Fixed array so that storm duration (D)            rmc
c              can now be up to 24h                              rmc
c   10/26/99   implemented the project file concept as in vfsm   jep
c    3/09/00   Version changed to 0.9, general program cleanup   rmc
c   16/06/00   Version changed to 1.0, erosion output organized  rmc
c   16/03/02   Version changed to 1.06 to couple with VFSMOD,    jep
c              author affiliation changed                        rmc
c    4/18/03   Fixed K - computed if we enter -1, other use      jep
c              entered value, also fixed dp output format        jep
c    4/19/03   dp now being read in                              jep
c    4/20/03   Runoff calculation for low CN revised             rmc
c    5/01/03   Added check for small runoff case to switch       rmc
c              to Williams sediment calculation that includes    jep
c              runoff.                                           jep
c   11/10/03   Reordered Erosion ieroty 1=Williams, 2=Gleams     jep
c                3=Foster to coincide with changes in Shell      jep
c   11/13/03   Fixed coef. on Type Ia - did not add new          rmc
c              hyet curves                                       rmc
c   01/10/05   Added changes suggested by U. of Guelph group     rmc
c              v2.4.1                                            rmc
c   09/15/11   Rewritten hydrograph calculation with convolution rmc
c              of excess rain steps, v3.0.0                      rmc
c   02/15/12   Added user table for 24-h hyetograph, v3.0.1      rmc
c   02/15/14   Fix user hyetograph option (jstype=5), v3.0.2     rmc
c   03/12/20   Fixed printing tp on the final hydrograph & added rmc
c              debug to print*,' full tabular hydrograph v3.0.3  rmc
c   8/10/20   Reordered Erosion ieroty 0=Foster; (none or) 1=    rmc
c              Williams, 2=Gleams to make Williams recommended,  rmc
C              3.0.4                                             rmc
c   24/11/20   a) Fixed i30 calculation, fixed typo in USLE K in rmc
c              structure factor data for vfSL texture; b) added  rmc
c              user option ("storm type"=5 in .INP file) to read rmc
C              actual cumulative hyetograph provided by user     rmc
c              hyetograph (hr,mm). v3.0.5                        rmc
c   15/06/21   a) Added check for duration of the storm (and     rmc
c              user normalized cumulative hyetograph) cannot     rmc
c              exceed 24 hr; b) simplified .iro, irn VFSMOD      rmc
c              outputs to  write 15-min hydro- and hydrographs   rmc
c              using moving average, mass balanced maintained).  rmc
C              v3.0.6                                            rmc
c   20/02/23   a) Added option to calculate sediment using MUSS  rmc
c              used in PRZM (ieroty=5 in .inp file); b) 5-min    rmc
c              hyeto-/hydrographs  are now written in .irn and   rmc
c              .iro VFSMOD files (option to select duration of   rmc
c              step is added as last input in first line of INP  rmc
c              but if not provided dincr=5min).  v3.0.7          rmc
c   20/02/23   Fixed processing issues with user hyetographs     rmc
c              (jstype=5,6) affecting runoff, rainfall and       rmc
c              sediment calculations and outfiles. Issue remains rmc
c              with Intel compiler integer/real conversions that rmc
c              affects array bounds when selecting compiler      rmc
c              optimization (-O or higher). The code must be     rmc   
c              compiled with option '-O -fp-model precise' to    rmc
c              ensure correct floating point preserve array      rmc
c              bounds. More in next version. v3.0.8              rmc
c   02/04/25   Fixed bug in mass balance in .iro and .irn where  rmc
c              the values for total rain and runoff reported on  rmc
c              the .out file do not match the integral under the rmc
c              irn and irn curves.    v3.0.9                     rmc
c                                                                rmc
c---------------------------------------------------------------
c    Compiling for Win32 and Unix environments9
c       1. The i/o for these operating systems is different.
c       2. Change the Unix/Win32 comments in the finput.f program 
c          to reflect your operating system.   3/9/00
c---------------------------------------------------------------
c common/hydgph:
c       rot(208), runoff time (units)
c       roq(208), runoff rate (m3/s)
c       u(208,2), unit hydrograph
c common/rain/:
c       rfix, maximum rain intensity (mm/h)
c       rti(200), rainfall time (hrs)
c       rfi(200), rainfall intensity (mm/h)
c       rcum(100,2), cumm rainfall (mm)
c       ref(100), excess rainfall intensity (mm/h)
c       ncum: number of steps if user hyetograph is read
c other:
c       nref=number of excess hyetograph steps
c       mref=number of unit hydrograph steps
c       nhyet=number of hyetograph steps
c       vol(m3), volro (mm)=runoff volume
C---------------------------------------------------------------
      implicit double precision (a-h, o-z)
      character*20 soilty
      CHARACTER*75 LISFIL(6)
      dimension sconc(6)

      common/hydgph/u(10000,2),qh(10000,2)
      common/rain/rfix,rti(10000),rfi(10000),rcum(10000,2),ref(10000,2),
     1 ncum

c------------------------------------------
c  get inputs and open files
c------------------------------------------
      call FINPUT(LISFIL)
      call getinp(P,CN,A,jstype,D,pL,Y,ek,cfact,pfact,soilty,
     1            ieroty,dp,om,dincr)
c--------------------------
c Calculate runoff volume by SCS method
c--------------------------
      call runoff(P,CN,xIa,Q)
      volro=Q
c--------------------------
c Calculate concentration time by SCS method
c--------------------------
      call calctc(pL,CN,Y,tc)
c--------------------------
c Calculate peak flow and time by SCS-TR55 method
c--------------------------
      call q_peak(A,Q,xIa,P,tc,jstype,qp,tp)
c--------------------------
c Output hydrology results
c--------------------------
      call results(P,CN,Q,A,tc,xIa,jstype,D,pL,Y,qp,tp,qdepth,
     1             ieroty)
c--------------------------
c Calculate SCS-unit hydrograph
c--------------------------
      call unit_hyd(Q,A,qp,tp,D,tc,qp5,tp5,mref)
      vol=volro*A*10000.d0/1000.d0

c--------------------------
c Calculate storm hyetograph from SCS storm type
c--------------------------
       call hyetgh(jstype,P,D,volro,qdepth,vol,qp,A,xIa,rtpeak,er,er1,
     1            erCoolm,ti,nref,tc,a1,b1,bigE,raimax30,nhyet)

c--------------------------
c Calculate storm hydrograph
c--------------------------
      call tab_hyd(Q,A,mref,nref,qp,nhyd)

c--------------------------
c do the modified usle to get erosion stuff
c--------------------------
      call musle(er,er1,erCoolm,ek,Y,pl,cfact,pfact,A,vol,tc,P,D,soilty,
     1           dp,ieroty,sconc,om,a1,b1,bigE,raimax30,qp)

c------------------------------------
c write vfsmod compatible input files
c------------------------------------
      call vfsout(dp,ieroty,sconc,A,pL,qp,tp,tc,D,ti,
     1           nhyet,nhyd,dincr)

      close(1)
      close(2)

      stop
      end
