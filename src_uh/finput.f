       SUBROUTINE FINPUT(LISFIL)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Create input and output file names from a command line input string    C
C     NOTE: Maximum lenght of command line string = 50                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*50 FILENM1
      CHARACTER*75 LISFIL(6)
      CHARACTER*3,SCOD(6)
      CHARACTER*1,DUMMY1
      character*200 linein
      character*1 slash
      DATA(SCOD(I),I=1,6)/'inp','out','hyt','iro','irn','isd'/

c*** Command line option to input filename
c*** Comment out the following depending for which system you compile

CWIN32*** Start of Win32 file i/o     ***
CWIN32      slash='\'
CWIN32      INARGS=NARGS()-1
CWIN32      IF (INARGS.EQ.1) THEN
CWIN32          CALL GETARG(1,FILENM1,IFSTATUS)
CWIN32*** End of Win32 file i/o       ***
cUNIX *** Start Unix file i/o         ***
CUNIX
      slash='/'
CUNIX
      INARGS=IARGC()
CUNIX
      IF (INARGS.EQ.1) THEN
CUNIX
         CALL GETARG(1,FILENM1)
cUNIX*** End of UNIX file i/o section ***
       ELSE
          WRITE(*,*)
          WRITE(*,105)
          WRITE(*,110)
          WRITE(*,130)
          WRITE(*,140)
          WRITE(*,150)
          WRITE(*,*)
          STOP
      ENDIF

c-----Write welcome message ---------------------------------------
      write(*,*)
      WRITE(*,160)
      WRITE(*,*)'                 @    @ @    @'
      WRITE(*,*)'                 @    @ @    @'
      WRITE(*,*)'                 @    @ @@@@@@'
      WRITE(*,*)'                 @    @ @    @'
      WRITE(*,*)'                  @@@@  @    @  Sept. 2023-v3.0.8'
      WRITE(*,160)
      WRITE(*,*)'        R.Munoz-Carpena              J.E. Parsons'
      WRITE(*,*)'          UFL - USA                   NCSU - USA'
      WRITE(*,*)'       carpena@ufl.edu'
      WRITE(*,160)
           WRITE(*,*)'      PROGRAM GENERATE RAINFALL, AND SOURCE AREA'
      WRITE(*,*)'      RUNOFF AND SEDIMENT INPUTS FOR VFSMOD.'
      WRITE(*,160)
c      WRITE(*,*)
C------ create I/O filenames from input string -------------------------
c------    or read filenames from a project file -----------------------

      ilstr=index(filenm1,'.')
      if (ilstr.gt.0) then
c      *** using project file (.prj or .lis) to read filenames
c      *** mods made 10/27/99, jep - push version to 1.0
c      ***  check to see if extension is .prj or .lis
          ilstr1=index(filenm1,'.prj')
          ilstr2=index(filenm1,'.lis')
          if ((ilstr1.gt.0).or.(ilstr2.gt.0)) then
c          *** fill filename array with safe names
             do 11 i=1,6
               dummy1=scod(i)
                   IF(DUMMY1.EQ.'i') THEN
                      WRITE(LISFIL(I),'(5A)')
     1            'inputs',slash,'dummy','.',SCOD(I)
                ELSE
                  WRITE(LISFIL(I),'(5A)')
     1            'output',slash,'dummy','.',SCOD(I)
               ENDIF
   11        continue
c
             open(unit=99,file=filenm1,status='old')
   12        read(99,'(a)',end=18) linein
             lpos=index(linein,'=')
             lstr=len(linein)
             if ((lpos.gt.0).and.(lstr.gt.0)) then
                do 14 jj=1,6
                   lpp = index(linein(1:lpos-1),scod(jj))
                   if (lpp.gt.0) lisfil(jj)=linein(lpos+1:)
   14           continue
             endif
             go to 12
c           ****** done
   18       continue
          else
            WRITE(*,*)
            WRITE(*,105)
            WRITE(*,110)
            WRITE(*,130)
            WRITE(*,140)
            WRITE(*,150)
            WRITE(*,*)
            STOP
          endif
      else
c     **** rafa's i/o scheme
        ILSTR=INDEX(FILENM1,' ')-1
        DO 101 I=1,6
            DUMMY1=SCOD(I)
            IF(DUMMY1.EQ.'i') THEN
                  WRITE(LISFIL(I),'(5A)')
     1        'inputs',slash,FILENM1(:ILSTR),'.',SCOD(I)
               ELSE
                  WRITE(LISFIL(I),'(5A)')
     1        'output',slash,FILENM1(:ILSTR),'.',SCOD(I)
            ENDIF
101     CONTINUE
      endif

       write(*,*)'*** Opening '
      DO 102 I=1,6
c            write(*,*)'*** Opening ',lisfil(i)
            write(*,'(70A)')lisfil(i)
102      CONTINUE
      WRITE(*,*)

C---------Open I/O files -------------------------------
      OPEN(1,FILE=LISFIL(1),STATUS='OLD')
      OPEN(2,FILE=LISFIL(2),STATUS='unknown')
      OPEN(10,FILE=LISFIL(3),STATUS='unknown')
      OPEN(12,FILE=LISFIL(4),STATUS='unknown')
      OPEN(14,FILE=LISFIL(5),STATUS='unknown')
      OPEN(15,FILE=LISFIL(6),STATUS='unknown')
      WRITE(2,220)LISFIL(2)
      WRITE(10,220)LISFIL(3)
      WRITE(2,*)

105   FORMAT('Name:    uh')
110   FORMAT(9x,'(Generate rainfall and runoff inputs for VFSMOD)')
130   FORMAT('Usage:   uh filename (max 8 characters or project name)')
CWIN32 identifier for the simulation
CWIN32 140   FORMAT('Version: v3.0.8 for Windows -09/2023')
cUNIX identifier for the simulation
cUNIX
140   FORMAT('Version: v3.0.8 for Unix -09/2023')
150   FORMAT('Authors: R.Munoz-Carpena & J.E.Parsons (UFL & NCSU)')
160   FORMAT(72('-'))
220   FORMAT('File: ',A40,9x,'UH v3.0.8, 09/2023')

      RETURN
      END
