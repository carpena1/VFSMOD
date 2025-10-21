      SUBROUTINE FINPUT(LISFIL,ISCR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Create input and output file names from a command line input string  C
C     NOTE: Maximum length of command line string = 25 characters      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C List of file units to open:                                          C
C     UNIT 1  LISFIL(1)   inputs/*.ikw                                 C
C     UNIT 2  LISFIL(2)   inputs/*.irn                                 C
C     UNIT 3  LISFIL(3)   inputs/*.iro                                 C
C     UNIT 7  LISFIL(4)   inputs/*.iso                                 C
C     UNIT 10 LISFIL(5)   output/*.osm                                 C
C     UNIT 11  LISFIL(6)  output/*.ohy                                 C
C     UNIT 12  LISFIL(7)  inputs/*.igr                                 C
C     UNIT 13  LISFIL(8)  output/*.og1                                 C
C     UNIT 14  LISFIL(9)  output/*.og2                                 C
C     UNIT 15 LISFIL(10)  output/*.osp                                 C
C     UNIT 16 LISFIL(11)  inputs/*.isd                                 C
C     UNIT 17 LISFIL(12)  inputs/*.iwq                                 C
C     UNIT 18 LISFIL(13)  output/*.owq                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*120 FILENM1
      CHARACTER*75 LISFIL(13)
      character*200 linein
      CHARACTER*3,SCOD(13)
      CHARACTER*1,DUMMY1,SISCR
      character*1 slash
      DATA(SCOD(I),I=1,13)/'ikw','irn','iro','iso','osm','ohy','igr',
     &  'og1','og2','osp','isd','iwq','owq'/

c*** Command line option to input filename
c*** Comment out the following depending for which system you complile
CDOS      slash='\'
CUNIX
      slash='/'
c***End of UNIX/DOS selection***********
CF2003*** New fortran 2003 standard for command line input ***
      iargs = command_argument_count()
      IF (IARGS.EQ.1) THEN
          CALL get_command_argument(1, FILENM1)
          ISCR=0
c-----Write welcome message ---------------------------------------
          WRITE(*,*)'   @     @ @@@@  @@@  @   @  @@@@  @@@'
          WRITE(*,*)'   @     @ @    @     @@ @@  @  @  @  @'
          WRITE(*,*)'    @   @  @@@   @@@  @ @ @  @  @  @   @'
          WRITE(*,*)'     @ @   @        @ @   @  @  @  @  @'
          WRITE(*,*)'      @    @     @@@  @   @  @@@@  @@ 10/2025 v4.6.1.1'
          WRITE(*,160)
          WRITE(*,*)'       R.Munoz-Carpena              J.E. Parsons'
          WRITE(*,*)'       U. of FL - USA               NCSU - USA'
          WRITE(*,*)'       carpena@ufl.edu        john_parsons@ncsu.edu'
          WRITE(*,*)
          WRITE(*,*)' Licensed under CC BY-ND 4.0'
          WRITE(*,*)' https://abe.ufl.edu/faculty/carpena/vfsmod/'
          WRITE(*,160)
          WRITE(*,*)'PROGRAM TO CALCULATE OVERLAND FLOW AND SEDIMENT',
     &      ' FILTRATION THROUGH A'
          WRITE(*,*)'VEGETATIVE FILTER STRIP OF AN INFLOW HYDROGRAPH',
     &      ' FROM AN ADJACENT FIELD,'
          WRITE(*,*)'DURING A STORM EVENT. VFSMOD HANDLES THE CASE OF',
     &      ' VARYING SURFACE COVER'
          WRITE(*,*)'AND SLOPES AT THE NODES AND TIME DEPENDENT',
     &      ' INFILTRATION FOR THE DOMAIN.'
          WRITE(*,160)
          WRITE(*,*)
       ELSEIF (IARGS.EQ.2) THEN
          CALL get_command_argument(1, FILENM1)
          CALL get_command_argument(2, SISCR)
          READ(SISCR, '(I)', IOSTAT=IOERR) ISCR
           IF (IOERR .NE. 0) THEN
            write(*,'(/,A50,/)')'ERROR: The second command argument is not a valid number'
            STOP
          ENDIF
          IF (ISCR.LT.1.OR.ISCR.GT.2) THEN
            write(*,'(/,A50,/)')'ERROR: The second command argument can only be 1-2'
            STOP
          ENDIF
       ELSE
          WRITE(*,*)
          WRITE(*,105)
          WRITE(*,110)
          WRITE(*,120)
          WRITE(*,130)
          WRITE(*,140)
          WRITE(*,150)
          WRITE(*,*)
          STOP
      ENDIF
c
C------ create I/O filenames from input string -------------------------
c------   or read filenames from project file  -------------------------

      ilstr=index(FILENM1,'.')
      if (ilstr.gt.0) then
c     *** using project file (.prj) to read filenames
c     ***   mods made 10/24/99, jep - push version to 1.04
c     ***  check for .prj or .lis to make sure this is our file
         ilstr1=index(FILENM1,'.prj')
         ilstr2=index(FILENM1,'.lis')
         if ((ilstr1.gt.0).or.(ilstr2.gt.0)) then
c     *** fill filename array with safe names
           do 11 i=1,13
                DUMMY1=SCOD(I)
                  if (dummy1.eq.'i') then
                      WRITE(LISFIL(I),'(5A)')
     1                   'inputs',slash,'dummy','.',SCOD(I)
                   else
                      WRITE(LISFIL(I),'(5A)')
     1                   'output',slash,'dummy','.',SCOD(I)
              endif
cdebug      print *,'Debug: i=',i,':',lisfil(i)
11               continue
c
            open (unit=99,file=FILENM1,iostat=istat,status='old')
            if (istat .ne. 0) then
               write(*,*)'ERROR: Project file cannot be opened:',FILENM1
               write(*,*)
               stop
            end if
12          read(99,'(a)',end=18) linein
            lpos = index(linein,'=')
            lstr = len(linein)
cdebug      print *,'lpos,lstr=',lpos,lstr,':',linein
            if ((lpos.gt.0).and.(lstr.gt.0)) then
                  do 14 jj=1,13
                   lpp=index(linein(1:lpos-1),scod(jj))
                   if (lpp.gt.0) lisfil(jj)=linein(lpos+1:)
14                continue
            endif
            go to 12
c       **** done
18      continue
cdebug        print *,(i,lisfil(i),i=1,13)
        else
          WRITE(*,*)
          WRITE(*,105)
          WRITE(*,110)
          WRITE(*,120)
          WRITE(*,130)
          WRITE(*,140)
          WRITE(*,150)
          WRITE(*,*)
          STOP
        endif

      else
c     **** rafa's new i/o scheme for filenames
        ILSTR=INDEX(FILENM1,' ')-1
        IF(ILSTR.LT.1)ILSTR=8
        DO 101 I=1,13
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
cdebug      print*,(i,lisfil(i),i=1,13)

C---------Open I/O files -------------------------------
c------------Inputs ----------
      OPEN(1,FILE=LISFIL(1),iostat=istat,STATUS='OLD')
      if (istat .ne. 0) then
            write(*,*)'ERROR: Project file cannot be opened:',LISFIL(1)
            write(*,*)
            stop
      end if
      OPEN(2,FILE=LISFIL(2),iostat=istat,STATUS='OLD')
      if (istat .ne. 0) then
            write(*,*)'ERROR: Project file cannot be opened:',LISFIL(2)
            write(*,*)
            stop
      end if
      OPEN(3,FILE=LISFIL(3),iostat=istat,STATUS='OLD')
      if (istat .ne. 0) then
            write(*,*)'ERROR: Project file cannot be opened:',LISFIL(3)
            write(*,*)
            stop
      end if
      OPEN(7,FILE=LISFIL(4),iostat=istat,STATUS='OLD')
      if (istat .ne. 0) then
            write(*,*)'ERROR: Project file cannot be opened:',LISFIL(4)
            write(*,*)
            stop
      end if
      OPEN(16,FILE=LISFIL(11),iostat=istat,STATUS='OLD')
      if (istat .ne. 0) then
            write(*,*)'ERROR: Project file cannot be opened:',LISFIL(11)
            write(*,*)
            stop
      end if
      OPEN(12,FILE=LISFIL(7),iostat=istat,STATUS='OLD')
      if (istat .ne. 0) then
            write(*,*)'ERROR: Project file cannot be opened:',LISFIL(7)
            write(*,*)
            stop
      end if
c------------Outputs ----------
      OPEN(10,FILE=LISFIL(5),STATUS='UNKNOWN')
      WRITE(10,220)LISFIL(5)
      OPEN(11,FILE=LISFIL(6),STATUS='UNKNOWN')
      WRITE(11,220)LISFIL(6)
      WRITE(11,*)
      OPEN(13,FILE=LISFIL(8),STATUS='UNKNOWN')
      WRITE(13,220)LISFIL(8)
      WRITE(13,*)
      OPEN(14,FILE=LISFIL(9),STATUS='UNKNOWN')
      WRITE(14,220)LISFIL(9)
      WRITE(14,*)
      OPEN(15,FILE=LISFIL(10),STATUS='UNKNOWN')
      WRITE(15,220)LISFIL(10)
      WRITE(15,*)

c*** in summary file, put the list of files for this run
      write(15,*) 'Files for this simulation'
      write(15,225) (i,scod(i),lisfil(i),i=1,11)

105   FORMAT('Name:    vfsm')
110   FORMAT(9x,'(VFSMOD model to calculate overland flow and pollutant')
120   FORMAT(9x,'trapping efficiency of Vegetative Filter Strips)')
130   FORMAT('Usage 1: vfsm projectname  (verbose execution)',/,
     &      '         projectname (max 25 chars.) is .prj project name (w/o extension)',
     &    /,'Usage 2: vfsm projectname 1 (quiet execution)',
     &    /,'Usage 3: vfsm projectname 2 (quiet execution and convergence control)')
C----Select unix/DOS versions for help message
Cunix
140   FORMAT('Version: v4.6.1.1 for Unix - 10/2025')
CDOS
c140   FORMAT('Version: v4.6.1.1 for Win64 - 10/2025')
150   FORMAT('Authors: R.Munoz-Carpena & J.E.Parsons (UFL & NCSU)')
160   FORMAT(72('-'))
220   FORMAT('File: ',A40,8x,'VFSMOD v4.6.1.1 10/2025')
225   format(3x,'File #=',i3,' code:',a3,'=',a)

      RETURN
      END
