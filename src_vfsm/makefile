# Fortran compiler given below
FC=ifort
#FC=f77
#FC=g77
#FC=xlf_r
#FC=xlf
#FC=gfortran

# LIBRARY NAMES GIVEN BELOW
LIBRARIES=

# LINKER GIVEN BELOW
LD=ifort
#LD=f77
#LD=g77
#LD=xlf_r
#LD=xlf
#LD=gfortran

# LINKER FLAGS GIVEN BELOW
# Intel Compiler: optimization with vectorization
LDFLAGS= -O3
#LDFLAGS= -fast
#LDFLAGS= -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
#LDFLAGS= -O0
#LDFLAGS= -qthreaded -O5 -qtune=auto -qarch=auto -qunroll=auto
#LDFLAGS= -qsmp -qthreaded -O5 -qtune=auto -qarch=auto -qunroll=auto
#LDFLAGS= -O -qmaxmem=-1
# ifort debug
#LDFLAGS= -g -traceback -check
# gfortran: optimization with vectorization
#LDFLAGS= -O3 -ftree-vectorize -ffixed-form -ffixed-line-length-none

# DEFINES
AR=ar
ARFLAGS=rv
CP=cp
RM=rm -f

# COMPILER FLAGS GIVEN BELOW
CFLAGS=$(LDFLAGS)
PFLAGS=$(CFLAGS)
FFLAGS=$(CFLAGS)


# INCLUDE SEARCH PATH GIVEN BELOW
INCLUDES=

.SUFFIXES: .uid .uil .o .c .f .for .F .FOR .hpf .mpl .p .pas .P .PAS .cob .cbl .COB .CBL .prj

# ALL SOURCES GIVEN BELOW
SRCS= vfsmod.f\
 assm.f\
 bca.f\
 cde.f\
 conver.f\
 einstein.f\
 elem.f\
 factor.f\
 finput.f\
 flow.f\
 forma.f\
 formb.f\
 gasub.f\
 gasubwt.f\
 grassed.f\
 grassin.f\
 ini.f\
 inputs.f\
 kwwrite.f\
 modify.f\
 ocf.f\
 outmass.f\
 points.f\
 quad.f\
 shapef.f\
 solve.f\
 step3.f\
 update.f\
 wqpest.f\
 wqsub.f

# OBJECTS GIVEN BELOW
OBJS= $(SRCS:.f=.o)

vfsm : $(OBJS)
	$(LD) -o $@ $(LDFLAGS) $(OBJS) $(LIBRARIES)

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $<

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDES) $<

.for.o:
	$(FC) -c $(FFLAGS) $(INCLUDES) $<

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDES) $<

.FOR.o:
	$(FC) -c $(FFLAGS) $(INCLUDES) $<

.f.a:
	$(FC) -c $(FFLAGS) $(INCLUDES) $<
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o

.for.a:
	$(FC) -c $(FFLAGS) $(INCLUDES) $<
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o

.F.a:
	$(FC) -c $(FFLAGS) $(INCLUDES) $<
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o

.FOR.a:
	$(FC) -c $(FFLAGS) $(INCLUDES) $<
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o

clean:
	$(RM) $(OBJS) $(PRJS) core vfsm

touchsrcs:
	touch $(SRCS)

# tell make how to create a .prj file from a .f file
.f.prj:
	ftnchek -nopretty -project -noextern -columns=132  $<

# set up macro PRJS containing project filenames
PRJS= $(OBJS:.o=.prj)

# "make check" will check everything that has been changed.
check: $(PRJS)
	ftnchek $(PRJS)

cleanck:
	$(RM) $(PRJS) core check
