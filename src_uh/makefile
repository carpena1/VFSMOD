### fortran make file for uh-2
# version as of 2/21/99 for $(FC)

# Fortran compiler given below
FC=ifort
#FC=gfortran
#FC=f77
#FC=g77

# LINKER GIVEN BELOW
LD=ifort
#LD=gfortran
#LD=f77
#LD=g77

OBJS = uh.o \
       finput.o \
       io.o \
       hydrograph.o \
       musle.o \
       hyetograph.o

SRCS = uh.f \
       finput.f \
       io.f \
       hydrograph.f \
       musle.f \
       hyetograph.f

# COMPILER FLAGS GIVEN BELOW

FTOPT = -O -fp-model precise
#FTOPT = -g -traceback -check
#FTOPT = -O2 -fp-model precise
#FTOPT = -fast -fp-model precise
#FTOPT = -O3 -fp-model precise
#FTOPT = -O -qmaxmem=-1 -fp-model precise
#FTOPT = -V -cross_reference -g -v -check_bounds
#FTOPT = -O
#FTOPT = -v -Qpath /ncsu/nagf90/lib -I/ncsu/nagf90/lib -L/ncsu/nagf90/lib
#FTNCOPT = -check_bounds
#DBGOPT = -check_bounds -g
OPT3 = -align records -align dcommons -O3
OPT4 = -v -align records -align dcommons -O4
OPT5 = -v -align records -align dcommons -O5
#LDOPT =  -Qpath /ncsu/nagf90/lib -L/ncsu/nagf90/lib

# INCLUDE SEARCH PATH GIVEN BELOW
INCLUDES=

.SUFFIXES: .uid .uil .o .c .f .for .F .FOR .hpf .mpl .p .pas .P .PAS .cob .cbl .COB .CBL .prj

uh: $(OBJS)
	$(FC) $(LDOPT) $(OBJS) -o uh

uh.o: uh.f
	$(FC) $(FTOPT) -c uh.f

finput.o: finput.f
	$(FC) $(FTOPT) -c finput.f

io.o: io.f
	$(FC) $(FTOPT) -c io.f

hydrograph.o: hydrograph.f
	$(FC) $(FTOPT) -c hydrograph.f

musle.o: musle.f
	$(FC)  $(FTOPT) -c musle.f

hyetograph.o: hyetograph.f
	$(FC)  $(FTOPT) -c hyetograph.f

debug: $(SRCS)
	$(FC) -v $(DBGOPT) $(SRCS) -o uh-db

opt3:  $(SRCS)
	$(FC) $(OPT3) $(SRCS) -o uh-3
	rm -f *.u

opt4:  $(SRCS)
	$(FC) $(OPT4) $(SRCS) -o uh-4
	rm -f *.u

opt5:  $(SRCS)
	$(FC) -O5  $(SRCS) -o uh-5
	rm -f *.u

clean:
	rm -f *.o core uh
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
