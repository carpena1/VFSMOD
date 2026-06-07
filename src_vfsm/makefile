# ─────────────────────────────────────────────────────────────────────────────
# Fortran compiler
# ─────────────────────────────────────────────────────────────────────────────
# RECOMMENDED: gfortran (macOS Apple Silicon / any Linux) — see LDFLAGS below
FC=gfortran
# ALTERNATIVE: ifort (Intel macOS/Linux only) — activate with matching LDFLAGS
#FC=ifort
#FC=f77
#FC=g77
#FC=xlf_r
#FC=xlf
#FC=nagfor

# LIBRARY NAMES GIVEN BELOW
LIBRARIES=

# ─────────────────────────────────────────────────────────────────────────────
# Linker  (keep in sync with FC above)
# ─────────────────────────────────────────────────────────────────────────────
LD=gfortran
#LD=ifort
#LD=f77
#LD=g77
#LD=xlf_r
#LD=xlf
#LD=nagfor

# ─────────────────────────────────────────────────────────────────────────────
# Compiler / linker flags
#
# BENCHMARKED: VFSMOD CR=0.6, KPG=1, optimal N, 11 L values [2-50 m]
#   2026-06-03  Intel Core i7-8569U @ 2.80 GHz  macOS Sonoma  (ifort 2021.10 / gfortran 15.2)
#   2026-06-06  Apple M5 Pro  ARM64  macOS Sequoia              (gfortran 15.2)
#
#   Machine             Build       Flags                                 med ms/run  Speedup
#   ─────────────────────────────────────────────────────────────────────────────────────────
#   Intel i7-8569U      ifort  -O3  -O3 -fp-model precise                  12.51 ms   1.00x  (ref)
#   Intel i7-8569U      ifort +xHst -O3 -xHost -fp-model precise           12.43 ms   1.01x
#   Intel i7-8569U      gfort  -O3  -O3 -march=native                      11.70 ms   1.07x  ← fastest on Intel
#   Intel i7-8569U      gfort +fmath -O3 -march=native -ffast-math         12.03 ms   1.04x
#   Apple M5 Pro        gfort  -O3  -O3 -march=native                      14.18 ms   0.88x  ← use on M5 Pro
#   Apple M5 Pro        gfort +fmath -O3 -march=native -ffast-math         14.19 ms   0.88x
#
# Note: at optimal N (5-13 nodes), wall time is dominated by process launch and
# file I/O (~13-15 ms); actual compute is <1 ms. The M5 Pro advantage will be
# much larger for compute-bound runs (large N, many iterations).
# All builds produce bit-identical output (max|Δq| = 0.00 m³/s) vs ifort ref.
# -ffast-math adds no benefit on either platform and is not recommended.
# ─────────────────────────────────────────────────────────────────────────────

# ── gfortran (Apple Silicon M-series / Linux) ─────────────────────────────────
# RECOMMENDED: -O3 -march=native  (safe, portable; fastest on Apple Silicon)
# -ffpe-summary=none  suppresses the end-of-run "IEEE_UNDERFLOW_FLAG" note that
#   gfortran's runtime prints when an underflow-to-zero occurred during the run.
#   The underflow itself is harmless (result == 0.0); ifort is silent about it.
LDFLAGS= -O3 -march=native -ffpe-summary=none -ffixed-line-length-132
# With loop unrolling (marginal benefit, no accuracy risk):
#LDFLAGS= -O3 -march=native -ffpe-summary=none -funroll-loops -ffixed-line-length-132
# Debug build (show ALL FPE signals for diagnostics):
#LDFLAGS= -g -fcheck=all -ffpe-summary=all -ffixed-line-length-132

# ── ifort (Intel macOS/Linux only — not available on Apple Silicon) ────────────
# RECOMMENDED: -O3 -fp-model precise  (accurate, safe, production)
# -Wl,-ld_classic  suppresses deprecation warnings on macOS Sonoma/Sequoia
#LDFLAGS= -O3 -fp-model precise -Wl,-ld_classic -extend-source 132
# With CPU-native vectorisation (+1% speed, same accuracy):
#LDFLAGS= -O3 -xHost -fp-model precise -Wl,-ld_classic -extend-source 132
# Debug build:
#LDFLAGS= -g -traceback -check -fp-model precise -Wl,-ld_classic -extend-source 132

# ── Other compilers (legacy) ──────────────────────────────────────────────────
# NAG fortran:
#LDFLAGS= -O -ieee=stop -132
# IBM xlf:
#LDFLAGS= -O -qmaxmem=-1 -fp-model precise

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
 multspcalc.f\
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

# ── Portable release build (gfortran) ─────────────────────────────────────────
# Links the gfortran runtime statically so the binary depends only on the system
# libc (libSystem on macOS) and needs no Homebrew/gcc runtime at the destination.
# This is the recipe used to produce executables/vfsm4621.   Usage:  make release
STATIC= -static-libgfortran -static-libquadmath -static-libgcc
release : $(OBJS)
	$(LD) -o vfsm $(LDFLAGS) $(STATIC) $(OBJS) $(LIBRARIES)

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
