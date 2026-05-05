              SUBROUTINE WQSUB(IWQ,TIME,N)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                             C
C   Water quality component skeleton for other pollututant processes. This    C
C   is provided to plug-in research modules (not available in the general     C
C   for management.                                                           C
C   See for example multi-reactive component based on TARSE:                  C
C   Pérez-Ovilla, O. 2010. A flexible numerical component to simulate         C
C     biogeochemical transport processes through vegetative filter strips.    C
C     Ph.D. dissertation. [Gainesville, Fla.]: University of Florida. URL     C
C     https://tinyurl.com/rx7hxrrh                                            C
C   Yu, C., R. Muñoz-Carpena, B. Gao and O. Perez-Ovilla. 2013. Effects of    C
C     ionic strength, particle size, flow rate, and vegetation type on        C
C     colloid transport through a dense vegetation saturated soil system:     C
C     Experiments and modeling. J. of Hydrology 499:316–323.                  C
C     doi:10.1016/j.jhydrol.2013.07.004.                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (MAXEQN=1001,MAXBND=40)
      IMPLICIT DOUBLE PRECISION   (A-H,O-Z)

C-------Read in main parameters of the program--------------

c      WRITE(18,*)IWQ,TIME,N

      RETURN
      END
