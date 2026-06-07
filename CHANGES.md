# VFSMOD — Changelog

_Markdown rendering of `src_vfsm/CHANGES.txt` (the authoritative source). Versions are listed newest first._

## 4.6.2.1

**(05/2026)**

- For the very small buffer sizes for SPOTMOD use, five inter-related fixes to the Green-Ampt infiltration and overland-flow ponding logic (gasub.f, vfsmod.f). These bugs were silent — producing zero outflow from the first segment in multi-segment SPOTMOD runs with moderate Ks values — and were diagnosed by tracing Green-Ampt state variables (NPOND, CU, NPFORCE, NEND) through debug output in GASUB and fixed as follows:
  - NEND reset per timestep (vfsmod.f). The end-of-runoff flag NEND was only ever set to 1 (never reset to 0). Once the condition (BCRO=0, XSCHK<= PZERO, NSTART=1) fired, NEND=1 was passed to every subsequent GASUB call. In GASUB the NEND=1 branch was setting NPOND=0 AND CU=-1 unconditionally, permanently forcing Case 1b (no ponding, R=0) for the remainder of the event. Fixed by resetting NEND=0 at the top of the timestep loop in vfsmod.f, so the flag only affects the single timestep when it is actually detected.
  - CU not overwritten in NEND block (gasub.f). Removed CU=-1 from the NEND=1 branch in GASUB. CU is only recomputed when the rainfall period changes (L<>LO). Forcing CU=-1 during the same period (L=LO) prevented the correct Case 1a (CU>0) branch from resuming after the NEND step, even though the soil remained genuinely ponded. With CU preserved, the ponding state is correctly re-evaluated from the value set at the start of the current rainfall period.
  - NPFORCE not reset at end-of-runoff (gasub.f). With NEND now reset to 0 every timestep, the one-timestep NEND=1 signal only cleared NPOND=0 but left NPFORCE=1 active. On every subsequent call, the line IF(NPFORCE.EQ.1.AND.NPOND.EQ.0) NPOND=1 immediately restored NPOND=1, keeping GASUB in Case 2 (ponded GA capacity) with no water available — causing infiltration and wetting-front depth to grow without bound after the storm ended. Fixed by adding NPFORCE=0 inside the NEND=1 block in GASUB, mirroring the reset already present in GASUBWT (line 90).
  - XSCHK spatial-average loop index (vfsmod.f). For the average-check mode (NCHK=-1), the accumulation loop used X(N) at every iteration instead of X(I), so XSCHK equalled N*X(N) regardless of conditions at upslope nodes. Because X(N) (downstream node) is the last to receive surface water, this biased XSCHK toward lower values and could delay the NPFORCE=1 trigger. Fixed by using X(I) to compute the true spatial average.
  - NPFORCE gated by TIME>=TP (vfsmod.f). After correcting the XSCHK average, the upstream BC node X(1)~BCRO>0 immediately raised XSCHK above PZERO at timestep 2, setting NPFORCE=1 before the soil reached its Green-Ampt ponding time TP. This prematurely switched GASUB to Case 2, where C1 = AGA*(TIME-TP+TPP) is negative for TIME<TP-TPP, causing RNEWTON to diverge and producing NaN in the overland-flow solution vector Q0. Fixed by adding the condition TIME>=TP to the NPFORCE trigger so that forced-ponding infiltration is only activated after the soil has reached its rainfall-driven ponding time.
- Updated Petrov-Galerkin (PG) stabilization parameters (inputs.f). For N>=15 nodes the old quartic polynomials in CRR (Munoz-Carpena et al. 1993) are replaced by improved bivariate fits f(CRR,N) selected from 16 candidate model forms using adjusted-R2 and BIC, calibrated from 120 digitized data points per parameter (Figs. 6-9, doi:10.1029/93WR00610):
```
PGPAR(1)=alpha_c: P3(CRR)+Q3(CRR)/N+R3(CRR)/N^2  R2=0.988  adjR2=0.987
PGPAR(2)=alpha_m: P2(CRR)+Q3(CRR)/N+R3(CRR)/N^2  R2=0.996  adjR2=0.996
PGPAR(3)=beta_c:  P3(CRR)+Q3(CRR)/N+R2(CRR)/N^2  R2=0.996  adjR2=0.995
PGPAR(4)=beta_m:  P2(CRR)+Q3(CRR)/N+R3(CRR)/N^2  R2=0.997  adjR2=0.997
```
These improve on the initial bivariate fits (R2 0.975-0.997 with 7-9 parameters) by using richer model forms (11-12 parameters) that better capture the nonlinear N-dependence at large CRR. The previous fits were constrained to Models B/F (Px+Qy/N or Px+Qy/N+Ry/N^2 with P2/Q2); the new search also considers P3 for f-infinity and Q3/R3 for all terms. For N<15 the original 1993 quartic fallback in CRR is retained unchanged.
- Bug fix: hydfilter divide-by-zero for near-step-function inflows (inputs.f). After the low-pass filter pass in hydfilter, if only a single non-zero inflow point survives, the triangular hydrograph width (Ttrend-Ttrini) collapses to zero, causing Qtrpeak = 2*QSUM1/(Ttrend-Ttrini) = Inf. This wrote Inf into BCROFF and silently zeroed all field inflow volume in the restarted simulation (output showed "Volume from up-field = 0"). Fixed by adding a guard that skips the triangular reconstruction (GOTO 65) when Ttrend-Ttrini<=0 or QSUM1<=0, preserving the already-filtered BCROFF.
- Input correction: sample.igr modified Manning n for sediment transport (VN) corrected from 0.02 to 0.012 s/cm^1/3 to match the intended sample scenario.
- Corrected option IDG=0 when several pesticides are used (IWQ>1 in .ikw file). Now it bypasses the degradation (multispc call) for all pesticides selected.
- Changes for gfortran compatibility (benchmarked 2026-06-03: gfortran -O3 -march=native is 7% faster than ifort -O3 -fp-model precise on Intel i7; both produce bit-identical output). Four non-standard ifort extensions removed from the source:
```
finput.f  line  70: READ(…,'(I)',…) bare-I format (no width) → list-directed *
inputs.f  line 653: ((a,b,c,d),j=2,n) nested implied-DO (ifort ext.) → flat (a,b,c,d,j=2,n)
wqpest.f  line 222: 'd' in column 1 (debug-line ifort ext.) → 'c' (standard comment)
wqpest.f  line 456: WRITE(18,'(A)'), LINE trailing comma (ifort ext.) → no comma
```
One additional runtime compatibility fix:
```
grassin.f line  30: ICO (implicit integer) was read directly from list-directed input
where many .igr files write it as 0.0/1.0 (real). ifort silently truncates; gfortran
rejects the token. Fixed by reading into DOUBLE PRECISION XTMP_ICO then ICO=NINT(XTMP_ICO).
```
Runtime note suppressed:
```
gfortran prints "IEEE_UNDERFLOW_FLAG" at exit when any underflow-to-zero occurred
during the run (harmless; ifort is silent). Suppressed with -ffpe-summary=none in
the recommended gfortran LDFLAGS (see makefile).
```
Shared DO-loop termination labels removed (2026-06-06, gfortran 15 warning on M5 Pro):
```
gfortran warns "Fortran 2018 deleted feature: Shared DO termination label" when two or
more nested DO loops share the same label on a single CONTINUE. Fixed in 5 files by
giving each inner loop its own unique label and an explicit CONTINUE statement:
  assm.f:   DO 10/DO 10 → DO 10/DO 11
  elem.f:   DO 10/DO 10 → DO 10/DO 11  (zeroing loop)
            DO 20/DO 20/DO 20 → DO 20/DO 21/DO 22  (integration loop)
  factor.f: DO 10/DO 10/DO 10 → DO 10/DO 11/DO 12
  ini.f:    DO 10/DO 10 → DO 10/DO 11
  solve.f:  DO 10/DO 10 → DO 10/DO 11
Build is now warning-free with gfortran 15 on Apple M5 Pro (ARM64).
```
- Sediment-trapped pesticide mass mfsed reformulated (wqpest.f, 2026-06). For the mechanistic mass-balance trapping (IWQPRO=3) the trapped mass is the phase split mf = dQ.mdi + dE.mpi, so the sediment-sorbed pool is now computed directly as the deposited sorbed fraction mfsed = dE.mpi = PDSED/100*DGMPI(JJ)*SAREA (clamped >=0), replacing the previous equilibrium re-partitioning estimate (Kd*Ctrap*SMtrap and its Freundlich analogue). The old estimator used a different physical model than the dP it was meant to decompose and under-counted the sorbed mass deposited with the trapped sediment for strongly sorbed compounds, triggering the >5% trapped-mass closure warning (e.g. sampleP3 compound 2, Kd=100: residual 16.2% -> 0.0%). The closure warning is retained as an honest CDE/trapping consistency diagnostic.
- CDE secant (chord) retardation for Freundlich isotherms (cde.f, 2026-06). The retardation factor R in the Huang & van Genuchten (1995) analytical leaching solution is now the secant slope of the isotherm, R = 1 + (rho/theta)*Kf*C^(N-1), instead of the tangent slope R = 1 + (rho/theta)*N*Kf*C^(N-1). For a self-sharpening loading front (favorable isotherm, N<1) the front travels at the chord slope by the Rankine-Hugoniot mass-balance shock condition, so the secant is the physically correct, mass-conserving retardation; it reduces exactly to the linear Kd at N=1 (bit-identical results for all linear-isotherm inputs). The fixed-point iteration on R/C1 is unchanged apart from dropping the N factor in both the initial estimate and the loop.
- Consistent sorbed-mass accounting in the CDE profile and mixing layer (cde.f, 2026-06). The nonlinear quadrature (qgausssorb) over the leaching profile is replaced by mfFp = (R-1)*mfFd for the soil profile and mfmlp = (R-1)*mfmld for the mixing layer, using the same secant R that governs transport. This makes the total infiltrated mass exact, mfF = R*mfFd = dQ.mdi, removing the +33% mfF over-count and the associated -6% to -15% closure error that the tangent R produced for samplePF (Freundlich, N=0.85); single-event closure for samplePF is now +0.82% (the irreducible residual from the linear-Kd split used by the dP trapping step). For Freundlich the two changes above must ship together: b1 mfsed with the tangent R would worsen closure; b1 mfsed with the secant R restores it.
- COMMON/WQ3/ extended with the array DGMPI(10), which carries the incoming sediment-sorbed mass mpi (DGPINp, mg/m2) from cde.f to wqpest.f for the mfsed = dE.mpi computation. The block declaration was made identical across all program units that use it (cde.f [declared twice], wqpest.f, multspcalc.f, inputs.f, outmass.f, vfsmod.f) to prevent block misalignment in future edits.
- Reference documentation added under docs/: mfsed_computation_approaches.docx (mfsed approaches, Freundlich secant retardation, three-state A/B/C comparison) with the two supporting figures secant_vs_tangent.png and leaching_profiles.png; the reproducible generator scripts are in docs/scripts/ (Python: secant_vs_tangent.py, leaching_profiles.py, make_mfsed_doc.py).

## 4.6.2

**(04/2026)**

- Added non-linear Freundlich sorption isotherm option (IKD=2) as an  alternative to the existing linear isotherms (IKD=0 direct Kd; IKD=1  via Koc and %OC). The Freundlich isotherm s=Kf*C^N is implemented  consistently across all sorption-dependent processes in the pesticide  component: (a) soil leaching CDE analytical solution (cde.f): the  concentration-dependent retardation factor R=1+N*Kf*C^(N-1)*rho_b/theta_s  is linearized at the average event flux concentration C1 via a  fixed-point iteration ensuring self-consistent R and C1; (b) incoming  pesticide phase partitioning (dissolved/sorbed) uses an effective  linearized Kd_eff=Kf*C^(N-1) evaluated at the estimated inlet  concentration; (c) outflow solid/liquid phase partitioning uses  Kd_eff at the estimated outflow concentration; (d) inter-event surface  residue porewater/sorbed split is solved via Newton-Raphson iteration  on the implicit equation Cd*Vw+Kf*Cd^N*DGMsoil=DGMRES, with the  linear solution as initial guess ensuring rapid convergence; (e) all  other partitioning uses linearized Kd_eff at the relevant local  concentration. For linear isotherms (IKD=0,1) VKF=VKD and VKN=1 are  set internally so all Freundlich branches reduce exactly to the  original linear code (backwards compatible). New input format for  IKD=2 in .iwq file: "2  Kf  N" where Kf is in L^N/kg and N is the  dimensionless Freundlich exponent (0<N<=1 for favorable isotherms).  New COMMON/WQ1/ variables VKF(10) and VKN(10) added to all 5 files  that declare this block (vfsmod.f, inputs.f, cde.f, wqpest.f,  outmass.f). Reference: Muñoz-Carpena, R., S. Reichenberger, R. Sur.  2026. Effect of linear and Freundlich isotherms in the prediction with  VFSMOD of pesticide residue remobilization from vegetative filter  strips. SETAC Europe 2026, Maastrich, Netherlands.
- Minor update to outmass.f: expanded COMMON/WQ1 declaration to include  the new VKF(10) and VKN(10) arrays for Freundlich parameters.
- Updated output format for sorption coefficients in .owq and CDE  output: linear isotherm now labelled "Linear sorption coeff. (Kd)";  Freundlich now printed as single line  "Non-linear Freundlich (Kf,Nf)= Kf L^N/Kg, N (-)".
- Bug fix: added VKF(10) and VKN(10) to COMMON/WQ1 in multspcalc.f (was missing from the initial v4.6.2 update, making this file inconsistent with the other 5 files that declare this block).
- Added comments throughout cde.f, inputs.f, wqpest.f documenting the  Freundlich vs linear branch points and key equations.
- Added leaching to equilibrium after runoff event and before the next  day degradation calculations. This improves the case of tracers and  highly mobile chemicals where some mixing layer mass can be lost  quickly after the event from soil saturation to field capacity.
- Added a `.owq` warning for pesticide trapped-mass closure only when the  metric `|mf-(mfF+mfsed)|/mi` exceeds 5%. This documents the remaining  approximate coupling between event trapping (`mf`) and CDE soil-profile  mass (`mfF`) while avoiding noise from negligible roundoff or small  discrepancies in realistic low-input cases.

## 4.6.1.1

**(01/2026)**

- Added test for acceptable limits of dispersivity provided by user, 0.005<DGLD(JJ)<0.40 m (0.5-40 cm). If value exceeds limits the closest bound is taken to ensure correct program execution. Note that DGLD(JJ)=0.05 m is recommended for pesticide calculations (PRZM default).

**(10/2025)**

- Fixed mass balance error for sediment that gave mistakenly 100%.
- Fixed the case with no degradation (IDG=0) that would not run. Note that because chemical leaching in the soil and remobilization are now required in all simulations, the IWQ inputs introduced in v4.6.1 are not compatible owith previous versions. These processes need new inputs that the program must have before running. However, The minimum IWQ now needed for no degradation (IDG=0) disregard lines after IMOB. Further, some inputs (ngday=1, dghalf=100., FC=0.25, dgLD= 0.05, dgT(i)=25., dgTheta(i)=0.25) although need to be provided in IWQ for reading, are not used as they are only needed for degradation. The values given above should be used in place (these would be repated when several compounds are selected (IWQ>1 in the IKW file), i.e. {dgHalf,dgPIN,dgLD,dgmres0}_j. Note also that the output in OWQ is also modified for the no degradation case, where the parts pertaining to degradation are not shown.

## 4.6.1

**(4/2025)**

- Modifications for SPOTMOD to ensure convergence even in segments <0.5 m. For this, command line argument "2" has been added that calls a new subroutine hydfilter (in inputs.f). When selected, the program checks internally for lack of convergence in overland flow solution (Q0 vector) and if NaN is detected it reformulates and reruns the problem with a) smaller Courant (CR=CR-0.2); b) checks for unstable inflow and smooths it to a triangular shape while ensuring mass balance; c) and changes the inflow check to the first node in the small segment (SCHK=0.).
- Added also exception to handle running pesticides with no inflow (rain remobilization).
- Added new makefile for new intel ifx fortran compiler (ifort fortran disconnued).
- Added example for 3 interactive pesticides (sampleP3.prj).
- Added optimized binary for Windows (/O3) for expected faster executions.
- (In preparation) Adding Freundlich sorption isotherm option to .iww file. To ensure backwards compatibiity, this is added as an additional option to the second line of iwq (IKD=2 values read are Kd and Freundlich exponent n; IKD=3 values read are Koc, OC, and Freundlich exponent n)

## 4.6.0

**(03/2025)**

- Added additional input checks if the file exits gracefully (i.e. missing project files).
- Fixed bug in multispecies where multiple independent species  (adjacency matrix=0, except diagonal -1) did not gave same results as the species run independently one by one.

**(01/2025)**

- Added multispecies deviation and remobilization.
- Added input error checking on a number of instances, and in particular around the new requirements when running multispecies.
- Checked compatibility with single species v4.5.2

## 4.5.1

**(02/2023)**

- Added option for automatic calculation of number of nodes (N). If the user selects N=-1 in the IKW file the program will internally calculate the number of nodes for the simulation based on the filter size (VL). The user can still select N if the number of nodes is specified in that file.

## 4.5.1

**(10/2022)**

- Fixed error in dynamic d50 calculation (NPART=8) where before it used the peak rainfall intensity from file .irn and it should be the average rainfall intensity from that file.

## 4.5.0

**(11/2021)**

- Mechanistic pesticide trapping: mass balance (default in place of empirical Sabbagh). Selected with option IWQPRO=3 in IWQ file.
- Internal sediment particle d50 dynamic calculation (avoids unrealistic d50 user parametrization). It considers not only source soil characteristics but hydrology and sediment characteristics for each individual event. Selected with option NPART=8 in ISD file.
- Improved surface residue fate component after event with new: physical pesticide soil leaching,  surface (sediment and mixing layer) mass balance and degradation, partial remobilization of surface porewater mass and  carry over. Selected with option IMOB=1 at the end of IWQ file.
- Fix printing bug fix on the hydrograph in .OHY output file that did not decline after storm.

## 4.4.3

**(8/2020)**

- Updated code to handle the updated soil (.ISO) file for the shallow water table (WT) case (see manual for details). Two changes were made to minimize and simplify the input requirements for WT conditions: a) the end bottom boundary condition (for t>tw) is now hard-set to the Dupuis-Forchheimer (DF) lateral draining stream (the code still contains the Vachaud vertical draining condition that is commented-out); b) the last line in the .iso file (if present) reads the heterogeneity ratio between horizontal and vertical conductivities for the lateral DF end bottom boundary condition, and if not given is assumed to 1 (homogenous); and c) the bubbling pressure is now only a parameter of Brooks-Corey soil characteristic curves, where the van Genuchten curves are applied to the full range from saturation.

- Updated model PDF documentation (online and in distribution package) to include shallow water table (description and command line and windows users), new pesticide equations and brief description of use for regulatory pesticide exposure assessments; new references.

## 4.4.2

**(9/2019)**

- Fixed bug for z (infiltration wetting front depth) written in output table in .OHY when using the shallow water table option, whereby z continued to grow after the wetting front intersected the water table.
- Added new pesticide remobilization schemes (options in IWQ file). This includes a new partial remobilization of just the (mixing layer) liquid phase residue, where the solid (sediment-bonded) residue stays on the surface and influences de next equilibrium partitioning for the next event in the series. This new scheme is tested for highly sorbed pesticides against field data, and its effect evaluated in long-term high tier regulatory environmental assessments.

## 4.4.1

**(8/2018)**

- For compatibility with gfortran, changed line 167 in input.f from "IF(ITHETATYPE.EQ.2) then hb=1/PARW(2)" to "IF(ITHETATYPE.EQ.2) hb=1/PARW(2)"
- Fixed plotting issue with decreasing "z" (wetting front) in .ohy output file

## 4.4.0

**(8/2018)**

- Coded alternative pesticide reduction (dP) equations including Sabbagh (2009), free fitting Sabbagh from Reichenberger et al. (2018), mass balance adapt. from Muñoz-Carpena (2015), and Chen (2016). These are now considered in the first input IWQPRO of the .iwq file.
- Fixed bug on tail of hydrograph that at times would not come down to 0 when sufficient time after inflow was allowed (additional time in IRN file with rain=0 (see manual).
- Changed format to scientific(E10.3) for Vi and Ei in IWQ file. Needed to check calculation of Fph.
- Additional SCHK option, where surface flooding (ponding) is triggered when the average water level in all nodes is > 0.

## 4.3.2

**(01/2016)**

- Change in output formats in IWQ to handle a case of no inflow and precipitation.

## 4.3.1

**(11/2015)**

- Change in output format for "Sediment inflow" in OWQ to handle large numbers.

## 4.3.0

**(10/2015)**

- Promoted to version 4.3.0 release for long-term pesticide assessments.
- Increased number of iterations in the infiltration Green-Ampt solution to improve convergence.
- Change in output format for FPH in OWQ to handle large numbers.

## 4.2.6

**(10/2015)**

- Fixed a rare output overflow on the "Runoff inflow reduction" line in OWQ when there is very little field inflow and most of the outflow comes from rain.
- Increased number of iterations in the infiltration Green-Ampt solution to improve convergence. This fixes negative infiltration issue encountered by some users.
- Improved integration of sediment outflow routine. This improves mass balance and fixes an issue with very small negative sediment outflows when the amount of sediment through the filter is very large.

## 4.2.5

**(08/2015)**

- Removed extra set of parenthesis on WRITE statement of line 524 of inputs.f. GFORTRAN did not handle this in the compilation.
- Added new dynamic time step (dt) calculation scheme of Jaber & Mohtar (J. Hyd. Eng. 2002, 7(1):3-11) to handle difficult problems (kinematic shock in sharp front problems). dt is calculated from N and estimated time of concentration in the filter and not directly from CR. CR is now internally calculated from the new dt values to obtain the Petrov-Galerkin coefficients.
- Changed format of Fph output in *.OWQ to scientific format to show cases when sediment reduction is close to 100% and Fph is large.
- Fixed output when incoming pesticide is very small (close to 0) so that residue partitioning and degradation does not produce NaN in OWQ.

IMPORTANT: Recommended input values
1. Number of nodes. In the file IKW, please change the number of node calculation based on VL (filter size) to the following N=odd(VL/3*11)/2)  (Note: half N of what we used to recommend and N must be odd at the end of the calc.)
2. To improve robustness for long-simulation pesticide exposure frameworks, in the file IKW, use a CR=0.5 for all sims. This is a conservative value that is likely not needed in most single simulations (normally CR=0.8).
3. Provide long enough simulation lenght to capture tail of hydrograph. For this it is recommended the the simulation is 25% longer than the rain (or inflow). For this add to lines to the end of the IRN value wit rain equals 0 as explained in manual. For example:

```
4 4.2593E-06 'Nrain, Rpeak (m/s)
0 4.2593E-06
10800 4.2593E-06 '(3 * 3600s)
10801 0    ' (next two lines set length of simulation to 25% longer than rain)
13500 0
```

## 4.2.4

**(04/2014)**

- Added error handling for no. of days in pest. degradation calculations (max 366, file IWQ), and steps in hydrograph (file IRO) and rainfall (file IRN) (max=200).
- Added error handling when no incoming flow or sediment is provided in IRO or IRN files so that it does not produce NaN in output files.

**(02/2014)**

- Added additional cases for the water table scenarios (gasubwt.f)
- Fixed calculation of pesticide degradation rate, Kref= Ln(2)/t_halflife (inputs.f)
- Fixed unit conversion for gas constant in pesticide degradation equation (outmass.f)
- Added a check for high sorption pesticides when all pesticide in filter is sediment-bonded (outmass.f)
- Various static improvements in IWQ output file (outmass.f)

## 4.2.3

**(08/2013)**

- Fixed bug in duration of the simulation handling (IRN file). If the user did not add a double end line in the IRN (i.e. end of simulation) to the IRN file, the program now sets the end of the simulation to the longest of the field inflow (IRO file) or rain series (IRN file).

## 4.2.1

**(08/2012)**

- Fixed minor bug for case when only lateral inflow is provided (no rain), like in some model testing laboratory scenarios (UCL, Belgium). Outflow hydrograph now ends in zero if sufficient time is provided).

## 4.2.0 (Pdegr)

**(05/2012)**

- New in-filter mass pesticide mass balance calculated at the end of the event. Prepared for integration in EU SWAN registration tool.
- Degradation subroutine to calculate degradation  of sediment bonded and mixing layer residue at the end of the event towards the beginning of the next event. It uses FOCUS equations that need daily surface soil temperature and moisture for every day betwwen two consecutive runoff events.
- Daily eair temperature (i.e. PRZM files or other source) can be used as mixing and surface soil moisture. Moisture content can be estimated based on running mass balance.
- For additional details see EU AIM and SWAN reports (http://abe.ufl.edu/carpena/vfsmod/FOCUSreports.shtml).

## 4.x.x (WT)

**(06/2011)**

- Added new subroutine to solve the soil infiltration problem for unsteady rain in the presence of a shallow water table using a modified Green-Ampt infiltration model as proposed by Salvucci and Entekhabi (1995), Chu (1996) and work by the authors of this program. The method was extended  to include mass balance on the surface as proposed by Skaggs (1982) and Khaleel (1982) in Hydrologic modeling of small watersheds, ASAE mon. no. 5, and Chu (1978, Water Resour. Res.). An extended soil input file (.iso) is required in this case. If an additional numeric parameter (WTD, water table depth, m) is found in the second line of the standard .iso file, then soil characteristic curve inputs are read  (below the second line in .iso) and the new subroutine (gasubwt.f) is called. Notice that SAV and OI values are ignored. Details of the structure of the extended .iso input file are provided in the user manual, and sample files are also the distribution package (see sampleWT.prj).
- New surface ponding forcing scheme (NPFORCE=1), when overland flow reached the check node NCHK (and zw>0, t<tw for the WT case) the infiltration is at capacity regardless of ponding or not from point excess calculation.

## 3.0.P (WQ-Pest.)

**(08/2008)**

- Error checking on input files. It will print error and exit gracefully if any of the project files are missing
- Restructured code to add water quality components while ensuring backwards compatibility. New IWQ flag added at the end of IKW file inputs. If "1" is present there it will expect a new input file with IWQ extension and produce a new output OWQ file. If no flag is present or a character other that "1" is there it will execute with no water quality component. The new program files WQSUB.f handles the processing of the water quality component when the component is selected. In addition, changes were made to VFSMOD.f, INPUTS.f, FINPUTS.f and OUTMASS.f.

## v2.4.5 since version 2.4.4

**(05/2007)**

-Variables used in the sediment deposition component (YT, XT, FI) now initialized to 0. This is needed by some versions of public domain compilers that do not do this initialization by default (g77, gfortran).
-Silent option now produces outputs in og1 in g or g/s as opposed to g/cm or g/cm.s as before. This is needed in the inverse optimization when optimizing FWIDTH.

## v2.4.4 since version 2.4.3

**(01/2007)**

-Checked for initial content OI< OS) change to exclude OI=OI case that would divide by 0 during the Green-Ampt calculation.

## v2.4.3 since version 2.4.2.a

**(01/2007)**

-Added silent mode (use: vfsm filename 1). This suppresses the welcome screen for batch simulations. This is used in the inverse calibration procedure.

## v2.4.2a since version 2.4.2

**(09/2006)**

-Instantaneous infiltration fi added to *.ohy output file to facilitate water quality calculations. Modified files kwwrite.f and inputs.f.

## v2.4.2 since version 2.2.1

**(12/2004)**

-Fixed seed for Newton-Raphson in step3.f:ln 82. Found convergence problem for large runoff events in an application in Michigan.
-Changed version to match GUI

## v2.2.1 since version 2.2.0

**(05/2003)**

- Fixed infiltration solution issue when initial F is 0. It was caused by the Newton-Raphson solution algorithm used in the Mein and Larson time implicit solution to the Green-Ampt equation. The solution will not take as a seed total infiltration F= 0.
- Added work-around for problem found when compiling vfsm with the unix g77 (gcc v3.1) fortran compiler. There is a bug in the compiler by which local variables are not treated properly. This gave different (incorrect) results for the cumulative flows calculated in kwwrite.f (written in .ohy file) between code compiled with g77 and other compilers (including those in PC). The cumulative values are now declared in a common block to force them to be static.

## v2.2.0 since version 2.01

**(05/2003)**

- Added exception check for when thetai=porosity that caused the model to fail when solving the Green-Ampt ecuation.
- Added runoff mass balance in .osp file and better handling of runoff in an impermeable plane (Ks=0).
- Increased decimals in osp file to improve resolution of Sensitivity Analysis graphs.
- Moved version up to match new public release of windows design version GUI (VFSMOD-W)

## v2.01 since version 1.06

**(04/2003)**

- Modified grassed.f, step3.f and outmass.f to handle small sediment events
- Several minor changes to integrate program with new design version GUI (VFSMOD-W)
- Moved version up to match new windows design version GUI (VFSMOD-W)

## v1.06 since version 1.05

**(03/2002)**

- Fixed problem when incoming sediment was 100% fine. This was caused by a line out of sequence in sediment filtration subroutine (step3.f). Fixed.
- Changed affiliation of one of the authors.

## v1.05 since version 1.04b

**(16/12/00)**

-New way of handling a sediment filled strip (grassed.f) . After this happens (NFUP=1), all the sediment inflow is routed to the outflow end (no deposition) and without stopping the simulation. This way the simulation summary (*.osm) now shows realistic trapping efficiencies and mass balance (outmass.f).

**(5/12/00)**

- Sediment trapping subroutine (step3.f) re-worked to avoid problems with equilibrium slope Se in some cases that led to NaNs.

**(20/10/00)**

- Output format for saturated hydraulic conductivity changed to exponential.

## v1.04b since version 1.04

**(28/02/00)**

- Windows GUI to handle UH and VFSMOD programs now available. Both fortran codes modified to produce outputs suitable for the GUI.
- Improved model outputs. New output files give more details on the simulations and a new output summary file summarizes filter performance (*.osp)
- Updated and expanded documentation. The manual has been revised and now contains numerous new sections on model use, reference tables for inputs to use in the model, and new documentation on the UH utility. An on-line version of the manual can be found on the web site and a PDF version is included in the distribution files.

**(05/02/2000)**

-Modified output format (from fixed to exponential) to allow for printing of long simulations and large sediment and flow quantities in files input.f, kwwrite.f and outmass.f.

## v1.04 since version 1.03

**(10/24/99)**

- Modified finput.f to accept project files (.lis or .prj). This enables the user to use any naming convention and directory structure for inputs and outputs. Documentation needs to be added to the User guide to illustrate this method of inputting the files.

**(10/16/99)**

- Expanded Rain array to 200 to match uh. In inputs.f, fixed totalrain calculation so the array will not try to access index=i-1 at i=1,

**(08/09/99)**

- Event for which strip fills up are now handled differently. The new procedure is, after filter strip is filled up, to stop sediment deposition and short-circuit incoming sediment to the filter exit while keeping routing water. At the end of the simulation the final sediment trapping efficiency is calculated based on the sediment captured before filling the filter and the amount that passed through after that time.

**(3/28/99)**

- Preprocessing program UH now implemented. The program prepares inputs for VFSmod based on readily available hydrological and soils data (NRCS method). This program is used in design approach.

## v1.03 since version 1.02

**(03/14/99)**

- Added Total rainfall (mm) for storm in the Storm parameters section in *.ohy and *.osm output files.

**(03/14/99)**

- Reworked finput subroutine to extend filename length to 25 characters
- Documentation (both Frame and PDF) updated.

**(03/10/99)**

- File *.igr splitted in two: *.igr (vegetation characteristics) and *.isd (incoming sediment characteristics) for interfacing with "uh"

**(02/28/99)**

- New output summary file summarizing filter performance (*.osp). New info on source area (width and length) is requested from the user (file *.iro).
- Minor fix on header of file *.osm (the version date was cut short).

## v1.02 since version 1.01

**(09/14/98)**

- The user does not need to specify the finite element mesh in the *.ikw files. The information now needed is: 1) x-position for each of the segments in which the filter is divided; and 2) surface properties (n, So) for each segment. The model now chooses the nodes where the changes occur based on the total number of nodes selected by the user.

- More detail information on surface properties is included in the output.

- Some parameters little used or confusing have been eliminated (IOUT and others).

- The documentation has been enlarged and revised.
