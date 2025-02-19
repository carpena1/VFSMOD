# VFSMOD

Vegetative Filter Strip Modeling System (VFSMOD) is a design-oriented vegetative filter strip modeling system. The MS-Windows VFSMOD-W graphical user interface (GUI) integrates the numerical model VFSMOD, a utility to generate source (upslope disturbed area) inputs for the model based on readily available NRCS site characteristics (UH), and advanced uncertainty and sensitivity analysis, inverse calibration and design menu-driven components.

Vegetative Filter Strip Modeling System (VFSMOD-W)

VFSMOD, the core of the modeling system, is a computer simulation model created to study hydrology, sediment and pollutant transport through vegetative filter strips (VFS). The model comprises the following modules:

- Finite element solution for the overland flow equations.
- Green-Ampt infiltration method for unsteady rainfall.
- University of Kentucky grass sediment deposition and filtration.
- Contaminant transport component (first implemented for pesticide reduction within the risk massessment and management regulatory process).

INSTALLATION INSTRUCTIONS FOR VFSmod v4.x.x
(See latest version changes in CHANGES.txt files in source code directories src_vfsm and scr_uh)

Obtaining VFSMOD
----------------------------------------------------------------
Clone this repository. You can find additional release information at https://abe.ufl.edu/carpena/vfsmod/

Installing and running VFSMOD
----------------------------------------------------------------
VFSMOD (and UH) source code is distributed both in DOS and UNIX versions along with make files and sample input and output files. The source code is written entirely in standard FORTRAN77 so that compilation should be straight forward following the included makefile and using the proper set of files for each platform (DOS/UNIX). Binaries for a few computer platforms can also be found at the internet sites. 


Installing in a DOS system (under MS-Windows '95 or later)
-----------------------------------------------------------------

a) Create a directory named VFSMOD

b) Expand the contents of the file vfsmodpc.zip. This should create the following directory structure
				VFSMOD
	SRC_VFSM	SRC_UH	DOCS	INPUTS		OUTPUT

c) The executable files VFSM.EXE (and UH.EXE) can be found in the parent directory VFSMOD.

d) Run the sample case named SAMPLE, by typing at the DOS prompt:
C:/VFSMOD-W> VFSM SAMPLE
or if you want "silent" execution (no banner) just type:
C:/VFSMOD-W> VFSM SAMPLE 1


Please note that the second part of the command issued (SAMPLE) refers to a set of files located in the subdirectory INPUTS. You could run a different problem by selecting a different set of input files with the condition that they are located in the subdirectory INPUTS. In this example, if you issue the DIR command within the INPUTS directory you should see the following files:

SAMPLE.IGR  SAMPLE.IKW  SAMPLE.IRN  SAMPLE.IRO  SAMPLE.ISD  SAMPLE.ISO

During the run a set of new files is created in the OUTPUT directory:

SAMPLE.OG1  SAMPLE.OG2  SAMPLE.OHY  SAMPLE.OSM  SAMPLE.OSP

The content of both input/output files is explained in detail in User's Manual.



UNIX system
-----------------

a) Create a directory named VFSMOD

	mkdir VFSMOD
	mv vfsm__ux.tar.gz VFSMOD
	cd VFSMOD

b) Expand the contents of the file vfsm20ux.tar.gz on the new directory. 

	gzcat vfsm__ux.tar.gz | tar xvf -

This should create the following directory structure
				VFSMOD
	src_vfsm	src_uh	docs	inputs		output

c) An installation script (setup) is included in the VFSMOD directory. To compile and install the program simply type setup. The script will compile the source code and copy the executable files (vfsm and uh) to the VFSMOD directory. If your FORTRAN compiler name is not f77 you will need to edit the makefile found in the src directories. You can also clean the executable and object files by typing setup clean.

d) Run the sample case named sample by typing at the UNIX prompt:
> vfsm sample
or if you want "silent" execution (no banner) just type:
> vfsm sample 1

Please note that the second part of the command issued (sample) refers to a set of files located in the subdirectory inputs. You could run a different problem by selecting a different set of input files with the condition that they are located in the subdirectory inputs. Note that you must have all the six input files in order to run the program. In our example, if you issue the ls command within the inputs directory you should see the following files:

sample.igr  sample.ikw  sample.irn  sample.iro  sample.isd  sample.iso

After you execute the command you should see a screen similar to the one given above. During the run a new set of files is created in the output directory:

sample.og1  sample.og2  sample.ohy  sample.osm  sample.osp

The content of both input/output files is explained in detail 
in the User's Manual.
