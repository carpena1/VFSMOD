# VFSMOD - Vegetative Filter Strip Modeling System
Written by Rafael Muñoz-Carpena, email: carpena@ufl.edu

VFSMOD (Vegetative Filter Strip Modeling System) is a computer simulation model created to study hydrology, sediment, and pollutant transport through vegetative filter strips (VFS). VFS are commonly used as a best management practice for runoff pollution mitigation from disturbed areas (agricultural, forest, urban, transportation lanes). VFSMOD, the core modeling engine, is a stand-alone command-line program written in Fortran with binaries distributed for different platforms. The stand-alone configuration allows for integration into other frameworks, including pesticide environmental risk assessment (ERA) tools in North America (VFS4PWC, PREM, VFSpipe, etc.) and Europe (FOCUS SWAN, pyVFSMOD). The model comprises the following modules:
- An advanced Petrov-Galerkin finite element solution for the overland flow equations.
- Green-Ampt infiltration method for unsteady rainfall.
- University of Kentucky grass sediment deposition and filtration.
- Contaminant transport component (first implemented for pesticide reduction within the risk assessment and management regulatory process).

Full details of the model can be found at https://abe.ufl.edu/carpena/vfsmod/

## Free Licensing 

VFSMOD by (c) Rafael Muñoz-Carpena is licensed under CC BY-ND 4.0 

The model is provided to you as an educational, research, and general application tool under the terms of the Creative Commons license, CC BY-ND 4.0 (Creative Commons Attribution-NoDerivatives 4.0 International). This license requires that reusers give credit to the creator. It allows reusers to copy and distribute the material in any medium or format in unadapted form only, even for commercial purposes.

## Changes

See latest version changes in CHANGES.txt files in the source code directories /src_vfsm and /scr_uh

## VFSMOD-GUI 

A VFSMOD platform-independent graphical user interface (GUI) was developed and maintained by Dr. Iñigo Berberana (inigo.barberena@unavarra.es) at the Public University of NAvarra (Spain), for his Ph.D. dissertation under the supervision of Drs. Muñoz-Carpena, Miguel Campo-Bescós and Javier Casalí. The GUI provides a complete modeling framework integrating:
- The numerical model VFSMOD, 
- A utility to generate source (upslope disturbed area) inputs for the model based on readily available NRCS site characteristics (UH)
- Advanced uncertainty and sensitivity analysis
- Inverse calibration
- Design under uncertainty
For portability (Windows, Linux, macOS) the GUI is entirely written in Python and packaged in a stand-alone installer (.exe for Windows, and .pkg for macOS that includes Intel and Apple silicon binaries). 

The GUI can be downloaded from Dr. Barberena's GitHub: https://github.com/Inigobarbe/VFSMOD-GUI. 

## Installing and running VFSMOD from the Command Line

Clone this repository for source code and example project files for running the command-line program. VFSMOD (and UH) source code is distributed both in Windows and Unix (macOS and Linux) versions, with source code, make files for GNU GFortran and Intel's ifort/ifx compilers, and sample input and output files. The source code is written entirely in FORTRAN so that compilation should be straightforward following the included makefile and using the proper set of files for each platform (Windows/UNIX). 

### Installing in Windows systems 

a) Create a directory named VFSMOD.

b) Expand the contents of the .zip file. This should create the following directory structure:

				      VFSMOD

	SRC_VFSM	SRC_UH	DOCS	INPUTS		OUTPUT

c) The executable files VFSM.EXE (and UH.EXE) can be found in the parent directory VFSMOD.

d) Run the sample case named SAMPLE, by typing at the DOS prompt:
> VFSM SAMPLE (or VFSM SAMPLE.PRJ)

or if you want "silent" execution (no banner) just type:

> VFSM SAMPLE 1 or VFSM SAMPLE.PRJ

If you just type _VFSM_ in the command line the program will provide a brief description of use.

Please note that the second part of the command issued (SAMPLE) refers to a set of files with the same name (with different extensions) located in the subdirectory INPUTS. You could run a different problem by selecting a different set of input files under the project file .PRJ with the condition that corresponding input files are located in the subdirectory INPUTS. See examples in the distribution directory and documentation online. In this example, if you issue the DIR command within the INPUTS directory you should see the following files:

_SAMPLE.IKW,  SAMPLE.IGR,  SAMPLE.IRN,  SAMPLE.IRO,  SAMPLE.ISD,  SAMPLE.ISO_

or for the case of a pesticide (see SAMPLEP.PRJ):

_SAMPLE.IKW,  SAMPLE.IGR,  SAMPLE.IRN,  SAMPLE.IRO,  SAMPLE.ISD,  SAMPLE.ISO, SAMPLEP.IWQ_

During the run, a set of new files is created in the OUTPUT directory:

_SAMPLE.OG1,  SAMPLE.OG2,  SAMPLE.OHY,  SAMPLE.OSM,  SAMPLE.OSP (SAMPLEP.OWQ)_

The content of both input/output files is explained in detail in the User's Manual online.


### Installing in UNIX systems

a) Create a directory named VFSMOD

	mkdir VFSMOD
	mv vfsm__ux.tar.gz VFSMOD
	cd VFSMOD

b) Expand the contents of the file vfsm20ux.tar.gz in the new directory. 

	gzcat vfsm__ux.tar.gz | tar xvf -

This should create the following directory structure:

				      VFSMOD

	src_vfsm	src_uh	docs	inputs		output

c) An installation script (setup) is included in the VFSMOD directory. To compile and install the program, simply type setup. The script will compile the source code and copy the executable files (vfsm and uh) to the VFSMOD directory. If your FORTRAN compiler name is not f77, you will need to edit the makefile found in the src directories. You can also clean the executable and object files by typing setup clean.

d) Run the sample case named sample by typing at the UNIX prompt:

> vfsm sample (or vfsm sample.prj)

or if you want "silent" execution (no banner), just type:

> vfsm sample 1 (or vfsm sample.prj 1)

If you just type _vfsm_ in the command line the program will provide a brief description of use.

Please note that the second part of the command issued (sample) refers to a set of files with the same name (with different extensions)located in the subdirectory inputs. You could run a different problem by selecting a different set of input files under the project file .prj with the condition that corresponding input files are located in the subdirectory INPUTS. See examples in the distribution directory and documentation online. In our example, if you issue the ls command within the inputs directory, you should see the following files:

_sample.igr,  sample.ikw,  sample.irn,  sample.iro,  sample.isd,  sample.iso_

0r for the case of with a pesticide (see sampleP.prj):

_sampleP.ikw,  sample.igr,  sample.irn,  sample.iro,  sample.isd,  sample.iso, sampleP.i wq_

After you execute the command, you should see a screen similar to the one given above. During the run, a new set of files is created in the output directory:

_sample.og1,  sample.og2,  sample.ohy,  sample.osm,  sample.osp (sampleP.owq)_

The content of both input/output files is explained in detail 
in the User's Manual online.
