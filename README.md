
# ISOTOPIA
ISOTOPIA is a software package for the prediction of radio-isotope production. It uses cross sections from the TENDL nuclear data library and calculates the radioactive yield as a function of irradiation time and other irradiation characteristics.

## Documentation and reference
A description of the code and its options can be found in the [ISOTOPIA User Manual (pdf)](https://github.com/arjankoning1/isotopia/blob/main/doc/isotopia.pdf).
The reference to be used for ISOTOPIA is

A.J. Koning, D. Rochman, J.-Ch. Sublet, N. Dzysiuk, M. Fleming, and S. van der Marck, *TENDL: Complete Nuclear Data Library for innovative Nuclear Science and Technology*, Nuclear Data Sheets 155,1 (2019).

## Installation

### Prerequisites:

The following are the prerequisites for compiling ISOTOPIA:
  - git (only if the package is downloaded via Github)
  - a recent Fortran compiler such as gcc (gfortran)

### Downloads:

To download ISOTOPIA, you can use one of the following options:
#### 1. Download the entire tar file:
```
https://nds.iaea.org/talys/isotopia.tar
tar zxf isotopia.tar
```

#### 2. Using git:
```
git clone https://github.com/arjankoning1/isotopia.git
```
The ISOTOPIA sample cases do not fall under the git repository. For that you need to download:
```
https://nds.iaea.org/talys/isotopiasamples.tar
tar zxf isotopiasamples.tar
```
The resulting *samples/* directory should be moved to the *isotopia/* directory.

### Installation instructions:

To install ISOTOPIA, you can use one of the following options:
#### 1. Using make:
```
cd isotopia/source
make
```
#### 2. Using the code_build script:
```
cd isotopia
code_build isotopia
```

The above will produce a *isotopia* executable in the *isotopia/bin* directory.
The compiler and its flags can be set in either the *source/Makefile* or in *code_build*.

## The ISOTOPIA package

The *isotopia/* directory contains the following directories and files:

+ `README.md` is this README file
+ `LICENSE` is the License file
+ `code_build` and `path_change` installation scripts
+ `source/` contains the Fortran source code of ISOTOPIA and the Makefile
+ `bin/` contains the executable after successful installation
+ `files/` contains abundace and decay data files needed for the calculations
+ `doc/` contains the tutorial in pdf format
+ `samples/` contains the input and output files of the sample cases, and the *verify* script for the user to run the sample cases

In total, you will need about 50 Mb of free disk space to install ISOTOPIA.

## Sample cases

A successful installation can be verified by running the sample cases. For each sample case, the results are written to a subdirectory *new/*, which can then be compared with the output provided in the ISOTOPIA package in the subdirectory *org/*. The entire sample set will take about 1 minute.
```
cd samples
./verify
```

ISOTOPIA works as follows:
```
isotopia < isotopia.inp > isotopia.out
```
assuming that *isotopia* is linked to the *isotopia/bin/isotopia* executable.

## License and Copyright
This software is distributed and copyrighted according to the [LICENSE](LICENSE) file.
