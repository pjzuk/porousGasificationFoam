#  porousGasificationMedia README
***

This folder contains files and programs created under
GNU GPL v3 License

by [Filip Mróz](https://github.com/Fafa87) & [Pawel Jan Zuk](https://github.com/pjzuk) 2025


if you use any part of this work please cite the scientific contribution:


OpenFOAM solver for thermal and chemical conversion in porous media
Pawel Jan Zuk, Bartosz Tużnik, Tadeusz Rymarz, Kamil Kwiatkowski,
Marek Dudyński, Flavio C. C. Galeazzo, Guenther C. Krieger Filho
Submitted to Computer Physics Communications

This is a OpenFOAM-v2406 build which can be also obtained form
https://github.com/pjzuk/porousGasificationFoam

The authors provide the equivalent (on the day of release)
implementation in OpenFOAM 8,
which can be obtained from repository:
https://github.com/btuznik/porousGasificationFoam


Table of content:
1. [Installation guide](#installation)
2. [Source guide](#source)
3. [Documentation](#doc)

<a name="installation"></a>
## Installation guide

The installation guide is prepared under the assumption 
that OpenFOAM is installed in standard location: `/opt/OpenFoam-8/`

To install the 'porousGasificationFoam':
1. Set the OpenFOAM environmental paths 

2. Check the environmental settings
    * run any OpenFOAM solver, e.g. icoFoam -help
    * type: `$ echo $WM_PROJECT_USER_DIR`

3. Optionally change the destination path. The default path is: `$WM_PROJECT_USER_DIR/`.
   To change destination path edit file porousGasificationMediaDirectories located
   in the library installation folder.

4. Set the package environment by typing from inside of porousGasificationFoam folder:

    `$ source porousGasificationMediaDirectories`

5. Run the install script:

   `$ . ./Allwmake`

6. Test the installation by running the solver:

    `$ porousGasificationFoam`

In case of errors:

1. Make sure the OpenFoam and all required packaches have been correctly installed,
   by testing one of the tutorial cases provided with OpenFOAM. 
   
2. Open new terminal and again set the nessesary paths:

    `$ source source source /<PATH TO OpenFOAM-v2406 installation folder>/OpenFOAM-v2406/etc/bashrc`\
    `$ source <PATH TO porousGasificationMediaDirectories>/porousGasificationMediaDirectiories`

3. Check the setting:

    run: `icoFoam -help` (output should be 'usage info')\
    `echo $WM_PROJECT_USER_DIR `\
    `echo $FOAM_HGS` \
    
    If missing set these variables manually.
    
4. Run instalaltion script and check for errors:

    `./Allwmake > log_install &`

<a name="source"></a>

## Source guide

###  Installation part

Files for installation and sourcing paths:

* `./README.md` -- readme file

* `./porousGasificationMediaDirectories` -- file with enviromental variables
                                            needed for the installation

* `./Allwmake` -- installation script

###  porousGasificationFoam -- solver

 porousGasificationFoam's main code, that uses porousGasificationMedia library.
 All calculations are scheduled here.

### porousGasificationMedia -- library

porousGasificationMedia library inculding four major parts:
1. pyrolysisModels -- classes that evaluate porous medium state and properties.
2. thermophysicalModels -- implementation of thermophysical and chemical properties od porous medium.
3. fieldPorosityModel -- implementation of mechanical properties of porous medium.
4. radiationModels -- heterogeneous radiation model

###  Utilities

1. setPorosity -- utility for creating porosity fields:
    *  porosityF -- porosity field
    *  Df -- Darcy porous resitance tensor
2. totalMassPorousGasificationFoam -- Integrate solid state mass over the whole
   computational domain. The calculation is performed for each stored time step.

<a name="doc"></a>
# Documentation

The documentation of the model can be generated with Doxygen. To build
the documentation the doxygen and the graphviz packages are required.
For Ubuntu users the packages can be obtained with the following command:

`$ sudo apt-get install doxygen graphviz`

To generate the documentation go to the `$POROUS_DOC_SRC/doc/Doxygen` directory
and type:

`$ ./Allwmake`

In order to view the documentation open the 
`$WM_PROJECT_DIR/doc/Doxygen/html/index.html` in a desired internet browser.
