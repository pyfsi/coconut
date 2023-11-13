### source_code_structure_wire_drawing ###
OpenFOAM code repository for wire drawing and wire rolling simulation  
  

### Important information ###

This code for structural calculation is a reduced version of the original one. The solver can be find in the coconut_plasticNonLinSolidFoam-directory 
The most important reduction is the absence to apply the layering technique which is described in the artikel. 
The layering technique is developed inside the Bekaert N.V. company and is classified. To let an FSI simulation run, 
use a wire which is longer without applying layering.
To get the complete code, please contact philip.cardiff@ucd.ie  
  

### Installation ###

First install and compile foam-extend-4.1/foam-extend on your system.  
Then compile source_code_structure_wire_drawing using the included `Allwmake` script.  
  

### Systems and Compilers ###

The compilaltion has been checked on with systems and compilers, for example:  

    - macOS 10.12.6 : gcc (GCC) 4.9.2 20141029 (prerelease)  
    - Ubuntu 16.04.3 LTS : gcc (Ubuntu 4.9.3-13ubuntu2) 4.9.3  
    - Ubuntu 16.04.3 LTS : icc (ICC) 18.0.2 20180210  
    - Red Hat Enterprise Linux Server release 6.9 (Santiago) : gcc (GCC) 4.9.2  
    - SUSE Linux Enterprise Server 11 (x86_64) : icc (ICC) 15.0.3 20150407
    - CO7: GCC/8.3.0   


