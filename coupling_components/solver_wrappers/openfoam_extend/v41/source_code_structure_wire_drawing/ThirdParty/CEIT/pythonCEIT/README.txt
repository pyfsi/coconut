This is a guide on how to use the FFT solver and its associated methods

Requirements:
There are some Python package requirements for the FFT solver
'numpy', 'scipy'
On some systems the Intel mkl system libraries can cause issues, if
errors appear when trying to run the solver, the additional package
'nomkl' should fix them

The ceit model 'mesostressAtIntPoint_v2.f' has been slightly modified
to allow it to be used with the Python utility f2py to compile it as
a Python module.

Steps:
1. Generate strain paths

The RVE requires as input the deformation gradient time series which it will
simulate. These are in the standard OpenFOAM tensor interpolation table format

Paths can be written by hand or extracted from previously run cases using the
post-processing utility 'extractDeformationGradients' (only for wire process lines).
This is the quickest way to get deformation paths but only at the written time steps

Alternatively there exists a function object 'deformationPath' which will generate
the deformation paths for each cell in a simulation

For wire process lines, the deformation paths at each pass (drawing1 etc)
need to be lined up and connected between each pass. A utility
'forward_append_process_line.py' exists for this. As the 'extractDeformationGradients'
utility will often generate deformation paths with too few data points (as it reads
written time-steps), a utility 'smoothPaths.ipynb' exists to interpolate the deformation
paths  between points while maintaining the deformation volume (J = det(F)). This
is a Jupyter notebook file

2. Run RVE solver
The 'ceit_fft.py' can then be run on the deformation paths to create an output
of stress vs deformation gradient. Some basic microstructural data is also recorded
such as mean lamellar spacing and the mean equivalent plastic strain in the ferrite
and cementite phases. Future work will include more microstructural data being considered
for the RVE.
