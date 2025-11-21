from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

import sys
if sys.version_info[0] == 2:
    # Abaqus 2023 and older work with Python 2.7; this version can't work with the symbols in the CoCoNuT logo
    # Therefore, make_surface has to be imported as a stand-alone module
    import imp

    coconut_path = imp.find_module('coconut')[1]
    fp, path_name, description = imp.find_module('make_surface',
                                                 [coconut_path + '/coupling_components/solver_wrappers/abaqus/extra/'])
    imp.load_module('make_surface', fp, path_name, description)
    from make_surface import *
else:
    # Abaqus 2024 and newer work with Python 3.10; a normal import is possible
    from coconut.coupling_components.solver_wrappers.abaqus.extra.make_surface import *

mdb = Mdb(pathName='case_tube3d.cae')
tubeModel = mdb.ModelFromInputFile(name='Model-1', inputFileName='mesh_tube3d.inp')
tubeMaterial = tubeModel.Material(name='Material')
tubeMaterial.Elastic(table=((300000.0, 0.3),))
tubeMaterial.Density(table=((1200.0,),))
tubeAssembly = tubeModel.rootAssembly
tubeInstance = tubeAssembly.instances['PART-1-1']
tubePart = tubeModel.parts['PART-1']
tubeModel.HomogeneousShellSection(integrationRule=SIMPSON, material='Material', name='ShellSection',
                                  nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, preIntegrate=OFF,
                                  temperature=GRADIENT, thickness=0.001, thicknessModulus=None, useDensity=OFF)
tubePart.SectionAssignment(offset=-0.5, region=Region(elements=tubePart.elements), sectionName='ShellSection')
tubeAssembly.regenerate()
tubeModel.DisplacementBC(amplitude=UNSET, createStepName='Initial', distributionType=UNIFORM, fieldName='',
                         localCsys=None, name='BC-1', region=tubeAssembly.sets['LEFTFIXED'], u1=SET, u2=SET, u3=SET,
                         ur1=UNSET, ur2=UNSET, ur3=UNSET)
tubeModel.DisplacementBC(amplitude=UNSET, createStepName='Initial', distributionType=UNIFORM, fieldName='',
                         localCsys=None, name='BC-2', region=tubeAssembly.sets['RIGHTFIXED'], u1=SET, u2=SET, u3=SET,
                         ur1=UNSET, ur2=UNSET, ur3=UNSET)
wallOutside = tubeAssembly.Surface(name='WALLOUTSIDE', side2Elements=tubeInstance.elements)
tubeAssembly.regenerate()
step1 = tubeModel.ImplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=0.01, nlgeom=ON, maxNumInc=100,
                                       haftol=1, initialInc=0.0001, minInc=0.0001, maxInc=0.01, amplitude=RAMP,
                                       noStop=OFF, nohaf=ON, initialConditions=OFF, timeIncrementationMethod=FIXED,
                                       application=QUASI_STATIC)
jobName = 'case_tube3d'
tubeJob = mdb.Job(name=jobName, model='Model-1', description='tube3d')
tubeJob.writeInput()
