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

e_modulus = 1.4e6
poisson_coefficient = 0.4
density = 1000
gravity = 0

mdb = Mdb(pathName='case_turek.cae')
beamModel = mdb.ModelFromInputFile(name='Model-1', inputFileName='mesh_turek.inp')
beamMaterial = beamModel.Material(name='Material')
beamMaterial.Elastic(table=((e_modulus, poisson_coefficient),))
beamMaterial.Density(table=((density,),))
beamAssembly = beamModel.rootAssembly
beamInstance = beamAssembly.instances['PART-1-1']
beamPart = beamModel.parts['PART-1']
beamModel.HomogeneousSolidSection(material='Material', name='BeamSection', thickness=1.0)
beamPart.SectionAssignment(offset=0.0, region=Region(elements=beamPart.elements), sectionName='BeamSection')
step1 = beamModel.StaticStep(name='Step-1', previous='Initial', timePeriod=1.0, nlgeom=ON, initialInc=0.01,
                             minInc=0.0001, maxNumInc=1000, amplitude=RAMP)
step1.Restart(frequency=99999, overlay=ON)
beamInsideMoving0 = SurfaceFromNodeSet(beamAssembly, beamInstance, 'BEAMINSIDEMOVING0', 'BEAMINSIDEMOVING0')
beamInsideMoving1 = SurfaceFromNodeSet(beamAssembly, beamInstance, 'BEAMINSIDEMOVING1', 'BEAMINSIDEMOVING1')
beamInsideMoving2 = SurfaceFromNodeSet(beamAssembly, beamInstance, 'BEAMINSIDEMOVING2', 'BEAMINSIDEMOVING2')
beamModel.Pressure(name='DistributedPressure0', createStepName='Step-1', distributionType=USER_DEFINED, field='',
                   magnitude=1.0, region=beamInsideMoving0)
beamModel.Pressure(name='DistributedPressure1', createStepName='Step-1', distributionType=USER_DEFINED, field='',
                   magnitude=1.0, region=beamInsideMoving1)
beamModel.Pressure(name='DistributedPressure2', createStepName='Step-1', distributionType=USER_DEFINED, field='',
                   magnitude=1.0, region=beamInsideMoving2)
beamModel.SurfaceTraction(name='DistributedShear0', createStepName='Step-1', region=beamInsideMoving0, magnitude=1,
                          traction=GENERAL, directionVector=((0, 0, 0), (1, 0, 0)), distributionType=USER_DEFINED,
                          follower=OFF)
beamModel.SurfaceTraction(name='DistributedShear1', createStepName='Step-1', region=beamInsideMoving1, magnitude=1,
                          traction=GENERAL, directionVector=((0, 0, 0), (1, 0, 0)), distributionType=USER_DEFINED,
                          follower=OFF)
beamModel.SurfaceTraction(name='DistributedShear2', createStepName='Step-1', region=beamInsideMoving2, magnitude=1,
                          traction=GENERAL, directionVector=((0, 0, 0), (1, 0, 0)), distributionType=USER_DEFINED,
                          follower=OFF)
if gravity > 0.0:
    beamModel.Gravity(name='Gravity', createStepName='Step-1', comp2=-gravity)
beamModel.DisplacementBC(name='FixedLeftEnd', createStepName='Step-1', region=beamAssembly.sets['BEAMINSIDEFIXED'],
                         u1=0, u2=0, ur3=0)
beamModel.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-1',
                             variables=('U', 'COORD'))
beamModel.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-2', variables=PRESELECT)
beamModel.HistoryOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='H-Output-1',
                               variables=PRESELECT)
jobName = 'case_turek'
beamJob = mdb.Job(name=jobName, model='Model-1', description='turek')
beamJob.writeInput()
