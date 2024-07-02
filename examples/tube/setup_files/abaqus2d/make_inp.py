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
import importlib.util
import sys

coconut_path = importlib.util.find_spec('coconut').submodule_search_locations[0]
spec = importlib.util.spec_from_file_location('make_surface',
                                              coconut_path +
                                              '/coupling_components/solver_wrappers/abaqus/extra/make_surface.py')
module = importlib.util.module_from_spec(spec)
sys.modules['make_surface'] = module
spec.loader.exec_module(module)
from make_surface import *

mdb = Mdb(pathName='case_tube2d.cae')
tubeModel = mdb.ModelFromInputFile(name='Model-1', inputFileName='mesh_tube2d.inp')
tubeMaterial = tubeModel.Material(name='Material')
tubeMaterial.Elastic(table=((300000.0, 0.3),))
tubeMaterial.Density(table=((1200.0,),))
tubeAssembly = tubeModel.rootAssembly
tubeInstance = tubeAssembly.instances['PART-1-1']
tubePart = tubeModel.parts['PART-1']
tubePart.setValues(space=AXISYMMETRIC, type=DEFORMABLE_BODY)
tubeModel.HomogeneousSolidSection(material='Material', name='TubeSection', thickness=1.0)
tubePart.SectionAssignment(offset=0.0, region=Region(elements=tubePart.elements), sectionName='TubeSection')
step1 = tubeModel.ImplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=0.0001, nlgeom=ON, maxNumInc=1,
                                       haftol=1, initialInc=0.0001, minInc=0.0001, maxInc=0.0001, amplitude=RAMP,
                                       noStop=OFF, nohaf=ON, initialConditions=OFF, timeIncrementationMethod=FIXED,
                                       application=QUASI_STATIC)
step1.Restart(frequency=99999, overlay=ON)
beamInsideMoving = SurfaceFromNodeSet(tubeAssembly, tubeInstance, 'BEAMINSIDEMOVING', 'BEAMINSIDEMOVING')
tubeModel.Pressure(name='DistributedPressure', createStepName='Step-1', distributionType=USER_DEFINED, field='',
                   magnitude=1, region=beamInsideMoving)
tubeModel.SurfaceTraction(name='DistributedShear', createStepName='Step-1', region=beamInsideMoving, magnitude=1,
                          traction=GENERAL, directionVector=((0, 0, 0), (1, 0, 0)), distributionType=USER_DEFINED)
tubeModel.DisplacementBC(name='FixedEnds', createStepName='Step-1', region=tubeAssembly.sets['BEAMINSIDEFIXED'], u1=0,
                         u2=0, ur3=UNSET)
tubeModel.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-1',
                             region=tubeAssembly.sets['BEAMINSIDEMOVING'], variables=('COORD', 'U'))
tubeModel.FieldOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='F-Output-2', variables=PRESELECT)
tubeModel.HistoryOutputRequest(createStepName='Step-1', frequency=LAST_INCREMENT, name='H-Output-1',
                               variables=PRESELECT)
jobName = 'case_tube2d'
tubeJob = mdb.Job(name=jobName, model='Model-1', description='tube2d')
tubeJob.writeInput()
