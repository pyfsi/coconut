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
import os
import shutil



mdb = Mdb(pathName='Test.cae')
yarnModel = mdb.ModelFromInputFile(name='Model-1',inputFileName='Yarn_Base.inp')
#yarnMaterial = yarnModel.Material(name = 'Material');
#yarnMaterial.Elastic(table= ((2.5e+09, 0.39), ))
#yarnMaterial.Density(table=((1140.0, ), ))
yarnAssembly = yarnModel.rootAssembly
yarnPart = yarnModel.parts['YARN']
yarnInstance = yarnAssembly.instances['INSTANCE-YARN']

yarnModel.CircularProfile(name='Profile-1', r=0.00036)
shear_modulus = 2.5e+09/(2*(1+0.39))
yarnModel.BeamSection(alphaDamping=0.0, beamShape=CONSTANT,
    betaDamping=0.0, centroid=(0.0, 0.0), compositeDamping=0.0,
    consistentMassMatrix=False, density=1140.0, dependencies=0, integration=
    BEFORE_ANALYSIS, name='Section-1', poissonRatio=0.39, profile='Profile-1',
    shearCenter=(0.0, 0.0), table=((2.5e+09, shear_modulus), ),
    temperatureDependency=OFF, thermalExpansion=OFF)
yarnPart.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=yarnPart.sets['ALL_ELEMENTS'], 
    sectionName = 'Section-1', thicknessAssignment=FROM_SECTION)

yarnAssembly.regenerate()


#
if (1.0> 0.0):
	rigidPart = mdb.models['Model-1'].Part(name='PART-2', dimensionality=THREE_D, type=ANALYTIC_RIGID_SURFACE)
	s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',sheetSize=0.5)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=STANDALONE)
	s.ConstructionLine(point1=(0.0, -0.25), point2=(0.0, 0.25))
	s.FixedConstraint(entity=g[2])
	s.Spot(point=(0.041-0.00036,0.0048))
	s.Spot(point=(0.001-0.00036,0.0048))
	s.Spot(point=(0.041-0.00036,-0.0352))
	s.Spot(point=(0.001-0.00036,0.0352))
	s.Spot(point=(0.00124-0.00036,0.04045))
	s.Spot(point=(0.00166-0.00036,0.283102382))
	s.ArcByCenterEnds(center=(0.041-0.00036,0.0048), direction=CLOCKWISE, point1=(0.041-0.00036,-0.0352), point2=(0.001-0.00036,0.0048))
	s.Line(point1=(0.001-0.00036,0.0048), point2=(0.001-0.00036,0.0352))
	s.Line(point1=(0.001-0.00036,0.0352), point2=(0.00124-0.00036,0.04045))
	s.Line(point1=(0.00124-0.00036,0.04045), point2=(0.00166-0.00036,0.283102382))
	rigidPart.AnalyticRigidSurfRevolve(sketch=s)
	s.unsetPrimaryObject()
	rigidPart = mdb.models['Model-1'].parts['PART-2']
	del mdb.models['Model-1'].sketches['__profile__']
	rigidPart.ReferencePoint(point=(0.0, 0.0, 0.0))
	yarnAssembly.Instance(name='PART-2-1', part=rigidPart, dependent=OFF)
	yarnAssembly.rotate(instanceList=('PART-2-1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle=270.0)
	rigidSurface = rigidPart.faces
	side1Faces = rigidSurface.getSequenceFromMask(mask=('[#7 ]', ), )
	rigidPart.Surface(side1Faces=side1Faces, name='RIGID')
	yarnAssembly.regenerate()
#

#

#
if (1.0>0.0):
	r1 = yarnAssembly.instances['PART-2-1'].referencePoints
	refPoints1=(r1[2], )
	region = yarnAssembly.Set(referencePoints=refPoints1, name='Set-8')
#

#
	yarnPart.Surface(circumElements = yarnPart.sets['ALL_ELEMENTS'].elements, name='Surf-1')

	yarnModel.ContactProperty('IntProp-1')
	yarnModel.interactionProperties['IntProp-1'].TangentialBehavior(dependencies=0, directionality=ISOTROPIC, elasticSlipStiffness=None,
    formulation=PENALTY, fraction=0.005, maximumElasticSlip=FRACTION, pressureDependency=OFF, shearStressLimit=None, slipRateDependency=OFF,
    table=((0.0, ), ), temperatureDependency=OFF)
	yarnModel.interactionProperties['IntProp-1'].NormalBehavior(allowSeparation=ON, constraintEnforcementMethod=DEFAULT, pressureOverclosure=HARD)
	yarnModel.interactionProperties['IntProp-1'].GeometricProperties(contactArea=1.0, padThickness=None)
	yarnModel.rootAssembly.regenerate()

#


#
#

#Need to assign a beam direction to be able to write a .inp file, this has to be replaced later 
yarnPart.assignBeamSectionOrientation(method= N1_COSINES, n1=(0.0, 0.0, -1.0), region=yarnPart.sets['ALL_ELEMENTS'])
##

jobName = 'Time0_Temp'
yarnJob = mdb.Job(name = jobName, model = 'Model-1', description = 'Yarn')
yarnJob.writeInput()
