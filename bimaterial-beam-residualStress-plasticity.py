#coding=utf-8

# This script performs a K_I analysis of a bimaterial beam 
# under three-point bending and thermal loading using Abaqus.
# Plasticity and residual stresses due to thermal expansion mismatch are considered.
# The script is structured into the following main parts:
# 1. Definition of parameters
# 2. Functions for model creation
#    - create_parts
#    - define_materials
#    - assemble_beam
#    - define_steps
#    - apply_bcs_and_loads
#    - define_crack
#    - mesh_assembly
# 3. Main execution block that calls the functions
#

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from abaqusConstants import Q_VECTORS
from abaqusConstants import SINGLE_NODE

from odbAccess import openOdb
from regionToolset import Region

import pickle
import json

# 1. PARAMETERS

#Script Settings
odbBase = 'C:/ABAQUS/' #Base directory in which the odb files are saved
disablePlasticity = False #If true, plasticity is disabled for both materials - useful for debugging

# Job Parameters
job_name = 'PlaneStrain'
num_cpus = 6
numDomains = 6

# Step Settings
initialIncrement = 0.001 # Initial increment size for the analysis steps
minIncrement = 1e-12 # Minimum increment size for the analysis steps
maxIncrement = 0.025 # Maximum increment size for the analysis steps
recordSteps = 0.1 # Step size for recording field outputs


# Geometric Parameters
beam_length = 40.0  # Total length of the beam
wre_height = 2.0    # Height of the WRe layer
tzm_height = 4.0    # Height of the TZM layer
support_span = 25.0 # Distance between the bending table supports

crack_length = 1.001  # Length of the crack - This parameter is redefined in the execution section to update it for the parametric study.
partition_length = 2.0*crack_length # Length of rectangular partition for structured meshing

plunger_radius = 1.5  # Radius of the plunger used for applying the load

# Material Properties
# Temperature influences on the elastic constants were considered, but showed only a small influence.
# Here room temperature values based on the following literature were used:
"""
TZM single crystal: J. Dickinson and P. Armstrong, Temperature dependence of the elastic constants of molybdenum, Journal of applied Physics, Vol. 38, No. 2, 602-606, (1967),
WRe single crystal: R. Ayres, G. Shannette, and D. Stein, Elastic constants of tungsten- rhenium alloys from 77 to 298K, Journal of Applied Physics, Vol. 46, No. 4, 1526-1530, (1975)
Conversion single to polycrystal: E. Kroener, Berechnung der elastischen konstanten des vielkristalls aus den konstanten des einkristalls, Zeitschrift für Physik, Vol. 151, No. 4, 504-518, (1958)
"""
# Thermal expansion coefficients were taken for Mo and W from:
"""
K. Wang and R. R. Reeber, “The role of defects on thermophysical properties:Thermal expansion of v, nb, ta, mo and w”, Materials Science and Engineering: R: Reports, vol. 23, no. 3, pp. 101–137, 1998.
"""

# Plasticity data was taken as temperature dependent perfectly plastic approximation from the publically available data of Plansee SE.
"""
TZM was fully recrystallized. Yield strength for this condition was taken from the data for sheet material found in: https://www.plansee.com/de/werkstoffe/molybdaen/eigenschaften.html
The yield strength of WRe was approximated from data for stress relieved tungsten sheet found in: https://www.plansee.com/de/werkstoffe/wolfram/eigenschaften.html
Links last accessed: 17.11.2025
"""

# TZM
tzm_poisson_ratio = 0.288302
tzm_elastic_modulus = 2*125010.0*(1+tzm_poisson_ratio)
tzm_expansion_coeff = 5.7e-06

tzm_plasticity = ((463.0, 0.0, 20.0), (365.0, 0.0, 200.0), (254.0, 0.0, 400.0), (175.0, 0.0, 600.0), (145.0, 0.0, 800.0), (130.0, 0.0, 1200.0), (75.0, 0.0, 1600.0))


# WRe
wre_poisson_ratio = 0.288302
wre_elastic_modulus = 2*161916.0*(1+wre_poisson_ratio)
wre_expansion_coeff = 4.7e-06

wre_plasticity = ((1222.0, 0.0, 200.0), (1080.0, 0.0, 300.0), (940.0, 0.0, 400.0), (806.0, 0.0, 500.0), (730.0, 0.0, 600.0), (680.0, 0.0, 700.0),(611.0, 0.0, 800.0), 
                 (550.0, 0.0, 900.0), (455.0, 0.0, 1000.0), (200.0, 0.0, 1200.0), (75.0, 0.0, 1400.0), (64.0, 0.0, 1600.0))

# Mesh Parameters
mesh_size = 0.1
wre_mesh_size = 0.1
seed_size = 0.01 # Mesh size of the edge along the crack front

elementType = CPE8 #CPE8 for plane strain / CPS8 for plane stress, append R for reduced integration

# Loading Parameters - Adapt to perform either three point bending or thermal stress analysis 
apply3PBending = True
applyThermalLoad = True

force_magnitude = 252.824179292929  # Magnitude of the force applied at the beam center for three point bending
T_initial = 1600. # Initial temperature, annealing treatment
T_coolDown = 0. # Final temperature after cooling, approx. room temperature
T_testing = 500. # Temperature during testing

# 2. MODEL CREATION FUNCTIONS

def create_parts(model, beam_length, wre_height, tzm_height, crack_length):
    """
    Creates the WRe and TZM parts of the bimaterial beam.
    Additionally create Plunger part to apply load and act as bending table supports.
    """
    # Create WRe part
    s = model.ConstrainedSketch(name='__profile__', sheetSize=50.0)
    s.rectangle(point1=(-beam_length / 2.0, 0.0), point2=(beam_length / 2.0, wre_height))
    part_wre = model.Part(name='WRe', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    part_wre.BaseShell(sketch=s)
    del model.sketches['__profile__']

    # Partition the WRe part for crack definition and meshing
    s = model.ConstrainedSketch(name='__profile__', sheetSize=80.09, transform=part_wre.MakeSketchTransform(
                                    sketchPlane=part_wre.faces[0], sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, wre_height / 2.0, 0.0)))
    part_wre.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=s)
    
    # Create crack line
    s.Line(point1=(0.0, wre_height/2.0), point2=(0.0, wre_height/2.0-crack_length))
    
    # Create lines for partitioning for structured mesh
    s.Line(point1=(-partition_length, wre_height/2.0), point2=(-partition_length, -wre_height/2.0))
    s.Line(point1=(partition_length, wre_height/2.0), point2=(partition_length, -wre_height/2.0))

    part_wre.PartitionFaceBySketch(faces=part_wre.faces.getSequenceFromMask(('[#1 ]',)), sketch=s)
    del model.sketches['__profile__']
    part_wre = model.parts['WRe']
    part_wre.PartitionEdgeByPoint(edge=part_wre.edges.findAt((-support_span/2.0, wre_height,0.0), ), point=(-support_span/2.0, wre_height,0.0))
    part_wre.PartitionEdgeByPoint(edge=part_wre.edges.findAt((support_span/2.0, wre_height,0.0), ), point=(support_span/2.0, wre_height,0.0))

    # Define left and right crack tip surfaces for initial contact (i.e. crack is introduced after annealing)
    crackEdge = part_wre.edges.findAt(((0.0, wre_height-crack_length/2.0,0.0),))    
    part_wre.Surface(name='Crack_right', side1Edges=crackEdge)
    part_wre.Surface(name='Crack_left', side2Edges=crackEdge)


    # Create TZM part
    s = model.ConstrainedSketch(name='__profile__', sheetSize=50.0)
    s.rectangle(point1=(-beam_length / 2.0, 0.0), point2=(beam_length / 2.0, -tzm_height))
    part_tzm = model.Part(name='TZM', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    part_tzm.BaseShell(sketch=s)
    del model.sketches['__profile__']

    #Partition the TZM part for applying load
    s = model.ConstrainedSketch(name='__edit__', objectToCopy=part_tzm.features['Shell planar-1'].sketch)
    part_tzm.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=s, upToFeature=part_tzm.features['Shell planar-1'])
    s.ConstructionLine(point1=(0.0, 1.25), point2=(0.0, -5.25))
    s.VerticalConstraint(addUndoState=False, entity=s.geometry[6])
    s.breakCurve(curve1=s.geometry[3], curve2=s.geometry[6], point1=(-12.375, -3.909), point2=(0.077, -6.248))
    part_tzm.features['Shell planar-1'].setValues(sketch=s)
    del model.sketches['__edit__']
    part_tzm.regenerate()


    # Create Plunger part
    s = model.ConstrainedSketch(name='__profile__', sheetSize=50.0)
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, plunger_radius))
    part_plunger = model.Part(name='Plunger', dimensionality=TWO_D_PLANAR, type=DISCRETE_RIGID_SURFACE)
    part_plunger.BaseWire(sketch=s)
    del model.sketches['__profile__']
    part_plunger.ReferencePoint(point=part_plunger.InterestingPoint(part_plunger.edges[0], CENTER))

def define_materials_and_sections(model, tzm_elastic_modulus, tzm_poisson_ratio, wre_elastic_modulus, wre_poisson_ratio):
    """
    Defines material properties and creates sections for WRe and TZM.
    """

    # Define materials
    model.Material(name='TZM')
    model.materials['TZM'].Elastic(dependencies=0, moduli=LONG_TERM, noCompression=OFF, noTension=OFF, table=((tzm_elastic_modulus, tzm_poisson_ratio), ), 
                                   temperatureDependency=OFF, type=ISOTROPIC)
    model.materials['TZM'].Expansion(type=ORTHOTROPIC, table=((tzm_expansion_coeff, tzm_expansion_coeff, 0.0), ))
    model.materials['TZM'].setValues(materialIdentifier='')
    model.materials['TZM'].setValues(description='')

    model.Material(name='WRe')
    model.materials['WRe'].Elastic(dependencies=0, moduli=LONG_TERM, noCompression=OFF, noTension=OFF, table=((wre_elastic_modulus, wre_poisson_ratio), ), 
                                   temperatureDependency=OFF, type=ISOTROPIC)
    model.materials['WRe'].Expansion(type=ORTHOTROPIC, table=((wre_expansion_coeff, wre_expansion_coeff, 0.0), ))
    model.materials['WRe'].setValues(materialIdentifier='')
    model.materials['WRe'].setValues(description='')

    if not disablePlasticity:
        model.materials['TZM'].Plastic(dataType=HALF_CYCLE, dependencies=0, hardening=ISOTROPIC, numBackstresses=1, rate=OFF, strainRangeDependency=OFF, table=tzm_plasticity, temperatureDependency=ON)
        model.materials['WRe'].Plastic(dataType=HALF_CYCLE, dependencies=0, hardening=ISOTROPIC, numBackstresses=1, rate=OFF, strainRangeDependency=OFF, table=wre_plasticity, temperatureDependency=ON)

    model.HomogeneousSolidSection(material='WRe', name='WRe', thickness=None)
    model.HomogeneousSolidSection(material='TZM', name='TZM', thickness=None)

    part_wre = model.parts['WRe']
    part_wre.Set(faces=part_wre.faces.getSequenceFromMask(('[#f ]',)), name='WRe-part')
    part_wre.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
                               region=part_wre.sets['WRe-part'], sectionName='WRe', 
                               thicknessAssignment=FROM_SECTION)
    
    part_tzm = model.parts['TZM']
    part_tzm.Set(faces=part_tzm.faces.getSequenceFromMask(('[#1 ]',)), name='TZM-Part')
    part_tzm.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
                               region=part_tzm.sets['TZM-Part'], sectionName='TZM', 
                               thicknessAssignment=FROM_SECTION)
    
    #Define Orientation for Orthotropic behaviour of WRe and TZM
    f = part_wre.faces
    faces = f.getSequenceFromMask(mask=('[#f ]', ), )
    region = Region(faces=faces)
    part_wre.MaterialOrientation(region=region, orientationType=GLOBAL, axis=AXIS_3, additionalRotationType=ROTATION_NONE, localCsys=None, fieldName='', stackDirection=STACK_3)
    
    f = part_tzm.faces
    faces = f.getSequenceFromMask(mask=('[#f ]', ), )
    region = Region(faces=faces)
    part_tzm.MaterialOrientation(region=region, orientationType=GLOBAL, axis=AXIS_3, additionalRotationType=ROTATION_NONE, localCsys=None, fieldName='', stackDirection=STACK_3)

def assemble_beam(model):
    """
    Assembles the beam by creating instances of the parts and tying them together.
    """
    a = model.rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    
    # Instance WRe
    part_wre = model.parts['WRe']
    inst_wre = a.Instance(name='WRe-1', part=part_wre, dependent=OFF)
    
    # Instance TZM
    part_tzm = model.parts['TZM']
    inst_tzm = a.Instance(name='TZM-1', part=part_tzm, dependent=OFF)

    # Place Plunger instances
    part_plunger = model.parts['Plunger']
    inst_plunger = a.Instance(name='Plunger', part=part_plunger, dependent=OFF)
    a.translate(instanceList=('Plunger', ), vector=(0.0, -tzm_height - plunger_radius, 0.0))

    inst_bendingeTable_left = a.Instance(name='BendingTableLeft', part=part_plunger, dependent=OFF)
    a.translate(instanceList=('BendingTableLeft', ), vector=(-support_span/2., wre_height + plunger_radius, 0.0))
    inst_bendingeTable_right = a.Instance(name='BendingTableRight', part=part_plunger, dependent=OFF)
    a.translate(instanceList=('BendingTableRight', ), vector=(support_span/2., wre_height + plunger_radius, 0.0))
    
    # Define surfaces for tie constraint
    bondingEdgesWRe = inst_wre.edges.findAt(((-beam_length/3.0, 0.0, 0.0),),((beam_length/3.0,0,0),),((0.0,0.0,0.0),))
    bondingEdgesTZM = inst_tzm.edges.findAt(((0.0,0.0,0.0),))

    a.Surface(name='WRe-Bond', side1Edges=bondingEdgesWRe)
    a.Surface(name='TZM-Bond', side1Edges=bondingEdgesTZM)
    
def define_step(model):
    """
    Defines the static analysis step.
    """
    model.StaticStep(initialInc=initialIncrement, minInc=minIncrement, maxInc=maxIncrement, maxNumInc=1000000, name='CoolDown', nlgeom=ON, previous='Initial')

    model.StaticStep(initialInc=initialIncrement, maxInc=maxIncrement, maxNumInc=100000000, minInc=minIncrement, 
                      name='HeatUp', previous='CoolDown', nlgeom=ON)

    model.StaticStep(adaptiveDampingRatio=0.05, continueDampingFactors=False, initialInc=initialIncrement, maxInc=maxIncrement, maxNumInc=100000000, minInc=minIncrement,
                     name='Establish Contact', previous='HeatUp', nlgeom=ON, stabilizationMagnitude=0.0002, stabilizationMethod=DISSIPATED_ENERGY_FRACTION)

    model.StaticStep(adaptiveDampingRatio=0.05, continueDampingFactors=False, initialInc=initialIncrement, maxInc=maxIncrement, maxNumInc=100000000, minInc=minIncrement,
                     name='Loading', previous='Establish Contact', nlgeom=ON, stabilizationMagnitude=0.0002, stabilizationMethod=DISSIPATED_ENERGY_FRACTION)
    
    
def apply_bcs_and_loads(model, support_span, force_magnitude):
    """
    Applies boundary conditions and loads to the assembly.
    """
    a = model.rootAssembly
    inst_wre = a.instances['WRe-1']
    inst_tzm = a.instances['TZM-1']

    # Define plunger and bending table as rigid bodies
    model.RigidBody(bodyRegion=Region(edges=a.instances['Plunger'].edges), name='Plunger', refPointRegion=Region(referencePoints= a.instances['Plunger'].referencePoints.values()))
    model.RigidBody(bodyRegion=Region(edges=a.instances['BendingTableLeft'].edges), name='BendingTableLeft', refPointRegion=Region(referencePoints= a.instances['BendingTableLeft'].referencePoints.values()))    
    model.RigidBody(bodyRegion=Region(edges=a.instances['BendingTableRight'].edges), name='BendingTableRight', refPointRegion=Region(referencePoints= a.instances['BendingTableRight'].referencePoints.values()))

    # Define support points for bending
    table_left=a.instances['BendingTableLeft']
    table_right=a.instances['BendingTableRight']
    support_edges = table_left.edges + table_right.edges
    
    bending_table=Region(edges=support_edges)

    # Apply displacement boundary condition to the support plungers to fix them
    model.DisplacementBC(name='BendingTable', createStepName='Initial', 
                         region=bending_table, u2=0.0, u1=0.0, ur3=0.0,
                         amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    
    # Apply boundary condition for beam at support points for thermal steps in y direction
    support_points_beam = inst_wre.vertices.findAt(((-support_span/2.0, wre_height,0.0),),((support_span/2.0, wre_height,0.0),))
    a.Set(name='BeamSupportPoints', vertices=support_points_beam)
    model.DisplacementBC(name='BeamSupportPoints', createStepName='Initial', 
                         region=a.sets['BeamSupportPoints'], u2=0.0, u1=UNSET, ur3=UNSET,
                         amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    
    model.boundaryConditions['BeamSupportPoints'].deactivate('Establish Contact') 
    # Keeping this boundary condition helps convergence in loading step.

    # Define the point for applying the load
    plungerEdges = a.instances['Plunger'].edges
    load_point = a.instances['Plunger'].vertices.findAt(((0.0, -tzm_height, 0.0),))
    a.Set(name='LoadPoint', vertices=load_point)
    a.Set(name='LoadPointTZM', vertices=a.instances['TZM-1'].vertices.findAt(((0.0, -tzm_height, 0.0),)))
    
    # Apply displacement boundary condition to the "loading" point to fix beam and plunger in x direction
    model.DisplacementBC(name='Plunger Fixed X', createStepName='Initial', 
                         region=Region(edges=plungerEdges), u1=0.0, u2=UNSET, ur3=0.0,
                         amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    model.DisplacementBC(name='LoadPointTZMthermal', createStepName='Initial', 
                         region=a.sets['LoadPointTZM'], u1=0.0, u2=UNSET, ur3=0.0,
                         amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    model.boundaryConditions['LoadPointTZMthermal'].deactivate('Establish Contact') 
    model.DisplacementBC(name='LoadPointTZMloading', createStepName='Establish Contact', 
                         region=a.sets['LoadPointTZM'], u1=0.0, u2=UNSET, ur3=UNSET,
                         amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

    # Apply force load for three point bending
    if apply3PBending:
        model.ConcentratedForce(name='Force/Thickness', createStepName='Establish Contact', 
            region=a.sets['LoadPoint'], cf2=force_magnitude/10.0,
            distributionType=UNIFORM, field='', localCsys=None)
        
        model.loads['Force/Thickness'].setValuesInStep(stepName='Loading', cf2=force_magnitude)
    
    # Apply temperature load to induce thermal stresses
    if applyThermalLoad:
        model.Temperature(createStepName='Initial', crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, 
                        distributionType=UNIFORM, magnitudes=(T_initial, ), name='Temp', region=Region(faces=inst_wre.faces.getSequenceFromMask(mask=('[#f ]', ), )+inst_tzm.faces.getSequenceFromMask(mask=('[#1 ]', ),),))
        model.predefinedFields['Temp'].setValuesInStep(magnitudes=(T_coolDown, ), stepName='CoolDown')
        model.predefinedFields['Temp'].setValuesInStep(magnitudes=(T_testing, ), stepName='HeatUp')

    """
    Considering the given beam in plane strain, for example in the middle of the beam, is not sensible with respect to the thermal stresses caused. 
    This is due to Abaqus assuming (typically correctly) that thermal expansion and contraction also happens in the thickness directiton, which can be assumed infinitely expanded for this case. This results in overestimated stresses in thickness direction.
    Thus, it is geometrically more sensible to assume, that the contraction in the long axis of the beam, its length, causes the relevant thermal stresses. Nonetheless, the stress field at the crack tip in the middle of the beam can still show plane strain character.
    This is best modelled by excluding the thickness thermal stresses using an orthotropic definition of the thermal expansion coefficients.
    """

    # Tie the WRe-TZM interface
    model.Tie(adjust=ON, master=a.surfaces['WRe-Bond'], name='Bond', 
              positionToleranceMethod=COMPUTED, slave=a.surfaces['TZM-Bond'], 
              thickness=ON, tieRotations=ON)
    
    # Define initial contact interaction for crack faces (i.e. bonded during annealing, crack introduced after cooling)
    model.ContactProperty('CrackFaceContact')
    model.interactionProperties['CrackFaceContact'].TangentialBehavior(formulation=ROUGH)
    model.interactionProperties['CrackFaceContact'].NormalBehavior(pressureOverclosure=HARD, allowSeparation=OFF, constraintEnforcementMethod=DEFAULT)
    model.SurfaceToSurfaceContactStd(name='CrackFaceContact', createStepName='Initial', 
        master=inst_wre.surfaces['Crack_right'], 
        slave=inst_wre.surfaces['Crack_left'], 
        sliding=SMALL, interactionProperty='CrackFaceContact', 
        initialClearance=OMIT, datumAxis=None, clearanceRegion=None)
    model.interactions['CrackFaceContact'].deactivate('HeatUp')

    # Define general contact property
    model.ContactProperty('GeneralContact')
    model.interactionProperties['GeneralContact'].TangentialBehavior(formulation=FRICTIONLESS)
    model.interactionProperties['GeneralContact'].NormalBehavior(pressureOverclosure=HARD, allowSeparation=OFF, constraintEnforcementMethod=DEFAULT)
    # As seperation of the surfaces is not expected at any time, allowSeparation is set to OFF for better convergence and to
    # create contact with the plunger at initialization. Thus no clearance is created at the plunger - beam interface in the thermal steps.

    # Define plunger - tzm contact
    model.SurfaceToSurfaceContactStd(name='PlungerContact', createStepName='Initial', 
        master=Region(side2Edges=a.instances['Plunger'].edges), 
        slave=Region(side2Edges=inst_tzm.edges.findAt(((-1.0, -tzm_height,0.0),))+inst_tzm.edges.findAt(((1.0, -tzm_height,0.0),))),
        sliding=FINITE, interactionProperty='GeneralContact', 
        initialClearance=OMIT, datumAxis=None, clearanceRegion=None, enforcement=NODE_TO_SURFACE)
    
    # Define bending table - WRe contact
    model.SurfaceToSurfaceContactStd(name='BendingTableContactLeft', createStepName='Establish Contact', 
        master=Region(side2Edges=a.instances['BendingTableLeft'].edges), 
        slave=Region(side1Edges=inst_wre.edges.findAt(((-support_span/2.0-1.0, wre_height,0.0),))+inst_wre.edges.findAt(((-support_span/2.0+1.0, wre_height,0.0),))),
        sliding=FINITE, interactionProperty='GeneralContact',
        initialClearance=OMIT, datumAxis=None, clearanceRegion=None, enforcement=NODE_TO_SURFACE)

    model.SurfaceToSurfaceContactStd(name='BendingTableContactRight', createStepName='Establish Contact', 
        master=Region(side2Edges=a.instances['BendingTableRight'].edges), 
        slave=Region(side1Edges=inst_wre.edges.findAt(((support_span/2.0-1.0, wre_height,0.0),))+inst_wre.edges.findAt(((support_span/2.0+1.0, wre_height,0.0),))),
        sliding=FINITE, interactionProperty='GeneralContact',
        initialClearance=OMIT, datumAxis=None, clearanceRegion=None, enforcement=NODE_TO_SURFACE)    

def define_crack(model):
    """
    Defines the crack, contour integral, and seam for the fracture analysis.
    """
    a = model.rootAssembly
    
    # Define crack front and tip
    crackFront=Region(vertices=model.rootAssembly.instances['WRe-1'].vertices.findAt(((0.0, wre_height-crack_length,0.0),)))
    crackTip=Region(vertices=model.rootAssembly.instances['WRe-1'].vertices.findAt(((0.0, wre_height-crack_length,0.0),)))
    # Define q-vectors for crack extension direction
    q_vectors=(((0.0, 0.0, 0.0),(0.0, -1.0, 0.0)),)
    
    # Create contour integral for crack
    # As shown, the crack is created without singular elements, this is due to the fact, that they would introduce an unphysical
    # singulartiy in the CoolDown (annealing) step, as the crack is not yet present at that time.
    a.engineeringFeatures.ContourIntegral(collapsedElementAtTip=NONE, crackFront=crackFront, crackTip=crackTip, extensionDirectionMethod=Q_VECTORS, 
    midNodePosition=0.5, name='Crack', qVectors=q_vectors, symmetric=OFF)

    # Define the crack seam
    crack_seam_edges = a.instances['WRe-1'].edges.findAt(((0.0, wre_height-crack_length+0.02,0.0),),((0.0, wre_height-0.25*crack_length,0.0),))
    a.Set(name='CrackSeam', edges=crack_seam_edges)
    a.engineeringFeatures.assignSeam(regions=a.sets['CrackSeam'])

def mesh_assembly(model, mesh_size):
    """
    Sets up mesh controls and generates the mesh for the assembly.
    """
    a = model.rootAssembly
    inst_wre = a.instances['WRe-1']
    # Seed the crack tip region
    edgeList = inst_wre.edges.findAt(((0.0, wre_height - crack_length/2.0,0.0),))+inst_wre.edges.findAt(((0.0, 0.0,0.0),))
    a.seedEdgeBySize(edges=edgeList, size=seed_size, constraint=FINER, deviationFactor=0.1)
    
    # Seed all part instances
    a.seedPartInstance(regions=(a.instances.values()), size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    a.seedPartInstance(regions=(a.instances['WRe-1'],), size=wre_mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    
    # Set mesh controls
    a.setMeshControls(elemShape=QUAD, regions=inst_wre.faces.findAt(((-support_span/2., wre_height-crack_length,0.0),))+inst_wre.faces.findAt(((support_span/2., wre_height-crack_length,0.0),)), technique=STRUCTURED)
    a.setMeshControls(algorithm=MEDIAL_AXIS, elemShape=QUAD, regions=inst_wre.faces.findAt(((0.0, wre_height-crack_length,0.0),)))

    # Set element types
    a.setElementType(elemTypes=(ElemType(elemCode=elementType, elemLibrary=STANDARD), 
                                ElemType(elemCode=elementType, elemLibrary=STANDARD)), 
                     regions=(inst_wre.faces.getSequenceFromMask(mask=('[#f ]',)) + \
                              a.instances['TZM-1'].faces.getSequenceFromMask(mask=('[#1 ]',)),))
    
    # Generate mesh for all instances
    a.generateMesh(a.instances.values())

def run_job(job, model, num_cpus):
    """
    Creates and submits the Abaqus job, and waits for its completion.
    """
    # Create history output request for K-factors
    model.HistoryOutputRequest(name='K_IC', createStepName='Loading', 
                                               contourIntegral='Crack', numberOfContours=20, 
                                               contourType=K_FACTORS, kFactorDirection=MERR,
                                                 rebar=EXCLUDE, sectionPoints=DEFAULT)
    # ABAQUS will give a warning that K-factors are requested for an elastoplastic material.
    # This can be ignored, as the K values are evaluated from the J-integral anyway. It is just convenient to compare all results based on K.
    
    model.fieldOutputRequests['F-Output-1'].setValues(timeInterval=recordSteps, variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 'NT', 'TEMP'))
    mdb.Job(name=job, model=model.name, description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=80, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, 
            echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, 
            userSubroutine='', scratch='', multiprocessingMode=DEFAULT, 
            numCpus=num_cpus, numDomains=numDomains)
            
    # Submit the job and wait for completion
    mdb.jobs[job].submit(consistencyChecking=OFF)
    mdb.jobs[job].waitForCompletion()
    KI = []
    
    odb_path = odbBase + job + '.odb'
    odb = openOdb(odb_path)
    # Get the historyOutputs dictionary for the region of interest
    hist_region = odb.steps['Loading'].historyRegions['ElementSet  ALL ELEMENTS']
    keys = hist_region.historyOutputs.keys()

    # look for keys that follow the pattern "K1 at K_IC_CRACK__<anything>_Contour_XX"
    prefix = 'K1 at K_IC_CRACK__'
    contour_marker = '_Contour_'

    for key in keys:
        if not key.startswith(prefix):
            continue
        if contour_marker not in key:
            continue
        # contour number is the substring after "_Contour_"
        contour_str = key.split(contour_marker)[-1]
        # guard: contour_str should be numeric (e.g. "01", "10")
        if not contour_str.isdigit():
            continue
        contour_no = int(contour_str)
        # read last available value for that history output
        try:
            val = hist_region.historyOutputs[key].data[-1][1]
        except Exception:
            val = None
        KI.append((contour_no, val))

    # sort by contour number to have ordered results
    KI.sort(key=lambda x: x[0])

    odb.close()

    return KI
    

# 3. MAIN EXECUTION BLOCK
if __name__ == '__main__':

    # Create a new model
    model = mdb.Model(name='Model-1')
    resultsSIFs = []

    # Call the functions to create and execute the model simulation
    create_parts(model, beam_length, wre_height, tzm_height, crack_length)
    define_materials_and_sections(model, tzm_elastic_modulus, tzm_poisson_ratio, wre_elastic_modulus, wre_poisson_ratio)
    assemble_beam(model)
    define_step(model)
    apply_bcs_and_loads(model, support_span, force_magnitude)
    define_crack(model)
    mesh_assembly(model, mesh_size)
    KI=run_job(job_name, model, num_cpus)
    print(KI)
    resultsSIFs.append((KI,crack_length,wre_height,tzm_height,force_magnitude,support_span,wre_poisson_ratio,wre_elastic_modulus/tzm_elastic_modulus))
    
    # Save results as pickle file
    with open('resultsSIFs.pkl', 'wb') as f:
        pickle.dump(resultsSIFs, f, protocol=pickle.HIGHEST_PROTOCOL)

    # Convert to JSON:
    json_ready = []
    for entry in resultsSIFs:
        KI, crack_length, wre_h, tzm_h, force, span, wre_nu, mod_ratio = entry
        json_ready.append({
            'KI': [[int(c), v] for c, v in KI],
            'crack_length': crack_length,
            'wre_height': wre_h,
            'tzm_height': tzm_h,
            'force_magnitude': force,
            'support_span': span,
            'wre_poisson_ratio': wre_nu,
            'modulus_ratio': mod_ratio
        })

    with open('resultsSIFs.json', 'w') as f:
        json.dump(json_ready, f, indent=2, sort_keys=True)