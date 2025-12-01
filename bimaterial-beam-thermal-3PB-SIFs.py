#coding=utf-8

# This script performs a parametric K_I analysis of a bimaterial beam 
# under three-point bending and/or thermal loading using Abaqus.
#
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
saveAll = True #If true, saves all models with different names. If false, overwrites the same model.
odbBase = 'C:/ABAQUS/' #Base directory in which the odb files are saved

# Geometric Parameters
beam_length = 30.0  # Total length of the beam
wre_height = 4.0    # Height of the WRe layer
tzm_height = 4.0    # Height of the TZM layer
support_span = 25.0 # Distance between the bending table supports

crack_length = 1.0  # Length of the crack - This parameter is redefined in the execution section to update it for the parametric study.
partition_length = 3.0*crack_length # Length of rectangular partition for structured meshing

# Material Properties
# Temperature influences on the elastic constants were considered, but showed only a small influence.
# Here room temperature values based on the following literature were used:
"""
TZM single crystal: J. Dickinson and P. Armstrong, Temperature dependence of the elastic constants of molybdenum, Journal of applied Physics, Vol. 38, No. 2, 602-606, (1967),
WRe single crystal: R. Ayres, G. Shannette, and D. Stein, Elastic constants of tungsten- rhenium alloys from 77 to 298K, Journal of Applied Physics, Vol. 46, No. 4, 1526-1530, (1975)
Conversion single to polycrystal: E. Kroener, Berechnung der elastischen konstanten des vielkristalls aus den konstanten des einkristalls, Zeitschrift für Physik, Vol. 151, No. 4, 504-518, (1958)
"""
# As the difference in thermal expansion coeff. is used for normalising the results, it is chosen as 0 for TZM and 1 for WRe

# TZM
tzm_poisson_ratio = 0.288302
tzm_elastic_modulus = 2*125010.0*(1+tzm_poisson_ratio)
tzm_expansion_coeff = 0


# WRe
wre_poisson_ratio = 0.288302
wre_elastic_modulus = 2*161916.0*(1+wre_poisson_ratio)
wre_expansion_coeff = 1

# Mesh Parameters
mesh_size = 0.2 # General mesh size
crack_tip_seed_number = 20 # Number of seeds along the crack front edge
crack_edge_bias_ratio = 5.0 # Bias ratio for edge seeding along the crack front
elementType = CPE8R #CPE8 for plane strain / CPS8 for plane stress


# Loading Parameters - Adapt to perform either three point bending or thermal stress analysis 
apply3PBending = False # If true, applies three point bending load
applyThermalLoad = True # If true, applies thermal load

force_magnitude = 160. # Magnitude of the force applied for three point bending (N/mm thickness)
T_initial = 1. # Initial temperature applied to the beam
T_step = 0. # Temperature in the loading step

# Job Parameters
job_name = 'K_IC_parametric'
num_cpus = 1

# 2. MODEL CREATION FUNCTIONS

def create_parts(model, beam_length, wre_height, tzm_height, crack_length):
    """
    Creates the WRe and TZM parts of the bimaterial beam.
    """

    # Create WRe part

    s = model.ConstrainedSketch(name='__profile__', sheetSize=50.0)
    s.rectangle(point1=(-beam_length / 2.0, 0.0), point2=(beam_length / 2.0, wre_height))
    part_wre = model.Part(name='WRe', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    part_wre.BaseShell(sketch=s)
    del model.sketches['__profile__']

    # Partition the WRe part for crack definition and meshing
    s = model.ConstrainedSketch(name='__profile__', sheetSize=80.09, 
                                transform=part_wre.MakeSketchTransform(
                                    sketchPlane=part_wre.faces[0], sketchPlaneSide=SIDE1, 
                                    sketchOrientation=RIGHT, origin=(0.0, wre_height / 2.0, 0.0)))
    part_wre.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=s)
    
    # Create crack line
    s.Line(point1=(0.0, wre_height/2.0), point2=(0.0, wre_height/2.0-crack_length))
    s.CircleByCenterPerimeter(center=(0.0, wre_height/2.0-crack_length), point1=(0.0, wre_height/2.0-crack_length+min(crack_length/2.0,(wre_height-crack_length)/2.0)))
    
    # Create lines for partitioning for structured mesh
    s.Line(point1=(-partition_length, wre_height/2.0), point2=(-partition_length, -wre_height/2.0))
    s.Line(point1=(partition_length, wre_height/2.0), point2=(partition_length, -wre_height/2.0))

    part_wre.PartitionFaceBySketch(faces=part_wre.faces.getSequenceFromMask(('[#1 ]',)), sketch=s)
    del model.sketches['__profile__']
    part_wre = model.parts['WRe']
    part_wre.PartitionEdgeByPoint(edge=part_wre.edges.findAt((-support_span/2.0, wre_height,0.0), ), point=(-support_span/2.0, wre_height,0.0))
    part_wre.PartitionEdgeByPoint(edge=part_wre.edges.findAt((support_span/2.0, wre_height,0.0), ), point=(support_span/2.0, wre_height,0.0))


    # Create TZM part

    s = model.ConstrainedSketch(name='__profile__', sheetSize=50.0)
    s.rectangle(point1=(-beam_length / 2.0, 0.0), point2=(beam_length / 2.0, -tzm_height))
    part_tzm = model.Part(name='TZM', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    part_tzm.BaseShell(sketch=s)
    del model.sketches['__profile__']

    # Partition the TZM part for applying load
    s = model.ConstrainedSketch(name='__edit__', objectToCopy=part_tzm.features['Shell planar-1'].sketch)
    part_tzm.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=s, 
                                           upToFeature=part_tzm.features['Shell planar-1'])
    s.ConstructionLine(point1=(0.0, 1.25), point2=(0.0, -5.25))
    s.VerticalConstraint(addUndoState=False, entity=s.geometry[6])
    s.breakCurve(curve1=s.geometry[3], curve2=s.geometry[6], 
                 point1=(-12.375, -3.909), point2=(0.077, -6.248))
    part_tzm.features['Shell planar-1'].setValues(sketch=s)
    del model.sketches['__edit__']
    part_tzm.regenerate()

def define_materials_and_sections(model, tzm_elastic_modulus, tzm_poisson_ratio, wre_elastic_modulus, wre_poisson_ratio):
    """
    Defines material properties and creates sections for WRe and TZM.
    """

    # Define materials
    model.Material(name='TZM')
    model.materials['TZM'].Elastic(table=((tzm_elastic_modulus, tzm_poisson_ratio),))
    
    model.Material(name='WRe')
    model.materials['WRe'].Elastic(table=((wre_elastic_modulus, wre_poisson_ratio),))

    # Define sections and assign them
    model.HomogeneousSolidSection(material='WRe', name='WRe', thickness=None)
    model.HomogeneousSolidSection(material='TZM', name='TZM', thickness=None)

    part_wre = model.parts['WRe']
    part_wre.Set(faces=part_wre.faces.getSequenceFromMask(('[#f ]',)), name='WRe-part')
    part_wre.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
                               region=part_wre.sets['WRe-part'], sectionName='WRe', 
                               thicknessAssignment=FROM_SECTION)
    
    #Define Orientation for Orthotropic behaviour of WRe
    f = part_wre.faces
    faces = f.getSequenceFromMask(mask=('[#f ]', ), )
    region = Region(faces=faces)
    part_wre.MaterialOrientation(region=region, orientationType=GLOBAL, axis=AXIS_3, additionalRotationType=ROTATION_NONE, localCsys=None, fieldName='', stackDirection=STACK_3)

    part_tzm = model.parts['TZM']
    part_tzm.Set(faces=part_tzm.faces.getSequenceFromMask(('[#1 ]',)), name='TZM-Part')
    part_tzm.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
                               region=part_tzm.sets['TZM-Part'], sectionName='TZM', 
                               thicknessAssignment=FROM_SECTION)

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
    
    # Define surfaces for tie constraint
    bondingEdgesWRe = inst_wre.edges.findAt(((-beam_length/3.0, 0.0, 0.0),),((beam_length/3.0,0,0),),((0.0,0.0,0.0),))
    bondingEdgesTZM = inst_tzm.edges.findAt(((0.0,0.0,0.0),))

    a.Surface(name='WRe-Bond', side1Edges=bondingEdgesWRe)
    a.Surface(name='TZM-Bond', side1Edges=bondingEdgesTZM)
    
    # Tie the two surfaces
    model.Tie(adjust=ON, master=a.surfaces['WRe-Bond'], name='Bond', 
              positionToleranceMethod=COMPUTED, slave=a.surfaces['TZM-Bond'], 
              thickness=ON, tieRotations=ON)

def define_step(model):
    """
    Defines the static analysis step.
    """
    model.StaticStep(maxNumInc=10000, name='Load', previous='Initial')

def apply_bcs_and_loads(model, support_span, force_magnitude):
    """
    Applies boundary conditions and loads to the assembly.
    """
    a = model.rootAssembly

    # Define support points for bending
    
    wre=a.instances['WRe-1']
    support_points = wre.vertices.findAt(((-support_span/2.0, wre_height,0.0),),((support_span/2.0, wre_height,0.0),))
    
    a.Set(name='BendingTablePoints', vertices=support_points)

    # Apply displacement boundary condition to the support points
    model.DisplacementBC(name='BendingTable', createStepName='Initial', 
                         region=a.sets['BendingTablePoints'], u2=0.0, u1=UNSET, ur3=UNSET,
                         amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    

    # Define the point for applying the load
    load_point_vert = a.instances['TZM-1'].vertices.findAt(((0.0,-tzm_height,0.0),))
    a.Set(name='LoadPoint', vertices=load_point_vert)
    
    # Apply displacement boundary condition to the "loading" point to fix beam in x direction
    model.DisplacementBC(name='LoadPoint', createStepName='Initial', 
                         region=a.sets['LoadPoint'], u1=0.0, u2=UNSET, ur3=UNSET,
                         amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    
    # Apply force load for three point bending
    if apply3PBending:
        model.ConcentratedForce(name='Force/Thickness', createStepName='Load', 
            region=a.sets['LoadPoint'], cf2=force_magnitude,
            distributionType=UNIFORM, field='', localCsys=None)
    
    # Apply temperature load to induce thermal stresses
    if applyThermalLoad:
        model.Temperature(createStepName='Initial', 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
        UNIFORM, magnitudes=(T_initial, ), name='Temp', region=Region(
        faces=a.instances['WRe-1'].faces.getSequenceFromMask(
        mask=('[#f ]', ), )+\
        a.instances['TZM-1'].faces.getSequenceFromMask(
        mask=('[#1 ]', ), ), 
        ))
        model.predefinedFields['Temp'].setValues(magnitudes=(T_step, ))
        model.predefinedFields['Temp'].setValuesInStep(magnitudes=(1.0, 
        ), stepName='Load')
        model.materials['TZM'].Expansion(table=((tzm_expansion_coeff, ), ))
        model.materials['WRe'].Expansion(type=ORTHOTROPIC, table=((wre_expansion_coeff, wre_expansion_coeff, 0.0), ))
    """
    Considering the given beam in plane strain, for example in the middle of the beam, is not sensible with respect to the thermal stresses caused. 
    This is due to Abaqus assuming (typically correctly) that thermal expansion and contraction also happens in the thickness directiton, which can be assumed infinitely expanded for this case. This results in overestimated stresses in thickness direction.
    Thus, it is geometrically more sensible to assume, that the contraction in the long axis of the beam, its length, causes the relevant thermal stresses. Nonetheless, the stress field at the crack tip in the middle of the beam can still show plane strain character.
    This is best modelled by excluding the thickness thermal stresses using an orthotropic definition of the thermal expansion coefficients.
    """

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
    model.rootAssembly.engineeringFeatures.ContourIntegral(
    collapsedElementAtTip=SINGLE_NODE, crackFront=crackFront, crackTip=crackTip, extensionDirectionMethod=Q_VECTORS, 
    midNodePosition=0.25, name='Crack', qVectors=q_vectors, 
    symmetric=OFF)

    # Define the crack seam
    crack_seam_edges = a.instances['WRe-1'].edges.findAt(((0.0, wre_height-crack_length+0.02,0.0),),((0.0, wre_height-0.25*crack_length,0.0),))
    a.Set(name='CrackSeam', edges=crack_seam_edges)
    a.engineeringFeatures.assignSeam(regions=a.sets['CrackSeam'])

def mesh_assembly(model, mesh_size, crack_tip_seed_number, crack_edge_bias_ratio):
    """
    Sets up mesh controls and generates the mesh for the assembly.
    """
    a = model.rootAssembly
    
    # Seed the crack tip region
    face=a.instances["WRe-1"].faces.findAt((0.0,0.1,0.0))
    edgesIdx=face.getEdges()
    edgeList=[]
    for edge in edgesIdx:
        edgeList.append(a.instances["WRe-1"].edges[edge])
    a.seedEdgeBySize(edges=edgeList, size=0.1)#min(crack_length,wre_height-crack_length)/2
    a.seedEdgeByNumber(constraint=FINER, edges=a.instances['WRe-1'].edges.findAt(((min(crack_length/2.0,(wre_height-crack_length)/2.0),wre_height-crack_length,0.0),)), 
                       number=crack_tip_seed_number)
    a.seedEdgeByBias(biasMethod=SINGLE, constraint=FINER, end2Edges=a.instances['WRe-1'].edges.findAt(((0.0,wre_height-crack_length+0.02,0.0),)), 
                     number=10, ratio=crack_edge_bias_ratio)
    
    # Seed the part instances
    a.seedPartInstance(regions=(a.instances['WRe-1'], a.instances['TZM-1']), size=mesh_size, 
                       deviationFactor=0.1, minSizeFactor=0.1)
    
    # Set mesh controls
    a.setMeshControls(elemShape=QUAD, regions=a.instances['WRe-1'].faces.getSequenceFromMask(mask=('[#3 ]',)) + \
                      a.instances['TZM-1'].faces.getSequenceFromMask(mask=('[#1 ]',)), technique=STRUCTURED)
    a.setMeshControls(elemShape=QUAD, regions=a.instances['WRe-1'].faces.getSequenceFromMask(mask=('[#4 ]',)))
    a.setMeshControls(regions=a.instances['WRe-1'].faces.getSequenceFromMask(('[#8 ]',)), technique=SWEEP)
    a.setMeshControls(algorithm=MEDIAL_AXIS, minTransition=ON, regions=a.instances['WRe-1'].faces.getSequenceFromMask(('[#4 ]',)))

    # Set element types
    a.setElementType(elemTypes=(ElemType(elemCode=elementType, elemLibrary=STANDARD), 
                                ElemType(elemCode=elementType, elemLibrary=STANDARD)), 
                     regions=(a.instances['WRe-1'].faces.getSequenceFromMask(mask=('[#f ]',)) + \
                              a.instances['TZM-1'].faces.getSequenceFromMask(mask=('[#1 ]',)),))
    
    # Generate mesh
    a.generateMesh(regions=(a.instances['WRe-1'], a.instances['TZM-1']))

def run_job(job, model, num_cpus):
    """
    Creates and submits the Abaqus job, and waits for its completion.
    """
    # Create history output request for K-factors
    model.HistoryOutputRequest(name='K_IC', createStepName='Load', 
                                               contourIntegral='Crack', numberOfContours=10, 
                                               contourType=K_FACTORS, kFactorDirection=MERR,
                                               frequency=LAST_INCREMENT, rebar=EXCLUDE, sectionPoints=DEFAULT)

    mdb.Job(name=job, model=model.name, description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=80, 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, 
            echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, 
            userSubroutine='', scratch='', multiprocessingMode=DEFAULT, 
            numCpus=num_cpus)
            
    # Submit the job and wait for completion
    mdb.jobs[job].submit(consistencyChecking=OFF)
    mdb.jobs[job].waitForCompletion()
    KI = []
    
    odb_path = odbBase + job + '.odb'
    odb = openOdb(odb_path)
    # Get the historyOutputs dictionary for the region of interest
    hist_region = odb.steps['Load'].historyRegions['ElementSet  ALL ELEMENTS']
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
    if not saveAll:
        model = mdb.Model(name='Model-1')
    resultsSIFs = []
    for i in range(1,3):  # Loop for parametric study; currently set to run once
        if saveAll:
            model = mdb.Model(name='Model_'+str(i))
        print('Running iteration '+str(i)+"\n")
        # Update parameters for each iteration if needed
        crack_length = i*0.1  # Example: vary crack length
        # Call the functions to create and execute the model simulation
        create_parts(model, beam_length, wre_height, tzm_height, crack_length)
        define_materials_and_sections(model, tzm_elastic_modulus, tzm_poisson_ratio, wre_elastic_modulus, wre_poisson_ratio)
        assemble_beam(model)
        define_step(model)
        apply_bcs_and_loads(model, support_span, force_magnitude)
        define_crack(model)
        mesh_assembly(model, mesh_size, crack_tip_seed_number, crack_edge_bias_ratio)
        KI=run_job(job_name+"_"+str(i), model, num_cpus)
        print(KI)
        resultsSIFs.append((KI,crack_length,wre_height,tzm_height,force_magnitude,support_span,wre_poisson_ratio,wre_elastic_modulus/tzm_elastic_modulus))
    with open('resultsSIFs.pkl', 'wb') as f:
        pickle.dump(resultsSIFs, f, protocol=pickle.HIGHEST_PROTOCOL)

    # Convert to JSON-serializable form: each entry -> dict, KI tuples -> lists
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