#> This is an example program which solves a weakly coupled 
#> FiniteElasticity-ALENavierStokes equation using OpenCMISS
#> calls. The 3d fluid domain has an obstacle in the middle. 
#>
#> By Chris Bradley
#>

#================================================================================================================================
#  Start Program
#================================================================================================================================

FLUID = 1
SOLID = 2
FSI = 3

#problemType = FSI
problemType = FLUID

LINEAR = 1
QUADRATIC = 2
CUBIC = 3
HERMITE = 4

NOTHING = 1
VELOCITY = 2
PRESSURE = 3
REFPRESSURE = 4

numberOfDimensions = 3

numberOfSolidXElements = 2
numberOfSolidYElements = 2
numberOfSolidZElements = 2
numberOfFluidX1Elements = 2
numberOfFluidX2Elements = 2
numberOfFluidY1Elements = 2
numberOfFluidY2Elements = 2
numberOfFluidZElements = 2

elementSize = 0.1
solidXSize = numberOfSolidXElements*elementSize
solidYSize = numberOfSolidYElements*elementSize
solidZSize = numberOfSolidZElements*elementSize
fluidX1Size = numberOfFluidX1Elements*elementSize
fluidX2Size = numberOfFluidX2Elements*elementSize
fluidY1Size = numberOfFluidY1Elements*elementSize
fluidY2Size = numberOfFluidY2Elements*elementSize
fluidZSize = numberOfFluidZElements*elementSize

uInterpolation = QUADRATIC
pInterpolation = LINEAR

#RBS = False
RBS = True

outputFrequency = 1 # Result output frequency

setupOutput = True
progressDiagnostics = True
debugLevel = 3

startTime = 0.0
stopTime  = 14.0
timeStep  = 0.1

# Inlet velocity parameters
A = 0.5
B = 2.0
C = -0.5

# Material properties
# NOTE: USE OF SI UNITS unless comment
# Low density fluid, rubber-like solid
fluidDynamicViscosity = 0.05  # kg / (m s)
fluidDensity  = 100           # kg m^-3
solidDensity  = 300           # kg m^-3
youngsModulus = 2.3E2         # Pa
poissonsRatio = 0.49          # [.]
# Neo-Hookean material law
shearModulus  = youngsModulus/(2.0*(1.0 + poissonsRatio))     # N / m^2
bulkModulus   = youngsModulus/(3.0*(1.0-2.0*poissonsRatio))
mooneyRivlin1 = 0.5*shearModulus                            # N / m^2
mooneyRivlin2 = 0.0
# Moving mesh
movingMeshKParameter   = 1.0       #default

solidPRef = 0.0
fluidPRef = 0.0

solidPInit = -mooneyRivlin1
fluidPInit = fluidPRef

# Set solver parameters
fsiDynamicSolverTheta    = [1.0]
nonlinearMaximumIterations      = 100000000 #default: 100000
nonlinearRelativeTolerance      = 1.0E-4    #default: 1.0E-05
nonlinearAbsoluteTolerance      = 1.0E-4    #default: 1.0E-10
nonlinearMaxFunctionEvaluations = 100000
nonlinearLinesearchAlpha        = 1.0
linearMaximumIterations      = 100000000 #default: 100000
linearRelativeTolerance      = 1.0E-4    #default: 1.0E-05
linearAbsoluteTolerance      = 1.0E-4    #default: 1.0E-10
linearDivergenceTolerance    = 1.0E5     #default: 1.0E5
linearRestartValue           = 30        #default: 30

#================================================================================================================================
#  Should not need to change anything below here.
#================================================================================================================================

if (uInterpolation == HERMITE):
    numberOfNodesXi = 2
else:
    numberOfNodesXi = uInterpolation+1

numberOfSolidXNodes = numberOfSolidXElements*(numberOfNodesXi-1)+1
numberOfSolidYNodes = numberOfSolidYElements*(numberOfNodesXi-1)+1
numberOfSolidZNodes = numberOfSolidZElements*(numberOfNodesXi-1)+1
numberOfSolidNodesPerZ = numberOfSolidXNodes*numberOfSolidYNodes
numberOfSolidNodes = numberOfSolidNodesPerZ*(numberOfSolidZElements*(numberOfNodesXi-1)+1)
numberOfFluidX1Nodes = numberOfFluidX1Elements*(numberOfNodesXi-1)+1
numberOfFluidX2Nodes = numberOfFluidX2Elements*(numberOfNodesXi-1)+1
numberOfFluidY1Nodes = numberOfFluidY1Elements*(numberOfNodesXi-1)+1
numberOfFluidY2Nodes = numberOfFluidY2Elements*(numberOfNodesXi-1)+1
numberOfFluidZNodes = numberOfFluidZElements*(numberOfNodesXi-1)+1
numberOfFluidXNodes1 = ((numberOfFluidX1Elements+numberOfSolidXElements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+1)
numberOfFluidXNodes2 = ((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)
numberOfFluidYNodes1 = ((numberOfFluidY1Elements+numberOfSolidYElements+numberOfFluidY2Elements)*(numberOfNodesXi-1)+1)
numberOfFluidYNodes2 = ((numberOfFluidY1Elements+numberOfFluidY2Elements)*(numberOfNodesXi-1)+2)
numberOfFluidNodesPerZ1 = ((numberOfFluidX1Elements+numberOfSolidXElements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+1)* \
                          ((numberOfFluidY1Elements+numberOfFluidY2Elements)*(numberOfNodesXi-1)+2)+ \
                          ((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)*(numberOfSolidYElements*(numberOfNodesXi-1)-1)
numberOfFluidNodesPerZ2 = ((numberOfFluidX1Elements+numberOfFluidX2Elements+numberOfSolidXElements)*(numberOfNodesXi-1)+1)* \
                          ((numberOfFluidY1Elements+numberOfFluidY2Elements+numberOfSolidYElements)*(numberOfNodesXi-1)+1)
numberOfFluidNodes = numberOfFluidNodesPerZ1*(numberOfSolidZNodes-1)+numberOfFluidNodesPerZ2*(numberOfFluidZNodes)
numberOfInterfaceNodesPerZ = 2*(numberOfSolidXElements*(numberOfNodesXi-1)+1)+2*(numberOfSolidYElements*(numberOfNodesXi-1)-1)
numberOfInterfaceNodes = 2*numberOfSolidXNodes*(numberOfSolidZNodes-1)+2*(numberOfSolidYNodes-2)*(numberOfSolidZNodes-1)+ \
                         numberOfSolidXNodes*numberOfSolidYNodes
numberOfSolidElementsPerZ = numberOfSolidXElements*numberOfSolidYElements
numberOfSolidElements = numberOfSolidElementsPerZ*numberOfSolidZElements
numberOfFluidXElements1 = numberOfFluidX1Elements+numberOfFluidX2Elements+numberOfSolidXElements
numberOfFluidXElements2 = numberOfFluidX1Elements+numberOfFluidX2Elements
numberOfFluidYElements1 = numberOfFluidY1Elements+numberOfFluidY2Elements+numberOfSolidYElements
numberOfFluidYElements2 = numberOfFluidY1Elements+numberOfFluidY2Elements
numberOfFluidElementsPerZ1 = numberOfFluidXElements1*numberOfFluidYElements1-numberOfSolidElementsPerZ
numberOfFluidElementsPerZ2 = numberOfFluidXElements1*numberOfFluidYElements1
numberOfFluidElements = numberOfFluidElementsPerZ1*numberOfSolidZElements+numberOfFluidElementsPerZ2*numberOfFluidZElements
numberOfInterfaceElementsPerZ = 2*numberOfSolidXElements+2*numberOfSolidYElements
numberOfInterfaceElements = 2*numberOfSolidXElements*numberOfSolidZElements+2*numberOfSolidYElements*numberOfSolidZElements+\
                            numberOfSolidXElements*numberOfSolidYElements
numberOfLocalNodes = numberOfNodesXi*numberOfNodesXi*numberOfNodesXi
numberOfLocalInterfaceNodes = numberOfNodesXi*numberOfNodesXi
localNodeIdx000 = 0
localNodeIdx100 = numberOfNodesXi-1
localNodeIdx010 = numberOfNodesXi*(numberOfNodesXi-1)
localNodeIdx110 = numberOfNodesXi*numberOfNodesXi-1
localNodeIdx001 = numberOfNodesXi*numberOfNodesXi*(numberOfNodesXi-1)
localNodeIdx101 = numberOfNodesXi-1+numberOfNodesXi*numberOfNodesXi*(numberOfNodesXi-1)
localNodeIdx011 = numberOfNodesXi*(numberOfNodesXi-1)+numberOfNodesXi*numberOfNodesXi*(numberOfNodesXi-1)
localNodeIdx111 = numberOfLocalNodes-1

contextUserNumber = 1

solidCoordinateSystemUserNumber     = 1
fluidCoordinateSystemUserNumber     = 2
interfaceCoordinateSystemUserNumber = 3
  
solidRegionUserNumber = 1
fluidRegionUserNumber = 2
interfaceUserNumber   = 3

uBasisUserNumber = 1
pBasisUserNumber = 2
interfaceBasisUserNumber = 3

solidMeshUserNumber     = 1
fluidMeshUserNumber     = 2
interfaceMeshUserNumber = 3
movingMeshUserNumber    = 4
  
solidDecompositionUserNumber     = 1
fluidDecompositionUserNumber     = 2
interfaceDecompositionUserNumber = 3
  
solidGeometricFieldUserNumber     = 11
solidFibreFieldUserNumber     = 12
solidEquationsSetFieldUserNumber = 13
solidDependentFieldUserNumber = 14
solidMaterialsFieldUserNumber = 15
solidSourceFieldUserNumber = 16

fluidGeometricFieldUserNumber     = 21
fluidEquationsSetFieldUserNumber = 22
fluidDependentFieldUserNumber = 23
fluidMaterialsFieldUserNumber = 24
fluidIndependentFieldUserNumber = 25
bcCellMLModelsFieldUserNumber = 26
bcCellMLStateFieldUserNumber = 27
bcCellMLParametersFieldUserNumber = 28
bcCellMLIntermediateFieldUserNumber = 29

movingMeshEquationsSetFieldUserNumber = 31
movingMeshDependentFieldUserNumber    = 32
movingMeshMaterialsFieldUserNumber    = 33
movingMeshIndependentFieldUserNumber  = 34

interfaceGeometricFieldUserNumber = 41
interfaceLagrangeFieldUserNumber  = 42
 
solidEquationsSetUserNumber  = 1
fluidEquationsSetUserNumber  = 2
movingMeshEquationsSetUserNumber = 3

bcCellMLUserNumber = 1

interfaceConditionUserNumber = 1
  
fsiProblemUserNumber = 1
  
DynamicSolverIndex = 1
LinearSolverMovingMeshIndex = 2
  
IndependentFieldMovingMeshUserNumberK = 1
LinearSolverMovingMeshEquationsUserNumber = 122

SolidEquationsSetIndex  = 1
FluidEquationsSetIndex  = 2
InterfaceConditionIndex = 1
SolidMeshIndex = 1
FluidMeshIndex = 2

#================================================================================================================================
#  Define functions
#================================================================================================================================

def GetElementNodes3D(elementNumber,localNodes3D,numberOfXNodes1,numberOfXNodes2,numberOfXNodes3,numberOfZNodes,zOffset,zDelta1,zDelta2,zDelta3):
    if (uInterpolation == LINEAR or uInterpolation == HERMITE):
        localNodes3D[localNodeIdx100] = localNodes3D[localNodeIdx000]+(numberOfNodesXi-1)
        localNodes3D[localNodeIdx010] = localNodes3D[localNodeIdx000]+numberOfXNodes1
        localNodes3D[localNodeIdx110] = localNodes3D[localNodeIdx010]+(numberOfNodesXi-1)
        localNodes3D[localNodeIdx001] = localNodes3D[localNodeIdx000]+(numberOfNodesXi-1)*numberOfZNodes+zOffset
        localNodes3D[localNodeIdx101] = localNodes3D[localNodeIdx100]+(numberOfNodesXi-1)*numberOfZNodes+zOffset
        localNodes3D[localNodeIdx011] = localNodes3D[localNodeIdx010]+(numberOfNodesXi-1)*numberOfZNodes+zOffset+zDelta3
        localNodes3D[localNodeIdx111] = localNodes3D[localNodeIdx110]+(numberOfNodesXi-1)*numberOfZNodes+zOffset+zDelta3
        uNodes3D = [localNodes3D[localNodeIdx000],localNodes3D[localNodeIdx100],localNodes3D[localNodeIdx010],localNodes3D[localNodeIdx110], \
                    localNodes3D[localNodeIdx001],localNodes3D[localNodeIdx101],localNodes3D[localNodeIdx011],localNodes3D[localNodeIdx111]]
    elif (uInterpolation == QUADRATIC):
        localNodes3D[localNodeIdx100] = localNodes3D[localNodeIdx000]+(numberOfNodesXi-1)
        localNodes3D[localNodeIdx010] = localNodes3D[localNodeIdx000]+numberOfXNodes1+numberOfXNodes3
        localNodes3D[localNodeIdx110] = localNodes3D[localNodeIdx010]+(numberOfNodesXi-1)
        localNodes3D[localNodeIdx001] = localNodes3D[localNodeIdx000]+(numberOfNodesXi-1)*numberOfZNodes+zOffset
        localNodes3D[localNodeIdx101] = localNodes3D[localNodeIdx100]+(numberOfNodesXi-1)*numberOfZNodes+zOffset
        localNodes3D[localNodeIdx011] = localNodes3D[localNodeIdx010]+(numberOfNodesXi-1)*numberOfZNodes+zOffset+zDelta1+zDelta3
        localNodes3D[localNodeIdx111] = localNodes3D[localNodeIdx110]+(numberOfNodesXi-1)*numberOfZNodes+zOffset+zDelta1+zDelta3
        localNodes3D[1] = localNodes3D[localNodeIdx000] + 1
        localNodes3D[3] = localNodes3D[localNodeIdx000] + numberOfXNodes1
        localNodes3D[4] = localNodes3D[3] + 1
        localNodes3D[5] = localNodes3D[4] + 1
        localNodes3D[7] = localNodes3D[localNodeIdx010] + 1
        localNodes3D[9] = localNodes3D[0]+numberOfZNodes
        localNodes3D[10] = localNodes3D[1]+numberOfZNodes
        localNodes3D[11] = localNodes3D[2]+numberOfZNodes
        localNodes3D[12] = localNodes3D[3]+numberOfZNodes
        localNodes3D[13] = localNodes3D[4]+numberOfZNodes
        localNodes3D[14] = localNodes3D[5]+numberOfZNodes
        localNodes3D[15] = localNodes3D[6]+numberOfZNodes
        localNodes3D[16] = localNodes3D[7]+numberOfZNodes
        localNodes3D[17] = localNodes3D[8]+numberOfZNodes
        localNodes3D[19] = localNodes3D[1]+2*numberOfZNodes+zOffset
        localNodes3D[21] = localNodes3D[3]+2*numberOfZNodes+zOffset+zDelta1
        localNodes3D[22] = localNodes3D[4]+2*numberOfZNodes+zOffset+zDelta1
        localNodes3D[23] = localNodes3D[5]+2*numberOfZNodes+zOffset+zDelta1
        localNodes3D[25] = localNodes3D[7]+2*numberOfZNodes+zOffset+zDelta1+zDelta3
        uNodes3D = localNodes3D
    elif (uInterpolation == CUBIC):
        localNodes3D[localNodeIdx100] = localNodes3D[localNodeIdx000]+(numberOfNodesXi-1)
        localNodes3D[localNodeIdx010] = localNodes3D[localNodeIdx000]+numberOfXNodes1+numberOfXNodes2+numberOfXNodes3
        localNodes3D[localNodeIdx110] = localNodes3D[localNodeIdx010]+(numberOfNodesXi-1)
        localNodes3D[localNodeIdx001] = localNodes3D[localNodeIdx000]+(numberOfNodesXi-1)*numberOfZNodes+zOffset
        localNodes3D[localNodeIdx101] = localNodes3D[localNodeIdx100]+(numberOfNodesXi-1)*numberOfZNodes+zOffset
        localNodes3D[localNodeIdx011] = localNodes3D[localNodeIdx010]+(numberOfNodesXi-1)*numberOfZNodes+zOffset+zDelta1+zDelta2+zDelta3
        localNodes3D[localNodeIdx111] = localNodes3D[localNodeIdx110]+(numberOfNodesXi-1)*numberOfZNodes+zOffset+zDelta1+zDelta2+zDelta3
        localNodes3D[1] = localNodes3D[localNodeIdx000] + 1
        localNodes3D[2] = localNodes3D[1] + 1
        localNodes3D[4] = localNodes3D[localNodeIdx000] + numberOfXNodes1
        localNodes3D[5] = localNodes3D[4] + 1
        localNodes3D[6] = localNodes3D[5] + 1
        localNodes3D[7] = localNodes3D[6] + 1
        localNodes3D[8] = localNodes3D[4] + numberOfXNodes2
        localNodes3D[9] = localNodes3D[8] + 1
        localNodes3D[10] = localNodes3D[9] + 1
        localNodes3D[11] = localNodes3D[10] + 1
        localNodes3D[13] = localNodes3D[localNodeIdx010] + 1
        localNodes3D[14] = localNodes3D[13] + 1
        localNodes3D[16] = localNodes3D[0]+numberOfZNodes
        localNodes3D[17] = localNodes3D[1]+numberOfZNodes
        localNodes3D[18] = localNodes3D[2]+numberOfZNodes
        localNodes3D[19] = localNodes3D[3]+numberOfZNodes
        localNodes3D[20] = localNodes3D[4]+numberOfZNodes
        localNodes3D[21] = localNodes3D[5]+numberOfZNodes
        localNodes3D[22] = localNodes3D[6]+numberOfZNodes
        localNodes3D[23] = localNodes3D[7]+numberOfZNodes
        localNodes3D[24] = localNodes3D[8]+numberOfZNodes
        localNodes3D[25] = localNodes3D[9]+numberOfZNodes
        localNodes3D[26] = localNodes3D[10]+numberOfZNodes
        localNodes3D[27] = localNodes3D[11]+numberOfZNodes
        localNodes3D[28] = localNodes3D[12]+numberOfZNodes
        localNodes3D[29] = localNodes3D[13]+numberOfZNodes
        localNodes3D[30] = localNodes3D[14]+numberOfZNodes
        localNodes3D[31] = localNodes3D[15]+numberOfZNodes
        localNodes3D[32] = localNodes3D[0]+2*numberOfZNodes
        localNodes3D[33] = localNodes3D[1]+2*numberOfZNodes
        localNodes3D[34] = localNodes3D[2]+2*numberOfZNodes
        localNodes3D[35] = localNodes3D[3]+2*numberOfZNodes
        localNodes3D[36] = localNodes3D[4]+2*numberOfZNodes
        localNodes3D[37] = localNodes3D[5]+2*numberOfZNodes
        localNodes3D[38] = localNodes3D[6]+2*numberOfZNodes
        localNodes3D[39] = localNodes3D[7]+2*numberOfZNodes
        localNodes3D[40] = localNodes3D[8]+2*numberOfZNodes
        localNodes3D[41] = localNodes3D[9]+2*numberOfZNodes
        localNodes3D[42] = localNodes3D[10]+2*numberOfZNodes
        localNodes3D[43] = localNodes3D[11]+2*numberOfZNodes
        localNodes3D[44] = localNodes3D[12]+2*numberOfZNodes
        localNodes3D[45] = localNodes3D[13]+2*numberOfZNodes
        localNodes3D[46] = localNodes3D[14]+2*numberOfZNodes
        localNodes3D[47] = localNodes3D[15]+2*numberOfZNodes
        localNodes3D[49] = localNodes3D[1]+3*numberOfZNodes+zOffset
        localNodes3D[50] = localNodes3D[2]+3*numberOfZNodes+zOffset
        localNodes3D[52] = localNodes3D[4]+3*numberOfZNodes+zOffset+zDelta1
        localNodes3D[53] = localNodes3D[5]+3*numberOfZNodes+zOffset+zDelta1
        localNodes3D[54] = localNodes3D[6]+3*numberOfZNodes+zOffset+zDelta1
        localNodes3D[55] = localNodes3D[7]+3*numberOfZNodes+zOffset+zDelta1
        localNodes3D[56] = localNodes3D[8]+3*numberOfZNodes+zOffset+zDelta1+zDelta2
        localNodes3D[57] = localNodes3D[9]+3*numberOfZNodes+zOffset+zDelta1+zDelta2
        localNodes3D[58] = localNodes3D[10]+3*numberOfZNodes+zOffset+zDelta1+zDelta2
        localNodes3D[59] = localNodes3D[11]+3*numberOfZNodes+zOffset+zDelta1+zDelta2
        localNodes3D[61] = localNodes3D[13]+3*numberOfZNodes+zOffset+zDelta1+zDelta2+zDelta3
        localNodes3D[62] = localNodes3D[14]+3*numberOfZNodes+zOffset+zDelta1+zDelta2+zDelta3
        uNodes3D = localNodes3D
    else:
        print('Invalid u interpolation')
        exit()                    
    pNodes3D = [localNodes3D[localNodeIdx000],localNodes3D[localNodeIdx100],localNodes3D[localNodeIdx010],localNodes3D[localNodeIdx110], \
              localNodes3D[localNodeIdx001],localNodes3D[localNodeIdx101],localNodes3D[localNodeIdx011],localNodes3D[localNodeIdx111]]
    if (debugLevel > 2):
        print('    Element %8d' %(elementNumber))
        print('      U Nodes: '+str(uNodes3D))
        print('      P Nodes: '+str(pNodes3D))
    return uNodes3D,pNodes3D

def GetElementNodes2D(elementNumber,localNodes2D,numberOfXNodes1,numberOfXNodes2,numberOfXNodes3):
    localNodes2D[localNodeIdx100] = localNodes2D[localNodeIdx000]+(numberOfNodesXi-1)
    if (uInterpolation == LINEAR or uInterpolation == HERMITE):
        localNodes2D[localNodeIdx010] = localNodes2D[localNodeIdx000]+numberOfXNodes1
    elif (uInterpolation == QUADRATIC):
        localNodes2D[localNodeIdx010] = localNodes2D[localNodeIdx000]+numberOfXNodes1+numberOfXNodes3
    elif (uInterpolation == CUBIC):
        localNodes2D[localNodeIdx010] = localNodes2D[localNodeIdx000]+numberOfXNodes1+numberOfXNodes2+numberOfXNodes3
    else:
        print('Invalid u interpolation')
        exit()                    
    localNodes2D[localNodeIdx110] = localNodes2D[localNodeIdx010]+(numberOfNodesXi-1)
    if (uInterpolation == LINEAR or uInterpolation == HERMITE):
        uNodes2D = [localNodes2D[localNodeIdx000],localNodes2D[localNodeIdx100],localNodes2D[localNodeIdx010],localNodes2D[localNodeIdx110]]
    elif (uInterpolation == QUADRATIC):
        localNodes2D[1] = localNodes2D[localNodeIdx000] + 1
        localNodes2D[3] = localNodes2D[localNodeIdx000] + numberOfXNodes1
        localNodes2D[4] = localNodes2D[3] + 1
        localNodes2D[5] = localNodes2D[4] + 1
        localNodes2D[7] = localNodes2D[localNodeIdx010] + 1
        uNodes2D = localNodes2D
    elif (uInterpolation == CUBIC):
        localNodes2D[1] = localNodes2D[localNodeIdx000] + 1
        localNodes2D[2] = localNodes2D[1] + 1
        localNodes2D[4] = localNodes2D[localNodeIdx000] + numberOfXNodes1
        localNodes2D[5] = localNodes2D[4] + 1
        localNodes2D[6] = localNodes2D[5] + 1
        localNodes2D[7] = localNodes2D[6] + 1
        localNodes2D[8] = localNodes2D[4] + numberOfXNodes2
        localNodes2D[9] = localNodes2D[8] + 1
        localNodes2D[10] = localNodes2D[9] + 1
        localNodes2D[11] = localNodes2D[10] + 1
        localNodes2D[13] = localNodes2D[localNodeIdx010] + 1
        localNodes2D[14] = localNodes2D[13] + 1
        uNodes2D = localNodes2D
    else:
        print('Invalid u interpolation')
        exit()                    
    pNodes2D = [localNodes2D[localNodeIdx000],localNodes2D[localNodeIdx100],localNodes2D[localNodeIdx010],localNodes2D[localNodeIdx110]]
    if (debugLevel > 2):
        print('    Element %8d' %(elementNumber))
        print('      U Nodes: '+str(uNodes2D))
        print('      P Nodes: '+str(pNodes2D))
    return uNodes2D,pNodes2D

def SetNodeParameters3D(nodeNumber,field,xPosition,yPosition,zPosition):
    field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,xPosition)
    field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,yPosition)
    field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,zPosition)
    if (debugLevel > 2):
        print('      Node        %d:' % (nodeNumber))
        print('         Position         = [ %.2f, %.2f, %.2f ]' % (xPosition,yPosition,zPosition))                 
    if (uInterpolation == HERMITE):
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,1.0)
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,0.0)
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,0.0)
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,0.0)
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,1.0)
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,0.0)
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,0.0)
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,0.0)
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,1.0)
        if (debugLevel > 2):
            print('        S1 derivative    = [ %.2f, %.2f, %.2f ]' % (1.0,0.0,0.0))                 
            print('        S2 derivative    = [ %.2f, %.2f, %.2f ]' % (0.0,1.0,0.0))                 
            print('        S3 derivative    = [ %.2f, %.2f, %.2f ]' % (0.0,0.0,0.0))                                 

#================================================================================================================================
#  Initialise OpenCMISS
#================================================================================================================================

# Import the libraries (OpenCMISS,python,numpy,scipy)
import numpy,csv,time,sys,os,pdb
from opencmiss.opencmiss import OpenCMISS_Python as oc

# Diagnostics
#oc.DiagnosticsSetOn(oc.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",[""])
#oc.ErrorHandlingModeSet(oc.ErrorHandlingModes.TRAP_ERROR)
oc.OutputSetOn("Testing")

context = oc.Context()
context.Create(contextUserNumber)

worldRegion = oc.Region()
context.WorldRegionGet(worldRegion)

# Get the computational nodes info
computationEnvironment = oc.ComputationEnvironment()
context.ComputationEnvironmentGet(computationEnvironment)
numberOfComputationalNodes = computationEnvironment.NumberOfWorldNodesGet()
computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()
        
#================================================================================================================================
#  Initial Data & Default Values
#================================================================================================================================

# (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
fluidEquationsSetOutputType = oc.EquationsSetOutputTypes.NONE
#fluidEquationsSetOutputType = oc.EquationsSetOutputTypes.PROGRESS
fluidEquationsOutputType = oc.EquationsOutputTypes.NONE
#fluidEquationsOutputType = oc.EquationsOutputTypes.TIMING
#fluidEquationsOutputType = oc.EquationsOutputTypes.MATRIX
#fluidEquationsOutputType = oc.EquationsOutputTypes.ELEMENT_MATRIX
solidEquationsSetOutputType = oc.EquationsSetOutputTypes.NONE
#solidEquationsSetOutputType = oc.EquationsSetOutputTypes.PROGRESS
solidEquationsOutputType = oc.EquationsOutputTypes.NONE
#solidEquationsOutputType = oc.EquationsOutputTypes.TIMING
#solidEquationsOutputType = oc.EquationsOutputTypes.MATRIX
#solidEquationsOutputType = oc.EquationsOutputTypes.ELEMENT_MATRIX
movingMeshEquationsSetOutputType = oc.EquationsSetOutputTypes.NONE
#movingMeshEquationsSetOutputType = oc.EquationsSetOutputTypes.PROGRESS
movingMeshEquationsOutputType = oc.EquationsOutputTypes.NONE
#movingMeshEquationsOutputType = oc.EquationsOutputTypes.TIMING
#movingMeshEquationsOutputType = oc.EquationsOutputTypes.MATRIX
#movingMeshEquationsOutputType = oc.EquationsOutputTypes.ELEMENT_MATRIX
interfaceConditionOutputType = oc.InterfaceConditionOutputTypes.NONE
#interfaceConditionOutputType = oc.InterfaceConditionOutputTypes.PROGRESS
interfaceEquationsOutputType = oc.EquationsOutputTypes.NONE
#interfaceEquationsOutputType = oc.EquationsOutputTypes.TIMING
#interfaceEquationsOutputType = oc.EquationsOutputTypes.PROGRESS
#interfaceEquationsOutputType = oc.EquationsOutputTypes.MATRIX
#interfaceEquationsOutputType = oc.EquationsOutputTypes.ELEMENT_MATRIX
# (NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
movingMeshLinearSolverOutputType = oc.SolverOutputTypes.NONE
#movingMeshLinearSolverOutputType = oc.SolverOutputTypes.PROGRESS
#movingMeshLinearSolverOutputType = oc.SolverOutputTypes.MATRIX
#fsiDynamicSolverOutputType = oc.SolverOutputTypes.NONE
fsiDynamicSolverOutputType = oc.SolverOutputTypes.PROGRESS
#fsiDynamicSolverOutputType = oc.SolverOutputTypes.MATRIX
#fsiNonlinearSolverOutputType = oc.SolverOutputTypes.NONE
fsiNonlinearSolverOutputType = oc.SolverOutputTypes.PROGRESS
#fsiNonlinearSolverOutputType = oc.SolverOutputTypes.MATRIX
#fsiLinearSolverOutputType = oc.SolverOutputTypes.NONE
fsiLinearSolverOutputType = oc.SolverOutputTypes.PROGRESS
#fsiLinearSolverOutputType = oc.SolverOutputTypes.MATRIX

if (setupOutput):
    print('SUMMARY')
    print('=======')
    print(' ')
    print('  Temporal parameters')
    print('  -------------------')
    print(' ')
    print('  Start time:     %.3f s' % (startTime))
    print('  Stop time:      %.3f s' % (stopTime))
    print('  Time increment: %.5f s' % (timeStep))
    print(' ')
    print('  Material parameters')
    print('  -------------------')
    print(' ')
    if (problemType != SOLID):
        print('    Fluid:')
        print('      Dynamic viscosity: {0:.3f} kg.m^-1.s^-1'.format(fluidDynamicViscosity))
        print('      Density: {0:.3f} kg.m^-3'.format(fluidDensity))
        print(' ')
    if (problemType != FLUID):
        print('    Solid:')
        print('      Density: {0:.3f} kg.m^-3'.format(solidDensity))
        print('      Young\'s modulus: {0:.3f} Pa'.format(youngsModulus))
        print('      Poisson\'s ratio: {0:.3f}'.format(poissonsRatio))
        print('      Neo-Hookean constant: {0:.3f}'.format(mooneyRivlin1))
        print(' ')
    print('  Mesh parameters')
    print('  -------------------')
    print(' ')
    print('    Number of dimensions: {0:d}'.format(numberOfDimensions))
    print('    uInterpolation: {}'.format(uInterpolation))
    print('    pInterpolation: {}'.format(pInterpolation))
    if (problemType != SOLID):
        print('    Fluid:')
        print('      Number of X1 elements: {0:d}'.format(numberOfFluidX1Elements))
        print('      Number of X2 elements: {0:d}'.format(numberOfFluidX2Elements))
        print('      Number of Y1 elements: {0:d}'.format(numberOfFluidY1Elements))
        print('      Number of Y2 elements: {0:d}'.format(numberOfFluidY2Elements))
        print('      Number of Z  elements: {0:d}'.format(numberOfFluidZElements))        
        print('      Number of nodes: {0:d}'.format(numberOfFluidNodes))
        print('      Number of elements: {0:d}'.format(numberOfFluidElements))
    if (problemType != FLUID):
        print('    Solid:')
        print('      Number of X elements: {0:d}'.format(numberOfSolidXElements))
        print('      Number of Y elements: {0:d}'.format(numberOfSolidYElements))
        print('      Number of Z elements: {0:d}'.format(numberOfSolidZElements))        
        print('      Number of nodes: {0:d}'.format(numberOfSolidNodes))
        print('      Number of elements: {0:d}'.format(numberOfSolidElements))
    if (problemType == FSI):
        print('    Interface:')
        print('      Number of nodes: {0:d}'.format(numberOfInterfaceNodes))
        print('      Number of elements: {0:d}'.format(numberOfInterfaceElements))
 
#================================================================================================================================
#  Coordinate Systems
#================================================================================================================================

if (progressDiagnostics):
    print(' ')
    print('Coordinate systems ...')

if (problemType != FLUID):
    # Create a RC coordinate system for the solid region
    solidCoordinateSystem = oc.CoordinateSystem()
    solidCoordinateSystem.CreateStart(solidCoordinateSystemUserNumber,context)
    solidCoordinateSystem.DimensionSet(numberOfDimensions)
    solidCoordinateSystem.CreateFinish()
if (problemType != SOLID):
    # Create a RC coordinate system for the fluid region
    fluidCoordinateSystem = oc.CoordinateSystem()
    fluidCoordinateSystem.CreateStart(fluidCoordinateSystemUserNumber,context)
    fluidCoordinateSystem.DimensionSet(numberOfDimensions)
    fluidCoordinateSystem.CreateFinish()
if (problemType == FSI):
    # Create a RC coordinate system for the interface region
    interfaceCoordinateSystem = oc.CoordinateSystem()
    interfaceCoordinateSystem.CreateStart(interfaceCoordinateSystemUserNumber,context)
    interfaceCoordinateSystem.DimensionSet(numberOfDimensions)
    interfaceCoordinateSystem.CreateFinish()

if (progressDiagnostics):
    print('Coordinate systems ... Done')
  
#================================================================================================================================
#  Regions
#================================================================================================================================

if (progressDiagnostics):
    print('Regions ...')

if (problemType != FLUID):
    # Create a solid region
    solidRegion = oc.Region()
    solidRegion.CreateStart(solidRegionUserNumber,worldRegion)
    solidRegion.label = 'SolidRegion'
    solidRegion.coordinateSystem = solidCoordinateSystem
    solidRegion.CreateFinish()
if (problemType != SOLID):
    # Create a fluid region
    fluidRegion = oc.Region()
    fluidRegion.CreateStart(fluidRegionUserNumber,worldRegion)
    fluidRegion.label = 'FluidRegion'
    fluidRegion.coordinateSystem = fluidCoordinateSystem
    fluidRegion.CreateFinish()
    
if (progressDiagnostics):
    print('Regions ... Done')

#================================================================================================================================
#  Bases
#================================================================================================================================

if (progressDiagnostics):
    print('Basis functions ...')

numberOfGaussXi =3
    
pBasis = oc.Basis()
pBasis.CreateStart(pBasisUserNumber,context)
pBasis.type = oc.BasisTypes.LAGRANGE_HERMITE_TP
pBasis.numberOfXi = 3
pBasis.interpolationXi = [oc.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
pBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
pBasis.CreateFinish()

uBasis = oc.Basis()
uBasis.CreateStart(uBasisUserNumber,context)
uBasis.type = oc.BasisTypes.LAGRANGE_HERMITE_TP
uBasis.numberOfXi = 3
if (uInterpolation == LINEAR):
    uBasis.interpolationXi = [oc.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
elif (uInterpolation == QUADRATIC):
    uBasis.interpolationXi = [oc.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*3
elif (uInterpolation == CUBIC):
    uBasis.interpolationXi = [oc.BasisInterpolationSpecifications.CUBIC_LAGRANGE]*3
elif (uInterpolation == HERMITE):
    uBasis.interpolationXi = [oc.BasisInterpolationSpecifications.CUBIC_HERMITE]*3
else:
    print('Invalid u interpolation')
    exit()
uBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
uBasis.CreateFinish()

if (problemType == FSI):
    interfaceBasis = oc.Basis()
    interfaceBasis.CreateStart(interfaceBasisUserNumber,context)
    interfaceBasis.type = oc.BasisTypes.LAGRANGE_HERMITE_TP
    interfaceBasis.numberOfXi = 2
    if (uInterpolation == LINEAR):
        interfaceBasis.interpolationXi = [oc.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*2
    elif (uInterpolation == QUADRATIC):
        interfaceBasis.interpolationXi = [oc.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*2
    elif (uInterpolation == CUBIC):
        interfaceBasis.interpolationXi = [oc.BasisInterpolationSpecifications.CUBIC_LAGRANGE]*2
    elif (uInterpolation == HERMITE):
        interfaceBasis.interpolationXi = [oc.BasisInterpolationSpecifications.CUBIC_HERMITE]*2
    else:
        print('Invalid u interpolation')
        exit()
    interfaceBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*2
    interfaceBasis.CreateFinish()
    
if (progressDiagnostics):
    print('Basis functions ... Done')
  
#================================================================================================================================
#  Mesh
#================================================================================================================================

if (progressDiagnostics):
    print('Meshes ...')    

pNodes3D = [0]*8
uNodes3D = [0]*numberOfLocalNodes
localNodes3D = [0]*numberOfLocalNodes

if (problemType != FLUID):
    solidNodes = oc.Nodes()
    solidNodes.CreateStart(solidRegion,numberOfSolidNodes)
    solidNodes.CreateFinish()

    solidMesh = oc.Mesh()
    solidMesh.CreateStart(solidMeshUserNumber,solidRegion,3)
    solidMesh.NumberOfElementsSet(numberOfSolidElements)
    solidMesh.NumberOfComponentsSet(2)

    solidUElements = oc.MeshElements()
    solidUElements.CreateStart(solidMesh,1,uBasis)
    
    solidPElements = oc.MeshElements()
    solidPElements.CreateStart(solidMesh,2,pBasis)
        
    # Solid mesh elements
    if (debugLevel > 2):
        print('  Solid Elements:')
    numberOfSolidXNodes = numberOfSolidXElements*(numberOfNodesXi-1)+1
    for zElementIdx in range(1,numberOfSolidZElements+1):
        for yElementIdx in range(1,numberOfSolidYElements+1):
            for xElementIdx in range(1,numberOfSolidXElements+1):
                elementNumber = xElementIdx+(yElementIdx-1)*numberOfSolidXElements+\
                                (zElementIdx-1)*numberOfSolidXElements*numberOfSolidYElements
                localNodes3D[localNodeIdx000] = (xElementIdx-1)*(numberOfNodesXi-1)+1+ \
                                              (yElementIdx-1)*(numberOfNodesXi-1)*numberOfSolidXNodes+\
                                              (zElementIdx-1)*(numberOfNodesXi-1)*numberOfSolidNodesPerZ
                [uNodes3D,pNodes3D] = GetElementNodes3D(elementNumber,localNodes3D,numberOfSolidXNodes,numberOfSolidXNodes,\
                                                        numberOfSolidXNodes,numberOfSolidNodesPerZ,0,0,0,0)
                solidUElements.NodesSet(elementNumber,uNodes3D)
                solidPElements.NodesSet(elementNumber,pNodes3D)
 
    solidUElements.CreateFinish()
    solidPElements.CreateFinish()

    solidMesh.CreateFinish()

if (problemType != SOLID):
    fluidNodes = oc.Nodes()
    fluidNodes.CreateStart(fluidRegion,numberOfFluidNodes)
    fluidNodes.CreateFinish()

    fluidMesh = oc.Mesh()
    fluidMesh.CreateStart(fluidMeshUserNumber,fluidRegion,3)
    fluidMesh.NumberOfElementsSet(numberOfFluidElements)
    fluidMesh.NumberOfComponentsSet(2)

    fluidUElements = oc.MeshElements()
    fluidUElements.CreateStart(fluidMesh,1,uBasis)
    
    fluidPElements = oc.MeshElements()
    fluidPElements.CreateStart(fluidMesh,2,pBasis)
                
    # Fluid mesh elements
    if (debugLevel > 2):
        print('  Fluid Elements:')
    # Elements around the obstacle
    for zElementIdx in range(1,numberOfSolidZElements+1):
        # Elements in front of the solid
        for yElementIdx in range(1,numberOfFluidY1Elements+1):
            for xElementIdx in range(1,numberOfFluidXElements1+1):
                elementNumber = xElementIdx+(yElementIdx-1)*numberOfFluidXElements1+(zElementIdx-1)*numberOfFluidElementsPerZ1
                localNodes3D[localNodeIdx000] = (xElementIdx-1)*(numberOfNodesXi-1)+1+ \
                                                (yElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidXNodes1+\
                                                (zElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidNodesPerZ1
                [uNodes3D,pNodes3D] = GetElementNodes3D(elementNumber,localNodes3D,numberOfFluidXNodes1,numberOfFluidXNodes1,\
                                                        numberOfFluidXNodes1,numberOfFluidNodesPerZ1,0,0,0,0)
                fluidUElements.NodesSet(elementNumber,uNodes3D)
                fluidPElements.NodesSet(elementNumber,pNodes3D)
        # Elements on the side of the solid
        elementOffset = numberOfFluidXElements1*numberOfFluidY1Elements
        nodeOffset = numberOfFluidXNodes1*(numberOfFluidY1Nodes-1)
        for yElementIdx in range(1,numberOfSolidYElements+1):
            # Elements to the left of the solid
            if (zElementIdx == numberOfSolidZElements):
                zOffset=(yElementIdx-2)*(numberOfNodesXi-1)*(numberOfSolidXNodes-2)+numberOfSolidXNodes-2
                zDelta1=numberOfSolidXNodes-2
                zDelta2=numberOfSolidXNodes-2
                zDelta3=numberOfSolidXNodes-2
            else:
                zOffset=0
                zDelta1=0
                zDelta2=0
                zDelta3=0
            for xElementIdx in range(1,numberOfFluidX1Elements+1):
                elementNumber = xElementIdx+elementOffset+(yElementIdx-1)*numberOfFluidXElements2+\
                                (zElementIdx-1)*numberOfFluidElementsPerZ1
                if (yElementIdx == 1): 
                    localNodes3D[localNodeIdx000] = (xElementIdx-1)*(numberOfNodesXi-1)+1+nodeOffset+\
                                                  (zElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidNodesPerZ1
                    [uNodes3D,pNodes3D] = GetElementNodes3D(elementNumber,localNodes3D,numberOfFluidXNodes1,numberOfFluidXNodes2,\
                                                            numberOfFluidXNodes2,numberOfFluidNodesPerZ1,0,0,zDelta2,zDelta3)
                else:
                    localNodes3D[localNodeIdx000] = (xElementIdx-1)*(numberOfNodesXi-1)+1+nodeOffset+ \
                                                  numberOfFluidXNodes1+(numberOfNodesXi-2)*numberOfFluidXNodes2+ \
                                                  (yElementIdx-2)*(numberOfNodesXi-1)*numberOfFluidXNodes2+ \
                                                  (zElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidNodesPerZ1
                    [uNodes3D,pNodes3D] = GetElementNodes3D(elementNumber,localNodes3D,numberOfFluidXNodes2,numberOfFluidXNodes2,\
                                                            numberOfFluidXNodes2,numberOfFluidNodesPerZ1,zOffset,zDelta1,zDelta2,zDelta3)
                fluidUElements.NodesSet(elementNumber,uNodes3D)
                fluidPElements.NodesSet(elementNumber,pNodes3D)
            # Elements to the right of the solid
            if (zElementIdx == numberOfSolidZElements):
                zOffset=(yElementIdx-1)*(numberOfNodesXi-1)*(numberOfSolidXNodes-2)
                zDelta1=numberOfSolidXNodes-2
                zDelta2=numberOfSolidXNodes-2
                zDelta3=numberOfSolidXNodes-2
            else:
               zOffset=0
               zDelta1=0
               zDelta2=0
               zDelta3=0
            for xElementIdx in range(1,numberOfFluidX2Elements+1):
                elementNumber = xElementIdx+elementOffset+numberOfFluidX1Elements+\
                                (yElementIdx-1)*numberOfFluidXElements2+\
                                (zElementIdx-1)*numberOfFluidElementsPerZ1
                localNodes3D[localNodeIdx000] = (xElementIdx-1)*(numberOfNodesXi-1)+\
                                                (numberOfFluidX1Elements+numberOfSolidXElements)*(numberOfNodesXi-1)+1+nodeOffset+ \
                                                (yElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidXNodes2+ \
                                                (zElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidNodesPerZ1
                if (yElementIdx == numberOfSolidYElements):
                    [uNodes3D,pNodes3D] = GetElementNodes3D(elementNumber,localNodes3D,numberOfFluidXNodes2,numberOfFluidXNodes2,\
                                                            numberOfFluidXNodes1,numberOfFluidNodesPerZ1,zOffset,zDelta1,zDelta2,0)
                else:
                    [uNodes3D,pNodes3D] = GetElementNodes3D(elementNumber,localNodes3D,numberOfFluidXNodes2,numberOfFluidXNodes2,\
                                                            numberOfFluidXNodes2,numberOfFluidNodesPerZ1,zOffset,zDelta1,zDelta2,zDelta3)
                fluidUElements.NodesSet(elementNumber,uNodes3D)
                fluidPElements.NodesSet(elementNumber,pNodes3D)                     
        # Elements behind the solid
        if (zElementIdx == numberOfSolidZElements):
            zOffset = (numberOfSolidXNodes-2)*(numberOfSolidYNodes-2)
        else:
            zOffset = 0
        elementOffset = numberOfFluidXElements1*numberOfFluidY1Elements+numberOfFluidXElements2*numberOfSolidYElements        
        nodeOffset = numberOfFluidXNodes1*numberOfFluidY1Nodes+numberOfFluidXNodes2*(numberOfSolidYNodes-2)
        for yElementIdx in range(1,numberOfFluidY2Elements+1):
            for xElementIdx in range(1,numberOfFluidX1Elements+numberOfFluidX2Elements+numberOfSolidXElements+1):
                elementNumber = xElementIdx+elementOffset+(yElementIdx-1)*(numberOfFluidXElements1)+\
                                (zElementIdx-1)*numberOfFluidElementsPerZ1
                localNodes3D[localNodeIdx000] = (xElementIdx-1)*(numberOfNodesXi-1)+1+nodeOffset+ \
                                              (yElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidXNodes1+\
                                              (zElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidNodesPerZ1
                [uNodes3D,pNodes3D] = GetElementNodes3D(elementNumber,localNodes3D,numberOfFluidXNodes1,numberOfFluidXNodes1,\
                                                        numberOfFluidXNodes1,numberOfFluidNodesPerZ1,zOffset,0,0,0)
                fluidUElements.NodesSet(elementNumber,uNodes3D)
                fluidPElements.NodesSet(elementNumber,pNodes3D)
    # Elements above the obstacle
    elementOffset = numberOfFluidElementsPerZ1*numberOfSolidZElements
    nodeOffset = numberOfFluidNodesPerZ1*numberOfSolidZElements*(numberOfNodesXi-1)
    for zElementIdx in range(1,numberOfFluidZElements+1):
        for yElementIdx in range(1,numberOfFluidYElements1+1):
            for xElementIdx in range(1,numberOfFluidXElements1+1):
                elementNumber = xElementIdx+elementOffset+(yElementIdx-1)*numberOfFluidXElements1+\
                                (zElementIdx-1)*numberOfFluidElementsPerZ2
                localNodes3D[localNodeIdx000] = (xElementIdx-1)*(numberOfNodesXi-1)+1+nodeOffset+ \
                                              (yElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidXNodes1+\
                                              (zElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidNodesPerZ2
                [uNodes3D,pNodes3D] = GetElementNodes3D(elementNumber,localNodes3D,numberOfFluidXNodes1,numberOfFluidXNodes1,\
                                                        numberOfFluidXNodes1,numberOfFluidNodesPerZ2,0,0,0,0)
                fluidUElements.NodesSet(elementNumber,uNodes3D)
                fluidPElements.NodesSet(elementNumber,pNodes3D)
    fluidUElements.CreateFinish()
    fluidPElements.CreateFinish()

    fluidMesh.CreateFinish()

if (progressDiagnostics):
    print('Meshes ... Done')    

#================================================================================================================================
#  Interface
#================================================================================================================================

if (problemType == FSI):
    if (progressDiagnostics):
        print('Interface ...')
    
        # Create an interface between the two meshes
        interface = oc.Interface()
        interface.CreateStart(interfaceUserNumber,worldRegion)
        interface.LabelSet('Interface')
        # Add in the two meshes
        solidMeshIndex = interface.MeshAdd(solidMesh)
        fluidMeshIndex = interface.MeshAdd(fluidMesh)
        interface.CoordinateSystemSet(interfaceCoordinateSystem)
        interface.CreateFinish()
        
    if (progressDiagnostics):
        print('Interface ... Done')
            
#================================================================================================================================
#  Interface Mesh
#================================================================================================================================

if (problemType == FSI):
    if (progressDiagnostics):
        print('Interface Mesh ...')
    
    pNodes2D = [0]*4
    uNodes2D = [0]*numberOfLocalInterfaceNodes
    localNodes2D = [0]*numberOfLocalInterfaceNodes

    # Create an interface mesh
    InterfaceNodes = oc.Nodes()
    InterfaceNodes.CreateStartInterface(interface,numberOfInterfaceNodes)
    InterfaceNodes.CreateFinish()
    
    interfaceMesh = oc.Mesh()
    interfaceMesh.CreateStartInterface(interfaceMeshUserNumber,interface,2)
    interfaceMesh.NumberOfElementsSet(numberOfInterfaceElements)
    interfaceMesh.NumberOfComponentsSet(1)
    
    interfaceElements = oc.MeshElements()
    interfaceElements.CreateStart(interfaceMesh,1,interfaceBasis)
        
    if (debugLevel > 2):
        print('  Interface Elements:')
    elementNumber = 0
    # Do elements on the sides on the obstacle
    for zElementIdx in range(1,numberOfSolidZElements+1):
        for xElementIdx in range(1,numberOfSolidXElements+1):
            elementNumber = elementNumber + 1
            localNodes2D[localNodeIdx000] = xElementIdx*(numberOfNodesXi-1)+(zElementIdx-1)*numberOfInterfaceNodesPerZ
            [uNodes2D,pNodes2D] = GetElementNodes2D(elementNumber,localNodes2D,numberOfFluidXNodes1,numberOfFluidXNodes1,\
                                                    numberOfFluidXNodes1,numberOfFluidNodesPerZ2,0,0,0,0)
            if (useHermite):
                interfaceHermiteElements.NodesSet(elementNumber,[localNodes1,localNodes3])
                if (debugLevel > 2):
                    print('    Element %8d; Nodes: %8d, %8d' % (elementNumber,localNodes1,localNodes3))
            else:
                localNodes2 = localNodes1 + 1
                interfaceQuadraticElements.NodesSet(elementNumber,[localNodes1,localNodes2,localNodes3])
                if (debugLevel > 2):
                    print('    Element %8d; Nodes: %8d, %8d, %8d' % (elementNumber,localNodes1,localNodes2,localNodes3))
                    
    if (useHermite):
        interfaceHermiteElements.CreateFinish()
    else:
        interfaceQuadraticElements.CreateFinish()
                
    interfaceMesh.CreateFinish()

    if (progressDiagnostics):
        print('Interface Mesh ... Done')
    

#================================================================================================================================
#  Mesh Connectivity
#================================================================================================================================

if (problemType == FSI):
    if (progressDiagnostics):
        print('Interface Mesh Connectivity ...')

    # Couple the interface meshes
    interfaceMeshConnectivity = oc.InterfaceMeshConnectivity()
    interfaceMeshConnectivity.CreateStart(interface,interfaceMesh)
    if (useHermite):
        interfaceMeshConnectivity.BasisSet(interfaceHermiteBasis)
    else:
        interfaceMeshConnectivity.BasisSet(interfaceQuadraticBasis)
        
    interfaceElementNumber = 0
    interfaceNodes = [0]*(numberOfInterfaceNodes)
    solidNodes = [0]*(numberOfInterfaceNodes)
    fluidNodes = [0]*(numberOfInterfaceNodes)
    localInterfaceNodes = [0]*numberOfNodesXi
    localSolidNodes = [0]*numberOfNodesXi
    localFluidNodes = [0]*numberOfNodesXi
    # Left edge of solid
    for interfaceElementIdx in range(1,numberOfSolidYElements+1):
        interfaceElementNumber = interfaceElementNumber + 1
        if (debugLevel > 2):
            print('  Interface Element %8d:' % (interfaceElementNumber))        
        solidElementNumber = (interfaceElementIdx - 1)*numberOfSolidXElements + 1
        fluidElementNumber = numberOfFluidX1Elements+(interfaceElementIdx - 1)*(numberOfFluidX1Elements + numberOfFluidX2Elements)
        # Map interface elements
        interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,solidMeshIndex,solidElementNumber)
        interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber)
        if (debugLevel > 2):
            print('    Solid Element %8d; Fluid Element %8d' % (solidElementNumber,fluidElementNumber))        
        localInterfaceNodes[0] = (interfaceElementIdx - 1)*(numberOfNodesXi - 1) + 1
        localSolidNodes[0] = (interfaceElementIdx - 1)*(numberOfNodesXi - 1)*(numberOfSolidXElements*(numberOfNodesXi - 1) + 1)+1
        localFluidNodes[0] = numberOfFluidX1Elements*(numberOfNodesXi-1) + 1 + \
                             (interfaceElementIdx - 1)*(numberOfNodesXi - 1)*\
                             ((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi - 1) + 2)
        if (not useHermite):
            localInterfaceNodes[1] = localInterfaceNodes[0]+1
            localSolidNodes[1] = localSolidNodes[0] + numberOfSolidXElements*(numberOfNodesXi-1)+1
            localFluidNodes[1] = localFluidNodes[0] + ((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)
        localInterfaceNodes[numberOfNodesXi-1] = localInterfaceNodes[0] + (numberOfNodesXi - 1)
        localSolidNodes[numberOfNodesXi-1] = localSolidNodes[0] + \
                                             (numberOfNodesXi - 1)*(numberOfSolidXElements*(numberOfNodesXi - 1) + 1)
        localFluidNodes[numberOfNodesXi-1] = localFluidNodes[0] + ((numberOfFluidX1Elements+numberOfFluidX2Elements)*\
                                                                   (numberOfNodesXi-1)+2)*(numberOfNodesXi - 1)
        # Map interface xi
        for localNodeIdx in range(0,numberOfNodesXi):
            xi=float(localNodeIdx)/float(numberOfNodesXi-1)
            solidXi = [0.0,xi]
            fluidXi = [1.0,xi]
            interfaceNodes[localInterfaceNodes[localNodeIdx]-1]=localInterfaceNodes[localNodeIdx]
            solidNodes[localInterfaceNodes[localNodeIdx]-1]=localSolidNodes[localNodeIdx]
            fluidNodes[localInterfaceNodes[localNodeIdx]-1]=localFluidNodes[localNodeIdx]
            interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,solidMeshIndex,solidElementNumber,localNodeIdx+1,1,solidXi)
            interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber,localNodeIdx+1,1,fluidXi)
            if (debugLevel > 2):
                print('    Local node    %8d:' % (localNodeIdx+1))        
                print('      Interface node    %8d:' % (localInterfaceNodes[localNodeIdx]))        
                print('      Solid node        %8d; Solid xi = [%.2f, %.2f ]' % (localSolidNodes[localNodeIdx],solidXi[0],solidXi[1]))
                print('      Fluid node        %8d; Fluid xi = [%.2f, %.2f ]' % (localFluidNodes[localNodeIdx],fluidXi[0],fluidXi[1]))
    # Top edge of solid
    for interfaceElementIdx in range(1,numberOfSolidXElements+1):
        interfaceElementNumber = interfaceElementNumber + 1
        if (debugLevel > 2):
            print('  Interface Element %8d:' % (interfaceElementNumber))        
        solidElementNumber = interfaceElementIdx + numberOfSolidXElements*(numberOfSolidYElements - 1)
        fluidElementNumber = interfaceElementIdx + numberOfFluidX1Elements + \
                             (numberOfFluidX1Elements+numberOfFluidX2Elements)*numberOfSolidYElements
        # Map interface elements
        interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,solidMeshIndex,solidElementNumber)
        interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber)
        if (debugLevel > 2):
            print('    Solid Element %8d; Fluid Element %8d' % (solidElementNumber,fluidElementNumber))        
        localInterfaceNodes[0] = (interfaceElementIdx - 1)*(numberOfNodesXi - 1) + 1 + (numberOfSolidYElements)*(numberOfNodesXi-1)
        localSolidNodes[0] = (interfaceElementIdx - 1)*(numberOfNodesXi - 1) + 1 + \
                                 (numberOfSolidXElements*(numberOfNodesXi-1)+1)*(numberOfSolidYElements)*(numberOfNodesXi-1)        
        localFluidNodes[0] = (interfaceElementIdx - 1)*(numberOfNodesXi - 1) + 1 + numberOfFluidX1Elements*(numberOfNodesXi-1) + \
                             ((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi - 1) + 2)*\
                             numberOfSolidYElements*(numberOfNodesXi-1)
        if (not useHermite):
            localInterfaceNodes[1] = localInterfaceNodes[0]+1
            localSolidNodes[1] = localSolidNodes[0] + 1
            localFluidNodes[1] = localFluidNodes[0] + 1
        localInterfaceNodes[numberOfNodesXi-1] = localInterfaceNodes[0] + (numberOfNodesXi - 1)
        localSolidNodes[numberOfNodesXi-1] = localSolidNodes[0] + (numberOfNodesXi - 1)
        localFluidNodes[numberOfNodesXi-1] = localFluidNodes[0] + (numberOfNodesXi - 1)        
        # Map interface xi
        for localNodeIdx in range(0,numberOfNodesXi):
            xi=float(localNodeIdx)/float(numberOfNodesXi-1)
            solidXi = [xi,1.0]
            fluidXi = [xi,0.0]
            interfaceNodes[localInterfaceNodes[localNodeIdx]-1]=localInterfaceNodes[localNodeIdx]
            solidNodes[localInterfaceNodes[localNodeIdx]-1]=localSolidNodes[localNodeIdx]
            fluidNodes[localInterfaceNodes[localNodeIdx]-1]=localFluidNodes[localNodeIdx]
            interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,solidMeshIndex,solidElementNumber,localNodeIdx+1,1,solidXi)
            interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber,localNodeIdx+1,1,fluidXi)
            if (debugLevel > 2):
                print('    Local node    %8d:' % (localNodeIdx+1))        
                print('      Interface node    %8d:' % (localInterfaceNodes[localNodeIdx]))        
                print('      Solid node        %8d; Solid xi = [%.2f, %.2f ]' % (localSolidNodes[localNodeIdx],solidXi[0],solidXi[1]))
                print('      Fluid node        %8d; Fluid xi = [%.2f, %.2f ]' % (localFluidNodes[localNodeIdx],fluidXi[0],fluidXi[1]))
    # right edge of solid
    for interfaceElementIdx in range(1,numberOfSolidYElements+1):
        if (interfaceElementIdx == 1):
            offset = numberOfSolidXElements*(numberOfNodesXi-1) 
        else:
            offset = 1
        interfaceElementNumber = interfaceElementNumber + 1
        if (debugLevel > 2):
            print('  Interface Element %8d:' % (interfaceElementNumber))
        solidElementNumber = (numberOfSolidYElements - interfaceElementIdx + 1)*numberOfSolidXElements 
        fluidElementNumber = (numberOfSolidYElements - interfaceElementIdx)*(numberOfFluidX1Elements + numberOfFluidX2Elements) \
                             + numberOfFluidX1Elements + 1
        # Map interface elements
        interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,solidMeshIndex,solidElementNumber)
        interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber)
        if (debugLevel > 2):
            print('    Solid Element %8d; Fluid Element %8d' % (solidElementNumber,fluidElementNumber))        
        localInterfaceNodes[0] = (interfaceElementIdx - 1)*(numberOfNodesXi - 1) + \
                                 (numberOfSolidXElements + numberOfSolidYElements)*(numberOfNodesXi - 1) + 1
        localSolidNodes[0] = (numberOfSolidXElements*(numberOfNodesXi - 1) + 1)* \
                             ((numberOfSolidYElements - interfaceElementIdx + 1)*(numberOfNodesXi-1)+1)
        localFluidNodes[0] = ((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi - 1) + 2)*\
                             (numberOfSolidYElements - interfaceElementIdx + 1)*(numberOfNodesXi - 1) + \
                             numberOfFluidX1Elements*(numberOfNodesXi-1)+ 1 + offset
        if (not useHermite):
            localInterfaceNodes[1] = localInterfaceNodes[0]+1
            localSolidNodes[1] = localSolidNodes[0] - numberOfSolidXElements*(numberOfNodesXi-1) - 1
            localFluidNodes[1] = localFluidNodes[0] - ((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2) \
                                 - offset + 1
        localInterfaceNodes[numberOfNodesXi-1] = localInterfaceNodes[0] + numberOfNodesXi - 1
        localSolidNodes[numberOfNodesXi-1] = localSolidNodes[0] - \
                                             (numberOfNodesXi - 1)*(numberOfSolidXElements*(numberOfNodesXi - 1) + 1)
        localFluidNodes[numberOfNodesXi-1] = localFluidNodes[0] - \
                                             ((numberOfFluidX1Elements+numberOfFluidX2Elements)*\
                                              (numberOfNodesXi-1)+2)*(numberOfNodesXi - 1) -offset+1
        # Map interface xi
        for localNodeIdx in range(0,numberOfNodesXi):
            xi=float(numberOfNodesXi-localNodeIdx-1)/float(numberOfNodesXi-1)
            solidXi = [1.0,xi]
            fluidXi = [0.0,xi]
            interfaceNodes[localInterfaceNodes[localNodeIdx]-1]=localInterfaceNodes[localNodeIdx]
            solidNodes[localInterfaceNodes[localNodeIdx]-1]=localSolidNodes[localNodeIdx]
            fluidNodes[localInterfaceNodes[localNodeIdx]-1]=localFluidNodes[localNodeIdx]
            interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,solidMeshIndex,solidElementNumber,localNodeIdx+1,1,solidXi)
            interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber,localNodeIdx+1,1,fluidXi)
            if (debugLevel > 2):
                print('    Local node    %8d:' % (localNodeIdx+1))        
                print('      Interface node    %8d:' % (localInterfaceNodes[localNodeIdx]))        
                print('      Solid node        %8d; Solid xi = [%.2f, %.2f ]' % (localSolidNodes[localNodeIdx],solidXi[0],solidXi[1]))
                print('      Fluid node        %8d; Fluid xi = [%.2f, %.2f ]' % (localFluidNodes[localNodeIdx],fluidXi[0],fluidXi[1]))
    # Map interface nodes
    interfaceMeshConnectivity.NodeNumberSet(interfaceNodes,solidMeshIndex,solidNodes,fluidMeshIndex,fluidNodes)        

    interfaceMeshConnectivity.CreateFinish()

    if (progressDiagnostics):
        print('Interface Mesh Connectivity ... Done')

#================================================================================================================================
#  Decomposition
#================================================================================================================================

if (progressDiagnostics):
    print('Decomposition ...')
    
if (problemType != FLUID):
    # Create a decomposition for the solid mesh
    solidDecomposition = oc.Decomposition()
    solidDecomposition.CreateStart(solidDecompositionUserNumber,solidMesh)
    solidDecomposition.TypeSet(oc.DecompositionTypes.CALCULATED)
    solidDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
    solidDecomposition.CalculateFacesSet(True)
    solidDecomposition.CreateFinish()

if (problemType != SOLID):
    # Create a decomposition for the fluid mesh
    fluidDecomposition = oc.Decomposition()
    fluidDecomposition.CreateStart(fluidDecompositionUserNumber,fluidMesh)
    fluidDecomposition.TypeSet(oc.DecompositionTypes.CALCULATED)
    fluidDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
    fluidDecomposition.CalculateFacesSet(True)
    fluidDecomposition.CreateFinish()

if (problemType == FSI):
    # Create a decomposition for the interface mesh
    interfaceDecomposition = oc.Decomposition()
    interfaceDecomposition.CreateStart(interfaceDecompositionUserNumber,interfaceMesh)
    interfaceDecomposition.TypeSet(oc.DecompositionTypes.CALCULATED)
    interfaceDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
    interfaceDecomposition.CreateFinish()

if (progressDiagnostics):
    print('Decomposition ... Done')
    
#================================================================================================================================
#  Geometric Field
#================================================================================================================================

if (progressDiagnostics):
    print('Geometric Field ...')

if (problemType != FLUID):    
    # Start to create a default (geometric) field on the solid region
    solidGeometricField = oc.Field()
    solidGeometricField.CreateStart(solidGeometricFieldUserNumber,solidRegion)
    # Set the decomposition to use
    solidGeometricField.DecompositionSet(solidDecomposition)
    # Set the scaling to use
    if (uInterpolation == HERMITE):
        solidGeometricField.ScalingTypeSet(oc.FieldScalingTypes.ARITHMETIC_MEAN)
    else:
        solidGeometricField.ScalingTypeSet(oc.FieldScalingTypes.NONE)
    solidGeometricField.VariableLabelSet(oc.FieldVariableTypes.U,'SolidGeometry')
    # Set the domain to be used by the field components.
    solidGeometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,1,1)
    solidGeometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,2,1)
    solidGeometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,3,1)
    # Finish creating the first field
    solidGeometricField.CreateFinish()

if (problemType != SOLID):
    # Start to create a default (geometric) field on the fluid region
    fluidGeometricField = oc.Field()
    fluidGeometricField.CreateStart(fluidGeometricFieldUserNumber,fluidRegion)
    # Set the decomposition to use
    if (uInterpolation == HERMITE):
        fluidGeometricField.ScalingTypeSet(oc.FieldScalingTypes.ARITHMETIC_MEAN)
    else:
        fluidGeometricField.DecompositionSet(fluidDecomposition)
    # Set the scaling to use
    fluidGeometricField.ScalingTypeSet(oc.FieldScalingTypes.NONE)
    fluidGeometricField.VariableLabelSet(oc.FieldVariableTypes.U,'FluidGeometry')
    # Set the domain to be used by the field components.
    fluidGeometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,1,1)
    fluidGeometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,2,1)
    fluidGeometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,3,1)
    # Finish creating the second field
    fluidGeometricField.CreateFinish()

if (problemType == FSI):
    # Start to create a default (geometric) field on the Interface
    interfaceGeometricField = oc.Field()
    interfaceGeometricField.CreateStartInterface(interfaceGeometricFieldUserNumber,interface)
    # Set the decomposition to use
    interfaceGeometricField.DecompositionSet(interfaceDecomposition)
    # Set the scaling to use
    if (uInterpolation == HERMITE):
        interfaceGeometricField.ScalingTypeSet(oc.FieldScalingTypes.ARITHMETIC_MEAN)
    else:
        interfaceGeometricField.ScalingTypeSet(oc.FieldScalingTypes.NONE)
    interfaceGeometricField.VariableLabelSet(oc.FieldVariableTypes.U,'InterfaceGeometry')
    # Set the domain to be used by the field components.
    interfaceGeometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,1,1)
    interfaceGeometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,2,1)
    interfaceGeometricField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,3,1)
    # Finish creating the first field
    interfaceGeometricField.CreateFinish()

if (progressDiagnostics):
    print('Geometric Field ... Done')
    
if (progressDiagnostics):
    print('Geometric Parameters ...')
    
if (problemType != FLUID):
    # Solid nodes
    if (debugLevel > 2):
        print('  Solid Nodes:')
    for zNodeIdx in range(1,numberOfSolidZNodes+1):
        for yNodeIdx in range(1,numberOfSolidYNodes+1):
            for xNodeIdx in range(1,numberOfSolidXNodes+1):
                nodeNumber = xNodeIdx+(yNodeIdx-1)*numberOfSolidXNodes+(zNodeIdx-1)*numberOfSolidNodesPerZ
                nodeDomain = solidDecomposition.NodeDomainGet(1,nodeNumber)
                if (nodeDomain == computationalNodeNumber):
                    xPosition = fluidX1Size + float(xNodeIdx-1)/float(numberOfSolidXNodes-1)*solidXSize
                    yPosition = fluidY1Size + float(yNodeIdx-1)/float(numberOfSolidYNodes-1)*solidYSize
                    zPosition = float(ZNodeIdx-1)/float(numberOfSolidZNodes-1)*solidZSize
                    SetNodeParameters3D(nodeNumber,solidGeometricField,xPosition,yPosition,zPosition)
    # Update fields            
    solidGeometricField.ParameterSetUpdateStart(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)
    solidGeometricField.ParameterSetUpdateFinish(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)

if (problemType != SOLID):                        
    if (debugLevel > 2):
        print('  Fluid Nodes:')
    # Nodes around of the solid
    for zNodeIdx in range(1,numberOfSolidZNodes+1):
        # Nodes in front of the solid
        for yNodeIdx in range(1,numberOfFluidY1Nodes+1):
            for xNodeIdx in range(1,numberOfFluidXNodes1+1):
                nodeNumber = xNodeIdx+(yNodeIdx-1)*numberOfFluidXNodes1+(zNodeIdx-1)*numberOfFluidNodesPerZ1
                nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
                if (nodeDomain == computationalNodeNumber):
                    xPosition = float(xNodeIdx-1)/float(numberOfFluidXNodes1-1)*(fluidX1Size+solidXSize+fluidX2Size)
                    yPosition = float(yNodeIdx-1)/float(numberOfFluidY1Nodes-1)*fluidY1Size
                    zPosition = float(zNodeIdx-1)/float(numberOfSolidZNodes-1)*solidZSize
                    SetNodeParameters3D(nodeNumber,fluidGeometricField,xPosition,yPosition,zPosition)
        # Nodes on the side of the solid
        for yNodeIdx in range(1,numberOfSolidYNodes-1):
            # Nodes to the left of the solid
            nodeOffset = numberOfFluidXNodes1*numberOfFluidY1Nodes
            for xNodeIdx in range(1,numberOfFluidX1Nodes+1):
                nodeNumber = xNodeIdx+nodeOffset+(yNodeIdx-1)*numberOfFluidXNodes2+(zNodeIdx-1)*numberOfFluidNodesPerZ1
                nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
                if (nodeDomain == computationalNodeNumber):
                    xPosition = float(xNodeIdx-1)/float(numberOfFluidX1Nodes-1)*fluidX1Size
                    yPosition = fluidY1Size+float(yNodeIdx)/float(numberOfSolidYNodes-1)*solidYSize
                    zPosition = float(zNodeIdx-1)/float(numberOfSolidZNodes-1)*solidZSize
                    SetNodeParameters3D(nodeNumber,fluidGeometricField,xPosition,yPosition,zPosition)
            # Nodes to the right of the solid
            nodeOffset = numberOfFluidX1Nodes+numberOfFluidXNodes1*numberOfFluidY1Nodes
            for xNodeIdx in range(1,numberOfFluidX2Nodes+1):
                nodeNumber = xNodeIdx+nodeOffset+(yNodeIdx-1)*numberOfFluidXNodes2+(zNodeIdx-1)*numberOfFluidNodesPerZ1
                nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
                if (nodeDomain == computationalNodeNumber):
                    xPosition = fluidX1Size+solidXSize+float(xNodeIdx-1)/float(numberOfFluidX2Nodes-1)*fluidX2Size
                    yPosition = fluidY1Size+float(yNodeIdx)/float(numberOfSolidYNodes-1)*solidYSize
                    zPosition = float(zNodeIdx-1)/float(numberOfSolidZNodes-1)*solidZSize
                    SetNodeParameters3D(nodeNumber,fluidGeometricField,xPosition,yPosition,zPosition)
        # Nodes behind the solid
        nodeOffset = numberOfFluidXNodes1*numberOfFluidY1Nodes+numberOfFluidXNodes2*(numberOfSolidYNodes-2)
        for yNodeIdx in range(1,numberOfFluidY2Nodes+1):
            for xNodeIdx in range(1,numberOfFluidXNodes1+1):
                nodeNumber = xNodeIdx+nodeOffset+(yNodeIdx-1)*numberOfFluidXNodes1+(zNodeIdx-1)*numberOfFluidNodesPerZ1
                nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
                if (nodeDomain == computationalNodeNumber):
                    xPosition = float(xNodeIdx-1)/float(numberOfFluidXNodes1-1)*(fluidX1Size+solidXSize+fluidX2Size)
                    yPosition = fluidY1Size+solidYSize+float(yNodeIdx-1)/float(numberOfFluidY2Nodes-1)*fluidY2Size
                    zPosition = float(zNodeIdx-1)/float(numberOfSolidZNodes-1)*solidZSize
                    SetNodeParameters3D(nodeNumber,fluidGeometricField,xPosition,yPosition,zPosition)
    # Nodes to the top of the solid
    nodeOffset = numberOfFluidNodesPerZ1*(numberOfSolidZNodes-1)
    for zNodeIdx in range(1,numberOfFluidZNodes+1):
        for yNodeIdx in range(1,numberOfFluidYNodes1+1):
            for xNodeIdx in range(1,numberOfFluidXNodes1+1):
                nodeNumber = xNodeIdx+nodeOffset+(yNodeIdx-1)*numberOfFluidXNodes1+(zNodeIdx-1)*numberOfFluidNodesPerZ2
                nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
                if (nodeDomain == computationalNodeNumber):
                    xPosition = float(xNodeIdx-1)/float(numberOfFluidXNodes1-1)*(fluidX1Size+solidXSize+fluidX2Size)
                    yPosition = float(yNodeIdx-1)/float(numberOfFluidYNodes1-1)*(fluidY1Size+solidYSize+fluidY2Size)
                    zPosition = solidZSize + float(zNodeIdx-1)/float(numberOfFluidZNodes-1)*fluidZSize
                    SetNodeParameters3D(nodeNumber,fluidGeometricField,xPosition,yPosition,zPosition)
    # Update fields            
    fluidGeometricField.ParameterSetUpdateStart(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)
    fluidGeometricField.ParameterSetUpdateFinish(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)

if (problemType == FSI):
    if (debugLevel > 2):
        print('  Interface Nodes:')
    # Left edge of interface nodes    
    for yNodeIdx in range(1,numberOfSolidYElements*(numberOfNodesXi-1)+1):
        nodeNumber = yNodeIdx
        #nodeDomain = interfaceDecomposition.NodeDomainGet(1,nodeNumber)
        nodeDomain = computationalNodeNumber
        if (nodeDomain == computationalNodeNumber):
            xPosition = fluidX1Size
            yPosition = float(yNodeIdx-1)/float(numberOfSolidYElements*(numberOfNodesXi-1))*solidYSize
            interfaceGeometricField.ParameterSetUpdateNodeDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                             1,oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,xPosition)
            interfaceGeometricField.ParameterSetUpdateNodeDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                             1,oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,yPosition)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Position         = [ %.2f, %.2f ]' % (xPosition,yPosition))                 
            if (useHermite):
                interfaceGeometricField.ParameterSetUpdateNodeDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                                 1,oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,0.0)
                interfaceGeometricField.ParameterSetUpdateNodeDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                                 1,oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,1.0)
                if (debugLevel > 2):
                    print('        S1 derivative    = [ %.2f, %.2f ]' % (1.0,0.0))                 
    # Top edge of interface nodes    
    for xNodeIdx in range(1,numberOfSolidXElements*(numberOfNodesXi-1)+2):
        nodeNumber = xNodeIdx+numberOfSolidYElements*(numberOfNodesXi-1)
        #nodeDomain = interfaceDecomposition.NodeDomainGet(1,nodeNumber)
        nodeDomain = computationalNodeNumber
        if (nodeDomain == computationalNodeNumber):
            xPosition = fluidX1Size+float(xNodeIdx-1)/float(numberOfSolidXElements*(numberOfNodesXi-1))*solidXSize
            yPosition = solidYSize
            interfaceGeometricField.ParameterSetUpdateNodeDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                             1,oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,xPosition)
            interfaceGeometricField.ParameterSetUpdateNodeDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                             1,oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,yPosition)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Position         = [ %.2f, %.2f ]' % (xPosition,yPosition))                 
            if (useHermite):
                interfaceGeometricField.ParameterSetUpdateNodeDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                                 1,oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,1.0)
                interfaceGeometricField.ParameterSetUpdateNodeDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                                 1,oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,0.0)
                if (debugLevel > 2):
                    print('        S1 derivative    = [ %.2f, %.2f ]' % (1.0,0.0))                 
    # Right edge of interface nodes    
    for yNodeIdx in range(1,numberOfSolidYElements*(numberOfNodesXi-1)+1):
        nodeNumber = yNodeIdx+(numberOfSolidYElements+numberOfSolidXElements)*(numberOfNodesXi-1)+1
        #nodeDomain = interfaceDecomposition.NodeDomainGet(1,nodeNumber)
        nodeDomain = computationalNodeNumber
        if (nodeDomain == computationalNodeNumber):
            xPosition = fluidX1Size+solidXSize
            yPosition = solidYSize-float(yNodeIdx)/float(numberOfSolidYElements*(numberOfNodesXi-1))*solidYSize
            interfaceGeometricField.ParameterSetUpdateNodeDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                             1,oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,xPosition)
            interfaceGeometricField.ParameterSetUpdateNodeDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                             1,oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,yPosition)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Position         = [ %.2f, %.2f ]' % (xPosition,yPosition))                 
            if (useHermite):
                interfaceGeometricField.ParameterSetUpdateNodeDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                                 1,oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,0.0)
                interfaceGeometricField.ParameterSetUpdateNodeDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                                 1,oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,-1.0)
                if (debugLevel > 2):
                    print('        S1 derivative    = [ %.2f, %.2f ]' % (1.0,0.0))                 

    # Update fields            
    interfaceGeometricField.ParameterSetUpdateStart(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)
    interfaceGeometricField.ParameterSetUpdateFinish(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Geometric Parameters ... Done')

# Export results
if (problemType != SOLID):
    fluidFields = oc.Fields()
    fluidFields.CreateRegion(fluidRegion)
    fluidFields.NodesExport("3DObstacleFluid","FORTRAN")
    fluidFields.ElementsExport("3DObstacleFluid","FORTRAN")
    fluidFields.Finalise()
if (problemType != FLUID):    
    solidFields = oc.Fields()
    solidFields.CreateRegion(solidRegion)
    solidFields.NodesExport("3DObstacleSolid","FORTRAN")
    solidFields.ElementsExport("3DObstacleSolid","FORTRAN")
    solidFields.Finalise()
if (problemType == FSI):
    interfaceFields = oc.Fields()
    interfaceFields.CreateInterface(interface)
    interfaceFields.NodesExport("3DObstacleInterface","FORTRAN")
    interfaceFields.ElementsExport("3DObstacleInterface","FORTRAN")
    interfaceFields.Finalise()

#================================================================================================================================
#  Equations Set
#================================================================================================================================

if (progressDiagnostics):
    print('Equations Sets ...')

if (problemType != FLUID):
    # Create the equations set for the solid region 
    solidEquationsSetField = oc.Field()
    solidEquationsSet = oc.EquationsSet()
    solidEquationsSetSpecification = [oc.EquationsSetClasses.ELASTICITY,
                                      oc.EquationsSetTypes.FINITE_ELASTICITY,
                                      oc.EquationsSetSubtypes.MOONEY_RIVLIN]
                                      #oc.EquationsSetSubtypes.MR_AND_GROWTH_LAW_IN_CELLML]
    solidEquationsSet.CreateStart(solidEquationsSetUserNumber,solidRegion,solidGeometricField,
                                  solidEquationsSetSpecification,solidEquationsSetFieldUserNumber,
                                  solidEquationsSetField)
    solidEquationsSet.OutputTypeSet(solidEquationsSetOutputType)
    solidEquationsSet.CreateFinish()
    
if (problemType != SOLID):
    # Create the equations set for the fluid region - ALE Navier-Stokes
    fluidEquationsSetField = oc.Field()
    fluidEquationsSet = oc.EquationsSet()
    if RBS:
        if (problemType == FSI):
            fluidEquationsSetSpecification = [oc.EquationsSetClasses.FLUID_MECHANICS,
                                              oc.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                              oc.EquationsSetSubtypes.ALE_RBS_NAVIER_STOKES]
        else:
            fluidEquationsSetSpecification = [oc.EquationsSetClasses.FLUID_MECHANICS,
                                              oc.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                              oc.EquationsSetSubtypes.TRANSIENT_RBS_NAVIER_STOKES]            
    else:
        if (problemType == FSI):
            fluidEquationsSetSpecification = [oc.EquationsSetClasses.FLUID_MECHANICS,
                                              oc.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                              oc.EquationsSetSubtypes.ALE_NAVIER_STOKES]
        else:
            fluidEquationsSetSpecification = [oc.EquationsSetClasses.FLUID_MECHANICS,
                                              oc.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                              oc.EquationsSetSubtypes.TRANSIENT_NAVIER_STOKES]
        
    fluidEquationsSet.CreateStart(fluidEquationsSetUserNumber,fluidRegion,fluidGeometricField,
                                  fluidEquationsSetSpecification,fluidEquationsSetFieldUserNumber,
                                  fluidEquationsSetField)
    fluidEquationsSet.OutputTypeSet(fluidEquationsSetOutputType)
    fluidEquationsSet.CreateFinish()

    if RBS:
        # Set boundary retrograde flow stabilisation scaling factor (default 0- do not use)
        fluidEquationsSetField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U1,
                                                           oc.FieldParameterSetTypes.VALUES,1,1.0)
        # Set max CFL number (default 1.0)
        fluidEquationsSetField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U1,
                                                           oc.FieldParameterSetTypes.VALUES,2,1.0E20)
        # Set time increment (default 0.0)
        fluidEquationsSetField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U1,
                                                           oc.FieldParameterSetTypes.VALUES,3,timeStep)
        # Set stabilisation type (default 1.0 = RBS)
        fluidEquationsSetField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U1,
                                                           oc.FieldParameterSetTypes.VALUES,4,1.0)
        
if (problemType == FSI):
    # Create the equations set for the moving mesh
    movingMeshEquationsSetField = oc.Field()
    movingMeshEquationsSet = oc.EquationsSet()
    movingMeshEquationsSetSpecification = [oc.EquationsSetClasses.CLASSICAL_FIELD,
                                           oc.EquationsSetTypes.LAPLACE_EQUATION,
                                           oc.EquationsSetSubtypes.MOVING_MESH_LAPLACE]
    movingMeshEquationsSet.CreateStart(movingMeshEquationsSetUserNumber,fluidRegion,fluidGeometricField,
                                       movingMeshEquationsSetSpecification,movingMeshEquationsSetFieldUserNumber,
                                       movingMeshEquationsSetField)
    movingMeshEquationsSet.OutputTypeSet(movingMeshEquationsSetOutputType)
    movingMeshEquationsSet.CreateFinish()
    
if (progressDiagnostics):
    print('Equations Sets ... Done')


#================================================================================================================================
#  Dependent Field
#================================================================================================================================

if (progressDiagnostics):
    print('Dependent Fields ...')

if (problemType != FLUID):
    # Create the equations set dependent field variables for the solid equations set
    solidDependentField = oc.Field()
    solidEquationsSet.DependentCreateStart(solidDependentFieldUserNumber,solidDependentField)
    solidDependentField.VariableLabelSet(oc.FieldVariableTypes.U,'SolidDependent')
    for componentIdx in range(1,numberOfDimensions+1):
        solidDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,componentIdx,1)
        solidDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.DELUDELN,componentIdx,1)
    solidDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,numberOfDimensions+1,2)
    solidDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.DELUDELN,numberOfDimensions+1,2)
    solidDependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U,numberOfDimensions+1,oc.FieldInterpolationTypes.NODE_BASED)
    solidDependentField.ComponentInterpolationSet(oc.FieldVariableTypes.DELUDELN,numberOfDimensions+1,oc.FieldInterpolationTypes.NODE_BASED)
    if (useHermite):
        solidDependentField.ScalingTypeSet(oc.FieldScalingTypes.ARITHMETIC_MEAN)
    else:
        solidDependentField.ScalingTypeSet(oc.FieldScalingTypes.NONE)
    solidEquationsSet.DependentCreateFinish()

    # Initialise the solid dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
    for componentIdx in range(1,numberOfDimensions+1):
        solidGeometricField.ParametersToFieldParametersComponentCopy(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,\
                                                                     componentIdx,solidDependentField,oc.FieldVariableTypes.U,
                                                                     oc.FieldParameterSetTypes.VALUES,componentIdx)
    solidDependentField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                    numberOfDimensions+1,solidPInit)
    
    solidDependentField.ParameterSetUpdateStart(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)
    solidDependentField.ParameterSetUpdateFinish(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)

if (problemType != FLUID):
    # Create the equations set dependent field variables for dynamic Navier-Stokes
    fluidDependentField = oc.Field()
    fluidEquationsSet.DependentCreateStart(fluidDependentFieldUserNumber,fluidDependentField)
    fluidDependentField.VariableLabelSet(oc.FieldVariableTypes.U,'FluidDependent')
    # Set the mesh component to be used by the field components.
    for componentIdx in range(1,numberOfDimensions+1):
        fluidDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,componentIdx,1)
        fluidDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.DELUDELN,componentIdx,1)
    fluidDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,numberOfDimensions+1,2)
    fluidDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.DELUDELN,numberOfDimensions+1,2)
    # fluidDependentField.ComponentInterpolationSet(oc.FieldVariableTypes.U,numberOfDimensions+1,oc.FieldInterpolationTypes.NODE_BASED)
    # fluidDependentField.ComponentInterpolationSet(oc.FieldVariableTypes.DELUDELN,numberOfDimensions+1,oc.FieldInterpolationTypes.NODE_BASED)
    # Finish the equations set dependent field variables
    fluidEquationsSet.DependentCreateFinish()

    # Initialise the fluid dependent field
    for componentIdx in range(1,numberOfDimensions+1):
        fluidDependentField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,componentIdx,0.0)
    # Initialise pressure component
    fluidDependentField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                    numberOfDimensions+1,fluidPInit)
    if RBS:
        fluidDependentField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.PRESSURE_VALUES,3,fluidPInit)
        
    fluidDependentField.ParameterSetUpdateStart(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)
    fluidDependentField.ParameterSetUpdateFinish(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)

if (problemType == FSI):        
    # Create the equations set dependent field variables for moving mesh
    movingMeshDependentField = oc.Field()
    movingMeshEquationsSet.DependentCreateStart(movingMeshDependentFieldUserNumber,movingMeshDependentField)
    movingMeshDependentField.VariableLabelSet(oc.FieldVariableTypes.U,'MovingMeshDependent')
    # Set the mesh component to be used by the field components.
    for componentIdx in range(1,numberOfDimensions+1):
        movingMeshDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,componentIdx,1)
        movingMeshDependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.DELUDELN,componentIdx,1)
    # Finish the equations set dependent field variables
    movingMeshEquationsSet.DependentCreateFinish()

    # Initialise dependent field moving mesh
    for ComponentIdx in range(1,numberOfDimensions+1):
        movingMeshDependentField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES, \
                                                             componentIdx,0.0)

    movingMeshDependentField.ParameterSetUpdateStart(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)
    movingMeshDependentField.ParameterSetUpdateFinish(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Dependent Fields ... Done')
     
#================================================================================================================================
#  Materials Field
#================================================================================================================================

if (progressDiagnostics):
    print('Materials Fields ...')

if (problemType != FLUID):
    # Create the solid materials field
    solidMaterialsField = oc.Field()
    solidEquationsSet.MaterialsCreateStart(solidMaterialsFieldUserNumber,solidMaterialsField)
    solidMaterialsField.VariableLabelSet(oc.FieldVariableTypes.U,'SolidMaterials')
    solidMaterialsField.VariableLabelSet(oc.FieldVariableTypes.V,'SolidDensity')
    solidEquationsSet.MaterialsCreateFinish()
    # Set Mooney-Rivlin constants c10 and c01 respectively
    solidMaterialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,mooneyRivlin1)
    solidMaterialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2,mooneyRivlin2)
    solidMaterialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.V,oc.FieldParameterSetTypes.VALUES,1,solidDensity)

if (problemType != SOLID):
    # Create the equations set materials field variables for dynamic Navier-Stokes
    fluidMaterialsField = oc.Field()
    fluidEquationsSet.MaterialsCreateStart(fluidMaterialsFieldUserNumber,fluidMaterialsField)
    # Finish the equations set materials field variables
    fluidEquationsSet.MaterialsCreateFinish()
    fluidMaterialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,fluidDynamicViscosity)
    fluidMaterialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,2,fluidDensity)
    
if (problemType == FSI):    
    # Create the equations set materials field variables for moving mesh
    movingMeshMaterialsField = oc.Field()
    movingMeshEquationsSet.MaterialsCreateStart(movingMeshMaterialsFieldUserNumber,movingMeshMaterialsField)
    # Finish the equations set materials field variables
    movingMeshEquationsSet.MaterialsCreateFinish()

    movingMeshMaterialsField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,\
                                                         movingMeshKParameter)
   
if (progressDiagnostics):
    print('Materials Fields ... Done')
    
#================================================================================================================================
#  Source Field
#================================================================================================================================

if (problemType != FLUID):
    if (gravityFlag):
        if (progressDiagnostics):
            print('Source Fields ...')
        #Create the source field with the gravity vector
        soidSourceField = oc.Field()
        solidEquationsSet.SourceCreateStart(solidSourceFieldUserNumber,solidSourceField)
        solidSourceField.ScalingTypeSet(oc.FieldScalingTypes.NONE)
        solidEquationsSet.SourceCreateFinish()
        for componentIdx in range(1,numberOfDimensions+1):
            solidSourceField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,component_Idx,
                                                         gravity[componentIdx-1])
        if (progressDiagnostics):
            print('Source Fields ... Done')
         
#================================================================================================================================
# Independent Field
#================================================================================================================================

if (problemType == FSI):
    if (progressDiagnostics):
        print('Independent Fields ...')

    # Create fluid mesh velocity independent field 
    fluidIndependentField = oc.Field()
    fluidEquationsSet.IndependentCreateStart(fluidIndependentFieldUserNumber,fluidIndependentField)
    fluidIndependentField.VariableLabelSet(oc.FieldVariableTypes.U,'FluidIndependent')
    # Set the mesh component to be used by the field components.
    for componentIdx in range(1,numberOfDimensions+1):
        fluidIndependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,componentIdx,1)
    # Finish the equations set independent field variables
    fluidEquationsSet.IndependentCreateFinish()
  
    # Create the moving mesh independent field 
    movingMeshIndependentField = oc.Field()
    movingMeshEquationsSet.IndependentCreateStart(movingMeshIndependentFieldUserNumber,movingMeshIndependentField)
    movingMeshIndependentField.VariableLabelSet(oc.FieldVariableTypes.U,'MovingMeshIndependent')
    # Set the mesh component to be used by the field components.
    for componentIdx in range(1,numberOfDimensions+1):
        movingMeshIndependentField.ComponentMeshComponentSet(oc.FieldVariableTypes.U,componentIdx,1)    
    # Finish the equations set independent field variables
    movingMeshEquationsSet.IndependentCreateFinish()

    # Initialise independent field moving mesh
    movingMeshIndependentField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,movingMeshKParameter)

    if (progressDiagnostics):
        print('Independent Fields ... Done')

#================================================================================================================================
#  Equations
#================================================================================================================================

if (progressDiagnostics):
    print('Equations ...')

if (problemType != FLUID):
    # Solid equations
    solidEquations = oc.Equations()
    solidEquationsSet.EquationsCreateStart(solidEquations)
    solidEquations.sparsityType = oc.EquationsSparsityTypes.SPARSE
    solidEquations.outputType = solidEquationsOutputType
    solidEquationsSet.EquationsCreateFinish()

if (problemType != SOLID):
    # Fluid equations 
    fluidEquations = oc.Equations()
    fluidEquationsSet.EquationsCreateStart(fluidEquations)
    fluidEquations.sparsityType = oc.EquationsSparsityTypes.SPARSE
    fluidEquations.outputType = fluidEquationsOutputType
    fluidEquationsSet.EquationsCreateFinish()

if (problemType == FSI):
    # Moving mesh equations
    movingMeshEquations = oc.Equations()
    movingMeshEquationsSet.EquationsCreateStart(movingMeshEquations)
    movingMeshEquations.sparsityType = oc.EquationsSparsityTypes.SPARSE
    movingMeshEquations.outputType = movingMeshEquationsOutputType
    movingMeshEquationsSet.EquationsCreateFinish()

if (progressDiagnostics):
    print('Equations ... Done')

#================================================================================================================================
#  CellML
#================================================================================================================================

if (progressDiagnostics):
    print('CellML ...')

if (problemType != SOLID):
    # Create CellML equations for the temporal fluid boundary conditions
    bcCellML = oc.CellML()
    bcCellML.CreateStart(bcCellMLUserNumber,fluidRegion)
    bcCellMLIdx = bcCellML.ModelImport("exponentialrampupinletbc.cellml")
    bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/A")
    bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/B")
    bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/C")
    bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/x")
    bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/y")
    bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/inletx")
    bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/inlety")
    bcCellML.CreateFinish()

    # Create CellML <--> OpenCMISS field maps
    bcCellML.FieldMapsCreateStart()
    # Map geometric field to x0 and y0
    bcCellML.CreateFieldToCellMLMap(fluidGeometricField,oc.FieldVariableTypes.U,1,oc.FieldParameterSetTypes.VALUES,
	                            bcCellMLIdx,"main/x",oc.FieldParameterSetTypes.VALUES)
    bcCellML.CreateFieldToCellMLMap(fluidGeometricField,oc.FieldVariableTypes.U,2,oc.FieldParameterSetTypes.VALUES,
	                            bcCellMLIdx,"main/y",oc.FieldParameterSetTypes.VALUES)
    # Map fluid velocity to ensure dependent field isn't cleared when the velocities are copied back
    bcCellML.CreateFieldToCellMLMap(fluidDependentField,oc.FieldVariableTypes.U,1,oc.FieldParameterSetTypes.VALUES,
	                            bcCellMLIdx,"main/inletx",oc.FieldParameterSetTypes.VALUES)
    bcCellML.CreateFieldToCellMLMap(fluidDependentField,oc.FieldVariableTypes.U,2,oc.FieldParameterSetTypes.VALUES,
	                            bcCellMLIdx,"main/inlety",oc.FieldParameterSetTypes.VALUES)
    # Map inletx and inlety to dependent field
    bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/inletx",oc.FieldParameterSetTypes.VALUES,
	                            fluidDependentField,oc.FieldVariableTypes.U,1,oc.FieldParameterSetTypes.VALUES)
    bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/inlety",oc.FieldParameterSetTypes.VALUES,
	                            fluidDependentField,oc.FieldVariableTypes.U,2,oc.FieldParameterSetTypes.VALUES)
    bcCellML.FieldMapsCreateFinish()


    # Create the CellML models field
    bcCellMLModelsField = oc.Field()
    bcCellML.ModelsFieldCreateStart(bcCellMLModelsFieldUserNumber,bcCellMLModelsField)
    bcCellMLModelsField.VariableLabelSet(oc.FieldVariableTypes.U,"BCModelMap")
    bcCellML.ModelsFieldCreateFinish()

    # Only evaluate BC on inlet nodes
    bcCellMLModelsField.ComponentValuesInitialiseIntg(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,1,0)
    if (debugLevel > 2):
        print('  CellML Boundary Conditions:')
        print('    Inlet Model Set:')
    for yNodeIdx in range(2,numberOfSolidYElements*(numberOfNodesXi-1)+1):
        nodeNumber = (yNodeIdx-1)*((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)+1
        nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            bcCellMLModelsField.ParameterSetUpdateNodeIntg(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                           1,oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,1)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
    for yNodeIdx in range(1,numberOfFluidYElements*(numberOfNodesXi-1)+1):
        nodeNumber = numberOfSolidYElements*(numberOfNodesXi-1)*((numberOfFluidX1Elements+numberOfFluidX2Elements)* \
                                                                 (numberOfNodesXi-1)+2) + (yNodeIdx-1)* \
                                                                 ((numberOfFluidX1Elements+ numberOfSolidXElements+ \
                                                                   numberOfFluidX2Elements)*(numberOfNodesXi-1)+1)+1
        nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            bcCellMLModelsField.ParameterSetUpdateNodeIntg(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,
                                                           1,oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,1)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))

    # Create the CellML state field
    bcCellMLStateField = oc.Field()
    bcCellML.StateFieldCreateStart(bcCellMLStateFieldUserNumber,bcCellMLStateField)
    bcCellMLStateField.VariableLabelSet(oc.FieldVariableTypes.U,"BCState")
    bcCellML.StateFieldCreateFinish()

    # Create the CellML parameters field
    bcCellMLParametersField = oc.Field()
    bcCellML.ParametersFieldCreateStart(bcCellMLParametersFieldUserNumber,bcCellMLParametersField)
    bcCellMLParametersField.VariableLabelSet(oc.FieldVariableTypes.U,"BCParameters")
    bcCellML.ParametersFieldCreateFinish()

    # Get the component numbers
    AComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,oc.CellMLFieldTypes.PARAMETERS,"main/A")
    BComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,oc.CellMLFieldTypes.PARAMETERS,"main/B")
    CComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,oc.CellMLFieldTypes.PARAMETERS,"main/C")
    # Set up the parameters field
    bcCellMLParametersField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,AComponentNumber,A)
    bcCellMLParametersField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,BComponentNumber,B)
    bcCellMLParametersField.ComponentValuesInitialiseDP(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,CComponentNumber,C)

    # Create the CELL intermediate field
    bcCellMLIntermediateField = oc.Field()
    bcCellML.IntermediateFieldCreateStart(bcCellMLIntermediateFieldUserNumber,bcCellMLIntermediateField)
    bcCellMLIntermediateField.VariableLabelSet(oc.FieldVariableTypes.U,"BCIntermediate")
    bcCellML.IntermediateFieldCreateFinish()

if (progressDiagnostics):
    print('CellML ... Done')

#================================================================================================================================
#  Interface Condition
#================================================================================================================================

if (problemType == FSI):
    if (progressDiagnostics):
        print('Interface Conditions ...')

    # Create an interface condition between the two meshes
    interfaceCondition = oc.InterfaceCondition()
    interfaceCondition.CreateStart(interfaceConditionUserNumber,interface,interfaceGeometricField)
    # Specify the method for the interface condition
    interfaceCondition.MethodSet(oc.InterfaceConditionMethods.LAGRANGE_MULTIPLIERS)
    # Specify the type of interface condition operator
    interfaceCondition.OperatorSet(oc.InterfaceConditionOperators.SOLID_FLUID)
    # Add in the dependent variables from the equations sets
    interfaceCondition.DependentVariableAdd(solidMeshIndex,solidEquationsSet,oc.FieldVariableTypes.U)
    interfaceCondition.DependentVariableAdd(fluidMeshIndex,fluidEquationsSet,oc.FieldVariableTypes.U)
    # Set the label
    interfaceCondition.LabelSet("FSI Interface Condition")
    # Set the output type
    interfaceCondition.OutputTypeSet(interfaceConditionOutputType)
    # Finish creating the interface condition
    interfaceCondition.CreateFinish()

    if (progressDiagnostics):
        print('Interface Conditions ... Done')

    if (progressDiagnostics):
        print('Interface Lagrange Field ...')
    
    # Create the Lagrange multipliers field
    interfaceLagrangeField = oc.Field()
    interfaceCondition.LagrangeFieldCreateStart(interfaceLagrangeFieldUserNumber,interfaceLagrangeField)
    interfaceLagrangeField.VariableLabelSet(oc.FieldVariableTypes.U,'InterfaceLagrange')
    # Finish the Lagrange multipliers field
    interfaceCondition.LagrangeFieldCreateFinish()
    
    for componentIdx in range(1,numberOfDimensions+1):
        interfaceLagrangeField.ComponentValuesInitialise(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES,componentIdx,0.0)

        interfaceLagrangeField.ParameterSetUpdateStart(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)
        interfaceLagrangeField.ParameterSetUpdateFinish(oc.FieldVariableTypes.U,oc.FieldParameterSetTypes.VALUES)

    if (progressDiagnostics):
        print('Interface Lagrange Field ... Done')

    if (progressDiagnostics):
        print('Interface Equations ...')

    # Create the interface condition equations
    interfaceEquations = oc.InterfaceEquations()
    interfaceCondition.EquationsCreateStart(interfaceEquations)
    # Set the interface equations sparsity
    interfaceEquations.sparsityType = oc.EquationsSparsityTypes.SPARSE
    # Set the interface equations output
    interfaceEquations.outputType = interfaceEquationsOutputType
    # Finish creating the interface equations
    interfaceCondition.EquationsCreateFinish()

    if (progressDiagnostics):
        print('Interface Equations ... Done')

#================================================================================================================================
#  Problem
#================================================================================================================================

if (progressDiagnostics):
    print('Problems ...')

# Create a FSI problem
fsiProblem = oc.Problem()
if (problemType == SOLID):
   fsiProblemSpecification = [oc.ProblemClasses.ELASTICITY,
                              oc.ProblemTypes.FINITE_ELASTICITY,
                              oc.ProblemSubtypes.QUASISTATIC_FINITE_ELASTICITY]
elif (problemType == FLUID):
    if RBS:
        fsiProblemSpecification = [oc.ProblemClasses.FLUID_MECHANICS,
                                   oc.ProblemTypes.NAVIER_STOKES_EQUATION,
                                   oc.ProblemSubtypes.TRANSIENT_RBS_NAVIER_STOKES]
    else:
        fsiProblemSpecification = [oc.ProblemClasses.FLUID_MECHANICS,
                                   oc.ProblemTypes.NAVIER_STOKES_EQUATION,
                                   oc.ProblemSubtypes.TRANSIENT_NAVIER_STOKES]
elif (problemType == FSI):
    if RBS:
        fsiProblemSpecification = [oc.ProblemClasses.MULTI_PHYSICS,
                                   oc.ProblemTypes.FINITE_ELASTICITY_NAVIER_STOKES,
                                   oc.ProblemSubtypes.FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE]
    else:
        fsiProblemSpecification = [oc.ProblemClasses.MULTI_PHYSICS,
                                   oc.ProblemTypes.FINITE_ELASTICITY_NAVIER_STOKES,
                                   oc.ProblemSubtypes.FINITE_ELASTICITY_NAVIER_STOKES_ALE]
        
fsiProblem.CreateStart(fsiProblemUserNumber,context,fsiProblemSpecification)
fsiProblem.CreateFinish()

if (progressDiagnostics):
    print('Problems ... Done')

#================================================================================================================================
#  Control Loop
#================================================================================================================================

if (progressDiagnostics):
    print('Control Loops ...')

# Create the fsi problem control loop
fsiControlLoop = oc.ControlLoop()
fsiProblem.ControlLoopCreateStart()
fsiProblem.ControlLoopGet([oc.ControlLoopIdentifiers.NODE],fsiControlLoop)
fsiControlLoop.LabelSet('TimeLoop')
fsiControlLoop.TimesSet(startTime,stopTime,timeStep)
fsiControlLoop.TimeOutputSet(outputFrequency)
fsiProblem.ControlLoopCreateFinish()

if (progressDiagnostics):
    print('Control Loops ... Done')

#================================================================================================================================
#  Solvers
#================================================================================================================================

if (progressDiagnostics):
    print('Solvers ...')

# Create the problem solver
bcCellMLEvaluationSolver = oc.Solver()
fsiDynamicSolver = oc.Solver()
fsiNonlinearSolver = oc.Solver()
fsiLinearSolver = oc.Solver()
movingMeshLinearSolver = oc.Solver()

fsiProblem.SolversCreateStart()
if (problemType == SOLID):
    # Solvers for growth Finite Elasticity problem
    # Get the BC CellML solver
    fsiProblem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,bcCellMLEvaluationSolver)
    bcCellMLEvaluationSolver.outputType = oc.SolverOutputTypes.PROGRESS
    # Get the nonlinear solver
    fsiProblem.SolverGet([oc.ControlLoopIdentifiers.NODE],2,fsiNonlinearSolver)
    fsiNonlinearSolver.NewtonLineSearchTypeSet(oc.NewtonLineSearchTypes.LINEAR)
    fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(oc.JacobianCalculationTypes.EQUATIONS) #(.FD/EQUATIONS)
    #fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(oc.JacobianCalculationTypes.FD) #(.FD/EQUATIONS)
    fsiNonlinearSolver.NewtonMaximumFunctionEvaluationsSet(nonlinearMaxFunctionEvaluations)
    fsiNonlinearSolver.OutputTypeSet(fsiNonlinearSolverOutputType)
    fsiNonlinearSolver.NewtonAbsoluteToleranceSet(nonlinearAbsoluteTolerance)
    fsiNonlinearSolver.NewtonMaximumIterationsSet(nonlinearMaximumIterations)
    fsiNonlinearSolver.NewtonRelativeToleranceSet(nonlinearRelativeTolerance)
    fsiNonlinearSolver.NewtonLineSearchAlphaSet(nonlinearLinesearchAlpha)
    # Get the dynamic nonlinear linear solver
    fsiNonlinearSolver.NewtonLinearSolverGet(fsiLinearSolver)
    #fsiLinearSolver.LinearTypeSet(oc.LinearSolverTypes.ITERATIVE)
    #fsiLinearSolver.LinearIterativeTypeSet(oc.IterativeLinearSolverTypes.GMRES)
    #fsiLinearSolver.LinearIterativeGMRESRestartSet(linearRestartValue)
    #fsiLinearSolver.LinearIterativeMaximumIterationsSet(linearMaximumIterations)
    #fsiLinearSolver.LinearIterativeDivergenceToleranceSet(linearDivergenceTolerance)
    #fsiLinearSolver.LinearIterativeRelativeToleranceSet(linearRelativeTolerance)
    #fsiLinearSolver.LinearIterativeAbsoluteToleranceSet(linearAbsoluteTolerance)
    fsiLinearSolver.LinearTypeSet(oc.LinearSolverTypes.DIRECT)
    fsiLinearSolver.OutputTypeSet(fsiLinearSolverOutputType)
elif (problemType == FLUID):
    # Solvers for coupled FiniteElasticity NavierStokes problem
    # Get the BC CellML solver
    fsiProblem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,bcCellMLEvaluationSolver)
    bcCellMLEvaluationSolver.outputType = oc.SolverOutputTypes.PROGRESS
    # Get the dynamic ALE solver
    fsiProblem.SolverGet([oc.ControlLoopIdentifiers.NODE],2,fsiDynamicSolver)
    fsiDynamicSolver.OutputTypeSet(fsiDynamicSolverOutputType)
    fsiDynamicSolver.DynamicThetaSet(fsiDynamicSolverTheta)
    # Get the dynamic nonlinear solver
    fsiDynamicSolver.DynamicNonlinearSolverGet(fsiNonlinearSolver)
    fsiNonlinearSolver.NewtonLineSearchTypeSet(oc.NewtonLineSearchTypes.LINEAR)
    fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(oc.JacobianCalculationTypes.EQUATIONS) #(.FD/EQUATIONS)
    #fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(oc.JacobianCalculationTypes.FD) #(.FD/EQUATIONS)
    fsiNonlinearSolver.NewtonMaximumFunctionEvaluationsSet(nonlinearMaxFunctionEvaluations)
    fsiNonlinearSolver.OutputTypeSet(fsiNonlinearSolverOutputType)
    fsiNonlinearSolver.NewtonAbsoluteToleranceSet(nonlinearAbsoluteTolerance)
    fsiNonlinearSolver.NewtonMaximumIterationsSet(nonlinearMaximumIterations)
    fsiNonlinearSolver.NewtonRelativeToleranceSet(nonlinearRelativeTolerance)
    fsiNonlinearSolver.NewtonLineSearchAlphaSet(nonlinearLinesearchAlpha)
    # Get the dynamic nonlinear linear solver
    fsiNonlinearSolver.NewtonLinearSolverGet(fsiLinearSolver)
    #fsiLinearSolver.LinearTypeSet(oc.LinearSolverTypes.ITERATIVE)
    #fsiLinearSolver.LinearIterativeMaximumIterationsSet(linearMaximumIterations)
    #fsiLinearSolver.LinearIterativeDivergenceToleranceSet(linearDivergenceTolerance)
    #fsiLinearSolver.LinearIterativeRelativeToleranceSet(linearRelativeTolerance)
    #fsiLinearSolver.LinearIterativeAbsoluteToleranceSet(linearAbsoluteTolerance)
    fsiLinearSolver.LinearTypeSet(oc.LinearSolverTypes.DIRECT)
    fsiLinearSolver.OutputTypeSet(fsiLinearSolverOutputType)
elif (problemType == FSI):
    # Solvers for coupled FiniteElasticity NavierStokes problem
    # Get the BC CellML solver
    fsiProblem.SolverGet([oc.ControlLoopIdentifiers.NODE],1,bcCellMLEvaluationSolver)
    bcCellMLEvaluationSolver.outputType = oc.SolverOutputTypes.PROGRESS
    # Get the dynamic ALE solver
    fsiProblem.SolverGet([oc.ControlLoopIdentifiers.NODE],2,fsiDynamicSolver)
    fsiDynamicSolver.OutputTypeSet(fsiDynamicSolverOutputType)
    fsiDynamicSolver.DynamicThetaSet(fsiDynamicSolverTheta)
    # Get the dynamic nonlinear solver
    fsiDynamicSolver.DynamicNonlinearSolverGet(fsiNonlinearSolver)
    fsiNonlinearSolver.NewtonLineSearchTypeSet(oc.NewtonLineSearchTypes.LINEAR)
    fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(oc.JacobianCalculationTypes.EQUATIONS) #(.FD/EQUATIONS)
    #fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(oc.JacobianCalculationTypes.FD) #(.FD/EQUATIONS)
    fsiNonlinearSolver.NewtonMaximumFunctionEvaluationsSet(nonlinearMaxFunctionEvaluations)
    fsiNonlinearSolver.OutputTypeSet(fsiNonlinearSolverOutputType)
    fsiNonlinearSolver.NewtonAbsoluteToleranceSet(nonlinearAbsoluteTolerance)
    fsiNonlinearSolver.NewtonMaximumIterationsSet(nonlinearMaximumIterations)
    fsiNonlinearSolver.NewtonRelativeToleranceSet(nonlinearRelativeTolerance)
    fsiNonlinearSolver.NewtonLineSearchAlphaSet(nonlinearLinesearchAlpha)
    # Get the dynamic nonlinear linear solver
    fsiNonlinearSolver.NewtonLinearSolverGet(fsiLinearSolver)
    #fsiLinearSolver.LinearTypeSet(oc.LinearSolverTypes.ITERATIVE)
    #fsiLinearSolver.LinearIterativeMaximumIterationsSet(linearMaximumIterations)
    #fsiLinearSolver.LinearIterativeDivergenceToleranceSet(linearDivergenceTolerance)
    #fsiLinearSolver.LinearIterativeRelativeToleranceSet(linearRelativeTolerance)
    #fsiLinearSolver.LinearIterativeAbsoluteToleranceSet(linearAbsoluteTolerance)
    fsiLinearSolver.LinearTypeSet(oc.LinearSolverTypes.DIRECT)
    fsiLinearSolver.OutputTypeSet(fsiLinearSolverOutputType)
    # Linear solver for moving mesh
    fsiProblem.SolverGet([oc.ControlLoopIdentifiers.NODE],3,movingMeshLinearSolver)
    movingMeshLinearSolver.OutputTypeSet(movingMeshLinearSolverOutputType)
# Finish the creation of the problem solver
fsiProblem.SolversCreateFinish()

if (progressDiagnostics):
    print('Solvers ... Done')

#================================================================================================================================
#  CellML Equations
#================================================================================================================================

if (progressDiagnostics):
    print('CellML Equations ...')

if (problemType != SOLID):
    # Create CellML equations and add BC equations to the solver
    bcEquations = oc.CellMLEquations()
    fsiProblem.CellMLEquationsCreateStart()
    bcCellMLEvaluationSolver.CellMLEquationsGet(bcEquations)
    bcEquationsIndex = bcEquations.CellMLAdd(bcCellML)
    fsiProblem.CellMLEquationsCreateFinish()

if (progressDiagnostics):
    print('CellML Equations ... Done')

#================================================================================================================================
#  Solver Equations
#================================================================================================================================

if (progressDiagnostics):
    print('Solver Equations ...')

# Start the creation of the fsi problem solver equations
fsiProblem.SolverEquationsCreateStart()
# Get the fsi dynamic solver equations
fsiSolverEquations = oc.SolverEquations()
if (problemType == SOLID):
    fsiNonlinearSolver.SolverEquationsGet(fsiSolverEquations)
else:
    fsiDynamicSolver.SolverEquationsGet(fsiSolverEquations)
fsiSolverEquations.sparsityType = oc.SolverEquationsSparsityTypes.SPARSE
if (problemType != FLUID):
    fsiSolidEquationsSetIndex = fsiSolverEquations.EquationsSetAdd(solidEquationsSet)
if (problemType != SOLID):
    fsiFluidEquationsSetIndex = fsiSolverEquations.EquationsSetAdd(fluidEquationsSet)
if (problemType == FSI):
    fsiInterfaceConditionIndex = fsiSolverEquations.InterfaceConditionAdd(interfaceCondition)
    # Set the time dependence of the interface matrix to determine the interface matrix coefficient in the solver matrix
    # (basiy position in big coupled matrix system)
    interfaceEquations.MatrixTimeDependenceTypeSet(fsiSolidEquationsSetIndex,True, \
                                                   [oc.InterfaceMatricesTimeDependenceTypes.STATIC,\
                                                    oc.InterfaceMatricesTimeDependenceTypes.FIRST_ORDER_DYNAMIC])
    interfaceEquations.MatrixTimeDependenceTypeSet(fsiFluidEquationsSetIndex,True, \
                                                   [oc.InterfaceMatricesTimeDependenceTypes.STATIC,\
                                                    oc.InterfaceMatricesTimeDependenceTypes.STATIC])
    
    # Create the moving mesh solver equations
    movingMeshSolverEquations = oc.SolverEquations()
    # Get the linear moving mesh solver equations
    movingMeshLinearSolver.SolverEquationsGet(movingMeshSolverEquations)
    movingMeshSolverEquations.sparsityType = oc.SolverEquationsSparsityTypes.SPARSE
    # Add in the equations set
    movingMeshEquationsSetIndex = movingMeshSolverEquations.EquationsSetAdd(movingMeshEquationsSet)
    
# Finish the creation of the fsi problem solver equations
fsiProblem.SolverEquationsCreateFinish()

if (progressDiagnostics):
    print('Solver Equations ...')

#================================================================================================================================
#  Boundary Conditions
#================================================================================================================================

if (progressDiagnostics):
    print('Boundary Conditions ...')

# Start the creation of the fsi boundary conditions
fsiBoundaryConditions = oc.BoundaryConditions()
fsiSolverEquations.BoundaryConditionsCreateStart(fsiBoundaryConditions)
if (problemType != FLUID):
    # Set no displacement boundary conditions on the bottom edge of the solid
    if (debugLevel > 2):
        print('  Solid Boundary Conditions:')
        print('    No Displacement Boundary conditions:')
    for xNodeIdx in range(1,numberOfSolidXElements*(numberOfNodesXi-1)+2):
        nodeNumber = xNodeIdx
        nodeDomain = solidDecomposition.NodeDomainGet(1,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.AddNode(solidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,1,oc.BoundaryConditionsTypes.FIXED,0.0)
            fsiBoundaryConditions.AddNode(solidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Displacement      = [ %.2f, %.2f ]' % (0.0,0.0))                 
            if (useHermite):
                fsiBoundaryConditions.AddNode(solidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(solidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(solidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(solidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(solidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(solidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
    if (debugLevel > 2):
        print('    Reference Solid Pressure Boundary Condition:')
        nodeNumber = numberOfSolidXElements*(numberOfNodesXi-1)+1
        nodeDomain = solidDecomposition.NodeDomainGet(1,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(solidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,3,oc.BoundaryConditionsTypes.FIXED,solidPRef)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Pressure         =   %.2f' % (solidPRef))

if (problemType != SOLID):                
    # Set inlet boundary conditions on the left hand edge
    if (debugLevel > 2):
        print('  Fluid Boundary Conditions:')
        print('    Inlet Boundary conditions:')
    for yNodeIdx in range(2,numberOfSolidYElements*(numberOfNodesXi-1)+1):
        nodeNumber = (yNodeIdx-1)*((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)+1
        nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,1,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
            fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,2,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Velocity         = [ %.2f, %.2f ]' % (0.0,0.0))                 
            if (useHermite):
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
    for yNodeIdx in range(1,numberOfFluidYElements*(numberOfNodesXi-1)+1):
        nodeNumber = numberOfSolidYElements*(numberOfNodesXi-1)*((numberOfFluidX1Elements+numberOfFluidX2Elements)* \
                                                                 (numberOfNodesXi-1)+2) + (yNodeIdx-1)* \
                                                                 ((numberOfFluidX1Elements+ numberOfSolidXElements+ \
                                                                   numberOfFluidX2Elements)*(numberOfNodesXi-1)+1)+1
        nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,1,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
            fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,2,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Velocity         = [ %.2f, %.2f ]' % (0.0,0.0))                 
            if (useHermite):
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED_INLET,0.0)
    # Set outlet boundary conditions on the right hand edge to have zero pressure
    if (debugLevel > 2):
        print('    Outlet Boundary conditions:')
    # Elements to the right of the solid
    for yElementIdx in range(2,numberOfSolidYElements+1):
        nodeNumber = (yElementIdx-1)*(numberOfNodesXi-1)*((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)+\
                     (numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2
        nodeDomain = fluidDecomposition.NodeDomainGet(2,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            if RBS:
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                              nodeNumber,numberOfDimensions+1,oc.BoundaryConditionsTypes.PRESSURE,fluidPRef)
            else:
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.DELUDELN,1, \
                                              oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                              nodeNumber,numberOfDimensions+1,oc.BoundaryConditionsTypes.FIXED,fluidPRef)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Pressure         =   %.2f' % (fluidPRef))                 
        if RBS:
            # Set the element normals for outlet stabilisation
            elementNumber = numberOfFluidX1Elements+numberOfFluidX2Elements+(yElementIdx-2)*(numberOfFluidX1Elements+numberOfFluidX2Elements)
            elementDomain = fluidDecomposition.ElementDomainGet(elementNumber)
            if (elementDomain == computationalNodeNumber):
                # Set the outflow normal to (0,0,+1)
                fluidEquationsSetField.ParameterSetUpdateElementDP(oc.FieldVariableTypes.V, \
                                                                   oc.FieldParameterSetTypes.VALUES, \
                                                                   elementNumber,5,+1.0)
                fluidEquationsSetField.ParameterSetUpdateElementDP(oc.FieldVariableTypes.V, \
                                                                   oc.FieldParameterSetTypes.VALUES, \
                                                                   elementNumber,6,0.0)
                # Set the boundary type
                fluidEquationsSetField.ParameterSetUpdateElementDP(oc.FieldVariableTypes.V, \
                                                                   oc.FieldParameterSetTypes.VALUES, \
                                                                   elementNumber,9, \
                                                                   oc.BoundaryConditionsTypes.PRESSURE)                                                
                if (debugLevel > 2):
                    print('      Element     %d:' % (elementNumber))
                    print('         Normal          = [ %.2f, %.2f ]' % (+1.0,0.0))
    # Elements above the solid
    for yElementIdx in range(1,numberOfFluidYElements+2):
        nodeNumber = numberOfSolidYElements*(numberOfNodesXi-1)*((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2) \
                     + ((yElementIdx-1)*(numberOfNodesXi-1)+1)*((numberOfFluidX1Elements+numberOfSolidXElements+numberOfFluidX2Elements)*\
                                                            (numberOfNodesXi-1)+1)
        nodeDomain = fluidDecomposition.NodeDomainGet(2,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            if RBS:
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                              nodeNumber,numberOfDimensions+1,oc.BoundaryConditionsTypes.PRESSURE,fluidPRef)
            else:
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.DELUDELN,1, \
                                              oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                              nodeNumber,numberOfDimensions+1,oc.BoundaryConditionsTypes.FIXED,fluidPRef)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Pressure         =   %.2f' % (fluidPRef))
        if RBS:
            # Set the element normals for outlet stabilisation
            elementNumber = (numberOfFluidX1Elements+numberOfFluidX2Elements)*numberOfSolidYElements+\
                            numberOfFluidX1Elements+numberOfFluidX2Elements+numberOfSolidXElements+\
                            (yElementIdx-2)*(numberOfFluidX1Elements+numberOfFluidX2Elements+numberOfSolidXElements)
            elementDomain = fluidDecomposition.ElementDomainGet(elementNumber)
            if (elementDomain == computationalNodeNumber):
                # Set the outflow normal to (0,0,+1)
                fluidEquationsSetField.ParameterSetUpdateElementDP(oc.FieldVariableTypes.V, \
                                                                   oc.FieldParameterSetTypes.VALUES, \
                                                                   elementNumber,5,+1.0)
                fluidEquationsSetField.ParameterSetUpdateElementDP(oc.FieldVariableTypes.V, \
                                                                   oc.FieldParameterSetTypes.VALUES, \
                                                                   elementNumber,6,0.0)
                # Set the boundary type
                fluidEquationsSetField.ParameterSetUpdateElementDP(oc.FieldVariableTypes.V, \
                                                                   oc.FieldParameterSetTypes.VALUES, \
                                                                   elementNumber,9, \
                                                                   oc.BoundaryConditionsTypes.PRESSURE)
                if (debugLevel > 2):
                        print('      Element     %d:' % (elementNumber))
                        print('         Normal          = [ %.2f, %.2f ]' % (+1.0,0.0))
                
    # Set no-slip boundary conditions on the bottom edge
    if (debugLevel > 2):
        print('    No-slip Boundary conditions:')
    for xNodeIdx in range(1,numberOfFluidX1Elements*(numberOfNodesXi-1)+2):
        nodeNumber = xNodeIdx
        nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,1,oc.BoundaryConditionsTypes.FIXED,0.0)
            fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
            if (useHermite):
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
    for xNodeIdx in range(1,numberOfFluidX2Elements*(numberOfNodesXi-1)+2):
        nodeNumber = numberOfFluidX1Elements*(numberOfNodesXi-1)+1+xNodeIdx
        nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,1,oc.BoundaryConditionsTypes.FIXED,0.0)
            fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
            if (useHermite):
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,1,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
    if (problemType == FLUID):
        # Set no slip around the solid
        # Left and right solid edge nodes
        for yNodeIdx in range(2,numberOfSolidYElements*(numberOfNodesXi-1)+1):
            nodeNumber1 = (yNodeIdx-1)*((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)+ \
                     numberOfFluidX1Elements*(numberOfNodesXi-1)+1
            nodeNumber2 = nodeNumber1+1
            nodeDomain1 = fluidDecomposition.NodeDomainGet(1,nodeNumber1)
            nodeDomain2 = fluidDecomposition.NodeDomainGet(1,nodeNumber2)
            if (nodeDomain1 == computationalNodeNumber):
                fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                              oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,1,
                                              oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                              oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,2,
                                              oc.BoundaryConditionsTypes.FIXED,0.0)
                if (debugLevel > 2):
                    print('      Node        %d:' % (nodeNumber1))
                if (useHermite):    
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,1,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,1,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,1,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,2,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,2,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,2,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
            if (nodeDomain2 == computationalNodeNumber):
                fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                              oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,1,
                                              oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                              oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,2,
                                              oc.BoundaryConditionsTypes.FIXED,0.0)
                if (debugLevel > 2):
                    print('      Node        %d:' % (nodeNumber2))
                if (useHermite):    
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,1,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,1,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,1,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,2,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,2,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,2,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
        # Top solid edge nodes
        for xNodeIdx in range(1,numberOfSolidXElements*(numberOfNodesXi-1)+2):
            nodeNumber = xNodeIdx+numberOfSolidYElements*(numberOfNodesXi-1)*((numberOfFluidX1Elements+numberOfFluidX2Elements)* \
                      (numberOfNodesXi-1)+2)+numberOfFluidX1Elements*(numberOfNodesXi-1)
            nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
            if (nodeDomain == computationalNodeNumber):
                fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                              oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                              oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                              oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                              oc.BoundaryConditionsTypes.FIXED,0.0)
                if (debugLevel > 2):
                    print('      Node        %d:' % (nodeNumber))
                if (useHermite):    
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,1,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,oc.FieldVariableTypes.U,1,
                                                  oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,2,
                                                  oc.BoundaryConditionsTypes.FIXED,0.0)
    # Set slip boundary conditions on the top edge
    if (debugLevel > 2):
        print('    Slip Boundary conditions:')
    for xNodeIdx in range(1,(numberOfFluidX1Elements+numberOfSolidXElements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2):
        nodeNumber = ((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)*numberOfSolidYElements*(numberOfNodesXi-1)+ \
                     ((numberOfFluidX1Elements+numberOfSolidXElements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+1)* \
                     numberOfFluidYElements*(numberOfNodesXi-1)+xNodeIdx
        nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,1,oc.BoundaryConditionsTypes.FIXED,0.0)
            fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
            if (useHermite):
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                              oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,2,oc.BoundaryConditionsTypes.FIXED,0.0)
    if (debugLevel > 2):
        print('    Reference Fluid Pressure Boundary Condition:')
        nodeNumber = (numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2
        nodeDomain = fluidDecomposition.NodeDomainGet(2,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(fluidDependentField,oc.FieldVariableTypes.U,1, \
                                          oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,3,oc.BoundaryConditionsTypes.FIXED,fluidPRef)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Pressure         =   %.2f' % (fluidPRef))

if (problemType == FSI):
    # Remove dof's at nodes where solid displacement and zero velocity is set (first n last interface node)
    if (debugLevel > 2):
        print('  Lagrange Boundary Conditions:')
        print('    Fixed Boundary conditions:')
    fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                  oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,1,1, \
                                  oc.BoundaryConditionsTypes.FIXED,0.0)
    fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                  oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,1,2, \
                                  oc.BoundaryConditionsTypes.FIXED,0.0)
    fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                  oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,numberOfInterfaceNodes,1, \
                                  oc.BoundaryConditionsTypes.FIXED,0.0)
    fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                  oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,numberOfInterfaceNodes,2, \
                                  oc.BoundaryConditionsTypes.FIXED,0.0)
    if (debugLevel > 2):
        print('      Node        %d:' % (1))
        print('      Node        %d:' % (numberOfInterfaceNodes))
    if (useHermite):
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                      oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,1,1, \
                                      oc.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                      oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,1,1, \
                                      oc.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                      oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,1,1, \
                                      oc.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                      oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,1,2, \
                                      oc.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                      oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,1,2, \
                                      oc.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                      oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,1,2, \
                                      oc.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                      oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,numberOfInterfaceNodes,1, \
                                      oc.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                      oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,numberOfInterfaceNodes,1, \
                                      oc.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                      oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,numberOfInterfaceNodes,1, \
                                      oc.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                      oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,numberOfInterfaceNodes,2, \
                                      oc.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                      oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,numberOfInterfaceNodes,2, \
                                      oc.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,oc.FieldVariableTypes.U,1, \
                                      oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,numberOfInterfaceNodes,2, \
                                      oc.BoundaryConditionsTypes.FIXED,0.0)               
# Finish FSI boundary conditions
fsiSolverEquations.BoundaryConditionsCreateFinish()

if (problemType == FSI):
    # Start the creation of the moving mesh boundary conditions
    movingMeshBoundaryConditions = oc.BoundaryConditions()
    movingMeshSolverEquations.BoundaryConditionsCreateStart(movingMeshBoundaryConditions)
    if (debugLevel > 2):
        print('  Moving Mesh Boundary Conditions:')
        print('    Fixed Wall Boundary conditions:')
    # Bottom edge nodes
    for xNodeIdx in range(1,(numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+3):
        nodeNumber = xNodeIdx
        nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                                oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                                oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
            if (useHermite):    
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
    # Side edges nodes
    for yNodeIdx in range(2,numberOfSolidYElements*(numberOfNodesXi-1)+1):
        nodeNumber1 = (yNodeIdx-1)*((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)+1
        nodeNumber2 = yNodeIdx*((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)
        nodeDomain1 = fluidDecomposition.NodeDomainGet(1,nodeNumber1)
        nodeDomain2 = fluidDecomposition.NodeDomainGet(1,nodeNumber2)
        if (nodeDomain1 == computationalNodeNumber):
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,1,
                                                oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,2,
                                                oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber1))
            if (useHermite):    
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
        if (nodeDomain2 == computationalNodeNumber):
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,1,
                                                oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,2,
                                                oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber2))
            if (useHermite):    
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
    for yNodeIdx in range(1,numberOfFluidYElements*(numberOfNodesXi-1)+1):
        nodeNumber1 = (yNodeIdx-1)*((numberOfFluidX1Elements+numberOfSolidXElements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+1)+1+\
                      ((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)*\
                      numberOfSolidYElements*(numberOfNodesXi-1)
        nodeNumber2 = yNodeIdx*((numberOfFluidX1Elements+numberOfSolidXElements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+1)+\
                      ((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)*\
                      numberOfSolidYElements*(numberOfNodesXi-1)
        nodeDomain1 = fluidDecomposition.NodeDomainGet(1,nodeNumber1)
        nodeDomain2 = fluidDecomposition.NodeDomainGet(1,nodeNumber2)
        if (nodeDomain1 == computationalNodeNumber):
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,1,
                                                oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,2,
                                                oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber1))
            if (useHermite):    
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
        if (nodeDomain2 == computationalNodeNumber):
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,1,
                                                oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,2,
                                                oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber2))
            if (useHermite):    
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
    # Top edge nodes
    for xNodeIdx in range(1,(numberOfFluidX1Elements+numberOfSolidXElements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2):
        nodeNumber = xNodeIdx +((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)* \
                     numberOfSolidYElements*(numberOfNodesXi-1)+numberOfFluidYElements*(numberOfNodesXi-1)* \
                     ((numberOfFluidX1Elements+ numberOfSolidXElements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+1)
        nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                                oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                                oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
            if (useHermite):    
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,1,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,2,
                                                    oc.BoundaryConditionsTypes.FIXED_WALL,0.0)
    if (debugLevel > 2):
        print('    Moving Wall Boundary conditions:')
    # Left and right solid edge nodes
    for yNodeIdx in range(2,numberOfSolidYElements*(numberOfNodesXi-1)+1):
        nodeNumber1 = (yNodeIdx-1)*((numberOfFluidX1Elements+numberOfFluidX2Elements)*(numberOfNodesXi-1)+2)+ \
                     numberOfFluidX1Elements*(numberOfNodesXi-1)+1
        nodeNumber2 = nodeNumber1+1
        nodeDomain1 = fluidDecomposition.NodeDomainGet(1,nodeNumber1)
        nodeDomain2 = fluidDecomposition.NodeDomainGet(1,nodeNumber2)
        if (nodeDomain1 == computationalNodeNumber):
            movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,1,
                                                oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
            movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,2,
                                                oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber1))
            if (useHermite):    
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,1,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,1,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,1,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,2,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,2,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,2,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
        if (nodeDomain2 == computationalNodeNumber):
            movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,1,
                                                oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
            movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,2,
                                                oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber2))
            if (useHermite):    
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,1,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,1,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,1,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,2,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,2,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,2,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
    # Top solid edge nodes
    for xNodeIdx in range(1,numberOfSolidXElements*(numberOfNodesXi-1)+2):
        nodeNumber = xNodeIdx+numberOfSolidYElements*(numberOfNodesXi-1)*((numberOfFluidX1Elements+numberOfFluidX2Elements)* \
                      (numberOfNodesXi-1)+2)+numberOfFluidX1Elements*(numberOfNodesXi-1)
        nodeDomain = fluidDecomposition.NodeDomainGet(1,nodeNumber)
        if (nodeDomain == computationalNodeNumber):
            movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                                oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
            movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                oc.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                                oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
            if (useHermite):    
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,1,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,oc.FieldVariableTypes.U,1,
                                                    oc.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,2,
                                                    oc.BoundaryConditionsTypes.MOVED_WALL,0.0)

    # Finish moving mesh boundary conditions
    movingMeshSolverEquations.BoundaryConditionsCreateFinish()

if (progressDiagnostics):
    print('Boundary Conditions ... Done')

#================================================================================================================================
#  Run Solvers
#================================================================================================================================

#quit()

# Solve the problem
print('Solving problem...')
start = time.time()
fsiProblem.Solve()
end = time.time()
elapsed = end - start
print('Calculation Time = %3.4f' %elapsed)
print('Problem solved!')
print('#')

#================================================================================================================================
#  Finish Program
#================================================================================================================================
