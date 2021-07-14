# script-version: 2.0
# Catalyst state generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [844, 539]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.0, -1.043081283569336e-06, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [2.443379419932841, 17.382832235131602, -20.055095079935207]
renderView1.CameraFocalPoint = [0.0, -1.043081283569336e-06, 0.0]
renderView1.CameraViewUp = [0.9500725701327958, 0.16904480929559032, 0.26227078360252315]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 10.09950424491012
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(844, 539)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'ExodusIIReader'
input = ExodusIIReader(registrationName='input', FileName=['/data/tjotaha/src/nalu_catalyst_5.9/Nalu/reg_tests/mesh/pipeTet.g.8.0', '/data/tjotaha/src/nalu_catalyst_5.9/Nalu/reg_tests/mesh/pipeTet.g.8.1', '/data/tjotaha/src/nalu_catalyst_5.9/Nalu/reg_tests/mesh/pipeTet.g.8.2', '/data/tjotaha/src/nalu_catalyst_5.9/Nalu/reg_tests/mesh/pipeTet.g.8.3', '/data/tjotaha/src/nalu_catalyst_5.9/Nalu/reg_tests/mesh/pipeTet.g.8.4', '/data/tjotaha/src/nalu_catalyst_5.9/Nalu/reg_tests/mesh/pipeTet.g.8.5', '/data/tjotaha/src/nalu_catalyst_5.9/Nalu/reg_tests/mesh/pipeTet.g.8.6', '/data/tjotaha/src/nalu_catalyst_5.9/Nalu/reg_tests/mesh/pipeTet.g.8.7'])
input.SideSetArrayStatus = []
input.ElementBlocks = ['Unnamed block ID: 1']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from input
inputDisplay = Show(input, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')
vtkBlockColorsLUT.InterpretValuesAsCategories = 1
vtkBlockColorsLUT.AnnotationsInitialized = 1
vtkBlockColorsLUT.Annotations = ['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10', '11', '11']
vtkBlockColorsLUT.ActiveAnnotatedValues = ['0', '1']
vtkBlockColorsLUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.63, 0.63, 1.0, 0.67, 0.5, 0.33, 1.0, 0.5, 0.75, 0.53, 0.35, 0.7, 1.0, 0.75, 0.5]

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# trace defaults for the display properties.
inputDisplay.Representation = 'Surface With Edges'
inputDisplay.ColorArrayName = ['FIELD', 'vtkBlockColors']
inputDisplay.LookupTable = vtkBlockColorsLUT
inputDisplay.SelectTCoordArray = 'None'
inputDisplay.SelectNormalArray = 'None'
inputDisplay.SelectTangentArray = 'None'
inputDisplay.OSPRayScaleArray = 'GlobalNodeId'
inputDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
inputDisplay.SelectOrientationVectors = 'None'
inputDisplay.ScaleFactor = 2.0
inputDisplay.SelectScaleArray = 'GlobalNodeId'
inputDisplay.GlyphType = 'Arrow'
inputDisplay.GlyphTableIndexArray = 'GlobalNodeId'
inputDisplay.GaussianRadius = 0.1
inputDisplay.SetScaleArray = ['POINTS', 'GlobalNodeId']
inputDisplay.ScaleTransferFunction = 'PiecewiseFunction'
inputDisplay.OpacityArray = ['POINTS', 'GlobalNodeId']
inputDisplay.OpacityTransferFunction = 'PiecewiseFunction'
inputDisplay.DataAxesGrid = 'GridAxesRepresentation'
inputDisplay.PolarAxes = 'PolarAxesRepresentation'
inputDisplay.ScalarOpacityFunction = vtkBlockColorsPWF
inputDisplay.ScalarOpacityUnitDistance = 0.2813007527538719
inputDisplay.OpacityArrayName = ['POINTS', 'GlobalNodeId']
inputDisplay.ExtractedBlockIndex = 2

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
inputDisplay.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 69056.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
inputDisplay.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 69056.0, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for vtkBlockColorsLUT in view renderView1
vtkBlockColorsLUTColorBar = GetScalarBar(vtkBlockColorsLUT, renderView1)
vtkBlockColorsLUTColorBar.Title = 'vtkBlockColors'
vtkBlockColorsLUTColorBar.ComponentTitle = ''

# set color bar visibility
vtkBlockColorsLUTColorBar.Visibility = 1

# show color legend
inputDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'CatalystTestImage_%.6ts%cm.png'
pNG1.Writer.ImageResolution = [844, 539]
pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.ExtractsOutputDirectory = 'catalyst_test_image_output'
options.GlobalTrigger = 'TimeStep'
options.CatalystLiveTrigger = 'TimeStep'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
