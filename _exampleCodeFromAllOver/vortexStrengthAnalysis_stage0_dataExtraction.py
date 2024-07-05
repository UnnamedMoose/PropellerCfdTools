# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 22:31:58 2020

@author: ALidtke
"""

import paraview.simple as pv
from paraview import numpy_support as ns
import numpy as np
import os
import shutil
import pandas
import re

# In the future, this might get updated somehow.
projectDir = os.getcwd()

# %% Constants.
filename = 'ORCA_solution.cgns'  # ReFRESCO output name.
guiderailFilename = os.path.join(projectDir, "pointCloud.ipcl")  # Name of file with centres of each slice (indexed point cloud from Rhino)
wallIndices = [15, 16, 17, 18, 19, 20]  # Which BC indices form the prop and hub (use Start trace -> Extract block to check).
surfaceTriangulationReduction = 0.95  # How much to coarsen stl of the geometry.
sliceRadiusByD = 0.08  # Size of each slice.
nStride = 1  # Use every n-th slice location from the file.

# %% Get constants from files, derive secondary properties, etc.
def getKeyword(key, controlFileString):
    return re.findall('<{}[\s?name="]?.*>.*</{}>'.format(key, key), controlFileString)[0].split(">")[1].split("</")[0]

with open(os.path.join(projectDir, "controls.xml")) as infile:
    controls = infile.read()
Dprop = float(getKeyword("referenceLength", controls))
Uref = float(getKeyword("referenceVelocity", controls))
rho = float(getKeyword("density", controls))
rClip = sliceRadiusByD * Dprop

# %% Read the guide rail.
pts = np.loadtxt(guiderailFilename, skiprows=1, usecols=[1,2,3])

# Compute normals.
normals = np.zeros(pts.shape)
normals[:-1,:] = pts[1:,:] - pts[:-1,:]
normals[-1,:] = pts[-1,:] - pts[-2,:]
normals /= np.linalg.norm(normals, axis=1)[:,np.newaxis]

# Subset to a manageable number of points.
pts = pts[::nStride,:]
normals = normals[::nStride,:]

# %% Read data and apply filters.
pv._DisableFirstRenderCameraReset()

# Read the file.
solution = pv.CGNSSeriesReader(registrationName=filename, FileNames=[os.path.join(projectDir, filename)])
solution.Bases = ['Base_Surface_elements', 'Base_Volume_elements']
solution.CellArrayStatus = ['Pressure', 'SkinFriction', 'Velocity', 'Vorticity', 'Yplus']

# create a new 'Extract Block'
extractBlock1 = pv.ExtractBlock(Input=solution)

# Properties modified on extractBlock1
extractBlock1.BlockIndices = wallIndices

# get active view
renderView1 = pv.GetActiveViewOrCreate('RenderView')

# get layout
layout1 = pv.GetLayout()

# show data in view
extractBlock1Display = pv.Show(extractBlock1, renderView1, 'UnstructuredGridRepresentation')

# hide data in view
pv.Hide(solution, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Merge Blocks'
mergeBlocks1 = pv.MergeBlocks(Input=extractBlock1)

# show data in view
mergeBlocks1Display = pv.Show(mergeBlocks1, renderView1, 'UnstructuredGridRepresentation')

# hide data in view
pv.Hide(extractBlock1, renderView1)

# update the view to ensure updated data information
# renderView1.Update()

# create a new 'Extract Surface'
extractSurface1 = pv.ExtractSurface(Input=mergeBlocks1)

# show data in view
extractSurface1Display = pv.Show(extractSurface1, renderView1, 'GeometryRepresentation')

# hide data in view
pv.Hide(mergeBlocks1, renderView1)

# update the view to ensure updated data information
# renderView1.Update()

# create a new 'Triangulate'
triangulate1 = pv.Triangulate(Input=extractSurface1)

# show data in view
triangulate1Display = pv.Show(triangulate1, renderView1, 'GeometryRepresentation')

# hide data in view
pv.Hide(extractSurface1, renderView1)

# update the view to ensure updated data information
# renderView1.Update()

# create a new 'Decimate'
decimate1 = pv.Decimate(Input=triangulate1)

# Properties modified on decimate1
decimate1.TargetReduction = surfaceTriangulationReduction

# show data in view
decimate1Display = pv.Show(decimate1, renderView1, 'GeometryRepresentation')

# hide data in view
pv.Hide(triangulate1, renderView1)

# save surface data
pv.SaveData(os.path.join(projectDir, "wallGrid.stl"), proxy=decimate1, FieldDataArrays=['ispatch'])
pv.SaveData(os.path.join(projectDir, "wallGrid.obj"), proxy=decimate1, FieldDataArrays=['ispatch'])

# %% Extract and save all the slices.

# Interpolate data onto grid vertices.
cellDatatoPointData1 = pv.CellDatatoPointData(Input=solution)

# Set inputs for slice extraction.
source = pv.FindSource('CellDatatoPointData1')
renderView1 = pv.GetActiveViewOrCreate('RenderView')

# Prepare outputs.
targetDir = os.path.join(projectDir, "sliceDataForVortexAnalysis")
if os.path.isdir(targetDir):
    shutil.rmtree(targetDir)
os.mkdir(targetDir)

# Act on each point.
calcs = []
for i in range(pts.shape[0]):
	x0 = pts[i,:]
	nHat = normals[i,:]

	# Create the slice.
	slays = pv.Slice(source)
	slays.SliceType = 'Plane'
	slays.SliceType.Origin = x0
	slays.SliceType.Normal = nHat

	# Subset.
	cleep = pv.Clip(Input=slays)
	cleep.ClipType = 'Cylinder'
	cleep.ClipType.Center = x0
	cleep.ClipType.Axis = nHat
	cleep.ClipType.Radius = rClip

	# Compute Cp.
	mb = pv.MergeBlocks(cleep)
	calcs.append(pv.Calculator(Input=mb))
	calcs[-1].ResultArrayName = 'Cp'
	calcs[-1].Function = 'Pressure / (0.5*{:.6e}^2*{:.6e})'.format(Uref, rho)

	# Save.
	pv.SaveData(os.path.join(targetDir, "slice_{:d}.csv".format(i)), proxy=calcs[-1])
	pv.SaveData(os.path.join(targetDir, "slice_{:d}.cgns".format(i)), proxy=calcs[-1])

# Set an orthogonal coordinate system on each slice.
# TODO figure out how to define the second direction such that it makes the most sense.
#   Now it's aligned with the global x-axis so that it's at least consisten for all slices,
#   but this doesn't have to hold true in all cases.
directionVectors_j = []
directionVectors_k = []
for i in range(pts.shape[0]):
    # Read the slice again.
    df = pandas.read_csv(os.path.join(targetDir, "slice_{:d}.csv".format(i)))
    # First axis aligned with normal
    iHat = normals[i,:]
    # Second aligned in the global x-direction (or as close as possible)
    # jHat = np.array([1., 0., 0.])
    jHat = df.loc[df["Points:0"].argmax(), ["Points:0", "Points:1", "Points:2"]].values - pts[i,:]
    jHat /= np.linalg.norm(jHat)
    kHat = np.cross(iHat, jHat)
    directionVectors_j.append(jHat)
    directionVectors_k.append(kHat)

    # Check coord system.
    with open(os.path.join(targetDir, "sliceCoordSystem_{:d}.obj".format(i)), "w") as f:
        xs = [
            pts[i, :],
            pts[i, :] + iHat*rClip,
            pts[i, :] + jHat*rClip,
            pts[i, :] + kHat*rClip,
        ]
        for x in xs:
            f.write("v {:.6e} {:.6e} {:.6e}\n".format(x[0], x[1], x[2]))
        f.write("l 1 2\n")
        f.write("l 1 3\n")
        f.write("l 1 4\n")

directionVectors_j = np.array(directionVectors_j)
directionVectors_k = np.array(directionVectors_k)

# Concatenate and save slice metadata as data frame.
df = pandas.DataFrame(data=np.hstack([pts, normals, normals, directionVectors_j, directionVectors_k]),
                      columns=[pattern.format(c) for pattern in ["{}", "n{}", "ei_{}", "ej_{}", "ek_{}"] for c in ["x", "y", "z"]])
df.to_csv(os.path.join(targetDir, "sliceProperties.csv"), index=False)

# Display.
# for calc in calcs:
# 	calc = pv.Show(calc, renderView1, 'UnstructuredGridRepresentation')

# Update the view.
# pv.Hide(source, renderView1)
# renderView1.Update()

