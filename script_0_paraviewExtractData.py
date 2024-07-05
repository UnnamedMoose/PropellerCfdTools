# trace generated using paraview version 5.11.0-RC1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
import os
import numpy as np
import pandas
import scipy.interpolate

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

Dprop = 0.34
radii = [r*Dprop/2. for r in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]]
outDir = "D:/git/propellerCfdTools"

# Read the prop surface.
propSurf = CGNSSeriesReader(
    registrationName='openWaterPropellerWithHub.cgns',
    FileNames=['D:/git/propellerCfdTools/openWaterPropellerWithHub.cgns'])
propSurf.Bases = ['Base_Surface_Elements']
propSurf.PointArrayStatus = ['Pressure', 'SkinFriction']

# Extract one blade, compute normals.
bladeSurfBlock = ExtractBlock(registrationName='bladeSurfBlock', Input=propSurf)
bladeSurfBlock.Selectors = ['/Root/Base_Surface_Elements/Zone_surface_2_2']
bladeSurf = ExtractSurface(registrationName='bladeSurf', Input=bladeSurfBlock)
bladeSurfNormals = GenerateSurfaceNormals(registrationName='bladeSurfNormals', Input=bladeSurf)
bladeSurfNormals.FlipNormals = 1

# Extract cuts along the blade.
# sliceData = []
for iRadius, radius in enumerate(radii):
    # Make a slice.
    slice = Slice(registrationName='Slice{:d}'.format(iRadius), Input=bladeSurfNormals)
    slice.SliceType = 'Cylinder'
    slice.SliceType.Axis = [1.0, 0.0, 0.0]
    slice.SliceType.Radius = radius

    # Save to a csv.
    sliceFile = 'sliceData_surface_r_{:.4f}.csv'.format(radius/(Dprop/2.))
    SaveData(os.path.join(outDir, sliceFile), proxy=slice,
        CellDataArrays=['Pressure', 'SkinFriction', 'Normals'],
        FieldDataArrays=['ispatch'],
        UseScientificNotation=1)

    # Load again.
    # sliceData.append(pandas.read_csv(os.path.join(outDir, sliceFile)))

"""
# # dataDir = "D:/propellerLPT_inception/blAnalysis/case_blAnalysis_J_0.600_rpm_411_f_0.500"
# dataDir = "D:/propellerLPT_inception/blAnalysis/case_blAnalysis_J_0.600_rpm_700_f_0.500"
# filename = "S5279prop-00000360_solution.cgns"
# outDir = "D:/propellerLPT_inception/blAnalysis/testData"
#
# # radii = [r*Dprop/2 for r in [0.7]]
# xcLocations = [0.03, 0.02, 0.01]
# propGeomIngesFile = "D:/propellerLPT_inception/geometry/S5279P.dat"
# dLine = 0.002


def readInges(filename):

    with open(filename, "r") as f:
        s = f.readlines()

    # Read global properties.
    _, D, scale, Dhub, nBlades = [float(v) for v in s[3].split()]
    # hubRatio = Dhub / D
    _, _, Pmean, BAR, _ = [float(v) for v in s[4].split()]
    nRadii, nSections = [int(v) for v in s[6].split()]

    # Read section data.
    sectionData = pandas.DataFrame(data=np.array([[float(v) for v in l.split()] for l in s[7:7+nRadii]]),
                                   columns=["r", "c", "skewLE", "skewMidChord", "pitch", "rake", "maxCamber",
                                            "xMaxCamber", "maxThickness", "xMaxThickness"])
    sectionData["r/R"] = sectionData["r"] / (D/2.)
    sectionData["c/D"] = sectionData["c"] / D
    sectionData["P/D"] = sectionData["pitch"] / D
    sectionData["rake/D"] = sectionData["rake"] / D
    sectionData["f/c"] = sectionData["maxCamber"] / np.clip(sectionData["c"], 1e-6, np.inf)
    sectionData["t/c"] = sectionData["maxThickness"] / np.clip(sectionData["c"], 1e-6, np.inf)
    # NOTE unsure about this one. It's close for the reference prop but not quite the same
    # as what came in the supposedly matching legacy ppg file.
    sectionData["skew"] = np.arctan(sectionData["skewMidChord"] / sectionData["r"]) / np.pi * 180

    # Read section shapes.
    sections = []
    for i in range(nRadii):
        i0 = 7+nRadii + i*nSections
        i1 = 7+nRadii + (i+1)*nSections
        sections.append(pandas.DataFrame(data=[[float(v) for v in l.split()] for l in s[i0:i1]],
                        columns=["x/c", "yBack/c", "yFace/c"]))
        sections[-1] /= max(1e-6, sectionData.loc[i, "c"])

    for f in sectionData:
        if f not in ["r/R", "c/D", "skew", "rake/D", "P/D", "t/c", "f/c"]:
            sectionData.drop(f, axis=1, inplace=True)

    return sectionData, sections


def formatVrmlHeader():
    return "#VRML V2.0 utf8\n\n"


def formatVrmlWireframe(pts, connections, colour=(0, 0, 0)):
    s = (
        "Shape {\n"
        "appearance Appearance {\n"
        "    material Material {\n"
        "        ambientIntensity 0\n"
        "        diffuseColor     " + "{:.4f} {:.4f} {:.4f}\n".format(colour[0], colour[1], colour[2]) +
        "        specularColor    1 0 0\n"
        "        emissiveColor    1 0 0\n"
        "        shininess        1\n"
        "        transparency     0\n"
        "    }\n"
        "}\n"
        "geometry IndexedLineSet {\n"
        "    coord Coordinate {\n"
        "        point [\n"
    )
    for pt in pts:
        s += "            {:.6e} {:.6e} {:.6e}\n".format(pt[0], pt[1], pt[2])
    s += (
        "        ]\n"
        "    }\n"
        "    coordIndex [\n"
    )
    for edge in connections:
        s += "            "
        s += " ".join(["{:d}".format(i) for i in edge]) + " -1,\n"
    s += (
        "        ]\n"
        "    }\n"
        "}\n"
    )
    return s

# create a new 'CGNS Series Reader'
s5279prop00000360_solutioncgns = CGNSSeriesReader(
    registrationName=filename,
    FileNames=[os.path.join(dataDir, filename)])
s5279prop00000360_solutioncgns.Bases = ['Base_Surface_elements', 'Base_Volume_elements']
s5279prop00000360_solutioncgns.CellArrayStatus = ['Pressure', 'SkinFriction', 'Velocity']

# Interpolate.
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=s5279prop00000360_solutioncgns)
cellDatatoPointData1.CellDataArraytoprocess = ['Pressure', 'SkinFriction', 'Velocity']

# create a new 'ExtractBlock' for the rotating cell zone.
extractBlock1 = ExtractBlock(registrationName='ExtractBlock1', Input=cellDatatoPointData1)
extractBlock1.Selectors = ['/Root/Base_Volume_elements/Interior_2']

# create a new 'ExtractBlock' for the blade.
extractBlock2 = ExtractBlock(registrationName='ExtractBlock2', Input=cellDatatoPointData1)
extractBlock2.Selectors = ['/Root/Base_Surface_elements/Zone_surface_2_2']
extractSurface1 = ExtractSurface(registrationName='ExtractSurface1', Input=extractBlock2)
generateSurfaceNormals1 = GenerateSurfaceNormals(registrationName='GenerateSurfaceNormals1', Input=extractSurface1)

# Extract cuts along the blade.
sliceFiles = []
for radius in radii:
    # slice1 = Slice(registrationName='Slice1', Input=extractBlock1)
    # slice1.SliceType = 'Cylinder'
    # slice1.SliceType.Axis = [1.0, 0.0, 0.0]
    # slice1.SliceType.Radius = radius
    # SaveData(os.path.join(outDir, 'sliceData_field_r_{:.3f}.cgns'.format(radius)), proxy=slice1)

    slice2 = Slice(registrationName='Slice2', Input=generateSurfaceNormals1)
    slice2.SliceType = 'Cylinder'
    slice2.SliceType.Axis = [1.0, 0.0, 0.0]
    slice2.SliceType.Radius = radius
    # SliceOffsetValues

    sliceFiles.append('sliceData_surface_r_{:.4f}.csv'.format(radius/(Dprop/2.)))
    SaveData(os.path.join(outDir, sliceFiles[-1]), proxy=slice2,
        CellDataArrays=['Pressure', 'SkinFriction', 'Velocity'],
        FieldDataArrays=['ispatch'],
        UseScientificNotation=1)

# Read the slices back to allow line probes to be created.

# Read the prop geometry in inges format
sectionData, sectionShapes = readInges(propGeomIngesFile)

# Vrml file for debugging and visualisation
outVrmlContents = formatVrmlHeader()

# start and end locations of the sampling lines
samplingLinesUpper = []
samplingLinesLower = []
# Coordinate system for each sampling point to resolve the flow into chordwise and spanwise.
lineCoordSystems = []

# Process each radius separately.
for iRadius, sliceFile in enumerate(sliceFiles):
    filename = os.path.join(outDir, sliceFile)

    # Read the slice data.
    surfData = pandas.read_csv(filename)

    # Get the radius of the slice.
    radius = np.mean(np.sqrt(surfData["Points:1"].values**2. + surfData["Points:2"].values**2.))
    r = radius / (Dprop/2.)

    # Find the leading and trailing edge. This is needed to sort the points correctly.
    coordFields = ["Points:0", "Points:1", "Points:2"]
    dMax = 0
    iEdges = []
    for i in range(surfData.shape[0]):
        dMaxTemp = 0.
        iEdgesTemp = None
        for j in range(surfData.shape[0]):
            d = np.linalg.norm(surfData[coordFields].values[i, :] - surfData[coordFields].values[j, :])
            if d > dMaxTemp:
                dMaxTemp = d
                iEdgesTemp = [i, j]
        if dMaxTemp > dMax:
            dMax = dMaxTemp
            iEdges = iEdgesTemp
    # Sort such that LE is first.
    # NOTE: in this case flow is in the +ve x-direction, so opposite to usual convention.
    if surfData.loc[iEdges[0], "Points:0"] > surfData.loc[iEdges[1], "Points:0"]:
        iEdges = (iEdges[1], iEdges[0])

    # Fit a plane to separate upper and lower surfaces.
    # NOTE: this will fail for large cambers but these are not typically used on propellers.
    x0 = np.array([(surfData["Points:0"].max() + surfData["Points:0"].min())/2, 0., 0.])
    v0 = surfData[coordFields].values[iEdges[0], :] - x0
    v1 = surfData[coordFields].values[iEdges[1], :] - x0
    vn = np.cross(v1, v0)
    vn /= np.linalg.norm(vn)

    surfData["isUpper"] = 0.
    for i in range(surfData.shape[0]):
        d = np.dot(vn, surfData[coordFields].values[i, :] - x0)
        if d >= 0.:
            surfData.loc[i, "isUpper"] = 1.

    # Get pper and lower surface points.
    iUpper = surfData[surfData["isUpper"] > 0].index.to_list()
    iLower = surfData[surfData["isUpper"] < 1e-3].index.to_list()
    for i in iEdges:
        if i not in iUpper:
            iUpper.append(i)
        if i not in iLower:
            iLower.append(i)

    # Sort along the curve. Start at the edge.
    def sortPoints(surfData, iEdges, indices):
        iSorted = [iEdges[0]]
        counter = 0
        while len(iSorted) < len(indices):
            d = np.linalg.norm(surfData[coordFields].values[indices, :] - surfData[coordFields].values[iSorted[-1], :], axis=1)
            for i in np.argsort(d):
                if indices[i] not in iSorted:
                    iSorted.append(indices[i])
                    break
            counter += 1
            if counter > len(indices)*10:
                raise RuntimeError(":*(")
        indices = iSorted
        return indices

    iUpper = sortPoints(surfData, iEdges, iUpper)
    iLower = sortPoints(surfData, iEdges, iLower)
    upperSurf = surfData.iloc[iUpper].copy().reset_index(drop=True)
    lowerSurf = surfData.iloc[iLower].copy().reset_index(drop=True)

    # Compute the section shape from prop geometry at this radius.
    for i in range(sectionData.shape[0]-1):
        if (r >= sectionData.loc[i, "r/R"]) and (r < sectionData.loc[i+1, "r/R"]):
            rGradMult = 1./(sectionData.loc[i+1, "r/R"] - sectionData.loc[i, "r/R"])*(r - sectionData.loc[i, "r/R"])
            # Get the section properties
            sectionProperties = (sectionData.iloc[i+1] - sectionData.iloc[i])*rGradMult + sectionData.iloc[i]
            # Get the section shape.
            # NOTE assuming x/c are the same!
            sectionShape = ((sectionShapes[i+1] - sectionShapes[i])*rGradMult + sectionShapes[i]).values
            sectionShape[:, 0] = sectionShapes[i]["x/c"]
            splineUpper = scipy.interpolate.InterpolatedUnivariateSpline(
                sectionShape[:, 0], sectionShape[:, 1], k=2)
            splineLower = scipy.interpolate.InterpolatedUnivariateSpline(
                sectionShape[:, 0], sectionShape[:, 2], k=2)
            # pitch angle in radians for the nose-tail line
            sectionProperties["thetaNt"] = np.arctan(sectionProperties["P/D"]/(sectionProperties["r/R"]*np.pi))
            break

    def pointOnProp(sectionProperties, xc, yc, Dprop):
        #Helper function which maps a 2D slice onto a 3D propeller radial offset
        #See John Carlton's Marine Propellers and Propulsion, eq. 3.17 on pp 47 (2nd ed.)
        # note that x/R is divided by 2 to yield x/D to be consistent with overall non-dimensionalisation scheme
        xcGenerator = 0.5
        x = -(sectionProperties["rake/D"] + sectionProperties["r/R"]/2*sectionProperties["skew"]*np.tan(sectionProperties["thetaNt"])) \
            + (xcGenerator-xc)*sectionProperties["c/D"]*np.sin(sectionProperties["thetaNt"]) \
            + yc*sectionProperties["c/D"]*np.cos(sectionProperties["thetaNt"])

        y = sectionProperties["r/R"]/2*np.sin(sectionProperties["skew"] \
            - ((xcGenerator - xc)*sectionProperties["c/D"]*np.cos(sectionProperties["thetaNt"]) \
            - yc*sectionProperties["c/D"]*np.sin(sectionProperties["thetaNt"]))/sectionProperties["r/R"]*2)

        z = sectionProperties["r/R"]/2*np.cos(sectionProperties["skew"] \
            - ((xcGenerator - xc)*sectionProperties["c/D"]*np.cos(sectionProperties["thetaNt"]) \
            - yc *sectionProperties["c/D"]*np.sin(sectionProperties["thetaNt"]))/sectionProperties["r/R"]*2)

        # NOTE the minus signs because the flow is in the +ve x-direction!
        return np.array([-x, -y, z])*Dprop

    # # NOTE the sides being reversed because the flow is in the +ve x-direction!
    # s = splineUpper
    # splineUpper = splineLower
    # splineLower = s

    # Test the offsets.
    xcRebuild = np.linspace(0, 1, 101)
    lowerSurfRemade = np.array([pointOnProp(sectionProperties, xc, splineLower(xc), Dprop) for xc in xcRebuild])
    upperSurfRemade = np.array([pointOnProp(sectionProperties, xc, splineUpper(xc), Dprop) for xc in xcRebuild])
    outVrmlContents += formatVrmlWireframe(lowerSurfRemade, [(i, i+1) for i in range(lowerSurfRemade.shape[0]-1)], colour=(0, 0, 1))
    outVrmlContents += formatVrmlWireframe(upperSurfRemade, [(i, i+1) for i in range(upperSurfRemade.shape[0]-1)], colour=(1, 0, 0))

    # TODO there is probably a better way to obtain the normal vectors but hey ho.
    #   Could compute a section at a small r offset and then get dx/dr, dy/dr, dz/dr and
    #   then also dx/dxc, dy/dxc, dz/dxc and somehow work from there.
    #   Alternatively, could compue 4 points around the point in question (+/- delta-xc, +/- delta_r)
    #   and compute four cross products, and then average the normal vectors
    def orientSamplingLine(xc, basePoint, tBlade, spline, surf, isUpper):
        global outVrmlContents

        # Find the nearest point in the grid to get the normal.
        iNearest = np.argmin(np.linalg.norm(surf[coordFields].values - basePoint, axis=1))
        normal = -1.*surf[["Normals:0", "Normals:1", "Normals:2"]].values[iNearest, :]
        # Offset base point along the normal to get the end point.
        endPoint = basePoint + normal*dLine
        # Start the line a tiny bit inwards to make sure it intersects the modelled blade surface.
        startPoint = basePoint - normal*0.2*tBlade

        # Compute the tangential direction.
        vTan = pointOnProp(sectionProperties, min(1., xc+1e-4), spline(xc), Dprop) \
            - pointOnProp(sectionProperties, max(0., xc-1e-4), spline(xc), Dprop)
        vTan /= np.linalg.norm(vTan)

        # Spanwise is orthogonal to the other two directions. Orient such that it
        # points towards the tip.
        if isUpper:
            vSpan = np.cross(normal, vTan)
        else:
            vSpan = np.cross(vTan, normal)

        # concatenate.
        line = np.vstack([startPoint, endPoint])

        # Add to visualisation.
        outVrmlContents += formatVrmlWireframe(line, [(0, 1)], colour=(1, 0.65, 0))
        outVrmlContents += formatVrmlWireframe([basePoint, basePoint+vTan*dLine], [(0, 1)], colour=(0.62, 0.12, 0.94))
        outVrmlContents += formatVrmlWireframe([basePoint, basePoint+vSpan*dLine], [(0, 1)], colour=(0, 1, 1))

        # Save the coordinate system for analysis of extracted data.
        lineCoordSystems.append({
            "r": r,
            "x/c": xc,
            "tBlade": tBlade,
            "isUpper": isUpper,
            "x": basePoint[0],
            "y": basePoint[1],
            "z": basePoint[2],
            "nx": normal[0],
            "ny": normal[1],
            "nz": normal[2],
            "tx": vTan[0],
            "ty": vTan[1],
            "tz": vTan[2],
            "sx": vSpan[0],
            "sy": vSpan[1],
            "sz": vSpan[2],
        })

        return line

    for iXc, xcSampling in enumerate(xcLocations):
        # Compute point on the blade where sampling is requested.
        basePointUpper = pointOnProp(sectionProperties, xcSampling, splineUpper(xcSampling), Dprop)
        basePointLower = pointOnProp(sectionProperties, xcSampling, splineLower(xcSampling), Dprop)
        localBladeThickness = np.linalg.norm(basePointLower - basePointUpper)

        # Orient the sampling lines and store them
        samplingLinesUpper.append(orientSamplingLine(
            xcSampling, basePointUpper, localBladeThickness, splineUpper, upperSurf, True))
        samplingLinesLower.append(orientSamplingLine(
            xcSampling, basePointLower, localBladeThickness, splineLower, lowerSurf, False))

        fieldsToSample = ['Pressure', 'Velocity', 'arc_length', 'vtkValidPointMask']

        # Create line monitors using Resample to Line
        resampleToLine1 = ResampleToLine(
            registrationName='ResampleToLine_{:d}'.format(len(xcLocations)*iRadius+iXc),
            Input=extractBlock1)
        resampleToLine1.Point1 = samplingLinesUpper[-1][0, :]
        resampleToLine1.Point2 = samplingLinesUpper[-1][1, :]
        resampleToLine1.Resolution = 100
        SaveData(os.path.join(outDir, "lineUpper_r_{:.4f}_xc_{:.4f}.csv".format(r, xcSampling)),
            proxy=resampleToLine1, PointDataArrays=fieldsToSample, UseScientificNotation=1)

        resampleToLine2 = ResampleToLine(
            registrationName='ResampleToLine_{:d}'.format(len(xcLocations)*iRadius+iXc),
            Input=extractBlock1)
        resampleToLine2.Point1 = samplingLinesLower[-1][0, :]
        resampleToLine2.Point2 = samplingLinesLower[-1][1, :]
        resampleToLine2.Resolution = 100
        SaveData(os.path.join(outDir, "lineLower_r_{:.4f}_xc_{:.4f}.csv".format(r, xcSampling)),
            proxy=resampleToLine2, PointDataArrays=fieldsToSample, UseScientificNotation=1)

        # Create the line monitors using Plot over line.
        # # samplingType = 'Sample At Cell Boundaries'  # Two points on either side of each face
        # samplingType = 'Sample At Segment Centers'  # One point in between two crossed faces.
        # # samplingType = 'Sample Uniformly'
        # plotOverLineUpper = PlotOverLine(
        #     registrationName='PlotOverLineUpper_{:d}'.format(len(xcLocations)*iRadius+iXc), Input=extractBlock1)
        # plotOverLineUpper.Point1 = samplingLinesUpper[-1][0, :]
        # plotOverLineUpper.Point2 = samplingLinesUpper[-1][1, :]
        # plotOverLineUpper.SamplingPattern = samplingType
        # if samplingType == "SampleUniformly":
        #     plotOverLineUpper.Resolution = 100

        # SaveData(os.path.join(outDir, "lineUpper_r_{:.4f}_xc_{:.4f}.csv".format(r, xcSampling)),
        #     proxy=plotOverLineUpper, PointDataArrays=fieldsToSample, UseScientificNotation=1)

        # plotOverLineLower = PlotOverLine(
        #     registrationName='PlotOverLineLower_{:d}'.format(len(xcLocations)*iRadius+iXc), Input=extractBlock1)
        # plotOverLineLower.Point1 = samplingLinesLower[-1][0, :]
        # plotOverLineLower.Point2 = samplingLinesLower[-1][1, :]
        # plotOverLineLower.SamplingPattern = samplingType
        # if samplingType == "SampleUniformly":
        #     plotOverLineLower.Resolution = 100
        # SaveData(os.path.join(outDir, "lineLower_r_{:.4f}_xc_{:.4f}.csv".format(r, xcSampling)),
        #     proxy=plotOverLineLower, PointDataArrays=fieldsToSample, UseScientificNotation=1)

# Save the vrml file.
with open(os.path.join(outDir, "analysisVis.vrml"), "w") as outf:
    outf.write(outVrmlContents)

# Save blade properties in human-readable format.
sectionData.to_csv(os.path.join(outDir, "bladeProperties.csv"), index=False)
# Save the coordinate systems for analysis.
lineCoordSystems = pandas.DataFrame(lineCoordSystems)
lineCoordSystems.to_csv(os.path.join(outDir, "lineCoordinateSystems.csv"), index=False)
"""
