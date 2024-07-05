# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 23:13:02 2020

@author: ALidtke
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas
import itertools
import os
import re
import sys
import scipy.interpolate
import shutil


# %% Function for picking up unique integer combinations.
def getUniqueEdges(faces):
    edgeVertexIndices = []
    faceEdgeIndices = np.zeros((len(faces), 3), dtype=np.int)
    for iFace in range(len(faces)):
        # Go over each edge of this face, assuming right-hand rule for ordering.
        for i in range(3):
            # Check which vertices making up this face are being paired together.
            iv0 = i
            iv1 = (i+1) % 3
            # Get actual vertex indices.
            iv0 = faces[iFace, iv0]
            iv1 = faces[iFace, iv1]
            # Check the direction of the edge. Always define it as (smaller, larger).
            direction = 1
            if iv0 > iv1:
                direction = -1
                iv0, iv1 = iv1, iv0
            edge = (iv0, iv1)
            # Check if this edge has been stored already. If not, add it to the edge list. If yes, just refer
            # this face to its location in the edge list, indicating direction.
            try:
                edgeIdx = edgeVertexIndices.index(edge)
            except ValueError:
                edgeVertexIndices.append(edge)
                edgeIdx = len(edgeVertexIndices) - 1
            faceEdgeIndices[iFace, i] = edgeIdx * direction
    nEdges = len(edgeVertexIndices)
    # Convert to array, no longer need a list of tuples for checking for duplicates.
    edgeVertexIndices = np.array(edgeVertexIndices)
    return nEdges, edgeVertexIndices, faceEdgeIndices


# %% Wall class used for IO and storing information about a triangulated geometry.
class Wall(object):
    def __init__(self, filename, scale=1):

        # - Read the vertices and triangulated faces from an .obj file. -
        vertices, faceVertexIndices = [], []

        with open(filename, "r") as f:
            for l in f:
                if re.match("v .*", l):
                    vertices.append([float(v) for v in l.replace("v", "").split()])
                elif re.match("f .*", l):
                    # Distinguish between simple .obj files and ones which explicityl define face no. (as in ReFRESCO).
                    if "//" not in l:
                        faceVertexIndices.append([(int(v)-1) for v in l.replace("f", "").split()])
                    else:
                        faceVertexIndices.append([(int(v)-1) for v in re.sub("//[0-9]+", "", l.replace("f", "")).split()])

        self.vertices = np.array(vertices) * scale
        self.faceVertexIndices = np.array(faceVertexIndices, dtype=np.int)
        self.nVertices = self.vertices.shape[0]
        self.nFaces = self.faceVertexIndices.shape[0]

        # - Allocate the bounding box. -
        self.bBox = np.zeros((2, 3))

        # - Identify point pairs making up edges. -
        self.nEdges, self.edgeVertexIndices, self.faceEdgeIndices = getUniqueEdges(self.faceVertexIndices)

        # - Allocate edge vectors and face properties. -
        self.edges = np.zeros((self.nEdges, 3))
        self.faceNormals = np.zeros((self.nFaces, 3))
        self.faceCentres = np.zeros((self.nFaces, 3))
        self.faceAreas = np.zeros(self.nFaces)

        # - Comptue geometric properties. -
        self.computeProperties()
        # TODO distributed normals?

    def computeProperties(self):
        # - Find the bounding box. -
        self.bBox[0, :] = np.min(self.vertices, axis=0)
        self.bBox[1, :] = np.max(self.vertices, axis=0)
        # Compute edge vectors.
        for iEdge in range(self.nEdges):
            self.edges = self.vertices[self.edgeVertexIndices[iEdge, 1], :] \
                - self.vertices[self.edgeVertexIndices[iEdge, 0], :]
        # Vectors along two edges in array form for the entire geometry in one go.
        e0 = self.vertices[self.faceVertexIndices[:, 1], :] - self.vertices[self.faceVertexIndices[:, 0], :]
        e1 = self.vertices[self.faceVertexIndices[:, 2], :] - self.vertices[self.faceVertexIndices[:, 0], :]
        # Compute normals using eo x e1
        self.faceNormals = np.cross(e0, e1, axis=1)
        # Compute area.
        self.faceAreas = 0.5 * scipy.linalg.norm(self.faceNormals, axis=1)
        # Normalise normals.
        for i in range(3):
            self.faceNormals[:, i] /= (2. * self.faceAreas)
        # Compute face centres.
        self.faceCentres = (
            self.vertices[self.faceVertexIndices[:, 0], :]
            + self.vertices[self.faceVertexIndices[:, 1], :]
            + self.vertices[self.faceVertexIndices[:, 2], :]
            ) / 3.

        self.centroid = np.mean(self.vertices, axis=0)

    def isPointInsideBoundingBox(self, point, offset):
        if (
                (point[0] < (self.bBox[0, 0]-offset)) or
                (point[1] < (self.bBox[0, 1]-offset)) or
                (point[2] < (self.bBox[0, 2]-offset)) or
                (point[0] > (self.bBox[1, 0]+offset)) or
                (point[1] > (self.bBox[1, 1]+offset)) or
                (point[2] > (self.bBox[1, 2]+offset))):
            return False
        else:
            return True

    def plot(self, ax, plotPoints=True, plotEdges=True, plotFaces=False, faceColour="red", faceAlpha=1.0,
             facesInclude=None, pointSize=6, edgesInclude=None, edgeWidth=1, pointEdgeColour="k"):

        if plotPoints:
            ax.plot(self.vertices[:, 0], self.vertices[:, 1], self.vertices[:, 2], ".", c=pointEdgeColour, ms=pointSize)

        if edgesInclude is None:
            edgesInclude = range(self.nEdges)
        if plotEdges:
            for iEdge in edgesInclude:
                pts = self.vertices[self.edgeVertexIndices[iEdge, :], :]
                ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], "-", lw=edgeWidth, c=pointEdgeColour)

        if facesInclude is None:
            facesInclude = range(self.nFaces)
        if plotFaces:
            for iFace in facesInclude:
                ax.add_collection3d(Poly3DCollection([list(self.vertices[self.faceVertexIndices[iFace, :], :])],
                                                     facecolors=faceColour, linewidths=0, edgecolors="None",
                                                     alpha=faceAlpha))


# %% Functions for dealing with triangulated geometries.

def checkIfOnFace(point, vertices, normal, area, normDistTol=1e-4, edgeTolNoProj=1e-2):
    e0 = vertices[1, :] - vertices[0, :]
    e1 = vertices[2, :] - vertices[0, :]

    # Check the projected distance from the face plane.
    dist = np.dot(point - vertices[0, :], normal)

    # Assume that face area is unchanged (rigid body motion).
    if ((np.abs(dist) / np.sqrt(area)) > normDistTol):
        return False, (0., 0.)

    # Use an algorithm that does not use projections, which should safeguard against degenerate 2D
    # triangles. Taken from: http://geomalgorithms.com/a06-_intersect-2.html

    # Compute all the vectors and dot products.
    uu = np.dot(e0, e0)
    uv = np.dot(e0, e1)
    vv = np.dot(e1, e1)
    w = point - vertices[0, :]
    wu = np.dot(w, e0)
    wv = np.dot(w, e1)
    D = uv * uv - uu * vv

    # Get and test parametric coords
    s = (uv * wv - vv * wu) / D
    if ((s >= -edgeTolNoProj) and (s <= (1.0 + edgeTolNoProj))):
        t = (uv * wu - uu * wv) / D
        if ((t >= -edgeTolNoProj) and ((s + t) <= (1.0 + edgeTolNoProj))):
            return True, (s, t)

    return False, (0., 0.)


def findRayIntersection(p0, p1, vertices, normal, area):
    v_P0_V0 = vertices[0, :] - p0
    v_P0_P1 = p1 - p0

    # See of the line is parrallel to the face.
    validContact = False
    contactPoint = np.zeros(3)
    fPlane = -1.0
    if (np.abs(np.dot(v_P0_P1/np.linalg.norm(v_P0_P1), normal)) > 1e-6):
        fPlane = np.dot(v_P0_V0, normal) / np.dot(v_P0_P1, normal)
        # If this is bound within <0, 1>, intersection occurs between p1 and p0.
        if ((fPlane >= 0.0) and (fPlane <= 1.0)):
            # Valid intersection with the plane that makes up this face, this occurrs at:
            contactPoint = p0 + fPlane*v_P0_P1
            # Test if the plane intersection point is inside the triangle.
            validPoint, triangleCoords = checkIfOnFace(contactPoint, vertices, normal, area)
            if validPoint:
                validContact = True

    return validContact, contactPoint, fPlane


def checkIfInside(point, vertices, faces, normals, areas, vertexTol=1e-8, pOut=None,
                  returnNoIntersections=False):
    # Pick a reference point outside.
    if pOut is None:
        pOut = 2*np.max(vertices, axis=0) - np.min(vertices, axis=0)

    # this will hold the last found interseciton of vPO and a face; multiple matches
    # may occur if the intersection point lands right on top of a vertex, but those
    # still should count as a single match in terms of determining whether or not
    # the enclosed point is inside the stl. Initialise as pOut to make sure the
    # first found intersection will not get discarded based on distance from
    # the saved validIntersection value.
    lastValidIntersection = pOut

    # initialise the key counter
    nIntersections = 0

    for iFace in range(faces.shape[0]):
        # Check if a ray crossing between the two points passes through this face.
        validIntersection, pIntersect, fIntersect = findRayIntersection(
            point, pOut, vertices[faces[iFace, :]], normals[iFace, :], areas[iFace])

        # If this intersection is the same as the previous one, ignore it;
        # this will happen when pIntersect is on top of a vertex shared by multiple faces.
        if ((validIntersection) and (np.linalg.norm(lastValidIntersection - pIntersect) < vertexTol)):
            validIntersection = False
        elif (validIntersection):
            lastValidIntersection = pIntersect
            nIntersections = nIntersections + 1

    if returnNoIntersections:
        return nIntersections
    else:
        if (nIntersections > 0) and ((nIntersections % 2) != 0):
            return True
        else:
            return False


# %% Constants.
#case = "./"
#case = "Orca1"
case = "X:\P31879.709\OptimisationCases\CompleteOptimisation\Optimisation_RANS\Orca9\prop\Simulation"

dataDir = os.path.join(case, "sliceDataForVortexAnalysis")
outDataDir = os.path.join(case, "processedSliceData")
saveFigs = True
plotOverviewFigs = True
plotDetailedFigs = True
includeCavitation = False  # Here be dragons - not tested with the new version.

interpMethod = "nearest"
approxMaxCoreRadiusByD = 0.05
nGridPointsSubsample = 101  # 51
nSectors = 13  # 13
nRadii = 201  # 101
vofThreshold = 0.05

# %% Read controls.
def getKeyword(key, controlFileString):
    return re.findall('<{}[\s?name="]?.*>.*</{}>'.format(key, key), controlFileString)[0].split(">")[1].split("</")[0]

with open(os.path.join(case, "controls.xml")) as infile:
    controls = infile.read()
Dprop = float(getKeyword("referenceLength", controls))
Uref = float(getKeyword("referenceVelocity", controls))
rho = float(getKeyword("density", controls))
approxMaxCoreRadius = approxMaxCoreRadiusByD * Dprop

# Make room for outputs.
if os.path.isdir(outDataDir):
    shutil.rmtree(outDataDir)
os.mkdir(outDataDir)

# %% Read data.
sliceData = []
sliceFiles = [f for f in os.listdir(dataDir) if re.match("slice_[0-9]+.csv", f)]
sliceNumbers = [int(f.split("_")[-1].replace(".csv", "")) for f in sliceFiles]
sliceFiles = np.array(sliceFiles)[np.argsort(sliceNumbers)]
for s in sliceFiles:
    sliceData.append(pandas.read_csv(os.path.join(dataDir, s)))
propStl = Wall(os.path.join(case, "wallGrid.obj"))

# %% Read slice properties.
sliceProperties = pandas.read_csv(os.path.join(dataDir, "sliceProperties.csv"))

# %% Resolve vectors onto the slice coordinate system.
for i in range(len(sliceData)):
    # First axis aligned with normal
    iHat = sliceProperties.loc[i, ["ei_x", "ei_y", "ei_z"]].values
    jHat = sliceProperties.loc[i, ["ej_x", "ej_y", "ej_z"]].values
    kHat = sliceProperties.loc[i, ["ek_x", "ek_y", "ek_z"]].values
    centre = sliceProperties.loc[i, ["x", "y", "z"]].values
    dirVecs = [iHat, jHat, kHat]
    # Resolve vectors onto the new coordinate system.
    for j, crd in enumerate(["Z", "X", "Y"]):
        sliceData[i]["planeCoord{}star".format(crd)] = \
            np.dot(sliceData[i][["Points:0", "Points:1", "Points:2"]].values - centre, dirVecs[j])
        sliceData[i]["planeVel{}star".format(crd)] = \
            np.dot(sliceData[i][["Velocity:0", "Velocity:1", "Velocity:2"]].values, dirVecs[j])
        sliceData[i]["planeVort{}star".format(crd)] = \
            np.dot(sliceData[i][["Vorticity:0", "Vorticity:1", "Vorticity:2"]].values, dirVecs[j])
    # Compute velocity and vorticity magnitudes.
    sliceData[i]["magVelPlane"] = \
        np.linalg.norm(sliceData[i][["planeVelXstar", "planeVelYstar"]].values, axis=1)
    sliceData[i]["magVort"] = \
        np.linalg.norm(sliceData[i][["Vorticity:0", "Vorticity:1", "Vorticity:2"]].values, axis=1)

# %% Plot each slice for figuring out what is going on.
if plotDetailedFigs:
    fields = ["planeVortZstar", "Cp", "magVort", "magVelPlane"]
    lims = [(-100, 100), (-1, 1), (0, 100), (0, 5)]
    labels = ["\omega_Z* [s^{-1}]", "C_P [-]", "|\omega| [1/s]", "|U_{plane}| [m/s]"]
    if includeCavitation:
        fields.append("VaporVolumeFraction")
        lims.append((0, 1))
        labels.append("\\alpha [-]")
    for i in range(len(sliceData)):
        for j in range(len(fields)):
            fig, ax = plt.subplots(1)
            ax.set_xlabel("x* [m]")
            ax.set_ylabel("y* [m]")
            ax.set_aspect("equal", "box")
            sct = ax.scatter(sliceData[i]["planeCoordXstar"].values, sliceData[i]["planeCoordYstar"].values,
                    c=sliceData[i][fields[j]], s=25, linewidths=0,
                    vmin=lims[j][0], vmax=lims[j][1])
            plt.colorbar(sct)
            if saveFigs:
                plt.savefig(os.path.join(outDataDir, "sliceData_{}_{:d}.png".format(fields[j], i)), dpi=200, bbox_inches="tight")

#%% Perform vortex identification.

vortexData = {k:[] for k in ["CoordinateX", "CoordinateY", "CoordinateZ",
              "coreRadius", "maxTangentialVelocity", "circulation", "cavityRadius"]}

# TODO this loop can be parallelised.
#for iSlice in [10]:
for iSlice in range(len(sliceData)):
    x = sliceData[iSlice]["planeCoordXstar"].values
    y = sliceData[iSlice]["planeCoordYstar"].values
    u = sliceData[iSlice]["planeVelXstar"].values
    v = sliceData[iSlice]["planeVelYstar"].values
    cp = sliceData[iSlice]["Cp"].values
    try:
        vof = sliceData[iSlice]["VaporVolumeFraction"].values
    except KeyError:
        vof = np.zeros(len(cp))
    magU = np.sqrt(u**2 + v**2)
    iMinCp = np.argmin(sliceData[iSlice]["Cp"].values)
    iMinVort = np.argmin(sliceData[iSlice]["planeVortZstar"].values)

    xg = np.linspace(np.min(x), np.max(x), nGridPointsSubsample)
    yg = np.linspace(np.min(y), np.max(y), nGridPointsSubsample)
    xg, yg = np.meshgrid(xg, yg)
    ug = scipy.interpolate.griddata(np.vstack((x, y)).T, u/magU, (xg, yg), method=interpMethod)
    vg = scipy.interpolate.griddata(np.vstack((x, y)).T, v/magU, (xg, yg), method=interpMethod)
    cpg = scipy.interpolate.griddata(np.vstack((x, y)).T, cp, (xg, yg), method=interpMethod)
    vofg = scipy.interpolate.griddata(np.vstack((x, y)).T, vof, (xg, yg), method=interpMethod)

    # Compute velocity angle.
    alphag = np.angle(ug + 1j*vg, deg=True)
    alphag[alphag<0] = alphag[alphag<0] + 360.

    # Perform vortex identification
    if not includeCavitation:
        # Divide into sectors.
        sectors = np.linspace(0, 360, nSectors)
        beta = np.zeros(alphag.shape)
        for j in range(1, len(sectors)):
            beta[(alphag < sectors[j]) & (alphag >= sectors[j-1])] = j

        # Find populated sectors.
        populatedSectors = np.unique(beta)
        populatedSectors = populatedSectors[np.abs(populatedSectors) > 1e-6]

        # For each point, find nearest points in other sectors.
        # lnValues = np.zeros(xg.shape)
        lnValuesV2_cw = np.zeros(xg.shape)
        # lnValuesV2_acw = np.zeros(xg.shape)

        for i in range(xg.shape[0]):
            for j in range(xg.shape[1]):
                # Distance between this point and all the rest.
                d = np.sqrt((xg[i,j] - xg)**2 + (yg[i,j] - yg)**2)

                # TODO need to tune this to allow switching between different techniques or using them in parallel.

                # # Vortfind
                # # Find the nearest point in each sector.
                # ijNearest = np.zeros((len(populatedSectors), 2), dtype=np.int)
                # for k, iSec in enumerate(populatedSectors):
                #     ijNearest[k,:] = np.unravel_index(np.argmin(d + np.abs(beta-iSec)*1e12), d.shape)
                #
                # # Compute the RMS distance.
                # lnValues[i,j] = np.sum((xg[i,j] - xg[ijNearest[:,0], ijNearest[:,1]])**2
                #    + (yg[i,j] - yg[ijNearest[:,0], ijNearest[:,1]])**2) / len(populatedSectors)

                # Vortfind v2.
                # Compute radii to the current point and radial velocity direction around it.
                r = np.dstack([xg - xg[i,j], yg - yg[i,j]])
                vel = np.dstack([ug, vg])
                omegaDir = np.cross(r, vel, axis=2) / np.clip(np.abs(np.cross(r, vel, axis=2)), 1e-12, 1e12)
                omegaDir = np.nan_to_num(omegaDir)#, nan=0)

                # Find the nearest point in each sector - anticlockwise.
                # ijNearest = np.zeros((len(populatedSectors), 2), dtype=np.int)
                # for k, iSec in enumerate(populatedSectors):
                #     ijNearest[k,:] = np.unravel_index(np.argmin(d + np.abs(beta-iSec)*1e12 + np.clip(omegaDir, -1, -0.1)*-1e12), d.shape)
                # # Compute the RMS distance.
                # lnValuesV2_acw[i,j] = np.sum((xg[i,j] - xg[ijNearest[:,0], ijNearest[:,1]])**2
                #   + (yg[i,j] - yg[ijNearest[:,0], ijNearest[:,1]])**2) / len(populatedSectors)

                # Repeat for clockwise rotation.
                ijNearest = np.zeros((len(populatedSectors), 2), dtype=np.int)
                for k, iSec in enumerate(populatedSectors):
                    ijNearest[k,:] = np.unravel_index(np.argmin(d + np.abs(beta-iSec)*1e12 + np.clip(omegaDir, 0.1, 1)*1e12), d.shape)
                lnValuesV2_cw[i,j] = np.sum((xg[i,j] - xg[ijNearest[:,0], ijNearest[:,1]])**2
                    + (yg[i,j] - yg[ijNearest[:,0], ijNearest[:,1]])**2) / len(populatedSectors)

        # ijCentre_vortfind = np.unravel_index(np.argmin(lnValues), xg.shape)
        ijCentre_vortfindV2_cw = np.unravel_index(np.argmin(lnValuesV2_cw), xg.shape)
        #ijCentre_vortfindV2_acw = np.unravel_index(np.argmin(lnValuesV2_acw), xg.shape)

        # B.S. Determine the vortex centre based on minCp
        #ijVortexCentre = ijCentre_vortfindV2_cw
        ijVortexCentre = np.unravel_index(np.argmin(cpg), xg.shape)

    # When considering cavitation, override vortex identification.
    elif includeCavitation:
        if np.sum(vofg>vofThreshold) > 0:
            xm = np.mean(xg[vofg>vofThreshold])
            ym = np.mean(yg[vofg>vofThreshold])
            d = np.sqrt((xm - xg)**2 + (ym - yg)**2)
            ijVortexCentre = np.unravel_index(np.argmin(d), d.shape)

#    # Compute omega for the identified vortex centre..
#    rg = np.dstack([xg - xg[ijCentre_vortfindV2_cw[0],ijCentre_vortfindV2_cw[1]], yg - yg[ijCentre_vortfindV2_cw[0],ijCentre_vortfindV2_cw[1]]])
#    velg = np.dstack([ug, vg])
#    omegag = np.cross(rg, velg, axis=2)

    # Compute tangential velocity - not the same as omega, of course.
    r = np.vstack([x - xg[ijVortexCentre[0],ijVortexCentre[1]],
                y - yg[ijVortexCentre[0],ijVortexCentre[1]],
                np.zeros(len(y))]).T
    tangentialDir = -1.0 * np.cross(r, np.array([0, 0, 1]))
    tangentialDir /= np.linalg.norm(tangentialDir, axis=1)[:,np.newaxis]
    vel = np.vstack([u, v, np.zeros(len(y))]).T
    tanVel = np.sum(vel*tangentialDir, axis=1)#np.tensordot(vel, tangentialDir, axes=0)
#    omega = np.cross(r, vel, axis=1)
#    alpha = np.angle(x + 1j*y, deg=True)
#    alpha[alpha<0] = alpha[alpha<0] + 360.

#    alphag = np.angle(xg + 1j*yg, deg=True)
#    alphag[alphag<0] = alphag[alphag<0] + 360.

    # Retrieve the coordinate system.
    iHat = sliceProperties.loc[iSlice, ["ei_x", "ei_y", "ei_z"]].values
    jHat = sliceProperties.loc[iSlice, ["ej_x", "ej_y", "ej_z"]].values
    kHat = sliceProperties.loc[iSlice, ["ek_x", "ek_y", "ek_z"]].values
    centre = sliceProperties.loc[iSlice, ["x", "y", "z"]].values

    # Interpolate flow velocity at different angles.
    radii = np.linspace(0, approxMaxCoreRadius, nRadii)
    tanVelTransects = np.zeros((len(sectors)-1, len(radii))) * np.nan
    vofTransects = np.zeros((len(sectors)-1, len(radii))) * np.nan
    ipt = 0
    with open(os.path.join(outDataDir, "sliceTransects_{:d}.obj".format(iSlice)), "w") as outfile:
        for i in range(len(sectors)-1):
            # Create points on each transect for this slice.
            xi = radii * np.cos(sectors[i]/180*np.pi) + xg[ijVortexCentre[0], ijVortexCentre[1]]
            yi = radii * np.sin(sectors[i]/180*np.pi) + yg[ijVortexCentre[0], ijVortexCentre[1]]

            # Compute line location in 3D.
            xyz0 = jHat * xi[0] + kHat * yi[0] + centre
            xyz1 = jHat * xi[-1] + kHat * yi[-1] + centre

            # Check if the transect cuts through the geometry. If so, do not use it.
            nint = checkIfInside(xyz0, propStl.vertices, propStl.faceVertexIndices, propStl.faceNormals, propStl.faceAreas,
                                  pOut=xyz1, returnNoIntersections=True)
            if nint > 0:
                continue

            # Save the transect geometry.
            outfile.write("v {:.6e} {:.6e} {:.6e}\n".format(xyz0[0], xyz0[1], xyz0[2]))
            outfile.write("v {:.6e} {:.6e} {:.6e}\n".format(xyz1[0], xyz1[1], xyz1[2]))
            outfile.write("l {:d} {:d}\n".format(ipt+1, ipt+2))
            ipt += 2

    #        omegaTransects[i,:] = scipy.interpolate.griddata(
    #                (np.nan_to_num(alpha.flatten(), nan=0), np.linalg.norm(r, axis=1).flatten()),
    #                omega.flatten(), np.vstack((np.ones(len(radii))*sectors[i], radii)).T, method='nearest')
            tanVelTransects[i,:] = np.nan_to_num(scipy.interpolate.griddata((x, y), tanVel, np.vstack((xi, yi)).T, method=interpMethod))
            if includeCavitation:
                vofTransects[i,:] = np.nan_to_num(scipy.interpolate.griddata((x, y), vof, np.vstack((xi, yi)).T, method=interpMethod))

    # Compute mean velocity profile and vortex properties.
    tanVelMean = np.mean(tanVelTransects[~np.isnan(tanVelTransects[:, 0]), :], axis=0)
    vofTransectsMean = np.mean(vofTransects[~np.isnan(tanVelTransects[:, 0]), :], axis=0)

    #vortexData["maxTangentialVelocity"].append(np.max(tanVelMean))

    # Core radius is where tangential velocity is at its peak.
    #vortexData["coreRadius"].append(radii[np.argmax(tanVelMean)])

    # B.S. Check whether the core radius corresponds to a radius equal to zero
    orderdTanVelMean = np.argsort(tanVelMean)[::-1]

    for index in orderdTanVelMean:
      if not radii[index] < 0.005:
         vortexData["coreRadius"].append(radii[index])
         vortexData["maxTangentialVelocity"].append(tanVelMean[index])
         break

    # Compute circulation from tangential velocity - Scully model.
    vortexData["circulation"].append(4*np.pi*vortexData["coreRadius"][-1] * vortexData["maxTangentialVelocity"][-1])
    if includeCavitation:
        try:
            vortexData["cavityRadius"].append(np.max(radii[vofTransectsMean > vofThreshold]))
        except ValueError:
            # If cavity no longer exists at this slice, set radius to nan
            vortexData["cavityRadius"].append(np.nan)
    else:
        vortexData["cavityRadius"] = 0.

    # Compute vortex position in 3D space.
    xStar = xg[ijVortexCentre[0], ijVortexCentre[1]]
    yStar = yg[ijVortexCentre[0], ijVortexCentre[1]]
    xyz = centre + xStar*jHat + yStar*kHat
    for i, c in enumerate(["X", "Y", "Z"]):
        vortexData["Coordinate"+c].append(xyz[i])

    # Plot figures depicting the solution found for each slice.
    if plotOverviewFigs:
        # ---
        # Plot the quiver plot and vortex identification criteria.
        fig, ax = plt.subplots(1)
        ax.set_xlabel("x* [m]")
        ax.set_ylabel("y* [m]")
        ax.set_aspect("equal", "box")

        for i in range(len(sectors)-1):
            xi = radii * np.cos(sectors[i]/180*np.pi) + xg[ijVortexCentre[0], ijVortexCentre[1]]
            yi = radii * np.sin(sectors[i]/180*np.pi) + yg[ijVortexCentre[0], ijVortexCentre[1]]
            ax.plot(xi, yi, "k-", lw=4, alpha=0.2, zorder=-10)

        sct = ax.scatter(x, y, c="g", s=25, linewidths=0)

        ax.quiver(xg, yg, ug, vg)

        markerSize = 20
        if includeCavitation:
            sct = ax.scatter(xg, yg, c=np.nan_to_num(vofg), vmin=0, vmax=1, s=markerSize, linewidths=0, zorder=-10)
        else:
            sct = ax.scatter(xg, yg, c=np.nan_to_num(cpg), vmin=-1, vmax=1, s=markerSize, linewidths=0, zorder=-10)
        plt.colorbar(sct)

        ax.plot(x[iMinCp], y[iMinCp], "o", ms=9, mfc="c", mew=0, mec="c", label="$Min. \; C_P$")
        ax.plot(x[iMinVort], y[iMinVort], "s", ms=9, mfc="m", mew=0, mec="m", label="$Min. \; \omega_Z*$")

        ax.plot(xg[ijVortexCentre[0], ijVortexCentre[1]], yg[ijVortexCentre[0], ijVortexCentre[1]],
             ">", ms=9, mfc="b", mew=0, mec="b", label="$Vortfind-CW$")

        ax.legend(loc="center left", prop={"size":18}, ncol=3, bbox_to_anchor=(0.01, 1.0))

        if saveFigs:
            plt.savefig(os.path.join(outDataDir, "sliceData_vortexIdentification_{:d}.png".format(iSlice)), dpi=200, bbox_inches="tight")

        # ---
        # Plot vortex cuts.
        fig, ax = plt.subplots(1)
        ax.set_xlabel("r [m]")
        ax.set_ylabel("U_theta [m/s]")
        for i in range(len(sectors)-1):
#            ax.plot(radii, tanVelTransects[i,:], "-", lw=2, c="r",
            ax.plot(radii, tanVelTransects[i,:], "-", lw=2, label="$\\theta^*={:.0f}^\circ$".format(sectors[i]), zorder=-10, alpha=0.7)
        ax.plot(radii, tanVelMean, "k-", lw=4, alpha=0.5, label="$Mean$", zorder=10)
        vTanModel = vortexData["circulation"][-1] / (2*np.pi) * (radii / (vortexData["coreRadius"][-1]**2 + radii**2))
        ax.plot(radii, vTanModel, "k--", lw=4, alpha=0.5, label="$Scully \; model$", zorder=15)
        ax.legend(loc="lower right", prop={"size":18}, ncol=3, bbox_to_anchor=(1.0, 1.0))
        if saveFigs:
            plt.savefig(os.path.join(outDataDir, "sliceData_vortexVelocities_{:d}.png".format(iSlice)), dpi=200, bbox_inches="tight")

# %% Save summary data.
print(vortexData)
vortexData = pandas.DataFrame(vortexData)
vortexData.to_csv(os.path.join(outDataDir, "vortexData.csv"), index=False)

with open(os.path.join("vortexResult.dat"), "w") as f:
    #f.write("{} {}".format(vortexData["circulation"].max(), 'circulation'))
    f.write("{} {}".format(vortexData["circulation"].iloc[-1], 'circulation'))

with open(os.path.join(outDataDir, "vortexTrajectory.obj"), "w") as f:
    for i in vortexData.index:
        f.write("v {:.6e} {:.6e} {:.6e}\n".format(vortexData.loc[i,"CoordinateX"], vortexData.loc[i,"CoordinateY"], vortexData.loc[i,"CoordinateZ"]))
    for i in range(1, len(vortexData)):
        f.write("l {:d} {:d}\n".format(i, i+1))

