# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 21:40:01 2022

@author: ALidtke
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas
import os
import re
import platform

font = {"family": "serif",
        "weight": "normal",
        "size": 16}

matplotlib.rc("font", **font)
matplotlib.rcParams["figure.figsize"] = (9, 6)

from everydaytools import nicePlots


def writeCoords(xyz, outfile, outDir):
    objStr = ""
    for i in range(xyz.shape[0]):
        objStr += "v {:.6f} {:.6f} {:.6f}\n".format(xyz[i, 0], xyz[i, 1], xyz[i, 2])

    with open(os.path.join(outDir, "{}.obj".format(outfile)), "w") as f:
        f.write(objStr)

    with open(os.path.join(outDir, "{}.dat".format(outfile)), "w") as f:
        for i in range(6):
            f.write("# comment\n")
        f.write("1 1 {:d}\n".format(xyz.shape[0]))
        f.write("\n")
        for i in range(xyz.shape[0]):
            f.write("{:.6e} {:.6e} {:.6e}\n".format(xyz[i, 0], xyz[i, 1], xyz[i, 2]))
            f.write("\n")
        f.write("\n")

    return objStr


def writePlane(name, inputFile, fields=["Velocity", "Pressure"], includeBoundary=False, interpolation="NEAREST_CELL_GRADIENT"):
    s = (
        '<monitor name="{}">'.format(name),
        '<MO_Plane>',
        '  <saveFrequency>1</saveFrequency>',
        '  <includeBoundaryData>{}</includeBoundaryData>'.format(str(includeBoundary).lower()),
        '  <interpolation><{}/></interpolation>'.format(interpolation),
        '  <fields>{}</fields>'.format(" ".join(fields)),
        '  <fileName>./data/{}</fileName>'.format(name),
        '  <fromFile>true</fromFile>',
        '  <inputFileName>{}.dat</inputFileName>'.format(inputFile),
        '  <forTecplot>true</forTecplot>',
        '</MO_Plane>',
        '</monitor>',
    )
    return "".join([l+"\n" for l in s])


dataDir = "testData"
dLine = 0.002
nProbes = 101
tMargin = 0.1
outDir = "sampling"

coords = pandas.read_csv(os.path.join(dataDir, "lineCoordinateSystems.csv"))

blProbes = []
allPoints = ""
controls = ""

for i in coords.index:
    if coords.loc[i, "isUpper"]:
        loc = "upper"
    else:
        loc = "lower"
    probes = "blProbe_{}_r_{:.4f}_xc_{:.4f}".format(loc, coords.loc[i, "r"], coords.loc[i, "x/c"])

    tBlade = coords.loc[i, "tBlade"]
    x0 = coords.loc[i, ["x", "y", "z"]].values
    n = coords.loc[i, ["nx", "ny", "nz"]].values
    xyz = x0 + n*np.linspace(-tMargin*tBlade, dLine, nProbes)[:, np.newaxis]

    allPoints += writeCoords(xyz, probes, outDir)
    controls += writePlane(probes, "./inputData/"+probes,
                            fields=["Velocity", "Pressure"],
                            interpolation="BARYCENTRIC") + "\n"
    blProbes.append(pandas.DataFrame(data=xyz, columns=["x", "y", "z"]))

# Save for visualisation and setting of controls.
with open(os.path.join(outDir, "allProbes.obj"), "w") as f:
    f.write(allPoints)

with open(os.path.join(outDir, "monitorControls.xml"), "w") as f:
    f.write(controls)
