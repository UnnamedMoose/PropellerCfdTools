# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 10:42:50 2022

@author: ALidtke
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas
import os
import re
import platform
from mpl_toolkits import mplot3d

font = {"family": "serif",
        "weight": "normal",
        "size": 16}

matplotlib.rc("font", **font)
matplotlib.rcParams["figure.figsize"] = (9, 6)

from everydaytools import nicePlots
from pymatt import dataReaders

# Where paraview-sampled data is
dataDir = "testData"

# Where the monitors are and where the processed data is written to.
outDir = "./case_blAnalysis_J_0.600_rpm_411_f_0.500"
case = "./case_blAnalysis_J_0.600_rpm_411_f_0.500/sampling_bary"
rpm = 411.
# outDir = "./case_blAnalysis_J_0.600_rpm_700_f_0.500"
# case = "./case_blAnalysis_J_0.600_rpm_700_f_0.500"
# rpm = 700.

coords = pandas.read_csv(os.path.join(dataDir, "lineCoordinateSystems.csv"))

# Only used for plotting.
filename_slice = os.path.join(dataDir, "sliceData_surface_r_0.7000.csv")

Uinf = 1.3974
Dprop = 0.34
omega = 2*np.pi*rpm/60

# lines = [
#     # "testData_sampledAtFaces/lineUpper_r_0.7000_xc_0.1000.csv",
#     # "testData_sampledAtCentres/lineUpper_r_0.7000_xc_0.1000.csv",
#     "testData/lineUpper_r_0.7000_xc_0.0300.csv",
#     # "testData_sampledAtFaces/lineLower_r_0.7000_xc_0.1000.csv",
#     # "testData_sampledAtCentres/lineLower_r_0.7000_xc_0.1000.csv",
#     "testData/lineLower_r_0.7000_xc_0.0300.csv",
# ]

# Read the lines extracted with Paraview
data = {}
lines = []
# lines += [f for f in os.listdir(dataDir) if f.startswith("lineUpper") or f.startswith("lineLower")]
# for iLine, line in enumerate(lines):
#     data[line] = pandas.read_csv(os.path.join(dataDir, line))

# Read ReFRESCO monitors - interpolation here is better than in paraview
monitorFiles = [f for f in os.listdir(os.path.join(case, "data")) if f.startswith("blProbe")]
for filename in monitorFiles:
    refrescoProbes = filename.replace(".dat", "")
    # "refresco_lineUpper_r_0.7000_xc_0.0300"
    data[refrescoProbes], fields, monitorData = dataReaders.read_refresco_probes(
        os.path.join(os.path.join(case, "data", filename)),
        "probes",
        experimentalFast=True)
    # Convert to Paraview-like format to reuse analysis code.
    lines.append(refrescoProbes)
    data[refrescoProbes].columns = ['Points:0', "Points:1", "Points:2"]
    data[refrescoProbes]["vtkValidPointMask"] = 1
    for i, field in enumerate(fields):
        f = field.replace("X", ":0").replace("Y", ":1").replace("Z", ":2")
        data[refrescoProbes][f] = monitorData[-1, :, i]

# lineLabels = ["upper", "lower", "upper-refresco"]

# Read the slice
surfData = pandas.read_csv(filename_slice)

# Process data.
coords["Vref"] = 0.
for iLine, line in enumerate(lines):

    if "upper" in line:
        isUpper = True
    else:
        isUpper = False

    r = float(re.findall("r_[0-9]{1}\.[0-9]{4}", line)[0].split("_")[1])
    xc = float(re.findall("xc_[0-9]{1}\.[0-9]{4}", line)[0].split("_")[1])

    Vref = np.sqrt((r*(Dprop/2.)*omega)**2. + Uinf**2.)

    coordSystem = coords[(np.abs(coords["r"]-r) < 1e-3) & (np.abs(coords["x/c"]-xc) < 1e-3) & (coords["isUpper"] == isUpper)]

    coords.loc[coordSystem.index.values[0], "Vref"] = Vref

    # Read the line monitor data.
    # Discard invalid points.
    data[line] = data[line][data[line]["vtkValidPointMask"] > 0.5].reset_index(drop=True)
    # Compute true distance from the wall (uses analytical description of the propeller blade
    # instead of the potentially crude grid).
    x0 = coordSystem[["x", "y", "z"]].values
    normal = coordSystem[["nx", "ny", "nz"]].values
    coordFields = ["Points:0", "Points:1", "Points:2"]
    # d = np.linalg.norm(data[line][coordFields].values - x0, axis=1)
    data[line]["d"] = np.sum(normal * (data[line][coordFields].values - x0), axis=1)
    data[line] = data[line][data[line]["d"] >= 0].reset_index(drop=True)
    # Add the effect of rotation
    # XXX NOTE for a typical prop set up this should be reversed
    data[line]["Velocity:1"] += -omega*data[line]["Points:2"]
    data[line]["Velocity:2"] += omega*data[line]["Points:1"]
    # Add BC at the wall (if the sampling is dense enough, the first
    # value should be almost zero anyway).
    data[line].loc[0, ["Velocity:0", "Velocity:1", "Velocity:2"]] = 0.
    data[line].loc[0, "d"] = 0.
    # Compute velocity magnitude.
    data[line]["Umag"] = np.sqrt(data[line]["Velocity:0"]**2. +
                                 data[line]["Velocity:1"]**2. +
                                 data[line]["Velocity:2"]**2.)
    # Resolve into tangential and spanwise components
    vt = coordSystem[["tx", "ty", "tz"]].values
    vs = coordSystem[["sx", "sy", "sz"]].values
    data[line]["Ut"] = np.sum(vt*data[line][["Velocity:0", "Velocity:1", "Velocity:2"]].values, axis=1)
    data[line]["Us"] = np.sum(vs*data[line][["Velocity:0", "Velocity:1", "Velocity:2"]].values, axis=1)
    data[line]["Ut/Uref"] = data[line]["Ut"] / Vref
    data[line]["Us/Uref"] = data[line]["Us"] / Vref

    data[line].to_csv(os.path.join(outDir, line+".csv"), index=False)

coords.to_csv(os.path.join(outDir, "lineCoordinateSystems.csv"), index=False)

# Situational awareness is key
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlabel('x')
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.scatter(surfData["Points:0"], surfData["Points:1"], surfData["Points:2"], c=surfData["Pressure"], s=5, cmap="jet")
for iLine, line in enumerate(lines):
    ax.plot(data[line]["Points:0"], data[line]["Points:1"], data[line]["Points:2"], ".", label=line)
# ax.legend()

for iTerm, term in enumerate(["x", "y", "z"]):
    fig, ax = nicePlots.niceFig([], [], "U_{} [m/s]".format(term), "Wall distance [m]")
    for iLine, line in enumerate(lines):
        ax.plot(data[line]["Velocity:{:d}".format(iTerm)], data[line]["d"], "-", lw=1, ms=5, label=line)
    # ax.legend()

for iTerm, term in enumerate(["s", "t"]):
    fig, ax = nicePlots.niceFig([], [], "U_{}/U_{{ref}}".format(term), "Wall distance [m]")
    for iLine, line in enumerate(lines):
        ax.plot(data[line]["U{}/Uref".format(term)], data[line]["d"], "-", lw=1, ms=5, label=line)
    # ax.legend()
    plt.savefig(os.path.join(outDir, "velocity_{}.png".format(term)), dpi=200, bbox_inches="tight")

fig, ax = nicePlots.niceFig([], [], "|U| [m/s]", "Wall distance [m]")
for iLine, line in enumerate(lines):
    ax.plot(data[line]["Umag"], data[line]["d"], "-", lw=1, ms=5, label=line)
# ax.legend()
plt.savefig(os.path.join(outDir, "velocityMagnitude.png"), dpi=200, bbox_inches="tight")
