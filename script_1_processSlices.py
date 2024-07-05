import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas
import os
from shapely.geometry import Point, Polygon

Dprop = 0.34
rpm = 600.
rhoInf = 998.9
muInf = 1.094e-3
radii = [r*Dprop/2. for r in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]]
outDir = "D:/git/propellerCfdTools"

# Process all the slices.
sliceData = []
for iRadius, radius in enumerate(radii):
    sliceFile = 'sliceData_surface_r_{:.4f}.csv'.format(radius/(Dprop/2.))
    slc = pandas.read_csv(os.path.join(outDir, sliceFile))

    # Compute the pressure coefficient.
    # NOTE: this used (nD)^2 as the reference velocity as is customary for propellers.
    slc["Cp"] = slc["Pressure"] / (rhoInf * (Dprop * rpm/60.)**2.)

    # Project onto the surface at a constant radius.
    slc["r"] = np.sqrt(slc["Points:1"]**2. + slc["Points:2"]**2.)
    slc["theta"] = -np.arctan2(slc["Points:2"], slc["Points:1"])
    slc["xPrime"] = slc["Points:0"]
    slc["yPrime"] = slc["r"]*slc["theta"]

    # Find the LE and TE point indices.
    chord = 0
    max_pair = None

    # Iterate over all pairs of points
    for i in range(len(slc)):
        for j in range(i+1, len(slc)):
            p1 = slc.loc[i, ["xPrime", "yPrime"]]
            p2 = slc.loc[j, ["xPrime", "yPrime"]]
            # Compute the distance between points p1 and p2
            distance = np.linalg.norm(p1 - p2)
            # Check if this is the largest distance found so far
            if distance > chord:
                chord = distance
                max_pair = [i, j]

    # Locate LE and TE
    if slc.loc[max_pair[0], "xPrime"] < slc.loc[max_pair[1], "xPrime"]:
        iLe, iTe = max_pair
    else:
        iTe, iLe = max_pair

    # Unit vector along the chord.
    vx = (slc.loc[iLe, ["xPrime", "yPrime"]] - slc.loc[iTe, ["xPrime", "yPrime"]]) / chord

    # Compute x/c.
    slc["x/c"] = 0.
    for i in range(slc.shape[0]):
        slc.loc[i, "x/c"] = 1. - \
            np.dot(vx, slc.loc[i, ["xPrime", "yPrime"]] - slc.loc[iTe, ["xPrime", "yPrime"]]) / chord

    # Find a point in the centre of the section, more or less.
    # This helps to deal with separating upper and lower surfaces when a lot of camber is present.
    xyCentre = slc[np.abs(slc["x/c"] - 0.5) < 0.1][["xPrime", "yPrime"]].mean().values

    # Separate upper and lower surfaces.
    polygon_lower = Polygon([
        slc.loc[iTe, ["xPrime", "yPrime"]],
        xyCentre,
        slc.loc[iLe, ["xPrime", "yPrime"]],
        [slc.loc[iTe, "xPrime"], slc.loc[iLe, "yPrime"]],
    ])

    slc["isUpper"] = 1.
    slc["isLower"] = 0.
    slc.loc[iTe, "isLower"] = 1.
    slc.loc[iLe, "isLower"] = 1.
    for i in range(slc.shape[0]):
        pt = Point(slc.loc[i, ["xPrime", "yPrime"]].values.copy())
        if polygon_lower.contains(pt):
            slc.loc[i, "isUpper"] = 0.
            slc.loc[i, "isLower"] = 1.

    # Save.
    sliceData.append(slc)
    slc.to_csv(os.path.join(outDir, "processed_"+sliceFile), index=False)

# Plot an example slice.
slc = sliceData[2]

# ===
fig, ax = plt.subplots()
ax.set_xlabel("y'")
ax.set_ylabel("x'")
ax.invert_yaxis()
cs = plt.scatter(slc["yPrime"], slc["xPrime"], c=slc["x/c"], s=15)
cbar = plt.colorbar(cs)
cbar.set_label("x/c")

df = slc[slc["isUpper"] > 0.5].copy().sort_values("xPrime")
ax.plot(df["yPrime"], df["xPrime"], "r--", zorder=-100, label="Upper")

df = slc[slc["isLower"] > 0.5].copy().sort_values("xPrime")
ax.plot(df["yPrime"], df["xPrime"], "m--", zorder=-100, label="Lower")

ax.legend()

# ===
fig, ax = plt.subplots()
ax.set_xlabel("x/c")
ax.set_ylabel("-C$_{P, \; n}$")
ax.plot(slc["x/c"], -slc["Cp"], ".")

plt.show()
