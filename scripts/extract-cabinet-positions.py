#!/usr/bin/env python

import numpy as np

from lofarantpos import geo
import pandas as pd

stationinfo_filename = "StationInfo.dat"

# extract column names
with open(stationinfo_filename) as stationinfo_file:
    for line in stationinfo_file:
        if line.startswith("# name "):
            column_names = line.split()[1:]

df = pd.read_csv(stationinfo_filename, delim_whitespace=True, comment='#', names=column_names)

# drop lines with no numbers
df = df.dropna()

# fix types
for column in df.columns:
    if column.startswith("nr"):
        df[column] = df[column].astype(int)
    elif column in ["HBAsplit", "LBAcal", "Aartfaac"]:
        df[column] = df[column].astype(bool)

df["ETRS-X"] = df["long"]
df["ETRS-Y"] = df["long"]
df["ETRS-Z"] = df["long"]

# Assume the cabinet locations are in ETRS (not verified)
df[["ETRS-X", "ETRS-Y", "ETRS-Z"]] = geo.xyz_from_geographic(
        np.deg2rad(df["long"]),
        np.deg2rad(df["lat"]),
        df["height"]).T

# Limit digits in ETRS positions, convert to string
for col in "ETRS-X", "ETRS-Y", "ETRS-Z":
    df[col] = df[col].map(lambda x: "{0:.3f}".format(x))

df.to_csv("../share/lofarantpos/stationinfo-tmp.csv", index=False)
