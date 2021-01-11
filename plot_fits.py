#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 09:26:45 2021

@author: miha
"""

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from parameters import d
from parameters import p






cStates = d["charge_states"]
filePaths = d["parsed_data_files"]

df = pd.DataFrame()
df["cState"] = cStates
df["fPath"] = filePaths


c = df["cState"] == 5
data = pd.read_csv(df[c]["fPath"].values[0])

ndf = pd.DataFrame(data)
ndf.columns = ["t", "i"]

t,i = ndf["t"],ndf["i"]

plt.plot(t,i)





