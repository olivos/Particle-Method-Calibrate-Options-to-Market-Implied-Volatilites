#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 12:27:06 2019

@author: oliv
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

paths = pd.read_csv("Paths.csv",header =None)
pnl = pd.read_csv("PNL30.csv",header =None)
pnl120 = pd.read_csv("PNL120.csv",header =None)


T = [0,0.00277778,0.00555556,0.00833333,0.0111111,0.0138889,0.0166667,0.0194444,0.0222222,0.025,0.0277778,0.0305556,0.0333333,0.0361111,0.0388889,0.0416667,0.0444444,0.0472222,0.05,0.0527778,0.0555556,0.0583333,0.0611111,0.0638889,0.0666667,0.0694444,0.0722222,0.075,0.0777778,0.0805556,0.0833333]

x = [0.3,0.35,0.4]
r30 = np.array([0.542144 , 0.575519 , 0.653081 ])
r120 = np.array([0.26941, 0.336617 , 0.458673])

#plt.hist(pnl.values[:,-1],bins = 50,label = "30 rehedgings",alpha = 0.8)
#plt.hist(pnl120.values[:,-1],bins = 50,label = "120 rehedgings",alpha = 0.8)
plt.plot(x,r30/r120)

plt.title("Reduction of std when hedging 4 times more often")
plt.legend()

