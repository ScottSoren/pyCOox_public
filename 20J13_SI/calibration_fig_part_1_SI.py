# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 09:26:45 2020

@author: scott
"""

from pyOER import all_measurements

for measurement in all_measurements():
    if "calibration" in measurement.category:
        ax = measurement.plot_experiment()
        ax[2].set_title(measurement.id)

