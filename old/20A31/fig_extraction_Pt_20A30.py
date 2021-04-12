#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 22:14:09 2020

@author: scott
"""
import numpy as np
from matplotlib import pyplot as plt
plt.close('all')
if True:
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')
    plt.rc('text', fontsize=8)
    plt.rc('lines', linewidth=0.5)
fontsize = 8
figwidth = 3.25 # half of a two-column figure
figheight = 2.75


from EC_MS import Dataset, CyclicVoltammagram, load_calibration_results


mdict = load_calibration_results('20A25_calibration_results.pkl') # 20A25 calibration is more trustworthy than 20A31!
CO2_M44, CO2_M46, CO2_M48 = mdict['CO2_M44'], mdict['CO2_M46'], mdict['CO2_M48']
O2_M32, O2_M34, O2_M36 = mdict['O2_M32'], mdict['O2_M34'], mdict['O2_M36']
H2, CO, He = mdict['H2'], mdict['CO'], mdict['He']
H2.F_cal = 2 # what it god damn should be. More accurate than that measured in calibration.py because of the tilt.
H2.cal_mat = {'M2':1/H2.F_cal, 'M4':-0.0007/H2.F_cal} # gets rid of the background due to He double-ionization

r = 0.0040

dataset = Dataset('./pickles/20A30_16O_05.pkl')


dataset.plot_experiment()

if True: # remove a bit of ugly noise right in the beginning
    from EC_MS import get_timecol
    tspan_shock = [295, 315]
    t = dataset.t
    mask = np.logical_or(t<tspan_shock[0], tspan_shock[1]<t)
    for col in dataset.data_cols:
        if get_timecol(col) == 'time/s':
            dataset.data[col] = dataset.data[col][mask]

dataset.set_background(t_bg=[1080, 1100])

if True: # correct calibration factor for M46 so that it follows the natural ratio when it should
    tspan_CO2_steady = [2950, 3000]
    x_44, y_44 = dataset.get_signal(mass='M44', tspan=tspan_CO2_steady)
    x_46, y_46 = dataset.get_signal(mass='M46', tspan=tspan_CO2_steady)
    r_obs = np.mean(y_46) / np.mean(y_44)
    CO2_M46.F_cal = CO2_M46.F_cal * r_obs/r

ax = dataset.plot_experiment(
        dataset,
        tspan=[100, 2800],
        mols=[[O2_M34, O2_M36, CO2_M46, CO2_M48], [O2_M32, CO2_M44]],
        logplot=False)
ax[-1].set_ylim([l/r for l in ax[0].get_ylim()])