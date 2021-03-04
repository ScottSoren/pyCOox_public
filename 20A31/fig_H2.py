#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 00:08:17 2020

@author: scott


The things I'm struggling with are pretty standard. I'm drifting. Some times
drifting is the most beautiful thing in the world. And I can't see it coming.
But luck strikes the prepared, and writing these articles will release a
heavy chain from my ankle and take a heavy weight off my shoulders so that
I can soon dance in the sun and catch a drifting flower.

"""

import pickle
from matplotlib import pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rc('text', fontsize=8)
plt.rc('lines', linewidth=0.5)

fontsize = 8
figwidth = 3.25 # half of a two-column figure
figheight = 2.75

from EC_MS import Dataset, CyclicVoltammagram, load_calibration_results

plt.close('all')

mdict = load_calibration_results('20A25_calibration_results.pkl') # 20A25 calibration is more trustworthy than 20A31!
CO2_M44, CO2_M46, CO2_M48 = mdict['CO2_M44'], mdict['CO2_M46'], mdict['CO2_M48']
O2_M32, O2_M34, O2_M36 = mdict['O2_M32'], mdict['O2_M34'], mdict['O2_M36']
H2, CO, He = mdict['H2'], mdict['CO'], mdict['He']
H2.F_cal = 2 # what it god damn should be. More accurate than that measured in calibration.py because of the tilt.
H2.cal_mat = {'M2':1/H2.F_cal, 'M4':-0.0007/H2.F_cal} # gets rid of the background due to He double-ionization

dataset = Dataset('./pickles/20A31_18O_01.pkl')

V_str, J_str = dataset.sync_metadata(RE_vs_RHE=0.715, A_el=0.196#*1e-3, J_str='J / [$\mu$A cm$^{-2}$]'
                             )


dataset.plot_experiment(mols=mdict)



if False: # the entire H2 thing
    dataset_H2 = dataset.cut(tspan=[17100, 20100], t_zero='start')
    dataset_H2.plot_experiment(mols=[H2, He])

if True: # RHE cal and HOR j_lim measurements
    dataset_H2 = dataset.cut(tspan=[18200, 19500], t_zero='start')
    axes = dataset_H2.plot_experiment(mols=[H2, He], logplot=False, unit='nmol/s')

    #axes[1].set_ylabel(J_str, fontsize=fontsize)
    axes[0].set_ylabel('cal. signal / [nmol s$^{-1}$]', fontsize=fontsize)
    axes[0].get_figure().set_figwidth(figwidth*0.9) # to compensate for other panel's right y-axis label
    axes[0].get_figure().set_figheight(figheight)
    axes[0].get_figure().tight_layout()

    plt.savefig('RHE_and_HOR_lim.png')
    plt.savefig('RHE_and_HOR_lim.svg')

    cvs_H2 = CyclicVoltammagram(dataset, tspan=[18000, 20000])
    cvs_H2.redefine_cycle(V=0.4, redox=1)

    cvs_H2.plot_experiment(J_str='cycle')


    cycle_He = cvs_H2[8]
    cycle_H2 = cvs_H2[5]

    diff = cycle_H2.subtract(cycle_He)

    ax = diff.plot(logplot=False)
    #ax[0].legend()




