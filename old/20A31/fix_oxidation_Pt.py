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

if True:
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')
    plt.rc('font', size=8)
    plt.rc('lines', linewidth=0.5)

fontsize = 8
figwidth = 3.25  # half of a two-column figure
figheight = 2.75

from EC_MS import Dataset, CyclicVoltammagram, load_calibration_results
from EC_MS import align_zero

plt.close('all')

mdict = load_calibration_results('20A25_calibration_results.pkl') # 20A25 calibration is more trustworthy than 20A31!
CO2_M44, CO2_M46, CO2_M48 = mdict['CO2_M44'], mdict['CO2_M46'], mdict['CO2_M48']
O2_M32, O2_M34, O2_M36 = mdict['O2_M32'], mdict['O2_M34'], mdict['O2_M36']
H2, CO, He = mdict['H2'], mdict['CO'], mdict['He']
H2.F_cal = 2 # what it god damn should be. More accurate than that measured in calibration.py because of the tilt.
H2.cal_mat = {'M2':1/H2.F_cal, 'M4':-0.0007/H2.F_cal} # gets rid of the background due to He double-ionization

dataset = Dataset('../pickles/20A31_18O_01.pkl')

V_str, J_str = dataset.sync_metadata(RE_vs_RHE=0.715, A_el=0.196#*1e-3, J_str='J / [$\mu$A cm$^{-2}$]'
                             )


dataset.plot_experiment(mols=mdict)

dataset_ox = dataset.cut(tspan=[14300, 17100], t_zero='start')


if True: # the entire oxidation experiment thing
    dataset_ox.set_background(t_bg=[430, 450])

    axes = dataset_ox.plot_experiment(
            mols=[[O2_M32, O2_M34, O2_M36, ], [H2, CO2_M44, CO2_M46, CO2_M48]],
            #mols=[O2_M32, O2_M34, O2_M36, H2, CO2_M44, CO2_M46, CO2_M48],
            logplot=False, t_bg=[430, 450]
            )
    #axes[1].set_ylabel(J_str, fontsize=fontsize)
    align_zero(axes[3], axes[0])
    axes[1].set_ylim([-0.2, 1.9])
    axes[2].set_ylim([-0.3, 0.4])
    axes[0].set_ylabel('cal. signal / [pmol s$^{-1}$]', fontsize=fontsize)
    axes[3].set_ylabel('cal. signal / [pmol s$^{-1}$]', fontsize=fontsize)
    axes[0].get_figure().set_figwidth(figwidth*1.1) # to compensate for other panel's right y-axis label
    axes[0].get_figure().set_figheight(figheight)
    for ax_i in axes:
        ax_i.set_xlabel("time / [s]")

    plt.savefig('fig_ox_1_Pt.png')
    plt.savefig('fig_ox_1_Pt.svg')


if True: # Cyclic voltammagram figure
    cv_ox = CyclicVoltammagram(dataset_ox)

    #cv_ox.reset()
    cv_ox.set_background(t_bg=[2380, 2400])

    red_1 = cv_ox[5]
    red_2 = cv_ox[6]

    ax = red_1.plot(mols=[CO2_M44, CO2_M46, CO2_M48, H2], logplot=False)

    red_2.plot(ax=ax, linestyle='--', mols=[CO2_M44, CO2_M46, CO2_M48, H2], logplot=False)

    dQ = cv_ox.get_difference(cycle_1=5, cycle_2=6, Vspan=[0.9, 0.58], redox=0,
                              ax=ax[1], color='c', alpha=0.5)

    ax[0].set_ylabel('cal. signal / [pmol s$^{-1}$]', fontsize=fontsize)
    ax[0].get_figure().set_figwidth(figwidth*0.9) # to compensate for other panel's right y-axis label
    ax[0].get_figure().set_figheight(figheight)
    for ax_i in ax:
        ax_i.set_xlim([0, 1.25])


    plt.savefig('fig_ox_2_Pt.png')
    plt.savefig('fig_ox_2_Pt.svg')

from EC_MS.Chem import NA, Far
import numpy as np

a0 = 3.92e-10
n_ML = 4/np.sqrt(3) / (a0 ** 2) / NA * (0.196e-4)

PtO_MLs = - dQ / (2*Far) / n_ML

tspan_ox = [1500, 2000] # well into the stable time period
x_36, y_36 = dataset_ox.get_flux(O2_M36, tspan=tspan_ox)
x_34, y_34 = dataset_ox.get_flux(O2_M34, tspan=tspan_ox)

r = np.mean(y_34) / np.mean(y_36)

x = r / (2 + r)

tspan_CO2_1 = [2425, 2475]
tspan_CO2_2 = [2475, 2550]


x_48_1, y_48_1 = dataset_ox.get_flux(CO2_M48, tspan=tspan_CO2_1, unit='pmol/s')
x_46_1, y_46_1 = dataset_ox.get_flux(CO2_M46, tspan=tspan_CO2_1, unit='pmol/s')
x_48_2, y_48_2 = dataset_ox.get_flux(CO2_M48, tspan=tspan_CO2_2, unit='pmol/s')
x_46_2, y_46_2 = dataset_ox.get_flux(CO2_M46, tspan=tspan_CO2_2, unit='pmol/s')

r_CO2_1 = np.mean(y_46_1) / np.mean(y_48_1)
r_CO2_2 = np.mean(y_46_2) / np.mean(y_48_2)

n_CO2_1 = np.trapz(y_48_1, x_48_1) + np.trapz(y_46_1, x_46_1)
n_CO2_2 = np.trapz(y_48_2, x_48_2) + np.trapz(y_46_2, x_46_2)