#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 03:55:27 2020

@author: scott
"""


from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt

from EC_MS.utils import Extraction, solve_carbonic_burst
from EC_MS.converters import mdict_from_SI2020_calibration_file
from EC_MS import standard_colors

forpublication = False
if forpublication:
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')
    plt.rc('font', size=8)
    plt.rc('lines', linewidth=0.5)

    fontsize = 8
    figwidth = 3.25 # half of a two-column figure
    figheight = 2.75


plt.close('all')

data_dir = Path("../../pickles").absolute()
data_file = data_dir / "20A31_18O_01.pkl"

calibration_dir = Path("../..").absolute()
calibration_file = calibration_dir / "20A25_sniffer_fixed.json"

mdict = mdict_from_SI2020_calibration_file(calibration_file)

experiment = Extraction(
    data_file=data_file,
    mdict=mdict,
    electrolyte="18O",
)
V_str, J_str = experiment.calibrate_EC(RE_vs_RHE=0.715, A_el=0.196)
experiment.plot_experiment()

#get water isotope ratio from O2 signals

t_bg = [14650, 14750]
tspan_OER = [15400, 15500]
alpha = experiment.get_alpha(tspan=tspan_OER, t_bg=t_bg)

burst = experiment.cut(tspan=[10200, 10750], t_zero=10278)
if True: # I don't like this, so I haven't bound it to dataset yet.
    from EC_MS import correct_shunt
    correct_shunt(burst.data, V_DL=[0.4, 0.6])

tspan_CO2 = [-10, 150]
#t_bg = [-25, -15]
t_bg = [400, 450]

plotit = True
if plotit:
    ax2 = burst.plot_experiment(mols=[], tspan=[-2, 100],
                          t_bg=t_bg,
                          logplot=False,)
    #plt.savefig('COox_faststrip_zoom_Pt.svg')



x = np.arange(1, 100, 0.25)
ys = []
for i, mass in [(0,'M44'), (1,'M46'),(2,'M48')]:
    molecule = experiment.mdict['CO2_' + mass]

    x0, y0 = burst.get_flux(molecule, tspan=tspan_CO2,
                            unit='pmol/s', t_bg=t_bg)
    # f = interp1d(x0, y0, kind='linear')
    # y = f(x)
    y = np.interp(x, x0, y0)
    ys += [y]

ys = np.array(ys)
ysum = np.sum(ys, axis=0)
yhats = ys / ysum

#fig3, ax3 = plt.subplots()
ax3 = ax2[0]
ax3.plot(x, yhats[0], color=standard_colors["M44"])
ax3.plot(x, yhats[1], color=standard_colors["M46"])
ax3.plot(x, yhats[2], color=standard_colors["M48"])
#ax3.set_yscale('log')
ax3.set_ylim([-0.03, 1.03])
ax3.tick_params(axis='y', left='on', right='on')
ax3.set_ylabel('partial signal')
ax3.set_xlabel('time / [s]')
#plt.savefig('COox_partial_signals_Pt.png')

#-------------- solve the model ------------
x_model = np.arange(0, 95, 0.5)

literaturek = True
if literaturek:
    SS = solve_carbonic_burst(k=0.037, alpha=alpha, tspan=x_model)
    ax3.plot(x_model, SS[:,0], '--', color=standard_colors["M44"])
    ax3.plot(x_model, SS[:,1], '--', color=standard_colors["M46"])
    ax3.plot(x_model, SS[:,2], '--', color=standard_colors["M48"])
    ax3.plot(0, alpha, '.', markersize=5, color=standard_colors["M44"])
    ax3.plot(0, 1-alpha, '.', markersize=5, color=standard_colors["M46"])
    ax3.plot(0, 0, '.', markersize=5, color=standard_colors["M48"])
    #plt.savefig('COox_literature_k_Pt.svg')


for ax in ax2:

    ax.set_xlim([-5, 100])

ax2[1].set_xlabel('time / [s]')
ax2[2].set_ylim([-0.1, 0.2])
ax2[2].set_ylabel(J_str)

if forpublication:
    ax3.get_figure().set_figwidth(figwidth*1.1)
    ax3.get_figure().set_figheight(figheight)
    #ax3.get_figure().tight_layout()

    plt.savefig('Pt_25C_CO_strip.png')
    plt.savefig('Pt_25C_CO_strip.svg')