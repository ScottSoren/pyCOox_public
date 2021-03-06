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
from EC_MS import standard_colors

forpublication = True
if forpublication:
    plt.rc("text", usetex=True)
    plt.rc("font", family="sans-serif")
    plt.rc("font", size=8)
    plt.rc("lines", linewidth=0.5)

    fontsize = 8
    figwidth = 3.25  # half of a two-column figure
    figheight = 2.75
else:
    plt.style.use("default")


plt.close("all")

data_dir = Path("../../pickles").absolute()
data_file = data_dir / "17J05_Pt_in_18O_02.pkl"

calibration_dir = Path("../..").absolute()
calibration_file = calibration_dir / "20A25_sniffer_fixed.json"
# ^ wrong but deson't matter

from spitze.quant import Calibration

calibration = Calibration.load(calibration_file)

experiment = Extraction(
    data_file=data_file,
    calibration_file=calibration_file,
    calibration=calibration,  # because I don't want to import spitze in EC_MS
    electrolyte="18O",
)
V_str, J_str = experiment.calibrate_EC(RE_vs_RHE=0.715, A_el=0.196)
experiment.plot_experiment(
    tspan=[16000, 20000],
)

t_bg = [18500, 18600]

experiment.set_background(t_bg=t_bg, masses=["M32", "M34", "M36", "M44", "M46", "M48"])
if True:  # experiment plot for SI
    ax = experiment.plot_experiment(
        mols=[
            ["CO", "He"],
            ["O2_M32", "O2_M34", "O2_M36", "CO2_M44", "CO2_M46", "CO2_M48"],
        ],
        unit_left="nmol/s",
        unit_right="pmol/s/cm^2",
        tspan=[18000, 19000],
        logplot=False,
    )
    ax[-1].set_ylim([-10, 300])
    ax[0].set_ylabel("sig. / [nmol s$^{-1}$]")
    ax[-1].set_ylabel("sig. / [pmol s$^{-1}$cm$^{-2}$]")
    if forpublication:
        fig = ax[0].get_figure()
        fig.set_figwidth(figwidth)
        fig.set_figheight(figheight)
        fig.savefig("Pt_35C_burst_experiment.png")
        fig.savefig("Pt_35C_burst_experiment.svg")
# a = b

# get water isotope ratio from O2 signals
tspan_OER = [18675, 18750]
alpha = experiment.get_alpha(tspan=tspan_OER, t_bg=t_bg)

burst = experiment.cut(tspan=[18500, 18800], t_zero=18620)
if False:  # I don't like this, so I haven't bound it to dataset yet.
    from EC_MS import correct_shunt

    correct_shunt(burst.data, V_DL=[0.4, 0.6])

tspan_CO2 = [-10, 150]
# t_bg = [-25, -15]
t_bg = [-50, -10]

plotit = True
if plotit:
    ax2 = burst.plot_experiment(
        mols=[],
        tspan=[-2, 100],
        t_bg=t_bg,
        logplot=False,
    )
    # plt.savefig('COox_faststrip_zoom_Pt.svg')

x = np.arange(3, 100, 0.25)
ys = []
for i, mass in [(0, "M44"), (1, "M46"), (2, "M48")]:
    molecule = experiment.mdict["CO2_" + mass]

    x0, y0 = burst.get_flux(molecule, tspan=tspan_CO2, unit="pmol/s", t_bg=t_bg)
    # from scipy.interpolate import interp1d
    # f = interp1d(x0, y0, kind='linear')
    # y = f(x)
    y = np.interp(x, x0, y0)  #
    ys += [y]

ys = np.array(ys)
ysum = np.sum(ys, axis=0)
yhats = ys / ysum

# fig3, ax3 = plt.subplots()
ax3 = ax2[0]
ax3.plot(x, yhats[0], color=standard_colors["M44"])
ax3.plot(x, yhats[1], color=standard_colors["M46"])
ax3.plot(x, yhats[2], color=standard_colors["M48"])
# ax3.set_yscale('log')
ax3.set_ylim([-0.03, 0.83])
ax3.tick_params(axis="y", left="on", right="on")
ax3.set_ylabel("partial signal")
ax3.set_xlabel("time / [s]")
# plt.savefig('COox_partial_signals_Pt.png')

# -------------- solve the model ------------
x_model = np.arange(1, 100, 0.5)

if True:  # 35 C K
    SS = solve_carbonic_burst(k=0.08, alpha=alpha, tspan=x_model)
    ax3.plot(x_model, SS[:, 0], "--", color=standard_colors["M44"])
    ax3.plot(x_model, SS[:, 1], "--", color=standard_colors["M46"])
    ax3.plot(x_model, SS[:, 2], "--", color=standard_colors["M48"])
    ax3.plot(0, alpha, ".", markersize=5, color=standard_colors["M44"])
    ax3.plot(0, 1 - alpha, ".", markersize=5, color=standard_colors["M46"])
    ax3.plot(0, 0, ".", markersize=5, color=standard_colors["M48"])
    # plt.savefig('COox_literature_k_Pt.svg')

if True:  # 25 C K
    SS = solve_carbonic_burst(k=0.037, alpha=alpha, tspan=x_model)
    ax3.plot(x_model, SS[:, 0], ":", color=standard_colors["M44"])
    ax3.plot(x_model, SS[:, 1], ":", color=standard_colors["M46"])
    ax3.plot(x_model, SS[:, 2], ":", color=standard_colors["M48"])

for ax in ax2:

    ax.set_xlim([-5, 70])

ax2[1].set_xlabel("time / [s]")
ax2[2].set_ylim([-0.1, 0.2])
ax2[2].set_ylabel(J_str)

if forpublication:
    ax3.get_figure().set_figwidth(figwidth * 1.1)
    ax3.get_figure().set_figheight(figheight)
    # ax3.get_figure().tight_layout()

    plt.savefig("Pt_35C_CO_strip.png")
    plt.savefig("Pt_35C_CO_strip.svg")
