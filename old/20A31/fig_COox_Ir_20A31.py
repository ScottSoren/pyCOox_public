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

forpublication = True
if forpublication:  # for the publication figure
    import matplotlib as mpl

    mpl.rcParams["figure.figsize"] = (3.25, 2.75)
    plt.rc("text", usetex=True)
    plt.rc("font", family="sans-serif")
    plt.rc("font", size=8)
    plt.rc("lines", linewidth=0.5)
else:
    plt.style.use("default")

fontsize = 8
figwidth = 3.25  # half of a two-column figure
figheight = 2.75

from EC_MS import plot_vs_potential, plot_experiment, sync_metadata, cut_dataset
from EC_MS import load_calibration_results, point_calibration, correct_shunt

plt.close("all")

mdict = load_calibration_results("20A25_calibration_results.pkl")
CO2_M44, CO2_M46, CO2_M48 = mdict["CO2_M44"], mdict["CO2_M46"], mdict["CO2_M48"]
O2_M32, O2_M34, O2_M36 = mdict["O2_M32"], mdict["O2_M34"], mdict["O2_M36"]
H2, CO, He = mdict["H2"], mdict["CO"], mdict["He"]
H2.F_cal = 2  # what it god damn should be. More accurate than that measured in calibration.py because of the tilt.
H2.cal_mat = {
    "M2": 1 / H2.F_cal,
    "M4": -0.0007 / H2.F_cal,
}  # gets rid of the background due to He double-ionization

with open("../pickles/20A31_18O_02.pkl", "rb") as f:
    data = pickle.load(f)

V_str, J_str = sync_metadata(
    data, RE_vs_RHE=0.715, A_el=0.196  # *1e-3, J_str='J / [$\mu$A cm$^{-2}$]'
)


plot_experiment(data, mols=mdict)

if True:  # plot CO stripping experiment

    tspan_CO_strip = [4450, 5200]
    t_bg_abs = [4930, 4950]

    t_bg = [
        t - tspan_CO_strip[0] for t in t_bg_abs
    ]  # it should be relative to the start of the interval.

    data_strip = cut_dataset(data, tspan=tspan_CO_strip, t_zero="start")

    # correct_shunt(data_ox, V_DL=[0.4, 0.6]) # not needed because this dataset is fucking beautiful.

    axes1 = plot_experiment(
        data_strip,
        mols=[[CO, He], [H2, CO2_M44, CO2_M46, CO2_M48]],
        removebackground="right",
        t_bg=t_bg,
        logplot=False,
        spec={"linewidth": 0.5},
    )

    for ax in axes1:
        ax.set_xlim([-50, 800])
        ax.set_xlabel("time / [s]", fontsize=fontsize)
    axes1[2].set_ylabel(J_str, fontsize=fontsize)
    axes1[1].set_ylabel(V_str, fontsize=fontsize)
    axes1[0].set_ylabel("cal. signal / [pmol s$^{-1}$]", fontsize=fontsize)
    axes1[-1].set_ylabel("cal. signal / [pmol s$^{-1}$]", fontsize=fontsize)
    axes1[0].get_figure().set_figwidth(
        figwidth * 1.1
    )  # to make room for right y-axis label
    axes1[0].get_figure().set_figheight(figheight)
    # axes1[0].get_figure().tight_layout()

    plt.savefig("Ir_CO_strip.png")
    plt.savefig("Ir_CO_strip.svg")

    # doesnt = exist

    tspan_c1 = [510, 610]
    tspan_c2 = [610, 710]

    axes2 = plot_vs_potential(
        data_strip,
        mols=[CO2_M44, CO2_M46, CO2_M48, H2],
        logplot=False,
        t_bg=t_bg,
        tspan=tspan_c1,
        spec={"linewidth": 0.5},
    )
    plot_vs_potential(
        data_strip,
        mols=[CO2_M44, CO2_M46, CO2_M48, H2],
        ax=axes2,
        logplot=False,
        t_bg=t_bg,
        tspan=tspan_c2,
        spec={"linestyle": "--", "linewidth": 0.5},
    )

    # for ax in axes2:
    # ax.set_xlabel(V_str, fontsize=fontsize)
    # ax.tick_params(labelsize=fontsize)
    axes2[0].set_ylim([-1.2, 12])
    axes2[1].set_ylabel(J_str, fontsize=fontsize)
    axes2[0].set_ylabel("cal. signal / [pmol s$^{-1}$]", fontsize=fontsize)
    axes2[0].get_figure().set_figwidth(
        figwidth * 0.9
    )  # to compensate for other panel's right y-axis label
    axes2[0].get_figure().set_figheight(figheight)
    # axes2[0].get_figure().tight_layout()

    plt.savefig("Ir_CO_strip_vs_potential.png")
    plt.savefig("Ir_CO_strip_vs_potential.svg")


if True:  # plot bulk CO oxidation experiment

    tspan_CO_ox = [5625, 6350]
    t_bg_abs = [5650, 5670]

    t_bg = [
        t - tspan_CO_ox[0] for t in t_bg_abs
    ]  # it should be relative to the start of the interval.

    data_ox = cut_dataset(data, tspan=tspan_CO_ox, t_zero="start")

    # correct_shunt(data_ox, V_DL=[0.4, 0.6]) # not needed because this dataset is fucking beautiful.

    axes1 = plot_experiment(
        data_ox,
        mols=[[CO, He], [H2, CO2_M44, CO2_M46, CO2_M48, O2_M32, O2_M34, O2_M36]],
        removebackground="right",
        t_bg=t_bg,
        logplot=False,
        spec={"linewidth": 0.5},
    )

    for ax in axes1:
        ax.set_xlim([-30, 760])
        ax.set_xlabel("time / [s]", fontsize=fontsize)
    axes1[2].set_ylabel(J_str, fontsize=fontsize)
    axes1[1].set_ylabel(V_str, fontsize=fontsize)
    axes1[0].set_ylabel("cal. signal / [pmol s$^{-1}$]", fontsize=fontsize)
    axes1[-1].set_ylabel("cal. signal / [pmol s$^{-1}$]", fontsize=fontsize)
    axes1[0].get_figure().set_figwidth(
        figwidth * 1.1
    )  # to make room for right y-axis label
    axes1[0].get_figure().set_figheight(figheight)
    # axes1[0].get_figure().tight_layout()

    plt.savefig("Ir_CO_ox.png")
    plt.savefig("Ir_CO_ox.svg")

    # doesnt = exist

    tspan_c1_abs = [5710, 5860]
    tspan_c2_abs = [6015, 6160]

    tspan_c1 = [
        t - tspan_CO_ox[0] for t in tspan_c1_abs
    ]  # it should be relative to the start of the interval.
    tspan_c2 = [
        t - tspan_CO_ox[0] for t in tspan_c2_abs
    ]  # it should be relative to the start of the interval.

    axes2 = plot_vs_potential(
        data_ox,
        mols=[CO2_M44, CO2_M46, CO2_M48, H2, O2_M32, O2_M34, O2_M36],
        logplot=False,
        t_bg=t_bg,
        tspan=tspan_c1,
        spec={"linewidth": 0.5},
    )
    plot_vs_potential(
        data_ox,
        mols=[CO2_M44, CO2_M46, CO2_M48, H2, O2_M32, O2_M34, O2_M36],
        ax=axes2,
        logplot=False,
        t_bg=t_bg,
        tspan=tspan_c2,
        spec={"linestyle": "--", "linewidth": 0.5},
    )

    # for ax in axes2:
    # ax.set_xlabel(V_str, fontsize=fontsize)
    # ax.tick_params(labelsize=fontsize)
    # axes2[0].set_ylim([-1, 20])
    axes2[1].set_ylabel(J_str, fontsize=fontsize)
    axes2[0].set_ylabel("cal. signal / [pmol s$^{-1}$]", fontsize=fontsize)
    axes2[0].get_figure().set_figwidth(
        figwidth * 0.9
    )  # to compensate for other panel's right y-axis label
    axes2[0].get_figure().set_figheight(figheight)
    # axes2[0].get_figure().tight_layout()

    plt.savefig("Ir_CO_ox_vs_potential.png")
    plt.savefig("Ir_CO_ox_vs_potential.svg")
