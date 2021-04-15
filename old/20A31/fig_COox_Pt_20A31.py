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
    figwidth = 7
    figheight = 5
    fontsize = 12

from EC_MS import plot_vs_potential, plot_experiment, sync_metadata, cut_dataset
from EC_MS import load_calibration_results, point_calibration, correct_shunt
from EC_MS.utils.extraction_class import get_EC_MS_mdict

# from siCalibration import Calibration


plt.close("all")

# calibration = Calibration.load('20A25_sniffer_fixed.json')
# mdict = get_EC_MS_mdict(calibration)
mdict = load_calibration_results("20A25_calibration_results.pkl")
# 20A25 calibration is more trustworthy than 20A31!
CO2_M44, CO2_M46, CO2_M48 = mdict["CO2_M44"], mdict["CO2_M46"], mdict["CO2_M48"]
O2_M32, O2_M34, O2_M36 = mdict["O2_M32"], mdict["O2_M34"], mdict["O2_M36"]
H2, CO, He = mdict["H2"], mdict["CO"], mdict["He"]
H2.F_cal = 2  # what it god damn should be. More accurate than that measured in calibration.py because of the tilt.
H2.cal_mat = {
    "M2": 1 / H2.F_cal,
    "M4": -0.0007 / H2.F_cal,
}  # gets rid of the background due to He double-ionization

with open("../pickles/20A31_18O_01.pkl", "rb") as f:
    data = pickle.load(f)

V_str, J_str = sync_metadata(
    data,
    RE_vs_RHE=0.715,
    A_el=0.196
    # *1e-3, J_str='J / [$\mu$A cm$^{-2}$]'
)

plot_experiment(data, mols=mdict)
# a = b
if True:  # plot CO stripping experiment

    tspan_CO_strip = [9700, 10500]
    t_bg_abs = [9900, 9920]

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
        unit_left="nmol/s",
        unit_right="pmol/s/cm^2",
    )

    for ax in axes1:
        ax.set_xlim([-30, 830])
        ax.set_xlabel("time / [s]", fontsize=fontsize)
    axes1[2].set_ylabel(J_str, fontsize=fontsize)
    axes1[1].set_ylabel(V_str, fontsize=fontsize)
    axes1[0].set_ylabel("cal. sig. / [nmol s$^{-1}$]", fontsize=fontsize)
    axes1[-1].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]", fontsize=fontsize)
    axes1[0].get_figure().set_figwidth(
        figwidth * 1.1
    )  # to make room for right y-axis label
    axes1[0].get_figure().set_figheight(figheight)
    # axes1[0].get_figure().tight_layout()

    if forpublication:
        plt.savefig("CO_strip.png")
        plt.savefig("CO_strip.svg")

    # doesnt = exist

    tspan_c1 = [530, 640]
    tspan_c2 = [640, 750]

    axes2 = plot_vs_potential(
        data_strip,
        mols=[CO2_M44, CO2_M46, CO2_M48, H2],
        logplot=False,
        t_bg=t_bg,
        tspan=tspan_c1,
        left_space=0.2,  # to avoid cutting the axis label
        spec={"linewidth": 0.5},
        unit="pmol/s/cm^2",
    )
    plot_vs_potential(
        data_strip,
        mols=[CO2_M44, CO2_M46, CO2_M48, H2],
        ax=axes2,
        logplot=False,
        t_bg=t_bg,
        tspan=tspan_c2,
        spec={"linestyle": "--", "linewidth": 0.5},
        unit="pmol/s/cm^2",
    )

    # for ax in axes2:
    # ax.set_xlabel(V_str, fontsize=fontsize)
    # ax.tick_params(labelsize=fontsize)
    axes2[0].set_ylim([-5, 100])
    axes2[1].set_ylabel(J_str, fontsize=fontsize)
    axes2[0].set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    axes2[1].set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    axes2[0].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]", fontsize=fontsize)
    axes2[0].get_figure().set_figwidth(
        figwidth * 0.9
    )  # to compensate for other panel's right y-axis label
    axes2[0].get_figure().set_figheight(figheight)
    # axes2[0].get_figure().tight_layout()

    if forpublication:
        plt.savefig("CO_strip_vs_potential.png")
        plt.savefig("CO_strip_vs_potential.svg")

if True:  # plot bulk CO oxidation experiment

    tspan_CO_ox = [11350, 12150]
    t_bg_abs = [11620, 11660]

    t_bg = [
        t - tspan_CO_ox[0] for t in t_bg_abs
    ]  # it should be relative to the start of the interval.

    data_ox = cut_dataset(data, tspan=tspan_CO_ox, t_zero="start")

    # correct_shunt(data_ox, V_DL=[0.4, 0.6]) # not needed because this dataset is
    # fucking beautiful.

    axes1 = plot_experiment(
        data_ox,
        mols=[[CO, He], [H2, CO2_M44, CO2_M46, CO2_M48, O2_M32, O2_M34, O2_M36]],
        removebackground="right",
        t_bg=t_bg,
        logplot=False,
        spec={"linewidth": 0.5},
        unit_left="nmol/s",
        unit_right="pmol/s/cm^2",
    )

    for ax in axes1:
        ax.set_xlim([-30, 830])
        ax.set_xlabel("time / [s]", fontsize=fontsize)
    axes1[2].set_ylabel(J_str, fontsize=fontsize)
    axes1[1].set_ylabel(V_str, fontsize=fontsize)
    axes1[0].set_ylabel("cal. sig. / [nmol s$^{-1}$]", fontsize=fontsize)
    axes1[-1].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]", fontsize=fontsize)
    axes1[0].get_figure().set_figwidth(
        figwidth * 1.1
    )  # to make room for right y-axis label
    axes1[0].get_figure().set_figheight(figheight)
    # axes1[0].get_figure().tight_layout()

    if forpublication:
        plt.savefig("CO_ox.png")
        plt.savefig("CO_ox.svg")

    # doesnt = exist

    tspan_c1_abs = [11440, 11610]
    tspan_c2_abs = [11790, 11960]

    tspan_c1 = [t - tspan_CO_ox[0] for t in tspan_c1_abs]
    # it should be relative to the start of the interval.
    tspan_c2 = [t - tspan_CO_ox[0] for t in tspan_c2_abs]
    # it should be relative to the start of the interval.

    axes2 = plot_vs_potential(
        data_ox,
        mols=[CO2_M44, CO2_M46, CO2_M48, H2, O2_M32, O2_M34, O2_M36],
        left_space=0.2,  # to avoid cutting the axis label
        logplot=False,
        t_bg=t_bg,
        tspan=tspan_c1,
        spec={"linewidth": 0.5},
        unit="pmol/s/cm^2",
    )
    plot_vs_potential(
        data_ox,
        mols=[CO2_M44, CO2_M46, CO2_M48, H2, O2_M32, O2_M34, O2_M36],
        ax=axes2,
        logplot=False,
        t_bg=t_bg,
        tspan=tspan_c2,
        spec={"linestyle": "--", "linewidth": 0.5},
        unit="pmol/s/cm^2",
    )

    # for ax in axes2:
    # ax.set_xlabel(V_str, fontsize=fontsize)
    # ax.tick_params(labelsize=fontsize)
    # axes2[0].set_ylim([-1, 20])
    axes2[1].set_ylabel(J_str, fontsize=fontsize)
    axes2[0].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]", fontsize=fontsize)
    axes2[0].get_figure().set_figwidth(figwidth * 0.9)
    axes2[0].get_figure().set_figheight(figheight)
    # axes2[0].get_figure().tight_layout()

    if forpublication:
        plt.savefig("CO_ox_vs_potential.png")
        plt.savefig("CO_ox_vs_potential.svg")


if True:  # get portion C^{16}O^{18}O lost to C^{18}O2
    import numpy as np
    from EC_MS import get_signal  # yuck, this is written the old way.

    tspan = [400, 800]
    x46, y46 = get_signal(data_ox, mass="M46", tspan=tspan)
    x48, y48 = get_signal(data_ox, mass="M48", tspan=tspan)

    Y46 = np.trapz(y46, x46)
    Y48 = np.trapz(y48, x48)

    x_lost = Y48 / (Y48 + Y46)
