#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 13:35:22 2020

@author: scott
"""

import pickle
import numpy as np
from matplotlib import pyplot as plt

plt.close("all")

from EC_MS import (
    load_calibration_results,
    Chem,
    Chip,
    Molecule,
    Dataset,
    CyclicVoltammagram,
    stagnant_operator,
    plot_operation,
    standard_colors,
)


plt.rc("text", usetex=True)
plt.rc("font", family="sans-serif")
plt.rc("font", size=8)
plt.rc("lines", linewidth=0.5)

fontsize = 8
figwidth = 3.25  # half of a two-column figure
figheight = 2.75

mdict = load_calibration_results("20A31_calibration_results.pkl")
CO2_M44, CO2_M46, CO2_M48 = mdict["CO2_M44"], mdict["CO2_M46"], mdict["CO2_M48"]
O2_M32, O2_M34, O2_M36 = mdict["O2_M32"], mdict["O2_M34"], mdict["O2_M36"]
H2, CO, He = mdict["H2"], mdict["CO"], mdict["He"]
H2.F_cal = 2  # what it god damn should be. More accurate than that measured in calibration.py because of the tilt.
H2.cal_mat = {
    "M2": 1 / H2.F_cal,
    "M4": -0.0007 / H2.F_cal,
}  # gets rid of the background due to He double-ionization

chip = Chip("DanChip-14-20A31")  # saved by calibration.py

dataset = Dataset("../pickles/20A31_18O_01.pkl")
A_el = 0.196

V_str, J_str = dataset.sync_metadata(
    RE_vs_RHE=0.715, A_el=A_el  # *1e-3, J_str='J / [$\mu$A cm$^{-2}$]'
)
sel_str = dataset.make_selector()


dataset.plot_experiment()

t_str = "time/s"

fitL = True
if fitL:
    tspan_HOR_steady = [19280, 19300]
    _, I = dataset.get_current(tspan_HOR_steady, unit="A")

    I_HOR = np.mean(I)

    # I_lim = j_lim*A_el = 2 * Chem.Far * A_el * H2.D/L * p0 / H2.kH
    T = 298.15  # 25 C in [K]
    kH = H2.kH * Chem.R * T  # in [bar/M] --> [Pa / (mol/m*3])]
    D = H2.D  # in [m^2/s]
    p0 = 1e5  # in Pa
    L = (
        p0 / kH * D * (2 * Chem.Far * 0.196e-4 / I_HOR)
    )  # [Pa / (Pa / (mol/m^3)) * m^2/s * C/mol * m^2 / (C/s)] = [m]
    # gives 130 um, woo-hoo! That seems about right.

modelCV = True
if modelCV:
    tspan_CVs = [11000, 13000]
    cv = CyclicVoltammagram(dataset, tspan=tspan_CVs)
    cv.get_sweeps()

    # cv.plot_experiment(J_str='sweep')
    # cv.plot_experiment(mols=[[He, CO], [CO2_M44, CO2_M46, CO2_M48, O2_M32, O2_M34, O2_M36]])

    # the cycles I want are close to 14 and 17, but the cycle switch is at an unfortunate place

    cv.redefine_cycle(redox=1, V=0.45)
    cv.plot_experiment(
        J_str="cycle",
        # J_str = 'scan_rate',
    )

    tspan_MS = [0, 180]
    tspan_EC = [0, 135]
    if True:  # complex version, keep MS tail
        cycle_CO = cv[5:7].cut(tspan=tspan_MS)
        cycle_He = cv[3:5].cut(tspan=tspan_MS)
    else:  # simple version, loose MS tail
        cycle_CO = cv[5]  # bloody fucking beautiful!
        cycle_He = cv[3]

    diff = cycle_CO.subtract(cycle_He)  # also bloody beautiful.

    plotit = True
    if plotit:

        # plt.savefig('COox_CVs.png')

        ax1 = cycle_CO.plot_vs_potential(
            mols=[H2, CO2_M44, CO2_M46, CO2_M48, O2_M32, O2_M34, O2_M36], logplot=False
        )
        cycle_He.plot_vs_potential(
            mols=[H2, CO2_M44, CO2_M46, CO2_M48, O2_M32, O2_M34, O2_M36],
            linestyle="--",
            ax=ax1,
        )
        # plt.savefig('COox_CVs_vs_potential.png')

    predictsignal = True
    if predictsignal:

        tspan_model = tspan_MS
        ax3 = diff.plot_experiment(
            mols=[CO2_M44, CO2_M46, CO2_M48],
            tspan=tspan_model,
            plotpotential=False,
            J_color="k",
            tspan_EC=tspan_EC,
            logplot=False,
        )
        ax3[1].set_xlim(ax3[0].get_xlim())
        plt.savefig("COox_CVs_predicted1.png")

        runsimulation = False
        if runsimulation:
            # L = 150e-6
            # chip.l_cap = 1.4e-3
            q0 = chip.capillary_flow(gas="CO") / Chem.NA
            t, j_el = diff.get_current(tspan=tspan_EC, unit="A/m^2")
            # need to fit L (above) to get a good model. Turns out we were way indented! L=160um
            results = stagnant_operator(
                tj=[t, j_el], tspan=tspan_model, L=L, mol=CO2_M44, n_el=2, q0=q0
            )

        ax3[0].plot(results["t"], results["j"] / 3, "k--")
        plt.savefig("COox_CVs_predicted_Pt.png")

        ax3[0].clear()
        for mass in ["M44", "M46", "M48"]:
            x, y = diff.get_signal(mass)
            y_norm = (y - min(y)) / (max(y) - min(y))
            ax3[0].plot(x, y_norm, standard_colors[mass])

        ax3[0].plot(results["t"], results["j"] / max(results["j"]), "k--")
        ax3[0].set_xlim(ax3[1].get_xlim())
        ax3[0].set_ylabel("normalized flux")
        plt.savefig("COox_CVs_predicted_normalized.png")


if True:  # make a 3-panel plot connecting j_diff to the predicted signal

    fig, ax = plt.subplots(nrows=3, sharex=True)

    ax[0].plot(t, j_el * 1e-1, "r-")
    plot_operation(results, ax=ax[1], colorbar=False)

    t_sim, j_sim = results["t"], results["j"]

    j_tot = np.zeros(t_sim.shape)
    ax[2].plot(t_sim, j_sim / A_el * 1e-3, "k--")
    ax[2].set_ylabel("$\dot{n}$ / [nmol s$^{-1}$cm$^{-2}$]")

    plotsignals = True  # plot individual norm. signals on right y-axis of bottom panel
    if plotsignals:
        ax = list(ax) + [ax[2].twinx()]
        ax[3].set_ylabel("norm. signal")
        ax[3].set_yticks([])  # no ticks.

    for m in [CO2_M44, CO2_M46, CO2_M48]:
        x, y = cycle_CO.get_flux(m, unit="pmol/s", removebackground=True)
        j_tot += np.interp(t_sim, x, y)
        y = y / max(y)
        if plotsignals:
            ax[3].plot(x, y, color=m.get_color())

    if plotsignals:
        ax[3].set_yticks([])  # no ticks.
    if False:  # then plot total
        ax[2].plot(t_sim, j_tot, "r")

    # ax[3].plot(t_sim, j_tot, 'k')

    ax[0].set_ylabel("$\Delta$J / [mA cm$^{-2}$]")
    ax[1].set_ylabel("y / [$\mu$m]")
    ax[2].set_xlabel("time / [s]")
    fig.subplots_adjust(hspace=0)
    for a in ax[0:2]:
        a.tick_params(axis="x", top="on", bottom="on", direction="in")
        a.tick_params(axis="y", left="on", right="on", direction="in")
    ax[0].xaxis.set_label_position("top")
    ax[0].set_xlabel("time / [s]")
    ax[0].tick_params(axis="x", labeltop="on")
    fig.set_figwidth(8)
    fig.set_figheight(6)

    ax_sig = ax[2]

    ax[0].get_figure().set_figwidth(figwidth * 1.1)
    ax[0].get_figure().set_figheight(figheight * 1.35)
    ax[0].get_figure().tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    if plotsignals:
        plt.savefig("COox_CVs_model_Pt.svg")
    else:
        plt.savefig("COox_CVs_model_total_Pt.svg")
