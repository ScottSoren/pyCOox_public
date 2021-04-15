#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 03:55:27 2020

@author: scott
"""


import pickle
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from scipy.integrate import odeint
from scipy.interpolate import interp1d

plt.rc("text", usetex=True)
plt.rc("font", family="sans-serif")
plt.rc("font", size=8)
plt.rc("lines", linewidth=0.5)

fontsize = 8
figwidth = 3.25  # half of a two-column figure
figheight = 2.75

from EC_MS import plot_vs_potential, plot_experiment, sync_metadata, cut_dataset
from EC_MS import load_calibration_results, point_calibration, correct_shunt
from EC_MS import get_signal

plt.close("all")

mdict = load_calibration_results("20A31_calibration_results.pkl")
CO2_M44, CO2_M46, CO2_M48 = mdict["CO2_M44"], mdict["CO2_M46"], mdict["CO2_M48"]
H2, CO, He = mdict["H2"], mdict["CO"], mdict["He"]
H2.F_cal = 2  # what it god damn should be. More accurate than that measured in calibration.py because of the tilt.
H2.cal_mat = {
    "M2": 1 / H2.F_cal,
    "M4": -0.0007 / H2.F_cal,
}  # gets rid of the background due to He double-ionization


with open("./pickles/20A31_18O_01.pkl", "rb") as f:
    data = pickle.load(f)
V_str, J_str = sync_metadata(
    data, RE_vs_RHE=0.715, A_el=0.196  # *1e-3, J_str='J / [$\mu$A cm$^{-1}$]'
)
t_str = "time/s"
plot_experiment(data)

# get water isotope ratio from O2 signals
getgamma = True
if getgamma:
    # t_bg = [14700, 14750]
    # tspan_O2 = [15000, 16000]
    t_bg = [14650, 14750]
    tspan_O2 = [15400, 15500]
    Y = {}
    for mass in ["M32", "M34", "M36"]:
        x, y = mdict["O2_" + mass].get_flux(
            data, tspan=tspan_O2, unit="nmol/s", t_bg=t_bg
        )
        Y[mass] = np.mean(y)
    if False:  # complex way to do it
        y_hat = np.array([Y["M32"], Y["M34"], Y["M36"]])
        y_hat = y_hat / np.sum(y_hat)

        # g is the H2(16)O / H2(18)O ratio, called gamma elsewhere
        def sqerror(g):
            return (
                (g ** 2 - y_hat[0]) ** 2
                + (2 * g * (1 - g) - y_hat[1]) ** 2
                + ((1 - g) ** 2 - y_hat[2]) ** 2
            )

        def testr(g):
            return np.array([g ** 2, 2 * g * (1 - g), (1 - g) ** 2])

        res = minimize(sqerror, 0.5)
        g = res.x[0]
    else:  # simple way to do it
        r = Y["M34"] / Y["M36"]
        g = r / (2 + r)


# a = b


def carbonic_ode(S, t, pars):
    """
    Equations from Scott's Thesis page 53
    """
    k = pars[0]  # rate constant / s^-1
    g = pars[1]  # H2(16)O / H2(18)O ratio

    S44, S46, S48 = S[0], S[1], S[2]

    dS44 = k * (-2 / 3 * (1 - g) * S44 + 1 / 3 * g * S46)
    dS46 = k * (2 / 3 * (1 - g) * S44 - 1 / 3 * S46 + 2 / 3 * g * S48)
    dS48 = k * (1 / 3 * (1 - g) * S46 - 2 / 3 * g * S48)

    return [dS44, dS46, dS48]


def solve_carbonic_burst(k=0.026, g=0.27, tspan=[0, 60]):
    """
    Returns the partial concentrations at M44, M46, and M48 following an
    initial burst of CO(16) oxidation given:
        g = the H2(18)O/H2(16)O ratio
        k = the rate constant for H2O + CO2 --> H2CO3 in s^-1
    """
    print("k = " + str(k))
    print("g = " + str(g))
    S0 = np.array([g, 1 - g, 0])
    pars = [k, g]
    if len(tspan) == 2:
        tspan = np.linspace(tspan[0], tspan[-1], 200)
    SS = odeint(carbonic_ode, S0, tspan, args=(pars,))
    return SS


data_burst = cut_dataset(data, tspan=[10200, 10750], t_zero=10278)
correct_shunt(data_burst, V_DL=[0.4, 0.6])
tspan_CO2 = [-10, 150]
# t_bg = [-25, -15]
t_bg = [400, 450]

plotit = True
if plotit:
    ax2 = plot_experiment(
        data_burst,
        mols=[],
        tspan=[-2, 100],
        t_bg=t_bg,
        logplot=False,
    )
    # plt.savefig('COox_faststrip_zoom_Pt.svg')


x = np.arange(1, 100, 0.25)
ys = []
for i, mass in [(0, "M44"), (1, "M46"), (2, "M48")]:
    x0, y0 = mdict["CO2_" + mass].get_flux(
        data_burst, tspan=tspan_CO2, unit="pmol/s", t_bg=t_bg
    )
    f = interp1d(x0, y0, kind="linear")
    ys += [f(x)]
ys = np.array(ys)
ysum = np.sum(ys, axis=0)
yhats = ys / ysum

# fig3, ax3 = plt.subplots()
ax3 = ax2[0]
ax3.plot(x, yhats[0], color=CO2_M44.get_color())
ax3.plot(x, yhats[1], color=CO2_M46.get_color())
ax3.plot(x, yhats[2], color=CO2_M48.get_color())
# ax3.set_yscale('log')
ax3.set_ylim([-0.03, 1.03])
ax3.tick_params(axis="y", left="on", right="on")
ax3.set_ylabel("partial signal")
ax3.set_xlabel("time / [s]")
# plt.savefig('COox_partial_signals_Pt.png')

# -------------- solve the model ------------
x_model = np.arange(0, 95, 0.5)

literaturek = True
if literaturek:
    SS = solve_carbonic_burst(k=0.037, g=g, tspan=x_model)
    ax3.plot(x_model, SS[:, 0], "--", color=CO2_M44.get_color())
    ax3.plot(x_model, SS[:, 1], "--", color=CO2_M46.get_color())
    ax3.plot(x_model, SS[:, 2], "--", color=CO2_M48.get_color())
    ax3.plot(0, g, ".", markersize=5, color=CO2_M44.get_color())
    ax3.plot(0, 1 - g, ".", markersize=5, color=CO2_M46.get_color())
    ax3.plot(0, 0, ".", markersize=5, color=CO2_M48.get_color())
    # plt.savefig('COox_literature_k_Pt.svg')


for ax in ax2:

    ax.set_xlim([-5, 100])

ax2[1].set_xlabel("time / [s]")
ax2[2].set_ylim([-0.1, 0.2])
ax2[2].set_ylabel(J_str)

ax3.get_figure().set_figwidth(figwidth * 1.1)
ax3.get_figure().set_figheight(figheight)
# ax3.get_figure().tight_layout()

plt.savefig("hydration_model_Pt.png")
plt.savefig("hydration_model_Pt.svg")
