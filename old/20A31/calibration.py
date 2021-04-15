#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 00:08:17 2020

@author: scott
"""
import pickle
import numpy as np
from matplotlib import pyplot as plt

from EC_MS import Dataset, standard_colors, is_time

from EC_MS import plot_experiment, cut_dataset
from EC_MS import point_calibration, calibration_curve
from EC_MS import recalibrate, chip_calibration, Molecule

from EC_MS import synchronize, sort_time, plot_signal
from EC_MS import save_calibration_results


plt.close("all")

if True:  # use the Dataset interface:
    dataset = Dataset("../pickles/20A31_18O_01.pkl")
else:
    with open("../pickles/20A31_18O_01.pkl", "rb") as f:
        data_1 = pickle.load(f)

with open("../pickles/20A31_18O_MS_data.pkl", "rb") as f:
    data_MS = pickle.load(f)

# plot_signal(data_MS)

mdict = {}

t_zero = 12250

dataset.tspan = dataset.tspan + t_zero
for col in dataset.data_cols:
    if is_time(col, dataset.data) or col.endswith("-x"):
        print(f"shifting {col}")
        dataset[col] = dataset[col] - t_zero

ax_0 = dataset.plot_experiment(
    unit="nA",
    tspan=[0, 7500],
    emphasis=None,
)

tspan_CO2 = [12589 - t_zero, 12832 - t_zero]  # [2100, 2300]
t_bg_CO2 = [1350, 1550]
F_cal_CO2 = 0
for mass in ["M44", "M46", "M48"]:

    m = dataset.point_calibration(mol="CO2", mass=mass, tspan=tspan_CO2, n_el=2)
    m.primary = mass
    F_cal_CO2 += m.F_cal
    m.name = "CO2_" + mass
    mdict["CO2_" + mass] = m

    x_bg, y_bg = dataset.get_signal(mass=mass, tspan=t_bg_CO2, unit="nA")
    x, y = dataset.get_signal(mass=mass, tspan=tspan_CO2, unit="nA")
    y_bg = np.ones(x.shape) * np.mean(y_bg)
    ax_0[0].fill_between(x, y_bg, y, color=standard_colors[mass], alpha=0.5)

for mass in ["M44", "M46", "M48"]:
    mdict["CO2_" + mass].F_cal = F_cal_CO2


tspan_O2 = [16000 - t_zero, 16600 - t_zero]
t_bg_O2 = [21500 - t_zero, 21600 - t_zero]
F_cal_O2 = 0
for mass in ["M32", "M34", "M36"]:
    m = dataset.point_calibration(
        mol="O2", mass=mass, tspan=tspan_O2, n_el=4, t_bg=t_bg_O2
    )
    m.primary = mass
    F_cal_O2 += m.F_cal
    mdict["O2_" + mass] = m

    x_bg, y_bg = dataset.get_signal(mass=mass, tspan=t_bg_O2, unit="nA")
    x, y = dataset.get_signal(mass=mass, tspan=tspan_O2, unit="nA")
    y_bg = np.ones(x.shape) * np.mean(y_bg)
    ax_0[0].fill_between(x, y_bg, y, color=standard_colors[mass], alpha=0.5)

for mass in ["M32", "M34", "M36"]:
    mdict["O2_" + mass].F_cal = F_cal_O2


# doesnt = exist
H2, ax = dataset.calibration_curve(
    ax=[ax_0, "new"],
    mol="H2",
    n_el=-2,
    t_bg=[17100 - t_zero, 17150 - t_zero],
    tspans=[
        [17250 - t_zero, 17300 - t_zero],
        [17510 - t_zero, 17530 - t_zero],
        [17725 - t_zero, 17750 - t_zero],
    ],
    out=["Molecule", "ax"],
    unit="nA",
)
ax_cal_curve_H2 = ax[1][1]

mdict["H2"] = H2

tspan_He = [58700, 58800]
tspan_air = [47000, 48000]
tspan_H2 = [35200, 35400]
tspan_CO = [17800, 18000]

chip = chip_calibration(
    data_MS, mol=mdict["O2_M32"], tspan=tspan_air, tspan_bg=tspan_He
)
chip.save("DanChip-14-20A31")

He = point_calibration(
    data_MS,
    mol="He",
    cal_type="external",
    chip=chip,
    tspan=tspan_He,
    tspan_bg=tspan_air,
)
mdict["He"] = He

CO = point_calibration(
    data_MS, mol="CO", cal_type="external", chip=chip, tspan=tspan_CO, tspan_bg=tspan_He
)
mdict["CO"] = CO

H2_2 = point_calibration(
    data_MS,
    mol="H2",
    cal_type="external",
    chip=chip,
    tspan=tspan_H2,
    tspan_bg=tspan_air,
)
# fuck that, for some reason that looks even worse.

save_calibration_results(mdict, "20A31_calibration_results.pkl")


mdict, ax_F_vs_f = recalibrate(
    internal=[mdict["CO2_M44"], mdict["O2_M32"], mdict["H2"]],
    external=[He, CO, H2_2],
)


fig_0 = ax_0[0].get_figure()
fig_0.set_figwidth(fig_0.get_figwidth() * 2.5)
fig_0.savefig("calibration_experiment.svg")
ax_cal_curve_H2.get_figure().savefig("H2 calibration curve.svg")
ax_F_vs_f.get_figure().savefig("F vs f.png")
