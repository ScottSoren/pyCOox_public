"""
This script makes each of the panels in Figure S1 of
Tracking Oxygen Atoms in Electrochemical CO Oxidation Part II...
(https://doi.org/10.1016/j.electacta.2021.137844)

It requires the relevant data file to be copied to ../data,
and for ixdat to be installed.
(see the README at https://github.com/ScottSoren/pyCOox_public)

This script requires ixdat 0.1+ (available as of April 15, 2021)

It also requires the calibration file Scott2021a_ElectrocmimActa_calibration.ix
included in the repository in ../paper_I_fig_S1

The attentive user may notice that the value of the calibrated signals are normalized to
electrode area here, but not in the publication.
"""
import numpy as np

from ixdat import Measurement
from ixdat.techniques.ec_ms import ECMSCalibration


calibration = ECMSCalibration.read(
    "../paper_I_fig_S1/Scott2021a_ElectrocmimActa_calibration.ix"
)
meas_Pt = Measurement.read(
    "../data/01_Pt_in_18O_electrolyte.pkl",
    reader="EC_MS",
    calibration_list=[calibration],
)
meas_Ir = Measurement.read(
    "../data/03_Ir_in_18O_electrolyte.pkl",
    reader="EC_MS",
    calibration_list=[calibration],
)

if False:  # plot entire measurements
    axes_Pt = meas_Pt.plot_measurement()
    axes_Ir = meas_Ir.plot_measurement()


#  ------------- Platinum reduction ---------------------------- #
meas_a = meas_Pt.cut(tspan=[14350, 17100])
meas_a.tstamp += meas_a.t[0]

meas_a.set_bg(
    tspan_bg=[0, 20],
    mass_list=["M2", "M32", "M34", "M36", "M44", "M46", "M48"],
)
axes_a = meas_a.plot_measurement(  # all on one axis here!
    mol_lists=[
        ["O2@M32", "O2@M34", "O2@M36", "CO2@M44", "CO2@M46", "CO2@M48"],
        [
            "H2",
        ],
    ],
    unit="pmol/s/cm^2",
    logplot=False,
    legend=False,
    remove_background=True,
)
axes_a[0].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]")

fig_a = axes_a[0].get_figure()
fig_a.savefig("paper_II_fig_S2a")

# for the CVs, the two parts to plot don't start at the same potential, so rather than
# using the CyclicVoltammagram selection tools, we just manually cut the right tspans:
meas_b1 = meas_a.cut(tspan=[2333, 2460]).as_cv()
meas_b2 = meas_a.cut(tspan=[2460, 2566]).as_cv()

# and plot these on shared axes:
axes_b = meas_b1.plot_vs_potential(
    mol_list=["H2", "CO2@M44", "CO2@M46", "CO2@M48"],
    unit="pmol/s/cm^2",
    legend=False,
    logplot=False,
    remove_background=True,
)
meas_b2.plot_vs_potential(
    mol_list=["H2", "CO2@M44", "CO2@M46", "CO2@M48"],
    unit="pmol/s/cm^2",
    legend=False,
    logplot=False,
    remove_background=True,
    linestyle="--",
    axes=axes_b,
)
axes_b[0].set_ylim([-5, 55])
axes_b[0].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]")

# what we're missing now is the highlight and the calculation of the charge passed.
# for this, we can use the select_sweep() method to grab the data:
vspan = [1.0, 0.5]
# note this selects cathodic sweep because 1.0 > 0.5
sweep_b1 = meas_b1.select_sweep(vspan=vspan)
sweep_b2 = meas_b2.select_sweep(vspan=vspan)

t1, v = sweep_b1.grab("potential")  # returns the time and potential for that Vspan
j1 = sweep_b1.grab_for_t("current", t1)  # returns norm. current in [mA/cm^2]

t2, v2 = sweep_b2.grab("potential")  # returns the time and potential for that Vspan
j2 = sweep_b2.grab_for_t("current", t2)  # returns norm. current in [mA/cm^2]

# for plotting, we need to ensure all vectors are the same length:
j2_interp = np.interp(-v, -v2, j2)  # j2_interp is same length as v and j1
axes_b[1].fill_between(v, j1, j2_interp, color="c")  # c for cyan

# for integrating, it's easier to use the raw values:
Q1 = np.trapz(j1, t1)  # charge passed in sweep 1 in [mC/cm^2]
Q2 = np.trapz(j2, t2)  # charge passed in sweep 2 in [mC/cm^2]

dQ = Q1 - Q2  # excess charge (negative for cathodic) passed:
print(f"reduction of PtOx involved a net charge of dQ={dQ} [mC/cm^2]")

fig_b = axes_b[0].get_figure()
fig_b.savefig("paper_II_fig_S2b")


#  ------------- Iridium reduction ---------------------------- #
meas_c = meas_Ir.cut(tspan=[7350, 10100])
meas_c.tstamp += meas_c.t[0]

meas_c.set_bg(
    tspan_bg=[0, 20], mass_list=["M2", "M32", "M34", "M36", "M44", "M46", "M48"]
)
axes_c = meas_c.plot_measurement(  # all on one axis here!
    mol_list=["H2", "O2@M32", "O2@M34", "O2@M36", "CO2@M44", "CO2@M46", "CO2@M48"],
    unit="pmol/s/cm^2",
    logplot=False,
    legend=False,
    remove_background=True,
)
axes_c[0].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]")

fig_c = axes_c[0].get_figure()
fig_c.savefig("paper_II_fig_S2c")

# for the CVs, the two parts to plot don't start at the same potential, so rather than
# using the CyclicVoltammagram selection tools, we just manually cut the right tspans:
meas_d1 = meas_c.cut(tspan=[2495, 2607]).as_cv()
meas_d2 = meas_c.cut(tspan=[2607, 2706]).as_cv()

# and plot these on shared axes:
axes_d = meas_d1.plot_vs_potential(
    mol_list=["H2", "CO2@M44", "CO2@M46", "CO2@M48"],
    unit="pmol/s/cm^2",
    legend=False,
    logplot=False,
    remove_background=True,
)
meas_d2.plot_vs_potential(
    mol_list=["H2", "CO2@M44", "CO2@M46", "CO2@M48"],
    unit="pmol/s/cm^2",
    legend=False,
    logplot=False,
    remove_background=True,
    linestyle="--",
    axes=axes_d,
)
axes_d[0].set_ylim([-5, 55])
axes_d[0].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]")

# what we're missing now is the highlight and the calculation of the charge passed.
# for this, we can use the select_sweep() method to grab the data:
vspan = [0.6, 0.05]
# note this selects cathodic sweep because 0.1 < 0.6
sweep_d1 = meas_d1.select_sweep(vspan=vspan)
sweep_d2 = meas_d2.select_sweep(vspan=vspan)

t1, v = sweep_d1.grab("potential")  # returns the time and potential for that Vspan
j1 = sweep_d1.grab_for_t("current", t1)  # returns norm. current in [mA/cm^2]

t2, v2 = sweep_d2.grab("potential")  # returns the time and potential for that Vspan
j2 = sweep_d2.grab_for_t("current", t2)  # returns norm. current in [mA/cm^2]

# for plotting, we need to ensure all vectors are the same length:
j2_interp = np.interp(-v, -v2, j2)  # j2_interp is same length as v and j1
axes_d[1].fill_between(v, j1, j2_interp, color="c")  # c for cyan

# for integrating, it's easier to use the raw values:
Q1 = np.trapz(j1, t1)  # charge passed in sweep 1 in [mC/cm^2]
Q2 = np.trapz(j2, t2)  # charge passed in sweep 2 in [mC/cm^2]

dQ = Q1 - Q2  # excess charge (negative for cathodic) passed:
print(f"reduction of IrOx involved a net charge of dQ={dQ} [mC/cm^2]")

fig_d = axes_d[0].get_figure()
fig_d.savefig("paper_II_fig_S2d")
