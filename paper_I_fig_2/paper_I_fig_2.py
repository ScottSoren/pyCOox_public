"""
This script makes each of the panels in Figure 2 of
Tracking Oxygen Atoms in Electrochemical CO Oxidation Part I...
(https://doi.org/10.1016/j.electacta.2021.137842)

It requires the relevant data file to be copied to ../data,
and for ixdat to be installed.
(see the README at https://github.com/ScottSoren/pyCOox_public)

This script requires ixdat 0.1+ (available as of April 15, 2021)

It also requires the calibration file Scott2021a_ElectrocmimActa_calibration.ix
included in the repository in ../paper_I_fig_S1

The attentive user may notice that the value of the calibrated signals are slightly
different here than in the publication. I think this was a mistake in the publication,
that I had used a different, less accurate, calibration. The mistake did not affect the
conclusions of the publication.
"""

from ixdat import Measurement
from ixdat.techniques.ec_ms import ECMSCalibration


# check which version of ixdat we're running to know the order of EC-MS axes
#    and how the isotopic molecules are named (with "_" or "@" between mol and mass):
import ixdat
if ixdat.__version__.startswith("0.1"):
    index_j_ax = 2
    sep = "_"
elif ixdat.__version__.startswith("0.2"):
    index_j_ax = 3
    sep = "@"
else:
    raise ImportError(f"This script doesn't run with ixdat version={ixdat.__version__}")


# ---------- load calibration and data, find data selections -------------------- #

meas = Measurement.read(
    "../data/01_Pt_in_18O_electrolyte.pkl", reader="EC_MS",
)
meas.calibration = ECMSCalibration.read(
    "../paper_I_fig_S1/Scott2021a_ElectrocmimActa_calibration.ix"
)

if False:  # plot the whole raw data to get an overview
    meas.plot_measurement()

tspan_CO_strip = [9700, 10500]
tspan_CO_ox = [11350, 12150]

if False:  # plot raw data
    for tspan in [tspan_CO_ox, tspan_CO_strip]:
        meas.plot_measurement(tspan=tspan)


# --- panel a: EC-MS plot with carrier gases on left y-axis and products on right --- #

meas_a = meas.cut(tspan=tspan_CO_ox)
meas_a.tstamp += meas_a.t[0]

axes_a = meas_a.plot_measurement(
    mol_lists=[
        ["He", "CO"],  # left y-axis
        ["H2", f"O2{sep}M32", f"O2{sep}M34", f"O2{sep}M36", f"CO2{sep}M44", f"CO2{sep}M46", f"CO2{sep}M48"],  # right
    ],
    tspan_bg=[
        None,
        [0, 20],
    ],  # [left, right] y-axes
    unit=["nmol/s", "pmol/s/cm^2"],  # [left, right] y-axes
    logplot=False,
    legend=False,
)
axes_a[0].set_ylabel("cal. sig. / [nmol s$^{-1}$]")
axes_a[-1].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]")

fig_a = axes_a[0].get_figure()
fig_a.savefig("paper_I_fig_2a.png")

# --- panel b: cyclic voltammatry MS plot of two of the cycles in (a) --- #

meas_b = meas_a.as_cv()
meas_b.redefine_cycle(start_potential=0.4, redox=1)
axes_b = meas_b[1].plot_vs_potential(
    mol_list=[f"H2", f"O2{sep}M32", f"O2{sep}M34", f"O2{sep}M36", f"CO2{sep}M44", f"CO2{sep}M46", f"CO2{sep}M48"],
    unit="pmol/s/cm^2",
    logplot=False,
    legend=False,
    removebackground=True,
)
meas_b[3].plot_vs_potential(
    mol_list=["H2", f"O2{sep}M32", f"O2{sep}M34", f"O2{sep}M36", f"CO2{sep}M44", f"CO2{sep}M46", f"CO2{sep}M48"],
    unit="pmol/s/cm^2",
    logplot=False,
    legend=False,
    removebackground=True,
    axes=axes_b,
    linestyle="--",
)
axes_b[0].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]")
axes_b[0].set_xlabel(meas_b.V_str)

fig_b = axes_b[0].get_figure()
fig_b.savefig("paper_I_fig_2b.png")


# --- panel c: EC-MS plot with carrier gases on left y-axis and products on right --- #

meas_c = meas.cut(tspan=tspan_CO_strip)
meas_c.tstamp += meas_c.t[0]

meas_c.set_bg(
    tspan_bg=[0, 20], mass_list=["M2", "M32", "M34", "M36", "M44", "M46", "M48"]
)  # these background values carry over to the ECMSCyclicVoltammagram meas_d :D
axes_c = meas_c.plot_measurement(
    mol_lists=[
        ["He", "CO"],  # left y-axis
        ["H2", f"O2{sep}M32", f"O2{sep}M34", f"O2{sep}M36", f"CO2{sep}M44", f"CO2{sep}M46", f"CO2{sep}M48"],  # right
    ],
    unit=["nmol/s", "pmol/s/cm^2"],  # [left, right] y-axes
    logplot=False,
    legend=False,
)
axes_c[0].set_ylabel("cal. sig. / [nmol s$^{-1}$]")
axes_c[-1].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]")
axes_c[-1].set_ylim(bottom=-5)

fig_c = axes_c[0].get_figure()
fig_c.savefig("paper_I_fig_2c.png")

# --- panel d: cyclic voltammatry MS plot of two of the cycles in (c) --- #

meas_d = meas_c.as_cv()
meas_d.redefine_cycle(start_potential=0.38, redox=0)
axes_d = meas_d[2].plot_vs_potential(
    mol_list=["H2", f"O2{sep}M32", f"O2{sep}M34", f"O2{sep}M36", f"CO2{sep}M44", f"CO2{sep}M46", f"CO2{sep}M48"],
    unit="pmol/s/cm^2",
    logplot=False,
    legend=False,
    removebackground=True,
)
meas_d[3].plot_vs_potential(
    mol_list=["H2", f"O2{sep}M32", f"O2{sep}M34", f"O2{sep}M36", f"CO2{sep}M44", f"CO2{sep}M46", f"CO2{sep}M48"],
    unit="pmol/s/cm^2",
    logplot=False,
    legend=False,
    removebackground=True,
    axes=axes_d,
    linestyle="--",
)
axes_d[0].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]")
axes_d[0].set_xlabel(meas_b.V_str)

fig_d = axes_d[0].get_figure()
fig_d.savefig("paper_I_fig_2d.png")
