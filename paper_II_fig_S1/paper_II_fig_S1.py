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
electrode area here, but not in the publication. They are also slightly different
than in the publication. I think this was a mistake in the publication,
that I had used a different, less accurate, calibration. The mistake did not affect the
conclusions of the publication.
"""

from ixdat import Measurement
from ixdat.techniques.ec_ms import ECMSCalibration


# ---------- load calibration and data, find data selections -------------------- #


calibration = ECMSCalibration.read(
    "../paper_I_fig_S1/Scott2021a_ElectrocmimActa_calibration.ix"
)

meas = Measurement.read(
    "../data/03_Ir_in_18O_electrolyte.pkl",
    reader="EC_MS",
    calibration_list=[calibration]
)
if False:  # plot the whole raw data to get an overview
    meas.plot_measurement()

meas.set_bg(  # define background just once for the whole script:
    tspan_bg=[4600, 4620],
    mass_list=["M2", "M32", "M34", "M36", "M44", "M46", "M48"],
)

tspan_CO_ox = [5680, 6300]
tspan_CO_strip = [4450, 5200]

if False:  # plot raw data
    for tspan in [tspan_CO_ox, tspan_CO_strip]:
        meas.plot_measurement(tspan=tspan)


# --- panel a: EC-MS plot with carrier gases on left y-axis and products on right --- #

meas_a = meas.cut(tspan=tspan_CO_ox)
meas_a.tstamp += meas_a.t[0]

axes_a = meas_a.plot_measurement(
    mol_lists=[
        ["He", "CO"],  # left y-axis
        ["H2", "O2@M32", "O2@M34", "O2@M36", "CO2@M44", "CO2@M46", "CO2@M48"],  # right
    ],
    unit=["nmol/s", "pmol/s/cm^2"],  # [left, right] y-axes,
    remove_background=[False, True],
    logplot=False,
    legend=False,
)
axes_a[0].set_ylabel("cal. sig. / [nmol s$^{-1}$]")
axes_a[-1].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]")

fig_a = axes_a[0].get_figure()
fig_a.savefig("paper_II_fig_S1a.png")

# --- panel b: cyclic voltammatry MS plot of two of the cycles in (a) --- #

meas_b = meas_a.as_cv()
meas_b.redefine_cycle(start_potential=0.2, redox=1)
axes_b = meas_b[1].plot_vs_potential(
    mol_list=["H2", "O2@M32", "O2@M34", "O2@M36", "CO2@M44", "CO2@M46", "CO2@M48"],
    unit="pmol/s/cm^2",
    logplot=False,
    legend=False,
    remove_background=True,
)
meas_b[3].plot_vs_potential(
    mol_list=["H2", "O2@M32", "O2@M34", "O2@M36", "CO2@M44", "CO2@M46", "CO2@M48"],
    unit="pmol/s/cm^2",
    logplot=False,
    legend=False,
    remove_background=True,
    axes=axes_b,
    linestyle="--",
)
axes_b[0].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]")
axes_b[0].set_xlabel(meas_b.v_name)

fig_b = axes_b[0].get_figure()
fig_b.savefig("paper_II_fig_S1b.png")


# --- panel c: EC-MS plot with carrier gases on left y-axis and products on right --- #

meas_c = meas.cut(tspan=tspan_CO_strip)
meas_c.tstamp += meas_c.t[0]

axes_c = meas_c.plot_measurement(
    mol_lists=[
        ["He", "CO"],  # left y-axis
        ["H2", "O2@M32", "O2@M34", "O2@M36", "CO2@M44", "CO2@M46", "CO2@M48"],  # right
    ],
    unit=["nmol/s", "pmol/s/cm^2"],  # [left, right] y-axes
    remove_background=[False, True],
    logplot=False,
    legend=False,
)
axes_c[0].set_ylabel("cal. sig. / [nmol s$^{-1}$]")
axes_c[-1].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]")
axes_c[-1].set_ylim(bottom=-5)

fig_c = axes_c[0].get_figure()
fig_c.savefig("paper_II_fig_S1c.png")

# --- panel d: cyclic voltammatry MS plot of two of the cycles in (c) --- #

meas_d = meas_c.as_cv()
meas_d.redefine_cycle(start_potential=0.15, redox=0)
axes_d = meas_d[2].plot_vs_potential(
    mol_list=["H2", "O2@M32", "O2@M34", "O2@M36", "CO2@M44", "CO2@M46", "CO2@M48"],
    unit="pmol/s/cm^2",
    logplot=False,
    legend=False,
    remove_background=True,
)
meas_d[3].plot_vs_potential(
    mol_list=["H2", "O2@M32", "O2@M34", "O2@M36", "CO2@M44", "CO2@M46", "CO2@M48"],
    unit="pmol/s/cm^2",
    logplot=False,
    legend=False,
    remove_background=True,
    axes=axes_d,
    linestyle="--",
)
axes_d[0].set_ylabel("cal. sig. / [pmol s$^{-1}$cm$^{-2}$]")
axes_d[0].set_xlabel(meas_b.v_name)
axes_d[0].set_ylim([-5, 50])

fig_d = axes_d[0].get_figure()
fig_d.savefig("paper_II_fig_S1d.png")
