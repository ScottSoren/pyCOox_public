"""
This script makes each of the panels in Figure S1 of
Tracking Oxygen Atoms in Electrochemical CO Oxidation Part II...
(https://doi.org/10.1016/j.electacta.2021.137844)

It requires the relevant data file to be copied to ../data,
and for ixdat to be installed.
(see the README at https://github.com/ScottSoren/pyCOox_public)

This script requires ixdat 0.1.1+ (available as of April 21, 2021)

It also requires the calibration file Scott2021a_ElectrocmimActa_calibration.ix
included in the repository in ../paper_I_fig_S1

The attentive user may notice that the calibratrion is corrected before use, using an
internal calibration. This is because this measurement was taken one week before the
calibration file, and MS sensitivity drifts. This is considered normal.

More concerning is the fact that the measured M34/M32 ratio after the sample should have
been reductively de-labeled seems a bit high compared to the nominal 0.0040, even though
the M46/M44 ratio immediately  seems right on. This is unexplained.
The nominal 0.0040 is used.
"""


import numpy as np
from ixdat import Measurement
from ixdat.techniques.ec_ms import ECMSCalibration
from ixdat.plotters.ms_plotter import STANDARD_COLORS


# check which version of ixdat we're running to know the order of EC-MS axes
#    and how the isotopic molecules are named (with "_" or "@" between mol and mass):
import ixdat
if ixdat.__version__.startswith("0.1"):
    index_j_ax = 2
    index_O16_ax = 3
    sep = "_"
elif ixdat.__version__.startswith("0.2"):
    index_j_ax = 3
    index_O16_ax = 2
    sep = "@"
else:
    raise ImportError(f"This script doesn't run with ixdat version={ixdat.__version__}")


# -------- load and select the data ---------- #
meas_full = Measurement.read("../data/04_Pt18Ox.pkl", reader="EC_MS")
meas = meas_full.cut(tspan=[100, 1900], t_zero="start")

# ------- fig 1a is a simple plot of the raw data. This is a one-liner in ixdat. -----
axes_a = meas.plot_measurement()
# and, to save it:
fig_a = axes_a[0].get_figure()
fig_a.savefig("paper_II_fig_1a.png")

# ------ fig 1b is an isotope experiment plot, with calibrated data on two axes ------ #

# First, we load, attach, and correct the calibration
calibration = ECMSCalibration.read(
    "../paper_I_fig_S1/Scott2021a_ElectrocmimActa_calibration.ix"
)
tspan_natural = [5000, 5200]  # steady OER from de-labeled catalyst
# tspan_bg = [3850, 3900]  # just before that OER period, signals at background level
tspan_bg = [8000, 8500]  # after that OER period, signals at background level
O2_cal_result = meas_full.ecms_calibration(
    mol="O2", mass="M32", n_el=4, tspan=tspan_natural, tspan_bg=tspan_bg
)
meas.calibration = calibration.scaled_to(O2_cal_result)

# -- First, we do this semi-manually to show the analysis explicitly:

# get the natural ratio from later in the experiment,
# where the 18O has already been reduced out of the Pt(18)Ox:
meas_full.plot_measurement()

if False:  # this seems to give too high beta :( ... can't tell why.
    t_M32, S_M32 = meas_full.grab_signal("M32", tspan=tspan_natural, t_bg=tspan_bg)
    t_M34, S_M34 = meas_full.grab_signal("M34", tspan=tspan_natural, t_bg=tspan_bg)
    beta = np.mean(S_M34) / np.mean(S_M32)  # The natural ratio (Equation S3)
else:
    beta = 0.0040  # the nominal natural ratio of 18O16O to 16O2

# then we set the background:
meas.set_bg(
    mass_list=["M2", "M32", "M34", "M36", "M44", "M46", "M48"], tspan_bg=[640, 680]
)

# make the plot:
axes_b = meas.plot_measurement(
    mol_lists=[
        [f"O2{sep}M34", f"O2{sep}M36", f"CO2{sep}M46", f"CO2{sep}M48"],  # 18-O on left y-axis
        [f"O2{sep}M32", f"CO2{sep}M44"],  # 16-O only on right -yaxis
    ],
    removebackground=True,
    unit="pmol/s/cm^2",
    logplot=False,
    legend=False,
)

# scale the axes according to the expected isotopic ratio:
ylim_O18 = np.array([-0.18, 18])
ylim_O16 = ylim_O18 / beta
axes_b[0].set_ylim(ylim_O18)
axes_b[index_O16_ax].set_ylim(ylim_O16)

# and make the highlights:
t, y_34_mol_s = meas.grab_flux(mol=f"O2{sep}M34", removebackground=True)
y_34 = y_34_mol_s * 1e12 / meas.A_el  # convert from [mol/s] to [pmol/cm^2]
y_32_mol_s = meas.grab_flux_for_t(mol=f"O2{sep}M32", t=t, removebackground=True)
y_32 = y_32_mol_s * 1e12 / meas.A_el  # convert from [mol/s] to [pmol/cm^2]
axes_b[0].fill_between(t, y_32 * beta, y_34, color="r", alpha=0.2)

y_46_mol_s = meas.grab_flux_for_t(mol=f"CO2{sep}M46", t=t, removebackground=True)
y_46 = y_46_mol_s * 1e12 / meas.A_el  # convert from [mol/s] to [pmol/cm^2]
y_44_mol_s = meas.grab_flux_for_t(mol=f"CO2{sep}M44", t=t, removebackground=True)
y_44 = y_44_mol_s * 1e12 / meas.A_el  # convert from [mol/s] to [pmol/cm^2]
axes_b[0].fill_between(t, y_44 * beta, y_46, color=STANDARD_COLORS["M46"], alpha=0.2)

# format and save the plot:
axes_b[0].set_ylabel("signal / [pmol s$^{-1}$cm$^{-2}$]")
axes_b[index_O16_ax].set_ylabel("signal / [pmol s$^{-1}$cm$^{-2}$]")
fig_b = axes_b[0].get_figure()
fig_b.savefig("paper_II_fig_1b.png")


# ---------- calculate the excess lattice O evolved amounts ---------- #

tspan_exchange = [100, 430]  # OER on Pt(18)Ox
tspan_extraction = [1050, 1300]  # CO ox while reducing Pt(18)Ox in CO

mask_exchange = np.logical_and(tspan_exchange[0] < t, t < tspan_exchange[-1])
mask_extraction = np.logical_and(tspan_extraction[0] < t, t < tspan_extraction[-1])

n_exchange = np.trapz(
    (y_34 - y_32 * beta)[mask_exchange], t[mask_exchange]
)  # in [pmol/cm^2]
n_extraction = np.trapz(
    (y_46 - y_34 * beta)[mask_extraction], t[mask_extraction]
)  # in [pmol/cm^2]
print(f"exces 18-O in O2 during OER: {n_exchange} [pmol/cm^2]")
print(f"exces 18-O in CO2 during COox: {n_extraction} [pmol/cm^2]")
