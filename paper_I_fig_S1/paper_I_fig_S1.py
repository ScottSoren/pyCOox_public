"""This script makes fig S1 of
Tracking Oxygen Atoms in Electrochemical CO Oxidation Part I...
(https://doi.org/10.1016/j.electacta.2021.137842)

It requires the relevant data file to be copied to ../data,
and for ixdat to be installed.
(see the README at https://github.com/ScottSoren/pyCOox_public)

The figure shows the calibration experiment in which the mass spec sensitivity factors,
reference electrode value, and working distance were obtained for use in both papers.

This script goes through all of the calculations and plotting, and can be used as a
tutorial on EC-MS calibration using ixdat. For the theory of quantitative mass
spectrometry employed here, see Chapter 2.2 of S. Scott, 2019:
https://orbit.dtu.dk/en/publications/isotope-labeling-studies-in-electrocatalysis-for-renewable-energy

The script produces an ixdat calibration file which is used by all the other scripts in
this repository: pyCOox_cal.ix, which is a JSON-formatted plain text file.

TODO:
    parts of the script require ixdat 0.1.0, which will be ready later today
"""
import numpy as np
from ixdat import Measurement
from ixdat.techniques.ms import Chip, MSCalResult
from ixdat.techniques.ec_ms import ECMSCalibration
from ixdat.constants import Faraday

# ---------------- load and plot the raw data ------------------------ #
# All the data used for calibration is here:
meas = Measurement.read("../data/01_Pt_in_18O_electrolyte.pkl", reader="EC_MS")

# We know the geometric area of the electrode, so we can normalize the current:
meas.normalize_current(A_el=0.196)  # its a 5mm diameter disk, area = 0.196 [cm^2]

# the part that we're interested in now starts rather late, so for convenience,
# lets shift t=0 to the right:
meas.tstamp += 12250

# subfigure a is a plot of the whole calibration experiment.
# The raw data is easily plotted with the plot_measurement function.
# See help(meas.plot_measurement) for customization.
axes_a = meas.plot_measurement(tspan=(0, 7500))
# A list of three axes for (i) ms, (ii) potential, and (iii) current is returned
# We will be adding to these axes to indicate the timespans used for the calibrations.

# It's nice to have this plot really wide:
fig_a = axes_a[0].get_figure()  # name the figure by getting it via the first axis.
fig_a.set_figwidth(
    fig_a.get_figwidth() * 2.5
)  # make it 2.5 times as wide as it is tall

if False:
    # save the figure as it is, without the highlights:
    fig_a.savefig("paper_I_fig_S1a_simple.png")

# --------------- The RE and working distance calibrations ---------------------- #

# The potential of the RHE with respect to the reference electrode is the steady raw
# potential ("Ewe/V" or, equivalently, "raw potential / [V]") while the electrolyte
# is saturated with hydrogen and there is zero current (OCP). This happens between
# 6200 s and 6300 s in the measurement. We use the "grab" function:
t_RHE, v_RHE = meas.grab("raw potential / [V]", tspan=(6200, 6300))
RE_vs_RHE = -np.mean(v_H2)

# Let's plot it as a thicker line to make it clear which value we've used:
axes_a[1].plot(t_RHE, v_RHE, color="black", linewidth=1)

# The working distance can be calculated by the mass-transport-limited HOR current.
# This is the raw current ("I/mA" or, equivalently, "raw current / [mA]") when the
# electrolyte is H2-saturated and the potential is in the double-layer region:
t_HOR, I_HOR_mA = meas.grab("raw current / [mA]", tspan=(7000, 7060))
I_lim = np.mean(I_HOR_mA) * 1e-3  # take the average and convert it from mA to [A]
J_lim = I_lim / (A_el * 1e-4)  # the geometric current density in [A/m^2]
# (the factor 1e-4 converts A_el from cm^2 to m^2)

# The working distance can be calculated from the HER limiting current as follows:
D_H2 = 4.5e-9  # diffusion constant of H2 in water, in [m^2/s]
c_sat_H2 = 0.78  # concentration of H2 in water saturated with 1 bar H2, in [mol/m^3]

# Then, we use eq. 10 from the SI to calculate the working distance:
L = 2 * Faraday / J_lim * D_H2 * C_sat_H2  # the working distance in [m]

# Let's plot it as a thicker line to make it clear which value we've used:
j_HOR = I_HOR_mA / A_el  # the current we used, here in [mA/cm^2] to match the plot.
axes_a[2].plot(t_HOR, j_HOR, color="red", linewidth=1)

# ---------------- Calibrate for H2 at m/z=2 -------------------------- #

# subfigure b is the only calibration curve in the dataset, for H2.
# For this we use the ixdat function which automates this:
cal_result_H2, ax_b, _ = meas.ecms_calibration_curve(
    mol="H2",
    mass="M2",
    n_el=-2,
    tspan=[(5000, 5050), (5260, 5280), (5475, 5500)],  # timespans with H2 evolution
    tspan_bg=(4850, 4900),
    ax="new",
    axes_measurement=axes_a,
)
# Because we asked for two types of plotting, three arguments are returned.
# The first, which we call cal_result_H2, is a MSCalResult for H2.
# The second, which we clal ax_H2, is the axis where it plots the calibration curve.
# The third is the axes where it plots the integrals, but we already have this as axes_a.

# The attribute cal_result_H2.F is the slope of the calibration curve,
# which is the sensitivity factor in [C/mol]:
print(cal_result_H2)  # prints: MSCalResult(name="H2_M2", mol="H2", mass="M2", F=1.9)

# This is the first of several MSCalResults, which we will collect in this list:
cal_results = [cal_result_H2]

# --------------- Calibrate for O2 at m/z=32, 34, and 36 --------------------- #

# Because we need to add the integrated signals at these three isotopes to calculate a
# sensitivity factor for O2, the quick-and-easy ECMSMeasurement.ecms_point_calibration
# function won't work. Instead, we do it manually:

S_int = 0  # we will add the integral at each mass to this
tspan_OER = (3750, 4350)  # a timespan in the experiment with steady OER
for mass in ["M32", "M34", "M36"]:
    S_int += meas.integrate_signal(
        mass=mass,
        tspan=tspan_OER,
        t_bg=(2400, 2500),  # a timespan in the experiment to use as a background
        ax=axes_a[0],  # the axis on which to indicate the integral
    )  # The integral is returned in [A * s] = [C]
Q = (
    meas.integrate(
        "raw_current / [mA]",
        tspan=tspan_OER,
        ax=axes_a[2],  # axes_a[2] is the current axis.
    )
    * 1e-3
)  # Charge passed in [C]. (factor 1e-3 converts mC to C)
n = Q / (4 * Faraday)  # amount of O2 produced by Faraday's law of electrolysis in [mol]

# The sensitivity factor is the ratio of the integrated signal to the amount of O2:
F_O2_Mx = S / n  # sensitivity factor in [C/mol]

# We assume that this sensitivity factor is the same for each isotope, and manually
# define the calibration results:
for mass in ["M32", "M34", "M36"]:
    cal = MSCalResult(
        mass=mass,
        mol=f"O2_{mass}",
        F=F_O2_Mx,
        cal_type="internal",
    )
    cal_results.append(cal)

# --------------- Calibrate for CO2 at m/z=44, 46, and 48 --------------------- #

# Same procedure as for O2 above, just different times and masses and 2 e- instead of 4

S_int = 0  # we will add the integral at each mass to this
tspan_COox = (339, 582)  # a timespan in the experiment with steady CO oxidation
for mass in ["M44", "M46", "M48"]:
    S_int += meas.integrate_signal(
        mass=mass,
        tspan=tspan_COox,
        t_bg=(2400, 2500),  # a timespan in the experiment to use as a background
        ax=axes_a[0],  # the axis on which to indicate the integral
    )  # The integral is returned in [A * s] = [C]
Q = (
    meas.integrate("raw_current", tspan=tspan_COox, ax=axes_a[2]) * 1e-3
)  # Charge passed in [C]. (factor 1e-3 converts mC to C)
n = Q / (
    2 * Faraday
)  # amount of CO2 produced by Faraday's law of electrolysis in [mol]

# The sensitivity factor is the ratio of the integrated signal to the amount of CO2:
F_O2_Mx = S / n  # sensitivity factor in [C/mol]

# We assume that this sensitivity factor is the same for each isotope, and manually
# define the calibration results:
for mass in ["M44", "M46", "M48"]:
    cal = MSCalResult(
        name=f"CO2_{mass}",
        mass=mass,
        mol="CO2",
        F=F_O2_Mx,
        cal_type="internal",
    )
    cal_results.append(cal)

# --------------- Calibrate for the pure gasses, CO, He, and H2 --------------------- #

chip = Chip()
# The chip object contains functions for calculating the molar flux of a pure gas (the
#   carrier gas to the capillary, and for using this to calculate a sensitivity factor.
#   In principle, the chip's capillary should be calibrated, but this usually only
#   changes the results by a few percent, and is not yet implemented in ixdat.

cal_CO = chip.carrier_gas_calibration(
    measurement=meas,  # the measurement with the calibration data
    mol="CO",  # the molecule to calibrate
    mass="M28",  # the mass to calibrate at
    tspan=(110, 140),  # the timespan to average the signal over
    ax=ax[0],  # the axis on which to indicate what signal is used with a thicker line
)
cal_He = chip.carrier_gas_calibration(
    measurement=meas, mol="He", mass="M4", tspan=(2200, 2400), ax=ax[0]
)
cal_H2_2 = chip.carrier_gas_calibration(
    measurement=meas, mol="H2", mass="M2", tspan=(6900, 7000), ax=ax[0]
)
cal_H2_2.name = "H2_M2_carrier"  # so that it doesn't clash with the EC-MS calibration.

cal_results += [cal_CO, cal_He, cal_H2_2]

# ------------------- Analyzing and saving the calibration ------------------- #

# Now we put all these results together in an ECMSCalibration object:
calibration = ECMSCalibration(
    name="Scott2021a_calibration",  # I'm naming this with the article's citation key.
    date="20A31",  # date of the calibration measurements (shorthand for 2020 January 31)
    setup="DTU sniffer",  # setup where calibration was made (the original chip EC-MS)
    ms_cal_results=cal_results,  # the mass spec calibrations
    RE_vs_RHE=RE_vs_RHE,  # the RE potential in [V]
    A_el=A_el,  # the geometric electrode area in [cm^2]
    L=L,  # the working distance in [m]
)

if False:  # Change this to True to make subfigure c
    # The code for the x-axis of subfigure c (the "F-vs-f" plot) is not yet implemented
    # in ixdat, so this requires the old EC_MS package
    from EC_MS import Molecule, recalibrate

    mdict = {}  # EC_MS works with a dictionary of Molecule objects.
    for cal in calibration.ms_cal_results:
        m = Molecule(cal.mol, primary=cal.mass, F_cal=cal.F)
        mdict[cal_point.name] = m

    _, ax_F_vs_f = recalibrate(
        internal=[mdict["CO2_M44"], mdict["O2_M32"], mdict["H2"]],
        external=[mdict["He"], CO, H2_2],
    )
    ax_F_vs_f.get_figure().savefig("F vs f.svg")


calibration.export()  # this creates Scott2021a_calibration.ix in the script's folder.
