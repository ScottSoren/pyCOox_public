
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from EC_MS.utils.extraction_class import Extraction
from EC_MS.converters import mdict_from_SI2020_calibration_file

plt.close("all")
if False:  # old wierd import during Spectro Inlet employment
    from spitze.quant import Calibration
    cal_dir = Path("..").absolute()
    c = Calibration.load("20A25_sniffer_fixed.json", cal_dir=cal_dir)
    c.F = c.F_0
    #  I have no fucking idea why this is necessary. Calibraiton.__init__ broke

calibration_dir = Path("../..").absolute()
calibration_file = calibration_dir / "20A25_sniffer_fixed.json"

mdict = mdict_from_SI2020_calibration_file(calibration_file)


forpublication = True
if forpublication:  # for the publication figure
    import matplotlib as mpl
    mpl.rcParams['figure.figsize'] = (3.25, 2.75)
    # plt.rc('text', usetex=True)  # crashingly slow
    plt.rc('font', family='sans-serif')
    plt.rc('font', size=6)
    plt.rc('lines', linewidth=0.5)
else:
    plt.style.use("default")


extraction_dir = Path(r"../extractions").absolute()
data_dir = Path(r"../../pickles").absolute()

file_A = extraction_dir / "scanned.json"
file_B = extraction_dir / "scanned_B.json"

extraction_A = Extraction.load(
        file_A, data_dir=data_dir, mdict=mdict,
) # Pt^{18}O in H2^{16}O
extraction_B = Extraction.load(
        file_B, data_dir=data_dir, mdict=mdict,
) # Pt^{18}O in H2^{16}O
extraction_B.sort_time()

extraction_A.timeshift(100)
extraction_B.timeshift(350)

tspan_experiment = [0, 1800]

t_str = "time / [s]"
 
if True:
    extraction_A.reset()
    axes = extraction_A.plot_experiment(
        masses=["M4", "M28", "M32", "M34", "M36", "M44", "M46", "M48", "M2"], 
        tspan=tspan_experiment, 
        unit="pA",
        logplot=True,
    )
    # axes[0].legend()
    extraction_A.set_background(t_bg=extraction_A.t_bg)
    
    axes[0].set_ylabel("signal / [pA]")
    axes[1].set_xlabel(t_str)
    axes[0].set_xlabel(t_str)
    axes[0].set_yticks([1e-1, 1, 1e1, 1e2, 1e3, 1e4])
    axes[1].set_yticks([0, 0.5, 1, 1.5])
    for ax in axes:
        ax.set_xticks([0, 250, 500, 750, 1000, 1250, 1500, 1750])
    fig = axes[0].get_figure()
    fig.savefig("fig_Pt_extraction_raw.png")
    fig.savefig("fig_Pt_extraction_raw.svg")

# --------- fig 5b, Pt^{18}O in H2^{16}O vs t ----------- "
if True:
    unit_right = "pmol/cm^2/s"
    axes = extraction_A.plot_exchange(
            mol="O2", tspan=tspan_experiment,
            unit_left="pmol/cm^2/s", unit_right=unit_right
    )
    extraction_A.plot_exchange(
            mol="CO2", axes=axes, tspan=tspan_experiment,
            unit_left="pmol/cm^2/s", unit_right=unit_right
    )
    axes[1].set_xlabel(t_str)
    axes[0].set_ylabel("signal / [pmol s$^{-1}$cm$^{-2}$]")
    if unit_right == "nmol/cm^2/s":
        axes[-1].set_ylabel("signal / [nmol s$^{-1}$cm$^{-2}$]")
    elif unit_right == "pmol/cm^2/s":
        axes[-1].set_ylabel("signal / [pmol s$^{-1}$cm$^{-2}$]")
    axes[-1].set_yticks([0, 1000, 2000, 3000, 4000])
    axes[1].set_yticks([0, 0.5, 1, 1.5])
    for ax in axes:
        ax.set_xticks([0, 250, 500, 750, 1000, 1250, 1500, 1750])
    fig = axes[0].get_figure()
    if forpublication:
        fig.savefig("fig_Pt_extraction_A.png")
        fig.savefig("fig_Pt_extraction_A.svg")

    if False:  # get excess O2
        excess_O2 = extraction_A.create_excess_mol("O2")
        x_ex, y_ex = extraction_A.get_flux(excess_O2, tspan=[100, 450])
        fig, ax = plt.subplots()
        ax.plot(x_ex, y_ex)
        Y_ex = np.trapz(y_ex, x_ex)

    if False:  # look at CO2 signals
        tspan_extraction = [1000, 1400]
        CO2_M44 = extraction_A.mdict["CO2_M44"]
        CO2_M46 = extraction_A.mdict["CO2_M46"]
        x1, y1 = extraction_A.get_flux(CO2_M44, tspan=tspan_extraction)
        x2, y2 = extraction_A.get_flux(CO2_M46, tspan=tspan_extraction)
        r_max = max(y2) / (max(y2) + max(y1))


# --------- fig 5c(?), Pt^{18}O in H2^{16}O vs U ----------- "

if False:
    axes = extraction_A.plot_extraction_vs_potential(mol="CO2", tspan=[1000, 1250])

# --------- fig 5b, Pt^{16}O in H2^{18}O vs t ----------- "
if False:
    axes = extraction_B.plot_exchange(mol="O2", tspan=tspan_experiment)
    extraction_B.plot_exchange(mol="CO2", axes=axes, tspan=tspan_experiment)
    axes[1].set_xlabel(t_str)
    fig = axes[0].get_figure()
    if forpublication:
        fig.savefig("fig_Pt_extraction_B.png")
        fig.savefig("fig_Pt_extraction_B.svg")

# --------- fig 5d(?), Pt^{16}O in H2^{18}O vs U ----------- "

if False:
    axes = extraction_B.plot_extraction_vs_potential(mol="CO2", tspan=[1000, 1250])

# --------- fig 6, Pt^{16}O in H2^{18}O vs U ----------- "
if False:
    tspan_extraction = [1000, 1250]
    axes = None
    extraction_A.linestyle = "-"
    extraction_B.linestyle = "--"
    x_lost = 0.222  # portion m/z=46 CO2 lost to m/z=48 CO2 in fig_COox_Pt_20A31.py
    # portion m/z=46 CO2 lost to m/z=44 in 16O electrolyte should be the same :
    extraction_A.x_lost = x_lost
    # portion m/z=44 CO2 lost to m/z=46 in 18O electrolyte should be the double-ish :
    extraction_B.x_lost = 2 * x_lost
    for extraction in [extraction_B,
                       extraction_A
                       ]:
        color = extraction.get_majors_and_minors(mol="CO2")[1][0].get_color()
        excess_CO2 = extraction.create_excess_mol("CO2")
        extraction.excess_CO2 = excess_CO2
        for mass in excess_CO2.cal_mat:
            excess_CO2.cal_mat[mass] *= 1/(1-extraction.x_lost)
        axes = extraction.plot_vs_potential(
            tspan=tspan_extraction, ax=axes, mols=[excess_CO2], logplot=False,
            color=color, linestyle=extraction.linestyle, t_bg=[1340, 1360],
        )


    axes[0].invert_xaxis()
    axes[1].invert_xaxis()
    current_ticks_uA = [-10, -5, 0, 5, 10]
    axes[1].set_yticks([y*1e-3 for y in current_ticks_uA])
    axes[1].set_yticklabels([str(y) for y in current_ticks_uA])
    axes[1].set_ylabel(r"J / [$\mu$A cm$^{-2}$]")
    axes[0].set_ylabel("lat. O in CO2 / [pmol/s]")

    fig = axes[0].get_figure()
    if forpublication:
        fig.savefig("fig_Pt_extraction_vs_U.png")
        fig.savefig("fig_Pt_extraction_vs_U.svg")

    if True: # calculate excess CO2
        x_A, y_A = extraction_A.get_flux(extraction_A.excess_CO2, tspan=[1100, 1250],
                                         unit="pmol/s", t_bg=[1340, 1360])
        Y_A = np.trapz(y_A, x_A)
        x_B, y_B = extraction_B.get_flux(extraction_B.excess_CO2, tspan=[1125, 1225],
                                         unit="pmol/s", t_bg=[1340, 1360])
        Y_B = np.trapz(y_B, x_B)
