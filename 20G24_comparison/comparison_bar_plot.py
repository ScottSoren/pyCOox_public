from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from EC_MS.utils.extraction_class import Extraction
from EC_MS.converters import mdict_from_SI2020_calibration_file

plt.close("all")

forpublication = False
if forpublication:  # for the publication figure
    import matplotlib as mpl
    figwidth, figheight = 3.25, 2.75
    mpl.rcParams['figure.figsize'] = (figwidth, figheight)
    # plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')
    plt.rc('font', size=8)
    plt.rc('lines', linewidth=0.5)
else:
    plt.style.use("default")

EXTRACTION_DIR = Path(__file__).parent.parent / "extractions"
DATA_DIR = Path(__file__).parent.parent.parent / "pickles"
CALIBRATION_DIR = Path(__file__).parent.parent

extraction_specs = {
    "PtOx":
        {"file":"scanned.json", "position":1, "t_start":210, "tspan_ratio":[490, 500], "fancy_name":"Pt$^{18}$O$_x$",},
    "IrO2":
        {"file":"Decade1C_all.json", "position":2, "t_start":457,"tspan_ratio":[800, 1050], "fancy_name":"Ir$^{18}$O$_2$",},
    "IrOx":
        {"file":"Jazz8_all.json", "position":3, "t_start":424, "fancy_name":"Ir$^{18}$O$_x$",},
    "IrOx_hyd":
        {"file":"Jazz8b_all.json", "position":4, "t_start":167, "fancy_name":"Ir$^{18}$O$_x\cdot y$H$_2$O",},
    }

calibration_dir = Path("..").absolute()
calibration_file = calibration_dir / "20A25_sniffer_fixed.json"
mdict = mdict_from_SI2020_calibration_file(calibration_file)


fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

extractions = {}

xticks = []
xticklabels = []

for name, spec in extraction_specs.items():
    path_to_extraction = EXTRACTION_DIR / spec["file"]

    extraction = Extraction.load(path_to_extraction, data_dir=DATA_DIR, mdict=mdict)

    if "tspan_ratio" in spec:
        extraction.get_alpha(tspan=spec["tspan_ratio"], ax=None)
        extraction.tspan_ratio = spec["tspan_ratio"]

    if True:  # correct calibration
        tspan_steady = extraction.tspan_ratio  # extraction.tspan_ratio
        O2 = extraction.point_calibration(mol="O2", mass="M32", tspan=tspan_steady, n_el=4)
        # since we use cal_mat but point_calibration gives F_cal, this is the comparison
        # to make. We should then adjust such that it makes sense.
        correction = O2.F_cal / extraction.mdict["O2_M32"].F_cal
        extraction.correction = correction
        for m, molecule in extraction.mdict.items():
            if hasattr(molecule, "cal_mat"):
                for mass, value in molecule.cal_mat.items():
                    molecule.cal_mat[mass] = value / correction
            if hasattr(molecule, "F_cal"):
                molecule.F_cal = molecule.F_cal * correction

    extractions[name] = extraction
    tspan_plot = extraction.tspan_exchange

    if True and not forpublication:
        ax = extraction.plot_exchange(tspan=tspan_plot, mol="O2")
        extraction.plot_exchange(tspan=tspan_plot, mol="CO2", axes=ax)

    O2 = extraction.mdict["O2_M32"]
    CO2 = extraction.mdict["CO2_M44"]

    O2_excess_M34 = extraction.create_excess_mol("O2")
    O2_M36 = extraction.mdict["O2_M36"]

    CO2_excess_M46 = extraction.create_excess_mol("CO2")
    CO2_M48 = extraction.mdict["CO2_M48"]

    t_int = 300

    tspan_int = [spec["t_start"], spec["t_start"]+t_int]
    position = spec["position"]

    totals = {}
    width = 0
    f_hyd_1 = 1 / ( 1 - 0.22)
    f_hyd_2 = 1 / ( 1 - 2*0.22)
    width = 1/8

    for (m, offset, ax, factor) in [
        (O2_excess_M34, -2, ax1, 1),
        (O2_M36, -2, ax1, 2),
        (CO2_excess_M46, -1, ax1, 1*f_hyd_1),
        (CO2_M48, -1, ax1, 2*f_hyd_2),
        (O2, 1, ax2, 1),
        (CO2, 2, ax2, 2),
    ]:
        x, y = extraction.get_flux(m, tspan=tspan_int, unit="pmol/s/cm^2")
        n = np.trapz(y, x) * factor
        n_dot = n / t_int  # in pmol/s

        if offset not in totals:
            totals[offset] = 0
        bottom = totals[offset]
        totals[offset] += n_dot

        color = m.get_color()
        pos = position + offset * width

        ax.bar(pos, n_dot, width=width, bottom=bottom, color=color)

    xticks += [position]
    xticklabels += [spec["fancy_name"]]

ylim2 = [0, 660]
ylim1 = [lim*0.02 for lim in ylim2]

ax1.set_xticks(xticks)
ax1.set_xticklabels(xticklabels)

ax1.set_ylim(ylim1)
ax2.set_ylim(ylim2)

ax1.set_ylabel("lattice $^{18}$O / [pmol s$^{-1}$cm$^{-2}$]")
ax2.set_ylabel("total product / [pmol s$^{-1}$cm$^{-2}]$")

ax1.set_title("average rates in first 300 s of electrolysis")

if forpublication:
    fig.set_figwidth(figwidth)
    fig = ax1.get_figure()
    fig.savefig("comparision.png")
    fig.savefig("comparision.svg")