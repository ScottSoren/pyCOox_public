import os
import numpy as np
from matplotlib import pyplot as plt
from EC_MS.utils.extraction_class import Extraction

forpublication = False

if forpublication:  # for the publication figure
    import matplotlib as mpl

    mpl.rcParams["figure.figsize"] = (3.25, 2.75)
    plt.rc("text", usetex=True)
    plt.rc("font", family="sans-serif")
    plt.rc("font", size=8)
    plt.rc("lines", linewidth=0.5)
else:
    plt.style.use("default")

extraction_dir = os.path.abspath(
    os.path.join(os.path.split(__file__)[0], "../extractions")
)
data_dir = os.path.abspath(os.path.join(os.path.split(__file__)[0], "../../pickles"))

file = os.path.join(extraction_dir, "Jazz8b_all.json")

extraction = Extraction.load(file, data_dir=data_dir)  # Ir^{18}O2 in H2^{16}O
if True:  # correct H2 signal and ^{16}O^{18}O signal for CO, and ratio
    extraction.mdict["H2"].cal_mat["M4"] = (
        -8.5e-4 * extraction.mdict["H2"].cal_mat["M2"]
    )
    extraction.mdict["O2_M34"].cal_mat["M28"] = (
        -2.3e-5 * extraction.mdict["H2"].cal_mat["M2"]
    )
    extraction.set_background(t_bg=[60, 100])
    # extraction.get_alpha(tspan=[750, 850], ax=None)


if True:  # correct calibration
    tspan_steady = extraction.tspan_ratio  # extraction.tspan_ratio
    O2 = extraction.point_calibration(mol="O2", mass="M32", tspan=tspan_steady, n_el=4)
    # since we use cal_mat but point_calibration gives F_cal, this is the comparison
    # to make. We should then adjust such that it makes sense.
    correction = O2.F_cal * extraction.mdict["O2_M32"].cal_mat["M32"]
    for m, molecule in extraction.mdict.items():
        for mass, value in molecule.cal_mat.items():
            molecule.cal_mat[mass] = value / correction
        if hasattr(molecule, "F_cal"):
            molecule.F_cal = molecule.F_cal * correction


if False:  # get an overview
    extraction.plot_experiment()

t_zero = 50
if True:  # fig10a
    # I have nothing at 5 mV/s :(
    tspan_experiment = [000, 2700]

    extraction.timeshift(t_zero=t_zero)
    axes = extraction.plot_exchange(
        mol="O2",
        tspan=tspan_experiment,
    )
    extraction.plot_flux(
        extraction.mdict["H2"],
        tspan=tspan_experiment,
        ax=axes[-1],
        unit="pmol/s",
        logplot=False,
    )
    axes = extraction.plot_exchange(mol="CO2", tspan=tspan_experiment, axes=axes)

    for ax in [axes[0], axes[-1]]:
        ax.set_ylim([lim / 3 for lim in ax.get_ylim()])

    if forpublication:
        fig = axes[0].get_figure()
        axes[0].set_xlabel("time / [s]")
        axes[1].set_xlabel("time / [s]")
        fig.savefig("fig_Ir_hydrous_a.png")
        fig.savefig("fig_Ir_hydrous_a.svg")

    extraction.timeshift(t_zero=-t_zero)

if True:  # figure 10b
    extraction.timeshift(t_zero=t_zero)
    tspan_c1 = np.array([1815, 1945]) - 50 + t_zero
    tspan_c2 = np.array([2095, 2245]) - 50 + t_zero
    tspan_c3 = np.array([2395, 2545]) - 50 + t_zero
    tspan_bg = np.array([1500, 1600]) - 50 + t_zero

    excess_O2 = extraction.create_excess_mol(mol="O2")
    excess_CO2 = extraction.create_excess_mol(mol="CO2")

    for mass, val in excess_CO2.cal_mat.items():
        excess_CO2.cal_mat[mass] = val * 1 / (1 - 0.4)

    H2 = extraction.mdict["H2"]

    axes = extraction.plot_vs_potential(
        mols=[[H2], [excess_O2, excess_CO2]],
        t_bg=tspan_bg,
        tspan=tspan_c1,
        logplot=False,
        linestyle="-",
        unit="pmol/s",
        right_space=0.85,
    )
    axes = extraction.plot_vs_potential(
        mols=[[H2], [excess_O2, excess_CO2]],
        t_bg=tspan_bg,
        tspan=tspan_c2,
        logplot=False,
        linestyle="--",
        unit="pmol/s",
        ax=axes,
    )
    axes = extraction.plot_vs_potential(
        mols=[[H2], [excess_O2, excess_CO2]],
        t_bg=tspan_bg,
        tspan=tspan_c3,
        logplot=False,
        linestyle=":",
        unit="pmol/s",
        ax=axes,
    )

    if True:  # fix some stuff
        from EC_MS import colorax, align_zero

        align_zero(axes[-1], axes[0])
        colorax(axes[0], "b")

    if forpublication:
        fig = axes[0].get_figure()
        axes[-1].set_ylabel("lat. O / [pmol/s]")
        fig.savefig("fig_Ir_hydrous_b.png")
        fig.savefig("fig_Ir_hydrous_b.svg")

if True:  # grab some numbers

    t_zero = 50
    tspan_experiment = [000, 1100]

    tspan = [t + 50 for t in tspan_experiment]

    x_O2, y_O2 = extraction.get_flux(excess_O2, tspan=tspan, unit="pmol/s")
    x_CO2, y_CO2 = extraction.get_flux(excess_CO2, tspan=tspan, unit="pmol/s")

    fig, ax = plt.subplots()

    ax.plot(x_O2, y_O2, "r")
    ax.plot(x_CO2, y_CO2, "purple")

    n_O2 = np.trapz(y_O2, x_O2)
    n_CO2 = np.trapz(y_CO2, x_CO2)
