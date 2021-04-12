import os
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from EC_MS.utils.extraction_class import Extraction
from EC_MS.converters import mdict_from_SI2020_calibration_file

plt.close("all")

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

calibration_dir = Path("../..").absolute()
calibration_file = calibration_dir / "20A25_sniffer_fixed.json"

mdict = mdict_from_SI2020_calibration_file(calibration_file)

extraction_dir = os.path.abspath(os.path.join(
    os.path.split(__file__)[0], "../extractions"
))
data_dir = os.path.abspath(os.path.join(
    os.path.split(__file__)[0], "../../pickles"
))

file = os.path.join(extraction_dir, "Decade1C_all.json")

extraction = Extraction.load(file, data_dir=data_dir, mdict=mdict)  # Ir^{18}O2 in H2^{16}O
if True:  # correct H2 signal and ratio
    extraction.mdict["H2"].cal_mat = {"M4": -8.5e-4 * 1/extraction.mdict["H2"].F_cal,
                                      "M2": 1/extraction.mdict["H2"].F_cal}
    # extraction.set_background(t_bg=[15325, 15375])
    # extraction.get_alpha(tspan=[14900, 15000], ax=None)

    tspan_bg = [1150, 1200]  # [150, 200]
    extraction.set_background(t_bg=tspan_bg)

    #extraction.get_alpha(tspan=[750, 850])
    extraction.alpha = 0.9980

if True:  # correct calibration
    tspan_steady = [750, 850]  # extraction.tspan_ratio
    O2 = extraction.point_calibration(
            mol="O2", mass="M32", tspan=tspan_steady, n_el=4
    )
    # since we use cal_mat but point_calibration gives F_cal, this is the comparison
    # to make. We should then adjust such that it makes sense.
    correction = O2.F_cal / extraction.mdict["O2_M32"].F_cal
    for mol, molecule in extraction.mdict.items():
        if mol in ["O2", "CO2"]:
            continue   # patch because stupid calibration is broken.
        if not hasattr(molecule, "cal_mat"):
            molecule.cal_mat = {molecule.primary: 1/molecule.F_cal}
        for mass, value in molecule.cal_mat.items():
            molecule.cal_mat[mass] = value / correction
        if hasattr(molecule, "F_cal"):
            molecule.F_cal = molecule.F_cal * correction
    
    # and, finally, for O2_M34, fix the background which increases in CO
    extraction.mdict["O2_M34"].cal_mat["M28"] = -8e-6 * extraction.mdict["O2_M34"].cal_mat["M34"]

t_zero = 200
extraction.timeshift(t_zero=t_zero)

if False:  # figure of data for SI

    if True:  # silence He after switch. It's wierd and distracting.
        extraction["M4-y"][extraction["M4-x"]>1900] = np.nan

    extraction.reset()
    axes = extraction.plot_experiment(
        masses=["M4", "M28", "M32", "M34", "M36", "M44", "M46", "M48", "M2"], 
        tspan=[000, 7500], 
        unit="pA",
        logplot=True,
    )
    # axes[0].legend()
    extraction.set_background(t_bg=[t - t_zero for t in extraction.t_bg])
    
    axes[0].set_ylabel("signal / [pA]")
    axes[1].set_xlabel("time / [s]")
    axes[0].set_xlabel("time / [s]")
    axes[0].set_yticks([1e-1, 1, 1e1, 1e2, 1e3, 1e4])
    
    if forpublication:
        fig = axes[0].get_figure()
        axes[0].set_xlabel("time / [s]")
        axes[1].set_xlabel("time / [s]")
        fig.set_figwidth(6.5)
        fig.savefig("fig_IrO2_extraction_a.png")
        fig.savefig("fig_IrO2_extraction_a.svg")

    fig = axes[0].get_figure()
    fig.savefig("fig_IrO2_extraction_raw.png")
    fig.savefig("fig_IrO2_extraction_raw.svg")


if True:  # main-text IrO2 fig a
    tspan_experiment = [000, 2965]

    axes = extraction.plot_exchange(mol="O2", tspan=tspan_experiment,
                                    unit_right="pmol/s/cm^2", unit_left="pmol/s/cm^2",
                                    )
    extraction.plot_flux(extraction.mdict["H2"], tspan=tspan_experiment,
                         ax=axes[-1], unit="pmol/s/cm^2", logplot=False)
    axes = extraction.plot_exchange(mol="CO2", tspan=tspan_experiment, axes=axes,
                                    unit_right="pmol/s/cm^2", unit_left="pmol/s/cm^2",)

    axes[1].set_xlabel("time / [s]")
    axes[0].set_ylabel("signal / [pmol s$^{-1}$cm$^{-2}$]")
    axes[-1].set_ylabel("signal / [pmol s$^{-1}$cm$^{-2}$]")


    if forpublication:
        fig = axes[0].get_figure()
        axes[0].set_xlabel("time / [s]")
        axes[1].set_xlabel("time / [s]")
        fig.savefig("r1_part_II_fig2a.png")
        fig.savefig("r1_part_II_fig2a.svg")

if True:  # main-text IrO2 fig b
    tspan_experiment = [3505, 5400]

    axes = extraction.plot_exchange(mol="O2", tspan=tspan_experiment,
                                    unit_right="pmol/s/cm^2", unit_left="pmol/s/cm^2",
                                    )
    extraction.plot_flux(extraction.mdict["H2"], tspan=tspan_experiment,
                         ax=axes[-1], unit="pmol/s/cm^2", logplot=False)
    axes = extraction.plot_exchange(mol="CO2", tspan=tspan_experiment, axes=axes,
                                    unit_right="pmol/s/cm^2", unit_left="pmol/s/cm^2",)

    axes[1].set_xlabel("time / [s]")
    axes[0].set_ylabel("signal / [pmol s$^{-1}$cm$^{-2}$]")
    axes[-1].set_ylabel("signal / [pmol s$^{-1}$cm$^{-2}$]")


    if forpublication:
        fig = axes[0].get_figure()
        axes[0].set_xlabel("time / [s]")
        axes[1].set_xlabel("time / [s]")
        fig.savefig("r1_part_II_fig2b.png")
        fig.savefig("r1_part_II_fig2b.svg")


if False:  # SI fig IrO2 a
    # I have nothing at 5 mV/s :(
    tspan_experiment = [000, 7500]

    axes = extraction.plot_exchange(mol="O2", tspan=tspan_experiment,
                                    unit_right="pmol/s/cm^2", unit_left="pmol/s/cm^2",
                                    right_space=0.9, left_space=0.1)
    extraction.plot_flux(extraction.mdict["H2"], tspan=tspan_experiment,
                         ax=axes[-1], unit="pmol/s/cm^2", logplot=False)
    axes = extraction.plot_exchange(mol="CO2", tspan=tspan_experiment, axes=axes,
                                    unit_right="pmol/s/cm^2", unit_left="pmol/s/cm^2",)

    axes[1].set_xlabel("time / [s]")
    axes[0].set_ylabel("signal / [pmol s$^{-1}$cm$^{-2}$]")
    axes[-1].set_ylabel("signal / [pmol s$^{-1}$cm$^{-2}$]")


    if forpublication:
        fig = axes[0].get_figure()
        axes[0].set_xlabel("time / [s]")
        axes[1].set_xlabel("time / [s]")
        fig.set_figwidth(6.5)
        fig.savefig("fig_IrO2_extraction_a.png")
        fig.savefig("fig_IrO2_extraction_a.svg")

extraction.timeshift(t_zero=-t_zero)

t_zero_CO_strip = 13300
if False:  # SI fig IrO2 b
    t_zero = t_zero_CO_strip
    extraction.timeshift(t_zero)
    tspan_experiment = [0, 1000]
    if True:  # assuming there's actuall been some O18 accumulation in the electrolyte
        tspan_bg = np.array([100, 200]) + 13000 - t_zero
        tspan_exchange = np.array([1900, 2000]) + 13000 - t_zero
        extraction.set_background(t_bg=tspan_bg)
        extraction.get_alpha(tspan=tspan_exchange, ax=None)
    axes = extraction.plot_exchange(mol="O2", tspan=tspan_experiment,
                                    unit_right="pmol/s/cm^2", unit_left="pmol/s/cm^2",
                                    )
    extraction.plot_flux(extraction.mdict["H2"], tspan=tspan_experiment,
                         ax=axes[-1], unit="pmol/s/cm^2", logplot=False)
    axes = extraction.plot_exchange(mol="CO2", tspan=tspan_experiment, axes=axes,
                                    unit_right="pmol/s/cm^2", unit_left="pmol/s/cm^2",)
    for ax in (axes[0], axes[-1]):
        ax.set_ylim([lim * 1.5 for lim in ax.get_ylim()])

    axes[1].set_xlabel("time / [s]")
    axes[0].set_ylabel("signal / [pmol s$^{-1}$cm$^{-2}$]")
    axes[-1].set_ylabel("signal / [pmol s$^{-1}$cm$^{-2}$]")

    if forpublication:
        fig = axes[0].get_figure()
        axes[0].set_xlabel("time / [s]")
        axes[1].set_xlabel("time / [s]")
        fig.savefig("fig_IrO2_extraction_b.png")
        fig.savefig("fig_IrO2_extraction_b.svg")

    extraction.timeshift(-t_zero)

if False:  # SI fig IrO2 c
    tspan_c1 = np.array([144, 244]) + t_zero_CO_strip
    tspan_c2 = np.array([444, 830]) + t_zero_CO_strip
    tspan_bg = np.array([100, 200]) + t_zero_CO_strip

    excess_O2 = extraction.create_excess_mol(mol="O2")
    excess_CO2 = extraction.create_excess_mol(mol="CO2")

    for mass, val in excess_CO2.cal_mat.items():
        excess_CO2.cal_mat[mass] = val * 1 / (1 - 0.4)

    H2 = extraction.mdict["H2"]

    axes = extraction.plot_vs_potential(
        mols=[[H2], [excess_O2, excess_CO2]], t_bg=tspan_bg,
        tspan=tspan_c1, logplot=False, linestyle="-", unit="pmol/s/cm^2",
        right_space=0.85
    )
    axes = extraction.plot_vs_potential(
        mols=[[H2], [excess_O2, excess_CO2]], t_bg=tspan_bg,
        tspan=tspan_c2, logplot=False, linestyle="--", unit="pmol/s/cm^2",
        ax=axes)
    axes[0].set_yticks([0, 250, 500, 750, 1000, 1250])
    axes[0].set_ylabel("   signal / [pmol s$^{-1}$cm$^{-2}$]")
    axes[-1].set_ylabel("   lat. O / [pmol s$^{-1}$cm$^{-2}$]")

    if True:  # fix some stuff
        from EC_MS import colorax, align_zero

        align_zero(axes[-1], axes[0])
        colorax(axes[0], 'b')

    if forpublication:
        fig = axes[0].get_figure()
        fig.savefig("fig_IrO2_extraction_c.png")
        fig.savefig("fig_IrO2_extraction_c.svg")


if True:  # calculate fluxes
    tspan_OER = [450, 1050]
    tspan_bg = [1150, 1200]

    excess_O2 = extraction.create_excess_mol(mol="O2")
    excess_CO2 = extraction.create_excess_mol(mol="CO2")
    for mass, val in excess_CO2.cal_mat.items():
        excess_CO2.cal_mat[mass] = val * 1 / (1 - 0.4)


    x_O2, y_O2 = extraction.get_flux(
        excess_O2, tspan=tspan_OER, t_bg=tspan_bg, unit="pmol/s"
    )
    x_CO2, y_CO2 = extraction.get_flux(
        excess_CO2, tspan=tspan_OER, t_bg=tspan_bg, unit="pmol/s"
    )

    if False:  # plot the excess
        fig, ax = plt.subplots()
        ax.plot(x_O2, y_O2, excess_O2.get_color())
        ax.plot(x_CO2, y_CO2, excess_CO2.get_color())

    n_ex_O2 = np.trapz(y_O2, x_O2)
    n_ex_CO2 = np.trapz(y_CO2, x_CO2)

    tspan_all_COox = [2000, 7500]
    CO2 = extraction.mdict["CO2_M44"]

    x_CO2_all_x, y_CO2_all_x = extraction.get_flux(
        excess_CO2, tspan=tspan_all_COox, t_bg=tspan_bg, unit="mol/s/cm^2"
    )
    x_CO2_all, y_CO2_all = extraction.get_flux(
        CO2, tspan=tspan_all_COox, t_bg=tspan_bg, unit="mol/s/cm^2"
    )
    n_CO2_all = np.trapz(y_CO2_all, x_CO2_all)
    n_ex_CO2_all = np.trapz(y_CO2_all_x, x_CO2_all_x)

    if False:
        ax.plot(x_CO2_all_x, y_CO2_all_x, color=excess_CO2.get_color())

    tspan_CO_strip = [13480, 13700]
    x, y = extraction.get_flux(CO2, tspan=tspan_CO_strip, unit="mol/s/cm^2")

    n_CO_strip = np.trapz(y, x)

    #if True:
    #    fig, ax = plt.subplots()
    #    ax.plot(x, y, color="brown")


    tspan_post_CO_OER = [13800, 14100]
    x_post_CO2, y_post_CO2 = extraction.get_flux(CO2, tspan=tspan_post_CO_OER, unit="mol/s/cm^2")
    x_post_O2, y_post_O2 = extraction.get_flux(O2, tspan=tspan_post_CO_OER, unit="mol/s/cm^2")
    x_post_x, y_post_x = extraction.get_flux(excess_CO2, tspan=tspan_post_CO_OER, unit="mol/s/cm^2")
    x_post_O2_x, y_post_O2_x = extraction.get_flux(excess_O2, tspan=tspan_post_CO_OER, unit="mol/s/cm^2")

    n_post_CO2 = np.trapz(y_post_CO2, x_post_CO2)
    n_post_O2 = np.trapz(y_post_O2, x_post_O2)
    n_post_CO2_x = np.trapz(y_post_x, x_post_x)
    n_post_O2_x = np.trapz(y_post_O2_x, x_post_O2_x)


    if True:
        fig, ax = plt.subplots()
        ax.plot(x_post_CO2, y_post_CO2, color="brown")
        ax.plot(x_post_O2, y_post_O2, color="black")
        ax.plot(x_post_x, y_post_x, color="purple")
        ax.plot(x_post_O2_x, y_post_O2_x, color="r")


    from EC_MS import Chem

    n_ML = 4/(np.sqrt(3)*(390e-12)**2) / Chem.NA * 1e-4


    print(f"CO2 in long extraction: {n_ex_CO2_all*1e9} nmol/cm^2 out of {n_CO2_all*1e9} nmol/cm^2")
    print(f"CO strip is {n_CO_strip*1e9} nmol/cm^2, compared to {n_ML*1e9} nmol/cm^2 surface Ir atoms")

    print(f"CO2 after strip: {n_post_CO2_x*1e12} pmol/cm^2 out of {n_post_CO2*1e9} nmol/cm^2")
    print(f"O2 after strip: {n_post_O2_x*1e12} pmol/cm^2 out of {n_post_O2*1e9} nmol/cm^2")



