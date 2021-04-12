import os
import numpy as np
from matplotlib import pyplot as plt
from EC_MS.utils.extraction_class import Extraction

forpublication = False
if forpublication:  # for the publication figure
    import matplotlib as mpl
    mpl.rcParams['figure.figsize'] = (3.25, 2.75)
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')
    plt.rc('font', size=8)
    plt.rc('lines', linewidth=0.5)
else:
    plt.style.use("default")

extraction_dir = os.path.abspath(os.path.join(
    os.path.split(__file__)[0], "../extractions"
))
data_dir = os.path.abspath(os.path.join(
    os.path.split(__file__)[0], "../../pickles"
))

file_A = os.path.join(extraction_dir, "Jazz8_all.json")
file_B = os.path.join(extraction_dir, "Jazz7_all.json")

extraction_A = Extraction.load(file_A, data_dir=data_dir)  # Ir^{18}O in H2^{16}O
extraction_B = Extraction.load(file_B, data_dir=data_dir)  # Ir^{18}O in H2^{16}O
for extraction in [extraction_A, extraction_B]:
    if True:  # correct H2 signal
        extraction.mdict["H2"].cal_mat["M4"] = -8.5e-4 * extraction.mdict["H2"].cal_mat["M2"]
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
    extraction_A.plot_experiment()
    extraction_B.plot_experiment()

if True:  # fig8a
    # I have nothing at 5 mV/s :(
    tspan_experiment = [000, 1800]

    extraction_A.timeshift(300)
    axes = extraction_A.plot_exchange(mol="O2", tspan=tspan_experiment)
    extraction_A.plot_flux(extraction_A.mdict["H2"], tspan=tspan_experiment,
                           ax=axes[-1], unit="pmol/s", logplot=False)
    axes = extraction_A.plot_exchange(mol="CO2", tspan=tspan_experiment, axes=axes)

    if forpublication:
        axes[0].set_xlabel("time / [s]")
        axes[1].set_xlabel("time / [s]")
        axes[0].get_figure().savefig("fig_Ir_extraction_a.png")
        axes[0].get_figure().savefig("fig_Ir_extraction_a.svg")


if True:  # fig8b, take 1
    tspan_experiment = [000, 450]

    extraction_B.timeshift(5450)
    axes = extraction_B.plot_exchange(mol="O2", tspan=tspan_experiment)
    extraction_B.plot_flux(extraction_B.mdict["H2"], tspan=tspan_experiment,
                           ax=axes[-1], unit="pmol/s", logplot=False)
    axes = extraction_B.plot_exchange(mol="CO2", tspan=tspan_experiment, axes=axes)


    if forpublication:
        axes[0].set_xlabel("time / [s]")
        axes[1].set_xlabel("time / [s]")
        axes[0].get_figure().savefig("fig_Ir_extraction_OER_in_CO.png")
        axes[0].get_figure().savefig("fig_Ir_extraction_OER_in_CO.svg")

extraction_A.name = "start"
extraction_B.name = "later"

if True:  # fig8b, take 2
    tspan_A = [124, 600]
    tspan_B = [55, 300]

    axes = "new"
    results = {}
    for extraction, tspan_0, linestyle in [
        (extraction_A, tspan_A, "-"), (extraction_B, tspan_B, "--")
    ]:
        t_zero = tspan_0[0]
        tspan = [t-t_zero for t in tspan_0]
        extraction.timeshift(t_zero)
        O2_excess = extraction.create_excess_mol(mol="O2")
        CO2_excess = extraction.create_excess_mol(mol="CO2")
        for mass, val in CO2_excess.cal_mat.items():
            CO2_excess.cal_mat[mass] = val * 1/(1-0.4)
        axes = extraction.plot_experiment(
            ax=axes,
            mols=[O2_excess, CO2_excess],
            logplot=False,
            tspan=tspan,
            unit="pmol/s",
            plotcurrent=False,
            spec=dict(linestyle=linestyle),
            verbose=False,
        )
        extraction.timeshift(-t_zero)
        for ax in axes:
            ax.set_xlim([t-20 for t in tspan])

        x_O2, y_O2 = extraction.get_flux(O2_excess, tspan=tspan)
        x_CO2, y_CO2 = extraction.get_flux(CO2_excess, tspan=tspan)

        fig, ax = plt.subplots()
        ax.plot(x_O2, y_O2, "r")
        ax.plot(x_CO2, y_CO2, "purple")

        n_O2 = np.trapz(y_O2, x_O2)
        n_CO2 = np.trapz(y_CO2, x_CO2)
        results[extraction.name] = {"n_O2": n_O2, "n_CO2": n_CO2}

        print(f"n_O2 = {n_O2}, n_CO2 = {n_CO2}")

    if forpublication:
        axes[0].set_xlabel("time / [s]")
        axes[1].set_xlabel("time / [s]")
        axes[0].get_figure().savefig("fig_Ir_extraction_OER_in_both.png")
        axes[0].get_figure().savefig("fig_Ir_extraction_OER_in_both.svg")
