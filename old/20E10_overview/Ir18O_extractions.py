# -*- coding: utf-8 -*-
"""
Created on Sun May  3 15:07:40 2020

@author: scott
"""
import os
from matplotlib import pyplot as plt

plt.close("all")

from EC_MS.utils.extraction_class import Extraction

data_dir = os.path.expanduser(
    "~/Dropbox (Spectro Inlets)/Soren_DTU/CO_and_lattice_O/Analysis/pickles"
)

extraction_specs = {
    "Jazz7_all": dict(
        data_file="20A30_16O_02.pkl",
        tspan_experiment=[0, 6000],
        tspan_exchange=[200, 550],
        tspan_extraction=[400, 1400],
        tspan_ratio=[9000, 9500],
        t_bg=[6500, 6600],
        electrolyte="16O",
        film="18O",
    ),
    "Decade1C_all": dict(
        data_files=["20A30_16O_03.pkl", "20A30_16O_04.pkl"],
        tspan_experiment=[0, 20000],
        tspan_exchange=[400, 1200],
        tspan_extraction=[2000, 6000],
        tspan_ratio=[14800, 15000],
        t_bg=[13300, 13400],
        electrolyte="16O",
        film="18O",
    ),
    "Jazz8_all": dict(
        data_file="20A30_16O_06.pkl",
        tspan_experiment=[0, 3500],
        tspan_exchange=[400, 1100],
        tspan_extraction=[1300, 2500],
        tspan_ratio=[950, 1000],
        t_bg=[3100, 3150],
        electrolyte="16O",
        film="18O",
    ),
    "Jazz9_all": dict(
        data_files=["20B02_16O_02.pkl"],  # "20B02_16O_01.pkl",
        tspan_experiment=[0, 12000],
        tspan_exchange=[2000, 5000],
        tspan_extraction=[5700, 6700],
        tspan_ratio=[11500, 11600],
        t_bg=[5500, 5600],
        electrolyte="16O",
        film="18O",
    ),
    "Jazz8b_all": dict(
        # b for "butterfly"
        data_file="20B02_16O_03.pkl",
        tspan_experiment=[0, 6000],
        tspan_exchange=[100, 1100],
        tspan_extraction=[1500, 4000],
        tspan_ratio=[5400, 5500],
        t_bg=[4200, 4400],
        electrolyte="16O",
        film="18O",
    ),
    "Jazz7_18O": dict(
        data_file="20A26_16O_02.pkl",
        tspan_experiment=[0, 8000],
        tspan_exchange=[100, 1100],
        tspan_extraction=[1500, 4000],
        tspan_ratio=[1100, 1200],
        t_bg=[5400, 5500],
        electrolyte="18O",
        film="16O",
    ),
    "Jazz6_18O": dict(
        # contaminated
        data_file="20A26_16O_03.pkl",
        tspan_experiment=[0, 5000],
        tspan_exchange=[900, 950],
        tspan_extraction=[1500, 4000],
        tspan_ratio=[900, 950],
        t_bg=[4800, 4900],
        electrolyte="18O",
        film="16O",
    ),
}
extractions = {}

# run_list = ["Decade1C_all", "Jazz8b_all"]
run_list = "all"
# run_list = ["Jazz7_all", "Jazz7_18O", "Jazz6_18O"]
calibration_file = "20A25_sniffer_fixed.json"
for name, spec in extraction_specs.items():
    if not run_list == "all":
        if not name in run_list:
            continue
    spec.update(data_dir=data_dir, calibration_file=calibration_file, name=name)
    print(f"---------- working on {name} ------------")
    extraction = Extraction(**spec)
    extraction.sync_metadata(RE_vs_RHE=0.715, A_el=0.196)
    # extraction.get_alpha(simple=False, ax=None)  # why was this here???

    axes = extraction.plot_exchange(mol="O2")
    if True:  # name in ["Decade1C_all", "Jazz8b_all"]:
        extraction.plot_flux(extraction.mdict["H2"], ax=axes[-1], unit="pmol/s")
        axes[-1].set_yscale("linear")
    extraction.plot_exchange(mol="CO2", axes=axes)
    axes[1].set_title(name)

    # plt.savefig(f"fig_Ir_extraction_{name}.png")
    # continue
    # extraction.plot_extraction_vs_potential(mol="CO2")
    # plt.savefig(f"fig_Pt_extraction_{name}_vs_U.png")
    extraction.quantify_extraction(mol="CO2")
    print(extraction.n_ex)

    json_name = name + ".json"
    extraction.name = name  # not sure why, but this gets reset by the dataset.
    extraction.save(json_name)
    extractions[name] = extraction
