# -*- coding: utf-8 -*-
"""
Created on Sun May  3 15:07:40 2020

@author: scott
"""
import os
from matplotlib import pyplot as plt

# plt.close("all")

from EC_MS.utils.extraction_class import Extraction

calibration_file = "20A25_sniffer.json"
extraction_dir = os.path.join(os.path.split(__file__)[0], "../extractions")

if False:  # fix calibration
    from siCalibration.calibration_class import Calibration

    calibration = Calibration.load(calibration_file)
    calibration.mol_list = [
        "H2", "He", "CO", "O2_M32", "O2_M34", "O2_M36", "CO2_M44", "CO2_M46", "CO2_M48"
    ]
    calibration.mass_list = [
        "M2", "M4", "M28", "M32", "M34", "M36", "M44", "M46", "M48"
    ]
    calibration.F["CO"]["M44"] = calibration.F["CO"]["M28"] * 1.2e-4
    calibration.F["CO"]["M34"] = calibration.F["CO"]["M28"] * 5.0e-6
    calibration_file = calibration_file.split(".")[0] + "_fixed.json"
    calibration.save(calibration_file)
else:
    calibration_file = calibration_file.split(".")[0] + "_fixed.json"

data_dir = os.path.expanduser(
    "~/Dropbox (Spectro Inlets)/Soren_DTU/CO_and_lattice_O/Analysis/pickles"
)

extraction_specs = {
    "scanned": dict(
        data_file="20A26_16O_01_fixed.pkl",
        tspan_experiment=[0, 2000],
        tspan_exchange=[200, 550],
        tspan_extraction=[1150, 1400],
        tspan_ratio=[5000, 5500],
        t_bg=[675, 700],
        electrolyte="16O",
        film="18O",
    ),
    "paused": dict(
        data_file="20A30_16O_05.pkl",
        tspan_experiment=[0, 7000],
        # tspan_exchange = [200, 550],
        tspan_extraction=[1700, 2100],
        tspan_ratio=[6000, 6300],
        t_bg=[1075, 1100],
        electrolyte="16O",
        film="18O",
    ),
    "drift": dict(
        data_file="20B01_16O_01.pkl",
        tspan_experiment=[0, 8000],
        # tspan_exchange = [200, 550],
        tspan_extraction=[3650, 3900],
        tspan_ratio=[5500, 6000],
        t_bg=[3500, 3600],
        electrolyte="16O",
        film="18O",
    ),
    "scanned_B": dict(
        data_files=["20A25_18O_01.pkl",
                    "20A25_18O_02.pkl"],
        tspan_experiment=[200, 2200],
        calibration_file=calibration_file,
        # tspan_exchange = [200, 550],
        tspan_extraction=[1400, 1600],
        tspan_ratio=[650, 700],
        t_bg=[1020, 1070],
        electrolyte="18O",
        film="16O",
    ),
    "contaminated_B": dict(
        data_file="20A26_16O_05.pkl",
        tspan_experiment=[0, 3000],
        # tspan_exchange = [200, 550],
        tspan_extraction=[1700, 2400],
        tspan_ratio=[800, 900],
        t_bg=[1100, 1150],
        electrolyte="18O",
        film="16O",
    ),
    "paused_B": dict(
        data_file="20A31_18O_01.pkl",
        tspan_experiment=[0, 2500],
        calibration_file=calibration_file,
        # tspan_exchange = [200, 550],
        tspan_extraction=[1600, 2050],
        tspan_ratio=[800, 1100],  # the ratio gets worse later
        t_bg=[3700, 3800],
        electrolyte="18O",
        film="16O",
    ),
}
extractions = {}

run_list = [
    "scanned", "paused", "drift",
#    "scanned_B", "contaminated_B", "paused_B"
]
run_list = "all"
for name in extraction_specs.keys():
    spec = extraction_specs[name]
    if not run_list == "all":
        if name not in run_list:
            continue
    spec.update(data_dir=data_dir, calibration_file=calibration_file)
    print(f"---------- working on {name} ------------")
    extraction = Extraction(**spec)

    extraction.plot_experiment(tspan="all")
    continue

    #continue # so that I can test just the Extraction initiation
    extraction.sync_metadata(RE_vs_RHE=0.715, A_el=0.196)
    #extraction.get_alpha(simple=False, ax=None)

    axes = extraction.plot_exchange(mol="O2")
    extraction.plot_exchange(mol="CO2", axes=axes)

    plt.savefig(f"fig_Pt_extraction_{name}.png")

    # extraction.plot_extraction_vs_potential(mol="CO2")
    # plt.savefig(f"fig_Pt_extraction_{name}_vs_U.png")

    extraction.quantify_extraction(mol="CO2")
    print(extraction.n_ex)

    json_name = os.path.join(extraction_dir, name + ".json")
    extraction.save(json_name)
    extractions[name] = extraction

print(extraction.name)
print(extraction.calibration_file)
print(extraction.calibration.F["CO"])
# empty last line so I can better select all with vim.

