#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 14:16:37 2019

@author: scott
"""

import os, re, pickle
from matplotlib import pyplot as plt

from EC_MS import download_cinfdata_set, load_EC_set, synchronize
from EC_MS import plot_experiment, plot_signal, trigger_cal

plt.close("all")


# data_dir = os.path.expanduser('~/o/FYSIK/list-SurfCat/setups/sniffer/Data/pc Cu') # linux
data_dir = r"C:\EC_data\scott\CO_and_lattice_O_project"  # Windows
# U:\FYSIK\list-SurfCat\setups\sniffer\Data\single_crystal_Cu\19J08_Ruben
MSset = {"comment": "20A30%"}  #'time':'2020-01-26%',
folder = "20A30_16O"


if not os.path.isdir("pickles"):
    os.mkdir("pickles")  # to put the synchronized data, as pickle files
if not os.path.isdir("overviews"):
    os.mkdir("overviews")  # to put the overview plots


# -------- first, we get all the MS data from the folder and combine it -------- #
if True:  # it's a bit slow...

    MS_data = download_cinfdata_set(**MSset)

    # now, as a sanity check, and to find our way around the experiment, we plot all of the MS data!
    ax = plot_signal(MS_data)
    ax.legend()
    plt.savefig(
        "./overviews/" + folder + "_MS_data.png"
    )  # and save it as an overview plot.

    # also, save the combined MS data as a pickle!
    with open(
        "./pickles/" + folder + "_MS_data.pkl", "wb"
    ) as pkl:  # defines the file, 'wb' means 'write binary'.
        pickle.dump(MS_data, pkl)  # save MS_data into the file


# -------- Then, we load each EC experiment and combining it with MS data -------- #

# figure out which files correspond to EC data:
data_folder = data_dir + os.sep + folder
data_list = os.listdir(data_folder)  # list of all the files in the data directory

EC_files = [
    f
    for f in data_list
    if ".mpt" in f and "timescan" not in f  # re.search(r'\A[0-9]+', f) and
]
# ^ you have to save the EC files starting with a number for this to work!
EC_files.sort()  # get them in an order that makes sense

tags = set([f[0:2] for f in EC_files])

datasets = {}

for tag in tags:  # loop over the EC file names.
    # if tag in ['01', '02', '03']: # these are from 19J17.
    #    continue
    print(
        "\n\n----------- working on " + tag + " -------------\n"
    )  # so you can follow the Console output
    EC_data = load_EC_set(
        data_folder,
        tag=tag,
        # fix_CP=True
        # ^ use Ece/V to work around an EC-Lab bug which writes 0's for <Ewe>/V in .mpt's for CP
    )
    # ^ load the EC data (no longer need to use rename_RGA_cols, parse_RGA_header, or timeshift)
    if EC_data["empty"]:
        print("WARNING: tag " + tag + " is empty!")
        continue

    dataset = synchronize([EC_data, MS_data], cut_buffer=120, append=False)
    # ^ combine the EC data with the MS data! To save memory, we don't keep all of the MS data:
    # Instead, we just keep the part that overlaps with the EC data, plus 120 s on either end.

    try:
        trigger_cal(dataset)
    except IndexError:
        print("WARNING!!! " + tag + " couldn't trigger_cal.")

    name = folder + "_" + tag
    with open("./pickles/" + name + ".pkl", "wb") as pkl:  # and save it as a pickle
        pickle.dump(dataset, pkl)

    # make the plot:
    ax = plot_experiment(
        dataset,
        tspan="all",  # including this makes the plot include the MS data on the edges
    )
    ax[0].legend()
    ax[1].set_title(name)  # give it a title, so we know what's what
    # ax[0].legend(loc='lower left')
    plt.savefig("./overviews/" + name + ".png")  # and save it in overviews.
