pyCOox_public
#############

Data, analysis, and plotting scripts for articles
=================================================

This repository is made for the curious reader of the two articles

- Soren B. Scott, Jakob Kinsgaard, Peter C. K. Vesborg, and Ib Chorkendorff.  Tracking oxygen atoms in electrochemical CO oxidation – Part I: Oxygen exchange via CO2 hydration. `Electrochimica Acta, 2021 <https://doi.org/10.1016/j.electacta.2021.137842>`_.

- Soren B. Scott, Jakob Kinsgaard, Peter C. K. Vesborg, and Ib Chorkendorff.  Tracking oxygen atoms in electrochemical CO oxidation - Part II: Lattice oxygen reactivity in oxides of Pt and Ir. `Electrochimica Acta, 2021 <https://doi.org/10.1016/j.electacta.2021.137844>`_.

It has all of the data analysis and plotting scripts. These scripts can serve as examples of electrochemistry - mass spectrometry (EC-MS) data analysis, especially for isotope-labeling experiments.

Setup
=====
To run the scripts:

1. Make sure you have python 3.6+ installed with the ``numpy``, ``scipy``, and ``matplotlib`` packages as a minimum. I recommend Anaconda python.

2. Install the latest version of `ixdat <https://ixdat.readthedocs.org>`_ by typing in your terminal or Anaconda prompt::

  $ pip install --upgrade ixdat

3. Clone or download this repository using git or github.

4. Copy all of the .pkl files from `here <https://www.dropbox.com/sh/trro30vogoy5k5p/AAAWf-Rs9bSNzcmPNnkzqyLJa?dl=0>`_ into the folder called data.

5. Run the scripts with your favorite pythin IDE. I recommend spyder for quick plotting/analysis or PyCharm for development.
A few are also available as .ipynb for use as tutorials with Jupyter Notebooks

The Scripts
===========

The following folders contain a script and some files it produces.
They are listed in an order logical for data flow.
Scripts may depend on files produced by other scripts above them in this list (if so this will be clear in the comments of the script.)

paper_I_fig_S1
--------------
Supplementary Figure S1 of `Tracking oxygen atoms in electrochemical CO oxidation – Part I: Oxygen exchange via CO2 hydration <https://doi.org/10.1016/j.electacta.2021.137842>`_
(Paper I) shows the derivation of the mass spec calibration (sensitivity factors), RHE calibration, and working distance calibration used throughout the articles.

The script paper_I_fig_S1.py works produces the three subfigures and functions as a tutorial on chip EC-MS calibration.
The folder also contains the calibration file (produced by the script) that other scripts in this repository read.

paper_I_fig_2
--------------
Figure 2 of Paper I is a repeat of the classic Figure 3 of `Trimarco et al, 2018 <https://doi.org/10.1016/j.electacta.2018.02.060>`_
but in 18-O labeled electrolyte. Two experiments are ploted: cyclic voltammatry in He-saturated and CO-saturated electrolyte (a and b),
and a CO stripping experiment (c and d).

The script makes the EC-MS plots of calibrated molecular fluxes together with electrohemical data both vs time (a and c) and vs potential (b and d).
It shows that manipulation of these EC-MS data sets becomes easy and (dare I admit it) fun with ixdat.

paper_I_fig_3
--------------
Figure 3 of Paper I

See ``old/20G07_hydration_model/Pt_25C_CO_strip.py``

paper_I_fig_4
--------------
Figure 4 of Paper I

See ``old/20G07_hydration_model/Pt_25C_CO_strip.py``

paper_I_fig_S2
--------------
Supplementary Figure S2 of Paper I

See ``old/20G07_hydration_model/Pt_35C_CO_burst.py``

paper_II_fig_S1
---------------
Figure S1 of
`Tracking oxygen atoms in electrochemical CO oxidation - Part II: Lattice oxygen reactivity in oxides of Pt and Ir. <https://doi.org/10.1016/j.electacta.2021.137844>`_.
(Paper II) is for iridium what Figure 2 of Paper I (itself a repeat of the classic Figure 3 of `Trimarco et al, 2018 <https://doi.org/10.1016/j.electacta.2018.02.060>`_) is for platinum:
Cyclic voltammatry in He-saturated and CO-saturated 1.0 M HClO4 electrolyte, and a CO stripping experiment.

The script makes the EC-MS plots of calibrated molecular fluxes together with electrohemical data both vs time and vs potential.
It shows that manipulation of these EC-MS data sets becomes easy and (dare I admit it) fun with ixdat.

Interesting, Ir behaves very similar to Pt with respect to hydrogen and carbon (HER and CO oxidation), but as the rest of the paper makes clear,
very differently when it comes to oxygen (OER).

paper_II_fig_S2
---------------
Supplementary Figure S2 of Paper II shows the formation of PtOx and IrOx during constant
anodic current on Pt and Ir, respectively, and their subsequent reduction to determine
the thickness of the oxide layer.

It shows how to calculated the difference in charge passed in corresponding parts of two
cyclic voltammograms using ixdat.

paper_II_fig_1
--------------
Figure 1 of Paper II

See ``old/20E16_Pt/fig_Pt_extraction.py``

paper_II_fig_S3
---------------
Supplementary Figure S3 of Paper II

See ``old/20E16_Pt/fig_Pt_extraction.py``

paper_II_fig_S4
---------------
Supplementary Figure S4 of Paper II

See ``old/20E16_Pt/fig_Pt_extraction.py``

paper_II_fig_S5
---------------
Supplementary Figure S5 of Paper II

See ``old/20E23_Ir/fig_Ir_extraction_sputtered_IrO2.py``

paper_II_fig_S6
---------------
Supplementary Figure S6 of Paper II

See ``old/20E23_Ir/fig_Ir_extraction_sputtered_IrO2.py``

paper_II_fig_3
--------------
Figure 3 of Paper II

See ``old/20E23_Ir/fig_Ir_extraction_sputtered_IrO2.py``

paper_II_fig_S7
---------------
Supplementary Figure S7 of Paper II

See ``old/20E23_Ir/fig_Ir_extraction_1.py``

paper_II_fig_S8
---------------
Supplementary Figure S8 of Paper II

See ``old/20E23_Ir/fig_Ir_extraction_butterfly_Ir.py``

paper_I_fig_5
-------------
Supplementary Figure 5 of Paper I

See ``old/20G24_comparison/comparison_bar_plot.py``


We're still working on it!
==========================

At present, not all of the scripts are reworked for use with ``ixdat``, and instead still require the legacy  `EC_MS <https://github.com/ScottSoren/EC_MS>`_ package.
Analysis and plotting which has not been converted is in the folder **old**, in the location indicated above.
The scripts in this folder are unfortunately not very well organized and readable. Please contact me if you need the script working for one of the figures before it is ready.

End
===
Enjoy, and if you find this useful, help us make ixdat even more useful for everyone:
https://ixdat.readthedocs.io/en/latest/introduction.html