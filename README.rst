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
(Paper I) shows the derivation of the mass spec calibration (sensitivity factors), RHE calibration, and working distance calibration used throughout the article.
See ``old/20J13_SI/calibration.py``

paper_I_fig_2
--------------
Figure 2 of Paper I

See ``old/20A31/fig_COox_Pt_20A31.py``


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
(Paper II)

See ``old/20A31/fig_COox_Ir_20A31.py``

paper_II_fig_S2
---------------
Supplementary Figure S2 of Paper II

See ``old/20A31/fix_oxidation_Pt.py`` and ``old/20A31/fix_oxidation_Ir.py``

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
https://ixdat.readthedocs.io/en/latest/introduction.html#