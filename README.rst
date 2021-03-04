pyCOox_public
=============
Data, analysis, and plotting scripts for articles
-------------------------------------------------

This repository is made for the curious reader of the two articles

- Soren B. Scott, Jakob Kinsgaard, Peter C. K. Vesborg, and Ib Chorkendorff.  Tracking oxygen atoms in electrochemical CO oxidation – Part I: Oxygen exchange via CO2 hydration. `Electrochimica Acta, 2021 <https://doi.org/10.1016/j.electacta.2021.137842>`_.
- Soren B. Scott, Jakob Kinsgaard, Peter C. K. Vesborg, and Ib Chorkendorff.  Tracking oxygen atoms in electrochemical CO oxidation - Part II: Lattice oxygen reactivity in oxides of Pt and Ir. `Electrochimica Acta, 2021 <https://doi.org/10.1016/j.electacta.2021.137844>`_.

It has all of the data analysis and plotting scripts. These scripts can serve as examples of electrochemistry - mass spectrometry (EC-MS) data analysis, especially for isotope-labeling experiments.

Right now, it is a mess, apologies! Over the next days, it will be reorganized so that it is clear which script makes which publication from the article.

The raw data used by these scripts is available `here. <https://www.dropbox.com/sh/owxna2hsocaw7vo/AADdQCNhZhvQ0uAD-xdCrno-a?dl=0>`_. This will also be reorganized in the next days.

Right now, this repository requires the legacy  `EC_MS <https://github.com/ScottSoren/EC_MS>`_ package.

In the coming weeks, it will be redone completely to instead work with The In-situ Experimental Tata Tool, `ixdat <https://github.com/ixdat/ixdat>`_, which is replacing ``EC_MS``.
