Copy raw data files to this folder!

A link giving access to the data files is in the repository README.

The raw data files are python pickle files (.pkl). When loaded, each is a large semi-homogeneous
dictionary of metadata and data. This is the primary data structure in the old ``EC_MS`` python package,
which was originally used to import the data from .mpt text files (Biologic's SP-150 potentiostat and EC-Lab software) and from
the `cinfdata <https://cinfdata-dababase-client.readthedocs.io/en/latest/index.html>`_ server
(Pfeiffer QMA 125 mass spectrometer and `PyExpLabSys <https://github.com/CINF/PyExpLabSys>`_).

``ixdat``'s "EC_MS" reader converts these dictionaries into an ``ECMSMeasurement`` in ixdat's much more
powerful relational data structure.

All of the scripts in this repository use a relative path pointing to the data in this folder. So once
the data is here, they should run without requiring any modification at all.

Enjoy, and if you find this useful, help us make ixdat even more useful for you and others:
https://ixdat.readthedocs.io/en/latest/introduction.html#