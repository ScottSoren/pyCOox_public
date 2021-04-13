from ixdat import Measurement

meas = Measurement.read("../data/04_Pt18Ox.pkl", reader="EC_MS")

axes_a = meas.plot_measurement(tspan=[0, 2000])
