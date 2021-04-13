from ixdat import Measurement

meas_A = Measurement.read("../data/04_Pt18Ox.pkl", reader="EC_MS")
axes_A = meas_A.plot_measurement(tspan=[0, 2000])

meas_B = Measurement.read("../data/05_Pt16Ox_in_18O_electrolyte.pkl", reader="EC_MS")
axes_B = meas_B.plot_measurement(tspan=[200, 2200])
