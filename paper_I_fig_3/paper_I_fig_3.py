from ixdat import Measurement

meas = Measurement.read("../data/01_Pt_in_18O_electrolyte.pkl", reader="EC_MS")

axes_a = meas.plot_measurement(tspan=[9700, 10500])
