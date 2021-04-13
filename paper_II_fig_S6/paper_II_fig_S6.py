from ixdat import Measurement

meas = Measurement.read("../data/06_Ir18O2.pkl", reader="EC_MS")

meas.plot_measurement(tspan=[000, 20000])
