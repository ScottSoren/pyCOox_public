from ixdat import Measurement

meas = Measurement.read("../data/09_Ir18Ox_yH2O.pkl", reader="EC_MS")

meas.plot_measurement()
