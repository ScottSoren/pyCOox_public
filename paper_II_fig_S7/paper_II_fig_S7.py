from ixdat import Measurement

meas_1 = Measurement.read("../data/07_Ir18Ox_1.pkl", reader="EC_MS")

axes_1 = meas_1.plot_measurement()


meas_2 = Measurement.read("../data/08_Ir18Ox_2.pkl", reader="EC_MS")

axes_2 = meas_2.plot_measurement()

