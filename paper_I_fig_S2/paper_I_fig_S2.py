from ixdat import Measurement

meas = Measurement.read("../data/02_Pt_in_18O_electrolyte_35C.pkl", reader="EC_MS")


meas.tstamp += 12250
axes_a = meas.plot_measurement(tspan=[0, 7500])
