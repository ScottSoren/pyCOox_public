from ixdat import Measurement

meas_Pt = Measurement.read("../data/01_Pt_in_18O_electrolyte.pkl", reader="EC_MS")
axes_Pt = meas_Pt.plot_measurement()

meas_Ir = Measurement.read("../data/03_Ir_in_18O_electrolyte.pkl", reader="EC_MS")
axes_Ir = meas_Ir.plot_measurement()
