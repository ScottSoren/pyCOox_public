from ixdat import Measurement

meas = Measurement.read("../data/01_Pt_in_18O_electrolyte.pkl", reader="EC_MS")

tspan_CO_strip = [9700, 10500]
tspan_CO_ox = [11350, 12150]

for tspan in [tspan_CO_ox, tspan_CO_strip]:

    axes_a = meas.plot_measurement(tspan=tspan)
