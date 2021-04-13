from ixdat import Measurement

meas = Measurement.read("../data/03_Ir_in_18O_electrolyte.pkl", reader="EC_MS")

tspan_CO_strip = [4450, 5200]
tspan_CO_ox = [5625, 6350]

for tspan in [tspan_CO_ox, tspan_CO_strip]:

    axes_a = meas.plot_measurement(tspan=tspan)