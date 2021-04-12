import os
import pickle
from EC_MS import Dataset

data_dir = os.path.abspath(os.path.join(os.path.split(__file__)[0], "../../pickles"))

dataset = Dataset(os.path.join(data_dir, "20B01_16O_01.pkl"))

example_set = dataset.cut(tspan=[8700, 9300], t_zero="start")

example_set["time/s"] = example_set["time/s*"]
example_set["t_str"] = "time/s"

if True:
    essential_cols = [
        "I/mA",
        '(Q-Qo)/C',
        # '<Ece>/V',
        # '<Ewe>/V',
        # '<I>/mA',
        'Ewe-Ece/V',
        'Ewe/V',
        'I/mA', 'M18-x',
        'M18-y',
        'M2-x',
        'M2-y',
        'M28-x',
        'M28-y',
        'M32-x',
        'M32-y',
        'M4-x',
        'M4-y',
        'M44-x',
        'M44-y',
        'control/V',
        'control/mA',
        'counter inc.',
        'cycle number',
        'dQ/C',
        'error',
        # 'file number_0',
        'time/s',
        'ox/red',
    ]

    for col in example_set.data_cols:
        if col not in essential_cols:
            del (example_set.data[col])

example_set.plot_experiment()

with open("example_set.pkl", "wb") as f:
    pickle.dump(example_set.data, f)
