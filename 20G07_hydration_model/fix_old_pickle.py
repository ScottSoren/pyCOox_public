# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 16:22:08 2020

@author: scott
"""
from pathlib import Path
import pickle

OLD_PICKLE_FOLDER = Path("../../19J23_Thesis_17K13_GroupMeeting").absolute()
NEW_PICKLE_FOLDER = Path("../../pickles").absolute()


with open(OLD_PICKLE_FOLDER / "datasets.pkl", "rb") as f:
    datas = pickle.load(f)

for key, data in datas.items():
    pkl_name = f"17J05_Pt_in_18O_{key}.pkl"
    with open(NEW_PICKLE_FOLDER / pkl_name, "wb") as f:
        pickle.dump(data, f)