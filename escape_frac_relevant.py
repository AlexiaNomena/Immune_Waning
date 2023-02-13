#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 10:44:40 2023

@author: raharinirina
"""

""" Extracting relevant escape fraction data """
import pandas as pd
import numpy as np
import pdb

Escape = pd.read_csv("Data/escape_data.csv") ### Too large to put in the repository (Dowloaded from https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps/blob/main/processed_data/escape_data.csv)

Condition = np.unique(Escape["condition"])

df = {}
df["condition"] = []
df["condition_type"] = []
df["condition_subtype"] = []
df["know_to_neutralize"] = []
df["mut. RBD sites"] = []
df["site_total_escape"] = []

for i in range(len(Condition)):
    test0 = Escape["condition"] == Condition[i]
    Type = Escape["condition_type"][test0]
    unique_type = np.unique(Type)
    for typ in unique_type:
        test1 = test0 & (Escape["condition_type"] == typ)
        subCond = Escape["condition_subtype"][test1]
        unique_subCond = np.unique(subCond)
        for sclass in unique_subCond:    
            test2 = test1 & (Escape["condition_subtype"] == sclass)
            virus_Aclass = Escape["known_to_neutralize"][test2]
            unique_vir = np.unique(virus_Aclass)
            for vir in unique_vir:
                test3 = test2 & (Escape["known_to_neutralize"] == vir)
                sites_Aclass = Escape["protein_site"][test3]
                unique_sites = np.unique(sites_Aclass)
                for prot_s in unique_sites:
                    test4 = test3 & (Escape["protein_site"] == prot_s)
                    escape_Aclass = Escape["site_total_escape"][test4]
                    escape_Aclass[escape_Aclass > 1] = 1  ### Check with the others if this data should be removed instead
                    unique_esc = np.unique(escape_Aclass)
                    for esc in unique_esc:
                        df["condition"].append(Condition[i])
                        df["condition_type"].append(typ)
                        df["condition_subtype"].append(sclass)
                        df["know_to_neutralize"].append(vir)
                        df["mut. RBD sites"].append(prot_s)
                        df["site_total_escape"].append(esc)

df = pd.DataFrame(data = df)
df.to_csv("Data/escape_fraction.csv")            