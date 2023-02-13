#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 10:44:40 2023

@author: raharinirina
"""

import pandas as pd
import numpy as np
import pdb

           

"""Extract antibody efficacy against wild type: IC50xx data for each relevant antibody and sites"""
Escape_data_site = pd.read_csv("Data/escape_data_site.csv")
Class_Connect = pd.read_csv("Data/antibody_classes.csv")
Escape_Fraction = pd.read_csv("Data/escape_fraction_per_site.csv")


unique_sites = np.array(np.unique(Escape_data_site["site"]))
IC50_all = np.array(Escape_data_site["IC50"])

relevant_AB = np.array(np.unique(Escape_Fraction["ab_class"]))
relevant_sites = np.array(np.unique(Escape_Fraction["site"]))

IC50xx_dic = {}
for i in range(len(relevant_sites)):
    site = unique_sites[i]
    if site in relevant_sites:
        where_site = Escape_data_site["site"] == site
        groups = Escape_data_site["group"][where_site]
        conditions = Escape_data_site["condition"][where_site]
        unique_groups = np.unique(groups)
        unique_conditions = np.unique(conditions)
        
        for group in unique_groups: 
            IC50xx_dic["(%s, %s)"%(site, group)] = []
            for cond in unique_conditions:
                if cond in np.array(Class_Connect["condition"][Class_Connect["group"] == group]):
                    ic50_cond = np.mean(IC50_all[where_site][(groups == group)&(conditions == cond)])
                    IC50xx_dic["(%s, %s)"%(site, group)].append(ic50_cond)


IC50xx_relevant = []   
for i in range(len(relevant_sites)):
    site = relevant_sites[i]
    for j in range(len(relevant_AB)):
        group_j = Class_Connect["group"][Class_Connect["condition_subtype"] == relevant_AB[j]]        
        IC50_j = []
        for group in group_j:
            try:
                ic50_group = np.mean(np.array(IC50xx_dic["(%s, %s)"%(site, group)]))
                IC50_j.append(ic50_group)
            except:
                pass
        IC50xx_relevant.append(np.mean(np.array(IC50_j)))

                    