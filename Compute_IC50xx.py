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
Escape_data_site = pd.read_csv("Data/escape_data_site.csv") # downloaded from https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps/blob/main/processed_data/escape_data_site.csv
Class_Connect = pd.read_csv("Data/antibody_classes.csv")
Escape_Fraction = pd.read_csv("Data/escape_fraction_per_site.csv")


unique_sites = np.array(np.unique(Escape_data_site["site"]))
IC50_all = np.array(Escape_data_site["IC50"])

relevant_sites = np.array(np.unique(Escape_Fraction["site"]))

IC50xx_dic = {}
for i in range(len(relevant_sites)):
    site = relevant_sites[i]
    where_site = Escape_data_site["site"] == site
    groups = Escape_data_site["group"][where_site]
    conditions = Escape_data_site["condition"][where_site]
    unique_groups = np.unique(groups)
    unique_conditions = np.unique(conditions)
    if len(where_site)!=0:
        for group in unique_groups: 
            IC50xx_dic["(%s, %s)"%(site, group)] = []
            for cond in unique_conditions:
                if cond in np.array(Class_Connect["condition"][Class_Connect["group"] == group]):
                    res = IC50_all[where_site][(groups == group)&(conditions == cond)]
                    ic50_cond = np.mean(res[~np.isnan(res)])
                    IC50xx_dic["(%s, %s)"%(site, group)].append(ic50_cond)

"""
len(unique_sites) < len(relevant_sites), therefore we can't compute the IC50 per relevant sites. 
We compute the IC50xx by averaging the IC50xx corresponding to the available sites that corresponds to the relevant_AB
"""
IC50xx_relevant = {}   
relevant_AB = np.array(np.unique(Escape_Fraction["ab_class"]))
for j in range(len(relevant_AB)):
    group_j = Class_Connect["group"][Class_Connect["condition_subtype"] == relevant_AB[j]]        
    IC50_j = []
    for group in group_j:
        IC50_sites = []
        for site in relevant_sites:
            try:
                if len(IC50xx_dic["(%s, %s)"%(site, group)])!= 0:
                    ic50_group = np.mean(np.array(IC50xx_dic["(%s, %s)"%(site, group)]))
                    IC50_sites.append(ic50_group)
            except:
                pass
        
        IC50_j.append(np.mean(np.array(IC50_sites)))
    IC50xx_relevant[relevant_AB[j]] = [np.mean(np.array(IC50_j))]

ic50_df = pd.DataFrame(IC50xx_relevant)
ic50_df.to_csv("Data/IC50xx_per_Ab_class.csv")

"""
# Compute an IC50xx per relevant sites and antibody class
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
"""                  