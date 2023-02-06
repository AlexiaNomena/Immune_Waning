#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 10:38:34 2023

@author: raharinirina
"""
from .PK import Antibody
from .VE import sqrt_diff, vaccine_efficacy_n_antibodies
from scipy.optimize import root

import numpy as np
import pdb

""" Projection of Infection probability as a funciton of COV19 Variant proportion """ 

def Antibody_infos(t, antibody_data, VE_data):
    t_max = antibody_data["t_max"]
    t_half = antibody_data["t_half"]
    params_dic = {"t_max":t_max, "t_half":t_half}
    c_t, c_dframe, ka, ke = Antibody(t = t, params_dic = params_dic, is_log = False, save = False)
    IC50 = root(sqrt_diff, 0.2, args = (VE_data["days"], VE_data["Efficacy"], len(t_max), c_dframe), method = "lm").x
    return IC50, c_dframe
    

def IE_Per_Variant(t, infected, total_population, antibody_data, VE_data_wild, variant_data):
    """ Antibody efficacy for wild-type and antibody dynamics against wild-type """
    ic50_wild, c_dframe = Antibody_infos(t, antibody_data, VE_data_wild)
    print("ic50 of wild-type", ic50_wild)
    "Antibody efficiacy dynamics against variant"
    res_dic = {}
    variant_prop = variant_data["proportion"]
    fold_res = variant_data["fold resistance"]
    variant_name = variant_data["name"]
    day_activation = antibody_data["day_activation"]
    
    for j in range(variant_prop.shape[1]):
        IC50_var = ic50_wild*fold_res[j]
        expect_ve = np.zeros((len(t), len(t)))
        for i in range(len(t)):
            for k in range(i + day_activation, len(t)):
                antibody_level = c_dframe.loc[k - (i + day_activation)][1:]
                expect_ve[i, k] = (infected[i]*variant_prop[k, j])*vaccine_efficacy_n_antibodies(antibody_level, IC50_var)
             
            expect_ve[i, :] = expect_ve[i, :]/(total_population[i])
        res_dic[variant_name[j]] = np.sum(expect_ve, axis = 0)
        
    return res_dic
        