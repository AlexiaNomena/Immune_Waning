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

""" Projection of Neutralization probability as a funciton of COV19 Variant proportion """ 

def Antibody_infos(t, antibody_data, VE_data):
    t_max = antibody_data["t_max"]
    t_half = antibody_data["t_half"]
    params_dic = {"t_max":t_max, "t_half":t_half}
    c_t, c_dframe, ka, ke = Antibody(t = t, params_dic = params_dic, is_log = False, save = False)
    IC50 = root(sqrt_diff, 0.2, args = (VE_data["days"], VE_data["Neut"], 1, c_dframe), method = "lm").x
    return IC50, c_dframe
    

def Neut_Total(t, infected, antibody_data, VE_data_wild, variant_data):
    """ Neutralization of wild-type and antibody dynamics against wild-type """
    ic50_wild, c_dframe = Antibody_infos(t, antibody_data, VE_data_wild)
    
    " Neutralization dynamics against variant"
    res_dic = {}
    variant_prop = variant_data["proportion"]
    fold_res = variant_data["fold resistance"]
    variant_name = variant_data["name"]
    for j in range(variant_prop.shape[1]):
        IC50_var = ic50_wild*fold_res[j]
        expect_ve = np.zeros(len(t))
        for i in range(len(t)):
            antibody_level = c_dframe.loc[i][1:]
            expect_ve[i] = (infected[i]*variant_prop[i, j])*vaccine_efficacy_n_antibodies(antibody_level, IC50_var)
        res_dic[variant_name[j]] = expect_ve
    
    return res_dic
        