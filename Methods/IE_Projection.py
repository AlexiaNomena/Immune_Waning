#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 10:38:34 2023

@author: raharinirina
"""
from .PK import Antibody
from .VE import sqrt_diff, vaccine_efficacy_n_antibodies, efficacy_n_antibodies
from scipy.optimize import root

import numpy as np
import pdb

""" Compute cross reactivity map """
mut_sites_library = {"Wuhan-Hu-1":[],
                        "Alpha": [501],
                         "Beta": [417, 484, 501],
                         "Gamma": [417, 484, 501],
                         "Delta": [452, 478],
                         "Lambda": [452, 490],
                         "Mu": [346, 484, 501],
                         "Omicron":[339, 371, 373, 375, 417, 440, 446, 477, 478, 484, 493, 496, 498, 501, 505],
                         "PMS20":[346, 417, 440, 445, 455, 475, 484, 501]}


def FR_xy(variant_1, variant_2, escape_per_sites, ab, mut_sites_per_variant = mut_sites_library, EF_func = ("MAX","MEAN"), GM = False):
    sites_1 = set(mut_sites_per_variant[variant_1])
    sites_2 = set(mut_sites_per_variant[variant_2])
    
    sites = list(sites_1.symmetric_difference(sites_2))
    if EF_func[1] == "MEAN":
        escape_data = escape_per_sites["mean_escape_fraction_per_site"].values
    elif EF_func[1] == "MAX":
        escape_data = escape_per_sites["max_escape_fraction_per_site"].values

    
    #### Bound Escape data####
    escape_data[escape_data >= 1.] = 0.99
    escape_data[escape_data <= 1./101] = 1./101
    
    Missed = []
    Greater_one = []
    prod_FR = 1
    for s in sites:
        S = (escape_per_sites["site"].values).astype(str) == str(s)
        AB = (escape_per_sites["group"].values == ab)
        where_s_ab = np.where(S & AB)[0]
        
        if len(where_s_ab != 0):
            if EF_func[0] == "MAX":
                EF = np.max(escape_data[where_s_ab].copy())
                
            elif EF_func[0] == "MEAN":
                EF = np.mean(escape_data[where_s_ab].copy())
            
            prod_FR *= 100*EF/(1.- EF)
            if EF/(1.- EF)<0:
                pdb.set_trace()
                
            if np.any(escape_data[where_s_ab]>1):
                cond_list = (escape_per_sites["condition"].values)[(where_s_ab) & (escape_data > 1)]
                for cond in cond_list:
                    Greater_one.append("Escape fraction = %.2f at site %s for AB %s and condition %s‚"%(escape_data[where_s_ab], s, ab, cond))
        else:
            Missed.append("No Escape fraction: at site %s and AB %s"%(s, ab))
    if GM:
        if len(sites) != 0:
            prod_FR = prod_FR**(1/len(sites))
    
    Greater_one.append("\nNum of EF > 1 = %d out of %d, Num of EF<1/101= %d out of %d"%(np.sum(escape_data > 1.), len(escape_data), np.sum(escape_data<=1./101), len(escape_data)))
    return prod_FR, Missed, Greater_one # ab needs to bind just one of the sites


def cross_reactivity(variant_name, escape_per_sites, Ab_classes, mut_sites_per_variant = mut_sites_library, EF_func = ("MAX","MEAN"), GM = False):
    FRxy = {}
    Missed = []
    Greater_one = []
    for ab in Ab_classes:
        FRxy_ab = np.ones((len(variant_name), len(variant_name)))
        for i in range(len(variant_name)):
            for j in range(len(variant_name)):
                if (i !=j) & (j>i):
                    tot_FR, missed, gOne = FR_xy(variant_name[i], variant_name[j], escape_per_sites, ab, mut_sites_per_variant, EF_func = EF_func, GM = GM)
                    if tot_FR != 0:
                        try:
                            FRxy_ab[i, j] = tot_FR
                        except:
                            pdb.set_trace()
                    
                    FRxy_ab[j, i] = FRxy_ab[i, j]
                    Missed +=missed
                    Greater_one += gOne
 
        FRxy[ab] = FRxy_ab
    
    Missed = np.unique(np.array(Missed))
    Greater_one = np.unique(np.array(Greater_one))           
    return FRxy, Missed, Greater_one


"""Expected Immunity Efficacy as a function of COVI19 variant proportions"""
def Immunity_one_present_variant(present_variant, t, PK_dframe, infection_data, tested_variant_list, 
                                 variant_name, variant_proportion, Ab_classes, IC50xx, Cross_react_dic):
    IM_res = np.zeros((len(infection_data), len(t)))
    x = list(variant_name).index(present_variant)
    
    for l in range(len(infection_data)):
        # expected num of people infected with variant x at time l
        infected_l = infection_data[l]*variant_proportion[l, x]
        
        for k in range(len(t)):
            if l <= k:
                antibody_level = PK_dframe.loc[k - l][1:]
                for j in range(len(tested_variant_list)):
                    y = list(variant_name).index(tested_variant_list[j])
                    
                    IC50xy = [Cross_react_dic[ab][x, y]*IC50xx[ab] for ab in Ab_classes]
                    Ab_Eff = efficacy_n_antibodies(antibody_level, IC50xy)
                    
                    # expected num of people exposed to variant y at time k given that they were infected with variant x at time l
                    IM_res[l, k] += infected_l*Ab_Eff
                    
    immunity_to_variant = np.sum(IM_res, axis = 0)
    
    return immunity_to_variant

""" Compute and plot several senarios """
import joblib as jb
from functools import partial
def Immunity_dynamics(t, PK_dframe, infection_data, present_variant_list, tested_variant_list, variant_name, variant_proportion, Ab_classes, IC50xx, Cross_react_dic, parallel = False):
    if not parallel:
        Expected_Immuned_list = []
        for present_variant in present_variant_list:
            immunity_to_variant = Immunity_one_present_variant(present_variant, t, PK_dframe, infection_data, 
                                                               tested_variant_list, variant_name, variant_proportion, 
                                                               Ab_classes, IC50xx, Cross_react_dic)
            
            Expected_Immuned_list.append(immunity_to_variant)
    else:
        pfunc = partial(Immunity_one_present_variant, t=t, PK_dframe = PK_dframe, infection_data = infection_data, 
                                                      tested_variant_list = tested_variant_list, variant_name = variant_name, 
                                                      variant_proportion = variant_proportion, Ab_classes = Ab_classes, 
                                                      IC50xx = IC50xx, Cross_react_dic = Cross_react_dic)
        
        Expected_Immuned_list = list(jb.Parallel(n_jobs = -1, prefer = "threads")(jb.delayed(pfunc)(present_variant_list[r]) for r in range(len(present_variant_list))))
    
    Expected_Immuned = np.sum(np.array(Expected_Immuned_list), axis = 0)
    return Expected_Immuned   
   

""" Projection of Infection probability as a funciton of COV19 Variant proportion """ 
def Antibody_infos(t, antibody_data, VE_data):
    t_max = antibody_data["t_max"]
    t_half = antibody_data["t_half"]
    params_dic = {"t_max":t_max, "t_half":t_half}
    c_t, c_dframe, ka, ke, cmax = Antibody(t = t, params_dic = params_dic, is_log = False, save = False)
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
""" First version of cross_reactivity map"""
def cross_reactivity_ver1(variant_name, escape_per_sites, Ab_classes, mut_sites_per_variant = mut_sites_library, EF_func = ("MAX","MEAN"), GM = False):
    FRxy = {}
    Missed = []
    Greater_one = []
    for ab in Ab_classes:
        FRxy_ab = np.ones((len(variant_name), len(variant_name)))
        for i in range(len(variant_name)):
            for j in range(len(variant_name)):
                if (i !=j) & (j>i):
                    tot_fub_xy, missed, gOne = FUB(variant_name[i], variant_name[j], escape_per_sites, ab, mut_sites_per_variant, EF_func, GM)
                    if tot_fub_xy != 0:
                        try:
                            FRxy_ab[i, j] = 100*tot_fub_xy/(1 - tot_fub_xy) #100/((1/tot_fub_xy) - 1)
                        except:
                            pdb.set_trace()
                    
                    FRxy_ab[j, i] = FRxy_ab[i, j]
                    Missed +=missed
                    Greater_one +=gOne
                    """
                    tot_fub_yx = FUB(variant_name[j], variant_name[i], escape_per_sites, ab, mut_sites_per_variant)
                    if tot_fub_yx != 0:
                        FRxy_ab[j, i] = 100/((1/tot_fub_yx) - 1)
                    """
        FRxy[ab] = FRxy_ab
    
    Missed = np.unique(np.array(Missed))
    Greater_one = np.unique(np.array(Greater_one))           
    return FRxy, Missed, Greater_one

def FUB(variant_1, variant_2, escape_per_sites, ab, mut_sites_per_variant = mut_sites_library, EF_func = ("MAX","MEAN"), GM = False):
    sites_1 = set(mut_sites_per_variant[variant_1])
    sites_2 = set(mut_sites_per_variant[variant_2])
    
    sites = list(sites_1.symmetric_difference(sites_2))
    if EF_func[1] == "MEAN":
        escape_data = escape_per_sites["mean_escape_fraction_per_site"].values
    elif EF_func[1] == "MAX":
        escape_data = escape_per_sites["max_escape_fraction_per_site"].values
    
    #### Bound Escape data####
    escape_data[escape_data >= 1.] = 0.99
    escape_data[escape_data <= 1./101] = 1./101
    
    Missed = []
    Greater_one = []
    all_bound = 1
    for s in sites:
        S = (escape_per_sites["site"].values).astype(str) == str(s)
        AB = (escape_per_sites["group"].values == ab)
        where_s_ab = np.where(S & AB)[0]
        try:
            if len(where_s_ab != 0):
                if EF_func[0] == "MAX":
                    EF = np.max(escape_data[where_s_ab].copy())
     
                elif EF_func[0] == "MEAN":
                    EF = np.mean(escape_data[where_s_ab].copy())
 
                all_bound *= (1 - min(0.99, EF))
                
                if np.any(escape_data[where_s_ab]>1):
                    cond_list = (escape_per_sites["condition"].values)[(where_s_ab) & (escape_data > 1)]
                    for cond in cond_list:
                        Greater_one.append("Escape fraction = %.2f at site %s for AB %s and condition %s "%(escape_data[where_s_ab], s, ab, cond))
            else:
                Missed.append("No Escape fraction: at site %s and AB %s"%(s, ab))
            
        except: # if there is no data for some mutation sites
            pdb.set_trace()
    return 1 - all_bound, Missed, Greater_one # ab needs to bind just one of the sites