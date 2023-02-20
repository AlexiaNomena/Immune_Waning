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


def FUB(variant_1, variant_2, escape_per_sites, ab, mut_sites_per_variant = mut_sites_library):
    sites_1 = set(mut_sites_per_variant[variant_1])
    sites_2 = set(mut_sites_per_variant[variant_2])
    
    sites = list(sites_1.symmetric_difference(sites_2))
    escape_data = escape_per_sites["mean_escape_fraction_per_site"].values
    Missed = []
    Greater_one = []
    all_bound = 1
    for s in sites:
        S = (escape_per_sites["site"].values).astype(str) == str(s)
        AB = (escape_per_sites["ab_class"].values == ab)
        where_s_ab = np.where(S & AB)[0]
        try:
            if len(where_s_ab != 0):
                all_bound *= (1 - min(0.99, escape_data[where_s_ab]))
                if escape_data[where_s_ab]>1:
                    Greater_one.append("Escape fraction = %.2f at site %s and AB %s"%(escape_data[where_s_ab], s, ab))
            else:
                Missed.append("No Escape fraction: at site %s and AB %s"%(s, ab))
            
        except: # if there is no data for some mutation sites
            pdb.set_trace()
    return 1 - all_bound, Missed, Greater_one # ab needs to bind just one of the sites

def cross_reactivity(variant_name, escape_per_sites, Ab_classes, mut_sites_per_variant = mut_sites_library):
    FRxy = {}
    Missed = []
    Greater_one = []
    for ab in Ab_classes:
        FRxy_ab = np.ones((len(variant_name), len(variant_name)))
        for i in range(len(variant_name)):
            for j in range(len(variant_name)):
                if (i !=j) & (j>i):
                    tot_fub_xy, missed, gOne = FUB(variant_name[i], variant_name[j], escape_per_sites, ab, mut_sites_per_variant)
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

"""Expected Immunity Efficacy as a function of COVI19 variant proportions"""
def Immunity_one_neutralized_variant(t, PK_dframe, infection_data, neutralized_variant, variant_list, variant_name, variant_proportion, Ab_classes, IC50xx, Cross_react_dic):
    IM_res = np.zeros((len(infection_data), len(t)))
    x = list(variant_name).index(neutralized_variant)
    
    for l in range(len(infection_data)):
        # expected num of people infected with variant x at time l
        infected_l = infection_data[l]
        
        for k in range(len(t)):
            if l <= k:
                antibody_level = PK_dframe.loc[k - l][1:]
                for j in range(len(variant_list)):
                    y = list(variant_name).index(variant_list[j])
                    
                    IC50xy = [Cross_react_dic[ab][x, y]*IC50xx[ab] for ab in Ab_classes]
                    Ab_Eff = efficacy_n_antibodies(antibody_level, IC50xy)
                    
                    # expected num of people exposed to variant y at time k given that they were infected with variant x at time l
                    infected_l_exposed_yk = variant_proportion[k, y]*infected_l
                    IM_res[l, k] += infected_l_exposed_yk*Ab_Eff
                    
    immunity_to_variant = np.sum(IM_res, axis = 0)
    
    return immunity_to_variant

""" Compute and plot several senarios """
def Immunity_dynamics(t, PK_dframe, infection_data, neutralized_variant_list, variant_list, variant_name, variant_proportion, Ab_classes, IC50xx, Cross_react_dic):
    Expected_Immuned_list = []
    for neutralized_variant in neutralized_variant_list:
        immunity_to_variant = Immunity_one_neutralized_variant(t, PK_dframe, infection_data, 
                                neutralized_variant, variant_list, 
                                variant_name, variant_proportion, 
                                Ab_classes, IC50xx, Cross_react_dic)
        
        Expected_Immuned_list.append(immunity_to_variant)
        
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
        