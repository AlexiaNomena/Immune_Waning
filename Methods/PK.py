#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 13:14:30 2023

@author: raharinirina
"""
import numpy as np
import pandas as pd

def Antibody(t, params_dic, is_log = False, save_to = "test"):
    """
    @brief: Compute Antibody Concentration as a function of time for N antibody classes
    
    Parameters
    ----------
    t : T time points ndarray (T) 
    params_dic : dictionary of parameters
                params_dic["t_max"] = t_{max}, ndarray (N, )
                params_dic["t_half"] = t_{1/2}, ndarray (N, ) 
    is_log : bool, optional
             True if return the log of concentration. The default is False.

    Returns
    -------
    Antibody concentration at time t.

    """
    t_max = params_dic["t_max"]
    t_half = params_dic["t_half"]
    
    # Antibody elimination rate
    k_e = np.log(2)/t_half
    
    # Antibody absorption rate
    k_a = np.ones(len(t_max)) # TBA
    
    c_max = (np.exp(- k_e*t_max) - np.exp(- k_a*t_max))
    
   
    
    c_t = (np.exp(- k_e[:, np.newaxis]*t) - np.exp(- k_a[:, np.newaxis]*t))/c_max[:, np.newaxis]
       
    # Build dataframe
    df = {}
    df["Days"] = t
    for i in range(len(t_max)):
        df["Ab_%d"%(i+1)] = c_t[i, :]
    
    
    df = pd.DataFrame(df)
    
    # save data
    df.to_csv(save_to)
         
    return c_t, df
        