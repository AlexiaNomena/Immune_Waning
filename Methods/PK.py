#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 13:14:30 2023

@author: raharinirina
"""
import numpy as np
import pandas as pd
from scipy.optimize import root
import warnings
warnings.filterwarnings("ignore")
import pdb

def Antibody(t, params_dic, is_log = False, save_to = "test", ka_solver = "lm"):
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
    save_to : string Path/filename (saved to csv)
    ka_solver : root solver method for finding absorbtion rate ka (see scipy.optimize.root)

    Returns
    -------
    Antibody concentration at time t.

    """
    t_max = params_dic["t_max"]
    t_half = params_dic["t_half"]
    
    # Antibody elimination rate
    ke = np.log(2)/t_half
    
    # Antibody absorption rate
    guess = np.ones(len(t_max))
    ka = root(ka_solve, guess, args = (ke, t_max), method = ka_solver).x
    print("\n k_a was found correctly:", np.all(np.isclose(ka_solve(ka, ke, t_max), np.zeros(len(t_max)))), "\n")

    # Compute Normalized Concentration
    c_max = (np.exp(- ke*t_max) - np.exp(- ka*t_max))
    c_t = (np.exp(- ke[:, np.newaxis]*t) - np.exp(- ka[:, np.newaxis]*t))/c_max[:, np.newaxis]
       
    # Build pandas dataframe‚
    df = {}
    df["Days"] = t
    for i in range(len(t_max)):
        df["Ab_%d"%(i+1)] = c_t[i, :]
        
    df = pd.DataFrame(df)
    
    # save data
    df.to_csv(save_to)
         
    return c_t, df

def ka_solve(ka, ke, t_max):
    if np.all(ka)>0:
        res = np.divide(t_max*(ka - ke) - (np.log(ka) - np.log(ke)), (ka - ke), out = np.ones(len(ke)), where = (ka - ke)!=0)
    else:
        res = 1
    return res
        