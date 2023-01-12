#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 15:27:56 2023

@author: raharinirina
"""
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages

#### Visualisation ###  
def PreFig(xsize = 12, ysize = 12):
    '''
    @brief: customize figure parameters
    '''
    matplotlib.rc('xtick', labelsize=xsize) 
    matplotlib.rc('ytick', labelsize=ysize)



def Display(t, c_t, is_log, labels, figsize = (7, 7), xysize = (15,15), labsize = 20, save_to = "test", xval = "x", yval = "f(x)"):
    PreFig(xsize = xysize[0], ysize = xysize[1])
    fig = plt.figure(figsize = figsize)
    ax = fig.add_subplot(1, 1, 1)
    
    for i in range(c_t.shape[0]):
        plt.plot(t, c_t[i, :], label = labels[i])
        
    if is_log:
        plt.ylabel("$\ln$ %s"%yval, fontsize = labsize)
    else:
        plt.ylabel("%s"%yval, fontsize = labsize)  
    
    plt.xlabel(xval, fontsize = labsize)
    plt.legend(fontsize = labsize)
    
    if save_to[:-3] == "pdf":
        ### save figure in pdf ###
        pdf = PdfPages(save_to +".pdf")
        pdf.savefig(fig, bbox_inches = "tight")
        pdf.close()
    else:
        plt.savefig(save_to)
    
    return fig, ax