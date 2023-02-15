#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 15:27:56 2023

@author: raharinirina
"""
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import pdb

#### Visualisation ###  
def PreFig(xsize = 12, ysize = 12):
    '''
    @brief: customize figure parameters
    '''
    matplotlib.rc('xtick', labelsize=xsize) 
    matplotlib.rc('ytick', labelsize=ysize)


def Display(t, Y, is_log, labels, figsize = (7, 7), xysize = (15,15), labsize = 20, save_to = "test", xval = "x", yval = "f(x)", linewidth = 3):
    PreFig(xsize = xysize[0], ysize = xysize[1])
    fig = plt.figure(figsize = figsize)
    ax = fig.add_subplot(1, 1, 1)
    for i in range(Y.shape[0]):
        plt.plot(t, Y[i, :], label = labels[i], linewidth = linewidth)
        
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

import seaborn as sns


def Heatmap(data_dic, row_labels, col_labels, annotsize = 12 ,colormap = None, save_to = "cross_reactivity", sub_fig_size = 6):
    dLabs = list(data_dic.keys())
    N = len(dLabs)
    if N%2 == 0:
        F = int(N/2) 
    else:
        F = (N//2) + 1
    PreFig()
    fig = plt.figure(figsize = (F*sub_fig_size, F*sub_fig_size))
    num = 1
    for i in range(N):
        ax = fig.add_subplot(F, F, num) 
        cMap = sns.heatmap(data = data_dic[dLabs[i]], cmap = colormap, xticklabels = row_labels, yticklabels = col_labels, 
                       cbar = False, annot = True, fmt = ".2f", annot_kws = {"size":annotsize})
        num +=1
        plt.title("FR to Ab "+dLabs[i], fontsize = 16)
    
    #plt.subplots_adjust(top=0.965, bottom=0.095, left=0.080, right = 0.75, hspace=0.175, wspace=0.175)
    plt.subplots_adjust(hspace=0.4, wspace=0.4)    
    pdf = PdfPages(save_to + ".pdf")
    pdf.savefig(fig, bbox_inches = "tight")
    pdf.close()
    return fig
    
