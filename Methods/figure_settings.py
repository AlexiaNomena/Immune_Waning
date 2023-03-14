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
import seaborn as sns

#### Visualisation ###  
def PreFig(xsize = 12, ysize = 12):
    '''
    @brief: customize figure parameters
    '''
    matplotlib.rc('xtick', labelsize=xsize) 
    matplotlib.rc('ytick', labelsize=ysize)


def Display(t, Y, is_log, labels, figsize = (7, 7), xysize = (15,15), labsize = 20, save_to = "test", xval = "x", yval = "f(x)", linewidth = 3, palette = None, linestyle = None):
    PreFig(xsize = xysize[0], ysize = xysize[1])
    fig = plt.figure(figsize = figsize)
    ax = fig.add_subplot(1, 1, 1)
    
    if linestyle is None:
        linestyle = ["-"]*Y.shape[0]
    
    if palette is None:
        for i in range(Y.shape[0]):
            plt.plot(t, Y[i, :], label = labels[i], linewidth = linewidth, linestyle = linestyle[i])
    else:
        col = sns.color_palette(palette, Y.shape[0])
        for i in range(Y.shape[0]):
            plt.plot(t, Y[i, :], label = labels[i], linewidth = linewidth, color = col[i], linestyle = linestyle[i])
        
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


def Heatmap(data_dic, row_labels, col_labels, annotsize = 12 , ticksize = None ,colormap = None, save_to = "cross_reactivity", 
            sub_fig_size = 6, num_row_col = None, cbar_kws = None, cbar_labsize = 12):
    
    dLabs = list(data_dic.keys())
    N = len(dLabs)
    
    pdf = PdfPages(save_to + ".pdf")    
    
    PreFig()
    

    for i in range(N):
        fig_sub = plt.figure(figsize = (sub_fig_size, sub_fig_size))
        
        if annotsize is not None:
            cMap = sns.heatmap(data = data_dic[dLabs[i]], cmap = colormap, xticklabels = row_labels, yticklabels = col_labels, 
                       cbar = False, annot = True, fmt = ".2f", annot_kws = {"size":annotsize})
        else:
            cMap = sns.heatmap(data = data_dic[dLabs[i]], cmap = colormap, xticklabels = row_labels, yticklabels = col_labels, 
                       cbar = True, annot = False, cbar_kws = cbar_kws)
            
            cMap.figure.axes[-1].yaxis.label.set_size(cbar_labsize)
            
        
        plt.title(dLabs[i], fontsize = cbar_labsize)
        if ticksize is None:
            plt.xticks(rotation = 90)
            plt.yticks(rotation = 0)
        else:
            plt.xticks(fontsize = ticksize[0], rotation = 90)
            plt.yticks(fontsize = ticksize[1], rotation = 0)
        
        try:
            pdf.savefig(cMap.figure, bbox_inches = "tight", dpi = 1200) 
        except:
            pdf.savefig(cMap.figure, bbox_inches = "tight") 
    
    if num_row_col is None:
        if N%2 == 0:
            F = int(N/2) 
        else:
            F = (N//2) + 1
        fig = plt.figure(figsize = (F*sub_fig_size, F*sub_fig_size))
    else:
        fig = plt.figure(figsize = (num_row_col[1]*sub_fig_size, num_row_col[0]*sub_fig_size))
    
    num = 1
    for i in range(N):
        if num_row_col is None:
            ax = fig.add_subplot(F, F, num)
        else:
            ax = fig.add_subplot(num_row_col[0], num_row_col[1], num)
        
        if annotsize is not None:
            cMap = sns.heatmap(data = data_dic[dLabs[i]], cmap = colormap, xticklabels = row_labels, yticklabels = col_labels, 
                       cbar = False, annot = True, fmt = ".2f", annot_kws = {"size":annotsize})
        else:
            cMap = sns.heatmap(data = data_dic[dLabs[i]], cmap = colormap, xticklabels = row_labels, yticklabels = col_labels, 
                       cbar = True, annot = False, cbar_kws = cbar_kws)
            
            cMap.figure.axes[-1].yaxis.label.set_size(cbar_labsize)
            
        num +=1
        
        plt.title(dLabs[i], fontsize = cbar_labsize)
        if ticksize is None:
            plt.xticks(rotation = 90)
            plt.yticks(rotation = 0)
        else:
            plt.xticks(fontsize = ticksize[0], rotation = 90)
            plt.yticks(fontsize = ticksize[1], rotation = 0)
    #plt.subplots_adjust(top=0.965, bottom=0.095, left=0.080, right = 0.75, hspace=0.175, wspace=0.175)
    plt.subplots_adjust(hspace=0.4, wspace=0.4)
    try:
        pdf.savefig(fig, bbox_inches = "tight", dpi = 1200)
    except:
        pdf.savefig(fig, bbox_inches = "tight")
    pdf.close()
    return fig
    
