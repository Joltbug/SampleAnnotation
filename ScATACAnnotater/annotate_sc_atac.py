#!/usr/bin/env python
# coding: utf-8

import os
import seaborn as sns

import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

import SortedNoDupeBedOverlap as bdO

import pandas as pd
import numpy as np
import csv
import scipy
import yaml
import getopt, sys
import scipy.stats as stats
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from multiprocessing import Pool

def parse_args(full_cmd_arguments):
    argument_list = full_cmd_arguments[1:]
    short_options = "i:r:b:f:m:t:vo"
    long_options = ["input=", "ref_folder=", "background=","figure_cols=","multiple_cores=","threshold=", "verbose", "output"]

    try:
        arguments, values = getopt.getopt(argument_list, short_options, long_options)
    except getopt.error as err:
        print (str(err))
        sys.exit(2)

    figure = False
    verbose = False
    inP = False
    outP = False
    ref = False
    bkg = False
    multiple = 1
    threshold = 50

    for a, v in arguments:
        if a in ("-v", "--verbose"):
            verbose = True
        elif a in ("-f", "--figure_cols"):
            figure = v
        elif current_argument in ("-m", "--multiple_cores"):
            multiple = v
        elif current_argument in ("-o", "--output"):
            outP = True
        elif current_argument in ("-i", "--input"):
            inP = v
        elif current_argument in ("-r", "--ref_folder"):
            ref = v
        elif current_argument in ("-t", "--threshold"):
            threshold = v
        
    configs = {'figure':figure,
              'verbose':verbose,
               'inP': inP,
               'outP': outP,
               'ref': ref,
               'bkg': bkg,
               'multiple': multiple,
               'threshold': threshold,
              }
        
    return configs

def check_input(configs):
    if (not configs['inP']):
        print("Missing input file")
        sys.exit(2)
    if (not configs['bkg']):
        print("Missing background file")
        sys.exit(2)
    if (not configs['ref']):
        print("Missing reference folder address")
        sys.exit(2)


def initialize_TSNE(mat1):
    mat = mat1.T / np.sum(mat1.T, axis = 1, keepdims=True)
    clf0 = PCA(n_components = 10)
    if(verbose): clf1 = TSNE(n_components = 2, verbose=1)
    else: clf1 = TSNE(n_components = 2)
    df_tsne = clf1.fit_transform(clf0.fit_transform(mat))
    return df_tsne

def create_plot(set_ref, figsize, tsne, score_array, verbose = False):
    if(verbose): print("Generating Plots")
    dist = len(set2_ref)/figure
    if (dist.is_integer()): 
        dist = int(dist)
        fig, axes = plt.subplots(dist, figure, figsize = ((figure *4), (dist*3)))
    else:
        fig, axes = plt.subplots(len(set_ref), 1, figsize = (3,(len(set2_ref) *4)))
        figure = 1
    for i, (c, _) in enumerate(set_ref):
        if(figure == 1 or dist == 1): ax = axes[i]
        else: ax = axes[i // figure, i % figure]
        vmax, vmin = np.percentile(scoresNP[:, i], 99), np.percentile(scoresNP[:, i], 1)
        ax.scatter(df_tsne[:,0], df_tsne[:,1], 
                s = 5, c = scoresNP[:, i], cmap = 'RdBu_r', vmax = vmax, vmin = vmin)
        ax.set_title(c)
        ax.set_xticks([])
        ax.set_yticks([])
    axes[-1,-1].axis('off')
    return plt  
    
if __name__ == "__main__":      
    configs = parse_args(sys.argv)
    check_input(configs)
    input_filename = configs['inP']
    read_input(input_filename, verbose = False)
    ref_subtype_peaks, ref_bk_peak = read_reference_data(ref_folder)
    
    if(configs['verbose']): print("Calculating Bed Overlaps")
    intersect = bdO.BedOverlap(id2peak, ref_bk_peak, configs['threshold'])
    #TODO: Add option to set overlap threshold

    set_ref = []
    for i, (ct, peaks) in enumerate(ref_subtype_peaks):
        peaks_intersect = bdO.BedOverlap(id2peak, peaks, configs['threshold'])
        ref_subtype_peaks[i] = (ct, peaks_intersect)
        set_ref.append((ct, peaks_intersect))

    if configs['verbose']: 
        print("Calculating Fischer Exact Scores")

    mat1 = data.drop("id2peak",axis=1).to_numpy()

    scores = compute_enrichment_score_parallel(mat, intersect, id2peak, set_ref, configs['multiple']
    
    if(configs['verbose']): print("Intializing TSNE")
    df_tsne = initialize_TSNE(mat1)
    
    plt = create_plot(set2_ref, configs['figure'], df_tsne, np.array(scores),configs['verbose'])
    if(configs['outP']): plt.savefig(input_filename + '.png')
    plt.show()

    if(configs['verbose'] and configs['outP']): print("Outputting")
    if(configs['outP']): output(configs['ref'], intersect, set_ref, scores)

    sys.exit(0)
