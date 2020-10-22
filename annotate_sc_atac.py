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


def get_fisher_exact(s1, s2, sbk):
    
    n1 = len(s1 & s2)
    n2 = len(s1 & (sbk - s2))
    n3 = len((sbk - s1) & s2)
    n4 = len((sbk - s1) & (sbk - s2))

    mat = [[n1, n2], [n3, n4]]

    oddsratio, pvalue = stats.fisher_exact(mat, 'greater')
    return oddsratio, pvalue

def scorefun(input):
    tmp = []
    for c, s2 in input[1]:
        s2 = set(s2[:,0])
        if not (input[0] & s2):
            tmp.append(0)
        else:
            tmp.append(get_fisher_exact(input[0], s2, input[2])[0])
    return tmp

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

def read_input(input_filename, verbose = False):
    if(verbose): print("Reading Input File")

    data = pd.read_csv(input_filename)

    data.rename( columns={'Unnamed: 0':'id2peak'}, inplace=True)
    id2peak = []
    for c in data['id2peak']:
        temp = c.split('_')
        id2peak.append((temp[0], int(temp[1]), int(temp[2])))

    id2peak.sort( key = lambda x: (x[0], x[1]))    
    return data, id2peak

def read_yaml(file):
    with open(file) as reffile: 
        info = yaml.load(reffile,Loader=yaml.FullLoader)
    return info

def read_reference_data(ref_folder, verbose = False):
    refinfo = read_yaml(os.path.join(ref_folder, 'info.yaml'))
    bkg_file = os.path.join(ref_folder, refinfo['background'])
    fn_target = [(c, os.path.join(ref_folder, refinfo['cell_type'][c])) for c in refinfo['cell_type']]

    if (verbose):
        print('Using background file: ',refinfo['background'])
        print('With cell types:')
        for cell in refinfo['cell_type']:
            print (cell,':',refinfo['cell_type'][cell])    
    
    ref_bk_peak = bdO.BedScan(bkg_file)
    
    ref_subtype_peaks = []
    for c, f1 in fn_target:
        scan = bdO.BedScan(f1)    
        ref_subtype_peaks.append((c, scan))
        
    return ref_subtype_peaks, ref_bk_peak

def compute_enrichment_score(mat, intersect, id2peak, set2_ref, num_cores):
    scores = []
    interset = set(intersect[:,0])
    if(num_cores >1):
        pool = Pool(num_cores)
        for i in range (mat.shape[1]):
            s1 = set([id2peak[j] for j in np.where(mat[:, i] > 0.5)[0]])
            input = [(s1, set2_ref, interset)]
            pool.apply_async(scorefun, input, callback=scores.append)
        pool.close()
        pool.join()
    else:
        mat1 = data.drop("id2peak",axis=1).to_numpy()
        for i in range (mat1.shape[1]):
            s1 = set([id2peak[j] for j in np.where(mat1[:, i] > 0.5)[0]])
            tmp = []
            for c, s2 in set2_ref:
                s2 = set(s2[:,0])
                if not (s1 & s2):
                    tmp.append(0)
                else:
                    tmp.append(get_fisher_exact(s1, s2, interset)[0])
            scores.append(tmp)
    return scores

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
    
def output(ref, intersect,set_ref,scores):
    #write output file for the intersection
    IntersectFile = open(fn_bk1+'Intersection.tsv', "w+")
    for i, (c, _) in enumerate(set_ref):
        IntersectFile.write(c +"\t")
    IntersectWriter= csv.writer(IntersectFile, delimiter='\t')
    for line in intersect: IntersectWriter.writerow(line)
    IntersectFile.close()

    #write output file for the intersection of cell types
    for c, s2 in set2_ref:
        setFile = open( c+")_Intersection.tsv" , "w+")
        setFile.write("Chr in Data\tChr in Bk\tOverlap percent\n")
        setWriter= csv.writer(setFile, delimiter='\t')
        for line in s2: setWriter.writerow(line)
        setFile.close()

     #write output file for fischer exact scores
    ScoreFile = open('DatafischerScores.tsv', "w+")
    ScoreWriter= csv.writer(ScoreFile, delimiter='\t')
    for line in scores: ScoreWriter.writerow(line)
    ScoreFile.close()    
    
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

    if(configs['verbose']): print("Calculating Fischer Exact Scores")

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
