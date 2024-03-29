#!/usr/bin/env python
# coding: utf-8

import os
import math
import seaborn as sns

import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')

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
    ## input: s1, s2, sbk, three sets
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
    ## TODO: remove background argument
    ## TODO: change to argparse?
    ## default values?
    short_options = "i:r:b:f:m:t:vo"
    long_options = ["input=", "ref_folder=","figure_cols=","multiple_cores=","threshold=", "verbose", "output"]

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
    multiple = 1
    threshold = 50

    for a, v in arguments:
        if a in ("-v", "--verbose"):
            verbose = True
        elif a in ("-f", "--figure_cols"):
            figure = int(v)
        elif a in ("-m", "--multiple_cores"):
            multiple = v
        elif a in ("-o", "--output"):
            outP = True
        elif a in ("-i", "--input"):
            inP = v
        elif a in ("-r", "--ref_folder"):
            ref = v
        elif a in ("-t", "--threshold"):
            threshold = v
        
    configs = {'figure':figure,
              'verbose':verbose,
               'inP': inP,
               'outP': outP,
               'ref': ref,
               'multiple': multiple,
               'threshold': threshold,
              }
        
    return configs

def check_input(configs):
    if (not configs['inP']):
        print("Missing input file")
        sys.exit(2)
    if (not configs['ref']):
        print("Missing reference folder address")
        sys.exit(2)

def intConvert(base):
    sublist = base.split('_')
    return [sublist[0], int(sublist[1]), int(sublist[2])]        

def read_input(input_filename, verbose = False):
    if(verbose): print("Reading Input File")
    
    data = pd.read_csv(inP, index_col = 0)
    id2peak = data.index.tolist()
    id2peak = list(map(intConvert, id2peak))
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
    
    ref_bk_peak = bdO.BedScanPd(bkg_file)
    
    ref_subtype_peaks = []
    for c, f1 in fn_target:
        scan = bdO.BedScanPd(f1)    
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
        ## TODO data is referenced externally, mat1 shoud be just mat       
        #mat = np.array(id2peak)
        for i in range (mat.shape[1]):
            s1 = set([id2peak[j] for j in np.where(mat[:, i] > 0.5)[0]])
            ## TODO change the following to be just
            ## tmp = scorefun((s1, set2_ref, interset))
            tmp = []
            for c, s2 in set2_ref:
                s2 = set(s2[:,0])
                if not (s1 & s2):
                    tmp.append(0)
                else:
                    tmp.append(get_fisher_exact(s1, s2, interset)[0])
            scores.append(tmp)
    return scores

def initialize_TSNE(mat1, verbose):
    mat = mat1.T / np.sum(mat1.T, axis = 1, keepdims=True)
    clf0 = PCA(n_components = 10)
    if(verbose): clf1 = TSNE(n_components = 2, verbose=1)
    else: clf1 = TSNE(n_components = 2)
    df_tsne = clf1.fit_transform(clf0.fit_transform(mat))
    return df_tsne

def create_plot(set_ref, figx, tsne, score_array, verbose = False):
    figure = 3
    dist = len(set2_ref)/figure
    if (dist.is_integer()): 
        dist = int(dist)
        fig, axes = plt.subplots(dist, figure, figsize = ((figure *4), (dist*3)))
    else:
        dist = math.ceil(dist)
        fig, axes = plt.subplots(dist, figure, figsize = ((figure *4), (dist*3)))
        #figure = 1
    
    x = 0
    y = 0
    for i, (c, _) in enumerate(set2_ref):
        if(figure == 1 or dist == 1): ax = axes[i]
        else: ax = axes[y][x]
        vmax, vmin = np.percentile(scoresNP[:, i], 99), np.percentile(scoresNP[:, i], 1)
        ax.scatter(df_tsne[:,0], df_tsne[:,1], 
                s = 5, c = scoresNP[:, i], cmap = 'RdBu_r', vmax = vmax, vmin = vmin)
        ax.set_title(c)
        ax.set_xticks([])
        ax.set_yticks([])
        x = x+1
        if (x == figure):
            x = 0
            y = y+ 1
    while(y != dist):
        ax = axes[y][x]
        ax.set_visible(False)
        x = x+1
        if (x == figure):
            x = 0
            y = y+ 1

    if(outP): plt.savefig(inP + '.png')
    plt.show()
    return plt
    
def output(ref, intersect,set_ref,scores):
    #write output file for the intersection
    IntersectFile = open(ref+'Intersection.tsv', "w+")
    for i, (c, _) in enumerate(set_ref):
        IntersectFile.write(c +"\t")
    IntersectWriter= csv.writer(IntersectFile, delimiter='\t')
    for line in intersect: IntersectWriter.writerow(line)
    IntersectFile.close()

    #write output file for the intersection of cell types
    for c, s2 in set_ref:
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
    
    ## TODO: add one or two example command
    configs = parse_args(sys.argv)
    check_input(configs)
    input_filename = configs['inP']
    data, id2peak = read_input(input_filename, verbose = False)
    ref_subtype_peaks, ref_bk_peak = read_reference_data(configs['ref'])
    
    if(configs['verbose']): print("Calculating Bed Overlaps")
    intersect = bdO.BedOverlap(id2peak, ref_bk_peak, configs['threshold'])
    #TODO: Add option to set overlap threshold

    set_ref = []
    for i, (ct, peaks) in enumerate(ref_subtype_peaks):
        peaks_intersect = bdO.BedOverlap(id2peak, peaks, configs['threshold'])
        ref_subtype_peaks[i] = (ct, peaks_intersect)
        set_ref.append((ct, peaks_intersect))

    if(configs['verbose']): print("Calculating Fischer Exact Scores")
    
    ## TODO if make changes to read_input, then
    ## mat = data.values
    mat = data.drop("id2peak",axis=1).to_numpy()

    scores = compute_enrichment_score(mat, intersect, id2peak, set_ref, configs['multiple'])
    
    ## TODO: spase between if and ( ?
    if(configs['verbose']): print("Intializing TSNE")
    df_tsne = initialize_TSNE(mat, configs['verbose'])
    
    plt = create_plot(set_ref, configs['figure'], df_tsne, np.array(scores),configs['verbose'])
    if(configs['outP']): plt.savefig(input_filename + '.png')
    plt.show()

    if(configs['verbose'] and configs['outP']): print("Outputting")
    if(configs['outP']): output(configs['ref'], intersect, set_ref, scores)

    sys.exit(0)
