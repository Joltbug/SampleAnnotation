import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import fisher_exact
from multiprocessing import Pool
from statsmodels.stats.multitest import multipletests
from scipy.io import mmread

def get_fisher_exact(var_in):
    s1, s2, sbk = var_in
    n1 = len(s1 & s2)
    n2 = len(s1 & (sbk - s2))
    n3 = len((sbk - s1) & s2)
    n4 = len((sbk - s1) & (sbk - s2))

    mat = [[n1, n2], [n3, n4]]

    oddsratio, pvalue = stats.fisher_exact(mat, 'greater')
    return oddsratio, pvalue

def process_one_sample(row2cols, id2cluster, id2peak, id2ct, pool = None):
    sbk = set([i for i in range(id2cluster.shape[0])])
    for ct in id2ct:
        set_cluster = set(np.where(id2cluster == ct)[0])
        inputs = [(cols, set_cluster, sbk) for r, cols in row2cols.items()]
        #print(inputs[:10])
        if pool is None:
            or_and_pvals = list(map(get_fisher_exact, inputs))
        else:
            or_and_pvals = pool.map(get_fisher_exact, inputs)
        #print(or_and_pvals )
        reject, adjp, _, _ = multipletests([p[1] for p in or_and_pvals])
        diffpeaks = [id2fea[i] for i, r in enumerate(reject) if r]
        pd.DataFrame(diffpeaks).to_csv('%s.bed' % ct, index = False, sep = '\t', header = False)
    return 
        
cl2ct = {4:'LMPP',
5:'CLP',
1:'HSC_MPP',
2:'MEP',
3:'CMP_BMP',
8:'GMP',
6:'Pro-B',
9:'MDP',
10:'pDC',
7:'Pre-B',
21:'Naive CD4 T1',
22:'Naive CD4 T2',
20:'Mature NK2',
17:'Basophil',
16:'Plasma cell',
19:'Mature NK1',
18:'Immature NK',
15:'Memory B',
11:'cDC',
13:'Monocyte 2',
12:'Monocyte 1',
14:'Naive B',
28:'Naive CD8 T3',
29:'Central memory CD8 T',
27:'Naive CD8 T2',
24:'Memory CD4 T',
23:'Naive Treg',
26:'Naive CD8 T1',
25:'Treg',
30:'Effector memory CD8 T',
31:'Gamma delta T'}

if __name__ == '__main__':
    
    df_barcode = pd.read_csv('./GSE129785_scATAC-Hematopoiesis-All.cell_barcodes.txt.gz', sep = '\t')
    id2fea = pd.read_csv('./GSE129785_scATAC-Hematopoiesis-All.peaks.txt.gz').Feature.values
    id2fea = [f.split('_') for f in id2fea]

    cluster = np.array([cl2ct[int(c[7:])] for c in df_barcode.Clusters.values])
    
    data_mat = mmread('./GSE129785_scATAC-Hematopoiesis-All.mtx')
    
    row2cols = {}
    for i, j in zip(data_mat.row, data_mat.col):
        if i not in row2cols:
            row2cols[i] = set()
        row2cols[i].add(j)
    id2ct = [cl2ct[i] for i in range(1, 32)]
    
    pool = Pool(30)
    process_one_sample(row2cols, cluster, id2fea, id2ct, pool)