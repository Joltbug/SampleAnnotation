import scipy.stats as stats
import numpy as np
from multiprocessing import Pool

def get_fisher_exact(s1, s2, sbk):
    n1, n2, n3, n4 = 0, 0, 0, 0
    for s in s1:
        if s not in sbk:
            continue
        if s in s2:
            n1 += 1
        else:
            n2 += 1
    for s in s2:
        if s not in sbk:
            continue    
        if s not in s1:
            n3 += 1
            
    n4 = len(sbk) - (n1 + n2 + n3)
    mat = [[n1, n2], [n3, n4]]

    oddsratio, pvalue = stats.fisher_exact(mat, 'greater')
    return oddsratio, pvalue

def list2id_dict(alist):
    return {p:i for i, p in enumerate(alist)}
    
def list2id(alist, id_dict, whitelist_set = None):
    if whitelist_set is None:
        return [id_dict[p] for p in alist]    
    else:
        return [id_dict[p] for p in alist if p in whitelist_set]  
    
def _scorefun(var_in):
    s1, ref_sets, bks = var_in
    tmp = []
    for c, s2 in ref_sets:
        if not (s1 & s2):
            tmp.append(0)
        else:
            tmp.append(get_fisher_exact(s1, s2, bks)[0])
    return tmp

def compute_enrichment_score(mat_in, bk_peaks, id2peak, set2ref_peak, num_cores = 1):
    ## mat: an numpy array of cells x peaks
    ## 
    
    col_to_use = np.array([i for i, p in enumerate(id2peak) if p in bk_peaks])
    peak_to_use = [p for i, p in enumerate(id2peak) if p in bk_peaks]
    peak2id = list2id_dict(peak_to_use)
    
    mat_in = mat_in[:, col_to_use]
    bk_ids = set(list2id(bk_peaks, peak2id))
    
    set2ref_id = []
    for c, peaks in set2ref_peak:
        set2ref_id.append((c, set(list2id(peaks, peak2id, bk_peaks))))
        
    if(num_cores > 1):
        pool = Pool(num_cores)
        map_fun = pool.map
    else:
        map_fun = map
        
    var_in = [] 
    for i in range(mat_in.shape[0]):
        s1 = set(mat_in[i].nonzero()[0])
        var_in.append((s1, set2ref_id, bk_ids))
    scores = map_fun(_scorefun, var_in)

    return scores
