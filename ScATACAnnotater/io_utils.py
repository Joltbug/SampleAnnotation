import SortedNoDupeBedOverlap as bdO
import os
import yaml
import pandas as pd

def read_bed_file(fileName, delimiter='\t', skiprows=0, header=None):
    ## default setting assumed input: tab separated, no header, no index
    ## assumed 0, 1, 2 columns are chromosome, start, end
    ## output: a list of peaks (chr:str, start:int, end:int)     
    
    data = pd.read_csv(fileName, delimiter = delimiter, 
                       skiprows = skiprows, header = header)
    output = data.values.tolist()
    out = list((row[0], int(row[1]), int(row[2])) for row in output)
    return out

def read_yaml(file):
    with open(file) as reffile: 
        info = yaml.load(reffile,Loader=yaml.FullLoader)
    return info

def read_reference_data(ref_folder, verbose = False):
    ## reading reference data
    
    refinfo = read_yaml(os.path.join(ref_folder, 'info.yaml'))
    bkg_file = os.path.join(ref_folder, refinfo['background'])
    fn_target = [(c, os.path.join(ref_folder, refinfo['cell_type'][c])) for c in refinfo['cell_type']]

    if (verbose):
        print('Using background file: ',refinfo['background'])
        print('With cell types:')
        for cell in refinfo['cell_type']:
            print (cell,':',refinfo['cell_type'][cell])    
    
    ref_bk_peak = read_bed_file(bkg_file)
    
    ref_subtype_peaks = []
    for c, f1 in fn_target:
        scan = set(read_bed_file(f1))    
        ref_subtype_peaks.append((c, scan))
        
    return ref_subtype_peaks, ref_bk_peak

def read_csv_input(input_filename, verbose = False):
    if(verbose): print("Reading Input File")

    data = pd.read_csv(input_filename)

    data.rename( columns={'Unnamed: 0':'id2peak'}, inplace=True)
    id2peak = []
    for c in data['id2peak']:
        temp = c.split('_')
        id2peak.append((temp[0], int(temp[1]), int(temp[2])))

    id2peak.sort( key = lambda x: (x[0], x[1]))    
    return data, id2peak

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