import getopt, sys

def parse_args():
    full_cmd_arguments = sys.argv
    argument_list = full_cmd_arguments[1:]
    short_options = "i:r:b:fvo"
    long_options = ["inpu=", "ref_folder=", "background=","figure", "verbose", "output"]

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

    for a, v in arguments:
        if a in ("-v", "--verbose"):
            verbose = True
        elif a in ("-f", "--figure"):
            figure = True
        elif current_argument in ("-o", "--output"):
            outP = True
        elif current_argument in ("-i", "--input"):
            inP = v
        elif current_argument in ("-r", "--ref_folder"):
            ref = v
        elif current_argument in ("-b", "--background="):
            bkg = v
    
    if (not inP):
        print("Missing input file")
        sys.exit(2)
    if (not bkg):
        print("Missing background file")
        sys.exit(2)
    
    return arguments

def read_reference_data(ref_folder):
    return

def read_sc_data():
    return

def write_output():
    return

if __name__ == '__main__':

    args = parse_args()
    #ref_peak_sets, bk_peaks = read_reference_data(args['ref_folder'], args['background'])
    #sc_data = read_sc_data(args['input'])
    #compute_common_peaks(ref_peak_sets, bk_peaks)
    #scores = compute_fischer_exact_scores(sc_data, ref_peak_sets)
    #write_output(scores, args['output'])
    