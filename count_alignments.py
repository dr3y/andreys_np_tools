import pandas as pd
import os
sys.path.append(os.path.join(".","nptools"))
from nanopore_scripts import *
import pickle

if(__name__=="__main__"):
    date = "210507"
    subfolder = "210526minimap3"
    genome_used = "0-3intPlasmid_short"
    data_folder = os.path.join("central","groups","murray-biocircuits","ashur","nanopore",date,subfolder)
    bclist = [f"BC{a:02d}" for a in range(6,11)]
    bcdict = {a:f"{a}_{genome_used}" for a in bclist}


    int_strings1 = allpaths([["p1","p2"],["p1","p2"]])
    int_strings2 = allpaths([["p1","p2"],["p1","p2"],["p1","p2"]])
    int_strings3 = allpaths([["p1","p2"],["p1","p2"],["p1","p2"],["p1","p2"]])

    genomelist = ["_".join(a) for a in [["p1"],["p2"]]+int_strings1+int_strings2+int_strings3]

    gdict = {a.replace("p","").replace("g",""):a for a in genomelist}
    outdict = countSequencing(bcdict,gdict,(429,466),data_folder,minlen=1000)
    print(outdict)
    pickle.dump(outdict,"outfile.pickle")