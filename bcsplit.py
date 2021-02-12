import sys
sys.path.append("./nptools")
import pandas as pd
import os
import importlib
import pickle
from nanopore_scripts import *
datasetfname = os.path.join("input","datasets.xlsx")
datapath = os.path.join("/","central","groups","murray-biocircuits","ashur","nanopore")

readsname = "allreads.fastq"
outname = "simprec"

df_data = pd.read_excel(datasetfname,sheet_name="alldata")
df_inducers = pd.read_excel(datasetfname,sheet_name="inducers",header=13)
df_seqs = pd.read_excel(datasetfname,sheet_name="sequences")


seqlist = list(df_data.date_sequenced.unique())

processreads = 100000
frontchecklength = 175
threshfrac = 0.3


for seqdataset in seqlist:
    
    subdf = df_data[df_data.date_sequenced==seqdataset]
    #prefix sequence(s)
    prefixes = list(subdf.prefix.unique())
    prefixseqslist = []
    for prefix in prefixes:
        if(pd.isna(prefix)):
            prefixlist = []
        elif("," in prefix):
            prefixlist = prefix.replace("[","").replace("]","").split(",")    
        else:
            prefixlist = [prefix]
        for pref in prefixlist:
            prefixseqslist += [df_seqs[df_seqs.name==pref].sequence.iloc[0]]
    #barcode sequences
    barcodes = list(subdf.barcode.unique())
    bcseqs = {}
    bcnames = []
    condlist = []
    for bc in barcodes:
        bcseq = df_seqs[df_seqs.name==bc].sequence.iloc[0]
        bcseqs[bc] = (bcseq,rc(bcseq))
        #CONDITIONS
        condstr = ""
        for i in range(4):
            cond = subdf[subdf.barcode==bc]["c"+str(i+1)].iloc[0]
            if(pd.isna(cond)):
                break
            if(condstr != ""):
                condstr += "_"+str(cond)
            else:
                condstr+=str(cond)
        condlist += [condstr]
    #plasmid barcodes
    plasbcs = []
    plasmidbarcodes = zip(list(subdf.variable1.unique()),list(subdf.variable2.unique()))
    for plasbc1,plasbc2 in plasmidbarcodes:
        if(pd.isna(plasbc1)):
            plasbcs += [[]]
        elif(pd.isna(plasbc2)):
            plasbcs += [df_seqs[df_seqs.name==plasbc1].sequence.iloc[0]]
        else:
            plasbcs += [df_seqs[df_seqs.name==plasbc1].sequence.iloc[0],\
                         df_seqs[df_seqs.name==plasbc2].sequence.iloc[0]]
    #postfix
    postfix = []
    postfixes = list(subdf.suffix.unique())
    for pfix in postfixes:
        pfixseq = df_seqs[df_seqs.name==pfix].sequence
        if(len(pfixseq)==0):
          postfix += [""]
        else:
            postfix+=[pfixseq.iloc[0]]
    if(any([pd.isna(bcseqs[a]) for a in bcseqs])):
        continue
    fastqfilename = os.path.join(datapath,str(seqdataset),readsname)
    print(f"date of seq is {seqdataset}")
    print(f"get reads from {fastqfilename}")
    print(f"prefixseqslist is {prefixseqslist}")
    print(f"plasbcs is {plasbcs}")
    print(f"postfixseq is {postfix}")
    print(f"bcseqs is {bcseqs}")
    print(f"condnames is {condlist}")
    print(f"put output into {os.path.join(datapath,str(seqdataset),str(seqdataset)+'_'+str(outname)+'.py')}")
    ''
    allseqDict,seqstats,unsorted=barcodeSplitAndCountRecords(fastqfilename,bcseqs,\
                                                barcode_detection_threshold=len(list(bcseqs.values())[0])*threshfrac,\
                                                end_threshold=len(postfix[0])*threshfrac,\
                                                processreads=processreads,\
                                                variable_sequences=plasbcs,\
                                                prefix_sequence=prefixseqslist,\
                                                postfix_sequence=postfix,\
                                               prefix_detection_threshold=len(prefixseqslist[0])*threshfrac,\
                                               variable_sequence_threshold=len(plasbcs[0])*threshfrac,\
                                                frontchecklength=frontchecklength,visualize=False,progressbar = False)
    #'''
    print("we had {} forward, {} reverse, and {} where we couldn't tell".format(seqstats[0],seqstats[1],seqstats[2]))
    allseqDict['conditions']=condlist+["none"]
    with open(os.path.join(datapath,str(seqdataset),str(seqdataset)+'_'+str(outname)+'.pickle'),'wb') as f:
        pickle.dump(allseqDict,f,pickle.HIGHEST_PROTOCOL)