import edlib
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import ipywidgets as widgets
from copy import deepcopy as dc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import dnaplotlib
def rc(seq):
    return str(Seq(seq).reverse_complement())
def countBarcodeStats(bcseqs,chopseqs='none'):
    """this function uses edlib to count the number of matches to given bcseqs. 
        chopseqs can be left, right, both, or none. This tells the program to 
        chop off one barcode from either the left, right, both, or none of the
        ends."""
    x=[]
    o1list = []
    o2list = []
    pcount = []
    jcount = []
    pjcount = []
    jpcount = []
    all_lists = {}
    switch_lists = {}
    run_lists = {}
    first_last = {}
    for bc in bcseqs:
        seqs = bcseqs[bc]
        #simpRecDF[simpRecDF.Barcode==bc].Sequence

        seqschop = []
        curpcount = 0
        curjcount = 0
        curjpcount = 0
        curpjcount = 0
        curbclist = []
        curswlist = []
        currunslist = []
        curfirstlast = [0,0,0]
        for a in seqs:
            anew = a
            if(chopseqs=='right'):
                anew = a[:-1]
            elif(chopseqs == 'left'):
                anew = a[1:]
            elif(chopseqs == 'both'):
                anew = a[1:-1]
            #if(len(anew)>0):
            seqschop+=[anew]
            pct = anew.count("P")
            jct = anew.count("J")
            curbclist+=[[pct,jct]]
            curpcount+=pct
            curjcount+=jct
            pjct = anew.count("PJ")
            jpct = anew.count("JP")
            curswlist += [[pjct,jpct]]
            curpjcount+=pjct
            curjpcount+=jpct
            currunslist += [longestRun(a,"PJ")]
            if(len(anew)>1):
                if(anew[0]=="J"):
                    curfirstlast[0]+=1 #J in the first position
                if(anew[-1]=="J"):
                    curfirstlast[1]+=1 #J in the last position
                curfirstlast[2]+=1 #this one counts all seqs
        first_last.update({bc:tuple(curfirstlast)})
        run_lists.update({bc:currunslist})
        all_lists.update({bc:curbclist})
        switch_lists.update({bc:curswlist})

        pcount+=[curpcount]
        jcount+=[curjcount]
        jpcount +=[curjpcount]
        pjcount +=[curpjcount]
    return all_lists,run_lists,switch_lists,first_last
def goodMatch(editDistance1,editDistance2=1000,thresh=10):
    if((editDistance1 < editDistance2) and (editDistance1 < thresh) and (editDistance1 >= 0)):
        #print("readforward")
        return 0
        #seqstats[0]+=1
    elif((editDistance2 < editDistance1) and (editDistance2 < thresh) and (editDistance2 >= 0)):
        #print("readreverse")
        return 1
    else:
        return -1
        #print("dunno")

def goodMatchGeneralized(matchlist,thresh=10):
    """This function looks through all the edit distances
    presented in the list 'matchlist', and returns the smallest one, or -1
    if they are all -1 or all above thresh"""
    bestmatch = -1
    minimatch = thresh+10
    for match_index in range(len(matchlist)):
        a =matchlist[match_index]
        if(a < thresh and a >=0):
            #at least one of the given matches fits in our threshold!
            if(a<minimatch):
                bestmatch = match_index
                minimatch = a
    return bestmatch
#"""
def rawcount(filename):
    f = open(filename, 'rb')
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.raw.read

    buf = read_f(buf_size)
    while buf:
        lines += buf.count(b'\n')
        buf = read_f(buf_size)

    return lines
def FastqMultiFileIterator(fastq_files_directory):
    """this function iterates through all the fastq files in a directory,
    returning one fastq record at a time"""
    fastq_files = glob.glob(os.path.join(fastq_files_directory,"*.fastq"))
    print(len(fastq_files))
    for fastq_file in fastq_files:
        with open(fastq_file) as open_fastq_file:
            for title, seq, qual in FastqGeneralIterator(open_fastq_file):
                yield title,seq,qual
def findBarcodes(sequence,bcorder,buffer=12,beginstate=1,thresh = 0.6,debug=False):
    """this function will go through the sequence and identify the bits and the order of them
    bcorder is formed like this: {1:[[bc1],[2]],
                                  2:[[bc2,bc3],[3]],
                                  3:[[bc4],[4]],
                                  4:[[bc5,bc6],[3,5]],
                                  5:[[bc7],[6]],
                                  6:[[bc8,bc9],["x"]]}
    bc is [sequence,bcname]
    each dictionary entry contains what possible barcodes could be in that spot, followed by what
    the next entry could be."""
    #we start by looking through the entire sequence to find where bc1 is. Then, we take the next chunk
    # (which is defined by the max of length between bc2 and bc3) and see if it contains bc2 or bc3. if NOT,
    # then move on to the subsequent chunk, which could contain bc4. if that is also not found, then go on
    # and look for bc5 or bc6. after that, we will search for bc4 OR bc7. If either of those are found
    # that determines what the next search is. if nothing is found, you go on along both paths.
    # but what if the length of sequence is different per path?  
    lastpos = 0
    accumulatedskip = 0 #how much we've already skipped because we couldnt find anything
    seqidentified = [] #identified sequence elements based on given barcodes
    
    possiblestates = [beginstate]
    searched = {}
    firsttime = 1
    while lastpos < len(sequence) and len(possiblestates)>0:
        #curstate = beginstate
        scores = []
        maxskip = 0 #how much we can skip based on max length
        possiblenewstates = []
        if(debug):
            print(possiblestates)
        for stateid in possiblestates:
            
            try:
                state = bcorder[stateid]
            except KeyError:
                #either we got "x", or something is wrong
                if(stateid == "x"):
                    continue
                else:
                    raise KeyError("couldn't find a state with name "+str(stateid))
            possiblenewstates += state[1]
            #we go through each option and see which is the most likely
            nextstate = state[1] #this is the next set of possible states if 
                                 #this state turns out to be correct
            #each possibly barcode though might have a different length!!!
            #if we failed to find ANYTHING, we need to advance forwards
            #ok so in that case let's keep the previous sequence. in case
            #we find anything it will advance forwards anyway.
            #but each state has its own set of possible barcodes so we
            #need to pick the longest of those to skip forward.
            #for example:
            # [  bc1  ][  bc2  ]----unk-------------------
            #                   [ last max ][ current max]
            #                   [   search here          ]
            #
            
            maxlen = max([len(bc[0]) for bc in state[0]])
            if(maxlen>maxskip):
                #keep track of what is the longest
                maxskip = maxlen
            
            bufferlen = maxlen+buffer
            if(firsttime):
                bufferlen+=50000
                firsttime = 0
            rightpos = lastpos+accumulatedskip+bufferlen
            if(rightpos > len(sequence)):
                rightpos = len(sequence)
            leftpos = lastpos-buffer
            if(leftpos < 0):
                leftpos = 0
            
            seqchunk = sequence[leftpos:rightpos]
            done = 1
            for bc in state[0]:
                bcseq = bc[0]
                bcid = bc[1]
                #going through the list of possible barcodes
                #print("aligning "+str(bcid)+" at "+str(leftpos)+":"+str(rightpos))
                ""
                if(bcid in searched):
                    #this is checking if we already tried searching for this
                    checklpos = searched[bcid][0]
                    checkrpos = searched[bcid][1]
                    #there are many scenarios where we are re-searching through
                    #the same area. But we are hoping that we can find something
                    #we didn't before because it spans the junction.
                    #however, if we are trying to search entirely within an area we already
                    #searched, then stop.
                    if(leftpos>=checklpos and rightpos<=checkrpos):
                        continue
                #"""
                done = 0
                #if the above variable is set to 0, we haven't searched for all possibilities
                #yet
                alignment = edlib.align(bcseq,seqchunk, mode="HW", task="path",k=-1)
                searched[bcid]=[leftpos,rightpos]
                alignmentidentity = 1.0-alignment['editDistance']/float(len(bcseq))
                #print("identity is "+str(alignmentidentity))
                if(alignmentidentity >= thresh):
                    #this means we found an alignment that is better than the threshold!
                    #there could be many of them....
                    #in this case all we are storing is how well it aligned, which barcode it was,
                    #and what the new end position would be.
                    if(debug):
                        print("found {} matching {} in {}:{}".format(bcid,alignmentidentity,leftpos,rightpos))
                    scores += [[alignmentidentity,bcid,lastpos+alignment['locations'][0][1],nextstate]]
            if(done):
                #this means there's nothing else to look for
                possiblenewstates = possiblenewstates[:-1]
                

        #now we have looked through all the possibilities
        if(len(scores)>0):
            #print("found!")
            #this happens if we found something

            sortedscores = sorted(scores)
            bestscore = sortedscores[-1]
            lastpos = bestscore[2]
            seqidentified += [bestscore[1]]
            if(bestscore[3]=="x"):
                #then we've reached the end!!!!
                break
            possiblestates = bestscore[3]
            accumulatedskip = 0
        else:
            #in this case we didn't find anything.
            accumulatedskip+=maxskip #make the search space longer
            seqidentified += ["unidentified"] #track that we didn't find anything
            possiblestates = possiblenewstates #the possibilities advance
    return seqidentified

def barcodeSplitAndCountRecords(fastq_files_directory,barcodes,processreads=10000,barcode_detection_threshold = 7,\
                                minimum_slice_length = 220,step_length = 80,overlap_length=30,end_threshold=12,\
                                prefix_sequence="CCCAGCAGGTATGATCCTGACGACGGAGCACGCCGTCGTCGACAAGCC",\
                                prefix_detection_threshold=12,\
                               variable_sequences=["CTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC","TTTCAATTTAATCATCCGGCTCGTATAATGTGTGGA"],\
                               postfix_sequence="CAAGCCCATTATTACCCTGTTATCCCTAGACACCAATCAGAGGCCACA",\
                               variable_sequence_threshold=12,
                               progressbar = True,frontchecklength=150,visualize=True):
    allseqs = []
    allseqDict = {}
    unsorted = []
    #taking care of the progress bar, below
    if(progressbar):
        pbar2 = widgets.FloatProgress(
            min=0,
            max=1.0
        )
        display(pbar2)
    if(type(fastq_files_directory)==str):
        #in_handle = open(fastq_file_name)
        fastq_files_list = glob.glob(os.path.join(fastq_files_directory,"*.fastq"))
        fastq_iterator = FastqMultiFileIterator(fastq_files_directory)
        total_number_of_reads = (rawcount(fastq_files_list[0])/4)* len(fastq_files_list)
    else:
        in_handle = fastq_files_directory
        fastq_iterator = [["_",a,"_"] for a in in_handle]
        total_number_of_reads = len(fastq_files_directory)

    rcprefix_sequence = rc(prefix_sequence)
    rcpostfix_sequence = rc(postfix_sequence)
    if((processreads > total_number_of_reads) or (processreads <= 0)):
        processreads = total_number_of_reads
    how_many_so_far = 0
    seqstats = [0,0,0]
    reverse_complement_variable_seqs = [rc(a) for a in variable_sequences]


    for title, seq, qual in fastq_iterator:
        #we use a function from BioPython which streams the file
        how_many_so_far+=1
        sequence_is_reverse=False
        #update the progress bar every 100 sequences we read

        if(how_many_so_far %100 == 0):
            if(pbar2 != None):
                pbar2.value = float(how_many_so_far/float(processreads))

        #chop off the beginning and the end for orientation checking
        front = seq[:frontchecklength]
        rev = seq[-frontchecklength:] #make sure this one is reverse complemented!
        #now, we are checking for the presence of the attB site. That marks the "front!"
        #actually that is the plasmid side therefor the newest added spacer
        fread = edlib.align(prefix_sequence,front, mode="HW", task="path",k=-1)["editDistance"]
        rread = edlib.align(rcpostfix_sequence,front, mode="HW", task="path",k=-1)["editDistance"]
        #rread = edlib.align(rcprefix_sequence,rev, mode="HW", task="path",k=-1)["editDistance"]
        curseq = seq #we're going to use this variable to store the sequence which will
        #subsequently be analyzed for having our barcodes
        #this next part compares the quality of match between the front and the back.
        #whichever has the best match will be deemed correct. If none of them have matches
        #that are good enough, then just toss them.
        if(fread >= 0 and rread >= 0):
            if(fread < rread and fread < prefix_detection_threshold):
                #if this happens then the read was already forwards
                seqstats[0]+=1
            elif(rread < fread and rread < prefix_detection_threshold):
                #if this happens then the read was reverse, and we reverse it
                seqstats[1]+=1
                sequence_is_reverse=True
                #curseq = rc(seq)
            else:
                #both sides scored the same on the attb site, so then we don't know
                #which direction is more correct. This seems very unlikely
                seqstats[2]+=1
        else:
            #if this happens then both sides failed to detect the attb site, and we
            #also don't know which direction is right. This seems more likely
            seqstats[2]+=1
        this_reads_barcode = "none"
        #next part finds barcodes in the edge sequences
        for bc in barcodes:
            #so we go through the list of all barcodes and align them one by one to ONLY the
            #"front" sequence, since that one is usually higher quality. It is possible though
            #that the front won't align, and in that case we could use the back.
            #also if the sequence is backwards, then it could be that whatever's in the beginning
            #doesn't belong to the read we are going to start aligning
            #for example:
            #--------[BC1]->-read1->-[rcBC1] #this is the optimal scenario
            #--------[BC1]-<-read1-<-[rcBC1] #this is also OK
            #--------[BC1]-<-read1-<-[rcBC1][BC2]-<-read2-<-[rcBC2] #in this case we will read the
            #                                                      sequence from read2 but assume the
            #                                                      barcode from read1. This is bad
            #--------[BC1]->-read1->-[rcBC1][BC2]-<-read2-<-[rcBC2] #this will have the same problem as above
            #--------[BC1]-<-read1-<-[rcBC1][BC2]->-read2->-[rcBC2] #this won't read anything
            #--------[BC1]->-read1->-[rcBC1][BC2]->-read2->-[rcBC2] #this one will work fine
            #
            #we analyze read1 no matter what direction it goes in, then chop off the sequence before
            #we get to read2, and put read2 back through the pipeline.
            is_front_barcode = edlib.align(barcodes[bc][0],front, mode="HW", task="path",k=-1)["editDistance"]
            if(goodMatchGeneralized([is_front_barcode],thresh=barcode_detection_threshold) == 0):
                this_reads_barcode = bc
                break
            #is_rev_barcode
        #if(this_reads_barcode == "none"):
            #so this means we didn't detect a barcode in the "front" of the sequence.
            #We can check the rear!! Only problem is if it's a multimer then this would be
            #wrong.

        ""
        i=0
        simplified_sequence = []

        while(True):
            #iterate through the read chopping out bit by bit
            slicing = [step_length*i-overlap_length*(i>0),step_length+step_length*i] #this decides where to chop
            i+=1
            subseq = curseq[slicing[0]:slicing[1]] #we have chopped out a section of the sequence
            if(len(curseq)-slicing[0]< (sum([len(a) for a in variable_sequences])/float(len(variable_sequences)))):
                #if the sequence we've got is shorter than any of the variable sequences we are looking for,
                #don't even bother checking it, because it won't match
                break
            if(not sequence_is_reverse):
                #this part aligns variable regions to the chunk we have
                matchlist = [edlib.align(a,subseq, mode="HW", task="path",k=-1)["editDistance"]\
                                                                     for a in variable_sequences]
                #maybe instead of a variable sequence we will have found the sequence that goes at the end.
                isitTheEnd = edlib.align(postfix_sequence,subseq, mode="HW", task="path",k=-1)["editDistance"]
            else:
                #this does the same thing but in reverse
                matchlist = [edlib.align(a,subseq, mode="HW", task="path",k=-1)["editDistance"]\
                                                        for a in reverse_complement_variable_seqs]
                #if the sequence is reverse, then we look for the beginning, and that's the end.
                isitTheEnd = edlib.align(rcprefix_sequence,subseq, mode="HW", task="path",k=-1)["editDistance"]
                

            #isitJ = edlib.align(Jonly,subseq, mode="HW", task="path",k=-1)["editDistance"]
            #isitP = edlib.align(Ponly,subseq, mode="HW", task="path",k=-1)["editDistance"]

            if(goodMatchGeneralized([isitTheEnd],thresh=end_threshold)==0):
                #this is true if we found the end
                if(sequence_is_reverse):
                    simplified_sequence+=['B']
                else:
                    simplified_sequence+=['E']
                
                nslice = step_length*(i+1)-overlap_length*(i+1>0)
                seqslice = curseq[nslice:]

                if(visualize):
                    print('.'.join([str(a) for a in simplified_sequence]))

                if(len(seqslice)>minimum_slice_length):
                    #if we found the end sequence before reaching the end of the actual sequence, then
                    #that means the second half of the sequence could be another read, so put it in
                    #the unsorted list.
                    unsorted += [curseq[nslice:]]
                break
            else:
                #we haven't found the end but instead found one of the barcodes.
                which = goodMatchGeneralized(matchlist,thresh = variable_sequence_threshold)
                #i wonder if this will work if we find NO barcodes? which has to return something i guess
                simplified_sequence+=[which]
        if(sequence_is_reverse):
            simplified_sequence = simplified_sequence[::-1]
        if(this_reads_barcode in allseqDict):
            allseqDict[this_reads_barcode]=allseqDict[this_reads_barcode]+[tuple(simplified_sequence)]
        else:
            allseqDict[this_reads_barcode]=[tuple(simplified_sequence)]
        #"""
        #processreads-=1
        if(how_many_so_far >= processreads):
            break
    return allseqDict,seqstats,unsorted
def bc2dhist(data,density=True,minlength = 2,maxlength = 10):
    """create a 2d histogram given a set of data such as
    [[p in read 1, j in read 1], [p in read 2,j in read 2],...]"""
    x = []
    y = []
    for a in data:
        if(a[0]+a[1]<=maxlength and a[0]+a[1]>=minlength):
            x+=[a[0]]
            y+=[a[1]]
    #x = [a[0] for a in data]
    #y = [a[1] for a in data]

    bchist = np.histogram2d(x,y,bins=[range(0,maxlength),\
                                range(0,maxlength)],density = density)
    labels = [a for a in range(0,maxlength)]
    return bchist, labels
def makeBCplot1(bc,crange,sqrange,experiment_dataframe,datalist,\
                              labs = ["P barcodes", "J barcodes"],makefig=True):
    if(makefig):
        plt.figure()
    bins = [int(sqrange[1])-int(sqrange[0])+1,int(sqrange[1])-int(sqrange[0])+1]
    h1,_,_,_ = plt.hist2d([a[1] for a in datalist[bc]],[a[0] for a in datalist[bc]],\
               cmin=crange[0],cmax=crange[2],bins=bins,range=[sqrange,sqrange],\
               normed=True,cmap = "Reds")

    #cbar = plt.colorbar()
    #cbar.solids.set_edgecolor("face")
    #plt.draw()
    h2,_,_,_ = plt.hist2d([a[1] for a in datalist[bc]],[a[0] for a in datalist[bc]],\
               cmin=crange[1],cmax=crange[2],bins=bins,range=[sqrange,sqrange],\
               normed=True,cmap = "viridis")

    #cbar = plt.colorbar()
    #cbar.solids.set_edgecolor("face")
    #plt.draw()
    #plt.clf()
    #plt.figure()
    sns.heatmap(h1,cmap = "Reds",vmin = crange[0],vmax=crange[1],annot=True)
    sns.heatmap(h2,cmap="viridis",vmin=crange[1],vmax=crange[2],annot=True)
    #print(bc)
    #print(bcs)
    #print(conditions)
    condstr = condNameBC(bc,experiment_dataframe)
    plt.title(condstr)
    plt.xlabel(labs[0])
    plt.ylabel(labs[1])
    return h1,h2
def condNameBC(bc,experiment_dataframe):
    """condense the name of a set of conditions to one string"""
    bcline = experiment_dataframe[experiment_dataframe["barcode"]==bc].iloc[0].fillna("")
    return "_".join([str(a) for a in bcline])
def diffPlot(control_hist,experimental_hist,cmap="RdBu",color_range = [-0.05,0.05],\
             labs=["P barcodes", "J barcodes"],control_barcode=None,exp_barcode=None,\
             experiment_dataframe=None,makefig=True,annot=True):
    diff_hist1 = experimental_hist-control_hist
    if(makefig):
        plt.figure()

    sns.heatmap(diff_hist1,cmap = "RdBu",vmin = color_range[0],vmax=color_range[1],annot=annot)
    if(type(experiment_dataframe)==type(None)):
        plt.title("difference")
    else:
        ctrlname = condNameBC(control_barcode,experiment_dataframe)
        expname = condNameBC(exp_barcode,experiment_dataframe)
        plt.title("{} - {}".format(expname,ctrlname))
    plt.xlabel(labs[0])
    plt.ylabel(labs[1])
def quantifyRecOrder(seqs,seqlimits=(0,100),letters="AB",rev=False,stretch=False,chop='none'):
    order1 = []
    order2 = []
    maxlen = seqlimits[1]
    frac = np.zeros([maxlen,2])
    #print(frac)
    for aposeq in seqs:
        seq = []
        for el in aposeq:
            if(el in letters):
                seq+=[el]
        if(chop=='right'):
            seq = seq[:-1]
        elif(chop=='left'):
            seq = seq[1:]
        elif(chop=='both'):
            seq = seq[1:-1]
        if(len(seq)< seqlimits[0] or len(seq) > seqlimits[1]):
            continue
            #if the sequence length is out of our desired window, skip!
        if(rev):
            #this makes it read the sequences backwards
            seq = seq[::-1]
        for lett in range(len(seq)):
            #this goes letter by letter
            #print(seq)
            #print(letters)
            #print(seq[lett])
            lind = letters.index(seq[lett])
            #this part tells us which column to put our count into,
            #so for example if the sequence at position 3 is "A", then
            #we add a '1' to the "A" column, etc
            if(stretch):
                #in this case we are assuming that each sequence recorded
                #the same length of time, but at a random resolution depending
                #on the length.
                #each letter then corresponds to a larger or smaller region
                #of the frac array.This means we need to extrapolate!!
                boxwidth = maxlen/len(seq)
                start = boxwidth*lett
                end = boxwidth*(lett+1)
                #print("{},{}".format(start,end))
                #all the whole numbers between start and end get added to column 1 or 0
                #fractional columns get a fraction of a 1.
                #example: start = 2.5 end = 7.5
                #2 gets 0.5,
                #3-6 gets 1
                #7 gets 0.5
                firstbox = int(start)
                lastbox = int(end)
                #possible scenarios: first and last are the same: in this case add the
                #total difference (which should be less than one) to the first box
                if(firstbox == lastbox):
                    frac[firstbox][lind]+=end-start
                else:
                    #in this case we add a fraction to the first box, a whole number to intervening
                    #boxes, and a fraction to the last box.
                    # 0  1  2  3  4  5  6
                    #[ ][ ][ ][ ][ ][ ][ ]
                    # |        |
                    #0.8      3.8
                    #first = 0
                    #last = 3
                    #in-between = 1,2
                    #0 gets 1-0.8 = 0.2
                    #3 gets 0.8
                    #1,2 gets 1
                    if(1.0-(start-firstbox)>.001):
                        frac[firstbox][lind]+=1.0-(start-firstbox)
                    if((end-lastbox)>.001):
                        frac[lastbox][lind]+=(end-lastbox)
                    for box in range(firstbox,lastbox)[1:]:
                        frac[box][lind]+=1.0


            else:
                #in this case we are assuming each sequence represents
                #an accurate record starting at the same time and ending
                #at a random time corresponding to the length.

                if(len(frac)>lett):
                    #this makes sure we don't run off the end of the array.
                    #if the current sequence is longer than the 'frac' array,
                    #then just don't count any more after you've put everything
                    #into frac.
                    frac[lett][lind]+=1
        if(len(seq)>1):
            order1+=[seq.count(letters)]
            order2+=[seq.count(letters[::-1])]

    realfrac = []
    for a in frac:
        realfrac+=[a[0]/(a[0]+a[1])]
    return order1,order2,realfrac
def longestRun(string,chars):
    """counts the number of iterations in the longest run of each character"""
    mode = 0
    current_run = 0
    best_run = [0]*len(chars)
    for val in string:
        if(val.lower() == chars[mode].lower()):
            #same character as we've seen! great
            current_run+=1
        elif(val.lower() in (chars[:mode]+chars[mode+1:]).lower()):
            #now we hit a character of the other type
            if(current_run > best_run[mode]):
                best_run[mode]=current_run
            current_run = 1
            mode = chars.index(val)
    if(current_run > best_run[mode]):
        #counting the last run in the string
        best_run[mode]=current_run
    return(best_run)
def diffPlotWrapper(ctrlBC,expBC,experiment_dataframe,dlist,labs=["Ps","Js"],crange=[0,.02,1],\
						sqrange=[-.5,6.5],annot=True,color_range=[-.05,.05]):
    plt.figure()
    control_hist1,control_hist2 = makeBCplot1(ctrlBC,crange,sqrange,\
                                    experiment_dataframe,dlist,labs,makefig=False)
    experiment_hist1,experiment_hist2 = makeBCplot1(expBC,crange,sqrange,\
                                    experiment_dataframe,dlist,labs,makefig=False)
    #print(experiment_hist1)
    #plot the raw data from each experiment, then REMOVE those plots, and plot`
    #the difference plot. This is because that's the only way to get all that data
    #in the right orientation.
    plt.clf()
    diffPlot(control_hist1,experiment_hist1,labs=labs,control_barcode=ctrlBC,\
             exp_barcode=expBC,experiment_dataframe=experiment_dataframe,annot=annot,\
			 color_range=color_range)
