

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec
import re
import os
import numpy as np
import pybedtools as bt


def plotTranscript(transcript,name,s,y=0):
    
    for i, (type_, pos) in enumerate(transcript):
        if type_ == "exon":
            lw = 5
        elif type_ in ("5UTR", "3UTR"):
            lw = 5
        elif type_ == "CDS":
            lw = 8
        else:
            lw = 1

        plt.plot([pos[1] ,pos[2]] ,[y,y], lw=lw, c="#888888")
    
    plt.ylim([-1,y+1])

    

def plotTranscripts(window, s):
    ys = []
    ls = []
        
    for i, (transcript, gtf_) in enumerate(window.items()):
        plotTranscript(gtf_,transcript,s,y=i*2)
        ys.append(i*2) 
        ls.append(transcript)
    plt.yticks(ys, ls,color="#333333", weight="demibold", ha="right", va="center")

    
    
    


                
                

def plotTacks(BED,N,
              outdir,gtf,
              samples,palette, windows=None):
    
    
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    
    
    if windows is None:
        print("Windows: Started!")
        windows = {}
        for i, row in BED.iterrows():
            Chr, Start,End, Name = row[0], row[1], row[2], row[3], 

            Start = int(Start)
            End = int(End)

            l = End - Start

            Start -= l
            End += l

            trns = Name.split(";")[1]

            print(Name, end=" :: ")
            winname = f"{Chr}:{Start}-{End}.{Name}"
            windows[winname] = {}

            
            if Start >= 0:

                gtf_ = bt.BedTool(gtf)

                bed = bt.BedTool(f"{Chr} {Start} {End}", from_string=True)

                gtf_ = bt.BedTool.intersect(gtf_,bed).to_dataframe()  


                for i, row in gtf_.iterrows():

                    Chr_, type_, Start_, End_, Strand_, meta = row[0], row[2], row[3], row[4], row[6], row[8]

                    genename = meta[meta.index("gene_name"):].split("\"")[1]
                    trnsname = meta[meta.index("transcript_id"):].split("\"")[1]
                    if  trnsname.find(trns) == -1:
                        continue
            
                    name = f"{genename}.{trnsname}"

                    if name not in windows[winname]:
                        windows[winname][name] = [(type_, (Chr_, int(Start_), int(End_), Strand_))]
                    else:
                        windows[winname][name].append((type_, (Chr_, int(Start_), int(End_), Strand_)))
                        
                        
        print("Windows: Done!")



        
    
    groups = {}
    for i, member in enumerate(samples):

        signal = member.split("-")[0]
        groups[i] = [
            i
            for i, member in enumerate(samples)
            if member.find(signal) != -1
        ] 




    print("Plot Windows: Started!")
    for ex, winname in enumerate(windows):

        Chr, s, e, name = re.split('[:,\-,.]',winname, 3)

        fig = plt.figure(figsize=[20,20+(0.3*len(windows[winname]))])
        gs = gridspec.GridSpec(
            nrows=len(samples)+1, ncols=1,
            height_ratios=[1 for _ in range(len(samples))] + [0.3*len(windows[winname])],
            top=0.95)


        Nbins = N.shape[-1]
        xs = np.linspace(int(s), int(e),Nbins)

        for i in range(len(samples)):
            fig.add_subplot(gs[i])

            plt.plot(xs, N[ex,i,:], c = palette[i])
            plt.xticks([])


            group = groups[i]

            plt.ylim([0,N[ex,group,:].max()])
            plt.yticks([N[ex,group,:].max(), N[ex,group,:].max()])

            plt.ylabel(samples[i], rotation=0, fontsize=20, ha="right")

            plt.xlim(xs[0], xs[-1])

        fig.add_subplot(gs[len(samples)])


        plotTranscripts(windows[winname],(s,e))
        sns.despine(top=True, right=True, left=True, bottom=True)

        fig.suptitle(f"{winname} (Width = {(int(e)-int(s)) // 1000} kb)",fontsize=30)

        plt.xlim(xs[0], xs[-1])

        fig.savefig(f"{outdir}/{winname.replace('/','_')}.pdf", bbox_inches="tight", pad_inches=1);
        
    print("Plot Windows: Done!")
    
    return windows

