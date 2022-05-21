
import os
import pickle
import multiprocessing 
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig as BW






def getIndex(l, nP):
    """
    This function splits list of regions into smaller chunks.
    """
    total = len(l)

    ppr = np.floor( total / nP)

    currents = []
    targets = []
    for pj in range(nP):
        if (pj+1) * ppr + ppr > total:
            currents.append(int(targets[-1]))
            targets.append(total)
            break
        currents.append(int(pj * ppr))
        targets.append(int((pj+1) * ppr))
    return currents, targets



def getSignal(poss, files, mi, Nbins,h, type_="mean", scaled=False):
    """
    This function gets signal.
    """
    S = np.zeros((len(poss), len(files), Nbins))
    for j,file in enumerate(files):
        print(f"{mi}: {file.split('/')[-1]}")
        bw = BW.open(file)
        for i,pos in enumerate(poss):
            if not scaled:

                center = (pos[1]+pos[2]) // 2
                if center-h < 0:
                    start = 1
                    startCrop = (Nbins//2) - (abs(center - start) // 20)
                else:
                    start = center -h
                    startCrop = 0
                if center +h > sizes.loc[pos[0], "Size"]:
                    end = sizes.loc[pos[0], "Size"]
                    endCrop = (Nbins//2) - (abs(center - end) // 20)
                else:
                    end = center +h
                    endCrop = 0
                bins = startCrop + endCrop
                
                try:
                    tmp = bw.stats(pos[0], start, end, nBins=Nbins-bins,  type=type_)
                except: 
                    print(pos)
            else:
                start = pos[1]
                end = pos[2]
                                
                startCrop = 0
                endCrop = 0
                
                hbin = Nbins // 4
                bins = startCrop + endCrop + 2*hbin
                
                tmp = bw.stats(pos[0], start, end, nBins=Nbins-bins,  type=type_)
                
                hbin1 = bw.stats(pos[0], start-h, start, nBins=hbin,  type=type_)
                hbin2 = bw.stats(pos[0], end, end+h, nBins=hbin,  type=type_)
                
                tmp = np.concatenate((hbin1, tmp,hbin2))            
            
            S[i,j,:] = np.nan_to_num(np.concatenate((np.zeros((startCrop)), tmp, np.zeros((endCrop)))))
    pickle.dump(S, open(f".tmp/{mi}.p", "wb"))




def readFile(path):
    poss = []
    with open(path, "r") as f:
        for line in f.readlines():
            row = line.strip().split("\t")
            poss.append((row[0], int(row[1]), int(row[2])))
    return poss




def concatSignal(out, nP):
    """
    This function concatenates multiple pickle files into one pickle file.
    """
    for i in range(nP):
        with open(f".tmp/{i}.p","rb") as f:
            tmp = pickle.load(f)
            if i == 0:
                A = tmp
            else:                
                A = np.concatenate((A, tmp), 0)
    print("Writing...")
    
    A = np.nan_to_num(A)
    with open(out, "wb") as f:
        pickle.dump(A, f)
        

def runSignal(
    BED, BWS, OUT,
    ref="hg19",
    scaled=False,Nbins=200,h=3000,
    type_="mean",
    nP=32
):
    """
    This function collects signal from BW files over multiple genomic regions using multiprocessing. According to genome build it fills extremes with 0.
    
    :param BED: Pandas DataFrame of regions.
    
    :param BWS: List of BW file locations.
    
    :param OUT: output location of the data.
    
    :param ref: Reference genome build. Default `hg19`
    
    :param scaled: If scaled, region of interest scaled to same size. Default `False`
    
    :param Nbins: Number of bins to cover the regions. Default `200`
    
    :param h: Range from center or TSS/TTS (BP). Default `3000` 
    
    :param type_: Type of collection of signal from BW. See pybigwig for details. Default `mean`
    
    :param nP: Number of processors. Default `32`
    """
    global sizes
    sizes = pd.read_table(f".ref/{ref}.chrom.sizes", names=["Chr", "Size"])
    print(ref)
    sizes = sizes.reset_index(drop=True).set_index("Chr")
    
    if type(BED) is list:
        poss = []
        BEDdf = pd.DataFrame()
        for path in BED:
            tmp = readFile(path)
            poss += tmp
            
            tmp = pd.DataFrame(tmp, names=["Chr", "Start", "End"])
            tmp["Set"] = Path(path).stem
            
            BEDdf = pd.concat([BEDdf, tmp])
            
        BED = BEDdf.sort_values(["Chr", "Start"])
        
        BEDdf = None
        
    else:
        poss = list(zip(BED["Chr"], BED["Start"], BED["End"]))
    
    samples = []
    for path in BWS:
        samples.append(Path(path).stem)
    
    
    currents, targets = getIndex(poss, nP)
    if not os.path.isdir(".tmp"):
        os.mkdir(".tmp")
        
    allprocesses = [
        multiprocessing.Process(
            target=getSignal,
            args=(
                poss[currents[i]:targets[i]],
                BWS,
                i,
                Nbins, h,
                type_,
                scaled
            )
        )
        for i in range(nP)
    ]
    for process in allprocesses:
        process.start()

    for process in allprocesses:
        process.join()

    concatSignal(OUT, nP)
    
    return BED, samples



 
    
        
    
    
