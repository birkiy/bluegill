
import os
import pickle
import multiprocessing 

import numpy as np
import pandas as pd
import pyBigWig as BW

                                                                                  


def getIndex(l, nP):
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

            else:
                start = pos[1]-h
                end = pos[2]+h
                startCrop = 0
                endCrop = 0
                bins = startCrop + endCrop

            try:
                tmp = bw.stats(pos[0], start, end, nBins=Nbins-bins,  type=type_)
            except:
                print(pos)

            S[i,j,:] = np.nan_to_num(np.concatenate((np.zeros((startCrop)), tmp, np.zeros((endCrop)))))
    pickle.dump(S, open(f".tmp/{mi}.p", "wb"))


                


def getSignal(poss, files, mi, Nbins,h, type_="mean", scaled=False):
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
                
                tmp = bw.stats(pos[0], start, end, nBins=Nbins-bins,  type=type_)
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









def concatSignal(out, nP):
    for i in range(nP):
        with open(f".tmp/{i}.p","rb") as f:
            tmp = pickle.load(f)
            if i == 0:
                A = tmp
            else:                
                A = np.concatenate((A, tmp), 0)
    print("Writing...")
    with open(out, "wb") as f:
        pickle.dump(A, f)
        

def runSignal(BED, BWS, OUT,
              ref="hg19", nP=32, Nbins=200,h=3000, type_="mean", scaled=False):
	global sizes
	sizes = pd.read_table(f"/groups/lackgrp/genomeAnnotations/{ref}/{ref}.chrom.sizes", names=["Chr", "Size"])
	print(ref)
	sizes = sizes.reset_index(drop=True).set_index("Chr")
	poss = list(zip(BED["Chr"], BED["Start"], BED["End"]))
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



