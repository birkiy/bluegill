
from ._utils.run import *


def runBed(BED):
    """
    This function organizes multiple BED files into one pandas DataFrame.
    
    :param BED: Has to be a dict keys are set names, and values are file locations.
    """
    poss = []
    BEDdf = pd.DataFrame()
    for set_, path in BED.items():
        tmp = readFile(path)
        poss += tmp

        tmp = pd.DataFrame(tmp, columns=["Chr", "Start", "End"])
        tmp["Set"] = set_

        BEDdf = pd.concat([BEDdf, tmp])

    BED = BEDdf.sort_values(["Chr", "Start"])
    
    return BED

def runSignal(
    BED, BWS, OUT,
    scaled=False,igv=False,Nbins=200,h=3000,
    type_="mean",
    nP=32
):
    """
    This function collects signal from BW files over multiple genomic regions using multiprocessing. According to genome build it fills extremes with 0.
    
    :param BED: Pandas DataFrame of regions.
    
    :param BWS: List of BW file locations.
    
    :param OUT: output location of the data.
    
    :param scaled: If scaled, region of interest scaled to same size. Default `False`
    
    :param Nbins: Number of bins to cover the regions. Default `200`
    
    :param h: Range from center or TSS/TTS (BP). Default `3000` 
    
    :param type_: Type of collection of signal from BW. See pybigwig for details. Default `mean`
    
    :param nP: Number of processors. Default `32`
    """
    if type(BED) is dict:
        poss = []
        BEDdf = pd.DataFrame()
        for set_, path in BED.items():
            tmp = readFile(path)
            poss += tmp
            
            tmp = pd.DataFrame(tmp, columns=["Chr", "Start", "End"])
            #tmp["Set"] = Path(path).stem
            tmp["Set"] = set_
            
            BEDdf = pd.concat([BEDdf, tmp])
            
        BED = BEDdf.sort_values(["Chr", "Start"])
        
        BEDdf = None
        
        poss = list(zip(BED["Chr"], BED["Start"], BED["End"]))
        
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
                scaled,
                igv
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



def runTMM(S):
    S = np.nan_to_num(S)
    R = 1 /(conorm.tmm_norm_factors(S.mean(2)) * S.mean(2).sum(0) / 1000000)
    
    return S * R[:,None]
    
        
    
    
