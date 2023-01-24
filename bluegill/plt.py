
from matplotlib import gridspec
from matplotlib.lines import Line2D

import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pybedtools as bt

import re
import os

def plotHeatmaps(
    N,BED,
    samples, palette, sets, colorPalette,
    ylim=(0,100),vmin=0,vmax=5,vmins=None, vmaxes=None,
    ylabrot=0, h=3000, clab="Signal",clabs=None,
    dpi=50, interpolation="antialiased",
    noSort=False, profile=True
):
    """
    Plots profiles and heatmaps of signal of genomic regions. 
    
    :param N: Numpy Signal matrix NxSxB, where number of regions is N, number of signal files is S, number of bins is B.
    
    :param BED: Pandas DataFrame of regions with the same order of N. Different sets are defined in `Set` column.
    
    :param samples: Names of signals. Dimention of S.
    
    :param palette: Colors of heatmaps. Default `Blues`
    
    :param sets: According to BED, which sets will be plotted.
        
    :param colorPalette: Colors of sets for profile.
    
    :param ylim: Ylim of profiles. Default `100`
    
    :param mmax: Multiplication coeficient to define the vmax of heatmap. Default `2`
    :param pow_: Power coeficient to define the vmax of heatmap. Default `0.5`
    
    :param ylabrot: Rotation angle of ylabels. Default `0`
    
    :param h: Range from center (BP). Default `3000`
    
    :param clab: Label of of colorbar. Default `TMM Signal`
    
    :param dpi: Dpi of figure. Default `50`
    
    :param interpolation: Normalization of heatmaps. Default `antialiased`
    
    :param noSort: Not sorting regions according to mean signal. Default `False`
    
    
    """
    
    if clabs is None:
        clabs = [clab]*len(samples)
    
    if vmaxes is None:
        vmaxes = [vmax]*len(samples)

    if vmins is None:
        vmins = [vmin]*len(samples)
    
    if palette is None:
        palette = len(samples) * ["Blues"]
    
    if "Set" not in BED.columns:
        BED["Set"] = "regions"
    

    if noSort:
        Nsorted = N
        sortedBED = BED
    else:
        orderRegion = np.argsort(N.mean((1,2)))[::-1]
        sortedBED = BED.iloc[orderRegion,:].reset_index()
        Nsorted = N[orderRegion, :,:]


    idxs = [sortedBED["Set"] == sets[i] for i in range(len(sets))]

    ratio = [sortedBED[idx].shape[0] for idx in idxs]

    if profile:
        ratio = [sum(ratio) // 4, *ratio]

    Nbins = N.shape[-1]
    
    fig = plt.figure(
        figsize=[len(samples)*3.5, 14], dpi=dpi
    )
    gs = gridspec.GridSpec(
        nrows=len(ratio), 
        ncols=len(samples),
        height_ratios=ratio
    )
    plt.subplots_adjust(
        hspace=0.05,
        wspace=0.3)

    
    for i, sample in enumerate(samples):
        if profile:
            j_plus = 1
            fig.add_subplot(gs[0,i])
            for j,idx in enumerate(idxs):
                NS = Nsorted[idx,i,:]
                plt.plot(NS.mean(0), c=colorPalette[sets[j]])
                plt.ylim(ylim)
                if i != 0:
                    plt.yticks([])
                else:
                    plt.ylabel(clab)
                plt.xticks([])


                if i == len(samples)-1:

                    handles = [  
                        Line2D([0], [0], marker='o', color=colorPalette[k], label=k, markersize=7, linestyle="None")
                        for k in colorPalette if k in sets
                    ]
                    labels = [k for k in colorPalette if k in sets]
                    plt.legend(handles, labels, frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))
            plt.title(sample)
        else:
            j_plus = 0

        for j,idx in enumerate(idxs):
            NS = Nsorted[idx,i,:]
            ax = fig.add_subplot(gs[j+j_plus,i])

            # interpolation = "None" if you need no normalization
            plt.imshow(NS,aspect="auto", cmap=palette[i], vmax=vmaxes[i], vmin=vmins[i] ,  interpolation=interpolation)
            if not profile and (j == 0):
                plt.title(sample)

            plt.xticks([]) 
            plt.yticks([])

            if i == 0:
                plt.ylabel(sets[j], rotation=ylabrot, ha="right", va="center")
            else:
                plt.ylabel("")

            if j == len(idxs)-1:
                
                plt.xticks([0, Nbins/2, Nbins], [f"-{h//1000}kb", "Center", f"+{h//1000}kb"], rotation=0)
                cax = fig.add_axes([ax.get_position().x0,ax.get_position().y0-0.05, ax.get_position().x1-ax.get_position().x0,0.01])
                cbar = plt.colorbar(cax=cax, orientation="horizontal")

                #cbar = plt.colorbar(location="bottom", pad=0.15)
                cbar.set_label(clabs[i]) 
            else:
                plt.xticks([])

            plt.yticks([])

            
    return  fig





def plotProfiles(
    N, BED,
    sets,colorPalette,
    nrows,ncols,
    ylim=15,ylab="Signal",
    dpi=100,
    h=3000
):

    idxs = [BED["Set"] == sets[i] for i in range(len(sets))]
    Nbins=N.shape[-1]

    
    fig = plt.figure(
        figsize=[ncols*3.5, nrows*3.5], dpi=dpi
    )
    
    gs = gridspec.GridSpec(
        nrows=nrows,ncols=ncols,
        hspace=0.4,wspace=0.2
    )




    for i,idx in enumerate(idxs):
        fig.add_subplot(gs[i])
        NS = N[idx,0,:]
        plt.plot(NS.mean(0), c=colorPalette[sets[i]],lw=2)
        plt.ylim([0,ylim])
        if i % ncols != 0:
            plt.yticks([])
        else:
            plt.ylabel(ylab)


        plt.xticks([])


        if i == (ncols*nrows)-1:

            handles = [  
                Line2D([0], [0], marker='o', color=colorPalette[k], label=k, markersize=7, linestyle="None")
                for k in colorPalette if k in sets
            ]
            labels = [k for k in colorPalette if k in sets]
            plt.legend(handles, labels, frameon=False, loc='center left', bbox_to_anchor=(1, nrows*0.5))  

        plt.title(sets[i], fontsize=25)

        plt.xticks([0, Nbins/2, Nbins], [f"-{h//1000}kb", "Center", f"+{h//1000}kb"], rotation=0)

    return fig


