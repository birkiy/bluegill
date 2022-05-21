import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D
import seaborn as sns

import numpy as np



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

    plt.rcParams["figure.dpi"] = dpi

    sns.set(font_scale=1.5, style="ticks")

    fig = plt.figure(
        figsize=[ncols*3.5, nrows*3.5]
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


