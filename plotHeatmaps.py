import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D
import seaborn as sns


def plotHeatmaps(samples, palette, colorPalette,
                 Nsorted,sets, ratio, idxs, ylabrot=0, ylim=100, h="3kb", legend=True, mmax=2,pow_=0.5, Nbins=200, dpi=50, interpolation="antialiased"):
    
    plt.rcParams["figure.dpi"] = dpi

    sns.set(font_scale=1.5, style="ticks")
    
    
    fig = plt.figure(
        figsize=[len(samples)*3.5, 14]
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
        fig.add_subplot(gs[0,i])
        for j,idx in enumerate(idxs):
            NS = Nsorted[idx,i,:]
            plt.plot(NS.mean(0), c=colorPalette[sets[j]])
            plt.ylim([0,ylim])
            if i != 0:
                plt.yticks([])
            else:
                plt.ylabel("TMM Signal")
            plt.xticks([])


            if i == len(samples)-1 and legend:

                handles = [  
                    Line2D([0], [0], marker='o', color=colorPalette[k], label=k, markersize=7, linestyle="None")
                    for k in colorPalette if k in sets
                ]
                labels = [k for k in colorPalette if k in sets]
                plt.legend(handles, labels, frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))  

        plt.title(sample)


        for j,idx in enumerate(idxs):
            NS = Nsorted[idx,i,:]
            ax = fig.add_subplot(gs[j+1,i])
            # interpolation = "None" if you need no normalization
            plt.imshow(NS,aspect="auto", cmap=palette[i], vmax=mmax*(Nsorted.max() ** pow_), vmin=0 ,  interpolation=interpolation)


            plt.xticks([]) 
            plt.yticks([])

            if i == 0:
                plt.ylabel(sets[j], rotation=ylabrot, ha="right")
            else:
                plt.ylabel("")

            if j == len(idxs)-1:
                plt.xticks([0, Nbins/2, Nbins], [f"-{h}", "Center", f"+{h}"], rotation=0)
                cax = fig.add_axes([ax.get_position().x0,ax.get_position().y0-0.05, ax.get_position().x1-ax.get_position().x0,0.01])
                cbar = plt.colorbar(cax=cax, orientation="horizontal")

                #cbar = plt.colorbar(location="bottom", pad=0.15)
                cbar.set_label("TMM Signal") 
            else:
                plt.xticks([])

            plt.yticks([])

            
    return  fig

