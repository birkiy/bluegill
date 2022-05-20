

def plotHeatmaps(samples, palette,
                 Nsorted,sets, ratio, idxs, ylabrot=0, ylim=100):
    
    plt.rcParams["figure.dpi"] = 50

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
            ax = plt.plot(NS.mean(0), c=colorPalette[sets[j]])
            #np.ceil(Nsorted.mean(0).max() / 100) *150
            plt.ylim([0,ylim])
            plt.yticks([])
            plt.xticks([])

            if i == len(samples)-1:

                handles = [  
                    Line2D([0], [0], marker='o', color=colorPalette[k], label=k, markersize=7, linestyle="None")
                    for k in colorPalette if k in sets
                ]
                labels = [k for k in colorPalette if k in sets]
                plt.legend(handles, labels, frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))  

        plt.title(sample)


        for j,idx in enumerate(idxs):
            NS = Nsorted[idx,i,:]
            fig.add_subplot(gs[j+1,i])
            plt.imshow(NS,aspect="auto", cmap=palette[i], vmax=2*(N.max() ** 0.5), vmin=0 )


            plt.xticks([]) 
            plt.yticks([])

            if i == 0:
                plt.ylabel(sets[j], rotation=ylabrot, ha="right")
            else:
                plt.ylabel("")

            if j == len(idxs)-1:
                plt.xticks([0, 100, 200], ["-3kb", "Center", "+3kb"], rotation=0)
            else:
                plt.xticks([])

            plt.yticks([])

            
    return  fig



    