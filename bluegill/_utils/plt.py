
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

    
    
    
