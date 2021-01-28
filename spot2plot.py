import glob,sys
import numpy as np
import matplotlib.pyplot as plt

def readspot(seq,spotd):
    fin=open(spotd,"r")
    resi_resn_prob_label=[]
    for x in fin:
        if not x.find("#")==0:
            (resi,resn,prob,label)=x.strip().split()
            print(resi,resn,prob,label)
            resi_resn_prob_label.append((resi,resn,prob,label))
    fin.close()
    return resi_resn_prob_label


if __name__=="__main__":
    seq=sys.argv[1]
    spotd=sys.argv[2] 
    prob_th=sys.argv[3]
    ymax=int(sys.argv[4])
    resi_resn_prob_label=readspot(seq,spotd)
    ########
    # Plot #
    ########
    fig=plt.figure()
    ax=fig.add_subplot(111)
    left=range(0,len(resi_resn_prob_label)) # xaxis
    height=[float(prob) for (resi,resn,prob,label) in  resi_resn_prob_label] # yaxis
    plt.bar(left,height,alpha=0.8,color="red")
    xticks=range(1,int(len(left)),int(500))
    plt.xticks(xticks)
    ax.yaxis.tick_right()
    plt.ylim([float(prob_th),ymax]) # 0.46 is criteria of disorder. 
    plt.title("")
    #plt.show()
    plt.savefig("Disorder_prediction.png",transparent=True)

    #########
    # usage #
    #########
    #python3.6 spot2plot.py Q9NR09.fasta s0.spotd "0.46" 1
    #SPOT-Disorder2: Protein disorder prediction(https://sparks-lab.org/server/spot-disorder2/)
