import glob,sys
import numpy as np
import matplotlib.pyplot as plt

def hrr2summary(hhr="sys.argv[1]",sequence="seq",divide=500,prob_th=30,ymax=30):
    total=np.zeros(len(seq))
    base=int(hhr.split(".")[0].split("_")[1])/int(divide)-1 # hhr name is like  "hhpred_2000.hhr"
    print("base:%s:%s-" % (base,int(base)*divide))
    fin=open(hhr,"r")
    flag=0
    counter=1
    header="No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM".split()
    for x in fin:
        if x.find("Prob ")>=0:
            flag+=1
        elif flag==1:
            try:
                #print(x.strip().split()[0],counter)
                if len(x.strip().split())>11 and int(x.strip().split()[0])==counter:
                    counter+=1
                    probability=float(x[35:40])
                    print("Prob:%s" % probability)
                    tmp=np.zeros(len(seq))
                    (queryHMMstart,queryHMMend)=tuple(x.strip().split()[-3].split("-"))
                    rebase=base*divide
                    print("rebase:%s" % rebase)
                    print(int(queryHMMstart)+rebase,int(queryHMMend)+rebase)
                    tmp[int(queryHMMstart)+int(rebase)-1:int(queryHMMend)+int(rebase)-1]=1 # Get pdb hit range in sequence
                    print("len:%s,sum:%s" % (len(tmp),np.sum(tmp,axis=0)))
                    if probability>prob_th:
                        total+=tmp
                elif flag==1 and x.find("Probab")>=0:# For correlation of Probab and other parameters 
                    print(x.strip())
                    Probab=x.strip().split()[0].split("=")[1]
                    Similarity=x.strip().split()[5].split("=")[1]
                    Evalue=x.strip().split()[1].split("=")[1]
                    print(Probab,Similarity,Evalue)
                    prob_simi_evalue_results.append((Probab,Similarity,Evalue))
                
            except:
                    continue


            
    fin.close()
    return(total,prob_simi_evalue_results)
    
def read_fasta(seq="sys.argv[1]"):
    fin=open(seq,"r")
    seqall=[]
    for x in fin:
        seqall.append(x.strip())
    fin.close()
    seq=[x for x in seqall if not x.find(">")==0 ] # Remove header
    seq=("").join(seq)
    lenres=len(("").join(seq))
    return seq

if __name__=="__main__":
    hhrs=glob.glob("*hhr")
    seq=read_fasta(seq=sys.argv[1])
    print(seq)
    alltotal=np.zeros(len(seq)) # Similar PDB frequency mapped to sequence.
    divide=sys.argv[2] # Divided range length of original fasta to HHpredict.
    prob_th=sys.argv[3]
    ymax=int(sys.argv[4])
    prob_simi_evalue_results=[]
    for hhr in hhrs:
        print("############")
        print(hhr)
        seq=read_fasta(seq=sys.argv[1])
        (total,simi_evalue_results)=hrr2summary(hhr=hhr,sequence=seq,divide=int(divide),prob_th=int(prob_th),ymax=int(ymax))
        alltotal+=total
    #print(total,np.sum(total,axis=0))
    print(alltotal,np.sum(alltotal,axis=0))
    ########
    # Plot
    ########
    # check residue colums have values
    #for (i,x) in enumerate(alltotal):
    #    print(i,x)
    left=range(0,len(seq))
    height=alltotal
    # Bar plot resi vi PDB frequence at Probability
    plt.bar(left,height,alpha=0.8)
    xticks=range(1,int(len(seq)),int(divide))
    plt.xticks(xticks)
    plt.ylim([0,ymax])
    plt.title("Probability>%s" % sys.argv[3])
    #plt.show()
    plt.savefig("Probability>%s.png" % sys.argv[3])
    plt.cla()
    prob=[float(p) for (p,s,e) in  prob_simi_evalue_results]
    print(prob)
    simi=[float(s) for (p,s,e) in  prob_simi_evalue_results]
    escore=[float(e) for (p,s,e) in  prob_simi_evalue_results]
    #plt.scatter(prob,simi)
    #plt.xticks(xxx)
    #plt.ylim([0,1])
    #plt.xlim([0,99])
    #plt.cla()

    # Scatter plot Probability vs Escore
    plt.scatter(prob,escore)
    plt.savefig("Prob_Escore.png")
    plt.cla()
    # Scatter plot Esocre<1
    plt.scatter(prob,escore)
    plt.ylim([0,1])
    plt.savefig("Prob_Escore_max1.png")
    plt.cla()
    #########
    # usage #
    #########
    #python3.6 extract_modify_hhr2summary.py Q9NR09.fasta 500 prob_th=0 ymax=50
