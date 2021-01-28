import sys

def divide_seq(seq="sys.argv[1]",resnumber=1000):
    fin=open(seq,"r")
    seqall=[]
    for x in fin:
        seqall.append(x.strip())
    fin.close()

    seq=[x for x in seqall if not x.find(">")==0 ]
    print(("").join(seq),len(("").join(seq)))
    seq=("").join(seq)
    lenres=len(("").join(seq))
    divide=int(lenres)/int(resnumber)

    for r in range(0,int(divide)+1):
        print(r)
        fout=open("%s.txt" % ((r+1)*int(resnumber)),"w")
        for (i,x) in enumerate(seq):
            #print(i,(r)*resnumber,(r+1)*resnumber)
            if i<(r+1)*int(resnumber) and i>=r*int(resnumber):
                print(i,(r)*resnumber,(r+1)*resnumber)
                fout.write("%s" % (x))
        fout.close()
    
    print(seq)
            
if __name__=="__main__":
    divide_seq(seq=sys.argv[1],resnumber=sys.argv[2])
    #python3.6 devide_seq.py Q9NR09.fasta 1000
    #python3.6 devide_seq.py Q9NR09.fasta 500

