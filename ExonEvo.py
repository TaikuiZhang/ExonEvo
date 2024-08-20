#!/usr/bin/env python3
import os,sys,re
import collections
import tempfile
import argparse
from EVOMODULE.Estimation import GENEUTM,GENEGFF,ARRSIM,EWIDTH,ARRIMG,MAPGENEG,MAPGENEYY,MAPS,LOOPMAP,Summary,Class_NR

try:
    import blosum as bl
except ImportError:
    os.system('pip3 install blosum')
    import blosum as bl
try:
    import numpy as np
except ImportError:
    os.system('pip3 install numpy')
    import numpy as np
try:
    from Bio import SeqIO
except ImportError:
    os.system('pip3 install Bio')
    from Bio import SeqIO

try:
    from Bio import Phylo
except ImportError:
    os.system('pip3 install Bio')
    from Bio import Phylo

try:
    from itertools import islice
except ImportError:
    os.system('pip3 install itertools')
    from itertools import islice

def handle():
    return None

def numpyNR(listB):
    NB=np.unique(listB)
    ND=[]
    for i in NB:
       ND.append(i)
    return ND
if __name__ == "__main__":
    print ("##ExonEvo v2023-Sep-28, a software for detecting conserved exons                                ##")
    print ("##developed by Taikui Zhang, PhD                                                                ##")
    print ("##Usage: ExonEvo.py --og HOG --gff HOG.GFF --pmsa HOG.fa --gtree HOG.tree --stree Ref.tree --csn 3 ##")
    parser = argparse.ArgumentParser()
parser.add_argument("--og", help="HOG0008435", required=True)
parser.add_argument("--gff", help="exon annotation", required=True)
parser.add_argument("--pmsa", help="protein sequence alignment", required=True)
parser.add_argument("--gtree", help="gene.tree", required=True)
parser.add_argument("--stree", help="Angiosperm.tree", required=True)
parser.add_argument("--csn",help="Minimum of continus sites",default=3,type=int)
args = parser.parse_args()

HOG=args.og
cmd1="mkdir "+str(HOG)
os.system(cmd1)
genetree=args.gtree
SpeciesTree=args.stree
Annofile=args.gff
PASTA=args.pmsa
MCSN=args.csn
tree = Phylo.read(genetree, "newick")
reftree = Phylo.read(SpeciesTree, "newick")
IDs=GENEUTM(tree,Annofile,PASTA)
GFFLOAD=GENEGFF(IDs,Annofile)
Y=len(IDs)
TotalW=EWIDTH(IDs,PASTA)
arrayYY=ARRSIM(IDs,PASTA)
arrayXYZ=ARRIMG(IDs,PASTA,TotalW,GFFLOAD)
arrayYY=MAPGENEG(arrayXYZ,arrayYY)

ExonLen={}
for line in GFFLOAD:
    line = line.strip()
    lines = line.strip().split("\t")
    ID=lines[0]
    Type=lines[7]
    Exon=lines[8]
    S=lines[9]
    E=lines[10]
    Len=abs(int(E)-int(S))+1
    LenK=ID+"\t"+Exon
    ExonLen[LenK]=Len

SPDB={}
SPDI=[]
GeneID={}
for CI in tree.get_terminals():
    Ct=CI.name
    Ct=Ct.strip()
    GeneID[Ct]=1
    Sp=Ct.strip().split(".")[0]
    SPDI.append(Sp)
    SPDI=numpyNR(SPDI)
    SPDB[Sp]=1
SpNHOG=len(SPDI)
SPDI=[]

cmd3="cp "+Annofile+" "+str(HOG)
os.system(cmd3)

#cmd4="cp "+PASTA+" "+str(HOG)
#os.system(cmd4)
#GeneNHOG=len(IDs)

DELETEEXON={}
EPDB=[]
DELETEDDB=[]
EPDB,DELETEDDB,DELETEEXON=MAPGENEYY(IDs,arrayXYZ,arrayYY,MCSN)
EPDBN=[]
SaveE={}
for pre in EPDB:
    GED1=pre.strip().split("\t")[1]
    GED2=pre.strip().split("\t")[-1]
    GED2=GED2.strip().split("'")[1]
    GED12=GED1+"\t"+GED2
    
    GED11=pre.strip().split("\t")[0]
    GED21=pre.strip().split("\t")[3]
    GED121=GED11+"\t"+GED21
    
    if DELETEEXON.get(GED12)==None and DELETEEXON.get(GED121)==None:
        EPDBN.append(pre)
        SaveE[GED12]=1
        SaveE[GED121]=1

EPDB=[]
EPDB=EPDBN
arrayXYZ=[]
arrayYY=[]

UEL={}
UTL={}
INPUT=[]
UEL,UTL,INPUT=MAPS(IDs,GFFLOAD,SaveE)

GFFLOAD=[]
MATCHEXON={}
for x in range(0,len(IDs)):
    GID=IDs[x]
    ED=UEL[GID]
    for y in range(0,len(ED)):
        exon=ED[y]
        SHARED=[]
        for z in range(0,len(EPDB)):
            IDz=EPDB[z].strip().split("\t")[0]
            EXONz=EPDB[z].strip().split("\t")[3]
            if GID==IDz and exon==EXONz:
                MATCHz=EPDB[z].strip().split("\t")[6]
                if int(MATCHz)>0:
                    SHARED.append(EPDB[z])
        if len(SHARED)>0:
            GE=GID+"\t"+exon
            if DELETEEXON.get(GE)==None:
                MATCHEXON[GE]=SHARED
LOOPMATCH={}
LMR={}
LOOPMATCH,c,LMR=LOOPMAP(INPUT,UTL,MATCHEXON,ExonLen,UEL)

Annobed=HOG+"/"+HOG+".summary"
fp = open(Annobed, "w")
Delete={}
for m in range(0,c):
    #print(m)
    DA=LOOPMATCH[str(m)]
    for n in range(0,c):
        if m!=n and Delete.get(m)==None and Delete.get(n)==None:
            DB=LOOPMATCH[str(n)]  
            Dcover=set(DA).intersection(set(DB))
            if len(Dcover)>0:
                merge=set(DA).union(set(DB))
                merge=list(merge)
                merge.sort()
                merge=deleteDuplicatedElement(merge)
                LOOPMATCH.pop(str(m), None)
                LOOPMATCH[str(m)]=merge
                Delete[n]=1
GeneNHOG=Y
for y in range(0,c):
    if Delete.get(y)==None:
        F=LOOPMATCH[str(y)]
        Class=LMR[y]
        FS=Summary(F,reftree)
        for fs in FS:
            out=HOG+"\t"+Class+"\t"+str(SpNHOG)+"\t"+str(GeneNHOG)+"\t"+fs
            print(out,file=fp)
fp.close()

Annobed=HOG+"/"+HOG+".delete"
fp = open(Annobed, "w")
for dog in DELETEDDB:
    dogbox=DELETEEXON[dog]
    outfile=dog+"\t"+Class_NR(dogbox)
    print(outfile,file=fp)
fp.close()

MATCHEXON={}
LOOPMATCH={}
UTL={}
UEL={}