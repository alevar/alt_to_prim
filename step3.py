#!/usr/bin/env python

from subprocess import call
import multiprocessing
import glob
import os

referenceFP="/ccb/salz3/avaraby/genomes/human/hg38/hg38_p8_genbank_chromosomes_all.fa"
chromosomeDir="/ccb/salz3/avaraby/genomes/human/hg38/genbank_chromosomes_prim"
gmapIdxDir="./gmapIDX"

# during this step the actual alignment is performed
# by using gmap and blat

outDir="./res_tmp"

childPIDs=[]

def child(gffFP):
    global referenceFP
    global chromosomeDir
    global gmapIdxDir

    chrom=gffFP.split("/")[-1].rstrip(".gff3")
    print(gffFP,chrom)
    
    gffreadArgs=["gffread","-E",gffFP,"-o-"]
    f=open(gffFP.rstrip(".gff3")+".gtf","wb")
    call(gffreadArgs,stdout=f)
    f.close()
#
    gffArgs=["gffread","-w",gffFP.rstrip(".gff3")+".fa","-g",referenceFP,gffFP.rstrip(".gff3")+".gtf"]
    call(gffArgs)

  #  third run gmap for each fasta against respective chromosome
    gmapArgs=["gmap","-D","./","-d",gmapIdxDir+"/"+chrom+"/"+chrom,"-t","2",gffFP.rstrip(".gff3")+".fa","-f","samse"]
    f=open(outDir+"/alt_alignments/"+chrom+".sam","wb")
    call(gmapArgs,stdout=f)
    f.close()

    blatArgs=["blat",chromosomeDir+"/"+chrom+".fa",gffFP.rstrip(".gff3")+".fa",outDir+"/alt_alignments/"+chrom+".psl"]
    call(blatArgs)
    
    os._exit(0)

def main():
    global childPIDs
    
    if not os.path.exists(outDir+"/alt_alignments"):
        os.mkdir(outDir+"/alt_alignments")
    for gffFP in glob.glob(outDir+"/alts/*gff3"):
        p=multiprocessing.Process(target=child,args=(gffFP,))
        childPIDs.append(p)
        p.start()

    while(len(childPIDs)>0):
        childPIDs[-1].join()
        childPIDs.remove(childPIDs[-1])

if __name__=="__main__":
    main()
