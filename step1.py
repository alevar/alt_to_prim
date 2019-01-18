#!/usr/bin/env python

import os
import pandas as pd
import subprocess

inputDir="/ccb/salz3/avaraby/genomes/human/hg38/genbank_chromosomes_prim"
gmapIdxDir="./gmapIDX"
assemblyReportFP="./data/GCF_000001405.38_GRCh38.p12_assembly_report.txt"
nomenclature="genbank"
chrMapCols=["name","role","molecule","type","genbank","rel","refseq","unit","seqLen","ucsc"]

k=10
q=2

# step 1 - preparing reference genome for alignment
# - download full reference enomes, including alternative scaffolds
# - download or extract primary scaffolds separately
# - build gmap index for the primary scaffolds
# - build blat index for the primary scaffolds

if not os.path.exists(gmapIdxDir+"/"):
    os.mkdir(gmapIdxDir+"/")
    
assemblyReport=pd.read_csv(assemblyReportFP,sep="\t",names=chrMapCols,comment="#")

print(set(assemblyReport[(assemblyReport['role']=="assembled-molecule")]["genbank"]))
for chrName in set(assemblyReport[(assemblyReport['role']=="assembled-molecule")]["genbank"]):
    curIdxDir=gmapIdxDir+"/"+chrName+"/"
    if not os.path.exists(curIdxDir):
        os.mkdir(curIdxDir)
    
    chromFaFP=inputDir+"/"+chrName+".fa"
    print(chromFaFP)
    assert os.path.exists(chromFaFP),"fasta file for chromosome "+chrName+" does not exist"
    
    # build gmap index
    subprocess.call(["gmap_build","--dir="+curIdxDir,"--db="+chrName,"-k",str(k),"-q",str(q),chromFaFP])

# Notes
# - should play a bit with the kmer-length value for gmap using -k flag
# - should also play a bit with the sampling interval for gmap using the -q flag


# THIS SHOULD BE INTERESTING AND USEFUL
# eventually, if there are any transcripts for which locations can not be accurately determined on the corresponding
# primary scaffolds, we should try aligning them against other primary scaffolds, and see if a better alignment is available

