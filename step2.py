#!/usr/bin/env python

from pybedtools import BedTool
from subprocess import call
import pandas as pd
import numpy as np
import glob
import csv
import os
import re

pd.set_option('display.max_columns',500)
gff3Cols=["seqid","source","type","start","end","score","strand","phase","attributes"]
chrMapCols=["name","role","molecule","type","genbank","rel","refseq","unit","seqLen","ucsc"]

chessFP="./data/chess2.1.2.gff"
assemblyReportFP="./data/GCF_000001405.38_GRCh38.p12_assembly_report.txt"

outDir="./res_tmp"
if not os.path.exists(outDir):
    os.mkdir(outDir)

assemblyReport=pd.read_csv(assemblyReportFP,sep="\t",names=chrMapCols,comment="#")

dfChess=pd.read_csv(chessFP,sep="\t",comment="#",names=gff3Cols)

# verify chromosome nomenclature and get the name of the nomenclature
setChromosomes=set(dfChess['seqid'])
curNomenclature=dict()
for chrID in setChromosomes:
    found=False
    for colName in ['name','molecule','genbank','refseq','ucsc']:
        if len(assemblyReport[assemblyReport[colName]==chrID])>0:
            curNomenclature.setdefault(colName,list())
            curNomenclature[colName].append(chrID)
            found=True
            break
    assert found==True,chrID+" not found"
assert len(curNomenclature)==1,"multiple nomenclatures identified: "+";".join(list(curNomenclature))

# now get names of all alternative and non-alternative chromosomes
setNonAltChrID=set(assemblyReport[((assemblyReport['unit']=="Primary Assembly")|\
                                       (assemblyReport['unit']=="non-nuclear"))][list(curNomenclature)[0]])
setAltChrID=set(assemblyReport[~((assemblyReport['unit']=="Primary Assembly")|\
                                     (assemblyReport['unit']=="non-nuclear"))][list(curNomenclature)[0]])

# create a version only with currently annotated primary scaffolds
dfGenePrim_Chess=dfChess[dfChess["seqid"].isin(setNonAltChrID)].reset_index(drop=True)
dfGenePrim_Chess.to_csv(outDir+"/chessPrim.gff",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)

# for each primary scaffold create a separate GFF
if not os.path.exists(outDir+"/prim_exons/"):
    os.mkdir(outDir+"/prim_exons/")
for chrom in set(dfGenePrim_Chess["seqid"]):
    dfGenePrim_Chess[(dfGenePrim_Chess["type"]=="exon")&(dfGenePrim_Chess["seqid"]==chrom)].reset_index(drop=True).to_csv(outDir+"/prim_exons/chessPrim_exon_"+chrom+".gff",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)

# get ID for all features
dfChess["ID"]=dfChess["attributes"].str.split("ID=",expand=True)[1].str.split(";",expand=True)[0]

#krser for attributes codes
codes=set()
def parseCodes(attr): #parses attributes and extracts codes into the global list
    global codes
    for i in attr.split(";"):
        codes.add(i.split("=")[0])

dfGene_Chess=dfChess[dfChess["type"]=="gene"].reset_index(drop=True)
dfGene_Chess.apply(lambda row: parseCodes(row["attributes"]),axis=1)
# extract each of the attributes into a separate column
for i in set(codes):
    dfGene_Chess[i]=dfGene_Chess["attributes"].str.split(i+"=",expand=True)[1].str.split(";",expand=True)[0]

# subset genes only on primary chromosomes
dfGenePrim_Chess=dfGene_Chess[dfGene_Chess['seqid'].isin(setNonAltChrID)].reset_index(drop=True)

# subset genes only on alternative scaffolds
dfGeneAlts_Chess=dfGene_Chess[dfGene_Chess["seqid"].isin(setAltChrID)].reset_index(drop=True)

#now need to identify the corresponding primary scaffold for each alternative
alt_to_prim_map=assemblyReport[assemblyReport["genbank"].isin(setAltChrID)][["molecule","genbank"]].merge(\
            assemblyReport[(assemblyReport["genbank"].isin(setNonAltChrID))\
                               &(assemblyReport["role"]=="assembled-molecule")][["molecule","genbank"]],\
                               on="molecule",how="left")
alt_to_prim_map.columns=["molecule","seqid_alt","seqid_prim"]
alt_to_prim_map.drop("molecule",axis=1,inplace=True)
dfGeneAlts_Chess=dfGeneAlts_Chess.merge(alt_to_prim_map,left_on='seqid',right_on='seqid_alt',how="left")
dfGeneAlts_Chess.drop("seqid_alt",axis=1,inplace=True)

# now let's build a transcript-level dataframe which will be merged onto the altsGene
dfT=dfChess[dfChess["type"]=="transcript"].reset_index(drop=True)
dfT['parent']=dfT.attributes.str.split("Parent=",expand=True)[1].str.split(";",expand=True)[0]
dfT=dfGeneAlts_Chess.merge(dfT[["ID",\
                                "start",\
                                "end",\
                                "source",\
                                "type",\
                                "score",\
                                "strand",\
                                "phase",\
                                "attributes",\
                                "parent"]], indicator=True, how='left', left_on="ID",right_on="parent")

dfT=dfT.rename(columns=({'ID_x':"geneID",'ID_y':"transcriptID",
                         'start_x':"geneStart",'start_y':"transcriptStart",
                         'end_x':"geneEnd",'end_y':"transcriptEnd",
                         'score_x':"geneScore",'score_y':"transcriptScore",
                         'strand_x':"geneStrand",'strand_y':"transcriptStrand",
                         'phase_x':"genePhase",'phase_y':"transcriptPhase",
                         'attributes_x':"geneAttributes",'attributes_y':"transcriptAttributes",
                         'source_x':"geneSource",'source_y':"transcriptSource",
                         'type_x':"geneType",'type_y':"transcriptType",
                         'seqid_x':"geneSeqid",'seqid_y':"transcriptSeqid"}))

# now need to enhance the data with exons
transcriptIDs=list(set(dfT["transcriptID"]))

dfExon=dfChess[(dfChess["type"]=="exon")&(dfChess["seqid"].isin(setAltChrID))].reset_index(drop=True)
dfExon["parent"]=dfExon["attributes"].str.split("Parent=",expand=True)[1].str.split(";",expand=True)[0]
dfExon=dfExon[dfExon["parent"].isin(transcriptIDs)].reset_index(drop=True)
dfExon=dfExon.merge(alt_to_prim_map,left_on='seqid',right_on='seqid_alt',how="left")
dfExon.drop("seqid_alt",axis=1,inplace=True)

if not os.path.exists(outDir+"/alts/"):
    os.mkdir(outDir+"/alts/")

def alignTranscript(case):
    #first build gffs for each chromosome
    chrs=list(set(dfT['seqid_prim']))
    for chrom in chrs:
        print(chrom)
        dfT[dfT['seqid_prim']==chrom][["seqid",\
               "transcriptSource",\
               "transcriptType",\
               "transcriptStart",\
               "transcriptEnd",\
               "transcriptScore",\
               "transcriptStrand",\
               "transcriptPhase",\
               "transcriptAttributes"]].to_csv(outDir+"/alts/"+chrom+\
                                               ".gff3",sep="\t",index=False,header=False,quoting=csv.QUOTE_NONE)
        # append exons to the transcript GFF
        dfExon[dfExon["seqid_prim"]==chrom][["seqid",\
                                        "source",\
                                        "type",\
                                        "start",\
                                        "end",\
                                        "score",\
                                        "strand",\
                                        "phase",\
                                        "attributes"]].to_csv(outDir+"/alts/"+chrom+".gff3",sep="\t",index=False,mode='a',header=False,quoting=csv.QUOTE_NONE)

alignTranscript(dfT)
