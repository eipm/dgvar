#!/usr/bin/env python
# select_candidate_variants.py
# select candidate germline variants
# 

import sys
from gzip import open
from string import atof
from os.path import exists

# function
def fix_splice_region(telements,thdic,cname):
	if cname == "Variant_Classification":
		return telements[thdic[cname]].replace("Splice_Region","Intron")
	return telements[thdic[cname]]

# main
if len(sys.argv) != 4:
	print "Usage: select_candidate_variants.py sample_id infile outfile"
	sys.exit(1)

sid = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

exaccut = 0.01
####caddcut = 15
####outcols = ["Category","#CHROM","POS","ID","REF","ALT","FILTER","VariantType","FS","GT","AlleleFreq","SnpEffEffect","SnpEffGeneName","SnpEffHGVS.c","SnpEffHGVS.p","exacAdjAF","clinvarHits","Hugo_Symbol","Chromosome","Start_Position","End_Position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","HGVSc","HGVSp","HGVSp_Short","Transcript_ID"]
outcols = ["Category","#CHROM","POS","ID","REF","ALT","FILTER","VariantType","GT","AlleleFreq","SnpEffEffect","SnpEffGeneName","SnpEffFeatureID","SnpEffHGVS.c","SnpEffHGVS.p","exacAdjAF","clinvarHits"]

# file exists?
if not exists(infile):
	print infile,": not available."
	sys.exit(1)
# read in combined annotation file
fin = open(infile,'rb')
variants = [[y.rstrip() for y in x.split("\t")] for x in fin.readlines()]
fin.close() 
#print sid,":",len(variants)-1,"variants in total,",

# screen candidates and output
print "PMid\tnumTotal\tnumHits\tnumClvP\tnumPass\tnumQC\tnumExAC\tnumClvB\tnumImpact\tnumInbreed"
# screen
numHits = numQC = numExAC = numClvP = numClvB = numImpact = numInbreed = numPass = 0
fout = open(outfile,'wb')
fout.write("Tumor_Sample_Barcode\t"+"\t".join(outcols)+"\n")
hdic = {}# header dictionary
for k,elements in enumerate(variants):
	if k == 0:# header line
		for j,y in enumerate(elements):
			hdic[y] = j
		continue
	mark = True
	category = elements[hdic["Category"]]
	gname = elements[hdic["Hugo_Symbol"]]
	exac = [x for x in elements[hdic["exacAdjAF"]].split(",") if x != "." and x != "NA"]
	#caddpass = elements[hdic["CADDGeneMatch"]].split(",")
	#caddphred = elements[hdic["CADDPhredScore"]].split(",")
	filt = elements[hdic["FILTER"]]
	impact = elements[hdic["SnpEffImpact"]].split(",")
	# filter
	if category == "C.F" or category == "F":
		mark = False
		numQC += 1
	####if gname.upper() not in genelist:
	####	mark = False
	####if category == "P.S":# protein truncating + exac < 1%
	####	mark = True
	elif len(exac) != 0 and min([atof(x) for x in exac]) > exaccut:
		mark = False
		numExAC += 1
	elif category == "P.C" or category == "LP.C":# annotated pathogenic/likely pathogenic variants by ClinVar
		mark = True
		numClvP += 1
	elif category == "B.C" or category == "LB.C":# annotated benign/likely benign variants by ClinVar
		mark = False
		numClvB += 1
	elif "HIGH" not in impact:
		mark = False
		numImpact += 1
	elif "InbreedingCoeff_Filter" in filt or "VQSRTrancheINDEL" in filt or "VQSRTrancheSNP" in filt:# keep AC_Adj0_Filter
		mark = False
		numInbreed += 1
	else:
		numPass += 1
	# output
	if mark:
		numHits += 1
		###fout.write(sid+"\t"+"\t".join([elements[hdic[x]] for x in outcols])+"\n")
		fout.write(sid+"\t"+"\t".join([fix_splice_region(elements,hdic,x) for x in outcols])+"\n")# fix Splice_Region annotation
fout.close()
# double check
if numHits != numClvP + numPass:
	print "Weird",sid,numHits,numClvP,numPass
	sys.exit(2)
# PMid\tnumTotal\tnumHits\tnumClvP\tnumPass\tnumQC\tnumExAC\tnumClvB\tnumImpact\tnumInbreed
print "%s\t%s"%(sid,"\t".join(["%d"%x for x in [len(variants)-1,numHits,numClvP,numPass,numQC,numExAC,numClvB,numImpact,numInbreed]]))

#print "Complete!"

