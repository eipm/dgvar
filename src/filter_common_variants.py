#!/usr/bin/env python
# filter_common_variants.py
# filter common variants in a cohort
# 

import sys
from collections import Counter
from gzip import open

if len(sys.argv) != 4:
	print "Usage: filter_common_variants.py infile freqtablefile outfile"
	sys.exit(1)

varfile = sys.argv[1]
freqtablefile = sys.argv[2]
outfile = sys.argv[3]

# read in variant frequency table
# do NOT consider genotype
fin = file(freqtablefile,'r')
ftable = [x.split() for x in fin.readlines()]
fin.close()
print len(ftable)-1,"variants loaded."

# variant dictionary
vdic = {}
for x in ftable[1:]:
	chrom,pos,ref,alt = x[0],x[1],x[2],x[3]
	vid = ":".join([chrom,pos,ref,alt])
	vdic[vid] = ''
print len(vdic),"common variants to filter."

# read in variants
fin = open(varfile,'r')
variants = [[y.strip() for y in x.split("\t")] for x in fin.readlines()]
fin.close()
print len(variants)-1,"variants loaded."
# header dictionary
hdic = {}
for k,x in enumerate(variants[0]):
	hdic[x] = k

# filter and ouptut
filt = Counter()
fout = open(outfile,'wb')
fout.write("\t".join(variants[0])+"\n")
for v in variants[1:]:# each variant
	chrom = v[hdic["#CHROM"]]
	pos = v[hdic["POS"]]
	ref = v[hdic["REF"]]
	alt = v[hdic["ALT"]]
	vid = ":".join([chrom,pos,ref,alt])
	if vid in vdic:
		filt[vid] += 1
	else:
		fout.write("\t".join(v)+"\n")
fout.close()

print "Filter out",len(filt),"unique variants,",sum([filt[v] for v in filt]),"coyies:"
print "".join(["%s\t%s\n"%(v,filt[v]) for v in filt])

print "Complete!"

