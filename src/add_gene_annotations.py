#!/usr/bin/env python
# add_gene_annotations.py
# 

import sys
from gzip import open

if len(sys.argv) != 4:
	print "Usage: add_gene_annotations.py in.vars.tsv ann.tsv out.vars.tsv"
	sys.exit(1)

infile = sys.argv[1]
annfile = sys.argv[2]
outfile = sys.argv[3]

# read in variants file
fin = open(infile,'r')
variants = [[y.strip() for y in x.split("\t")] for x in fin.readlines()]
fin.close()
print len(variants)-1,"variants loaded."
# header dictionary
vhdic = {}
for k,x in enumerate(variants[0]):
	vhdic[x] = k

# read in gene annotations
fin = file(annfile,'r')
ann = [x.split() for x in fin.readlines()]
fin.close()
print len(ann)-1,"gene annotations loaded."
# gene dictionary
ahdic = {}
for k,a in enumerate(ann):
	if k == 0:
		continue
	ahdic[a[0]] = "\t".join(a[1:])

# merge and output
fout = open(outfile,'w')
fout.write("\t".join(variants[0])+"\t"+"\t".join(ann[0][1:])+"\n")
for v in variants[1:]:
	gene = v[vhdic["SnpEffGeneName"]].upper()
	fout.write("\t".join(v))
	if gene in ahdic:
		fout.write("\t"+ahdic[gene]+"\n")
	else:
		fout.write("\t"+"\t".join(["no" for k in xrange(len(ann[0])-1)])+"\n")
fout.close()

print "Complete!"

