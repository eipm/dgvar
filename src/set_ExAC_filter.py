#!/usr/bin/env python
# set_ExAC_filter.py
# assign 'ExAC filter' to a given set of variants
# 

import sys
from gzip import open
from re import search

if len(sys.argv) != 4:
	print "Usage: set_ExAC_filter.py in.vcf.gz ExAC.filtered.vcf.gz out.vcf.gz"
	sys.exit(1)

infile = sys.argv[1]
exacfile = sys.argv[2]
outfile = sys.argv[3]

# load ExAC non-PASS entries
print "Loading ExAC non-PASS entries...",
sys.stdout.flush()
fin = open(exacfile,'rb')
exac = [x.rstrip().split("\t") for x in fin.readlines() if not search("^#",x)]
fin.close()
print "ok.",len(exac),"variants detected."
sys.stdout.flush()

# build variant dictionary
vardic = {}
for entry in exac:# each variant
	chrom = entry[0]
	pos = entry[1]
	ref = entry[3]
	alt = entry[4]
	filt = entry[6]
	eid = ":".join([chrom,pos,ref,alt])
	if eid not in vardic:
		vardic[eid] = filt
	else:
		print "Error: duplicate variant entry!",eid
		sys.exit(2)

# read in input vcf
print "Loading vcf file to annotate...",
sys.stdout.flush()
fin = open(infile,'rb')
vcf = fin.readlines()
fin.close()
print "ok."
sys.stdout.flush()

# output
numHits = 0
print "Annotating and writing to file...",
fout = open(outfile,'wb')
for line in vcf:
	if search("^#",line):# header line
		fout.write(line)
	else:
		entry = line.split("\t")
		chrom = entry[0]
		pos = entry[1]
		ref = entry[3]
		alts = entry[4].split(",")
		filts = [x for x in entry[6].split(";") if x != "PASS"]
		mark = False
		for alt in alts:# each alternative allele
			eid = ":".join([chrom,pos,ref,alt])
			if eid in vardic:
				filts.append(vardic[eid])
				mark = True
		if mark == True:
			numHits += 1
		reportfilts = []
		if len(filts) == 0:
			reportfilts.append("PASS")
		else:
			reportfilts = list(set(filts))# remove duplicates if any
		fout.write("\t".join(entry[:6]+[";".join(reportfilts)]+entry[7:]))
fout.close()

print "ok.",numHits,"hits assigned new filters."

print "Complete!"

