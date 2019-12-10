#!/home/taz2008/softwares/Python-2.7.6/python
# screen_germline_variants.v5.py
# screen germline variants for reporting
# AUTHOR: Tuo Zhang
# DATE: 9/22/2018
# VERSION: 4.0
# NEW: support new ClinVar vcf annotations
# PREV: add support to ExAC non-PASS filter
# PREV: modified to assign category for all genes
# PREV.PREV: bug fix: empty list issue when calling 'zip' in function 'runFilterExAC'
#      0) skip snpEff protein_protein_contact annotations; 1) extract AF info in AFR,AMR,EAS,FIN,NFE,OTH,SAS populations; 2) add a population based AF filter; 3) add a filter to remove heterozygous variants in MUTYH gene; 4) add a filter to remove snpEff "HIGH" impact predictions based on evidence "protein_protein_contact" (obselete; no use any more since such annotations have been skipped from the very beginning)
# NOTE:
# 

import sys
from string import atoi
from string import atof
from re import search
from os.path import exists
from os.path import dirname
from gzip import open

class GermlineVars:
	"Class for processing germline variants from a VCF file"
	#def __init__(self, sampleID, vcfFile, ACMGGeneFile, BROCAGeneFile, CancerGeneFile, outPrefix, logFile):
	def __init__(self, sampleID, vcfFile, outPrefix, logFile):
		self.sampleID = sampleID
		self.vcfFile = vcfFile
		#self.ACMGGeneFile = ACMGGeneFile
		#self.BROCAGeneFile = BROCAGeneFile
		#self.CancerGeneFile = CancerGeneFile
		self.outPrefix = outPrefix
		self.logFile = logFile
		self.gsep = "/"
		self.depthCutoff = 10
		self.alleleFreqCutoff = 0.35
		self.snpFSCutoff = 60.0
		self.indelFSCutoff = 200.0
		self.dbNSFPMetaIdxCutoff = 5
		self.ExACCutoff = 0.01# AF > 0.01 indicates common variant
		self.ExACAdjANCutoff = 100# at least this many people used in the AF calculation
		self.predMinCutoff = 3# at least this many available damage predictors
		self.predOppCutoff = 1# allowed number of opponent damage predictors
		self.numVars = {}# number of total variants and variants in different categories.
		self.ACMGGeneList = []# list of ACMG genes
		self.BROCAGeneList = []# list of BROCA genes
		self.CancerGeneList = []# list of cancer genes, based on COSMIC cancer gene consensus
		self.header = []# VCF header lines
		self.rawInfo = []# raw VCF data: a list, each element is a dictionary
		self.annInfo = []# annotation data extracted from INFO column: a list, same as above
		self.filters = {}# filter markers: a dictionary, each element is a list (one filter) that stores TRUE or FALSE indicating pass or not pass a filter for each variant entry
		self.numValidDamagePreds = []# per variant: number of available damage predictors
		self.category = []# category assigned to each variant: a list
		self.columnsHeader = []# headers in output file
		self.columnsBasic = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER"]# basic columns for output
		self.columnsSnpEffANN = ["SnpEffEffect","SnpEffImpact","SnpEffGeneName","SnpEffGeneID","SnpEffFeatureType","SnpEffFeatureID","SnpEffBiotype","SnpEffRank","SnpEffHGVS.c","SnpEffHGVS.p","SnpEffcDNAPos","SnpEffCDSPos","SnpEffProtPos","SnpEffDist","SnpEffMsg"]
		self.columnsFilter = []# Filter IDs for output header
		self.categoryIDs = ["P.C","P.S","LP.C","LP.P","LB.C","LB.P","B.C","VUS","O","D","R","C.F","C.E","F"]
		self.filterIDSetPredictors = ["Polyphen2_HVAR","Polyphen2_HDIV","SIFT","MutationTaster","MetaSVM","MetaLR"]# selected damage predictor filters for making consensus
		self.fieldSnpEffANN = ["SnpEffAlt","SnpEffEffect","SnpEffImpact","SnpEffGeneName","SnpEffGeneID","SnpEffFeatureType","SnpEffFeatureID","SnpEffBiotype","SnpEffRank","SnpEffHGVS.c","SnpEffHGVS.p","SnpEffcDNAPos","SnpEffCDSPos","SnpEffProtPos","SnpEffDist","SnpEffMsg"]# 16 fields in the SnpEff ANN
		self.fielddbNSFP = ["dbNSFP_genename","dbNSFP_cds_strand","dbNSFP_Ensembl_geneid","dbNSFP_Ensembl_transcriptid","dbNSFP_SIFT_pred","dbNSFP_Polyphen2_HDIV_pred","dbNSFP_Polyphen2_HVAR_pred","dbNSFP_MutationTaster_pred","dbNSFP_MetaSVM_pred","dbNSFP_MetaLR_pred","dbNSFP_Reliability_index","dbNSFP_ExAC_AF"]# dbNSFP_ExAC_AF not used
		self.fieldImpactPred = ["dbNSFP_SIFT_pred","dbNSFP_Polyphen2_HDIV_pred","dbNSFP_Polyphen2_HVAR_pred","dbNSFP_MutationTaster_pred","dbNSFP_MetaSVM_pred","dbNSFP_MetaLR_pred","dbNSFP_Reliability_index","dbNSFP_ExAC_AF"]# dbNSFP_ExAC_AF not used
		self.fieldExAC = ["exacAF","exacAdjAF","exacAF_AFR","exacAF_AMR","exacAF_EAS","exacAF_FIN","exacAF_NFE","exacAF_SAS","exacAF_OTH","exacAN_Adj","exacAF_MALE","exacAF_FEMALE"]
		self.fieldClinVarBasic = ["clinvarCLNVARID","clinvarRS","clinvarCLNSIG","clinvarCLNSIGCONF","clinvarCLNREVSTAT","clinvarGENEINFO","clinvarDBVARID","clinvarORIGIN","clinvarMC","clinvarSSR"]
		self.fieldClinVarMulti = ["clinvarCLNDN","clinvarCLNDISDB","clinvarCLNVI"]

	def verifyFiles(self):
		if not exists(self.vcfFile):
			self.printMsg("Error: cannot find the input VCF file: %s"%(self.vcfFile), 100)
		if not exists(dirname(self.logFile)):
			self.printMsg("Error: cannot find folder for logfile: %s"%(self.logFile), 101)
		if not exists(dirname(self.outPrefix)):
			self.printMsg("Error: cannot find output folder: %s"%(dirname(self.outPrefix)), 102)
		#if not exists(self.ACMGGeneFile):
		#	self.printMsg("Error: cannot find the ACMG gene list file: %s"%(self.ACMGGeneFile), 103)
		#if not exists(self.BROCAGeneFile):
		#	self.printMsg("Error: cannot find the BROCA gene list file: %s"%(self.BROCAGeneFile), 104)
		#if not exists(self.CancerGeneFile):
		#	self.printMsg("Error: cannot find the Cancer gene list file: %s"%(self.CancerGeneFile), 105)

	def readVCF(self):
		# Load the query VCF
		self.printMsg("Loading VCF file...",0)
		# read in VCF file
		fin = open(self.vcfFile,'rb')
		vcf = fin.readlines()
		fin.close()
		# record header lines
		self.header.extend([x for x in vcf if search("^#",x)])
		# correct data header line: ID --> sample ID
		theader = self.header[-1].rstrip().split("\t")
		theader[-1] = self.sampleID
		# extract raw info
		for v in vcf:
			if search("^#",v):# ignore headers
				continue
			tdic = {}# create a dictionary
			for k,element in enumerate(v.rstrip().split("\t")):
				tdic[theader[k]] = element
			self.rawInfo.append(tdic)
			# initialize self.annInfo
			self.annInfo.append({})
			# initialize self.category
			self.category.append("X")# un-categorized
			# initialize self.numValidDamagePreds
			self.numValidDamagePreds.append(0)
		self.printMsg("ok.\n", 0)
		self.numVars["Total"] = len(self.rawInfo)
		self.printMsg("%d variants loaded.\n"%(self.numVars["Total"]), 0)

	def getCancerGenes(self, tcancer):
		# Extract cancer genes
		# dictionary of items in header
		tcdic = {}
		for k,x in enumerate(tcancer[0]):
			tcdic[x.rstrip()] = k
		tcgenes = [x[tcdic["Gene Symbol"]].upper() for x in tcancer[1:]]
		tnumgenes = len(tcgenes)
		for x in tcancer[1:]:
			tcgenes.extend([y.strip().upper() for y in x[tcdic["Synonyms"]].rstrip().replace('"','').split(",") if y != ""])
		return (tnumgenes,tcgenes)

	def readGeneList(self):
		# Load ACMG/BROCA gene list as well as Cancer gene list
		self.printMsg("Loading ACMG/BROCA/Cancer gene lists...", 0)
		# read in ACMG gene list file
		fin = file(self.ACMGGeneFile,'r')
		acmg = [x.split() for x in fin.readlines()]
		fin.close()
		# read in BROCA gene list file
		fin = file(self.BROCAGeneFile,'r')
		broca = [x.split() for x in fin.readlines()]
		fin.close()
		# read in Cancer gene list file
		fin = file(self.CancerGeneFile,'r')
		cancer = [x.split("\t") for x in fin.readlines()]
		fin.close()
		# extract ACMG genes
		self.ACMGGeneList.extend([x[0].upper() for x in acmg])
		# extract BROCA genes
		self.BROCAGeneList.extend([x[0].upper() for x in broca])
		# extract Cancer genes
		####tcancergenes = [x[0].upper() for x in acmg if x[1] == "yes"] + [x[0].upper() for x in broca if x[1] == "yes"]
		####self.CancerGeneList.extend(list(set(tcancergenes)))
		tnum, tcancergenes = self.getCancerGenes(cancer)
		self.CancerGeneList.extend(tcancergenes)
		self.printMsg("ok.\n", 0)
		self.printMsg("%d ACMG genes.\n"%(len(self.ACMGGeneList)), 0)
		self.printMsg("%d BROCA genes.\n"%(len(self.BROCAGeneList)), 0)
		self.printMsg("%d Cancer genes.\n"%(tnum), 0)

	def printMsg(self, tmsg, texit):
		# print a regular message; print an error message and exit; or print a warning message and continue
		if texit > 0:# error
			print tmsg.rstrip()
			sys.exit(texit)
		elif texit == 0:# regular
			print tmsg,
		else:# warning
			flog = file(self.logFile,'a')
			flog.write(tmsg.rstrip()+"\n")
			flog.close()

	def splitInfo(self, tinfostr):
		# For one variant entry:
		# split Info string and convert into dictionary
		tinfodic = {}
		for tentry in tinfostr.split(";"):# each entry
			tpat = search("(.*)=(.*)",tentry)
			if tpat:# descriptor = value
				tkey,tval = tpat.groups()
				if tkey in tinfodic:# unique descriptor ID?
					self.printMsg("Error: non-unique descriptor ID detected in VCF INFO column - %s"%tkey,99)
				tinfodic[tkey] = tval
			else:# descriptor without a value
				tkey = tentry
				if tkey in tinfodic:# unique descriptor ID?
					self.printMsg("Error: non-unique descriptor ID detected in VCF INFO column - %s"%tkey,99)
				tinfodic[tkey] = True
		return tinfodic

	def splitSnpEffAnn(self, tsnpeffstr):
		# For one variant entry:
		# split SnpEff Annotations and convert into a dictionary (ALT allele as key)
		# keep the original annotation order for each key (ALT allele)
		# skip "protein_protein_contact|HIGH" annotations if any
		tanndic = {}
		for tann in tsnpeffstr.split(","):# each annotation
			telements = tann.split("|")
			# in case Compound variants: two or more variants affecting the annotations
			# (e.g. two consecutive SNPs conforming a MNP, two consecutive frame_shift variants that "recover" the frame).
			# in this case, the Allele field should include a reference to the other variant/s included in the annotation;
			talt = telements[0].split("-")[0]# extract the right alelle 
			teff = telements[1]
			timp = telements[2]
			if teff == "protein_protein_contact" and timp == "HIGH":# skip
				continue
			if talt not in tanndic:
				tanndic[talt] = [telements]
			else:
				tanndic[talt].append(telements)
		return tanndic

	def selectSnpEffAnn(self, tainfos, talts):
		# For one variant entry:
		# select one SnpEff annotation for the given ALT allele
		# current strategy: use the first annotation, i.e. accept the default snpEff ordering (the most highest impact)
		# Effect sort order. When multiple effects are reported, SnpEff sorts the effects the following way:
		# - Putative impact: Effects having higher putative impact are first.
		# - Effect type: Effects assumed to be more deleterious effects first.
		# - Canonical trancript before non-canonical.
		# - Marker genomic coordinates (e.g. genes starting before first).
		retdic = {}
		if "ANN" not in tainfos:# no snpEff annotations
			for i,tkey in enumerate(self.fieldSnpEffANN):
				retdic[tkey] = ""
		else:
			tsnpeffstr = tainfos["ANN"]
			tanndic = self.splitSnpEffAnn(tsnpeffstr)
			if not tanndic:# empty dictionary, i.e. no annotations
				for i,tkey in enumerate(self.fieldSnpEffANN):
					retdic[tkey] = ""
			else:
				for talt in talts:# each ALT
					if talt not in tanndic:
						self.printMsg("Error: no matched ALT allele in the snpEff Annotations - %s"%talt, 98)
				# dictionary: 16 fields, field id as key
				for i,tkey in enumerate(self.fieldSnpEffANN):
					retdic[tkey] = ','.join([tanndic[talt][0][i] for talt in talts])# use the first available annotation
		return retdic

	def addAnnInfo(self, sn, tadic, tkeys):
		# For one variant entry: sn - serial number of the variant, e.g. 0 indicates the first variant entry
		# add certain annotations (from a dictionary) to self.annInfo at a given position
		for tkey in tkeys:
			if tkey in self.annInfo[sn]:# non-unique key
				self.printMsg("Error: non-unique descriptor ID detected in snpEff ANN column - %s"%tkey,97)
			elif tkey in tadic:# if detect an annotation
				self.annInfo[sn][tkey] = tadic[tkey]

	def cleanAnndbNSFP(self, tainfos):
		# For one variant entry:
		# clean dbNSFP annotations
		retdic = {}
		for tkey in self.fielddbNSFP:
			if tkey in self.fieldImpactPred:# convert to list
				if tkey in tainfos:# available
					retdic[tkey] = tainfos[tkey].split(",")# list of predictions
				else:
					retdic[tkey] = []# no available predictitons
			else:# string
				if tkey in tainfos:# available
					retdic[tkey] = tainfos[tkey]
				else:
					retdic[tkey] = ""
		return retdic

	def cleanAnnExAC(self, tainfos, talts):
		# For one variant entry:
		# clean ExAC annotations
		retdic = {}
		for tkey in self.fieldExAC:# convert to list
			if tkey in tainfos:# available
				retdic[tkey] = tainfos[tkey].split(",")# list of predictions
			else:
				retdic[tkey] = ["." for ta in talts]# no available predictitons
		return retdic

	def cleanAnnClinVar(self, tainfos):
		# For one variant entry:
		# clean ClinVar annotations
		retdic = {}
		# basic info
		for tkey in self.fieldClinVarBasic:
			if tkey in tainfos:# available
				retdic[tkey] = tainfos[tkey]
			else:
				retdic[tkey] = ""
		# disease info
		if "clinvarCLNDN" in tainfos:# annotation available
			tCLNHits = []#(CLNDN[i],CLNDISDB[i])
			# split the annotations
			tCLNDNs = tainfos["clinvarCLNDN"].split("|")
			tCLNDISDBs = ["" for j in xrange(len(tCLNDNs))]
			if "clinvarCLNDISDB" in tainfos:
				tainfos["clinvarCLNDISDB"].split("|")
			# double check number of annotations
			if not len(tCLNDNs) == len(tCLNDISDBs):
				self.printMsg("Error: number of ClinVar annotations does not equal: %s,%d,%d"%(tainfos["clinvarCLNVARID"],len(tCLNDNs),len(tCLNDISDBs)),96)
			# pair annotations
			for i in xrange(len(tCLNDNs)):
				tCLNHits.append((tCLNDNs[i],tCLNDISDBs[i]))
			retdic["clinvarHits"] = tCLNHits
		else:# no annotation
			retdic["clinvarHits"] = []
		# omim annotation
		if "clinvarCLNVI" in tainfos:# annotation available
			tCLNVIs = [[x[:x.rfind(":")], x[x.rfind(":")+1:]] for x in tainfos["clinvarCLNVI"].split("|")]# GeneDx:GDX:1685905
			tCLNOMIMs = []
			for clnviid,clnvival in tCLNVIs:
				if clnviid == "OMIM_Allelic_Variant":# OMIM annotation available
					tCLNOMIMs.append(clnvival)
			retdic["clinvarOMIM"] = tCLNOMIMs
		else:
			retdic["clinvarOMIM"] = []
		return retdic

	def getReadStat(self, tainfos, tcall, tformat):
		# For one variant entry:
		# calculate allele frequency, strand bias (FS), genotype
		retdic = {}
		# variant type
		retdic["VariantType"] = tainfos["VariantType"]
		# strand bias (FS)
		retdic["FS"] = atof(tainfos["FS"])
		# fetch genotype info
		tfs = tformat.split(":")
		tcs = tcall.split(":")
		# double check
		if not len(tfs) == len(tcs):
			self.printMsg("Error: number of elements does not match in FORMAT and Call: %s; %s"%(tcall,tformat), 94)
		tdic = {}
		for k,x in enumerate(tfs):
			tdic[x] = tcs[k]
		retdic["GT"] = tdic["GT"]
		# allele frequency
		tfreqs = [atoi(x) for x in tdic["AD"].split(",")]
		# double check
		if len(tfreqs) != 2:# allowed, e.g. genotype '1/2' - 0,12,12
			self.printMsg("Warning: number of elements in AD != 2: %s"%tcall, -1)
		retdic["AlleleFreq"] = tfreqs
		return retdic

	def processAnnotations(self):
		# For all variant entries:
		# process annotations and extract useful info.
		self.printMsg("Processing annotations...", 0)
		for i,rinfos in enumerate(self.rawInfo):# each variant entry
			ainfos = self.splitInfo(rinfos["INFO"])# ?descriptor=?value;
			# snpEff annotations
			self.addAnnInfo(i, self.selectSnpEffAnn(ainfos, rinfos["ALT"].split(",")), self.fieldSnpEffANN)
			# dbNSFP annotations: impact predictions + ExAC_AF(not used)
			self.addAnnInfo(i, self.cleanAnndbNSFP(ainfos), self.fielddbNSFP)
			# ExAC annotations
			self.addAnnInfo(i, self.cleanAnnExAC(ainfos, rinfos["ALT"].split(",")), self.fieldExAC)
			# ClinVar annotations
			self.addAnnInfo(i, self.cleanAnnClinVar(ainfos), self.fieldClinVarBasic+["clinvarHits","clinvarOMIM"])
			# read statistics
			self.addAnnInfo(i, self.getReadStat(ainfos, rinfos[self.sampleID], rinfos["FORMAT"]), ["VariantType", "FS","GT","AlleleFreq"])
		self.printMsg("ok.\n", 0)

	def addFilters(self, tfilterset):
		# For a set of filters:
		# add each filter (list of TRUE/FALSE) to the big table: self.filters
		# input 'tfilterset' - [(tkey1, tfilter1), (tkey2, tfilter2), ...]
		for tkey,tfilter in tfilterset:
			# double check list size
			if len(tfilter) != len(self.rawInfo):
				self.printMsg("Error: number of elements in filter does not match the number of variant entries: %d,%d."%(len(tfilter),len(self.rawInfo)), 93)
			# add to self.filters
			self.filters[tkey] = tfilter

	def runFilterGATK(self):
		# a single filter based on the 'FILTER' column in VCF
		tfilter = []
		for x in self.rawInfo:# each variant
			#if x["FILTER"] == "PASS":
			ty = x["FILTER"].split(";")
			if "LowCoverage" not in ty and "LowQual" not in ty and "SnpCluster" not in ty:# PASS GATK filter, not necessary ExAC filter
				tfilter.append(True)
			else:
				tfilter.append(False)
		self.columnsFilter.append("GATK")
		return [("GATK",tfilter)]

	def runFilterVAF(self):
		# a single filter based on variant allele frequency
		tfilter = []
		for x in self.annInfo:# each variant
			# sum AF of genotypes other than '0', then divided by total read depth
			if sum(x["AlleleFreq"][1:]) >= float(sum(x["AlleleFreq"])) * self.alleleFreqCutoff:
				tfilter.append(True)
			else:
				tfilter.append(False)
		self.columnsFilter.append("VAF")
		return [("VAF",tfilter)]

	def runFilterDP(self):
		# a single filter based on read depth
		tfilter = []
		for x in self.annInfo:# each variant
			if sum(x["AlleleFreq"]) >= self.depthCutoff:
				tfilter.append(True)
			else:
				tfilter.append(False)
		self.columnsFilter.append("DP")
		return [("DP",tfilter)]

	def runFilterFS(self):
		# a single filter based on 'FS' column in VCF (strand bias)
		tfilter = []
		# filter a SNP with FS > 60.0 and an indel with FS > 200.0
		for x in self.annInfo:# each variant
			if x["VariantType"] == "SNP":# SNP
				if x["FS"] <= self.snpFSCutoff:
					tfilter.append(True)
				else:
					tfilter.append(False)
			else:# Indel and other
				if x["FS"] <= self.indelFSCutoff:
					tfilter.append(True)
				else:
					tfilter.append(False)
		self.columnsFilter.append("FS")
		return [("FS",tfilter)]

	def runFilterACMG(self):
		# a single filter based on whether a variant lies in a ACMG gene
		tfilter = []
		for x in self.annInfo:# each variant
			# use the first gene name annotation in Multi-allelic case
			if x["SnpEffGeneName"].split(",")[0].upper() in self.ACMGGeneList:
				tfilter.append(True)
			else:
				tfilter.append(False)
		self.columnsFilter.append("ACMG")
		return [("ACMG",tfilter)]

	def runFilterBROCA(self):
		# a single filter based on whether a variant lies in a BROCA gene
		tfilter = []
		for x in self.annInfo:# each variant
			# use the first gene name annotation in Multi-allelic case
			if x["SnpEffGeneName"].split(",")[0].upper() in self.BROCAGeneList:
				tfilter.append(True)
			else:
				tfilter.append(False)
		self.columnsFilter.append("BROCA")
		return [("BROCA",tfilter)]

	def runFilterCancer(self):
		# a single filter based on whether a variant lies in a cancer gene
		# currently only for ACMG+BROCA genes
		tfilter = []
		for x in self.annInfo:# each variant
			# use the first gene name annotation in Multi-allelic case
			if x["SnpEffGeneName"].split(",")[0].upper() in self.CancerGeneList:
				tfilter.append(True)
			else:
				tfilter.append(False)
		self.columnsFilter.append("CANCER")
		return [("CANCER",tfilter)]

	def runFilterSnpEffImpact(self):
		# four filters based on variant's snpEff impact (HIGH, MODERATE, LOW, MODIFIER)
		timpacts = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
		tfilters = []
		for timpact in timpacts:
			tfilters.append([])
			for x in self.annInfo:# each variant
				if timpact in x["SnpEffImpact"].split(","):# consider multi-allelic case
					tfilters[-1].append(True)
				else:
					tfilters[-1].append(False)
		self.columnsFilter.extend(["SnpEff%s"%x for x in timpacts])
		return zip(["SnpEff%s"%x for x in timpacts], tfilters)

	def runFilterSnpEffEffect(self):
		# a single filter for variant of SnpEff HIGH impact but weird SnpEffEffect (e.g. protein_protein_contact)
		tfilter = []
		for x in self.annInfo:# each variant
			if "HIGH" in x["SnpEffImpact"].split(",") and "protein_protein_contact" in x["SnpEffEffect"].split(","):# consider multi-allelic case; this may introduce false negative
				tfilter.append(False)
			else:
				tfilter.append(True)
		self.columnsFilter.append("SnpEffClean")
		return [("SnpEffClean",tfilter)]

	def runFilterPredDamage(self):
		# six filters based on six damage predictors (Polyphen2_HVAR, Polyphen2_HDIV, SIFT, MutationTaster, MetaSVM, MetaLR)
		# update self.numValidDamagePreds: a number per variant recording available damage predictions
		tfilterset = []
		# Filter - Polyphen2_HVAR
		# "D":"probably damaging"; "P":"possibly damaging"; "B":"benign"
		tfilter = []
		for k,x in enumerate(self.annInfo):# each variant
			tmark = False# benign, should be filtered out
			for y in x["dbNSFP_Polyphen2_HVAR_pred"]:# each prediction
				# double check
				if y not in ["D","P","B","."]:
					self.printMsg("Error: invalid Polyphen2_HVAR_pred prediction: %s"%y, 92)
				if y in ["D","P"]:
					tmark = True
					break
			tfilter.append(tmark)
			if len(x["dbNSFP_Polyphen2_HVAR_pred"]) > 0:# has predictions
				self.numValidDamagePreds[k] += 1
		self.columnsFilter.append("Polyphen2_HVAR")
		tfilterset.append(("Polyphen2_HVAR",tfilter))
		# Filter - Polyphen2_HDIV
		# "D":"porobably damaging"; "P":"possibly damaging"; "B":"benign"
		tfilter = []
		for k,x in enumerate(self.annInfo):# each variant
			tmark = False# benign, should be filtered out
			for y in x["dbNSFP_Polyphen2_HDIV_pred"]:# each prediction
				# double check
				if y not in ["D","P","B","."]:
					self.printMsg("Error: invalid Polyphen2_HDIV_pred prediction: %s"%y, 91)
				if y in ["D","P"]:
					tmark = True
					break
			tfilter.append(tmark)
			if len(x["dbNSFP_Polyphen2_HDIV_pred"]) > 0:# has predictions
				self.numValidDamagePreds[k] += 1
		self.columnsFilter.append("Polyphen2_HDIV")
		tfilterset.append(("Polyphen2_HDIV",tfilter))
		# Filter - SIFT
		# "D(amaging)"; "T(olerated)"
		tfilter = []
		for k,x in enumerate(self.annInfo):# each variant
			tmark = False# benign, should be filtered out
			for y in x["dbNSFP_SIFT_pred"]:# each prediction
				# double check
				if y not in ["D","T","."]:
					self.printMsg("Error: invalid SIFT_pred prediction: %s"%y, 90)
				if y in ["D"]:
					tmark = True
					break
			tfilter.append(tmark)
			if len(x["dbNSFP_SIFT_pred"]) > 0:# has predictions
				self.numValidDamagePreds[k] += 1
		self.columnsFilter.append("SIFT")
		tfilterset.append(("SIFT",tfilter))
		# Filter - MutationTaster
		# "A":"disease_causing_automatic"; "D":"disease_causing"; "N":"polymorphism"; "P":"polymorphism_automatic"
		tfilter = []
		for k,x in enumerate(self.annInfo):# each variant
			tmark = False# benign, should be filtered out
			for y in x["dbNSFP_MutationTaster_pred"]:# each prediction
				# double check
				if y not in ["A","D","N","P","."]:
					self.printMsg("Error: invalid MutationTaster_pred prediction: %s"%y, 89)
				if y in ["A","D"]:
					tmark = True
					break
			tfilter.append(tmark)
			if len(x["dbNSFP_MutationTaster_pred"]) > 0:# has predicitons
				self.numValidDamagePreds[k] += 1
		self.columnsFilter.append("MutationTaster")
		tfilterset.append(("MutationTaster",tfilter))
		# Filter - MetaSVM
		# "T(olerated)" or "D(amaging)"
		# Reliability_index: Number of observed component scores, Ranges from 1 to 10
		tfilter = []
		for k,x in enumerate(self.annInfo):# each variant
			tmark = False# benign, should be filtered out
			# check reliability score
			trscores = [atoi(y) for y in x["dbNSFP_Reliability_index"] if not search("\D+",y)]
			if len(trscores) != 0 and min(trscores) > self.dbNSFPMetaIdxCutoff:# good
				for y in x["dbNSFP_MetaSVM_pred"]:# each prediction
					# double check
					if y not in ["D","T","."]:
						self.printMsg("Error invalid MetaSVM_pred prediction: %s"%y, 88)
					if y in ["D"]:
						tmark = True
						break
				self.numValidDamagePreds[k] += 1
			else:# not good
				pass# should filtered out, not enough available predictors
			tfilter.append(tmark)
		self.columnsFilter.append("MetaSVM")
		tfilterset.append(("MetaSVM",tfilter))
		# Filter - MetaLR
		# "T(olerated)" or "D(amaging)"
		# Reliability_index: Number of observed component scores, Ranges from 1 to 10
		tfilter = []
		for k,x in enumerate(self.annInfo):# each variant
			tmark = False# benign, should be filtered out
			# check reliability score
			trscores = [atoi(y) for y in x["dbNSFP_Reliability_index"] if not search("\D+",y)]
			if len(trscores) != 0 and min(trscores) > self.dbNSFPMetaIdxCutoff:# good
				for y in x["dbNSFP_MetaLR_pred"]:# each prediction
					# double check
					if y not in ["D","T","."]:
						self.printMsg("Error: invalid MetaSVM_pred prediction: %s"%y, 87)
					if y in ["D"]:
						tmark = True
						break
				self.numValidDamagePreds[k] += 1
			else:# not good
				pass# should filtered out, not enough available predictors
			tfilter.append(tmark)
		self.columnsFilter.append("MetaLR")
		tfilterset.append(("MetaLR",tfilter))
		return tfilterset

	def runFilterClinVar(self):
		# 19 possible single values for CLNSIG:
		# 'protective', '_drug_response', 'Pathogenic', 'Affects', 'Likely_pathogenic',
		# 'Conflicting_interpretations_of_pathogenicity', 'risk_factor', 'drug_response',
		# '_association', 'Benign', 'other', '_risk_factor', 'Uncertain_significance',
		# 'not_provided', '_other', '_Affects', '_protective', 'Likely_benign', 'association'
		# can be any combinations of values above
		# 13 filters based on whether a variant has a specific CLNSIG
		# (Pathogenic, Likely_pathogenic, Likely_benign, Benign, drug_response, risk_factor, Uncertain_significance,
		# association, other, not_provided, Affects, protective, Conflicting_interpretations_of_pathogenicity)
		tsigs = ["Pathogenic","Likely_pathogenic","Likely_benign","Benign","drug_response","risk_factor","Uncertain_significance","association","other","not_provided","Affects","protective","Conflicting_interpretations_of_pathogenicity"]
		tlabels = ["P","LP","LB","B","DR","RF","US","AS","O","NP","AF","PO","CIP"]
		tfilters = []
		for k,tsig in enumerate(tsigs):
			tfilters.append([])
			for x in self.annInfo:# each variant
				# debug: output CLNSIG value
				#if k == 0:
				#	print x["clinvarCLNSIG"]
				# assign filter
				if tsig in x["clinvarCLNSIG"]:
					tfilters[-1].append(True)
				else:
					tfilters[-1].append(False)
		self.columnsFilter.extend(["CLNSIG_%s"%(tlabels[k]) for k,x in enumerate(tsigs)])
		return zip(["CLNSIG_%s"%(tlabels[k]) for k,x in enumerate(tsigs)], tfilters)

#	def runFilterClinVarConflict(self):
#		# a single filter based on whether a variant has conflict ClinVar annotations (pathogenic and benign)
#		# filter return True if no conflict is detected; otherwise False
#		tfilter = []
#		for k in xrange(self.numVars["Total"]):# each variant
#			if (self.filters["CLNSIG5"][k] or self.filters["CLNSIG4"][k]) and (self.filters["CLNSIG3"][k] or self.filters["CLNSIG2"][k]):# conflict
#				tfilter.append(False)
#			else:# ok
#				tfilter.append(True)
#		self.columnsFilter.append("ClinVarConflictAnn")
#		return [("ClinVarConflictAnn",tfilter)]

	def runFilterExAC(self):
		# three filters based on the frequency (AF, AdjAF, and AF in sub-population) of a given variant in the ExAC database
		# check AN_Adj to make sure there is enough people for AF calculation
		# frequency > 0.01 --> common SNPs
		tfilterAF = []
		tfilterAdjAF = []
		tfilterSubAF = []
		for x in self.annInfo:# each variant
			# choose the minimum ExAC freq in case multiple values are detected
			# get AN_Adj
			tadjan = 99999# dummy number
			if len(x["exacAN_Adj"]) != 0 and search("\d+",x["exacAN_Adj"][0]):# available AN_Adj
				tadjan = atoi(x["exacAN_Adj"][0])
			# AF
			####taf = x["exacAF"]
			taf = [atof(y) for y in x["exacAF"] if search("\d+",y)]
			if len(taf) == 0:# not available, this is a novel variant
				tfilterAF.append(True)
			elif min(taf) <= self.ExACCutoff:# at least one ALT variant is rare; may introduce false positive call
				tfilterAF.append(True)
			else:# this is a common variant
				tfilterAF.append(False)
			# AdjAF
			tadjaf = [atof(y) for y in x["exacAdjAF"] if search("\d+",y)]
			if tadjan < self.ExACAdjANCutoff:# not enough people(samples) in AdjAF calculation
				tfilterAdjAF.append(True)
			elif len(tadjaf) == 0:# not available, this is a novel variant
				tfilterAdjAF.append(True)
			elif min(tadjaf) <= self.ExACCutoff:# at least one ALT variant is rare; may introduce false positive call
				tfilterAdjAF.append(True)
			else:# this is a common variant
				tfilterAdjAF.append(False)
			# AF in sub-population
			# pair AFs for each ALT allele
			tafSubAFsRaw = zip(x["exacAF_AFR"],x["exacAF_AMR"],x["exacAF_EAS"],x["exacAF_FIN"],x["exacAF_NFE"],x["exacAF_SAS"])
			# AF: string to float
			tafSubAFs = [[atof(y) for y in v if search("\d+",y)] for v in tafSubAFsRaw]
			# add a dummy low number to avoid empty list error from using 'max' function
			for v in tafSubAFs:
				v.append(-99)
			# get the maximum AF across all sub-populations for each ALT allele
			tafSubAFsMax = [max(v) for v in tafSubAFs]
			# perform filtering
			if tadjan < self.ExACAdjANCutoff:# not enough people(samples) in AdjAF calculation
				tfilterSubAF.append(True)
			elif len(tafSubAFsMax) == 0:# not available, this is a novel variant
				tfilterSubAF.append(True)
			elif min(tafSubAFsMax) <= self.ExACCutoff:# at least one ALT variant is rare; may introduce false positive call
				tfilterSubAF.append(True)
			else:
				tfilterSubAF.append(False)
		self.columnsFilter.extend(["ExACAF","ExACAdjAF","ExACSubAF"])
		return [("ExACAF",tfilterAF),("ExACAdjAF",tfilterAdjAF),("ExACSubAF",tfilterSubAF)]

#	def runFilterClinvarCOMMON(self):
#		# a single filter based on the COMMON entry by ClinVar
#		tfilter = []
#		for x in self.annInfo:# each variant
#			if x["clinvarCOMMON"] == "1":
#				tfilter.append(False)
#			else:
#				tfilter.append(True)
#		self.columnsFilter.append("CLINVARCOMMON")
#		return [("CLINVARCOMMON",tfilter)]

	def runFilterPredDamageConsensus(self):
		# two filters based on consensus of the six damage predictors
		tfilterT = []# consensus (T)olerance
		tfilterD = []# consensus (D)amage
		for k in xrange(self.numVars["Total"]):# each variant
			if self.numValidDamagePreds[k] < self.predMinCutoff:# not enough available damage predictors, for 'benign'
				tfilterT.append(False)
				tfilterD.append(False)
				continue
			# enough damage predictors
			tpreds = [self.filters[x][k] for x in self.filterIDSetPredictors]
			tnumDamagePreds = len([x for x in tpreds if x])
			if tnumDamagePreds >= len(self.filterIDSetPredictors) - self.predOppCutoff:
				tfilterT.append(False)
				tfilterD.append(True)
			elif tnumDamagePreds <= self.predOppCutoff:
				tfilterT.append(True)
				tfilterD.append(False)
			else:
				tfilterT.append(False)
				tfilterD.append(False)
		self.columnsFilter.extend(["PredT","PredD"])
		return [("PredT",tfilterT),("PredD",tfilterD)]

	def runFilterMUTYH(self):
		# a single filter based on whether a variant is a heterozygous variant in MUTYH gene
		tfilter = []
		for x in self.annInfo:# each variant
			# use the first gene name annotation in Multi-allelic case
			if x["SnpEffGeneName"].split(",")[0].upper() == "MUTYH":
				if "0" not in x["GT"]:# homozygous variant in MUTYH
					tfilter.append(True)
				else:# heterozygous variant in MUTYH:
					tfilter.append(False)
			else:# variant not in MUTYH
				tfilter.append(True)
		self.columnsFilter.append("MUTYH_HOM")
		return [("MUTYH_HOM",tfilter)]

	def processFilters(self):
		# For all filters:
		# process filtering on the variants set
		self.printMsg("Performing filters...", 0)
		# 1) Filter 'GATK' - 'FILTER' column in VCF
		self.addFilters(self.runFilterGATK())
		# 2) Filter 'VAF' - variant allele frequency
		self.addFilters(self.runFilterVAF())
		# 3) Filter 'DP' - read depth
		self.addFilters(self.runFilterDP())
		# 4) Filter 'FS' - strand bias
		self.addFilters(self.runFilterFS())
		# 5) Filter 'ACMG' - ACMG gene?
		self.addFilters(self.runFilterACMG())
		# 6) Filter 'BROCA' - BROCA gene?
		self.addFilters(self.runFilterBROCA())
		# 7) Filter 'CANCER' - Cancer gene?
		self.addFilters(self.runFilterCancer())
		# 8-11) Filters 'SnpEffHIGH', 'SnpEffMODERATE', 'SnpEffLOW', 'SnpEffMODIFIER'
		self.addFilters(self.runFilterSnpEffImpact())
		# 12-17) Filters 'Polyphen2_HVAR', 'Polyphen2_HDIV', 'SIFT', 'MutationTaster', 'MetaSVM', 'MetaLR'
		self.addFilters(self.runFilterPredDamage())
		# 18-26) Filters 'CLNSIG0', 'CLNSIG1', 'CLNSIG2', 'CLNSIG3', 'CLNSIG4', 'CLNSIG5', 'CLNSIG6', 'CLNSIG7', 'CLNSIG255'
		# 18-30) Fitlers 'CLNSIG_P', 'CLNSIG_LP', etc
		self.addFilters(self.runFilterClinVar())
		# 27-29) Filters 'ExACAF', 'ExACAdjAF', 'ExACSubAF' - common variant?
		self.addFilters(self.runFilterExAC())
		# 30) Filter 'CLINVARCOMMON'
		#self.addFilters(self.runFilterClinvarCOMMON())
		# 31-32) Filter 'PredT', 'PredD'
		self.addFilters(self.runFilterPredDamageConsensus())
		# 33) Filter 'ClinVarConflictAnn'
		#self.addFilters(self.runFilterClinVarConflict())
		# 34) Filter 'MUTYH_HOM'
		self.addFilters(self.runFilterMUTYH())
		# 35) Filter 'SnpEffClean'
		self.addFilters(self.runFilterSnpEffEffect())
		self.printMsg("ok.\n", 0)

	def statFilter(self, tfiltName):
		# statistic on a given filter: how many variants have passed that filter, and how many do not?
		if tfiltName not in self.filters:
			self.printMsg("Error: invalid filter name for statFilter(): %s"%tfiltName, 85)
		tnumPASS = len([x for x in self.filters[tfiltName] if x])
		tnumFail = len([x for x in self.filters[tfiltName] if not x])
		self.printMsg("Filter %s: %d variants passed and %d failed.\n"%(tfiltName,tnumPASS,tnumFail), 0)

	def statFilterSet(self, tfiltNameSet):
		# statistic on each filter in a filter set
		for tfiltName in tfiltNameSet:
			self.statFilter(tfiltName)

	def assignCategory(self):
		# assign category for each variant
		# main-category
		# P - pathogenic; LP - likely pathogenic; LB - likely benign; B - benigh; VUS - vairants of unknown significance; F - filtered variants; C - currently filtered out but may need manually check; D - drug response; H - histocompatibility; O - other (confers sensitivity; risk factor; association; protective; Affects)
		# sub-category
		# P: P.C - pathogenic inferred from ClinVar
		# P: P.S - pathogenic inferred from snpEff 'HIGH' impact annotation
		# LP: LP.C - likely pathogenic inferred from ClinVar
		# LP: LP.P - likely pathogenic inferred from consensus disease predictions
		# LB: LB.C - likely benign inferred from ClinVar
		# LB: LB.P - likely benign inferred from consensus disease predictions
		# B: B.C - benign inferred from ClinVar
		# C: C.F - check inferred by low VAF/DP or high FS
		# C: C.E - check inferred by inconsistency (P/LP in ClinVar but high frequency in ExAC; both benign and pathogenic by ClinVar); heterozygous variant in MUTYH gene; weird SnpEffEffect annotations (e.g. protein_protein_contact)
		self.printMsg("Assigning category...", 0)
		for k in xrange(self.numVars["Total"]):# each variant
			if not self.filters["GATK"][k]:# not PASS GATK filtering
				self.category[k] = "F"
			#elif not self.filters["ACMG"][k] and not self.filters["BROCA"][k]:# neither ACMG nor BROCA genes
			#	self.category[k] = "F"
			#elif not self.filters["VAF"][k] or not self.filters["DP"][k] or not self.filters["FS"][k]:# low VAF/DP or high FS
			elif not self.filters["VAF"][k] or not self.filters["DP"][k]:# low VAF/DP
				self.category[k] = "C.F"
			elif self.filters["CLNSIG_CIP"][k]:# conflict annotations from ClinVar
				self.category[k] = "C.E"
			elif self.filters["CLNSIG_P"][k]:# clinvar pathogenic
				if self.filters["ExACAdjAF"][k] and self.filters["MUTYH_HOM"][k]:# ExAC rare + non-heterozygous in MUTYH
					self.category[k] = "P.C"
				else:# ExAC common
					self.category[k] = "C.E"
			elif self.filters["CLNSIG_LP"][k]:# clinvar likely pathogenic
				if self.filters["ExACAdjAF"][k] and self.filters["MUTYH_HOM"][k]:# ExAC rare + non-heterozygous in MUTYH
					self.category[k] = "LP.C"
				else:# ExAC common
					self.category[k] = "C.E"
			elif self.filters["CLNSIG_B"][k]:# clinvar benign
				self.category[k] = "B.C"
			elif self.filters["CLNSIG_LB"][k]:# clinvar likely benign
				self.category[k] = "LB.C"
			elif self.filters["SnpEffHIGH"][k]:# snpEff HIGH impact
				if self.filters["ExACAdjAF"][k] and self.filters["ExACSubAF"][k] and self.filters["MUTYH_HOM"][k] and self.filters["SnpEffClean"][k]:# ExAC rare + non-heterozygous in MUTYH + no weird SnpEffEffect annotation (e.g. protein_protein_contact)
					self.category[k] = "P.S"
				else:# ExAC common
					self.category[k] = "VUS"
			elif self.filters["PredD"][k]:# consensus damage prediction
				if self.filters["ExACAdjAF"][k] and self.filters["ExACSubAF"][k] and self.filters["MUTYH_HOM"][k]:# ExAC rare + non-heterozygous in MUTYH
					self.category[k] = "LP.P"
				else:# ExAC common
					self.category[k] = "VUS"
			elif self.filters["PredT"][k]:# consensus tolerance prediction
				if self.filters["ExACAdjAF"][k] and self.filters["ExACSubAF"][k]:# ExAC rare
					self.category[k] = "VUS"
				else:# ExAC common
					self.category[k] = "LB.P"
			elif self.filters["CLNSIG_DR"][k]:# clinvar drug response
				self.category[k] = "D"
			elif self.filters["CLNSIG_RF"][k]:# clinvar risk factor
				self.category[k] = "R"
			elif self.filters["CLNSIG_O"][k]:# clinvar other
				self.category[k] = "O"
			else:
				self.category[k] = "VUS"
		self.printMsg("ok.\n", 0)

	def statCategory(self):
		# statistic on assigned categories
		# initialize entry counts for each category
		for x in self.categoryIDs:
			self.numVars["Cate%s"%x] = 0
		# go over each variant's category
		for x in self.category:
			self.numVars["Cate%s"%x] += 1
		# output
		for x in self.categoryIDs:
			self.printMsg("Category %s = %d\n"%(x,self.numVars["Cate%s"%x]), 0)

	def writeBasicInfoToStr(self, ti):
		# for a single variant entry (at ith location):
		# collect basic info and convert to list of strings for output
		tStrList = []
		tStrList.extend([self.rawInfo[ti][x] for x in self.columnsBasic])
		return tStrList

	def writeSnpEffAnnToStr(self, ti):
		# for a single variant entry (at ith location):
		# collect snpEff annotations and convert to list of strings for output
		tStrList = []
		tStrList.extend([self.annInfo[ti][x] for x in self.columnsSnpEffANN])
		return tStrList

	def writeDbNSFPAnnToStr(self, ti):
		# for a single variant entry (at ith location):
		# collect dbNSFP annotations and convert to list of strings for output
		tStrList = []
		tStrList.extend([",".join(self.annInfo[ti][x]) for x in self.fieldImpactPred])
		return tStrList

	def writeExACAnnToStr(self, ti):
		# for a single variant entry (at ith location):
		# collect ExAC annotations and convert to list of strings for output
		tStrList = []
		tStrList.extend([",".join(self.annInfo[ti][x]) for x in self.fieldExAC])
		return tStrList

	def getOMIMLink(self, tomimvarid):
		# generate a link to variant description on OMIM website
		if "." not in tomimvarid:
			self.printMsg("Warning: weird OMIM_Allelic_Variant: %s"%tomimvarid, -1)
			return tomimvarid
		return "http://www.omim.org/entry/%s#%s"%(tomimvarid.split('.')[0],tomimvarid.split('.')[1])

	def writeClinvarAnnToStr(self, ti):
		# for a single variant entry (at ith location):
		# collect ClinVar annotations and convert to list of strings for output
		tStrList = []
		# add clinvar basic info
		for x in self.fieldClinVarBasic:
			tStrList.append(self.annInfo[ti][x])
		## add clinvarPM: Variant is Precious(Clinical,Pubmed Cited)
		#if self.annInfo[ti]["clinvarPM"]:# True
		#	tStrList.append("Y")
		#else:
		#	tStrList.append("N")
		# add clinvarHits: (CLNDNs,CLNDISDBs)
		tStrList.append(",".join(["%s|%s"%x for x in self.annInfo[ti]["clinvarHits"]]))
		# add clinvarOMIM
		tStrList.append(",".join([self.getOMIMLink(x) for x in self.annInfo[ti]["clinvarOMIM"]]))
		return tStrList

	def writeExtendInfoToStr(self, ti):
		# for a single variant entry (at ith location):
		# collect extended info and convert to list of strings for output
		tStrList = []
		# add VariantType: SNP/DEL*/INS*/MULTI*
		tStrList.append(self.annInfo[ti]["VariantType"])
		# add FS: strand bias
		tStrList.append("%.3f"%(self.annInfo[ti]["FS"]))
		# add GT: 0/1; 1/1; 1/2
		tStrList.append(self.annInfo[ti]["GT"])
		# add AlleleFreq: 0/1; 0/1/2
		tStrList.append("/".join(["%d"%x for x in self.annInfo[ti]["AlleleFreq"]]))
		return tStrList

	def writeCategoryToStr(self, ti):
		# for a single variant entry (at ith location):
		# collect category and convert to list of strings for output
		return [self.category[ti]]

	def writeFiltersToStr(self, ti):
		# for a single variant entry (at ith location):
		# collect filter and convert to list of strings for output
		tStrList = []
		for x in self.columnsFilter:
			if self.filters[x][ti]:# pass filter
				tStrList.append("Y")
			else:# not pass filter
				tStrList.append("N")
		return tStrList

	def writeEntryToStr(self, ti):
		# for a single variant entry (at ith location):
		# collect basic/annotation info and convert to list of string for output
		toutStrList = []
		# category
		toutStrList.extend(self.writeCategoryToStr(ti))
		# basic info
		toutStrList.extend(self.writeBasicInfoToStr(ti))
		# extended info
		toutStrList.extend(self.writeExtendInfoToStr(ti))
		# snpEff annotation info
		toutStrList.extend(self.writeSnpEffAnnToStr(ti))
		# dbNSFP annotation info (raw predictions + ExAC)
		toutStrList.extend(self.writeDbNSFPAnnToStr(ti))
		# ExAC annotation
		toutStrList.extend(self.writeExACAnnToStr(ti))
		# ClinVar annotation info
		toutStrList.extend(self.writeClinvarAnnToStr(ti))
		# filters
		toutStrList.extend(self.writeFiltersToStr(ti))
		return toutStrList

	def writeHeaderToStr(self):
		# for all variant entries:
		# should be in the same order as in self.writeEntryToStr
		self.columnsHeader.append("Category")
		self.columnsHeader.extend(self.columnsBasic)
		self.columnsHeader.extend(["VariantType","FS","GT","AlleleFreq"])
		self.columnsHeader.extend(self.columnsSnpEffANN)
		self.columnsHeader.extend(self.fieldImpactPred)
		self.columnsHeader.extend(self.fieldExAC)
		self.columnsHeader.extend(self.fieldClinVarBasic)
		self.columnsHeader.extend(["clinvarHits","clinvarOMIM"])
		self.columnsHeader.extend(self.columnsFilter)

	def writeDataTofile(self):
		# for all variant entries:
		# write results to file
		self.printMsg("Writing results to file...", 0)
		# prepare file header line
		self.writeHeaderToStr()
		# output
		fall = open(self.outPrefix+self.sampleID+".info.all.txt.gz",'wb')
		#frep = file(self.outPrefix+self.sampleID+".info.report.txt",'w')
		fall.write("\t".join(self.columnsHeader)+"\n")
		#frep.write("\t".join(self.columnsHeader)+"\n")
		# output
		for k in xrange(self.numVars["Total"]):
			fall.write("\t".join(self.writeEntryToStr(k))+"\n")
			#if self.category[k] != "F":# anything that are not filtered
			#	frep.write("\t".join(self.writeEntryToStr(k))+"\n")
		fall.close()
		#frep.close()
		self.printMsg("ok.\n", 0)

	def run(self):
		# run main program
		self.printMsg("######  Start processing %s  ######\n"%(self.sampleID), 0)
		# variant screening
		self.verifyFiles()
		#self.readGeneList()
		self.readVCF()
		self.processAnnotations()
		self.processFilters()
		self.assignCategory()
		self.writeDataTofile()
		# statistics on results
		self.printMsg("\nStats\n", 0)
		self.statFilterSet(self.columnsFilter)
		self.statCategory()
		self.printMsg("###### Complete processing %s ######\n"%(self.sampleID), 0)

def help():
	print "Usage: screen_germline_variants.py sampleID inputFile outPrefix logFile ACMGGeneFile BROCAGeneFile CancerGeneFile"
	print "\tsampleID: PM sample id"
	print "\tinputFile: the input VCF file."
	print "\toutPrefix: the output prefix (outFile = outPrefix+sampleID+.info.all.txt)."
	print "\tlogFile: the log file (in appending mode, no deleting)."
	print "\tACMGGeneFile: the list of ACMG genes."
	print "\tBROCAGeneFile: the list of BROCA genes."
	print "\tCancerGeneFile: the list of consensus cancer genes by COSMIC."

def main(targs):
	#sampleID, inputFile, outPrefix, logFile, ACMGGeneFile, BROCAGeneFile, CancerGeneFile = None, None, None, None, None, None, None
	sampleID, inputFile, outPrefix, logFile = None, None, None, None
	try:
		sampleID = targs[1]
		inputFile = targs[2]
		outPrefix = targs[3]
		logFile = targs[4]
		#ACMGGeneFile = targs[5]
		#BROCAGeneFile = targs[6]
		#CancerGeneFile = targs[7]
	except IndexError:
		help()
		sys.exit(1)

	#var = GermlineVars(sampleID, inputFile, ACMGGeneFile, BROCAGeneFile, CancerGeneFile, outPrefix, logFile)
	var = GermlineVars(sampleID, inputFile, outPrefix, logFile)
	var.run()

if __name__ == "__main__":
	main(sys.argv)

