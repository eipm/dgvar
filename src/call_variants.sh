#!/usr/bin/env sh
# call_variants.sh
# identify high-confidence deleterious germline variants (DGVs)
# Version 1.5
# Author: Tuo Zhang
# Date: 10/2/2018
# 
#
# modified based on IPM_germline_pipeline.dev.sh
# call variants and add preliminary annotations for TCGA bams
# currently only focus on HaloPlex target regions
# test: using new ClinVar annotations
# 

if [ $# -ne 2 ]
then
	echo Usage: call_variants.sh bamfile sample_id
	exit 1
fi

# BAM file used for variant calling
bamfile=$1
# Sample id (arbitrary name without any space)
sid=$2

# folders
basedir=/athena/ipm/scratch/users/taz2008/redteam2
pipedir=${basedir}/germline_variant_v2
workdir=${basedir}/request/Bishoy/171101/updated/info
srcdir=${workdir}/../code
##codedir=${pipedir}/standalone/src
dbdir=${pipedir}/standalone/db
tooldir=${pipedir}/standalone/tools
vardir=${workdir}/raw_variants/${sid}
varlogdir=${vardir}/logs
vartmpdir=${vardir}/snp_tmp
anndir=${vardir}/snpeff_gatk_v2
annlogdir=${anndir}/logs
anntmpdir=${anndir}/snpeff_tmp
filtdir=${anndir}/filt
filtlogdir=${filtdir}/logs
#cadddir=${workdir}/cadd/${sid}
#caddlogdir=${cadddir}/logs
#mafdir=${workdir}/../../info/maf
#maflogdir=${mafdir}/logs
####ipmannsuitedir=${vardir}/ipmannsuite
####ipmannsuitelogdir=${ipmannsuitedir}/logs

# tools
samtools=${tooldir}/samtools-0.1.19/samtools
java=/softlib/exe/x86_64/bin/java
java7=/home/taz2008/softwares/jre1.7.0_51/bin/java
python=/home/taz2008/softwares/Python-2.7.6/python
#perl=/home/taz2008/perl/perl-5.22.2/bin/perl
gatk=${tooldir}/gatk/2.5.2/GenomeAnalysisTK.jar
snpeff=/home/taz2008/softwares/snpEff_v4.2/snpEff.jar
snpsift=/home/taz2008/softwares/snpEff_v4.2/SnpSift.jar
#vcf2maf=/athena/ipm/scratch/users/taz2008/softwares/mskcc-vcf2maf-2f82fa4/vcf2maf.pl
####ipmannsuite=/home/kwe2001/scratch/IPM_Annotation_Suite.pl

# configuration
config=/home/taz2008/softwares/snpEff_v4.2/snpEff.config
mem=10g
anndb=GRCh37.75
#flag="-no-upstream -no-downstream -no-intergenic -formatEff -v -o gatk"
flag="-no-upstream -no-downstream -no-intergenic -v"
#dbNSFPcols="genename,Uniprot_acc,Uniprot_id,Uniprot_aapos,cds_strand,Ensembl_geneid,Ensembl_transcriptid,aapos_SIFT,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,Reliability_index,ExAC_AC,ExAC_AF,clinvar_rs,clinvar_clnsig,clinvar_trait"
ExACcols="AF,AdjAF,AN_Adj,AF_AFR,AF_AMR,AF_EAS,AF_FIN,AF_NFE,AF_SAS,AF_OTH,AF_MALE,AF_FEMALE"
CLNcols="CLNVARID,CLNDN,CLNDISDB,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNVI,DBVARID,ORIGIN,GENEINFO,MC,RS,SSR"

# databases
genomedir=${dbdir}/genome
refseq=${genomedir}/human_g1k_b37.fasta
snpdb=${dbdir}/dbsnp/dbsnp_137.b37.vcf
####clinvar=${dbdir}/clinvar/clinvar_20160531.vcf
####clinvar=${dbdir}/clinvar/clinvar_20160531.multiallelic.expanded.sorted.vcf.gz
clinvar=${dbdir}/clinvar.20180805/clinvar_20180805.fix.vcf.gz
#dbnsfp=${dbdir}/dbNSFP_v2.9.1/dbNSFP2.9.1.txt.gz
####exac=${dbdir}/ExAC/ExAC.r0.3.1.sites.vep.vcf.gz
####exac=${dbdir}/ExAC/ExAC.r0.3.1.sites.vep.multiallelic.expanded.sorted.PASS.vcf.gz
exac=${dbdir}/ExAC/ExAC.r0.3.1.sites.vep.multiallelic.expanded.v2.sorted.PASS.vcf.gz
exacnopass=${dbdir}/ExAC/ExAC.r0.3.1.sites.vep.multiallelic.expanded.v2.sorted.nonPASS.vcf.gz
target=${dbdir}/target/HaloPlex.b37.bed
logofile=${dbdir}/logo/WCM_2Line_CoBrand_RGB.jpg
#acmgfile=${dbdir}/GeneSet/ACMG_genes.txt
#brocafile=${dbdir}/GeneSet/BROCA_genes.txt
#cancerfile=${dbdir}/GeneSet/Census_allThu_Jun_16_21-46-18_2016.tsv
#veppath=/athena/ipm/scratch/users/taz2008/softwares/VEP/vep
#vepdata=/athena/ipm/scratch/users/taz2008/softwares/VEP/.vep
#veprefseq=${vepdata}/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
#vepexac=${vepdata}/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz

freqtablefile=${pipedir}/patch/info/processed.3.v2/frequency.no_genotype.no_FS.txt

####################### functions #######################
# check whether the previous command was successful
# if not, exit with a pre-defined error message
mycheck(){
	if [ $? -ne 0 ]
	then
		echo -e "\nError: $1"
		exit 99
	fi
}
#########################################################

echo "Start to process Sample ${sid}."

# check output folder
if [ -d ${vardir} ]
then
	echo "Output folder already exists, nothing is done."
	exit 0
else
	mkdir ${vardir}
	mkdir ${varlogdir}
	mkdir ${vartmpdir}
	mkdir ${anndir}
	mkdir ${annlogdir}
	mkdir ${anntmpdir}
	mkdir ${filtdir}
	mkdir ${filtlogdir}
	#mkdir ${cadddir}
	#mkdir ${caddlogdir}
	#mkdir ${mafdir}
	#mkdir ${maflogdir}
	####mkdir ${ipmannsuitedir}
	####mkdir ${ipmannsuitelogdir}
fi

# check source bam
echo -n "Checking source bam file..."
if [ -e ${bamfile} ]
then
        echo "ok."
else
        echo -e "\nCannot locate the source bam file: ${bamfile}"
	exit 1
fi

# copy source bam (and .bai) file to local directory
echo -n "Copying bam and preparing idx bai..."
cd ${vardir}
mycheck "Failed to change to directory ${vardir}"

bam=Sample_${sid}.bam
cp ${bamfile} ${bam}
mycheck "Failed to copy source bam ${bamfile} to the local directory `pwd`"

bamprefix=${bamfile%.bam}
if [ -e ${bamprefix}.bai ]
then
	cp ${bamprefix}.bai  Sample_${sid}.bai
elif [ -e ${bamfile}.bai ]
then
	cp ${bamfile}.bai  ${bam}.bai
else
	${samtools} index ${bam}
fi
mycheck "Failed to get index for bam ${bamfile}"
echo "ok."

# produce raw variant calls
echo -n "Make raw variant calls using GATK UnifiedGenotyper..."
${java} -Xmx${mem} -Djava.io.tmpdir=${vartmpdir} -jar ${gatk} -glm BOTH -R ${refseq} -T UnifiedGenotyper -D ${snpdb} -o ${vardir}/${sid}.UG.raw.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -minIndelFrac 0.1 -dcov 5000 -A AlleleBalance -A BaseCounts -A GCContent -A LowMQ -A SampleList -A VariantType -L ${target} -nt 8 -nct 4 \
-I ${bam} \
>${varlogdir}/${sid}.GATK.UnifiedGenotyper.log 2>&1
mycheck "Failed to make variant calls using GATK UnifiedGenotyper"
echo "ok."

# filter raw variant calls
echo -n "Filter raw variants using GATK VariantFiltration..."
${java} -Xmx${mem} -Djava.io.tmpdir=${vartmpdir} -jar ${gatk} -R ${refseq} -T VariantFiltration -o ${vardir}/${sid}.UG.filtered.vcf --variant ${vardir}/${sid}.UG.raw.vcf --clusterWindowSize 10 --clusterSize 3  --filterExpression "(DP - MQ0) < 10" --filterName "LowCoverage" \
>${varlogdir}/${sid}.GATK.VariantFiltration.log 2>&1
mycheck "Failed to pre-filter variant calls using GATK VariantFiltration"
echo "ok."

# annotate variants using snpeff
echo -n "Annotate variants using snpEff..."
${java7} -Xmx${mem} -Djava.io.tmpdir=${anntmpdir} -jar ${snpeff} ann -c ${config} -s ${anndir}/${sid}.summary.snpeff.UG.html ${flag} ${anndb} ${vardir}/${sid}.UG.filtered.vcf >${anndir}/${sid}.UG.filtered.snpeff.vcf 2>${annlogdir}/${sid}.snpeff.ann.log
mycheck "Failed to annotate variants using snpEff"
echo "ok."

# add clinVar annotations using snpsift
echo -n "Add ClinVar annotations using snpSift..."
${java7} -Xmx${mem} -Djava.io.tmpdir=${anntmpdir} -jar ${snpsift} annotate -v -noId -name clinvar -info ${CLNcols} ${clinvar} ${anndir}/${sid}.UG.filtered.snpeff.vcf >${anndir}/${sid}.UG.filtered.snpeff.clinvar.vcf 2>${annlogdir}/${sid}.snpsift.clinvar.log
mycheck "Failed to add ClinVar annotations"
echo "ok."

# add ExAC AF info using snpsift
echo -n "Add ExAC AF (allele frequency) annotations by snpSift..."
${java7} -Xmx${mem} -Djava.io.tmpdir=${anntmpdir} -jar ${snpsift} annotate -v -noId -name exac -info ${ExACcols} ${exac} ${anndir}/${sid}.UG.filtered.snpeff.clinvar.vcf >${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.vcf 2>${annlogdir}/${sid}.snpsift.ExAC.log
mycheck "Failed to add ExAC AF annotations"
echo "ok."

# add ExAC filterings using snpsift
echo -n "Add ExAC filterings by snpSift..."
${java7} -Xmx${mem} -Djava.io.tmpdir=${anntmpdir} -jar ${snpsift} annotate -v -noId -name exac -info ${ExACcols} ${exacnopass} ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.vcf >${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.extend.vcf 2>${annlogdir}/${sid}.snpsift.ExAC.extend.log
mycheck "Failed to add ExAC filterings"
echo "ok."
gzip ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.extend.vcf
mycheck "Failed to gzip vcf ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.extend.vcf"

# update FILTER column
echo -n "Set ExAC filtering..."
${python} ${srcdir}/set_ExAC_filter.py ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.extend.vcf.gz ${exacnopass} ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.patch.vcf.gz >${annlogdir}/${sid}.set_ExAC_filter.log
mycheck "Failed to set ExAC filtering in vcf"
echo "ok."

# extract info and screen
echo -n "Pre-categorize variants..."
${python} ${srcdir}/screen_germline_variants.py ${sid} ${anndir}/${sid}.UG.filtered.snpeff.clinvar.ExAC.patch.vcf.gz ${filtdir}/ ${filtdir}/${sid}.screen_variants.warning.txt >${filtlogdir}/${sid}.screen_germline_variants.log 2>&1
mycheck "Failed to pre-categorize variants"
echo "ok."

# select candidates
echo -n "Selecting candidate variants..."
${python} ${srcdir}/select_candidate_variants.v4.counts.no_FS.per_sample.py ${sid} ${filtdir}/${sid}.info.all.VEP.v2.txt.gz ${filtdir}/${sid}.candidates.v2.no_FS.txt.gz >${filtlogdir}/${sid}.select_candidate_variants.v4.counts.no_FS.per_sample.log 2>&1
mycheck "Failed to select candidate variants"
echo "ok."

# filter common variants in the IPM_cohort
${python} ${srcdir}/remove_artifects.no_FS.py ${filtdir}/${sid}.candidates.v2.no_FS.txt.gz ${freqtablefile} ${filtdir}/${sid}.candidates.v2.no_FS.filter_common.txt >${filtlogdir}/${sid}.remove_artifects.no_FS.log 2>&1
gzip ${filtdir}/${sid}.candidates.v2.no_FS.filter_common.txt

# search for potential false positive variants
${python} ${srcdir}/remove_artifects.mis_align.no_FS.py ${filtdir}/${sid}.candidates.v2.no_FS.filter_common.txt.gz ${filtdir}/${sid}.info.all.txt.gz ${filtdir}/${sid}.candidates.v2.no_FS.filter_common.clean.txt.gz ${filtdir}/${sid}.candidates.v2.no_FS.filter_common.manual_review.txt.gz >${filtlogdir}/${sid}.remove_artifects.mis_align.no_FS.log 2>&1

# cleaning
rm -rf ${vartmpdir}
rm -rf ${anntmpdir}
#rm ${vardir}/${bam}
rm ${anndir}/${sid}.UG.filtered.snpeff.vcf
gzip ${anndir}/${sid}.UG.filtered.snpeff.dbnsfp.vcf
rm ${vardir}/${sid}.UG.filtered.vcf
#gzip ${vardir}/${sid}*vcf
gzip ${anndir}/${sid}*vcf
#gzip ${anndir}/*.txt
gzip ${anndir}/*.html
gzip ${filtdir}/*.txt

echo "Complete processing Sample ${sid}."

