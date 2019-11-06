#!/usr/bin/env sh
# IPM_germline_pipeline.dev.sh
# IPM germline pipeline (re-designed)
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

####bamfile=/zenodotus/elementolab/scratch/kwe2001/PM_Cases/Sample_PM198_EBC1_1_Ctrl_HALO/Sample_PM198_EBC1_1_Ctrl_HALO.md.filt.bam
bamfile=$1
####sid=PM198_EBC1_1_Ctrl_HALO
sid=$2

# folders
basedir=/athena/ipm/scratch/users/taz2008/redteam2
pipedir=${basedir}/germline_variant_v2
workdir=${basedir}/request/Bishoy/171101/updated/info
srcdir=${workdir}/../code
codedir=${pipedir}/standalone/src
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
mafdir=${workdir}/../../info/maf
maflogdir=${mafdir}/logs
####ipmannsuitedir=${vardir}/ipmannsuite
####ipmannsuitelogdir=${ipmannsuitedir}/logs

# tools
samtools=${tooldir}/samtools-0.1.19/samtools
java=/softlib/exe/x86_64/bin/java
java7=/home/taz2008/softwares/jre1.7.0_51/bin/java
python=/home/taz2008/softwares/Python-2.7.6/python
perl=/home/taz2008/perl/perl-5.22.2/bin/perl
gatk=${tooldir}/gatk/2.5.2/GenomeAnalysisTK.jar
snpeff=/home/taz2008/softwares/snpEff_v4.2/snpEff.jar
snpsift=/home/taz2008/softwares/snpEff_v4.2/SnpSift.jar
vcf2maf=/athena/ipm/scratch/users/taz2008/softwares/mskcc-vcf2maf-2f82fa4/vcf2maf.pl
####ipmannsuite=/home/kwe2001/scratch/IPM_Annotation_Suite.pl

# configuration
config=/home/taz2008/softwares/snpEff_v4.2/snpEff.config
mem=10g
anndb=GRCh37.75
#flag="-no-upstream -no-downstream -no-intergenic -formatEff -v -o gatk"
flag="-no-upstream -no-downstream -no-intergenic -v"
dbNSFPcols="genename,Uniprot_acc,Uniprot_id,Uniprot_aapos,cds_strand,Ensembl_geneid,Ensembl_transcriptid,aapos_SIFT,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,MutationTaster_score,MutationTaster_pred,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,Reliability_index,ExAC_AC,ExAC_AF,clinvar_rs,clinvar_clnsig,clinvar_trait"
ExACcols="AF,AdjAF,AN_Adj,AF_AFR,AF_AMR,AF_EAS,AF_FIN,AF_NFE,AF_SAS,AF_OTH,AF_MALE,AF_FEMALE"
CLNcols="CLNVARID,CLNDN,CLNDISDB,CLNREVSTAT,CLNSIG,CLNSIGCONF,CLNVI,DBVARID,ORIGIN,GENEINFO,MC,RS,SSR"

# databases
genomedir=${dbdir}/genome
refseq=${genomedir}/human_g1k_b37.fasta
snpdb=${dbdir}/dbsnp/dbsnp_137.b37.vcf
####clinvar=${dbdir}/clinvar/clinvar_20160531.vcf
####clinvar=${dbdir}/clinvar/clinvar_20160531.multiallelic.expanded.sorted.vcf.gz
clinvar=${dbdir}/clinvar.20180805/clinvar_20180805.fix.vcf.gz
dbnsfp=${dbdir}/dbNSFP_v2.9.1/dbNSFP2.9.1.txt.gz
####exac=${dbdir}/ExAC/ExAC.r0.3.1.sites.vep.vcf.gz
####exac=${dbdir}/ExAC/ExAC.r0.3.1.sites.vep.multiallelic.expanded.sorted.PASS.vcf.gz
exac=${dbdir}/ExAC/ExAC.r0.3.1.sites.vep.multiallelic.expanded.v2.sorted.PASS.vcf.gz
exacnopass=${dbdir}/ExAC/ExAC.r0.3.1.sites.vep.multiallelic.expanded.v2.sorted.nonPASS.vcf.gz
target=${dbdir}/target/HaloPlex.b37.bed
logofile=${dbdir}/logo/WCM_2Line_CoBrand_RGB.jpg
acmgfile=${dbdir}/GeneSet/ACMG_genes.txt
brocafile=${dbdir}/GeneSet/BROCA_genes.txt
cancerfile=${dbdir}/GeneSet/Census_allThu_Jun_16_21-46-18_2016.tsv
veppath=/athena/ipm/scratch/users/taz2008/softwares/VEP/vep
vepdata=/athena/ipm/scratch/users/taz2008/softwares/VEP/.vep
veprefseq=${vepdata}/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
vepexac=${vepdata}/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
freqtablefile=${pipedir}/patch/info/processed.3.v2/frequency.no_genotype.no_FS.txt

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
:<<TEST
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
#:<<TEST
bam=Sample_${sid}.bam
cp ${bamfile} ${bam}
if [ $? -ne 0 ]
then
	echo -e "\nFailed to copy source bam ${bamfile} to the local directory `pwd`"
	exit 2
fi
${samtools} index ${bam}
echo "ok."

# produce raw variant calls
echo -n "Make raw variant calls by GATK UnifiedGenotyper..."
${java} -Xmx${mem} -Djava.io.tmpdir=${vartmpdir} -jar ${gatk} -glm BOTH -R ${refseq} -T UnifiedGenotyper -D ${snpdb} -o ${vardir}/${sid}.UG.raw.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -minIndelFrac 0.1 -dcov 5000 -A AlleleBalance -A BaseCounts -A GCContent -A LowMQ -A SampleList -A VariantType -L ${target} -nt 8 -nct 4 \
-I ${bam} \
>${varlogdir}/${sid}.GATK.UnifiedGenotyper.log 2>&1
echo "ok."

# filter raw variant calls
echo -n "Filter raw variants by GATK VariantFiltration..."
${java} -Xmx${mem} -Djava.io.tmpdir=${vartmpdir} -jar ${gatk} -R ${refseq} -T VariantFiltration -o ${vardir}/${sid}.UG.filtered.vcf --variant ${vardir}/${sid}.UG.raw.vcf --clusterWindowSize 10 --clusterSize 3  --filterExpression "(DP - MQ0) < 10" --filterName "LowCoverage" \
>${varlogdir}/${sid}.GATK.VariantFiltration.log 2>&1
echo "ok."

# evaluate variant calls
#echo -n "Evaluate variants by GATK VariantEval..."
#${java} -Xmx${mem} -Djava.io.tmpdir=${vartmpdir} -jar ${gatk} -R ${refseq} -T VariantEval -D ${snpdb} -o ${vardir}/${sid}.UG.filtered.eval.gatkreport -eval ${vardir}/${sid}.UG.filtered.vcf -nt 4 -L ${target} \
#>${varlogdir}/${sid}.GATK.VariantEval.log 2>&1
#echo "ok."
TEST

# cp vcf files
oldvcffile2=${pipedir}/variants_pipeline_result/${sid}/${sid}.UG.filtered.vcf.gz
if [ ! -e ${oldvcffile2} ]
then
        echo "Cannot find vcf file: ${oldvcffile2}"
        exit 100
fi
zcat ${oldvcffile2} >${vardir}/${sid}.UG.filtered.vcf

# annotate variants using snpeff
echo -n "Annotate variants by snpEff..."
${java7} -Xmx${mem} -Djava.io.tmpdir=${anntmpdir} -jar ${snpeff} ann -c ${config} -s ${anndir}/${sid}.summary.snpeff.UG.html ${flag} ${anndb} ${vardir}/${sid}.UG.filtered.vcf >${anndir}/${sid}.UG.filtered.snpeff.vcf 2>${annlogdir}/${sid}.snpeff.ann.log
echo "ok."

# add dbNSFP annotations using snpsift
echo -n "Add dbNSFP annotations by snpSift..."
${java7} -Xmx${mem} -Djava.io.tmpdir=${anntmpdir} -jar ${snpsift} dbnsfp -v -db ${dbnsfp} -f ${dbNSFPcols} ${anndir}/${sid}.UG.filtered.snpeff.vcf >${anndir}/${sid}.UG.filtered.snpeff.dbnsfp.vcf 2>${annlogdir}/${sid}.snpsift.dbnsfp.log
echo "ok."

:<<TEST2
# cp vcf files
oldvcffile1=${pipedir}/variants_pipeline_result/${sid}/snpeff_gatk_v2/${sid}.UG.filtered.snpeff.dbnsfp.vcf.gz
oldvcffile2=${pipedir}/variants_pipeline_result/${sid}/${sid}.UG.filtered.vcf.gz
if [ ! -e ${oldvcffile1} ]
then
	echo "Cannot find vcf file: ${oldvcffile1}"
	exit 99
fi
if [ ! -e ${oldvcffile2} ]
then
	echo "Cannot find vcf file: ${oldvcffile2}"
	exit 100
fi
zcat ${oldvcffile1} >${anndir}/${sid}.UG.filtered.snpeff.dbnsfp.vcf
zcat ${oldvcffile2} >${vardir}/${sid}.UG.filtered.vcf
TEST2

# add clinVar annotations using snpsift
echo -n "Add ClinVar annotations by snpSift..."
${java7} -Xmx${mem} -Djava.io.tmpdir=${anntmpdir} -jar ${snpsift} annotate -v -noId -name clinvar -info ${CLNcols} ${clinvar} ${anndir}/${sid}.UG.filtered.snpeff.dbnsfp.vcf >${anndir}/${sid}.UG.filtered.snpeff.dbnsfp.clinvar.vcf 2>${annlogdir}/${sid}.snpsift.clinvar.log
echo "ok."

# add ExAC AF info using snpsift
echo -n "Add ExAC AF annotations by snpSift..."
${java7} -Xmx${mem} -Djava.io.tmpdir=${anntmpdir} -jar ${snpsift} annotate -v -noId -name exac -info ${ExACcols} ${exac} ${anndir}/${sid}.UG.filtered.snpeff.dbnsfp.clinvar.vcf >${anndir}/${sid}.UG.filtered.snpeff.dbnsfp.clinvar.ExAC.vcf 2>${annlogdir}/${sid}.snpsift.ExAC.log
echo "ok."

# add ExAC AF info using snpsift
echo -n "Add ExAC AF no-Pass annotations by snpSift..."
${java7} -Xmx${mem} -Djava.io.tmpdir=${anntmpdir} -jar ${snpsift} annotate -v -noId -name exac -info ${ExACcols} ${exacnopass} ${anndir}/${sid}.UG.filtered.snpeff.dbnsfp.clinvar.ExAC.vcf >${anndir}/${sid}.UG.filtered.snpeff.dbnsfp.clinvar.ExAC.extend.vcf 2>${annlogdir}/${sid}.snpsift.ExAC.extend.log
echo "ok."
gzip ${anndir}/${sid}.UG.filtered.snpeff.dbnsfp.clinvar.ExAC.extend.vcf

# update FILTER column
echo -n "Set ExAC filter..."
${python} ${srcdir}/set_ExAC_filter.py ${anndir}/${sid}.UG.filtered.snpeff.dbnsfp.clinvar.ExAC.extend.vcf.gz ${exacnopass} ${anndir}/${sid}.UG.filtered.snpeff.dbnsfp.clinvar.ExAC.patch.vcf.gz >${annlogdir}/${sid}.set_ExAC_filter.log
echo "ok."

# extract info and screen
echo -n "Categorize variants..."
${srcdir}/screen_germline_variants.v5.py ${sid} ${anndir}/${sid}.UG.filtered.snpeff.dbnsfp.clinvar.ExAC.patch.vcf.gz ${filtdir}/ ${filtdir}/${sid}.screen_variants.warning.txt ${acmgfile} ${brocafile} ${cancerfile} >${filtlogdir}/${sid}.screen_germline_variants.v5.log 2>&1
echo "ok."

#### get CADD for variants
####echo -n "Requesting for CADD scores..."
####${python} ${srcdir}/simplify_vcf.py ${vardir}/${sid}.UG.filtered.vcf ${cadddir}/${sid}.for_CADD.vcf >${caddlogdir}/simplify_vcf.${sid}.log 2>&1
####${python} ${srcdir}/run_CADD.v2.py ${cadddir}/${sid}.for_CADD.vcf ${cadddir}/${sid}.webpage.html.gz ${cadddir}/${sid}.CADD.gz >${caddlogdir}/run_CADD.${sid}.log 2>&1
###### !!! need to wait until server generates results
####${python} ${srcdir}/download_CADD.py ${cadddir}/${sid}.webpage.html.gz ${cadddir}/${sid}.CADD.gz >${caddlogdir}/download_CADD.${sid}.log 2>&1
####sleep 20
####${python} ${srcdir}/add_CADD.py ${filtdir}/${sid}.info.all.txt ${cadddir}/${sid}.CADD.gz ${filtdir}/${sid}.info.all.CADD.v2.txt.gz >${filtlogdir}/${sid}.add_CADD.log 2>&1
####echo "ok."

# get MAF for variants
echo -n "Generate MAF..."
${perl} ${vcf2maf}  --input-vcf  ${vardir}/${sid}.UG.filtered.vcf  --output-maf  ${mafdir}/${sid}.vep.maf --vep-path ${veppath}  --vep-data ${vepdata}  --ref-fasta ${veprefseq}  --filter-vcf ${vepexac}  --vcf-normal-id ${sid} >${maflogdir}/vcf2maf.${sid}.log 2>&1
####gzip ${mafdir}/${sid}.UG.filtered.vep.vcf
gzip ${filtdir}/${sid}.info.all.txt
### !!! wait for CADD
###${python} ${srcdir}/add_VEP.py ${filtdir}/${sid}.info.all.CADD.v2.txt.gz ${mafdir}/${sid}.vep.maf ${filtdir}/${sid}.info.all.CADD.VEP.v2.txt.gz >${filtlogdir}/${sid}.add_VEP.log 2>&1
${python} ${srcdir}/add_VEP.py ${filtdir}/${sid}.info.all.txt.gz ${mafdir}/${sid}.vep.maf ${filtdir}/${sid}.info.all.VEP.v2.txt.gz >${filtlogdir}/${sid}.add_VEP.log 2>&1
echo "ok."

# select candidates
${python} ${srcdir}/select_candidate_variants.v4.counts.no_FS.per_sample.py ${sid} ${filtdir}/${sid}.info.all.VEP.v2.txt.gz ${filtdir}/${sid}.candidates.v2.no_FS.txt.gz >${filtlogdir}/${sid}.select_candidate_variants.v4.counts.no_FS.per_sample.log 2>&1

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

