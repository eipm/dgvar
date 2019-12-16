#!/usr/bin/env sh
# eg.sh
# run a quick test on an example bam file
# 

workdir=`pwd`

if [ ! -d ${workdir}/../results ]
then
	mkdir ${workdir}/../results
fi

sh ${workdir}/../call_variants.sh  ${workdir}/test.bam  test  ${workdir}/test.bed

