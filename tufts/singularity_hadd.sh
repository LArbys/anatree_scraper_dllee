#!/bin/bash

CONTAINER=/cluster/kappa/90-days-archive/wongjiradlab/larbys/images/dllee_unified/singularity-dllee-unified-080717.img

LARBYSDIR=/cluster/kappa/90-days-archive/wongjiradlab
OUTPUTDIR=${LARBYSDIR}/larbys/data/comparison_samples/1e1p/out_week080717/scraped_anatree
MERGEDFILE=${LARBYSDIR}/larbys/data/comparison_samples/1e1p/out_week080717/scrapedana_merged_comparison_1e1p.root

module load singularity
singularity exec ${CONTAINER} bash -c "source /usr/local/bin/thisroot.sh && hadd -f ${MERGEDFILE} ${OUTPUTDIR}/*.root"

