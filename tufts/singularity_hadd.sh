#!/bin/bash

CONTAINER=/cluster/kappa/90-days-archive/wongjiradlab/larbys/images/dllee_unified/singularity-dllee-unified-071017.img

LARBYSDIR=/cluster/kappa/90-days-archive/wongjiradlab
OUTPUTDIR=${LARBYSDIR}/larbys/data/mcc8.1/bnbnu_only_anatree/scraped
MERGEDFILE=${LARBYSDIR}/larbys/data/mcc8.1/bnbnu_only_anatree/anatree_scraped_merged.root

module load singularity
singularity exec ${CONTAINER} bash -c "source /usr/local/bin/thisroot.sh && hadd -f ${MERGEDFILE} ${OUTPUTDIR}/*.root"

