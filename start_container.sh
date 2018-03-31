#!/bin/bash

cd /cluster/kappa/90-days-archive/wongjiradlab/twongj01/anatree_scraper_dllee/

CONTAINER=/cluster/kappa/90-days-archive/wongjiradlab/larbys/images/dllee_unified/singularity-dllee-unified-071017.img

module load singularity
singularity shell ${CONTAINER}

