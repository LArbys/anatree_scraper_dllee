#!/bin/bash
#
#SBATCH --job-name=anatree
#SBATCH --output=log_anatree.log
#SBATCH --ntasks=50
#SBATCH --time=30:00
#SBATCH --mem-per-cpu=4000


LARBYSDIR=/cluster/kappa/90-days-archive/wongjiradlab

CONTAINER=/cluster/kappa/90-days-archive/wongjiradlab/larbys/images/dllee_unified/singularity-dllee-unified-071017.img
WORKDIR=${LARBYSDIR}/twongj01/anatree_scraper_dllee/tufts
INPUTLISTDIR=${WORKDIR}/inputlist
JOBLIST=${WORKDIR}/jobidlist.txt
BINARY=${LARBYSDIR}/twongj01/anatree_scraper_dllee/scrape_anatree
OUTPUTDIR=${WORKDIR}/output/

mkdir -p ${OUTPUTDIR}
module load singularity
srun singularity exec ${CONTAINER} bash -c "cd ${WORKDIR} && source run_grid_scraper.sh ${WORKDIR} ${INPUTLISTDIR} ${JOBLIST} ${BINARY} ${OUTPUTDIR}"

