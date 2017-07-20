#!/bin/bash
#
#SBATCH --job-name=anatree
#SBATCH --output=log_anatree.log
#SBATCH --ntasks=50
#SBATCH --time=20:00
#SBATCH --mem-per-cpu=4000


LARBYSDIR=/cluster/kappa/90-days-archive/wongjiradlab

CONTAINER=/cluster/kappa/90-days-archive/wongjiradlab/larbys/images/dllee_unified/singularity-dllee-unified-071017.img
WORKDIR=${LARBYSDIR}/grid_jobs/anatree_scraper_dllee/tufts
INPUTLISTDIR=${WORKDIR}/inputlist
JOBLIST=${WORKDIR}/jobidlist.txt
BINARY=${LARBYSDIR}/grid_jobs/anatree_scraper_dllee/scrape_anatree
OUTPUTDIR=${LARBYSDIR}/larbys/data/mcc8.1/bnbnu_only_anatree/scraped

mkdir -p ${OUTPUTDIR}
module load singularity
srun singularity exec ${CONTAINER} bash -c "cd ${WORKDIR} && source run_grid_scraper.sh ${WORKDIR} ${INPUTLISTDIR} ${JOBLIST} ${BINARY} ${OUTPUTDIR}"

