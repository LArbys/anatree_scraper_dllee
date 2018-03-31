#!/bin/bash

source /usr/local/bin/thisroot.sh

workdir=$1
inputlistdir=$2
joblist=$3
binary=$4
outputdir=$5

cd ${workdir}

# GET JOBID
let "proc_line=${SLURM_PROCID}+1"
let JOBID=`sed -n ${proc_line}p ${joblist}`

echo "proc_line: ${proc_line}"
echo "JOBID: ${JOBID}"

# GET INPUT LIST
inputlist_fullpath=`printf %s/filelist_%03d.txt ${inputlistdir} ${JOBID}`
echo $inputlist_fullpath

# make jobid
jobdir=`printf slurm_job%04d ${JOBID}`
mkdir ${jobdir}
cd ${jobdir}

# SCRAPE
export PATH=.:${PATH}
cp ${binary} .
./scrape_anatree -l ${inputlist_fullpath}


# RENAME OUTPUT FILE
outfilename=`printf anatree_scraped_%03d.root ${JOBID}`
mv output_scrapedtree.root ${outfilename}
cp ${outfilename} ${outputdir}/${outfilename}

cd ../
rm -r ${jobdir}



