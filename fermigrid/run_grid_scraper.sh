#!/bin/bash

echo "PATH VARIABLES"
echo $PATH
echo $PATH > pathvars.txt

# Get JOB ID FILE
ifdh cp /pnfs/uboone/persistent/users/tmw/anatree_scraper/bnb_nu/jobids.txt jobids.txt

echo "copied jobids.txt"
tail -n 10 jobids.txt

# GET JOBID
let "proc_line=${PROCESS}+1"
let JOBID=`sed -n ${proc_line}p jobids.txt`

echo "proc_line: ${proc_line}"
echo "JOBID: ${JOBID}"

# GET INPUT LIST
inputlist_fullpath=`printf /pnfs/uboone/persistent/users/tmw/anatree_scraper/bnb_nu/inputlist/filelist_%03d.txt ${JOBID}`
inputlist_name=`printf filelist_%03d.txt ${JOBID}`
echo "Input lists:"
echo $inputlist_fullpath
echo $inputlist_name
ifdh cp ${inputlist_fullpath} ${inputlist_name}

# GET INPUT FILES
ifdh cp -f ${inputlist_name}

# CHANGE INPUT TO USE LOCAL PATHS
cat ${inputlist_name} | awk -F/ '{ print $NF }' > input_list.txt

# SCRAPE
${MRB_INSTALL}/uboonecode/v06_34_00/slf6.x86_64.e14.prof/bin/./scrape_anatree -l input_list.txt

# RENAME OUTPUT FILE
outfilename=`printf anatree_scraped_%03d.root ${JOBID}`
mv output_scrapedtree.root ${outfilename}
ls
ifdh cp ${outfilename} /pnfs/uboone/persistent/users/tmw/anatree_scraper/bnb_nu/scraped_anatree_out/${outfilename}

# CLEAN UP INPUT, SO IT DOESN'T TRANSFER
rm ana_hist_*.root



