Fermigrid scripts for running scrape_anatree

* build the binary
* for convenience, make a link to the larsoft code. I'll assume the link is named `larsoft`
* copy the binary into `larsoft/localProducts_larsoft_vXX_YY_ZZ_eNN_prof/uboonecode/vXX_YY_ZZ/slf6.x86_64.eNN.prof/bin/`
* in the larsoft directory, run `make_tar_uboone.sh larsoft.tar`
* make an output directory for the files on `/pnfs`
* get the list of anatree files from samweb and put them into a file: `samweb list-definition-files [anatree def] | sort > samplefilenames.txt`
* get the location of the one of the files: `samweb locate-file [filename]`
* make another list with the full path: `cat samplefilenames.txt | awk '{ print [anatree dir]/$0 }' > anatree_pathlist.txt`
* use make_inputlists.py to make input filelists. be sure to adjust name of text file with path lists and/or the number of jobs
* make a jobids.txt: `seq 0 N`, where N is the number of jobs that matches the value in make_inputlists.py
* copy over the `inputlist` files onto pnfs: `mkdir [pnfs folder]/inputlist`; ifdh cp -r inputlist/* [pnfs folder]/inputlist/`
* copy over the `jobids.txt` file onto pngs: `ifdh cp jobids.txt [pnfs folder]/jobids.txt`
* edit `run_grid_scraper.sh` to have the right locations of the 1) jobids.txt file 2) the inputlist dir and 3) output folder for the scraped root files
* double-check the last step. the calls to ifdh better be right, else the job will spin it's wheels for HOURS.
* edit the project xml file
* launch a test job by setting <numjobs>1</numjobs> in the project xml
* if fails, use `jobsub_fetchlog --jobid=` to understand what happened. You can get the jobid from `[pnfs folder]/log/vXX_YY_ZZ/jobids.list`
* if succeed, change the <numjobs>N</numjobs> in the project xml
* Note: the jobs are quick, 5-10 mins. If longer, an ifdh call is probably WRONG. Unfortunately, you need for the calls to time out to get error logs

