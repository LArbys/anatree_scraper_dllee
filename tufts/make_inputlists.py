import os,sys

# This script splits the input anatree filelists into chunks defined by the number of jobs
# Each file is about 6 MB. And each worker node is alloted only 40 GB of disk space.
# One should only aim for using 35 GB of space
# there a total of 20567 files, which amounts to a total of 123 GB of data
# scraped output files are in the 100 of kB, making the shrinking of this information worth it

ANATREE_SOURCEDIR="/cluster/kappa/90-days-archive/wongjiradlab/larbys/data/mcc8.1/bnbnu_only_anatree/anatree"

anafiles = os.listdir( ANATREE_SOURCEDIR )
analist = []
for f in anafiles:
    if ".root" not in f:
        continue
    analist.append( ANATREE_SOURCEDIR+"/"+f.strip() )

nfiles = len(analist)

num_jobs = 50
nfiles_per_job = nfiles/num_jobs
mem_per_job =  nfiles_per_job*6/1e3


print "files per job: ",nfiles_per_job
print "total file size per job: ",mem_per_job
if mem_per_job>35.0:
    print "too many files per job. increase number of jobs."
    sys.exit(-1)


jobid=0
nfiles_in_job = 0
os.system("mkdir -p inputlist")
jobflist = open("inputlist/filelist_%03d.txt"%(jobid),'w')
jobidlist = open("jobidlist.txt", 'w')
for fline in analist:
    fline = fline.strip()
    print >> jobflist,fline
    nfiles_in_job += 1
    if nfiles_in_job>=nfiles_per_job:
        jobflist.close()
        print >> jobidlist,jobid
        jobid += 1
        print "Making filelist for jobid #",jobid
        jobflist = open("inputlist/filelist_%03d.txt"%(jobid),'w')
        nfiles_in_job = 0
jobflist.close()
jobidlist.close()
