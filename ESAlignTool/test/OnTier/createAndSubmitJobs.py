#!/usr/bin/env python

import sys, os, shutil, re, subprocess, time
from optparse import OptionParser



cfi_template = """FileNames = [
INPUT_FILES
]
"""


bash_template = """#!/bin/bash
BATCHDIR=${PWD}
#export SCRAM_ARCH=slc5_amd64_gcc462
#export SCRAM_ARCH=slc6_amd64_gcc481
#export SCRAM_ARCH=slc6_amd64_gcc491
export SCRAM_ARCH=slc6_amd64_gcc530
cd MAIN_WORKDIR
eval `scram runtime -sh`
cp -v MAIN_WORKDIR/CMSSW_cfg.py $BATCHDIR/CMSSW_cfg.py
cp -v MAIN_WORKDIR/inputMatrixElements_cfi.py $BATCHDIR/
cp -v DATASET_WORKDIR/input/inputFiles_JOB_NUMBER_cfi.py $BATCHDIR/inputFiles_cfi.py
if [ -e MAIN_WORKDIR/inputDB.db ]
then
    cp -v MAIN_WORKDIR/inputDB.db $BATCHDIR/
fi
if [ -e MAIN_WORKDIR/myJSON.txt ]
then
    cp -v MAIN_WORKDIR/myJSON.txt $BATCHDIR/
fi
cd $BATCHDIR
echo "Running CMSSW job"
cmsRun CMSSW_cfg.py CFG_PARAMETERS OutFilename=OUTPUT_FILENAME_JOB_NUMBER.root >>DATASET_WORKDIR/output/OUTPUT_FILENAME_JOB_NUMBER.txt 2>DATASET_WORKDIR/output/OUTPUT_FILENAME_JOB_NUMBER.error.txt
exitcode=$?
 
echo "Copying file  OUTPUT_FILENAME.root to DATASET_WORKDIR/output/OUTPUT_FILENAME_JOB_NUMBER.root"
cp -v OUTPUT_FILENAME_JOB_NUMBER.root DATASET_WORKDIR/output/OUTPUT_FILENAME_JOB_NUMBER.root

exit $exitcode
"""


# usage description
usage = """Usage: ./createAndSubmitJobsOnCMSDATA.py [options]\n
Example: ./createAndSubmitJobsOnCMSDATA.py -w LXBatch_Jobs -d datasetList.txt -c btagvalidation_cfg.py\n
For more help: ./createAndSubmitJobsOnCMSDATA.py --help
"""

def main():
  # input parameters
  parser = OptionParser(usage=usage)

  parser.add_option("-w", "--main_workdir", dest="main_workdir", action='store', help="Main working directory", metavar="MAIN_WORKDIR")
 # parser.add_option("-d", "--dataset_list", dest="dataset_list", action='store', help="Text file containing a list of datasets to be processed", metavar="DATASET_LIST")
  parser.add_option("-d", "--dataset_list", dest="dataset_list", action='store', help="Text file containing a list of datasets to be processed", metavar="DATASET_LIST")
  parser.add_option("-l", "--file_list", dest="file_list", action='store', help="Text file containing a list of datasets to be processed", metavar="FILE_LIST")
  parser.add_option("-D", "--DBFileInput",  dest="DBFileInput",  action='store', help="Input DB file", default='', metavar="DBFILEINPUT")
  parser.add_option("-J", "--myJSONfile", dest="myJSONfile", action='store', help="JSON file for selecting event", default='', metavar="MYJSONFILE")
  parser.add_option("-o", "--output_filename", dest="output_filename", action='store', default='AlignmentFile', help="Output ROOT filename (Default set to AlignmentFile)", metavar="OUTPUT_FILENAME")
  parser.add_option("-E", "--eos_path", dest="eos_path", action='store', help="EOS path to copy output files to (This parameter is optional)", metavar="EOS_PATH")
  parser.add_option('-m', '--match', dest="match", action='store', help='Only files containing the MATCH string in their names will be considered (This parameter is optional)', metavar='MATCH')
  parser.add_option("-c", "--cmssw_cfg", dest="cmssw_cfg", action='store', help="CMSSW configuration file", metavar="CMSSW_CFG")
  parser.add_option('-f', '--fraction', dest='fraction', action='store', default='1.0', help='Fraction of files to be processed. Default value is 1 (This parameter is optional)', metavar='FRACTION')
  parser.add_option("-q", "--queue", dest="queue", action='store', default='1nh', help="LXBatch queue (choose among cmst3 8nm 1nh 8nh 1nd 1nw). Default is '1nh' (This parameter is optional)", metavar="QUEUE")
  parser.add_option("-n", "--no_submission", dest="no_submission", action='store_true', default=False, help="Create the necessary configuration files and skip the job submission (This parameter is optional)")

  (options, args) = parser.parse_args()

  # make sure all necessary input parameters are provided
  if not (options.main_workdir and options.dataset_list and options.cmssw_cfg):
    print usage
    sys.exit()

  main_workdir = options.main_workdir
  dataset_list = options.dataset_list
  cmssw_cfg = options.cmssw_cfg
  file_list = options.file_list

  # redefine main_workdir as an absolute path (if not defined in such form already)
  if not re.search("^/", main_workdir):
    main_workdir = os.path.join(os.getcwd(),main_workdir)

  # create the main working directory
  if not os.path.exists(main_workdir) : 
    os.mkdir(main_workdir)
  else :
    print 'Warning ' + main_workdir + ' exists. '

  # copy the dataset list file to the main_workdir
  shutil.copyfile(dataset_list,os.path.join(main_workdir,'datasetList.txt'))
  shutil.copyfile(file_list,os.path.join(main_workdir,'filelist.txt'))

  # copy the CMSSW cfg file to the cfg_files_dir
  shutil.copyfile(cmssw_cfg,os.path.join(main_workdir,'CMSSW_cfg.py'))

  if options.myJSONfile != '':
    shutil.copyfile( options.myJSONfile, os.path.join(main_workdir,'myJSON.txt'))
  if options.DBFileInput != '':
    shutil.copyfile( options.DBFileInput, os.path.join(main_workdir,'inputDB.db'))

  # look for pileup distribution files and copy them into main_workdir
  cfg_dirname = os.path.dirname(cmssw_cfg)
  if cfg_dirname=='':
    cfg_dirname = os.getcwd()
  for filename in os.listdir(cfg_dirname):
    if not os.path.isfile(os.path.join(cfg_dirname,filename)):
      continue
    if re.search("inputMatrixElements_cfi.py", filename):
      shutil.copy(os.path.join(cfg_dirname,filename),main_workdir)

  # open and read the dataset_list file
  dataset_list_file = open(dataset_list,"r")
  dataset_list_lines = dataset_list_file.readlines()

  # loop over datasets
  for line in dataset_list_lines:
    line_elements = line.split()
    if (len(line_elements)==0 or line_elements[0][0]=='#'): continue

    output_filename = options.output_filename
    cfg_parameters = ''
    if( len(line_elements)>3 ):
      cfg_parameters = line_elements[3]
      if( line_elements[3].split('=')[0]=='outFilename' ): output_filename = line_elements[3].split('=')[1].replace('.root','')
      for par in range(4,len(line_elements)):
        cfg_parameters = cfg_parameters + ' ' + line_elements[par]
        if( line_elements[par].split('=')[0]=='outFilename' ): output_filename = line_elements[par].split('=')[1].replace('.root','')

    dataset = line_elements[0].lstrip('/').replace('/','__')
    print 'Processing ' + line_elements[0]

    dataset_workdir = os.path.join(main_workdir,dataset)

    # create the dataset working directory
    os.mkdir(dataset_workdir)
    os.mkdir(os.path.join(dataset_workdir,'input'))
    os.mkdir(os.path.join(dataset_workdir,'output'))

    filelist = []
   # process_input_dir(line_elements[2], options.match, filelist)
    filelist_file = open(file_list,"r")
    filelist = filelist_file.readlines()


    ##################################################
    njobs = line_elements[1]
    numfiles = int(len(filelist)*float(options.fraction))
    ijobmax=int(njobs)
    if ijobmax > numfiles:
        ijobmax=numfiles
        print '  Number of jobs requested exceeds the total number of input files.\n  The number of jobs set to ' + str(ijobmax) + ' to match the number of input files'
    filesperjob = int(numfiles/ijobmax)
    if numfiles%ijobmax!=0:
        filesperjob = filesperjob + 1
        ijobmax = int(numfiles/filesperjob)
        if numfiles%filesperjob!=0:
            ijobmax = ijobmax + 1
        if ijobmax != int(njobs):
            print '  Could not create the exact number of jobs requested.\n  For optimal job splitting, the number of jobs set to ' + str(ijobmax)
    #################################################
    ifile = 0
    for ijob in range(ijobmax):
        # prepare the list of input files
        #input_files = '    \'root://eoscms.cern.ch/' + filelist[ifile] + '\''
        #input_files = '    \'root://se01.grid.nchc.org.tw//' + filelist[ifile] + '\''
        input_files = '    \'root://se01.grid.nchc.org.tw//' + filelist[ifile] 
        for i in range(filesperjob-1):
            if ifile>(numfiles-2):
                break
            ifile = ifile + 1
            #input_files += (',\n    \'root://eoscms.cern.ch/' + filelist[ifile] + '\'')
            #input_files += (',\n    \'root://se01.grid.nchc.org.tw//' + filelist[ifile] + '\'')
            input_files += (',\n    \'root://se01.grid.nchc.org.tw//' + filelist[ifile] )
        ifile = ifile + 1

        ## create input cfi file
        input_files_cfi = open(os.path.join(dataset_workdir,'input','inputFiles_' + str(ijob) + '_cfi.py'),'w')
        input_files_cfi.write(re.sub('INPUT_FILES',input_files,cfi_template))
        input_files_cfi.close()

        ## create Bash script
        bash_script = open(os.path.join(dataset_workdir,'input','job_' + str(ijob) + '.sh'),'w')
        bash_script_content = re.sub('MAIN_WORKDIR',main_workdir,bash_template)
        bash_script_content = re.sub('DATASET_WORKDIR',dataset_workdir,bash_script_content)
        bash_script_content = re.sub('JOB_NUMBER',str(ijob),bash_script_content)
        bash_script_content = re.sub('CFG_PARAMETERS',cfg_parameters,bash_script_content)
        bash_script_content = re.sub('OUTPUT_FILENAME',output_filename,bash_script_content)
#        bash_script_content = re.sub('EOS_PATH',eos_path,bash_script_content)
        bash_script.write(bash_script_content)
        bash_script.close()

        if(not options.no_submission):
          #time.sleep(2)
            cmd = 'qsub -q cms -V -N job_'+str(ijob) +' -o ' + os.path.join(dataset_workdir,'output','job_' + str(ijob) + '.log') +' '+ os.path.join(dataset_workdir,'input','job_' + str(ijob) + '.sh')
            os.system(cmd)
#            print cmd
  # close all open files
  dataset_list_file.close()


if __name__ == "__main__":
  main()
