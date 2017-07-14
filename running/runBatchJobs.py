#!/usr/bin/python
import os
import sys
import re
import argparse
import subprocess

parser = argparse.ArgumentParser(description='Prepare and submit ntupling jobs')
parser.add_argument("-m", "--macro",         dest="macro", default="runSomething.C", help="file to be run. [Default: runSomething.C]")
parser.add_argument("-b", "--runBatch",      dest="runBatch", default=False, help="Should we setup a condor job? [Default: False]")
parser.add_argument("-w", "--computeWeight", dest="compW", default=False, help="Should we include the parameters to calculate weight? [Default: False]")
parser.add_argument("-i", "--input",         dest="input", default="procdatasets.conf", help="input config or directory [Default: procdatasets.conf]")
parser.add_argument("-o", "--outputDir",     dest="outdir", default="", help="Output directory for ntuples. [Default: \"\"]")
parser.add_argument("-j", "--jobdir"       , dest="jobdir", default="jobs", help="Directory for job files  [Default: jobs]")
parser.add_argument("-r", "--runningDir"    , dest="runningDir", default="running/", help="Where to find helper files for running  [Default: running/]")
parser.add_argument("-s", "--runScript"    , dest="runScript", default="running/runBatchJob.sh", help="File that tells condor how to run  [Default: running/runBatchJob.sh]")
parser.add_argument("-n", "--numFiles",      dest="numFiles", default=5, help="Number of files per job if no config [Default: 5]")
parser.add_argument("-t", "--treeInt",       dest="treeInt" , default=1, help="treeInt if no config [Default: 1]")
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()


def prepareSampleJob(libName,outList, name, filelist, nFilesPerJob, treeInt, weightJob = False,  xsec = 1, numE = 1):
    nFiles = len(filelist)
    iF = 0
    iF2 = 0
    iJ = 0
    while iF < nFiles:
        if(iF2 == 0):
            jobfiles = open("{0}/files_{1}_{2}.txt".format(args.jobdir,name, iJ), "w")    
        if filelist[iF].startswith("/eos/uscms/store/user") :
            jobfiles.write("root://cmseos:1094/%s" % (re.match(r'/eos/uscms(.*)',filelist[iF]).group(1) ) )        
        else :
            jobfiles.write(filelist[iF])        
        jobfiles.write("\n")
        iF2 += 1
        iF += 1
        if(iF2 == int(nFilesPerJob)  or  iF ==  nFiles):
            jobfiles.close()
            
            inputF = "files_{0}_{1}.txt".format(name,iJ)
            outputF = "out_{0}_{1}.root".format(name,iJ)
            if not args.runBatch: 
                inputF  = os.path.normpath(os.path.join(os.path.join(os.getcwd(), args.jobdir),inputF))
                outputF = os.path.normpath(os.path.join(os.path.join(os.getcwd(), args.outdir),outputF))
            if weightJob  :
                CMD = "root -l -b -q \'{cfg}+(\"{INF}\",{TreeInt},\"{OUTF}\",{xs},{nE})\'".format( cfg=args.macro,INF=inputF,TreeInt=treeInt,OUTF=outputF,xs=xsec,nE=numE)
            else :
                CMD = "root -l -b -q \'{cfg}+(\"{INF}\",{TreeInt},\"{OUTF}\")\'".format( cfg=args.macro,INF=inputF,TreeInt=treeInt,OUTF=outputF)
            
            if args.runBatch :
                jobscript = open("{0}/submit_{1}_{2}.sh".format(args.jobdir,name,iJ), "w")
                jobscript.write("""
cat > submit.cmd <<EOF
universe                = vanilla
Requirements            = (Arch == "X86_64") && (OpSys == "LINUX")
Executable              = {runscript}
Arguments               = "{cmd}" {outputdir}
Output                  = logs/job_{samp}_{num}.out
Error                   = logs/job_{samp}_{num}.err
Log                     = logs/job_{samp}_{num}.log
use_x509userproxy       = true
initialdir              = {jobdir}
Should_Transfer_Files   = YES
transfer_input_files    = {cfgdir},{libdir},{workdir}/files_{samp}_{num}.txt
WhenToTransferOutput    = ON_EXIT
Queue
EOF

condor_submit submit.cmd;
rm submit.cmd""".format(
                runscript=script, cmd=CMD, outputdir=args.outdir,  workdir=os.path.normpath(os.path.join(os.getcwd(), args.jobdir)), num=iJ, jobdir=args.jobdir,
                cfgdir=os.path.normpath(os.path.join(os.getcwd(), args.macro)),libdir=os.path.normpath(os.path.join(os.getcwd(),  libName)),samp=name
                ))
                jobscript.close()
                os.system("chmod +x %s/submit_%s_%d.sh" % (args.jobdir, name, iJ))        
                outList.append("./{jobdir}/submit_{samp}_{num}.sh".format(jobdir=args.jobdir,samp=name,num=ijob))
            else :
                outList.append(CMD + " &")
            iF2 = 0
            iJ += 1

def getFileList(inputDir,name="") :    
    if name == "":
        grepCmd = "egrep '.*root'"
    else :
        grepCmd = "egrep '.*%s(-ext[0-9]*|)_[0-9]*.root'" % ( name)
        
    if inputDir.startswith("/eos/uscms/store/user") or inputDir.startswith("/store/user") :
        cmd = (" eos root://cmseos.fnal.gov find -f %s | %s" % ( args.input,grepCmd))
        prefix = "root://cmseos:1094/"
    else:
        cmd = ("find -f %s | %s" % (args.input,grepCmd))
        prefix = ""
    
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result = ps.communicate()
    return result[0].rstrip('\n').split('\n')                



os.system("mkdir -p %s" % args.jobdir)
os.system("mkdir -p %s/logs" % args.jobdir)
if args.outdir.startswith("/eos/cms/store/user") or args.outdir.startswith("/store/user") :
    os.system("eosmkdir -p %s" % (args.outdir))
else :
    os.system("mkdir -p %s" % args.outdir)
    
compCmd = ("root -l -b -q \'%s/compileMacro.C(\"%s\")\'" % (args.runningDir, args.macro))
ps = subprocess.Popen(compCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
print "Compiling macro:"
output = ps.communicate()            
libName = re.sub(r'(.+)\.(\w+)', r'\1_\2.so', args.macro)
if not os.path.isfile(libName) :
    print "could not find " + libName
    print output[0]
    print output[1]
    exit()            
else:
    print "Compiled!"

outList = []
if os.path.isdir(args.input) :
    fileList = getFileList(args.input,"")
    prepareSampleJob(libName,outList,"all", fileList, args.numFiles,args.treeInt)
else :
    inputData = open(args.inputData, "r")
    for line in inputData:
        if re.match("\s*#.*", line) : 
            continue
        match = re.match("^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$", line)
        if not match : 
            print "Do not understand:"
            print line
            continue
        fileList = getFileList(match.group(6),match.group(1))
        prepareSampleJob(libName,outList,match.group(1), fileList, match.group(5), match.group(2), 
                         args.compW, match.group(3),match.group(4))
        
if args.runBatch:
    subscript = open("submitall.sh", "w")
    subscript.write("#!/bin/bash")
    for line in outList :
        subscript.write("%s\n" & line)
else:
    for line in outList :
        print line

print "Done!"
exit()
