import argparse
from subprocess import check_call
import os, sys, glob, shutil
from os.path import join as pjoin


parser = argparse.ArgumentParser()
parser.add_argument("--test-archive", type=str)
parser.add_argument("--install-prefix", type=str)
parser.add_argument("--test-source-dir", type=str)
args = parser.parse_args()

sys.path.insert(0,args.test_source_dir)

import util

run_dir = pjoin(os.getcwd(),"run")
if os.path.exists(run_dir):
    shutil.rmtree(run_dir,ignore_errors=True)
inp_dir = util.tar_extractall_safe_single_dir(args.test_archive,path=run_dir)
#inp_dir = pjoin(run_dir,"Bacillus.anthracis.sterne.PNNL.chunk.4sp")
out_dir = inp_dir+".out"
run_exe = pjoin(args.install_prefix,"bin","pgp_htc")
run_cmd = [run_exe,inp_dir,out_dir]

try:
    check_call(run_cmd)
    #pass
except:
    print >> sys.stderr, "Running the pipeline resulted in error"
    raise

out_dir_data = pjoin(out_dir,"run")
gff_dir_out = pjoin(out_dir_data,"DerivedData","GFFs")
assert os.path.isdir(gff_dir_out),"Output GFF dir does not exist"
gff_dir_exp = pjoin(inp_dir,"DerivedData.expected","GFFs")
if os.path.isdir(gff_dir_exp):
    gff_files_exp = glob.glob(pjoin(gff_dir_exp,"*.gff"))
    gff_files_exp = dict((os.path.basename(f),f) for f in gff_files_exp)
    gff_files_out = glob.glob(pjoin(gff_dir_out,"*.gff"))
    gff_files_out = dict((os.path.basename(f),f) for f in gff_files_out)
    assert set(gff_files_exp) == set(gff_files_out),"Different number of GFF files"
    for (f_exp,path_exp) in gff_files_exp.items():
        path_out = gff_files_out[f_exp]
        recs_exp = open(path_exp,"r").readlines()
        recs_out = open(path_out,"r").readlines()
        assert len(recs_exp) == len(recs_out),\
                "Different # of lines in {} ({}) and {} ({}).".\
                format(path_exp,len(recs_exp),path_out,len(recs_out))
        print "Found {} output GFF records".format(len(recs_out))
print
print "Pipeline command was: {}".format(" ".join(run_cmd))
print
print "Pipeline output is in: {}".format(out_dir_data)
print

