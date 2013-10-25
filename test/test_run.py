import argparse
from subprocess import check_call
import os, sys, glob, shutil, shlex
from os.path import join as pjoin


parser = argparse.ArgumentParser(description="Test driver program.",
        epilog="Any remaining options will be passed directly to pgp_run program.")
parser.add_argument("test_archive", type=os.path.abspath, 
    help="Path to structured input directory. "+\
            "It can be either archive in 'tar.gz' format or an actual directory")
parser.add_argument("install_prefix", type=os.path.abspath,
    help="Path to --prefix used when installing the package")

args,args_other = parser.parse_known_args()

if len(args_other) == 1:
    args_other = shlex.split(args_other[0])

sys.path.insert(0,os.path.dirname(os.path.abspath(sys.argv[0])))

import util

run_dir = pjoin(os.getcwd(),"run_"+os.path.basename(args.test_archive))
if os.path.exists(run_dir):
    shutil.rmtree(run_dir,ignore_errors=True)
os.makedirs(run_dir)
if os.path.isfile(args.test_archive):
    inp_dir = util.tar_extractall_safe_single_dir(args.test_archive,path=run_dir)
    out_dir = inp_dir+".out"
elif os.path.isdir(args.test_archive):
    inp_dir = args.test_archive
    out_dir = pjoin(run_dir,os.path.basename(inp_dir)+".out")
else:
    raise ValueError("test_archive is neither file nor directory: {}".\
            format(args.test_archive))
#inp_dir = pjoin(run_dir,"Bacillus.anthracis.sterne.PNNL.chunk.4sp")
run_exe = pjoin(args.install_prefix,"bin","pgp_run")
run_cmd = [run_exe,inp_dir,out_dir] + args_other
#run_cmd = [run_exe,inp_dir,out_dir, "-T", "sge", "-B", '-l fast -P 0534 -b n -S /bin/bash']

run_cmd_str = " ".join([ '"{}"'.format(c) for c in run_cmd ])

print "Pipeline command is: {}".format(run_cmd_str)

try:
    check_call(run_cmd)
    #pass
except:
    print >> sys.stderr, "Running the pipeline resulted in error"
    raise

out_dir_data = out_dir
gff_dir_out = pjoin(out_dir_data,"DerivedData","GFFs")
assert os.path.isdir(gff_dir_out),"Output GFF dir does not exist"
assert glob.glob(pjoin(out_dir_data,"DerivedData","*.txt")),\
        "Output report files were not created"
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
        len_recs_exp = len(recs_exp)
        len_recs_out = len(recs_out)
        max_len_recs = max(len_recs_exp,len_recs_out)
        #we give 5% slack for numerical stability across runs and architectures
        #when comparing the number of annotated peptides with the reference run
        assert max_len_recs == 0 or \
                float(abs(len_recs_exp - len_recs_out))/max_len_recs <= 0.1,\
                "Substantially different # of peptides in {} ({}) and {} ({}).".\
                format(path_exp,len_recs_exp,path_out,len_recs_out)
        print "Found {} output GFF records".format(len_recs_out)
print
print "Pipeline command was: {}".format(run_cmd_str)
print
print "Pipeline output is in: {}".format(out_dir_data)
print

