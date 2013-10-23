import argparse
from subprocess import check_call
import os, sys, glob, shutil
from os.path import join as pjoin
import ConfigParser

parser = argparse.ArgumentParser(description="""Execute entire PGP pipeline, possibly 
running different stages with different cluster job requirements or locally on
submit node.""",
epilog="""You can also include any combination of Makeflow options, which will be passed directly to the Makeflow.
 In particular, if running on a cluster, you will need to provide the Makeflow
 arguments described below.""")

#parser.parse_args(['--', '-f'])
parser.add_argument("--pgp-config-file", type=os.path.abspath, 
        help="PGP config file")
parser.add_argument("input_dir", type=os.path.abspath, 
        help="Input structured data directory")
parser.add_argument("output_dir", type=os.path.abspath, 
        help="Output and work directory")
parser.add_argument("--pgp-top-makefile", type=os.path.abspath, 
        help="Top level Makeflow makefile")
parser.add_argument("--do-not-run", action="store_true", 
        help="Do not run the pipeline, only generate the workflow")
parser.add_argument("-T","--batch-type", type=str, 
        help="Makeflow: Batch system type [%(default)s]",
        default="local")
parser.add_argument("-B","--batch-options", type=str, 
        help="Makeflow: Add these options to all batch submit files. Quote in '' if spaces are present inside.")
parser.add_argument("--pgp-batch-options-extra-prepare", type=str, default="",
        help="Makeflow: Add these options to --batch-options for 'prepare' (serial) stage. Quote in '' if spaces are present inside.")
parser.add_argument("--pgp-batch-options-extra-process", type=str, default="",
        help="Makeflow: Add these options to --batch-options for 'process' (parallel) stage. Quote in '' if spaces are present inside.")
parser.add_argument("--pgp-local-prepare", action="store_true", 
        help="Run the 'prepare' stage on a local machine regardless of --batch-type setting")
parser.add_argument("--pgp-process-mpi", action="store_true", 
        help="Use MPI backend for the 'process' stage")
args,makeflow_args_other = parser.parse_known_args()

#print args
#print makeflow_args_other
    
assert os.path.isfile(args.pgp_config_file), "Config file does not exist: {}".format(args.config_file)

config = ConfigParser.SafeConfigParser()
config.read(args.pgp_config_file)

ini_section = "main"

default_top_makefile = config.get(ini_section,"make_top_in")
makeflow_exe = config.get(ini_section,"makeflow")
env_file = config.get(ini_section,"task_env")
command_file = "this_pgp_run.sh"

if not args.pgp_top_makefile:
    args.pgp_top_makefile = default_top_makefile

assert not os.path.exists(args.output_dir),\
    "Output directory should not already exist: {}".format(args.output_dir)

os.makedirs(args.output_dir)

top_makefile_out = pjoin(args.output_dir,"pgp_run.mkf")

shutil.copy(args.pgp_top_makefile,top_makefile_out)

top_makefile_out = os.path.basename(top_makefile_out)

run_cmd = [makeflow_exe] + makeflow_args_other + [top_makefile_out]

os.chdir(args.output_dir)

def _evar(out,name,val):
    out.write('export {}="{}"\n'.format(name,val))

with open(command_file,"w") as out:
    out.write("""#!/bin/bash
source {env_file}
""".format(env_file=env_file))
    _evar(out,"PGP_INPUT_DIR",args.input_dir)
    _evar(out,"PGP_OUTPUT_DIR",args.output_dir)
    _evar(out,"MAKEFLOW_BATCH_QUEUE_TYPE",args.batch_type)
    if not args.batch_options:
        batch_options = ""
    else:
        batch_options = args.batch_options
    _evar(out,"BATCH_OPTIONS",batch_options)
    _evar(out,"PGP_BATCH_OPTIONS_PREPARE", " ".join((batch_options,args.pgp_batch_options_extra_prepare)).strip())
    _evar(out,"PGP_BATCH_OPTIONS_PROCESS", " ".join((batch_options,args.pgp_batch_options_extra_process)).strip())
    _evar(out,"PGP_BATCH_LOCAL_PREPARE", 1 if args.pgp_local_prepare else 0)
    if args.pgp_process_mpi:
        local_nested = 0
        wrapper_nested = config.get(ini_section,"wrapper_nested_mpi")
        assert os.path.isfile(wrapper_nested),\
                "Wrapper script for nested MPI makeflows does not exist: {}".format(wrapper_nested)
    else:
        local_nested = 1
        wrapper_nested = ""
    _evar(out,"PGP_BATCH_LOCAL_PROCESS", local_nested)
    _evar(out,"PGP_NESTED_WRAPPER_PROCESS", wrapper_nested)

    out.write(" ".join(['"{}"'.format(c) for c in run_cmd])+"\n")

os.chmod(command_file,os.stat(command_file).st_mode|0755)

print "PGP: Makeflow command is saved to a file if you want to execute it later as ./{}".\
    format(command_file)

if not args.do_not_run:
    print "PGP: starting the run"
    check_call(["bash",command_file])

