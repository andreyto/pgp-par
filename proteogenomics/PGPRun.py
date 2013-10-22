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

if args.batch_type == "local":
    nested_makeflow_command = "MAKEFLOW"
else:
    nested_makeflow_command = "{} --batch-type {}".\
            format(makeflow_exe,args.batch_type)
    if args.batch_options:
        assert not "'" in args.batch_options, \
                "Use double quotes inside of --batch-options value instead of single quotes"
        nested_makeflow_command += " --batch-options '{}'".\
                format(args.batch_options)

with open(args.pgp_top_makefile,"r") as inp:
    text = inp.read()
text = text.replace("##PGP_INPUT_DIR##",args.input_dir)
text = text.replace("##PGP_OUTPUT_DIR##",args.output_dir)
text = text.replace("##PGP_NESTED_MAKEFLOW_COMMAND##",nested_makeflow_command)

top_makefile_out = pjoin(args.output_dir,"pgp_run.mkf")

with open(top_makefile_out,"w") as out:
    out.write(text)

top_makefile_out = os.path.basename(top_makefile_out)

makeflow_args = ["--batch-type",args.batch_type] + makeflow_args_other
if args.batch_options:
    makeflow_args += ["--batch-options",args.batch_options]

run_cmd = [makeflow_exe] + makeflow_args + [top_makefile_out]

os.chdir(args.output_dir)
with open(command_file,"w") as out:
    out.write("""#!/bin/bash
source {env_file}
{run_cmd}""".format(env_file=env_file,run_cmd=" ".join(['"{}"'.format(c) for c in run_cmd])))
#os.chmod(command_file,0x755)
os.chmod(command_file,os.stat(command_file).st_mode|0755)

print "PGP: Makeflow command is saved to a file if you want to execute it later as ./{}".\
    format(command_file)

if not args.do_not_run:
    print "PGP: starting the run"
    check_call(run_cmd)

