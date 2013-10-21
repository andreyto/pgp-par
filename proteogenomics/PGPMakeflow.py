### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   Copyright J. Craig Venter Institute 2012
#   See docs.PGP/COPYING file distributed along with the proteogenomics 
#   package for the copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

"""Implementation of the PGP workflow that uses CCTools Makeflow.
Makeflow can be found here: http://nd.edu/~ccl/software/makeflow/
"""

import ConfigParser, exceptions, logging, optparse, os, re, string, sys, glob

from MGT.Util import *

def basedir(path):
    """Return right-most component of the directory path.
    This gives identical results both for /some/base and /some/base/"""
    dirn,basen = os.path.split(path)
    if not basen:
        basen = os.path.basename(dirn)
    return basen


def assert_log(cond,msg,logger):
    """If condition is False, log error and raise AssertionError"""
    if not cond:
        logger.error(msg)
        raise AssertionError(msg)

log = logging.getLogger(__name__)

config_file_def = "config/pgp_makeflow.ini"
ini_section = "main"

makeflow_rule_tpl = """\
%(targets)s: %(inputs)s
    %(cmd)s
"""

prep_script_tpl = """\
#!/bin/bash
source %(env)s
set -e
rm -f %(prep_flag_ok)s
prep_script=$PGP_HOME/ClusterSub.py
$PGP_PYTHON $prep_script -t %(inp_dir)s -s %(run_dir)s
#copy gbk files into run_dir too
cp %(gbk_src)s %(gbk_dest)s
touch %(prep_flag_ok)s
"""

inspect_script_tpl = """\
#!/bin/bash
source %(env)s
set -e
#cleaning up inspect_res at the start of the rule, otherwise bzip2 fails 
#inside inspect_wrapper
[ -n "$1" ] || exit 1
[ -n "$2" ] || exit 1
#rsync ../1739749481324218539/test.vics/ResultsX/$(basename $1) ResultsX/
#rsync ../1739749481324218539/test.vics/ResultsX/pepnovo/$(basename $3) ResultsX/pepnovo/
rm -f $1
SGE_TASK_ID=$2 $PGP_HOME/sgeInspect.sh 
"""

pvalue_script_tpl = """\
#!/bin/bash
source %(env)s
set -e
pvalue_script=$PGP_HOME/PValue.py
pvalue_args="-w . -p 0.1 -S 0.5 -1 -H"
#rsync -a ../1739749481324218539/test.vics/ResultsX/pvalue10H/ ResultsX/pvalue10H/
cd %(cwd)s
$PGP_PYTHON $pvalue_script -r %(pepnovo_dir)s $pvalue_args
"""

msgf_script_tpl = """\
#!/bin/bash
source %(env)s
set -e
#rsync ../1739749481324218539/test.vics/ResultsX/msgfOfPepnovo/$(basename $3) ResultsX/msgfOfPepnovo/
%(msgf_java)s %(msgf_java_args)s -jar %(msgf_jar)s -inspect $1 -d $2 > $3
"""

postproc_script_tpl = """\
#!/bin/bash
source %(env)s
set -e
rm -f %(postproc_flag_ok)s
postproc_script=$PGP_HOME/ProteogenomicsPostProcessing.py
$PGP_PYTHON $postproc_script \
-r %(msgf_dir)s \
-o %(_6frame_trie)s \
-b %(gbk)s \
-G %(peptides_gff)s \
-w %(other_out)s \
-p 1e-10
touch %(postproc_flag_ok)s
"""

class pgp_makeflow(object):
    """Generator of Makeflow input for the PGP pipeline."""

    def __init__(self,pgp_home,config,inp_dir):
        self.pgp_home = pgp_home
        self.config = config
        self.inp_dir = inp_dir

    def gen_all(self):
        config = self.config
        
        # create "prepare data" makefile
        make_prep = config.get(ini_section,"make_prep")
        self.mf_out = open(make_prep,"w")
        tasks_prep = self.gen_prep()
        self.mf_out.close()

        # execute "prepare data" makefile which is needed
        # to generate "process data" makefile
        run([config.get(ini_section,"makeflow"),make_prep])

        # create "process data" makefile
        start_dir = os.getcwd()
        try:
            os.chdir(tasks_prep[0]["proj_dir"])
            self.mf_out = open(config.get(ini_section,"make_proc"),"w")
            tasks_inspect = self.gen_inspect()
            tasks_pvalue = self.gen_pvalue(tasks_inspect=tasks_inspect)
            tasks_msgf = self.gen_msgf(tasks_pvalue=tasks_pvalue)
            tasks_postproc = self.gen_postproc(tasks_msgf=tasks_msgf)
            self.mf_out.close()
        finally:
            os.chdir(start_dir)

    def gen_prep(self):
        config = self.config
        mf_out = self.mf_out
        task_env = config.get(ini_section,"task_env")
        prep_flag_ok = config.get(ini_section,"prep_flag_ok")
        inp_dir = self.inp_dir
        run_dir = config.get(ini_section,"run_dir")
        work_dir = os.getcwd()
        db_genomic_dir = config.get(ini_section,"db_genomic_dir")
        gbk_ext = config.get(ini_section,"gbk_ext")
        gbk_src = pjoin(inp_dir,db_genomic_dir,"*"+gbk_ext)
        gbk_dest = pjoin(work_dir,run_dir,db_genomic_dir)
        script = prep_script_tpl % dict(
                env=task_env,
                inp_dir=inp_dir,
                run_dir=pjoin(work_dir,run_dir),
                prep_flag_ok=prep_flag_ok,
                gbk_src=gbk_src,
                gbk_dest=gbk_dest
                )
        script_file = config.get(ini_section,"prep_wrapper")
        strToFile(script,script_file)
        cmd = "bash %s" % (script_file,)
        tasks = [ dict(proj_dir=run_dir,
            prep_res=prep_flag_ok) ]
        prep_res_all = " ".join([ task["prep_res"] \
                    for task in tasks ])
        mf_out.write(makeflow_rule_tpl %\
                dict(targets=prep_res_all,
                    inputs=inp_dir,
                    cmd=cmd)
                )
        return tasks

    def gen_inspect(self):
        config = self.config
        mf_out = self.mf_out
        task_env = config.get(ini_section,"task_env")
        tasks_dir = config.get(ini_section,"inspect_tasks_dir")
        assert_log(os.path.isdir(tasks_dir),"Inspect tasks directory is missing",log)
        tasks_files = glob.glob(pjoin(tasks_dir,"*.in"))
        assert_log(len(tasks_files),"No job files found in inspect jobs directory",log)
        script = inspect_script_tpl % dict(
                env=task_env,
                )
        script_file = config.get(ini_section,"inspect_wrapper")
        strToFile(script,script_file)

        tasks = [ dict(inspect_task_id=stripSfx(os.path.basename(j)),
            inspect_task_file=j) \
                    for j in tasks_files ]
        
        for task in tasks:
            task["inspect_res"] = pjoin(config.get(ini_section,"inspect_dir"),
                    task["inspect_task_id"]+".txt.bz2")
            task["pepnovo_res"] = pjoin(config.get(ini_section,"pepnovo_dir"),
                    task["inspect_task_id"]+".res")
        
        for task in tasks:
            cmd = "bash %s %s %s %s" % (script_file,
                    task["inspect_res"],
                    task["inspect_task_id"],
                    task["pepnovo_res"]
                    )
            mf_out.write(makeflow_rule_tpl %\
                    dict(targets="%s %s" % (task["inspect_res"],task["pepnovo_res"]),
                        inputs=task["inspect_task_file"],
                        cmd=cmd)
                    )
            
        return tasks

    def gen_pvalue(self,tasks_inspect):
        config = self.config
        mf_out = self.mf_out
        task_env = config.get(ini_section,"task_env")
        pepnovo_dir = config.get(ini_section,"pepnovo_dir")
        pvalue_dir = config.get(ini_section,"pvalue_dir")
        makedir(pvalue_dir)
        script = pvalue_script_tpl % dict(
                env=task_env,
                cwd=pvalue_dir,
                pepnovo_dir=os.path.abspath(pepnovo_dir)
                )
        script_file = config.get(ini_section,"pvalue_wrapper")
        strToFile(script,script_file)
        cmd = "bash %s" % (script_file,)
        pepnovo_res_all = " ".join([ task["pepnovo_res"] for task in tasks_inspect ])
        tasks = [ dict(pvalue_task_id=task["inspect_task_id"],
            pvalue_res=pjoin(pvalue_dir,task["inspect_task_id"]+".res")) \
                    for task in tasks_inspect ]
        pvalue_res_all = " ".join([ task["pvalue_res"] \
                    for task in tasks ])
        mf_out.write(makeflow_rule_tpl %\
                dict(targets=pvalue_res_all,
                    inputs=pepnovo_res_all,
                    cmd=cmd)
                )
        return tasks

    def gen_msgf(self,tasks_pvalue):
        config = self.config
        mf_out = self.mf_out
        task_env = config.get(ini_section,"task_env")
        pvalue_dir = config.get(ini_section,"pvalue_dir")
        mzxml_dir = config.get(ini_section,"mzxml_dir")
        msgf_dir = config.get(ini_section,"msgf_dir")
        #somehow on TACC Ranger JAVA_HOME defined by 'module add'
        #does not transfer to compute nodes as it should with 
        #qsub -V
        makedir(msgf_dir)
        script = msgf_script_tpl % dict(
                env=task_env,
                msgf_java=config.get(ini_section,"msgf_java"),
                msgf_java_args=config.get(ini_section,"msgf_java_args"),
                msgf_jar=config.get(ini_section,"msgf_jar"),
                )
        script_file = config.get(ini_section,"msgf_wrapper")
        strToFile(script,script_file)
        tasks = []
        for task_pvalue in tasks_pvalue:
            task = dict(msgf_task_id=task_pvalue["pvalue_task_id"],
                    msgf_res=pjoin(msgf_dir,task_pvalue["pvalue_task_id"]+".msgf"))
            cmd = "bash %s %s %s %s" % (script_file,
                task_pvalue["pvalue_res"],
                mzxml_dir,
                task["msgf_res"])
            mf_out.write(makeflow_rule_tpl %\
                    dict(targets=task["msgf_res"],
                        inputs=task_pvalue["pvalue_res"],
                        cmd=cmd)
                    )
            tasks.append(task)
        return tasks

    def gen_postproc(self,tasks_msgf):
        config = self.config
        mf_out = self.mf_out
        task_env = config.get(ini_section,"task_env")
        msgf_dir = config.get(ini_section,"msgf_dir")
        postproc_flag_ok_root = config.get(ini_section,"postproc_flag_ok")
        db_genomic_dir = config.get(ini_section,"db_genomic_dir")
        gbk_ext = config.get(ini_section,"gbk_ext")
        gbk_glob = pjoin(db_genomic_dir,"*"+gbk_ext)
        postproc_results_dir = config.get(ini_section,"postproc_results_dir")
        makedir(postproc_results_dir)
        gff_dir = config.get(ini_section,"gff_dir")
        makedir(gff_dir)
        """
        -r /path/to/archive/PSM.pipeline \
        -o /path/to/archive/Databases/Genomic/NC_000911.6frame.trie \
        -b /path/to/archive/Databases/Genomic/NC_000911.gbk \
        -G /path/to/archive/DerivedData/GFFs/NC_000911.peptides.gff \
        -w /path/to/archive/DerivedData/NC_000911.outputstub.txt \
        -p 1e-10
        """
        msgf_res_all = " ".join([ task["msgf_res"] for task in tasks_msgf ])
        tasks = []
        gbk_files = glob.glob(gbk_glob)
        for gbk in gbk_files:
            acc = os.path.basename(gbk).rsplit(gbk_ext,1)[0].strip()
            _6frame_trie = pjoin(db_genomic_dir,
                    acc+config.get(ini_section,"_6frame_trie_ext"))
            peptides_gff = pjoin(gff_dir,
                    acc+config.get(ini_section,"peptides_gff_ext"))
            other_out = pjoin(postproc_results_dir,
                    acc+config.get(ini_section,"other_out_ext"))
            postproc_flag_ok = postproc_flag_ok_root+"."+acc
            script = postproc_script_tpl % dict(
                    env=task_env,
                    msgf_dir=msgf_dir,
                    gbk=gbk,
                    _6frame_trie=_6frame_trie,
                    peptides_gff=peptides_gff,
                    other_out=other_out,
                    postproc_flag_ok=postproc_flag_ok
                    )
            script_file = config.get(ini_section,"postproc_wrapper")+"."+acc
            strToFile(script,script_file)
            task = dict(postproc_task_id=acc,
                    postproc_res=postproc_flag_ok)
            cmd = "bash %s" % (script_file,)
            mf_out.write(makeflow_rule_tpl %\
                    dict(targets=task["postproc_res"],
                        inputs=msgf_res_all,
                        cmd=cmd)
                    )
            tasks.append(task)
        all_postproc_flag_ok = " ".join([ task["postproc_res"] for task in tasks ])
        mf_out.write(makeflow_rule_tpl %\
                dict(targets=postproc_flag_ok_root,
                    inputs=all_postproc_flag_ok,
                    cmd="LOCAL echo ok > {}".format(postproc_flag_ok_root))
                )
        return [dict(postproc_res=postproc_flag_ok_root)]

def getProgOptions():
    from optparse import OptionParser, make_option
    option_list = [
        make_option(None, "--config-file",
        action="store", 
        type="string",
        help="Config file [pgp_home/%s]" % (config_file_def,),
        dest="config_file"),
        make_option(None, "--inp-dir",
        action="store", 
        type="string",
        help="Directory with input data, accessible to submit node",
        dest="inp_dir"),
        make_option(None, "--work-dir",
        action="store", 
        type="string",
        help="Working directory, accessible to both submit and compute nodes",
        dest="work_dir"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args()

    return options,args

if __name__ == "__main__":
    pgp_home = os.path.dirname(sys.argv[0])
    opt,args = getProgOptions()
    assert opt.inp_dir is not None,"--inp-dir is mandatory argument"
    opt.inp_dir = os.path.abspath(opt.inp_dir)
    if opt.config_file:
        config_file = opt.config_file
    else:
        config_file = pjoin(pgp_home,config_file_def)
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    pgp_home = config.get(ini_section,"pgp_home")
    if opt.work_dir is None:
        opt.work_dir=os.getcwd()
    else:
        makedir(opt.work_dir)
    start_dir = os.getcwd()
    try:
        os.chdir(opt.work_dir)
        mkf = pgp_makeflow(config=config,pgp_home=pgp_home,inp_dir=opt.inp_dir)
        mkf.gen_all()
    finally:
        os.chdir(start_dir)

