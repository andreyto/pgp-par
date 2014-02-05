PGP: Parallel Proteogenomics Pipeline for MPI clusters, high-throughput batch clusters and multicore workstations
================================================================================================================

The open source, automated proteogenomics pipeline (PGP) described here
is freely available under [GPLv3](http://www.gnu.org/copyleft/gpl.html)
license. It is designed to run in a range of different parallel Linux
computing environments:

-   HPC (high-performance computing) clusters that are set up to
    efficiently schedule only large (100s+ of cores) parallel MPI jobs
    under a control of batch queuing system such as (formerly) Sun Grid
    Engine (SGE), SLURM or PBS/Torque. Our primary targets for this use
    case were compute clusters of XSEDE (<https://www.xsede.org/>), the
    federation of supercomputers supported by the US National Science
    Foundation. XSEDE allocates its resources to outside researchers
    through a peer reviewed proposal system. The biologists can now use
    our software on this major computational resource.
-   HTC (high-throughput computing) clusters widely used as local
    bioinformatics computing resources. These clusters are configured to
    efficiently schedule large numbers of serial jobs under a control of
    batch queuing system.
-   A single multi-core workstation without a batch queuing system
    (including a case of single core machine).

Parallel execution ability is important for proteogenomic annotation
software due to a high volume of required computations (order of 100
CPU\*hrs for a typical bacterial genome).

### Citations
Please cite this article if you use this software
[[PMID:24470574](http://bioinformatics.oxfordjournals.org/content/early/2014/01/27/bioinformatics.btu051)].

### The significance of proteogenomic annotation

Our pipeline is a tool for improving the existing genomic annotations
from available proteomics mass spectrometry data. As most genome
annotation pipelines consist of automated gene finding, they lack
experimental validation of primary structure
[[PMC2265698](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2265698/),
[PMC2238897](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2238897/)],
having to rely on DNA centric sources of data such as sequence homology,
transcriptome mapping, codon frequency, etc. By incorporating the
orthogonal set of data, proteogenomics is able to discover novel genes,
post-translational modifications (PTMs) and correct the erroneous
primary sequence annotations.

Quick Start
-----------

To test the pipeline, you can run it on a single workstation (preferably
multicore one).

Make sure that your are in a directory where you have write access and
enough free space (2G to process the included small test dataset).

Follow the Installation instructions from this manual, downloading the
source code, editing the environment scripts and CMake configuration
variables if necessary, and then running the installation command:

    SOURCE/install --target-env htc --prefix ../PGP

where `SOURCE` stands for a path to a directory where you have
downloaded or checked out the source code. The command above will
install the pipeline into a directory `PGP` in your current working
directory.

The installation procedure will automatically test the pipeline in a
local (multicore) execution mode on (severely) subsampled input data
included in the package distribution as
`SOURCE/proteogenomics/TestSuite/Bacillus.anthracis.sterne.PNNL.chunk.4sp.tar.gz`
and report the test outcome. This small test will take about 10 min to
run on a quad core workstation. To reduce the run time, the test dataset
contains only four MS runs and only the first 460 Kbp of the reference
chromosome.

The test invokes `PGP/bin/pgp_run` script that is the top level
executable of the PGP pipeline.

The test will also print the pipeline command line that was executed,
which will include the locations of the input and output directories.
You can copy that command line, modify the options and try running it
again from the shell prompt. Pass `--help` argument to `pgp_run` for a
description of available options.

From under the output directory, you can load the GFF3 files from
`DerivedData/GFFs` and the corresponding reference from
`Databases/Genomic` into any of the genomic browsers such as [NCBI
Genome workbench](http://www.ncbi.nlm.nih.gov/tools/gbench/), Artemis
[[PMC3278759](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3278759/)] or
[CLC Genomics
workbench](http://www.clcbio.com/products/clc-genomics-workbench/) for a
visual inspection of the newly created annotations.

See the description of `PGP_TEST_PGP_RUN_EXTRA_OPTIONS` CMake variable
below if you want the automated test to be executed with a distributed
cluster backend instead of the default local multi-core backend.

I). The PGP Algorithm
---------------------

The full protocol of our pipeline is described in detail in
[[PMC3219674](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3219674/)]. In
that study, we have applied the pipeline to 46 genomes spanning eight
bacterial and archaeal phyla across the tree of life. These diverse
datasets facilitated the development of a robust approach for
proteogenomics that is functional across genomes varying in %GC, gene
content, proteomic sampling depth, phylogeny, and genome size. In
addition to finding evidence for 682 novel proteins, 1336 new start
sites, and numerous likely mis-annotations, we discovered sites of
post-translational maturation in the form of proteolytic cleavage of
1175 signal peptides. The output files from this study are available at
(<http://omics.pnl.gov/pgp/overview.php>). These results have been
submitted back to the NCBI, where they are being used for generating
improved updated RefSeq files. The NCBI RefSeq updates can be seen in
the Genbank flat files (.gbk) of the corresponding genomes wherever the
proteomics data are listed as an experimental evidence. One example is
the [Mycobacterium tuberculosis H37Rv
genome](ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Mycobacterium_tuberculosis_H37Rv_uid57777/NC_000962.gbk)
containing the CDS attributes /experiment="EXISTENCE: identified in
proteomics study".

The pipeline combines several open source proteomics tools from
(<http://proteomics.ucsd.edu>) with our own post-processing steps.

Briefly, tandem mass spectra are searched by Inspect
[[PMID:16013882](http://www.ncbi.nlm.nih.gov/pubmed/16013882)] against a
translation of the genome and subsequently rescored with PepNovo
[[PMID:15858974](http://www.ncbi.nlm.nih.gov/pubmed/15858974)] and MSGF
[[PMC2689316](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2689316/)].
The pipeline translates the DNA sequence in all six frames to generate a
reference protein database.

Each stop to stop open reading frame (ORF) is included regardless of
coding potential. The concatenated decoy records for estimating the
statistical significance of the predicted peptides are generated by
shuffling each ORF. Significant peptide/spectrum matches (PSM) are those
with a ***p-value*** of 1.e-10 or better. All confident peptides are
mapped onto their genomic location (nucleotide coordinates) and grouped
into sets within an ORF.

The last stage of the pipeline employs five ORF filters. First, we
remove low complexity peptides and peptides which are more than 750 bp
from the next in-frame peptide. We remove ORFs which lack a uniquely
mapping peptide or which lack a fully tryptic peptide. Finally, we
require two peptides per protein.

Peptides from ORFs which meet all of the above criteria are saved in the
output GFF3 files in coordinates relative to the reference genome. The
user can load these files into genomic browsers such as NCBI Genome
workbench or Artemis for subsequent analysis.

### Example of interpreting the output

Peptides mapping to regions of the genome which lack a protein
annotation represent either novel genes, or 5' extensions of current
genes. For example, at 3.376 MB in the B. anthracis Sterne genome lies
the BAS3403 gene, the small spore protein Tlp. Many peptides from this
protein were discovered by the pipeline. In the open reading frame
directly upstream, the protocol detected other peptides which were not
part of any currently annotated protein. Blasting this ORF revealed
homology to another spore coat protein from B. cereus, not annotated in
any anthrax genome.

![Image of Novel Anthrax Sporecoat
Protein](https://bitbucket.org/andreyto/proteogenomics/raw/master/proteogenomics/docs.PGP/anthrax_sporecoat.png)

II). Parallelization strategy
-----------------------------

Our goal was to build a pipeline that would be portable across different
parallel execution environments that users might encounter. The overall
algorithm is embarrassingly parallel for the most part. Specifically, it
is possible to process each spectrum file independently throughout all
computationally intensive stages of the algorithm. There is a global
synchronization point in a middle of the workflow when the algorithm
needs to build a histogram of all scores for p-value computation. We
have selected a distributed workflow model where multiple serial
processes are executed concurrently following a dependency graph defined
by required input and output files. This model exploits inherent
parallelism of the problem sufficiently well while being compatible with
a wide variety of execution environments such as standalone multicore
machines, high-throughput compute clusters and MPI clusters (more on the
latter below).

We achieve the portability across different execution environments by
using the Makeflow workflow execution engine
(<http://nd.edu/~ccl/software/makeflow/>). Makeflow also provides a high
degree of fault tolerance against compute node failures (in non-MPI
execution mode) and restart capability in case of master node failure.

Our installation procedure builds its own local copy of the CCTools
package that contains Makeflow and associated backend executables.

Given a directory with input data, our pipeline's code generates a
description of the workflow (in a language that is based on the `make`
file syntax). Makeflow can then execute this workflow using different
parallel "backends" selected by the user at run-time. While its backend
implementations are straightforward on high-throughput clusters (they
submit each task to the batch queuing system) and multicore workstations
(they run multiple subprocesses), the architecture is more complicated
on HPC clusters that often have scheduling policies tuned to efficiently
allocate only large MPI jobs. On such MPI clusters, Makeflow supports
execution of complex workflows composed of many interdependent serial
jobs through a "glide-in" mechanism. Specifically, the user submits the
MPI backend executable of Makeflow as a single parallel MPI job (e.g.
with 100 MPI process ranks). Then, the user starts the Makeflow master
on a single node. The individual MPI process ranks connect to the
master, and the master farms out single serial tasks to the ranks
following the order of workflow dependencies. A given backend rank
executes each incoming task as a separate subprocess.

For both HTC and MPI clusters, the pipeline assumes that there is a
shared file system mounted in the same location on all compute nodes as
well as on the job submission ("login") node.

III). Installation
------------------

### A). Dependencies that have to be present on your system to build and run the pipeline

The integrated installation procedure needs:

-   [Git](http://git-scm.com/) version control system to checkout the
    source from the public repository on BitBucket
    (<https://bitbucket.org/andreyto/proteogenomics>). If you instead
    download the package archive from BitBucket, you will not need Git.
-   [CMake](http://www.cmake.org/) configuration and build utility
    (version 2.8 or higher)
-   *BASH* shell
-   C++ compiler (gcc or MPI wrappers, depending on the targeted
    execution environment)
-   Python interpreter (\>=2.7) and Python development libraries
-   GNU Make
-   Internet connection. The installation procedure might try to
    download some Python dependencies from Python Package Index (PyPi).

All software components listed above are available as standard packages
in major Linux distributions.

For run-time, you will need:

-   Python (\>=2.7) and Java (\>=1.6) runtimes
-   MPI environment if MPI backend is used for execution
-   Some batch queuing system supported by Makeflow (e.g. SGE or any of
    its clones, PBS, Condor) if a distributed backend will be used for
    execution

### B). Getting the sources

    git clone git@bitbucket.org:andreyto/proteogenomics.git

this will put the sources into a subdirectory `proteogenomics` under your 
current working directory.

**or** download an archive from the BitBucket repository
(<https://bitbucket.org/andreyto/proteogenomics>).
On BitBucket, you should choose the download archive under `Downloads` 
menu -> tab `Branches` -> Branch `master` -> link to one of the 
available archive formats. For example, this is the link to the 
["gz" archive](https://bitbucket.org/andreyto/proteogenomics/get/master.tar.gz).
Expand the downloaded "gz" archive with `tar -zxf master.tar.gz` 
and rename the resulting folder to`proteogenomics`.

### C). Building and installing

A shell script `install` is a wrapper that activates the necessary shell
environment and then calls CMake's configuration, build, install and
test stages in a single run. Our CMake build procedure attempts to
automatically install the pipeline itself as well as some of its
dependencies included in our repository. When you are executing this
script, your current working directory (CWD) must be **outside** of and
**not** directly above the `proteogenomics` source directory. For
example, in a directory where you have issued the `git clone` command
described above, you can do:

    mkdir build; cd build; ../proteogenomics/install --prefix PGP_ROOT --target-env TARGET_ENV

where the placeholders `PGP_ROOT` and `TARGET_ENV` must be replaced with
the actual values as described below.

-   `--target-env TARGET_ENV` Sets the computing environment name
    -   `ranger` is the environment configured for XSEDE TACC Ranger
        cluster (that was our primary testing site for the MPI execution
        mode).
    -   `htc` covers a generic case of HTC clusters and multi-core workstations.
    -   `htc_bu_scc` is for BU SCC cluster used in HTC mode.

    The computing environments are discussed in detail below under
    "Customizing the build procedure". `htc` is assumed by default.

-   `--prefix PGP_ROOT` Sets the installation directory. Everything will
    be installed under that directory, which will be created if it does
    not already exist. Wrapper scripts for a typical pipeline invocation
    will be under `PGP_ROOT/bin` and configurations files under
    `PGP_ROOT/config`.

#### Customizing the build procedure

Although the CMake that is invoked by `install` will try to figure out
the locations of necessary libraries and executables on your system
automatically, it might still need some explicit instructions that you
can provide by editing specific files under the checked out **source**
directory (designated by `SOURCE` placeholder further in the text)
**before** you invoke `install`. The files are located under
`SOURCE/config` and `SOURCE/bin`. You are more likely to need editing
them on compute clusters, where the administrators often build multiple
customized versions of various packages and place them under
non-standard locations (such as `/usr/local/packages`). There are also
local differences in how the MPI programs are invoked if MPI execution
backend is used.

Under `SOURCE/config` and `SOURCE/bin` directories, we have
subdirectories named after each of the target execution environments for
which we have tuned and tested building and execution of the pipeline
(designated by `TARGET_ENV` placeholder here). Currently, those are
`ranger` and `htc`. There is also a subdirectory `noarch` for files that
are common to all environments. You typically will not have to edit
anything under `noarch`. To tune the build and run time for your
specific environment, you might have to edit two or three files under
the environment-specific subdirectory that is the closest to your target
system, and then pass the corresponding environment name to the
`install` such as `--target-env ranger`. Those files are:

-   `SOURCE/config/TARGET_ENV/pgp_login.rc` This file will be "sourced"
    by install and job submission scripts at least on the login node
    (the machine where you both run the `install` and submit the
    pipeline for execution). You should put where any Bash shell
    commands and variables needed to properly configure your environment
    (such as setting `LD_LIBRARY_PATH` and `PYTHONPATH` variables to
    make sure that your GCC compiler and Python interpreter work if they
    are installed in non-standard locations).

    For the `ranger` environment, the `pgp_login.rc` uses "module"
    commands to make available the proper dependencies such as MPI
    compilers for building the Makeflow. "Module" command is a
    standardized user environment management script used by XSEDE
    clusters (<https://www.xsede.org/software-environments>).
    "module" commands are also used for `htc_bu_scc` environment.

    On XSEDE systems, the users will have to modify the details of
    specific package versions activated by the "module" command in order
    to adapt the building and execution environment to their specific
    cluster.

-   `SOURCE/config/TARGET_ENV/toolchain.cmake` We customize CMake
    variables through this file. The variable that you most likely want
    to modify where is called `PGP_TEST_PGP_RUN_EXTRA_OPTIONS`. You can
    edit it to specify any extra options to `PGP_ROOT/bin/pgp_run` to be
    used when the automated test is executed at the end of the
    installation procedure. In particular, you can tell the pipeline to
    run the test using distributed cluster backend, and describe the
    necessary batch system parameters. The existing `toolchain.cmake`
    files contain commented out examples of setting this variable. There
    is a good chance that you would not need to modify anything else in
    this file. If you do, you can look at the toolchain file under the
    `ranger` or `htc` subdirectory as an example, and consult the CMake
    [documentation](http://cmake.org/cmake/help/v2.8.8/cmake.html).

-   `SOURCE/bin/TARGET_ENV/pgp_wrapper_nested_mpi.in` When the user
    selects MPI backend for execution, this wrapper script is submitted
    as a parallel batch job (a batch job spanning multiple nodes). This
    script is responsible for starting Makeflow master with a given make
    file and starting MPI worker daemons. The existing file under
    `ranger` subdirectory contains detailed comments explaining its
    inner workings. You will likely only need to modify the MPI command
    line because that might differ from site to site even if the same
    MPI library is used. This script is not installed if
    `install --target-env htc` is used.

You can also create copies of either `ranger` or `htc` subdirectories
under `SOURCE/config` and `SOURCE/bin`, using a name that is appropriate
for your target environment (e.g.
`cp -a SOURCE/bin/ranger SOURCE/bin/my_cluster;` and
`cp -a SOURCE/config/ranger SOURCE/config/my_cluster`); edit the files
in the new subdirectories as you see fit; and use this name as
`install --target-env my_cluster`.

During installation, files from the `TARGET_ENV` subdirectories are
processed by CMake with variable substitution for any template
parameters escaped by the `@` symbol on both sides, and then copied into
`PGP_ROOT/config` and `PGP_ROOT/bin` directories, along with any files
from `noarch` subdirectories.

The generic `htc` configuration has been tested on Ubuntu 12.04
workstation with `local` execution backend and on CentOS 6 JCVI HTC 
cluster with both `local` and `sge` backends.

IV). Running the pipeline
-------------------------

At the high level, the pipeline execution consists of two stages. The
first, serial stage prepares the input data and generates the workflow
file for the second stage. The second stage processes the input data in
multiple parallel tasks (with some internal barrier synchronization
steps, for example, for computing the *p-values*).

Because scheduling policies and node availability are typically vastly
different between serial and parallel jobs on MPI clusters like XSEDE
systems, we have given to the user a fine grained control over running
these two stages ("prepare" and "process") in such environments as
described below. Even for the HTC cluster environments, it might be
useful to force the "prepare" stage to execute locally on the submit
node if the file system where the input data resides is not mounted on
the compute nodes.

The users can always modify the parameters of Makeflow execution and
backend job submission by following the Makeflow User Manual
(<http://www3.nd.edu/~ccl/software/manuals/makeflow.html>), possibly
using the parameters of our testing runs as the starting points (see
also description of `PGP_TEST_PGP_RUN_EXTRA_OPTIONS` CMake variable in
this document).

Some site-specific tuning of job submission options will likely to be
required on any computational cluster, considering the multitude of
customizations of both the operating and batch systems that cluster
administrators typically employ.

The following steps must be taken to execute the proteogenomics
analysis:

A). Data Requirements: Genomic & Proteomic Data
-----------------------------------------------

1.  Annotated nucleotide (mandatory) and protein (optional) sequence
    files (e.g. **gbk**, **fna** and **faa** RefSeq files from NCBI
    <http://www.ncbi.nlm.nih.gov/refseq/>).
2.  Mass spectrometric data (spectra) in mzxml format (mzXML conversion
    tools can be found at
    <http://tools.proteomecenter.org/wiki/index.php?title=Formats:mzXML>)

B). Structured Input Data Directory
-----------------------------------

1.  The input data mentioned under **'A)'** should be organized in a
    specific directory structure as shown by the listings of the
    directories below. This can be easily accomplished by uncompressing
    the `template.directory.tar.gz` provided under
    `PGP_ROOT/proteomics/data`. Assuming the annotated genomic data
    comes from speciesX, uncompress the template directory inside
    speciesX folder:

        user@assembly [\~/speciesX]% ls -l

        total 32

        drwxr-x--- 5 user tigr 4096 Aug 30 12:01 Databases

        drwxr-x--- 3 user tigr 4096 Aug 30 12:01 DerivedData

        drwxr-x--- 2 user tigr 4096 Aug 30 12:01 OrthologyClusters

        drwxr-x--- 3 user tigr 4096 Aug 30 12:01 PeptidomeSubmission

        drwxr-x--- 2 user tigr 4096 Aug 30 12:01 Publications

        drwxr-x--- 2 user tigr 4096 Aug 30 12:01 RawZips

        drwxr-x--- 4 user tigr 4096 Aug 30 12:01 RefSeqSubmission

        drwxr-x--- 2 user tigr 4096 Aug 30 12:01 mzxml

2.  Put your RefSeq fna and gbk files for speciesX into
    `Databases/Genomic` folder and faa file(s) - into
    `Databases/Predictions` folder. Any additional protein databases
    including common contaminants such as trypsin and keratin (but not
    the DNA), which also need to be searched with Inspect, go into the
    `Databases/Proteomic` folder. The template directory already has
    some unpacked Common.RS files in this directory (as shown below).

        user@assembly [\~/tmp/db/Databases]: [12:07:22]% ls -l

        total 12

        drwxr-x--- 2 user tigr 4096 Feb 19 2010 Genomic

        drwxr-x--- 2 user tigr 4096 Feb 19 2010 Predictions

        drwxr-x--- 2 user tigr 4096 Feb 19 2010 Proteomic

        ./Genomic:

        total 0

        ./Predictions:

        total 0

        ./Proteomic:

        total 16

        -rw-r----- 1 user tigr 1472 Feb 19 2010 Common.RS.index

        -rw-r----- 1 user tigr 8254 Feb 19 2010 Common.RS.trie

3.  Your spectra mzxml file(s) should be compressed as `tar.gz` format
    (e.g. by `tar -czf spectra.tar.gz *.mzXML` and placed in the `mzxml`
    directory.

4.  The remaining empty directories (some optional) serve as
    placeholders to be populated by the data generated during the
    analysis. Some of them are the suggested locations for the last,
    manual curation stages of the proteogenomic re-annotation process
    aimed to culminate in the submission of results to NCBI RefSeq as
    described in the pipeline protocol publication
    [[PMC3219674](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3219674/)].
    If the user is only interested in the results of the automatic
    annotation but not in the NCBI submission, those directories can be
    ignored.

-   **DerivedData** - copied from output directory (see below). Specific
    files or results that should be produced with every dataset
-   **GFFs** - a sub-directory copied from output directory (see below)
    -   contains the gff file with peptide annotations for the genome

-   **OrthologyClusters** - folder with MSA (multiple sequence
    alignments) of orthologous clusters from related organisms
-   **Publications** - if there are relevant publications related to the
    data, they go here.
-   **RawZips** - holds the raw Inspect and PepNovo rescored data in bz2
    files
-   **PeptidomeSubmission** â€“ the Peptidome submission needs meta data
    file, PSM results file, spectra & the results file - as explained as
    below.
-   **Peptidome** data generation: meta data generated by the user
-   **PepXML** - these are pepXML conversions of Inspect's output (PSMs)
-   **mzxml** - spectra files (mzXML.tar.gz)
-   **ProteinInference** - the protein inference as made by the PGP
    script
-   **RefseqSubmission** - originals submitted - contains new
    annotations generated by the PGP

C). Execution on HTC clusters or multicore workstations
-------------------------------------------------------

Execute:

    PGP_ROOT/bin/pgp_run [options] INPUT_DIR OUTPUT_DIR

where `OUTPUT_DIR` should be on a shared file system if you are running
on a cluster. The arguments in square brackets `[...]` are optional.
They can be given before or after positional arguments.

Run `PGP_ROOT/bin/pgp_run --help` for a description of all options.

The options marked as `Makeflow:` in the help output will be passed to
the `makeflow` program that executes the generated workflow graph. If
the `--local-prep` flag is provided, the data preparation stage will be
executed locally on the submit node. You might have to use that flag if
your compute nodes do not have access to the file system where the
`INPUT_DIR` resides (e.g. it is on some kind of archival storage only
available on the login node of your cluster). `OUTPUT_DIR` should be on
a shared file system accessible from cluster compute nodes.

On a single workstation without a batch system, the Makeflow arguments
can be omitted. In that case, Makeflow will run the number of concurrent
subprocesses that equals the number of cores on the machine.

On SGE cluster, the command might look like

    PGP_ROOT/bin/pgp_run INPUT_DIR OUTPUT_DIR -T sge -B 'SGE options'

where the string `SGE options` (must be in quotes as shown above) is
passed verbatim to the SGE when the Makeflow submits each task in the
workflow. This string must contain the options that you would pass to a
`qsub` in order to run a single job on your system. For example, on JCVI
internal cluster the full command would look like

    PGP_ROOT/bin/pgp_run INPUT_DIR OUTPUT_DIR -T sge -B '-P 0000 -b n -S /bin/bash'

You can look at the Makeflow manual for possible command line arguments
if you want to use other backends such as PBS or Condor, or modify the
behavior of Makeflow (e.g. to constrain the number of concurrent jobs).

You can execute `PGP_ROOT/vendor/bin/makeflow --help` to have Makeflow
print a short description of options that it can take. In particular,
the description of the `-T` option contains the list of batch systems
that Makeflow supports directly. Yet other types of batch systems can be
used through work-queue backend by submitting their workers through
mechanisms native to each batch system (see Makeflow manual).

In the commands above, the controlling process of the Makeflow (the
"master") and, optionally, the data preparation stage will run on the
login node, and the script `pgp_run` will be in a running state until
the entire workflow has finished. This script itself will use very
little resources (outside of the local `prepare` stage), with Makeflow 
master merely farming  out tasks to the batch system. 

In case you do not want to wait for`pgp_run`
on the login node, you can also submit the command above itself to your
batch system through a standard (e.g. qsub) mechanism. This would
require that your compute nodes are allowed to submit new jobs
themselves, because Makeflow will be submitting new jobs from the node
where the master script is executing.

`pgp_run` generates a shell script `this_pgp_run.sh` under `output_dir`
that will be used to run the generated workflow. Instead of having
`pgp_run` to execute that workflow script immediately by default, you
can use `--do-not-run` option to `pgp_run`, and execute
`./this_pgp_run.sh` later from shell after making `output_dir` your
current working directory. You can also submit this script to your batch
system (assuming your compute nodes are allowed to submit jobs too, and
they can access `INPUT_DIR`).

When the pipeline finishes successfully, it prints to standard output
"Nothing left to do". If it fails, it prints "Makeflow failed".

While it runs, you can use
`PGP_ROOT/vendor/bin/makeflow_monitor OUTPUT_DIR/pgp_proc.mkf.makeflowlog`
to monitor the progress of the `process` stage. Once it finishes, you
can use
`PGP_ROOT/vendor/bin/makeflow_log_parser OUTPUT_DIR/pgp_proc.mkf.makeflowlog`
to print a summary of CPU time used by the tasks.

D). Execution on XSEDE MPI clusters
-----------------------------------

The procedure is the same as described above for HTC clusters. The
difference is in providing a few extra options to `pgp_run`:

1.  `--pgp-batch-options-extra-prepare` should provide additional batch
    system options that would request a serial job for the `prepare`
    stage. These option string will be concatenated with the common
    batch option string provided by -B|--batch-options argument. This
    concatenated value will be used by the Makeflow when submitting the
    `prepare` job to the batch system (see the example below).

2.  `--pgp-batch-options-extra-process` should provide additional batch
    system options that would request a parallel MPI job for the
    `process` stage (see the example below).

3.  `--pgp-process-mpi` activates the use of MPI workflow backend during
    the `process` stage.

**Example:** On TACC Ranger, the command would look like this (note the
shell line continuation backlash symbol `\`):

    PGP_ROOT/bin/pgp_run --pgp-process-mpi -T sge \
    -B '-V -b n -S /bin/bash' \
    --pgp-batch-options-extra-prepare '-pe 1way 16 -q serial -l h_rt=01:30:00' \
    --pgp-batch-options-extra-process '-pe 16way 64 -q normal -l h_rt=12:00:00' \
    INPUT_DIR OUTPUT_DIR

The following will happen when that command line is processed:

-   Makeflow will submit jobs to SGE batch system
-   The `prepare` stage will be submitted to a queue called `serial`
    with `qsub` options
    `-V -b n -S /bin/bash -pe 1way 16 -q serial -l h_rt=01:30:00`. The
    `-V` option tells SGE to propagate environment variables from login
    (submit) node to the compute nodes. This switch **has** to be used
    on XSEDE clusters that have `module` system for configuring the user
    environment. The `-pe 1way 16` requests the minimum number of CPU
    cores possible on that specific cluster (a single node). The
    `-l h_rt=01:30:00` requests 1.5 hours of run-time limit.
-   The `process` stage will be submitted to a queue called `normal` as
    a single parallel MPI job running on 64 CPU cores (4 nodes) with a
    (wall clock) run time limit of 12 hours. The Makeflow will farm out
    serial tasks to its "glide-in" daemons within this MPI job.

The users should consult their cluster manual to find out how to specify
batch parameters for serial and parallel jobs. The number of CPU cores
in the parallel job should not exceed one third of the number of MS runs
in the input directory in order to ensure reasonable load balancing and
CPU utilization within the MPI process.

V). Results
-----------

**Structured Output Data Directory.** When the PGP analysis is complete,
in addition to some log files, the following directories are created in
`OUTPUT_DIR` directory.

-   **Databases:** the 6-frame translated refseq genomic data copied
    from input directory.
-   **mzxml:** MS spectral data copied from input directory.
-   **DerivedData:** contains the proteogenomics analysis results and
    various reports for each input dataset. **That directory contains
    the main output of the pipeline.**
-   **DerivedData/GFFs:** a subdirectory with the GFF3 files, one per
    chromosome, containing peptide annotations of the genome generated
    from mapping the proteomics spectra.

    See
    [[PMC3219674](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3219674/)]
    regarding the interpretation and downstream use of these
    annotations.

-   **ResultsX:** the directory contains results from Inspect analysis.
-   **pepnovo:** the directory contains results from PepNovo analysis.
-   **pvalue10H:** the sub-directory with the output of PValue.py (with
    -p 0.1, -H).
-   **msgfOfPepnovo:** the sub-directory with the output of MS\_GF
    validation.
-   **jobs, Done, output:** other directories with pipeline working data
    (could be empty).
