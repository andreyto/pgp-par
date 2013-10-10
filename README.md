Automated Proteogenomics Pipeline for MPI clusters, high-througput batch clusters and multicore workstations
========================================================================

The open source, automated proteogenomics pipeline (PGP) described here is freely accessible to the scientific community. It is designed to run in different kinds of parallel Linux computing environments:

- HPC (high-performance computing) clusters that are set up to efficiently schedule only large (100s+ of cores) parallel MPI jobs under a control of batch queuing system such as (formerly) Sun Grid Engine (SGE), SLURM or PBS/Torque. Our primary targets for this use case were compute clusters of XSEDE ([https://www.xsede.org/][1]), the federation of supercomputers supported by the US National Science Foundation. XSEDE allocates its resources to outside researchers through a peer reviewed proposal system. We considered it as an important requirement that the biologists could use our software on this major computational resource.

- HTC (high-throughput computing) clusters widely used as local bioinformaics computing resources. These clusters are configured to efficiently schedule large numbers of serial jobs under a control of batch queuing system.

- A single multi-core workstation without a batch queuing system (including a case of single core machine).

Parallel execution ability is important for proteogenomic annotation software due to a high volume of required computations \(order of 100 CPU*hrs for a typical bacterial genome\).

###The significance of proteogenomic annotation
Our pipeline is a tool for improving the
existing genomic annotations from available proteomics mass
spectrometry data. As most genome annotation pipelines consist of automated gene finding, they lack experimental validation of primary structure [[PMC2265698](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2265698/), [PMC2238897](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2238897/)], having to rely on DNA centric sources of data such as sequence homology, transcriptome mapping, codon frequency, etc. By incorporating the orthogonal set of data, proteogenomics is able to discover novel genes, post-translational modifications (PTMs) and correct the erroneous primary sequence annotations.

Quick Start
-----------
To test the pipeline, you can run it on a single workstation (preferably multicore one).

Make sure that your are in a directory where you have write access and enough free space (20G).

Follow the Installation instructions below, first editing environment scripts and CMake configuration variables if necessary, and then running the installation command: `proteogenomics/config/installPGP.sh -e htc -r PGP`

This will install the pipeline into a directory `PGP` in your working directory.

The installation procedure will automatically test the pipeline in a local (multicore) execution mode on reduced (subsampled) spectra data included in the package distribution and report the test status.

You can also download the full size version of the testing dataset (8.3GB in size) as:

    **URL** 

Unpack it in your current working directory:

    tar -zxf Cyanobacterium.synechocystis.PCC6803.tar.gz

And execute:
    
    PGP/proteogenomics/config/run.PGP.htc.sh Cyanobacterium.synechocystis.PCC6803  Cyanobacterium.synechocystis.PCC6803.results

This will take about 100 CPU*hrs. Following completion, check Cyanobacterium.synechocystis.PCC6803.results/GFFs to make sure that there are output GFF3 files where. You can load the GFF3 files and corresponding reference from Cyanobacterium.synechocystis.PCC6803/Databases/Genomic into any of the genomic browsers such as [NCBI Genome workbench](http://www.ncbi.nlm.nih.gov/tools/gbench/), Artemis [[PMC3278759](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3278759/)] or [CLC Genomics workbench](http://www.clcbio.com/products/clc-genomics-workbench/) for a visual inspection of the newly created annotations.

I). The PGP Algorithm
----------------------

The full protocol of our pipeline is described in detail in  [[PMC3219674](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3219674/)]. In that study, we have applied the pipeline to 46 genomes spanning eight bacterial and archaeal phyla across the tree of life. These diverse datasets facilitated the development of a robust approach for proteogenomics that is functional across genomes varying in %GC, gene content, proteomic sampling depth, phylogeny, and genome size. In addition to finding evidence for 682 novel proteins, 1336 new start sites, and numerous likely misannotations, we discovered sites of post-translational maturation in the form of proteolytic cleavage of 1175 signal peptides. The output files from this study are available at ([http://omics.pnl.gov/pgp/overview.php][2]).

The pipeline combines several open source proteomics tools from ([http://proteomics.ucsd.edu][3]) with our own post-processing steps.

Briefly, tandem mass spectra are searched by Inspect [[PMID:16013882](http://www.ncbi.nlm.nih.gov/pubmed/16013882)] against a translation of the genome and subsequently rescored with PepNovo [[PMID:15858974](http://www.ncbi.nlm.nih.gov/pubmed/15858974)] and MSGF [[PMC2689316](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2689316/)]. The pipeline  translates the DNA sequence in all six frames to generate a reference protein database.

Each stop to stop open reading frame (ORF) is included regardless of
coding potential. The concatenated decoy records for estimating the statistical significance of the predicted peptides are generated by shuffling each ORF.
Significant peptide/spectrum matches (PSM) are those with a pvalue of
e-10 or better, which leads to a peptide level false discovery rate of ~0.3%. All confident peptides are mapped onto their genomic location
(nucleotide coordinates) and grouped into sets within an ORF. 

The last stage of the pipeline employs five ORF filters. First, we remove low complexity peptides and peptides which are more than 750 bp from the next in-frame peptide. We remove ORFs which lack a uniquely mapping peptide or which lack a fully tryptic
peptide. Finally, we require two peptides per protein.

Peptides from ORFs which meet all of the above criteria are saved in the output GFF3 files in coordinates relative to the reference genome. The user can load these files into genomic browsers such as NCBI Genome workbench or Artemis for subsequent analysis.

### Example of interpreting the output

Peptides mapping to regions of the genome which lack a protein
annotation represent either novel genes, or 5' extensions of current
genes. For example, at 3.376 MB in the B. anthracis Sterne genome lies
the BAS3403 gene, the small spore protein Tlp. Many peptides from this
protein were discovered by the pipeline. In the open reading frame
directly upstream, the protocol detected other peptides which were not part of any currently annotated protein. Blasting this ORF revealed homology to another spore coat protein from B. cereus, not annotated
in any anthrax genome.

![Image of Novel Anthrax Sporecoat Protein](https://bitbucket.org/andreyto/proteogenomics/raw/master/docs.PGP/anthrax_sporecoat.png)

II). Parallelization strategy
-----------------
Our goal was to build a pipeline that would be portable across different parallel execution environments that users might encounter. The overall algorithm is embarassingly parallel for the most part. Specifically, it is possible to process each spectrum file independently throughout all computationally intensive stages of the algorithm. There is a global synchronization point in a middle of the workflow when the algorithm needs to build a histogram of all scores for p-value computation. We have selected a distributed workflow model where multiple serial processes are executed cuncurrently following a dependency graph defined by required input and output files. This model exploits inherent parallelism of the problem sufficiently well while being compatible with a wide variety of execution environments such as standalone multicore machines, high-throughput compute clusters and MPI clusters (more on the latter below).

We achieve the portability across different execution environments by using the Makeflow workflow execution engine ([http://nd.edu/~ccl/software/makeflow/][4]). Makeflow also provides a high degree of fault tolerance against compute node failures (in non-MPI execution mode) and restart capability in case of master node failure. 

Our installation procedure builds its own local copy of the CCTools package that contains Makeflow and associated backend executables.

Given a directory with input data, our pipeline's code generates a description of the workflow (in a language that is based on the `make` file syntax). Makeflow can then execute this workflow using different parallel "backends" selected by the user at run-time. While its backend implementations are straigtforward on high-throughput clusters (they submit each task to the batch queuing system) and multicore workstations (they run multiple subprocesses), the architecture is more complicated on HPC clusters that often have scheduling policies tuned to efficiently allocate only large MPI jobs. On such MPI clusters,  Makeflow supports execution of complex workflows composed of many interdependent serial jobs through a "glide-in" mechanism. Specifically, the user submits the MPI backend executable of Makeflow as a single parallel MPI job (e.g. with 100 MPI process ranks). Then, the user starts the Makeflow master on a single node. The individual MPI process ranks connect to the master, and the master farms out single serial tasks to the ranks following the order of workflow dependencies. A given backend rank executes each incoming task as a separate subprocess.

For both HTC and MPI clusters, the pipeline assumes that there is a shared file system mounted in the same location on all compute nodes as well as on the job submition ("login") node.

III). Installation
-----------------

### A). Dependencies that have to be present on your system to build and run the pipeline

The integrated installation procedure needs: 

> - [Git](http://git-scm.com/) version control system to checkout the source from the public repository on BitBucket ([https://bitbucket.org/andreyto/proteogenomics][5]). If you instead downlaod the package archive from BitBucket, you will not need Git
> - [CMake](http://www.cmake.org/) configuration and build utility (version 2.8 or higher)
> - *Wget* to download dependencies
> - *BASH* shell
> - C++ compiler (gcc or MPI wrappers, depending on the targeted execution environment)
> - Python interpreter and Python development libraries
> - [Boost Python](http://www.boost.org/) interface library
> - GNU Make
> - Internet connection. The installation procedure will try to download an archive with some additional dependencies (such as UCSD proteomics tools) from our BitBucket repository

All software components listed above are available as standard packages in major Linux distributions.

For run-time, you will need:

> - Java, Python, Boost Python shared libraries
> - MPI environment if MPI backend if used for execution

### B). Getting the sources
    
    git clone git@bitbucket.org:andreyto/proteogenomics.git

### C). Building and installing

A shell script `proteogenomics/config/installPGP.sh` is a wrapper that calls CMake configuration, build, test and install stages in a single run. Our CMake build procedure attempts to automatically install the pipeline itself as well as some of its dependencies. When you are executing this script, your current working directory (CWD) must be **outside** of the `proteogenomics` source directory. For example, your CWD can be a directory where you have issued the `git clone` command described above.

If the full build procedure does not work for you right away, you can temporarily comment out the latter stages in the `installPGP.sh` script in order to debug the configuration step of CMake.

The directory in which the software is installed by the `installPGP.sh` is
referred to as **PGP\_ROOT** in this document and the scripts. The PGP
code itself will be in `proteogenomics` subfolder of `PGP_ROOT` (referred to as  **PGP\_HOME**). All the
configuration files needed to set the environment for running the
pipeline will reside in the `config` subdirectory (referred to as **config\_dir**) under
PGP\_HOME.

**Usage:** ***installPGP.sh -h***

The options defining the cluster environment and the path to the
installation directory must be provided as indicated below.

    OPTIONS:

    -h help message

    -e computing environment name ("ranger" is the environment configured for XSEDE Ranger cluster; "htc" covers both HTC clusters and multi-core workstations). The computing environments are discussed in detail in README under "Customizing the build procedure".

    -r the installation directory (PGP\_ROOT)

**example:**

    installPGP.sh -e ranger -r <path to PGP_ROOT>

When the script exits successfully after installing the software, it will print messages to standard output, informing the user what to do next for executing the pipeline.

####Customizing the build procedure
Although CMake invoked by `installPGP.sh` will try to figure out the locations of necessary libraries on your system automatically, it might still need some explicit instructions that you can provide by editing some files under checked out **source** directory `proteogenomics/config` **before** `installPGP.sh` is invoked. This is especially likely needed on compute clusters, where the administrators often build multiple customized versions of libraries such as Boost and place them in non-standard locations (such as /usr/local/packages). Pointing CMake to the right Boost Python bindings library is typically the most challenging part of the build procedure on such heavily modified systems. There has to be a match between the versions of the selected Python interpreter, C++ compiler (or mpiCC wrapper) and Boost Python library. Additionally, the pipeline environment file might have to be modified by appending the directory where Boost shared libraries are located to the `LD_LIBRARY_PATH` environment variable.

We customize CMake variables through the so called "toolchain" files in `proteogenomics/config/`, with files named after a specific target computing environment. For example, there is `toolchain.ranger.cmake` that sets necessary CMake variables for XSEDE Ranger cluster. TACC Ranger was the XSEDE cluster where we have tested the MPI execution mode.

When the user selects "ranger" installation environment, our pipeline instantiates several scripts that make it more simple for the user to later execute the pipeline under the MPI Makeflow backend described above. At build time, the "ranger" option also uses "module" commands to make available the proper dependencies such as MPI compilers for building Makeflow as well as Boost libraries. "Module" command is a standardised user environment management script used by XSEDE clusters ([https://www.xsede.org/software-environments][7]). 

On XSEDE systems, the users will have to modify the details of specific package versions activated my the "module" command in order to adapt the building and execution environment to their specific cluster. If they decide to call their cluster `my_cluster`, they can copy every file that has a name like `*.ranger.*` in the source `proteogenomics/config` directory to a file that has the word `ranger` replaced with the word `my_cluster`, modify the necessary content in the new files, and then run the installation procedure as `installPGP.sh -e my_cluster -r <path to PGP_ROOT>`. Alternatively, the users can simply edit the `*.ranger.*` files and use `installPGP.sh -e ranger -r <path to PGP_ROOT>`.

Using other Makeflow backends such as SGE serial job submission or multicore workstation is much more simple, and those cases are covered by the sample execution scripts instantiated when "htc" environment is selected during installation. 

To customize the build procedure on HTC systems, the users can edit `*.htc.*` files.

IV). Running the pipeline
--------------------------

At high level, the pipeline execution consists of two stages. The first, serial stage prepares the input data and generates the workflow file for the second stage. The second stage performs the analysis in multiple parallel processes (with some internal barrier synchronisation steps, for example, for computing the *p-values*). 

Because scheduling policies and node availability are typically vastly different between serial and parallel jobs on MPI clusters like XSEDE systems, we have given to the user the control over launching these two stages in these environments as described below. For HTC environments, these two stages are executed consecutively by a single script.

The users can always modify the parameters of Makeflow execution and backend job submission by following the Makeflow User Manual ([http://www3.nd.edu/~ccl/software/manuals/makeflow.html][8]), possibly using our template scipts as the starting points. Our scripts have extensive annotations in inline comments. Some site-specific tuning of job submission scripts will likely to be required on any computational cluster, considering the multitude of customizations of both the operating and batch systems that cluster administrators typically employ.

The following steps must be taken to execute the proteogenomics analysis:

**A). Data Requirements: Genomic & Proteomic Data.**

1). Annotated nucleotide (mandatory) and protein (optional) sequence
files (e.g. **gbk**, **fna** and **faa** RefSeq files from NCBI ([http://www.ncbi.nlm.nih.gov/refseq/][9])).

2). Mass spectrometric data (spectra) in mzxml format (mzXML conversion tools can be found at
[http://tools.proteomecenter.org/wiki/index.php?title=Formats:mzXML][10]),

**B). Structured Input Data Directory:**

1). The input data mentioned under **'A'** should be organized in a
specific directory structure as shown by the listings of the directories below. This can be easily accomplished by uncomressing the `template.directory.tar.gz` provided in the data subdirectory of the
proteogenomics folder (PGP\_HOME) downloaded by the pipeline. Assuming the annotated genomic data comes from speciesX, uncompress the template directory inside speciesX folder:

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

2). Put your RefSeq fna and gbk files for speciesX into
`Databases/Genomic` folder and faa file(s) - into
`Databases/Predictions` folder. Any additional protein databases including common contaminants such as trypsin and keratin (but not the DNA), which also need to be searched with Inspect, go into the `Databases/Proteomic` folder. The template directory already has some
unpacked Common.RS files in this directory (as shown below).

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

3). Your spectra mzxml file(s) should be compressed in gz tar format and placed in the `mzxml` directory.

4). The remaining empty directories (some optional) serve as
placeholders to be populated by the data generated during the analysis. Some of them are the suggested locations for the last, manual curation stages of the proteogenomic re-annotation process aimed to culminate in the submission of results to NCBI Peptidome and RefSeq as described in the pipeline protocol publication []. If the user is only interested in the results of the automatic annotation but not in the NCBI submission, those directories can be ignored.

**DerivedData** - copied from output directory (see below). Specific
files or results that should be produced with every dataset

**GFFs** - a sub-directory copied from output directory (see below) -
contains the gff file with peptide annotations for the genome

**OrthologyClusters** - folder with MSA (multiple sequence alignments)
of orthologous clusters from related organisms

**Publications** - if there are relevant publications related to the
data, they go here.

**RawZips** - holds the raw Inspect and PepNovo rescored data in bz2
files

**PeptidomeSubmission** – the Peptidome submission needs meta data file, PSM results file, spectra & the results file - as explained as below.

**Peptidome** data generation: meta data generated by the user

**PepXML** - these are pepXML conversions of Inspect's output (PSMs)

**mzxml** - spectra files (mzXML.tar.gz)

**ProteinInference** - the protein inference as made by the PGP script

**RefseqSubmission** - originals submitted - contains new annotations generated by the PGP

**C). Execution on XSEDE MPI clusters**

**i. Preparation of Data & the Workflow**

In this step, the input data for the pipeline are prepared and a workflow specifying the various PGP analyses to be executed on the target distributed computed environment is generated.

Copy a config_dir/prepPGPdata.ranger.qsub template script to your working directory and edit it to reflect your actual data location and batch system parameters as per the instructions below.

Make sure to do the following before executing the qsub script (substitute actual location for `config_dir` below):

a) `source config_dir/pgp_makeflow_env_master.sh`
    
b) For example, SGE-specific options `-A (account name), -pe (number of nodes and cores), -l h_rt (requested
run-time)` have to be properly modified in the qsub script on a XSEDE cluster. A "serial" queue is specified in the script - 
edit the script if such queue does not exist on your cluster. This serial job needs only the minimum number of cores allowed on the particular cluster. 

c) values for INPUT (--input-dir path for input directory in the format described above) and OUT
(--work-dir path for output results) are set in the qsub script. These are MANDATORY.

**Usage**:

Submit the job to your batch system. In case of SGE or PBS, it will look like:

    qsub prepPGPdata.ranger.qsub
    
Alternatively, you can execute this script locally on your login node if this is tolerated by the cluster user policy. In that case you do not have to edit any batch system related options inside the script, but still need to edit the data locations:

    bash prepPGPdata.ranger.qsub

-   After successful execution of the qsub script, the **log** file
    should contain the instructions (at the end) needed for running the next qsub script (runPGP.ranger.qsub) for PGP analysis.   
      

**ii. Running the Analysis:** 

This qsub script performs the actual proteogenomics analysis and the post-processing of
the resulting data according to the workflow generated in the previous
step. Copy the template script `config_dir/runPGP.ranger.qsub` to your working directory and modify the options inside the script to fit your execution environment.

**Note:** Make sure to do the following before submitting **runPGP.ranger.qsub** script:

a) `source config_dir/pgp_makeflow_env_master.sh`. The reason for this to be done on the login node of the XSEDE cluster before submitting the script is because the XSEDE `module` commands must be executed only on the login node. The batch system is configured to propagate the resulting environment variables to the compute nodes. This process might have some quirks depending on your specific cluster. You might have to tune this script accordingly.

b) SGE-specific options `-A (account name), -pe (number of nodes and cores), -l h_rt (run-time
requested)` are properly modified in the qsub script. These are MANDATORY. Also edit the queue name if necessary. Replace the options if your cluster uses a different queuing system, such as PBS or SLURM. Also make sure that the commands for starting MPI jobs inside the script match your environment (e.g. XSEDE clusters that have Infiniband interconnect use `ibrun` instead of `mpirun`).

c) change directory (‘cd’) to `--work-dir/<speciesX>` (as defined in the previous qsub script `prepPGPdata.ranger.qsub`. 

**Usage:**

`qsub <edited script location>/runPGP.ranger.qsub`

-   if the script is executed successfully, the results of the analysis
    are written to --work-dir/speciesX ( --work-dir = OUT variable as
    defined in prepPGPdata.ranger.qsub).

**D). Execution on HTC clusters or multicore workstations**

In this execution environments, there is a single script that runs the serial data preparation stage followed by the parallel processing stage. There is no need to edit the template script. Rather, you should execute the script `config_dir/runPGP.htc.sh`, directly passing to it the necessary command line arguments described below.

**Usage:**

    config_dir/run.PGP.htc.sh INPUT OUT [Makeflow arguments]
    
where `INPUT` is your prepared directory with input data, and `OUT` is the intended location of the output. `OUT` should be on a shared file system if you are running on a cluster. The optional Makeflow arguments are passed directly to the Makeflow instance that will be running the data processing stage.

On a single workstation without a batch system, the command can be simply

    config_dir/run.PGP.htc.sh INPUT OUT

In that case, Makeflow will run the number of concurrent subprocesses that equals the number of cores on the machine.

On SGE cluster, the command might look like

    config_dir/run.PGP.htc.sh INPUT OUT -T sge -B 'SGE options'

where the string `SGE options` (must be in quotes as shown above) is passed verbatim to the SGE when the Makeflow submits each task in the workflow. This string must contain the options that you would pass to a `qsub` in order to run a single job on your system. For example, on JCVI internal cluster the full command would look like

    config_dir/run.PGP.htc.sh INPUT OUT -T sge -B '-P 0000 -b n -S /bin/bash'

You can look at the Makeflow manual for possible command line arguments if you want to use other backends such as Condor or modify the behavior of Makeflow (e.g. to constrain the number of concurrent jobs).

In the commands above, both the data preparation as well as the controlling process of the Makeflow (the "master") will run on the login node, and the script `run.PGP.htc.sh` will keep running until the entire workflow has finished. It will use very little resources during the data processing stage, merely farming out tasks to the batch system. You can also submit the command above itself to your batch system through a standard (e.g. qsub) mechanism in case you do not want to wait for it to finish on the login node. This would require that your compute nodes are allowed to submit new jobs themselves, because Makeflow will be submitting new jobs from the node where the master script is executing.

V). Results
------------

**Structured Output Data Directory:** When the PGP analysis is complete,
in addition to some log files, the following directories are created in
-–work-dir/speciesX directory.

**Databases:** the 6-frame translated refseq genomic data copied from
input directory.

**mzxml:** MS spectral data copied from input directory.

**DerivedData:** contains the proteogenomics analysis results and
various reports for each input dataset.

**GFFs:** a sub-directory with the gff file containing peptide
annotations of the genome generated from mapping the proteomics spectra. **This can be considered the main output of the pipeline.** See [[PMC3219674](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3219674/)] regarding the interpretation and downstream use of these annotations.

**ResultsX:** the directory contains results from Inspect analysis.

**pepnovo:** the directory contains results from PepNovo analysis.

**pvalue10H:** the sub-directory with the output of PValue.py (with -p
0.1, -H).

**msgfOfPepnovo:** the sub-directory with the output of MS\_GF
validation.

**jobs, Done, output:** other directories with pipeline working data
(could be empty).

VI). Testing
-----------

We have created a sample archive of an input directory completely populated with spectra and genomic data for *Cyanobacterium synechocystis PCC6803*. You can download it **wget URL**, unpack with `tar -xzf` and supply the path of the resulting directory as input to the pipeline. Note that this is a real size dataset, and it will take about 100 CPU*hrs to process it.

There is also the input archive for the same genome (as referenced in the Quick Start section), with spectra data drastically subsampled in order to provide a quick way of testing that the pipeline can run to completion after the installation. Due to low coverage in this artificially constructed sample, the output should not be used for any biological interpretation.


  [1]: https://www.xsede.org/
  [2]: http://omics.pnl.gov/pgp/overview.php
  [3]: http://proteomics.ucsd.edu
  [4]: http://nd.edu/~ccl/software/makeflow/
  [5]: https://bitbucket.org/andreyto/proteogenomics
  [6]: Customizing%20the%20build%20procedure
  [7]: https://www.xsede.org/software-environments
  [8]: http://www3.nd.edu/~ccl/software/manuals/makeflow.html
  [9]: http://www.ncbi.nlm.nih.gov/refseq/
  [10]: http://tools.proteomecenter.org/wiki/index.php?title=Formats:mzXML