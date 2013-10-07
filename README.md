<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<HTML>
<HEAD>
    
<META HTTP-EQUIV="CONTENT-TYPE" CONTENT="text/html; charset=windows-1252">
<TITLE></TITLE>
    
<META NAME="GENERATOR" CONTENT="OpenOffice.org 3.3  (Win32)">
<META NAME="AUTHOR" CONTENT="Venepally, Pratap">
<META NAME="CREATED" CONTENT="20120831;20490000">
<META NAME="CHANGEDBY" CONTENT="Venepally, Pratap">
<META NAME="CHANGED" CONTENT="20120904;19490000">
<META NAME="AppVersion" CONTENT="14.0000">
<META NAME="Company" CONTENT="J. Craig Venter Institute">
<META NAME="DocSecurity" CONTENT="0">
<META NAME="HyperlinksChanged" CONTENT="false">
<META NAME="LinksUpToDate" CONTENT="false">
<META NAME="ScaleCrop" CONTENT="false">
<META NAME="ShareDoc" CONTENT="false">
<STYLE TYPE="text/css">
    <!--
        @page { margin-left: 1in; margin-right: 1in; margin-top: 1in; margin-bottom: 0.5in }
        @page:first { }
        P { margin-bottom: 0.08in; direction: ltr; color: #000000; widows: 2; orphans: 2 }
        A:link { color: #0000ff; so-language: zxx }
    -->
    </style>
</HEAD>
<BODY LANG="en-US" TEXT="#000000" LINK="#0000ff" DIR="LTR">
<P ALIGN=CENTER STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=5 STYLE="font-size: 20pt"><B>Automated Proteogenomics
Pipeline<BR></B></FONT><FONT SIZE=4 STYLE="font-size: 16pt"><B>Sam Payne
& Andrey Tovchigrechko, JCVI.</B></FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>The semi-automated proteogenomics pipeline (PGP) described
here is freely accessible to the scientific community and is intended to
be used on TeraGrid, the National Science Foundation-supported Sun
Constellation Linux Cluster (Ranger), which was recently integrated into
XSEDE project (https://www.xsede.org/). In addition, the pipeline can
also be run on local SGE clusters or even on single high-memory,
multi-CPU Unix systems.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B><BR>Use:</B></FONT><FONT SIZE=3> The pipeline is
primarily intended as a tool for improving the existing genomic
annotations from available proteomics and mass spectrometry data. The
resulting improvements might include discovering novel genes,
post-translational modifications (PTMs) and correction of the erroneous
primary sequence annotations.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=4><B>I). The PGP Process:</B></FONT><FONT SIZE=4> The
pipeline runs the following steps, developed at
http://proteomics.ucsd.edu, in sequence:<BR> </FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>1). InsPecT:</B></FONT><FONT SIZE=3> Inspect generates
peptide/spectrum matches (PSM) by comparing MS/MS data with 6-frame
translated amino acid sequence database [1]. It uses peptide sequence
tags (PSTs) to filter the database and is designed to detect
post-translational modifications (PTM) [2, 3].<BR></FONT><BR>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>2). PepNovo:</B></FONT><FONT SIZE=3> In addition to
serving as a de novo sequencing program based on probabilistic network
modeling [4], PepNovo can also verify the database search results
produced by tools such as InsPecT in the previous step. PepNovo rescores
the raw results from InsPecT, replacing the MQScore and delta score
fields with the scores obtained with the rank models.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>3). Validation:</B></FONT><FONT SIZE=3> Statistical
validation of PSMs, to filter out false matches, are performed by
p-value estimation using two different algorithms: Initially, the PSMs
identified by Inspect are assigned p-value based on the score
distribution of hits to the target and decoy databases [5, 6] and those
with p \< 0.01 are retained. These results are then processed by a
second algorithm (MS\_GF) [7], which provides its own p-value
estimations for identifying false hits.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>4). Post-processing:</B></FONT><FONT SIZE=3> The final
stage performs proteogenomics analyses and saves the annotation results.
Some of the analyses include: A). mapping peptides to the genome, B).
clustering and filtering of the peptides, C). identifying mis-predicted
and novel proteins, D). protein inference & E). Cleavage Analysis,
etc.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=4><B>II). Installation:</B></FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>A). The pipeline requires the following Computer
resources & Software: </B></FONT><FONT SIZE=3>Unix BASH shell
environment with optional access to Xsede (TeraGrid (Ranger) or other
SGE cluster;\
Git utility for downloading proteogenomics package from Bitbucket
repository
(</FONT><A HREF="https://bitbucket.org/">https://bitbucket.org/</A><FONT SIZE=3>);
Perl, JAVA & Python languages. </FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>B).</B></FONT><FONT SIZE=3>
</FONT><FONT SIZE=3><B>installPGP.sh:</B></FONT><FONT SIZE=3> : To
install the package and other dependencies required by the pipeline, i).
download the installPGP.sh script
(</FONT><FONT COLOR="#1f497d"><FONT SIZE=3><B>wget
</B></FONT></FONT><A HREF="https://bitbucket.org/andreyto/proteogenomics/raw/320bd5e412f0/config/installPGP.sh"><FONT SIZE=3><B>https://bitbucket.org/andreyto/proteogenomics/raw/320bd5e412f0/config/installPGP.sh</B></FONT></A><FONT COLOR="#1f497d"><FONT SIZE=3><B>),</B></FONT></FONT><FONT SIZE=3>
ii). make it executable and iii). run it from a directory of your choice
as shown below under
</FONT><FONT COLOR="#1f497d"><FONT SIZE=3><B>Usage</B></FONT></FONT><FONT SIZE=3>.
</FONT><FONT COLOR="#ff0000"><FONT SIZE=3>It only needs to be run
once.</FONT></FONT><FONT SIZE=3> The directory in which the software is
installed by the</FONT><FONT SIZE=3><B>
installPGP.sh</B></FONT><FONT SIZE=3> is referred to as
</FONT><FONT SIZE=3><B>PGP\_ROOT</B></FONT><FONT SIZE=3> in the document
and the script. The PGP software is in
</FONT><FONT COLOR="#ff0000"><FONT SIZE=3><B>'proteogenomics'</B></FONT></FONT><FONT COLOR="#ff0000"><FONT SIZE=3>
</FONT></FONT><FONT SIZE=3>subfolder
(</FONT><FONT SIZE=3><B>PGP\_HOME</B></FONT><FONT SIZE=3>). All the
configuration files needed to set the environment for running the
pipeline reside in the
</FONT><FONT COLOR="#ff0000"><FONT SIZE=3><B>config</B></FONT></FONT><FONT COLOR="#ff0000"><FONT SIZE=3>
</FONT></FONT><FONT SIZE=3>subdirectory
(</FONT><FONT SIZE=3><B>config\_dir</B></FONT><FONT SIZE=3>) under
PGP\_HOME. </FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT><FONT COLOR="#1f497d"><FONT SIZE=3><B>Usage:
</B></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3>installPGP.sh
\<options\></FONT></FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>The options defining the cluster environment and the path
to the installation directory must be provided as indicated below.
</FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> OPTIONS:</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">       
<FONT SIZE=3> -h help message</FONT>
</P>
<P STYLE="margin-left: 1.88in; text-indent: -0.38in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>-e computing environment name (‘ranger’ is the currently
supported environment, corresponding to TACC Ranger on XSEDE
site)</FONT>
</P>
<P STYLE="margin-left: 1.81in; text-indent: -0.31in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>-r the installation directory (PGP\_ROOT)</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT><FONT SIZE=3><B>example:</B></FONT><FONT SIZE=3>
</FONT>
</P>
<P STYLE="margin-left: 0.5in; text-indent: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT COLOR="#c00000"><FONT SIZE=3><B>installPGP.sh -e ranger -r \<path
to PGP\_ROOT\> &</B></FONT></FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<UL>
    <LI><P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT COLOR="#1f497d"><FONT SIZE=3>When the script exits successfully
after installing the software, it will print messages to STDOUT,
informing what to do to next for executing the pipeline. </FONT></FONT>
</P>
</UL>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=4><B>III). Running the pipeline</B></FONT><FONT SIZE=4>:
Since the run-environment and the job submission requirements might
differ on different Unix clusters, a semi-automated pipeline with
minimal user input is created. This consists of two separate batch jobs.
The prerequisites and the detailed instructions are described
below.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT><FONT SIZE=3><B>A). Data Requirements: Genomic &
Proteomic Data.</B></FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>1). Annotated nucleotide (mandatory) and protein (optional)
sequence files (e.g. </FONT><FONT SIZE=3><B>gbk</B></FONT><FONT SIZE=3>,
</FONT><FONT SIZE=3><B>fna</B></FONT><FONT SIZE=3> and
</FONT><FONT SIZE=3><B>faa</B></FONT><FONT SIZE=3> refseq files from
NCBI). </FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>2). Mass spectrometric data (spectra) in mzxml format
(mzXML conversion tools can be found at
http://tools.proteomecenter.org/wiki/index.php?title=Formats:mzXML),</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="text-indent: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>B). Organized Input Data Directory:</B></FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>1). The input data mentioned under
</FONT><FONT SIZE=3><B>‘A</B></FONT><FONT SIZE=3>’ should be organized
in a specific directory structure as shown by the listings of the
directories below. </FONT><FONT COLOR="#1f497d"><FONT SIZE=3><B>This can
be easily accomplished by uncomressing the 'template.directory.tar.gz'
provided in the data subdirectory of the proteogenomics folder
(PGP\_HOME) downloaded by the pipeline.</B></FONT></FONT><FONT SIZE=3>
</FONT><FONT COLOR="#1f497d"><FONT SIZE=3><U>Assuming the annotated
genomic data comes from speciesX, uncompress the template directory
inside speciesX folder:</U></FONT></FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> user@assembly [\~/speciesX]% ls -l</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> total 32</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> drwxr-x--- 5 user tigr 4096 Aug 30 12:01 Databases</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> drwxr-x--- 3 user tigr 4096 Aug 30 12:01
DerivedData</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> drwxr-x--- 2 user tigr 4096 Aug 30 12:01
OrthologyClusters</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> drwxr-x--- 3 user tigr 4096 Aug 30 12:01
PeptidomeSubmission</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> drwxr-x--- 2 user tigr 4096 Aug 30 12:01
Publications</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> drwxr-x--- 2 user tigr 4096 Aug 30 12:01 RawZips</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> drwxr-x--- 4 user tigr 4096 Aug 30 12:01
RefSeqSubmission</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> drwxr-x--- 2 user tigr 4096 Aug 30 12:01 mzxml</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>2). The available refseq fna and gbk records should be
saved to Database/Genomic folder and faa file(s) should be placed in
Database/Predictions folder. Any additional protein databases including
common contaminants such as trypsin and keratin (not DNA to be put into
six frames), which need to be searched with Inspect, go into the
Databases/Proteomic folder. The template directory already has some
unpacked Common.RS files in this directory (as shown below).</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> user@assembly [\~/tmp/db/Databases]: [12:07:22]% ls
-l</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> .:</FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>total 12</FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>drwxr-x--- 2 user tigr 4096 Feb 19 2010 Genomic</FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>drwxr-x--- 2 user tigr 4096 Feb 19 2010 Predictions</FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>drwxr-x--- 2 user tigr 4096 Feb 19 2010 Proteomic</FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>./Genomic:</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> total 0</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> ./Predictions:</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> total 0</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> ./Proteomic:</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> total 16</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> -rw-r----- 1 user tigr 1472 Feb 19 2010
Common.RS.index</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> -rw-r----- 1 user tigr 8254 Feb 19 2010
Common.RS.trie</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>3). The spectra mzxml file(s) should be compressed in gz
tar format and placed in the
</FONT><FONT SIZE=3><B>mzxml</B></FONT><FONT SIZE=3> directory.</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>4). The remaining empty directories (some optional) serve
as placeholders to be populated by the data generated by the analysis. A
brief description of each is as follows:</FONT>
</P>
<P STYLE="margin-left: 1.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>DerivedData</B></FONT><FONT SIZE=3> - copied from output
directory (see below). Specific files or results that should be produced
with every dataset</FONT>
</P>
<P STYLE="margin-left: 2in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>GFFs - a sub-directory copied from output directory (see
below) - contains the gff file with peptide annotations for the
genome</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-left: 1.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>OrthologyClusters</B></FONT><FONT SIZE=3> - folder with
MSA (multiple sequence alignments) of orthologous clusters from related
organisms</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-left: 1.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>Publications</B></FONT><FONT SIZE=3> - if there are
relevant publications related to the data, they go here.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT><FONT SIZE=3><B>RawZips</B></FONT><FONT SIZE=3> -
holds the raw Inspect and PepNovo rescored data in bz2 files</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-left: 1.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>PeptidomeSubmission</B></FONT><FONT SIZE=3> – the
submission needs meta data file, PSM results file, spectra & the results
file - as explained as below.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT><FONT SIZE=3><B>Peptidome</B></FONT><FONT SIZE=3>
data generation: meta data generated by the user</FONT>
</P>
<P STYLE="margin-left: 2in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>PepXML</B></FONT><FONT SIZE=3> - these are pepXML
conversions of Inspect's output (PSMs)</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT><FONT SIZE=3><B>mzxml</B></FONT><FONT SIZE=3> -
spectra files (mzXML.tar.gz)</FONT>
</P>
<P STYLE="margin-right: -0.25in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>
</FONT><FONT SIZE=3><B>ProteinInference</B></FONT><FONT SIZE=3> - the
protein inference as made by the PGP script</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT><FONT SIZE=3><B>RefseqSubmission</B></FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> originals</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> submitted - contains new annotations generated by the
PGP</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>C). Preparation of Data & the Workflow
(prepPGPdata.ranger.qsub):</B></FONT><FONT SIZE=3> In this step, the
input data for the pipeline are prepared and a workflow specifying the
various PGP analyses to be executed on the distributed computed
environment (clusters) is generated. This is accomplished by using
</FONT><FONT COLOR="#1f497d"><FONT SIZE=3><B>makeflow program
</B></FONT></FONT><FONT SIZE=3>which enables analysis to be performed on
“different systems including single multicore machines, Condor and SGE
batch systems or
</FONT><FONT COLOR="#000000"><FONT FACE="Times, serif">bundled</FONT></FONT><FONT COLOR="#000000"><FONT FACE="Times, serif"> </FONT></FONT><A HREF="http://nd.edu/~ccl/software/workqueue"><FONT COLOR="#000080"><FONT FACE="Times, serif">Work
Queue</FONT></FONT></A><FONT COLOR="#000000"><FONT FACE="Times, serif"> </FONT></FONT><FONT COLOR="#000000"><FONT FACE="Times, serif">system”
</FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3><B>(</B></FONT></FONT><A HREF="http://nd.edu/~ccl/software/makeflow/">http://nd.edu/\~ccl/software/makeflow/</A><FONT COLOR="#1f497d"><FONT SIZE=3><B>).
</B></FONT></FONT><FONT SIZE=3> The
</FONT><FONT COLOR="#1f497d"><FONT SIZE=3><B>makeflow</B></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3>
</FONT></FONT><FONT SIZE=3>program comes as part of the
</FONT><FONT COLOR="#1f497d"><FONT SIZE=3><B>cctools</B></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3>
</FONT></FONT><FONT SIZE=3>package that is already installed in the
previous step
(</FONT><FONT COLOR="#1f497d"><FONT SIZE=3><B>installPGP.sh</B></FONT></FONT><FONT SIZE=3>).
</FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT COLOR="#1f497d"><FONT SIZE=3><I><B>Note:</B></I></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3><I>
the default options in this qsub script (prepPGPdata.ranger.qsub) are
'ranger-specific'. If the cluster environment is not Ranger, modify
these options as required by the user's own environment.
</I></FONT></FONT><FONT COLOR="#ff0000"><FONT SIZE=3><I><B>Make sure to
do the following before executing the qsub script:
</B></I></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3><I><B>a).
source config\_dir/pgp\_makeflow\_env\_master.sh b). SGE-specific
options -A (account name), -pe (number of nodes and cores), & -l h\_rt
(requested run-time) are properly modified
</B></I></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3><I><U><B>in the
qsub
script</B></U></I></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3><I><B>;
c). values for DB (--input-dir path for annotated genome) and OUT
(--work-dir path for output results) are set
</B></I></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3><I><U><B>in the
qsub
script</B></U></I></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3><I><B>.
These are MANDATORY.</B></I></FONT></FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT COLOR="#1f497d"><FONT SIZE=3><B>Usage</B></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3>:
</FONT></FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT COLOR="#c00000"><FONT SIZE=3><B>qsub
config\_dir/prepPGPdata.ranger.qsub</B></FONT></FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<UL>
    <LI><P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT COLOR="#1f497d"><FONT SIZE=3>After successful execution of the
qsub script, the
</FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3><B>log</B></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3>
file should end with the instructions needed for running the next qsub
script (runPGP.ranger.qsub) for PGP analysis.</FONT></FONT><FONT SIZE=3>
<BR></FONT><BR>
</P>
</UL>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>D). Running the Analysis
(runPGP.ranger.qsub):</B></FONT><FONT SIZE=3> This qsub script performs
the various proteogenomics analyses and the post-processing of the
resulting data according to the workflow generated in the previous step.
The default options in this qsub script are 'ranger-specific'. if the
environment is not Ranger cluster, modify these options as required by
the user's own environment.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT COLOR="#1f497d"><FONT SIZE=3><I><B>Note:
</B></I></FONT></FONT><FONT COLOR="#ff0000"><FONT SIZE=3><I><B>Make sure
to do the following before executing
</B></I></FONT></FONT><FONT COLOR="#ff0000"><FONT SIZE=3><B>runPGP.ranger.qsub</B></FONT></FONT><FONT COLOR="#ff0000"><FONT SIZE=3><I><B>
script:
</B></I></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3><I><B>a).
source config\_dir/pgp\_makeflow\_env\_master.sh b). SGE-specific
options -A (account name), -pe (number of nodes and cores), & -l h\_rt
(run-time requested) are properly modified
</B></I></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3><I><U><B>in the
qsub
script</B></U></I></FONT></FONT><FONT COLOR="#1f497d"><FONT SIZE=3><I><B>;
c). change directory (‘cd’) to --work-dir/speciesX' (as defined in the
previous qsub script ‘ prepPGPdata.ranger.qsub’. These are
MANDATORY.</B></I></FONT></FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT COLOR="#1f497d"><FONT SIZE=3><B>Usage:</B></FONT></FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT COLOR="#c00000"><FONT SIZE=3><B>qsub
config\_dir/runPGP.ranger.qsub</B></FONT></FONT><FONT COLOR="#c00000"><FONT SIZE=3>
</FONT></FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<UL>
    <LI><P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT COLOR="#1f497d"><FONT SIZE=3>if the script is executed
successfully, the results of the analysis are written to
--work-dir/speciesX ( --work-dir = OUT variable as defined in
prepPGPdata.ranger.qsub).</FONT></FONT>
</P>
</UL>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=4><B>IV). Results:</B></FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-left: 0.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>Organized Output Data Directory:</B></FONT><FONT SIZE=3>
When the PGP analysis is complete, in addition to some log files, the
following directories are created in –work-dir/speciesX
directory.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>Databases:</B></FONT><FONT SIZE=3> the 6-frame
translated refseq genomic data copied from input directory.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT><FONT SIZE=3><B>mzxml:</B></FONT><FONT SIZE=3> MS
spectral data copied from input directory.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-left: 1in; margin-bottom: 0in; line-height: 100%">
<A NAME="_GoBack"></A>
<FONT SIZE=3><B>DerivedData:</B></FONT><FONT SIZE=3> contains the
proteogenomics analysis results and various reports for each input
dataset.</FONT>
</P>
<P STYLE="margin-left: 1.5in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>GFFs</B></FONT><FONT SIZE=3> - a sub-directory with the
gff file containing peptide annotations of the genome.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT><FONT SIZE=3><B>ResultsX:</B></FONT><FONT SIZE=3>
the directory contains results from Inspect analysis. </FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT><FONT SIZE=3><B>pepnovo:</B></FONT><FONT SIZE=3>
the directory contains results from PepNovo analysis.</FONT>
</P>
<P STYLE="margin-right: -0.19in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3> </FONT><FONT SIZE=3><B>pvalue10H:</B></FONT><FONT SIZE=3>
the sub-directory with the output of PValue.py (with -p 0.1, -H).</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>
</FONT><FONT SIZE=3><B>msgfOfPepnovo:</B></FONT><FONT SIZE=3> the
sub-directory with the output of MS\_GF validation.</FONT>
</P>
<P STYLE="margin-left: 2.38in; margin-right: -0.13in; text-indent: -1.38in; margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3><B>jobs, Done, output:</B></FONT><FONT SIZE=3> other
directories with analysis tacking data (could be empty). </FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=4><B>V). Citations:</B></FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>[1]. S. Tanner, H. Shu, A. Frank, L.Wang, E. Zandi, M.
Mumby, P.A. Pevzner, and V. Bafna. Inspect: Fast and accurate
identification of post-translationally modified peptides from tandem
mass spectra. Anal. Chem., 77(14):4626–4639, 2005.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>[2]. Identification of Post-translational Modifications via
Blind Search of Mass-Spectra. Dekel Tsur, Stephen Tanner, Ebrahim Zandi,
Vineet Bafna, Pavel A. Pevzner. Nature Biotechnology 23, 1562-2567
(2005). </FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>[3]. Frank, A., Tanner, S., Bafna, V. and Pevzner, P.
"Peptide sequence tags for fast database search in mass-spectrometry",
J. Proteome Res. 2005 Jul-Aug;4(4):1287-95.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>[4]. Frank, A. and Pevzner, P. "PepNovo: De Novo Peptide
Sequencing via Probabilistic Network Modeling", Analytical Chemistry
77:964-973, 2005.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>[5]. Elias J.E., Gygi S.P. Target-decoy search strategy for
increased confidence in large-scale protein identifications by mass
spectrometry. Nat Methods 4,(3):207-14 (2007)</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>[6]. Castellana NE, Payne SH, Shen Z, Stanke M, Bafna V, et
al. Discovery and revision of Arabidopsis genes by proteogenomics. Proc
Natl Acad Sci U S A. 2008;105:21034–21038.</FONT>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<BR>
</P>
<P STYLE="margin-bottom: 0in; line-height: 100%">
<FONT SIZE=3>[7]. Sangtae Kim, Nitin Gupta and Pavel Pevzner. Spectral
Probabilities and Generating Functions of Tandem Mass Spectra: A Strike
against Decoy Databases. J. Proteome Res., 7 (8), 3354-3363,
2008.</FONT>
</P>
<P STYLE="margin-bottom: 0in">
<BR>
</P>
</BODY>
</HTML>


