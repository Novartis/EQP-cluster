Exon quantification pipeline (EQP-cluster)
==========================================

EQP-cluster is a Unix-based RNA-seq quantification pipeline which
takes a set of sample Fastq files as input, aligns them against
reference files, and generates files with
the gene, exon, or junction counts for each sample.

It provides scripts to facilitate the distributed execution of
alignment and quantification operations for each of the samples and
is designed to support the Univa\* Grid Engine (UGE)
batch-queuing/scheduling system (job submission via `qsub`).

If only the quantification step is of interest and your Fastq files are
pre-aligned, then you can also use  
[https://github.com/novartis/EQP-QM](https://github.com/novartis/EQP-QM)  
which is a Unix based RNA-seq quantification module; it uses SAM/BAM
genome alignment files as input and creates gene, exon, and junctions counts.

--------------
\*: Univa® is a registered trademark of Univa Corporation

Version
-------

This document corresponds to EQP-cluster version 2.1.0


Installation of EQP-cluster
---------------------------

After downloading and unpacking the files from GitHub, the download
directory should contain the following entries:  
`LICENSE.txt              bin/                     exon-pipeline-scripts/`
`README.md                exon-pipeline-files.tgz  samples.tgz`  

To complete the installation of EQP, go to the `bin` directory and
call the script `create-bowtie2-indices.sh`:  
`cd bin`  
`create-bowtie2-indices.sh`

Note that `create-bowtie2-indices.sh` requires that the programs
`samtools` and `bowtie2-build` are installed. See below for the
version requirements.

The script `create-bowtie2-indices.sh` extracts the files from the
tar file `exon-pipeline-files.tgz` into the directory `exon-pipeline-files`
and launches three calls to `bowtie2-build` in the background to construct
the Bowtie2 indices of the genome Fasta file, the transcript Fasta file,
and a custom junction Fasta file. The computation of the Bowtie2 indices
will take several hours.

**Limitations**:  
Note that currently EQP-cluster makes use of human Ensembl76 reference files
(on which the files contained in the directory `exon-pipeline-files` are based)
and is, thus, configured to quantify only human data; furthermore, due to
the way the custom junction Fasta file is created EQP-cluster
works best with reads of length 101bp or less.

Other reference file package will be made available in the future.

Once the Bowtie2 indices are computed, EQP-cluster is ready to be used.
The directory containing the subdirectories `bin`, `exon-pipeline-scripts`,
and `exon-pipeline-files` will be called `<project directory>` in the following.


Dependencies
------------

* Python (>= version 2.6.5, imported libraries: copy, gettext, gzip,
  numpy, os, re, sets, sys, textwrap, warnings)
* Java (>= version 1.6)
* samtools (>= version 0.1.17)
* bedtools (>= version 2.24.0)
* bowtie2  (>= version 2.0.5)
* Univa® Grid Engine


Running EQP-cluster
-------------------

For the processing of samples EQP relies one a fairly rigid directory
structure in which reference files are located at predefined locations
in the `exon-pipeline-files` subdirectory, scripts and programs are
located at predefined locations in the `exon-pipeline-scripts`
subdirectory, and Fastq files are located at predefined locations in the
`samples` subdirectory. In following we describe how to setup a project
specific directory structure in the samples subdirectory.

### Creating a samples directory
For each of your samples you need to create a subdirectory in `<project
directory>/samples`. For each sample directory <sample> EQP assumes that
it contains one or two Fastq files (or links to Fastq files) called
either <sample>.fastq.gz or <sample>\_1.fastq.gz and <sample>\_2.fastq.gz
depending on whether you want to process single-read or paired-end
reads. See the document `manual/EQP-manual.pdf` for more information on
how create samples directories in a more automated way using the
convenience script `create-fastq-links.sh`.

After the creation of the sample directories, the content of the
<project directory> should be:  
`LICENSE.txt  exon-pipeline-files/ samples/`  
`README.md    exon-pipeline-scripts/`    
`bin/         manual/`              


### Generating the counts 
After the creation of the sample directories EQP is ready for execution
using the script `process-samples.sh` in the subdirectory `bin`. A
typical (and complete) processing of the Fastq files in the samples
subdirectory is accomplished by the following calls:  
1. `process-samples.sh split`  
Split the input Fastq files into chunks of 25M reads and filter out
reads that originate completely from the poly-A tails. Please note that
this step needs to be executed even if the input Fastq files contain
less than 25M reads as EQP changes the reads identifiers in this step
to a format of which is required in the `compute-counts` step below.  
All sub-processes must finish before running the next command.  
2. `process-samples.sh align-all`  
Align all chunk Fastq files against the transcript, genome, and junction
Fasta reference file.  
3. `process-samples.sh compute-counts`  
Generate the gene, exon, and junction counts for each chunk. This step
can be executed directly after the alignment step without the need to
wait for the completion of the jobs.  
Again all sub-processes must finish before the next call.  
4. `process-samples.sh combine-counts`  
Combine the counts of all chunks of one sample into an aggregate count
and combine the counts for the samples in a file with has a matrix
structure with the gene, exon, or junction identifiers as the rows and
the counts for each sample as the columns.

Please read the document `manual/EQP-manual.pdf` for a much detailed
description of the different operations and various options, for
instance, to check the status of the submitted jobs.

### Location of the results
The results can then be found in the files:
`<process directory>/samples/count-files/<process directory base>-bowtie2-pe-<count type>.cnt`
where

*  `<process directory base>` is the base name of <process directory>, that is,
just the name of the directory without the complete path
*  `<count type>` is one of gene, exon, or junction.


EQP-cluster and UGE
-------------------

EQP is designed to work with the Univa® Grid Engine (UGE, formerly Sun
Grid Engine, SGE). Most of the processing steps are submitted to the UGE
scheduler via the `qsub` command (with the notable exception of
`process-samples.sh combine-counts` which is executed on the machine on
which the command is issued). Once submitted, a job is then executed
asynchronously on a node of the cluster managed by UGE. Below we show which
resources are used by the `qsub` call of EQP. It is also possible to run
EQP without UGE. In this case the option:  
`  -noqsub <number of cores>`  
needs to be supplied to all commands. With this option EQP launches
processes in the background until the number of cores given by the `-noqsub`
option is exhausted. Usually, the different commands require more than one
core:

*  `split` requires one core
*  `align-all` six cores
*  `compute-counts` four cores

Note that the number of cores used by EQP will be at
most the number given after `-noqsub`. If this number is smaller than the
minimum number of required cores for one process, then, nevertheless, at
least one process will be launched using the number of cores given by the
`-noqsub` option. 

The format of the qsub command used by EQP is as follows:  
`qsub -V -l m_mem_free=<mem> -l h_rt=<run time> -N <job name> -hold_jid <dependence jobs>`  
`  -o <log file> -pe smp <num cores> <command>`  
The parameters `-V, -N, -hold_jid, -o`, and `-pe smp` are generally
accepted by UGE `qsub` whereas the parameters `-l m_mem_free=<mem> -l
h_rt=<run time>` are dependent on the configuration of the UGE
instance. Here, `<mem>` is the memory requirement per core and `<run time>`
is either given as the number of seconds (14400 for four hours, default
value) or as hours (12:00:00 for compute-counts as it sometimes may
take more than four hours).

If your UGE instance does not accept the parameters `m_mem_free` and `h_rt`
or requires other parameters (e.g. a queue name), the recommendation is
to edit the script `process-samples.sh` and to adapt the calls to `qsub`
in the function `submitJob` which is the only place where `qsub`
jobs are actually submitted.

Testing EQP-cluster
-------------------
In order to test EQP just unpack the file `samples.tgz`:
`tar xfz samples.tgz`  
Then go to the `bin` directory, and execute run `run-test.sh`:
`cd bin`  
`run-test.sh`  
This executes the sequence of `split`, `align-all`, `compute-counts`, and
`combine-counts` for the two samples `sample1` and `sample2` (which are
contained in the `samples` directory) using the option `-noqsub`.

If you want to test whether the submission of jobs via UGE works, just
execute the four commands  
`process-samples.sh split`  
`process-samples.sh align-all`  
`process-samples.sh compute-counts`  
`process-samples.sh combine-counts`  
as described above and compare the newly computed counts in the
directory `samples/count-files` with the previously computed counts in the
directory `samples/comparison-counts` by calling (in the `samples`
subdirectory):  
`run-test.sh compare-only`


License
=======

Copyright 2015 Novartis Institutes for Biomedical Research

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
You may obtain a copy of the License at

`http://www.apache.org/licenses/LICENSE-2.0`

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed
on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for
the specific language governing permissions and limitations under the License.
