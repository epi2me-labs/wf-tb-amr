# Mycobacterium tuberculosis workflow

Mycobacterium tuberculosis workflow for multiplexed Nanopore sequencing data.



## Introduction

<!---This section of documentation typically contains a list of things the workflow can perform also any other intro.--->

`wf-tb-amr` is a workflow for determining the antibiotic resistance of
Mycobacterium tuberculosis targeted sequencing samples. The workflow handles
multiplexed sequencing runs and provides clear and simple reports summarising
the predicted resistance profile of each sample according to genetic variants
discovered.



## Compute requirements

Recommended requirements:

+ CPUs = 16
+ Memory = 32GB

Minimum requirements:

+ CPUs = 8
+ Memory = 16GB

Approximate run time: 5 minutes per sample

ARM processor support: True




## Install and run

<!---Nextflow text remains the same across workflows, update example cmd and demo data sections.--->

These are instructions to install and run the workflow on command line. You can also access the workflow via the [EPI2ME application](https://labs.epi2me.io/downloads/).

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and software resources, therefore nextflow will need to be installed before attempting to run the workflow. 

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or
[singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of 
the required software. Both methods are automated out-of-the-box provided 
either docker or singularity is installed. This is controlled by the [`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles) parameter as exemplified below. 

It is not required to clone or download the git repository in order to run the workflow. 
More information on running EPI2ME workflows can be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow. This will pull the repository in to the assets folder of nextflow and provide a list of all parameters available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-tb-amr --help
```
A demo dataset is provided for testing of the workflow. It can be downloaded using: 
```
wget <https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-tb-amr/wf-tb-amr-demo.tar.gz> 
tar -xzvf wf-tb-amr-demo.tar.gz 
```
The workflow can be run with the demo data using: 
```
nextflow run epi2me-labs/wf-tb-amr \ 
--fastq wf-tb-amr-demo/fastq \
--sample_sheet wf-tb-amr-demo/sample_sheet.csv \
-profile standard 
```
For further information about running a workflow on the cmd line see https://labs.epi2me.io/wfquickstart/  



## Related protocols

<!---Hyperlinks to any related protocols that are directly related to this workflow, check the community for any such protocols.--->

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Inputs

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for this workflow and should have one of the following values; `test_sample`, `positive_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Reference Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| reference | string | NCBI accession for reference genome. | By default the workflow uses NC_000962.3. WARNING: If you change this parameter but don't alter the variant database, Genbank file, and the amplicon BED (all generated for NC_000962.3) to match, then the behaviour of the workflow is unlikely to be as expected. |  |
| amplicons_bed | string | The location of the amplicons for the assay. | A BED file describing the location of the amplicons used to generate the data to be processed by this workflow. Based on NC_000962.3. |  |
| variant_db | string | Variant database in VCF format. | A list of variants, in VCF format, that this assay will genotype. SNPs only. |  |
| genbank | string | Genbank file for organism of interest. | Genbank file used for variant annotation, defaults to NC000962.3 annotations. |  |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| align_threads | number | Number of CPU threads to use per alignment task. | The total CPU resource used by the workflow is constrained by the executor configuration. | 4 |
| mpileup_threads | number | Number of CPU threads to use per mpileup task. | The total CPU resource used by the workflow is constrained by the executor configuration. | 4 |
| maf | number | Minimum mutant allele frequency to consider. | By default the workflow will filter any variant which is present at less than 10% allele frequency (0.1). Change this parameter to alter this filtering behaviour. Minimum is set at 1%. | 0.1 |
| downsample | integer | Number of reads to downsample to in each direction, leave empty for no downsampling. | Downsampling can help with run times without significantly impacting the result of your analysis. By default no downsampling is performed, but you can specify the coverage to downsample to in each direction here, to a minimum of 100. |  |
| minimum_read_support | integer | The minimum number of reads to consider for a variant call on each strand. | By default the workflow expects a minimum of 5 reads on each strand supporting a variant. This is to ensure that when using the `maf`, very low coverage regions do not contribute to potentially false positive variant calls. | 5 |
| ntc_threshold | string | Comma separated string of x,y - where x is the read count threshold and y is the number of amplicons i.e. 20,3 - fail is more than 20 reads in more than 3 amplicons. |  | 20,3 |
| sample_threshold | string | Comma separated string of x,y - where x is the read count threshold and y is the number of amplicons. For example 20,8 means a sample will pass only if at least 20 reads were detected in each of at least 8 amplicons. |  | -20,8 |
| positive_threshold | string | Comma separated string of x,y - where x is the read count threshold and y is the number of amplicons i.e. 20,2 - fail is less than 20 reads in less than 2 amplicons. |  | -20,2 |
| strand_bias | integer | Set a threshold for strand bias filtering. | Strand bias is represented as a Phred scaled p-value from a Fisher's exact test, with a value close to 0 being preferable. | 1000 |
| report_config | string | Report configuration file. | The report can be configured to help with translation. See `report_config.eng.json` in the primer scheme directory. Here you can provide a path to your own report configuration file. |  |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| disable_ping | boolean | Enable to prevent sending a workflow ping. |  | False |






## Outputs

Outputs files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | ./wf-tb-amr-report.html | The report for all samples in the workflow run | aggregated |
| Workflow CSV results summary | ./wf-tb-amr-report.csv | The CSV summary of the results of the workflow | aggregated |
| Alignment | ./{{ alias }}.bam | Aligned reads for the sample in BAM format | per-sample |
| Alignment index | ./{{ alias }}.bam.bai | An index file for the alignment in BAI format | per-sample |
| Variants | ./{{ alias }}.final.vcf | Called, annotated variants for the sample. | per-sample |
| Per sample report | ./{{ alias }}_report.html | The report for a single sample | per-sample |




## Pipeline overview

<!---High level numbered list of main steps of the workflow and hyperlink to any tools used. If multiple workflows/different modes perhaps have subheadings and numbered steps. Use nested numbering or bullets where required.--->
### 1. Concatenates input files and generate per read stats.

The [fastcat/bamstats](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow. It will also output per read stats including average read lengths and qualities.

### 2. Align reads to NC_000962.3 reference genome

[minimap2](https://github.com/lh3/minimap2) is used to align reads from the samples to the Mycobacterium tuberculosis reference genome FASTA (NC_000962.3). This step also discards unmapped reads and generates statistics from the resulting BAM file.

### 3. Downsample reads

Downsampling is off by default. This optional step downsamples the data to a specified number of reads.

### 4. Call variants using mpileup

The [bcftools](https://github.com/samtools/bcftools) mpileup tool is used to determine base composition of pre-defined variants.

### 5. Phase variants

The [whatshap](https://github.com/whatshap/whatshap) tool is used to phase variants. Codon numbers are then added using [vcf-annotator](https://github.com/rpetit3/vcf-annotator), and the results processed to ensure that phased variants are correctly annotated.

### 6. Report results

The workflow outputs an HTML report with overall results for all samples in the run, indivdual sample HTML reports, and a summary CSV file.



## Troubleshooting

<!---Any additional tips.--->
+ If the workflow fails please run it with the demo dataset to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).



## FAQ's

<!---Frequently asked questions, pose any known limitations as FAQ's.--->

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-tb-amr/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).



## Related blog posts

<!---Any other sections that are relevant specifically to this workflow and may be useful to users eg. ## Related blog posts. ## Learning center links.--->
+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



