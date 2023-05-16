# wf-tb-amr | Mycobacterium Tuberculosis Antimicrobial Resistance Identification

This repository contains a [nextflow](https://www.nextflow.io/) workflow for
the identification of variants causing anti-microbial resistance in Mycobacterium
tuberculosis targeted sequencing data.




## Introduction

`wf-tb-amr` is a workflow for determining the antibiotic resistance of
Mycobacterium tuberculosis targeted sequencing samples. The workflow handles
multiplexed sequencing runs and provides clear and simple reports summarising
the predicted resistance profile of each sample according to genetic variants
discovered.


### Workflow details


1. Align reads to NC_000962.3 reference genome (minmap2)
2. Use mpileup to determine base composition of pre-defined variants, "genotyping" (bcftools)
3. Phase variants (whatshap)
4. Report results


### Required inputs


Like all of our workflows you need to specify your input data. In this case data
must be provided in `.fastq(.gz)` format and it must be demultiplexed.

This workflow also requires a sample sheet which identifies which samples on the
run are test samples and which are controls. The sample sheet must have three
columns `barcode`, `alias`, and `type`:

* `barcode` - the barcode of the sample (i.e. barcode02).
* `alias` - the unique identifier you wish to use to refer to the sample.
* `type` - the type of sample, this can be:
  * `test_sample` - the samples for which you wish to identify antimicrobial resistance.
  * `positive_control` - a sample known to be positive for MTB.
  * `no_template_control` - a PCR blank.

The controls are assessed for performance to determine the validity of the assay.




## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[singularity]((https://docs.sylabs.io/guides/latest/user-guide/) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-tb-amr --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* a HTML document detailing QC metrics and the primary findings of the workflow.
* a CSV file containing a machine readable version of the results.
* a [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file for each sample.
* a [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) file for each sample.




## Useful links

* [WHO Variant Catalogue](https://www.who.int/publications-detail-redirect/9789240028173)
* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [singularity](https://docs.sylabs.io/guides/latest/user-guide/)
