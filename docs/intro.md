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
