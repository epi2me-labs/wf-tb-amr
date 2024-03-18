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


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


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


