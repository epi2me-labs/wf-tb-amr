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