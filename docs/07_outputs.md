Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | ./wf-tb-amr-report.html | The report for all samples in the workflow run | aggregated |
| Workflow CSV results summary | ./wf-tb-amr-report.csv | The CSV summary of the results of the workflow | aggregated |
| Alignment | ./{{ alias }}.bam | Aligned reads for the sample in BAM format | per-sample |
| Alignment index | ./{{ alias }}.bam.bai | An index file for the alignment in BAI format | per-sample |
| Variants | ./{{ alias }}.final.vcf | Called, annotated variants for the sample. | per-sample |
| Per sample report | ./{{ alias }}_report.html | The report for a single sample | per-sample |
