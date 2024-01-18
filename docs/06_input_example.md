<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts FASTQ files as input.

The FASTQ input parameter for this workflow accepts the path to a directory containing one level of sub-directories which in turn contain FASTQ files. The data is assumed to be multiplexed with the names of the sub-directories as barcodes. A sample sheet can be provided with `--sample_sheet`.

```   
input_directory
├── barcode01
│   ├── reads0.fastq
│   └── reads1.fastq
├── barcode02
│   ├── reads0.fastq
│   ├── reads1.fastq
│   └── reads2.fastq
└── barcode03
└── reads0.fastq
```