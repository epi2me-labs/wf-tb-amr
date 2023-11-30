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