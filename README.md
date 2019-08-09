

- [Nextflow installation](#nextflow-installation)
- [Software environment](#software-environment)
- [Nextflow version](#nextflow-version)
- [Before running the pipeline](#before-running-the-pipeline)
  - [Download example input data](#download-example-input-data)
- [Running the pipeline - execution profiles](#running-the-pipeline---execution-profiles)
  - [Local/server profiles](#localserver-profiles)
  - [HPC (SLURM) profiles](#hpc-slurm-profiles)
- [Running with pre-configured software environment](#running-with-pre-configured-software-environment)
- [Other runtime parameters](#other-runtime-parameters)

##  Nextflow installation

Check if you have Java 8 or newer (`java -version`, should be 1.8 or newer), and if you do, run

```
curl -s https://get.nextflow.io | bash
```

This will place the executable in your working directory and you should be able to run it

```
./nextflow
```

It probably makes sense to move the executable to a [directory accessible via `$PATH`](https://askubuntu.com/questions/60218/how-to-add-a-directory-to-the-path), just to be able to run `nextflow` rather than having to remember to type the full path to nextflow each time you want to run it.

Depending on the system this may suffice:

```
mkdir -p ~/bin
mv ./nextflow ~/bin
```

For consistency of how nextflow output is presented in the terminal, you should use the specific version of nextflow, e.g. `NXF_VER=19.01.0 nextflow run ...`. Alternatively, use `-ansi-log false` with later versions of nextflow.


## Software environment

Software environment required by the workflow is captured in Docker and Singularity container images we make available on Docker Hub and Singularity Hub. Alternatively you can use Conda environment as specified in [`conf/conda.yaml`](conf/conda.yaml).
You could also "simply" make sure that all required software is available in the execution environment (not recommended).


* Singularity [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2468)


* Docker  [![Docker Pulls](https://img.shields.io/docker/pulls/rsuchecki/nextflow-embl-abr-webinar.svg)](https://hub.docker.com/r/rsuchecki/nextflow-embl-abr-webinar)
  [![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/rsuchecki/nextflow-embl-abr-webinar.svg)](https://hub.docker.com/r/rsuchecki/nextflow-embl-abr-webinar)
  [![Docker Cloud build](https://img.shields.io/docker/cloud/build/rsuchecki/nextflow-embl-abr-webinar.svg)](https://hub.docker.com/r/rsuchecki/nextflow-embl-abr-webinar)


## Nextflow version

 [![Nextflow](https://img.shields.io/badge/nextflow-19.01.0-orange.svg)](https://www.nextflow.io/)



## Before running the pipeline

Although not strictly necessary for running the pipeline, it makes sense
to start by cloning this repo and moving to the directory

```
git clone https://github.com/csiro-crop-informatics/nextflow-embl-abr-webinar.git
cd nextflow-embl-abr-webinar
git checkout workshop
```

### Download example input data

```
nextflow run setup_data.nf
```

## Running the pipeline - execution profiles

The intended use is with one of the execution profiles defined in [nextflow.config](nextflow.config). Profiles define
* how we access the software environment
* what compute environment is used

### Local/server profiles

```
nextflow run main.nf -profile conda
nextflow run main.nf -profile docker
nextflow run main.nf -profile singularity
```

### HPC (SLURM) profiles

```
nextflow run main.nf -profile slurm,conda
nextflow run main.nf -profile slurm,singularity
```

## Running with pre-configured software environment

If you happen to already have all the required tools installed (available on `$PATH`), you should be able to simply run


```
nextflow run main.nf
```

or

```
nextflow run main.nf -profile slurm
```

## Other runtime parameters

* By default, only one of the 16 wheat accessions will be processed. To process a subset of say, 5 use `--take 5` or to process all use `--take all`.



