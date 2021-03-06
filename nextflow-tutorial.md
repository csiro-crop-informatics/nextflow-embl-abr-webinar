This workflow re-implements in Nextflow the logic of the Snakemake workflow presented in [this webinar](https://github.com/UofABioinformaticsHub/2019_EMBL-ABR_Snakemake_webinar#data-for-the-webinar) using the same data sets.


- [Nextflow installation](#Nextflow-installation)
- [Software environment](#Software-environment)
- [Nextflow version](#Nextflow-version)
- [Running the example pipeline](#Running-the-example-pipeline)
  - [Local/server profiles](#Localserver-profiles)
  - [HPC (SLURM) profiles](#HPC-SLURM-profiles)
  - [Cloud profile](#Cloud-profile)
- [Contributing](#Contributing)

##  Nextflow installation

You can install Nextflow using conda

```
conda install nextflow
conda update nextflow
```

Alternatively, check if you have Java 8 or newer (`java -version`, should be 1.8 or newer), and if you do, run

```
curl -s https://get.nextflow.io | bash
```

This will place the executable in your working directory and you should be able to run it

```
./nextflow
```

It probably makes sense to move the executable to a [directory accessible via `$PATH`](https://askubuntu.com/questions/60218/how-to-add-a-directory-to-the-path), just to be able to run `nextflow` rather than having to remember to type the full path to nextflow each time you want to run it.

For consistency of how nextflow output is presented in the terminal, you should use the specific version of nextflow, e.g. `NXF_VER=19.01.0 nextflow run ...`. Alternatively, use `-ansi-log false` with later versions of nextflow.


## Software environment


* Conda or either one of Docker or Singularity to run the appropriate container with all the remaining required software.


  [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2468)


  [![Docker Pulls](https://img.shields.io/docker/pulls/rsuchecki/nextflow-embl-abr-webinar.svg)](https://hub.docker.com/r/rsuchecki/nextflow-embl-abr-webinar)
  [![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/rsuchecki/nextflow-embl-abr-webinar.svg)](https://hub.docker.com/r/rsuchecki/nextflow-embl-abr-webinar)
  [![Docker Cloud build](https://img.shields.io/docker/cloud/build/rsuchecki/nextflow-embl-abr-webinar.svg)](https://hub.docker.com/r/rsuchecki/nextflow-embl-abr-webinar)


## Nextflow version

 [![Nextflow](https://img.shields.io/badge/nextflow-19.01.0-orange.svg)](https://www.nextflow.io/)


## Running the example pipeline

If you happen to already have all the required tools installed (available on `$PATH`), you should be able to simply

```
nextflow run csiro-crop-informatics/nextflow-embl-abr-webinar
```

By default, only one of the 16 wheat accessions will be processed. To process a subset of say, 5 use `--take 5` or to process all use `--take all`.


The intended use is with one of the execution [profiles](https://github.com/csiro-crop-informatics/nextflow-embl-abr-webinar/blob/13bf8b1a041541ad8bc1e30d5bdef23e2b37b67f/nextflow.config#L44-L94).

### Local/server profiles

```
nextflow run csiro-crop-informatics/nextflow-embl-abr-webinar -profile conda
nextflow run csiro-crop-informatics/nextflow-embl-abr-webinar -profile docker
nextflow run csiro-crop-informatics/nextflow-embl-abr-webinar -profile singularity
```


### HPC (SLURM) profiles

```
nextflow run csiro-crop-informatics/nextflow-embl-abr-webinar -profile slurm,conda
nextflow run csiro-crop-informatics/nextflow-embl-abr-webinar -profile slurm,singularity,singularitymodule
```


### Cloud profile

If you are new to AWS batch and/or nextflow, follow [this blog post](https://antunderwood.gitlab.io/bioinformant-blog/posts/running_nextflow_on_aws_batch/), once you are done, or you already use AWS batch, you may be able to simply run

```
nextflow run csiro-crop-informatics/nextflow-embl-abr-webinar -profile awsbatch \
  -work-dir s3://your_s3_bucket/work --outdir s3://your_s3_bucket/results
```

after replacing `your_s3_bucket` with a bucket you have created on S3.

## Contributing

* Please submit [GitHub issues](https://github.com/csiro-crop-informatics/nextflow-embl-abr-webinar/issues/new) to indicate what needs to be clarified or explained.
* Pull requests are very welcome.
