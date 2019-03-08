This tutorial re-implements in Nextflow the logic of the Snakemake workflow presented in [this webinar](https://github.com/UofABioinformaticsHub/2019_EMBL-ABR_Snakemake_webinar#data-for-the-webinar).


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



# Software environment


* Conda or either one of Docker or Singularity to run the appropriate container with all the remaining required software


  [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2468)


  ![Docker Pulls](https://img.shields.io/docker/pulls/rsuchecki/nextflow-embl-abr-webinar.svg)





# Notes

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.02.0--edge-orange.svg)](https://www.nextflow.io/)