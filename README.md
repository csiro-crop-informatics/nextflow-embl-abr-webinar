# Nextflow: Scalable, Sharable and Reproducible Computational Workflows across Clouds and Clusters (EMBL-ABR Nextflow webinar)

>Large analysis workflows are fragile ecosystems of software tools, scripts and dependencies. This complexity commonly makes these workflows not only irreproducible but sometimes even not re-runnable outside their original development environment. Nextflow is a reactive workflow framework and a domain specific programming language which follows the dataflow paradigm and offers an alternative, and arguably superior, approach to developing, executing and sharing pipelines.

>In this webinar we will follow the steps required for developing sharable, version controlled, container-backed workflows, which can be seamlessly executed across different environments from a laptop to cluster to cloud. We will do this by leveraging Nextflow’s integration with code and container image hosting services such as GitHub and Docker Hub, and out of the box support for various HPC cluster schedulers and the Amazon AWS cloud.

Date/time: Thursday 14 March 2019 13:00-14:00 AEDT /12:00-13:00 AEST

Presenter: Radoslaw Suchecki, CSIRO

Registration: https://attendee.gotowebinar.com/register/8408436403729692931)

# Software environment

You will need
* Java 8 or newer (check with `java -version`, should be 1.8 or newer)

  and

* [![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.02.0--edge-orange.svg)](https://www.nextflow.io/)
  ```
  curl -s https://get.nextflow.io | bash
  ```
* Conda or either one of Docker or Singularity to run the appropriate container with all the remaining required software


  [![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/2468)


  ![Docker Pulls](https://img.shields.io/docker/pulls/rsuchecki/nextflow-embl-abr-webinar.svg)