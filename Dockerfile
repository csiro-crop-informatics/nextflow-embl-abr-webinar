FROM rsuchecki/miniconda3:4.5.12_d42c6c234cbabb3737a145df6a52230cf2841923

LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

RUN conda install --override-channels -c conda-forge -c bioconda -c default \
  fastqc=0.11.8 \
  multiqc=1.7 \
  trimmomatic=0.36 \
  pigz=2.3.4 \
  bwa=0.7.17 \
  samtools=1.9 \
  htslib=1.9 \
  unzip=6.0 \
  tabix=0.2.6 \
  gnu-wget=1.18