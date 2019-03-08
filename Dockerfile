FROM rsuchecki/miniconda3:4.5.12

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

LABEL maintainer="Rad Suchecki <rad.suchecki@csiro.au>"
SHELL ["/bin/bash", "-c"]

COPY conf/conda.yaml /
RUN conda env create -f /conda.yaml && conda clean -a
ENV PATH /opt/conda/envs/tutorial/bin:$PATH