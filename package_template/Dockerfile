# from continuumio/miniconda3
from quay.io/nyuad_cgsb/easybuild

##When run on travis run this from 1 level up of the pkg cache
##  - "/opt/anaconda3/pkgs"
##So it is run from /opt/anaconda3/

ARG CONDA_PKG=local
ENV CONDA_PKG ${CONDA_PKG}

RUN echo ${CONDA_PKG}
RUN which conda
RUN conda info
RUN mkdir -p /opt/conda/conda-bld/linux-64

## We have to add the pkgs and the build-system
COPY pkgs/* /opt/conda/pkgs/

RUN conda config --set always_yes yes --set changeps1 no
RUN conda config --add channels r
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda config --add channels nyuad-cgsb

RUN ls -lah /opt/conda/conda-bld/linux-64
RUN conda config --add channels file://opt/conda/conda-bld/linux-64
RUN conda index /opt/conda/conda-bld/linux-64
RUN conda install -y --use-local gencore_rnaseq
