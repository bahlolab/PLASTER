FROM nfcore/base:1.14

LABEL \
  author="Jacob Munro" \
  description="Container for PLASTER pipeline" \
  maintainer="Bahlo Lab"

# set the conda env name
ARG NAME='PLASTER_v0.0'

# Install the conda environment
COPY inst/environment.yml /
RUN conda env create -f /environment.yml \
    && conda clean -a -y \
    && conda env export --name $NAME > $NAME.yml

# Install R packages from source
COPY inst/install_packages.R inst/github_packages.txt /
RUN /opt/conda/envs/$NAME/bin/Rscript --vanilla install_packages.R GITHUB:github_packages.txt

# set env variables
ENV TZ=Etc/UTC \
    R_HOME=/usr/local/lib/R/ \
    R_ENVIRON=/usr/local/lib/R/etc/Renviron R_LIBS_USER=/usr/local/lib/R/site-library \
    PATH="/opt/conda/envs/$NAME/bin:/opt/conda/bin:${PATH}"
