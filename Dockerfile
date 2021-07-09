FROM nfcore/base:1.14

LABEL \
  author="Jacob Munro" \
  description="Container for PLASTER pipeline" \
  maintainer="Bahlo Lab"

# set the conda env name
ARG NAME='PLASTER_v0.0'

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml \
    && conda clean -a -y \
    && conda env export --name $NAME > $NAME.yml

# set path
ENV PATH="/opt/conda/envs/$NAME/bin:/opt/conda/bin:${PATH}"
