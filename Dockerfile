# Generated by: Neurodocker version 0.7.0+0.gdc97516.dirty
# Latest release: Neurodocker version 0.7.0
# Timestamp: 2020/08/14 22:31:57 UTC
# 
# Thank you for using Neurodocker. If you discover any issues
# or ways to improve this software, please submit an issue or
# pull request on our GitHub repository:
# 
#     https://github.com/ReproNim/neurodocker

FROM debian:buster

USER root

ARG DEBIAN_FRONTEND="noninteractive"

ENV LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8" \
    ND_ENTRYPOINT="/neurodocker/startup.sh"
RUN export ND_ENTRYPOINT="/neurodocker/startup.sh" \
    && apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           apt-utils \
           bzip2 \
           ca-certificates \
           curl \
           locales \
           unzip \
	   python3-pip \
	   git \
	   gcc \
	   vim \
	   nano \
	   ssh-client \
	   xvfb \
	   x11-utils \
	   ssh \
	   libx11-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
    && dpkg-reconfigure --frontend=noninteractive locales \
    && update-locale LANG="en_US.UTF-8" \
    && chmod 777 /opt && chmod a+s /opt \
    && mkdir -p /neurodocker \
    && if [ ! -f "$ND_ENTRYPOINT" ]; then \
         echo '#!/usr/bin/env bash' >> "$ND_ENTRYPOINT" \
    &&   echo 'set -e' >> "$ND_ENTRYPOINT" \
    &&   echo 'export USER="${USER:=`whoami`}"' >> "$ND_ENTRYPOINT" \
    &&   echo 'if [ -n "$1" ]; then "$@"; else /usr/bin/env bash; fi' >> "$ND_ENTRYPOINT"; \
    fi \
    && chmod -R 777 /neurodocker && chmod a+s /neurodocker
RUN test "$(getent passwd neuro)" || useradd --no-user-group --create-home --shell /bin/bash neuro
USER neuro
WORKDIR /home/neuro
ENTRYPOINT ["/neurodocker/startup.sh"]

ENV CONDA_DIR="/opt/miniconda-latest" \
    PATH="/opt/miniconda-latest/bin:$PATH"
RUN export PATH="/opt/miniconda-latest/bin:$PATH" \
    && echo "Downloading Miniconda installer ..." \
    && conda_installer="/tmp/miniconda.sh" \
    && curl -fsSL --retry 5 -o "$conda_installer" https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash "$conda_installer" -b -p /opt/miniconda-latest \
    && rm -f "$conda_installer" \
    && conda update -yq -nbase conda \
    && conda config --system --prepend channels conda-forge \
    && conda config --system --set auto_update_conda false \
    && conda config --system --set show_channel_urls true \
    && sync && conda clean -y --all && sync \
    && conda create -y -q --name neuro \
    && conda install -y -q --name neuro \
           "python=3.6" \
	   "pytest" \
           "pandas" \
           "scipy" \
           "matplotlib" \
    && sync && conda clean -y --all && sync \
    && bash -c "source activate neuro" \
    && pip3 install boto \
    && pip3 install h5py \
    && pip3 install nose \
    && pip3 install sklearn \
    && pip3 install scipy \
    && pip3 install pillow \
    && pip3 install xvfbwrapper \
    && pip3 install meshio  \
    && pip3 install nibabel \
    && pip3 install vtk \
    && pip3 install numpy \
    && pip3 install jupyterlab \
    && pip3 install notebook \
    && rm -rf ~/.cache/pip* \
    && sync

RUN bash -c 'source activate neuro'

USER root

RUN mkdir /data && chmod 777 /data && chmod a+s /data

RUN mkdir /output && chmod 777 /output && chmod a+s /output

RUN mkdir /home/neuro/repo && chmod 777 /home/neuro/repo && chmod a+s /home/neuro/repo

USER neuro 

RUN bash -c 'source activate neuro'

ARG SSH_KEY
ENV SSH_KEY=$SSH_KEY
RUN mkdir /home/neuro/.ssh/
RUN echo "$SSH_KEY" > /home/neuro/.ssh/id_rsa
RUN chmod 600 /home/neuro/.ssh/id_rsa
RUN touch /home/neuro/.ssh/known_hosts
RUN ssh-keyscan github.com >> /home/neuro/.ssh/known_hosts

CMD git clone git@github.com:hptaylor/connectome_harmonic_core.git /home/neuro/repo ;'bash'

USER root

RUN rm -rf /opt/conda/pkgs/*

USER neuro 

WORKDIR /home/neuro
