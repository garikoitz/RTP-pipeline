
# Create Docker container that can run Matlab (mrDiffusion and afq analysis), ANTs, FSL, mrTrix.

# Start with the Matlab r2017a runtime container
FROM  flywheel/matlab-mcr:v92.1
MAINTAINER Michael Perry <lmperry@stanford.edu>

ENV FLYWHEEL /flywheel/v0
WORKDIR ${FLYWHEEL}

# Because we're coming in from a Matlab-MCR we need to unset LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH ""

###########################
# Configure neurodebian (https://github.com/neurodebian/dockerfiles/blob/master/dockerfiles/xenial-non-free/Dockerfile)
RUN set -x \
	&& apt-get update \
	&& { \
		which gpg \
		|| apt-get install -y --no-install-recommends gnupg \
	; } \
	&& { \
		gpg --version | grep -q '^gpg (GnuPG) 1\.' \
		|| apt-get install -y --no-install-recommends dirmngr \
	; } \
	&& rm -rf /var/lib/apt/lists/*

RUN set -x \
	&& export GNUPGHOME="$(mktemp -d)" \
	&& gpg --keyserver ha.pool.sks-keyservers.net --recv-keys DD95CC430502E37EF840ACEEA5D32F012649A5A9 \
	&& gpg --export DD95CC430502E37EF840ACEEA5D32F012649A5A9 > /etc/apt/trusted.gpg.d/neurodebian.gpg \
	&& rm -rf "$GNUPGHOME" \
	&& apt-key list | grep neurodebian

RUN { \
	echo 'deb http://neuro.debian.net/debian xenial main'; \
	echo 'deb http://neuro.debian.net/debian data main'; \
	echo '#deb-src http://neuro.debian.net/debian-devel xenial main'; \
} > /etc/apt/sources.list.d/neurodebian.sources.list

RUN sed -i -e 's,main *$,main contrib non-free,g' /etc/apt/sources.list.d/neurodebian.sources.list; grep -q 'deb .* multiverse$' /etc/apt/sources.list || sed -i -e 's,universe *$,universe multiverse,g' /etc/apt/sources.list


############################
# Install dependencies
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --force-yes \
    xvfb \
    xfonts-100dpi \
    xfonts-75dpi \
    xfonts-cyrillic \
    zip \
    unzip \
    python \
    imagemagick \
    wget \
    subversion \
    fsl-5.0-core \
    fsl-first-data \
    jq \
    ants


############################
# MRTRIX 3

# Here we download and build MRTRIX3 from source.
RUN apt-get install -y \
    git \
    g++ \
    bsdtar \
    python \
    python-numpy \
    libeigen3-dev \
    zlib1g-dev \
    libqt4-opengl-dev \
    libgl1-mesa-dev \
    libfftw3-dev \
    libtiff5-dev

ENV mrtrix3COMMIT=8cef83213c4dcce7be1296849bda2b097004dd0c
RUN curl -#L  https://github.com/MRtrix3/mrtrix3/archive/$mrtrix3COMMIT.zip | bsdtar -xf- -C /usr/lib
WORKDIR /usr/lib/
RUN mv mrtrix3-${mrtrix3COMMIT} mrtrix3
RUN chmod -R +rwx /usr/lib/mrtrix3
WORKDIR /usr/lib/mrtrix3
RUN  ./configure && \
    ./build && \
    ./set_path

ENV PATH /usr/lib/mrtrix3/bin:$PATH


############################
# DTIINIT

# ADD the dtiInit Matlab Stand-Alone (MSA) into the container.
COPY dtiinit/source/bin/dtiInitStandAloneWrapper /usr/local/bin/dtiInit

# Add bet2 (FSL) to the container
ADD https://github.com/vistalab/vistasoft/raw/97aa8a8/mrAnatomy/Segment/bet2 /usr/local/bin/

# ADD AFQ and mrD templates via svn hackery
ENV TEMPLATES /templates
RUN mkdir $TEMPLATES
RUN svn export --force https://github.com/yeatmanlab/AFQ.git/trunk/templates/ $TEMPLATES
RUN svn export --force https://github.com/vistalab/vistasoft.git/trunk/mrDiffusion/templates/ $TEMPLATES

# Add the MNI_EPI template and JSON schema files to the container
ADD https://github.com/vistalab/vistasoft/raw/97aa8a8/mrDiffusion/templates/MNI_EPI.nii.gz $TEMPLATES
ADD https://github.com/vistalab/vistasoft/raw/97aa8a8/mrDiffusion/dtiInit/standalone/dtiInitStandAloneJsonSchema.json $TEMPLATES

# Copy the help text to display when no args are passed in.
COPY dtiinit/source/doc/help.txt /opt/help.txt

# Ensure that the executable files are executable
RUN chmod +x /usr/local/bin/bet2 /usr/local/bin/dtiInit

# Configure environment variables for bet2
ENV FSLOUTPUTTYPE NIFTI_GZ

# Copy and configure code
WORKDIR ${FLYWHEEL}
COPY dtiinit/source/run ${FLYWHEEL}/run_dtiinit
COPY dtiinit/source/parse_config.py ${FLYWHEEL}/dtiinit_parse_config.py


############################
# AFQ Browser
# It built the docker container but then it failed in the execution. Just removing it, never used. 
# RUN apt-get update -qq \
#     && apt-get install -y \
#     git \
#     python-dev \
#     libblas-dev \
#     liblapack-dev \
#     libatlas-base-dev \
#     gfortran \
#     python-numpy \
#     python-pandas \
#     python-scipy \
#     python-pip

# We need to start by upgrading setuptools, or run into https://github.com/yeatmanlab/AFQ-Browser/issues/101
# RUN pip install --upgrade setuptools

# Bust the cache to force the next steps:
# ENV BUSTCACHE 11

# Install AFQ-Browser from my branch:
# RUN pip install git+https://github.com/arokem/AFQ-Browser.git@zip

# Copy AFQ Browser run
# COPY afq-browser/run ${FLYWHEEL}/run_afq-browser.py


############################
# AFQ


# THIS IS IF WE ARE USING THE GIT LFS VERSION
# Install git-lfs
# RUN apt-get install -y software-properties-common && \
#     add-apt-repository -y ppa:git-core/ppa && \
#     curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash && \
#     apt-get install -y git-lfs && \
#     git lfs install

# Download binary from git-lfs
# WORKDIR /tmp
# ENV COMMIT 9a8c9c9
# RUN git clone https://github.com/scitran-apps/afq-pipeline.git && \
#     cd /tmp/afq-pipeline && \
#     git reset ${COMMIT} && \
#     git lfs pull && \
#     cp /tmp/afq-pipeline/afq/source/bin/compiled/AFQ_StandAlone_QMR /usr/local/bin/AFQ


# AND THIS IF WE ARE COPYING THE BINARY FROM LOCAL
# ADD the source Code and Binary to the container
COPY afq/source/bin/compiled/RTP /usr/local/bin/AFQ
COPY afq/run ${FLYWHEEL}/run_afq
COPY afq/source/parse_config.py ${FLYWHEEL}/afq_parse_config.py
RUN chmod +x /usr/local/bin/AFQ ${FLYWHEEL}/afq_parse_config.py


# AND THIS WAS COMMON
# ADD the source Code and Binary to the container
COPY afq/run ${FLYWHEEL}/run_afq
COPY afq/source/parse_config.py ${FLYWHEEL}/afq_parse_config.py
RUN chmod +x /usr/local/bin/AFQ ${FLYWHEEL}/afq_parse_config.py


# Set the diplay env variable for xvfb
ENV DISPLAY :1.0

############################

# Configure entrypoint
COPY run ${FLYWHEEL}/run
RUN chmod +x ${FLYWHEEL}/*
ENTRYPOINT ["/flywheel/v0/run"]
COPY manifest.json ${FLYWHEEL}/manifest.json

