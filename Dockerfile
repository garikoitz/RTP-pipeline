
# Create Docker container that can run Matlab (mrDiffusion and afq analysis), ANTs, FSL, mrTrix.

# Start with the Matlab r2020b 
FROM ubuntu:xenial
# install mrc dependencies
RUN apt-get -qq update && apt-get -qq install -y \ 
    unzip \
    xorg \
    wget \
    curl && \
    mkdir /mcr-install && \
    mkdir /opt/mcr && \
    cd /mcr-install && \
    wget https://ssd.mathworks.com/supportfiles/downloads/R2020b/Release/4/deployement_files/installer/complete/glnxa64/MATLAB_Runtime_R2020b_Update_4_glnxa64.zip && \
    cd /mcr-install && \
    unzip -q MATLAB_Runtime_R2020b_Update_4_glnxa64.zip && \ 
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf mcr-install

# Configure environment variables for MCR
ENV LD_LIBRARY_PATH /opt/mcr/v99/runtime/glnxa64:/opt/mcr/v99/bin/glnxa64:/opt/mcr/v99/sys/os/glnxa64:/opt/mcr/v99/extern/bin/glnxa64

# Start with the Matlab r2018b runtime container
# FROM  flywheel/matlab-mcr:v95
MAINTAINER Garkoitz Lerma-Usabiaga  <garikoitz@gmail.com>

ENV FLYWHEEL /flywheel/v0
WORKDIR ${FLYWHEEL}

# Because we're coming in from a Matlab-MCR we need to unset LD_LIBRARY_PATH
# ENV LD_LIBRARY_PATH ""

RUN apt-get update --fix-missing \
 && apt-get install -y wget bzip2 ca-certificates \
      libglib2.0-0 libxext6 libsm6 libxrender1 \
      git mercurial subversion curl grep sed dpkg gcc
RUN apt-get install -y libxt6 libxcomposite1 libfontconfig1 libasound2

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
## install mrtrix3 requirements
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

# Use tag instead of commit number, it is much more clear. 
# ENV mrtrix3COMMIT=8cef83213c4dcce7be1296849bda2b097004dd0c
# RUN curl -#L  https://github.com/MRtrix3/mrtrix3/archive/$mrtrix3COMMIT.zip | bsdtar -xf- -C /usr/lib
# WORKDIR /usr/lib/
# RUN mv mrtrix3-${mrtrix3COMMIT} mrtrix3
# RUN chmod -R +rwx /usr/lib/mrtrix3
# WORKDIR /usr/lib/mrtrix3
# RUN  ./configure && \
#     ./build && \
#     ./set_path


## install and compile mrtrix3
RUN git clone https://github.com/MRtrix3/mrtrix3.git
RUN cd mrtrix3 && git fetch --tags && git checkout tags/3.0.1 && ./configure -nogui && ./build

## manually add to path
ENV PATH=$PATH:/flywheel/v0//mrtrix3/bin


#https://wiki.ubuntu.com/DashAsBinSh
RUN rm /bin/sh && ln -s /bin/bash /bin/sh



# DTIINIT: DELETE THIS: RIGHT NOW,  just commenting

# ADD the dtiInit Matlab Stand-Alone (MSA) into the container.
# COPY dtiinit/source/bin/dtiInitStandAloneWrapper /usr/local/bin/dtiInit

# Add bet2 (FSL) to the container
# Maintain this here until we are sure we are not using it anymore
ADD https://github.com/vistalab/vistasoft/raw/97aa8a8/mrAnatomy/Segment/bet2 /usr/local/bin/
# Ensure that the executable files are executable
RUN chmod +x /usr/local/bin/bet2 
# Configure environment variables for bet2
ENV FSLOUTPUTTYPE NIFTI_GZ

# ADD AFQ and mrD templates via svn hackery
# this stopped working at some point, it is better to just use the templates we know we need and add them in the matlab compilation compilation
# Maintaining the folder just in case, maybe we can edit the code so that copies templates directly there, I guess this will be empty
ENV TEMPLATES /templates
RUN mkdir $TEMPLATES
# RUN svn export --force https://github.com/yeatmanlab/AFQ.git/trunk/templates/ $TEMPLATES
# RUN svn export --force https://github.com/vistalab/vistasoft.git/trunk/mrDiffusion/templates/ $TEMPLATES

# Add the MNI_EPI template and JSON schema files to the container
# ADD https://github.com/vistalab/vistasoft/raw/97aa8a8/mrDiffusion/templates/MNI_EPI.nii.gz $TEMPLATES
# ADD https://github.com/vistalab/vistasoft/raw/97aa8a8/mrDiffusion/dtiInit/standalone/dtiInitStandAloneJsonSchema.json $TEMPLATES

# Copy the help text to display when no args are passed in.
# COPY dtiinit/source/doc/help.txt /opt/help.txt



# Copy and configure code
WORKDIR ${FLYWHEEL}


# This is what we dont want to happen anymore
# COPY dtiinit/source/run ${FLYWHEEL}/run_dtiinit
# COPY dtiinit/source/parse_config.py ${FLYWHEEL}/dtiinit_parse_config.py

############################
# AFQ


# AND THIS IF WE ARE COPYING THE BINARY FROM LOCAL
# ADD the source Code and Binary to the container
COPY afq/source/bin/compiled/RTP /usr/local/bin/RTP
COPY afq/run ${FLYWHEEL}/run
COPY afq/source/parse_config.py ${FLYWHEEL}/parse_config.py
COPY afq/includeFiles/ply2obj.py ${FLYWHEEL}/ply2obj.py
RUN chmod +x /usr/local/bin/RTP ${FLYWHEEL}/parse_config.py

# Set the diplay env variable for xvfb
ENV DISPLAY :1.0

############################

# Configure entrypoint
RUN chmod +x ${FLYWHEEL}/*
ENTRYPOINT ["/flywheel/v0/run"]
COPY manifest.json ${FLYWHEEL}/manifest.json

RUN ldconfig
