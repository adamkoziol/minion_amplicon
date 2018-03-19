# Dockerfile for COWBAT
FROM ubuntu:16.04

MAINTAINER Dr. Adam G. Koziol <adam.koziol@inspection.gc.ca>

ENV DEBIAN_FRONTEND noninteractive



# Install packages
RUN apt-get update -y -qq && apt-get install -y \
	python-dev \
	git \
	curl \
	wget \
	python3-pip \
	ttf-dejavu \
	nano  

ENV PATH /usr/sbin:$PATH
RUN useradd -ms /bin/bash/ ubuntu
USER ubuntu

WORKDIR HOME
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /home/ubuntu/miniconda.sh
RUN bash /home/ubuntu/miniconda.sh -b -p /home/ubuntu/miniconda
ENV PATH /home/ubuntu/miniconda/bin:$PATH
RUN echo $PATH
	    # && rm -rf miniconda.sh \
RUN conda install -y python=3 \
	    && conda update conda	

# Add miniconda to the PATH
# ENV PATH $HOME/miniconda/bin:$PATH

# Upgrade pip
RUN pip3 install --upgrade pip

#------

# Install the pipeline
WORKDIR /home/ubuntu/
ENV PATH /home/ubuntu/minion_amplicon:$PATH
RUN git clone https://github.com/adamkoziol/minion_amplicon.git
WORKDIR /home/ubuntu/minion_amplicon
RUN conda env create


# Set the language to use utf-8 encoding - encountered issues parsing accented characters in Mash database
ENV LANG C.UTF-8