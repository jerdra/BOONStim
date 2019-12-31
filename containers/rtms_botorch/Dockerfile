FROM ubuntu:16.04

RUN apt-get update && \
	apt-get install -y --no-install-recommends \
	software-properties-common \
	gcc \
	g++ \
	cmake \
	make \
	git \
	libboost-all-dev \
	libblas-dev \
	liblapack-dev \
	libglu1 \
	libxrender-dev \
	libxcursor1 \
	libxft2 \
	libxinerama1 \
	curl \
	gfortran \
	&& add-apt-repository ppa:deadsnakes/ppa \
	&& apt-get update \
	&& apt-get install -y --no-install-recommends python3.7 python3.7-dev \
	&& curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
	&& python3.7 /get-pip.py && rm /get-pip.py \
	&& rm /usr/bin/python3 /usr/bin/python \
	&& ln -s /usr/bin/python3.7 /usr/bin/python3 \
	&& ln -s /usr/bin/python3.7 /usr/bin/python

#Get GMSH and Simnibs modules
COPY ./gmsh-sdk /gmsh-sdk

#Set up SimNIBS python package and botorch using Ax
RUN	pip install -f https://github.com/simnibs/simnibs/releases/tag/v3.0 simnibs
RUN	pip install ax-platform jupyter numba matplotlib docopt

#Set up Field Optimization package
#Force a git clone
RUN	mkdir -p /field_pkgs/
ADD	https://api.github.com/repos/jerdra/fieldopt/git/refs/heads/master version.json
RUN	git clone -b master https://github.com/jerdra/fieldopt.git /field_pkgs/fieldopt

ENV	PYTHONPATH=$PYTHONPATH:/gmsh-sdk/lib/:/field_pkgs

RUN	apt-get update && apt-get install -y --no-install-recommends vim

ENTRYPOINT /bin/bash

