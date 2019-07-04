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
	gfortran && \
	git clone https://github.com/wujian16/Cornell-MOE.git \
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

#Set up environment for installing Cornell-MOE
ENV	MOE_CC_PATH=/usr/bin/gcc \
	MOE_CXX_PATH=/usr/bin/g++ \
	MOE_CMAKE_OPTS="-D MOE_PYTHON_INCLUDE_DIR=/usr/include/python3.7 -D MOE_PYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython3.7m.so.1.0" 
	#PYTHONPATH=$PYTHONPATH:/simnibs/Python_modules/src/:/gmsh-sdk/lib/

#Set up SimNIBS python package and compile Cornell-MOE
RUN	pip install -f https://github.com/simnibs/simnibs/releases/latest simnibs 
RUN	cd /Cornell-MOE && python ./setup.py install \
	&& pip install jupyter numba matplotlib


#Set up Field Optimization package
RUN	mkdir -p /field_pkgs/ \
	&& git clone https://github.com/jerdra/fieldopt /field_pkgs/

ENV	PYTHONPATH=$PYTHONPATH:/Cornell-MOE:/gmsh-sdk/lib/:/field_pkgs

RUN	apt-get update && apt-get install -y --no-install-recommends vim
	
ENTRYPOINT /bin/bash

