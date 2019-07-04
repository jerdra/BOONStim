FROM neurodebian:xenial

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	connectome-workbench \
	&& mkdir /msm
COPY ./msm /msm
