FROM fieldopt:dev

RUN	apt-get update \
	&& apt-get install -y --no-install-recommends xvfb
RUN	pip install pyvista

ENTRYPOINT ["/bin/bash", "-c"]
