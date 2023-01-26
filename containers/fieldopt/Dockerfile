FROM fieldopt:dev

RUN	apt-get update \
	&& apt-get install -y --no-install-recommends xvfb
RUN	pip install pyvista pythreejs

ENTRYPOINT ["/bin/bash", "-c"]
