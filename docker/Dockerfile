FROM continuumio/anaconda3:latest

# RUN /opt/conda/bin/pip3 install --no-cache-dir pyscf[geomopt]
RUN apt-get install -y --no-install-recommends libxrender1 libxext6 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


COPY ./ /workspace/neoqm
WORKDIR /workspace/neoqm
RUN /opt/conda/bin/pip3 install -e .
