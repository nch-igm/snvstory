FROM continuumio/miniconda3:4.8.2

LABEL container.base.image="continuumio/miniconda3:4.8.2" \
      software.name="IGM Churchill Ancestry" \
      software.version="3.0.0"

# install system requirements so these can potentially be cached
ENV NUMBA_VERSION="0.48" PATH="/opt/conda/bin:$PATH"
RUN conda install -y \
      -c conda-forge \
      matplotlib=3.2.1 xgboost=1.0.2  \
      pandas=1.0.1 numpy=1.18.1 \
      scipy=1.4.1 umap-learn=0.4.2 \
      scikit-learn=0.24.1 \
      bokeh=2.4.3 && \
      # clean up
      conda clean --all

# install the application requirements to cache
ARG SERVICE_NAME=Ancestry
WORKDIR /opt/${SERVICE_NAME}
COPY ./requirements.txt ./
ARG CODEARTIFACT_AUTH_TOKEN
RUN pip install -r requirements.txt --extra-index-url https://aws:$CODEARTIFACT_AUTH_TOKEN@igm-257995316808.d.codeartifact.us-east-2.amazonaws.com/pypi/igm-pypi-store/simple/
RUN rm requirements.txt

# copy source code
COPY ./igm_churchill_ancestry ./igm_churchill_ancestry

ENTRYPOINT ["python3", "-m", "igm_churchill_ancestry"]

ENV TMP_DIR=/data
