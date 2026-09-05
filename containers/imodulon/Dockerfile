FROM python:3.12.8-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
COPY containers/imodulon/requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt
COPY bin/imodulon-analysis /usr/local/bin/imodulon-analysis
COPY bin/imodulon_analysis /usr/local/bin/imodulon_analysis
COPY bin/workflow_glue/__init__.py bin/workflow_glue/util.py bin/workflow_glue/transcript_biotypes.py /usr/local/bin/workflow_glue/
RUN imodulon-analysis --help
