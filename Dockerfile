FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:latest

# Install metamist from the specified commit URL
RUN git clone https://github.com/populationgenomics/metamist.git \
  && cd metamist \
  && python regenerate_api.py \
  && git checkout 679defc3232a351160c9af3eb6504c4d3d7676bb \
  && pip install .
# RUN pip install git+https://github.com/populationgenomics/metamist@679defc3232a351160c9af3eb6504c4d3d7676bb

# Copy other files as usual
COPY README.md .
COPY setup.py .
COPY cpg_workflows cpg_workflows
COPY seqr-loading-pipelines/hail_scripts hail_scripts
COPY gnomad_methods/gnomad gnomad
RUN pip install .
