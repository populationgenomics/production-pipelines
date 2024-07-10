FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:latest

RUN pip install --no-install-compiled --upgrade -r requirements.txt

# Install metamist from the specified commit URL
RUN git clone https://github.com/populationgenomics/metamist.git \
  && cd metamist \
  && git checkout e2dd77588ff2271e970b5d43108d451afd828227 \
  && pip install .

# Copy other files as usual
COPY README.md .
COPY setup.py .
COPY cpg_workflows cpg_workflows
COPY seqr-loading-pipelines/hail_scripts hail_scripts
COPY gnomad_methods/gnomad gnomad
RUN pip install .


# e2dd77588ff2271e970b5d43108d451afd828227
# FCF 679defc3232a351160c9af3eb6504c4d3d7676bb
