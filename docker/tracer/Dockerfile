FROM python:3

# Metadata
LABEL container.base.image="python:3"
LABEL software.name="fit-trace"
LABEL software.version="0.1.1"
LABEL software.description="Simulates mutation-selection balance"
LABEL software.website="https://github.com/EndyLab/codon-tables"
LABEL software.documentation="https://github.com/EndyLab/codon-tables/README.md"
LABEL tags="thundeflask"

# get repository and set up any files necessary
RUN git clone "https://github.com/EndyLab/codon-tables.git"
RUN cd codon-tables && git checkout docker
RUN pip3 install -r /codon-tables/res/requirements.txt

# get environmental variable for parameters
ENV DATA_DIR=""
ENV PARAM_FILE=""

CMD ["./codon-tables/scripts/run_docker.sh"]
