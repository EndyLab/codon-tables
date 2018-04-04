FROM python:3

# Metadata
LABEL container.base.image="python:3"
LABEL software.name="fit-trace"
LABEL software.version="0.1.1"
LABEL software.description="Simulates genetic codes competing for resources"
LABEL software.website="https://github.com/EndyLab/codon-tables"
LABEL software.documentation="https://github.com/EndyLab/codon-tables/README.md"
LABEL tags="thundeflask"

# get repository and set up any files necessary
RUN apt update && apt install nano
RUN git clone "https://github.com/EndyLab/codon-tables.git"
RUN git config --global user.email "<>"
RUN git config --global user.name "dummy"
RUN git -C codon-tables pull origin docker
RUN pip3 install -r /codon-tables/res/requirements.txt

# get environmental variable for parameters
ENV DATA_DIR=""
ENV PARAM_FILE=""

CMD ["./codon-tables/scripts/run_docker_competition.sh"]
