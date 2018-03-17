#!/bin/bash

git -C codon-tables pull
pip3 install -r ./codon-tables/scripts/requirements.txt
python3 ./codon-tables/scripts/docker_test.py
