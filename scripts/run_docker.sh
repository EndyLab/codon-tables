#!/bin/bash

echo "Updating from git"
git -C codon-tables pull

echo "Updating requirements"
pip3 install -r ./codon-tables/scripts/requirements.txt

echo "Running!"
python3 ./codon-tables/scripts/docker_test.py
