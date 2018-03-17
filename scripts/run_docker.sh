#!/bin/bash

git -C codon-tables pull
python3 ./codon-tables/scripts/docker_test.py
