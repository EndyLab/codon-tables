#!/bin/bash

AWS_ACCESS_KEY_ID=`cat ~/.aws/credentials | grep access_key_id  | sed 's/.*=\(.*\)/\1/'`
AWS_SECRET_ACCESS_KEY=`cat ~/.aws/credentials | grep secret_access_key  | sed 's/.*=\(.*\)/\1/'`

docker run -it endylab-codon-tables-tracer bash $@
