#!/bin/bash

AWS_ACCESS_KEY_ID=`cat ~/.aws/credentials | grep access_key_id  | sed 's/.*=\(.*\)/\1/'`
AWS_SECRET_ACCESS_KEY=`cat ~/.aws/credentials | grep secret_access_key  | sed 's/.*=\(.*\)/\1/'`

docker run --mount type=bind,source=`pwd`/data,target=/data \
-e AWS_BUCKET="endylab-codon-table-simulations" \
-e DATA_DIR=test-simulation \
-e PARAM_FILE=3-15_StandardCode_0_0_params.pickle \
-e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
-e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY codon-table $@
