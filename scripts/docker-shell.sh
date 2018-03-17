#!/bin/bash

echo $AWS_ACCESS_KEY_ID

docker run -it --mount type=bind,source=`pwd`/data,target=/data \
-e AWS_BUCKET="endylab-codon-table-simulations" \
-e DATA_DIR=test-simulation \
-e PARAM_FILE=3-15_StandardCode_0_0_params.pickle \
-e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
-e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY codon-table bash $@
