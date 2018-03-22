#!/bin/bash

AWS_ACCESS_KEY_ID=`cat ~/.aws/credentials | grep access_key_id  | sed 's/.*=\(.*\)/\1/'`
AWS_SECRET_ACCESS_KEY=`cat ~/.aws/credentials | grep secret_access_key  | sed 's/.*=\(.*\)/\1/'`

docker run --mount type=bind,source='/home/jonathan/Dropbox/Lab/ATD/codon-tables/data/squeeze-test/squeeze-test-4/',target=/data \
-e DATA_DIR='/data/' \
-e PARAM_FILE='2018-03-22_Standard Code_0_0_params.pickle' \
-e AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
-e AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY codon-table $@
