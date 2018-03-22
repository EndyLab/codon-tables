import logging
logging.basicConfig(level=logging.INFO)

import pickle
from tqdm import tqdm
import requests
import os
import os.path
import boto3

# get environmental varaiables
datapath = os.environ['DATA_DIR'] + 'output/'
awsbucket = os.environ['AWS_BUCKET'] if 'AWS_BUCKET' in os.environ else ''

## DRAFTING: WORK WITH LOCAL FILES FOR NOW ##
datapath = ~/Lab/ATD/codon-tables/data/squeeze-test-1/output/
# prepare boto3 client
s3 = boto3.resource('s3', region_name='us-west-1')
## recursively walk through files to upload
# for root, dirs, files in os.walk(datapath):
    # for filename in files:
        ## get local and remote paths
        # local_path = os.path.join(root, filename)
        # infile = s3_upload_path + filename
        ## write to s3
        # with open(local_path, 'rb') as handle:
            # s3.Bucket(bucket).put_object(Key=outfile, Body=handle)
