import boto3
import logging

logging.basicConfig(level=logging.DEBUG)
s3 = boto3.resource('s3', region_name='us-west-1')
s3.Bucket('endylab-codon-table-simulations').download_file('test-simulation/params/3-15_StandardCode_0_0_params.pickle', './test.pickle')
