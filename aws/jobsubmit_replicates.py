import logging
logging.basicConfig(level=logging.INFO)

import os
import boto3

# import environmental variables from .config_replicates.py
from config_replicates import *

# define some path names
S3_PARAM_DIR = S3_HOME_DIR + 'params/'
S3_OUTPUT_DIR = S3_HOME_DIR + 'output/'

## DO NOT TOUCH: JOB INVARIANT PARAMETERS ##
arrayProperties = {
    'size' : NUM_CORES
}
containerOverrides = {
    'vcpus' : NUM_CORES,
    'memory' : RAM, # in MB
    'environment' : [
        {
            'name' : 'DATA_DIR',
            'value' : S3_HOME_DIR
        },
        {
            'name' : 'PARAM_FILE',
            'value' : GENERALIZED_FILENAME
        },
        {
            'name' : 'AWS_BUCKET',
            'value' : BUCKETNAME
        }
    ],
}
retryStrategy={
    'attempts' : ATTEMPTS # or however many you want
}
############
# Do stuff #
############

logging.info("Submitting {0} to {1}".format(JOBNAME, JOBQUEUE))
batchClient = boto3.client('batch')
response = batchClient.submit_job(
    jobName=JOBNAME, jobQueue=JOBQUEUE, arrayProperties=ARRAYPROPERTIES,
    dependsOn=DEPENDS_ON, jobDefinition=JOBDEFINITION, parameters=PARAMETERS,
    containerOverrides=CONTAINEROVERRIDES, retryStrategy=RETRYSTRATEGY
)
success_string = (
    "Job Submitted! Check https://{0}.console.aws.amazon.com/batch/home?region={0}#/dashboard"
    + " to visually track simulation status."
).format(S3_REGION)
logging.info(success_string)
