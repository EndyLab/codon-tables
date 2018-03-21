import logging
logging.basicConfig(level=logging.INFO)

import os
import boto3

##############
# Prep stuff #
##############

# s3 variables
num_cores = 16
bucketname = 'endylab-codon-table-simulations'
s3_home_dir = 'test-simulation-4/python-jobput-test/
s3_region = 'us-west-1'

# define some path names
s3_param_dir = s3_home_dir + 'params/'
s3_output_dir = s3_home_dir + 'output/'

# DO TOUCH : JOB SPECIFIC PARAMETERS #
jobName = 'python-generated-job'
# dependsOn=[
#     {
#         'jobId' : 'name of dependent job'
#         'type' : 'N_TO_N'/'SEQUENTIAL'
#     }
# ]
# parameters = {
#     'string' : 'string'
# }
containerOverrides = {
    'vcpus' : num_cores,
    'memory' : 1000, # in MB
    'environment' : [
        {
            'name' : 'DATA_DIR', # don't change
            'value' : 'path/to/datafile/data'
        }
        {
            'name' : 'PARAM_FILE', # don't change
            'value' : 'generic_filename_with_batchnum_replaced_with{0}'
        }
        {
            'name' : 'AWS_BUCKET', # don't change
            'value' : 'endylab-codon-table-simulations' # use this one generally
        }
    ],
}
retryStrategy={
    'attempts' : 2 # or however many you want
}

## DO NOT TOUCH: JOB INVARIANT PARAMETERS ##
jobQueue = 'arn:aws:batch:us-west-1:508804160265:job-queue/codon-simulation-jobqueue'
arrayProperties = {
    'size' : num_cores
}
jobDefinition = 'arn:aws:batch:us-west-1:508804160265:job-definition/codon-sim-replicates:1'
