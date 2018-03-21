import logging
logging.basicConfig(level=logging.INFO)

import os
import boto3

##############
# Prep stuff #
##############

# s3 variables
bucketname = 'endylab-codon-table-simulations'
s3_home_dir = 'test-simulation-4/python-jobput-test/
s3_region = 'us-west-1'
generalized_filename = '2018-03-19_Standard Code_0_{0}_params.pickle'

# AWS Batch variables
num_cores = 16
RAM = 2000 # in MB
attempts = 2 # per instance
jobName = 'python-generated-job'
jobDefinition = 'arn:aws:batch:us-west-1:508804160265:job-definition/codon-sim-replicates:1'
jobQueue = 'arn:aws:batch:us-west-1:508804160265:job-queue/codon-simulation-jobqueue'

# define some path names
s3_param_dir = s3_home_dir + 'params/'
s3_output_dir = s3_home_dir + 'output/'

# job specific parameters
dependsOn=[]
# dependsOn=[
#     {
#         'jobId' : 'name of dependent job'
#         'type' : 'N_TO_N'/'SEQUENTIAL'
#     },
#     {
#         'jobId' : 'name of dependent job'
#         'type' : 'N_TO_N'/'SEQUENTIAL'
#     }
# ]
parameters = {}
# parameters = {
#     'string' : 'string',
#     'string' : 'string'
# }
## DO NOT TOUCH: JOB INVARIANT PARAMETERS ##
arrayProperties = {
    'size' : num_cores
}
containerOverrides = {
    'vcpus' : num_cores,
    'memory' : RAM, # in MB
    'environment' : [
        {
            'name' : 'DATA_DIR',
            'value' : s3_home_dir
        },
        {
            'name' : 'PARAM_FILE'
            'value' : generalized_filename
        },
        {
            'name' : 'AWS_BUCKET',
            'value' : bucketname
        }
    ],
}
retryStrategy={
    'attempts' : attempts # or however many you want
}
############
# Do stuff #
############

logging.info("Submitting {0} to {1}".format(jobName, jobQueue))
batchClient = boto3.client('batch')
response = batchClient(
    jobName=jobName, jobQueue=jobQueue, arrayProperties=arrayProperties,
    dependsOn=dependsOn, jobDefinition=jobDefinition, parameters=parameters,
    containerOverrides=containerOverrides, retryStrategy=retryStrategy
)
success_string = (
    "Job Submitted! Check https://{0}.console.aws.amazon.com/batch/home?region={0}#/dashboard"
    + " to visually track simulation status."
).format(s3_region)
logging.info(success_string)
