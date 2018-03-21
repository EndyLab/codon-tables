import logging
logging.basicConfig(level=logging.INFO)

from datetime import date
from src.bacteria import strain
from src.parameterizer import *
import pickle
import os

##############
# Prep stuff #
##############

# absolute path to codon-table directory
PATH='/home/jonathan/Lab/ATD/codon-tables/'

# population variables
N_pop = 1e6
sim_num = 0
batch_num = 0
strains = [strain(N_pop=N_pop)]
T_0 = 0
T_sim = 1000
dt = 0.1
t_extra = 5
N_sims = 32
mut_param = [1, 2]
date = str(date.today())
code = 'Standard Code'
filepath = 'replicate-test-2/' # this is the name for the local directory for this sim

# s3 variables
num_cores = 16
bucketname = 'endylab-codon-table-simulations'
s3_upload_dir = 'test-simulation-4/'
s3_region = 'us-west-1'

############
# Do stuff #
############

# define local and remote paths
pickle_path = PATH + 'data/' + filepath + 'params/'
s3_upload_path = s3_upload_dir + filepath + 'params/'

# generate base parameter dictionary
logging.info("Generating Base Parameter Dictionary")
params = genParamDict(
    sim_num, batch_num, strains, N_pop,
    T_0, T_sim, dt, t_extra, N_sims, mut_param,
    date, code, filepath
)

# prepare set of parameter dictionaries
logging.info("Generating Parameter Batch")
paramDicts = batcher(params, num_cores)

# inform user of the load sharing scheme
for params in paramDicts:
    batch_num = params['batch_num']
    n_sim = params['N_sims']
    print('Batch Num: {0}; n_sims = {1}'.format(batch_num, n_sim))

# store parameter files in pickle format
logging.info("Saving Parameter Files Locally to {0}".format(pickle_path))
paramPickler(paramDicts, pickle_path)

# upload files to s3
logging.info("Uploading Parameter Directory {0} to {1}:{2}".format(
    pickle_path, bucketname, s3_upload_dir
    )
)
paramUpload(pickle_path, bucketname, s3_upload_path, s3_region)
success_string = (
    "Success! Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'"
    + " to see parameter files."
).format(s3_region)
logging.info(success_string)
