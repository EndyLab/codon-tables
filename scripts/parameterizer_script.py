from datetime import date
from src.bacteria import strain
from src.parameterizer import *
import pickle
import os
import boto3

##############
# Prep stuff #
##############

# population variables
N_pop = 3e5
sim_num = 0
batch_num = 0
strains = [strain(N_pop=N_pop)]
T_0 = 0
T_sim = 1000
dt = 0.1
t_extra = 5
N_sims = 100
mut_param = [1, 2]
date = str(date.today())
code = 'Standard Code'
filepath = 'upload-test/' # this is the name for the local directory for this sim

# s3 variables
num_cores = 32
bucketname = 'endylab-codon-table-simulations'
s3_upload_dir = 'test-upload/'
s3_region = 'us-west-1'

# local pickle options
pickle_path = '../data/' + filepath 'params/'

############
# Do stuff #
############

# generate base parameter dictionary
params = genParamDict(
    sim_num, batch_num, strains, N_pop,
    T_0, T_sim, dt, t_extra, N_sims, mut_param,
    date, code, filepath
)

# prepare set of parameter dictionaries
paramDicts = batcher(params, num_cores)
# inform user of the load sharing scheme
for params in paramDicts:
    batch_num = params['batch_num']
    n_sim = params['N_sims']
    print('Batch Num: {0}; n_sims = {1}'.format(batch_num, n_sim))

# pickle_file = '{0}_{1}_{2}_{3}_params.pickle'.format(date, code, sim_num, batch_num)

# output parameter dictionaries to pickle files
paramPickler(paramDicts, pickle_path, pickle_file)

# prepare boto3 client
s3 = boto3.resource('s3', region_name=s3_region)

# recursively loop through files to upload
for root, dirs, files in os.walk(pickle_path):
    for filename in files:
        # get local path
        local_path = os.path.join(root, filename)
        # massage filename
        outfile = s3_upload_dir + 'data/' + filepath + 'params/' + pickle_file
