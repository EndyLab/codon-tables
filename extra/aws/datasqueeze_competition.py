import logging
logging.basicConfig(level=logging.INFO)

import pickle
from tqdm import tqdm
import requests
import os
import os.path
import shutil
import pandas as pd
import boto3

# import environmental variables from .config_replicates.py
from config_competition import *

# define local small output file path
s3_path = S3_HOME_DIR + 'output/'
pickle_path = PATH + 'data/' + s3_path
# prepare boto3 client
s3 = boto3.client('s3', region_name=S3_REGION)
# get all files in bucket for simulation and download locally
logging.info("Getting Output Files From {0}:{1}".format(BUCKETNAME, s3_path))
filenames = [
    dict['Key'] for dict in s3.list_objects_v2(
        Bucket=BUCKETNAME,
        Prefix=s3_path,
        Delimiter='/'
    )['Contents']
]
# download files locally
logging.info("Writing Batch Output Files to {0}".format(pickle_path + 'batch/'))
os.makedirs(pickle_path + 'batch/', exist_ok=True)

pbar = tqdm(filenames)
local_filenames = []
for s3_filename in pbar:
    pbar.set_description('Saving {0}'.format(s3_filename))
    local_filename = PATH + 'data/' + s3_filename
    basepath = os.path.basename(s3_filename)
    s3.download_file(BUCKETNAME, s3_filename, local_filename)
    # move to ./batch/ subdirectory
    new_local = pickle_path + 'batch/' + basepath
    shutil.move(local_filename, new_local)
    local_filenames.append(new_local)

# loop through files and grab dataframe
logging.info("Concatenating Dataframes")
df_list = []
sim_increment = 0
pbar = tqdm(local_filenames)
for filename in pbar:
    pbar.set_description('Concatenating {0}'.format(filename))
    with open(filename, 'rb') as handle:
        df = pickle.load(handle)
    # label individual simulation runs
    n_sims = max(df['sim'])
    df['sim'] += sim_increment
    df_list.append(df)
    # increment sim_increment counter
    sim_increment += (n_sims + 1)

# concatenate dataframes and save to output file
summary_file = (
    PATH + 'data/' +
    '{0}_{1}_{2}_concatenated.pickle'.format(DATE, CODE, SIM_NUM)
)
logging.info("Writing Concatenated Data to {0}".format(summary_file))
dfs = pd.concat(df_list)
summary_file = (
    pickle_path + '{0}_{1}_{2}_concatenated.pickle'.format(DATE, CODE, SIM_NUM)
)
with open(summary_file, 'wb') as handle:
    pickle.dump(dfs, handle)
logging.info("Success! Simulation results saved to {0}".format(summary_file))
