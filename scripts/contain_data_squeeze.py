import logging
logging.basicConfig(level=logging.INFO)

# import necessary modules
import numpy as np
import pandas as pd
from tqdm import tqdm
import pickle
import os
import os.path
import boto3
from copy import deepcopy as copy

# set S3 variables
bucketname = 'endylab-codon-table-simulations'
s3_path = 'manuscript/data/contours/'
s3_region = 'us-west-1'
s3 = boto3.client('s3', region_name=s3_region)
filenames = [
    dict['Key'] for dict in s3.list_objects_v2(
        Bucket=bucketname,
        Prefix=s3_path,
        Delimiter='/'
    )['Contents']
]
# download competition simulations locally
logging.info("Writing Batch Output Files Locally")
pbar = tqdm(filenames[1:])
local_filenames = []
for s3_filename in pbar:
    pbar.set_description('Saving {0}'.format(s3_filename))
    basepath = os.path.basename(s3_filename)
    local_filename = 'res/{0}'.format(basepath)
    s3.download_file(bucketname, s3_filename, local_filename)
    local_filenames.append(local_filename)
logging.info("Download Successful!")

# concatenate dataframes into one big dataframe
logging.info("Unpacking Dataframes")
dfs = []
pbar = tqdm(local_filenames)
for file in pbar:
    pbar.set_description('Unpacking {0}'.format(file))
    with open(file, 'rb') as handle:
        dfs.append(pickle.load(handle).loc[1000])
logging.info("Concatenating Dataframes")
DF = pd.concat(dfs, copy=False)


# massage dataframes into proper format
N_0 = list(set(DF.loc[DF['code'] == 'RED20']['N_0']))
N_0.sort()
print(N_0)
num_reps = len(DF.loc[(DF['N_0']== N_0[0]) & (DF['code'] == 'RED20')])
codes = [code for code in colordict.keys()]
for code in tqdm(codes, desc='codes'):
    logging.info("Processing Data for {0}".format(code))
    for n_0 in tqdm(N_0, desc='initial conditions'):
        DF.loc[(DF['code'] == code)&(DF['N_0'] == n_0), 'sim'] = np.arange(num_reps)

DF.loc[:,'popfrac'] = (DF.loc[:,'popfrac'] == 0)
DF.loc[:,'N_0'] /= 1e6

# write compressed DF to file
logging.info("Writing Concatenated DataFrame to File")
s3_outpath = 'manuscript/data/'
base_filename = 'lin_contour_data.pickle'
pickle_path = '../res/'
# pickle first
with open(pickle_path+base_filename, 'wb') as handle:
    pickle.dump(DF, handle)
# write to s3
with open(pickle_path+base_filename, 'rb') as data:
    s3.upload_fileobj(data, bucketname, s3_outpath+base_filename)
success_string = (
    "Success! Inset saved to {0}:{1}. ".format(bucketname, s3_outpath+base_filename)
    + "Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'".format(s3_region)
    + " to see output figure file."
)
logging.info(success_string)
