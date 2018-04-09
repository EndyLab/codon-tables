import logging
logging.basicConfig(level=logging.INFO)

# import necessary modules
import numpy as np
import pandas as pd
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from scipy import stats
from scipy.special import erfc
from scipy.signal import gaussian
from scipy.ndimage import convolve1d
import peakutils
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import random
from random import shuffle
import pickle
import os
import os.path
import boto3
from copy import deepcopy as copy
from src.codonTable import codonTable
from src.codonUtils import utils
from src.thunderflask import thunderflask
from src.bacteria import strain

sns.set_context("paper")
sns.set_style('white')
sns.set_style('ticks')

# set variables
bucketname = 'endylab-codon-table-simulations'
s3_path = ''
s3_region = 'us-west-1'

# get pickle files and concatenate
s3 = boto3.client('s3', region_name=s3_region)
logging.info("Getting Output Files From {0}:{1}".format(bucketname, s3_path))
filenames = [
    dict['Key'] for dict in s3.list_objects_v2(
        Bucket=bucketname,
        Prefix=s3_path,
        Delimiter='/'
    )['Contents']
]
# download files locally
logging.info("Writing Batch Output Files Locally")
pbar = tqdm(filenames[1:])
local_filenames = []
for s3_filename in pbar:
    pbar.set_description('Saving {0}'.format(s3_filename))
    basepath = os.path.basename(s3_filename)
    local_filename = 'res/{0}'.format(basepath)
    # s3.download_file(bucketname, s3_filename, local_filename)
    local_filenames.append(local_filename)
logging.info("Download Successful!")

# concatenate dataframes into one big dataframe
logging.info("Unpacking Dataframes")
dfs = []
pbar = tqdm(local_filenames)
for file in pbar:
    pbar.set_description('Unpacking {0}'.format(file))
    with open(file, 'rb') as handle:
        dfs.append(pickle.load(handle))
logging.info("Concatenating Dataframes")
DF = pd.concat(dfs, copy=False)

# extract dataframe for figure 3
codes_3b = ['Colorado', 'Standard Code']
codes_3c = list(set(DF['code']) - set(codes_3b))
f = lambda code: code in codes_3b
g = lambda code: not f(code)
DF_3b = DF.loc[DF['code'].map(f)]
DF_3c = DF.loc[DF['code'].map(g)]

logging.info("Plotting 3B: Mean Fitness Traces")
colordict = {
    'Standard Code' : 'blue',
    'Colorado' : 'red',
    'FF20' : 'green',
    'FF16' : 'orange',
    'Reductionist20' : 'purple',
    'Reductionist15' : 'brown'
}
plt.figure()
ax1 = sns.tsplot(
    data=DF_3b,
    time='time',
    value='fitness',
    unit='sim',
    condition='code',
    color=colordict,
    ci='sd'
)
# format plot
logging.info("Formatting 3B")
sns.despine()
plt.xlim([0, 1000])
plt.ylim([0, 1.3])
plt.legend()
plt.title('Mean Fitness vs Time (1000 Simulations)')
plt.xlabel('Time (in generations)')
plt.ylabel('Mean Fitness')

# save output
logging.info("Saving Figure 3B to S3")
figure_basename = '3b_vector.svg'
figure_path = '/home/ubuntu/' + figure_basename
figure_s3path = s3_path + figure_basename
plt.savefig(figure_path)
with open(figure_path, 'rb') as data:
    s3.upload_fileobj(data, bucketname, figure_s3path)
success_string = (
    "Success! Figure saved to {0}:{1}. ".format(bucketname, figure_s3path)
    + "Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'".format(s3_region)
    + " to see output figure file."
)
logging.info(success_string)

# move on to 3C:
logging.info("Plotting 3C: Mean Fitness Traces")
plt.figure()
ax1 = sns.tsplot(
    data=DF_3c,
    time='time',
    value='fitness',
    unit='sim',
    condition='code',
    color=colordict,
    ci='sd'
)
# format plot
logging.info("Formatting 3C")
sns.despine()
plt.xlim([0, 1000])
# plt.ylim([0, 1.3])
plt.legend()
plt.title('Mean Fitness vs Time (1000 Simulations)')
plt.xlabel('Time (in generations)')
plt.ylabel('Mean Fitness')

# save output
logging.info("Saving Figure 3C and inset to S3")
figure_basename = '3c_vector.svg'
figure_path = '/home/ubuntu/' + figure_basename
figure_s3path = s3_path + figure_basename
plt.savefig(figure_path)
with open(figure_path, 'rb') as data:
    s3.upload_fileobj(data, bucketname, figure_s3path)
success_string = (
    "Success! Figure saved to {0}:{1}. ".format(bucketname, figure_s3path)
    + "Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'".format(s3_region)
    + " to see output figure file."
)
logging.info(success_string)

ax1 = sns.tsplot(
    data=DF_3b.loc[DF_3b['code'] == 'Standard Code'],
    time='time',
    value='fitness',
    unit='sim',
    condition='code',
    color=colordict,
    ci='sd'
)
figure_basename = '3c_inset_vector.svg'
figure_path = '/home/ubuntu/' + figure_basename
figure_s3path = s3_path + figure_basename
plt.savefig(figure_path)
with open(figure_path, 'rb') as data:
    s3.upload_fileobj(data, bucketname, figure_s3path)
success_string = (
    "Success! Inset saved to {0}:{1}. ".format(bucketname, figure_s3path)
    + "Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'".format(s3_region)
    + " to see output figure file."
)
logging.info(success_string)
