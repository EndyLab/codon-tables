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

# formatting
sns.set_context("paper")
sns.set_style('white')
sns.set_style('ticks')

labelsize=16
width = 4
height = width / 1.618

plt.rc('font', family='serif')
plt.rc('xtick', labelsize=labelsize)
plt.rc('ytick', labelsize=labelsize)
plt.rc('axes', labelsize=labelsize)

# set variables
bucketname = 'endylab-codon-table-simulations'
s3_path = 'manuscript/fig5/'
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

# extract dataframe for figure 5
dashcodes = ['PROMISC20', 'PROMISC14']
f = lambda code: code not in dashcodes

logging.info("Plotting 5b-1: solid line traces")
colordict = {
    'Standard Code' : 'blue',
    'Colorado' : 'red',
    'FF20' : 'green',
    'FF16' : 'orange',
    'RED20' : 'purple',
    'RED14' : 'brown',
    'PROMISC20': 'purple',
    'PROMISC14': 'brown'
}
plt.figure()
ax1 = sns.tsplot(
    data=DF.loc[(DF['code'].map(f)) & (DF['code'] != 'Standard Code')],
    time='time',
    value='fitness',
    unit='sim',
    condition='code',
    color=colordict,
    ci='sd',
    linestyle='-'
)
logging.info("Plotting 5b-2: dashed line traces")
ax2 = sns.tsplot(
    data=DF.loc[(DF['code'].map(lambda code: not f(code)))],
    time='time',
    value='fitness',
    unit='sim',
    condition='code',
    color=colordict,
    ci='sd',
    linestyle='--'
)
# format plot
logging.info("Formatting 5B_inset")
sns.despine()
plt.xlim([0, 1000])
plt.xticks([0, 500, 1000])
plt.ylim([-0.05, 0.6])
plt.yticks([0, 0.3, 0.6])
plt.legend()
plt.title('Mean Fitness vs Time (1000 Simulations)', fontsize=labelsize)
plt.xlabel('Time (in generations)')
plt.ylabel('Mean Fitness')

# save output
logging.info("Saving Figure 5B inset to S3")
figure_basename = '5b_inset_vector.svg'
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

# move on to 5b main:
logging.info("Plotting 5B Main: Mean Fitness Traces")
ax3 = sns.tsplot(
    data=DF.loc[DF['code'] == 'Standard Code'],
    time='time',
    value='fitness',
    unit='sim',
    condition='code',
    color=colordict,
    ci='sd',
    linestyle='-'
)
# format plot
logging.info("Formatting 5B Main")
sns.despine()
plt.xlim([0, 1000])
plt.xticks([0, 200, 400, 600, 800, 1000])
plt.ylim([-0.05, 1.05])
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
plt.legend()
plt.title('Mean Fitness vs Time (1000 Simulations)', fontsize=labelsize)
plt.xlabel('Time (in generations)')
plt.ylabel('Mean Fitness')

figure_basename = '5b_vector.svg'
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
