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

labelsize=16
width = 4
height = width / 1.618

plt.rc('font', family='serif')
plt.rc('xtick', labelsize=labelsize)
plt.rc('ytick', labelsize=labelsize)
plt.rc('axes', labelsize=labelsize)

# set variables
bucketname = 'endylab-codon-table-simulations'
s3_path = 'manuscript/contour/'
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
        dfs.append(pickle.load(handle))
logging.info("Concatenating Dataframes")
DF = pd.concat(dfs, copy=False)

# extract dataframe for figure 3

## NOTE: Use this block if splitting evolutionary rates into two subfigures
# codes_3b = ['Colorado', 'Standard Code']
# codes_3c = list(set(DF['code']) - set(codes_3b))
# f = lambda code: code in codes_3b
# g = lambda code: not f(code)
# DF_3b = DF.loc[DF['code'].map(f)]
# DF_3c = DF.loc[DF['code'].map(g)]

# define functions for calculating containment probabilities
def contain_probability(DF, code):
    df = DF.loc[DF['code'] == code]
    N_0 = list(set(df['N_0']))
    N_0.sort()
    num_sims = len(df.loc[df['N_0'] == N_0[0]].loc[0])
    t = df.loc[df['sim'] == 0]['time']
    contain_probability = np.zeros((len(N_0), len(t)))
    for i, n_0 in enumerate(tqdm(N_0, desc='Processing Initial Conditions: ', leave=False)):
        lildf = df.loc[df['N_0'] == n_0]
        for j in tqdm(range(len(t)), desc='Processing sims: '):
            weedf = lildf.loc[j]
            contain_probability[i, j] = sum(weedf['popfrac'] == 0) / num_sims

    return contain_probability, t, N_0

def endpoint_contain(DF, code):
    df = DF.loc[DF['code'] == code]
    N_0 = list(set(df['N_0']))
    N_0.sort()
    num_sims = len(df.loc[df['N_0'] == N_0[0]].loc[0])
    timedf = df.loc[df['sim'] == 0]['time']
    endt = df.iloc[-1]
    endind = df.index[-1]
    df = df.loc[endind]
    contain = np.zeros(len(N_0))
    for i, n_0 in enumerate(N_0):
        lildf = df.loc[df['N_0'] == n_0]
        contain[i] = sum(lildf['popfrac'] == 0) / num_sims
    return contain, N_0

# calculate containment probabilities
Standard_contour, t, N_0 = contain_probability(DF, 'Standard')
Colorado_contour, t, N_0 = contain_probability(DF, 'Colorado')
FF20_contour, t, N_0 = contain_probability(DF, 'FF20')
FF16_contour, t, N_0 = contain_probability(DF, 'FF16')
RED20_contour, t, N_0 = contain_probability(DF, 'RED20')
RED15_contour, t, N_0 = contain_probability(DF, 'RED15')
PROMISC20_contour, t, N_0 = contain_probability(DF, 'PROMISC20')
PROMISC15_contour, t, N_0 = contain_probability(DF, 'PROMISC15')

# save output
file_basename = 'contour_caching_lin.pickle'
file_path = '/home/ubuntu/' + file_basename
file_s3path = s3_path + file_basename
plt.savefig(file_path)
with open(file_path, 'rb') as data:
    s3.upload_fileobj(data, bucketname, file_s3path)
success_string = (
    "Success! Inset saved to {0}:{1}. ".format(bucketname, file_s3path)
    + "Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'".format(s3_region)
    + " to see output figure file."
)
logging.info(success_string)
