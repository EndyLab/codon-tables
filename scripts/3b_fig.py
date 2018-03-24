import logging
logging.basicConfig(level=logging.INFO)

# import necessary modules
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from scipy import stats
from scipy.special import erfc
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

# set variables
bucketname = 'endylab-codon-table-simulations'
s3_path = 'manuscript/fig3b/'
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

# begin plotting
def plotter(df, code, color):
    local_df = copy(df[df['code'] == code])
    for i in tqdm(range(max(df['sim'])+1), desc='Plotting individual traces for {0}'.format(code)):
        t = np.array(local_df[local_df['sim'] == i]['time'])
        f = np.array(local_df[local_df['sim'] == i]['fitness'])
        plt.plot(t, f, color=color, alpha=0.1)
        del t, f
    df_mean = local_df.groupby('time').mean()['fitness']
    t_mean = np.array(local_df_mean.index)
    f_mean = np.array(local_df_mean.values)
    mean_handle = plt.plot(t_mean, f_mean, color=color, alpha=0.7, label='{0} (mean)'.format(code))

logging.info("Plotting Standard Code")
plotter(DF, 'Standard Code', 'gray')
logging.info("Plotting Colorado Code")
plotter(DF, 'Colorado', 'red')
logging.info("Plotting FF20")
plotter(DF, 'FF20', 'green')
# logging.info("Plotting FF16")
# plotter(DF, 'FF16', 'gray')
# logging.info("Plotting RED20")
# plotter(DF, 'RED20', 'gray')
# logging.info("Plotting Standard Code")
# plotter(DF, 'Standard Code', 'gray')

# format plot
sns.set_style('white')
sns.set_style('ticks')
sns.despine()
sns.plt.ylim(0,1.3)
