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
pbar = tqdm(filenames[1:7])
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
def replicate_plotter(df, code, color):
    local_df = copy(df[df['code'] == code])
    for i in tqdm(range(max(df['sim'])+1), desc='Plotting individual traces for {0}'.format(code)):
        t = np.array(local_df[local_df['sim'] == i]['time'])
        f = np.array(local_df[local_df['sim'] == i]['fitness'])
        plt.plot(t, f, color=color, alpha=0.03)
        del t, f
    del local_df

def interleave_plotter(DF, colordict):
    for i in tqdm(range(max(DF['sim'])+1), desc='Plotting individual traces (interleaved)'):
        for code in colordict.keys():
            local_DF = DF[DF['code'] == code]
            t = np.array(local_DF[local_DF['sim'] == i]['time'])
            f = np.array(local_DF[local_DF['sim'] == i]['fitness'])
            plt.plot(t, f, color=colordict[code], alpha=0.03)
            del t, f

def mean_plotter(df, code, color):
    local_df = copy(df[df['code'] == code])
    df_mean = local_df.groupby('time').mean()['fitness']
    t_mean = np.array(df_mean.index)
    f_mean = np.array(df_mean.values)
    mean_handle = plt.plot(t_mean, f_mean, color=color, alpha=1, linewidth=3, label='{0}'.format(code))
    del local_df

logging.info("Plotting Interleaved Replicates")
colordict = {
    'Standard Code' : 'blue',
    'Colorado' : 'red',
    'FF20' : 'green'
}
interleave_plotter(DF, colordict)
# logging.info("Plotting Replicates (Standard Code)")
# replicate_plotter(DF, 'Standard Code', 'gray')
# logging.info("Plotting Replicates (Colorado Code)")
# replicate_plotter(DF, 'Colorado', 'red')
# logging.info("Plotting Replicates (FF20)")
# replicate_plotter(DF, 'FF20', 'green')
logging.info("Plotting Mean (Standard Code)")
mean_plotter(DF, 'Standard Code', 'blue')
logging.info("Plotting Mean (Colorado Code)")
mean_plotter(DF, 'Colorado', 'red')
logging.info("Plotting Mean (FF20)")
mean_plotter(DF, 'FF20', 'green')

# logging.info("Plotting FF16")
# plotter(DF, 'FF16', 'gray')
# logging.info("Plotting RED20")
# plotter(DF, 'RED20', 'gray')
# logging.info("Plotting Standard Code")
# plotter(DF, 'Standard Code', 'gray')

# format plot
logging.info("Formatting Figure")
sns.despine()
plt.xlim([0, 1000])
plt.ylim([0, 1.3])

# save output
logging.info("Saving Figure to S3")
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
