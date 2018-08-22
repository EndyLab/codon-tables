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
# set up colors
color_palette = sns.color_palette("Paired", 10, desat=0.75).as_hex()

colordict = {
    'Standard Code' : color_palette[1],
    'Colorado' : color_palette[5],
    'FF20' : color_palette[3],
    'FF16' : color_palette[2],
    'RED20' : color_palette[7],
    'RED15' : color_palette[9],
    'PROMISC20' : color_palette[7],
    'PROMISC15' : color_palette[9]
}
# create mini dataframes

labelsize=16
width = 8 / 1.5
height = 6 / 1.5 #width / 1.618

plt.rc('font', family='serif')
plt.rc('xtick', labelsize=labelsize)
plt.rc('ytick', labelsize=labelsize)
plt.rc('axes', labelsize=labelsize)

# set variables
bucketname = 'endylab-codon-table-simulations'
s3_path = 'manuscript/fig5c/'
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
pbar = tqdm(filenames[1:-1])
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
N_0 = list(set(DF.loc[DF['code'] == 'Standard']['N_0']))
N_0.sort()
num_reps = len(DF.loc[(DF['N_0']== N_0[0]) & (DF['code'] == 'Standard')])
codes = [code for code in colordict.keys() if code not in ['FF20', 'FF16', 'Colorado', 'FFQUAD']]
for code in tqdm(codes, desc='codes'):
    logging.info("Processing Data for {0}".format(code))
    for n_0 in tqdm(N_0, desc='initial conditions'):
        DF.loc[(DF['code'] == code)&(DF['N_0'] == n_0), 'sim'] = np.arange(num_reps)

DF.loc[:,'popfrac'] = (DF.loc[:,'popfrac'] == 0)
DF.loc[:,'N_0'] /= 1e6
# extract dataframe for figure 5
logging.info("Plotting 5C")

df = DF.loc[DF['code'].map(lambda code: code not in ['Colorado', 'PROMISC20', 'PROMISC15'])]
df_2 =  DF.loc[DF['code'].map(lambda code: code in ['PROMISC20', 'PROMISC15'])]
# plot solid and dashed tsplots
sns.tsplot(data=df, time='N_0', value='popfrac', unit='sim', condition='code', err_style='boot_traces', n_boot=100, color=colordict)
sns.tsplot(data=df_2, time='N_0', value='popfrac', unit='sim', condition='code',
           err_style='boot_traces', n_boot=100, color=colordict, linestyle='--')

# format plot
logging.info("Formatting 5C Main")
sns.despine()
plt.ylabel('Containment Probability')
plt.xlabel('Invasive Pop. Fraction')
fig = plt.gcf()
fig.set_size_inches(width, height)
plt.legend()
plt.show()

figure_basename = '5c_vector.svg'
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
