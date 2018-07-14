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
pbar = tqdm(filenames[1:])
local_filenames = []
for s3_filename in pbar:
    pbar.set_description('Saving {0}'.format(s3_filename))
    basepath = os.path.basename(s3_filename)
    local_filename = 'res/{0}'.format(basepath)
    # s3.download_file(bucketname, s4_filename, local_filename)
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
logging.info("Plotting 5C")

color_palette = sns.color_palette("Paired", 10, desat=0.75).as_hex()

colordict = {
    'Standard' : color_palette[1],
    'Colorado' : color_palette[5],
    'FF20' : color_palette[3],
    'FF16' : color_palette[2],
    'RED20' : color_palette[7],
    'RED15' : color_palette[6],
    'PROMISC20' : color_palette[9],
    'PROMISC15' : color_palette[8]
}
# create mini dataframes
df = DF.loc[DF['code'].map(lambda code: code not in ['Colorado', 'PROMISC20', 'RED20'])]
df_2 =  DF.loc[DF['code'].map(lambda code: code in ['PROMISC20', 'RED20'])]
# plot solid and dashed tsplots
sns.tsplot(data=df, time='N_0', value='popfrac', unit='sim', condition='code', err_style='boot_traces', n_boot=100, color=colordict)
sns.tsplot(data=df_2, time='N_0', value='popfrac', unit='sim', condition='code',
           err_style='boot_traces', n_boot=100, color=colordict, linestyle='--')

ax = plt.gca()
# plt.xlim([0.7, 1])
# plt.ylim([])
# plt.yticks([0, 0.5, 1])
# plt.xticks([])

path = '/home/jon/Lab/Fast Fail/Figures/Figure 5/'
# plt.title('Containment Probability vs Invasive Pop. Fraction'.format(code), fontsize=labelsize)
sns.despine()
plt.ylabel('Containment Probability')
plt.xlabel('Invasive Pop. Fraction')
fig = plt.gcf()
fig.set_size_inches(width, height)
plt.savefig(path+'5c.svg')
plt.show()
# format plot
logging.info("Formatting 5B Main")
sns.despine()
plt.xlim([0, 1000])
plt.xticks([0, 200, 400, 600, 800, 1000])
plt.ylim([-0.05, 1.05])
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
plt.legend()
# plt.title('Mean Fitness vs Time (1000 Simulations)', fontsize=labelsize)
plt.xlabel('Time (in generations)')
plt.ylabel('Mean Fitness')
fig = plt.gcf()
fig.set_size_inches(width, height)

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