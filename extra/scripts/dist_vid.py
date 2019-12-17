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
from matplotlib import animation
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
wanted_codes = ['Colorado', 'Standard Code', 'FF20']
f = lambda code: code in wanted_codes
DF_3b = DF.loc[DF['code'].map(f)]

##############################################
# Sup Video: Distribution Evolving Over Time #
##############################################

logging.info("Creating Movie")
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 1.6))

# define general parameters
sims = set(DF_3b['sim'])
codes = set(DF_3b['code'])
colordict = {
    'Colorado' : 'red',
    'Standard Code' : 'blue',
    'FF20': 'green'
}

# define video parameters
fps = 30
bumper = 30
skip = 20
frames = int( (len(DF_3b.loc[(DF_3b['code'] == 'FF20') & (DF_3b['sim'] == 1)]['time']) - bumper) / skip )
dpi = 100

def framer(nFrame):
    plt.cla()
    # adjust frame with offset
    framenum = int((nFrame*skip) + bumper - 1)
    # get current fitness from simulations
    data = DF_3b.loc[framenum]

    # plot distribution
    for code in codes:
        ax = sns.distplot(data.loc[data['code'] == code]['fitness'], kde=True, hist=True, rug=False, norm_hist=True, color=colordict[code], label=code)
    plt.xlim([0,1.6])
    plt.yticks(visible=False)
    sns.despine(left=True)
    ax.axes.get_yaxis().set_visible(False)
    plt.xlabel('Mean Fitness')
    plt.ylabel('Probability')
    t = data['time'].iloc[0]
    t_before_decimal = int(t)
    t_after_decimal = t - t_before_decimal
    t_string = str(t_before_decimal) + str(t_after_decimal)[1:3]
    plt.title('Distribution of Mean Fitnesses Across Replicates (t={0})'.format(t_string))
    plt.legend()


anim = animation.FuncAnimation(fig, framer, frames=frames)
figure_basename = 'test_vid.gif'
figure_path = '/home/ubuntu/' + figure_basename
figure_s3path = s3_path + figure_basename
anim.save(figure_path, writer='imagemagick', dpi=dpi, fps=fps);
fig = plt.figure()
ax = plt.axes(xlim=(0, 1.6))

logging.info("Uploading Movie to S3")
with open(figure_path, 'rb') as data:
    s3.upload_fileobj(data, bucketname, figure_s3path)
success_string = (
    "Success! Gif saved to {0}:{1}. ".format(bucketname, figure_s3path)
    + "Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'".format(s3_region)
    + " to see output figure file."
)
logging.info(success_string)
