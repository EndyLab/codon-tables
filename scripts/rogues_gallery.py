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
s3_path = 'BIOE_retreat/rogues_gallery/'
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

## NOTE: Use this block if splitting evolutionary rates into two subfigures
# codes_3b = ['Colorado', 'Standard Code']
# codes_3c = list(set(DF['code']) - set(codes_3b))
# f = lambda code: code in codes_3b
# g = lambda code: not f(code)
# DF_3b = DF.loc[DF['code'].map(f)]
# DF_3c = DF.loc[DF['code'].map(g)]

DF_3b = DF
logging.info("Plotting 3B: Mean Fitness Traces")
colors = sns.palplot(sns.cubehelix_palette(6, start=0, rot=1, reverse=False)).as_hex()
colordict = {
    'Standard Code' : colors[0],
    'FF20' : colors[2],
    'FF16' : colors[3],
    'RED20' : colors[1],
    'RED15' : colors[4],
    'FFQUAD' : colors[5]
}
plt.figure()
ax1 = sns.tsplot(
    data=DF_3b,
    time='time',
    value='fitness',
    unit='sim',
    condition='code',
    color=colordict,
    ci='sd',
)
# format plot
logging.info("Formatting 3B")
sns.despine()
plt.xlim([0, 1000])
plt.ylim([-0.05, 1.3])
plt.legend()
# plt.title('Hyper-evolvable', fontsize=labelsize)
plt.xlabel('Time (in generations)')
plt.ylabel('Mean Fitness')

for l in ax1.ax.lines:
    print(l.get_linewidth())
    plt.setp(l,linewidth=3)
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

# # move on to 3C:
# logging.info("Plotting 3C: Mean Fitness Traces")
# plt.figure()
# ax1 = sns.tsplot(
#     data=DF_3c,
#     time='time',
#     value='fitness',
#     unit='sim',
#     condition='code',
#     color=colordict,
#     ci='sd'
# )
# # format plot
# logging.info("Formatting 3C")
# sns.despine()
# plt.xlim([0, 1000])
# # plt.ylim([0, 1.3])
# plt.xticks([0, 500, 1000])
# plt.yticks([0, 0.15, 0.3])
# plt.legend()
# plt.title('Hypo-evolvable', fontsize=labelsize)
# plt.xlabel('Time (in generations)')
# plt.ylabel('Mean Fitness')
#
# # save output
# logging.info("Saving Figure 3C and inset to S3")
# figure_basename = '3c_inset_vector.svg'
# figure_path = '/home/ubuntu/' + figure_basename
# figure_s3path = s3_path + figure_basename
# plt.savefig(figure_path)
# with open(figure_path, 'rb') as data:
#     s3.upload_fileobj(data, bucketname, figure_s3path)
# success_string = (
#     "Success! Figure saved to {0}:{1}. ".format(bucketname, figure_s3path)
#     + "Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'".format(s3_region)
#     + " to see output figure file."
# )
# logging.info(success_string)
#
# ax1 = sns.tsplot(
#     data=DF_3b.loc[DF_3b['code'] == 'Standard Code'],
#     time='time',
#     value='fitness',
#     unit='sim',
#     condition='code',
#     color=colordict,
#     ci='sd'
# )
# plt.xticks([0, 200, 400, 600, 800, 1000])
# plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
# figure_basename = '3c_vector.svg'
# figure_path = '/home/ubuntu/' + figure_basename
# figure_s3path = s3_path + figure_basename
# plt.savefig(figure_path)
# with open(figure_path, 'rb') as data:
#     s3.upload_fileobj(data, bucketname, figure_s3path)
# success_string = (
#     "Success! Inset saved to {0}:{1}. ".format(bucketname, figure_s3path)
#     + "Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'".format(s3_region)
#     + " to see output figure file."
# )
# logging.info(success_string)

# # continue on to figure 3c and 3d
# logging.info('Calculating Lag Times and Growth Rates')
# sims = set(DF_3b['sim'])
# codes = set(DF_3b['code'])
# lags = []
# rates = []
# # loop over codes
# for code in tqdm(codes, desc='Looping over Codes'):
#     # declare storage variables
#     t_lag = np.zeros(len(sims))
#     rate = np.zeros(len(sims))
#     DF = DF_3b.loc[DF_3b['code'] == code]
#     for i, sim in enumerate(tqdm(sims, desc='Looping over Sims')):
#         # extract data for this sim
#         data = DF.loc[DF['sim'] == sim]
#         t = data['time'].values
#         f = data['fitness'].values
#         # smooth with gaussian filter
#         gaussian_filter = gaussian(30, 10)
#         filtered_signal = convolve1d(f, gaussian_filter/gaussian_filter.sum())
#         # calculate first derivative
#         delt = np.diff(t)
#         t_avg = (t[1:]+t[:-1])/2
#         filt_grad = np.diff(filtered_signal)/delt
#         # find peaks
#         peak_ind = peakutils.indexes(filt_grad, thres=0.05, min_dist=int(30/delt.mean()))
#         # get timestamp for this point
#         t_lag[i] = t_avg[peak_ind[0]]
#         t_ind = int(peak_ind[0])
#         # get estimate for evolutionary rate
#         dt = t[-1]  - t[t_ind]
#         dx = f[-1] - f[t_ind]
#         rate[i] = dx/dt
#     # store arrays in list
#     lags.append(t_lag)
#     rates.append(rate)
#
# # collate data into a dataframe
# logging.info('Collating Data into Dataframe')
# dfs = []
# for (lag, rate, code) in zip(lags, rates, codes):
#     d = pd.DataFrame({
#         'lag' : lag,
#         'rate' : rate,
#         'code' : [code for i in range(len(lag))]
#     })
#     dfs.append(d)
# DF_3cd = pd.concat(dfs)
# # plot 3c and save
# logging.info("Plotting 3C: Lag Time Distributions")
# # plot violin plots for lag times
# plt.figure()
# ax2 = sns.violinplot(
#     x='lag',
#     y='code',
#     data=DF_3cd,
#     palette=colordict,
#     inner='box',
#     order=wanted_codes
# )
# plt.title('Distribution of Lag Times (N=1000)')
# plt.xlabel('Lag Time (in generations)')
# sns.despine(trim=True)
# # save output
# logging.info("Saving Figure 3C to S3")
# figure_basename = '3c_vector.svg'
# figure_path = '/home/ubuntu/' + figure_basename
# figure_s3path = s3_path + figure_basename
# plt.savefig(figure_path)
# with open(figure_path, 'rb') as data:
#     s3.upload_fileobj(data, bucketname, figure_s3path)
# success_string = (
#     "Success! Figure saved to {0}:{1}. ".format(bucketname, figure_s3path)
#     + "Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'".format(s3_region)
#     + " to see output figure file."
# )
# logging.info(success_string)
#
# # plot 3D and save
# logging.info("Plotting 3D: Evolutionary Rate Distributions")
# # plot violin plots for lag times
# plt.figure()
# ax3 = sns.violinplot(
#     x='rate',
#     y='code',
#     data=DF_3cd,
#     palette=colordict,
#     inner='box',
#     order=wanted_codes
# )
# plt.title('Distribution of Evolutionary Rates (N=1000)')
# plt.xlabel(r'Evolutionary Rates (in \frac{1}{gen^2})')
# sns.despine(trim=True)
#
# # save output
# logging.info("Saving Figure 3D to S3")
# figure_basename = '3d_vector.svg'
# figure_path = '/home/ubuntu/' + figure_basename
# figure_s3path = s3_path + figure_basename
# plt.savefig(figure_path)
# with open(figure_path, 'rb') as data:
#     s3.upload_fileobj(data, bucketname, figure_s3path)
# success_string = (
#     "Success! Figure saved to {0}:{1}. ".format(bucketname, figure_s3path)
#     + "Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'".format(s3_region)
#     + " to see output figure file."
# )
# logging.info(success_string)
#
# # plot 3E and save
# logging.info("Plotting 3E: Endpoint Fitness Distributions")
# # get endpoint fitness from simulations
# endpoints = {}
# for code in tqdm(codes, desc='Looping through codes'):
#     endpoints[code] = []
#     df = DF_3b.loc[DF_3b['code'] == code]
#     for sim in tqdm(sims, desc='Looping through sims'):
#         endpoints[code].append(df.loc[df['sim'] == sim,'fitness'].iloc[-1])
# DF_endtimes = pd.DataFrame.from_dict(endpoints)
#
# # plot distribution
# plt.figure()
# for code in codes:
#     sns.distplot(
#         DF_endtimes[code],
#         kde=True,
#         hist=True,
#         rug=False,
#         norm_hist=True,
#         color=colordict[code],
#         label=code,
#         rug_kws={"alpha" : 0.03}
#     )
# sns.despine(trim=True)
# plt.xlabel('Endpoint Fitness')
# plt.ylabel('Probability')
# plt.legend()
# plt.title('Distribution of Endpoint Fitnesses')
# # save output
# logging.info("Saving Figure 3E to S3")
# figure_basename = '3e_vector.svg'
# figure_path = '/home/ubuntu/' + figure_basename
# figure_s3path = s3_path + figure_basename
# plt.savefig(figure_path)
# with open(figure_path, 'rb') as data:
#     s3.upload_fileobj(data, bucketname, figure_s3path)
# success_string = (
#     "Success! Figure saved to {0}:{1}. ".format(bucketname, figure_s3path)
#     + "Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'".format(s3_region)
#     + " to see output figure file."
# )
# logging.info(success_string)
