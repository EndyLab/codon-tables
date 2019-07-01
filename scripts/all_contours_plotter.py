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

# define plotting parameters
figure_param = {
    'labelsize':16,
	'width': 4,
    'height' : 4 / 1.618,
    'font':'serif'
}

# set variables
bucketname = 'endylab-codon-table-simulations'
s3_path = 'manuscript/sup_figs/contour_plots/'
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
pbar = tqdm(filenames[1:9])
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

# repackage data
concat_file = 'lin_contour_concat.pickle'
logging.info("Pickling Dataframes")
with open('res/'+ concat_file, 'wb') as handle:
    pickle.dump(DF, handle, protocol=4)
logging.info("Uploading concatenated datafile")
with open('res/'+concat_file, 'rb') as data:
    s3.upload_fileobj(data, bucketname, s3_path+concat_file)
logging.info("Success! Concatenated datafile uploaded")

# # define useful helper functions
# def contain_probability(DF, code):
#     df = DF.loc[DF['code'] == code]
#     N_0 = list(set(df['N_0']))
#     N_0.sort()
#     num_sims = len(df.loc[df['N_0'] == N_0[0]].loc[0])
#     t = df.loc[df['sim'] == 0]['time']
#     contain_probability = np.zeros((len(N_0), len(t)))
#     for i, n_0 in enumerate(tqdm(N_0, desc='Processing Initial Conditions: ', leave=False)):
#         lildf = df.loc[df['N_0'] == n_0]
#         for j in tqdm(range(len(t)), desc='Processing sims: '):
#             weedf = lildf.loc[j]
#             contain_probability[i, j] = sum(weedf['popfrac'] == 0) / num_sims
#
#     return contain_probability, t, N_0
#
# def contour_plotter(DF, code, path, figure_param):
#     # unpack figure parameters
#     labelsize=figure_param['labelsize']
#     width = figure_param['width']
#     height = figure_param['height']
#
#     plt.rc('font', family=figure_param['font'])
#     plt.rc('xtick', labelsize=figure_param['labelsize'])
#     plt.rc('ytick', labelsize=figure_param['labelsize'])
#     plt.rc('axes', labelsize=figure_param['labelsize'])
#     # call contain_probability to get relevant arrays
#     contour, t, N_0 = contain_probability(DF, code)
#     # plot resulting contour
#     X, Y = np.meshgrid(t, np.array(N_0)/1e6)
#     CS = plt.contour(X, Y, contour, 20, cmap=plt.cm.viridis_r, vmin=0, vmax=1)
#     plt.clabel(CS, inline=1, fontsize=10)#, colors='black')
#     # format figure
#     ax = plt.gca()
#     plt.xticks([i*200 for i in range(6)])
#     plt.clim(0,1)
#     plt.title('{0}'.format(code if code != 'Standard' else 'Standard Code'), fontsize=labelsize)
#     plt.xlabel('Time (generations)')
#     plt.ylabel('Invasive Pop. Fraction')
#     sns.despine()
#     fig = plt.gcf()
#     fig.set_size_inches(width, height)
#     plt.savefig(path+'{0}_contour_lines.svg'.format(code))
# # plot contours
# path = '/home/ubuntu/'
# codes = set(DF['code'])
# for code in tqdm(codes, desc='Plotting Codes'):
#     # plot and format
#     logging.info("Plotting Contours: {0}".format(code))
#     plt.figure()
#     contour_plotter(DF, code, path, figure_param)
#     # upload to s3
#     figure_basename = '{0}_contour_lines.svg'.format(code)
#     figure_path = path + figure_basename
#     figure_s3path = s3_path + figure_basename
#     # testing save functionality
#     plt.savefig(figure_path)
#     logging.info("Uploading Plot: {0}".format(code))
#     with open(figure_path, 'rb') as data:
#         s3.upload_fileobj(data, bucketname, figure_s3path)
