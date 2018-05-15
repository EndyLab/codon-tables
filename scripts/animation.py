import logging
logging.basicConfig(level=logging.INFO)

# import necessary modules
import numpy as np
import pandas as pd
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import matplotlib.animation as animation
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
import copy
import boto3
from src.codonTable import codonTable
from src.codonUtils import utils
from src.thunderflask import thunderflask
from src.bacteria import strain
from src.codonOptimizer import tableOptimizer

# simulation video

# control aesthetics
labelsize=20
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=labelsize)
plt.rc('ytick', labelsize=labelsize)
plt.rc('axes', labelsize=labelsize)

# populate sim
LUCA = strain(N_pop=1e6, fitness=0, mu=2e-5)
sim = thunderflask(LUCA)
# initialize some variables
T_curr = 0
mut_param = [1, 2]
dt = 0.1
T_sim = 500

# run simulation
logging.info("Starting Simulation Run")
sim.simulate(T_sim, dt, T_curr, mut_param, save_all=True, prune_strains=True,
             show_progress=True)
logging.info("Simulation Complete")

# make animation
logging.info("Animating Results")
# set up list of strains
n = len(sim.allStrains)
colors = pl.cm.viridis(np.linspace(0,1,n))
strainlist = [(time_ind, bact) for time_ind, bact in enumerate(sim.allStrains)]

sorting_list = []
for time_ind, bact in tqdm(strainlist, desc='Looping through all strains'):
    try:
        pop_size = max(bact.poptrace)
    except:
        pop_size = 0
    sorting_list.append((pop_size, time_ind, bact))
sortedlist = sorted(sorting_list, key=lambda x: x[0])
endlist = [(pop_ind, time_ind, bact) for pop_ind, (pop_size, time_ind, bact) in enumerate(reversed(sortedlist))]
shuffle(endlist)

# set plot characteristics
labelsize=20
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=labelsize)
plt.rc('ytick', labelsize=labelsize)
plt.rc('axes', labelsize=labelsize)

width = 6
height = 1.75* width
# set video characteristics
fps = 30
vid_length = 5 # in sec
tot_frames = fps*vid_length

def framer(nFrame):
    # calculate curr_time
    curr_time = nFrame*(500/tot_frames)
    # plot pop traces
    fig, axarr = plt.subplots(2, sharex=True)
    print('Frame {0}/{1}'.format(nFrame, tot_frames))
    for pop_ind, time_ind, bact in endlist:
        if (pop_ind < 30) or (pop_ind % 20 ==0) :
            t = bact.timepoints
            pop = bact.poptrace
            stopind = next((int(ind) for ind, time in enumerate(t) if time >= curr_time), -1)
            axarr[0].semilogy(t[:stopind], pop[:stopind], color=colors[time_ind])
    axarr[0].set_yticks([])
    axarr[0].set_xticks([])
    # plot mean fitness
    t = np.array(sim.f_avgtrace['timepoints'])
    f = np.array(sim.f_avgtrace['f_avg'])
    stopind = next((int(ind) for ind, time in enumerate(t) if time >= curr_time), -1)
    axarr[1].plot(t[:stopind], f[:stopind], 'k')
    plt.xlabel('Time (gen)')
    axarr[0].set_ylabel('Pop. Frac')
    axarr[1].set_ylabel('Fitness (1/gen)')
    plt.suptitle('Simulating Population...', fontsize=labelsize)
    plt.xlim([0, 500])
    axarr[0].set_ylim([1e0, 10**(6.2)])
    axarr[1].set_ylim([0, 0.4])
    fig.set_size_inches(width, height)
    ax = plt.gca()
    sns.despine()
    plt.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        left=False,
        labelbottom=False,
        labelleft=False
    )
# set S3 variables
bucketname = 'endylab-codon-table-simulations'
s3_path = 'quals/simulation_vid/'
s3_region = 'us-west-1'
fig = plt.gcf()
figure_basename = 'sim_vid.gif'
figure_path = '/home/ubuntu/' + figure_basename
figure_s3path = s3_path + figure_basename
anim = animation.FuncAnimation(fig, framer, frames=tot_frames)
anim.save(figure_path, writer='imagemagick', dpi=dpi, fps=fps);
fig = plt.figure()

logging.info("Uploading Movie to S3")
with open(figure_path, 'rb') as data:
    s3.upload_fileobj(data, bucketname, figure_s3path)
success_string = (
    "Success! Gif saved to {0}:{1}. ".format(bucketname, figure_s3path)
    + "Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'".format(s3_region)
    + " to see output figure file."
)
logging.info(success_string)
