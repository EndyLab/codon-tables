import logging
logging.basicConfig(level=logging.INFO)

import numpy as np
from scipy.interpolate import interp1d as interp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
from tqdm import tqdm
from src.thunderflask import thunderflask
from src.bacteria import strain
from src.ffgen import ffgen
from src.codonTable import codonTable
from src.codonUtils import utils
from src.bacteria import strain
import requests
import os
import os.path
import boto3

# get environmental variables and define filepaths
datapath = os.environ['DATA_DIR']
paramfile = datapath + "params/" + os.environ['PARAM_FILE']
awsbucket = os.environ['AWS_BUCKET'] if 'AWS_BUCKET' in os.environ else ''

# make file directories if necessary
os.makedirs(os.path.dirname(paramfile), exist_ok=True)
logging.info("Starting simulation run")

# optionally download parameter file from S3 (for AWS Batch)
if awsbucket != "":
    s3 = boto3.resource('s3', region_name="us-west-1")
    if 'AWS_BATCH_JOB_ARRAY_INDEX' in os.environ:
        logging.info("Batch job: using updating param with index")
        paramfile = paramfile.format(os.environ['AWS_BATCH_JOB_ARRAY_INDEX'])

    logging.info("Downloading params from S3: {}/{}".format(awsbucket, paramfile))
    s3.Bucket(awsbucket).download_file(paramfile, paramfile)
    logging.info("Download complete")

# load parameters
logging.info("Loading param file {}".format(paramfile))
with open(paramfile, 'rb') as handle:
    param = pickle.load(handle)

# initialize simulation
logging.info("Initializing simulation")

sim_num = param['sim_num']
batch_num = param['batch_num']
strains = param['strains']
N_pop = param['N_pop']
T_0 = param['T_0']
T_sim = param['T_sim']
dt = param['dt']
t_extra = param['t_extra']
N_sims = param['N_sims']
mut_param = param['mut_param']
date = param['date']
code = param['code']
filepath = param['filepath']

# get output filename and make output filepath
filename = (
    '{0}_{1}_sim={2}_batch={3}_favg_traces_'
    'N_pop={4}e{5}=T={6}_N={7}_b={8}_l={9}.pickle'.format(
        date, code, sim_num, batch_num,
        str(N_pop)[0], int(np.log10(N_pop)),
        T_sim, N_sims, mut_param[0], mut_param[1]
    )
)

outfile = datapath + '/output/' + filepath + filename
os.makedirs(os.path.dirname(outfile), exist_ok=True)

logging.info("Running simulation")

# initialize list of dataframes and interpolation times
dataframes = []
newtimes = np.linspace(0, T_sim, int((T_sim)/dt))

# run N simulations
for i in tqdm(range(N_sims), desc='Simulation Number: '):
    sim = thunderflask(strains)
    sim.simulate(T_sim+t_extra, dt, T_0, mut_param)
    t = sim.f_avgtrace['timepoints']
    f_avg = sim.f_avgtrace['f_avg']
    interp_fxn = interp(t, f_avg)
    newf = interp_fxn(newtimes)
    df = pd.DataFrame({
        'time' : newtimes,
        'fitness' : newf,
        'sim' : [i for j in range(len(newf))],
        'code' : [code for j in range(len(newf))]
    })
    dataframes.append(df)

# package data into pandas dataframe
df_sc = pd.concat(dataframes)

# pickle results
with open(outfile, 'wb') as handle:
    pickle.dump(df_sc, handle)

# optionally upload output to S3 (for AWS Batch)
if awsbucket != "":
    logging.info("Uploading output to S3: {}/{}".format(awsbucket, outfile))
    s3 = boto3.resource('s3', region_name="us-west-1")

    with open(outfile, 'rb') as f:
        s3.Bucket(awsbucket).put_object(Key=outfile, Body=f)

    logging.info("Upload complete")
