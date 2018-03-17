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
import logging

logging.basicConfig(filename="/simulation.log",level=logging.DEBUG)

# get environmental variables and define filepaths
datapath = os.environ['DATA_DIR']
paramfile = datapath + "/params/" + os.environ['PARAM_FILE']
awsbucket = os.environ['AWS_BUCKET'] if 'AWS_BUCKET' in os.environ else ''

os.makedirs(os.path.dirname(paramfile), exist_ok=True)
logging.info("Starting simulation run")

if awsbucket != "":
    logging.info("Downloading params from S3: {}/{}".format(awsbucket, paramfile))
    s3 = boto3.resource('s3')
    s3.Bucket(awsbucket).download_file(paramfile, paramfile)
    logging.info("Download complete")

logging.info("Loading param file", paramfile)
with open(paramfile, 'rb') as handle:
    param = pickle.load(handle)

# response = requests.get(paramfile)
# param = pickle.load(response.data)

logging.info("Initializing simulation")

# initialize variables
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
# get output filename
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

# initialize list of dictionaries of arrays (i know, it's too much)
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
        'value' : newf,
        'sim' : [i for j in range(len(newf))],
        'code' : [code for j in range(len(newf))]
    })
    dataframes.append(df)
# package data into pandas dataframe
df_sc = pd.concat(dataframes)
# pickle results
with open(outfile, 'wb') as handle:
    pickle.dump(df_sc, handle)
