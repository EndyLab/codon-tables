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
from sys import argv

# get appropriate table
table = {
    'Standard Code' : utils.standardTable,
    'Fast Fail' : ffgen.triplet(),
    'Colorado' : utils.coloradoTable
}
# get dictionary of parameters from passed pickle file
filename = argv[1]
with open(filename, 'rb') as handle:
    param = pickle.load(handle)
# initialize variables
sim_num = param['sim_num']
batch_num = param['batch_num']
T_curr = param['T_curr']
mut_param = param['mut_param']
dt = param['dt']
N_sims = param['N_sims']
T_sim = param['T_sim']
t_extra = param['t_extra']
date = param['date']
code = param['code']
filepath = param['filepath']
filename = '{0}_{1}_sim={2}_batch={3}_favg_traces_T={2}_N={3}_b={4}_l={5}.pickle'.format(date,
                                                                              code,
                                                                              sim_num,
                                                                              batch_num,
                                                                              T_sim,
                                                                              N_sims,
                                                                              mut_param[0], 
                                                                              mut_param[1]) 
# initialize list of dictionaries of arrays (i know, it's too much) 
dataframes = []
newtimes = np.linspace(0, T_sim, int((T_sim)/dt))
# run N simulations
for i in tqdm(range(N_sims), desc='Simulation Number: '):
    LUCA = strain(N_pop=1e6, table=table[code], fitness=0, mu=2e-5)
    sim = thunderflask(LUCA)
    sim.simulate(T_sim+t_extra, dt, T_curr, mut_param)
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
with open('{0}/{1}'.format(filepath, filename), 'wb') as handle:
    pickle.dump(df_sc, handle)
