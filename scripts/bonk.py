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

# get appropriate table
table = {
    'Standard Code' : utils.standardTable,
    'Fast Fail' : ffgen.triplet(),
    'Colorado' : utils.coloradoTable
}
# initialize some variables
T_curr = 0
mut_param = [0.5, 0.33]
dt = 0.1
N_sims = 3
T_sim = 100
t_extra = 5
date = '2-21'
code = 'Fast Fail'
filepath = 'res/2-14 Traces'
filename = '{0}_{1}_favg_traces_T={2}_N={3}_b={4}_l={5}.pickle'.format(date,
                                                                       code,
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
# plot results
ax = sns.tsplot(data=df_sc, time='time', value='value', unit='sim')
plt.title('{0}: <F> vs time ({1} Replicates)'.format(code, N_sims))
plt.xlabel('Time (in generations)')
plt.ylabel('Mean Fitness')
plt.show()
print('done')
# fig, axarr = plt.subplots(2, sharex=True)
# axarr[0].plot(t, f)
# axarr[0].set_title('Fast Fail Code: Mean Fitness vs Time')
# dt = np.diff(t)
# t_avg = (t[1:]+t[:-1])/2
# gradf = np.diff(f)/dt
# axarr[1].plot(t_avg, gradf)
# axarr[1].set_title(r'Fast Fail Code: $\frac{df}{dt}$ vs Time')
# plt.xlabel('Time (in generations)')
# plt.ylabel('Fitness (%)')
# plt.show()
