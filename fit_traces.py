import numpy as np
from scipy.interpolate import interp1d as interp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
from tqdm import tqdm
from thunderflask import thunderflask
from bacteria import strain
from ffgen import ffgen
from codonUtils import utils
from codonTable import codonTable
from bacteria import strain

# set up dictionary of tables
tables = {
    'Standard Code' : utils.standardTable,
    'Fast Fail' : ffgen.triplet(),
    'Colorado' : utils.coloradoTable
}
# initialize some variables
T_curr = 0
mut_param = [1, 2]
dt = 0.1
N_sims = 10
T_sim = 250
t_extra = 5
date = '2-15'
code = 'Standard Code'
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
newtimes2 = np.linspace(1, T_sim, len(newtimes))
print(len(newtimes) == len(newtimes2))
# run N simulations
for i in tqdm(range(N_sims), desc='Simulation Number: '):
    LUCA = strain(N_pop=1e6, table=tables[code], fitness=0, mu=2e-5)
    sim = thunderflask(LUCA)
    sim.simulate(T_sim+t_extra, dt, T_curr, mut_param)
    t = np.array(sim.f_avgtrace['timepoints'])
    f_avg = np.array(sim.f_avgtrace['f_avg'])
    interp_fxn = interp(t, f_avg)
    newf = interp_fxn(newtimes)
    dts = np.diff(t)
    t_avg = (t[1:]+t[:-1])/2
    gradf = np.diff(f_avg)/dts
    grad_fxn = interp(t_avg, gradf)
    newgrad = grad_fxn(newtimes2)
    df = pd.DataFrame({
        'time' : newtimes, 
        'f_avg' : newf, 
        'dfdt' : newgrad,
        'sim' : [i for j in range(len(newf))],
        'code' : [code for j in range(len(newf))]
    })
    dataframes.append(df)
# package data into pandas dataframe
df_sc = pd.concat(dataframes)
# pickle results
with open('{0}/{1}'.format(filepath, filename), 'wb') as handle:
    pickle.dump(df_sc, handle)
# plot results
ax = sns.tsplot(data=df_sc, time='time', value='dfdt', unit='sim')
plt.title('{0}: df/dt vs time ({1} Replicates)'.format(code, N_sims))
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

