from src.bacteria import strain
from src.parameterizer import *
import pickle

# initialize LUCA strain
N_pop = 3e5
LUCA = strain(N_pop=N_pop)

# initialize some variables
sim_num = 0
batch_num = 0
strains = [LUCA]
T_0 = 0
T_sim = 100
dt = 0.1
t_extra = 5
N_sims = 100
mut_param = [1, 2]
date = '3-15'
code = 'Standard Code'
filepath = 'test1'

# pickle that shiznit
num_cores = 32
params = genParamDict(sim_num, batch_num, strains, N_pop, T_0, T_sim, dt, t_extra,
                      N_sims, mut_param, date, code, filepath)
paramDicts = batcher(params, num_cores)
for params in paramDicts:
    batch_num = params['batch_num']
    n_sim = params['N_sims']
    print('Batch Num: {0}; n_sims = {1}'.format(batch_num, n_sim))

pickle_path = '../data/params/'
pickle_file = '{0}_{1}_{2}_{3}_params.pickle'.format(date, code, sim_num, batch_num)
# import ipdb; ipdb.set_trace()
paramPickler(paramDicts, pickle_path, pickle_file)
