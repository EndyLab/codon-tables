import pickle

def genParamDict(sim_num, batch_num, T_curr, T_sim, N_sims, mut_param, dt, 
                 t_extra, date, code, filepath):
    '''A function that is used to generate a dictionary of parameters, given a
    set of inputs. Used to pass parameters to fit_traces.py and derivative
    scripts via pickle file'''
    return locals()


# initialize some variables
sim_num = 0
batch_num = 0
T_curr = 0
mut_param = [1, 2]
dt = 0.1
N_sims = 10
T_sim = 10
t_extra = 5
date = '3-13'
code = 'Colorado'
filepath = '/home/jonathan/Lab/Fast Fail/Traces'

# pickle that shiznit
params = genParamDict(sim_num, batch_num, T_curr, T_sim, N_sims, mut_param, 
                      dt, t_extra, date, code, filepath)
pickle_path = '../'
pickle_file = '{0}_{1}_{2}_{3}_params.pickle'.format(date, code, sim_num,
                                                     batch_num)
with open(pickle_path+pickle_file, 'wb') as handle:
    pickle.dump(params, handle)
