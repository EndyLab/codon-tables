# import necessary modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from thunderflask import thunderflask
from bacteria import strain
from tqdm import tqdm
import cProfile

# populate sim
LUCA = strain(N_pop=1e6, fitness=0, mu=2e-5)
sim = thunderflask(LUCA)
# initialize some variables
T_curr = 0
mut_param = [1, 2]
dt = 1

#profile sim
cProfile.run('sim.simulate(300, dt, T_curr, mut_param, save_all=True, prune_strains=False)')
