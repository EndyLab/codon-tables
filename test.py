import numpy as np
import matplotlib.pyplot as plt
from thunderflask import thunderflask
from bacteria import strain
from tqdm import tqdm

# run several iterations
# populate sim
LUCA = strain(N_pop=1e6, fitness=0, mu=2e-5)
sim = thunderflask(LUCA)
# initialize some variables
T_curr = 0
mut_param = [1, 2]
dt = 1

# def iterate(sim, iteration, dt, T_curr, f_avg, mut_param):
#     # run mutation sim
#     sim.mutationSim(dt, T_curr, mut_param)
#     # run stochastic sim
#     T_next, taus = sim.stochSim(T_approx=1, T_curr=T_curr, f_avg=f_avg)
#     # do analytic simulation
#     sim.analyticSim(T_curr=T_curr, taus=taus, f_avg=f_avg)
#
#     # shuffle strains and look at changes
#     sim.strainShuffle(T_curr=T_next, f_avg=f_avg)
#     # update f_avg
#     f_avg = sim.updateF_avg()
#
#     return T_next, f_avg, sim
#
# for i in tqdm(range(1,331)):
#     T_curr, f_avg, sim = iterate(sim, i, dt, T_curr, f_avg, mut_param)

# run simulation
# import ipdb; ipdb.set_trace()
sim.simulate(500, dt, T_curr, mut_param)

for bact in sim.estStrains:
    t = bact.timepoints
    pop = bact.poptrace
    plt.semilogy(t, pop)
plt.title('Analytic Simulation of Big Strains')
plt.show()
