import numpy as np
import matplotlib.pyplot as plt
from thunderflask import thunderflask
from bacteria import strain
from tqdm import tqdm
from ffgen import ffgen
from codonTable import codonTable
from bacteria import strain

# populate sim with a fast fail organism
table = ffgen.triplet()
LUCA = strain(N_pop=1e6, table=table, fitness=0, mu=2e-5)
sim = thunderflask(LUCA)
# initialize some variables
T_curr = 0
mut_param = [2, 3.5]
dt = 1

# run simulation
sim.simulate(1000, dt, T_curr, mut_param)
t = np.array(sim.f_trace['timepoints'])
f = np.array(sim.f_trace['fitnesses'])
fig, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(t, f)
axarr[0].set_title('Fast Fail Code: Mean Fitness vs Time')
dt = np.diff(t)
t_avg = (t[1:]+t[:-1])/2
gradf = np.diff(f)/dt
axarr[1].plot(t_avg, gradf)
axarr[1].set_title(r'Fast Fail Code: $\frac{df}{dt}$ vs Time')
plt.xlabel('Time (in generations)')
plt.ylabel('Fitness (%)')
plt.show()

for i, bact in enumerate(sim.estStrains):
    # if i % 10 == 0:
    t = bact.timepoints
    pop = bact.poptrace
    plt.semilogy(t, pop)
plt.title('Fast Fail Code: Population Traces for Established Strains')
plt.show()