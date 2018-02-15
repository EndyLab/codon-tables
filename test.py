import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
# import seaborn as sns
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
dt = 0.3

# run simulation
sim.simulate(100, dt, T_curr, mut_param, save_all=True, prune_strains=True,
             show_progress=False)
t = np.array(sim.f_avgtrace['timepoints'])
f = np.array(sim.f_avgtrace['f_avg'])
fig, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(t, f)
axarr[0].set_title('Standard Code: Mean Fitness vs Time')
dt = np.diff(t)
t_avg = (t[1:]+t[:-1])/2
gradf = np.diff(f)/dt
axarr[1].plot(t_avg, gradf)
axarr[1].set_title(r'Standard Code: $\frac{df}{dt}$ vs Time')
plt.xlabel('Time (in generations)')
plt.ylabel('Fitness (%)')
plt.show()

# t = []
# f = []
# for i, time in tqdm(enumerate(sim.f_trace['timepoints'])):
#     for fitness in sim.f_trace['fitnesses'][i]:
#         t.append(time)
#         f.append(fitness)
# t = np.array(t)
# f = np.array(f)
# with sns.axes_style("white"):
#     sns.jointplot(x=t, y=f, kind="hex", color='k')
# plt.xlabel('Time (in generations)')
# plt.ylabel('Fitness (%)')
# plt.show()

n = len(sim.allStrains)
colors = pl.cm.viridis(np.linspace(0,1,n))
for i, bact in enumerate(sim.allStrains):
    t = bact.timepoints
    pop = bact.poptrace
    plt.semilogy(t, pop, color=colors[i])
# for i, (key, (times, pops)) in enumerate(sim.poptrace.items()):
#     if i % 10  == 0:
#         plt.semilogy(times, pops, color=colors[i])
plt.xlabel('Time (in generations)')
plt.ylabel('Population Size')
plt.title('Standard Code: Population Traces for Established Strains')
plt.show()
