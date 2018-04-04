import numpy as np
from scipy.interpolate import interp1d as interp
import matplotlib.pyplot as plt
from src.codonTable import codonTable
from src.codonUtils import utils
from src.thunderflask import thunderflask
from src.bacteria import strain
from src.ffgen import ffgen

# populate sim with two strains
fftable = ffgen.triplet()
failer = strain(N_pop=0.5*1e6, table=fftable, code='FF20', fitness=0, mu=2e-5)
WT = strain(N_pop=0.5*1e6, fitness=0, mu=2e-5)
strains = [failer, WT]
sim = thunderflask(strains=strains)
# initialize some variables
T_curr = 0
mut_param = [1, 2]
dt = 0.1
# run simulation
sim.simulate(T=1000, dt=dt, mut_param=mut_param, save_all=True, competition=True)
t = np.array(sim.f_avgtrace['timepoints'])
x1 = np.array(sim.popfrac['Standard Code'])
x2 = np.array(sim.popfrac['FF20'])
plt.plot(t, x1, label='Standard Code')
plt.plot(t, x2, label='FF20')
plt.title('Fast Fail Code: Mean Fitness vs Time')
plt.xlabel('Time (in generations)')
plt.ylabel('Population Fraction')
plt.legend()
plt.show()

sim.f_avgtrace['timepoints'].append(1000)
alt_t = np.array(sim.f_avgtrace['timepoints'])
sim.popfrac['FF20'].append(0)
alt_x2 = np.array(sim.popfrac['FF20'])
fxn = interp(alt_t, alt_x2)
newt = np.linspace(0, 1000, int(1000/dt) + 1)
newx2 = fxn(newt)

plt.plot(t, x2, label='FF20')
plt.plot(newt, newx2, label='interpolated')
plt.legend()
plt.show()
