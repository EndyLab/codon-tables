import numpy as np
from scipy.interpolate import interp1d as interp
import matplotlib.pyplot as plt
from src.codonTable import codonTable
from src.codonUtils import utils
from src.thunderflask import thunderflask
from src.bacteria import strain
from src.ffgen import ffgen
from src.parameterizer import *

# calculate number of replicates per core
N_sims = 500
num_cores = 256
N_per = int(N_sims/num_cores)
# handle case where num_cores > N_sims
if N_per == 0:
    for batch_num in range(N_sims):
        pass
# handle case where num_cores < N_sims
else:
    modulo = N_sims % num_cores
    for batch_num in range(num_cores):
        pass

shorthand = min(N_sims, num_cores)
print(batch_num == (shorthand - 1))
