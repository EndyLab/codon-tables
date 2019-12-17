
#import dependencies
import numpy as np
import matplotlib.pyplot as plt
from src.codonTable import codonTable
from src.codonUtils import utils
from src.codonoptimizer import MonteCarlo
from tqdm import tqdm

# Declare arrays to store statistics here

# Initialize simulator and initial codon table (default to Standard Table)
sim = MonteCarlo()
# Generate N tables and store statistics
N = 100
for i in tqdm(range(N)):
    table, Ws, Es = sim.GDA()
