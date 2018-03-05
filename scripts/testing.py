# import necessary modules
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import pickle
from src.codonTable import codonTable
from src.codonUtils import utils
from src.thunderflask import thunderflask
from src.bacteria import strain

sns.set_context("paper")
# get fast fail table
with open('res/ff_table_manuscript.pickle', 'rb') as handle:
        fftable = codonTable(table=pickle.load(handle))
# get random table
with open('res/random_table_manuscript.pickle', 'rb') as handle:
    rand = codonTable(table=pickle.load(handle))
# get standard code
sc = codonTable()  
# get colorado code
col = codonTable(table=utils.coloradoTable)

# define colors
orange = '#ef6c00'
blue = '#2196f3'
green = '#008000'

# plot and save some shiznit
path = '/home/jonathan/Lab/Fast Fail/Figures/Figure 1/'
rand.plotGraph(filename=path+'rand_graph.svg')
#sc.plotGraph(filename=path+'sc_graph.svg')
#col.plotGraph(filename=path+'col_graph.svg')
