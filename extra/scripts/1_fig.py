import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import pickle
from src.codonTable import codonTable
from src.codonUtils import utils

def fig1a_c():
    # get fast fail table
    with open('res/ff_table_manuscript.pickle', 'rb') as handle:
        fftable = codonTable(table=pickle.load(handle))

    # get standard code
    sc = codonTable()  

    # get colorado code
    col = codonTable(table=utils.coloradoTable)

    # plot and save some shiznit
    path = '/home/jonathan/Lab/Fast Fail/Figures/Figure 1/'
    fftable.plotGraph(filename=path+'ff_graph.svg')
    sc.plotGraph(filename=path+'sc_graph.svg')
    col.plotGraph(filename=path+'col_graph.svg')
    # return
    return

def fig1d():
    return

