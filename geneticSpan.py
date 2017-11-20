# import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Bio.pairwise2 as distance # this is used to align two strings
from codonTable import codonTable
from codonUtils import utils

class geneticSpan:
    '''A class designed to quantify the 'genetic spanning rate' of a given
    codon table. The genetic spanning rate given a set of proteins is defined
    as the sum of indels/point mutations in the DNA sequence (i.e. distance)
    required to go from one member to another, over all pairs. This can be
    normalized by the spanning rate of the Standard Table for the same set.
    '''
    def __init__(self):
