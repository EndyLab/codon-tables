#import dependencies
import numpy as np
import matplotlib.pyplot as plt
from uuid import uuid4 as uuid
from codonUtils import utils

class strain():
    '''A class used to represent a bacterial strain in culture'''
    def __init__(self, N_pop=100, fitness=0, mu=1e-3,
        t_large=0, codonTable=None, ID=None, lineage=[]):
        '''The init function for the strain class.

        Parameters
        ----------
        - int/float N_pop: current population size of strain
        - float fitness: the absolute fitness of an individual strain
        - float mu: mutation rate of bacterial strain (1/genome-gen)
        - float t_large: timestep in which a strain becomes large enough to
            simulate analytically
        - dict codonTable: genetic code for bacterial strain
        - str ID: a unique string identifier for bacterial strain
        - list <str> lineage: a list tracking the lineage of bacterial strain

        Returns
        -------
        strain obj: returns handle for strain instance
        '''
        # generate values for optional parameters
        if codonTable == None:
            codonTable = utils.standardTable
        if ID == None:
            ID = uuid()
        # store parameters in class attributes
        self.N_pop = N_pop
        self.fitness = fitness
        self.mu = mu
        self.t_large = t_large
        self.codonTable = codonTable
        self.ID = str(ID)
        self.lineage = lineage
        self.timepoints = []
        self.poptrace = []

# run as script
if __name__ == '__main__':
    print(help(strain))
    test_strain = strain()
    print('all good!')
