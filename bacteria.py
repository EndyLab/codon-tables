#import dependencies
import numpy as np
import matplotlib.pyplot as plt
from uuid import uuid4 as uuid
from codonUtils import utils

class strain():
    '''A class used to represent a bacterial strain in culture'''
    def __init__(self, N_pop=100, pd=0.5, pm=0.5, kd=0.3, km=0.2, td=10, tm=20,
        mu=1e-3, t_large=0, codonTable=None, ID=None, lineage=[]):
        '''The init function for the strain class.

        Parameters
        ----------
        - int/float N_pop: current population size of strain
        - float pd, pm: maximum probability of cell division/death respectively
        - float kd, km: time dependence of cell division/death respectively
        - float td, tm: lag time for cell division/death respectively
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
        self.growParam = [pd, pm, kd, km, td, tm]
        self.t_large = t_large
        self.mu = mu
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
