#import dependencies
import numpy as np
import matplotlib.pyplot as plt
import random
from codonTable import codonTable
from codonUtils import utils
from bacteria import strain

# define simulator class
class thunderflask():
    '''Welcome to the Thunderdome of E. coli. A class used to simulate the
    genetic diversity of a batch culture of e. coli given a genetic code'''

    def __init__(self, strains=None, T=36000, dT=600, n_stochthresh=1e4,
        simflag=True):
        '''
        Initializer for thunderflask class.

        Parameters
        ----------
        -[bacteria.strain obj] strains: a single strain or list of strains to
            simulate
        - int T: the approximate total simulation time (in sec)
        - int dT: the approximate simulation time (in sec) for each module
        - int n_stochthresh: the population size past which the code treats a
            strain using the numerical integrator rather than the stochastic sim
        - bool simflag: an optional flag to run simulation upon initialization

        Returns
        -------
        thunderflask obj
        '''
        # if no strains given, initialize wild type E. coli with 100 CFU
        if strains == None:
            strains = [strain(100)]
        # if one strain given, package into list
        elif type(strains) != list:
            strains = [strains]
        # declare attributes for numerical and stochastic regimes for strains
        self.bigStrains = []
        self.smallStrains = strains
        # partition initial strain list appropriately
        self.strainShuffle()
        # 

    def simulate(self):
        '''
        Main method of thunderflask class. Used to run genetic diversity
        simulation. There are four main modules in this simulation:

        (1) Stochastic -> (2) Numerical -> (3) Pop. Management -> (4) Mutations

        Consult documentation of each method for specific information. In brief,
        small and large population strains are handled separately in (1) and
        (2), and new mutations are generated in (4).

        Parameters
        ----------

        Returns
        -------
        '''
        # set up main loop

        # main loop

            # stochastic simulation

            # numerical simulation

            # strain swapping (if low pop strains get large or vice versa)

            # mutation simulation

        # return results of simulation
        return

    def stochSim(self):
        '''
        '''
        return

    def numericalSim(self):
        '''
        '''
        return

    def strainShuffle(self):
        '''
        '''
        return

    def mutationSim(self):
        '''
        '''
        return
