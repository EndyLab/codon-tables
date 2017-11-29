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

    def __init__(self, strains=None):
        '''
        Initializer for thunderflask class.

        Parameters
        ----------
        [bacteria.strain obj] strains: a single strain or list of strains to
            simulate

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
        # declare attributes
        self.bigStrains = []
        self.smallStrains = strains
        self.poptrace = {}
        # partition initial strain list appropriately
        self.strainShuffle()

    def simulate(self, T=36000, dT=600, timestep=1, n_stochthresh=1e4):
        '''
        Main method of thunderflask class. Used to run genetic diversity
        simulation. There are four main modules in this simulation:

        (1) Stochastic -> (2) Numerical -> (3) Pop. Management -> (4) Mutations

        Consult documentation of each method for specific information. In brief,
        small and large population strains are handled separately in (1) and
        (2), and new mutations are generated in (4).

        Parameters
        ----------
        - float T: the approximate total simulation time (in sec)
        - float dT: the approximate simulation time (in sec) for each module
        - float timestep: the temporal resolution (in sec) for population traces
        - float n_stochthresh: the population size past which the code treats a
            strain using the numerical integrator rather than the stochastic sim
        Returns
        -------
        '''
        # initialize population trace arrays for initial strains
        for strain in self.bigStrains + self.smallStrains:
            numSteps = np.ceil(T/timestep)
            self.poptrace[strain.ID] = np.zeros(numSteps)
            self.poptrace[strain.ID][0] = strain.N_pop
        # main loop
        T_tot = 0
        count = 0
        while T_tot < T:
            # stochastic simulation

            # numerical simulation

            # strain swapping (if low pop strains get large or vice versa)

            # mutation simulation

        # return results of simulation
        return

    def stochSim(self, T_approx, dt, startind):
        '''
        Method used to perform stochastic simulation of bacterial growth
        dynamics. Used to handle strains whose total population size is less
        than a given threshold.

        Parameters
        ----------
        - float T_approx: the approximate simulation time (in sec)
        - float dt: the time resolution of population traces (in sec)
        - int startind: the index of the simulation start time in pop traces

        Returns
        ------
        - float T_elapsed: the actual simulation time (in sec)
        - int index: the number of timesteps performed
        '''
        # declare numpy array of traces for this epoch with some buffer
        numIter = np.ceil(T_approx/dt) + T_approx*0.1 # 10% buffer for safety
        numStrains = len(self.smallStrains)
        trace = np.zeros(numStrains, numIter)
        # declare initial population sizes
        for i, strain in enumerate(self.smallStrains):
            trace[i, 0] = strain.N_pop
        # declare array of reaction propensities
        a_i = np.zeros(2*numStrains)
        # declare iteration counter variables
        T_elapsed = 0
        index = 0
        # loop while T_approx has not been reached
        while (T_elapsed < T_approx):
            # calculate reaction propensities
            for i, strain in enumerate(self.smallStrains):
                '''
                growprop = ????
                deathprop = ?????
                a_i[2*i] = growprop*strain.N_pop # growth propensity
                a_i[2*i+1] = deathprop*strain.N_pop # death propensity
                '''
            # calculate time to next reaction
            a_cumsum = np.cumsum(a_i)
            a_tot = a_cumsum[-1]
            tau = np.log(1/np.random.rand())/a_tot
            # choose next reaction
            rxnval = a_tot * np.random.rand()
            for i in range(numStrains):
                if a_cumsum[i] < rxnval:
                    continue
                else:
                    break
            strainind = i//2
            # even rxns are growth, odd are death
            popchange = (-1)**i
            # propagate populations
            dIndex = np.ceil(tau/dt)
            '''
            make this a tiling thing, does not automatically broadcast
            traces[:, index:index+dIndex+1] = traces[:, index]
            '''
            trace[strainind, index+dIndex+1] += popchange
            # update T_elapsed and index
            T_elapsed += tau
            index += dIndex

        # update current population size for each strain and package traces
        for i, strain in enumerate(self.smallStrains):
            strain.N_pop = trace[i, index]
            self.poptrace[strain.ID][startindex:index+1] = trace[i,:]
        # return T_elapsed and current index
        return T_elapsed, index

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
