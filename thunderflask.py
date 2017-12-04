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
        self.deadStrains = []
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
        # main loop
        T_tot = 0
        index = 0
        #while T_tot < T:
            # stochastic simulation
            #T_stoch, ind_stoch = self.stochSim(dT, timestep, index)
            # numerical simulation

            # strain swapping (if low pop strains get large or vice versa)

            # mutation simulation

        # return results of simulation
        return

    def stochSim(self, T_approx, T_0):
        '''
        Method used to perform stochastic simulation of bacterial growth
        dynamics. Used to handle strains whose total population size is less
        than a given threshold.

        Parameters
        ----------
        - float T_approx: the approximate simulation time (in sec)
        - float T_0: the current simulation time (in sec)

        Returns
        ------
        - float T_elapsed: the actual simulation time (in sec)
        - np.array taus: the time steps between reactions (in sec)
        '''
        ###############################################
        # Initialize Stochastic Simulation attributes #
        ###############################################
        # Declare number of strains
        numStrains = len(self.smallStrains)
        # calculate a_tot for initial conditions to predict average timestep
        a_i = self.__rxnpropensity(numStrains, T_0)
        a_tot = a_i.sum()
        # declare number of iterations to perform
        numIter = int(np.ceil(T_approx*a_tot)) # dt ~ 1/atot --> iter = T/dt
        # construct dictionary of possible reactions
        rxndict = self.__rxndict(numStrains)
        #####################
        # Declare variables #
        #####################
        # declare numpy array of traces for this epoch
        trace = np.zeros((numStrains, numIter))
        # declare numpy array for storing time steps
        taus = np.zeros(numIter)
        # declare initial population sizes
        for i, strain in enumerate(self.smallStrains):
            trace[i, 0] = strain.N_pop
        # declare iteration counter variables
        T_elapsed = 0
        #############
        # Main Loop #
        #############
        for timeind in range(1, numIter):
            # calculate reaction propensities
            a_i = self.__rxnpropensity(numStrains, T_0 + T_elapsed)
            # calculate time to next reaction
            a_cumsum = np.cumsum(a_i)
            a_tot = a_cumsum[-1]
            tau = np.log(1/np.random.rand())/a_tot
            taus[timeind] = tau
            # choose next reaction
            rxnval = a_tot * np.random.rand()
            i = np.argmax(a_cumsum > rxnval)
            # update population sizes
            trace[:,timeind] = trace[:,timeind-1] + rxndict[i]
            # update T_elapsed
            T_elapsed += tau

        # update population size and timepoints for each strain
        for i, strain in enumerate(self.smallStrains):
            strain.N_pop = trace[i, -1]
            strain.timepoints += (np.cumsum(taus) + T_0).tolist()
            strain.poptrace += trace[i,:].tolist()
        # return T_elapsed and time intervals
        return T_elapsed, taus

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

    #######################
    #   Private Methods   #
    #######################
    def __rxnpropensity(self, numStrains, t):
        '''A private method used to calculate the reaction propensities
        of the stochastic regime strains.

        Parameters
        ----------
        - int numStrains: the number of strains to simulate
        - float t: the current simulation time

        Returns
        -------
        np.array a_i: numpy array of reaction propensities
        '''
        # declare array for reaction propensities
        a_i = np.zeros(2*numStrains)
        # loop through strains and calculate the birth/death propensities
        for i, strain in enumerate(self.smallStrains):
            pd, pm, kd, km, td, tm = strain.growParam
            growprop = pd/(1 + np.exp(kd*(td - t)))
            deathprop = pm/(1 + np.exp(km*(tm - t)))
            a_i[2*i] = growprop*strain.N_pop # growth propensity
            a_i[2*i+1] = deathprop*strain.N_pop # death propensity
        # return results
        return a_i

    def __rxndict(self, numStrains):
        '''A private method used to generate a reaction dictionary for the
        stochastic simulator. Each item is a numpy array of length numStrains
        which describes the change in each species count

        Parameters
        ----------
        int numStrains: the number of strains to simulate

        Returns
        -------
        dict rxndict: python dict representing the possible reactions
        '''
        # declare rxndict to return
        rxndict = {}
        # generate two arrays for each strain (grow and death) and store
        for i in range(numStrains):
            growarray = np.zeros(numStrains)
            # set ith element to be 1
            growarray[i] = 1
            # even elements are growth, odd are death
            rxndict[2*i] = growarray
            rxndict[2*i+1] = -growarray
        # return resulting dict
        return rxndict

# debug script
if __name__ == '__main__':
    # generate some strains
    strains = []
    for i in range(10):
        strains.append(strain())
    sim = thunderflask(strains)
    T_elapsed, taus = sim.stochSim(600,5)
    print('T_elapsed = {0}'.format(T_elapsed))
    for i in range(5):
        t = sim.smallStrains[i].timepoints
        trace = sim.smallStrains[i].poptrace
        plt.plot(t, trace, label='Strain {0}'.format(i))
    plt.legend()
    plt.show()
