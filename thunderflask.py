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
            numSteps = int(np.ceil(T/timestep))
            self.poptrace[strain.ID] = np.zeros(numSteps)
            self.poptrace[strain.ID][0] = strain.N_pop
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

    def stochSim(self, T_approx):
        '''
        Method used to perform stochastic simulation of bacterial growth
        dynamics. Used to handle strains whose total population size is less
        than a given threshold.

        Parameters
        ----------
        float T_approx: the approximate simulation time (in sec)

        Returns
        ------
        - float T_elapsed: the actual simulation time (in sec)
        - np.array taus: the time steps between reactions (in sec)
        '''
        ###############################################
        # Initialize Stochastic Simulation attributes #
        ###############################################
        # calculate a_tot for initial conditions to predict average timestep
        a_i = self.__rxnpropensity()
        a_tot = a_i.sum()
        # declare number of iterations to perform and number of strains
        numIter = int(np.ceil(T_approx*a_tot)) # dt ~ 1/atot --> iter = T/dt
        numStrains = len(self.smallStrains)
        # construct dictionary of possible reactions
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
        # DEBUG: track reaction choices
        ivals = []
        #############
        # Main Loop #
        #############
        for timeind in range(numIter):
            # calculate reaction propensities
            a_i = self.__rxnpropensity()
            # calculate time to next reaction
            a_cumsum = np.cumsum(a_i)
            a_tot = a_cumsum[-1]
            tau = np.log(1/np.random.rand())/a_tot
            taus.append(tau)
            # choose next reaction
            rxnval = a_tot * np.random.rand()
            i = np.argmax(a_cumsum > rxnval)
            ivals.append(i)
            strainind = i//2
            # even rxns are growth, odd are death
            popchange = (-1)**i
            # propagate populations
            dIndex = int(np.ceil(tau/dt))
            trace[:, index:index+dIndex+1] = np.broadcast_to(trace[:, index],
                (dIndex+1, numStrains)).T
            trace[strainind, index+dIndex+1] += popchange
            # update T_elapsed and index
            T_elapsed += tau
            index += dIndex

        # update current population size for each strain and package traces
        for i, strain in enumerate(self.smallStrains):
            strain.N_pop = trace[i, index]
            self.poptrace[strain.ID][startind:startind+index+1] = trace[i,:index+1]
        # return T_elapsed and current index
        return T_elapsed, index, ivals

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
    def __rxnpropensity(self):
        '''A private method used to calculate the reaction propensities
        of the stochastic regime strains.

        Parameters
        ----------
        None

        Returns
        -------
        np.array a_i: numpy array of reaction propensities
        '''
        # declare array for reaction propensities
        a_i = np.zeros(2*numStrains)
        # loop through strains and calculate the birth/death propensities
        for i, strain in enumerate(self.smallStrains):
            # NOTE these are temporary reaction propensities for debug
            growprop = strain.k_grow
            deathprop = strain.k_death
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
    sim = thunderflask()
    sim.simulate(T=1000, dT=600, timestep=0.1)
    T_elapsed, index, ivals = sim.stochSim(600,0.1,0)
    print('T_elapsed = {0}; index = {1}'.format(T_elapsed, index))
    print('i values: {0}'.format(ivals[:5]))
    trace = sim.poptrace[sim.smallStrains[0].ID]
    plt.plot(trace[0:index])
    plt.show()
