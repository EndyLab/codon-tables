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

    def simulate(self):
        '''
        Main method of thunderflask class. Used to run genetic diversity
        simulation. There are four main modules in this simulation:

        (1) Mutations-> (2) Numerical -> (3) Stochastic -> (4) Pop. Management

        Consult documentation of each method for specific information. In brief,
        small and large population strains are handled separately in (2) and
        (3), and new mutations are generated in (1).

        Loosely follows the design outlined in Desai and Fisher 2007. The model
        makes the following assumptions
        - (on average) fixed population size N
        - only beneficial mutations can establish (survive drift) and
            contribute significantly to genetic diversity
        - fixed beneficial mutation rate Ub; number of mutations drawn from
            Poisson distribution
        - does NOT assume mutations all have the same effect; drawn from
            exponential distribution
        - assumes populations reaching size of at least 3/[f - <f>] have
            established

        Parameters
        ----------
        None

        Returns
        -------
        None
        '''
        # main loop
        # while not done:
            # numerical simulation

            # stochastic simulation

            # strain swapping (if low pop strains get large or vice versa)

            # mutation simulation

        # return results of simulation
        return

    def mutationSim(self):
        '''
        '''
        return

    def stochSim(self, T_approx, T_0):
        '''
        Method used to perform stochastic simulation of bacterial growth
        dynamics. Used to handle strains whose total population size is less
        than a given threshold.

        Parameters
        ----------
        - float T_approx: the approximate simulation time (in generations)
        - float T_0: the current simulation time (in generations)

        Returns
        ------
        - float T_elapsed: the actual simulation time (in generations)
        - np.array taus: the time steps between reactions (in generations)
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
            # update T_elapsed
            T_elapsed += tau
            # break simulation if too much time has elapsed
            if T_elapsed > T_approx:
                T_elapsed -= tau
                break
            # otherwise, update simulation
            else:
                # append time interval
                taus[timeind] = tau
                # choose next reaction
                rxnval = a_tot * np.random.rand()
                i = np.argmax(a_cumsum > rxnval)
                # update population sizes
                trace[:,timeind] = trace[:,timeind-1] + rxndict[i]

        # update population size and timepoints for each strain
        for i, strain in enumerate(self.smallStrains):
            strain.N_pop = trace[i, -1]
            strain.timepoints += (np.cumsum(taus) + T_0).tolist()
            strain.poptrace += trace[i,:].tolist()
        # return T_elapsed and time intervals
        return T_elapsed, taus

        def analyticSim(self, T_0, taus):
            ''' Method used to simulate large population strains analytically. Loosely follows the design outlined in Desai and Fisher 2007.
            Parameters
            ----------
            - float T_0: the current simulation time (in generations)
            - np.array taus: the time steps between reactions (in generations)

            Returns
            -------
            None
            '''
            # calculate array of timepoints
            t_array = np.cumsum(taus) + T_0
            # calculate population traces for each strain
            for strain in self.bigStrains:
                # unpack parameters
                N_0 = strain.N_pop
                pd, pm, kd, km, td, tm = strain.growParam
                t_offset = strain.t_large
                # adjust time parameters by time offset
                td -= t_offset
                tm -= t_offset
                t_adjusted = t_array - t_offset
                # simulate
                # update strain attributes
                strain.N_pop = trace[-1]
                strain.timepoints += t_array.tolist()
                strain.poptrace += trace.tolist()

            # return from the method
            return

    def strainShuffle(self, T_curr, threshold):
        ''' A method used handle exchanging strains between small and large
        population groups.

        Parameters
        ----------
        - float T_curr: current time in the simulation
        - float threshold: the threshold population number that differentiates
            small and large population strains

        Returns
        -------
        None
        '''
        # loop through small strains
        for i, strain in enumerate(self.smallStrains):
            # move strain to bigStrains if above the threshold
            if strain.N_pop > threshold:
                strain.t_large = T_curr
                strain.bigStrains.append(strain)
                __ = strain.smallStrains.pop(i)

        # loop through large strains
        for i, strain in enumerate(self.bigStrains):
            # move strain to smallStrains if below the threshold
            if strain.N_pop <= threshold:
                strain.smallStrains.append(strain)
                __ = strain.bigStrains.pop(i)

        # return from method
        return

    #######################
    #   Private Methods   #
    #######################
    def __rxnpropensity(self, f_avg):
        '''A private method used to calculate the reaction propensities
        of the stochastic regime strains.

        Parameters
        ----------
        float f_avg: the current average fitness in the population

        Returns
        -------
        np.array a_i: numpy array of reaction propensities
        '''
        # declare array for reaction propensities
        a_i = np.zeros(2*len(self.smallStrains))
        # loop through strains and calculate the birth/death propensities
        for i, strain in enumerate(self.smallStrains):
            f = strain.fitness
            growprop = 1 + (f - f_avg)
            deathprop = 1
            a_i[2*i] = growprop*strain.N_pop # growth propensity
            a_i[2*i+1] = deathprop*strain.N_pop # death propensity
        # return results
        return a_i

    @staticmethod
    def __rxndict(numStrains):
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

    @staticmethod
    def __horowitzContinuous(N_0, param, t_array):
        ''' A static method used to calculate the growth rate of strains in the
        continuous regime. Uses the birth and death probabilities from Horowitz
        et al. 2010 and the closed form integration of these probabilities.

        Parameters
        ----------
        - float N_0: strain size at initial time
        - list param: strain specific parameters for simulation
        - np.array t_array: array of timepoints over which to simulate

        Returns
        -------
        np.array trace: array representing the population level at given times
        '''
        # unpack simulation parameters
        pd, pm, kd, km, td, tm = param
        # calculate curve and return population trace
        trace = ( N_0 * np.exp((pd - pm)*t_array)
            * ( (1 + np.exp(kd*(td - t_array)))**(pd/kd)
            / (1 + np.exp(km*(tm - t_array)))**(pm/km)) )
        return trace

# debug script
if __name__ == '__main__':
    # generate some strains
    strains = []
    for i in range(10):
        strains.append(strain())
    sim = thunderflask(strains)
    sim.bigStrains = sim.smallStrains
    taus = np.ones(5000) * 0.1
    sim.analyticSim(10, taus)
    for i in range(5):
        t = sim.bigStrains[i].timepoints
        trace = sim.bigStrains[i].poptrace
        plt.plot(t, trace, label='Strain {0}'.format(i))
    print(t[-1])
    plt.legend()
    plt.show()
