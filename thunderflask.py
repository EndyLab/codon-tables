#import dependencies
import numpy as np
from scipy.stats import halfgennorm
import matplotlib.pyplot as plt
import random
from tqdm import tqdm
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
        self.allStrains = strains
        self.estStrains = []
        self.f_avgtrace = {
            'timepoints' : [],
            'f_avg' : []
        }
        self.f_trace = {
            'timepoints' : [],
            'fitnesses' : [],
        }
        self.poptrace = {}
        # partition initial strain list appropriately
        self.strainShuffle(T_curr=0, f_avg=0, save_established=True)
        # populate poptrace key/value pairs
        for bact in strains:
            self.poptrace[bact.ID] = ([], [])

    def simulate(self, T=500, dt=1, T_0=0, mut_param=[1,2], twiddle=3,
            save_established=False, save_dead=False, save_all=False,
            prune_strains=True):
        '''
        Main method of thunderflask class. Used to run genetic diversity
        simulation. There are four main modules in this simulation:

        (1) Stochastic -> (2) Numerical -> (3) Pop. Management -> (4) Mutations

        Consult documentation of each method for specific information. In brief,
        small and large population strains are handled separately in (1) and
        (2), and new mutations are generated in (4).

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
        - float T: the total time over which to simulate (in generations)
        - float dt: the time epoch over which to run each epoch (in generations)
        - float T_0: the initial time for the simulation (in generations)
        - list<float> mut_param: a list of parameters to pass to mutation module
        - float twiddle: a user defined number to adjust the thresholding
        - bool save_established: tells the simulation whether or not to store
            established species
        - bool save_dead: tells the simulation whether or not to store dead
            species
        - bool save_all: tells simulation whether or not to to save all species
        - bool prune_strains: tells the simulation whether or not to prune
            small strains below mean fitness

        Returns
        -------
        None
        '''
        #initialize current time variable
        T_curr = T_0
        # loop until time is reached
        # while T_curr < T:
        for i in tqdm(range(int(T/dt))):
            # update average fitness
            f_avg, fs = self.updateF_avg()
            # run stochastic simulation
            T_next, __ = self.stochSim(dt, T_curr, f_avg)
            # run numerical simulation
            res = 10
            taus = np.ones(res)/res
            self.analyticSim(T_curr=T_curr, dt=dt, taus=taus, f_avg=f_avg)
            # run mutation simulation
            self.mutationSim(T_curr=T_next, dt=dt, mut_param=mut_param,
                save_all=save_all)
            # shuffle strains
            self.strainShuffle(T_curr=T_next, f_avg=f_avg, twiddle=twiddle,
                save_established=save_established, save_dead=save_dead)
            # update current time
            T_curr = T_next
            # update traces
            self.tracer(T_curr, f_avg, fs)
        # return when completed
        return

    def stochSim(self, T_approx, T_curr, f_avg, prune_strains=True):
        '''
        Method used to perform stochastic simulation of bacterial growth
        dynamics. Used to handle strains whose total population size is less
        than a given threshold.

        Parameters
        ----------
        - float T_approx: the approximate simulation time (in generations)
        - float T_curr: the current simulation time (in generations)
        - float f_avg: the current average fitness in the population
        - bool prune_strains: tells the simulation whether or not to prune
            small strains below mean fitness

        Returns
        ------
        - float T_curr: time at end of simulation (in generations)
        - np.array taus: the time steps between reactions (in generations)
        '''
        # Optionally Prune strains to remove any less fit than the mean
        if prune_strains:
            self.smallStrains = [bact for bact in self.smallStrains
                if bact.fitness >= f_avg]
        ###############################################
        # Initialize Stochastic Simulation attributes #
        ###############################################
        # Declare number of strains
        numStrains = len(self.smallStrains)
        # Exit sim if the number of strains is zero
        if numStrains == 0:
            return T_curr, []
        # construct dictionary of possible reactions
        rxndict = self.__rxndict(numStrains)
        #####################
        # Declare variables #
        #####################
        # declare numpy array of traces for this epoch
        trace = np.zeros(numStrains)
        # declare list for storing time steps
        taus = []
        # declare initial population sizes
        for i, bacteria in enumerate(self.smallStrains):
            trace[i] = bacteria.N_pop
        # declare end condition
        T_end = T_curr + T_approx
        #############
        # Main Loop #
        #############
        while T_curr < T_end:
            # calculate reaction propensities
            a_i = self.__rxnpropensity(f_avg, trace)
            # break
            # calculate time to next reaction (break if rxnprop go to 0)
            a_cumsum = np.cumsum(a_i)
            a_tot = a_cumsum[-1]
            if a_tot <= 0:
                timeind -= 1
                break
            tau = np.log(1/np.random.rand())/a_tot
            # update T_curr
            T_curr += tau
            # append time interval
            taus.append(tau)
            # choose next reaction
            rxnval = a_tot * np.random.rand()
            i = np.argmax(a_cumsum > rxnval)
            # update population sizes
            trace += rxndict[i]
            # update population trace for each strain
            ind = int(i/2)
            self.smallStrains[ind].timepoints.append(T_curr)
            self.smallStrains[ind].poptrace.append(trace[ind])

        # update population size for each strain
        for i, bacteria in enumerate(self.smallStrains):
            bacteria.N_pop = trace[i]
        # return T_elapsed and time intervals
        return T_curr, taus

    def analyticSim(self, T_curr, dt, taus, f_avg):
        ''' Method used to simulate large population strains analytically. Loosely follows the design outlined in Desai and Fisher 2007.
        Parameters
        ----------
        - float T_curr: the current simulation time (in generations)
        - float dt: the time over which to simulate (in generations)
        - np.array taus: the time steps between reactions (in generations)
        - float f_avg: the current average fitness in the population

        Returns
        -------
        None
        '''
        # if taus is empty, generate an evenly spaced array of timepoints
        if len(taus) == 0:
            taus = np.linspace(0, dt, 10)
        # calculate array of timepoints
        t = np.cumsum(taus)
        t_array = t + T_curr
        # calculate population traces for each strain
        for bacteria in self.bigStrains:
            # unpack parameters
            N_0 = bacteria.N_pop
            f = bacteria.fitness - f_avg
            # calculate exponential growth of strain
            trace =  N_0 * np.exp(f*t)
            trace[trace < 0] = 0
            # update strain attributes
            bacteria.N_pop = trace[-1]
            bacteria.timepoints += t_array.tolist()
            bacteria.poptrace += trace.tolist()

        # return from the method
        return

    def strainShuffle(self, T_curr, f_avg, twiddle=3,
            min_threshold=1e2, max_threshold=1e4,
            save_established=False, save_dead=False,
            prune_strains=True):
        ''' A method used to handle exchanging strains between small and large
        population groups. Has one 'magic number' parameter to allow the user
        to alter the establishment threshold. Enforces a minimum threshold for
        analytic simulation as well.

        Parameters
        ----------
        - float T_curr: current time in the simulation
        - float threshold: the threshold population number that differentiates
            small and large population strains
        - float f_avg: the current average fitness in the population
        - float twiddle: a user defined number to adjust the thresholding
        - float min_threshold: the minimum population allowed for thresholding
        - float max_threshold: the maximum population allowed for thresholding
        - bool save_established: tells the simulation whether or not to store
            established species
        - bool save_dead: tells the simulation whether or not to store dead
            species
        - bool prune_strains: tells the simulation whether or not to prune
            small strains below mean fitness

        Returns
        -------
        None
        '''
        # loop through small strains
        small_toRemove = []
        for i, bacteria in enumerate(self.smallStrains):
            # calculate 1/(f-<f>)
            f = bacteria.fitness
            if f - f_avg <= 0.0:
                threshold = np.infty
            else:
                threshold = twiddle/(f - f_avg)
            # cap threshold between minimum and maximum allowed values
            threshold = max(min_threshold, threshold)
            threshold = min(max_threshold, threshold)
            # move strain to bigStrains if above the threshold
            if bacteria.N_pop > threshold:
                bacteria.t_est = T_curr
                self.bigStrains.append(bacteria)
                small_toRemove.append(i)
                # save established strains if requested
                if save_established:
                    self.estStrains.append(bacteria)
            # kill dead strains
            elif bacteria.N_pop < 1:
                small_toRemove.append(i)
                # save dead strains if requested
                if save_dead:
                    self.deadStrains.append(bacteria)
        # remove small strains that are large or dead
        smallStrains = [bact for i, bact in enumerate(self.smallStrains)
            if i not in small_toRemove]
        self.smallStrains = smallStrains

        # loop through large strains
        big_toRemove = []
        for i, bacteria in enumerate(self.bigStrains):
            # check thresholding if prune_strains=True
            if not prune_strains:
                # calculate 1/(f-<f>)
                f = bacteria.fitness
                if f - f_avg <= 0.0:
                    threshold = np.infty
                else:
                    threshold = twiddle/(f - f_avg)
                # cap threshold between minimum and maximum allowed values
                threshold = max(min_threshold, threshold)
                threshold = min(max_threshold, threshold)
                # move strain to smallStrains if below the threshold
                if bacteria.N_pop <= threshold:
                    self.smallStrains.append(bacteria)
                    big_toRemove.append(i)
                # check for death
                elif bacteria.N_pop < 1:
                    # kill strain
                    big_toRemove.append(i)
                    # save dead strains if requested
                    if save_dead:
                        self.deadStrains.append(bacteria)
            # if not prune strains, just check for death
            elif bacteria.N_pop < 1:
                # kill strain
                big_toRemove.append(i)
                # save dead strains if requested
                if save_dead:
                    self.deadStrains.append(bacteria)
        # remove small strains that are large or dead
        bigStrains = [bact for i, bact in enumerate(self.bigStrains)
            if i not in big_toRemove]
        self.bigStrains = bigStrains

        # return from method
        return

    def mutationSim(self, T_curr, dt, mut_param, save_all=False):
        '''
        Method used to determine the number of new strains to generate in a
        given period of time. Also handles generation and storage of these new
        strains via helper function(s). Relies on helper functions __mutNum()
        to choose the number of mutants to generate, and __mutStrength() to
        choose the strength of those mutations.

        Parameters
        ----------
        - float dt: the time over which to generate mutants (in generations)
        - float T_curr: current time in the simulation
        - list<float> mut_param: a list of floats to pass to __mutStrength()
        - bool save_all: tells simulation to save in self.allStrains

        Returns
        -------
        None
        '''
        # loop through established strains
        for bacteria in self.bigStrains:
            # randomly draw the number of mutants
            n_mut = self.__mutNum(bacteria, dt)
            # generate a vector of mutation effects using __mutStrength()
            dfs = self.__mutStrength(n_mut, mut_param)
            # generate new strains from ancestor using __mutate()
            self.__mutate(bacteria, dfs, T_curr, save_all)

        # return from function
        return

    def updateF_avg(self):
        '''A method used to calculate the current average fitness of the
        system. Weights the average by the population size of each strain. Also
        returns the fitnesses of every living strain as an array

        Parameters
        ----------
        None

        Returns
        -------
        float f_avg: updated average fitness of the population
        np.array fs: array of fitnesses
        '''
        # declare lists of population sizes and fitnesses
        pops = []
        fs = []
        # populate lists
        for bact in self.smallStrains + self.bigStrains:
            pops.append(bact.N_pop)
            fs.append(bact.fitness)
        # convert to numpy arrays and calculate fitness*population
        pops = np.array(pops)
        fs = np.array(fs)
        f_weighted = fs*pops
        # calculate average and return appropriate variables
        f_avg = f_weighted.sum()/pops.sum()
        return f_avg, fs

    def tracer(self, T_curr, f_avg, fs):
        '''A method used to update appropriate metric traces

        Parameters
        ----------
        - float T_curr: the current simulation time (in generations)
        - float f_avg: the current average fitness in the population
        - np.array fs: array of fitnesses

        Returns
        -------
        None
        '''
        # update f_avgtrace
        self.f_avgtrace['timepoints'].append(T_curr)
        self.f_avgtrace['f_avg'].append(f_avg)
        # update poptrace
        for bact in self.smallStrains + self.bigStrains:
            self.poptrace[bact.ID][0].append(T_curr)
            self.poptrace[bact.ID][1].append(bact.N_pop)
        # # update f_trace
        # self.f_trace['timepoints'].append(T_curr)
        # self.f_trace['fitnesses'].append(fs)
        # return from function
        return

    #######################
    #   Private Methods   #
    #######################
    def __rxnpropensity(self, f_avg, populations=[]):
        '''A private method used to calculate the reaction propensities
        of the stochastic regime strains.

        Parameters
        ----------
        - float f_avg: the current average fitness in the population
        - np.array populations: array of populations sizes for each strain; if
            empty, pulls population size from sim.smallStrains[i].N_pop

        Returns
        -------
        np.array a_i: numpy array of reaction propensities
        '''
        # declare array for unweighted propensities
        props = np.ones(2*len(self.smallStrains))
        # if populations is None, generate from sim object
        if len(populations) == 0:
            populations = np.zeros(len(self.smallStrains))
            for i, bacteria in enumerate(self.smallStrains):
                populations[i] = bacteria.N_pop
        # copy population sizes (to have two elements per strain)
        populations = np.repeat(populations, 2)
        # loop through strains and calculate the birth propensities
        for i, bacteria in enumerate(self.smallStrains):
            # update birth propensity (death is always 1)
            f = bacteria.fitness
            props[2*i] = 1 + (f - f_avg)
        # return results
        return props*populations

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

    def __mutNum(self, bacteria, dt):
        ''' A private method used to generate a number of mutations to
        introduce in this generation. Assumes the generation of novel mutations
        follows a Poisson Process whose expectation value is:

                l_mut = (N*Ub)dt

        Where N is the current strain population size, Ub is the per
        genome-generation beneficial mutation rate, and dt is the time period
        over which to generate mutants (in generations).

        Parameters
        ----------
        - bacteria.strain bacteria: bacterial strain to mutate
        - float dt: time period over which to mutate (in generations)

        Returns
        -------
        int n_mut: number of mutants to introduce this generation
        '''
        # extract strain specific values and calculate expectation value
        N = bacteria.N_pop
        Ub = bacteria.mu
        l_mut = N*Ub*dt
        # randomly draw the number of mutants from a Poisson distribution
        n_mut = np.random.poisson(l_mut)
        return n_mut

    def __mutStrength(self, n_mut, mut_param):
        ''' A private method used to generate a vector of mutation effects
        given a number of mutants to generate. Draws these mutation effects
        from a generalized halfnormal distribution of the form:

                         beta*lam     -(lam*s)^beta
                P(s) = ------------- e
                       Gamma(1/beta)

        where s is the fitness strength, b is a shape parameter and l is a
        scaling parameter. When b = 1, the distribution is exponential. When b
        = 2, the distribution is normal. Automatically scales mutation strength
        to a percentage (divides output by 100)

        Parameters
        ----------
        - int n_mut: number of mutants to introduce this generation
        - list<float> mut_param: the shape and scaling parameters for the halfgennorm
            distribution

        Returns
        -------
        np.array dfs: the strengths of the mutations generated
        '''
        # unpack mutation distribution parameters
        beta, lam = mut_param
        # make n_mut draws from a halfgennorm distribution
        dfs = halfgennorm.rvs(beta, scale=lam, size=n_mut)/100
        return dfs

    def __mutate(self, bacteria, dfs, T_curr, save_all=False):
        ''' A private method used to generate new strains and package them,
        given an ancestral strain and a vector of mutation effects. Also
        creates a key/value pair in self.poptrace dict. Optionally saves to all
        strains

        Parameters
        ----------
        - bacteria.strain bacteria: bacterial strain to mutate
        - np.array dfs: the strengths of the mutations to apply
        - float T_curr: current time in the simulation
        - bool save_all: tells simulation to save in self.allStrains

        Returns
        -------
        None
        '''
        # extract cell lineage and fitness from ancestor strain
        lineage = bacteria.lineage
        f = bacteria.fitness
        # loop through mutation strengths
        for df in dfs:
            # generate a new strain
            mutant = strain(N_pop=1, table=bacteria.table,
                t_0=T_curr, fitness=f+df, lineage=lineage)
            # package resulting strain back into simulation
            self.smallStrains.append(mutant)
            if save_all:
                self.allStrains.append(mutant)
            # create a key value pair in self.poptrace
            self.poptrace[mutant.ID] = ([], [])

        # return from method
        return

