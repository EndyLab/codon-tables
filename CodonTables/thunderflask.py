#import dependencies
import numpy as np
from copy import deepcopy as copy
from scipy.stats import halfgennorm, binom
import matplotlib.pyplot as plt
import random
from CodonTables.table import CodonTable
from CodonTables.utils import utils
from CodonTables.bacteria import Strain

# define simulator class
class Thunderflask():
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
        self.big_strains = []
        self.small_strains = strains
        self.dead_strains = []
        self.all_strains = strains
        self.est_strains = []
        self.f_avgtrace = {
            'timepoints' : [],
            'f_avg' : []
        }
        self.f_trace = {
            'timepoints' : [],
            'fitnesses' : [],
        }
        self.poptrace = {}
        self.popfrac = {}
        # partition initial strain list appropriately
        self.strain_shuffle(T_curr=0, f_avg=0, save_dead=False, save_established=True,
                           prune_strains=False)
        # populate poptrace key/value pairs
        for bact in strains:
            self.poptrace[bact.ID] = ([], [])

        # initialize population fraction dictionary
        codes = set([bact.code for bact in strains])
        for code in codes: self.popfrac[code] = []

    def simulate(self, T=500, dt=1, T_0=0, mut_param=[1,2], twiddle=3,
            save_established=False, save_dead=False, save_all=False,
            prune_strains=True, show_progress=True,
            competition=False):
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
        - bool competition: tells simulation to prematurely terminate if
            population fraction becomes 1 for any given genetic code

        Returns
        -------
        None
        '''
        #initialize current time variable
        T_curr = T_0
        # update average fitness and store initial conditions
        f_avg, fs = self.update_f_avg()
        self.tracer(T_curr=T_curr, f_avg=f_avg, fs=fs)
        # loop until time is reached, either with or without tqdm
        if show_progress == True or type(show_progress) == str:
            # import either tqdm or tqdm_notebook depending on input
            if show_progress == 'notebook':
                from tqdm import tqdm_notebook as tqdm
            else:
                from tqdm import tqdm
            # run simulation
            for i in tqdm(range(int(T/dt)), desc='Iteration Number: '):
                T_curr = self.iterate(T_curr=T_curr, dt=dt,
                                      mut_param=mut_param, twiddle=twiddle,
                                      save_established=save_established,
                                      save_dead=save_dead, save_all=save_all,
                                      prune_strains=prune_strains)
                # break if popfraction is met
                if self.__check_domination(competition):
                    print('Population Sweep at time t={0}'.format(T_curr))
                    break
        else:
            while T_curr < T:
                T_curr = self.iterate(T_curr=T_curr, dt=dt,
                                      mut_param=mut_param, twiddle=twiddle,
                                      save_established=save_established,
                                      save_dead=save_dead, save_all=save_all,
                                      prune_strains=prune_strains)
                # break if popfraction is met
                if self.__check_domination(competition):
                    print('Population Sweep at time t={0}'.format(T_curr))
                    break
        return

    def iterate(self, T_curr, dt, mut_param, twiddle,
                save_established, save_dead, save_all, prune_strains):
        '''
        Method used to perform one iteration of the simulation. This allows
        simulate() to optionally give tqdm progress bar.

        Parameters
        ----------
        - float T_curr: the current simulation time (in generations)
        - float dt: the time epoch over which to run each epoch (in generations)
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
        float T_curr: the current simulation time (in generations) after one
            iteration
        '''
        # update average fitness
        f_avg, fs = self.update_f_avg()
        # run stochastic simulation
        T_next, __ = self.stoch_sim(T_approx=dt, T_curr=T_curr, f_avg=f_avg,
                                   prune_strains=prune_strains)
        # handle when no stoch sim occured for numerical simulation (set to dt)
        if T_curr == T_next:
            T_next += dt
        # run numerical simulation
        res = 10
        taus = np.ones(res)*(T_next - T_curr)/res # generates res timepoints
        self.analytic_sim(T_curr=T_curr, dt=dt, taus=taus, f_avg=f_avg)
        # run mutation simulation
        self.mutation_sim(T_curr=T_next, dt=dt, mut_param=mut_param, save_all=save_all)
        # shuffle strains
        self.strain_shuffle(T_curr=T_next, f_avg=f_avg, twiddle=twiddle,
                           save_established=save_established, save_dead=save_dead,
                           prune_strains=prune_strains)
        # update current time
        T_curr = T_next
        # update traces
        self.tracer(T_curr=T_curr, f_avg=f_avg, fs=fs)
        # return current time
        return T_curr

    def stoch_sim(self, T_approx, T_curr, f_avg, prune_strains):
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
            self.small_strains = [bact for bact in self.small_strains
                if bact.fitness >= f_avg]
        ###############################################
        # Initialize Stochastic Simulation attributes #
        ###############################################
        # Declare number of strains
        num_strains = len(self.small_strains)
        # Exit sim if the number of strains is zero
        if num_strains == 0:
            return T_curr, []
        # construct dictionary of possible reactions
        rxndict = self.__rxndict(num_strains)
        #####################
        # Declare variables #
        #####################
        # declare numpy array of traces for this epoch
        trace = np.zeros(num_strains)
        # declare list for storing time steps
        taus = []
        # declare initial population sizes
        for i, bacteria in enumerate(self.small_strains):
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
            self.small_strains[ind].timepoints.append(T_curr)
            self.small_strains[ind].poptrace.append(trace[ind])

        # update population size for each strain
        for i, bacteria in enumerate(self.small_strains):
            bacteria.N_pop = trace[i]
        # return T_elapsed and time intervals
        return T_curr, taus

    def analytic_sim(self, T_curr, dt, taus, f_avg):
        ''' Method used to simulate large population strains analytically.
        Loosely follows the design outlined in Desai and Fischer 2007.

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
        for bacteria in self.big_strains:
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

    def strain_shuffle(self, T_curr, f_avg, save_established, save_dead, prune_strains,
                      twiddle=3,min_threshold=1e2, max_threshold=1e4):
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
        small_to_remove = []
        for i, bacteria in enumerate(self.small_strains):
            # calculate 1/(f-<f>)
            f = bacteria.fitness
            if f - f_avg <= 0.0:
                threshold = np.infty
            else:
                threshold = twiddle/(f - f_avg)
            # cap threshold between minimum and maximum allowed values
            threshold = max(min_threshold, threshold)
            threshold = min(max_threshold, threshold)
            # move strain to big_strains if above the threshold
            if bacteria.N_pop > threshold:
                bacteria.t_est = T_curr
                self.big_strains.append(bacteria)
                small_to_remove.append(i)
                # save established strains if requested
                if save_established:
                    self.est_strains.append(bacteria)
            # kill dead strains
            elif bacteria.N_pop < 1:
                small_to_remove.append(i)
                # save dead strains if requested
                if save_dead:
                    self.dead_strains.append(bacteria)
        # remove small strains that are large or dead
        small_strains = [bact for i, bact in enumerate(self.small_strains)
            if i not in small_to_remove]
        self.small_strains = small_strains

        # loop through large strains
        big_to_remove = []
        for i, bacteria in enumerate(self.big_strains):
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
                # move strain to small_strains if below the threshold
                if bacteria.N_pop <= threshold:
                    self.small_strains.append(bacteria)
                    big_to_remove.append(i)
                # check for death
                elif bacteria.N_pop < 1:
                    # kill strain
                    big_to_remove.append(i)
                    # save dead strains if requested
                    if save_dead:
                        self.dead_strains.append(bacteria)
            # if not prune strains, just check for death
            elif bacteria.N_pop < 1:
                # kill strain
                big_to_remove.append(i)
                # save dead strains if requested
                if save_dead:
                    self.dead_strains.append(bacteria)
        # remove small strains that are large or dead
        big_strains = [bact for i, bact in enumerate(self.big_strains)
            if i not in big_to_remove]
        self.big_strains = big_strains

        # return from method
        return

    def mutation_sim(self, T_curr, dt, mut_param, save_all):
        '''
        Method used to determine the number of new strains to generate in a
        given period of time. Also handles generation and storage of these new
        strains via helper function(s). Relies on helper functions __mut_num()
        to choose the number of mutants to generate, and __mut_strength() to
        choose the strength of those mutations.

        Parameters
        ----------
        - float dt: the time over which to generate mutants (in generations)
        - float T_curr: current time in the simulation
        - list<float> mut_param: a list of floats to pass to __mut_strength()
        - bool save_all: tells simulation to save in self.all_strains

        Returns
        -------
        None
        '''
        # loop through established strains
        for bacteria in self.big_strains:
            # randomly draw the number of mutants
            n_mut = self.__mut_num(bacteria, dt)
            # generate a vector of mutation effects using __mut_strength()
            dfs = self.__mut_strength(n_mut, mut_param)
            # generate new strains from ancestor using __mutate()
            self.__mutate(bacteria, dfs, T_curr, save_all)

        # return from function
        return

    def update_f_avg(self):
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
        for bact in self.small_strains + self.big_strains:
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
        # update poptrace and calculate population fractions
        N_tot = 0
        n_code = {};
        for code in self.popfrac.keys(): n_code[code] = 0
        for bact in self.small_strains + self.big_strains:
            # update poptrace
            self.poptrace[bact.ID][0].append(T_curr)
            self.poptrace[bact.ID][1].append(bact.N_pop)
            # update N_tot and n_code counters
            N_tot += bact.N_pop
            n_code[bact.code] += bact.N_pop
        # update popfrac
        for code, fractrace in self.popfrac.items():
            fractrace.append(n_code[code]/N_tot)
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
            empty, pulls population size from sim.small_strains[i].N_pop

        Returns
        -------
        np.array a_i: numpy array of reaction propensities
        '''
        # declare array for unweighted propensities
        props = np.ones(2*len(self.small_strains))
        # if populations is None, generate from sim object
        if len(populations) == 0:
            populations = np.zeros(len(self.small_strains))
            for i, bacteria in enumerate(self.small_strains):
                populations[i] = bacteria.N_pop
        # copy population sizes (to have two elements per strain)
        populations = np.repeat(populations, 2)
        # loop through strains and calculate the birth propensities
        for i, bacteria in enumerate(self.small_strains):
            # update birth propensity (death is always 1)
            f = bacteria.fitness
            props[2*i] = 1 + (f - f_avg)
        # return results
        return props*populations

    @staticmethod
    def __rxndict(num_strains):
        '''A private method used to generate a reaction dictionary for the
        stochastic simulator. Each item is a numpy array of length num_strains
        which describes the change in each species count

        Parameters
        ----------
        int num_strains: the number of strains to simulate

        Returns
        -------
        dict rxndict: python dict representing the possible reactions
        '''
        # declare rxndict to return
        rxndict = {}
        # generate two arrays for each strain (grow and death) and store
        for i in range(num_strains):
            growarray = np.zeros(num_strains)
            # set ith element to be 1
            growarray[i] = 1
            # even elements are growth, odd are death
            rxndict[2*i] = growarray
            rxndict[2*i+1] = -growarray
        # return resulting dict
        return rxndict

    @staticmethod
    def __horowitz_continuous(N_0, param, t_array):
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

    def __mut_num(self, bacteria, dt, threshold=1):
        ''' A private method used to generate a number of mutations to
        introduce in this generation. Assumes the generation of novel mutations
        follows a Poisson Process whose expectation value is:

                l_mut = (N*Ub)dt

        Where N is the current strain population size, Ub is the per
        genome-generation beneficial mutation rate, and dt is the time period
        over which to generate mutants (in generations). Automatically
        determines when the Poisson assumptions break down and switches to N
        independent draws from a binomial distribution in this case.

        Parameters
        ----------
        - bacteria.strain bacteria: bacterial strain to mutate
        - float dt: time period over which to mutate (in generations)
        - float threshold: the cutoff expectation value after which the
            simulation will switch to independent draws from binomial
            distribution

        Returns
        -------
        int n_mut: number of mutants to introduce this generation
        '''
        # extract strain specific values and calculate expectation value
        N = bacteria.N_pop
        Ub = bacteria.mu
        l_mut = N*Ub*dt
        # determine whether or not to use binomial or Poisson distribution
        if l_mut < threshold:
            # randomly draw from binomial distribution
            n_mut = np.random.binomial(N, Ub*dt)
        else:
            # randomly draw the number of mutants from a Poisson distribution
            n_mut = np.random.poisson(l_mut)
        return n_mut

    def __mut_strength(self, n_mut, mut_param):
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

    def __mutate(self, bacteria, dfs, T_curr, save_all):
        ''' A private method used to generate new strains and package them,
        given an ancestral strain and a vector of mutation effects. Also
        creates a key/value pair in self.poptrace dict. Optionally saves to all
        strains

        Parameters
        ----------
        - bacteria.strain bacteria: bacterial strain to mutate
        - np.array dfs: the strengths of the mutations to apply
        - float T_curr: current time in the simulation
        - bool save_all: tells simulation to save in self.all_strains

        Returns
        -------
        None
        '''
        # extract cell characteristics from ancestor strain
        lineage = bacteria.lineage
        f = bacteria.fitness
        code = bacteria.code
        table = bacteria.table
        # loop through mutation strengths
        for df in dfs:
            # generate a new strain
            mutant = strain(
                N_pop=1, code=code, table=table,
                t_0=T_curr, fitness=f+df, lineage=copy(lineage)
            )
            # package resulting strain back into simulation
            self.small_strains.append(mutant)
            if save_all:
                self.all_strains.append(mutant)
            # create a key value pair in self.poptrace
            self.poptrace[mutant.ID] = ([], [])

        # return from method
        return

    def __check_domination(self, competition):
        '''A private method used to check whether or not one genetic code has
        swept a population. Returns a boolean telling self.simulate whether or
        not to terminate

        Parameters
        ----------
        bool competition: tells simulation to prematurely terminate if
            population fraction becomes 1 for any given genetic code

        Returns
        -------
        bool decision: continuation condition for simulation
        '''
        # immediately declare continuation if competition == false
        if competition == False:
            decision = False
        # otherwise, check continuation condition
        else:
            # loop through popfractions to see if any one code dominates
            decision = 1 in [popfrac[-1] for popfrac in self.popfrac.values()]
        return decision
