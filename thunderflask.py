#mport dependencies
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
        self.estStrains = []
        self.f_trace = {
            'timepoints' : [],
            'fitnesses' : []
        }
        # partition initial strain list appropriately
        self.strainShuffle(T_curr=0, f_avg=0)

    def simulate(self, T=500, dt=1, T_0=0, mut_param=[1,2],
            save_established=False, save_dead=False):
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
        - bool save_established: tells the simulation whether or not to store
            established species
        - bool save_dead: tells the simulation whether or not to store dead
            species

        Returns
        -------
        None
        '''
        #initialize current time variable
        T_curr = T_0
        # loop until time is reached
        # while T_curr < T:
        for i in tqdm(range(int(T/dt))):
            # update average fitness and store
            f_avg = self.updateF_avg()
            self.f_trace['timepoints'].append(T_curr)
            self.f_trace['fitnesses'].append(f_avg)
            # run stochastic simulation
            T_next, taus = self.stochSim(dt, T_curr, f_avg)
            # run numerical simulation
            self.analyticSim(T_curr, dt, taus, f_avg)
            # run mutation simulation
            self.mutationSim(dt, T_next, mut_param)
            # shuffle strains
            self.strainShuffle(T_next, f_avg, save_established, save_dead)
            # update current time
            T_curr = T_next
        # return when completed
        return

    def stochSim(self, T_approx, T_curr, f_avg):
        '''
        Method used to perform stochastic simulation of bacterial growth
        dynamics. Used to handle strains whose total population size is less
        than a given threshold.

        Parameters
        ----------
        - float T_approx: the approximate simulation time (in generations)
        - float T_curr: the current simulation time (in generations)
        - float f_avg: the current average fitness in the population

        Returns
        ------
        - float T_curr: time at end of simulation (in generations)
        - np.array taus: the time steps between reactions (in generations)
        '''
        ######################################################
        # Prune strains to remove any less fit than the mean #
        ######################################################
        self.smallStrains = [bact for bact in self.smallStrains
            if bact.fitness >= f_avg]
        ###############################################
        # Initialize Stochastic Simulation attributes #
        ###############################################
        # Declare number of strains
        numStrains = len(self.smallStrains)
        # calculate a_tot for initial conditions to predict average timestep
        a_i = self.__rxnpropensity(f_avg)
        a_tot = a_i.sum()
        # declare number of iterations to perform; expect pop to double in 1 gen
        numIter = int(2*np.ceil(T_approx*a_tot)) # dt ~ 1/atot --> iter = T/dt
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
        for i, bacteria in enumerate(self.smallStrains):
            trace[i, 0] = bacteria.N_pop
        # declare iteration counter variables
        T_elapsed = 0
        #############
        # Main Loop #
        #############
        for timeind in range(1, numIter):
            # calculate reaction propensities
            a_i = self.__rxnpropensity(f_avg, trace[:, timeind-1])
            # calculate time to next reaction (break if rxnprop go to 0)
            a_cumsum = np.cumsum(a_i)
            a_tot = a_cumsum[-1]
            if a_tot <= 0:
                timeind -= 1
                break
            tau = np.log(1/np.random.rand())/a_tot
            # update T_elapsed
            T_elapsed += tau
            # break simulation if too much time has elapsed
            if T_elapsed > T_approx:
                timeind -= 1
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
        for i, bacteria in enumerate(self.smallStrains):
            bacteria.N_pop = trace[i, timeind]
            bacteria.timepoints += (np.cumsum(taus[:timeind+1]) +
                T_curr).tolist()
            bacteria.poptrace += trace[i,:timeind+1].tolist()
        # return T_elapsed and time intervals
        T_curr = T_curr + T_elapsed
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

    def strainShuffle(self, T_curr, f_avg, twiddle=2,
            min_threshold=1e2, max_threshold=1e4,
            save_established=False, save_dead=False):
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
            # note: commented block allows for full stochastic simulation, start to finish. However in the current iteration, strains are pruned when they have below average fitness. Thus traces are truncated. With this code commented out, and the if statement below changed to an elif statement, the full simulation can be run

            # # calculate 1/(f-<f>)
            # f = bacteria.fitness
            # if f - f_avg <= 0.0:
            #     threshold = np.infty
            # else:
            #     threshold = twiddle/(f - f_avg)
            # # cap threshold between minimum and maximum allowed values
            # threshold = max(min_threshold, threshold)
            # threshold = min(max_threshold, threshold)
            # # move strain to smallStrains if below the threshold
            # if bacteria.N_pop <= threshold:
            #     self.smallStrains.append(bacteria)
            #     big_toRemove.append(i)
            # move dead strains to deadStrains
            if bacteria.N_pop < 1:
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

    def mutationSim(self, dt, T_curr, mut_param):
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
            self.__mutate(bacteria, dfs, T_curr)

        # return from function
        return

    def updateF_avg(self):
        '''A method used to calculate the current average fitness of the system. Weights the average by the population size of each strain.

        Parameters
        ----------
        None

        Returns
        -------
        float f_avg: updated average fitness of the population
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
        # calculate average and return
        f_avg = f_weighted.sum()/pops.sum()
        return f_avg

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
        # declare array for reaction propensities
        a_i = np.zeros(2*len(self.smallStrains))
        # if populations is None, generate from sim object
        if len(populations) == 0:
            populations = np.zeros(len(self.smallStrains))
            for i, bacteria in enumerate(self.smallStrains):
                populations[i] = bacteria.N_pop
        # loop through strains and calculate the birth/death propensities
        for i, bacteria in enumerate(self.smallStrains):
            f = bacteria.fitness
            growprop = 1 + (f - f_avg)
            deathprop = 1
            a_i[2*i] = growprop * populations[i] # growth propensity
            a_i[2*i+1] = deathprop * populations[i] # death propensity
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

    def __mutate(self, bacteria, dfs, T_curr):
        ''' A private method used to generate new strains and package them, given an ancestral strain and a vector of mutation effects.

        Parameters
        ----------
        - bacteria.strain bacteria: bacterial strain to mutate
        - np.array dfs: the strengths of the mutations to apply
        - float T_curr: current time in the simulation

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

        # return from method
        return

# debug script
if __name__ == '__main__':
    # generate a single strain to begin with
    LUCA = strain(N_pop = 1e6)
    # create simulator object
    sim = thunderflask(LUCA)
    # artificially increase LUCA's fitness by setting f_avg = 0.1
    f_avg = 1.0
    # lets do some mutations!!
    mut_sim = [1, 1]
    sim.mutationSim(dt=1, T_curr=0, mut_param=mut_sim)
    # # begin debugging
    # run a round of stochastic simulation on these new mutants
    runtime = 50 # in generations
    T_curr = 0
    sim.stochSim(runtime, T_curr, f_avg)

    # # debug zombie strains
    # bonk = strain(N_pop=1, fitness=0)
    # sim.smallStrains.append(bonk)
    # runtime = 5 # in generations
    # T_curr = 0
    # sim.stochSim(runtime, T_curr, f_avg)
    # plot some of the trajectories of the stochastic strains
    for i in range(10):
        bact = sim.smallStrains[i]
        t = bact.timepoints
        pop = bact.poptrace
        plt.plot(t, pop)
    plt.show()
    # # test analytic component
    # taus = np.random.rand(1000)*0.01
    # sim.analyticSim(T_curr=0, taus=taus, f_avg=f_avg)
