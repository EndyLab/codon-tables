#import dependencies
import numpy as np
import matplotlib.pyplot as plt
from numba import jit
from codonTable import codonTable
from codonUtils import utils

#define Monte Carlo Simulation class
class MonteCarlo:
    '''MonteCarlo: A class designed to optimize a codon table given an
    arbitrary objective function to minimize/maximize
    '''
    def __init__(self, table=utils.standardTable,
                costfxn=None, wobble_rule = 'standard'):
        '''MonteCarlo.__init__: the init function for the MonteCarlo class.
        Optionally allows the user to specify the starting codon table,
        associated objective function, and which wobble_rules should be
        followed. Default values are supplied if the user does not specify
        them (i.e. Standard Code, 47 block wobble rule)

        Parameters
        ----------
        - dict table = utils.standardTable: a python dict representing the
            codon table. Also optionally accepts a codonTable object.
        - func costfxn = None: a function that takes a python dict as an
            input and outputs a cost associated with that table. If no
            funciton is supplied, defaults to maxMutMinVar.
        - str wobble_rule = 'standard': a string telling the simulator which
            wobble rules to follow for accepting new tablesi

            Acceptable inputs:
            - 'standard' : 47 blocks, 2 stop codons
            - 'preserveBlock' : maintain same block structure as standard
                table
            - 'unrestricted' : 63 open blocks, at least 0 of every AA and
                stop.
        Returns
        -------
        MonteCarlo obj: returns an instance of the MonteCarlo object
        '''
        # handle codonTable objects being passed
        if type(table) == codonTable:
            table = codonTable.codonTable
        # handle costfxn input
        if costfxn == None:
            costfxn = self.maxMutMinVar
        # calculate some attributes
        resiCounts = utils.getAAcounts(table)
        connectivity = utils.getCodonConnectivity(table)
        # assign attributes
        self.table = table
        self.costfxn = costfxn
        self.resiCounts = resiCounts
        self.connectivity = connectivity
        self.utils = utils
        self.wobble_rule = wobble_rule

    ######################################
    ##          Public  Methods         ##
    ######################################


    def run(self, objective = None, algorithm = None, minflag = False):
        '''MonteCarlo.run: the method used to handle searching through codon space using a defined search algorithm and objective function.

        Parameters
        ----------

        Returns
        -------
        '''
        # handle default selections
        if objective == None:
            objective = self.maxMutMinVar
        if algorithm == None:
            algorithm = self.GDA
        return False

    def GDA(self):
        #ToDo: fill out GDA
        return False


    def tableShuffle(self, table):
        '''Takes a codon table as an input and shuffles it to produce a new,
            similar table in order to traverse table space.

        Parameters
        ----------
        dict table: a python dict representing a codon table starting point

        Returns
        -------
        dict newTable: a python dict representing the next codon table
        '''
        
        return False



    @staticmethod
    def maxMutMinVar(table, counts, connectivity):
        '''MonteCarlo.maxMutMinVar: the default cost function for class.
        Implements a version of the cost function used by Novozhilov et al.
        2007, but optimizes for maximizing mutability while minimizing variance
        per mutation

        Parameters
        ----------
        - dict table: a python dictionary representing the codon table to
            evaluate
        - dict counts: a python dictionary representing the number of codons
            encoding for a particular amino acid in the given table
        - dict connectivity: a python dictionary mapping str codon --> list of
            (str codon, int dist) representing which codons can be reached from
            the input codon without tracing intermediate residues, along with
            the number of mutations to that codon

        Returns
        -------
        float metric: the absolute number score of the given table
        '''
        # initialize cost variable, declare polar requirement scale
        metric = 0
        PRS = utils.PRS
        # loop over source codons
        for c_1 in table.keys():
            # calculate AA degeneracy; skip stop codons
            AA_1 = table[c_1]
            if AA_1 == '*':         # skips stops
                continue
            f_c = 1/counts[AA_1]
            # loop over sink codons
            for (c_2, dist) in connectivity[c_1]:
                # get second AA, skip stop codons
                AA_2 = table[c_2]
                if AA_2 == '*':
                    continue
                # calculate components of objective function
                d_c12 = 1/(1 + (PRS[AA_1] - PRS[AA_2])**2)
                P_c12 = (1/12)**dist
                # calculate contribution to metric from this pair
                metric += f_c * P_c12 * d_c12
        # return resulting metric
        return metric

    ######################################
    ##          Private Methods         ##
    ######################################

# Debugging
sim = MonteCarlo()
table = sim.table
counts = sim.resiCounts
connect = sim.connectivity

sim.maxMutMinVar(table, counts, connect)
