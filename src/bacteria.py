#import dependencies
import numpy as np
import matplotlib.pyplot as plt
from uuid import uuid4 as uuid
from src.codonUtils import utils
from src.codonTable import codonTable

class strain():
    '''A class used to represent a bacterial strain in culture'''
    def __init__(self, code=None, N_pop=100, fitness=0, mu=2e-5, t_0=0,
        t_est=0, table=None, ID=None, lineage=[]):
        '''The init function for the strain class.

        Parameters
        ----------
        - str code: a human readable version of the genetic code used by strain
        - int/float N_pop: current population size of strain
        - float fitness: the absolute fitness of an individual strain
        - float mu: mutation rate of bacterial strain (1/genome-gen) assuming
            Standard Code (automatically adjusts based on genetic code)
        - float t_0: timestep in which a strain first arrives
        - float t_est: timestep in which a strain becomes large enough to
            simulate analytically (established)
        - dict table: genetic code for bacterial strain
        - str ID: a unique string identifier for bacterial strain
        - list <str> lineage: a list tracking the lineage of bacterial strain

        Returns
        -------
        strain obj: returns handle for strain instance
        '''
        # generate values for optional parameters
        if table == None:
            code = 'Standard Code'
            table = codonTable(utils.standardTable)
        elif type(table) == dict:
            table = codonTable(table=table)

        if ID == None:
            ID = str(uuid())
        # update lineage to include own ID
        lineage.append(ID)
        # store parameters in class attributes
        self.code = code
        self.N_pop = N_pop
        self.fitness = fitness
        self.t_0 = t_0
        self.t_est = t_est
        self.table = table
        self.mu = self.getMutRate(mu)
        self.ID = str(ID)
        self.lineage = lineage
        self.timepoints = []
        self.poptrace = []

    def getMutRate(self, mu):
        '''A private method used to calculate what the mutation rate should be for a particular strain given the rate for cells with the Standard Code.

        Parameters
        ----------
        float mu: the mutation rate of this strain given the Standard Code

        Returns
        -------
        float mu_adj: the adjusted mutation rate given the gentic code of the
            strain
        '''
        # get dict form of codon table and of standard code
        table = self.table.codonDict
        sc = utils.standardTable
        # if table is identical to the standard table, return mu
        if table == sc:
            mu_adj = mu
        # if not, calculate how many substitutions are available to each code
        else:
            # initialize counters for the standard code and strain code
            count_sc = 0
            count = 0
            # calculate mutational availability metric for standard code
            for (c1, c2) in utils.tripletMutPairs:
                # increment counters if substitution is nonsynonymous AND
                # neither codon encodes for a STOP
                count_sc += not (
                    (sc[c1] == '*' or sc[c2] == '*') or
                    (sc[c1] == sc[c2]))
            # calculate corresponding availability metric for strain code
            mut_pairs = utils.getMutPairs(table)
            for (c1, c2) in mut_pairs
                count += not (
                    (table[c1] == '*' or table[c2] == '*') or
                    (table[c1] == table[c2]))
            # adjust mutation rate by ratio of availability metrics
            f_sc = count_sc / len(utils.tripletMutPairs)
            f_table = count / len(mut_pairs)
            mu_adj = mu*(f_table/t_sc)

        # return the adjusted mutation rate
        return mu_adj

# run as script
if __name__ == '__main__':
    import ipdb; ipdb.set_trace()
    x = codonTable()
    y = strain(table=x)
    print('done')
