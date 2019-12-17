#import dependencies
import numpy as np
import matplotlib.pyplot as plt
import random
from CodonTables.table import CodonTable
from CodonTables.utils import utils
from copy import copy

#define Monte Carlo Simulation class
class TableOptimizer:
    '''A class designed to optimize a codon table given an arbitrary objective
    function to minimize/maximize
    '''
    def __init__(self, table=utils.standardTable,
                costfxn=None, wobble_rule = 'standard',
                debug=False):
        '''the init function for the TableOptimizer class. Optionally allows the
        user to specify the starting codon table, associated objective
        function, and which wobble_rules should be followed. Default values are
        supplied if the user does not specify them (i.e. Standard Code, 47
        block wobble rule)

        Parameters
        ----------
        - dict table = utils.standardTable: a python dict representing the
            codon table. Also optionally accepts a CodonTable object.
        - func costfxn = None: a function that takes a python dict as an
            input and outputs a cost associated with that table. If no
            funciton is supplied, defaults to maxMutMinVar.
        - str wobble_rule = 'standard': a string telling the simulator which
            wobble rules to follow for accepting new tables

            Acceptable inputs:
            - 'standard' : 48 blocks, 2 stop blocks
            - 'preserveBlock' : maintain same block structure as standard
                table
            - 'unrestricted' : 63 open blocks, at least 0 of every AA and
                stop.

        - bool debug: optional boolean telling the initializer whether or not
            to run self.__debug()

        Returns
        -------
        TableOptimizer obj: returns an instance of the TableOptimizer object
        '''
        # handle CodonTable objects being passed
        if type(table) == CodonTable:
            table = CodonTable.CodonTable
        # handle costfxn input
        if costfxn == None:
            costfxn = self.maxMutMinVar
        # calculate some attributes
        resiCounts = utils.getAAcounts(table)
        connectivity = utils.getCodonConnectivity(table)
        # determine block structure based on wobble rule
        blockChoices = {
            'standard' : utils.standardBlock,
            'preserveBlock' : utils.naturalBlock,
            'unrestricted' : utils.unrestrictedBlock
        }
        # assign attributes
        self.table = table
        self.blockStruct = copy(blockChoices[wobble_rule])
        self.costfxn = costfxn
        self.resiCounts = resiCounts
        self.connectivity = connectivity
        self.utils = utils
        self.wobble_rule = wobble_rule

        # optionally run debugger
        if debug:
            self.__debugGDA()

    ######################################
    ##          Public  Methods         ##
    ######################################

    def GDA(self, dW=None, W=0, W_stop=float("inf"), maxIter=1000,
            preserveBlock = False, preserveStop = False, subFunc='minPRS'):
        '''The Great Deluge Algorithm for optimizing an objective function over
        codon table space

        Parameters
        ----------
        - float dW=None: represents the rate of change (rain flux) of min
            allowable energy (water level); if None, initialized to 1% of
            cost(self.table)
        - float W=0: represents the initial minimal acceptable energy level
            (water level); defaults to 0
        - float W_stop=inf: represents the maximum allowable water level; stops
            algorithm when W > W_stop; defaults to an infinte value
        - int maxIter=1000: represents the maximum number of iterations to
            perform; defaults to 1000
        - bool preserveBlock=False: a bool that tells the function whether or
            not to preserve block structure when shuffling the table
        - bool preserveStop=False: a bool that tells the function whether or not
            to shuffle blocks encoding for STOP
        - str subFunc: an optional parameter specifying the substitution metric
            to be used in maximization. Defaults to __sub_minPRS().

            Acceptable inputs:
            - 'minPRS' --> __sub_minPRS()
            - 'Gilis' --> __sub_Gilis()
            - 'SCV' --> __sub_SCV()

        Returns
        -------
        - dict table: a python dict representing an optimized codon table
        - np.array Ws: a numpy array representing water level vs iteration
        - np.array Es: a numpy array representing cost(table) vs iteration
        '''
        # define dictionary of acceptable substitution metrics
        subs = {
            'minPRS' : self.__sub_minPRS,
            'Gilis' : self.__sub_Gilis,
            'SCV' : self.__sub_SCV
        }
        # declare data structures to return
        table = self.table
        Ws = np.zeros(maxIter)
        Es = np.zeros(maxIter)
        # initialize iteration counter and Energy value
        count = 0
        E = self.costfxn(table, subs[subFunc])
        # if standard table is below initial water level, adjust W
        if (E < W):
            W = E
        # populate first elements of arrays
        Ws[0] = W
        Es[0] = E
        # if not specified, calculate dW to be 1% of initial Es
        if dW is None:
            dW = abs(E/100)
        # loop through algorithm until an end condition is met
        while (count < maxIter - 1):
            # increment counter and make a tentative move
            count += 1
            newTable = self.tableShuffle(table, preserveBlock, preserveStop)
            E_new = self.costfxn(newTable, subs[subFunc])
            # determine whether or not to accept change
            if (E_new > W):
                # if so, update table, energy level and water level
                table = newTable
                E = E_new
                W += dW
            # store values in arrays
            Ws[count] = W
            Es[count] = E
            # determine whether or not to end algorithm based on water level
            if W > W_stop:
                # pad arrays with current value to end
                if count + 1 < maxIter:
                    Ws[count+1:] = W
                    Es[count+1:] = E
                # break loop
                break
        # return table and arrays
        return table, Ws, Es

    def tableShuffle(self, table, preserveBlock = False, preserveStop = False):
        '''Takes a codon table as an input and shuffles it to produce a new,
            similar table in order to traverse table space.

        Parameters
        ----------
        - dict table: a python dict representing a codon table starting point
        - bool preserveBlock=False: a bool that tells the function whether or
            not to preserve block structure when shuffling the table
        - bool preserveStop=False: a bool that tells the function whether or not
            to shuffle blocks encoding for STOP
        Returns
        -------
        dict newTable: a python dict representing the next codon table
        '''
        # convert table to block form and get residue degeneracy
        blockForm = utils.tableToBlocks(table, self.blockStruct)
        blockCounts = utils.getBlockCounts(blockForm)
        # Case: switch two blocks
        if (random.random() < 0.5) | preserveBlock:
            switchable = [block for block, AA in blockForm.items()
                            if not (preserveStop and AA == '*')]
            toBlock = random.choice(switchable)
            fromBlock = random.choice([block for block in switchable
                                        if not block == toBlock])
            # get AA for toBlock, switch with fromBlock, then set fromBlock
            tempAA = blockForm[toBlock]
            blockForm[toBlock] = blockForm[fromBlock]
            blockForm[fromBlock] = tempAA
        # Case: change identity of a block
        else:
            # get list of blocks that can be mutated
            mutable = [block for block, AA in blockForm.items()
                    if blockCounts[AA] > 1
                    and not (preserveStop and AA == '*')]
            # get list of AA that can be inserted into block
            listAA = [AA for AA in utils.residues
                    if not (preserveStop and AA == '*')]
            # pick which block to switch and perform switch
            changeBlock = random.choice(mutable)
            blockForm[changeBlock] = random.choice(listAA)
        # convert blockForm back to table form and return
        return utils.blocksToTable(blockForm, self.blockStruct)

    ######################################
    ##          Cost  Functions         ##
    ##                                  ##
    ## NOTE: all cost functions should  ##
    ##  take the same inputs!!          ##
    ######################################

    @staticmethod
    def maxMutMinVar(table, subFunc):
        '''the default cost function for class. Implements a version of the
        cost function used by Novozhilov et al. 2007, but optimizes for
        maximizing mutability while minimizing variance per mutation

        Parameters
        ----------
        - dict table: a python dictionary representing the codon table to
            evaluate
        - func subFunc(str AA_1, str AA_2): a function representing the penalty
            for substituting two amino acids

        Returns
        -------
        float metric: the absolute number score of the given table
        '''
        # initialize cost variable, declare polar requirement scale
        metric = 0
        PRS = utils.PRS
        # generate dictionary of AA counts and codon connectivity
        counts = utils.getAAcounts(table)
        connectivity = utils.getCodonConnectivity(table)
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
                d_c12 = subFunc(AA_1, AA_2)
                ## this version of d is good for increasing diff b/w AA
                # d_c12 = (1 + (PRS[AA_1] - PRS[AA_2])**2)
                P_c12 = (1/12)**dist
                # calculate contribution to metric from this pair
                metric += f_c * P_c12 * d_c12
        # return resulting metric
        return metric

    ######################################
    ##          Private Methods         ##
    ######################################

    # substitution function
    @staticmethod
    def __sub_minPRS(AA_1, AA_2):
        '''a private method used to calculate the cost of substituting two amino
        acids using the polar requirement scale. Higher values are given to subs
        that minimize the change in PRS. Output is bounded by [0,1].

        Parameters
        ----------
        - str AA_1: a string representing the "from" amino acid
        - str AA_2: a string representing the "to" amino acid

        Returns
        -------
        float d_c12: a float representing the cost of this substitution
        '''
        # unpack PRS dict and return calculated cost
        PRS = utils.PRS
        return 1/(1 + (PRS[AA_1] - PRS[AA_2])**2)

    @staticmethod
    def __sub_Gilis(AA_1, AA_2):
        '''a private method used to calculate the cost of substituting two amino
        acids using the Gilis substitution matrix. Essentially just packages
        the input for accessing the utils.Gilis dict. This matrix returns
        larger values for substitutions that preserve physicochemistry and
        smaller ones otherwise. Returns int values in range [-6, 7]

        Parameters
        ----------
        - str AA_1: a string representing the "from" amino acid
        - str AA_2: a string representing the "to" amino acid

        Returns
        -------
        int d_c12: an int representing the cost of this substitution
        '''
        # unpack Gilis dict and return  cost
        Gilis = utils.Gilis
        return Gilis[(AA_1, AA_2)]

    @staticmethod
    def __sub_SCV(AA_1, AA_2):
        '''a private method used to calculate the cost of substituting two amino
        acids using the SCV substitution matrix. Essentially just packages the
        input for accessing the utils.SCV dict. This matrix returns larger
        values for substitutions that preserve physicochemistry and smaller
        ones otherwise. Returns int values in range [-15, 12]

        Parameters
        ----------
        - str AA_1: a string representing the "from" amino acid
        - str AA_2: a string representing the "to" amino acid

        Returns
        -------
        int d_c12: an int representing the cost of this substitution
        '''
        # unpack SCV dict and return cost
        SCV = utils.SCV
        return SCV[(AA_1, AA_2)]

    def __debugGDA(self):
        '''a private method used to debug GDA
        '''
        table, Ws, Es = self.GDA(maxIter=3000, preserveBlock=True, preserveStop=True)
        # plot results
        fig, axArray = plt.subplots(2, sharex=True)
        axArray[0].plot(Es)
        axArray[0].set_title('Simulation Metrics vs Iteration')
        axArray[0].set_ylabel('Energy')
        axArray[1].plot(Ws)
        axArray[1].set_ylabel('Water')
        plt.show()
        # represent resulting codon tables
        Table = CodonTable(table)
        fig2 = Table.plot3d('Optimized Codon Table: Node Color=Hydropathy')
        fig3 = Table.plotGraph('Optimized Codon Graph: Node Color=Residue Degeneracy', node_val='count')

# Debugging
if __name__ == '__main__':
    sim = TableOptimizer(debug=True)
