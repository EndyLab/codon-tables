# import necessary modules
import numpy as np
from scipy.special import comb
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Bio.Data.CodonTable
import networkx as nx
import random
from copy import copy
import pickle
import os

# define path to file directory
path = os.path.dirname(os.path.abspath(__file__))

class utils:
    ''' A module used to package useful definitions and functions when
    manipulating codon tables. Class attributes are listed below.

    Class Attributes
    ----------------
    - dict standardTable: a dict representing the Standard Table
    - list(str) dNTPs: a list of strings representing the DNA NTPs
    - list(str) rNTPs: a list of strings representing the RNA NTPs
    - list(str) tripletCodons: a list of string representing the set of rNTP
        codons
    - set tuples(str) tripletMutPairs: a set of tuples of string representing
        pairs of triplet rNTP codons one mutation away
    - dict PRS: a dict representing the Polar Requirement Scale
    - dict kdHydrophobicity: a dict representing Kyte Doolittle hydrophobicity
    - dict unrestrictedBlock: a dict representing a completely unfettered block
        structure; used in simulations where wobble rule is ignored
    - dict standardBlock: a dict representing the most biologically permissive
        block structure (48 blocks)
    - dict naturalBlock: a dict representing the block structure of the Standard
        Table (25 blocks). More restrictive than other defined block structures
    '''
    ###########################
    # define class properties #
    ###########################

    # define standard table
    standardTable = Bio.Data.CodonTable.standard_rna_table.forward_table
    standardTable['UAA'] = '*'
    standardTable['UAG'] = '*'
    standardTable['UGA'] = '*'
    # unpickle additional class properties
    with open(path+'/res/utilsDefinitions.pickle', 'rb') as handle:
        unPickled = pickle.load(handle)
    [dNTPs, rNTPs, residues, tripletCodons, tripletMutPairs,
     PRS, kdHydrophobicity, Gilis, SCV,
     unrestrictedBlock, standardBlock, naturalBlock] = unPickled

    @staticmethod
    def getAAcounts(table):
        ''' A function that takes a Codon Table and finds the counts of each
        AA. Returns a dictionary mapping AA to their respective counts.

        Parameters
        ----------
        dict table: a python dict representing the codon table

        Returns
        -------
        dict AA_count: a python dict mapping amino acids to degeneracy
        '''
        # declare dictionary of AA counts
        AA_count = {}
        # iterate over key and value pairs in self.table
        for codon, AA in table.items():
            # handle case where AA is previously uncounted
            if AA not in AA_count:
                # add AA to AA_count and initialize count value to 1
                AA_count[AA] = 1
            # else, increment AA count
            else:
                AA_count[AA] += 1
        # return AA_count dictionary
        return AA_count

    @staticmethod
    def getBlockCounts(blocks):
        ''' A function that takes a Codon Table represented in block structure
        form and finds the number of blocks encoding each AA. Returns a
        dictionary mapping AA to their respective counts.

        Parameters
        ----------
        dict blocks: a python dict representing the codon table in block form

        Returns
        -------
        dict blockCounts: a python dict mapping amino acids to degeneracy
        '''
        # initialize dict of counts and populate keys
        blockCounts = {}
        for AA in utils.residues:
            blockCounts[AA] = 0
        # increment counts
        for AA in blocks.values():
            blockCounts[AA] += 1
        # return blockCounts
        return blockCounts

    @staticmethod
    def getCodonConnectivity(table):
        '''getCodonConnectivity(dict table) a function that takes a codon table
        and finds the graph distance between codon pairs. Connectivity is
        defined as follows: two codons c and c' are as connected if a series of
        point mutations can convert c to c' changing which amino acid it
        encodes for, until the last mutation (i.e. only) one AA change is
        allowed per path.

        Outputs a dict of str --> list of (str, int) tuples representing a list
        of the connected codons and their distance. Implemented as breadth
        first search.

        Parameters
        ----------
        dict table: a python dict representing the codon table

        Returns
        -------
        dict distDict: a python dictionary representing the adjacency matrix of
            a codon table with respect to codon neighbors.
        '''
        # declare dictionary of distances
        distDict = {}
        # loop through all possible codons
        for codon in table.keys():
            # declare temporary cache to hold discovered codons; store first codon in there
            cache = set()
            cache.add(codon)
            # declare neighbors list
            neighbors = []
            # use connectRecurse to map connectivity
            distDict[codon] = utils.__connectRecurse(codon, 1, table,
                neighbors, cache);
        # return codonDist
        return distDict

    @staticmethod
    def __connectRecurse(codon, level, table, neighbors, cache):
        ''' A recursive helper function that finds all of a codon's nearest
        neighbors and how far away they are. Returns a list of tuples
        representing the codons and their distances away.

        Codons are said to be connected if going from c --> c' converts the
        decoded AA from A --> A' without an intermediate residue.

        Parameters
        ----------
        - str codon: a string representing the input codon
        - int level: the current number of mutations away from the start codon
        - dict table: a python dict representing the codon table
        - list neighbors: the current list of the base codon's nearest neighbors
        - set cache: memoization set to store previously visited codons
        Returns
        -------
        list neighbors: returns updated neighbors list
        '''
        # define list of codons to recursively look through after searching this level
        recurse_list = []
        # loop through every codon one mutation away
        for i, base in enumerate(codon):
            for nt in utils.rNTPs:
                # handle if nt is the same as base
                if nt == base:
                    continue
                # if not, generate new codon
                c_new = codon[:i] + nt + codon[i+1:]
                # Base case: c_new already found
                if c_new in cache:
                    continue
                # Base case: found terminus
                elif table[c_new] != table[codon]:
                    # add distance to neighbors list
                    neighbors.append((c_new, level))
                    # add c_new to cache of found codons
                    cache.add(c_new)
                # Recursive case
                else:
                    # add c_new to cache of found codons
                    cache.add(c_new)
                    # append c_new to list of codons to recurse through
                    recurse_list.append(c_new)
        # iterate over codons to recursively search for connectivity
        for c in recurse_list:
            # append results to neighbors list
            neighbors = utils.__connectRecurse(c, level + 1,
                                            table, neighbors, cache)
        # return resulting list
        return neighbors

    @staticmethod
    def getResiConnectivity(table):
        ''' getResiConnectivity(dict table): a function that takes a dictionary
        representing a codon table and outputs a dictionary mapping amino acids
        to their respective neighbors, along with number of mutations away.
        '''
        # call getCodonConnectivity
        codonDistDict = utils.getCodonConnectivity(table)
        # declare dict to return
        resiDistDict = {}
        # loop over codons
        for c1, codonNeighbors in codonDistDict.items():
            # extract amino acid for c1 and declare neighbors list
            A1 = table[c1]
            aaNeighbors = []
            # loop over elements of neighbors list
            for (c2, level) in codonNeighbors:
                # convert neighbors to residues and store
                A2 = table[c2]
                aaNeighbors.append((A2, level))
            # store resulting list in resiDistDict
            resiDistDict[A1] = aaNeighbors
        # return dictionary
        return resiDistDict

    @staticmethod
    def getCodonNeighbors(codon):
        '''A static method used to get all codons one mutation away from the given codon.

        Parameters
        ----------
        str codon: the codon whose neighbors will be returned

        Returns
        -------
        list<str> neighbors: a list of codons one mutation away
        '''
        # declare list of neighbors
        neighbors = []
        # generate nearest neighbors by looping over codon positions
        for i, base in enumerate(codon):
            for nt in utils.rNTPs:
                # handle if nt is the same as base
                if nt == base:
                    continue
                # if not, generate new codon
                c_new = codon[:i] + nt + codon[i+1:]
                # store new codon in neighbors
                neighbors.append(c_new)
        # return resulting list
        return neighbors

    @staticmethod
    def tableToBlocks(table, blockStruct):
        '''A static method that takes a codon table and returns the
        representation as blocks of codons (individual tRNAs) as opposed to
        individual codons.

        Parameters
        ----------
        - dict table: a python dict representing the codon table
        - dict blockStruct: a python dict representing the table block structure

        Returns
        -------
        - dict blocks: a python dict representing the codon table in block form
        - bool False: an "exception" if input table does not match blockStruct
        '''
        # run checkBlock to confirm proper block structure, returns False if not
        if utils.checkBlock(table, blockStruct) != True:
            return False
        # declare dictionary to return
        blocks = {}
        # loop over blockStruct and populate blocks
        for blockInd, codonList in blockStruct.items():
            blocks[blockInd] = table[codonList[0]]
        # return populated blocks dict
        return blocks

    @staticmethod
    def blocksToTable(blocks, blockStruct):
        '''A static method that takes a codon table represented in block
        structure form and returns the representation as a traditional codon
        table

        Parameters
        ----------
        - dict blocks: a python dict representing the codon table in block form
        - dict blockStruct: a python dict representing the table block structure

        Returns
        -------
        - dict table: a python dict representing the codon table
        - bool False: an "exception" if input table does not match blockStruct
        '''
        # declare table to return
        table = {}
        # loop over blocks in blockStruct and assign to table using blocks
        for blockInd, codonList in blockStruct.items():
            blockAA = blocks[blockInd]
            for codon in codonList:
                table[codon] = blockAA
        # return filled codon table
        return table

    @staticmethod
    def checkBlock(table, blockStruct):
        '''A static method used to check whether a given codon table conforms
        to the given block structure

        Parameters
        ----------
        - dict table: a python dict representing the codon table
        - dict blockStruct: a python dict representing the table block structure

        Returns
        -------
        bool valid: true->table conforms to block structure; false otherwise
        '''
        # loop over codons in each block; return false if they code for
        # different residues
        for codonList in blockStruct.values():
            #initialize set of residues that a block codes for and populate
            blockResidues = set()
            for codon in codonList:
                blockResidues.add(table[codon])
            # return false if the set is more than one element long
            if len(blockResidues) > 1:
                return False
        # if function reaches this point, return True
        return True

    @staticmethod
    def randomTable(wobble_rule = 'standard'):
        '''A static method used to generate a random codon table, optionally
        defining the block structure. Will guarantee each amino acid be
        represented by at least one block in the table.

        Parameters
        ----------
        str wobble_rule = 'standard': a string telling the simulator which
            wobble rules to follow for accepting new tables

        Acceptable inputs:
        - 'standard' : 48 blocks
        - 'preserveBlock' : maintain same block structure as standard
            table
        - 'unrestricted' : 63 open blocks, at least 1 of every AA and
            stop.

        Returns
        -------
        dict table: a python dict representing the codon table to return
        '''
        # determine block structure based on wobble rule
        blockChoices = {
            'standard' : utils.standardBlock,
            'preserveBlock' : utils.naturalBlock,
            'unrestricted' : utils.unrestrictedBlock
        }
        blockStruct = copy(blockChoices[wobble_rule])
        # get blocks to assign
        blocks = list(blockStruct.keys())
        random.shuffle(blocks)
        # randomly assign one block to each residue
        residues = copy(utils.residues)
        for i, AA in enumerate(residues):
            block = blocks[-(i+1)]
            blockStruct[block] = AA
        # truncate to ignore assigned blocks
        blocks = blocks[:-(i+1)]
        # randomly assign values to the remaining blocks
        for block in blocks:
            AA = random.choice(residues)
            blockStruct[block] = AA
        # convert blockStruct to table and return
        return utils.blocksToTable(blockStruct, blockChoices[wobble_rule])

    @staticmethod
    def numTables(l_aa, b):
        '''A static method used to calculate the number of codon tables
        realizable given a number of amino acids to include, length of the
        codon, and number of blocks. Relies on an inclusion/exclusion criterion
        (i.e. count the total number of codon tables, minus the number that do
        not include one AA, plus the number that do not include two AAs...)

        l_aa = length of amino acid alphabet (20 + 1 stop)
        b = number of blocks to assign (triplet most permissive = 48, quadruplet most permissive = 192)

        n = l_aa^b + Sum_i^(l_aa-1) [(-1)^i * binomial(l_aa, i) * (l_aa - i)^b]

        Parameters
        ----------
        - int l_aa: the number of amino acids + Stop to encode
        - int b: the number of blocks in the codon table

        Returns
        -------
        - int n: the number of possible tables
        - str num: n, represented in scientific notation as a string
        '''
        # handle string processing
        mag = -1
        tempN = n
        while (tempN > 0):
            # increment mag for each order of magnitude
            tempN = tempN // 10
            mag += 1
        # create string representing n in scientific notation
        strN = str(n)[:3]
        num = '{0}.{1}E{2}'.format(strN[0], strN[1:], mag)
        return n, num

    @staticmethod
    def silencicity(table):
        '''A static method used to calculate the silencicity of a codon table.
        Silencicity is a lab defined metric calculating the fraction of all
        mutations that are synonymous out of all possible ones.

        Parameters
        ----------
        dict table: a python dict representing the codon table to analyze

        Returns
        -------
        float silencicity: a float representing the silencicity metric
        '''
        # initialize counter and get number of possible mutation pairs
        synMut = 0
        totalMut = len(utils.tripletMutPairs)
        # loop over mutation pairs and increment for synonymous mutations
        for (c1, c2) in utils.tripletMutPairs:
            if(table[c1] == table[c2]):
                synMut += 1
        # return fraction of synonymous mutations
        return synMut/totalMut

if __name__ == '__main__':
    table = utils.randomTable()
    from codonTable import codonTable
    test = codonTable(table)
    test.plotGraph()
