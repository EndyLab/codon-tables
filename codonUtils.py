# import necessary modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Bio.Data.CodonTable
import networkx as nx
import pickle

class utils():
    ###########################
    # define class properties #
    ###########################

    # define standard table
    standardTable = Bio.Data.CodonTable.standard_rna_table.forward_table
    standardTable['UAA'] = '*'
    standardTable['UAG'] = '*'
    standardTable['UGA'] = '*'
    # unpickle additional class properties
    with open('res/utilsDefinitions.pickle', 'rb') as handle:
        unPickled = pickle.load(handle)
    [dNTPs, rNTPs, PRS, kdHydrophobicity, unrestrictedBlock,
        standardBlock, naturalBlock] = unPickled

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
            neighbors = utils.connectRecurse(c, level + 1, table, neighbors, cache)

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
