# import necessary modules
import numpy as np
import pandas as pd
from src.codonUtils import utils
from src.codonTable import codonTable
from queue import Queue as queue
from random import choice

class ffgen():
    ''' A class used to generate fast fail codon tables
    '''

    def __init__(self):
        return

    @staticmethod
    def triplet():
        '''
        A static method used to generate triplet decoding, fast fail genetic codes
        using a rational, 'top down' approach. Snakes along codon table to fill
        16 maximally distant codons first (all greater than one mutation from
        each other), then randomly chooses 4 codons to place the remaining
        amino acids.

        Parameters
        ----------
        None

        Returns
        ----------
        dict table: a triplet fast fail table
        '''
        # return table built by reducedTriplet, with no knockouts
        return ffgen.reducedTriplet(knockout=0)

    @staticmethod
    def reducedTriplet(knockout=0):
        '''
        A static method used to generate triplet decoding, fast fail genetic
        codes with an arbitrary number of amino acids randomly removed from the
        set of encoded amino acids. Uses a rational, 'top down' approach.
        Snakes along codon table to fill 16 maximally distant codons first (all
        greater than one mutation from each other), then randomly chooses 4 - n
        codons to place the remaining amino acids. If the number of knockouts
        is greater than 4, it skips the second step and only includes 16 - n in
        the first step.

        Parameters
        ----------
        int knockout: number of amino acids to omit

        Returns
        ----------
        dict table: a triplet fast fail table
        '''

        ############################
        # fill out first 16 codons #
        ############################

        # declare caching variables
        usedCodons = set()
        unusedAA = set(utils.residues[:-1])
        # declare budding table
        table = {}
        # randomly permute rNTPs for positions 1 and 2, store 3 unpermuted
        pos1 = np.random.permutation(utils.rNTPs)
        pos2 = np.random.permutation(utils.rNTPs)
        pos3 = utils.rNTPs
        # store nucleotides for positions 1 and 2 in queues
        queue1 = queue()
        queue2 = queue()
        for nt1, nt2 in zip(pos1, pos2):
            queue1.put(nt1)
            queue2.put(nt2)

        # handle knockout cases
        if knockout >= 20: 
            # if all amino acids are to be knocked out, return a table of all
            # stops
            for codon in utils.tripletCodons:
                table[codon] = '*'
        else:
            # otherwise, populate first 16 elements of table (unless knockout
            # is reached)
            count = 0
            for i in range(queue1.qsize()):
                # get first nucleotide of next codon
                nt1 = queue1.get()
                for j in range(queue2.qsize()):
                    # break if knockout is reached
                    if count > (20 - knockout): break
                    # get second nucleotide of next codon
                    nt2 = queue2.get()
                    # get third nucleotide
                    nt3 = pos3[j]
                    # build codon, assign it a residue
                    codon = nt1 + nt2 + nt3
                    AA = choice(tuple(unusedAA))
                    table[codon] = AA
                    # update caching variables
                    usedCodons.add(codon)
                    unusedAA.remove(AA)
                    # re-enqueue second nucleotide
                    queue2.put(nt2)
                    # increment counter
                    count += 1
                # re-enqueue first nucleotide
                queue1.put(nt1)
                # dequeue and re-enqueue second position to shift array
                nt2 = queue2.get()
                queue2.put(nt2)

            # place the last 4 elements, if necessary 
            numLeft = 4 - knockout
            if numLeft > 0:
                availableCodons = set(utils.tripletCodons) - usedCodons
                for i in range(numLeft):
                    # pick a codon from the available set and assign it an amino acid
                    codon = choice(tuple(availableCodons))
                    AA = choice(tuple(unusedAA))
                    table[codon] = AA
                    # update cache variables
                    availableCodons = ffgen.updateAvailable3(codon, availableCodons)
                    usedCodons.add(codon)
                    unusedAA.remove(AA)

            # assign unused codons to STOP
            remainingCodons = set(utils.tripletCodons) - usedCodons
            for codon in remainingCodons:
                table[codon] = '*'

        # return built table
        return table

    @staticmethod
    def updateAvailable3(newCodon, availableSet):
        ''' A static method used to update the set of codons that can be used
        for triplet decoding fast fail code, given that a new codon is
        occupied.''' 

        # iterate over remaining codons
        copySet = list(availableSet)
        for codon in copySet:
            # remove codons that have two nucleotide overlaps
            count = 0
            for i in range(len(newCodon)):
                count += (codon[i] == newCodon[i])
            if count >= 2:
                availableSet.remove(codon)
        # return updated set
        return availableSet

if __name__ == '__main__':
    # availableCodons = set(utils.tripletCodons)
    # count = 0
    # while len(availableCodons) > 0:
    #     codon = random.choice(tuple(availableCodons))
    #     availableCodons = ffgen.updateAvailable3(codon, availableCodons)
    #     count += 1
    print(ffgen.triplet())
