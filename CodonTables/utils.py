# import necessary modules
import numpy as np
from scipy.special import comb
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Bio.Data.CodonTable
import networkx as nx
from collections import deque
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
    - dict standard_table: a dict representing the Standard Table
    - dict colorado_table: a dict representing the optimized code from Pines et
        al 2017
    - dict RED20: a dict representing RED20
    - dict RED15: a dict representing RED15
    - list(str) dNTPs: a list of strings representing the DNA NTPs
    - list(str) rNTPs: a list of strings representing the RNA NTPs
    - list(str) triplet_codons: a list of string representing the set of rNTP
        codons
    - set tuples(str) triplet_mut_pairs: a set of tuples of string representing
        pairs of triplet rNTP codons one mutation away
    - list(str) quadruplet_codons: a list of string representing the set of rNTP
        codons
    - set tuples(str) quadruplet_mut_pairs: a set of tuples of string representing
        pairs of quadruplet rNTP codons one mutation away
    - dict PRS: a dict representing the Polar Requirement Scale
    - dict kd_hydropathy: a dict representing Kyte Doolittle hydropath
    - dict unrestricted_block: a dict representing a completely unfettered block
        structure; used in simulations where wobble rule is ignored
    - dict standard_block: a dict representing the most biologically permissive
        block structure (48 blocks)
    - dict natural_block: a dict representing the block structure of the Standard
        Table (25 blocks). More restrictive than other defined block structures
    '''
    ###########################
    # define class properties #
    ###########################

    # define standard table
    standard_table = Bio.Data.CodonTable.standard_rna_table.forward_table
    standard_table['UAA'] = '*'
    standard_table['UAG'] = '*'
    standard_table['UGA'] = '*'
    # unpickle additional class properties
    with open(path+'/res/utils_definitions.pickle', 'rb') as handle:
        un_pickled = pickle.load(handle)
    [dNTPs, rNTPs, residues, triplet_codons, triplet_mut_pairs,
     quadruplet_codons, quadruplet_mut_pairs,
     PRS, kd_hydropathy, Gilis, SCV,
     unrestricted_block, standard_block, natural_block,
     basepair_WC, wobble_WC, colorado_table] = un_pickled
    # unpickle RED15 and RED20 tables
    with open(path+'/res/RED20.pickle', 'rb') as handle:
         RED20 = pickle.load(handle)
    with open(path+'/res/RED15.pickle', 'rb') as handle:
         RED15 = pickle.load(handle)

    @staticmethod
    def get_aa_counts(table):
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
    def get_block_counts(blocks):
        ''' A function that takes a Codon Table represented in block structure
        form and finds the number of blocks encoding each AA. Returns a
        dictionary mapping AA to their respective counts.

        Parameters
        ----------
        dict blocks: a python dict representing the codon table in block form

        Returns
        -------
        dict block_counts: a python dict mapping amino acids to degeneracy
        '''
        # initialize dict of counts and populate keys
        block_counts = {}
        for AA in utils.residues:
            block_counts[AA] = 0
        # increment counts
        for AA in blocks.values():
            block_counts[AA] += 1
        # return block_counts
        return block_counts

    @staticmethod
    def is_ambiguous(table):
        '''A staticmethod that takes a codon table as a dictionary and returns True if it is ambiguous and False if not.

        Parameters
        ----------
        dict table: a python dict representing the codon table

        Returns
        -------
        bool ambiguous: boolean representing the ambiguity of the table
        '''
        # use utils.promiscuity method to determine ambiguity
        try:
            __ = utils.promiscuity(table, allow_ambiguous=False) #fails if ambiguous
            ambiguous = False
        except:
            ambiguous =True
        return ambiguous

    @staticmethod
    def is_promiscuous(table):
        '''A staticmethod that takes a codon table as a dictionary and returns True if it represents a promiscuous table and False if not.

        Parameters
        ----------
        dict table: a python dict representing the codon table

        Returns
        -------
        bool ambiguous: boolean representing the promiscuity of the table
        '''
        # this is a one liner, but a tad obfuscated. Checks to see if each codon encodes for only one AA (thus is type str).
        # returns true if any of the elements are not strings
        return sum(type(AA) != str for AA in table.values()) > 0

    @staticmethod
    def is_one_to_one(table):
        '''A staticmethod that takes a codon table as a dictionary and returns
            True if it represents a One-To-One genetic code and False otherwise.

            A one-to-one code is defined as a code in which every amino acid is
            represented with exactly one codon. This defines an unambiguous
            mapping of protein sequence to corresponding DNA sequence.

        Parameters
        ----------
        dict table: a python dict representing the codon table

        Returns
        -------
        bool one2one: boolean; True if One-To-One, and False otherwise
        '''
        # declare storage dict to count amino acid number
        aa_set = set(aa for aa in table.values())
        aa_counts = {aa:0 for aa in aa_set}
        # count number of amino acids
        for aa in table.values():
            aa_counts[aa] += 1
        # iterate through dictionary and check counts
        one2one = True
        for aa, count in aa_counts.items():
            #skip stop and null signals:
            if aa in {'*', '0'}:
                continue
            elif count > 1:
                one2one = False
                break
        return one2one

    @staticmethod
    def get_codon_connectivity(table):
        '''get_codon_connectivity(dict table) a function that takes a codon table
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
        dict dist_dict: a python dictionary representing the adjacency matrix of
            a codon table with respect to codon neighbors.
        '''
        # declare dictionary of distances
        dist_dict = {}
        # loop through all possible codons
        for codon in table.keys():
            # declare temporary cache to hold discovered codons; store first codon in there
            cache = set(codon)
            # declare queue of codons to check
            codon_deque = deque()
            # declare neighbors list
            neighbors = []
            # use connect_recurse to map connectivity
            dist_dict[codon] = utils.__connect_recurse(codon, 1, table,
                neighbors, codon_deque, cache);
        # return codon_dist
        return dist_dict

    @staticmethod
    def __connect_recurse(codon, level, table, neighbors, codon_deque, cache):
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
        - deque codon_deque: the queue of codons to search recursively
        - set cache: memoization set to store previously visited codons
        Returns
        -------
        list neighbors: returns updated neighbors list
        '''
        # import ipdb; ipdb.set_trace()
        # loop through every codon one mutation away
        for i, base in enumerate(codon):
            for nt in utils.rNTPs:
                # handle if nt is the same as base
                if nt == base: continue
                # if not, generate new codon
                c_new = codon[:i] + nt + codon[i+1:]
                # Base case: c_new already found
                if c_new in cache: continue
                # Base case: found terminus
                elif table[c_new] != table[codon]:
                    # add distance to neighbors list
                    neighbors.append((c_new, level))
                    # add c_new to cache of found codons
                    cache.add(str(c_new))
                # Recursive case
                else:
                    # add c_new to cache of found codons
                    cache.add(c_new)
                    # append c_new to queue of codons to recurse through
                    codon_deque.appendleft((c_new, level))
        # iterate over codons to recursively search for connectivity
        while not len(codon_deque) == 0:
            # get next codon to search
            c, newlevel = codon_deque.pop()
            # append results to neighbors list
            neighbors = utils.__connect_recurse(c, newlevel + 1, table,
                                               neighbors, codon_deque, cache)
        # return resulting list
        return neighbors

    @staticmethod
    def get_resi_connectivity(table):
        ''' get_resi_connectivity(dict table): a function that takes a dictionary
        representing a codon table and outputs a dictionary mapping amino acids
        to their respective neighbors, along with number of mutations away.
        '''
        # call get_codon_connectivity
        codon_dist_dict = utils.get_codon_connectivity(table)
        # declare dict to return
        resi_dist_dict = {}
        # loop over codons
        for c1, codon_neighbors in codon_dist_dict.items():
            # extract amino acid for c1 and declare neighbors list
            A1 = table[c1]
            aa_neighbors = []
            # loop over elements of neighbors list
            for (c2, level) in codon_neighbors:
                # convert neighbors to residues and store
                A2 = table[c2]
                aa_neighbors.append((A2, level))
            # store resulting list in resi_dist_dict
            if A1 not in resi_dist_dict:
                resi_dist_dict[A1] = aa_neighbors
            else:
                resi_dist_dict[A1] += aa_neighbors
        # return dictionary
        return resi_dist_dict

    @staticmethod
    def get_codon_neighbors(codon):
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
    def table_to_blocks(table, block_struct):
        '''A static method that takes a codon table and returns the
        representation as blocks of codons (individual tRNAs) as opposed to
        individual codons.

        Parameters
        ----------
        - dict table: a python dict representing the codon table
        - dict block_struct: a python dict representing the table block structure

        Returns
        -------
        - dict blocks: a python dict representing the codon table in block form
        - bool False: an "exception" if input table does not match block_struct
        '''
        # run check_block to confirm proper block structure, returns False if not
        if utils.check_block(table, block_struct) != True:
            return False
        # declare dictionary to return
        blocks = {}
        # loop over block_struct and populate blocks
        for block_ind, codon_list in block_struct.items():
            blocks[block_ind] = table[codon_list[0]]
        # return populated blocks dict
        return blocks

    @staticmethod
    def blocks_to_table(blocks, block_struct):
        '''A static method that takes a codon table represented in block
        structure form and returns the representation as a traditional codon
        table

        Parameters
        ----------
        - dict blocks: a python dict representing the codon table in block form
        - dict block_struct: a python dict representing the table block structure

        Returns
        -------
        - dict table: a python dict representing the codon table
        - bool False: an "exception" if input table does not match block_struct
        '''
        # declare table to return
        table = {}
        # loop over blocks in block_struct and assign to table using blocks
        for block_ind, codon_list in block_struct.items():
            block_aa = blocks[block_ind]
            for codon in codon_list:
                table[codon] = block_aa
        # return filled codon table
        return table

    @staticmethod
    def check_block(table, block_struct):
        '''A static method used to check whether a given codon table conforms
        to the given block structure

        Parameters
        ----------
        - dict table: a python dict representing the codon table
        - dict block_struct: a python dict representing the table block structure

        Returns
        -------
        bool valid: true->table conforms to block structure; false otherwise
        '''
        # loop over codons in each block; return false if they code for
        # different residues
        for codon_list in block_struct.values():
            #initialize set of residues that a block codes for and populate
            block_residues = set()
            for codon in codon_list:
                block_residues.add(table[codon])
            # return false if the set is more than one element long
            if len(block_residues) > 1:
                return False
        # if function reaches this point, return True
        return True

    @staticmethod
    def random_table(wobble_rule = 'standard'):
        '''A static method used to generate a random codon table, optionally
        defining the block structure. Will guarantee each amino acid be
        represented by at least one block in the table.

        Parameters
        ----------
        str wobble_rule = 'standard': a string telling the simulator which
            wobble rules to follow for accepting new tables

        Acceptable inputs:
        - 'standard' : 48 blocks
        - 'preserve_block' : maintain same block structure as standard
            table
        - 'unrestricted' : 63 open blocks, at least 1 of every AA and
            stop.

        Returns
        -------
        dict table: a python dict representing the codon table to return
        '''
        # determine block structure based on wobble rule
        block_choices = {
            'standard' : utils.standard_block,
            'preserve_block' : utils.natural_block,
            'unrestricted' : utils.unrestricted_block
        }
        block_struct = copy(block_choices[wobble_rule])
        # get blocks to assign
        blocks = list(block_struct.keys())
        random.shuffle(blocks)
        # randomly assign one block to each residue
        residues = copy(utils.residues)
        for i, AA in enumerate(residues):
            block = blocks[-(i+1)]
            block_struct[block] = AA
        # truncate to ignore assigned blocks
        blocks = blocks[:-(i+1)]
        # randomly assign values to the remaining blocks
        for block in blocks:
            AA = random.choice(residues)
            block_struct[block] = AA
        # convert block_struct to table and return
        return utils.blocks_to_table(block_struct, block_choices[wobble_rule])

    @staticmethod
    def num_tables(l_aa, b):
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
        temp_n = n
        while (temp_n > 0):
            # increment mag for each order of magnitude
            temp_n = temp_n // 10
            mag += 1
        # create string representing n in scientific notation
        str_n = str(n)[:3]
        num = '{0}.{1}E{2}'.format(str_n[0], str_n[1:], mag)
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
        # initialize counter and get mutation pairs
        syn_mut = 0
        mut_pairs = utils.get_mut_pairs(table)
        total_mut = len(mut_pairs)
        # loop over mutation pairs and increment for synonymous mutations
        for (c1, c2) in mut_pairs:
            if(table[c1] == table[c2]):
                syn_mut += 1
        # return fraction of synonymous mutations
        return syn_mut/total_mut

    @staticmethod
    def mutability(table):
        '''A static method used to calculate the average chemical variability
        of single point mutations in a given genetic code. For each
        nonsynonymous single point mutation, it calculates the chemical
        distance between the previously encoded amino acid and its replacement
        after mutation. The mean of these values is then returned.

        Parameters
        ----------
        dict table: a python dict representing the codon table to analyze

        Returns
        -------
        float mutability: a float representing the silencicity metric
        '''
        # initialize counter and running metric, and get mutation pairs
        nonsyn_mut = 0
        metric = 0
        mut_pairs = utils.get_mut_pairs(table)
        total_mut = len(mut_pairs)
        # get Kyte-Doolittle hydropathy metric
        kd = utils.kd_hydropathy
        # loop over mutation pairs
        for (c1, c2) in mut_pairs:
            # increment counter and metric if nonsynonymous
            if not (table[c1] == table[c2]):
                # increment counter
                nonsyn_mut += 1
                # increment metric
                aa1 = table[c1]
                aa2 = table[c2]
                metric += np.abs(kd[aa1] - kd[aa2])
        # if there are no nonsynonymous mutations, return 0
        if nonsyn_mut == 0:
            mutability = 0
        # else, return the average dKD per mutation
        else:
            mutability = metric / nonsyn_mut
        return mutability

    @staticmethod
    def promiscuity(table, allow_ambiguous=False):
        '''A static method used to generate the genetic code resulting from
        considering tRNA promiscuity. Uses Crick Wobble Hypothesis. Raises an
        exception if the table generated is ambiguous (more than one signal
        acceptable for a given codon)

        Parameters
        ----------
        - dict table: the codon table to promsicuitize
        - bool allow_ambiguous: flag telling code whether to accept ambiguity

        Returns
        -------
        dict promsicuous: the resulting table when considering tRNA promiscuity
        '''
        # handle type errors for input table
        if type(table) != dict:
            try:
                table = table.codon_dict # attempt to convert to dict
            except:
                raise ValueError("Input table is not type dict or CodonTable")
        # declare table to return
        promiscuous = {}
        for codon in utils.triplet_codons:
            promiscuous[codon] = '0'
        # loop over codons to reassign
        for codon, AA in table.items():
            # skip assignments to STOP
            if AA == '0': continue
            # get codons that would be decoded in reality
            wobble = utils.wobble_WC[utils.basepair_WC[codon[-1]]]
            codons = [codon[:2]+nt3 for nt3 in wobble]
            # determine if there is ambiguity
            acceptable = [AA, '0']
            for c in codons:
                if promiscuous[c] not in acceptable:
                    # raise error if allow_ambiguous = False
                    if not allow_ambiguous:
                        raise ValueError('input code generates ambiguous code upon promiscuization')
                    else:
                        # else, package all nonstop codons as tuple
                        AAs = tuple(
                            [aa for aa in promiscuous[c] if aa != '0'] +
                            [AA]
                        )
                        promiscuous[c] = AAs
                # otherwise, package as simple str --> str mapping
                else:
                    promiscuous[c] = AA
        return promiscuous

    @staticmethod
    def mut_pair_num(table):
        '''
        A static method that calculates the number of pairs of codons one
        mutation away from each other. Treats mutations with directionality. In
        general, the number of mutational pairs is equal to the number of
        codons in a table multiplied by the number of unique codons within one
        mutation. Let a = alphabet length (generally 4), l = codon length
        (generally 3)

                n = (a^l) * l(a-1)

        Parameters
        ----------
        dict table: the codon table to analyze

        Returns
        -------
        int mut_num: the number of distinct mutational pairs.
        '''
        # get list of all codons in table
        codon_list = list(table)
        # get alphabet size
        alphabet = set()
        for codon in codon_list:
            for nt in codon:
                alphabet.add(nt)
        a = len(alphabet)
        # get codon length
        l = len(codon_list[0])
        # calculate mut_num and return
        return (a**l) * l * (a-1)

    @staticmethod
    def get_mut_pairs(table):
        '''
        A static method used to generate the set of all pairs of codons one
        mutation away given a codon table.

        Parameters
        ----------
        dict table: the codon table to analyze

        Returns
        -------
        set<(str, str)> mut_pairs: a set of distinct mutational pairs.
        '''
        # declare set of mutational pairs
        mut_pairs = set()
        # get list of codons and iterate over them
        codon_list = list(table)
        for codon in codon_list:
            # iterate over each base in the codon
            for i, base in enumerate(codon):
                for nt in utils.rNTPs:
                    # handle if nt is the same as base
                    if nt == base:
                        continue
                    # if not, generate new codon
                    c_new = codon[:i] + nt + codon[i+1:]
                    # add to set
                    mut_pairs.add((codon, c_new))
        return mut_pairs

    @staticmethod
    def order_NTPs(sortable, nucleic_acid='RNA'):
        '''A static method used to sort iterables by standard order of NTPs. For RNA, U-C-A-G. For DNA, T-C-A-G. Returns sorted object.

        Parameters
        ----------
        - iterable sortable: the object to sort
        - str nucleic_acid: the type of nucleic acid considered

        Returns
        -------
        iterable sorted_obj: the sorted object
        '''
        # define ordering dictionary
        orderdict = {
            'RNA' : ['U', 'C', 'A', 'G'],
            'DNA' : ['T', 'C', 'A', 'G']
        }
        # raise error if nucleic_acid flag invalid
        if nucleic_acid.upper() not in orderdict:
            raise ValueError('nucleic_acid flag set to invalid option (use DNA or RNA)')
        # attempt sorting
        try:
            order = orderdict[nucleic_acid.upper()]
            sorted_obj = sorted(
                sortable, key=lambda word: [order.index(nt) for nt in word]
            )
        except ValueError:
            print('Variable to sort broke the code :/')
            # raise error
            sorted_obj = False
        return sorted_obj


if __name__ == '__main__':
    table = {
        'UUU' : 'F',
        'UCA' : 'S',
        'UCG' : 'L',
        'AUG' : 'M',
    }
    newtable = utils.promiscuity(table, allow_ambiguous=True)
    from CodonTables.table import CodonTable
    new_table = CodonTable(newtable)
    new_table.codon_dict
