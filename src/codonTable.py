# import necessary modules
import numpy as np
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as LSC
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import networkx as nx
from src.codonUtils import utils

class codonTable:
    ''' A class used to handle codon table objects.
    '''
    def __init__(self, table=None, ordering=None, norm=True):
        '''Automatically loads object with a codonTable and a
        comparison function between amino acids "ordering". bool norm is used
        to tell dictToGraph whether or not to set node values based on
        absolute value of metric or just the ordering. This nomenclature is
        consistent throughout the object.

        Parameters
        ----------
        - dict table=None: a python dict representing a codon table, optionally
            takes a string 'random' to assign a random table
        - dict ordering=None: represents the mapping of residues to some metric
        - bool norm=True: true->use residue ordering; false->use metric
            absolute value
        Returns
        -------
        codonTable obj: returns a codonTable object
        '''
        # optionally generate a random table, or load a preset table
        if type(table) == str:
            table_options = {
                'STANDARD':utils.standardTable,
                'COLORADO':utils.coloradoTable,
                'RED20':utils.RED20,
                'RED15':utils.RED15,
                'RANDOM':utils.randomTable()
            }
            try:
                table = table_options[table.upper()]
            except:
                raise ValueError('Table string not recognized. Use one of the following options: {0}'.format(set(table_options.keys())))
        # default to standard table if not specified or wrong datatype
        elif table == None or not type(table) == dict:
            table = utils.standardTable

        # Define AA color ordering if not defined
        if ordering == None:
            # default to Kyte and Doolittle ordering
            ordering = utils.kdHydrophobicity

        # determine table ambiguity
        ambiguous = utils.isPromiscuous(table)

        # determine if table is one-to=one
        oneToOne = utils.isOneToOne(table)

        # Assign assign instance attributes
        self.utils = utils
        self.ordering = ordering
        self.codonDict = table
        self.codonTable = self.dictToTable(table)
        self.ambiguous = ambiguous
        self.oneToOne = oneToOne
        # assign attributes for unambiguous tables
        if not ambiguous:
            self.codonGraph = self.dictToGraph(table, norm)
            self.codonSparseMat = nx.adjacency_matrix(self.codonGraph)
            self.codonAdjMat = np.array(self.codonSparseMat.todense())

    def sortOrdering(self):
        '''
        Returns the order of amino acids as defined by the mapping dictionary
        as a numpy array

        Parameters
        ----------
        None

        Returns
        -------
        np.array AA_array: an array of residues sorted by a metric
        '''
        # declare ordering array
        hydr_array = np.zeros(len(self.ordering))
        # declare list of AA
        AA_list = []
        # iterate through dictionary and store
        for i, (AA, hydr) in enumerate(self.ordering.items()):
            AA_list.append(AA)
            hydr_array[i] = hydr

        # convert AA_list to a np.array
        AA_array = np.array(AA_list)
        # sort along hydr_array and return indices
        indices = np.argsort(hydr_array)
        # return sorted AA_array
        return AA_array[indices]

    def getVal(self, AA, norm=True):
        ''' Returns the "value" of an amino acid, using the self.ordering
        dictionary. The norm flag is most relevant to this function, as it
        determines whether to return the absolute value of self.ordering[AA] or
        the amino acid's relative index when ordered.
        '''
            # extract value to store
        if norm == True:
            # use AA order instead of absolute metric value
            order = self.sortOrdering()
            return np.where(order == AA)[0]
        else:
            # use metric value instead of order
            metric = self.ordering
            return metric[AA]

    def getCube(self, norm=True):
        ''' Takes the codon table (in dictionary form) and reporting function
        and returns a cube of 64 locations (represented as a 4x4x4 cube with
        values). For levity and brevity, the returned cube is called Borg. The
        normalization flag (norm=True) tells the method to use the order of AA
        as opposed to the exact value of the ordering metric for visualization

        Parameters
        ----------
        bool norm=True: true->use residue ordering; false->use metric
            absolute value

        Returns
        -------
        np.array Borg: a 4x4x4 np.array representing the codon table in 3D
        '''
            # declare Borg cube
        Borg = np.zeros([4,4,4])
        # define coordinate mappings for codon nt
        codonToInt = {
            'A' : 0,
            'U' : 1,
            'C' : 2,
            'G' : 3,
        }
        # loop over codonTable items
        for codon, AA in self.codonDict.items():
            # extract x, y and z values
            x = codonToInt[codon[0]]
            y = codonToInt[codon[1]]
            z = codonToInt[codon[2]]
            # store value in Borg cube
            Borg[x][y][z] = self.getVal(AA, norm)
        # return the Borg cube
        return Borg

    def getScatterData(self, norm=True):
        ''' Conerts the codon table to a format that is easily parsed by
        scatter(). returns 64 element long arrays representing the x, y, z and
        hydrophobicity values.

        Parameters
        ----------
        bool norm=True: true->use residue ordering; false->use metric
            absolute value

        Returns
        -------
        - np.array xs: represents the x values of position for each codon
        - np.array ys: represents the y values of position for each codon
        - np.array zx: represents the z values of position for each codon
        - np.array vals: represents the metric value for each codon (color)
        '''
        # declare arrays to return
        xs = np.zeros(64)
        ys = np.zeros(64)
        zs = np.zeros(64)
        vals = np.zeros(64)
        stops = []
        # define coordinate mappings for codon nt
        codonToInt = {
            'U' : 0,
            'C' : 1,
            'A' : 2,
            'G' : 3,
        }
        # loop over codonTable items
        for i, (codon, AA) in enumerate(self.codonDict.items()):
            # package stop codons separately into a list of tuples
            if AA == '*':
                pos = (codonToInt[codon[0]], codonToInt[codon[1]],
                    codonToInt[codon[2]])
                stops.append(pos)
                continue
            # extract x, y and z values
            xs[i] = codonToInt[codon[0]]
            ys[i] = codonToInt[codon[1]]
            zs[i] = codonToInt[codon[2]]
            # extract value to store
            vals[i] = self.getVal(AA, norm)
        # return arrays
        return xs, ys, zs, vals, stops

    def dictToTable(self, table=None, suffix=''):
        '''A method used to convert a dict representing a codon table to a
        pandas DataFrame representing the same table in a human readable form

        Paramters
        ---------
        - dict table=None: python dict representing a codon table
        - str suffix: string to append to each codon (used for quadruplet
            tables and higher)

        Returns
        -------
        pd.DataFrame df: pandas DataFrame representing a codon table in a more
            traditional format
        '''
        # handle default options for table
        if table == None:
            table = self.codonDict
        # get list of rNTPs
        NTP = utils.rNTPs
        # recurse if quadruplet code or higher dimmensional
        if (len(list(table)[0]) > 3) and (suffix == ''):
            # declare storage variables
            suffixes = set()
            dfs = []
            # get suffixes
            for codon, AA in table.items():
                suffixes.add(codon[3:])
            # for each suffix, get subtable and convert to dataframe
            suffixes = utils.orderNTPs(list(suffixes))
            subtables = {}
            for suffix in suffixes:
                subdict = {
                    codon:table[codon] for codon in table if codon[3:] == suffix
                }
                subtables[suffix] = subdict
                df = self.dictToTable(table=table, suffix=suffix)
                dfs.append(df)
            # get the concatenation of the resulting dataframes
            DF = pd.concat(dfs, keys=suffixes, copy=False)
        else:
            # declare list of dataframes
            dfs = []
            # for each first position, define new dict to become a dataframe
            for c1 in NTP:
                dfDict = {}
                # for each second position, generate a new list of codon/AA pairs
                for c2 in NTP:
                    row = []
                    # for each third position, populate row list
                    for c3 in NTP:
                        codon = c1+c2+c3+suffix
                        '''NOTE: this is where things get funky'''
                        element = '{0} : {1}'.format(codon, table[codon])
                        row.append(element)
                    # add row to dict under label of second position
                    dfDict[c2] = row
                # create new data frame and add to list of data frames
                df = pd.DataFrame(dfDict, index=NTP)
                dfs.append(df[NTP])
            # get the concatenation of the resulting dataframes
            DF = pd.concat(dfs, keys=NTP, copy=False)
        return DF

    def dictToGraph(self, table=None, norm=True):
        ''' Takes a dictionary representing a codon table as an input and
        returns a networkx graph representing that table. If no table is given,
        the function will default to using self.codonDict.

        Parameters
        ----------
        - dict table=None: dict representing a codon table
        - bool norm=True: true->use residue ordering; false->use metric
            absolute value

        Returns
        -------
        nx.Graph codonGraph: a networkx graph representing the codon table
        '''
        # define Prob of single point mutations
        p_mut = 1/12
        # handle default table values
        if table == None:
            table = self.codonDict
        # declare graph
        codonGraph = nx.Graph()
        # call utils.getAAcounts to get degeneracy
        aaCounts = self.utils.getAAcounts(table)
        # call utils.getResiConnectivity to get dict representing graph
        graphDict = self.utils.getResiConnectivity(table)
        # initialize all nodes
        for resi in graphDict.keys():
            # add node for current residue
            count = aaCounts[resi]
            kd = self.getVal(resi, norm)
            codonGraph.add_node(resi, count=count, kd=kd)

        # parse edges from graphDict into networkx graph
        for A1, neighbors in graphDict.items():
            # declare dictionary to hold edge weights
            edgeDict = {}
            # populate edge dictionary
            for (A2, level) in neighbors:
                if A2 not in edgeDict:
                    edgeDict[A2] = p_mut**level
                else:
                    edgeDict[A2] += p_mut**level
            # convert edge dictionary into edges
            for A2, weight in edgeDict.items():
                codonGraph.add_edge(A1, A2, weight=weight)

        #return network
        return codonGraph

    def plot3d(self, title="", color="viridis", ctitle="", norm=True):
        ''' Represents self.table in 3D space. Returns the figure handle of the
        visualization. Optionally puts a title in the figure.

        Parameters
        ----------
        - str title="": an optional input to define the title of the plot
        - str color="viridis": an optional input used to set the color of the
            color. Accepted inputs:
                - red
                - green
                - blue
                - purple
                - grey
                - viridis
        - str ctitle="": an optional input to define title of color bar
        - bool norm=True: true->use residue ordering; false->use metric
            absolute value

        Returns
        -------
        plt.figure fig: matplotlib figure handle for the resulting plot
        '''
        # raise error if table is higher dimensional than triplet
        dimension = len(list(self.codonDict)[0])
        if (dimension > 3):
            raise ValueError(
                'Cannot plot 3d representation of {0}D genetic code'.format(dimension)
            )
        # call getScatterData to extract arrays from codon table
        xs, ys, zs, vals, stops = self.getScatterData(norm)
        # initialize figure and axes objects
        fig = plt.figure()
        ax = Axes3D(fig)
        # pick colormap
        if (color == "viridis"):
            softened = False
        else:
            softened = True
        cmap = self.colormap(color, softened)
        # plot data
        scat = ax.scatter(xs, ys, zs, s=1200, cmap=cmap, c=vals)
        # plot colorbar
        cbar = plt.colorbar(scat, shrink=0.8, location="left")
        cbar.set_alpha(0.5)
        cbar.outline.set_visible(False)
        cbar.set_ticks([])
        cbar.set_label(ctitle)
        # plot stops as grey balls
        xstop = []
        ystop = []
        zstop = []
        for (x, y, z) in stops:
            xstop.append(x)
            ystop.append(y)
            zstop.append(z)
        ax.scatter(xstop,ystop,zstop,s=1200,c='grey')
        # set axes labels and ticks
        ticks = np.arange(0,4)
        tick_labels = ['U', 'C', 'A', 'G']

        ax.set_xlabel('Position 1')
        ax.set_ylabel('Position 2')
        ax.set_zlabel('Position 3')

        ax.set_xticks(ticks)
        ax.set_xticklabels(tick_labels)

        ax.set_yticks(ticks)
        ax.set_yticklabels(tick_labels)

        ax.set_zticks(ticks)
        ax.set_zticklabels(tick_labels)
        # format title
        if title != "":
            plt.title(title)
        #show figure
        ax.view_init(elev=8, azim=-11)
        plt.show()
        #return figure
        return fig

    def plotGraph(self, title="", ctitle="",
                    color='viridis', norm=True,
                    nodeSize='count', nodeColor='kd',
                    weightingFxn=None,
                    filename=None):
        ''' Represents self.codonDict as a network capturing the adjacency of
        the amino acids. Two amino acids are defined as adjacent if a codon
        representing AA_1 can be mutated to represent AA_2 without representing
        an intermediate residue. The edge weight is defined by the probability
        of mutating to another AA:

            P(a1, a2) = \Sigma_i p^{n_i}

        nodeSize tells the plotter what values to use in representing node ize.
        nodeColor tells the plotter what should be represented using the node
        color value (defaults to kdHydropathy). Returns figure handle for
        plotted figure.

        Parameters
        ----------
        - str title="": an optional input to define the title of the plot
        - str ctitle="": an optional input to define title of color bar
        - bool norm=True: true->use residue ordering; false->use metric
            absolute value
        - str color="viridis": an optional input used to set the color of the
            network. Accepted inputs:
                - red
                - green
                - blue
                - purple
                - grey
                - viridis
        - str nodeSize='degeneracy': an optional input to specify the
            node size based on a metric. Defaults to degeneracy
        - str nodeColor='kd': an optional input to specify the
            node size based on a metric. Defaults to degeneracy
        - fxn(array<floats>) weightingFxn: a function used to remap edge
            weights. Optional
        - str filename: optionally save figure to file with filename

        Returns
        -------
        plt.figure fig: matplotlib figure handle for the resulting plot

        Known Bugs
        ----------
        Edge weight may be changing between iterations! confirm this and fix
        '''
        # pick colormap
        if (color == "viridis"):
            softened = False
        else:
            softened = True
        cmap = self.colormap(color, softened)
        # unpack graph
        G = self.codonGraph
        # unpack and normalize edge weights and node values for visualization
        weights = np.array(
            [edge['weight'] for (a1, a2, edge) in G.edges(data=True)]
        )
        # Optionally weigh edges nonlinearly for easier visualization
        if weightingFxn != None:
            weights = weightingFxn(weights)
        weights /= np.mean(weights)
        node_size = [data[nodeSize]*200 for (node, data) in G.nodes(data=True)]
        node_color = np.array(
            [data[nodeColor] for (node, data) in G.nodes(data=True)]
        ).ravel()
        # set up layout
        positions = nx.spring_layout(G, iterations=100)
        # draw graph
        fig = plt.figure()
        plt.axis('off')
        nodes = nx.draw_networkx_nodes(G, positions, cmap=cmap,
                node_color=node_color, node_size=node_size)
        edges = nx.draw_networkx_edges(G, positions, width=weights)
        labels = nx.draw_networkx_labels(G, positions)
        # plot colorbar
        cbar = plt.colorbar(nodes)
        cbar.outline.set_visible(False)
        cbar.set_ticks([])
        cbar.set_label(ctitle)
        # recolor stop and null codons to white, grey, respectively
        stops = []; nulls = [];
        for i, AA in enumerate(G.nodes()):
            if AA == '*':
                stops.append(int(i))
            elif AA == '0':
                nulls.append(int(i))
        stops = np.array(stops, dtype=int)
        nulls = np.array(nulls, dtype=int)
        nx.draw_networkx_nodes(G, positions,
            nodelist=['*'], node_size=node_size[stops[0]],
            node_color='white')
        nx.draw_networkx_nodes(G, positions,
            nodelist=['0'], node_size=node_size[nulls[0]],
            node_color='grey')
        #format graph
        if title != "": plt.title(title)
        # optionally save figure
        if type(filename) == str: plt.savefig(filename)
        plt.show()
        return fig

    @staticmethod
    def colormap(color, softened=False):
        '''a static method used to generate a color bar for plots.
        Optionally cuts off low end of the color scheme for visibility

        Parameters
        ----------
        - str color: a string representing the color of color bar to use.
            Accepted inputs:
                - red
                - green
                - blue
                - purple
                - grey
                - viridis
        - bool softened=False: a bool telling the function whether or not to
            clip lighter colors for visibility

        Returns
        -------
        cmap: the resulting color bar
        '''
        # pick colormap from options
        options = {
            'red' : cm.Reds,
            'blue' : cm.Blues,
            'green' : cm.Greens,
            'purple' : cm.Purples,
            'grey' : cm.Greys,
            'viridis' : cm.viridis
        }
        cmap = options[color]
        # optionally clip colormap
        if softened:
            start = 0.3
            stop = 1
            colors = cmap(np.linspace(start, stop, cmap.N))
            cmap = LSC.from_list('Upper Half', colors)
        return cmap

    def reverseTranslate(self, prot_seq, stop_codon='UGA'):
        '''a method used to translate a given protein sequence, assuming the
        codon table is one-to-one. Raises an error otherwise

        Parameters
        ----------
        - str prot_seq: string representing amino acid sequence to translate
        - str stop_codon: stop codon to use for GOI (default to UGA)

        Returns
        -------
        str gene: string representing translated amino acid sequence
        '''
        # handle error if table is not one-to-one
        if not self.oneToOne:
            raise TypeError('Cannot translate sequence. Codon table is not One-To-One')
        # otherwise, create reverse translation dictionary
        rev_dict = {aa:codon for codon, aa in self.codonDict.items()}
        rev_dict['*'] = stop_codon
        # translate gene and return
        gene = ''
        for aa in prot_seq:
            gene += rev_dict[aa]
        return gene

### Test script
if __name__ == '__main__':
    test = codonTable('RED20')
    test_peptide ='MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*'
    gene = test.reverseTranslate(test_peptide)
    print('Reverse Translated Gene {0}'.format(gene))
