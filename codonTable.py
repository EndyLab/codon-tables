# import necessary modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
from codonUtils import utils

class codonTable:
    ''' A class used to handle codon table objects.
    '''
    def __init__(self, codonTable=None, ordering=None, norm=True):
        '''Automatically loads object with a codonTable and a
        comparison function between amino acids "ordering". bool norm is used
        to tell tableToGraph whether or not to set node values based on
        absolute value of metric or just the ordering. This nomenclature is
        consistent throughout the object.

        Parameters
        ----------
        - dict codonTable=None: represents codon table
        - dict ordering=None: represents the mapping of residues to some metric
        - bool norm=True: true->use residue ordering; false->use metric
            absolute value
        Returns
        -------
        codonTable obj: returns a codonTable object
        '''
        # Define codonTable if not user defined
        if codonTable == None:
            # load standard codon table using Biopython
            codonTable = utils.standardTable
        # Define AA color ordering if not defined
        if ordering == None:
            # default to Kyte and Doolittle ordering
            ordering = utils.kdHydrophobicity

        # Assign assign instance attributes
        self.utils = utils
        self.codonTable = codonTable
        self.ordering = ordering
        self.codonGraph = self.tableToGraph(codonTable, norm)
        self.sparseMat = nx.adjacency_matrix(self.codonGraph)
        self.adjMat = np.array(self.sparseMat.todense())

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
        for codon, AA in self.codonTable.items():
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
        # define coordinate mappings for codon nt
        codonToInt = {
            'A' : 0,
            'U' : 1,
            'C' : 2,
            'G' : 3,
        }
        # loop over codonTable items
        for i, (codon, AA) in enumerate(self.codonTable.items()):
            # skip stop codons for now
            if AA == '*':
                continue
            # extract x, y and z values
            xs[i] = codonToInt[codon[0]]
            ys[i] = codonToInt[codon[1]]
            zs[i] = codonToInt[codon[2]]
            # extract value to store
            vals[i] = self.getVal(AA, norm)
        # return arrays
        return xs, ys, zs, vals

    def tableToGraph(self, table=None, norm=True):
        ''' Takes a dictionary representing a codon table as an input and
        returns a networkx graph representing that table. If no table is given,
        the function will default to using self.codonTable.

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
            table = self.codonTable
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

    def plot3d(self, title="", norm=True):
        ''' Represents self.table in 3D space. Returns the figure handle of the
        visualization. Optionally puts a title in the figure.

        Parameters
        ----------
        - str title="": an optional input to define the title of the plot
        - bool norm=True: true->use residue ordering; false->use metric
            absolute value

        Returns
        -------
        plt.figure fig: matplotlib figure handle for the resulting plot
        '''
        # call getScatterData to extract arrays from codon table
        xs, ys, zs, vals = self.getScatterData(norm)
        # initialize figure and axes objects
        fig = plt.figure()
        ax = Axes3D(fig)
        # plot data
        scat = ax.scatter(xs,ys,zs,s=1000,c=vals)
        cbar = plt.colorbar(scat)
        # set axes labels and ticks
        ticks = np.arange(0,4)
        tick_labels = ['A', 'U', 'C', 'G']

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
        plt.show()
        #return figure
        return fig

    def plotGraph(self, title="", norm=True,
                    nodeSize='count', nodeColor='kd'):
        ''' Represents self.codonTable as a network capturing the adjacency of
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
        - bool norm=True: true->use residue ordering; false->use metric
            absolute value
        - str nodeSize='degeneracy': an optional input to specify the
            node size based on a metric. Defaults to degeneracy
        - str nodeColor='kd': an optional input to specify the
            node size based on a metric. Defaults to degeneracy

        Returns
        -------
        plt.figure fig: matplotlib figure handle for the resulting plot

        Known Bugs
        ----------
        Edge weight may be changing between iterations! confirm this and fix
        '''
        # unpack graph
        G = self.codonGraph
        # unpack and normalize edge weights and node values for visualization
        weights = np.array(
            [edge['weight'] for (a1, a2, edge) in G.edges(data=True)]
        )
        weights /= np.mean(weights)
        node_size = [data[nodeSize]*200 for (node, data) in G.nodes(data=True)]
        node_color = [data[nodeColor] for (node, data) in G.nodes(data=True)]
        # set up layout
        positions = nx.spring_layout(G, iterations=100)
        # draw graph
        fig = plt.figure()
        plt.axis('off')
        nodes = nx.draw_networkx_nodes(G, positions,
                node_color=node_color, node_size=node_size)
        edges = nx.draw_networkx_edges(G, positions, width=weights)
        labels = nx.draw_networkx_labels(G, positions)
        # draw color bar
        cbar = plt.colorbar(nodes)
        # format graph
        if title != "":
            plt.title(title)
        plt.show()
        return fig

### Test script
if __name__ == '__main__':
    test = codonTable()
    fig = test.plot3d('Standard Codon Table: Node Color=Hydropathy')
    fig2 = test.plotGraph('Standard Codon Table: Node Color=Hydropathy Degeneracy', nodeSize='count', nodeColor='kd')
