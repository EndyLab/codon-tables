from codonTable import codonTable
from codonUtils import utils

table = utils.standardTable
connectivity = utils.getResiConnectivity(table)
print(len(connectivity))
