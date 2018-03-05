import pickle
from src.codonUtils import utils

# initialize table and set all signals to stops
fftable = {}
for codon in utils.tripletCodons:
    fftable[codon] = '*'

# populate individual codons
fftable['UUC'] = 'C'
fftable['CUU'] = 'A'
fftable['AUU'] = 'Q'
fftable['AUA'] = 'G'
fftable['GUG'] = 'N'
fftable['UCA'] = 'W'
fftable['CCU'] = 'P'
fftable['CCC'] = 'D'
fftable['ACG'] = 'R'
fftable['GCU'] = 'H'
fftable['GCA'] = 'M'
fftable['UAU'] = 'Y'
fftable['UAG'] = 'E'
fftable['CAG'] = 'S'
fftable['AAC'] = 'T'
fftable['GAA'] = 'F'
fftable['UGG'] = 'V'
fftable['CGA'] = 'K'
fftable['AGU'] = 'I'
fftable['GGC'] = 'L'

# pickle data
filename = 'ff_table_manuscript.pickle'
with open('res/'+filename, 'wb') as handle:
    pickle.dump(fftable, handle)
