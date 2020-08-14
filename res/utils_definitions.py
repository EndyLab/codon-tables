import pickle
import pandas as pd

# define NTPs list
dNTPs = ['T', 'C', 'A', 'G']
rNTPs = ['U', 'C', 'A', 'G']
# define list of natural amino acids, including stop
residues = [
    'G', 'A', 'V', 'L', 'I', 'P', 'M', 'C', 'S', 'F', 'Y',
    'W', 'T', 'N', 'Q', 'D', 'E', 'R', 'H', 'K', '*'
];
# define polar requirement scale
PRS = {
    'F' : 5.0,
    'L' : 4.9,
    'I' : 4.9,
    'M' : 5.3,
    'V' : 5.6,
    'S' : 7.5,
    'P' : 6.6,
    'T' : 6.6,
    'A' : 7.0,
    'Y' : 5.7,
    'H' : 8.4,
    'Q' : 8.6,
    'N' : 10.0,
    'K' : 10.1,
    'D' : 13.0,
    'E' : 12.5,
    'C' : 11.5,
    'W' : 5.3,
    'R' : 9.1,
    'G' : 7.9,
}

# define default hydrophobicity metric
kdHydrophobicity = {
    'I' : 4.5,
    'V' : 4.2,
    'L' : 3.8,
    'F' : 2.8,
    'C' : 2.5,
    'M' : 1.9,
    'A' : 1.8,
    '*' : 0, #stop codon's given neutral value for visualization
    '0' : 0, #null codon's given neutral value for visualization
    'G' : -0.4,
    'T' : -0.7,
    'W' : -0.9,
    'S' : -0.8,
    'Y' : -1.3,
    'P' : -1.6,
    'H' : -3.2,
    'E' : -3.5,
    'Q' : -3.5,
    'D' : -3.5,
    'N' : -3.5,
    'K' : -3.9,
    'R' : -4.5,
}

# define block structures and codon set
unrestrictedBlock = {}
tripletCodons = []
count = 0;
for c1 in rNTPs:
    for c2 in rNTPs:
        for c3 in rNTPs:
            tripletCodons.append(c1+c2+c3)
            unrestrictedBlock[count] = [c1+c2+c3]
            count+=1

standardBlock = {}
for i in range(48):
    standardBlock[i] = []

count = 0
for c1 in rNTPs:
    for c2 in rNTPs:
        for c3 in rNTPs:
            standardBlock[count].append(c1+c2+c3)
            if (c3 != 'U'):
                count+=1

# define quadruplet codon set
quadrupletCodons = []
for nt1 in rNTPs:
    for nt2 in rNTPs:
        for nt3 in rNTPs:
            for nt4 in rNTPs:
                quadrupletCodons.append(nt1+nt2+nt3+nt4)

naturalBlock = {
    0 : ['UUU', 'UUC'],
    1 : ['UUA', 'UUG'],
    2 : ['CUU', 'CUC', 'CUA', 'CUG'],
    3 : ['AUU', 'AUC', 'AUA'],
    4 : ['AUG'],
    5 : ['GUU', 'GUC', 'GUA', 'GUG'],
    6 : ['UCU', 'UCC', 'UCA', 'UCG'],
    7 : ['CCU', 'CCC', 'CCA', 'CCG'],
    8 : ['ACU', 'ACC', 'ACA', 'ACG'],
    9 : ['GCU', 'GCC', 'GCA', 'GCG'],
    10: ['UAU', 'UAC'],
    11: ['UAA', 'UAG'],
    12: ['CAU', 'CAC'],
    13: ['CAA', 'CAG'],
    14: ['AAU', 'AAC'],
    15: ['AAA', 'AAG'],
    16: ['GAU', 'GAC'],
    17: ['GAA', 'GAG'],
    18: ['UGU', 'UGC'],
    19: ['UGA'],
    20: ['UGG'],
    21: ['CGU', 'CGC', 'CGA', 'CGG'],
    22: ['AGU', 'AGC'],
    23: ['AGA', 'AGG'],
    24: ['GGU', 'GGC', 'GGA', 'GGG']
}

# define Watson Crick Wobbling Rules
basepairWC = {
    'U' : 'A',
    'C' : 'G',
    'A' : 'U',
    'G' : 'C'
}
wobbleWC = {
    'U' : ['A', 'G'],
    'C' : ['G'],
    'A' : ['U', 'C'],
    'G' : ['U', 'C'],
    'I' : ['A', 'C', 'U']
}

#define all pairs of codons 1 mutation away
tripletMutPairs = set()
for codon in tripletCodons:
    for i, base in enumerate(codon):
        for nt in rNTPs:
            # handle if nt is the same as base
            if nt == base:
                continue
            # if not, generate new codon
            c_new = codon[:i] + nt + codon[i+1:]
            # add to set
            tripletMutPairs.add((codon, c_new))

# define all pairs of quadruplet codons 1 mutation away
quadrupletMutPairs = set()
for codon in quadrupletCodons:
    for i, base in enumerate(codon):
        for nt in rNTPs:
            # handle if nt is the same as base
            if nt == base:
                continue
            # if not, generate new codon
            c_new = codon[:i] + nt + codon[i+1:]
            # add to set
            quadrupletMutPairs.add((codon, c_new))
# define refactored [sic] code from Pines et al 2017 (aka Colorado code)
coloradoTable = {
    'GAA' : 'V',
    'UCG' : 'V',
    'CGU' : 'V',
    'UGA' : 'L',
    'AAU' : 'L',
    'CUC' : 'L',
    'CCA' : 'I',
    'GGG' : 'I',
    'UUU' : 'I',
    'UAC' : 'I',
    'CAG' : 'A',
    'AUA' : 'A',
    'GCU' : 'A',
    'AGC' : 'A',
    'GAU' : 'E',
    'ACA' : 'E',
    'UUC' : 'E',
    'CGG' : 'E',
    'UGU' : 'D',
    'AAC' : 'D',
    'GUG' : 'D',
    'UAA' : '*',
    'UCU' : 'P',
    'AUG' : 'P',
    'GUC' : 'P',
    'CAA' : 'P',
    'GAC' : 'T',
    'UCA' : 'T',
    'CCC' : 'S',
    'AGG' : 'S',
    'AUU' : 'Q',
    'GGA' : 'Q',
    'UGC' : 'N',
    'CAU' : 'N',
    'GCG' : 'M',
    'CUA' : 'M',
    'AAA' : 'C',
    'UUG' : 'C',
    'GGU' : 'C',
    'CUU' : 'G',
    'AGU' : 'G',
    'ACC' : 'G',
    'UAG' : 'G',
    'UGG' : 'R',
    'GCA' : 'R',
    'CAC' : 'R',
    'GGC' : 'H',
    'CCG' : 'H',
    'UUA' : 'H',
    'ACU' : 'H',
    'CGA' : 'K',
    'UCC' : 'K',
    'GUU' : 'K',
    'AAG' : 'K',
    'CCU' : 'Y',
    'GAG' : 'Y',
    'AUC' : 'Y',
    'CGC' : 'W',
    'ACG' : 'W',
    'GUA' : 'W',
    'UAU' : 'W',
    'GCC' : 'F',
    'CUG' : 'F',
    'AGA' : 'F',
}
# massage Gilis.csv files into the proper format
df = pd.read_csv('res/Gilis.csv')
indices = df['Unnamed: 0'].tolist()
df.pop('Unnamed: 0')
df.index = indices
# populate dictionary representing the substitution matrix
Gilis = {}

for AA1 in indices:
    for AA2 in indices:
        Gilis[(AA1, AA2)] = df[AA1][AA2]

# repeat for scv.csv
df = pd.read_csv('res/scv.csv')
indices = df['Unnamed: 0'].tolist()
df.pop('Unnamed: 0')
df.index = indices

SCV = {}

for AA1 in indices:
    for AA2 in indices:
        SCV[(AA1, AA2)] = df[AA1][AA2]

if __name__ == '__main__':
    # time to pickle!
    toDump = [dNTPs, rNTPs, residues, tripletCodons, tripletMutPairs,
                quadrupletCodons, quadrupletMutPairs,
                PRS, kdHydrophobicity, Gilis, SCV,
                unrestrictedBlock, standardBlock, naturalBlock,
                basepairWC, wobbleWC,
                coloradoTable]
    with open('res/utilsDefinitions.pickle', 'wb') as handle:
        pickle.dump(toDump, handle)

    # test the pickle
    with open('res/utilsDefinitions.pickle', 'rb') as handle:
        unDumped = pickle.load(handle)
    # taste the pickle
    if (toDump == unDumped) :
        print("Mmmm, that's a tasty pickle")
    else:
        print("Hmmm, something's up with this pickle")
