import pickle
if __name__ == '__main__':
    # define NTPs list
    dNTPs = ['T', 'C', 'A', 'G']
    rNTPs = ['U', 'C', 'A', 'G']
    # define list of natural amino acids, including stop
    residues = ['G', 'A', 'V', 'L', 'I', 'P', 'M', 'C', 'S', 'F', 'Y', 'W', 'T', 'N', 'Q', 'D', 'E', 'R', 'H', 'K', '*'];
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

    # define block structures
    unrestrictedBlock = {}
    count = 0;
    for c1 in rNTPs:
        for c2 in rNTPs:
            for c3 in rNTPs:
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

    # time to pickle!
    toDump = [dNTPs, rNTPs, residues, PRS, kdHydrophobicity, unrestrictedBlock, standardBlock, naturalBlock]
    with open('utilsDefinitions.pickle', 'wb') as handle:
        pickle.dump(toDump, handle)

    # test the pickle
    with open('utilsDefinitions.pickle', 'rb') as handle:
        unDumped = pickle.load(handle)
    # taste the pickle
    if (toDump == unDumped) :
        print("Mmmm, that's a tasty pickle")
    else:
        print("Hmmm, something's up with this pickle")
