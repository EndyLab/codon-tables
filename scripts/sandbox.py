# define physical RED15 and RED20 tables
red15 = {
    codon:'*' for codon in utils.tripletCodons
}
red15['UUC'] = 'F'
red15['CUA'] = 'L'
red15['AUG'] = 'M'
red15['GUC'] = 'V'
red15['UCA'] = 'S'
red15['CCG'] = 'P'
red15['ACC'] = 'T'
red15['GCA'] = 'A'
red15['UAC'] = 'Y'
red15['CAG'] = 'Q'
red15['AAA'] = 'K'
red15['GAC'] = 'D'
red15['UGG'] = 'W'
red15['CGU'] = 'R'
red15['GGA'] = 'G'

red20 = copy.deepcopy(red15)
red20['AUC'] = 'I'
red20['CAC'] = 'H'
red20['AAC'] = 'N'
red20['GAA'] = 'E'
red20['UGC'] = 'C'
