#!/usr/bin/env python3

import math
import numpy as np

class property_entropy:
    def __init__(self, aln, matrix, gap_penalty = 0, pseudocount = 0.0000001,
                 z_score = False, property_partition = 'Mirny',
                 **kwargs):
        self.matrix = matrix
        self.aln = aln
        self.pseudocount = pseudocount
        self.z_score = z_score

        self.gap_penalty = gap_penalty
        blosum_background_distr = [0.078, 0.051, 0.041,
                                   0.052, 0.024, 0.034,
                                   0.059, 0.083, 0.025,
                                   0.062, 0.092, 0.056,
                                   0.024, 0.044, 0.043,
                                   0.059, 0.055, 0.014,
                                   0.034, 0.072]

        self.property_partition = property_partition

        self.amino_acids = ['A', 'R', 'N', 'D',
                            'C', 'Q', 'E', 'G',
                            'H', 'I', 'L', 'K',
                            'M', 'F', 'P', 'S',
                            'T', 'W', 'Y', 'V',
                            '-']

        #dictionary to map from amino acid to its row/column in a similarity matrix
        self.aa_to_index = {}
        for i, aa in enumerate(self.amino_acids):
            self.aa_to_index[aa] = i

    def run(self,
            **kwargs):
        """
        Calculate the entropy of a column col relative to a partition of the
        amino acids. Similar to Mirny '99. sim_matrix and bg_distr are ignored, but
        could be used to define the sets.

        # Source code modified from: Capra JA and Singh M.
        # Predicting functionally important residues from sequence
        # conservation. Bioinformatics. 23(15): 1875-1882, 2007.
        """
        from rotifer.seq.core.core import conservation

        # Mirny and Shakn. '99'
        # Williamson '95
        property_partition_dc = {'Mirny': [['A','V','L','I','M','C'],
                                    ['F','W','Y','H'],
                                    ['S','T','N','Q'],
                                    ['K','R'],
                                    ['D', 'E'],
                                    ['G', 'P'],
                                    ['-']],
                                 'Williamson': [['V','L', 'I','M'],
                                         ['F','W','Y'],
                                         ['S','T'],
                                         ['N','Q'],
                                         ['H','K','R'],
                                         ['D','E'],
                                         ['A','G'],
                                         ['P'],
                                         ['C'],
                                         ['-']] }


        property_partition = property_partition_dc[self.property_partition]

        con = conservation(self.aln,self.matrix, self.pseudocount)
        fc = con.weighted_freq_count_pseudocount()
        res = []

        # sum the aa frequencies to get the property frequencies
        for col in range(self.matrix.shape[1]):
            prop_fc = [0.] * len(property_partition)

            for p in range(len(property_partition)):
                for aa in property_partition[p]:
                    prop_fc[p] += fc.loc[aa,col]

            h = 0.
            for i in range(len(prop_fc)):
                if prop_fc[i] != 0:
                    h += prop_fc[i] * math.log(prop_fc[i])

            h /= math.log(min(len(property_partition), self.matrix.shape[0]))

            inf_score = 1 - (-1 * h)

            if self.gap_penalty == 1:
                res.append(inf_score * con._weighted_gap_penalty(col))
            else:
                res.append(inf_score)

        if self.z_score:
            return np.array(con._zscore(res))

        else:
            return np.array(res)
