#!/usr/bin/env python3
import math
import numpy as np
import pandas as pd
##
class bit_score:
    def __init__(self,aln, matrix, gap_penalty = 0, pseudocount = 0.0000001,
                 z_score = False,
                 **kwargs):
        self.matrix = matrix
        self.aln = aln
        self.pseudocount = pseudocount
        self.z_score = z_score
        # self.normalized = normalized

        self.gap_penalty = gap_penalty
        blosum_background_distr = [0.078, 0.051, 0.041,
                                   0.052, 0.024, 0.034,
                                   0.059, 0.083, 0.025,
                                   0.062, 0.092, 0.056,
                                   0.024, 0.044, 0.043,
                                   0.059, 0.055, 0.014,
                                   0.034, 0.072]


        # Ri = log2(20) -(Hi - en)
        # i: position in the aln, Hi (Shannon entropy), en:The approximation for the small-sample correction
        # Hi = - (sum (f(b,i) * log2(f(b,i))) b is the aa/base i is the position
        # en = (1/ln2) * (s-1)/(2n) where s = 4 for DNA/RNA 20 for aa n the number of sequences
        # the height for a given aa is:
        # f(b,i) * Ri where the height will be for base b.

    def run(self,
             **kwargs):
        """
        col is each column

        Calculates the Shannon entropy of the column col. sim_matrix  and
        bg_distr are ignored. If gap_penalty == 1, then gaps are penalized. The
        entropy will be between zero and one because of its base. See p.13 of
        Valdar 02 for details. The information score 1 - h is returned for the sake
        of consistency with other scores.

        # Source code modified from: Capra JA and Singh M.
        # Predicting functionally important residues from sequence
        # conservation. Bioinformatics. 23(15): 1875-1882, 2007.
        """
        from rotifer.seq.core.core import conservation
        con = conservation(self.aln,self.matrix, self.pseudocount)
        fc = con.weighted_freq_count_pseudocount()
        res = []

        en_1 = (1/math.log(2))
        en_2 = ( 19 / 2 * self.aln.shape[1])
        en = en_1/en_2
        fc = fc.drop('-')

        df = pd.DataFrame(index = fc.index)
        for col in range(self.matrix.shape[1]):
            Hi = fc[col].apply(lambda x: x * math.log2(x)).sum()
            Hi = (-1) * Hi
            corrected_Hi = Hi+en
            Ri = (math.log2(20)) - (corrected_Hi)
            bits = fc[col].apply(lambda x: x * Ri)
            df = df.join(bits)

        return df





