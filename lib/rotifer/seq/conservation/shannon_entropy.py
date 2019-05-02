#!/usr/bin/env python3
import math
import numpy as np
##
class shannon_entropy:
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

        for col in range(self.matrix.shape[1]):

            h = 0.
            for i in range(len(fc[col])):
                if fc[col][i] != 0:
                    h += fc[col][i] * math.log(fc[col][i])

            #h /= math.log(len(fc))
            h /= math.log(min(len(fc[col]), self.matrix.shape[0]))

            inf_score = 1 - (-1 * h)

            if self.gap_penalty == 1:
                res.append(inf_score * con._weighted_gap_penalty(col))
            else:
                res.append(inf_score)

        if self.z_score:
            return np.array(con._zscore(res))

        else:
            return np.array(res)


