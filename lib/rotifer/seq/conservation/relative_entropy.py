#!/usr/bin/env python3

import math
import numpy as np

class relative_entropy:
    def __init__(self,aln, matrix, gap_penalty = 0, pseudocount = 0.0000001,
                z_score = False,
                 **kwargs):
        self.matrix = matrix
        self.aln = aln
        self.pseudocount = pseudocount
        self.z_score = z_score

        self.gap_penalty = gap_penalty
        self.blosum_background_distr = [0.078, 0.051, 0.041,
                                   0.052, 0.024, 0.034,
                                   0.059, 0.083, 0.025,
                                   0.062, 0.092, 0.056,
                                   0.024, 0.044, 0.043,
                                   0.059, 0.055, 0.014,
                                   0.034, 0.072]

        self.amino_acids = ['A', 'R', 'N', 'D',
                            'C', 'Q', 'E', 'G',
                            'H', 'I', 'L', 'K',
                            'M', 'F', 'P', 'S',
                            'T', 'W', 'Y', 'V',
                            '-']

    def run(self,
            **kwargs):
        """Calculate the relative entropy of the column distribution with a
        background distribution specified in bg_distr. This is similar to the
        approach proposed in Wang and Samudrala 06. sim_matrix is ignored."""

        from rotifer.seq.core.core import conservation

        con = conservation(self.aln,self.matrix, self.pseudocount)
        fc = con.weighted_freq_count_pseudocount()
        fc = fc.reindex(self.amino_acids)
        res = []

        #remove gap count
        for col in range(self.matrix.shape[1]):
            distr = self.blosum_background_distr[:]

            if len(distr) == 20:
                new_fc = fc[col][:-1]
                s = sum(new_fc)

                for i in range(len(new_fc)):
                    new_fc[i] = new_fc[i] / s
                fc2 = new_fc

            # if len(fc) != len(distr): return -1

            d = 0.

            for i in range(len(fc2)):
                if distr[i] != 0.0:
                    d += fc2[i] * math.log(fc2[i]/distr[i])

            d /= math.log(len(fc2))

            if self.gap_penalty == 1:
                res.append(d * con._weighted_gap_penalty(col))

            else:
                res.append(d)

        if self.z_score:
            return np.array(con._zscore(res))

        else:
            return np.array(res)

