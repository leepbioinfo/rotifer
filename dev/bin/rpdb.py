#!/usr/bin/env python3
import sys
from Bio.PDB import *

# Get distance between all pairs

# Get get secondary structure

import rotifer.core.cli as corecli

__version__ = 0.001
__authors__ = 'Gilberto Kaihami'

def parse_cli():
    parser = corecli.parser(description = 'Get distance from pdb file')

    parser.add( long_arg = '--accession',
                short_arg = '-a',
                dest = 'accession',
                helper = 'Input PDB file',
                action = "store"
                )

    parser.add( long_arg = '--output',
                short_arg = '-o',
                dest = 'output',
                helper = 'Output (distance / buried)',
                action = "store"
                )
    # Add another options here

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = parse_cli()
    pdb1 = args.accession

    parser=PDBParser(PERMISSIVE=1)

    structure = parser.get_structure('1UBQ',pdb1)

#    rd = ResidueDepth(model, pdb_file)
    model = structure[0]

    chain = model['A']

    # residue distance all vs all
    if args.output == 'distance':
        for res in chain:
            for res2 in chain:
                try:
                    print(res.id[1],res.id[2], res2.id[1], res2, res['CA'] - res2['CA'], sep = '\t')
                except:
                    pass

    # Residue depth

    else:
        rd = ResidueDepth(model)

        for res in chain:
            try:
                print(rd[('A', res.id)]) # residue_depth, ca_depth
            except:
                pass




