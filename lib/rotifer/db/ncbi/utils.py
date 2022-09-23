import types

def best_ipgs(ipgs):
    best = ipgs.sort_values(['id','order'], ascending=[True,True])
    best = best.drop_duplicates(['id'], keep='first')
    return best

def ipgs_to_dicts(ipgs):
    if len(ipgs) == 0:
        return dict(), dict()

    # Split
    assemblies = ipgs[ipgs.assembly.notna()].id.unique().tolist()
    nucleotides = ipgs[~ipgs.id.isin(assemblies)]
    assemblies = ipgs[ipgs.id.isin(assemblies)]

    # By assembly
    if len(assemblies) > 0:
        assemblies = assemblies.assembly.unique().tolist()
        assemblies = ipgs[ipgs.assembly.isin(assemblies)]
        assemblies = assemblies.filter(['assembly','pid','representative'])
        assemblies = assemblies.drop_duplicates(ignore_index=True)
        assemblies = assemblies.groupby('assembly').apply(lambda x: x.set_index('pid').representative.to_dict())
        assemblies = assemblies.to_dict()
    else:
        assemblies = dict()

    # By nucleotide
    if len(nucleotides) > 0:
        nucleotides = nucleotides.nucleotide.unique().tolist()
        nucleotides = ipgs[ipgs.nucleotide.isin(nucleotides)]
        nucleotides = nucleotides.filter(['nucleotide','pid','representative'])
        nucleotides = nucleotides.drop_duplicates(ignore_index=True)
        nucleotides = nucleotides.groupby('nucleotide').apply(lambda x: x.set_index('pid').representative.to_dict())
        nucleotides = nucleotides.to_dict()
    else:
        nucleotides = dict()

    return assemblies, nucleotides


