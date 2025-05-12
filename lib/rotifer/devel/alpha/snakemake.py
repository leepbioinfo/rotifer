from rotifer.db import ncbi
import pandas as pd

def pid2selectdipg(pid,
                   after=3,
                   before=3,
                   strand=False,
                   eukaryotes=False):
    ic = ncbi.IPGCursor()
    gnc = ncbi.GeneNeighborhoodCursor(after=after, before=before, strand=strand, eukaryotes=eukaryotes)
    ipgtable = ic.fetchall(pid)
    i1 = ipgtable.query('pid in @pid')
    missing = list(set(pid) - set(i1.pid))
    i2 = ipgtable.query('representative in @missing')
    i2.pid = i2.representative
    ipg = pd.concat([i1,i2]).drop_duplicates('pid')
    a = ncbi.assemblies(taxonomy=False)
    a.rename({'#assembly_accession': 'assembly_accession'}, axis=1, inplace=True)
    a = a.query('assembly_accession in @ipg.assembly.unique().tolist()')
    tc =   ncbi.TaxonomyCursor()
    taxonomies = tc.fetchall(a.taxid.tolist())
    a = a[['assembly_accession','taxid']]
    result = a.merge(taxonomies,
                     how='left').merge(ipg[['pid','nucleotide','start','stop', 'assembly']],
                                      left_on='assembly_accession',
                                       right_on='assembly',

                                       how='left'
                                       )
    result['position'] = result.nucleotide + ":" + result.start.astype(str) + '-' +  + result.stop.astype(str)
    return result
