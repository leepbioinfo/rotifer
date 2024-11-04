
def extract_by_hmm(seqobj1, seqobj2, suffixes=('','_YyYyYy')):
    '''
    Use a HMM derived from seqobj2 to find similar regions
    in the sequences in seqobj1.

    If more than one match is found in a protein, two 
    different sequences will be returned with the same
    identifier and can only be distinguished using the
    start and end columns.

    Returns
    -------
    Unaligned sequence set object
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import sys
    import pandas as pd
    import rotifer
    #from rotifer.core import loadpath
    
    hmmer2table = f'{rotifer.config["base"]}/bin/hmmer2table'
    domain2architecture = f'{rotifer.config["base"]}/bin/domain2architecture'
    architecture2table = f'{rotifer.config["base"]}/bin/architecture2table'

    sharedColumns = set()
    nseqobj = seqobj1.copy()
    with tempfile.TemporaryDirectory() as tmpdirname:
        nseqobj.df.sequence = nseqobj.df.sequence.str.replace('-','')
        nseqobj.to_file(f'{tmpdirname}/seqfile') 

        if isinstance(seqobj2, rotifer.devel.beta.sequence.sequence):
            seqobj2.to_file(f'{tmpdirname}/model')
            Popen(f'hmmbuild {tmpdirname}/model.hmm {tmpdirname}/model',
                  stdout=PIPE,
                  shell=True).communicate()
            model = f'{tmpdirname}/model.hmm'

        Popen(f'hmmsearch {model} {tmpdirname}/seqfile > {tmpdirname}/hhsearch_r',
              stdout=PIPE,
              shell=True).communicate()
        exe  = f'{hmmer2table} {tmpdirname}/hhsearch_r'
        exe += f' | {domain2architecture} | {architecture2table} > {tmpdirname}/hhsearch_table'
        Popen(exe, stdout=PIPE, shell=True).communicate()
        t = pd.read_csv(f'{tmpdirname}/hhsearch_table', sep="\t")
        t.rename({'ID':'id'}, axis=1, inplace=True)
        sharedColumns = set(t.columns).intersection(set(seqobj1.df.columns))

    if "start" in sharedColumns:
        start = f'start{suffixes[1]}'
    else:
        start = "start"
    if "end" in sharedColumns:
        end = f'end{suffixes[1]}'
    else:
        end = "end"

    nseqobj.df = nseqobj.df.merge(t, on='id', how='left', suffixes=suffixes)
    nseqobj.df[end] = nseqobj.df[end].fillna(nseqobj.df.length)
    nseqobj.df[start] = nseqobj.df[start].fillna(1)
    nseqobj.df.sequence = nseqobj.df.apply(lambda x: x.sequence[(int(x[start])-1):int(x[end])], axis=1)
    nseqobj = nseqobj.filter('domain == "model"')
    return (nseqobj)
