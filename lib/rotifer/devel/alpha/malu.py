def findFilesByPrefixSuffix(
        prefix,
        suffixes=[".scan.arch","_cluster.tsv"],
        path=".",
        order=['c100i100','c80i70','c80i0','pfam','aravind']):
    import os
    from glob import glob
    import pandas as pd

    # Finding...
    tables = []
    for fname in glob(f'''{path}/{prefix}.*'''):
        name = os.path.basename(fname)
        name = name[len(prefix):]
        if name[0] == ".":
            name = name[1:]
        for pattern in suffixes:
            if name[-len(pattern):] == pattern:
                name = name[:len(name)-len(pattern)]
                if name[-1] == ".":
                    name = name[:-1]
                tables.append([ name, prefix, pattern, fname ])
    tables = pd.DataFrame(tables, columns=['target','prefix','suffix','path'])
    if order:
        myorder = { order[i]: i for i in range(0,len(order)) }
        tables['idx'] = tables['target'].map(myorder)
        tables.sort_values('idx', inplace=True)
    else:
        tables['idx'] = tables.index.tolist()
    return tables

def read_tables(
        prefix,
        suffixes=[".scan.arch","_cluster.tsv"],
        path=".",
        order=['c100i100','c80i70','c80i0','pfam','aravind'],
        concat=[],
        source=[],
        merge=[],
        colnames={'c100i100':['c100i100','pid'], 'c80i70':['c80i70','c100i100'], 'c80i0':['c80i0','c80i70']},
        rename={'pfam':{'ID':'c100i100','architecture':'pfam'}, 'aravind':{'ID':'c100i100','architecture':'aravind'}},
        filter_columns={'pfam':['c100i100','pfam'], 'aravind':['c100i100','aravind']},
        sep="\t"):
    """
    This function loads, process and either join or concatenate a set of tables.

    The list of tables is derived iby matching prefixes and suffixes to a target directory.
    """
    import pandas as pd

    # Loading...
    mergeIndex = 0
    concatIndex = 0
    df = pd.DataFrame()
    tables = findFilesByPrefixSuffix(prefix, suffixes, path, order)
    print('this is the findfilesfunc result')
    print(tables)
    for idx, row in tables.iterrows():
        print(f"this is {idx}")
        name = row['target']
        print(name)

        # Set column names
        columns = None
        if name in colnames:
            columns = colnames[name]
            print(columns)
        # Read
        tmp = pd.read_csv(row['path'], sep=sep, names=columns)
        print('this is before adjust loop')
        print(tmp.head(2))

        # Adjust
        if name in rename:
            tmp.rename(rename[name], inplace=True,axis=1)
        if name in filter_columns:
            tmp = tmp.filter(filter_columns[name])
        if name in source:
            tmp['source'] = name
        
        print('this is after passing adjust loop')
        print(tmp.head(2))
        # Operate
        print('Reaching Operate loop')
        print(f"df.empty: {df.empty}")
        if df.empty:
            print(f"Setting df to: {name}")
            df = tmp
            print('this is printed if df was empty')
            print(df.head(2))
        else:
            print(f"Trying to merge or concat: {name}")
            if merge and merge[mergeIndex] == name:
                print('merge loop')
                print(f'{merge} and {merge[mergeIndex]}')
                df = df.merge(tmp, how='left')
                mergeIndex = mergeIndex + 1
                print(df.head(2))
            elif concat and concat[concatIndex] == name:
                print('concat loop')
                print(f'{concat} and {concat[concatIndex]}')
                df = pd.concat([ df, tmp ],ignore_index=True)
                concatIndex = concatIndex + 1
                print(df.head(2))
        print(f"{name} was not merged or concatenated. Skipping.")

    # return dataframe
    return df

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
