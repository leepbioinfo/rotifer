#!/usr/bin/env python3
import pandas as pd

def select_col(df_col,col_ls):
    '''
    Return a list of columns.
    This function accepts a list of strings parse it to return the corresponding columns in the df.
    Use this to drop or keep columns
    ------------------------
    INPUT:
    df_col = an array containing the columns in a pandas DF (df.columns)
    col_ls = a list of strings in the format (column_number; start..end; column_name)
    ------------------------
    OUTPUT:
    cols: a list with columns name
    ------------------------
    EXAMPLES:
    df:
       A B C D E F
    a  1 0 1 1 0 1
    b  2 1 1 0 1 0

    col_ls = ['0', '1..3', 'F']
    cols = select_col(df.columns, col_ls)
    print(cols)
    >> ['A','B','C','D', 'F']

    Here the first number selected the first column (0)
    The 1..3 selected the columns 1 to 3 (columns 1,2,3)
    The last item in the list ('F') selected the column by name.
    '''
    cols = []
    for col in col_ls:
        try:
            col = int(col)
        except:
            col = col
        if isinstance(col, int):
            cols.append(df_col[col])
        else:
            if '..' in col:
                _ = col.split('..')
                start = int(_[0])
                end = int(_[1])
                for x in range(start, end +1):
                    cols.append(df_col[x])
            else:
                cols.append(col)
    return cols

def fasta2df(fi, split_seq = False):
    tmp = {}
    splitted = {}
    y = 0
    for x in range(0, len(fi)):
        if fi[x].startswith('>'):
            header = fi[x]
            header = fi[x].replace('>', '')
            tmp[header] = ''
        else:
            tmp[header] += fi[x]
            if split_seq:
                for e in fi[x]:
                    splitted[str(y)] = e
                    y+=1
    if split_seq:
        print([x.split(' ')[0] for x in tmp.keys()])
        df = pd.DataFrame({'Header':list(tmp.keys()),
                           'ID': [x.split(' ')[0] for x in tmp.keys()],
                           'Seq': list(tmp.values())})
        df2 = df['Seq'].apply(lambda x: pd.Series(list(x)))
        df2['Header'] = df['Header']
        df2['ID'] = df['ID']
        return df2
    else:
        df = pd.DataFrame({'Header':list(tmp.keys()),
                           'ID': [x.split(' ')[0] for x in tmp.keys()],
                           'Seq': list(tmp.values())})
        df = df.astype(str)
    return df
