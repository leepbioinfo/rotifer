#!/usr/bin/env python3

import sys
import pandas as pd
from collections import OrderedDict

def pad(data, side='left', fillchar=' '):
    data = data.astype(str)
    if isinstance(data,pd.DataFrame):
        side = { x: side[x] if x in side else 'left' for x in data.columns }
        data = data.apply(lambda x: x.str.pad(x.str.len().max(), side=side[x.name], fillchar=fillchar))
    elif isinstance(data,pd.Series):
        data = data.str.pad(data.str.len().max(), side=side, fillchar=fillchar)
    return data

def to_color(data, mode='fg', colors=None, padding=None, fillchar=' '):
    """
    Convert a Series to DType str and colorize it.
    """
    import yaml
    from rotifer.core import functions as rcf

    if colors == None:
        colors = yaml.load(open(rcf.findDataFiles(":colors/aminoacids.yml")[0]), Loader=yaml.Loader)

    def color_res(s, cs):
        if s in colors:
            return cs(s, colors[s])
        else:
            return cs(s, "000")

    def color_bg(s, color = ''):
        '''
        s: String
        '''
        if color == "000":
            color = f'49;5;000'
        else:
            color = f'48;5;{color}'
        return f'\033[{color}m{s}\033[0m'

    def color_fg(s, color = ''):
        if color == "000":
            color = f'39;5;000'
        else:
            color = f'38;5;{color}'
        return f'\033[{color}m{s}\033[m'

    # Choose the coloring function
    color_switch = {'background':color_bg, 'bg':color_bg, 'foreground':color_fg, 'fg':color_fg}

    # Color the data
    if padding:
        data = pad(data, side=padding, fillchar=fillchar)
    if isinstance(data,pd.DataFrame):
        data = data.applymap(lambda x: ''.join([color_res(y, color_switch[mode]) for y in x]))
    elif isinstance(data,pd.Series):
        data = data.map(lambda x: ''.join([color_res(y, color_switch[mode]) for y in x]))

    return data

def optimize_memory_usage(df):
    '''
    Optimize memory usage of a Pandas DataFrame.
    
    Usage:
      from rotifer.pandas.functions import optimize_df
      df = optimize_memory_usage(df)

    Returns:
      A memory optimized version of the input Pandas
      DataFrame with the exact same content as the input.

    Parameters:
      df : a Pandas DataFrame
    '''

    # Make a copy of the input DataFrame
    df = df.copy()

    # int optimization
    converted_int = df.select_dtypes(include=['int'])
    converted_int = converted_int.apply(pd.to_numeric, downcast = "integer")

    # Float
    converted_float = df.select_dtypes(include=['float'])
    converted_float = converted_float.apply(pd.to_numeric,downcast='float')

    # Object
    df_obj = df.select_dtypes(include=['object'])
    converted_obj = pd.DataFrame()

    for col in df_obj.columns:
        num_unique_values = len(df_obj[col].unique())
        num_total_values = len(df_obj[col])
        if num_unique_values / num_total_values <= 0.5:
            converted_obj.loc[:,col] = df_obj[col].astype('category')
        else:
            converted_obj.loc[:,col] = df_obj[col]

    # Adding converted objects to a better df
    df[converted_int.columns] = converted_int
    df[converted_float.columns] = converted_float
    df[converted_obj.columns] = converted_obj
    return df

def print_everything(max_rows=10000000000, max_columns=1000000000, max_colwidth=1000000000, width=100000, verbose=False):
    """
    A function to set pandas display options to print entire
    DataFrames instead of truncating output to a few columns
    and rows.

    Usage:
      import rotifer.pandas.functions as rpf
      rpf.print_everything()

    Parameters:
      max_rows     : maximum number of DataFrame rows to print
      max_columns  : maximum number of DataFrame columns to show
      max_colwidth : maximum length of columns
                     All columns wider than this will be truncated
      width        : maximum length for rows
      verbose      : warn on STDERR when changing options
    
    WARNING:
      Very large Pandas DataFrames could take too long to print!!!
    """

    pd.options.display.max_rows = max_rows
    pd.options.display.max_columns = max_columns
    pd.options.display.max_colwidth = max_colwidth
    pd.options.display.width = width

    if verbose:
        print('Pandas options have been reset!', file=sys.stderr)

def show_display_options():
    '''
    List all Pandas display options that control printed output
    formatting (see pandas.options.display).
    '''
    return pd.DataFrame([ (x,getattr(pd.options.display,x)) for x in dir(pd.options.display) ], columns=('option','value'))
