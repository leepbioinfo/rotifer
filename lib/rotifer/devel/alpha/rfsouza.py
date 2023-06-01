#!/usr/bin/env python3

import os
import sys
import Bio
import argcomplete
import rotifer.core.cli as corecli
import pandas as pd
import gspread
from oauth2client.service_account import ServiceAccountCredentials

__version__ = '0.1'
__authors__ = 'Robson Francisco de Souza'

# Subroutines
def parse_cli():
    configdir = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(sys.argv[0]))), "etc")
    jsonfile = os.path.basename(sys.argv[0]) + '.json'
    parser = corecli.parser(description = 'Multiple sequence alignment to identity or similarity matrix')
    parser.add(long_arg = '--credentials',
               short_arg = '-c',
               dest = 'credentials',
               helper = "JSON file with authorized Google Developer's Service account key",
               default = f'{os.path.join(configdir,jsonfile)}')
    parser.add(long_arg = '--outfile',
               short_arg = '-o',
               dest = 'outfile',
               helper = "Output file name",
               default = sys.stdout)
    parser.add(long_arg = '--outformat',
               short_arg = '-of',
               dest = 'outformat',
               helper = "Output file format. Could be any of the formats supported by pandas.",
               default = 'csv')
    parser.add(long_arg = '--sheet',
               short_arg = '-s',
               dest = 'sheet',
               helper = "Name of the output sheet in the Excel file.",
               default = 'excel')
    #action = corecli.action.autoopen)
    argcomplete.autocomplete(parser)
    args = parser.parse_args()
    return args

def main():
    args = parse_cli()
    gc = connect(args.credentials)
    sheet = get_sheet(gc, args.worksheet, args.sheet)
    resultado.to_csv(args.outfile, sheet_name=args.sheet, index=False)

def connect(jsonfile):
    scope = ['https://spreadsheets.google.com/feeds', 'https://www.googleapis.com/auth/drive']
    credentials = ServiceAccountCredentials.from_json_keyfile_name(jsonfile, scope)
    gc = gspread.authorize(credentials) # gspread.Client
    return gc

def get_sheet(gc, name, index):
    sheet = gc.open(name)
    sheet = sheet.get_worksheet(index).get_all_values()
    h = sheet.pop(0)
    sheet = pd.DataFrame(sheet, columns=h)
    return sheet

if __name__ == '__main__':
    main()
