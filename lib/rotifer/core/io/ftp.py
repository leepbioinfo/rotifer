#!/usr/bin/env python3

from ftplib import FTP
from os.path import expanduser
import yaml
import os
import sys

__version__ = 0.1
__authors__ = 'Robson F. Souza'

def open(url = 'ftp.ncbi.nlm.nih.gov', user='', password = 'nicastro@iq.usp.br', unpack=False):
    '''
    Fetch a filehandle to a file located on an FTP site.
    '''

    ftp = FTP(url)
    ftp.login(user = user, passwd = password)
    for x,y in genomes.iterrows():
        aa = y.assembly_accession
        path = y.ftp_path.split('nlm.nih.gov')[1]
        genome_file = path.split('/')[-1]
        ftp.cwd(path)
        ftp.retrbinary("RETR "+ f'{genome_file}_protein.faa.gz', open(f'{genome_file}.faa.gz', 'wb').write)
        ftp.retrbinary("RETR "+ f'{genome_file}_feature_table.txt.gz', open(f'{genome_file}_feature_table.txt.gz', 'wb').write)
