#!/usr/bin/env python3

import os
import sys
import subprocess
import time


def create_logo(hmm, outpath = ''):
    '''
    Get sequence logo from a HMM.
    Uses skylign to do it
    return a png figure
    hmm => HMM file
    outpath => the output path (no file name)
    Returns a PNG figure, with name FILE.png at local path or custom path (use outpath = '')
    '''
    if outpath:
        out = os.path.join(outpath, hmm.split('/')[-1])
    else:
        out = hmm

    tries = 0
    while tries < 4:
        get_url_logo = subprocess.Popen(["curl -H 'Accept:application/json' -F file='@{0}' -F processing=hmm http://skylign.org".format(hmm)], 
                             shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        url_path = get_url_logo.communicate()[0].decode('utf-8')
        if 'http:' in url_path:
            url = url_path.split('"')[3]
            os.system("curl -H 'Accept:image/png' {0} > {1}.png 2> /dev/null".format(url,out))
            break
        else:
            time.sleep(0.5)
            if tries == 3:
                sys.stderr.write('No logo output {0}\nDue Server error\nTry again later\n '.format(hmm))
            tries +=1

if __name__ == '__main__':
    sys.stderr.write('Testing logo:\n')
    create_logo(sys.argv[1], outpath = '.')
