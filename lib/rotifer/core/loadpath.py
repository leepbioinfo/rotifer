#!/usr/bin/env python3

from os.path import expanduser
import os
class path:
    def __init__(self):
        home = expanduser('~')
        local_conf = os.path.join(home, '.rotifer/config')
        self.localpath = ''
        try:
            self.localpath(local_conf+'path.yaml')
        except:
            pass

        self.globalpath = os.path.abspath(__file__).split('rotifer')[0]


if __name__ == '__main__':
    a = path().global_path
    print(path().localpath)
