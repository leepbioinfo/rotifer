#!/usr/bin/env python3
from log import log
import os

aa = 'aaaa'

def f(x):
    log({1: 'This infoAAA',
              2: 'This is Warning',
              3: 'DEBUG'},
             level = 1,
        name = __name__
        )
