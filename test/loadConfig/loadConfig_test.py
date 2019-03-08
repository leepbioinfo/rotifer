#!/usr/bin/env python3

import rotifer.core.functions as rcf
import os
import sys

# Create examples
a = rcf.loadConfig(':dir1.dir2.test.key1', user_path = os.getcwd())
