#!/usr/bin/env python
from __future__ import division
import argparse
import sys
import os
import re
try:
    from itertools import izip as zip
except ImportError:
    pass

if len(sys.argv) <2:
    print("Please start command with a relevant tool:\n\tcallMUT\tTransloCapture")
    sys.exit()

if sys.argv[1].upper()=="CALLMUT":
    from . import callMUT as mod

elif sys.argv[1].upper()=="TRANSLOCAPTURE":
    from . import TransloCapture as mod

def main(argus):
    mod.main_process(argus)

args = mod.argypargy()

if __name__=="__main__":
    main(args)