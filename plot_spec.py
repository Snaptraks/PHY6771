import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pdb

specData=open('NWKZKA.spr')

for line in specData:

    lines=line.strip()
    columns=line.split()

    print columns[0],columns[1]
