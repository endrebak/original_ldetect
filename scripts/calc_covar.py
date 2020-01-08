#!/usr/bin/env python3

import sys
import pandas as pd
# import numpy as np
# import pandas as pd

# calculate Wen/Stephens shrinkage LD estimate
gmapfile = sys.argv[1] # genetic map
indfile = sys.argv[2] #list of individuals
# NE = 11418.0
# CUTOFF = 1e-7
# outfile = sys.argv[5] # outfile file

NE = float(sys.argv[3])
CUTOFF = float(sys.argv[4])

from ldetect2.src.calc_covar import calc_covar

df = pd.read_table(sys.stdin, header=None, index_col=None)

# allpos = df.pop(0) #.tolist()
# allrs = df.pop(1) #.tolist()

pos2gpos = pd.read_table(gmapfile, usecols=[1, 2], index_col=0, header=None, sep=" ", squeeze=True)

# print(pos2gpos,)

# print(allpos, file=sys.stderr)
# print(allrs, file=sys.stderr)
# print(pos2gpos[allpos.values], file=sys.stderr)
# raise
# print(pos2gpos, file=sys.stderr)
# print(pos2gpos.reindex(allpos.values))
# missing = set(allpos) - set(pos2gpos)
# print(len(allpos), file=sys.stderr)
# print(len(pos2gpos), file=sys.stderr)
# print(sorted(list(missing))[:10], file=sys.stderr)
# raise

calc_covar(df, pos2gpos, indfile, NE, CUTOFF)
