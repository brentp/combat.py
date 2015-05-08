#!/usr/bin/env python

import pandas as pd

p = pd.read_table('py-batch.txt', index_col=0)
r = pd.read_table('r-batch.txt', index_col=0)

assert (p - r).max().max() < 1e-4