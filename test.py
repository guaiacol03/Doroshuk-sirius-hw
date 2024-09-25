import numpy as np
import pandas as pd
PADDING = b'0'
VALID_NUCLEOTIDES = np.array(["A", "T", "G", "C", "a", "t", "g", "c"], dtype="S1")

a = np.array([['a', 't', 'g', 'c', 'a', 't', 'g', 'c'], ['a', 'a', 'g', 'c', 'a', 't', 'c', 'c']], dtype='S1')
adapters = [['a', 't', 'c'], ['a', 't', 'g', 'c']]


t = a[:, :1]
print(1)

