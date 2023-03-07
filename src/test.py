import sys
import os
sys.path.append(os.getcwd() + '/..')


from base import *

inp = r"""icntl = 7
emin(2) = 1.3e3
"""
print(Parameters().definition())
