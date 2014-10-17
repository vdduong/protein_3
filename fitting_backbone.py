# fitting a triplet of library residues into the backbone

import os
import re
import random
from routines import *
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D
from collections import deque
from auction_algo import *




dict_cs = {}
dict_cs['H_2'] = 8.57 # 0.02
dict_cs['HB2_2'] = 2.21
dict_cs['HA_1'] = 4.16
dict_cs['QB_1'] = 1.33
dict_cs['QG2_41'] = 0.85
dict_cs['QE_44'] = 7.41

dict_cs['H_3'] = 8.72
dict_cs['HG2_3'] = 2.1
dict_cs['HE22_2'] = 6.84
dict_cs['HA_2'] = 4.52
dict_cs['HA_3'] = 5.35
dict_cs['HB3_3'] = 1.88
dict_cs['HB_68'] = 1.69
dict_cs['HB2_3'] = 1.96
dict_cs['HE21_2'] = 7.95
dict_cs['QG2_67'] = 0.78

dict_cs['H_4'] = 9.15
dict_cs['QG1_41'] = 0.92
dict_cs['QB_40'] = 1.05
dict_cs['HB3_4'] = 2.77
dict_cs['HA_41'] = 4.89
dict_cs['HB2_4'] = 2.55
dict_cs['HA_39'] = 5.62
dict_cs['QG2_42'] = 0.74
dict_cs['HA_4'] = 5.16
dict_cs['QD_4'] = 7.07

# 2: GLN, 3: GLU, 4: PHE

