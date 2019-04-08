import sympy
from sympy.matrices.expressions.blockmatrix import block_collapse
from sympy.matrices.expressions.matexpr import matrix_symbols
from sympy.utilities.autowrap import ufuncify

import numpy as np
import time
from numba import jit

# @jit

a = np.array([[0,1],[-2,3]])
t = time.time()

from abstractions import MatlabAbstraction, MatlabTA

HWC = MatlabAbstraction('Abstraction Data/HWC_Q_M4.mat')
ST = MatlabAbstraction('Abstraction Data/STEERING_Q_M4.mat')
mTA_HWC = MatlabTA(HWC)
mTA_ST = MatlabTA(ST)

from ControlLoop import ControlLoop
cl1 = ControlLoop(mTA_HWC, name="cl1")
cl2 = ControlLoop(mTA_ST, name="cl2")
#
from Network import Network
net = Network(1, 5)


from timedautomata import NTA
ntga = NTA(cl1, cl2, net)

from datetime import datetime
with open(f'xml/ntga_{datetime.now()}.xml', 'w') as file:
    file.write(ntga.to_xml())