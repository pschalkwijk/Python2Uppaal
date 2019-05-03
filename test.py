import numpy as np
import time
import shortuuid
import subprocess

# location of your uppaal installation
VERIFYTA = "/home/pschalkwijk/Documents/uppaal/bin-Linux/verifyta.sh"

# system matrices, not yet used
a = np.array([[0,1],[-2,3]])
b = np.array([[0],[1]])

t = time.time()

from ta import MatlabAbstraction, MatlabTA
CL1_M_20 = MatlabAbstraction('data/CL1_M_20.mat')
CL2_M_20 = MatlabAbstraction('data/CL2_M_20.mat')
mTA_CL1 = MatlabTA(CL1_M_20)
mTA_CL2 = MatlabTA(CL2_M_20)

from ControlLoop import ControlLoop
loops = []
for i in range(2):
    loops.append(ControlLoop(mTA_CL1, f'cl{i}'))

for i in range(2):
    loops.append(ControlLoop(mTA_CL2, f'cl{2+i}'))
#
from Network import Network
net = Network(1, 1) # the network is 0.005/0.001 seconds active
loops.append(net)

# NOTE: To be more specific, this could be described with a NTGA.
# Uppaal doesn't care as long as edges are marked controllable or uncontrollable.
from ta import NTA
ntga = NTA(*loops)

from datetime import datetime
filename = f'ntga_{datetime.now()}'
with open(f'xml/{filename}.xml', 'w') as file:
    file.write(ntga.to_xml())

strategy_name = shortuuid.uuid()[:6]
with open(f'queries/{filename}.q', 'w') as file:
    file.write(f"strategy {strategy_name} = control: A[] not ({net.name}{net.index}.Bad)")


arg_list = [VERIFYTA, '-u','-s','--print-strategies','strat', f'xml/{filename}.xml', f'queries/{filename}.q']
print(' '.join(arg_list))
verify_ta = subprocess.run(arg_list, stdout=subprocess.PIPE)
print(verify_ta.stdout.decode('utf-8'))
print(f"strategy written to strat/{strategy_name}")
