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
for i in range(1):
    loops.append(ControlLoop(mTA_CL1, f'cl{i}'))

for i in range(1):
    loops.append(ControlLoop(mTA_CL2, f'cl{2+i}'))
#
from Network import Network
net = Network(1, 5) # the network is 0.005/0.001 seconds active
loops.append(net)

# NOTE: To be more specific, this could be described with a NTGA.
# Uppaal doesn't care as long as edges are marked controllable or uncontrollable.
from ta import NTA
ntga = NTA(*loops)

from datetime import datetime
strategy_name = shortuuid.uuid()[:6]
filename = f'ntga_{datetime.now()}'
with open(f'xml/{strategy_name}.xml', 'w') as file:
    file.write(ntga.to_xml())


with open(f'queries/{strategy_name}.q', 'w') as file:
    file.write(f"strategy {strategy_name} = control: A[] not ({net.name}{net.index}.Bad)")


arg_list = [VERIFYTA, '-u','-s','--generate-strategy','2','--print-strategies','strat', f'xml/{strategy_name}.xml', f'queries/{strategy_name}.q']
print(' '.join(arg_list))
verify_ta = subprocess.run(arg_list, stdout=subprocess.PIPE)
print(f"strategy written to strat/{strategy_name}")
result = verify_ta.stdout.decode('utf-8')
with open(f'results/{strategy_name}.txt', 'w') as file:
    file.write(result)
print(result)
print(f"strategy written to results/{strategy_name}.txt")
