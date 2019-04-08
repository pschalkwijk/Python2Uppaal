import numpy as np
import time
import shortuuid
import subprocess

# location of your uppaal installation
VERIFYTA = "/media/sf_Code/uppaal64-4.1.20-stratego-5/bin-Linux/verifyta.sh"

# system matrices, not yet used
a = np.array([[0,1],[-2,3]])
b = np.array([[0],[1]])

t = time.time()

from ta import MatlabAbstraction, MatlabTA
HWC = MatlabAbstraction('data/HWC_Q_M4.mat')
ST = MatlabAbstraction('data/STEERING_Q_M4.mat')
mTA_HWC = MatlabTA(HWC)
mTA_ST = MatlabTA(ST)

from ControlLoop import ControlLoop
cl1 = ControlLoop(mTA_HWC, name="cl1")
cl2 = ControlLoop(mTA_ST, name="cl2")
#
from Network import Network
net = Network(1, 5) # the network is 0.005/0.001 seconds active

# NOTE: To be more specific, this could be described with a NTGA.
# Uppaal doesn't care as long as edges are marked controllable or uncontrollable.
from ta import NTA
ntga = NTA(cl1, cl2, net)

from datetime import datetime
filename = f'ntga_{datetime.now()}'
with open(f'xml/{filename}.xml', 'w') as file:
    file.write(ntga.to_xml())

strategy_name = shortuuid.uuid()[:6]
with open(f'queries/{filename}.q', 'w') as file:
    file.write(f"strategy {strategy_name} = control: A[] not ({net.name}{net.index}.Bad)")


arg_list = [VERIFYTA, '-t', '0', '--print-strategies', 'strat', f'xml/{filename}.xml', f'queries/{filename}.q']
verify_ta = subprocess.run(arg_list, stdout=subprocess.PIPE)
print(verify_ta.stdout.decode('utf-8'))
print(f"strategy written to strat/{strategy_name}")
