import numpy as np
import time
import shortuuid
import subprocess

# location of your uppaal installation
VERIFYTA = "/home/pschalkwijk/Documents/uppaal/bin-Linux/verifyta.sh"

# t = time.time()

# Starting with a number sometimes causes Uppaal to crash, so we use letters only
shortuuid.set_alphabet("ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz")

from ta import MatlabAbstraction, MatlabTA
CL1_M_20 = MatlabAbstraction('data/CL1_M_8.mat')
CL2_M_20 = MatlabAbstraction('data/CL2_M_8.mat')
mTA_CL1 = MatlabTA(CL1_M_20)
mTA_CL2 = MatlabTA(CL2_M_20)

from ControlLoop import ControlLoop
from Network import Network
from ta import NTA

loops = []
loops.append(ControlLoop(mTA_CL1, f'cl1'))
loops.append(ControlLoop(mTA_CL2, f'cl2'))
net = Network(1,5)
loops.append(net)
for loop in loops:
    loop.template.layout(auto_nails=True)

# NOTE: To be more specific, this could be described with a NTGA.
# Uppaal doesn't care as long as edges are marked controllable or uncontrollable.
ntga = NTA(*loops)
ntga.template.declaration = f'{ntga.template.declaration}\nint EarNum;\nint EarMax = 4;'

# We will use the following folders, and they should exist already:
# - ./xml/
# - ./queries/
# - ./strat/
# - ./results/

# Our strategy is named simple and a short id is added for uniqueness
strategy_name = f'simple_{shortuuid.uuid()[:4]}'

# Write our model to xml/{strategy_name}.xml
with open(f'xml/{strategy_name}.xml', 'w') as file:
    file.write(ntga.to_xml())

# Write a query to queries/{strategy_name}.q
with open(f'queries/{strategy_name}.q', 'w') as file:
    file.write(f"strategy {strategy_name} = control: A[] not ({net.name}{net.index}.Bad)")

# Execute verifyta with the files just created. Output target directory 'strat'
arg_list = [VERIFYTA, '-u','-s','--generate-strategy','2','--print-strategies','strat', f'xml/{strategy_name}.xml', f'queries/{strategy_name}.q']
print(' '.join(arg_list))
verify_ta = subprocess.run(arg_list, stdout=subprocess.PIPE)

# catch the results (summary) and write to file
result = verify_ta.stdout.decode('utf-8')
with open(f'results/{strategy_name}.txt', 'w+') as file:
    file.writelines([f'{strategy_name}\n'])
    file.write(result)
print(result)
print(f"strategy written to strat/{strategy_name}")
print(f"results written to results/{strategy_name}.txt")
