import numpy as np
import time
import shortuuid
import subprocess

# location of your uppaal installation
VERIFYTA = "/home/pschalkwijk/Documents/uppaal/bin-Linux/verifyta.sh"

# system matrices, not yet used
# a = np.array([[0,1],[-2,3]])
# b = np.array([[0],[1]])

# t = time.time()

# Starting with a number sometimes causes Uppaal to crash, so we use letters only
shortuuid.set_alphabet("ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz")

from ta import MatlabAbstraction, MatlabTA
CL1_M_20 = MatlabAbstraction('data/CL1_M_20.mat')
CL2_M_20 = MatlabAbstraction('data/CL2_M_20.mat')
mTA_CL1 = MatlabTA(CL1_M_20)
mTA_CL2 = MatlabTA(CL2_M_20)

from ControlLoop import ControlLoop
from Network import Network
from ta import NTA

for t in range(10):
    loops = []
    cl1 = 0
    cl2 = 0
    for i in range(t):
        if i % 2:
            loops.append(ControlLoop(mTA_CL1, f'cl{i}'))
            cl1 += 1
        else:
            loops.append(ControlLoop(mTA_CL2, f'cl{i}'))
            cl2 += 1

    #
    net = Network(1, 5) # the network is 0.005/0.001 seconds active
    loops.append(net)
    for loop in loops:
        loop.template.layout(auto_nails=True)
    # NOTE: To be more specific, this could be described with a NTGA.
    # Uppaal doesn't care as long as edges are marked controllable or uncontrollable.
    ntga = NTA(*loops)


    strategy_name = f'{shortuuid.uuid()[:4]}_CL1_{cl1}_CL2_{cl2}'
    with open(f'xml/{strategy_name}.xml', 'w') as file:
        file.write(ntga.to_xml())


    with open(f'queries/{strategy_name}.q', 'w') as file:
        file.write(f"strategy {strategy_name} = control: A[] not ({net.name}{net.index}.Bad)")

    arg_list = [VERIFYTA, '-u','-s','--generate-strategy','1','--print-strategies','strat', f'xml/{strategy_name}.xml', f'queries/{strategy_name}.q']
    print(' '.join(arg_list))
    verify_ta = subprocess.run(arg_list, stdout=subprocess.PIPE)
    result = verify_ta.stdout.decode('utf-8')
    with open(f'results/{strategy_name}.txt', 'w+') as file:
        file.writelines([f'{strategy_name}\n'])
        file.write(result)
    print(result)
    print(f"strategy written to strat/{strategy_name}")
    print(f"results written to results/{strategy_name}.txt")
