import shortuuid
import subprocess

"""
Constants and environment variables
"""
# location of your uppaal installation
VERIFYTA = "/home/pschalkwijk/Documents/uppaal/bin-Linux/verifyta.sh"

# Starting with a number sometimes causes Uppaal to crash, so we use letters only
shortuuid.set_alphabet("ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz")

run = shortuuid.uuid()[:4]
for i in range(10):
    """
    Abstractions
    """
    from ta import MatlabAbstraction, MatlabTA
    # Import the abstractions made in matlab and parse into a timed automaton
    CL1 = MatlabAbstraction(f'data/CL1_M_{20*(i+1)}.mat')
    CL2 = MatlabAbstraction(f'data/CL2_M_{20*(i+1)}.mat')
    mTA_CL1 = MatlabTA(CL1)
    mTA_CL2 = MatlabTA(CL2)

    """
    Control loops and network
    """
    from ControlLoop import ControlLoop
    from Network import Network
    from ta import NTA

    # add the contraints as used the control loops in Dieky's work
    cl1 = ControlLoop(mTA_CL1, f'cl1')
    # create a second loop form the same system with a different initial location
    cl2 = ControlLoop(mTA_CL2, f'cl2')

    # create a network as described by Dieky
    net = Network(1, 5) # the network is 0.005/0.001 seconds active

    # Use the auto-layout function to create a graph that's humanly readable
    # this uses "dot" and is slow for really large graphs. If you don't care about readability, skip this
    cl1.template.layout(auto_nails=True)
    cl2.template.layout(auto_nails=True)
    net.template.layout(auto_nails=True)

    # NOTE: To be more specific, this could be described with a NTGA.
    # Uppaal doesn't care as long as edges are marked controllable or uncontrollable.
    # Combine the two control loops and the network into a NT(G)A, and add the declaration of earlier update parameter
    ntga = NTA(cl1, cl2, net)
    ntga.template.declaration = f'{ntga.template.declaration}\nint EarNum;\nint EarMax = 4;'

    """
    Create Strategy
    """
    # we generate a unique ID for the strategy, so we don't overwrite anything
    strategy_name = f'Scalability_{run}_M_{20*(i+1)}'
    # write to a file, using the built in to_xml function
    with open(f'scalablity/{strategy_name}.xml', 'w') as file:
        file.write(ntga.to_xml())

    # create the query is our scheduler objective: control such that we don't reach network.Bad state
    with open(f'scalability/{strategy_name}.q', 'w') as file:
        file.write(f"strategy {strategy_name} = control: A[] not ({net.name}{net.index}.Bad)")

    arg_list = [VERIFYTA, '-u','-s','--generate-strategy','2','--print-strategies','scalability', f'scalability/{strategy_name}.xml', f'scalability/{strategy_name}.q']
    print(' '.join(arg_list))
    verify_ta = subprocess.run(arg_list, stdout=subprocess.PIPE)
    result = verify_ta.stdout.decode('utf-8')
    with open(f'scalability/{strategy_name}.txt', 'w+') as file:
        file.writelines([f'{strategy_name}\n'])
        file.write(result)
    print(result)
    print(f"strategy written to scalability/{strategy_name}")
    print(f"results written to scalability/{strategy_name}.txt")


    # """
    # Parse results
    # """
    # from ta import parser
    # strategy = parser()
    #
    # strat_ta = None
    # with open(f'demo/{strategy_name}', 'r') as strat:
    #     strat_ta = strategy.parse(strat.read())
    # print(strat_ta.locations)
