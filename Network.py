
from TA import TGA
import pprint
import pyuppaal
last_location_id = 0
class Network:
    """
    Class creating a communication network
    """

    edges = set()
    actions_u = set()
    actions_c = set()

    def __init__(self, channels, delta, sync='up'):
        self.channels = range(channels)
        self.delta = delta

        #self.locations = {self.int_to_name(chan, channels) for chan in range(pow(2, channels))}
        # create pyuppaal locations
        # all locations in a network are urgent
        # locations = {self.int_to_name(chan, channels) for chan in range(pow(2, channels))}
        self.locations = {loc: pyuppaal.Location(urgent=True, name=loc, id=chan)
                          for loc, chan in [(self.int_to_name(chan, channels),chan)
                                            for chan in range(pow(2, channels))]}

        self.clocks = {f'c{channel}' for channel in self.channels}
        self.actions_u = {f'{sync}?'}

        #self.edges.update({*[edge for i in range(pow(2, channels))
        edges = {*[edge for i in range(pow(2, channels))
                            for k in range(channels)
                            for edge in self.create_double_edge(i, pow(2, k) | i, channels, delta)
                            if (pow(2, k) | i) is not i]}
        # Add 'Bad state' and edges to bad afterwards
        self.locations.update({'Bad': pyuppaal.Location(urgent=True, name='Bad')})

        edges.update({(self.int_to_name(pow(2,channels)-1,channels),
                        True,
                        frozenset(self.actions_u),
                        frozenset(),
                        frozenset(),
                        'Bad')})
        edges.update({('Bad',
                        True,
                        frozenset(self.actions_u),
                        frozenset(),
                        frozenset(),
                        'Bad')})
        self.edges = edges
        self.transitions = self.generate_transitions()
        self.initlocation = self.locations[self.int_to_name(0,channels)]

    @staticmethod
    def int_to_name(i, n):
        """
        Convert an integer location to an bitwise On/Off location,
        where On denotes an active bit, and Off an inactive bit, including leading zeros
        i should always be smaller than pow(2,channels)

        :param i: int
        :return: str
        """
        assert i < pow(2,n)
        return ''.join(['On' if int(c) else 'Off' for c in format(i, f'0{n}b')])

    def create_double_edge(self, lower, higher, channels, delta):
        """
        Create a transition (lower, guard, contr, uncontr, resets, higher)
        and a transition (higher, guard, contr, uncontr, resets, lower)
        :param lower:
        :param higher:
        :param channels:
        :return:
        """
        # freeze sets for hashing
        empty = frozenset()
        actions_u = frozenset(self.actions_u)
        actions_c = frozenset(self.actions_c)
        clock = f'c{(lower ^ higher) .bit_length() - 1}'
        clockset = frozenset({clock})
        # taking an extra channel into work resets clock
        up = (self.int_to_name(lower, channels),True, False, actions_u, clockset,self.int_to_name(higher, channels))
        # channel is released if clock equals delta, no clock resets or actions
        # TODO: keep delta a variable?
        down = (self.int_to_name(higher, channels),f'{clock}=={delta}',True,False,empty, self.int_to_name(lower, channels))
        return up, down

    def generate_transitions(self):
        transitions = []
        for (source, guard, actions_c, actions_u, resets, target) in self.edges:
            props = {}
            if guard:
                props.update({'guard': str(guard).lower() if type(guard) is bool else guard})
            if actions_u:
                props.update({'synchronisation': ','.join(actions_u)})
            if resets:
                props.update({'assignment': ','.join([f'{clock}=0' for clock in resets])})

            transitions.append(pyuppaal.Transition(self.locations[source],
                                                   self.locations[target],
                                                   **props))
        return transitions
