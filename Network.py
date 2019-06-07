import shortuuid

from ta import TGA, pyuppaal


class Network(TGA):
    """
    Class creating a communication network
    """

    def __init__(self, channels, delta, name='Network', sync='up'):
        super().__init__()
        channel_range = range(channels)
        self.delta = delta
        self.name = name
        self.index = shortuuid.uuid()[:6]
        # all locations in a network are urgent
        self.clocks = {f'c{channel}' for channel in channel_range}
        # invariant = ' && '.join([f'{clock}<={delta}' for clock in self.clocks])
        self.invariants = {self.int_to_name(i, channels): self.int_to_invariant(i, channels, delta)
                           for i in range(1, pow(2, channels))}

        # pprint.pprint(self.invariants)
        self.locations = {(self.int_to_name(chan, channels)) for chan in range(pow(2, channels))}

        self.actions_u = {f'{sync}?'}

        edges = {*[edge for i in range(pow(2, channels))
                   for k in range(channels)
                   for edge in self.create_double_edge(i, (pow(2, k) | i), channels, delta)
                   if (pow(2, k) | i) is not i]}

        # Add 'Bad state' and edges to bad afterwards
        self.locations.add('Bad')

        edges.update({(self.int_to_name(pow(2, channels)-1, channels),
                       True,
                       frozenset(),
                       frozenset(self.act_u),
                       frozenset(),
                       'Bad')})
        edges.update({('Bad',
                       True,
                       frozenset(),
                       frozenset(self.act_u),
                       frozenset(),
                       'Bad')})
        self.edges = edges
        self.l0 = self.int_to_name(0, channels)

    @staticmethod
    def int_to_name(i, n):
        """
        Convert an integer location to an bitwise On/Off location,
        where On denotes an active bit, and Off an inactive bit, including leading zeros
        i should always be smaller than pow(2,channels)

        :param int i: the number to convert
        :param int n: the number of channels (bits)
        :return str: On for 1, Off for 0
        """
        assert i < pow(2,n)
        name = ''.join(['On' if int(c) else 'Off' for c in reversed(format(i, f'0{n}b'))])
        return name

    @staticmethod
    def int_to_invariant(i, n, delta):
        assert i < pow(2, n)
        invariant = ' && '.join([f'c{index}<={delta}' for index, c in enumerate(reversed(format(i, f'0{n}b'))) if int(c)])
        return invariant

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
        actions_u = frozenset(self.act_u)
        actions_c = frozenset(self.act_c)
        clock = f'c{(lower ^ higher) .bit_length() - 1}'
        # print(clock, format(lower, f'0{channels}b'), format(higher, f'0{channels}b'))
        clockset = frozenset({clock})
        # taking an extra channel into work resets clock
        up = (self.int_to_name(lower, channels),True, False, actions_u, clockset,self.int_to_name(higher, channels))
        # channel is released if clock equals delta, no clock resets or actions
        # TODO: keep delta a variable?
        down = (self.int_to_name(higher, channels), f'{clock}=={delta}', True, False, empty, self.int_to_name(lower, channels))
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
            props.update({'controllable': False})
            transitions.append(pyuppaal.Transition(source,
                                                   target,
                                                   **props))
        return transitions

    def to_xml(self):
        template = self.template
        # template.layout(auto_nails=True)
        return template.to_xml()
