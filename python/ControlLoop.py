import shortuuid

from ta import TGA, pyuppaal


class ControlLoop(TGA):
    """
    Create a CL system to use with a Network as in the paper by Dieky,
    to be exported to Uppaal (single sigma)
    """

    def __init__(self, nta, name='ControlLoop', initial_location=None, sync='up', d=5):
        super().__init__()
        self.name = name
        self.sync = sync
        self.index = shortuuid.uuid()[:6]

        self.locations.update({f'R{location}' for location in nta.locations})
        self.actions_u.update(nta.actions)
        self.clocks.update(nta.clocks)
        self.edges = {self.uncontrollable(edge) for edge in nta.edges}
        self.edges.update({self.early(edge) for edge in nta.edges})
        self.invariants.update(nta.invariants)

        # create inital location

        # create early trigger locations
        self.urgent = {f"Ear{loc}" for loc in nta.locations}
        self.locations.update(self.urgent)
        # Add urgent invariant to all early locations; Or don't; yk, we
        # self.invariants.update({location: 'urgent' for location in early.keys()})

        # Add edge from Ri to Eari
        self.edges.update([(f'R{location}', f'c < {nta.abstraction.limits[location][0]} &&'
                                            f'{nta.abstraction.limits[location][0]} - {d} <= c && EarNum < EarMax', False, False, frozenset(),
                            f'Ear{location}')
                           for location in nta.locations])

        if initial_location is None:
            self.locations.update({f'R0'})
            self.edges.update([(f'R0', False, False, False, frozenset(),f'R{location}') for location in nta.locations])
        else:
            self.locations.update({f'R0'})
            self.edges.update([(f'R0', False, False, False, frozenset(),f'R{location}') for location in initial_location])
        self.urgent.update({"R0"})
        self.l0 = 'R0'

    def early(self, edge):
        """
        Convert an edge from (l,g,a,c,l') -> (Ear(l),g,a_c,a_u,c,l')
        :param edge:
        :return:
        """
        (start, guard, assignment, clocks, end) = edge
        ia = set(assignment)
        ia.update({'EarNum = EarNum + 1'})
        assignment = frozenset(ia)
        return f'Ear{start}', guard, assignment, frozenset({f'{self.sync}!'}), clocks, f'R{end}'

    def uncontrollable(self, edge):
        """
        Convert an edge from (l,g,a,c,l') -> (l,g,a_c,a_u,c,l')
        :param edge:
        :return:
        """
        (start, guard, internal_action, clocks, end) = edge
        ia = set(internal_action)
        ia.update({'EarNum = 0'})
        internal_action = frozenset(ia)
        return f'R{start}', guard, internal_action, frozenset({f'{self.sync}!'}), clocks, f'R{end}'

    def controllable(self, edge):
        """
        Convert an edge from (l,g,a,c,l') -> (l,g,a_c,a_u,c,l')
        :param edge:
        :return:
        """
        (start, guard, internal_action, clocks, end) = edge
        return f'R{start}', guard, internal_action, False, clocks, f'R{end}'

    def generate_locations(self):
        locations = [pyuppaal.Location(invariant=self.invariants.get(loc), name=loc, urgent=(loc in self.urgent))
                     for loc in self.locations]
        return locations

    def generate_transitions(self):
        """ convert edges to transitions """
        transitions = []
        for (source, guard, actions_c, actions_u, resets, target) in self.edges:
            props = {}
            if guard:
                props.update({'guard': str(guard).lower() if type(guard) is bool else guard})
            if actions_u:
                props.update({'synchronisation': '\n'.join(actions_u)})
                props.update({'controllable': False})
            if resets:
                props.update({'assignment': ', '.join([f'{clock}=0' for clock in resets])})
            if actions_c:
                clock_assignments = props.get('assignment', False)
                action_assignments = ', '.join([action for action in actions_c])
                if clock_assignments:
                    action_assignments = ', '.join([clock_assignments, action_assignments])
                props.update({'assignment': action_assignments})
            if target not in self.locations:
                print(f'{target} is not found in set of locations')
            else:
                transitions.append(pyuppaal.Transition(source, target, **props))
        return transitions

    def generate_clocks(self):
        return f"clock {', '.join(self.clocks)};"

    def generate_declarations(self):
        return f'{self.generate_clocks()}\n'

    def to_xml(self, layout=False):
        template = self.template
        if layout:
            template.layout(auto_nails=True)
        return template.to_xml()


class sigmaControlLoop(TGA):
    def __init__(self, *ntas, name='ControlLoop', initial_location=None, sync='up', d=5):
        super().__init__()
        self.name = name
        self.sync = sync
        self.index = shortuuid.uuid()[:6]
        self.ntas = ntas
        for index, nta in enumerate(ntas):
            self.locations.update({f'R{location}_s{index}' for location in nta.locations})
            self.actions_u.update(nta.actions)
            self.clocks.update(nta.clocks)
            self.edges = {self.uncontrollable(edge, index) for edge in nta.edges}
            self.edges.update({self.early(edge, index) for edge in nta.edges})
            self.invariants.update({f'{key}_s{index}': value} for (key, value) in nta.invariants.items())

            # create early trigger locations and make choosing sigma urgent
            self.urgent = {f"Ear{loc}_s{index}" for loc in nta.locations}
            self.urgent.update({f'R{location}' for location in nta.locations})
            self.locations.update(self.urgent)

            # Add edge from Ri to Eari
            self.edges.update([(f'R{location}_s{index}', f'c < {nta.abstraction.limits[location][0]} &&'
                                                f'{nta.abstraction.limits[location][0]} - {d} <= c && EarNum < EarMax', False, False, frozenset(),
                                f'Ear{location}_s{index}')
                               for location in nta.locations])
            self.edges.update([(f'R{location}',True,False,False,frozenset(),f'R{location}_s{index}')]
                              for location in nta.locations)
            if initial_location is None:
                self.locations.update({f'R0'})
                self.edges.update([(f'R0', False, False, False, frozenset(),f'R{location}')
                                   for location in nta.locations])
            else:
                self.locations.update({f'R0'})
                self.edges.update([(f'R0', False, False, False, frozenset(),f'R{location}')
                                   for location in initial_location])
        self.urgent.update({"R0"})
        self.l0 = 'R0'

    def early(self, edge, index):
        """
        Convert an edge from (l,g,a,c,l') -> (Ear(l),g,a_c,a_u,c,l')
        :param edge:
        :return:
        """
        (start, guard, assignment, clocks, end) = edge
        ia = set(assignment)
        ia.update({'EarNum = EarNum + 1'})
        assignment = frozenset(ia)
        return f'Ear{start}_s{index}', guard, assignment, frozenset({f'{self.sync}!'}), clocks, f'R{end}'

    def uncontrollable(self, edge, index):
        """
        Convert an edge from (l,g,a,c,l') -> (l,g,a_c,a_u,c,l')
        :param edge:
        :return:
        """
        (start, guard, internal_action, clocks, end) = edge
        ia = set(internal_action)
        ia.update({'EarNum = 0'})
        internal_action = frozenset(ia)
        return f'R{start}', guard, internal_action, frozenset({f'{self.sync}!'}), clocks, f'R{end}'

    def controllable(self, edge):
        """
        Convert an edge from (l,g,a,c,l') -> (l,g,a_c,a_u,c,l')
        :param edge:
        :return:
        """
        (start, guard, internal_action, clocks, end) = edge
        return f'R{start}', guard, internal_action, False, clocks, f'R{end}'

    def generate_locations(self):
        locations = [pyuppaal.Location(invariant=self.invariants.get(loc), name=loc, urgent=(loc in self.urgent))
                     for loc in self.locations]
        return locations

    def generate_transitions(self):
        """ convert edges to transitions """
        transitions = []
        for (source, guard, actions_c, actions_u, resets, target) in self.edges:
            props = {}
            if guard:
                props.update({'guard': str(guard).lower() if type(guard) is bool else guard})
            if actions_u:
                props.update({'synchronisation': '\n'.join(actions_u)})
                props.update({'controllable': False})
            if resets:
                props.update({'assignment': ', '.join([f'{clock}=0' for clock in resets])})
            if actions_c:
                clock_assignments = props.get('assignment', False)
                action_assignments = ', '.join([action for action in actions_c])
                if clock_assignments:
                    action_assignments = ', '.join([clock_assignments, action_assignments])
                props.update({'assignment': action_assignments})
            if target not in self.locations:
                print(f'{target} is not found in set of locations')
            else:
                transitions.append(pyuppaal.Transition(source, target, **props))
        return transitions

    def generate_clocks(self):
        return f"clock {', '.join(self.clocks)};"

    def generate_declarations(self):
        return f'{self.generate_clocks()}\n'

    def to_xml(self, layout=False):
        template = self.template
        if layout:
            template.layout(auto_nails=True)
        return template.to_xml()
