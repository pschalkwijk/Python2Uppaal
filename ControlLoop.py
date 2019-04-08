from timedautomata import TGA
import shortuuid
import pyuppaal


class ControlLoop(TGA):
    """
    Create a CL system to use with a Network as in the paper by Dieky,
    to be exported to Uppaal
    """

    def __init__(self, nta, name='ControlLoop', initial_location='R1', sync='up'):
        super().__init__()
        self.name = name
        self.sync = sync
        self.index = shortuuid.uuid()[:6]

        # Use the sets of the nta as a base
        # the base locations R1 are urgent and don't have invariants
        # self.locations = {loc: pyuppaal.Location(urgent=True, name=loc, id=loc)
        #                   for loc in nta.base_locations}
        # self.locations.update({loc: pyuppaal.Location(invariant=nta.invariants[loc], name=loc, id=loc)
        #                        for loc in nta.sigma_locations})

        self.locations.update({f'R{location}' for location in nta.locations})
        self.actions_u.update(nta.actions)
        self.clocks.update(nta.clocks)
        self.edges = {self.uncontrollable(edge) for edge in nta.edges}
        self.invariants.update(nta.invariants)

        # TODO: create inital location

        # create early trigger locations
        self.urgent = {f"Ear{loc}" for loc in nta.locations}
        self.locations.update(self.urgent)
        # Add urgent invariant to all early locations; Or don't; yk, we
        # self.invariants.update({location: 'urgent' for location in early.keys()})

        # Add edge from Ri to Eari
        # TODO: add guard from MIET to tau_min (use tolerance better)
        self.edges.update([(f'R{location}', f'{nta.abstraction.limits[location][0]-5}>=c &&'
                                            f' {nta.abstraction.limits[location][0]} <= c', False, False, frozenset(),
                            f'Ear{location}')
                           for location in nta.locations])
        # self.sigma = 's'
        # guard_map = {}
        # for sigma, ta in nta.TA.items():
        #     for (Ri_s,_,_,_,Ri) in ta.edges:
        #         edge = (f"EarR{Ri_s}",f"R{Ri}")
        #         sigma_set = guard_map.get(edge, list())
        #         sigma_set.append(f" s=={sigma} ")
        #         guard_map[edge] = sigma_set
        #
        # self.edges.update([(Ear, '||'.join(sigmas), False, frozenset(self.actions_u.union({f'{self.sync}!'})),
        #                     frozenset(self.clocks), Ri)
        #                   for (Ear, Ri), sigmas in guard_map.items()])

        # TODO: Add a decent initial location
        self.l0 = initial_location


    def uncontrollable(self, edge):
        """
        Convert an edge from (l,g,a,c,l') -> (l,g,a_c,a_u,c,l')
        :param edge:
        :return:
        """
        (start, guard, internal_action, clocks, end) = edge
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
                props.update({'synchronisation': ','.join(actions_u)})
                props.update({'controllable': False})
            if resets:
                props.update({'assignment': ', '.join([f'{clock}=0' for clock in resets])})
            if actions_c:
                clock_assignments = props.get('assignment',False)
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

    # def generate_template(self):
    #     return pyuppaal.Template(self.name,
    #                              declaration=self.generate_declarations(),
    #                              locations=[value for value in self.locations.values()],
    #                              transitions=self.transitions,
    #                              initlocation=self.initlocation)

    def to_xml(self):
        template = self.template
        template.layout(auto_nails=True)
        return template.to_xml()
