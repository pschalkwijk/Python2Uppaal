from abstractions import SigmaTA
from TA import TA
from pyuppaal import NTA

class NTGA(NTA):
    def __init__(self, *TGA, synchronisation='up'):
        system = ''
        for index, tga in enumerate(TGA):
            system += f'{tga.name}{index} = {tga.name}();\n'
            system += f'system {tga.name}{index};\n'
        super(NTGA, self).__init__(
            declaration=f'broadcast chan {synchronisation};',
            templates=TGA,
            system=system
        )

class SigmaNTA(TA):
    base_locations = set()
    sigma_locations = set()

    def __init__(self, *ta, **kwargs):
        """
        Init network with tuple of timed-automata or abstractions
        :param ta: tuple
        :param kwargs:
        """
        self.TA = {}
        assert type(ta) == tuple
        for index, item in enumerate(ta):
            if isinstance(item, SigmaTA):
                self.TA[index] = item
            else:
                self.TA[index] = SigmaTA(item)

        # combine the clocks
        self.clocks.update(*[ta.clocks for ta in self.TA.values()])

        # Add all the base locations to self.locations
        self.base_locations.update([f"R{location}" for ta in self.TA.values() for location in ta.locations])
        # Base locations are urgent
        self.invariants.update({location: 'urgent' for location in self.locations})
        # Add a transition from Ri to Ri_sigma for each Ri, sigma assigning the index for sigma to s
        self.edges.update([(location, True, f's={sigma}', frozenset(), f'{location}_s{sigma}')
                           for location in self.locations for sigma in self.TA.keys()])
        # Add R_s locations for each sigma
        self.sigma_locations.update([f"R{location}_s{sigma}" for sigma, ta in self.TA.items() for location in ta.locations])

        # Add all invariants index with sigma
        self.invariants.update({f"R{location}_s{sigma}": constraint for sigma, ta in self.TA.items()
                                for location, constraint in ta.invariants.items()})
        # Update the edges in the TA's from (Ri -> Rj) to (Ri_sigma -> Rj)
        self.edges.update([self.edge_to_sigma_edge(edge, sigma)
                           for sigma, ta in self.TA.items() for edge in ta.edges])

    @property
    def locations(self):
        return self.base_locations.union(self.sigma_locations)

    def edge_to_sigma_edge(self, edge, sigma):
        """
        Change the starting location from i to Ri_sSigmag
        Change the target location from i to Ri
        :param sigma: int
        :param edge: (l, g, {a},{c}, l')
        :return: 
        """
        (start, guard, actions, clocks, target) = edge
        return f'R{start}_s{sigma}', guard, actions, clocks, f'R{target}'
