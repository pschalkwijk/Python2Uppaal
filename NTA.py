from abstractions import SigmaTA
from TA import TA


class SigmaNTA(TA):
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

        # Add all the base locations to self.locations
        self.locations.add("R%s" % location for ta in self.TA.values() for location in ta.locations)
        # Add R_s locations for each sigma
        self.locations.add("R%s_s%s" % (location, sigma) for sigma, ta in self.TA.items() for location in ta.locations)
        # combine the clocks
        self.clocks.add("%s%s" %(clock, sigma) for sigma, ta in self.TA.items() for clock in ta.clocks)
        # Add all invariants index with sigma
        self.invariants.update({"R%s_s%s" % (location, sigma) : constraint for sigma, ta in self.TA.items()
                                for location, constraint in ta.invariants.items()})
        # TODO: urgent locations, constrain as 'urgent' or 'c <= 0'?
        # TODO: edge transitions will go from Ri -> Ri_sj urgently, and from Ri_sj to Ri with guards and actions from sj
        # TODO: How to parse the actions?
