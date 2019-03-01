from itertools import groupby, count

from TA import TA


class AbstractedTA:
    def __init__(self, abstraction):
        self.abstraction = abstraction
        self.locations = self.transitions_to_locations(abstraction.transition)
        self.clocks = {'c'}
        self.actions = {'*'}
        self.edges = self.transitions_to_edges(abstraction.transition)
        self.invariants = self.transitions_to_invariants(abstraction.transition) # map invariants to locations

    @staticmethod
    def transitions_to_locations(transitions):
        """
        Create a list of unique locations that are actually reachable
        :param transitions: dict
        :return: set
        """
        locations = set()
        for (start, step), end in transitions.items():
            locations.add(start)
            locations.update(end)
        return locations

    @staticmethod
    def transitions_to_edgemap(transitions):
        edge_map = {}
        for (start, step), targets in transitions.items():
            for target in targets:
                if (start, target) in edge_map:
                    edge_map[(start, target)].append(step)
                    edge_map[(start, target)].sort()
                else:
                    edge_map[(start, target)] = [step]
        return edge_map

    def transitions_to_edges(self, transitions):
        """
        Create guards for a set of transitions
        :param transitions: dict
        :return: set
        """
        def as_range(it):
            l = list(it)
            return (l[0],l[-1])

        edge_map = self.transitions_to_edgemap(transitions)

        edges = set()
        action_set = frozenset(self.actions)
        clock_set = frozenset(self.clocks)
        for (start, end), value in edge_map.items():
            # TODO: difference between delay and discrete transitions
            intervals = [as_range(g) for _,g in groupby(value, key=lambda n, c=count(): n-next(c))]
            edges.add(tuple(val for i in intervals for guard in self.interval_to_guard(i) for val in [start, guard, action_set, clock_set, end]))
        return edges

    def interval_to_guard(self, interval):
        """
        Convert an interval to a guard
        :param interval: tuple
        :return: string
        """
        assert type(interval) is tuple
        assert len(interval) == 2
        guard = set() # FIXME: should this be a set, a single evaluation or a string?
        if interval[0] == interval[1]:
            for clock in self.clocks:
                guard.add('%s==%s;' % (str(interval[0]),clock))
        else:
            for clock in self.clocks:
                guard.add('%s<=%s<=%s;' % (str(interval[0]),clock, str(interval[1])))
        return guard

    @staticmethod
    def transitions_to_invariants(transitions):
        """
        Create the mapping of invariants to locations.
        Each location is upper bounded by the final time step found in the transition table
        :param transitions: dict
        :return: dict
        """
        upper_bound = {}
        for (start, step) in transitions.keys():
            upper_bound[start] = max(upper_bound.get(start, 0), step)
        # FIXME: variable set of clocks? immutable, hashable table instead of dict?
        invariants = {location: "c<=%s" % str(final_step) for location, final_step in upper_bound.items()}
        return invariants


class SigmaTA(AbstractedTA):
    def __init__(self, abstraction):
        if isinstance(abstraction, AbstractedTA):
            abstraction = abstraction.abstraction
        super(SigmaTA, self).__init__(abstraction)
        self.sigma = abstraction.trigger.sigma
