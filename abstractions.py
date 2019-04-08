from itertools import groupby, count
import scipy.io as sio
import numpy as np

from timedautomata import timed_automaton


@timed_automaton
class TA:
    """
    Base function for abstracting TA's. Is decorated as a TA
    """
    def __init__(self, abstraction, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.clocks = {'c'}
        self.abstraction = abstraction
        self.parse_abstraction(abstraction)
        print(abstraction)

    def parse_abstraction(self, abstraction):
        pass

    def interval_to_guard(self, interval):
        """
        Convert an interval to a guard
        :param interval: tuple
        :return: string
        """
        assert type(interval) is tuple
        assert len(interval) == 2
        lower, upper = interval
        guard = set() # FIXME: should this be a set, a single evaluation or a string?
        if lower == 0:
            if lower == upper:
                lower = upper = 1
            else:
                lower = 1
        if lower == upper:
            for clock in self.clocks:
                guard.add(f'{clock}=={lower}')
        else:
            for clock in self.clocks:
                guard.add(f'{clock}>={lower} && {clock}<={upper}')
        return guard


class ETCTimeTA(TA):
    def parse_abstraction(self, abstraction):
        super().parse_abstraction(abstraction)
        self.locations = self.transitions_to_locations(abstraction.transition)
        self.edges = self.transitions_to_edges(abstraction.transition)
        self.invariants = self.transitions_to_invariants(abstraction.transition)  # map invariants to locations

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
            """
            The PETC from etctime is modelled as such that k starts at 1, i.e. [0,2] => k=1, k=3
            """
            l = list(it)
            return (l[0]+1,l[-1]+1)

        edge_map = self.transitions_to_edgemap(transitions)

        edges = set()
        action_set = frozenset(self.actions)
        clock_set = frozenset(self.clocks)
        for (start, end), value in edge_map.items():
            # TODO: difference between delay and discrete transitions
            intervals = [as_range(g) for _,g in groupby(value, key=lambda n, c=count(): n-next(c))]
            edges.add(tuple(val for i in intervals for guard in self.interval_to_guard(i) for val
                            in [start, guard, action_set, clock_set, end]))
        return edges


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
        invariants = {location: f"c<={final_step}" for location, final_step in upper_bound.items()}
        return invariants


class MatlabAbstraction:
    def __init__(self, filename, tol=0.001):
        mat = sio.loadmat(filename, squeeze_me=True)
        self.transitions = mat.get('Reachable_regions')
        self.scale_factor = np.int64(1/tol)
        lower_limits = (mat.get('Tau_s_opt')/tol).astype(np.int64)
        upper_limits = (mat.get('Tau_s_max')/tol).astype(np.int64)
        self.regions = np.arange(1, len(self.transitions)+1)
        self.limits = {i: tuple([lower_limits[i-1], upper_limits[i-1]]) for i in self.regions}

        class trig:
            def __init__(self, sigma):
                self.sigma = sigma
        self.trigger = trig(mat.get('alpha'))


class MatlabTA(TA):
    def parse_abstraction(self, abstraction):
        super().parse_abstraction(abstraction)
        self.locations = set(abstraction.regions)
        self.edges = self.transitions_to_edges(abstraction.transitions)
        self.invariants = self.transitions_to_invariants(abstraction.limits)

    def transitions_to_edges(self, transitions):
        """
        In matlab the transitions dict is the reachable set in form:
        {start: [end1, end2, etc]}
        :param transitions:
        :return:
        """
        edges = set()
        action_set = frozenset(self.actions)
        clock_set = frozenset(self.clocks)
        for index, endpoints in enumerate(transitions):
            start = index+1
            guard_set = self.interval_to_guard(self.abstraction.limits[start])
            # Our guard_set only has a single guard by definition
            guard = guard_set.pop()
            for end in endpoints:
                edges.add(tuple(val for val in [start, guard, action_set, clock_set, end]))
        return edges

    @staticmethod
    def transitions_to_invariants(limits):
        """
        Create the mapping of invariants to locations.
        Each location is upper bounded by the final time step found in the transition table
        :param limits: dict
        :return: dict
        """
        upper_bound = {}
        for location, (lower, upper) in limits.items():
            upper_bound[location] = upper
        # FIXME: variable set of clocks? immutable, hashable table instead of dict?
        invariants = {location: f"c<={final_step}" for location, final_step in upper_bound.items()}
        return invariants
