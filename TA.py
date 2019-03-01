"""
CLass modeling a TA following Dieky
* Creates locations and edges/transitions based on TrafficAbstraction.transitions
* all actions are internal (*). For use with network later on.
* On each transition the clock (c) is reset.
* a guards are created for each interval for which a transition from l -> l' is possible.
"""

class TA:
    locations = set()
    clocks = set()
    actions = set()
    edges = set()
    invariants = dict()

class TGA(TA):
    actions_c = set()
    actions_u = set()

    @property
    def actions(self):
        return self.actions_c.union(self.actions_u)