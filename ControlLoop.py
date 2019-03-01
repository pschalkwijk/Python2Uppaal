from TA import TA

class ControlLoop(TGA):
    """
    Create a CL system to use with a Network as in the paper by Dieky
    """
    def __init__(self, nta, sync='up!'):
        # Use the sets of the nta as a base
        self.locations.update(nta.locations)
        self.actions_u.update(nta.actions)
        self.clocks.update(nta.clocks)
        self.edges.update(nta.edges)
        self.invariants.update(nta.invariants)

        # create early trigger locations
        early = {"Ear%s_s%s" % (location, sigma) for sigma, ta in nta.TA.items() for location in ta.locations}
        self.locations.update(early)
        # Add urgent invariant to all early locations
        self.invariants.update({location: 'urgent' for location in early})
        # TODO: create edges for early updates, with sync as uncontrollable action
        # TODO: add sync as uncontrollable action to edges (discrete transitions only)