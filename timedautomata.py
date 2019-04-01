"""
CLass modeling a TA following Dieky
* Creates locations and edges/transitions based on TrafficAbstraction.transitions
* all actions are internal (*). For use with network later on.
* On each transition the clock (c) is reset.
* a guards are created for each interval for which a transition from l -> l' is possible.
"""

import pyuppaal
import logging

from functools import wraps


class ta_base:
    """
    A single timed automaton consists of:
     - a set of finite locations (L)
     - an initial location (l0)
     - a set of finite actions (Act)
     - a set of finitely many real-valued clocks (C)
     - a set of edges E, where E is a subset of L x Act x B(C) x 2^C x L
     - A mapping (Inv -> B(C)) assigning invariants to locations
    """
    def __init__(self):
        self.locations = set()
        self.l0 = ""
        self.clocks = set()
        self.actions = set()
        self.edges = set()
        self.invariants = dict()


def outdate_cache(fn):
    """
    Decorator, use before setting a property that should be in the TA
    Before setting the new value we set the class template to 'dirty'
    :param fn:
    :return:
    """
    @wraps(fn)
    def wrapped(self, *args, **kwargs):
        self._template_cached = False
        return fn(self, *args, **kwargs)

    return wrapped


def update_cache(fn):
    """
    Decorator, use to update the pyuppaal template before returning it.
    :param fn:
    :return:
    """
    @wraps(fn)
    def wrapped(self, *args, **kwargs):
        if not self._template_cached:
            self._pyuppaal = self.create_template()
            self._template_cached = True
        return fn(self, *args, **kwargs)

    return wrapped


def timed_automaton(cls):
    """
    Decorator to turn a class into a timed automaton.
    Internal value _ta is 'read-only', should only be edited from inside this decorator
    PyUppaal Template will be cached if possible, otherwise created again

    :param cls:
    :return:
    """
    class TimedAutomaton(cls):
        def __init__(self):
            super().__init__()
            self._ta = ta_base()
            self._template_cached = False
            self._pyuppaal = pyuppaal.Template(cls.__name__)

        def generate_declarations(self):
            """
            Overload this function with a more detailed variant in your TA
            :return:
            """
            return f"clock {', '.join(self.clocks)};"

        def generate_locations(self):
            """
            Overload this function with a more detailed variant in your TA
            :return:
            """
            locations = [pyuppaal.Location(invariant=self.invariants.get(loc), name=loc) for loc in self.locations]
            return locations

        def generate_transitions(self):
            """
            Overload this function with a more detailed variant in your TA
            :return:
            """
            transitions = [pyuppaal.Transition(source, target, guard=guard) for
                           (source, guard, action, select, target) in self.edges]
            return transitions

        def assign_initial_location(self, template):
            """
            Overload this function with a more detailed variant in your TA
            :return:
            """
            try:
                template.initlocation = template.get_location_by_name(self.l0)
            except AssertionError as a:
                logging.debug(f'No initial location matching {self.l0} found in current template')

        def create_template(self):
            """
            overwrite this function in with a more detailed function
            :return:
            """
            locations = self.generate_locations()
            transitions = self.generate_transitions()
            declarations = self.generate_declarations()
            template = pyuppaal.Template(self._pyuppaal.name, declaration=declarations, locations=locations,
                                         transitions=transitions)
            self.assign_initial_location(template)
            return template

        @property
        def locations(self):
            return self._ta.locations

        @locations.setter
        @outdate_cache
        def locations(self, locations):
            if len(locations) is 0:
                self._ta.locations = set()
            else:
                self._ta.locations.update(locations)

        @property
        def l0(self):
            return self._ta.l0

        @l0.setter
        @outdate_cache
        def l0(self, initial_location):
            self._ta.l0 = initial_location

        @property
        def actions(self):
            return self._ta.actions

        @actions.setter
        @outdate_cache
        def actions(self, actions):
            if len(actions) is 0:
                self._ta.actions = set()
            else:
                self._ta.actions.update(actions)

        @property
        def clocks(self):
            return self._ta.clocks

        @clocks.setter
        @outdate_cache
        def clocks(self, clocks):
            if len(clocks) is 0:
                self._ta.clocks = set()
            else:
                self._ta.clocks.update(clocks)

        @property
        def edges(self):
            return self._ta.edges

        @edges.setter
        @outdate_cache
        def edges(self, edges):
            if len(edges) is 0:
                self._ta.edges = set()
            else:
                self._ta.edges.update(edges)

        @property
        def invariants(self):
            return self._ta.invariants

        @invariants.setter
        @outdate_cache
        def invariants(self, invariants):
            if len(invariants) is 0:
                self._ta.invariants = dict()
            else:
                self._ta.invariants.update(invariants)

        @property
        def ta(self):
            return self._ta

        @property
        def name(self):
            return self._pyuppaal.name

        @name.setter
        def name(self, name):
            self._pyuppaal.name = name

        @property
        @update_cache
        def template(self):
            return self._pyuppaal

    return TimedAutomaton


def game_automaton(cls):
    """
    Decorator to turn a class into a game automaton.

    :param cls:
    :return:
    """
    class GameAutomaton(cls):
        def __init__(self):
            super().__init__()
            self.act_c = set()
            self.act_u = set()

        @property
        def actions_c(cls):
            return cls.act_c

        @actions_c.setter
        @outdate_cache
        def actions_c(cls, actions):
            if len(actions) is 0:
                cls.act_c = set()
            else:
                cls.act_c.update(actions)

        @property
        def actions_u(cls):
            return cls.act_c

        @actions_u.setter
        @outdate_cache
        def actions_u(cls, actions):
            if len(actions) is 0:
                cls.act_u = set()
            else:
                cls.act_u.update(actions)

        @property
        def actions(cls):
            return cls.actions_c.union(cls.actions_u)
    return GameAutomaton


def timed_game_automaton(cls):
    @game_automaton
    @timed_automaton
    class TimedGameAutomaton(cls):
        pass

    return TimedGameAutomaton


def network_timed_automata(cls):
    """
    Decorator to turn a class into a network of timed automata.
    Internal value _tas is 'read-only' and should only be edited from inside this decorator
    :param cls:
    :return:
    """

    base = timed_automaton(cls)

    def nta_to_ta(fn):
        def wrapper(self, properties, *args, **kwargs):
            for property in properties:
                ta, prop = property.split(".")
                if ta not in self._tas.keys():
                    self._tas[ta] = self._TA()
                    self._tas[ta].name = ta
                setattr(self._tas[ta], fn.__name__, {prop})
                # print(self._tas[ta])
            return fn(self, properties, *args, **kwargs)
        return wrapper

    class NetworkTimedAutomata(base):
        def __init__(self):
            super().__init__()
            self._tas = dict()
            self._TA = TA

        def generate_system(self):
            system_ass = ""
            systems = []
            for index, ta in enumerate(self._tas.keys()):
                system_ass += f'{ta}{index} = {ta}();\n'
                systems .append(f'{ta}{index}')
            return system_ass + f"system {', '.join(systems)};"

        def generate_declarations(self):
            if len(self.actions) is 0:
                return ""
            else:
                return f'broadcast chan {", ".join(self.actions)};'

        def create_template(self):
            return pyuppaal.NTA(templates=[ta.template for ta in self._tas.values()],
                                declaration=self.generate_declarations(),
                                system=self.generate_system())

        @property
        def locations(self):
            return base.locations.fget(self)

        @locations.setter
        @nta_to_ta
        def locations(self, locations):
            base.locations.fset(self, locations)

        @property
        def clocks(self):
            return base.clocks.fget(self)

        @clocks.setter
        @nta_to_ta
        def clocks(self, clocks):
            base.clocks.fset(self, clocks)

        @property
        def l0(self):
            return base.l0.fget(self)

        @l0.setter
        @nta_to_ta
        def l0(self, inital_locations):
            base.l0.fset(self, inital_locations)

        @property
        def actions(self):
            return base.actions.fget(self)

        @actions.setter
        @nta_to_ta
        def actions(self, actions):
            base.actions.fset(self, actions)

        @property
        def edges(self):
            return base.edges.fget(self)

        @edges.setter
        @nta_to_ta
        def edges(self, edges):
            base.edges.fset(self, edges)

        @property
        def invariants(self):
            return base.invariants.fget(self)

        @invariants.setter
        @nta_to_ta
        def invariants(self, invariants):
            base.invariants.fset(self, invariants)

        @property
        def ta(self):
            raise AttributeError(f"{self.__class__} object has no attribute 'ta'")

        @property
        def nta(self):
            return self._ta

    return NetworkTimedAutomata


def network_timed_game_automata(cls):
    base = network_timed_automata(cls)

    @game_automaton
    class NetworkTimedGameAutomata(base):
        def __init__(self):
            super().__init__()
            self._TA = TGA

    return NetworkTimedGameAutomata


@timed_automaton
class TA:
    pass


@timed_game_automaton
class TGA:
    pass


@network_timed_automata
class NTA:
    pass


@network_timed_game_automata
class NTGA:
    pass

@timed_automaton
class PTA:
    pass


@timed_game_automaton
class SPTGA:
    pass

