from parsimonious.grammar import Grammar
from parsimonious.nodes import NodeVisitor
from .timedautomata import network_timed_game_automata

TiGaGrammar = Grammar(
    r"""
    strategy        = in_state states newline
    
    in_state        = in_open locations newline? in_condition newline in_close
    in_condition    = "(" space_del_text+ ")"
    in_open         = "Initial state:" newline
    in_close        = "Note: The 'strategy' is not guaranteed to be a strategy." newline+ in_text
    in_text         = "Strategy to win:" / "Strategy to avoid losing"
    
    states          = state+
    state           = st_open locations delay* move*
    st_open         = newline+ "State:" ws*
    
    location        =  ws text+
    loc_addendum    = ws* text* ws*
    locations       = "(" location+ ws ")" loc_addendum
    
    delays          = delay+
    delay           = newline* del_open invariants del_close
    del_open        = "While you are in" tab
    del_close       = ", wait."
    
    moves           = move+
    move            = newline* mv_open invariants mv_close transition
    mv_open         = "When you are in" ws
    mv_close        = ","
       
    transition      = tr_open trans conditions
    tr_open         = ws* "take transition" ws*
    tr_close        = ws
    trans           = start to end
    to              = "->"
    start           = text+
    end             = text+
    
    cond            = ","? ws space_del_text+
    conditions      = ws "{" cond+ "}"
           
    inv             = "&&"? ws? inv_text+ ws?
    inv_text        = ~"[A-z0-9.=><\-+']*"i
    invariant       = " || "? "(" inv+ ")"
    invariants      = invariant+
    
    tab             = ~"\t"
    ws              = ~"\s*"i
    newline         = ~"\n*"
    text            = ~"[A-z0-9.&=><:#']*"i
    space_del_text  = ws* text+
    """
)


@network_timed_game_automata
class parser(NodeVisitor):
    grammar = TiGaGrammar

    def visit_strategy(self, node, visited_children):
        """
        strategy = in_state states newline

        :param node:
        :param visited_children:
        :return:
        """
        in_state, states, newline = visited_children
        return self.nta

    def visit_state(self, node, visited_children):
        """
        state = st_open locations delays? moves?
        :param node:
        :param visited_children:
        :return:
        """
        st_open, locations, delays, moves = visited_children
        self.locations.add(locations)
        self.invariants.update({locations: delays})
        for move in moves:
            source, guards, target = move
            print(source, guards, target)
        #add edges, (start, end,
        # state = {name: child for (name, child) in visited_children if child}
        return {"state": locations}

    def visit_location(self, node, visited_children):
        """
        location =  ws text+
        :param node:
        :param visited_children:
        :return:
        """
        ws, text = visited_children
        if len(text) > 0:
            return "".join(text)


    def visit_locations(self, node, visited_children):
        """
        locations = "(" location+ " )" loc_addendum
        :param node:
        :param visited_children:
        :return:
        """
        _, locations, _, ws, addendum = visited_children
        return frozenset(location for location in locations if location)

    def visit_invariant(self, node, visited_children):
        """
        invariant = " || "? "(" inv+ ")"
        :param node:
        :param visited_children:
        :return:
        """
        _or, _, invariant, _ = visited_children
        # print(invariant)
        return frozenset(inv for inv in invariant if len(inv) > 0)

    def visit_move(self, node, visited_children):
        """
        move = newline* mv_open invariants mv_close transition
        :param node:
        :param visited_children:
        :return:
        """
        newline, mv_open, guards, close, transition = visited_children
        source, target = transition
        return source, guards, target

    def visit_delay(self, node, visited_children):
        """
        delay = newline* del_open invariants del_close
        :param node:
        :param visited_children:
        :return:
        """
        newline, del_open, invariants, del_close = visited_children
        return set(invariants)

    def visit_transition(self, node, visited_children):
        """
        transition = tr_open trans conditions
        :param node:
        :param visited_children:
        :return:
        """
        tr_open, trans, conditions = visited_children
        return trans

    def visit_trans(self, node, visited_children):
        """
        trans = start to end
        :param node:
        :param visited_children:
        :return:
        """
        start, to, end = visited_children
        return start, end

    def visit_start(self, node, *args):
        """
        start = text+
        :param node:
        :param visited_children:
        :return:
        """
        return node.text

    def visit_end(self, node, *args):
        """
        end = text+
        :param node:
        :param visited_children:
        :return:
        """
        return node.text

    def visit_in_state(self, node, visited_children):
        """
        in_state = in_open locations newline? in_condition newline in_close
        :param node:
        :param visited_children:
        :return:
        """
        in_open, locations, newline, in_condition, newline, in_close = visited_children
        return visited_children

    def visit_in_condition(self, node, visited_children):
        """
        in_condition = "(" space_del_text+ ")"

        :param node:
        :param visited_children:
        :return:
        """
        _, initial_condition, _ = visited_children
        return node.expr_name, node.text

    def visit_inv(self, node, visited_children):
        """
        inv = "&&"? ~"[A-z0-9.=><\-+\s']*"i

        :param node:
        :param args:
        :return:
        """
        _and, ws, inv, ws = visited_children
        if len(inv) > 0:
            return "".join(inv)

    def visit_text(self, node, *args):
        return node.text

    def visit_inv_text(self, node, *args):
        return node.text

    def visit_space_del_text(self, node, *args):
        """
        space_del_text = ws* text+
        :param node:
        :param visited_children:
        :return:
        """
        return node.text

    def generic_visit(self, node, visited_children):
        # print(f'skipping {node.expr_name}')
        return visited_children

#
# with open('strat/test', 'r') as strat:
#         SV = TiGaParser()
#         SV.parse(strat.read())
