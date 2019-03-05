"""
A traffic-model abstraction for event-triggered control systems,
as done by C. Hop based on the work of A. Kolarijani as found in
"A Formal Traffic Characterization of LTI Event-triggered Control Systems, A.S. Kolarijani et al"
"""
from etctime import TabuadaPETC, LinearPlant, LinearController, LinearETC
import numpy as np
import scipy as sc
from scipy.linalg import expm

class TrafficModel:
    """
    This is a base class for a traffic model.

    Attributes:
        plant (LinearPlant): The plant to be controlled
        controller (LinearController): The controller controlling the plant
        trigger (LinearETC): The event-triggered system
    """
    plant = None
    controller = None
    trigger = None


class TrafficModelPETC(TrafficModel):
    """
    This is a periodic-ETC version of the base traffic model
    Attributes:
        M (list): Transition matrices to state
        N (list): Transition matrices to input and output
        Q (list): Quadratic forms for the cones
    """
    M = []
    N = []
    Q = []

    def __init__(self, trigger, N_conv=5, sigma_bar=1, l=100, m=4):
        """
        Initialize the traffic model with the triggering combination

        :param TabuadaPETC trigger: Combination of plant, controller and triggering coefficient
        :param int N_conv: order of the taylor approximation of Phi. Should be larger than 2.
        :param float sigma_bar: Upper limit for the global lower bound on sampling
        :param int l: number of subdivisions in the interval [0,sigma_bar]
        :param int m: Number of considered subdivisions for the angle(s) theta (ranging from 0 to pi),
         must be an even number!
        """
        # We check if N_conv > 2
        try:
            assert (N_conv > 2)
        except AssertionError as e:
            print(f"Order of taylor approximation (N_conv) should be larger than 2, not {N_conv}")
            raise e
        # We check if m is an even number
        try:
            assert (m % 2 == 0)
        except AssertionError as e:
            print(f"Number of angle subdivisions (m) should be an even number, not {m}")
            raise e

        self.trigger = trigger
        self.plant = trigger.plant
        self.controller = trigger.controller

        # Find the state space dimension of the plant (n <= 3), and use it to calculate
        # the number of regions we use for half the state space (q)
        (n, _) = trigger.plant.A.size()

        try:
            assert (n <= 3)
        except AssertionError as e:
            print(f"State space dimensions larger than 3 are not supported")
            raise e

        q = pow(m, n-1)
        self._line_search(0, 0.001, sigma_bar)

    def _line_search(self, tau_min, d_tau, sigma_bar, nu=None):
        """
        Find tau_opt using a (constrained) line-search

        :param float tau_min: time to start the line search
        :param float d_tau: step size in the line search
        :param float sigma_bar: upper limit for the global lower bound
        :param list nu:
        :return float: return the found optimal value
        """

        # Number of steps is end - start / step size
        steps = np.floor((sigma_bar-tau_min)/d_tau)
        tau_opt = 0

        A = self.plant.A
        (n ,_) = self.plant.A.size()

        def M_fun(sigma_prima):
            eAt = lambda t: expm(A*t)
            return eAt
        return tau_min

