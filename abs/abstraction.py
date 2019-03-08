"""
A traffic-model abstraction for event-triggered control systems,
as done by C. Hop based on the work of A. Kolarijani as found in
"A Formal Traffic Characterization of LTI Event-triggered Control Systems, A.S. Kolarijani et al"
"""
from etctime import TabuadaPETC, LinearPlant, LinearController, LinearETC
import numpy as np
import scipy as sc

from numpy import matmul
from numpy.linalg import matrix_power
from scipy.linalg import expm, inv
from scipy.special import factorial

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
        (n, _) = trigger.plant.A.shape

        try:
            assert (n <= 3)
        except AssertionError as e:
            print(f"State space dimensions larger than 3 are not supported")
            raise e

        q = pow(m, n-1)
        self.l = l # number of partitions
        self._line_search(0, 0.001, sigma_bar)

    def _line_search(self, tau_min, d_tau, sigma_bar, nu=None, N_conv=5):
        """
        Find tau_opt using a (constrained) line-search

        :param float tau_min: time to start the line search
        :param float d_tau: step size in the line search
        :param float sigma_bar: upper limit for the global lower bound
        :param list nu:
        :param int N_conv: order of the taylor approximation
        :return float: return the found optimal value
        """

        # Number of steps is end - start / step size
        steps = int(np.floor((sigma_bar-tau_min)/d_tau))+1
        tau_opt = 0

        A = self.plant.A
        (n, _) = self.plant.A.shape
        B = self.plant.B
        K = self.controller.K
        alpha = self.trigger.sigma

        L = self._L(A, B, K, n, N_conv, self.l, alpha, sigma_bar)

        for i in range(steps):
            tau_s = tau_min+ i*d_tau



        return 1

    @staticmethod
    def _L_sum(A, n, Nc):
        Ai = np.eye(n)
        M = np.zeros((Nc * n, Nc * n))
        Mn = np.zeros((n, Nc * n))
        f = 1
        for k in range(Nc):
            f *= (k + 1)
            M += np.kron(np.eye(Nc, k=-k), Ai / f)
            Mn[:, k * n:(k + 1) * n] = Ai / f
            Ai = Ai @ A

        return M @ Mn.T, Mn

    @staticmethod
    def _L(A, B, K, n, N_conv, l, alpha, sigma_bar):
        invA = inv(A)
        def M_fun(s):
            # FIXME: distinguis between singular and not singular matrices: can't use inv on singular matrix)
            return np.matmul(invA, (expm(A*s)-np.eye(n)))
        L = np.zeros((n*N_conv+n, n*l))
        Lsum, Afac = TrafficModelPETC._L_sum(A, n, N_conv-1)
        dsigma_prime = sigma_bar / l
        sigma_prime = 0
        In = np.eye(n)
        abk = (A - B @ K)
        for j in range(l):
            M = M_fun(sigma_prime)
            P1 = In @ M @ abk
            N = A @ M + In
            P2 = N @ abk
            Li = np.zeros(((N_conv+1)*n, n))
            Li[:n, :] = In - P1 - P1.T + (1 - alpha)*P1@P1.T
            Li[n:2*n, :] = ((1 - alpha)*P1.T)@P2 + P2.T @ ((1 - alpha)*P1.T)
            for i in range(N_conv-1):
                Li[(i+2)*n:(i+3)*n, :] = ((1 - alpha) * P1.T) @ Afac[:n, n*i:n*(i+1)] @ P2 \
                                        + P2.T @ Afac[:n, n*i:n*(i+1)].T @ ((1 - alpha) * P1) \
                                        + (1 - alpha) * P2.T @ Lsum[n*i:n*(i+1), :n] @ P2
            sigma_prime += dsigma_prime
            L[:, j*n:(j+1)*n] = Li
        return L
