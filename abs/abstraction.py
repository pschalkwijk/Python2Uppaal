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
        self.phi = None

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
        print(self._line_search(0, 0.001, sigma_bar))

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
        self.L = L
        for i in range(steps):
            tau_s = tau_min+i*d_tau
            # try:
            phi, max_eig = self._phi(tau_s, self.l, sigma_bar, N_conv, n, L)
            if max_eig > 0:
                break
            self.phi = phi
            tau_opt = tau_s
        return tau_opt

    @staticmethod
    def _phi(tau_s, l, sigma_bar, N_conv, n, L):
        dsigma_prime = sigma_bar/l
        i_max = int(np.floor(tau_s/dsigma_prime))
        sigma_prime = dsigma_prime**np.arange(0, N_conv+1)
        sigma_prime = sigma_prime.reshape(1, N_conv+1, 1, 1)
        phi = np.zeros((i_max+1, N_conv+1, n, n))
        for i in range(i_max):
            phi[i, :, :, :] = (sigma_prime * L[i, :, :, :]).cumsum(axis=1)

        if tau_s < sigma_bar:
            sigma_prime = (tau_s - i_max*dsigma_prime)**np.arange(N_conv+1)
            sigma_prime = sigma_prime.reshape(1, N_conv+1, 1, 1)
            phi[i_max, :, :, :] = (sigma_prime * L[i_max, :, :, :]).cumsum(axis=1)
        eig, eiv = np.linalg.eig(phi)
        max_eig = np.amax(eig.reshape(np.prod(eig.shape)))

        return phi, max_eig

    @staticmethod
    def _L_sum(a, n, Nconv):
        Nc = Nconv - 1
        k = np.zeros((Nc, n, n))
        m = np.zeros((Nc, n, n))
        b = np.eye(n) @ a
        k[0, :, :] = np.eye(n)
        f = 2
        for i in range(n, Nconv):
            m[i - n, :, :] = b / f
            t1 = m[:i - n, :, :]
            t2 = np.transpose(m[:i - n, :, :][::-1], (0, 2, 1))
            L1 = (b + b.T) / f
            k[i - 1, :, :] = (t2 @ t1).sum(axis=0) + L1

            f *= i + 1
            b = b @ a
        m[-1,:,:] = m[-2,:,:] @ a
        return k, m

    @staticmethod
    def _L(A, B, K, n, N_conv, l, alpha, sigma_bar):
        invA = inv(A)
        def M_fun(s):
            # FIXME: distinguis between singular and not singular matrices: can't use inv on singular matrix)
            return np.matmul(invA, (expm(A*s)-np.eye(n)))
        L = np.zeros((l, N_conv+1, n, n))
        Lsum, Afac = TrafficModelPETC._L_sum(A, n, N_conv)
        dsigma_prime = sigma_bar / l
        sigma_prime = 0
        In = np.eye(n)
        abk = (A - B @ K)
        for j in range(l):
            M = M_fun(sigma_prime)
            P1 = In + M @ abk
            N = A @ M + In
            P2 = N @ abk
            Li = np.zeros((N_conv+1, n, n))
            Li[0, :, :] = In - P1 - P1.T + (1 - alpha)*P1.T @ P1
            Li[1, :, :] = ((1 - alpha)*P1.T - In)@P2 + P2.T @ ((1 - alpha)*P1 - In)
            for i in range(N_conv-1):
                Li[i + 2, :, :] = ((1 - alpha) * P1.T - In) @ Afac[i, :,:] @ P2 \
                                        + P2.T @ Afac[i, :, :].T @ ((1 - alpha) * P1 - In) \
                                        + (1 - alpha) * P2.T @ Lsum[i, :, :] @ P2
            sigma_prime += dsigma_prime
            L[j, :, :, :] = Li
        return L
