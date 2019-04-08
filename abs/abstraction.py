"""
A traffic-model abstraction for event-triggered control systems,
as done by C. Hop based on the work of A. Kolarijani as found in
"A Formal Traffic Characterization of LTI Event-triggered Control Systems, A.S. Kolarijani et al"
"""
from etctime import TabuadaPETC, LinearPlant, LinearController, LinearETC
import numpy as np
from scipy.linalg import expm, inv, block_diag
import picos as pic
import cvxpy as cvx


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

        # q = pow(m, n-1)
        # self.l = l # number of partitions
        # tau_opt = self._line_search(0, 0.001, sigma_bar)
        # nu = self._nu(self.plant.A, self.plant.B, self.controller.K, N_conv, l, self.trigger.sigma, sigma_bar, tau_opt)
        # tau_opt_nu = self._tau_opt_nu(tau_opt, l, sigma_bar, N_conv, n, nu)
        # print(tau_opt_nu)
        # self.lower_bounds(n, m, sigma_bar, tau_opt, nu, l, N_conv)

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
        sigma_prime = dsigma_prime**np.arange(N_conv+1)
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
    def _nu(A, B, K, N_Conv, l, alpha, sigma_bar, tau_opt, del_sig=0.001, L=None):
        i_max = int(np.floor((sigma_bar/l)/del_sig) + 1)
        (n, _) = A.shape
        if L is None:
            L = TrafficModelPETC._L(A, B, K, n, N_Conv, l, alpha, sigma_bar)
        invA = inv(A)
        In = np.eye(n)
        abk = A - B@K

        def int_eAs(s):
            # FIXME: distinguis between singular and not singular matrices: can't use inv on singular matrix)
            return invA @ (expm(A*s)-In)

        def _phi(s):
            Lambda = In + int_eAs(s)@abk
            return (In - Lambda.T) @ (In - Lambda) - alpha*Lambda.T@Lambda

        max_eig = np.zeros((i_max*l, 2))
        for i in range(i_max):
            subval = np.zeros((l, n, n))
            sig_prim = i*del_sig
            sigma_prime = (sig_prim)**np.arange(N_Conv+1)
            sigma_prime = sigma_prime.reshape(1,N_Conv+1,1,1)
            phi_hat = (sigma_prime*L[:, :, :, :]).sum(axis=1)
            for r in range(l):
                phi_c = _phi(sig_prim+r*sigma_bar/l)
                subval[r,:,:] = phi_c - phi_hat[r,:,:]
                eig = np.linalg.eigvalsh(subval[r,:,:])
                max_eig[i*l+r, :] = eig

        nu = np.amax(np.abs(max_eig))
        return nu

    def _tau_opt_nu(self, tau_opt, l, sigma_bar, N_Conv, n, nu, d_tau=0.001, L=None):
        if L is None:
            L = self._L(self.plant.A, self.plant.B, self.controller.K, n, N_Conv, l, self.trigger.sigma, sigma_bar)
        In = np.eye(n)
        j_max = int(np.floor(tau_opt/d_tau))
        for j in range(j_max):
            tau = tau_opt - j*d_tau
            phi, max_eig = TrafficModelPETC._phi(tau, l, sigma_bar, N_Conv, n, L)
            phi = phi + np.reshape(In*nu, (1, 1, n, n))
            eig, eiv = np.linalg.eig(phi)
            max_eig = np.amax(eig.reshape(np.prod(eig.shape)))
            if max_eig <= 0:
                return tau

    def lower_bounds(self, n, m, sigma_bar, tau_opt, nu, l, N_conv, d_tau=0.001, L=None, e_tol=0.1, sedumi_eps=10^(-5)):
        """
        Find the lower time bounds on the sampling time for each region and for an nD system,
        using projections of each region onto 2D planes and their corresponding 2x2 Q-matrices
        :return:
        """

        if L is None:
            L = self._L(self.plant.A, self.plant.B, self.controller.K, n, N_conv, l, self.trigger.sigma, sigma_bar)

        dimensions = [np.arange(1,m+1) for i in range(1, n-1)]
        dimensions.append(np.arange(1,2*m+1))
        # we can take int(m/2) since m is always an even number
        q = pow(m, n - 1)
        Q = self.isotropic_covering_2d(m)
        if n is 3:
            regions = np.stack(np.meshgrid(*dimensions)).transpose().reshape(q, n-1, 1)
        elif n is 2:
            regions = np.stack(*dimensions).reshape(2*m,1)

        for mm in range(q):
            print(f'{100*mm/(q+1)}%')
            Q_sec = regions[mm,:] #np.array(*regions[mm,:])
            f = int(np.floor((sigma_bar - tau_opt)/d_tau))
            for ff in range(f):
                tau_s = tau_opt + d_tau
                phi, _ = self._phi(tau_s, l, sigma_bar, N_conv, n, L)
                # if f is 1:
                for y in range(int(np.floor(tau_s*l/sigma_bar))):
                    for z in range(N_conv):
                        prob = pic.Problem()
                        e = prob.add_variable('e',n-1)
                        # create variables
                        Q_found = np.empty((n-1, n, n))
                        Q_LMI = np.empty((n - 1, n, n))  # n in {2,3} probably
                        for i in range(n-1):
                            index_Qfound = Q_sec[i]-1
                            Q_found[i] = Q[index_Qfound]
                        inEq = 0
                        for q in range(n-1):
                            Q2D = Q_found[q]
                            z1 = (n-2)*q
                            z2 = (n-2)*(n-1-q)
                            Q_nD = block_diag(*[np.zeros((z1, z1)), Q2D, np.zeros((z2, z2))])
                            # Q_LMI[i] = Q_nD
                            inEq += e[q]*Q_nD
                        con_e = [e >= e_tol, phi[y,z] + nu*np.eye(n) + inEq <= sedumi_eps]
                        prob.add_list_of_constraints(con_e)
                        prob.solve(solver='cvxopt', verbose=0)
        return regions

    def isotropic_covering_2d(self, m):
        Q = np.empty((m, 2, 2))
        thet_region = np.arange(1, m+1)
        thet_min = (thet_region-1)*np.pi/m
        thet_max = thet_region*np.pi/m
        for i in range(m):
            a1 = np.array([[-np.sin(thet_min[i])],[np.cos(thet_min[i])]])
            a2 = np.array([[np.sin(thet_max[i])],[-np.cos(thet_max[i])]])
            Q[i] = a1.T@a2 + a2.T@a1
        return Q


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
        m[-1,: ,:] = m[-2,:,:] @ a
        return k, m

    @staticmethod
    def _L(A, B, K, n, N_conv, l, alpha, sigma_bar):
        invA = inv(A)
        def M_fun(s):
            # FIXME: distinguis between singular and not singular matrices: can't use inv on singular matrix)
            return invA @ (expm(A*s)-np.eye(n))
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
