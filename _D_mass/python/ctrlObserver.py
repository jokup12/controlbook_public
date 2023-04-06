#full state feedback
import numpy as np
import control as cnt
import massParam as P



class ctrlObserver:
    def __init__(self):
        tr = 1
        zeta = .707
        integrator_pole = .35




        #State Space Matrices
        self.A = np.array([[0, 1], [-P.k/P.m, -P.b/P.m]])
        self.B = np.array([[0], [1/P.m]])
        self.C = np.array([[1, 0]])
        D = np.array([0])

        A1 = np.vstack((np.hstack((self.A, np.zeros((np.size(self.A,1),1)))),
                        np.hstack((-self.C, np.array([[0]]))) ))
        B1 = np.vstack( (self.B, 0) )

        #gain calculations
        wn = np.pi / (2 * tr * np.sqrt(1 - zeta**2))
        des_char_poly = np.convolve([1, 2 * zeta * wn, wn ** 2],
                                    [1, -integrator_pole])
        des_poles = np.roots(des_char_poly)

        #CC = cnt.ctrb(A, B)
        CC1 = cnt.ctrb(A1, B1)

        # check for controllability
        if np.linalg.matrix_rank(CC1) != np.size(A1, 1):
            print("The system is not controllable")
        else:
            self.K1 = (cnt.place(A1, B1, des_poles))
            self.K = self.K1[0][0:2]
            self.ki = self.K1[0][2]
            #self.kr = -1.0 / (C @ np.linalg.inv(A - B @ self.K) @ B)

        #observer
        tr_obs = tr/10
        zeta_obs = .707
        wn_obs = np.pi / (2 * tr_obs * np.sqrt(1 - zeta_obs**2))
        des_obsv_char_poly = [1, 2 * zeta_obs * wn_obs, wn_obs ** 2]
        des_obsv_poles = np.roots(des_obsv_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(self.A.T, self.C.T)) != 2:
            print("The system is not observerable")
        else:
            self.L = cnt.acker(self.A.T, self.C.T, des_obsv_poles).T
        print('K: ', self.K)
        print('ki ', self.ki)
        print('L^T: ', self.L.T)

        self.errordot = 0.0 # estimated derivative of error
        self.errord1 = 0.0 # Error delayed by one sample
        self.integrator = 0.0 # integrator
        self.z = 0.0
        self.zdot = 0.0
        self.zd1 = 0.0
        # variables to implement integrator
        self.integrator = 0.0  # integrator
        self.error_d1 = 0.0  # error signal delayed by 1 sample
        self.x_hat = np.array([
            [0.0],  # theta_hat_0
            [0.0],  # thetadot_hat_0
        ])
        self.F_d1 = 0.0  # control force, delayed 1 sample
        self.F = 0.0


    def update(self, ref, state):

        z_hat = self.update_observer(state)
        #self.z = state[0][0]
        #self.zdot = x[1][0]


        error = ref - z_hat
        self.integrator = self.integrator + (P.Ts/2) * (error + self.errord1)

        #self.zdot = ((2 * P.sigma - P.Ts) / (2 * P.sigma + P.Ts)) * self.zdot \
        #            + (2 / (2 * P.sigma + P.Ts)) * (self.z - self.zd1)

        #F = self.kp * (self.ref - self.z) - self.kd * self.zdot
        self.F = -self.K @ z_hat + self.ki*self.integrator
        self.F = self.F.item(0)
        self.F = self.saturate(self.F, P.F_max)

        self.zd1 = self.z
        self.errord1 = error
        return self.F, z_hat

    def update_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.x_hat, y_m)
        F2 = self.observer_f(self.x_hat + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f(self.x_hat + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f(self.x_hat + P.Ts * F3, y_m)
        self.x_hat += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        return self.x_hat

    def observer_f(self, x_hat, y_m):
        # compute feedback linearizing force
        z_hat = x_hat[0][0]
        Fe = z_hat * P.k

        # xhatdot = A*(xhat-xe) + B*(u-ue) + L(y-C*xhat)
        zhat_dot = self.A @ x_hat\
                   + self.B * (self.F)\
                   + self.L * (y_m - self.C @ x_hat)
        return zhat_dot


    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u