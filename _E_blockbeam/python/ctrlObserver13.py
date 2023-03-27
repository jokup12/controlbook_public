#full state feedback
import numpy as np
import control as cnt
import blockbeamParam as P



class ctrlObserver:
    def __init__(self):
        #tuning params
        trTheta = 2.0
        trZ = 1.0
        zetaZ = .7
        zetaTheta = .7
        integrator_pole = -5.0
        wnTheta = np.pi / (2 * trTheta * np.sqrt(1 - zetaTheta**2))
        wnZ = np.pi / (2 * trZ * np.sqrt(1 - zetaZ ** 2))


        #observer poles
        tr_z_obs = trZ / 10.0  # rise time for position
        tr_theta_obs = trTheta / 10.0  # rise time for angle
        zeta_obs = .707
        dist_obs_pole = -5.0
        wn_z_obs = np.pi / (2 * tr_z_obs * np.sqrt(1 - zeta_obs ** 2))
        wn_theta_obs = np.pi / (2 * tr_theta_obs * np.sqrt(1 - zeta_obs ** 2))


        #State Space Matrices
        coefficientA = P.m2*P.length**2/3 + P.m1*P.ze**2
        self.A = np.array([[0.0, 0.0, 1.0, 0.0],
                           [0.0, 0.0, 0.0, 1.0],
                           [0.0, -P.g, 0.0, 0.0],
                           [-P.m1*P.g/coefficientA, 0.0, 0.0, 0.0]])
        self.B = np.array([[0.0],
                           [0.0],
                           [0.0],
                           [P.length/coefficientA]])
        self.C = np.array([[1.0, 0.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0, 0.0]])
        self.D = np.array([[0.0], [0.0]])

        self.Cr = self.C[0]
        #print('Cr: ', self.Cr)

        A1 = np.vstack((np.hstack((self.A, np.zeros((np.size(self.A, 1), 1)))),
                        np.hstack((-self.Cr, np.zeros((1))))))
        #print('A1:', A1)
        B1 = np.vstack( (self.B, np.zeros((1,1))) )
        #print('B1: ', B1)

        CC = cnt.ctrb(self.A, self.B)
        CC1 = cnt.ctrb(A1, B1)

        #gain calculations
        des_char_poly = np.convolve(
            np.convolve([1, 2 * zetaZ * wnZ, wnZ ** 2],
                        [1, 2 * zetaTheta * wnTheta, wnTheta ** 2]),
            [1, -integrator_pole])
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(CC1) != 5:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, B1, des_poles)
            self.K = K1[0][0:4]
            self.ki = K1[0][4]

        #print('K: ', self.K)
        #print('ki: ', self.ki)
        #print('poles: ', des_poles)

        #observer
        #A2 = np.concatenate((
        #        np.concatenate((self.A, self.B), axis=1),
        #        np.zeros((1,5))), axis=0)
        #C2 = np.concatenate((self.C, np.zeros((2,1))), axis=1)

        des_obs_char_poly = np.convolve(
                [1, 2 * zetaZ * wn_z_obs, wn_z_obs**2],
                [1, 2 * zetaTheta * wn_theta_obs, wn_theta_obs**2])
        des_obs_poles = np.roots(des_obs_char_poly)
        print('observer poles: ', des_obs_poles)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(self.A.T, self.C.T)) != 4:
            print("The system is not observerable")
            print('rank = ', np.linalg.matrix_rank(cnt.ctrb(self.A.T, self.C.T)))
        else:
            #print('A.T ', self.A.T)
            #print('C.T ', self.C.T)
            self.L = cnt.place(self.A.T, self.C.T, des_obs_poles).T

        #self.observer_state = np.array([[0],[0],[0],[0],[0]])
        #print('L^T: ', self.L.T)


        self.z = 0.0
        self.zdot = 0.0
        self.zd1 = 0.0
        self.errordotZ = 0.0
        self.errorZd1 = 0.0
        self.integratorZ = 0.0 # integrator
        self.theta = 0.0
        self.thetadot = 0.0
        self.thetad1 = 0.0
        self.errordotTheta = 0.0 # estimated derivative of error
        self.errord1Theta = 0.0 # Error delayed by one sample
        self.integratorTheta = 0.0 # integrator
        self.thetaR = 0.0
        self.Fd1 = 0.
        self.F = 0.
        self.x_hat = np.array([[0.], [0.], [0.], [0.]])




    def update(self, ref, x):

        x_hat = self.update_observer(x)
        z_hat = self.Cr @ x_hat

        errorZ = ref - z_hat
        self.integratorZ = self.integratorZ + (P.Ts/2)*(errorZ + self.errorZd1)
        self.errorZd1 = errorZ


        xe = np.array([[P.ze], [0.], [0.], [0.]])
        x_tilde = x_hat - xe

        F_tilde = -self.K @ x_tilde - self.ki*self.integratorZ
        F_unsat = F_tilde[0] + P.Fe
        self.F = self.saturate(F_unsat.item(0), P.Fmax)
        self.integratorAntiWindup(self.F, F_unsat)
        self.Fd1 = self.F

        return self.F, x_hat

    def update_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.x_hat, y_m)
        F2 = self.observer_f(self.x_hat + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f(self.x_hat + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f(self.x_hat + P.Ts * F3, y_m)
        self.x_hat += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        return self.x_hat

    def observer_f(self, x_hat, y_m):
        xe = np.array([[P.ze], [0.], [0.], [0.]])
        xhat_dot = self.A @ x_hat \
                   + self.B * (self.F-P.Fe) \
                   + self.L @ (y_m-self.C @ x_hat)
        return xhat_dot

    def integratorAntiWindup(self, F, F_unsat):
        if self.ki != 0:
            self.integratorZ = self.integratorZ + P.Ts/self.ki*(F-F_unsat)

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u