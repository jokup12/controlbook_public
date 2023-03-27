#full state feedback
import numpy as np
import control as cnt
import massParam as P



class ctrlStateFeedback:
    def __init__(self):
        tr = 1
        zeta = .707
        integrator_pole = -.35


        #State Space Matrices
        A = np.array([[0, 1], [-P.k/P.m, -P.b/P.m]])
        B = np.array([[0], [1/P.m]])
        C = np.array([[1, 0]])
        D = np.array([0])

        A1 = np.vstack((np.hstack((A, np.zeros((np.size(A,1),1)))),
                        np.hstack((-C, np.array([[0]]))) ))
        B1 = np.vstack( (B, 0) )

        #gain calculations
        wn = np.pi / (2 * tr * np.sqrt(1 - zeta**2))
        des_char_poly = np.convolve([1, 2 * zeta * wn, wn ** 2],
                                    [1, -integrator_pole])
        des_poles = np.roots(des_char_poly)

        CC = cnt.ctrb(A, B)
        CC1 = cnt.ctrb(A1, B1)

        # check for controllability
        if np.linalg.matrix_rank(CC1) != np.size(A1, 1):
            print("The system is not controllable")
        else:
            self.K1 = (cnt.place(A1, B1, des_poles))
            self.K = self.K1[0][0:2]
            self.ki = self.K1[0][2]
            #self.kr = -1.0 / (C @ np.linalg.inv(A - B @ self.K) @ B)
        print('K1: ', self. K1)
        print('K: ', self.K)
        print('ki: ', self.ki)
        print('poles: ', des_poles)


        #self.ki = 1.6
        self.errordot = 0.0 # estimated derivative of error
        self.errord1 = 0.0 # Error delayed by one sample
        self.integrator = 0.0 # integrator
        self.z = 0.0
        self.zdot = 0.0
        self.zd1 = 0.0


    def update(self, ref, x):
        self.z = x[0][0]
        #self.zdot = x[1][0]
        error = ref - self.z
        self.integrator = self.integrator + (P.Ts/2) * (error + self.errord1)


        self.zdot = ((2 * P.sigma - P.Ts) / (2 * P.sigma + P.Ts)) * self.zdot \
                    + (2 / (2 * P.sigma + P.Ts)) * (self.z - self.zd1)

        #F = self.kp * (self.ref - self.z) - self.kd * self.zdot
        F = -self.K @ x - self.ki*self.integrator
        F = F.item(0)
        F = self.saturate(F, P.F_max)

        self.zd1 = self.z
        self.errord1 = error
        return F

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u