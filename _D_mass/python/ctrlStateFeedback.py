#full state feedback
import numpy as np
import control as cnt
import massParam as P



class ctrlStateFeedback:
    def __init__(self):
        tr = 2.0
        zeta = .7

        #State Space Matrices
        A = np.array([[0, 1], [-P.k/P.m, -P.b/P.m]])
        B = np.array([[0], [1/P.m]])
        C = np.array([[1, 0]])
        D = np.array([0])

        #gain calculations
        wn = np.pi / (2 * tr * np.sqrt(1 - zeta**2))
        des_char_poly = (1, 2 * zeta * wn, wn ** 2)
        des_poles = np.roots(des_char_poly)

        CC = cnt.ctrb(A, B)

        # check for controllability
        if np.linalg.matrix_rank(CC) == np.size(A, 1):
            print('the system is controllable')
        else:
            print('the system is not controllable')

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 2:
            print("The system is not controllable")
        else:
            self.K = (cnt.place(A, B, des_poles))
            self.kr = -1.0 / (C @ np.linalg.inv(A - B @ self.K) @ B)
        print('K: ', self.K)
        print('kr: ', self.kr)
        print('poles: ', des_poles)


        self.z = 0.0
        self.zdot = 0.0
        self.zd1 = 0.0




    def update(self, ref, x):
        self.z = x[0][0]
        #self.zdot = x[1][0]
        self.zdot = ((2 * P.sigma - P.Ts) / (2 * P.sigma + P.Ts)) * self.zdot \
                    #+ (2 / (2 * P.sigma + P.Ts)) * (self.z - self.zd1)

        #F = self.kp * (self.ref - self.z) - self.kd * self.zdot
        F = -self.K @ x + self.kr*ref
        F = F.item(0)
        F = self.saturate(F, P.F_max)


        self.zd1 = self.z
        return F

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u