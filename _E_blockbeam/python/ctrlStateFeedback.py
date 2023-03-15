#full state feedback
import numpy as np
import control as cnt
import blockbeamParam as P



class ctrlStateFeedback:
    def __init__(self):
        trTheta = .2
        trZ = 10.0*trTheta
        zeta = .707

        wnTheta = np.pi / (2 * trTheta * np.sqrt(1 - zeta**2))
        des_char_poly_theta = (1, 2 * zeta * wnTheta, wnTheta ** 2)
        des_poles_theta = np.roots(des_char_poly_theta)

        wnZ = np.pi / (2 * trZ * np.sqrt(1 - zeta ** 2))
        des_char_poly_Z = (1, 2 * zeta * wnZ, wnZ ** 2)
        des_poles_Z = np.roots(des_char_poly_Z)

        des_poles = np.array([des_poles_theta[0], des_poles_theta[1], des_poles_Z[0], des_poles_Z[1]])

        #State Space Matrices
        coefficientA = P.m2*P.length**2/3 + P.m1*P.ze**2
        A = np.array([[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0], [0.0, -P.g, 0.0, 0.0], [-P.m1*P.g/coefficientA, 0.0, 0.0, 0.0]])
        B = np.array([[0.0], [0.0], [0.0], [P.length/coefficientA]])
        C = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
        D = np.array([[0.0], [0.0]])

        CC = cnt.ctrb(A, B)
        print(np.linalg.matrix_rank(CC))

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 4:
            print("The system is not controllable")
        else:
            self.K = (cnt.place(A, B, des_poles))
            Cr = C[0]
            self.kr = -1.0 / (Cr @ np.linalg.inv(A - B @ self.K) @ B)
        print('K: ', self.K)
        print('kr: ', self.kr)
        print('poles: ', des_poles)


        self.z = 0.0
        self.zdot = 0.0
        self.zd1 = 0.0
        self.theta = 0.0
        self.thetadot = 0.0
        self.thetad1 = 0.0



    def update(self, ref, x):
        self.z = x[0][0]
        self.theta = x[1][0]
        #self.zdot = x[2][0]
        #self.thetadot = x[3][0]
        x_tilde = x - np.array([[P.z0], [0], [0], [0]])
        zr_tilde = ref - P.z0
        self.zdot = ((2 * P.sigma - P.Ts) / (2 * P.sigma + P.Ts)) * self.zdot \
                    + (2 / (2 * P.sigma + P.Ts)) * (self.z - self.zd1)
        self.thetadot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.thetadot \
                        + (2/(2*P.sigma + P.Ts))*(self.theta - self.thetad1)

        #F = self.kp * (self.ref - self.z) - self.kd * self.zdot
        #F = -self.K @ x + self.kr*ref
        F_tilde = -self.K @ x_tilde + self.kr*zr_tilde
        F = F_tilde[0] + P.Fe
        F = F.item(0)
        F = self.saturate(F, P.Fmax)


        self.zd1 = self.z
        self.thetad1 = self.theta
        return F

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u