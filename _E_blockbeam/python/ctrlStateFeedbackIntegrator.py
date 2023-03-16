#full state feedback
import numpy as np
import control as cnt
import blockbeamParam as P



class ctrlStateFeedback:
    def __init__(self):
        trTheta = .1
        trZ = 10.0*trTheta
        zetaZ = .4
        zetaTheta = .4
        integrator_pole = .28

        wnTheta = np.pi / (2 * trTheta * np.sqrt(1 - zetaTheta**2))
        des_char_poly_theta = (1, 2 * zetaTheta * wnTheta, wnTheta ** 2)
        des_poles_theta = np.roots(des_char_poly_theta)

        wnZ = np.pi / (2 * trZ * np.sqrt(1 - zetaZ ** 2))
        des_char_poly_Z = (1, 2 * zetaZ * wnZ, wnZ ** 2)
        des_poles_Z = np.roots(des_char_poly_Z)

        des_char_poly = np.convolve(
            np.convolve([1, 2 * zetaZ * wnZ, wnZ ** 2],
                        [1, 2 * zetaTheta * wnTheta, wnTheta ** 2]),
                        np.poly([integrator_pole]))
        des_poles = np.roots(des_char_poly)

        #State Space Matrices
        coefficientA = P.m2*P.length**2/3 + P.m1*P.ze**2
        A = np.array([[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0], [0.0, -P.g, 0.0, 0.0], [-P.m1*P.g/coefficientA, 0.0, 0.0, 0.0]])
        B = np.array([[0.0], [0.0], [0.0], [P.length/coefficientA]])
        C = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
        D = np.array([[0.0], [0.0]])

        Cr = C[0]
        print('Cr: ', Cr)
        B1 = np.vstack( (B, np.zeros((1,1))) )
        print('B1: ', B1)

        A1 = np.vstack((np.hstack((A, np.zeros((np.size(A, 1), 1)))),
                        np.hstack((-Cr, np.zeros((1))))))
        print('A1:', A1)

        CC = cnt.ctrb(A, B)
        CC1 = cnt.ctrb(A1, B1)
        #print(np.linalg.matrix_rank(CC))

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(CC1) != 5:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, B1, des_poles)
            self.K = K1[0][0:4]
            self.ki = K1[0][4]

        print('K: ', self.K)
        print('ki: ', self.ki)
        print('poles: ', des_poles)


        self.z = 0.0
        self.zdot = 0.0
        self.zd1 = 0.0
        self.errordotZ = 0.0
        self.errord1Z = 0.0
        self.integratorZ = 0.0 # integrator
        self.theta = 0.0
        self.thetadot = 0.0
        self.thetad1 = 0.0
        self.errordotTheta = 0.0 # estimated derivative of error
        self.errord1Theta = 0.0 # Error delayed by one sample
        self.integratorTheta = 0.0 # integrator
        self.thetaR = 0.0



    def update(self, ref, x):
        self.z = x[0][0]
        self.theta = x[1][0]
        errorZ = ref - self.z
        errorTheta = self.thetaR - self.theta
        #self.zdot = x[2][0]
        #self.thetadot = x[3][0]
        self.integratorTheta = self.integratorTheta + (P.Ts/2)*(errorTheta + self.errord1Theta)
        self.integratorZ = self.integratorZ + (P.Ts/2)*(errorZ + self.errord1Z)


        self.zdot = ((2 * P.sigma - P.Ts) / (2 * P.sigma + P.Ts)) * self.zdot \
                    + (2 / (2 * P.sigma + P.Ts)) * (self.z - self.zd1)
        self.thetadot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.thetadot \
                        + (2/(2*P.sigma + P.Ts))*(self.theta - self.thetad1)

        F_tilde = -self.K @ x + self.ki*self.integratorZ
        F = F_tilde[0] + P.Fe
        F = F.item(0)
        F = self.saturate(F, P.Fmax)


        self.errord1Z = errorZ
        self.errord1Theta = errorTheta
        self.zd1 = self.z
        self.thetad1 = self.theta

        return F

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u