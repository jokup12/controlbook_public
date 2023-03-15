import numpy as np
import blockbeamParam as P


class ctrlPD:
    def __init__(self):

        A = (P.m2*P.length**2/3) + (P.m1*P.length**2/4)
        # inner loop
        # tr = 0.8 # part (a)
        trTheta = .24 #rise time
        zeta = 0.707
        # desired natural frequency
        wnTheta = np.pi/(2*trTheta*np.sqrt(1-zeta**2))

        #inner loop
        b0Theta = P.length/A
        a1Theta = 0.0
        a0Theta = 0.0
        alpha1Theta = 2.0 * zeta * wnTheta
        alpha0Theta = wnTheta**2
        # compute PD gains
        self.kdTheta = (alpha1Theta-a1Theta) / b0Theta
        self.kpTheta = (alpha0Theta - a0Theta) / b0Theta

        # outer loop
        trZ = 10.0*trTheta
        zeta = .707
        wnZ = np.pi/(2*trZ*np.sqrt(1-zeta**2))
        b0Z = -P.g
        a0Z = 0.0
        a1Z = 0.0
        alpha1Z = 2.0 * zeta * wnZ
        alpha0Z = wnZ ** 2
        # compute PD gains for outer loop
        self.kdZ = (alpha1Z-a1Z) / b0Z
        self.kpZ = (alpha0Z - a0Z) / b0Z

        print('kpTheta: ', self.kpTheta)
        print('kdTheta: ', self.kdTheta)
        print('kpZ: ', self.kpZ)
        print('kdZ: ', self.kdZ)

    def update(self, r ,state):
        z = state[0][0]
        theta = state[1][0]
        zdot = state[2][0]
        thetadot = state[3][0]


        thetaR = self.kpZ*(r-z) - self.kdZ*zdot   #theta reference
        FPD = self.kpTheta*(thetaR-theta) - self.kdTheta*thetadot   #z reference
        Fe = (P.m1*P.g*z)/P.length + P.m2*P.g/2
        F = FPD + Fe

        F = self.saturate(F, P.Fmax)
        return F


    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
            print('limit reached')
        return u








