import numpy as np
import VTOLParam as P


class ctrlPD:
    def __init__(self):

        self.M = P.mc + 2*P.mr

        #height loop
        trH = 4  # rise time
        zeta = 0.707
        # desired natural frequency
        wnH = np.pi / (2 * trH * np.sqrt(1 - zeta ** 2))
        b0H = 1/self.M
        a1H = 0.0
        a0H = 0.0
        alpha1H = 2.0 * zeta * wnH
        alpha0H = wnH ** 2
        # compute PD gains
        self.kdH = (alpha1H - a1H) / b0H
        self.kpH = (alpha0H - a0H) / b0H



        # inner theta loop
        # tr = 0.8 # part (a)
        trTheta = .4 #rise time
        zeta = 0.707
        # desired natural frequency
        wnTheta = np.pi/(2*trTheta*np.sqrt(1-zeta**2))
        b0Theta = 1/(P.Jc +2*P.mr*P.d)
        a1Theta = 0.0
        a0Theta = 0.0
        alpha1Theta = 2.0 * zeta * wnTheta
        alpha0Theta = wnTheta**2
        # compute PD gains
        self.kdTheta = (alpha1Theta - a1Theta) / b0Theta
        self.kpTheta = (alpha0Theta - a0Theta) / b0Theta



        # outer z loop
        trZ = 10.0*trTheta
        zeta = .707
        wnZ = np.pi/(2*trZ*np.sqrt(1-zeta**2))
        b0Z = -P.g
        a0Z = 0.0
        a1Z = P.mu/self.M
        alpha1Z = 2.0 * zeta * wnZ
        alpha0Z = wnZ ** 2
        # compute PD gains for outer loop
        self.kdZ = (alpha1Z-a1Z) / b0Z
        self.kpZ = (alpha0Z - a0Z) / b0Z

        print('kpH: ', self.kpH)
        print('kdH: ', self.kdH)
        print('kpTheta: ', self.kpTheta)
        print('kdTheta: ', self.kdTheta)
        print('kpZ: ', self.kpZ)
        print('kdZ: ', self.kdZ)

    def update(self, r ,state):
        z = state[0][0]
        h = state[1][0]
        theta = state[2][0]
        zdot = state[3][0]
        hdot = state[4][0]
        thetadot = state[5][0]


        FPD = self.kpH*(r-h) - self.kdH*hdot     #height reference
        Fe = self.M*P.g
        F = FPD + Fe

        thetaR = self.kpZ*(r-z) - self.kdZ*zdot   #theta reference
        TPD = self.kpTheta*(thetaR-theta) - self.kdTheta*thetadot   #z reference
        Te = 0
        T = TPD + Te



        #print(u.shape)


        u = np.array([[F], [T]])
        v = P.mixing @ np.array([[F], [T]]) #changes to fr and fl
        self.fr = v[0][0]
        self.fl = v[1][0]
        self.fr = self.saturate(self.fr, P.Fmax)
        self.fl = self.saturate(self.fl, P.Fmax)
        v = np.array([[self.fr], [self.fl]])
        u = np.array([[self.fr+self.fl], [P.d*(self.fr-self.fl)]])

        return u, v


    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
            print('limit reached')
        return u








