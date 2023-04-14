import numpy as np
import VTOLParam as P


class ctrlPID:
    def __init__(self):
        self.M = P.mc + 2 * P.mr

        # height loop
        trH = 5  # rise time
        zeta = 0.707
        # desired natural frequency
        wnH = np.pi / (2 * trH * np.sqrt(1 - zeta ** 2))
        b0H = 1 / self.M
        a1H = 0.0
        a0H = 0.0
        alpha1H = 2.0 * zeta * wnH
        alpha0H = wnH ** 2
        # compute PD gains
        self.kdH = (alpha1H - a1H) / b0H
        self.kpH = (alpha0H - a0H) / b0H
        self.kiH = .2
        self.h = 0.0
        self.errordotH = 0.0
        self.errord1H = 0.0
        self.integratorH = 0.0
        self.hd1 = 0.0
        self.hdot = 0.0

        # inner theta loop
        # tr = 0.8 # part (a)
        trTheta = .5  # rise time
        zeta = 0.707
        # desired natural frequency
        wnTheta = np.pi / (2 * trTheta * np.sqrt(1 - zeta ** 2))
        b0Theta = 1 / (P.Jc + 2 * P.mr * P.d)
        a1Theta = 0.0
        a0Theta = 0.0
        alpha1Theta = 2.0 * zeta * wnTheta
        alpha0Theta = wnTheta ** 2
        # compute PD gains
        self.kdTheta = (alpha1Theta - a1Theta) / b0Theta
        self.kpTheta = (alpha0Theta - a0Theta) / b0Theta
        self.theta = 0.0
        self.errordotTheta = 0.0
        self.errord1Theta = 0.0
        self.thetaref = 0.0
        self.thetad1 = 0.0
        self.thetadot = 0.0

        # outer z loop
        trZ = 10.0 * trTheta
        zeta = .707
        wnZ = np.pi / (2 * trZ * np.sqrt(1 - zeta ** 2))
        b0Z = -P.g
        a0Z = 0.0
        a1Z = P.mu / self.M
        alpha1Z = 2.0 * zeta * wnZ
        alpha0Z = wnZ ** 2
        # compute PD gains for outer loop
        self.kdZ = (alpha1Z - a1Z) / b0Z
        self.kpZ = (alpha0Z - a0Z) / b0Z
        self.kiZ = -0.001
        self.z = 0.0
        self.errordotZ = 0.0
        self.errord1Z = 0.0
        self.integratorZ = 0.0
        self.zd1 = 0.0
        self.zdot = 0.0

        ################
        #temporary values for testing
        self.kpH = 4.218
        self.kdH = 4.779
        self.kiH = 1.0
        self.kpTheta = 7.099
        self.kdTheta = 1.064
        self.kiTheta = 0
        self.kpZ = -.147
        self.kdZ = -.213
        self.kiZ = 0


        print('Longitude gains    kpH: ', self.kpH, '  kdH: ', self.kdH, '  kiH: ', self.kiH)
        print('Inner Loop gains   kpTheta: ', self.kpTheta, '  kdTheta: ', self.kdTheta)
        print('Outer Loop gains   kpZ: ', self.kpZ, '  kdZ: ', self.kdZ, '   kiZ: ', self.kiZ)


    def update(self, zref, href, state):
        #print(state.shape)

        self.z = state[0][0]
        self.h = state[1][0]
        self.theta = state[2][0]
        errorZ = zref - self.z
        errorH = href - self.h
        errorTheta = self.thetaref - self.theta
        self.integratorH = self.integratorH + (P.Ts/2)*(errorH + self.errord1H)
        self.integratorZ = self.integratorZ + (P.Ts/2)*(errorZ + self.errord1Z)
        if(errorH < .1):
            self.integratorH = 0
        if(errorZ < .1):
            self.integratorZ = 0


        self.zdot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.zdot \
                    + (2/(2*P.sigma + P.Ts))*(self.z - self.zd1)
        self.hdot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.hdot \
                    + (2/(2*P.sigma + P.Ts))*(self.h - self.hd1)
        self.thetadot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.thetadot \
                    + (2/(2*P.sigma + P.Ts))*(self.theta - self.thetad1)

        FPD = self.kpH * (href - self.h) - self.kdH * self.hdot + self.kiH*self.integratorH   #height reference
        Fe = self.M * P.g
        F = FPD + Fe

        self.thetaref = self.kpZ * (zref - self.z) - self.kdZ * self.zdot + self.kiZ*self.integratorZ #outer loop, thetaR
        TPD = self.kpTheta * (self.thetaref - self.theta) - self.kdTheta * self.thetadot  #inner loop, Z
        Te = 0
        T = TPD + Te

        u = np.array([[F], [T]])
        v = P.mixing @ np.array([[F], [T]])  # changes to fr and fl
        self.fr = v[0][0]
        self.fl = v[1][0]
        #self.fr = self.saturate(self.fr, P.Fmax)
        #self.fl = self.saturate(self.fl, P.Fmax)
        v = np.array([[self.fr], [self.fl]]) #[fr, fl]T
        u = np.array([[self.fr + self.fl], [P.d * (self.fr - self.fl)]]) #[F, torque]T


        self.errord1H = errorH
        self.errord1Z = errorZ
        self.errord1Theta = errorTheta
        self.hd1 = self.h
        self.zd1 = self.z
        self.thetad1 = self.theta

        return u, v

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
            print('limit reached')
        return u








