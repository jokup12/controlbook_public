import numpy as np
import blockbeamParam as P


class ctrlPID:
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
        #self.kiTheta = -.05

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
        self.kiZ = -.02

        print('kpTheta: ', self.kpTheta)
        print('kdTheta: ', self.kdTheta)
        print('kpZ: ', self.kpZ)
        print('kdZ: ', self.kdZ)

        self.errordotTheta = 0.0 # estimated derivative of error
        self.errord1Theta = 0.0 # Error delayed by one sample
        self.errordotZ = 0
        self.errord1Z = 0
        self.integratorTheta = 0.0 # integrator
        self.integratorZ = 0.0 # integrator
        self.theta = 0
        self.thetadot = 0
        self.thetad1 = 0
        self.z = 0.0      #estimated position
        self.zdot = 0.0   #estimated derivative
        self.zd1 = 0.0    #position delayed 1 sample
        self.thetaR = 0.0

    def update(self, r ,state):
        self.z = state[0][0]
        self.theta = state[1][0]
        errorZ = r - self.z
        errorTheta = self.thetaR - self.theta
        self.integratorTheta = self.integratorTheta + (P.Ts/2)*(errorTheta + self.errord1Theta)

        self.zdot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.zdot \
                    + (2/(2*P.sigma + P.Ts))*(self.z - self.zd1)
        self.thetadot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.thetadot \
                        + (2/(2*P.sigma + P.Ts))*(self.theta - self.thetad1)


        self.thetaR = self.kpZ*(r-self.z) - self.kdZ*self.zdot #+ self.kiZ*self.integratorZ  #outer loop
        if(np.abs(self.zdot) < .05):
            self.integratorZ = self.integratorZ + (P.Ts / 2) * (errorZ + self.errord1Z)
            self.thetaR = self.thetaR + self.kiZ*self.integratorZ
            print('integratorZ = ', self.integratorZ)
        else:
            print('zdot = ', self.zdot)
            self.integratorZ = 0        #resets the integrator

        FPD = self.kpTheta*(self.thetaR-self.theta) - self.kdTheta*self.thetadot  #inner loop
        Fe = (P.m1*P.g*self.z)/P.length + P.m2*P.g/2
        F = FPD + Fe

        self.errord1Z = errorZ
        self.errord1Theta = errorTheta
        self.zd1 = self.z
        self.thetad1 = self.theta
        #F = self.saturate(F, P.Fmax)
        return F


    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
            print('limit reached')
        return u








