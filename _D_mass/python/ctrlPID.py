import numpy as np
import massParam as P


class ctrlPID:
    def __init__(self):
        #  tuning parameters
        #tr = 0.8 # part (a)
        tr = 1.9 #rise time
        zeta = 0.7
        # desired natural frequency
        wn = np.pi/(2*tr*np.sqrt(1-zeta**2))
        alpha1 = 2.0 * zeta * wn
        alpha0 = wn**2
        # compute PD gains
        self.kd = P.m*(alpha1-P.b/P.m)
        self.kp = P.m*(alpha0-P.k/P.m)
        self.ki = 1.6

        self.errordot = 0.0 # estimated derivative of error
        self.errord1 = 0.0 # Error delayed by one sample
        self.integrator = 0.0 # integrator
        self.z = 0.0      #estimated position
        self.zdot = 0.0   #estimated derivative
        self.zd1 = 0.0    #position delayed 1 sample


        print('kp: ', self.kp)
        print('kd: ', self.kd)

    def update(self, r, state):
        self.z = state[0][0]
        error = r - self.z
        self.integrator = self.integrator + (P.Ts/2) * (error + self.errord1)

        self.zdot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts))*self.zdot \
                    + (2/(2*P.sigma + P.Ts))*(self.z - self.zd1)

        #solve for force with kp, kd, and ki
        F = self.kp*(r-self.z) - self.kd*self.zdot + self.ki*self.integrator
        F = saturate(F, P.F_max)
        self.errord1 = error
        self.zd1 = self.z
        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
        print('limit reached')
    return u








