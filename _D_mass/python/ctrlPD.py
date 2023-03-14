import numpy as np
import massParam as P


class ctrlPD:
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
        print('kp: ', self.kp)
        print('kd: ', self.kd)        

    def update(self, r ,state):
        z = state[0][0]
        zdot = state[1][0]

        F = self.kp*(r-z)-self.kd*zdot
        F = saturate(F, P.F_max)
        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
        print('limit reached')
    return u








