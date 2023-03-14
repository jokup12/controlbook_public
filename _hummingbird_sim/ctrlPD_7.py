import numpy as np
import hummingbirdParam as P


class ctrlPD:
    def __init__(self):

        self.M = P.m1 + P.m2 + P.m3

        #height loop
        trTheta = 1.5  # rise time
        zeta = 0.707
        # desired natural frequency
        wnTheta = np.pi / (2 * trTheta * np.sqrt(1 - zeta ** 2))
        b0Theta = P.ellT / (P.m1*P.ell1**2 + P.m2*P.ell2**2 + P.J1y + P.J2y)
        a1Theta = 0.0
        a0Theta = 0.0
        alpha1Theta = 2.0 * zeta * wnTheta
        alpha0Theta = wnTheta ** 2
        # compute PD gains
        self.kdTheta = (alpha1Theta - a1Theta) / b0Theta
        self.kpTheta = (alpha0Theta - a0Theta) / b0Theta
        self.thetadot = 0.0
        self.theta1 = 0.0


        # inner theta loop  IGNORE, CHANGE LATER
        # tr = 0.8 # part (a)
        #trTheta = .4 #rise time
        #zeta = 0.707
        # desired natural frequency
        #wnTheta = np.pi/(2*trTheta*np.sqrt(1-zeta**2))
        #b0Theta = 1#1/(P.Jc +2*P.mr*P.d)
        a1Theta = 0.0
        a0Theta = 0.0
        alpha1Theta = 2.0 * zeta * wnTheta
        alpha0Theta = wnTheta**2
        # compute PD gains
        #self.kdTheta = (alpha1Theta - a1Theta) / b0Theta
        #self.kpTheta = (alpha0Theta - a0Theta) / b0Theta

        # outer z loop  IGNORE, CHANGE LATER
        trZ = 10.0*trTheta
        zeta = .707
        wnZ = np.pi/(2*trZ*np.sqrt(1-zeta**2))
        b0Z = -P.g
        a0Z = 0.0
        a1Z = 0
        alpha1Z = 2.0 * zeta * wnZ
        alpha0Z = wnZ ** 2
        # compute PD gains for outer loop
        self.kdZ = (alpha1Z-a1Z) / b0Z
        self.kpZ = (alpha0Z - a0Z) / b0Z

        print('B0H: ', b0Theta)
        print('kpTheta: ', self.kpTheta)
        print('kdTheta: ', self.kdTheta)
        #print('kpTheta: ', self.kpTheta)
        #print('kdTheta: ', self.kdTheta)
        #print('kpZ: ', self.kpZ)
        #print('kdZ: ', self.kdZ)

    def update(self, ref ,state):
        phi = state[0][0]
        theta = state[1][0]
        psi = state[2][0]
        #phidot = state[3][0]
        #thetadot = state[4][0]
        #psidot = state[5][0]

        #dirty derivative and error for theta
        self.thetadot = ((2*P.sigma - P.Ts)/(2*P.sigma + P.Ts)) * self.thetadot \
                + (2/(2*P.sigma + P.Ts))*(theta-self.theta1)
        self.theta1 = theta
        #if (abs(ref[1][0] - theta) < .1)


        #height reference PD control
        FPD = self.kpTheta*(ref[1][0]-theta) - self.kdTheta*self.thetadot
        Fe = np.cos(theta)*(P.m1*P.ell1 + P.m2*P.ell2)*P.g/P.ellT
        F = FPD + Fe


        T = 0.0
        

        u = np.array([[F], [T]])
        v = np.linalg.inv(np.array([[P.d, -P.d], [1.0, 1.0]])) @ np.array([[T], [F]]) #changes to fl and fr
        v = v/P.km

        #saturate fl and fr
        fl = v[0][0]
        fr = v[1][0]
        v[0][0] = self.saturate(fl, P.PWM_min, P.PWM_max)
        v[1][0] = self.saturate(fr, P.PWM_min, P.PWM_max)
        u = np.array([[fl+fr], [T]])   #[P.d*(fl-fr)])
        v[0][0] = fl
        v[1][0] = fr

        return u, v


    def saturate(self, u, minlimit, maxlimit):
        if abs(u) < minlimit:
            u = minlimit
            #print('min reached')
        if abs(u) > maxlimit:
            u = maxlimit
            print('max reached')
        return u








