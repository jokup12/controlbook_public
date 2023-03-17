import numpy as np
import VTOLParam as P


class VTOLDynamics:
    def __init__(self, alpha):
        # Initial state conditions
        self.state = np.array([
            [P.z0],  # z initial position
            [P.h0],  # h initial height
            [P.theta0],  # Theta initial orientation
            [P.zdot0],  # zdot initial velocity
            [P.hdot0],  # hdot initial velocity
            [P.thetadot0],  # Thetadot initial velocity
        ])
        # Mass of the core, kg
        self.mc = P.mc * (1. + alpha * (2. * np.random.rand() - 1.))
        # mass of the right and left rotors, kg
        self.mr = P.mr * (1. + alpha * (2. * np.random.rand() - 1.))
        # self.ml = P.ml * (1. + alpha * (2. * np.random.rand() - 1.))
        # length between core and rotors, m
        self.d = P.d * (1. + alpha * (2. * np.random.rand() - 1.))
        # moment of inertia, kg m2
        self.Jc = P.Jc * (1. + alpha * (2. * np.random.rand() - 1.))
        # drag mu, kg/s
        self.mu = P.mu * (1. + alpha * (2. * np.random.rand() - 1.))
        # the gravity constant is well known, so we don't change it.
        self.g = P.g
        # sample rate at which the dynamics are propagated
        self.Ts = P.Ts

        print('Mass =', self.mc + 2*self.mr)


    def update(self, u):
        # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        # saturate the input torque

        self.rk4_step(u)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output
        return y

    def f(self, state, u):
        # Return xdot = f(x,u), the system state update equations
        # re-label states for readability
        fr = u[0][0]
        fl = u[1][0]
        z = state[0][0]
        h = state[1][0]
        theta = state[2][0]
        zdot = state[3][0]
        hdot = state[4][0]
        thetadot = state[5][0]

        #F = fr + fl

        zddot = (-(fr+fl)*np.sin(theta)-self.mu*zdot + .1)/(self.mc+2*self.mr)  #the .1 is wind force
        hddot = (-(self.mc+2*self.mr)*self.g + (fr+fl)*np.cos(theta))/(self.mc+2*self.mr)
        thetaddot = (self.d*(fr - fl))/(self.Jc + 2*self.mr*self.d*self.d)

        xdot = np.array([[zdot], [hdot], [thetadot], [zddot], [hddot], [thetaddot]])
        return xdot

    def h(self):
        # return the output equations
        # could also use input u if needed
        theta = self.state[0][0]
        y = np.array([[theta]])
        return y

    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

    def saturate(u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
            print('limit reached')
        return u