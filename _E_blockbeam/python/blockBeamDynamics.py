import numpy as np
import blockbeamParam as P


class blockBeamDynamics:
    def __init__(self, alpha):
        # Initial state conditions
        self.state = np.array([
            [P.z0],  # z initial position
            [P.theta0],  # Theta initial orientation
            [P.zdot0],  # zdot initial velocity
            [P.thetadot0],  # Thetadot initial velocity
        ])
        self.m1 = P.m1 * (1. + alpha * (2. * np.random.rand() - alpha))
        self.m2 = P.m2 * (1. + alpha * (2. * np.random.rand() - alpha))
        self.length = P.length * (1. + alpha * (2. * np.random.rand() - alpha))
        # the gravity constant is well known, so we don't change it.
        self.g = P.g
        # sample rate at which the dynamics are propagated
        self.Ts = P.Ts


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
        z = state[0][0]
        theta = state[1][0]
        zdot = state[2][0]
        thetadot = state[3][0]

        zddot = (z*thetadot**2 - self.g*np.sin(theta))
        thetaddot = (u*self.length*np.cos(theta)-2*self.m1*z*zdot*thetadot
                     -self.m1*self.g*z*np.cos(theta)
                     -self.m2*self.g*self.length*np.cos(theta)/2) \
                    / ((self.m2*self.length**2)/3 + self.m1*z**2)

        xdot = np.array([[zdot], [thetadot], [zddot], [thetaddot]])
        return xdot

    def h(self):
        # return the output equations
        # could also use input u if needed
        z = self.state[0][0]
        theta = self.state.item(1)
        y = np.array([[z], [theta]])
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
    return u